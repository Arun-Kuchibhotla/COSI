"""Four confidence-interval constructors from already computed estimates."""
import math
import numpy as np
from scipy.special import ndtri
from cosi_calibration import EPS, _control, _real_scalar

# Validate and return a nonempty vector of finite batch estimates.
def _validate_estimators(W):
    raw = np.asarray(W)
    if raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf":
        raise ValueError("W must be a nonempty numeric vector of finite batch estimates")
    W = raw.astype(float)
    if not np.all(np.isfinite(W)):
        raise ValueError("W must contain finite batch estimates")
    return W


# Evaluate and validate h(0),...,h(B), or compute the fixed normal scores.
def _score_values(B, method="symmetry", h=None):
    if method == "normality":
        values = ndtri(.5+np.arange(B+1)/(2*(B+1)))
    elif h is None:
        values = np.arange(B+1,dtype=float)
    elif callable(h):
        values = [h(i) for i in range(B+1)]
        if not all(_real_scalar(v) and math.isfinite(v) for v in values):
            raise ValueError("h must return a finite real scalar at each integer 0,...,B")
    else:
        values = h
    raw = np.asarray(values)
    if raw.shape != (B+1,) or raw.dtype.kind not in "iuf":
        raise ValueError("h must give B+1 numeric scores")
    values = raw.astype(float)
    if (not np.all(np.isfinite(values)) or np.any(values < 0)
            or np.any(np.diff(values) < 0) or not np.isfinite(values.sum())):
        raise ValueError("h must give finite, nonnegative, nondecreasing scores")
    if values[0] != 0:
        raise ValueError("The stated Walsh weights and N_h require h(0)=0")
    return values


# Validate the cutoff for binomial order statistics or signed-rank inversion.
def _check_cutoff(c, B, N=None):
    if not _real_scalar(c) or not math.isfinite(c):
        raise ValueError("c must be one finite numeric cutoff")
    if N is None:
        if c != int(c) or not -1 <= c <= B:
            raise ValueError("Binomial c must be an integer in [-1,B]")
    elif not (c == -1 or 0 <= c <= N/2):
        raise ValueError("Signed-rank c must be -1 or lie in [0,N_h/2]; it denotes kappa")


# Return supplied tie-breakers, or draw nonzero U(-1,1) values using a given RNG.
def _tie_breakers(B, U=None, rng=None):
    if U is None:
        uniform = np.random.uniform if rng is None else rng.uniform
        U = uniform(-1,1,B)
        while np.any(U == 0):
            U[U == 0] = uniform(-1,1,np.count_nonzero(U == 0))
    U = np.asarray(U,dtype=float)
    if U.shape != (B,) or not np.all(np.isfinite(U)) or np.any((U < -1)|(U > 1)|(U == 0)):
        raise ValueError("U must contain B nonzero values in [-1,1]")
    return U.copy()


# Construct [W_(B-c),W_(c+1)] from supplied estimates and binomial c* at alpha/2.
def cosi_median(W, c):
    W = _validate_estimators(W)
    _check_cutoff(c,len(W))
    return _interval(_order_value(W,len(W)-c),_order_value(W,c+1))


# Construct a weighted signed-rank interval. c denotes kappa_{B,alpha}^h,
# not a binomial cutoff; h may be a function or the B+1 scores h(0),...,h(B).
def cosi_symmetry(W, c, h=None, U=None, control=None):
    W = _validate_estimators(W)
    B = len(W)
    control = _control(control)
    scores = _score_values(B,"symmetry",h)
    _check_cutoff(c,B,math.fsum(scores))
    # With sign scores only diagonal Walsh averages have positive weight.
    # Use the closed median convention, including ties; no U or Walsh pairs
    # are needed. Calibrated kappa is integer; floor also defines other cutoffs.
    if np.all(scores[1:] == 1):
        binomial_c = B if c == -1 else B-math.floor(c)-1
        return cosi_median(W,binomial_c)
    if B*(B+1)//2 > control["max_walsh"]:
        raise ValueError("Number of Walsh pairs exceeds max_walsh")
    return _walsh_interval(W,_tie_breakers(B,U),scores,c,control)


# Construct the signed-rank interval with fixed normal scores and supplied kappa.
def cosi_normality(W, c, U=None, control=None):
    W = _validate_estimators(W)
    return cosi_symmetry(W,c,h=_score_values(len(W),"normality"),U=U,control=control)


# Construct the unimodal-star interval using independent anchor a and c* at alpha.
def cosi_unimodal_star(W, a, c):
    W = _validate_estimators(W)
    if not _real_scalar(a) or not math.isfinite(a):
        raise ValueError("a must be one finite anchor estimate")
    _check_cutoff(c,len(W))
    distance = (1+math.sqrt(2))*np.abs(W-a)
    lower,upper = np.minimum(W,a)-distance,np.maximum(W,a)+distance
    result = _interval(_order_value(lower,len(W)-c),_order_value(upper,c+1))
    result.update(L=lower,U=upper)
    return result

# Represent real endpoints, endpoint inclusion, and whether the resulting interval is empty.
def _interval(lower, upper, lower_closed=True, upper_closed=True):
    lower, upper = float(lower), float(upper)
    if math.isnan(lower) or math.isnan(upper):
        raise ArithmeticError("Undefined interval endpoint")
    closed = dict(lower=bool(math.isfinite(lower) and lower_closed),
                  upper=bool(math.isfinite(upper) and upper_closed))
    return dict(CI=dict(lower=lower, upper=upper), closed=closed,
                empty=bool(lower > upper or (lower == upper and not all(closed.values()))))


# Return a one-based order statistic with extended values outside the observed range.
def _order_value(x, k):
    if k <= 0:
        return -math.inf
    if k > len(x):
        return math.inf
    return float(np.partition(x, int(k)-1)[int(k)-1])


# Invert weighted augmented Walsh quantiles, retaining lexicographic endpoint inclusion.
def _walsh_interval(W, U, h_values, kappa, control):
    W, U, h_values = map(lambda x: np.asarray(x,dtype=float), (W,U,h_values))
    B, N = len(W), math.fsum(h_values[1:])
    order = np.lexsort((U,W))
    W, U = W[order], U[order]
    count = B*(B+1)//2
    av, au, weights = (np.empty(count) for _ in range(3))
    increments = np.diff(h_values)
    at = 0
    for i in range(B):
        length = B-i
        target = slice(at, at+length)
        av[target] = W[i]/2+W[i:]/2
        au[target] = U[i]/2+U[i:]/2
        weights[target] = increments[:length]
        at += length
    order = np.lexsort((au,av))
    av, au, weights = av[order], au[order], weights[order]
    # Accumulating in extended precision approaches R's cumsum behavior.
    cumulative = np.cumsum(weights, dtype=np.longdouble).astype(float)
    tol = 0.0 if np.all(h_values == np.floor(h_values)) and N <= 2**52 else control["mass_tol"]*max(1,N)
    if abs(cumulative[-1]-N) > max(tol, 8*EPS*max(1,N)):
        raise ArithmeticError("Total Walsh weight differs from sum h(i); check h(0)=0")
    cumulative[-1] = N

    # Find a weighted lexicographic quantile, optionally its right-hand limit.
    def quantile(gamma, right=False):
        if (not right and gamma <= 0) or (right and gamma < 0):
            return (-math.inf, 0.0)
        if gamma > N or (right and gamma >= N):
            return (math.inf, 0.0)
        if right and gamma == 0:
            indices = np.flatnonzero(cumulative > 0)
        elif right:
            indices = np.flatnonzero(cumulative > gamma+tol)
        else:
            indices = np.flatnonzero((cumulative > 0) & (cumulative >= gamma-tol))
        if not indices.size:
            return (math.inf, 0.0)
        i = indices[0]
        return (float(av[i]), float(au[i]))

    lower, upper = quantile(kappa, True), quantile(N-kappa)
    result = _interval(lower[0], upper[0], lower[1] < 0, upper[1] > 0)
    result.update(augmented_endpoints=np.asarray([lower,upper]), mass_tolerance=tol)
    return result


# Test real-valued interval membership, respecting open endpoints and empty intervals.
def contains_cosi_ci(object, t):
    """Membership including open endpoints, infinite endpoints and empty sets."""
    if not isinstance(object,dict) or not {"CI","closed","empty"} <= object.keys():
        raise TypeError("object must be a cosi_ci result")
    values = np.asarray(t,dtype=float)
    if object["empty"]:
        return np.zeros(values.shape,dtype=bool)
    lo, hi = object["CI"]["lower"], object["CI"]["upper"]
    lower = values >= lo if object["closed"]["lower"] else values > lo
    upper = values <= hi if object["closed"]["upper"] else values < hi
    return lower & upper & np.isfinite(values)
