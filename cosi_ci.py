"""COSI confidence intervals matching the accompanying base-R implementation.

Requires NumPy and SciPy. See cosi_ci_README.md for assumptions and examples.
Use a shared ``randomness`` dictionary for R/Python comparisons: equal native
R and NumPy seeds do not imply equal draws. Row identifiers in returned
batch_rows and in shared permutations are ONE-based, matching R.
"""
from __future__ import annotations

import math
from functools import reduce
from numbers import Real

import numpy as np
from scipy.optimize import brentq
from scipy.special import log_ndtr, ndtri, ndtri_exp

from cosi_smallest_B import smallest_B_distribution, binomial_cutoff_distribution

EPS = np.finfo(float).eps


def _real_scalar(value):
    return isinstance(value, Real) and not isinstance(value, (bool, np.bool_))


def _control(control=None):
    defaults = dict(max_states=2**20, max_walsh=2_000_000, max_dp_work=50_000_000,
                    mass_tol=64*EPS, max_iterations=256, concentration_B=40)
    if control is None:
        control = {}
    if not isinstance(control, dict) or set(control)-set(defaults):
        raise ValueError("control must use the documented option names")
    defaults.update(control)
    for name, x in defaults.items():
        if (not _real_scalar(x) or not math.isfinite(x) or x < 0
                or (name != "mass_tol" and (x < 1 or x != math.floor(x)))):
            raise ValueError(f"Invalid control[{name!r}]")
        if name != "mass_tol":
            defaults[name] = int(x)
    return defaults


def _draw_law(law, randomized, conservative, rng, u=None):
    support = np.asarray(law["support"])
    if not randomized:
        return support.max().item() if conservative == "max" else support.min().item()
    if support.size == 1:
        return support[0].item()
    if u is None:
        u = rng.random()
    cumulative = np.cumsum(law["probabilities"])
    cumulative[-1] = 1.0
    return support[np.searchsorted(cumulative, u, side="right")].item()


def make_cosi_randomness(n, B, seed=None):
    """Make a data-independent plan to reuse in Python and R.

    Convert arrays to lists before JSON export. A plan has a one-based row
    permutation, a uniform cutoff draw, and B U(-1,1) tie-breakers.
    ``seed`` accepts the usual NumPy default_rng arguments, including a Generator.
    """
    for x in (n, B):
        if not _real_scalar(x) or not math.isfinite(x) or x < 1 or x != int(x):
            raise ValueError("n and B must be positive integers")
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(int(n)) + 1
    cutoff_uniform = float(rng.random())
    uniforms = rng.uniform(-1, 1, int(B))
    while np.any(uniforms == 0):
        uniforms[uniforms == 0] = rng.uniform(-1, 1, np.count_nonzero(uniforms == 0))
    return dict(permutation=permutation, cutoff_uniform=cutoff_uniform, uniforms=uniforms)


def _check_randomness(randomness, n, B, signed_rank):
    if randomness is None:
        return
    if not isinstance(randomness, dict) or set(randomness) != {
            "permutation", "cutoff_uniform", "uniforms"}:
        raise ValueError("randomness needs permutation, cutoff_uniform, uniforms")
    p = np.asarray(randomness["permutation"])
    if p.ndim != 1 or p.size != n or not np.array_equal(np.sort(p), np.arange(1, n+1)):
        raise ValueError("randomness['permutation'] must permute 1,...,n")
    u = randomness["cutoff_uniform"]
    if not _real_scalar(u) or not math.isfinite(u) or not 0 <= u < 1:
        raise ValueError("randomness['cutoff_uniform'] must lie in [0,1)")
    if signed_rank:
        u = np.asarray(randomness["uniforms"], dtype=float)
        if (u.shape != (B,) or not np.all(np.isfinite(u))
                or np.any((u < -1) | (u > 1) | (u == 0))):
            raise ValueError("randomness['uniforms'] needs B nonzero values in [-1,1]")


def _kappa_law(cm, cp, pm, pp, alpha, engine):
    # Larger kappa rejects more, hence the conservative endpoint is the MINIMUM.
    if pm >= alpha or pp < alpha or pp <= pm:
        raise ArithmeticError("Could not resolve the signed-rank critical-value bracket")
    tau = (alpha-pm)/(pp-pm)
    probabilities = np.asarray([1-tau, tau])
    keep = probabilities > 0
    return dict(support=np.asarray([cm, cp])[keep], probabilities=probabilities[keep],
                rejection_probabilities=np.asarray([pm, pp])[keep],
                c_minus=cm, c_plus=cp, P_minus=pm, P_plus=pp, tau=tau,
                engine=engine, exact=True)


def _log_cosh(x):
    x = np.abs(np.asarray(x, dtype=float))
    out = x + np.log1p(np.exp(-2*x)) - math.log(2)
    small = x < 1
    out[small] = np.log1p(2*np.sinh(x[small]/2)**2)
    return out


def _concentration_calibration(scores, alpha, control, reason="general real scores"):
    """Tightest computed radius from four finite-sample-valid inequalities.

    Work with T=(G-N/2)/(N/2), whose coefficients sum to one. Every finite
    Chernoff parameter is valid, even if numerical optimization does not finish.
    """
    scores = np.asarray(scores, dtype=float)
    N = math.fsum(scores)
    if N == 0:
        return _kappa_law(-1, 0, 0, 1, alpha, "zero-scores")
    mu = N/2
    w = scores[scores > 0]/N
    variance = math.fsum(w*w)
    sd = math.sqrt(variance)
    log_alpha = math.log(alpha)
    L = math.log(2)-log_alpha
    # Bentkus--Dzindzalieta (2015), Theorem 1.1, optimal constant-factor
    # Gaussian comparison. This is an inequality, not a normal approximation.
    log_C = -math.log(4)-float(log_ndtr(-math.sqrt(2)))
    r_hoeffding = sd*math.sqrt(2*L)
    r_gaussian = -sd*float(ndtri_exp(log_alpha-math.log(2)-log_C))
    r_chebyshev = sd*math.exp(-log_alpha/2)
    # Mathematical tie implied by C_*; avoid a library-rounding tie-break.
    if alpha == .5:
        r_gaussian = r_chebyshev = sd*math.sqrt(2)

    def slope(lam):
        z = lam*w
        part = np.empty_like(z)
        small = z < 20
        part[small] = z[small]*np.tanh(z[small])-_log_cosh(z[small])
        zz = z[~small]
        e = np.exp(-2*zz)
        part[~small] = math.log(2)-np.log1p(e)-2*zz*e/(1+e)
        return math.fsum(part)-L

    lambda_H = math.sqrt(2*L/variance)
    candidates = [lambda_H]
    hi = lambda_H
    f_hi = slope(hi)
    for _ in range(control["max_iterations"]):
        if f_hi >= 0 or hi > np.finfo(float).max/2:
            break
        hi *= 2
        f_hi = slope(hi)
    candidates.append(hi)
    if math.isfinite(f_hi) and f_hi >= 0:
        try:
            root = brentq(slope, 0, hi, xtol=math.sqrt(EPS),
                          maxiter=control["max_iterations"])
            if root > 0 and math.isfinite(root):
                candidates.append(root)
        except (ValueError, RuntimeError):
            pass  # Keep the best finite parameter already evaluated.
    candidate_radii = [L/lam + math.fsum(_log_cosh(lam*w)/lam) for lam in candidates]
    chosen = int(np.argmin(candidate_radii))
    lam = candidates[chosen]
    K = math.fsum(_log_cosh(lam*w))
    radii = dict(chernoff=candidate_radii[chosen], gaussian=r_gaussian,
                 hoeffding=r_hoeffding, chebyshev=r_chebyshev)
    best = min(radii, key=radii.get)
    cushion = (control["mass_tol"]+256*EPS)*max(1, N)
    raw_kappa = mu*(1-radii[best])
    kappa = raw_kappa-cushion if math.isfinite(raw_kappa) and raw_kappa > cushion else -1.0
    log_bound = -math.inf
    if kappa >= 0:
        effective_kappa = kappa+control["mass_tol"]*max(1, N)
        r = (mu-effective_kappa)/mu
        bounds = [math.log(2)+K-lam*r,
                  math.log(2)+log_C+float(log_ndtr(-r/sd)),
                  math.log(2)-r*r/(2*variance), math.log(variance)-2*math.log(r)]
        log_bound = min(0, *bounds)
        if math.isnan(log_bound) or log_bound > log_alpha:
            kappa, log_bound = -1.0, -math.inf
    result = dict(support=np.asarray([kappa]), probabilities=np.asarray([1.0]),
                rejection_probabilities=np.asarray([np.nan]), engine="concentration",
                exact=False, bound_used=best, rejection_probability_bound=math.exp(log_bound),
                log_rejection_probability_bound=log_bound, radius=mu*radii[best],
                radius_fraction=radii[best], candidate_radii={k: mu*v for k,v in radii.items()},
                numerical_cushion=cushion, fallback_reason=reason)
    result["lambda"] = lam
    return result


def _law_from_pmf(values, probabilities, N, alpha, engine):
    cdf = np.cumsum(probabilities)
    k = int(np.searchsorted(cdf, alpha/2))
    if k >= len(values):
        raise ArithmeticError("Failed to locate a signed-rank quantile")
    cm = -1 if k == 0 else float(values[k-1])
    pm = 0 if k == 0 else min(1.0, 2*float(cdf[k-1]))
    return _kappa_law(cm, float(values[k]), pm, min(1.0, 2*float(cdf[k])), alpha, engine)


def _subset_sums(scores):
    sums = np.asarray([0.0])
    for s in scores:
        sums = np.concatenate((sums, sums+s))
    return np.unique(sums, return_counts=True)


def _mitm_calibration(scores, N, alpha, control):
    n = len(scores)
    m = n//2
    if 2**((n+1)//2) > control["max_states"] or n > 52:
        return _concentration_calibration(scores, alpha, control,
                                          "meet-in-the-middle state/counting limit")
    a, ac = _subset_sums(scores[:m])
    b, bc = _subset_sums(scores[m:])
    cumulative_b = np.concatenate(([0], np.cumsum(bc)))

    def count_pairs(x, strict=False, positions=False):
        j = np.searchsorted(b, x-a, side="right")
        # Addition-based corrections keep support-boundary comparisons
        # consistent despite subtraction rounding in x-a.
        while True:
            ii = np.flatnonzero(j > 0)
            sums = a[ii]+b[j[ii]-1]
            ii = ii[sums >= x if strict else sums > x]
            if not ii.size:
                break
            j[ii] -= 1
        while True:
            ii = np.flatnonzero(j < len(b))
            sums = a[ii]+b[j[ii]]
            ii = ii[sums < x if strict else sums <= x]
            if not ii.size:
                break
            j[ii] += 1
        count = int(np.sum(ac*cumulative_b[j]))
        return (count, j) if positions else count

    total = 2**n
    rank = max(1, math.ceil(alpha*2**(n-1)))
    if count_pairs(0) >= rank:
        cp = 0.0
    else:
        lo, hi = 0.0, float(a[-1]+b[-1])
        converged = False
        for _ in range(control["max_iterations"]):
            mid = lo+(hi-lo)/2
            if mid == lo or mid == hi:
                converged = True
                break
            count, j = count_pairs(mid, positions=True)
            if count < rank:
                lo = mid
            else:
                ii = np.flatnonzero(j > 0)
                hi = float(np.max(a[ii]+b[j[ii]-1]))
        if not converged:
            return _concentration_calibration(scores, alpha, control,
                                              "pair-sum selection iteration limit")
        _, j = count_pairs(lo, positions=True)
        ii = np.flatnonzero(j < len(b))
        cp = float(np.min(a[ii]+b[j[ii]]))
    count, j = count_pairs(cp, strict=True, positions=True)
    ii = np.flatnonzero(j > 0)
    cm = float(np.max(a[ii]+b[j[ii]-1])) if ii.size else -1
    return _kappa_law(cm, cp, 2*count/total, min(1, 2*count_pairs(cp)/total),
                      alpha, "meet-in-the-middle")


def _lattice_pmf(weights):
    # Convolution with [1/2 at 0, 1/2 at weight], using vector slices.
    pmf = np.asarray([1.0])
    for w in weights:
        w = int(w)
        old = pmf/2
        pmf = np.zeros(len(old)+w)
        pmf[:len(old)] = old
        pmf[w:w+len(old)] += old
    return pmf


def _sr_calibration(scores, alpha, control):
    scores = np.asarray(scores, dtype=float)
    N, B = math.fsum(scores), len(scores)
    positive = scores[scores > 0]
    if not positive.size:
        return _kappa_law(-1, 0, 0, 1, alpha, "zero-scores")
    p_zero = 2.0**(1-len(positive))
    if alpha <= p_zero:
        return _kappa_law(-1, 0, 0, p_zero, alpha, "zero-threshold")
    if alpha/2 == 0:
        return _concentration_calibration(scores, alpha, control,
                                          "discrete tail probability underflow")
    if B <= 1023 and scores[0] > 0 and np.all(scores == scores[0]*np.arange(1,B+1)):
        # SciPy has no public qsignrank/psignrank equivalents. Compute the
        # same Wilcoxon distribution by an exact lattice recurrence instead.
        pmf = _lattice_pmf(range(1, B+1))
        keep = np.flatnonzero(pmf > 0)
        return _law_from_pmf(keep*scores[0], pmf[keep], N, alpha, "signed-rank")
    step = float(positive.min())
    ratios = positive/step
    if np.all(positive == np.floor(positive)) and positive.max() <= 2**52:
        step = float(reduce(math.gcd, map(int,positive)))
        ratios = positive/step
    if np.all(ratios == np.floor(ratios)) and np.all(step*ratios == positive):
        states = math.fsum(ratios)+1
        if states <= control["max_states"] and len(ratios)*states <= control["max_dp_work"]:
            pmf = _lattice_pmf(ratios)
            keep = np.flatnonzero(pmf > 0)
            return _law_from_pmf(keep*step, pmf[keep], N, alpha, "lattice-dynamic-programming")
    if len(positive) <= min(14, math.log2(control["max_states"])):
        values, counts = _subset_sums(positive)
        return _law_from_pmf(values, counts/2**len(positive), N, alpha, "enumeration")
    if B >= control["concentration_B"]:
        return _concentration_calibration(scores, alpha, control, "B reached concentration_B")
    return _mitm_calibration(positive, N, alpha, control)


def _interval(lower, upper, lower_closed=True, upper_closed=True):
    lower, upper = float(lower), float(upper)
    if math.isnan(lower) or math.isnan(upper):
        raise ArithmeticError("Undefined interval endpoint")
    closed = dict(lower=bool(math.isfinite(lower) and lower_closed),
                  upper=bool(math.isfinite(upper) and upper_closed))
    return dict(CI=dict(lower=lower, upper=upper), closed=closed,
                empty=bool(lower > upper or (lower == upper and not all(closed.values()))))


def _order_value(x, k):
    if k <= 0:
        return -math.inf
    if k > len(x):
        return math.inf
    return float(np.partition(x, int(k)-1)[int(k)-1])


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


class CosiResult(dict):
    """Dictionary with the same named result fields as the R list."""
    def __repr__(self):
        text = f"COSI {self['method']} interval; B={self['B']}, batches={self['batches']}, alpha={self['alpha']:g}\n"
        law = self["kappa_distribution"]
        if law and law["engine"] == "concentration":
            text += (f"Conservative concentration calibration ({law['bound_used']}); "
                     f"rejection bound <= {law['rejection_probability_bound']:.6g}\n")
        if self["empty"]:
            return text+"empty set"
        left, right = ("[" if self["closed"]["lower"] else "("), ("]" if self["closed"]["upper"] else ")")
        return text+f"{left}{self['CI']['lower']:.7g}, {self['CI']['upper']:.7g}{right}"


def cosi_ci(data, estimate, alpha=.05, gamma=1/3, method="median", randomized=True,
            h=None, control=None, randomness=None, *, seed=None):
    """Construct a COSI interval from balanced, disjoint batch estimates.

    ``estimate`` receives a 2-D row subset and must return a finite real scalar.
    Methods: median, symmetry (default h(i)=i), normality, unimodal-star.
    ``randomized`` controls cutoff draws only; splits and tie-breakers remain
    random. ``randomness`` overrides all random inputs, enabling R/Python
    agreement. Otherwise ``seed`` controls a NumPy Generator. Do not provide
    both. Concentration cutoffs are deterministic for either randomized flag.
    """
    if not callable(estimate):
        raise TypeError("estimate must be a function")
    if not _real_scalar(alpha) or not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0,1)")
    if not _real_scalar(gamma) or not math.isfinite(gamma) or not 0 < gamma < .5:
        raise ValueError("gamma must lie strictly between 0 and 1/2")
    if not isinstance(randomized, (bool,np.bool_)):
        raise ValueError("randomized must be True or False")
    method = {"symmetric":"symmetry", "unimodal":"unimodal-star"}.get(method,method)
    if method not in ("median","symmetry","normality","unimodal-star"):
        raise ValueError("Unknown method")
    control = _control(control)
    data = np.asarray(data)
    if data.ndim == 1:
        data = data[:,None]
    if data.ndim != 2 or min(data.shape) == 0:
        raise ValueError("data must be a nonempty vector or matrix")
    n = data.shape[0]
    split_calibration = smallest_B_distribution(alpha,gamma)
    B = split_calibration["B"]
    batches = B+int(method == "unimodal-star")
    if batches > n:
        raise ValueError(f"Not enough samples for splitting: need {batches} batches; sample size is {n}")
    signed_rank = method in ("symmetry","normality")
    _check_randomness(randomness,n,B,signed_rank)
    if randomness is not None and seed is not None:
        raise ValueError("Supply randomness or seed, not both")
    rng = None if randomness is not None else np.random.default_rng(seed)
    cutoff_uniform = None if randomness is None else randomness["cutoff_uniform"]
    cutoff = cutoff_law = kappa = kappa_law = h_values = None
    if method == "median":
        cutoff_law = split_calibration
        cutoff = _draw_law(cutoff_law,randomized,"max",rng,cutoff_uniform)
    elif method == "unimodal-star":
        cutoff_law = binomial_cutoff_distribution(B,alpha)
        cutoff = _draw_law(cutoff_law,randomized,"max",rng,cutoff_uniform)
    else:
        if B*(B+1)//2 > control["max_walsh"]:
            raise ValueError(f"Need {B*(B+1)//2} Walsh pairs; max_walsh is {control['max_walsh']}")
        if method == "normality":
            h_values = ndtri(.5+np.arange(B+1)/(2*(B+1)))
        else:
            h = (lambda i: i) if h is None else h
            if not callable(h):
                raise TypeError("h must be a function for method='symmetry'")
            values = [h(i) for i in range(B+1)]
            if not all(_real_scalar(v) and math.isfinite(v) for v in values):
                raise ValueError("h must return a finite real scalar at each of 0,...,B")
            h_values = np.asarray(values,dtype=float)
        if h_values[0] != 0:
            raise ValueError("The stated Walsh weights and N_h require h(0)=0")
        if np.any(h_values < 0) or np.any(np.diff(h_values) < 0) or not np.isfinite(h_values.sum()):
            raise ValueError("h values must be nonnegative, nondecreasing, with finite sum")
        kappa_law = _sr_calibration(h_values[1:],alpha,control)
        kappa = _draw_law(kappa_law,randomized,"min",rng,cutoff_uniform)

    permutation = (rng.permutation(n)+1 if randomness is None else
                   np.asarray(randomness["permutation"],dtype=int))
    shuffled = data[permutation-1,:]
    # This matches R's split(1:n, sort((1:n) %% batches)), including which
    # labeled batches receive the additional rows when n is not divisible.
    labels = np.sort(np.arange(1,n+1) % batches)
    groups = [np.flatnonzero(labels == j) for j in range(batches)]
    all_estimates = np.empty(batches)
    for j, indices in enumerate(groups):
        try:
            value = estimate(shuffled[indices,:])
        except Exception as error:
            raise ValueError(f"estimate failed on batch {j+1}: {error}") from error
        if not _real_scalar(value) or not math.isfinite(value):
            raise ValueError(f"estimate must return a finite real scalar; failed on batch {j+1}")
        all_estimates[j] = value

    anchor = uniforms = None
    if method == "median":
        W = all_estimates
        interval = _interval(_order_value(W,B-cutoff),_order_value(W,cutoff+1))
    elif method == "unimodal-star":
        anchor = float(all_estimates[0])
        W = all_estimates[1:]
        distance = (1+math.sqrt(2))*np.abs(W-anchor)
        lower, upper = np.minimum(W,anchor)-distance, np.maximum(W,anchor)+distance
        interval = _interval(_order_value(lower,B-cutoff),_order_value(upper,cutoff+1))
        interval.update(L=lower,U=upper)
    else:
        W = all_estimates
        uniforms = (rng.uniform(-1,1,B) if randomness is None else
                    np.asarray(randomness["uniforms"],dtype=float).copy())
        while np.any(uniforms == 0):
            uniforms[uniforms == 0] = rng.uniform(-1,1,np.count_nonzero(uniforms == 0))
        interval = _walsh_interval(W,uniforms,h_values,kappa,control)
    interval.update(B=B,batches=batches,ci_est=W,all_estimates=all_estimates,
                    anchor=anchor,uniforms=uniforms,c=cutoff,cutoff_distribution=cutoff_law,
                    kappa=kappa,kappa_distribution=kappa_law,h_values=h_values,
                    N_h=None if h_values is None else math.fsum(h_values),
                    split_calibration=split_calibration,
                    batch_sizes={str(j):len(i) for j,i in enumerate(groups)},
                    batch_rows={str(j):permutation[i] for j,i in enumerate(groups)},alpha=float(alpha),
                    gamma=float(gamma),method=method,randomized=bool(randomized))
    return CosiResult(interval)


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
