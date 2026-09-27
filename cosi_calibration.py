"""Critical-value calibration; no data splitting or interval construction."""
import math
import random
from contextlib import contextmanager, nullcontext
from functools import reduce
from numbers import Real
import numpy as np
from scipy.optimize import brentq
from scipy.special import log_ndtr, ndtri_exp
from scipy.stats import binom
from cosi_smallest_B import _cutoff_probabilities, _logsf_series

EPS = np.finfo(float).eps


# Evaluate the two-sided binomial rejection probability without tail underflow.
def _sign_rejection_probability(B, k):
    if k < 0:
        return 0.0
    if k >= B/2:
        return 1.0
    if k == 0:
        return math.ldexp(1.0,1-B)
    tail = float(binom.cdf(k,B,.5))
    return min(1.0,2*tail if tail >= 1e-300 else
               math.exp(math.log(2)+_logsf_series(B,B-int(k)-1)))


# Compute the sign-score cutoff law by transforming the binomial law through
# kappa=B-c_star-1; retain the binomial law to couple draws with the median method.
def _sign_calibration(B, alpha):
    q, p_lower, p_upper = _cutoff_probabilities(B,float(alpha))
    bin_support = np.asarray([q-1,q],dtype=int)
    bin_prob = np.asarray([p_lower,p_upper])
    keep = bin_prob > 0
    bin_law = dict(support=bin_support[keep],probabilities=bin_prob[keep])
    cm, cp = B-q-1, B-q
    pm = _sign_rejection_probability(B,cm)
    pp = _sign_rejection_probability(B,cp)
    probabilities = bin_prob[::-1].copy()
    keep = probabilities > 0
    support = np.asarray([cm,cp])[keep]
    rejection = np.asarray([pm,pp])[keep]
    level_exact = bool(cp != B/2 or p_lower == 0)
    tau = p_lower
    # At an exact CDF level the singleton is c_plus, not c_minus.
    if p_lower == 0:
        cp,cm,pp = cm,cm-1,pm
        pm,tau = _sign_rejection_probability(B,cm),1.0
    law = dict(support=support,probabilities=probabilities[keep],
               rejection_probabilities=rejection,
               c_minus=cm,c_plus=cp,P_minus=pm,P_plus=pp,tau=tau,
               engine="sign-binomial",exact=True,median_compatible=True,
               nominal_level_exact=level_exact,
               rejection_probability_bound=float(np.dot(rejection,probabilities[keep])))
    return law,bin_law


# Calibrate and draw kappa for h=[h(1),...,h(B)], excluding h(0)=0.
# Return the chosen kappa and its law. A seed is local; a supplied Generator
# lets cosi_ci use its existing stream. u optionally supplies the cutoff uniform.
# Unit scores default to median-compatible randomization, conservative at some
# central cutoffs. Set median_compatible=False for the canonical P_h mixture.
def cosi_kappa(h, alpha=0.05, seed=None, *, randomized=True, control=None, u=None,
               median_compatible=True):
    scope = nullcontext(seed) if isinstance(seed,np.random.Generator) else _seed_scope(seed)
    with scope as rng:
        raw = np.asarray(h)
        if raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf":
            raise ValueError("h must be a nonempty numeric vector h(1),...,h(B), without h(0)")
        scores = raw.astype(float)
        if (not np.all(np.isfinite(scores)) or np.any(scores < 0) or
                np.any(np.diff(scores) < 0) or not np.isfinite(scores.sum())):
            raise ValueError("h must be finite, nonnegative, and nondecreasing with finite sum")
        if not _real_scalar(alpha) or not math.isfinite(alpha) or not 0 < alpha < 1:
            raise ValueError("alpha must lie strictly between 0 and 1")
        if not isinstance(randomized,(bool,np.bool_)):
            raise ValueError("randomized must be True or False")
        if not isinstance(median_compatible,(bool,np.bool_)):
            raise ValueError("median_compatible must be True or False")
        if u is not None and (not _real_scalar(u) or not math.isfinite(u) or not 0 <= u < 1):
            raise ValueError("u must lie in [0,1)")
        control = _control(control)
        if median_compatible and np.all(scores == 1):
            law,bin_law = _sign_calibration(len(scores),alpha)
            c_star = _draw_law(bin_law,randomized,"max",rng,u)
            kappa = len(scores)-c_star-1
        else:
            law = _sr_calibration(scores,alpha,control)
            kappa = _draw_law(law,randomized,"min",rng,u)
        return dict(kappa=kappa,distribution=law)

# Identify a real numeric scalar, excluding logical values.
def _real_scalar(value):
    return isinstance(value, Real) and not isinstance(value, (bool, np.bool_))


# Supply a local generator and seed callback randomness from numpy.random and
# random; restore both global RNG states even if validation or estimation fails.
@contextmanager
def _seed_scope(seed):
    if seed is None:
        yield np.random.default_rng()
        return
    if (not _real_scalar(seed) or not math.isfinite(seed) or seed != int(seed)
            or not 0 <= seed <= 2147483647):
        raise ValueError("seed must be None or an integer in [0,2147483647]")
    seed = int(seed)
    numpy_state, python_state = np.random.get_state(), random.getstate()
    try:
        np.random.seed(seed)
        random.seed(seed)
        yield np.random.default_rng(seed)
    finally:
        np.random.set_state(numpy_state)
        random.setstate(python_state)


# Fill in and validate resource limits and numerical tolerances.
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


# Draw a cutoff from its law, or select the requested conservative support value.
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


# Compute the two-point signed-rank cutoff law achieving the requested rejection level.
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


# Compute log(cosh(x)) stably for small and large arguments.
def _log_cosh(x):
    x = np.abs(np.asarray(x, dtype=float))
    out = x + np.log1p(np.exp(-2*x)) - math.log(2)
    small = x < 1
    out[small] = np.log1p(2*np.sinh(x[small]/2)**2)
    return out


# Compute a conservative signed-rank cutoff using the smallest of four valid tail radii.
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

    # Evaluate the derivative equation used to minimize the Chernoff radius.
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


# Convert a discrete score-sum distribution into the randomized signed-rank cutoff law.
def _law_from_pmf(values, probabilities, N, alpha, engine):
    cdf = np.cumsum(probabilities)
    k = int(np.searchsorted(cdf, alpha/2))
    if k >= len(values):
        raise ArithmeticError("Failed to locate a signed-rank quantile")
    cm = -1 if k == 0 else float(values[k-1])
    pm = 0 if k == 0 else min(1.0, 2*float(cdf[k-1]))
    return _kappa_law(cm, float(values[k]), pm, min(1.0, 2*float(cdf[k])), alpha, engine)


# Enumerate distinct subset sums and their multiplicities for an exact Bernoulli-sum law.
def _subset_sums(scores):
    sums = np.asarray([0.0])
    for s in scores:
        sums = np.concatenate((sums, sums+s))
    return np.unique(sums, return_counts=True)


# Select signed-rank cutoff support by counting sums of two half-distributions.
def _mitm_calibration(scores, N, alpha, control):
    n = len(scores)
    m = n//2
    if 2**((n+1)//2) > control["max_states"] or n > 52:
        return _concentration_calibration(scores, alpha, control,
                                          "meet-in-the-middle state/counting limit")
    a, ac = _subset_sums(scores[:m])
    b, bc = _subset_sums(scores[m:])
    cumulative_b = np.concatenate(([0], np.cumsum(bc)))

    # Count weighted pairs of half-sums below a threshold, optionally returning boundary positions.
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


# Compute an integer-weight Bernoulli-sum probability mass function by dynamic programming.
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


# Select an exact score-law algorithm when feasible, otherwise use concentration calibration.
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
