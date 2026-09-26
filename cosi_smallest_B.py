"""Find the smallest B for the randomized binomial cutoff in COSI.

Requires NumPy and SciPy.  The image's finite-sample bounds eliminate an
impossible interval; the remaining integers are checked in increasing order.
No monotonicity of the feasibility condition in B is assumed.

The randomization is independent of the binomial variable.  At an exact CDF
level it can degenerate to one integer.  Its largest support value is
q_B = min{k: P(Binom(B, 1/2) > k) <= alpha/2}.
The revised search also requires q_B < B.
smallest_B returns B, a randomized c_star, and its full distribution.
"""

import math
import numpy as np
from scipy.special import ndtri_exp
from scipy.stats import binom


def _check_alpha(alpha):
    alpha = float(alpha)
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be a finite number in (0, 1)")
    return alpha


def _exact_tail_leq(n, k, alpha):
    """Integer-arithmetic comparison, used only for small, near-equality cases."""
    # Symmetry: P(Y > k) = 2**(-n) * sum_{j=0}^{n-k-1} choose(n, j).
    coefficient, numerator = 1, 0
    for j in range(n - k):
        numerator += coefficient
        coefficient = coefficient * (n - j) // (j + 1)
    a_num, a_den = alpha.as_integer_ratio()
    return 2 * numerator * a_den <= a_num * (1 << n)


def _logsf_series(n, k):
    """Upper-tail log probability when SciPy's ordinary tail underflows.

    All calls have k >= n/2.  Starting at k+1, successive binomial masses
    decrease, so their ratios give a stable, positive series.
    """
    if k >= n:
        return -math.inf
    j = k + 1
    log_first = float(binom.logpmf(j, n, 0.5))
    term = total = 1.0
    eps = np.finfo(float).eps
    while j < n:
        ratio = (n - j) / (j + 1)
        # Later ratios are smaller: this bounds all omitted terms.
        if term * ratio / (1.0 - ratio) <= eps * total:
            break
        term *= ratio
        total += term
        j += 1
    return log_first + math.log(total)


def _tail_leq(n, k, alpha):
    """Vectorized test P(Binom(n, 1/2) > k) <= alpha/2."""
    n, k = np.broadcast_arrays(np.asarray(n, dtype=np.int64),
                               np.asarray(k, dtype=np.int64))
    shape = n.shape
    n, k = n.ravel(), k.ravel()
    tail = binom.sf(k, n, 0.5)  # Do not form 1 - binom.cdf(...).
    target = alpha / 2.0
    passed = tail <= target

    # Positive subnormal probabilities have poor relative precision.  Recompute
    # those comparisons in log space; this also handles alpha/2 underflow.
    if target < 1e-300:
        for i in np.flatnonzero((tail < 1e-300) & (k < n)):
            passed[i] = _logsf_series(int(n[i]), int(k[i])) <= (
                math.log(alpha) - math.log(2.0))

    # Inclusive equality matters.  Resolve small near-ties using exact integers
    # instead of changing the mathematical condition by adding a tolerance.
    near = np.abs(tail - target) <= 64 * np.finfo(float).eps * np.maximum(tail, target)
    for i in np.flatnonzero(near & (n <= 52)):
        passed[i] = _exact_tail_leq(int(n[i]), int(k[i]), alpha)
    return passed.reshape(shape)


def largest_support(B, alpha):
    """Largest support value of c*_{B, alpha/2}, including degenerate cases."""
    alpha = _check_alpha(alpha)
    if isinstance(B, bool) or not np.isscalar(B) or not math.isfinite(B) or B < 1 or B != int(B):
        raise ValueError("B must be a positive integer")
    if B > 2**52:
        raise ValueError("B exceeds the supported double-precision integer range")
    B = int(B)
    # Binary search is valid in the cutoff k for fixed B (not in B itself).
    lo, hi = (B + 1) // 2 - 1, B
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bool(_tail_leq(B, mid, alpha)):
            hi = mid
        else:
            lo = mid
    return hi


def _smallest_B_value(alpha, gamma, chunk_size=65536):
    """Find the smallest B with q_B/B - 1/2 <= gamma AND q_B < B.

    alpha : number in (0, 1); the cutoff uses alpha/2 internally.
    gamma : number in (0, 1/2].
    chunk_size : maximum number of candidate B's evaluated at once.

    q_B < B is equivalent to 2**(1-B) <= alpha.  For gamma < 1/2 it
    follows from the other condition; for gamma = 1/2 it determines B.

    Uses the finite-sample bounds in the supplied image, not a normal
    approximation to the final answer.  With alpha fixed, the number of
    candidates is O(1/gamma) as gamma tends to zero, instead of O(1/gamma**2).
    Probabilities and noninteger inputs use double precision.  Near probability
    ties with B > 52 are subject to the numerical accuracy of SciPy.
    """
    alpha = _check_alpha(alpha)
    gamma = float(gamma)
    if not math.isfinite(gamma) or not 0 < gamma <= 0.5:
        raise ValueError("gamma must be a finite number in (0, 1/2]")
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, (int, np.integer)) or chunk_size < 1:
        raise ValueError("chunk_size must be a positive integer")

    # q_B < B iff P(Binom(B, 1/2) > B-1) = 2**(-B) <= alpha/2.
    # Thus B >= ceil(log2(2/alpha)).  Compute this exactly for the input
    # float: alpha = m*2**e with 1/2 <= m < 1 implies that ceiling is 2-e.
    # This avoids log rounding at powers of two and alpha/2 underflow.
    first = 2 - math.frexp(alpha)[1]

    # Exact edge cases also avoid huge search bounds when gamma is tiny.
    if gamma == 0.5:
        return first
    if alpha >= 0.5:
        return 2

    log_tail = math.log(alpha) - math.log(2.0)
    z = -float(ndtri_exp(log_tail))  # Phi^{-1}(1-alpha/2), without cancellation.

    # Image: c <= z*sqrt(B)/2 + 1, and q_B = ceil(B/2) + c.
    # Thus q_B - B/2 <= z*sqrt(B)/2 + 3/2.
    root = (z + math.sqrt(z*z + 24*gamma)) / (4*gamma)
    if not math.isfinite(root) or root > math.sqrt(2**52 - 64):
        raise ValueError("Search bound exceeds the supported double-precision integer range")
    upper_real = root * root
    pad = 32 * np.finfo(float).eps * max(1.0, upper_real)
    last = max(first + 1, math.ceil(upper_real + pad))
    intervals = [(first, last)]

    # Image's lower bound implies q_B-B/2 >= z*sqrt(B)/2 - D.
    # Hence feasibility requires 2*gamma*t^2 - z*t + 2*D >= 0, t=sqrt(B).
    A = max(z**3 / 5, math.sqrt(1 + (2*math.log(2)-1) / first**2))
    D = 1.5 + A / (2*math.sqrt(first))
    discriminant = z*z - 16*gamma*D
    if discriminant > 0:
        s = math.sqrt(discriminant)
        # Stable quadratic roots: avoid subtracting two nearly equal numbers.
        left = (4*D / (z+s))**2
        right = ((z+s) / (4*gamma))**2
        # Keep both boundary neighborhoods to protect against root roundoff.
        end1 = min(last, math.ceil(left + pad))
        start2 = max(first, end1 + 1, math.floor(right - pad))
        intervals = [(first, end1), (start2, last)]

    for start, stop in intervals:
        for low in range(start, stop + 1, int(chunk_size)):
            B = np.arange(low, min(stop + 1, low + chunk_size), dtype=np.int64)
            # q_B <= floor(B/2 + gamma*B).  This avoids subtracting q_B/B - 1/2.
            # Separate the parity term so adding 1/2 does not erase tiny gamma.
            k = B // 2 + np.floor(gamma*B + 0.5*(B % 2)).astype(np.int64)
            # Both requirements together are q_B <= min(k, B-1).
            k = np.minimum(k, B - 1)
            hits = np.flatnonzero(_tail_leq(B, k, alpha))
            if hits.size:
                return int(B[hits[0]])
    raise ArithmeticError("No solution within the image's upper bound; check numerical accuracy")


def _cutoff_probabilities(B, alpha):
    """Return q_B and the probabilities at q_B-1 and q_B.

    If S(k)=P(Binom(B,1/2)>k) and m=P(Binom(B,1/2)=q_B), then
    p_lower=(alpha/2-S(q_B))/m and p_upper=(S(q_B-1)-alpha/2)/m.
    Hence the randomized upper-tail probability is exactly alpha/2 in
    exact arithmetic.  No data enter this randomization.
    """
    first = 2 - math.frexp(alpha)[1]
    if B == first:
        # Here q_B=B-1.  This formula stays accurate even for subnormal alpha
        # and is exactly zero at a degenerate power-of-two boundary.
        q = B - 1
        p_lower = (math.ldexp(alpha, B - 1) - 1.0) / B
        p_upper = 1.0 - p_lower
    else:
        q = largest_support(B, alpha)
        if B <= 52:
            # Exact integer numerators handle small degenerate CDF levels.
            tail_num = sum(math.comb(B, j) for j in range(B - q))
            a_num, a_den = alpha.as_integer_ratio()
            denominator = 2 * a_den * math.comb(B, q)
            numerator = a_num * (1 << B) - 2 * a_den * tail_num
            p_lower = numerator / denominator
            p_upper = (denominator - numerator) / denominator
        else:
            tails = binom.sf([q, q - 1], B, 0.5)
            mass = float(binom.pmf(q, B, 0.5))
            target = alpha / 2.0
            if target >= 1e-300 and mass >= 1e-300:
                p_lower = (target - tails[0]) / mass
                p_upper = (tails[1] - target) / mass
            else:
                # Log differences with expm1 avoid underflow and cancellation.
                logs = [math.log(t) if t >= 1e-300 else _logsf_series(B, k)
                        for t, k in zip(tails, [q, q - 1])]
                log_target = math.log(alpha) - math.log(2.0)
                log_mass = float(binom.logpmf(q, B, 0.5))
                p_lower = math.exp(log_target - log_mass) * (-math.expm1(logs[0] - log_target))
                p_upper = math.exp(logs[1] - log_mass) * (-math.expm1(log_target - logs[1]))

    # Roundoff can move an endpoint slightly outside [0,1].  Reject material
    # violations; compute the larger probability as a complement of the
    # smaller one to preserve accuracy and make the probabilities sum to one.
    if not all(math.isfinite(p) and -1e-10 <= p <= 1 + 1e-10
               for p in (p_lower, p_upper)):
        raise ArithmeticError("Could not resolve cutoff probabilities in double precision")
    p_lower, p_upper = (float(min(1.0, max(0.0, p))) for p in (p_lower, p_upper))
    if p_lower <= p_upper:
        p_upper = 1.0 - p_lower
    else:
        p_lower = 1.0 - p_upper
    if p_upper <= 0:
        raise ArithmeticError("Upper support probability is numerically unresolved")
    return q, p_lower, p_upper


def smallest_B(alpha, gamma, chunk_size=65536, *, seed=None):
    """Return the smallest feasible B together with c*_{B, alpha/2}.

    Parameters
    ----------
    alpha, gamma, chunk_size : as in the original search.
        The two requirements concern q_B, the largest support value:
        q_B/B - 1/2 <= gamma and q_B < B.
    seed : optional nonnegative integer or numpy.random.Generator.
        Controls the independent randomization only; B is deterministic.
        With seed=None, NumPy creates a fresh generator.  Use a Generator
        for successive draws from one reproducible random stream.

    Returns
    -------
    dict with B, c_star, largest_support, support, probabilities.
        c_star is one draw, not the largest support value.  The probabilities
        correspond to the support entries in order; zero-probability support
        points are omitted.  To obtain only the integer B, use result['B'].
    """
    alpha = _check_alpha(alpha)
    rng = np.random.default_rng(seed)
    B = _smallest_B_value(alpha, gamma, chunk_size)
    q, p_lower, p_upper = _cutoff_probabilities(B, alpha)
    if p_lower == 0:
        c_star, support, probabilities = q, (q,), (1.0,)
    else:
        # Sample the less likely outcome directly for numerical stability.
        c_star = (q - int(rng.binomial(1, p_lower)) if p_lower <= p_upper
                  else q - 1 + int(rng.binomial(1, p_upper)))
        support, probabilities = (q - 1, q), (p_lower, p_upper)
    return {"B": B, "c_star": c_star, "largest_support": q,
            "support": support, "probabilities": probabilities}


if __name__ == "__main__":
    for gamma in (1/2, 1/3, 1/4, 1/5, 1/10, 1/20, 1/100):
        result = smallest_B(0.05, gamma, seed=42)
        print(f"alpha=0.05, gamma={gamma:.8g}: {result}")
