"""Find the smallest B for the randomized binomial cutoff in COSI.

Requires NumPy and SciPy.  The image's finite-sample bounds eliminate an
impossible interval; the remaining integers are checked in increasing order.
No monotonicity of the feasibility condition in B is assumed.

The randomization is independent of the binomial variable.  At an exact CDF
level it can degenerate to one integer.  Its largest support value is
q_B = min{k: P(Binom(B, 1/2) > k) <= alpha/2}.
The revised search also requires q_B < B.
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


def smallest_B(alpha, gamma, chunk_size=65536):
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


if __name__ == "__main__":
    for gamma in (1/2, 1/3, 1/4, 1/5, 1/10, 1/20, 1/100):
        B = smallest_B(0.05, gamma)
        print(f"alpha=0.05, gamma={gamma:.8g}: B={B}, q_B={largest_support(B, 0.05)}")
