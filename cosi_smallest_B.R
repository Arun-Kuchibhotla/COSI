# Find the smallest B for the randomized binomial cutoff in COSI.
# Uses only base R (the stats package).  Source this file, then call
# smallest_B(alpha = 0.05, gamma = 1/3).
#
# Independent randomization is assumed.  The largest support value is
# q_B = min{k: P(Binom(B, 1/2) > k) <= alpha/2}.
# The randomization may degenerate to a single integer at an exact CDF level.
# The revised search also requires q_B < B.

.cosi_check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a finite number in (0, 1)")
  }
}

.cosi_tail_leq <- function(B, k, alpha) {
  # R evaluates log upper tails directly.  Neither 1-alpha/2 nor 1-CDF is
  # formed, and alpha/2 need not be representable as a positive number.
  log_tail <- pbinom(k, size = B, prob = 0.5,
                     lower.tail = FALSE, log.p = TRUE)
  log_target <- log(alpha) - log(2)
  passed <- log_tail <= log_target

  # Resolve small, nearly equal probabilities without a permissive tolerance.
  # Pascal's triangle has exact integer entries and sums in double precision
  # for B <= 52, and division by 2^B is exact in binary arithmetic.
  near <- which(B <= 52 & is.finite(log_tail) &
                abs(log_tail - log_target) <=
                  64 * .Machine$double.eps * max(1, abs(log_target)))
  for (i in near) {
    n <- B[i]
    row <- 1
    for (j in seq_len(n)) {
      row <- c(row, 0) + c(0, row)
    }
    # By symmetry, P(Y > k) = P(Y <= n-k-1).
    m <- n - k[i]
    tail <- if (m <= 0) 0 else sum(row[seq_len(m)]) / 2^n
    passed[i] <- tail <= alpha / 2
  }
  passed
}

largest_support <- function(B, alpha) {
  # Largest support value of c*_{B, alpha/2}; binary search is valid in k.
  .cosi_check_alpha(alpha)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > 2^52) {
    stop("B must be a positive integer no larger than 2^52")
  }
  lo <- ceiling(B/2) - 1
  hi <- B
  while (hi - lo > 1) {
    mid <- floor(lo + (hi - lo)/2)
    if (.cosi_tail_leq(B, mid, alpha)) hi <- mid else lo <- mid
  }
  hi
}

smallest_B <- function(alpha, gamma, chunk_size = 65536) {
  # Find the smallest B with q_B/B - 1/2 <= gamma AND q_B < B.
  # alpha is the overall level: the alpha/2 adjustment is made internally.
  # chunk_size bounds memory consumption during the vectorized search.
  # q_B < B is equivalent to 2^(1-B) <= alpha.  For gamma < 1/2 it
  # follows from the other condition; for gamma = 1/2 it determines B.
  #
  # With alpha fixed, the image reduces the candidate count to O(1/gamma)
  # as gamma tends to zero, instead of O(1/gamma^2).
  # Near probability ties with B > 52 depend on pbinom's numerical accuracy.
  .cosi_check_alpha(alpha)
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
      gamma <= 0 || gamma > 0.5) {
    stop("gamma must be a finite number in (0, 1/2]")
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 ||
      chunk_size != floor(chunk_size)) {
    stop("chunk_size must be a positive integer")
  }

  # q_B < B iff P(Binom(B, 1/2) > B-1) = 2^(-B) <= alpha/2.
  # Start at ceiling(log2(2/alpha)), calculated without forming 2/alpha.
  first <- max(2, ceiling(1 - log2(alpha)))
  # Correct possible log rounding at powers of two using exact binary powers.
  # Comparing 2^(1-B) with alpha also avoids underflow of alpha/2.
  while (2^(1-first) > alpha) first <- first + 1
  while (first > 2 && 2^(2-first) <= alpha) first <- first - 1

  # Exact edge cases avoid unnecessary work (and huge bounds for tiny gamma).
  if (gamma == 0.5) return(first)
  if (alpha >= 0.5) return(2)

  log_tail <- log(alpha) - log(2)
  z <- qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)

  # Image: c <= z*sqrt(B)/2 + 1 and q_B = ceiling(B/2) + c.
  # Therefore q_B-B/2 <= z*sqrt(B)/2 + 3/2.
  root <- (z + sqrt(z*z + 24*gamma)) / (4*gamma)
  if (!is.finite(root) || root > sqrt(2^52 - 64)) {
    stop("Search bound exceeds the supported double-precision integer range")
  }
  upper_real <- root^2
  pad <- 32 * .Machine$double.eps * max(1, upper_real)
  last <- max(first + 1, ceiling(upper_real + pad))
  intervals <- list(c(first, last))

  # The image also gives q_B-B/2 >= z*sqrt(B)/2-D, with this constant D.
  # A feasible B must satisfy 2*gamma*t^2-z*t+2*D >= 0 for t=sqrt(B).
  A <- max(z^3/5, sqrt(1 + (2*log(2)-1)/first^2))
  D <- 1.5 + A/(2*sqrt(first))
  discriminant <- z*z - 16*gamma*D
  if (discriminant > 0) {
    s <- sqrt(discriminant)
    # Stable quadratic formulas avoid cancellation at the smaller root.
    left <- (4*D/(z+s))^2
    right <- ((z+s)/(4*gamma))^2
    # The open interval (left, right) is impossible.  Retain a small margin
    # around each endpoint to guard against floating-point rounding.
    end1 <- min(last, ceiling(left + pad))
    start2 <- max(first, end1 + 1, floor(right - pad))
    intervals <- list(c(first, end1), c(start2, last))
  }

  # Search retained intervals in increasing order; feasibility is not monotone
  # in B, so an ordinary binary search in B would be incorrect.
  for (interval in intervals) {
    low <- interval[1]
    stop_at <- interval[2]
    while (low <= stop_at) {
      high <- min(stop_at, low + chunk_size - 1)
      B <- seq(from = low, to = high, by = 1)
      # q_B <= floor(B/2 + gamma*B).  Keep the parity term separate and avoid
      # the potentially misleading subtraction q_B/B - 1/2.
      k <- floor(B/2) + floor(gamma*B + 0.5*(B %% 2))
      # Both requirements together are q_B <= min(k, B-1).
      k <- pmin(k, B - 1)
      hits <- which(.cosi_tail_leq(B, k, alpha))
      if (length(hits)) return(B[hits[1]])
      low <- high + 1
    }
  }
  stop("No solution within the image's upper bound; check numerical accuracy")
}

# Examples:
# smallest_B(0.05, 1/2)       # 6 (previously 1 without the strict condition)
# smallest_B(0.05, 1/3)       # 6
# smallest_B(0.05, 0.10)      # 90
# smallest_B(0.05, 0.01)      # 9551
# largest_support(6, 0.05)   # 5
