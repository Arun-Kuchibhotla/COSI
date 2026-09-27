# Find the smallest B for the randomized binomial cutoff in COSI.
# Uses only base R (the stats package).  Source this file, then call
# smallest_B(alpha = 0.05, gamma = 1/3).
#
# Independent randomization is assumed.  The largest support value is
# q_B = min{k: P(Binom(B, 1/2) > k) <= alpha/2}.
# The randomization may degenerate to a single integer at an exact CDF level.
# The revised search also requires q_B < B.
# smallest_B returns B, a randomized c_star, and its full distribution.
# smallest_B_distribution returns the deterministic calibration without a draw.

# Validate that the significance level is finite and strictly between zero and one.
.cosi_check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a finite number in (0, 1)")
  }
}

# Test the inclusive binomial upper-tail bound used in cutoff selection.
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

# Find the largest positive-probability support value of the binomial cutoff at alpha/2.
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

# Find the smallest feasible B using the image bounds to prune candidates, then exact tail checks.
.cosi_smallest_B_value <- function(alpha, gamma, chunk_size = 65536) {
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

# Compute the mixing probabilities for the two adjacent binomial cutoff values at alpha/2.
.cosi_cutoff_probabilities <- function(B, alpha) {
  # If S(k)=P(Binom(B,1/2)>k), m=P(Binom(B,1/2)=q), and q=q_B, then
  # p_lower=(alpha/2-S(q))/m and p_upper=(S(q-1)-alpha/2)/m.
  # Their mixture has upper-tail probability alpha/2 in exact arithmetic.
  if (2^(1-B) <= alpha && (B == 2 || 2^(2-B) > alpha)) {
    # At the smallest B allowed by q_B<B, we know q_B=B-1.  This formula
    # avoids alpha/2 underflow, even for the smallest positive double.
    q <- B - 1
    p_lower <- (alpha / 2^(1-B) - 1) / B
    p_upper <- 1 - p_lower
  } else {
    q <- largest_support(B, alpha)
    if (B <= 52) {
      # Pascal's triangle yields exact dyadic masses at small B, so exact
      # CDF levels give a genuinely degenerate cutoff instead of a tiny atom.
      row <- 1
      for (j in seq_len(B)) row <- c(row, 0) + c(0, row)
      mass <- row[q + 1] / 2^B
      m <- B - q
      tail <- if (m <= 0) 0 else sum(row[seq_len(m)]) / 2^B
      p_lower <- (alpha/2 - tail) / mass
      p_upper <- (tail + mass - alpha/2) / mass
    } else {
      # Log upper tails avoid subtracting CDFs from one.  expm1 evaluates
      # exp(x)-1 accurately when its argument is close to zero.
      log_tail <- pbinom(c(q, q - 1), B, 0.5,
                        lower.tail = FALSE, log.p = TRUE)
      log_target <- log(alpha) - log(2)
      log_mass <- dbinom(q, B, 0.5, log = TRUE)
      p_lower <- exp(log_target - log_mass) * (-expm1(log_tail[1] - log_target))
      p_upper <- exp(log_tail[2] - log_mass) * (-expm1(log_target - log_tail[2]))
    }
  }

  # Reject material errors; trim roundoff at endpoints and compute the larger
  # probability as a complement to keep the probabilities summing to one.
  p <- c(p_lower, p_upper)
  if (any(!is.finite(p)) || any(p < -1e-10 | p > 1 + 1e-10)) {
    stop("Could not resolve cutoff probabilities in double precision")
  }
  p <- pmin(1, pmax(0, p))
  if (p[1] <= p[2]) p[2] <- 1 - p[1] else p[1] <- 1 - p[2]
  if (p[2] <= 0) stop("Upper support probability is numerically unresolved")
  list(q = q, probabilities = p)
}

# Return the smallest feasible B, its cutoff distribution, and one randomized cutoff draw.
smallest_B <- function(alpha, gamma, chunk_size = 65536, seed = NULL) {
  # Return a list with B, c_star, largest_support, support, probabilities.
  # c_star is ONE randomized draw; B and the distribution are deterministic.
  # The constraints concern the largest support value, not an individual draw.
  # Use result$B to obtain the integer returned by the previous version.
  #
  # seed=NULL uses R's current random stream (so set.seed() works normally).
  # A supplied seed makes a reproducible local draw and restores the caller's
  # RNG state.  For successive draws, call set.seed() once, then omit seed.
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max)) {
    stop("seed must be NULL or an integer between 0 and .Machine$integer.max")
  }
  B <- .cosi_smallest_B_value(alpha, gamma, chunk_size)
  calibration <- .cosi_cutoff_probabilities(B, alpha)
  q <- calibration$q
  p <- calibration$probabilities

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  if (p[1] == 0) {
    c_star <- q
    support <- q
    probabilities <- 1
  } else {
    # Sample the less likely outcome directly for numerical stability.
    c_star <- if (p[1] <= p[2]) q - rbinom(1, 1, p[1]) else
                                q - 1 + rbinom(1, 1, p[2])
    support <- c(q - 1, q)
    probabilities <- p
  }
  list(B = B, c_star = c_star, largest_support = q,
       support = support, probabilities = probabilities)
}

# Examples:
# result <- smallest_B(0.05, 1/3, seed = 42)
# result$B                   # 6
# result$c_star              # a draw in {4,5}
# result$support             # c(4,5)
# result$probabilities       # c(0.1,0.9)
# smallest_B(0.05, 0.10)$B   # 90
# largest_support(6, 0.05)   # 5

# Return the smallest feasible B and its cutoff support and probabilities without drawing a cutoff.
smallest_B_distribution <- function(alpha, gamma, chunk_size = 65536) {
  # RNG-free interface for confidence-interval construction.  Keeping the
  # distribution separate lets a caller choose randomization or conservatism.
  B <- .cosi_smallest_B_value(alpha, gamma, chunk_size)
  calibration <- .cosi_cutoff_probabilities(B, alpha)
  q <- calibration$q
  p <- calibration$probabilities
  keep <- p > 0
  list(B = B, largest_support = q,
       support = c(q - 1, q)[keep], probabilities = p[keep])
}

# Return cutoff support and probabilities at a supplied binomial tail level.
binomial_cutoff_distribution <- function(B, level) {
  # Distribution of c*_{B,level} for ANY level in (0,1), with independent
  # randomization satisfying P(Binom(B,1/2) <= c*) = 1-level.
  # In particular, level=alpha is needed by the unimodal-star construction.
  # Its support can include -1 when level is large; this is intentional.
  .cosi_check_alpha(level)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > 2^52) {
    stop("B must be a positive integer no larger than 2^52")
  }
  if (B <= 52) {
    # Exact dyadic probabilities, including equality at a CDF jump.
    row <- 1
    for (j in seq_len(B)) row <- c(row, 0) + c(0, row)
    # Compute the binomial upper-tail probability at a candidate cutoff.
    tail_at <- function(k) {
      if (k < 0) return(1)
      if (k >= B) return(0)
      sum(row[seq_len(B - k)]) / 2^B
    }
    # Test the inclusive binomial upper-tail bound used in cutoff selection.
    tail_leq <- function(k) tail_at(k) <= level
  } else {
    # Test the inclusive binomial upper-tail bound used in cutoff selection.
    tail_leq <- function(k) pbinom(k, B, 0.5, lower.tail = FALSE,
                                  log.p = TRUE) <= log(level)
  }
  # Search in k; unlike a search in B, this predicate is monotone.
  lo <- -1
  hi <- B
  while (hi - lo > 1) {
    mid <- floor(lo + (hi - lo)/2)
    if (tail_leq(mid)) hi <- mid else lo <- mid
  }
  q <- hi
  if (B <= 52) {
    mass <- row[q + 1] / 2^B
    p <- c((level - tail_at(q))/mass,
           (tail_at(q - 1) - level)/mass)
  } else {
    lt <- pbinom(c(q, q - 1), B, 0.5, lower.tail = FALSE, log.p = TRUE)
    lm <- dbinom(q, B, 0.5, log = TRUE)
    la <- log(level)
    p <- c(exp(la - lm) * (-expm1(lt[1] - la)),
           exp(lt[2] - lm) * (-expm1(la - lt[2])))
  }
  if (any(!is.finite(p)) || any(p < -1e-10 | p > 1 + 1e-10)) {
    stop("Could not resolve binomial cutoff probabilities in double precision")
  }
  p <- pmin(1, pmax(0, p))
  if (p[1] <= p[2]) p[2] <- 1 - p[1] else p[1] <- 1 - p[2]
  if (p[2] <= 0) stop("Upper binomial support probability is unresolved")
  keep <- p > 0
  list(B = B, level = level, largest_support = q,
       support = c(q - 1, q)[keep], probabilities = p[keep])
}
