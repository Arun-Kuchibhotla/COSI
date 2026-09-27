# COSI intervals from independent batch estimates.
# First source("cosi_smallest_B.R"), then source this file.
# Uses base R only.  Use set.seed(...) for reproducible splitting and cutoffs.

.cosi_ci_control <- function(control) {
  defaults <- list(max_states = 2^20, max_walsh = 2e6,
                   max_dp_work = 5e7, mass_tol = 64*.Machine$double.eps,
                   max_iterations = 256, concentration_B = 40)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% names(defaults))))) {
    stop("control must be a named list using documented control options")
  }
  defaults[names(control)] <- control
  for (nm in names(defaults)) {
    z <- defaults[[nm]]
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
        z < 0 || (nm != "mass_tol" && (z < 1 || z != floor(z)))) {
      stop(paste0("Invalid control$", nm))
    }
  }
  defaults
}

.cosi_draw_law <- function(law, randomized, conservative = c("max", "min")) {
  conservative <- match.arg(conservative)
  if (!randomized) {
    return(if (conservative == "max") max(law$support) else min(law$support))
  }
  if (length(law$support) == 1L) return(law$support)
  law$support[sample.int(length(law$support), 1L, prob = law$probabilities)]
}

.cosi_kappa_law <- function(c_minus, c_plus, P_minus, P_plus, alpha, engine) {
  # Larger kappa rejects more often, so its conservative choice is the
  # SMALLEST support value.  This direction differs from the binomial cutoff.
  if (P_minus >= alpha || P_plus < alpha || P_plus <= P_minus) {
    stop("Could not resolve the signed-rank critical-value bracket numerically")
  }
  tau <- (alpha - P_minus)/(P_plus - P_minus)
  probabilities <- c(1 - tau, tau)
  keep <- probabilities > 0
  list(support = c(c_minus, c_plus)[keep], probabilities = probabilities[keep],
       rejection_probabilities = c(P_minus, P_plus)[keep],
       c_minus = c_minus, c_plus = c_plus, P_minus = P_minus, P_plus = P_plus,
       tau = tau, engine = engine, exact = TRUE)
}

.cosi_log_cosh <- function(x) {
  # log(cosh(x)) without overflow, or cancellation when x is small.
  x <- abs(x)
  out <- x + log1p(exp(-2*x)) - log(2)
  small <- x < 1
  out[small] <- log1p(2*sinh(x[small]/2)^2)
  out
}

.cosi_concentration_calibration <- function(scores, alpha, control,
                                            reason = "general real scores") {
  # Conservative calibration, with no enumeration or Monte Carlo sampling.
  # Normalize so T=(G-N/2)/(N/2)=sum w_i*Rademacher_i. This prevents
  # overflow in sum(scores^2) and makes the optimization scale invariant.
  N <- sum(scores)
  if (N == 0) return(.cosi_kappa_law(-1, 0, 0, 1, alpha, "zero-scores"))
  mu <- N/2
  w <- scores[scores > 0]/N
  variance <- sum(w^2)
  sd <- sqrt(variance)
  log_alpha <- log(alpha)
  L <- log(2) - log_alpha

  # Four rigorous radii for P(|T| >= radius) <= alpha. The Gaussian
  # comparison is Bentkus--Dzindzalieta, Theorem 1.1 (2015), with its
  # optimal universal multiplicative constant; it is NOT a normal CLT.
  # https://arxiv.org/abs/1307.3451
  log_C <- -log(4) - pnorm(sqrt(2), lower.tail = FALSE, log.p = TRUE)
  r_hoeffding <- sd*sqrt(2*L)
  r_gaussian <- sd*qnorm(log_alpha - log(2) - log_C,
                         lower.tail = FALSE, log.p = TRUE)
  r_chebyshev <- sd*exp(-log_alpha/2)

  # Chernoff radius: [L + sum log cosh(lambda*w_i)]/lambda.
  # Its derivative changes sign where lambda*K'(lambda)-K(lambda)=L.
  # For large z, evaluate z*tanh(z)-log(cosh(z)) without subtracting
  # two large, nearly equal numbers. Zero normalized weights are harmless.
  slope <- function(lambda) {
    z <- lambda*w
    part <- numeric(length(z))
    small <- z < 20
    part[small] <- z[small]*tanh(z[small]) - .cosi_log_cosh(z[small])
    zz <- z[!small]
    e <- exp(-2*zz)
    part[!small] <- log(2) - log1p(e) - 2*zz*e/(1 + e)
    sum(part) - L
  }
  lambda_H <- sqrt(2*L/variance)
  candidates <- lambda_H
  hi <- lambda_H
  f_hi <- slope(hi)
  for (iteration in seq_len(control$max_iterations)) {
    if (f_hi >= 0 || hi > .Machine$double.xmax/2) break
    hi <- 2*hi
    f_hi <- slope(hi)
  }
  candidates <- c(candidates, hi)
  if (is.finite(f_hi) && f_hi >= 0) {
    root <- tryCatch(uniroot(slope, c(0, hi), f.lower = -L, f.upper = f_hi,
                            tol = sqrt(.Machine$double.eps),
                            maxiter = control$max_iterations,
                            check.conv = TRUE)$root,
                     error = function(e) NA_real_)
    if (is.finite(root) && root > 0) candidates <- c(candidates, root)
  }
  # Every FINITE lambda gives a valid radius, even if the optimizer does
  # not converge. Do not replace lambda by infinity: at the edge of a
  # discrete support, a limiting radius can change an inclusive tail event.
  candidate_radii <- vapply(candidates, function(lambda)
    L/lambda + sum(.cosi_log_cosh(lambda*w)/lambda), numeric(1))
  chosen <- which.min(candidate_radii)
  lambda <- candidates[chosen]
  K <- sum(.cosi_log_cosh(lambda*w))
  radii <- c(chernoff = candidate_radii[chosen], gaussian = r_gaussian,
             hoeffding = r_hoeffding, chebyshev = r_chebyshev)
  best <- which.min(radii)

  # A smaller kappa gives a wider interval. Subtract a numerical cushion
  # that also covers the Walsh cumulative-weight tolerance used below.
  # This is a double-precision safeguard, not an arbitrary-precision proof.
  cushion <- (control$mass_tol + 256*.Machine$double.eps)*max(1, N)
  raw_kappa <- mu*(1 - radii[best])
  kappa <- if (is.finite(raw_kappa) && raw_kappa > cushion)
             unname(raw_kappa - cushion) else -1

  # Report a bound for the rejection event of the IMPLEMENTED interval.
  # Account for a possibly inward Walsh mass tolerance when checking it.
  log_bound <- -Inf
  if (kappa >= 0) {
    effective_kappa <- kappa + control$mass_tol*max(1, N)
    r <- (mu - effective_kappa)/mu
    bounds <- c(chernoff = log(2) + K - lambda*r,
                gaussian = log(2) + log_C +
                  pnorm(r/sd, lower.tail = FALSE, log.p = TRUE),
                hoeffding = log(2) - r^2/(2*variance),
                chebyshev = log(variance) - 2*log(r))
    log_bound <- min(0, bounds)
    # If extreme floating-point arithmetic defeats the check, the whole
    # real line remains a valid deterministic fallback.
    if (is.na(log_bound) || log_bound > log_alpha) {
      kappa <- -1
      log_bound <- -Inf
    }
  }
  list(support = kappa, probabilities = 1, rejection_probabilities = NA_real_,
       engine = "concentration", exact = FALSE, bound_used = names(radii)[best],
       rejection_probability_bound = exp(log_bound),
       log_rejection_probability_bound = log_bound,
       radius = unname(mu*radii[best]), radius_fraction = unname(radii[best]),
       candidate_radii = radii*mu, lambda = lambda,
       numerical_cushion = cushion, fallback_reason = reason)
}

.cosi_law_from_pmf <- function(values, probabilities, N, alpha, engine) {
  cdf <- cumsum(probabilities)
  k <- which(cdf >= alpha/2)[1L]
  if (is.na(k)) stop("Failed to locate a signed-rank quantile")
  c_plus <- values[k]
  c_minus <- if (k == 1L) -1 else values[k - 1L]
  # The two tails overlap at N/2: P_h(N/2) is 1, not 2*F_G(N/2).
  P_plus <- min(1, 2*cdf[k])
  P_minus <- if (k == 1L) 0 else min(1, 2*cdf[k - 1L])
  .cosi_kappa_law(c_minus, c_plus, P_minus, P_plus, alpha, engine)
}

.cosi_subset_sums <- function(scores) {
  sums <- 0
  for (s in scores) sums <- c(sums, sums + s)
  runs <- rle(sort(sums))
  list(values = runs$values, counts = as.numeric(runs$lengths))
}

.cosi_mitm_calibration <- function(scores, N, alpha, control) {
  # Split the Bernoulli sum into two independent halves.  Only the two
  # half-distributions are stored, rather than all 2^length(scores) sums.
  n <- length(scores)
  m <- floor(n/2)
  needed <- 2^ceiling(n/2)
  if (needed > control$max_states || n > 52) {
    return(.cosi_concentration_calibration(scores, alpha, control,
                    "meet-in-the-middle state/counting limit"))
  }
  left <- .cosi_subset_sums(scores[seq_len(m)])
  right <- .cosi_subset_sums(scores[seq.int(m + 1L, n)])
  a <- left$values
  b <- right$values
  rb <- c(0, cumsum(right$counts))
  nb <- length(b)

  count_pairs <- function(x, strict = FALSE, positions = FALSE) {
    j <- findInterval(x - a, b)
    # Correct subtraction-rounding at a pair-sum boundary using addition.
    # Keeping distinct right values makes this cheap even with repeated sums.
    repeat {
      ii <- which(j > 0)
      bad <- if (strict) a[ii] + b[j[ii]] >= x else a[ii] + b[j[ii]] > x
      ii <- ii[bad]
      if (!length(ii)) break
      j[ii] <- j[ii] - 1L
    }
    repeat {
      ii <- which(j < nb)
      good <- if (strict) a[ii] + b[j[ii] + 1L] < x else a[ii] + b[j[ii] + 1L] <= x
      ii <- ii[good]
      if (!length(ii)) break
      j[ii] <- j[ii] + 1L
    }
    count <- sum(left$counts * rb[j + 1L])
    if (positions) list(count = count, j = j) else count
  }

  total <- 2^n
  rank <- max(1, ceiling(alpha * 2^(n - 1L)))
  if (count_pairs(0) >= rank) {
    c_plus <- 0
  } else {
    lo <- 0                         # count(lo) < rank
    hi <- max(a) + max(b)            # count(hi) >= rank
    converged <- FALSE
    for (iteration in seq_len(control$max_iterations)) {
      mid <- lo + (hi - lo)/2
      if (mid == lo || mid == hi) {
        converged <- TRUE
        break
      }
      z <- count_pairs(mid, positions = TRUE)
      if (z$count < rank) {
        lo <- mid
      } else {
        ii <- which(z$j > 0)
        hi <- max(a[ii] + b[z$j[ii]]) # tighten to an actual support value
      }
    }
    if (!converged) return(.cosi_concentration_calibration(scores, alpha, control,
                                            "pair-sum selection iteration limit"))
    z <- count_pairs(lo, positions = TRUE)
    ii <- which(z$j < nb)
    c_plus <- min(a[ii] + b[z$j[ii] + 1L])
  }
  z <- count_pairs(c_plus, strict = TRUE, positions = TRUE)
  ii <- which(z$j > 0)
  c_minus <- if (!length(ii)) -1 else max(a[ii] + b[z$j[ii]])
  P_minus <- 2*z$count/total
  P_plus <- min(1, 2*count_pairs(c_plus)/total)
  .cosi_kappa_law(c_minus, c_plus, P_minus, P_plus, alpha, "meet-in-the-middle")
}

.cosi_sr_calibration <- function(scores, alpha, control) {
  # Distribution-based calibration of G=sum scores[i]*Bernoulli(1/2).
  # Prefer fast exact discrete algorithms. Large general problems use
  # rigorous concentration bounds, evaluated in floating-point arithmetic.
  N <- sum(scores)
  B <- length(scores)
  positive <- scores[scores > 0]
  if (!length(positive)) {
    return(.cosi_kappa_law(-1, 0, 0, 1, alpha, "zero-scores"))
  }
  # This endpoint probability has a closed form, even when alpha/2 would
  # underflow.  Zero scores do not affect the number of sign configurations.
  P_zero <- 2^(1 - length(positive))
  if (alpha <= P_zero) {
    return(.cosi_kappa_law(-1, 0, 0, P_zero, alpha, "zero-threshold"))
  }
  if (alpha/2 == 0) {
    return(.cosi_concentration_calibration(scores, alpha, control,
                                          "discrete tail probability underflow"))
  }
  if (B <= 1023 && scores[1L] > 0 && all(scores == scores[1L]*seq_len(B))) {
    # h(i)=a*i: use R's optimized Wilcoxon signed-rank distribution.
    k <- qsignrank(alpha/2, n = B)
    # Quantile routines may fuzz probabilities at a jump.  Restore the
    # inclusive definition by checking the neighboring CDF values directly.
    P_at <- function(j) if (j < 0) 0 else min(1, 2*psignrank(j, n = B))
    while (P_at(k) < alpha) k <- k + 1
    while (k > 0 && P_at(k - 1) >= alpha) k <- k - 1
    cp <- scores[1L]*k
    cm <- if (k == 0) -1 else scores[1L]*(k - 1)
    pm <- if (k == 0) 0 else 2*psignrank(k - 1, n = B)
    pp <- min(1, 2*psignrank(k, n = B))
    return(.cosi_kappa_law(cm, cp, pm, pp, alpha, "signed-rank"))
  }
  lattice_step <- min(positive)
  ratios <- positive/lattice_step
  if (all(positive == floor(positive)) && max(positive) <= 2^52) {
    gcd <- function(a, b) {
      while (b != 0) { r <- a %% b; a <- b; b <- r }
      a
    }
    lattice_step <- Reduce(gcd, positive)
    ratios <- positive/lattice_step
  }
  # Also recognize real multiples such as h(i)=i^2/10, without approximate
  # rounding of an arbitrary score function onto a lattice.
  if (all(ratios == floor(ratios)) && all(lattice_step*ratios == positive)) {
    step <- lattice_step
    weights <- ratios
    states <- sum(weights) + 1
    if (states <= control$max_states &&
        length(weights)*states <= control$max_dp_work) {
      # Dynamic programming on an integer lattice, allowing gaps in support.
      pmf <- 1
      for (w in weights) {
        old <- pmf/2
        len <- length(old)
        pmf <- numeric(len + w)
        pmf[seq_len(len)] <- old
        ii <- seq_len(len) + w
        pmf[ii] <- pmf[ii] + old
      }
      keep <- which(pmf > 0)
      return(.cosi_law_from_pmf((keep - 1)*step, pmf[keep], N, alpha,
                                "lattice-dynamic-programming"))
    }
  }
  if (length(positive) <= min(14, log2(control$max_states))) {
    law <- .cosi_subset_sums(positive)
    return(.cosi_law_from_pmf(law$values, law$counts/2^length(positive),
                              N, alpha, "enumeration"))
  }
  if (B >= control$concentration_B) {
    return(.cosi_concentration_calibration(scores, alpha, control,
                                          "B reached concentration_B"))
  }
  .cosi_mitm_calibration(positive, N, alpha, control)
}

.cosi_interval <- function(lower, upper, lower_closed = TRUE, upper_closed = TRUE) {
  # Closed/open endpoints matter after projecting lexicographic inequalities
  # onto the real line.  Retain the endpoints even if they describe the empty set.
  if (anyNA(c(lower, upper))) stop("Undefined interval endpoint")
  closed <- c(lower = is.finite(lower) && lower_closed,
              upper = is.finite(upper) && upper_closed)
  empty <- lower > upper || (lower == upper && !all(closed))
  list(CI = c(lower = lower, upper = upper), closed = closed, empty = empty)
}

.cosi_order_value <- function(x, k) {
  if (k <= 0) return(-Inf)
  if (k > length(x)) return(Inf)
  sort(x, partial = k)[k]
}

.cosi_walsh_interval <- function(W, U, h_values, kappa, control) {
  B <- length(W)
  N <- sum(h_values[-1L])
  # Augmented order is lexicographic, including when estimates tie.
  ord <- order(W, U, method = "radix")
  W <- W[ord]
  U <- U[ord]
  count <- B*(B + 1)/2
  av <- au <- weights <- numeric(count)
  increment <- diff(h_values)
  at <- 0
  for (i in seq_len(B)) {
    j <- seq.int(i, B)
    ii <- at + seq_along(j)
    # Dividing before addition avoids unnecessary overflow of W_i+W_j.
    av[ii] <- W[i]/2 + W[j]/2
    au[ii] <- U[i]/2 + U[j]/2
    weights[ii] <- increment[j - i + 1L]
    at <- at + length(j)
  }
  ord <- order(av, au, method = "radix")
  av <- av[ord]
  au <- au[ord]
  weights <- weights[ord]
  cumulative <- cumsum(weights)
  # Integers have exact cumulative masses in the supported range.  For real
  # scores, allow a small mass tolerance for telescoping-sum roundoff; never
  # use a tolerance to merge or reorder the (average, tie-breaker) pairs.
  tol <- if (all(h_values == floor(h_values)) && N <= 2^52) 0 else
           control$mass_tol * max(1, N)
  if (abs(tail(cumulative, 1L) - N) > max(tol, 8*.Machine$double.eps*max(1, N))) {
    stop("Total Walsh weight differs from sum h(i); check h(0)=0")
  }
  cumulative[length(cumulative)] <- N

  Q <- function(gamma, right = FALSE) {
    if ((!right && gamma <= 0) || (right && gamma < 0)) return(c(-Inf, 0))
    if (gamma > N || (right && gamma >= N)) return(c(Inf, 0))
    # Q(gamma+) uses strictly greater cumulative weight.  It is not generally
    # Q(gamma+1); that shortcut is invalid for the normal-score method.
    if (right && gamma == 0) {
      ii <- which(cumulative > 0)
    } else if (right) {
      ii <- which(cumulative > gamma + tol)
    } else {
      ii <- which(cumulative > 0 & cumulative >= gamma - tol)
    }
    if (!length(ii)) return(c(Inf, 0))
    c(av[ii[1L]], au[ii[1L]])
  }

  lower <- Q(kappa, right = TRUE)
  upper <- Q(N - kappa)
  # (l,u_l) < (t,0) includes t=l iff u_l<0.  The upper endpoint is
  # included iff u_u>0.  Equal real endpoints can therefore form a singleton.
  result <- .cosi_interval(lower[1L], upper[1L], lower[2L] < 0, upper[2L] > 0)
  result$augmented_endpoints <- rbind(lower = lower, upper = upper)
  colnames(result$augmented_endpoints) <- c("value", "tie_breaker")
  result$mass_tolerance <- tol
  result
}

cosi_ci <- function(data, estimate, alpha = 0.05, gamma = 1/3,
                    method = c("median", "symmetry", "normality", "unimodal-star"),
                    randomized = TRUE, h = function(x) x, control = list()) {
  # Main interface.  Rows are observations; estimate always receives a matrix.
  # randomized toggles critical-value randomization ONLY.  Splitting and the
  # U(-1,1) tie-breakers remain random, including when randomized=FALSE.
  # Concentration calibration uses one deterministic kappa for either setting.
  if (!exists("smallest_B_distribution", mode = "function") ||
      !exists("binomial_cutoff_distribution", mode = "function")) {
    stop("First source the accompanying revised cosi_smallest_B.R")
  }
  if (!is.function(estimate)) stop("estimate must be a function")
  .cosi_check_alpha(alpha)
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
      gamma <= 0 || gamma >= 0.5) stop("gamma must lie strictly between 0 and 1/2")
  if (!is.logical(randomized) || length(randomized) != 1L || is.na(randomized)) {
    stop("randomized must be TRUE or FALSE")
  }
  if (identical(method, "symmetric")) method <- "symmetry"
  if (identical(method, "unimodal")) method <- "unimodal-star"
  method <- match.arg(method, c("median", "symmetry", "normality", "unimodal-star"))
  control <- .cosi_ci_control(control)
  data <- as.matrix(data)
  nn <- nrow(data)
  if (is.null(nn) || nn < 1 || ncol(data) < 1) stop("data must have rows and columns")

  # This interface consumes no random numbers.  In particular, setting
  # randomized=FALSE does not first generate and discard a randomized cutoff.
  split_calibration <- smallest_B_distribution(alpha, gamma)
  B <- split_calibration$B
  batches <- B + as.integer(method == "unimodal-star")
  if (batches > nn) {
    stop(sprintf("Not enough samples for splitting: need %.0f batches; sample size is %d (alpha=%g, gamma=%g)",
                 batches, nn, alpha, gamma))
  }

  cutoff <- NULL
  cutoff_law <- NULL
  kappa <- NULL
  kappa_law <- NULL
  h_values <- NULL
  if (method == "median") {
    cutoff_law <- split_calibration
    cutoff <- .cosi_draw_law(cutoff_law, randomized, "max")
  } else if (method == "unimodal-star") {
    # Follow the user's final formula: c*_{B,alpha}, NOT c*_{B,alpha/2}.
    # B itself is still chosen using the alpha/2 split-calibration constraint.
    cutoff_law <- binomial_cutoff_distribution(B, level = alpha)
    cutoff <- .cosi_draw_law(cutoff_law, randomized, "max")
  } else {
    if (B*(B + 1)/2 > control$max_walsh) {
      stop(sprintf("Need %.0f Walsh pairs; control$max_walsh is %.0f",
                   B*(B + 1)/2, control$max_walsh))
    }
    if (method == "normality") {
      h_values <- qnorm(1/2 + (0:B)/(2*(B + 1)))
    } else {
      if (!is.function(h)) stop("h must be a function for method='symmetry'")
      # h need not be vectorized; only its values at 0,...,B are needed.
      h_values <- vapply(0:B, function(i) {
        z <- h(i)
        if (!is.numeric(z) || length(z) != 1L || !is.finite(z)) {
          stop("h must return one finite numeric value at each of 0,...,B")
        }
        as.numeric(z)
      }, numeric(1))
    }
    if (h_values[1L] != 0) {
      stop("The stated Walsh weights and N_h require h(0)=0")
    }
    if (any(h_values < 0) || any(diff(h_values) < 0) || !is.finite(sum(h_values))) {
      stop("h(0),...,h(B) must be nonnegative, nondecreasing, with finite sum")
    }
    kappa_law <- .cosi_sr_calibration(h_values[-1L], alpha, control)
    kappa <- .cosi_draw_law(kappa_law, randomized, "min")
  }

  # Shuffle once and create balanced, disjoint batches exactly as specified.
  permutation <- sample.int(nn)
  shuffled <- data[permutation, , drop = FALSE]
  TMP <- split(seq_len(nn), sort(seq_len(nn) %% batches))
  all_estimates <- vapply(seq_along(TMP), function(idx) {
    value <- tryCatch(estimate(shuffled[TMP[[idx]], , drop = FALSE]),
                      error = function(e) stop(sprintf("estimate failed on batch %d: %s", idx, conditionMessage(e))))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop(sprintf("estimate must return one finite numeric value; failed on batch %d", idx))
    }
    as.numeric(value)
  }, numeric(1))

  anchor <- NULL
  uniforms <- NULL
  if (method == "median") {
    W <- all_estimates
    interval <- .cosi_interval(.cosi_order_value(W, B - cutoff),
                               .cosi_order_value(W, cutoff + 1))
  } else if (method == "unimodal-star") {
    anchor <- all_estimates[1L]
    W <- all_estimates[-1L]
    d <- (1 + sqrt(2))*abs(W - anchor)
    L <- pmin(W, anchor) - d
    U <- pmax(W, anchor) + d
    interval <- .cosi_interval(.cosi_order_value(L, B - cutoff),
                               .cosi_order_value(U, cutoff + 1))
    interval$L <- L
    interval$U <- U
  } else {
    W <- all_estimates
    uniforms <- runif(B, -1, 1)
    # Exact zero has probability zero in the mathematical construction.
    while (any(uniforms == 0)) {
      ii <- which(uniforms == 0)
      uniforms[ii] <- runif(length(ii), -1, 1)
    }
    interval <- .cosi_walsh_interval(W, uniforms, h_values, kappa, control)
  }

  result <- c(interval, list(B = B, batches = batches, ci_est = W,
                 all_estimates = all_estimates, anchor = anchor, uniforms = uniforms,
                 c = cutoff, cutoff_distribution = cutoff_law,
                 kappa = kappa, kappa_distribution = kappa_law,
                 h_values = h_values, N_h = if (is.null(h_values)) NULL else sum(h_values),
                 split_calibration = split_calibration,
                 batch_sizes = lengths(TMP),
                 batch_rows = lapply(TMP, function(ii) permutation[ii]),
                 alpha = alpha, gamma = gamma, method = method, randomized = randomized))
  class(result) <- "cosi_ci"
  result
}

contains_cosi_ci <- function(object, t) {
  # Membership respects open/closed endpoints and the empty-set flag.
  if (!inherits(object, "cosi_ci")) stop("object must be a cosi_ci result")
  if (!is.numeric(t)) stop("t must be numeric")
  if (object$empty) return(rep(FALSE, length(t)))
  low <- if (object$closed[1L]) t >= object$CI[1L] else t > object$CI[1L]
  high <- if (object$closed[2L]) t <= object$CI[2L] else t < object$CI[2L]
  low & high & is.finite(t)
}

print.cosi_ci <- function(x, ...) {
  cat(sprintf("COSI %s interval; B=%g, batches=%g, alpha=%g\n",
              x$method, x$B, x$batches, x$alpha))
  if (identical(x$kappa_distribution$engine, "concentration")) {
    cat(sprintf("Conservative concentration calibration (%s); rejection bound <= %.6g\n",
                x$kappa_distribution$bound_used,
                x$kappa_distribution$rejection_probability_bound))
  }
  if (x$empty) {
    cat("empty set\n")
  } else {
    cat(if (x$closed[1L]) "[" else "(",
        format(x$CI[1L]), ", ", format(x$CI[2L]),
        if (x$closed[2L]) "]\n" else ")\n", sep = "")
  }
  invisible(x)
}
