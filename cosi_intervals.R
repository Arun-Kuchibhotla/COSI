# Four confidence-interval constructors from supplied batch estimates.
# Load through source("cosi_ci.R"), or source cosi_calibration.R first.

# Validate a nonempty vector of finite scalar batch estimates; return numeric W.
.cosi_validate_estimators <- function(W) {
  if (!is.numeric(W) || !is.null(dim(W)) || !length(W) || any(!is.finite(W))) {
    stop("W must be a nonempty numeric vector of finite batch estimates")
  }
  as.numeric(W)
}

# Evaluate/validate h(0),...,h(B), or return the prescribed normal scores.
.cosi_score_values <- function(B, method = "symmetry", h = function(i) i) {
  values <- if (method == "normality") {
    qnorm(0.5 + (0:B)/(2*(B + 1)))
  } else if (is.function(h)) {
    vapply(0:B, function(i) {
      z <- h(i)
      if (!is.numeric(z) || length(z) != 1L || !is.finite(z)) {
        stop("h must return one finite numeric value at each integer 0,...,B")
      }
      as.numeric(z)
    }, numeric(1))
  } else h
  if (!is.numeric(values) || !is.null(dim(values)) || length(values) != B + 1 ||
      any(!is.finite(values)) || any(values < 0) || any(diff(values) < 0) ||
      !is.finite(sum(values))) stop("h must give B+1 finite, nonnegative, nondecreasing scores")
  if (values[1L] != 0) stop("The stated Walsh weights and N_h require h(0)=0")
  as.numeric(values)
}

# Validate c as a binomial order-statistic cutoff or signed-rank threshold.
.cosi_check_cutoff <- function(cutoff, B, N = NULL) {
  if (!is.numeric(cutoff) || length(cutoff) != 1 || !is.finite(cutoff)) {
    stop("c must be one finite numeric cutoff")
  }
  if (is.null(N)) {
    if (cutoff != floor(cutoff) || cutoff < -1 || cutoff > B) stop("Binomial c must be an integer in [-1,B]")
  } else if (!(cutoff == -1 || (cutoff >= 0 && cutoff <= N/2))) {
    stop("Signed-rank c must be -1 or lie in [0,N_h/2]; it denotes kappa")
  }
  invisible(NULL)
}

# Return supplied tie-breakers or generate B independent nonzero U(-1,1) draws.
.cosi_tie_breakers <- function(B, U = NULL) {
  if (is.null(U)) {
    U <- runif(B, -1, 1)
    while (any(U == 0)) U[U == 0] <- runif(sum(U == 0), -1, 1)
  }
  if (!is.numeric(U) || length(U) != B || any(!is.finite(U)) ||
      any(U < -1 | U > 1 | U == 0)) stop("U must contain B nonzero values in [-1,1]")
  as.numeric(U)
}

# Construct [W_(B-c), W_(c+1)] using c drawn from c*_{B,alpha/2}.
# This function only uses supplied estimates and cutoff; it does no calibration.
cosi_median <- function(W, c) {
  W <- .cosi_validate_estimators(W)
  B <- length(W)
  .cosi_check_cutoff(c, B)
  .cosi_interval(.cosi_order_value(W, B-c), .cosi_order_value(W, c+1))
}

# Construct the weighted signed-rank interval from W, h, and optional U.
# Here c MUST be the signed-rank cutoff kappa_{B,alpha}^h, not a binomial c*.
cosi_symmetry <- function(W, c, h = function(i) i, U = NULL, control = list()) {
  W <- .cosi_validate_estimators(W)
  B <- length(W)
  control <- .cosi_ci_control(control)
  h_values <- .cosi_score_values(B, "symmetry", h)
  .cosi_check_cutoff(c, B, sum(h_values))
  # Sign scores give positive weight only to diagonal Walsh averages. Use
  # closed median endpoints, including ties; U and Walsh pairs are unnecessary.
  if (all(h_values[-1L] == 1)) {
    binomial_c <- if (c == -1) B else B-floor(c)-1
    return(cosi_median(W,binomial_c))
  }
  if (B*(B+1)/2 > control$max_walsh) stop("Number of Walsh pairs exceeds control$max_walsh")
  .cosi_walsh_interval(W, .cosi_tie_breakers(B,U), h_values, c, control)
}

# Construct the normal-score signed-rank interval from W and its calibrated
# kappa (passed as c); the normal scores are determined by length(W).
cosi_normality <- function(W, c, U = NULL, control = list()) {
  W <- .cosi_validate_estimators(W)
  scores <- .cosi_score_values(length(W), "normality")
  cosi_symmetry(W, c, h = scores, U = U, control = control)
}

# Construct the unimodal-star interval using W, independent anchor a, and
# c drawn from c*_{B,alpha}; also return the expanded per-estimate L and U.
cosi_unimodal_star <- function(W, a, c) {
  W <- .cosi_validate_estimators(W)
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a)) stop("a must be one finite anchor estimate")
  B <- length(W)
  .cosi_check_cutoff(c, B)
  d <- (1 + sqrt(2))*abs(W-a)
  L <- pmin(W,a)-d
  U <- pmax(W,a)+d
  result <- .cosi_interval(.cosi_order_value(L,B-c), .cosi_order_value(U,c+1))
  result$L <- L
  result$U <- U
  result
}

# Represent real endpoints, endpoint inclusion, and whether the resulting interval is empty.
.cosi_interval <- function(lower, upper, lower_closed = TRUE, upper_closed = TRUE) {
  # Closed/open endpoints matter after projecting lexicographic inequalities
  # onto the real line.  Retain the endpoints even if they describe the empty set.
  if (anyNA(c(lower, upper))) stop("Undefined interval endpoint")
  closed <- c(lower = is.finite(lower) && lower_closed,
              upper = is.finite(upper) && upper_closed)
  empty <- lower > upper || (lower == upper && !all(closed))
  list(CI = c(lower = lower, upper = upper), closed = closed, empty = empty)
}

# Return a one-based order statistic with extended values outside the observed range.
.cosi_order_value <- function(x, k) {
  if (k <= 0) return(-Inf)
  if (k > length(x)) return(Inf)
  sort(x, partial = k)[k]
}

# Invert weighted augmented Walsh quantiles, retaining lexicographic endpoint inclusion.
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

  # Find a weighted lexicographic quantile, optionally its right-hand limit.
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

# Test membership in a constructor or main-function result, respecting endpoints.
contains_cosi_ci <- function(object, t) {
  if (!is.list(object) || !all(c("CI","closed","empty") %in% names(object))) {
    stop("object must be a COSI interval result")
  }
  if (!is.numeric(t)) stop("t must be numeric")
  if (object$empty) return(rep(FALSE,length(t)))
  low <- if (object$closed[1L]) t >= object$CI[1L] else t > object$CI[1L]
  high <- if (object$closed[2L]) t <= object$CI[2L] else t < object$CI[2L]
  low & high & is.finite(t)
}
