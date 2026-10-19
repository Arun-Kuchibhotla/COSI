# Find the smallest B for the randomized binomial cutoff in COSI.
# Uses only base R (the stats package).  Source this file, then call
# smallest_B(alpha = 0.05, delta = 1/3).
#
# The (1-\alpha/2)-th quantile of binomial distribution is
# q_B = min{k: P(Binom(B, 1/2) > k) <= alpha/2}.
# c*_{B,\alpha/2} used in COSI intervals is randomized between
# q_B and q_B - 1 so that P(Binom(B, 1/2) \le c^*_{B,\alpha/2}) = 1 - \alpha/2.
# The randomization may degenerate to a single integer at an exact CDF level.
# The following function smallest_B finds the smallest B such that q_B \le B(1/2 + delta), 
# given alpha and delta. 
smallest_B <- function(alpha = 0.05, delta = 1/3) {
  stopifnot(
    is.numeric(alpha), length(alpha) == 1L,
    is.finite(alpha), alpha > 0, alpha < 1,
    is.numeric(delta), length(delta) == 1L,
    is.finite(delta), delta > 0
  )
  if(delta >= 0.5){
  	return(1L)
  }

  B0 <- ceiling(1 - log2(alpha))
  candidates <- c(B0, B0 + 1)  # earliest remaining odd and even candidates
  log_alpha <- log(alpha) - log(2)

  repeat {
    j <- which.min(candidates)
    B <- candidates[j]

    if (!is.finite(B))
      stop("The search requires B beyond the supported range B < 2^53.")

    k <- qbinom(log_alpha, size = B, prob = 0.5,
                lower.tail = FALSE, log.p = TRUE)

    if (!is.finite(k))
      stop("The binomial quantile evaluation failed.")

    excess <- k - B / 2

    # Equivalent to k <= B * (0.5 + delta).
    if (excess <= B * delta)
      return(B)

    # Skip values of the same parity that cannot satisfy the condition.
    # This works because if q_B satisfies q_{B+2} \ge q_B + 1 so that
    # q_B - B/2 is separately monotonic for even and odd B.
    next_B <- B + 2 * max(1, floor((excess / delta - B) / 2))

    candidates[j] <-
      if (is.finite(next_B) && next_B < 2^53) next_B else Inf
  }
}

## The following function computes augmented order
## statistics
augmented_os <- function(x, k){
	c(-Inf, x, Inf)[k + 1L]
}

## The following function finds p
## such that 
## p*P(Binom(B, 1/2) > q_B) + (1-p)P(Binom(B, 1/2) > q_B - 1) = \alpha/2
randomization_quantile_dist <- function(B, alpha = 0.05){
	log_tail <- log(alpha) - log(2)
	q_B <- qbinom(log_tail, size = B, prob = 0.5,
		        lower.tail = FALSE, log.p = TRUE)
	# log_S[1] = log P(Z > q_B)
	# log_S[2] = log P(Z > q_B - 1)
	log_S <- pbinom(c(q_B, q_B - 1), size = B, prob = 0.5,
	               lower.tail = FALSE, log.p = TRUE)

	if (is.na(log_S[1]) || !is.finite(log_S[2]) ||
	  log_S[1] >= log_S[2])
	stop("Could not evaluate the two adjacent binomial tails.")

	# Allow only rounding-scale disagreement at a CDF jump.
	tol <- 64 * .Machine$double.eps * max(1, abs(log_tail), abs(log_S[2]))

	if (log_S[1] > log_tail + tol || log_S[2] < log_tail - tol)
	stop("The binomial quantile and tail evaluations disagree.")

	# Stable evaluation of
	# p = (P(Z > q_B - 1) - alpha/2) / P(Z = q_B).
	p <- expm1(log_tail - log_S[2]) / expm1(log_S[1] - log_S[2])
	p <- min(1, max(0, p))  # limit final rounding to [0, 1]

	return(list(B = B, q_B = q_B, p = p))
}


## The following function cosi_median computes the generalized HulC
## confidence intervals. These are also the COSI intervals when
## the scale invariance condition is zero median.
## It takes as input a vector of B real-numbers and the
## quantile value c that lies between 0 and B. If c is not
## given, then it is computed from B and alpha. Return a vector
## containing the lower confidence bound and upper confidence bound.
cosi_median <- function(w, c = NULL, alpha = 0.05, seed = NULL){
	w <- sort(as.vector(w))
	B <- length(w)
	if(!is.null(seed)){
		set.seed(seed)
	}
	if(!is.null(c)){
		return(list(CI = c(augmented_os(w, B-c), augmented_os(w, c+1))))
	} else {
		U <- runif(1)
		if(B == ceiling(1 - log2(alpha))){
			if(U <= 1 + (1 - exp(log(alpha) + (B-1)*log(2)))/B){
				return(list(CI = c(w[1], w[B])))
			} else {
				return(list(CI = c(w[2], w[B-1])))
			}
		} else {
			tmp <- randomization_quantile_dist(B, alpha = alpha)
			p <- tmp$p
			q_B <- tmp$q_B
			q_star <- q_B - (U > tmp$p)
			return(list(CI = c(augmented_os(w, B-q_star), augmented_os(w, q_star + 1))))
		}
	}
}

## The following function cosi_sym computes the COSI
## intervals when the scale invariance condition is 
## symmetry around zero. It takes as input a vector of B real numbers,
## a vector h of B non-negative non-decreasing numbers, 
## and the quantile value c that lies between -1 and sum(h)/2. If c is not
## given, then it is computed from B and alpha. 
# remotes::install_github(
#   "Arun-Kuchibhotla/Weighted-Bernoulli",
#   subdir = "weightedBernoulli",
#   force = TRUE
# )
cosi_sym <- function(w, h, c = NULL, alpha = 0.05, seed = NULL) {
	if (!is.null(seed)){
		set.seed(seed)
	} 
	B <- length(w)
	stopifnot(B > 0, length(h) == B, all(is.finite(w)), all(is.finite(h)),
				all(h >= 0), all(diff(h) >= 0), alpha > 0, alpha < 1)
	cal <- list(kappa = c, exact = NA)
	if (is.null(c)){
		cal <- weightedBernoulli::qweighted_bernoulli(
				p = alpha, weights = h, prob = 0.5, lower.tail = FALSE,
				randomized = TRUE, two.sided = TRUE, draw = TRUE)
	} 
	c <- unname(cal$kappa)
	u <- runif(B, -1, 1) 
	o <- order(w, u) 
	x <- cbind(w, u)[o, , drop = FALSE]
	ij <- which(upper.tri(matrix(0, B, B), diag = TRUE), arr.ind = TRUE)
	A <- (x[ij[, 1], , drop = FALSE] + x[ij[, 2], , drop = FALSE]) / 2
	weight <- diff(c(0, h))[ij[, 2] - ij[, 1] + 1]
	o <- order(A[, 1], A[, 2])
	A <- rbind(c(-Inf, 0), A[o, , drop = FALSE], c(Inf, 0))
	weight <- c(0, weight[o], 0)
	L <- if (c < 0) A[1, ] else A[c(which(cumsum(weight) > c), nrow(A))[1], ]
	H <- if (c < 0) A[nrow(A), ] else A[tail(c(1, which(rev(cumsum(rev(weight))) > c)), 1), ]
	closed <- c(L[2] < 0, H[2] > 0)
	list(CI = unname(c(L[1], H[1])), closed = unname(closed),
			empty = unname(L[1] > H[1] || (L[1] == H[1] && !all(closed))),
			c = c, exact = cal$exact, calibration = cal)
}


## The following function cosi_normal computes the
## COSI intervals when the scale invariance condition is
## zero mean normality. It takes the same inputs as the cosi_sym
## function, except for h.
cosi_normal <- function(w, c = NULL, alpha = 0.05, seed = NULL){
	B <- length(as.vector(w))
	h <- qnorm(1/2 + {1:B}/{2*(B+1)})
	return(cosi_sym(w, h, c = c, alpha = alpha, seed = seed))
}

## The following function cosi_unimodal computed the
## COSI interval when the scale-invariance property is
## unimodal*, i.e., P(\eta*W \in (a, 2a)) \le 1/2 for all a > 0,
## and \eta in {+1, -1}. It takes same inputs as cosi_median with an
## addition of pilot. If pilot = NULL, the function uses ALL B estimates and intersects
## over all real auxiliary choices. It does not choose an observed pilot.
## The rejection cutoff is c = 1 + c*_{B-1, alpha/B}, so each consecutive
## cluster contains c + 1 = 2 + c*_{B-1, alpha/B} observations.
cosi_unimodal <- function(w, c = NULL, alpha = 0.05, seed = NULL, pilot = NULL){
	w <- as.vector(w)
	B <- length(w)
	# Validate before sorting, which otherwise removes missing values.
  stopifnot(
      is.numeric(w), !is.complex(w), B >= 1L,
      all(is.finite(w)),
      is.numeric(alpha), length(alpha) == 1L,
      is.finite(alpha), alpha > 0, alpha < 1
  )
  w <- as.numeric(w)
	if(!is.null(seed)){
		set.seed(seed)
	}
	U <- runif(1)
	if(!is.null(pilot)){
		tmp <- abs(w - pilot)
		lower <- sort(pmin(w, pilot) - (1 + sqrt(2)) * gap)
    upper <- sort(pmax(w, pilot) + (1 + sqrt(2)) * gap)
		if(!is.null(c)){
			return(list(CI = c(augmented_os(lower, B-c), augmented_os(upper, c+1))))
		} else {
			tmp <- randomization_quantile_dist(B, alpha = 2*alpha)
			p <- tmp$p
			q_B <- tmp$q_B
			q_star <- q_B - (U > tmp$p)
			return(list(CI = c(augmented_os(lower, B-q_star), augmented_os(upper, q_star+1))))
		}
	} else {
		stopifnot(is.numeric(w), B >= 1L, all(is.finite(w)),
			length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
		tmp <- randomization_quantile_dist(B - 1L, alpha = 2*alpha/B)
		c_star <- tmp$q_B - (U > tmp$p)
		c <- 1L + c_star
		# One cutoff realization is shared across the entire intersection.
		# c = B means there are no clusters and the interval is R.
		if(c == B) return(list(CI = c(-Inf, Inf), c = c))
		w <- sort(w)
		i <- seq_len(B - c)
		gap <- w[i + c] - w[i]
		return(list(CI = c(max(w[i] - gap,
			min(w[i+c] + gap)), c = c))
	}
}