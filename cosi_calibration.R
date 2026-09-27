# Critical-value distributions and conservative concentration calibration.
# This module never splits data or constructs confidence intervals.
# Source cosi_smallest_B.R before this file, or load everything via cosi_ci.R.

# Evaluate the two-sided binomial rejection probability without tail underflow.
.cosi_sign_rejection_probability <- function(B, k) {
  if (k < 0) return(0)
  if (k >= B/2) return(1)
  if (k == 0) return(2^(1-B))
  min(1,exp(log(2)+pbinom(k,B,.5,log.p=TRUE)))
}

# Compute the sign-score kappa law using kappa=B-c_star-1; keep the binomial
# law so identical random inputs select the same interval as the median method.
.cosi_sign_calibration <- function(B, alpha) {
  bin <- .cosi_cutoff_probabilities(B,alpha)
  bin_support <- c(bin$q-1,bin$q)
  keep <- bin$probabilities > 0
  bin_law <- list(support=bin_support[keep],probabilities=bin$probabilities[keep])
  cm <- B-bin$q-1
  cp <- B-bin$q
  pm <- .cosi_sign_rejection_probability(B,cm)
  pp <- .cosi_sign_rejection_probability(B,cp)
  probabilities <- rev(bin$probabilities)
  keep <- probabilities > 0
  support <- c(cm,cp)[keep]
  rejection <- c(pm,pp)[keep]
  level_exact <- cp != B/2 || bin$probabilities[1L] == 0
  tau <- bin$probabilities[1L]
  # At an exact CDF level the singleton is c_plus, not c_minus.
  if (bin$probabilities[1L] == 0) {
    cp <- cm
    cm <- cp-1
    pp <- pm
    pm <- .cosi_sign_rejection_probability(B,cm)
    tau <- 1
  }
  law <- list(support=support,probabilities=probabilities[keep],
    rejection_probabilities=rejection,c_minus=cm,c_plus=cp,
    P_minus=pm,P_plus=pp,tau=tau,engine="sign-binomial",exact=TRUE,
    median_compatible=TRUE,nominal_level_exact=level_exact,
    rejection_probability_bound=sum(rejection*probabilities[keep]))
  list(kappa_law=law,binomial_law=bin_law)
}

# Calibrate and draw kappa for h=c(h(1),...,h(B)), excluding h(0)=0.
# Return kappa and its law. Seed the call locally; u may supply the cutoff draw.
# Unit scores use median-compatible randomization, conservative at some central
# cutoffs. median_compatible=FALSE selects the canonical P_h mixture instead.
cosi_kappa <- function(h, alpha = 0.05, seed = NULL, randomized = TRUE,
                       control = list(), u = NULL, median_compatible = TRUE) {
  restore <- .cosi_seed_scope(seed)
  on.exit(restore(),add=TRUE)
  if (!is.numeric(h) || !is.null(dim(h)) || !length(h) || any(!is.finite(h)) ||
      any(h < 0) || any(diff(h) < 0) || !is.finite(sum(h))) {
    stop("h must be the finite, nonnegative, nondecreasing vector h(1),...,h(B), without h(0)")
  }
  .cosi_check_alpha(alpha)
  if (!is.logical(randomized) || length(randomized) != 1 || is.na(randomized)) {
    stop("randomized must be TRUE or FALSE")
  }
  if (!is.logical(median_compatible) || length(median_compatible) != 1 || is.na(median_compatible)) {
    stop("median_compatible must be TRUE or FALSE")
  }
  if (!is.null(u) && (!is.numeric(u) || length(u) != 1 || !is.finite(u) || u < 0 || u >= 1)) {
    stop("u must lie in [0,1)")
  }
  control <- .cosi_ci_control(control)
  if (median_compatible && all(h == 1)) {
    calibrated <- .cosi_sign_calibration(length(h),alpha)
    law <- calibrated$kappa_law
    c_star <- .cosi_draw_law(calibrated$binomial_law,randomized,"max",u)
    kappa <- length(h)-c_star-1
  } else {
    law <- .cosi_sr_calibration(h,alpha,control)
    kappa <- .cosi_draw_law(law,randomized,"min",u)
  }
  list(kappa=kappa,distribution=law)
}

# Seed all R randomness for this call and return a function restoring the
# caller's RNG state, including when an estimator or validation step errors.
.cosi_seed_scope <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be NULL or an integer in [0,2147483647]")
  }
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_seed <- if (had_seed) get(".Random.seed",envir=.GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv)
    else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) {
      rm(".Random.seed",envir=.GlobalEnv)
    }
    invisible(NULL)
  }
}

# Fill in and validate resource limits and numerical tolerances.
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

# Draw a cutoff from its law, or select the requested conservative support value.
.cosi_draw_law <- function(law, randomized, conservative = c("max", "min"), u = NULL) {
  conservative <- match.arg(conservative)
  if (!randomized) {
    return(if (conservative == "max") max(law$support) else min(law$support))
  }
  if (length(law$support) == 1L) return(law$support)
  if (!is.null(u)) {
    cumulative <- cumsum(law$probabilities)
    cumulative[length(cumulative)] <- 1
    return(law$support[which(u < cumulative)[1L]])
  }
  law$support[sample.int(length(law$support), 1L, prob = law$probabilities)]
}

# Compute the two-point signed-rank cutoff law achieving the requested rejection level.
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

# Compute log(cosh(x)) stably for small and large arguments.
.cosi_log_cosh <- function(x) {
  # log(cosh(x)) without overflow, or cancellation when x is small.
  x <- abs(x)
  out <- x + log1p(exp(-2*x)) - log(2)
  small <- x < 1
  out[small] <- log1p(2*sinh(x[small]/2)^2)
  out
}

# Compute the log MGF, lambda*psi'(lambda)-psi(lambda), and psi''(lambda).
# Nonnegative weights and stable exponential formulas avoid cosh overflow.
.cosi_rademacher_cgf <- function(weights, lambda) {
  z <- lambda*weights
  e <- exp(-2*z)
  log_mgf <- sum(.cosi_log_cosh(z))
  entropy <- log(2)-log1p(e)-2*z*e/(1+e)
  small <- z < 1
  entropy[small] <- z[small]*tanh(z[small])-.cosi_log_cosh(z[small])
  curvature <- sum(weights^2*4*e/(1+e)^2)
  list(log_mgf=log_mgf, entropy=sum(entropy), curvature=curvature)
}

# Compute a conservative Chernoff radius using an analytical starting value
# and three fixed refinements when the optimizer is finite; no root search.
# gap_bound controls the error relative to the optimal Chernoff radius in
# exact arithmetic, with an extra floating-point allowance in the result.
.cosi_finite_chernoff_radius <- function(weights, log_level) {
  variance <- sum(weights^2)
  m <- sum(weights > 0)
  L <- log(2)-log_level
  lambda <- sqrt(2*L/variance)
  finite_optimum <- L < m*log(2)
  if (finite_optimum) {
    t <- tan(L/(m*sqrt(3)))
    lambda <- sqrt((m/variance)*(2*sqrt(3)*t)/(1-sqrt(3)*t))
  }
  candidates <- numeric(0)
  for (j in seq_len(if (finite_optimum) 4L else 1L)) {
    values <- .cosi_rademacher_cgf(weights, lambda)
    K <- values$log_mgf
    H <- values$entropy
    curvature <- values$curvature
    candidates <- c(candidates,lambda)
    radius <- (L+K)/lambda
    gap <- if (finite_optimum) max(0,L-H)/lambda+
      1024*.Machine$double.eps*max(1,radius) else Inf
    if (j == 4L || !finite_optimum || curvature <= 0) break
    square <- lambda^2+2*max(0,L-H)/curvature
    if (!is.finite(square)) break
    lambda <- sqrt(square)
  }
  list(radius=radius, lambda=tail(candidates,1L), log_mgf=K,
       gap_bound=gap, lambdas=candidates, refinements=length(candidates)-1L)
}

# Invert the binomial positive-part second-moment bound by scanning its cells.
# Return the quantile, minimizing shift, and guarded log second moment.
# Log probabilities and a conditional-variance recurrence avoid underflow
# and cancellation from subtracting nearly equal raw second moments.
.cosi_binomial_p2_quantile <- function(n, p, log_beta) {
  guard <- 256*.Machine$double.eps*(n+1)
  if (log_beta <= n*log(p)+guard) return(NULL)
  log_pmf <- dbinom(0:n,n,p,log=TRUE)
  tail_log <- -Inf
  delta <- variance <- 0
  best <- NULL
  for (k in seq.int(n,0)) {
    atom_log <- log_pmf[k+1L]
    if (k == n) {
      tail_log <- atom_log
    } else {
      combined <- max(atom_log,tail_log)+log1p(exp(-abs(atom_log-tail_log)))
      tail_weight <- exp(tail_log-combined)
      atom_weight <- exp(atom_log-combined)
      previous_mean <- delta+1
      variance <- tail_weight*variance+tail_weight*atom_weight*previous_mean^2
      delta <- tail_weight*previous_mean
      tail_log <- combined
    }
    # delta is E[Y-k | Y >= k]; variance is Var(Y | Y >= k).
    log_ratio <- tail_log+guard-log_beta
    if (log_ratio > 0 && variance > 0) {
      log_denominator <- if (log_ratio < 50) log(expm1(log_ratio)) else
        log_ratio+log1p(-exp(-log_ratio))
      log_gap <- (log(variance)-log_denominator)/2
      if (k > 0 && log_gap >= log(delta+1)) {
        offset <- -1
      } else {
        offset <- min(0,delta-exp(log_gap))
      }
    } else if (k > 0) {
      offset <- -1
    } else next
    moment <- (delta-offset)^2+variance
    if (moment <= 0 || !is.finite(moment)) next
    log_moment <- tail_log+log(moment)+guard
    log_radius <- (log_moment-log_beta)/2
    if (log_radius >= log(.Machine$double.xmax)) next
    shift <- k+offset
    quantile <- shift+exp(log_radius)
    if (is.null(best) || quantile < best$quantile) {
      best <- list(quantile=quantile,shift=shift,log_moment=log_moment)
    }
  }
  best
}

# Construct the matching-variance binomial comparison and two-sided radius.
# Retain the shift and moment to check the final implemented rejection bound.
.cosi_binomial_p2_radius <- function(weights, log_level) {
  weights <- weights[weights > 0]
  m <- length(weights)
  variance <- sum(weights^2)
  y <- max(weights)
  ratio <- variance/(m*y^2)
  p <- ratio/(1+ratio)
  step <- y*(1+ratio)
  inverse <- .cosi_binomial_p2_quantile(m,p,log_level-log(2))
  if (is.null(inverse)) return(list(radius=Inf,shift=0,log_moment=Inf))
  list(radius=step*(inverse$quantile-m*p),shift=step*(inverse$shift-m*p),
       log_moment=2*log(step)+inverse$log_moment)
}

# Calibrate a deterministic conservative cutoff from analytical concentration
# radii, three fixed Chernoff refinements, and a binomial P2 comparison.
# Preserve existing output fields and add diagnostics for the finite checks.
.cosi_concentration_calibration <- function(scores, alpha, control,
                                             reason="general real scores") {
  N <- sum(scores)
  if (N == 0) return(.cosi_kappa_law(-1,0,0,1,alpha,"zero-scores"))
  mu <- N/2
  w <- scores[scores > 0]/N
  variance <- sum(w^2)
  sd <- sqrt(variance)
  log_alpha <- log(alpha)
  L <- log(2)-log_alpha
  chernoff <- .cosi_finite_chernoff_radius(w,log_alpha)
  comparison <- .cosi_binomial_p2_radius(w,log_alpha)

  # Sharp constant-factor Gaussian comparison, not a normal approximation.
  log_C <- -log(4)-pnorm(sqrt(2),lower.tail=FALSE,log.p=TRUE)
  r_gaussian <- sd*qnorm(log_alpha-log(2)-log_C,lower.tail=FALSE,log.p=TRUE)
  r_hoeffding <- sd*sqrt(2*L)
  r_chebyshev <- sd*exp(-log_alpha/2)
  if (alpha == .5) r_gaussian <- r_chebyshev <- sd*sqrt(2)

  # Explicit conservative inversion of Pinelis's Q bound via Mills's ratio.
  C_P <- 5*sqrt(2*pi*exp(1))*(2*pnorm(1)-1)
  # A lower bound for z avoids cancellation as alpha approaches one.
  z <- if (alpha > .9) sqrt(pi/2)*(1-alpha) else
    qnorm(log_alpha-log(2),lower.tail=FALSE,log.p=TRUE)
  A <- 1+C_P*(z+1/z)/(9+z^2)
  r_pinelis <- sd*qnorm(log_alpha-log(2)-log(A),lower.tail=FALSE,log.p=TRUE)
  # Tomaszewski controls a strict tail at one standard deviation.
  r_tomaszewski <- if (alpha >= .5) sd*(1+512*.Machine$double.eps) else Inf
  radii <- c(chernoff=chernoff$radius,binomial_p2=comparison$radius,
             gaussian=r_gaussian,pinelis=r_pinelis,hoeffding=r_hoeffding,
             chebyshev=r_chebyshev,tomaszewski=r_tomaszewski)
  best <- which.min(radii)
  cushion <- (control$mass_tol+256*.Machine$double.eps)*max(1,N)
  raw_kappa <- mu*(1-radii[best])
  kappa <- if (is.finite(raw_kappa) && raw_kappa > cushion)
    unname(raw_kappa-cushion) else -1
  log_bound <- -Inf
  if (kappa >= 0) {
    effective_kappa <- kappa+control$mass_tol*max(1,N)
    r <- (mu-effective_kappa)/mu
    x <- r/sd
    q_parts <- c(pnorm(x,lower.tail=FALSE,log.p=TRUE),
                 log(C_P)-x^2/2-log(2*pi)/2-log(9+x^2))
    log_Q <- max(q_parts)+log1p(exp(-abs(diff(q_parts))))
    bounds <- c(log(2)+chernoff$log_mgf-chernoff$lambda*r,
                log(2)+log_C+pnorm(x,lower.tail=FALSE,log.p=TRUE),
                log(2)+log_Q,log(2)-r^2/(2*variance),log(variance)-2*log(r))
    if (is.finite(comparison$log_moment) && r > comparison$shift) {
      bounds <- c(bounds,log(2)+comparison$log_moment-2*log(r-comparison$shift))
    }
    if (r > sd) bounds <- c(bounds,-log(2))
    log_bound <- min(0,bounds)
    if (is.na(log_bound) || log_bound > log_alpha) {
      kappa <- -1
      log_bound <- -Inf
    }
  }
  list(support=kappa,probabilities=1,rejection_probabilities=NA_real_,
       engine="concentration",exact=FALSE,bound_used=names(radii)[best],
       rejection_probability_bound=exp(log_bound),log_rejection_probability_bound=log_bound,
       radius=unname(mu*radii[best]),radius_fraction=unname(radii[best]),
       candidate_radii=radii*mu,numerical_cushion=cushion,fallback_reason=reason,
       lambda=chernoff$lambda,lambda_candidates=chernoff$lambdas,
       chernoff_refinements=chernoff$refinements,chernoff_gap_bound=mu*chernoff$gap_bound)
}

# Convert a discrete score-sum distribution into the randomized signed-rank cutoff law.
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

# Enumerate distinct subset sums and their multiplicities for an exact Bernoulli-sum law.
.cosi_subset_sums <- function(scores) {
  sums <- 0
  for (s in scores) sums <- c(sums, sums + s)
  runs <- rle(sort(sums))
  list(values = runs$values, counts = as.numeric(runs$lengths))
}

# Select signed-rank cutoff support by counting sums of two half-distributions.
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

  # Count weighted pairs of half-sums below a threshold, optionally returning boundary positions.
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

# Select an exact score-law algorithm when feasible, otherwise use concentration calibration.
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
    # Compute the greatest common divisor used to identify an exact score lattice.
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
