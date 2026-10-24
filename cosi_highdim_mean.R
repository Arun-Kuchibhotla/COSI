# Core high-dimensional COSI confidence sets.
# Source this file to define cosi_hd(). No simulation runs when sourced.
# A call such as fit <- cosi_hd(X) uses approximately half of X for the pilot.
# pilot='empirical' uses the sample mean and ordinary sample covariance.
# pilot='mom' uses coordinatewise MOM and entrywise MOM covariance estimates.
# pilot='both' returns both choices using IDENTICAL splits and cutoff draws.
#
# fit$median$contains(theta): membership of the first (nonrectangular) set.
# fit$mp$lower, fit$mp$upper: coordinate bounds of the second set.
# fit$median$volume_root and fit$mp$volume_root: Vol(C)^(1/d).
# fit$Sigma_hat(rows,cols): requested block of the pilot covariance estimate.
# fit$Sigma_hat(full=TRUE): materialize the full covariance matrix if desired.


# Reuse the univariate calibration functions; source their file only if needed.
if (!exists("smallest_B", mode="function") ||
    !exists("randomization_quantile_dist", mode="function"))
  source("cosi_univariate.R", local=TRUE)

# Smallest odd K with scalar MOM bad-block bound <=alpha, independent of d.
# Chebyshev gives bad-block probability <=1/8 at threshold
# sqrt(8)*sigma_k/sqrt(floor(n_pilot/K)). At alpha=.05 this gives K=3.
# Uniform MOM bounds can put dimension in the expected supremum rather than K;
# see Minsker, arXiv:1812.03523, Eq. (2.7). Its constants are unspecified, so
# this scalar calibration is not a simultaneous 1-alpha pilot accuracy claim.
# Conditional COSI coverage permits any independent pilot and positive scales.
.hd_mom_groups <- function(alpha=.05,max_groups=Inf) {
  K <- 1L
  while (pbinom(K%/%2L,K,1/8,lower.tail=FALSE)>alpha) K <- K+2L
  if (K>max_groups) stop("Pilot is too small for MOM tuning: need ",
    "at least ",2*K," pilot observations for location and paired-square scales.")
  K
}

.hd_group_means <- function(X,group,K) rowsum(X,group,reorder=FALSE)/tabulate(group,K)
.hd_col_median <- function(X) {
  if (requireNamespace("matrixStats",quietly=TRUE)) matrixStats::colMedians(X)
  else apply(X,2L,stats::median)
}

.hd_cutoff <- function(B,alpha,randomize=TRUE,one_sided=FALSE) {
  # The unchanged helper divides its input alpha by TWO internally.
  z <- randomization_quantile_dist(B,alpha=if(one_sided) 2*alpha else alpha)
  z$q_B-as.integer(randomize && runif(1)>z$p)
}

.hd_order_row <- function(X,k) {
  if (k<=0L) return(rep(-Inf,ncol(X)))
  if (k>nrow(X)) return(rep(Inf,ncol(X)))
  # At the default B=9 this vectorized insertion network is inexpensive.
  Y <- X
  if (nrow(Y)>1L) for (i in 2:nrow(Y)) for (j in seq.int(i-1L,1L)) {
    lo <- pmin(Y[j,],Y[j+1L,]); hi <- pmax(Y[j,],Y[j+1L,])
    Y[j,] <- lo; Y[j+1L,] <- hi
  }
  unname(Y[k,])
}

# Log[(hi/scale)^d - ((hi-width)/scale)^d]. Using the original
# band width avoids cancellation when a narrow interval is far from the pilot.
.hd_log_cone_band <- function(hi, width, scale, d) {
  out <- rep(-Inf, length(hi))
  keep <- hi > 0 & width > 0
  if (!any(keep)) return(out)
  h <- hi[keep]; w <- width[keep]
  fraction <- pmin(1, w / h)
  correction <- log(-expm1(d * log1p(-fraction)))
  # Only reachable for extreme floating-point scales: use 1-(1-u)^d ~ d*u.
  tiny <- fraction == 0
  correction[tiny] <- log(d) + log(w[tiny]) - log(h[tiny])
  out[keep] <- d * (log(h) - log(scale[keep])) + correction
  out
}

.hd_log_volume_selected <- function(lower, upper, pilot, scale) {
  d <- length(lower)
  if (!d || length(upper) != d || length(pilot) != d || length(scale) != d)
    stop("Bounds, pilot, and scales must have the same positive dimension.")
  if (anyNA(lower) || anyNA(upper) || any(!is.finite(pilot)) ||
      any(!is.finite(scale)) || any(scale <= 0))
    stop("Bounds cannot be missing; pilots and positive scales must be finite.")
  nonempty_band <- lower < upper
  if (!any(nonempty_band)) return(-Inf)
  if (any(nonempty_band & (!is.finite(lower) | !is.finite(upper))))
    return(Inf)
  if (d == 1L) return(log(upper - lower))

  # Each signed selected-coordinate cone contributes
  # prod(scale) * 2^(d-1)/d * (hi_standardized^d - lo_standardized^d).
  lp <- ln <- rep(-Inf, d)
  positive <- nonempty_band & upper > pilot
  if (any(positive)) {
    hi <- upper[positive] - pilot[positive]
    width <- upper[positive] - pmax(lower[positive], pilot[positive])
    lp[positive] <- .hd_log_cone_band(hi, width, scale[positive], d)
  }
  negative <- nonempty_band & lower < pilot
  if (any(negative)) {
    hi <- pilot[negative] - lower[negative]
    width <- pmin(upper[negative], pilot[negative]) - lower[negative]
    ln[negative] <- .hd_log_cone_band(hi, width, scale[negative], d)
  }
  terms <- c(lp, ln)
  largest <- max(terms)
  if (is.infinite(largest)) return(largest)
  log_sum <- largest + log(sum(exp(terms - largest)))
  sum(log(scale)) + (d - 1) * log(2) - log(d) + log_sum
}


# alpha tunes the MOM group count without any dependence on dimension.
# The covariance entries use the SAME K and alpha input. They are MOMs of
# paired outer products, which have expectation Sigma. A covariance error bound
# requires additional product/fourth moments; finite covariance alone does not
# provide one. The entrywise MOM estimate need not be positive semidefinite.

.hd_membership_function <- function(lower,upper,pilot=NULL,scale=NULL) {
  force(lower); force(upper); force(pilot); force(scale)
  function(theta) {
    if (!is.numeric(theta) || length(theta)!=length(lower) || any(!is.finite(theta)))
      stop("theta must be a finite vector with one entry per coordinate.")
    if (is.null(pilot)) return(all(lower<=theta & theta<=upper))
    v <- (pilot-theta)/scale
    signed_k <- which.min(c(v,-v))
    k <- (signed_k-1L) %% length(theta)+1L
    lower[k]<=theta[k] && theta[k]<=upper[k]
  }
}

# The returned function stores only pilot residuals or paired differences,
# not the complete data or any d-by-d array. Covariance blocks are evaluated
# deterministically and never consume random numbers.
.hd_covariance_function <- function(basis,method,groups=NULL,K=NULL) {
  force(basis); force(method); force(groups); force(K)
  d <- ncol(basis)
  function(rows=NULL,cols=NULL,full=FALSE) {
    if (is.null(rows)) {
      if (!isTRUE(full)) stop("Supply rows and columns, or use full=TRUE.")
      rows <- seq_len(d)
    }
    if (is.null(cols)) cols <- rows
    valid <- function(i) is.numeric(i) && length(i)>0L &&
      all(is.finite(i) & i==floor(i) & i>=1L & i<=d)
    if (!valid(rows) || !valid(cols)) stop("Covariance indices must be in 1,...,d.")
    if (method=="empirical")
      return(crossprod(basis[,rows,drop=FALSE],basis[,cols,drop=FALSE])/(nrow(basis)-1L))
    # Chunk BOTH axes so that entrywise MOM never allocates K*d*d numbers.
    ans <- matrix(NA_real_,length(rows),length(cols))
    for (a in seq.int(1L,length(rows),by=128L)) {
      ia <- seq.int(a,min(length(rows),a+127L))
      for (b in seq.int(1L,length(cols),by=128L)) {
        ib <- seq.int(b,min(length(cols),b+127L))
        estimates <- matrix(NA_real_,K,length(ia)*length(ib))
        for (g in seq_len(K)) {
          ii <- which(groups==g)
          estimates[g,] <- as.vector(crossprod(
            basis[ii,rows[ia],drop=FALSE],basis[ii,cols[ib],drop=FALSE])/length(ii))
        }
        ans[ia,ib] <- matrix(.hd_col_median(estimates),length(ia),length(ib))
      }
    }
    ans
  }
}

.hd_estimate_pilot <- function(X0,method,alpha,scale_floor_relative=1e-8) {
  n0 <- nrow(X0); d <- ncol(X0)
  if (method=="empirical") {
    vartheta <- colMeans(X0)
    residuals <- sweep(X0,2L,vartheta,"-")
    variance_hat <- colSums(residuals^2)/(n0-1L)
    Sigma_hat <- .hd_covariance_function(residuals,"empirical")
    K <- NA_integer_
  } else {
    pairs <- n0%/%2L
    K <- .hd_mom_groups(alpha,max_groups=pairs)
    groups <- rep(seq_len(K),length.out=n0)
    vartheta <- .hd_col_median(.hd_group_means(X0,groups,K))
    odd <- seq.int(1L,2L*pairs,by=2L)
    differences <- (X0[odd,,drop=FALSE]-X0[odd+1L,,drop=FALSE])/sqrt(2)
    pair_groups <- rep(seq_len(K),length.out=pairs)
    variance_hat <- .hd_col_median(.hd_group_means(differences^2,pair_groups,K))
    Sigma_hat <- .hd_covariance_function(differences,"mom",pair_groups,K)
  }
  # Sigma_hat and variance_hat retain the raw covariance estimates. A tiny
  # relative floor is used ONLY for the scales entering the confidence sets.
  positive <- variance_hat[is.finite(variance_hat) & variance_hat>0]
  if (!length(positive)) stop("The pilot supplied no positive marginal variances.")
  floor_sd <- scale_floor_relative*sqrt(median(positive))
  regularized <- !is.finite(variance_hat) | variance_hat<floor_sd^2
  variance_for_scale <- variance_hat
  variance_for_scale[regularized] <- floor_sd^2
  list(vartheta=unname(vartheta),variance_hat=unname(variance_hat),
    sigma_hat=sqrt(unname(variance_for_scale)),Sigma_hat=Sigma_hat,
    groups=K,regularized=sum(regularized),scale_floor=floor_sd)
}

.hd_fit_sets <- function(W,pilot,cm,cp,S0,blocks,alpha,delta,method) {
  B <- nrow(W); d <- ncol(W)
  lower <- .hd_order_row(W,B-cm); upper <- .hd_order_row(W,cm+1L)
  log_volume <- .hd_log_volume_selected(lower,upper,pilot$vartheta,pilot$sigma_hat)
  median_set <- list(contains=.hd_membership_function(lower,upper,pilot$vartheta,pilot$sigma_hat),
    volume_root=exp(log_volume/d),log_volume=log_volume,cutoff=cm,
    coordinate_lower=lower,coordinate_upper=upper)
  distances <- sweep(abs(sweep(W,2L,pilot$vartheta,"-")),2L,pilot$sigma_hat,"/")
  # Maximum OF coordinatewise order statistics, exactly as in the definition.
  radius <- max(.hd_order_row(distances,cp+1L))
  lower <- pilot$vartheta-radius*pilot$sigma_hat
  upper <- pilot$vartheta+radius*pilot$sigma_hat
  log_volume <- sum(log(2*radius*pilot$sigma_hat))
  mp_set <- list(lower=lower,upper=upper,
    contains=.hd_membership_function(lower,upper),
    volume_root=exp(log_volume/d),log_volume=log_volume,cutoff=cp,radius=radius)
  list(median=median_set,mp=mp_set,
    volume_root=c("cosi-m"=median_set$volume_root,"cosi-m+"=mp_set$volume_root),
    vartheta=pilot$vartheta,variance_hat=pilot$variance_hat,sigma_hat=pilot$sigma_hat,
    Sigma_hat=pilot$Sigma_hat,pilot_method=method,pilot_groups=pilot$groups,
    scale_regularized=pilot$regularized,scale_floor=pilot$scale_floor,
    B=B,n_pilot=length(S0),S0=S0,blocks=blocks,block_means=W,
    alpha=alpha,delta=delta)
}

cosi_hd <- function(X,alpha=.05,delta=.30,pilot=c("mom","empirical","both"),
    n_pilot=floor(nrow(X)/2),randomize=TRUE,seed=2026) {
  pilot <- match.arg(pilot)
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X) || ncol(X)<1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix, with observations in rows.")
  n <- nrow(X); B <- smallest_B(alpha,delta)
  if (length(n_pilot)!=1L || !is.finite(n_pilot) || n_pilot<2L ||
      n_pilot!=floor(n_pilot) || n-n_pilot<B)
    stop("Need an integer pilot size >=2 and at least B inference observations.")
  # An explicit seed makes this call reproducible WITHOUT changing the caller's
  # random-number stream. seed=NULL uses the caller's stream, for simulations.
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    if (had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
    on.exit(if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv)
      else rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
    set.seed(seed)
  }
  permutation <- sample.int(n)
  S0 <- permutation[seq_len(n_pilot)]
  inference <- permutation[-seq_len(n_pilot)]
  groups <- rep(seq_len(B),length.out=length(inference))
  blocks <- split(inference,groups)
  W <- .hd_group_means(X[inference,,drop=FALSE],groups,B)
  X0 <- X[S0,,drop=FALSE]
  methods <- if(pilot=="both") c("empirical","mom") else pilot
  pilots <- lapply(methods,function(method) .hd_estimate_pilot(X0,method,alpha))
  # One common cutoff realization per set, shared across coordinates AND pilot
  # choices, makes the two pilot comparisons use exactly the same randomization.
  cm <- .hd_cutoff(B,alpha,randomize,one_sided=FALSE)
  cp <- .hd_cutoff(B,alpha,randomize,one_sided=TRUE)
  fits <- lapply(seq_along(methods),function(i)
    .hd_fit_sets(W,pilots[[i]],cm,cp,S0,blocks,alpha,delta,methods[i]))
  names(fits) <- methods
  if(pilot=="both") fits else fits[[1L]]
}
