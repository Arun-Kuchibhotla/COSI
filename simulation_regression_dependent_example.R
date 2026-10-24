# Four regression examples with serial dependence and gapped COSI blocks.
# Keep this script beside the UNCHANGED cosi_univariate.R and cosi1d.R.
# source("simulation_regression_dependent_example.R") runs the complete study.
# Set RUN_SIMULATION=FALSE at the end to load functions without running it.
# Required packages: weightedBernoulli and ggplot2 (>=3.4.0).
# Competitors come first; then generators, gapped COSI, simulations, and plots.

load_dependent_regression_functions <- function(core_dir=".",envir=.GlobalEnv) {
  for(pkg in c("weightedBernoulli","ggplot2"))
    if(!requireNamespace(pkg,quietly=TRUE)) stop("Install '",pkg,"' first.")
  core <- file.path(core_dir,c("cosi_univariate.R","cosi_univariate(1).R"))
  core <- core[file.exists(core)]
  if(!length(core)) stop("Cannot find cosi_univariate.R.")
  files <- c(core[1],file.path(core_dir,"cosi1d.R"))
  # Load their function definitions without running package installation lines.
  # The supplied core files themselves are not edited.
  for(file in files) {
    if(!file.exists(file)) stop("Cannot find ",file)
    for(expr in parse(file))
      if(is.call(expr) && is.symbol(expr[[1]]) &&
         as.character(expr[[1]]) %in% c("<-","=") && is.call(expr[[3]]) &&
         identical(expr[[3]][[1]],as.name("function"))) eval(expr,envir=envir)
  }
  if(!exists("cosi1D_dependent",envir=envir,mode="function"))
    stop("The supplied cosi1d.R must define cosi1D_dependent().")
  invisible(files)
}

# Dependent-regression competitors. Data columns: response, then covariates.
# Every regression includes an intercept; coefficient=1 selects the first slope.

dep_rank_error <- function() {
  stop(structure(list(message="Rank-deficient regression design.", call=NULL),
    class=c("dep_rank_error", "error", "condition")))
}

dep_ols_coefficient <- function(data, coefficient=1L) {
  Z <- cbind(1, data[, -1L, drop=FALSE])
  fit <- .lm.fit(Z, data[, 1L])
  if (fit$rank < ncol(Z)) dep_rank_error()
  unname(fit$coefficients[coefficient+1L])
}

# The prefix sums let a consecutive window contribute its X'X and X'Y
# without repeatedly copying its rows. Summing these same quantities over
# resampled blocks also gives the OLS fit of the complete bootstrap sample.
dep_regression_workspace <- function(data, coefficient=1L) {
  Z <- cbind(1, data[, -1L, drop=FALSE]); Y <- data[, 1L]
  n <- nrow(Z); q <- ncol(Z)
  if (coefficient < 1L || coefficient >= q)
    stop("coefficient must select one of the supplied slopes.")
  fit <- .lm.fit(Z, Y)
  if (fit$rank < q) dep_rank_error()
  products <- matrix(0, n, q*q)
  for (j in seq_len(q)) for (k in seq_len(q))
    products[, j+(k-1L)*q] <- Z[,j]*Z[,k]
  moments <- cbind(products, Z*Y)
  prefix <- rbind(rep(0, ncol(moments)), apply(moments, 2L, cumsum))
  list(n=n, q=q, coefficient=coefficient, theta_hat=unname(fit$coefficients[coefficient+1L]),
    Z=Z, residuals=fit$residuals, prefix=prefix)
}

dep_coefficient_from_moments <- function(moments, q, coefficient) {
  gram <- matrix(moments[seq_len(q*q)], q, q)
  rhs <- moments[q*q+seq_len(q)]
  upper <- tryCatch(chol(gram), error=function(e) NULL)
  if (is.null(upper)) dep_rank_error()
  beta <- backsolve(upper, forwardsolve(t(upper), rhs))
  unname(beta[coefficient+1L])
}

dep_window_moments <- function(work, starts, size) {
  # Row s of prefix sums the observations preceding start s.
  work$prefix[starts+size, , drop=FALSE] - work$prefix[starts, , drop=FALSE]
}

dep_window_coefficients <- function(work, size, starts) {
  distinct <- unique(starts)
  moments <- dep_window_moments(work, distinct, size)
  values <- vapply(seq_along(distinct), function(i)
    dep_coefficient_from_moments(moments[i,], work$q, work$coefficient), numeric(1))
  values[match(starts, distinct)]
}

# Bartlett-kernel Newey-West covariance, without prewhitening. The finite-
# sample factor agrees with adjust=TRUE: multiply the covariance by n/(n-q).
# The maximum lag is floor(4*(n/100)^(2/9)); kernel bandwidth is lag+1.
# This is a fixed sample-size rule, not automatic plug-in bandwidth selection.
# References: Newey and West (1987), doi:10.2307/1913610;
# statsmodels cov_hac documentation for the fixed lag rule.
dep_wald_hac <- function(data, alpha=.05, coefficient=1L,
    lag=floor(4*(NROW(data)/100)^(2/9)), work=NULL) {
  if (is.null(work)) work <- dep_regression_workspace(data, coefficient)
  n <- work$n
  lag <- min(as.integer(lag), n-1L)
  if (lag < 0L) stop("lag must be nonnegative.")
  bread <- chol2inv(chol(crossprod(work$Z)))
  influence <- drop(work$Z %*% bread[, coefficient+1L])*work$residuals
  variance <- sum(influence^2)
  if (lag > 0L) for (h in seq_len(lag))
    variance <- variance + 2*(1-h/(lag+1))*
      sum(influence[(h+1L):n]*influence[seq_len(n-h)])
  variance <- variance*n/(n-work$q)
  se <- sqrt(max(0, variance))
  list(CI=work$theta_hat+c(-1,1)*qnorm(1-alpha/2)*se,
    standard_error=se, max_lag=lag, kernel_bandwidth=lag+1L)
}

# Basic pairs bootstrap using overlapping, noncircular moving blocks.
# Select ceil(n/ell) block starts with replacement and truncate the last block
# to obtain exactly n rows. Type-6 quantiles use the (R_boot+1)*p convention.
dep_moving_block_bootstrap <- function(data, alpha=.05, R_boot=500L,
    block_length=max(2L, round(NROW(data)^(1/3))), coefficient=1L, work=NULL) {
  if (is.null(work)) work <- dep_regression_workspace(data, coefficient)
  n <- work$n; ell <- as.integer(block_length)
  if (ell < 1L || ell > n) stop("block_length must lie between 1 and n.")
  blocks <- ceiling(n/ell)
  starts <- matrix(sample.int(n-ell+1L, blocks*R_boot, replace=TRUE),
    nrow=R_boot, ncol=blocks)
  moments <- matrix(0, R_boot, ncol(work$prefix))
  for (j in seq_len(blocks)) {
    size <- min(ell, n-(j-1L)*ell)
    moments <- moments + dep_window_moments(work, starts[,j], size)
  }
  estimates <- vapply(seq_len(R_boot), function(r)
    dep_coefficient_from_moments(moments[r,], work$q, coefficient), numeric(1))
  quantiles <- quantile(estimates, c(alpha/2,1-alpha/2), type=6, names=FALSE)
  list(CI=unname(2*work$theta_hat-quantiles[c(2,1)]),
    block_length=ell, R_boot=R_boot)
}

# Estimated-rate subsampling with consecutive time windows. This retains
# Larry Wasserman's 100-point inverse-CDF grid, interquantile spacings, and
# log-size regression. Only the sampling units change: consecutive windows
# replace unordered subsets. Sampling starts with replacement gives N draws
# from the empirical distribution of all overlapping windows of each size.
dep_subsampling_regression <- function(data, alpha=.05,
    b=round(NROW(data)^(2/3)), N=500L, beta=(9:1)/10,
    min_rows=3L*NCOL(data), coefficient=1L, work=NULL) {
  if (is.null(work)) work <- dep_regression_workspace(data, coefficient)
  n <- work$n
  sizes <- round(n^beta)
  keep <- sizes >= min_rows & sizes < n & !duplicated(sizes)
  sizes <- sizes[keep]
  if (length(sizes) < 2L || b < min_rows || b >= n)
    stop("Need two feasible rate sizes and min_rows <= b < n.")

  subsample_quantiles <- function(size, probabilities, gamma=0) {
    starts <- sample.int(n-size+1L, N, replace=TRUE)
    estimates <- dep_window_coefficients(work, size, starts)
    root <- size^gamma*(estimates-work$theta_hat)
    grid <- seq(min(root), max(root), length.out=100L)
    G <- vapply(grid, function(z) mean(root <= z), numeric(1))
    vapply(probabilities, function(p) min(grid[G >= p]), numeric(1))
  }

  probabilities <- c(.1,.6,.2,.7,.3,.8,.4,.9)
  y <- vapply(sizes, function(size) {
    q <- subsample_quantiles(size, probabilities)
    mean(log(q[c(2,4,6,8)]-q[c(1,3,5,7)]))
  }, numeric(1))
  log_sizes <- log(sizes)
  gamma_hat <- -sum((y-mean(y))*(log_sizes-mean(log_sizes)))/
    sum((log_sizes-mean(log_sizes))^2)
  if (is.nan(gamma_hat)) gamma_hat <- 0  # Original fallback, without clipping.
  q <- subsample_quantiles(b, c(alpha/2,1-alpha/2), gamma=gamma_hat)
  list(CI=unname(work$theta_hat-q[c(2,1)]/n^gamma_hat),
    gamma_hat=unname(gamma_hat), b=b, rate_sizes=sizes, N=N)
}

COSI_METHODS <- c("median","sym","normal","unimodal")
DEPENDENT_SCENARIOS <- c("1. Correlated Gaussian design", "2. Skewed design and errors",
  "3. Nonlinear mean, Student t(3) errors", "4. Nonlinear mean, skewed Pareto errors")
DEPENDENT_LABELS <- c(HAC="Wald (HAC)",
  Bootstrap="Moving-block bootstrap",Subsampling="Subsampling (estimated rate)",
  median="COSI median",sym="COSI symmetric",normal="COSI normal scores",
  unimodal="COSI unimodal*")
DEPENDENT_COLORS <- c(HAC="#59636E",Bootstrap="#E69F00",
  Subsampling="#CC79A7",median="#0072B2",sym="#009E73",normal="#D55E00",unimodal="#7B2CBF")

# A positive gap preserves time order in the supplied wrapper. For fixed B,
# this sequence tends to infinity while discarding an asymptotically vanishing
# fraction of rows. At our n grid the gaps are 2,2,2,3,3,3.
dependent_gap <- function(n) max(1L,floor(log(n)/2))
dependent_regression_truth <- function(idx) c(1,1,1.75,1.75)[idx]

# Exactly stationary Gaussian AR(1), with marginal covariance Sigma.
# The first draw already has the stationary law, so no burn-in is needed.
stationary_gaussian_ar <- function(n,p=1L,phi=.6,Sigma=diag(p)) {
  innovations <- matrix(rnorm(n*p),n,p) %*% chol(Sigma)
  out <- innovations
  if(n>1L) for(i in 2:n)
    out[i,] <- phi*out[i-1L,]+sqrt(1-phi^2)*innovations[i,]
  out
}

# Response is first, followed by five covariates. Independent latent G and H
# paths each have AR coefficient .6. G has Toeplitz cross-coordinate covariance
# rho^|j-k|. Gaussian-copula transforms give exact stationary non-Gaussian
# margins while preserving geometric mixing of the latent process.
# Cases3/4 retain the preceding nonlinear mean. Their first projection slope
# is1.75, and their population intercept is1.5+.5*rho (1.85 at rho=.7).
# Pareto shape2.5, rather than1.5, keeps the squared-error objective finite.
gen_data_regression_dependent <- function(n,idx,p=5L,rho=.7,phi_x=.6,phi_e=.6) {
  if(!idx %in% 1:4 || p<2L || abs(rho)>=1 || abs(phi_x)>=1 || abs(phi_e)>=1)
    stop("Use idx=1:4, p>=2 and correlations strictly between -1 and 1.")
  Sigma <- rho^abs(outer(seq_len(p),seq_len(p),"-"))
  G <- stationary_gaussian_ar(n,p,phi_x,Sigma)
  H <- drop(stationary_gaussian_ar(n,1L,phi_e))
  X <- if(idx==2L) exp(G/2)-exp(1/8) else G
  beta <- rep(c(1,-1,.5,.5,-.5),length.out=p)
  mu <- 1+drop(X %*% beta)
  if(idx %in% 3:4) mu <- mu+.5*X[,1]^2+.25*X[,1]^3+.5*X[,1]*X[,2]
  sigma <- if(idx==1L) rep(1,n) else .5+.5*abs(X[,1])
  log_tail <- pnorm(H,lower.tail=FALSE,log.p=TRUE)
  error <- switch(as.character(idx),
    `1`=H,
    `2`=-log_tail-1,
    `3`=qt(log_tail,df=3,lower.tail=FALSE,log.p=TRUE)/sqrt(3),
    `4`={a <- 2.5; sqrt(a*(a-2))*((a-1)/a*exp(-log_tail/a)-1)})
  colnames(X) <- paste0("x",seq_len(p))
  cbind(y=mu+sigma*error,X)
}

# A failed fit is retained as an unbounded interval; do not redraw its data.
dependent_safe <- function(expr) tryCatch(expr,dep_rank_error=function(e)
  list(CI=c(-Inf,Inf),status="rank_deficient",gamma_hat=NA_real_))

dependent_record <- function(ans, truth) {
  ci <- ans$CI
  closed <- if (is.null(ans$closed)) c(TRUE,TRUE) else ans$closed
  empty <- isTRUE(ans$empty) || ci[1]>ci[2] || (ci[1]==ci[2] && !all(closed))
  covered <- !empty && (truth>ci[1] || (truth==ci[1] && closed[1])) &&
    (truth<ci[2] || (truth==ci[2] && closed[2]))
  data.frame(lower=ci[1], upper=ci[2], width=if(empty) 0 else ci[2]-ci[1],
    covered=covered, closed_lower=closed[1], closed_upper=closed[2], empty=empty,
    cutoff=if(is.null(ans[["c"]])) NA_real_ else ans[["c"]],
    calibration_exact=if(is.null(ans$calibration_exact)) NA else ans$calibration_exact,
    status=if(is.null(ans$status)) "ok" else ans$status)
}

# Call the existing gap-based wrapper for every COSI method. A shared seed
# gives identical contiguous inference blocks to median/sym/normal. Unimodal*
# uses one additional separated block for its pilot, exactly as before.
dependent_cosi_fits <- function(data,alpha=.05,delta=.30,gap=dependent_gap(NROW(data)),
    randomize=TRUE,split_seed=NULL,coefficient=1L) {
  estimate <- function(d) dep_ols_coefficient(d,coefficient)
  fits <- lapply(COSI_METHODS,function(method)
    dependent_safe(cosi1D_dependent(data,estimate,alpha=alpha,delta=delta,
      gap=gap,randomize=randomize,seed=split_seed,method=method)))
  names(fits) <- COSI_METHODS
  for(method in COSI_METHODS) if(!is.null(fits[[method]]$status)) {
    B_ci <- smallest_B(alpha,delta)
    B <- B_ci+as.integer(method=="unimodal")
    fits[[method]]$B <- B
    fits[[method]]$B_ci <- B_ci
    fits[[method]]$n_used <- NROW(data)-(B-1L)*gap
  }
  if(all(vapply(fits[1:3],function(z) is.null(z$status),logical(1))))
    stopifnot(identical(fits$median$blocks,fits$sym$blocks),
      identical(fits$median$blocks,fits$normal$blocks))
  fits
}

dependent_regression_design <- function(n_grid=c(100,200,300,500,750,1000),
    alpha=.05,delta=.30,p=5L,gap_function=dependent_gap) {
  B <- smallest_B(alpha,delta)
  n_grid <- sort(unique(as.integer(n_grid)))
  gap <- vapply(n_grid,gap_function,numeric(1))
  if(any(gap<1 | gap!=floor(gap))) stop("Use positive integer gaps to preserve time order.")
  used <- n_grid-(B-1L)*gap
  used_u <- n_grid-B*gap
  minimum <- p+3L
  if(any(floor(used_u/(B+1L))<minimum))
    stop("Too few rows per regression block after gaps; increase n or delta.")
  data.frame(n=n_grid,delta=delta,gap=gap,B=B,unimodal_B=B+1L,
    n_used=used,unimodal_n_used=used_u,smallest_block=floor(used/B),
    unimodal_smallest_block=floor(used_u/(B+1L)),minimum_required=minimum)
}

simulate_regression_dependent <- function(R=500,n_grid=c(100,200,300,500,750,1000),
    alpha=.05,delta=.30,randomize=TRUE,seed=100,R_boot=500,N_sub=500,
    p=5L,rho=.7,phi_x=.6,phi_e=.6,gap_function=dependent_gap,
    workers=1L,checkpoint_dir=NULL) {
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  design <- dependent_regression_design(n_grid,alpha,delta,p,gap_function)
  jobs <- expand.grid(idx=1:4,n=design$n)
  if(!is.null(checkpoint_dir)) dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  one_cell <- function(j) {
    idx <- jobs$idx[j]; n <- jobs$n[j]
    cell_seed <- as.integer(seed+1000000*idx+n)
    set.seed(cell_seed)
    truth <- dependent_regression_truth(idx)
    gap <- design$gap[match(n,design$n)]
    out <- vector("list",R)
    for(r in seq_len(R)) {
      data <- gen_data_regression_dependent(n,idx,p,rho,phi_x,phi_e)
      split_seed <- sample.int(.Machine$integer.max,1)
      work <- dep_regression_workspace(data)
      sub <- dependent_safe(dep_subsampling_regression(data,alpha,N=N_sub,work=work))
      baseline <- list(HAC=dep_wald_hac(data,alpha,work=work),
        Bootstrap=dependent_safe(dep_moving_block_bootstrap(data,alpha,R_boot,work=work)),
        Subsampling=sub)
      # Isolate wrapper seed resets from the next simulated dataset.
      simulation_rng <- get(".Random.seed",envir=.GlobalEnv)
      fits <- c(baseline,dependent_cosi_fits(data,alpha,delta,gap,randomize,split_seed))
      assign(".Random.seed",simulation_rng,envir=.GlobalEnv)
      records <- do.call(rbind,lapply(fits,dependent_record,truth=truth))
      out[[r]] <- cbind(data.frame(idx=idx,scenario=DEPENDENT_SCENARIOS[idx],
        n=n,replicate=r,method=names(fits),truth=truth,alpha=alpha,delta=delta,
        n_covariates=p,rho=rho,phi_x=phi_x,phi_e=phi_e,gap=gap,
        B=vapply(fits,function(z) if(is.null(z$B)) NA_real_ else z$B,numeric(1)),
        B_ci=vapply(fits,function(z) if(is.null(z$B_ci)) NA_real_ else z$B_ci,numeric(1)),
        n_used=vapply(fits,function(z) if(is.null(z$n_used)) n else z$n_used,numeric(1)),
        pilot=vapply(fits,function(z) if(is.null(z$pilot)) NA_real_ else z$pilot,numeric(1)),
        gamma_hat=c(NA,NA,sub$gamma_hat,rep(NA_real_,4)),
        subsample_size=c(NA,NA,round(n^(2/3)),rep(NA_real_,4)),
        bootstrap_block_length=c(NA,round(n^(1/3)),NA,rep(NA_real_,4)),
        hac_max_lag=c(baseline$HAC$max_lag,rep(NA_real_,6)),
        hac_bandwidth=c(baseline$HAC$kernel_bandwidth,rep(NA_real_,6)),
        split_seed=split_seed,cell_seed=cell_seed),records)
      if(r%%100L==0L) message("Scenario ",idx,", n=",n,": ",r,"/",R)
    }
    ans <- do.call(rbind,out); rownames(ans) <- NULL
    reference <- ans$width[ans$method=="HAC"]
    stopifnot(length(reference)==R,all(is.finite(reference) & reference>0))
    ans$wald_hac_width <- rep(reference,each=length(DEPENDENT_LABELS))
    ans$width_ratio <- ans$width/ans$wald_hac_width
    if(anyNA(ans$width_ratio)) stop("Undefined paired Wald width ratio.")
    if(!is.null(checkpoint_dir)) saveRDS(ans,file.path(checkpoint_dir,
      paste0("scenario_",idx,"_n_",n,".rds")))
    message("Completed scenario ",idx,", n=",n,".")
    ans
  }
  if(workers>1L && .Platform$OS.type!="windows") {
    pieces <- parallel::mclapply(seq_len(nrow(jobs)),one_cell,mc.cores=workers,
      mc.set.seed=FALSE,mc.preschedule=FALSE)
    if(any(vapply(pieces,inherits,logical(1),"try-error"))) stop("A worker failed.")
  } else pieces <- lapply(seq_len(nrow(jobs)),one_cell)
  results <- do.call(rbind,pieces)
  results <- results[order(results$idx,results$n,results$replicate,
    match(results$method,names(DEPENDENT_LABELS))),]; rownames(results) <- NULL
  list(results=results,design=design,settings=list(R=R,n_grid=design$n,alpha=alpha,
    delta=delta,randomize=randomize,seed=seed,R_boot=R_boot,N_sub=N_sub,
    p=p,rho=rho,phi_x=phi_x,phi_e=phi_e,gap_rule="max(1,floor(log(n)/2))",
    width_reference="HAC",rng_kind=RNGkind()))
}

read_dependent_results <- function(csv_file,scenarios=NULL) {
  d <- read.csv(csv_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!is.null(scenarios)) d <- d[d$idx %in% scenarios,]
  if (!nrow(d)) stop("No rows for the requested scenarios.")
  if (length(unique(d$delta))!=1L || length(unique(d$alpha))!=1L)
    stop("Select one delta and alpha before plotting.")
  stopifnot(!anyNA(d$covered),!anyNA(d$width_ratio))
  d$scenario <- factor(d$scenario,levels=unique(d$scenario[order(d$idx)]))
  d
}

dependent_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, color = "#25313D",
        margin = ggplot2::margin(b = 13)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5EAF0", linewidth = .35),
      panel.border = ggplot2::element_rect(fill = NA, color = "#D9E0E7", linewidth = .4),
      panel.spacing = grid::unit(1.1, "lines"),
      strip.background = ggplot2::element_rect(fill = "#F1F4F7", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 12, color = "#25313D",
        margin = ggplot2::margin(t = 7, b = 7)),
      axis.text = ggplot2::element_text(color = "#42505D", size = 11),
      axis.title = ggplot2::element_text(color = "#25313D", size = 12),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.size = grid::unit(.42, "cm"), legend.spacing.y = grid::unit(.12, "cm"),
      legend.margin = ggplot2::margin(t = 8), plot.margin = ggplot2::margin(12, 14, 8, 10))
}



plot_coverage_n_dependent <- function(csv_file,scenarios=1:4,ylim=NULL) {
  d <- read_dependent_results(csv_file,scenarios)
  s <- aggregate(covered~scenario+n+method,d,mean); names(s)[4] <- "coverage"
  if (is.null(ylim)) ylim <- c(min(.8,floor(min(s$coverage)*20)/20),1.005)
  order <- c("HAC","Bootstrap","Subsampling","spacer",COSI_METHODS)
  ggplot2::ggplot(s,ggplot2::aes(n,coverage,color=method,group=method))+
    ggplot2::geom_hline(yintercept=1-unique(d$alpha),linetype="dashed",color="#7F8993")+
    ggplot2::geom_line(linewidth=.75)+ggplot2::geom_point(size=1.8)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_color_manual(values=c(DEPENDENT_COLORS,spacer=NA_character_),
      limits=order,breaks=order,labels=c(DEPENDENT_LABELS,spacer="")[order],drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=sort(unique(d$n)))+
    ggplot2::scale_y_continuous(breaks=seq(0,1,.05),labels=function(z) paste0(round(100*z),"%"))+
    ggplot2::coord_cartesian(ylim=ylim)+
    ggplot2::labs(title="Empirical coverage",x="Sample size",y="Empirical coverage",
      color=NULL,subtitle=NULL,caption=NULL)+dependent_plot_theme()+
    ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE,
      override.aes=list(alpha=c(1,1,1,0,1,1,1,1))))
}

# Compute hinges and whiskers before plotting, INCLUDING unbounded intervals.
# Inf values remain in the empirical distribution and in the saved medians.
# This avoids geom_boxplot() silently removing nonfinite observations.
dependent_boxplot_data <- function(d) {
  groups <- split(d,interaction(d$idx,d$n,d$method,drop=TRUE))
  out <- lapply(groups,function(z) {
    q <- grDevices::boxplot.stats(z$width_ratio)$stats
    cbind(z[1,c("idx","scenario","n","method")],
      data.frame(ymin=q[1],lower=q[2],middle=q[3],upper=q[4],ymax=q[5],
        count=NROW(z),infinite=sum(is.infinite(z$width_ratio))))
  })
  ans <- do.call(rbind,out); rownames(ans) <- NULL; ans
}

plot_width_ratios_dependent <- function(csv_file,scenarios=1:4,ylim=c(1/16,4)) {
  d <- read_dependent_results(csv_file,scenarios)
  labels <- DEPENDENT_LABELS[names(DEPENDENT_LABELS)!="unimodal"]
  d <- d[d$method %in% names(labels),]
  d$method <- factor(d$method,levels=names(labels))
  b <- dependent_boxplot_data(d)
  b$n_label <- factor(b$n,levels=sort(unique(d$n)))
  ggplot2::ggplot(b,ggplot2::aes(x=n_label,fill=method,
      ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax))+
    ggplot2::geom_hline(yintercept=1,linetype="longdash",color="#7F8993",linewidth=.5)+
    ggplot2::geom_boxplot(stat="identity",width=.80,
      position=ggplot2::position_dodge2(width=.86,padding=.14),
      outlier.shape=NA,color="#28333D",linewidth=.3,fatten=1.8)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_fill_manual(values=DEPENDENT_COLORS[names(labels)],
      breaks=names(labels),labels=labels,drop=FALSE)+
    ggplot2::scale_y_continuous(breaks=c(1/16,1/8,1/4,1/2,1,2,4),
      labels=c("1/16","1/8","1/4","1/2","1","2","4"),expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_trans(y="log10",ylim=ylim)+
    ggplot2::labs(title="Interval width relative to Wald (HAC)",x="Sample size",
      y="Width / Wald (HAC) width (log scale)",fill=NULL,subtitle=NULL,caption=NULL)+
    dependent_plot_theme()+ggplot2::guides(fill=ggplot2::guide_legend(nrow=2,byrow=TRUE))
}

# Redraw the SAME two 2-by-2 figures from the saved CSV, with no simulation.
# PNG and PDF versions use the same ggplot object.
save_dependent_plots <- function(csv_file,output_dir=dirname(csv_file)) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  plots <- list(empirical_coverage=plot_coverage_n_dependent(csv_file),
    width_ratios=plot_width_ratios_dependent(csv_file))
  for (name in names(plots)) for (ext in c("pdf","png"))
    ggplot2::ggsave(file.path(output_dir,paste0(name,".",ext)),plots[[name]],
      device=if (ext=="pdf") grDevices::cairo_pdf else "png",
      width=13,height=8.5,dpi=200,bg="white")
  invisible(plots)
}

run_dependent_regression_study <- function(core_dir=".",
    output_dir="dependent_regression_results",R=500,n_grid=c(100,200,300,500,750,1000),
    alpha=.05,delta=.30,randomize=TRUE,seed=100,R_boot=500,N_sub=500,
    p=5L,rho=.7,phi_x=.6,phi_e=.6,workers=1L) {
  load_dependent_regression_functions(core_dir)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  sim <- simulate_regression_dependent(R,n_grid,alpha,delta,randomize,seed,
    R_boot,N_sub,p,rho,phi_x,phi_e,workers=workers,
    checkpoint_dir=file.path(output_dir,"checkpoints"))
  csv <- file.path(output_dir,"dependent_regression_intervals.csv")
  saveRDS(sim,file.path(output_dir,"simulation.rds"))
  write.csv(sim$results,csv,row.names=FALSE)
  d <- read_dependent_results(csv)
  coverage <- aggregate(covered~scenario+n+method,d,mean); names(coverage)[4] <- "coverage"
  widths <- aggregate(width_ratio~scenario+n+method,d,median)
  status <- aggregate(replicate~scenario+n+method+status,d,length); names(status)[5] <- "count"
  infinite <- aggregate(width_ratio~scenario+n+method,d,function(z) mean(is.infinite(z)))
  names(infinite)[4] <- "fraction_infinite"
  write.csv(coverage,file.path(output_dir,"empirical_coverage.csv"),row.names=FALSE)
  write.csv(widths,file.path(output_dir,"median_width_ratios.csv"),row.names=FALSE)
  write.csv(widths[widths$method=="unimodal",],
    file.path(output_dir,"unimodal_median_width_ratios.csv"),row.names=FALSE)
  write.csv(status,file.path(output_dir,"fit_status_counts.csv"),row.names=FALSE)
  write.csv(infinite,file.path(output_dir,"infinite_interval_rates.csv"),row.names=FALSE)
  write.csv(sim$design,file.path(output_dir,"block_counts.csv"),row.names=FALSE)
  sub <- d[d$method=="Subsampling",]
  rates <- aggregate(gamma_hat~scenario+n,sub,function(x)
    c(median=median(x),q25=quantile(x,.25),q75=quantile(x,.75)))
  write.csv(rates,file.path(output_dir,"subsampling_rate_summary.csv"),row.names=FALSE)
  plots <- save_dependent_plots(csv,output_dir)
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))
  invisible(list(csv_file=csv,coverage=coverage,median_width_ratios=widths,
    design=sim$design,plots=plots))
}

# ----------------------- EDIT SETTINGS AND RUN -----------------------
# delta=.30 gives9 inference estimates, with a tenth block for the unimodal*
# pilot. Positive gaps mean no permutation. All five covariates and an
# intercept are fitted in every block. At n100, the smallest block has8 rows.
# Width boxplots omit unimodal*, use a log axis limited to[1/16,4], and save
# its median ratios separately. Both legends have exactly two rows; all COSI
# entries appear together on row two. No subtitles, captions, or error bars.
# Every cell has its own seed, giving serial/parallel reproducibility.
RUN_SIMULATION <- TRUE
CORE_DIR <- "."
OUTPUT_DIR <- "dependent_regression_results"
N_GRID <- c(100,200,300,500,750,1000)
REPLICATIONS <- 500
ALPHA <- .05
DELTA <- .30
RANDOMIZE <- TRUE
SEED <- 100
BOOTSTRAP_REPLICATIONS <- 500
SUBSAMPLING_REPLICATIONS <- 500
WORKERS <- 1L  # Increase on Linux/macOS;1 also works on Windows.

if(RUN_SIMULATION) study <- run_dependent_regression_study(core_dir=CORE_DIR,
  output_dir=OUTPUT_DIR,R=REPLICATIONS,n_grid=N_GRID,alpha=ALPHA,delta=DELTA,
  randomize=RANDOMIZE,seed=SEED,R_boot=BOOTSTRAP_REPLICATIONS,
  N_sub=SUBSAMPLING_REPLICATIONS,workers=WORKERS)
