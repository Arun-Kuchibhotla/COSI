# Four regression examples at fixed delta=.30, increasing sample size.
# Editable settings and the call running the entire study are at the END.
# Required packages: weightedBernoulli, boot, ggplot2 (>= 4.0.0).
# Both PDF/PNG figures and all summaries are produced from the saved CSV.

load_regression_functions <- function(core_dir = ".", envir = .GlobalEnv) {
  for (pkg in c("weightedBernoulli", "boot", "ggplot2"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install the package '", pkg, "' before running the study.")
  core <- file.path(core_dir, c("cosi_univariate(1).R", "cosi_univariate.R"))
  core <- core[file.exists(core)]
  if (!length(core)) stop("Cannot find cosi_univariate(1).R or cosi_univariate.R.")
  files <- c(core[1], file.path(core_dir, c("cosi1d.R", "subsampling_CI.R")))
  for (file in files) {
    if (!file.exists(file)) stop("Cannot find ", file)
    for (expr in parse(file)) {
      if (is.call(expr) && is.symbol(expr[[1]]) &&
          as.character(expr[[1]]) %in% c("<-", "=") &&
          is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function")))
        eval(expr, envir = envir)
    }
  }
  invisible(files)
}


COSI_METHODS <- c("median", "sym", "normal", "unimodal")
REGRESSION_SCENARIO_LABELS <- c(
  "1. Correlated Gaussian design", "2. Skewed design and errors",
  "3. Nonlinear mean, Student t(3) errors", "4. Nonlinear mean, skewed Pareto errors")
REGRESSION_METHOD_LABELS <- c(Wald="Wald", Bootstrap="Basic bootstrap",
  Subsampling="Subsampling (estimated rate)", median="COSI median",
  sym="COSI symmetric", normal="COSI normal scores", unimodal="COSI unimodal*")
REGRESSION_COLORS <- c(Wald="#59636E", Bootstrap="#E69F00", Subsampling="#CC79A7",
  median="#0072B2", sym="#009E73", normal="#D55E00", unimodal="#7B2CBF")

# Target: first slope in E[ZZ']^{-1} E[ZY], Z=(1,X'). This agrees with the
# squared-error projection when E[Y^2] is finite. Scenario 4 has infinite error
# variance; its target is defined by the population normal equations.
regression_truth <- function(idx) c(1, 1, 1.75, 1.75)[idx]

# Response FIRST, then p covariates. Rows are independent. Within each row,
# G ~ N(0,Sigma), Sigma[j,k]=rho^|j-k|. X=G, except in scenario 2 where
# X=exp(G/2)-exp(1/8), coordinatewise. beta=(1,-1,.5,.5,-.5).
# Y=mu(X)+sigma(X)*error, with errors independent of G and mean zero.
# 1: mu=1+X'beta, sigma=1, error ~ N(0,1).
# 2: mu=1+X'beta, sigma=.5+.5|X1|, error ~ Exp(1)-1.
# 3: mu=1+X'beta+.5 X1^2+.25 X1^3+.5 X1 X2,
#    sigma=.5+.5|X1|, error ~ t(3)/sqrt(3).
# 4: same nonlinear mean and sigma; error=(1/3)U^(-2/3)-1, U~Unif(0,1).
# Gaussian moments give first population slope 1+3/4=1.75 in scenarios 3-4.
gen_data_regression <- function(n, idx, p=5L, rho=.7) {
  if (!idx %in% 1:4) stop("idx must be in 1:4.")
  if (p < 2L || abs(rho) >= 1) stop("Use p >= 2 and abs(rho) < 1.")
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n*p), n, p) %*% chol(Sigma)
  if (idx == 2L) X <- exp(X/2) - exp(1/8)
  colnames(X) <- paste0("x", seq_len(p))
  beta <- rep(c(1, -1, .5, .5, -.5), length.out=p)
  mu <- 1 + drop(X %*% beta)
  if (idx %in% 3:4)
    mu <- mu + .5*X[,1]^2 + .25*X[,1]^3 + .5*X[,1]*X[,2]
  sigma <- if (idx == 1L) rep(1, n) else .5 + .5*abs(X[,1])
  error <- switch(as.character(idx),
    `1`=rnorm(n), `2`=rexp(n)-1, `3`=rt(n, df=3)/sqrt(3),
    `4`=(1/3)*runif(n)^(-2/3)-1)
  cbind(y=mu + sigma*error, X)
}

# OLS always includes an intercept and ALL supplied covariates.
# A rank-deficient fit is explicitly reported; no covariate is silently dropped.
regression_rank_error <- function() {
  stop(structure(list(message="Rank-deficient regression design.", call=NULL),
    class=c("regression_rank_error", "error", "condition")))
}
ols_coefficient <- function(data, coefficient=1L) {
  Z <- cbind(1, data[, -1, drop=FALSE])
  fit <- .lm.fit(Z, data[,1])
  if (fit$rank < ncol(Z)) regression_rank_error()
  unname(fit$coefficients[coefficient + 1L])
}

# HC0 sandwich covariance, matching vcovHC(..., type="HC") in HulC's notebook.
# Pairs bootstrap, rather than residual bootstrap, accommodates misspecification.
wald_regression <- function(data, alpha=.05, coefficient=1L) {
  Z <- cbind(1, data[, -1, drop=FALSE])
  fit <- .lm.fit(Z, data[,1])
  if (fit$rank < ncol(Z)) regression_rank_error()
  bread <- chol2inv(chol(crossprod(Z)))
  influence <- drop(Z %*% bread[, coefficient + 1L])*fit$residuals
  se <- sqrt(sum(influence^2))
  list(CI=unname(fit$coefficients[coefficient + 1L] +
    c(-1,1)*qnorm(1-alpha/2)*se), standard_error=se)
}

basic_bootstrap_regression <- function(data, alpha=.05, R_boot=500, coefficient=1L) {
  fit <- boot::boot(data, function(d, i)
    ols_coefficient(d[i, , drop=FALSE], coefficient), R=R_boot)
  list(CI=unname(boot::boot.ci(fit, conf=1-alpha, type="basic")$basic[4:5]))
}

# Use your unchanged rate-estimated subsampling_ci(), with regression-feasible
# sizes only. Retain the original 100-point inverse-CDF grid and rate estimate.
subsampling_regression <- function(data, alpha=.05, b=round(NROW(data)^(2/3)),
    N=500, beta=(9:1)/10, min_rows=3L*NCOL(data), coefficient=1L) {
  n <- NROW(data)
  sizes <- round(n^beta)
  keep <- sizes >= min_rows & sizes < n & !duplicated(sizes)
  beta <- beta[keep]
  if (length(beta) < 2L || b < min_rows || b >= n)
    stop("Need at least two feasible rate-estimation sizes and min_rows <= b < n.")
  fit <- subsampling_ci(data, estimate=function(d) ols_coefficient(d, coefficient),
    alpha=alpha, b=b, N=N, beta=beta)
  if (anyNA(fit$CI) || !is.finite(fit$gamma_hat))
    stop("Nonfinite estimated subsampling rate.")
  fit$rate_sizes <- round(n^beta)
  fit
}

# Conservative fallback for the rare singular random design: keep the
# replication, return R, and record its status. Never condition on a successful
# split or redraw an unsuccessful split or subsample.
regression_safe <- function(expr) {
  tryCatch(expr, regression_rank_error=function(e)
    list(CI=c(-Inf, Inf), status="rank_deficient", gamma_hat=NA_real_))
}
regression_record <- function(ans, truth) {
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

regression_calibrations <- function(B_values, alpha) {
  ans <- lapply(B_values, function(B) {
    out <- lapply(c("sym", "normal"), function(method) {
      h <- if (method == "sym") seq_len(B) else
        qnorm(1/2 + seq_len(B)/(2 * (B + 1)))
      z <- weightedBernoulli::qweighted_bernoulli(
        p = alpha, weights = h, prob = .5, lower.tail = FALSE,
        randomized = TRUE, two.sided = TRUE, draw = FALSE)
      z$h <- h
      z
    })
    setNames(out, c("sym", "normal"))
  })
  setNames(ans, as.character(B_values))
}

# Same weighted-method calculation and RNG draws as cosi1D(), using the
# partition estimates returned by its median call and the cached calibration.
regression_weighted_fit <- function(w, n, split_seed, method, cal,
                              alpha, delta, randomize) {
  set.seed(split_seed)
  invisible(sample.int(n))  # Advance past cosi1D()'s identical permutation.
  cutoff <- if (randomize && runif(1) < cal$tau) cal$c_plus else cal$c_minus
  ans <- if (method == "normal") {
    cosi_normal(w, c = cutoff, alpha = alpha)
  } else {
    cosi_sym(w, h = cal$h, c = cutoff, alpha = alpha)
  }
  ans$B <- length(w)
  ans$B_ci <- length(w)
  ans$delta <- delta
  ans$method <- method
  ans$w <- unname(w)
  ans$calibration_exact <- isTRUE(cal$exact)
  ans
}


# Same B and partition for median/sym/normal. Your unchanged unimodal wrapper
# uses B+1 total blocks, reserves the first for an independent pilot, and uses
# B estimates in the interval. It therefore has a different partition.
regression_cosi_fits <- function(data, alpha, delta, randomize, split_seed, cal,
    coefficient=1L) {
  estimate <- function(d) ols_coefficient(d, coefficient)
  med <- regression_safe(cosi1D(data, estimate, alpha=alpha, delta=delta,
    randomize=randomize, seed=split_seed, method="median"))
  B <- smallest_B(alpha, delta)
  if (is.null(med$status)) {
    med$calibration_exact <- TRUE
    sym <- regression_weighted_fit(med$w, NROW(data), split_seed, "sym", cal$sym,
      alpha, delta, randomize)
    normal <- regression_weighted_fit(med$w, NROW(data), split_seed, "normal", cal$normal,
      alpha, delta, randomize)
  } else {
    sym <- normal <- med
    med$B <- sym$B <- normal$B <- B
    med$B_ci <- sym$B_ci <- normal$B_ci <- B
  }
  uni <- regression_safe(cosi1D(data, estimate, alpha=alpha, delta=delta,
    randomize=randomize, seed=split_seed, method="unimodal"))
  uni$B <- B + 1L; uni$B_ci <- B
  uni$calibration_exact <- TRUE
  list(median=med, sym=sym, normal=normal, unimodal=uni)
}


# Fixed cell streams make results independent of parallel scheduling.
regression_cell_seed <- function(seed,idx,n) as.integer(seed + 1000000*idx + n)

regression_design <- function(n_grid=c(100,200,300,500,750,1000), alpha=.05,
    delta=.30, scenarios=1:4, p=5L, min_block=NULL) {
  n_grid <- sort(unique(n_grid))
  if (any(n_grid<1 | n_grid!=floor(n_grid))) stop("n_grid must contain positive integers.")
  max_d <- p+1L
  # Require at least two residual degrees of freedom in the smallest block.
  if (is.null(min_block)) min_block <- max_d+2L
  if (min_block<max_d) stop("min_block must allow all regression coefficients.")
  B <- smallest_B(alpha,delta)
  smallest <- floor(n_grid/(B+1L))
  if (any(smallest<min_block)) stop("Some sample sizes have too few rows per block. ",
    "Increase n or delta; B is not capped and no sample size is silently omitted.")
  data.frame(n=n_grid,delta=delta,B=B,unimodal_B=B+1L,
    median_smallest_block=floor(n_grid/B),unimodal_smallest_block=smallest,
    minimum_required=min_block)
}

illustrate_regression <- function(n=1000,alpha=.05,delta=.30,
    randomize=TRUE,seed=100,scenarios=1:4,p=5L,rho=.7) {
  design <- regression_design(n,alpha,delta,scenarios,p)
  cal <- regression_calibrations(design$B,alpha)[[1L]]
  labels <- REGRESSION_SCENARIO_LABELS
  out <- lapply(scenarios,function(idx) {
    set.seed(regression_cell_seed(seed,idx,n))
    data <- gen_data_regression(n,idx,p,rho)
    split_seed <- sample.int(.Machine$integer.max,1)
    fits <- c(list(Wald=regression_safe(wald_regression(data,alpha))),
      regression_cosi_fits(data,alpha,delta,randomize,split_seed,cal))
    cbind(data.frame(idx=idx,scenario=labels[idx],method=names(fits),n=n,
      delta=delta,truth=regression_truth(idx),
      B=vapply(fits,function(f) if(is.null(f$B)) NA_real_ else f$B,numeric(1))),
      do.call(rbind,lapply(fits,regression_record,truth=regression_truth(idx))))
  })
  ans <- do.call(rbind,out); rownames(ans) <- NULL; ans
}

simulate_coverage_n_regression <- function(R=500,
    n_grid=c(100,200,300,500,750,1000),alpha=.05,delta=.30,
    randomize=TRUE,seed=100,R_boot=500,N_sub=500,subsample_exp=2/3,
    scenarios=1:4,p=5L,rho=.7,min_block=NULL,
    checkpoint_dir=NULL,workers=1L) {
  design <- regression_design(n_grid,alpha,delta,scenarios,p,min_block)
  cal <- regression_calibrations(unique(design$B),alpha)[[1L]]
  labels <- REGRESSION_SCENARIO_LABELS
  jobs <- expand.grid(idx=scenarios,n=design$n)
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  one_cell <- function(j) {
    idx <- jobs$idx[j]; n <- jobs$n[j]
    cell_seed <- regression_cell_seed(seed,idx,n)
    set.seed(cell_seed)
    truth <- regression_truth(idx)
    b <- round(n^subsample_exp)
    out <- vector("list",R)
    for (r in seq_len(R)) {
      data <- gen_data_regression(n,idx,p,rho)
      split_seed <- sample.int(.Machine$integer.max,1)
      sub <- regression_safe(subsampling_regression(data,alpha,b=b,N=N_sub))
      baseline <- list(Wald=regression_safe(wald_regression(data,alpha)),
        Bootstrap=regression_safe(basic_bootstrap_regression(data,alpha,R_boot)),
        Subsampling=sub)
      simulation_rng <- get(".Random.seed",envir=.GlobalEnv)
      fits <- c(baseline,regression_cosi_fits(data,alpha,delta,randomize,split_seed,cal))
      records <- do.call(rbind,lapply(fits,regression_record,truth=truth))
      out[[r]] <- cbind(data.frame(idx=idx,scenario=labels[idx],replicate=r,
        n=n,n_covariates=NCOL(data)-1L,truth=truth,alpha=alpha,delta=delta,
        rho=rho,method=names(fits),split_seed=split_seed,
        B=vapply(fits,function(a) if(is.null(a$B)) NA_real_ else a$B,numeric(1)),
        B_ci=vapply(fits,function(a) if(is.null(a$B_ci)) NA_real_ else a$B_ci,numeric(1)),
        pilot=vapply(fits,function(a) if(is.null(a$pilot)) NA_real_ else a$pilot,numeric(1)),
        gamma_hat=c(NA,NA,sub$gamma_hat,rep(NA_real_,4)),
        b=c(NA,NA,b,rep(NA_real_,4))),records)
      assign(".Random.seed",simulation_rng,envir=.GlobalEnv)
      if (r %% 250L==0L) message("Scenario ",idx,", n=",n,": ",r,"/",R)
    }
    ans <- do.call(rbind,out)
    keys <- c("idx","n","replicate","delta")
    wald <- ans[ans$method=="Wald",c(keys,"width")]
    names(wald)[names(wald)=="width"] <- "wald_width"
    ans <- merge(ans,wald,by=keys,all.x=TRUE,sort=FALSE)
    stopifnot(all(is.finite(ans$wald_width) & ans$wald_width>0))
    ans$width_ratio <- ans$width/ans$wald_width
    if (anyNA(ans$width_ratio)) stop("Undefined paired Wald width ratio.")
    ans$cell_seed <- cell_seed
    if (!is.null(checkpoint_dir)) saveRDS(ans,file.path(checkpoint_dir,
      paste0("scenario_",idx,"_n_",n,".rds")))
    message("Completed scenario ",idx,", n=",n,".")
    ans
  }
  if (workers>1L && .Platform$OS.type!="windows") {
    pieces <- parallel::mclapply(seq_len(nrow(jobs)),one_cell,mc.cores=workers,
      mc.set.seed=FALSE,mc.preschedule=FALSE)
    if (any(vapply(pieces,inherits,logical(1),"try-error"))) stop("A simulation worker failed.")
  } else pieces <- lapply(seq_len(nrow(jobs)),one_cell)
  columns <- names(pieces[[1L]])
  results <- do.call(rbind,lapply(pieces,function(z) z[,columns]))
  results <- results[order(results$idx,results$n,results$replicate,
    match(results$method,names(REGRESSION_METHOD_LABELS))),]
  rownames(results) <- NULL
  list(results=results,design=design,alpha=alpha,
    settings=list(R=R,n_grid=design$n,delta=delta,seed=seed,R_boot=R_boot,N_sub=N_sub,
      p=p,rho=rho,randomize=randomize,scenarios=scenarios,subsample_exp=subsample_exp))
}

read_regression_results <- function(csv_file,scenarios=NULL) {
  d <- read.csv(csv_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!is.null(scenarios)) d <- d[d$idx %in% scenarios,]
  if (!nrow(d)) stop("No rows for the requested scenarios.")
  if (length(unique(d$delta))!=1L || length(unique(d$alpha))!=1L)
    stop("Select one delta and alpha before plotting.")
  stopifnot(!anyNA(d$covered),!anyNA(d$width_ratio))
  d$scenario <- factor(d$scenario,levels=unique(d$scenario[order(d$idx)]))
  d
}

regression_plot_theme <- function() {
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



plot_coverage_n_regression <- function(csv_file,scenarios=1:4,ylim=NULL) {
  d <- read_regression_results(csv_file,scenarios)
  s <- aggregate(covered~scenario+n+method,d,mean); names(s)[4] <- "coverage"
  if (is.null(ylim)) ylim <- c(min(.8,floor(min(s$coverage)*20)/20),1.005)
  order <- c("Wald","Bootstrap","Subsampling","spacer",COSI_METHODS)
  ggplot2::ggplot(s,ggplot2::aes(n,coverage,color=method,group=method))+
    ggplot2::geom_hline(yintercept=1-unique(d$alpha),linetype="dashed",color="#7F8993")+
    ggplot2::geom_line(linewidth=.75)+ggplot2::geom_point(size=1.8)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_color_manual(values=c(REGRESSION_COLORS,spacer=NA_character_),
      limits=order,breaks=order,labels=c(REGRESSION_METHOD_LABELS,spacer="")[order],drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=sort(unique(d$n)))+
    ggplot2::scale_y_continuous(breaks=seq(0,1,.05),labels=function(z) paste0(round(100*z),"%"))+
    ggplot2::coord_cartesian(ylim=ylim)+
    ggplot2::labs(title="Empirical coverage",x="Sample size",y="Empirical coverage",
      color=NULL,subtitle=NULL,caption=NULL)+regression_plot_theme()+
    ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE,
      override.aes=list(alpha=c(1,1,1,0,1,1,1,1))))
}

# Compute hinges and whiskers before plotting, INCLUDING unbounded intervals.
# Inf values remain in the empirical distribution and in the saved medians.
# This avoids geom_boxplot() silently removing nonfinite observations.
regression_boxplot_data <- function(d) {
  groups <- split(d,interaction(d$idx,d$n,d$method,drop=TRUE))
  out <- lapply(groups,function(z) {
    q <- grDevices::boxplot.stats(z$width_ratio)$stats
    cbind(z[1,c("idx","scenario","n","method")],
      data.frame(ymin=q[1],lower=q[2],middle=q[3],upper=q[4],ymax=q[5],
        count=NROW(z),infinite=sum(is.infinite(z$width_ratio))))
  })
  ans <- do.call(rbind,out); rownames(ans) <- NULL; ans
}

plot_width_ratios_regression <- function(csv_file,scenarios=1:4,ylim=c(1/16,4)) {
  d <- read_regression_results(csv_file,scenarios)
  labels <- REGRESSION_METHOD_LABELS[names(REGRESSION_METHOD_LABELS)!="unimodal"]
  d <- d[d$method %in% names(labels),]
  d$method <- factor(d$method,levels=names(labels))
  b <- regression_boxplot_data(d)
  b$n_label <- factor(b$n,levels=sort(unique(d$n)))
  ggplot2::ggplot(b,ggplot2::aes(x=n_label,fill=method,
      ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax))+
    ggplot2::geom_hline(yintercept=1,linetype="longdash",color="#7F8993",linewidth=.5)+
    ggplot2::geom_boxplot(stat="identity",width=.80,
      position=ggplot2::position_dodge2(width=.86,padding=.14),
      outlier.shape=NA,color="#28333D",linewidth=.3,median.linewidth=.55,staplewidth=.35)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_fill_manual(values=REGRESSION_COLORS[names(labels)],
      breaks=names(labels),labels=labels,drop=FALSE)+
    ggplot2::scale_y_continuous(breaks=c(1/16,1/8,1/4,1/2,1,2,4),
      labels=c("1/16","1/8","1/4","1/2","1","2","4"),expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_transform(y="log10",ylim=ylim)+
    ggplot2::labs(title="Interval width relative to Wald",x="Sample size",
      y="Width / Wald width (log scale)",fill=NULL,subtitle=NULL,caption=NULL)+
    regression_plot_theme()+ggplot2::guides(fill=ggplot2::guide_legend(nrow=2,byrow=TRUE))
}

# Redraw the SAME two 2-by-2 figures from the saved CSV, with no simulation.
# PNG and PDF versions use the same ggplot object.
save_regression_plots <- function(csv_file,output_dir=dirname(csv_file)) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  plots <- list(empirical_coverage=plot_coverage_n_regression(csv_file),
    width_ratios=plot_width_ratios_regression(csv_file))
  for (name in names(plots)) for (ext in c("pdf","png"))
    ggplot2::ggsave(file.path(output_dir,paste0(name,".",ext)),plots[[name]],
      device=if (ext=="pdf") grDevices::cairo_pdf else "png",
      width=13,height=8.5,dpi=200,bg="white")
  invisible(plots)
}

run_regression_study <- function(core_dir=".",output_dir="regression_results",
    R=500,n_grid=c(100,200,300,500,750,1000),alpha=.05,delta=.30,
    randomize=TRUE,seed=100,R_boot=500,N_sub=500,subsample_exp=2/3,
    scenarios=1:4,p=5L,rho=.7,min_block=NULL,workers=1L) {
  load_regression_functions(core_dir)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  sim <- simulate_coverage_n_regression(R,n_grid,alpha,delta,randomize,seed,
    R_boot,N_sub,subsample_exp,scenarios,p,rho,min_block,
    checkpoint_dir=file.path(output_dir,"checkpoints"),workers=workers)
  csv <- file.path(output_dir,"regression_sample_size_intervals.csv")
  write.csv(sim$results,csv,row.names=FALSE)
  d <- read_regression_results(csv)
  coverage <- aggregate(covered~scenario+n+method,d,mean); names(coverage)[4] <- "coverage"
  widths <- aggregate(width_ratio~scenario+n+method,d,median)
  status <- aggregate(replicate~scenario+n+method+status,d,length); names(status)[5] <- "count"
  infinite <- aggregate(width_ratio~scenario+n+method,d,function(z) mean(is.infinite(z)))
  names(infinite)[4] <- "fraction_infinite"
  write.csv(coverage,file.path(output_dir,"regression_empirical_coverage.csv"),row.names=FALSE)
  write.csv(widths,file.path(output_dir,"regression_median_width_ratios.csv"),row.names=FALSE)
  unimodal <- reshape(widths[widths$method=="unimodal",c("scenario","n","width_ratio")],
    idvar="scenario",timevar="n",direction="wide")
  names(unimodal) <- sub("width_ratio[.]","n_",names(unimodal))
  write.csv(unimodal,file.path(output_dir,"regression_unimodal_median_ratios.csv"),row.names=FALSE)
  write.csv(status,file.path(output_dir,"regression_fit_status_counts.csv"),row.names=FALSE)
  write.csv(infinite,file.path(output_dir,"regression_infinite_interval_rates.csv"),row.names=FALSE)
  write.csv(sim$design,file.path(output_dir,"regression_block_counts.csv"),row.names=FALSE)
  saveRDS(sim,file.path(output_dir,"simulation.rds"))
  plots <- save_regression_plots(csv,output_dir)
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))
  invisible(list(csv_file=csv,coverage=coverage,unimodal_medians=unimodal,
    design=sim$design,plots=plots))
}

# ----------------------- EDIT SETTINGS AND RUN -----------------------
# delta=.30 had coverage similar to .25 in the mean study, with 9 rather
# than 12 COSI blocks. At n=100 this leaves about 10-11 rows per block for
# six regression coefficients; unimodal* uses 10 blocks including its pilot.
# Width display is 1/16 to 4 (log scale). Unimodal* is omitted from width
# plots, and its median width ratios are saved separately for your caption.
# Both legends have two rows, with all displayed COSI methods on row two.
# No coverage error bars, subtitles, or bottom captions.
# workers=1 is portable to Windows; increase on Linux/macOS if desired.
CORE_DIR <- "."
OUTPUT_DIR <- "regression_results"
N_GRID <- c(100, 200, 300, 500, 750, 1000)
REPLICATIONS <- 500
ALPHA <- .05
DELTA <- .30
RANDOMIZE <- TRUE
SEED <- 100
BOOTSTRAP_REPLICATIONS <- 500
SUBSAMPLING_REPLICATIONS <- 500
WORKERS <- 1L

study <- run_regression_study(core_dir=CORE_DIR,output_dir=OUTPUT_DIR,
  R=REPLICATIONS,n_grid=N_GRID,alpha=ALPHA,delta=DELTA,randomize=RANDOMIZE,
  seed=SEED,R_boot=BOOTSTRAP_REPLICATIONS,N_sub=SUBSAMPLING_REPLICATIONS,
  workers=WORKERS)
