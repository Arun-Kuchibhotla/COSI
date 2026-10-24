# High-dimensional mean simulations and plots.
# Keep this file beside cosi_highdim_mean.R and use that folder as R's
# working directory. Source this file to run the study with settings below.
# Competitive methods are defined first, followed by the data generators,
# simulation loop, CSV summaries, and ggplot figures. Core COSI calculations
# live in cosi_highdim_mean.R, which reuses cosi_univariate.R calibration.
source("cosi_highdim_mean.R")

HD_SCENARIOS <- c("1. Correlated Gaussian", "2. Symmetric heavy tails",
  "3. Rare symmetric shocks", "4. Rare asymmetric shocks")
HD_METHODS <- c("Gaussian","Empirical","Multiplier",
  "median_empirical","mp_empirical","median_mom","mp_mom")
HD_LABELS <- c(Gaussian="Gaussian (known covariance)",Empirical="Empirical bootstrap",
  Multiplier="Gaussian multiplier",median_empirical="cosi-m (empirical)",
  mp_empirical="cosi-m+ (empirical)",median_mom="cosi-m (MOM)",mp_mom="cosi-m+ (MOM)")
HD_COLORS <- c(Gaussian="#59636E",Empirical="#E69F00",Multiplier="#CC79A7",
  median_empirical="#0072B2",mp_empirical="#56B4E9",median_mom="#009E73",mp_mom="#7B2CBF")

# Oracle Gaussian box: a diagnostic for Gaussian approximation failure.
# Covariance is KNOWN for this extra benchmark, not estimated from the data.
# Equicorrelation .5 in scenario1, identity in scenarios2-4.
hd_gaussian_cdf <- function(t,d,rho=0) {
  if (t<=0) return(0)
  if (!is.finite(t)) return(1)
  if (rho==0) return(exp(d*log1p(-2*pnorm(t,lower.tail=FALSE))))
  integrate(function(z) {
    mass <- pnorm((t-sqrt(rho)*z)/sqrt(1-rho))-
      pnorm((-t-sqrt(rho)*z)/sqrt(1-rho))
    exp(d*log(pmax(mass,0))+dnorm(z,log=TRUE))
  },-9,9,rel.tol=1e-8,abs.tol=1e-10,subdivisions=300L)$value
}
hd_gaussian_critical <- function(d,alpha=.05,rho=0) {
  independent <- qnorm(-expm1(log1p(-alpha)/d)/2,lower.tail=FALSE)
  if (rho==0) return(independent)
  uniroot(function(t) hd_gaussian_cdf(t,d,rho)-(1-alpha),
    c(0,independent+1),tol=1e-9)$root
}

# Bootstrap helpers for simultaneous mean confidence boxes.
# These helpers use supplied, fixed coordinate scales; they do not recompute
# standard deviations inside individual resamples.
hd_rowmax_abs <- function(M) {
  A <- abs(M)
  A[cbind(seq_len(nrow(A)), max.col(A, ties.method="first"))]
}

# Columns 1:R_boot are empirical-bootstrap weights; remaining columns are
# centered Gaussian multiplier weights. Every column sums to zero up to
# floating-point rounding, so the weights can multiply the original data.
hd_bootstrap_weights <- function(n, R_boot=500L) {
  if (length(n)!=1L || !is.finite(n) || n<2 || n!=floor(n))
    stop("n must be an integer at least two.")
  if (length(R_boot)!=1L || !is.finite(R_boot) || R_boot<1 || R_boot!=floor(R_boot))
    stop("R_boot must be a positive integer.")
  empirical <- (stats::rmultinom(R_boot, size=n, prob=rep(1/n,n))-1)/sqrt(n)
  multiplier <- matrix(stats::rnorm(n*R_boot), nrow=n, ncol=R_boot)
  multiplier <- sweep(multiplier, 2L, colMeans(multiplier), "-")/sqrt(n)
  weights <- cbind(empirical, multiplier)
  attr(weights,"R_boot") <- as.integer(R_boot)
  weights
}

# Maxima are accumulated over column chunks. One combined cross-product per
# chunk evaluates both bootstrap methods without constructing a d x d matrix.
hd_bootstrap_critical <- function(X, weights, scale, alpha=.05, column_chunk=512L) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X)<2L || ncol(X)<1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix with at least two observations.")
  if (!is.matrix(weights) || !is.numeric(weights) || nrow(weights)!=nrow(X) ||
      ncol(weights)<2L || ncol(weights)%%2L!=0L || any(!is.finite(weights)))
    stop("weights must have nrow(X) rows and 2*R_boot finite numeric columns.")
  if (!is.numeric(scale) || length(scale)!=ncol(X) || any(!is.finite(scale) | scale<=0))
    stop("scale must contain one positive finite standard deviation per coordinate.")
  if (length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("alpha must lie in (0,1).")
  if (length(column_chunk)!=1L || !is.finite(column_chunk) ||
      column_chunk<1L || column_chunk!=floor(column_chunk))
    stop("column_chunk must be a positive integer.")
  R_boot <- ncol(weights)%/%2L
  maxima <- numeric(2L*R_boot)
  for (first in seq.int(1L,ncol(X),by=column_chunk)) {
    columns <- seq.int(first,min(ncol(X),first+column_chunk-1L))
    values <- crossprod(weights,X[,columns,drop=FALSE])
    values <- sweep(values,2L,scale[columns],"/")
    maxima <- pmax(maxima,hd_rowmax_abs(values))
  }
  empirical_draws <- maxima[seq_len(R_boot)]
  multiplier_draws <- maxima[R_boot+seq_len(R_boot)]
  quantile_index <- min(R_boot,ceiling((R_boot+1)*(1-alpha)))
  list(empirical=unname(sort(empirical_draws,partial=quantile_index)[quantile_index]),
    multiplier=unname(sort(multiplier_draws,partial=quantile_index)[quantile_index]),
    empirical_draws=empirical_draws,multiplier_draws=multiplier_draws,
    quantile_index=quantile_index,R_boot=R_boot)
}


hd_rectangle <- function(center,halfwidth,critical=NA_real_) {
  lower <- center-halfwidth; upper <- center+halfwidth
  log_volume <- sum(log(2*halfwidth))
  list(lower=lower,upper=upper,contains=.hd_membership_function(lower,upper),
    critical=critical,log_volume=log_volume,volume_root=exp(log_volume/length(center)))
}

# All methods use the SAME data. Both COSI pilot choices also use the SAME
# sample split and cutoff draws. The full-data bootstrap benchmarks are
# unchanged, and their original coordinate SDs are fixed across resamples.
hd_compare <- function(X,alpha=.05,delta=.30,n_pilot=floor(nrow(X)/2),
    randomize=TRUE,R_boot=500,seed=NULL,gaussian_rho=0,column_chunk=512L,
    gaussian_critical=NULL) {
  if (!is.null(seed)) set.seed(seed)
  cs <- cosi_hd(X,alpha=alpha,delta=delta,pilot="both",n_pilot=n_pilot,
    randomize=randomize,seed=NULL)
  n <- nrow(X); d <- ncol(X); center <- colMeans(X)
  residuals <- sweep(X,2L,center,"-")
  sample_sd <- sqrt(colSums(residuals^2)/(n-1))
  if (any(sample_sd<=0)) stop("A bootstrap coordinate has zero sample SD.")
  boot <- hd_bootstrap_critical(residuals,hd_bootstrap_weights(n,R_boot),sample_sd,
    alpha,column_chunk)
  q <- if(is.null(gaussian_critical)) hd_gaussian_critical(d,alpha,gaussian_rho)
    else gaussian_critical
  fits <- list(Gaussian=hd_rectangle(center,rep(q/sqrt(n),d),q),
    Empirical=hd_rectangle(center,boot$empirical*sample_sd/sqrt(n),boot$empirical),
    Multiplier=hd_rectangle(center,boot$multiplier*sample_sd/sqrt(n),boot$multiplier),
    median_empirical=cs$empirical$median,mp_empirical=cs$empirical$mp,
    median_mom=cs$mom$median,mp_mom=cs$mom$mp)
  list(fits=fits,B=cs$mom$B,n_pilot=n_pilot,mom_groups=cs$mom$pilot_groups,
    pilot_methods=c("none","none","none","empirical","empirical","mom","mom"),
    pilot_groups=c(rep(NA_integer_,5),rep(cs$mom$pilot_groups,2)),
    scale_floored=c(rep(0L,3),rep(cs$empirical$scale_regularized,2),rep(cs$mom$scale_regularized,2)),
    max_mean=sqrt(n)*max(abs(center)),gaussian_critical=q)
}

hd_record <- function(fit,truth) data.frame(covered=fit$contains(truth),
  log_volume=fit$log_volume,volume_root=fit$volume_root,
  critical=if(is.null(fit$critical)) NA_real_ else fit$critical,
  cutoff=if(is.null(fit$cutoff)) NA_real_ else fit$cutoff)

# All four settings have mean zero and UNIT marginal variances.
# 1: X_ik=(Z_ik+F_i)/sqrt(2): equicorrelation .5.
# 2: independent t(2.5)/sqrt(5): finite variance, infinite 3rd/4th moments.
# 3: X_ik=(Z_ik+sqrt(n*d)*R_ik*A_ik)/sqrt(2),
#    A_ik~Bernoulli(1/(n*d)), R_ik~Rademacher. Covariance=I_d.
#    This is explicitly a dimension-dependent triangular array, not prefixes
#    of one fixed law. Exactly one shock has probability about exp(-1).
# 4: X_ik=(Z_ik+30*(A_ik-.01))/sqrt(1+900*.01*.99), A~Bernoulli(.01).
#    This also stress-tests the marginal-median assumption of COSI.
gen_data_hd_mean <- function(n,d,idx) {
  shocks <- 0L
  if (idx==1L) {
    X <- (matrix(rnorm(n*d),n,d)+rnorm(n))/sqrt(2)
  } else if (idx==2L) {
    X <- matrix(rt(n*d,df=2.5),n,d)/sqrt(5)
  } else if (idx==3L) {
    X <- matrix(rnorm(n*d),n,d)
    hit <- which(runif(n*d)<1/(n*d)); shocks <- length(hit)
    if (shocks) X[hit] <- X[hit]+sqrt(n*d)*sample(c(-1,1),shocks,replace=TRUE)
    X <- X/sqrt(2)
  } else if (idx==4L) {
    X <- matrix(rnorm(n*d),n,d)-.30
    hit <- which(runif(n*d)<.01); shocks <- length(hit)
    if (shocks) X[hit] <- X[hit]+30
    X <- X/sqrt(1+900*.01*.99)
  } else stop("idx must be in 1:4.")
  attr(X,"shock_count") <- shocks
  X
}

simulate_hd_mean <- function(R=500,n=100,d_grid=c(10,30,100,300,1000,3000,10000),
    alpha=.05,delta=.30,n_pilot=floor(n/2),randomize=TRUE,seed=100,
    R_boot=500,column_chunk=512L,workers=1L,checkpoint_dir=NULL) {
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  d_grid <- sort(unique(as.integer(d_grid)))
  # MOM tuning depends on alpha only; use the same K at every dimension.
  .hd_mom_groups(alpha,max_groups=floor(n_pilot/2))
  if (n-n_pilot<smallest_B(alpha,delta)) stop("Too few inference observations.")
  jobs <- expand.grid(idx=1:4,d=d_grid)
  q_ind <- vapply(d_grid,hd_gaussian_critical,numeric(1),alpha=alpha,rho=0)
  q_cor <- vapply(d_grid,hd_gaussian_critical,numeric(1),alpha=alpha,rho=.5)
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  one_cell <- function(j) {
    idx <- jobs$idx[j]; d <- jobs$d[j]
    cell_seed <- as.integer(seed+100000*idx+d); set.seed(cell_seed)
    q <- if(idx==1L) q_cor[match(d,d_grid)] else q_ind[match(d,d_grid)]
    out <- vector("list",R); example <- NULL
    for (r in seq_len(R)) {
      X <- gen_data_hd_mean(n,d,idx)
      fit <- hd_compare(X,alpha,delta,n_pilot,randomize,R_boot,
        gaussian_rho=if(idx==1L) .5 else 0,column_chunk=column_chunk,gaussian_critical=q)
      if (r==1L) example <- fit$fits
      records <- do.call(rbind,lapply(fit$fits,hd_record,truth=numeric(d)))
      out[[r]] <- cbind(data.frame(idx=idx,scenario=HD_SCENARIOS[idx],d=d,n=n,
        replicate=r,method=names(fit$fits),alpha=alpha,delta=delta,
        B=fit$B,n_pilot=n_pilot,pilot_method=fit$pilot_methods,
        pilot_groups=fit$pilot_groups,mom_groups=fit$mom_groups,
        scale_floored=fit$scale_floored,shock_count=attr(X,"shock_count"),
        oracle_max=fit$max_mean,oracle_critical=q,cell_seed=cell_seed),records)
      if (r%%100L==0L) message("Scenario ",idx,", d=",d,": ",r,"/",R)
    }
    ans <- do.call(rbind,out); rownames(ans) <- NULL
    ans$width_ratio <- ans$volume_root/(2*q/sqrt(n))
    if (!is.null(checkpoint_dir)) saveRDS(list(results=ans,example=example),
      file.path(checkpoint_dir,paste0("scenario_",idx,"_d_",d,".rds")))
    message("Completed scenario ",idx,", d=",d,".")
    list(results=ans,example=example)
  }
  if (workers>1L && .Platform$OS.type!="windows") {
    pieces <- parallel::mclapply(seq_len(nrow(jobs)),one_cell,mc.cores=workers,
      mc.set.seed=FALSE,mc.preschedule=FALSE)
    if (any(vapply(pieces,inherits,logical(1),"try-error"))) stop("A worker failed.")
  } else pieces <- lapply(seq_len(nrow(jobs)),one_cell)
  results <- do.call(rbind,lapply(pieces,`[[`,"results"))
  results <- results[order(results$idx,results$d,results$replicate,
    match(results$method,HD_METHODS)),]; rownames(results) <- NULL
  examples <- lapply(pieces,`[[`,"example")
  names(examples) <- paste0("scenario_",jobs$idx,"_d_",jobs$d)
  list(results=results,examples=examples,settings=list(R=R,n=n,d_grid=d_grid,
    alpha=alpha,delta=delta,n_pilot=n_pilot,randomize=randomize,seed=seed,R_boot=R_boot,
    pilot_methods=c("empirical","mom"),shared_split=TRUE,shared_cutoffs=TRUE,
    mom_groups=.hd_mom_groups(alpha),mom_tuning="scalar_binomial_alpha_only",
    rng_kind=RNGkind()))
}


hd_plot_theme <- function() {
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





read_hd_results <- function(csv_file) {
  d <- read.csv(csv_file,stringsAsFactors=FALSE)
  if (length(unique(d$alpha))!=1L || length(unique(d$delta))!=1L)
    stop("Select one alpha and delta before plotting.")
  d$scenario <- factor(d$scenario,levels=HD_SCENARIOS)
  d$method <- factor(d$method,levels=HD_METHODS)
  d
}

plot_hd_coverage <- function(csv_file) {
  d <- read_hd_results(csv_file)
  s <- aggregate(covered~scenario+d+method,d,mean); names(s)[4] <- "coverage"
  ggplot2::ggplot(s,ggplot2::aes(d,coverage,color=method,group=method))+
    ggplot2::geom_hline(yintercept=1-unique(d$alpha),linetype="dashed",color="#7F8993")+
    ggplot2::geom_line(linewidth=.75)+ggplot2::geom_point(size=1.8)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_color_manual(values=c(HD_COLORS,spacer=NA_character_),
      limits=c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7]),
      breaks=c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7]),
      labels=c(HD_LABELS,spacer="")[c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7])],drop=FALSE)+
    ggplot2::scale_x_log10(breaks=sort(unique(d$d)),labels=sort(unique(d$d)))+
    ggplot2::scale_y_continuous(breaks=seq(0,1,.1),labels=function(z) paste0(round(100*z),"%"))+
    ggplot2::coord_cartesian(ylim=c(0,1.005))+
    ggplot2::labs(title="Empirical simultaneous coverage",x="Dimension",y="Empirical coverage",
      color=NULL,subtitle=NULL,caption=NULL)+hd_plot_theme()+
    ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE,
      override.aes=list(alpha=c(1,1,1,0,1,1,1,1))))
}

plot_hd_width_ratios <- function(csv_file,ylim=NULL) {
  d <- read_hd_results(csv_file)
  if (is.null(ylim)) {
    q <- aggregate(width_ratio~scenario+d+method,d,function(x) quantile(x,.9))
    upper <- max(4,2^ceiling(log2(max(q$width_ratio))))
    ylim <- c(1/16,upper)
  }
  d$dimension <- factor(d$d,levels=sort(unique(d$d)))
  groups <- split(d,interaction(d$idx,d$d,d$method,drop=TRUE))
  b <- do.call(rbind,lapply(groups,function(z) {
    q <- grDevices::boxplot.stats(z$width_ratio)$stats
    cbind(z[1,c("scenario","dimension","method")],
      data.frame(ymin=q[1],lower=q[2],middle=q[3],upper=q[4],ymax=q[5]))
  }))
  breaks <- 2^seq(floor(log2(ylim[1])),ceiling(log2(ylim[2])))
  labels <- ifelse(breaks<1,paste0("1/",round(1/breaks)),as.character(breaks))
  ggplot2::ggplot(b,ggplot2::aes(dimension,fill=method,ymin=ymin,lower=lower,
      middle=middle,upper=upper,ymax=ymax))+
    ggplot2::geom_hline(yintercept=1,linetype="longdash",color="#7F8993",linewidth=.5)+
    ggplot2::geom_boxplot(stat="identity",width=.82,
      position=ggplot2::position_dodge2(width=.86,padding=.14),
      outlier.shape=NA,color="#28333D",linewidth=.3,median.linewidth=.55,staplewidth=.35)+
    ggplot2::facet_wrap(~scenario,ncol=2)+
    ggplot2::scale_fill_manual(values=c(HD_COLORS,spacer=NA_character_),
      limits=c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7]),
      breaks=c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7]),
      labels=c(HD_LABELS,spacer="")[c(HD_METHODS[1:3],"spacer",HD_METHODS[4:7])],drop=FALSE)+
    ggplot2::scale_y_continuous(breaks=breaks,labels=labels,expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_transform(y="log10",ylim=ylim)+
    ggplot2::labs(title="Set size relative to the Gaussian reference",x="Dimension",
      y="Volume-equivalent width ratio (log scale)",fill=NULL,subtitle=NULL,caption=NULL)+
    hd_plot_theme()+ggplot2::guides(fill=ggplot2::guide_legend(nrow=2,byrow=TRUE,
      override.aes=list(alpha=c(1,1,1,0,1,1,1,1))))
}

save_hd_plots <- function(csv_file,output_dir=dirname(csv_file),width_ylim=NULL) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  plots <- list(empirical_coverage=plot_hd_coverage(csv_file),
    width_ratios=plot_hd_width_ratios(csv_file,width_ylim))
  for (name in names(plots)) for (ext in c("pdf","png"))
    ggplot2::ggsave(file.path(output_dir,paste0(name,".",ext)),plots[[name]],
      device=if(ext=="pdf") grDevices::cairo_pdf else "png",
      width=13,height=8.5,dpi=200,bg="white")
  invisible(plots)
}

# Summaries and figures are recomputed from the SAVED replication-level CSV.
run_hd_study <- function(output_dir="highdim_mean_results",R=500,n=100,
    d_grid=c(10,30,100,300,1000,3000,10000),alpha=.05,delta=.30,
    n_pilot=floor(n/2),randomize=TRUE,seed=100,R_boot=500,
    column_chunk=512L,workers=1L,width_ylim=NULL) {
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 first.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  sim <- simulate_hd_mean(R,n,d_grid,alpha,delta,n_pilot,randomize,seed,
    R_boot,column_chunk,workers,checkpoint_dir=file.path(output_dir,"checkpoints"))
  csv <- file.path(output_dir,"highdim_mean_intervals.csv")
  saveRDS(sim$results,file.path(output_dir,"simulation_records.rds"))
  write.csv(sim$results,csv,row.names=FALSE)
  d <- read_hd_results(csv)
  coverage <- aggregate(covered~scenario+d+method,d,mean); names(coverage)[4] <- "coverage"
  widths <- aggregate(width_ratio~scenario+d+method,d,median)
  write.csv(coverage,file.path(output_dir,"empirical_coverage.csv"),row.names=FALSE)
  write.csv(widths,file.path(output_dir,"median_width_ratios.csv"),row.names=FALSE)
  diag <- d[d$method=="Gaussian",]
  cells <- split(diag,interaction(diag$idx,diag$d,drop=TRUE))
  diagnostic <- do.call(rbind,lapply(cells,function(z) {
    t <- sort(z$oracle_max); R <- length(t)
    cdf <- vapply(t,hd_gaussian_cdf,numeric(1),d=z$d[1],rho=if(z$idx[1]==1L) .5 else 0)
    data.frame(idx=z$idx[1],scenario=z$scenario[1],d=z$d[1],
      gaussian_coverage=mean(z$covered),
      gaussian_cdf_KS=max(seq_len(R)/R-cdf,cdf-(seq_len(R)-1)/R),
      probability_any_shock=mean(z$shock_count>0),
      mom_groups=z$mom_groups[1],scale_floored=max(z$scale_floored))
  }))
  diagnostic <- diagnostic[order(diagnostic$idx,diagnostic$d),]
  write.csv(diagnostic,file.path(output_dir,"gaussian_approximation_diagnostics.csv"),row.names=FALSE)
  saveRDS(sim$examples,file.path(output_dir,"example_sets.rds"))
  saveRDS(sim$settings,file.path(output_dir,"settings.rds"))
  plots <- save_hd_plots(csv,output_dir,width_ylim)
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))
  invisible(list(csv_file=csv,coverage=coverage,median_widths=widths,
    diagnostics=diagnostic,examples=sim$examples,plots=plots))
}

# ----------------------- EDIT SETTINGS AND RUN -----------------------
# The TOTAL n includes the pilot half. At n=100, the other50 rows supply
# B=9 means from blocks of sizes5/6. Both pilots share all inference blocks.
# The MOM mean and covariance entries use the same alpha-tuned group count,
# independent of d: K=3 at alpha=.05 for every dimension.
# Exact set volume, membership, and covariance estimation live in the CORE file.
# Set RUN_SIMULATION=FALSE to load definitions and redraw existing CSV results.
RUN_SIMULATION <- TRUE
OUTPUT_DIR <- "highdim_mean_results"
REPLICATIONS <- 500
N <- 100
D_GRID <- c(10,30,100,300,1000,3000,10000)
ALPHA <- .05
DELTA <- .30
N_PILOT <- 50
RANDOMIZE <- TRUE
SEED <- 100
BOOTSTRAP_REPLICATIONS <- 500
WORKERS <- 1L  # Portable to Windows; increase on Linux/macOS if desired.
WIDTH_YLIM <- NULL

if (RUN_SIMULATION) study <- run_hd_study(output_dir=OUTPUT_DIR,
  R=REPLICATIONS,n=N,d_grid=D_GRID,alpha=ALPHA,delta=DELTA,n_pilot=N_PILOT,
  randomize=RANDOMIZE,seed=SEED,R_boot=BOOTSTRAP_REPLICATIONS,
  workers=WORKERS,width_ylim=WIDTH_YLIM)
