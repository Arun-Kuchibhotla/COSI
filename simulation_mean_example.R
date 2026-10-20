# Mean confidence interval simulations and plots from the saved CSV.
# Keep this file beside your three supplied core R files.
# Run once:
#   source("simulation_mean_example.R")
#   study <- run_mean_study(core_dir = ".", output_dir = "mean_simulation_results")
# To redraw figures later without simulations or loading core functions:
#   source("simulation_mean_example.R")
#   plot_width_ratios("mean_simulation_results/mean_simulation_intervals.csv")
#   plot_coverage_delta("mean_simulation_results/mean_simulation_intervals.csv")

# Load function definitions unchanged. Skip the core files' top-level package
# installation commands and their hard-coded source("cosi_univariate.R").
load_mean_functions <- function(core_dir = ".", envir = .GlobalEnv) {
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
COSI_SCENARIO_LABELS <- c("1. Gaussian", "2. Student t (df = 1.5)",
                          "3. Four dominant shocks", "4. Heteroscedastic Gaussian")
COSI_METHOD_LABELS <- c(Wald="Wald", median="COSI median", sym="COSI symmetric",
                        normal="COSI normal scores", unimodal="COSI unimodal*")
COSI_COLORS <- c(Wald="#424B54", median="#0072B2", sym="#009E73",
                 normal="#D55E00", unimodal="#9B59B6")
# (1) Four independent common-mean-zero designs.
# Scenario 3 is a triangular array. Scenario 4 has a bounded variance profile:
# approximately 10% of observations have SD 8, all remaining observations SD 1.
gen_data_mean <- function(n,idx) {
  switch(idx,
         rnorm(n),
         rt(n,df=1.5),
         {x <- rnorm(n); k <- seq_len(min(4L,n));
         x[k] <- x[k]+n*sample(c(-1,1),length(k),replace=TRUE); x},
         {sigma <- rep(1,n); sigma[seq.int(1L,n,by=10L)] <- 8; rnorm(n,sd=sigma)})
}

# (2) One data realization per scenario; all four COSI methods use the same
# partition. Seeded data generation and calls leave the caller's RNG intact.
.cosi_record <- function(ans,truth=0) {
  ci <- ans$CI
  closed <- if (is.null(ans$closed)) c(TRUE,TRUE) else ans$closed
  empty <- isTRUE(ans$empty) || ci[1]>ci[2] || (ci[1]==ci[2] && !all(closed))
  covered <- !empty && (truth>ci[1] || (truth==ci[1] && closed[1])) &&
    (truth<ci[2] || (truth==ci[2] && closed[2]))
  data.frame(lower=ci[1],upper=ci[2],width=if(empty) 0 else ci[2]-ci[1],
             covered=covered,closed_lower=closed[1],closed_upper=closed[2],empty=empty,
             cutoff=if(is.null(ans[["c"]])) NA_real_ else ans[["c"]],
             calibration_exact=if(is.null(ans$calibration_exact)) NA else ans$calibration_exact)
}

illustrate_mean <- function(n=100,alpha=.05,delta=1/3,randomize=TRUE,seed=100) {
  ans <- NULL
  set.seed(seed)
  for (idx in 1:4) {
    x <- gen_data_mean(n, idx)
    split_seed <- sample.int(.Machine$integer.max, 1)
    
    for (m in c("Wald", COSI_METHODS)) {
      fit <- if (m == "Wald") {
        list(CI=mean(x)+c(-1,1)*qnorm(1-alpha/2)*sd(x)/sqrt(length(x)),
             closed=c(TRUE,TRUE),empty=FALSE)
      } else {
        cosi1D(x, mean, alpha, delta, randomize,
             seed = split_seed, method = m)
      }
      
      row <- data.frame(
        idx = idx,
        scenario = COSI_SCENARIO_LABELS[idx],
        method = m,
        B = if (m == "Wald") NA else fit$B,
        delta = delta
      )
      
      ans <- rbind(ans, cbind(row, .cosi_record(fit)))
    }
  }
  
  rownames(ans) <- NULL
  ans
}

# ggplot functions return ggplot objects: append + theme(...), + labs(...), etc.
.cosi_plot_data <- function(d) {
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 first.")
  d$scenario <- factor(d$scenario,levels=COSI_SCENARIO_LABELS)
  d$method <- factor(d$method,levels=names(COSI_METHOD_LABELS))
  d
}
.cosi_theme <- function() ggplot2::theme_minimal(base_size=12)+
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                 strip.text=ggplot2::element_text(face="bold"),legend.position="bottom",
                 legend.title=ggplot2::element_blank(),plot.title.position="plot",
                 plot.title=ggplot2::element_text(face="bold"),plot.caption=ggplot2::element_text(hjust=0))
.cosi_colors <- function(fill=FALSE) {
  if(fill) ggplot2::scale_fill_manual(values=COSI_COLORS,labels=COSI_METHOD_LABELS,drop=FALSE) else
    ggplot2::scale_color_manual(values=COSI_COLORS,labels=COSI_METHOD_LABELS,drop=FALSE)
}

plot_illustration_mean <- function(illustration) {
  d <- .cosi_plot_data(illustration)
  good <- d[!d$empty & is.finite(d$lower) & is.finite(d$upper),]
  spans <- lapply(split(good,good$scenario,drop=TRUE),function(z) {
    ran <- range(c(z$lower,z$upper,0));pad <- max(diff(ran)*.1,.01)
    data.frame(scenario=z$scenario[1],left=ran[1]-pad,right=ran[2]+pad)
  })
  bounds <- do.call(rbind,spans)
  whole <- merge(d[!d$empty & is.infinite(d$lower) & is.infinite(d$upper),],bounds,by="scenario")
  ggplot2::ggplot(good,ggplot2::aes(y=method,color=method))+
    ggplot2::geom_vline(xintercept=0,linetype=2,color="grey55")+
    ggplot2::geom_segment(ggplot2::aes(x=lower,xend=upper,yend=method),linewidth=.8)+
    ggplot2::geom_point(ggplot2::aes(x=lower),size=1.8)+
    ggplot2::geom_point(ggplot2::aes(x=upper),size=1.8)+
    ggplot2::geom_segment(data=whole,ggplot2::aes(x=left,xend=right,yend=method),
                          arrow=grid::arrow(ends="both",length=grid::unit(.07,"inches")),linewidth=.8)+
    ggplot2::facet_wrap(~scenario,ncol=2,scales="free_x")+.cosi_colors()+
    ggplot2::scale_y_discrete(limits=rev(names(COSI_METHOD_LABELS)),
                              labels=rev(c("Wald","Median","Sym","Normal","Unimodal*")),drop=FALSE)+
    ggplot2::labs(title="One realization of each confidence interval",x="Parameter value",y=NULL)+
    .cosi_theme()+ggplot2::theme(legend.position="none")
}

# The following code plots one replication of confidence intervals.
# tmp <- illustrate_mean(delta = 0.3, seed = 1200)
# plot_illustration_mean(tmp)

# Compute each data-independent weighted calibration once. The supplied
# cosi1D() and low-level interval functions are not modified.
mean_calibrations <- function(B_values, alpha) {
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
mean_weighted_fit <- function(w, n, split_seed, method, cal,
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

add_width_ratios <- function(results) {
  wald <- results[results$method == "Wald",
    c("scenario", "replicate", "delta", "width")]
  names(wald)[4] <- "wald_width"
  stopifnot(!anyDuplicated(wald[c("scenario", "replicate", "delta")]))
  out <- merge(results, wald,
    by = c("scenario", "replicate", "delta"), all.x = TRUE, sort = FALSE)
  stopifnot(nrow(out) == nrow(results),
    all(is.finite(out$wald_width) & out$wald_width > 0))
  out$width_ratio <- out$width / out$wald_width
  out <- out[order(out$idx, out$replicate, out$delta,
    match(out$method, c("Wald", "Bootstrap", "Subsampling", COSI_METHODS))), ]
  rownames(out) <- NULL
  out
}

simulate_coverage_delta <- function(
    R = 500, n = 100, alpha = .05,
    deltas = c(.10, .15, .20, .25, .30, 1/3),
    randomize = TRUE, seed = 100,
    R_boot = 500, N_sub = 500, subsample_exp = 2/3,
    checkpoint_dir = NULL) {
  if (!is.null(seed)) set.seed(seed)
  B_grid <- vapply(deltas, function(d) smallest_B(alpha, delta = d), numeric(1))
  # Your unimodal wrapper uses one extra block for its independent pilot.
  keep <- B_grid + 1L <= n
  if (any(!keep)) warning("Dropping delta values requiring too many blocks, including the pilot.")
  deltas <- deltas[keep]
  B_grid <- B_grid[keep]
  stopifnot(length(deltas) > 0)
  calibrations <- mean_calibrations(unique(B_grid), alpha)
  b <- round(n^subsample_exp)
  out <- vector("list", 4L * R * length(deltas))
  at <- 0L
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  for (idx in 1:4) {
    for (r in seq_len(R)) {
      x <- gen_data_mean(n, idx)
      split_seed <- sample.int(.Machine$integer.max, 1)
      sub <- subsampling_ci(x, estimate = mean, b = b, N = N_sub, alpha = alpha)
      if (any(!is.finite(c(sub$CI, sub$gamma_hat))))
        stop("Nonfinite subsampling result in scenario ", idx, ", replication ", r)
      baseline <- list(
        Wald = list(CI = mean(x) + c(-1, 1) * qnorm(1 - alpha/2) * sd(x)/sqrt(n)),
        Bootstrap = basic_bootstrap_mean(x, alpha, R_boot),
        Subsampling = sub)

      # cosi1D() calls set.seed(). Restore the simulation RNG after its calls
      # so that future samples do not depend on the last method/delta called.
      simulation_rng <- get(".Random.seed", envir = .GlobalEnv)
      for (g in seq_along(deltas)) {
        d <- deltas[g]
        med <- cosi1D(x, mean, alpha = alpha, delta = d,
          randomize = randomize, seed = split_seed, method = "median")
        med$calibration_exact <- TRUE
        cal <- calibrations[[as.character(B_grid[g])]]
        fits <- list(
          median = med,
          sym = mean_weighted_fit(med$w, n, split_seed, "sym", cal$sym,
            alpha, d, randomize),
          normal = mean_weighted_fit(med$w, n, split_seed, "normal", cal$normal,
            alpha, d, randomize),
          unimodal = cosi1D(x, mean, alpha = alpha, delta = d,
            randomize = randomize, seed = split_seed, method = "unimodal"))
        fits$unimodal$calibration_exact <- TRUE
        fits <- c(baseline, fits)
        records <- do.call(rbind, lapply(fits, .cosi_record))
        at <- at + 1L
        out[[at]] <- cbind(data.frame(
          idx = idx, scenario = COSI_SCENARIO_LABELS[idx], replicate = r,
          n = n, alpha = alpha, delta = d, method = names(fits),
          B = vapply(fits, function(a) if (is.null(a$B)) NA_real_ else a$B, numeric(1)),
          B_ci = vapply(fits, function(a) if (is.null(a$B_ci)) NA_real_ else a$B_ci, numeric(1)),
          pilot = vapply(fits, function(a) if (is.null(a$pilot)) NA_real_ else a$pilot, numeric(1)),
          gamma_hat = c(NA, NA, sub$gamma_hat, rep(NA_real_, 4))), records)
      }
      assign(".Random.seed", simulation_rng, envir = .GlobalEnv)
      if (r %% 100L == 0L) message("Scenario ", idx, ": ", r, "/", R, " replications.")
    }
    if (!is.null(checkpoint_dir)) saveRDS(do.call(rbind, out[seq_len(at)]),
      file.path(checkpoint_dir, paste0("scenario_", idx, "_checkpoint.rds")))
    message("Completed scenario ", idx, ": ", R, " replications.")
  }
  results <- add_width_ratios(do.call(rbind, out))
  coverage <- aggregate(covered ~ scenario + delta + method, results, mean)
  names(coverage)[4] <- "coverage"
  list(results = results, coverage = coverage,
    design = data.frame(delta = deltas, B = B_grid, unimodal_B = B_grid + 1L),
    alpha = alpha, subsample_size = b)
}

MEAN_METHOD_LABELS <- c(Wald = "Wald", Bootstrap = "Basic bootstrap",
  Subsampling = "Subsampling (estimated rate)", median = "COSI median",
  sym = "COSI symmetric", normal = "COSI normal scores", unimodal = "COSI unimodal*")
MEAN_COLORS <- c(Wald = "#59636E", Bootstrap = "#E69F00", Subsampling = "#CC79A7",
  median = "#0072B2", sym = "#009E73", normal = "#D55E00", unimodal = "#7B2CBF")

read_mean_results <- function(csv_file) {
  d <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("scenario", "replicate", "delta", "method", "width_ratio", "covered") %in% names(d)),
    !anyNA(d$covered), !anyNA(d$width_ratio))
  d$scenario <- factor(d$scenario, levels = COSI_SCENARIO_LABELS)
  d
}

mean_plot_theme <- function() {
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

plot_coverage_delta <- function(csv_file, ylim = c(.80, 1.005)) {
  d <- read_mean_results(csv_file)
  alpha <- unique(d$alpha)
  stopifnot(length(alpha) == 1L)
  s <- aggregate(covered ~ scenario + delta + method, d, mean)
  names(s)[4] <- "coverage"
  order <- c("Wald", "Bootstrap", "Subsampling", "spacer", COSI_METHODS)
  s$method <- factor(s$method, levels = names(MEAN_METHOD_LABELS))
  ticks <- sort(unique(s$delta))
  tick_labels <- vapply(ticks, function(z)
    if (abs(z - 1/3) < 1e-9) "1/3" else sprintf("%.2f", z), character(1))
  ggplot2::ggplot(s, ggplot2::aes(delta, coverage, color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "#7F8993") +
    ggplot2::geom_line(linewidth = .75) + ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::scale_color_manual(values = c(MEAN_COLORS, spacer = NA_character_),
      limits = order, breaks = order, labels = c(MEAN_METHOD_LABELS, spacer = "")[order], drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = ticks, labels = tick_labels) +
    ggplot2::scale_y_continuous(breaks = seq(.8, 1, .05),
      labels = function(z) paste0(round(100 * z), "%")) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(title = "Empirical coverage", x = expression(delta), y = "Empirical coverage",
      color = NULL, subtitle = NULL, caption = NULL) + mean_plot_theme() +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
      override.aes = list(alpha = c(1, 1, 1, 0, 1, 1, 1, 1))))
}

plot_width_ratios <- function(csv_file, ylim = c(1/16, 4)) {
  d <- read_mean_results(csv_file)
  labels <- MEAN_METHOD_LABELS[names(MEAN_METHOD_LABELS) != "unimodal"]
  d <- d[d$method %in% names(labels), ]
  d$method <- factor(d$method, levels = names(labels))
  ticks <- sort(unique(d$delta))
  tick_labels <- vapply(ticks, function(z)
    if (abs(z - 1/3) < 1e-9) "1/3" else sprintf("%.2f", z), character(1))
  d$delta_label <- factor(d$delta, levels = ticks, labels = tick_labels)
  ggplot2::ggplot(d, ggplot2::aes(delta_label, width_ratio, fill = method)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "longdash", color = "#7F8993", linewidth = .5) +
    ggplot2::geom_boxplot(width = .80,
      position = ggplot2::position_dodge2(width = .86, padding = .14),
      outlier.shape = NA, color = "#28333D", linewidth = .30,
      median.linewidth = .55, staplewidth = .35) +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::scale_fill_manual(values = MEAN_COLORS[names(labels)],
      breaks = names(labels), labels = labels, drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4),
      labels = c("1/16", "1/8", "1/4", "1/2", "1", "2", "4"),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_transform(y = "log10", ylim = ylim) +
    ggplot2::labs(title = "Interval width relative to Wald", x = expression(delta),
      y = "Width / Wald width (log scale)", fill = NULL, subtitle = NULL, caption = NULL) +
    mean_plot_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}

run_mean_study <- function(core_dir = ".", output_dir = "mean_simulation_results",
    R = 500, n = 100, alpha = .05, deltas = c(.10, .15, .20, .25, .30, 1/3),
    randomize = TRUE, seed = 100, R_boot = 500, N_sub = 500, subsample_exp = 2/3) {
  load_mean_functions(core_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sim <- simulate_coverage_delta(R, n, alpha, deltas, randomize, seed,
    R_boot, N_sub, subsample_exp, checkpoint_dir = file.path(output_dir, "checkpoints"))
  csv_file <- file.path(output_dir, "mean_simulation_intervals.csv")
  write.csv(sim$results, csv_file, row.names = FALSE)

  # Both figures below are generated by reading the saved CSV, not sim.
  d <- read_mean_results(csv_file)
  coverage <- aggregate(covered ~ scenario + delta + method, d, mean)
  names(coverage)[4] <- "coverage"
  write.csv(coverage, file.path(output_dir, "empirical_coverage.csv"), row.names = FALSE)
  widths <- aggregate(width_ratio ~ scenario + delta + method, d, median)
  write.csv(widths, file.path(output_dir, "median_width_ratios.csv"), row.names = FALSE)
  unimodal <- reshape(widths[widths$method == "unimodal", c("scenario", "delta", "width_ratio")],
    idvar = "scenario", timevar = "delta", direction = "wide")
  write.csv(unimodal, file.path(output_dir, "unimodal_median_ratios.csv"), row.names = FALSE)
  write.csv(sim$design, file.path(output_dir, "block_counts.csv"), row.names = FALSE)
  saveRDS(sim, file.path(output_dir, "simulation.rds"))

  coverage_plot <- plot_coverage_delta(csv_file)
  width_plot <- plot_width_ratios(csv_file)
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(output_dir, paste0("empirical_coverage.", ext)),
      coverage_plot, width = 13, height = 8.5, dpi = 200, bg = "white")
    ggplot2::ggsave(file.path(output_dir, paste0("width_ratios.", ext)),
      width_plot, width = 13, height = 8.5, dpi = 200, bg = "white")
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(list(csv_file = csv_file, coverage = coverage,
    unimodal_medians = unimodal, coverage_plot = coverage_plot, width_plot = width_plot))
}
