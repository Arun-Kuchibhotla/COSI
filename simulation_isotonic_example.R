# Monotone regression: simulation study and CSV-based ggplot figures.
# Keep this file beside core_isotonic.R, cosi_univariate.R, isotonic_fast.c.
# To run the complete study (500 Monte Carlo replications per cell):
#   source("simulation_isotonic_example.R")
#   study <- run_isotonic_study(core_dir=".", output_dir="isotonic_results")
# To redraw existing results without rerunning simulations:
#   export_isotonic_plots("isotonic_results/isotonic_simulation_records.csv",
#     illustrations_file="isotonic_results/illustration_bands.csv",
#     observations_file="isotonic_results/illustration_observations.csv")
# Native acceleration is optional. A base-R fallback is always available.
# source() defines functions; it does not start the simulation.

# Optional native acceleration for repeated isotonic fits. The R fallback
# calls iso_fit() and is always available. No random numbers are drawn here.

iso_load_accelerator <- function(directory = ".", compile = TRUE, quiet = TRUE) {
  if (is.loaded("iso_estimates_many_native")) return(invisible(TRUE))
  source_file <- file.path(directory, "isotonic_fast.c")
  library_file <- file.path(directory,
    paste0("isotonic_fast", .Platform$dynlib.ext))
  needs_compile <- !file.exists(library_file) ||
    (file.exists(source_file) &&
      file.info(source_file)$mtime > file.info(library_file)$mtime)
  if (compile && needs_compile && file.exists(source_file)) {
    result <- tryCatch(system2(file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", shQuote(normalizePath(source_file)),
        "-o", shQuote(normalizePath(library_file, mustWork = FALSE))),
      stdout = TRUE, stderr = TRUE), error = function(e) NULL)
    status <- attr(result, "status")
    if (!quiet && (is.null(result) || (!is.null(status) && status != 0L)))
      message("Native compilation was unavailable; using the base-R fits.")
  }
  if (file.exists(library_file)) {
    tryCatch(dyn.load(library_file), error = function(e) {
      if (!quiet) message("Native library could not be loaded; using base R.")
      NULL
    })
  }
  invisible(is.loaded("iso_estimates_many_native"))
}

iso_estimates_many <- function(x_sorted, y_sorted, indices, grid) {
  if (!is.matrix(indices)) indices <- matrix(indices, ncol = 1L)
  storage.mode(indices) <- "integer"
  K <- length(grid)
  if (isTRUE(getOption("isotonic.native", TRUE)) &&
      is.loaded("iso_estimates_many_native")) {
    return(.Call("iso_estimates_many_native", as.double(x_sorted),
      as.double(y_sorted), indices, as.double(grid)))
  }
  estimates <- vapply(seq_len(ncol(indices)), function(r) {
    i <- indices[, r]
    iso_fit(x_sorted[i], y_sorted[i], grid)$estimate
  }, numeric(K))
  matrix(estimates, nrow = K, ncol = ncol(indices))
}


# Data-generating schemes and figures for the isotonic regression study.
# This file is joined to the simulation driver when the public script is built.

ISO_SCENARIO_LABELS <- c(
  "1. Linear, Gaussian",
  "2. Curved exponential",
  "3. Cubic, flat at the center",
  "4. Heteroscedastic, beta design",
  "5. Student t with 3 degrees of freedom",
  "6. Skewed Pareto errors",
  "7. Constant function",
  "8. Plateau followed by a rise",
  "9. Unequal slopes at the center",
  "10. Increasing function with a jump",
  "11. Sparse design near the center",
  "12. Student t with 1.5 degrees of freedom"
)

ISO_METHODS <- c("Bootstrap", "Subsampling", "Subsampling_est", "HulC",
                 "cosi-m", "cosi-m+")
ISO_METHOD_LABELS <- c(
  Bootstrap = "Pairs bootstrap",
  Subsampling = "Subsampling (1/3)",
  Subsampling_est = "Subsampling (estimated rate)",
  HulC = "HulC (Bonferroni)",
  "cosi-m" = "cosi-m (rectangle)",
  "cosi-m+" = "cosi-m+"
)
ISO_METHOD_COLORS <- c(
  Bootstrap = "#666E77", Subsampling = "#DDAA33",
  Subsampling_est = "#AA6699", HulC = "#CC6677",
  "cosi-m" = "#0072B2", "cosi-m+" = "#009E73"
)

iso_truth <- function(x, scenario) {
  if (length(scenario) != 1L || is.na(scenario) ||
      !scenario %in% seq_along(ISO_SCENARIO_LABELS))
    stop("scenario must be an integer from 1 to 12.")
  switch(as.integer(scenario),
    x,
    (exp(2 * x) - 1) / (exp(2) - 1),
    .5 + 4 * (x - .5)^3,
    x,
    x,
    x,
    rep(.5, length(x)),
    4 * pmax(x - .5, 0)^2,
    .2 + ifelse(x <= .5, .4 * (x - .5), 1.6 * (x - .5)),
    .4 * x + .6 * (x >= .5),
    x,
    x
  )
}

gen_data_isotonic <- function(n, scenario) {
  if (length(n) != 1L || !is.finite(n) || n < 2L || n != floor(n))
    stop("n must be an integer of at least 2.")
  if (length(scenario) != 1L || is.na(scenario) ||
      !scenario %in% seq_along(ISO_SCENARIO_LABELS))
    stop("scenario must be an integer from 1 to 12.")
  x <- if (scenario == 4L) {
    stats::rbeta(n, 2, 3)
  } else if (scenario == 11L) {
    .5 + sample(c(-1, 1), n, replace = TRUE) * .5 * sqrt(stats::runif(n))
  } else {
    stats::runif(n)
  }
  noise <- if (scenario == 4L) {
    (.1 + .4 * x) * stats::rnorm(n)
  } else if (scenario == 5L) {
    .2 * stats::rt(n, df = 3) / sqrt(3)
  } else if (scenario == 6L) {
    a <- 2.5
    # exp(Exp(rate=a)) has Pareto tail P(Z>z)=z^(-a), z>=1.
    # This centering and scaling gives mean zero and variance one.
    z <- sqrt(a * (a - 2)) * ((a - 1) / a * exp(stats::rexp(n, a)) - 1)
    .2 * z
  } else if (scenario == 12L) {
    # Match the interquartile range of the .2*N(0,1) reference errors.
    .2 * stats::rt(n, df = 1.5) * stats::qnorm(.75) / stats::qt(.75, 1.5)
  } else {
    .2 * stats::rnorm(n)
  }
  truth <- iso_truth(x, scenario)
  list(x = x, y = truth + noise, truth = truth)
}

isotonic_scenario_table <- function() {
  data.frame(
    scenario = seq_along(ISO_SCENARIO_LABELS),
    scenario_label = ISO_SCENARIO_LABELS,
    design = c(rep("Uniform(0,1)", 3), "Beta(2,3)",
      rep("Uniform(0,1)", 6), "0.5 + 0.5*S*sqrt(U)", "Uniform(0,1)"),
    mean_function = c("x", "(exp(2*x)-1)/(exp(2)-1)",
      "0.5 + 4*(x-0.5)^3", "x", "x", "x", "0.5",
      "4*max(x-0.5,0)^2", "0.2 + 0.4*(x-0.5) left; 0.2 + 1.6*(x-0.5) right",
      "0.4*x + 0.6*I(x>=0.5)", "x", "x"),
    noise = c(rep("0.2*N(0,1)", 3), "(0.1+0.4*x)*N(0,1)",
      "0.2*t(3)/sqrt(3)", "0.2*standardized Pareto(shape=2.5)",
      rep("0.2*N(0,1)", 5), "0.2*t(1.5)*qnorm(.75)/qt(.75,1.5)"),
    feature = c("Positive constant slope", "Varying positive slope",
      "Zero first and second derivatives at the center", "Nonuniform design and changing variance",
      "Symmetric heavy tails; finite variance", "Asymmetric heavy tails; finite variance",
      "Globally flat mean", "A flat region adjacent to an increasing region",
      "Different left and right slopes", "A discontinuity at the center",
      "Design density vanishes at the center", "Finite mean but infinite variance"),
    stringsAsFactors = FALSE
  )
}

isotonic_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E3E7EB", linewidth = .25),
      panel.border = ggplot2::element_rect(fill = NA, color = "#C5CDD5", linewidth = .3),
      panel.spacing = grid::unit(.8, "lines"),
      strip.background = ggplot2::element_rect(fill = "#F1F4F7", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 9,
        color = "#25313D", margin = ggplot2::margin(t = 5, b = 5)),
      axis.text = ggplot2::element_text(color = "#35424E"),
      axis.title = ggplot2::element_text(color = "#25313D"),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8.5),
      legend.key.width = grid::unit(.36, "cm"),
      legend.key.height = grid::unit(.32, "cm"),
      legend.spacing.x = grid::unit(.08, "cm"),
      legend.spacing.y = grid::unit(.06, "cm"),
      plot.margin = ggplot2::margin(6, 7, 4, 6))
}

.isotonic_plot_data <- function(records, scenarios) {
  d <- records[records$scenario %in% scenarios & records$method %in% ISO_METHODS, ]
  if (!nrow(d)) stop("No simulation records match these scenarios and methods.")
  d$scenario_label <- factor(ISO_SCENARIO_LABELS[d$scenario],
    levels = ISO_SCENARIO_LABELS[scenarios])
  d$method <- factor(d$method, levels = ISO_METHODS)
  d
}

.isotonic_legend_order <- function() c(ISO_METHODS, "spacer1", "spacer2")
.isotonic_legend_labels <- function() {
  c(ISO_METHOD_LABELS, spacer1 = "", spacer2 = "")[.isotonic_legend_order()]
}
.isotonic_legend_colors <- function() {
  c(ISO_METHOD_COLORS, spacer1 = NA_character_, spacer2 = NA_character_)
}

.isotonic_log_coordinate <- function(ylim = NULL) {
  name <- if ("coord_transform" %in% getNamespaceExports("ggplot2"))
    "coord_transform" else "coord_trans"
  getExportedValue("ggplot2", name)(y = "log10", ylim = ylim)
}

plot_isotonic_coverage <- function(records, scenarios = 1:4, alpha = .05,
    outcome = "covered", ylim = NULL) {
  d <- .isotonic_plot_data(records, scenarios)
  if (!outcome %in% names(d)) stop("Cannot find coverage column: ", outcome)
  d$coverage_value <- d[[outcome]]
  s <- aggregate(coverage_value ~ scenario_label + n + method, d, mean)
  # Keep every empirical value, with a tighter vertical range in each panel.
  # Supplying ylim = c(lower, upper) instead gives all panels a common range.
  panel_limits <- aggregate(coverage_value ~ scenario_label, s, min)
  panel_limits$lower <- pmax(0, pmin(1 - alpha - .05,
    .05 * floor((panel_limits$coverage_value - .02 + 1e-10) / .05)))
  # A small margin below zero keeps zero-coverage markers inside the frame.
  panel_limits$lower[panel_limits$lower == 0] <- -.02
  anchors <- rbind(
    data.frame(scenario_label = panel_limits$scenario_label,
      n = min(d$n), coverage_value = panel_limits$lower),
    data.frame(scenario_label = panel_limits$scenario_label,
      n = min(d$n), coverage_value = 1.01))
  coverage_breaks <- function(limits) {
    step <- if (diff(limits) > .7) .25 else if (diff(limits) > .35) .10 else .05
    sort(unique(round(c(seq(0, 1, by = step), 1 - alpha, 1), 10)))
  }
  legend_order <- .isotonic_legend_order()
  ggplot2::ggplot(s, ggplot2::aes(n, coverage_value, color = method, group = method)) +
    ggplot2::geom_blank(data = anchors, ggplot2::aes(n, coverage_value),
      inherit.aes = FALSE) +
    ggplot2::geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "#7F8993") +
    ggplot2::geom_line(linewidth = .45) +
    ggplot2::geom_point(size = 1.3, show.legend = TRUE) +
    ggplot2::facet_wrap(~ scenario_label, ncol = 2,
      scales = if (is.null(ylim)) "free_y" else "fixed") +
    ggplot2::scale_color_manual(values = .isotonic_legend_colors(),
      limits = legend_order, breaks = legend_order, labels = .isotonic_legend_labels(),
      drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$n))) +
    ggplot2::scale_y_continuous(breaks = coverage_breaks,
      labels = function(z) paste0(round(100 * z), "%"),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = "Sample size", y = "Simultaneous coverage",
      title = NULL, subtitle = NULL, caption = NULL) +
    isotonic_plot_theme() +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
      override.aes = list(alpha = c(rep(1, 6), 0, 0))))
}

plot_isotonic_width <- function(records, scenarios = 1:4, relative = FALSE,
    ratio_ylim = c(1/16, 8)) {
  d <- .isotonic_plot_data(records, scenarios)
  d$n_label <- factor(d$n, levels = sort(unique(d$n)))
  d$width_value <- if (relative) d$width_ratio else d$avg_width
  legend_order <- .isotonic_legend_order()
  new_boxplots <- utils::packageVersion("ggplot2") >= "4.0.0"
  boxplot_style <- if (new_boxplots) {
    list(median.linewidth = .4, staplewidth = .35)
  } else {
    list(fatten = 1.75)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(n_label, width_value, fill = method)) +
    do.call(ggplot2::geom_boxplot, c(list(width = .82,
      position = ggplot2::position_dodge2(width = .88, padding = .12),
      outlier.shape = NA, color = "#28333D", linewidth = .23, show.legend = TRUE),
      boxplot_style)) +
    ggplot2::facet_wrap(~ scenario_label, ncol = 2,
      scales = if (relative) "fixed" else "free_y") +
    ggplot2::scale_fill_manual(values = .isotonic_legend_colors(),
      limits = legend_order, breaks = legend_order, labels = .isotonic_legend_labels(),
      drop = FALSE) +
    ggplot2::labs(x = "Sample size",
      y = if (relative) "Band width / HulC width (log scale)" else "Mean band width (log scale)",
      title = NULL, subtitle = NULL, caption = NULL) +
    isotonic_plot_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
      override.aes = list(alpha = c(rep(1, 6), 0, 0),
        colour = c(rep("#28333D", 6), NA, NA))))
  if (relative) {
    p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "longdash",
      color = "#7F8993", linewidth = .35) +
      ggplot2::scale_y_continuous(breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8),
        labels = c("1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8"),
        expand = ggplot2::expansion(mult = 0)) +
      .isotonic_log_coordinate(ratio_ylim)
  } else {
    p <- p + .isotonic_log_coordinate()
  }
  p
}

# Expand right-continuous step bounds before drawing a ribbon. This avoids
# visually replacing the monotone envelope with linear interpolation.
.isotonic_step_ribbon <- function(d) {
  groups <- split(d, interaction(d$n, d$method, drop = TRUE))
  out <- lapply(groups, function(z) {
    z <- z[order(z$x), ]
    if (nrow(z) < 2L) return(z)
    index <- rep(seq_len(nrow(z)), each = 2L)
    a <- z[index, ]
    a$x <- c(z$x[1], rep(z$x[-1], each = 2L), z$x[nrow(z)])
    a
  })
  do.call(rbind, out)
}

plot_isotonic_illustration <- function(bands, scenario, observations = NULL,
    method = "cosi-m+", ylim = NULL) {
  d <- bands[bands$scenario == scenario, ]
  if (!nrow(d) || !method %in% d$method)
    stop("Illustration data are missing the requested scenario or method.")
  selected <- d[d$method == method, ]
  selected <- selected[order(selected$n, selected$x), ]
  selected$n_label <- factor(paste0("n = ", selected$n),
    levels = paste0("n = ", sort(unique(selected$n))))
  ribbon <- .isotonic_step_ribbon(selected)
  reference <- d[d$method == "HulC", ]
  reference$n_label <- factor(paste0("n = ", reference$n), levels = levels(selected$n_label))
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = ribbon,
      ggplot2::aes(x, ymin = lower, ymax = upper, fill = method), alpha = .25) +
    ggplot2::geom_line(data = selected,
      ggplot2::aes(x, truth, color = "Truth"), linewidth = .6) +
    ggplot2::geom_step(data = selected,
      ggplot2::aes(x, estimate, color = "Isotonic fit"), linewidth = .4, direction = "hv") +
    ggplot2::facet_wrap(~ n_label, ncol = 2) +
    ggplot2::scale_fill_manual(values = ISO_METHOD_COLORS[method],
      labels = ISO_METHOD_LABELS[method]) +
    ggplot2::scale_color_manual(values = c("Truth" = "#111111",
      "Isotonic fit" = "#0072B2", "HulC (Bonferroni)" = "#A45161")) +
    ggplot2::labs(x = "Covariate", y = expression(f[0](x)),
      title = ISO_SCENARIO_LABELS[scenario], subtitle = NULL, caption = NULL) + isotonic_plot_theme()
  if (nrow(reference)) {
    p <- p + ggplot2::geom_step(data = reference,
      ggplot2::aes(x, lower, color = "HulC (Bonferroni)"),
      direction = "hv", linewidth = .35, linetype = "dashed") +
      ggplot2::geom_step(data = reference,
        ggplot2::aes(x, upper, color = "HulC (Bonferroni)"),
        direction = "hv", linewidth = .35, linetype = "dashed")
  }
  if (!is.null(observations)) {
    points <- observations[observations$scenario == scenario, ]
    points$n_label <- factor(paste0("n = ", points$n), levels = levels(selected$n_label))
    p <- p + ggplot2::geom_point(data = points, ggplot2::aes(x, y),
      size = .35, alpha = .15, color = "#657280", inherit.aes = FALSE)
  }
  p + ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1),
      color = ggplot2::guide_legend(order = 2))
}


# Full-sample comparisons for the isotonic regression bands.
# The core file supplies iso_fit() and iso_monotone_band().

# Preserve the 100-point inverse-CDF calculation used in the HulC code.
iso_rate_quantiles <- function(values, probs) {
  xgrid <- seq(min(values), max(values), length.out = 100L)
  cdf <- findInterval(xgrid, sort(values)) / length(values)
  vapply(probs, function(p) xgrid[which(cdf >= p)[1L]], numeric(1))
}

iso_benchmarks <- function(x, y, grid, scale, alpha = .05,
    R_boot = 500L, N_sub = 500L, N_rate = 200L,
    subsample_exp = 2/3, beta = c(.9, .75, .6, .45)) {
  n <- length(x)
  K <- length(grid)
  stopifnot(length(y) == n, length(scale) == K,
    all(is.finite(scale) & scale > 0),
    length(alpha) == 1L, alpha > 0, alpha < 1,
    R_boot >= 2L, N_sub >= 2L, N_rate >= 2L)

  # Sort once. Sorted sampled indices then preserve the covariate ordering.
  ordering <- order(x)
  x <- x[ordering]
  y <- y[ordering]
  estimates_at <- function(indices)
    iso_estimates_many(x, y, indices, grid)

  center <- iso_fit(x, y, grid)$estimate
  standardized_errors <- function(estimates)
    abs(sweep(sweep(estimates, 1L, center, "-"), 1L, scale, "/"))
  sup_distances <- function(estimates)
    apply(standardized_errors(estimates), 2L, max)
  quantile_root <- function(root)
    unname(stats::quantile(root, probs = 1 - alpha, type = 6))

  make_band <- function(critical, exponent, sample_size) {
    radius <- critical * (sample_size / n)^exponent
    band <- iso_monotone_band(center - radius * scale,
      center + radius * scale, grid)
    band$center <- center
    band$critical <- critical
    band$radius <- radius
    band$rate <- exponent
    band$subsample_size <- sample_size
    band
  }

  # Ordinary pairs bootstrap. The sample contains n pairs, including repeats.
  bootstrap_indices <- replicate(R_boot, sort(sample.int(n, n, replace = TRUE)))
  bootstrap_root <- sup_distances(estimates_at(bootstrap_indices))
  bootstrap <- make_band(quantile_root(bootstrap_root), 0, n)

  # Use the same construction subsamples for the two subsampling methods.
  b <- as.integer(round(n^subsample_exp))
  if (b < 8L || b >= n)
    stop("The construction subsample size must be at least 8 and less than n.")
  subsample_indices <- replicate(N_sub, sort(sample.int(n, b, replace = FALSE)))
  subsample_errors <- standardized_errors(estimates_at(subsample_indices))
  subsample_root <- apply(subsample_errors, 2L, max)
  critical <- quantile_root(subsample_root)
  subsampling <- make_band(critical, 1/3, b)

  # Apply the scalar interquantile-range regression separately at every anchor.
  # Thus a flat or cubic region does not impose its rate on the whole curve.
  rate_sizes <- sort(unique(as.integer(round(n^beta))))
  rate_sizes <- rate_sizes[rate_sizes >= 8L & rate_sizes < n]
  rate_coordinate <- which.min(abs(grid - .5))
  probs <- c(.1, .6, .2, .7, .3, .8, .4, .9)
  odd <- c(1L, 3L, 5L, 7L)
  even <- c(2L, 4L, 6L, 8L)
  mean_log_width <- matrix(NA_real_, nrow = K, ncol = length(rate_sizes),
    dimnames = list(as.character(grid), as.character(rate_sizes)))

  for (j in seq_along(rate_sizes)) {
    size <- rate_sizes[j]
    rate_indices <- replicate(N_rate, sort(sample.int(n, size, replace = FALSE)))
    values <- sweep(estimates_at(rate_indices), 1L, center, "-")
    for (k in seq_len(K)) {
      quantiles <- iso_rate_quantiles(values[k, ], probs)
      widths <- quantiles[even] - quantiles[odd]
      if (all(is.finite(widths) & widths > 0))
        mean_log_width[k, j] <- mean(log(widths))
    }
  }

  gammahat_vector <- rep(0, K)
  rate_status_vector <- rep("ok", K)
  if (length(rate_sizes) < 2L) {
    rate_status_vector[] <- "insufficient_rate_sizes"
  } else {
    log_size <- log(rate_sizes)
    for (k in seq_len(K)) {
      if (any(!is.finite(mean_log_width[k, ]))) {
        rate_status_vector[k] <- "degenerate_interquantile_widths"
      } else {
        gammahat_vector[k] <- -sum((mean_log_width[k, ] -
          mean(mean_log_width[k, ])) * (log_size - mean(log_size))) /
          sum((log_size - mean(log_size))^2)
        if (!is.finite(gammahat_vector[k])) {
          gammahat_vector[k] <- 0
          rate_status_vector[k] <- "nonfinite_rate"
        }
      }
    }
  }
  gammahat <- gammahat_vector[rate_coordinate]
  rate_status <- if (all(rate_status_vector == "ok")) "ok" else "coordinate_fallback"

  # Rescale each coordinate before forming the maximum. Finite estimated
  # exponents are retained as obtained, without truncation.
  estimated_errors <- sweep(subsample_errors, 1L, (b/n)^gammahat_vector, "*")
  estimated_errors[subsample_errors == 0] <- 0
  estimated_root <- apply(estimated_errors, 2L, max)
  critical_est <- quantile_root(estimated_root)
  subsampling_est <- make_band(critical_est, 0, n)
  subsampling_est$rate <- gammahat_vector
  subsampling_est$subsample_size <- b
  subsampling_est$rate_status <- rate_status
  subsampling_est$rate_status_vector <- rate_status_vector

  # Coordinatewise HulC with one common full-sample partition.
  bH <- as.integer(ceiling(log2(2 * K / alpha)))
  if (bH > n) stop("Too few observations for the HulC partition.")
  permutation <- sample.int(n)
  blocks <- split(permutation, rep(seq_len(bH), length.out = n))
  block_estimates <- vapply(blocks, function(indices)
    iso_fit(x[sort(indices)], y[sort(indices)], grid)$estimate, numeric(K))
  hulc <- iso_monotone_band(apply(block_estimates, 1L, min),
    apply(block_estimates, 1L, max), grid)
  hulc$B <- bH

  list(fits = list(Bootstrap = bootstrap, Subsampling = subsampling,
      Subsampling_est = subsampling_est, HulC = hulc),
    gammahat = gammahat, gammahat_vector = gammahat_vector,
    gammahat_min = min(gammahat_vector), gammahat_median = median(gammahat_vector),
    gammahat_max = max(gammahat_vector),
    rate_status = rate_status, rate_status_vector = rate_status_vector,
    rate_coordinate = rate_coordinate, rate_t = grid[rate_coordinate],
    rate_sizes = rate_sizes, rate_mean_log_width = mean_log_width,
    b = b, bH = bH)
}


# Simulation, saved records and figures for the isotonic regression study.
# Sourcing the public script defines functions; run_isotonic_study() runs it.

.iso_source_hashes <- function(files) {
  files <- files[file.exists(files)]
  hashes <- tools::md5sum(files)
  names(hashes) <- basename(names(hashes))
  hashes
}

.iso_band_covered <- function(band, truth) {
  !band$empty && all(band$lower <= truth & truth <= band$upper)
}

.iso_cell <- function(scenario, n, settings, output_dir, resume) {
  checkpoint <- file.path(output_dir, "checkpoints",
    sprintf("scenario_%02d_n_%04d.rds", scenario, n))
  signature <- list(scenario = scenario, n = n, settings = settings)
  if (resume && file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (!identical(saved$signature, signature))
      stop("Checkpoint settings differ: ", checkpoint,
        ". Use a new output directory or resume=FALSE.")
    message("Reused scenario ", scenario, ", n=", n, ".")
    return(saved)
  }

  grid <- isotonic_grid(n)
  K <- length(grid)
  center <- which.min(abs(grid - .5))
  truth <- iso_truth(grid, scenario)
  records <- vector("list", settings$R)
  grid_counts <- matrix(0, K, 11L, dimnames = list(NULL,
    c("split_below", "split_above", "split_equal", "pilot_below",
      "pilot_above", "pilot_equal", "full_below", "full_above",
      "full_equal", "pilot_clamped", "inference_clamped")))
  pilot_counts <- matrix(NA_real_, settings$R, K)
  rate_estimates <- matrix(NA_real_, settings$R, K)
  rate_fallbacks <- matrix(FALSE, settings$R, K)
  example <- NULL

  for (replicate in seq_len(settings$R)) {
    rep_seed <- settings$seed + 1000000 * scenario + 1000 * n + replicate
    set.seed(rep_seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
      sample.kind = "Rejection")
    dat <- gen_data_isotonic(n, scenario)
    cosi <- cosi_isotonic(dat$x, dat$y, grid = grid,
      alpha = settings$alpha, delta = settings$delta,
      randomize = settings$randomize, seed = rep_seed + 50000000)

    # Benchmarks share the data and pilot scales; their resampling stream is
    # separate from both data generation and the randomized COSI partition.
    set.seed(rep_seed + 100000000)
    benchmark <- iso_benchmarks(dat$x, dat$y, grid, cosi$scale,
      alpha = settings$alpha, R_boot = settings$R_boot,
      N_sub = settings$N_sub, N_rate = settings$N_rate,
      subsample_exp = settings$subsample_exp, beta = settings$beta)
    rate_estimates[replicate, ] <- benchmark$gammahat_vector
    rate_fallbacks[replicate, ] <- benchmark$rate_status_vector != "ok"
    fits <- c(benchmark$fits, cosi$fits)[ISO_METHODS]
    if (any(vapply(fits, is.null, logical(1))))
      stop("A requested method is missing from the fitted bands.")
    reference_width <- fits$HulC$avg_width
    if (!is.finite(reference_width) || reference_width <= 0)
      stop("The paired HulC reference width is not finite and positive.")
    full <- iso_fit(dat$x, dat$y, grid)

    records[[replicate]] <- do.call(rbind, lapply(ISO_METHODS, function(method) {
      band <- fits[[method]]
      data.frame(scenario = scenario,
        scenario_label = ISO_SCENARIO_LABELS[scenario], n = n,
        replicate = replicate, method = method, alpha = settings$alpha,
        delta = settings$delta, K = K,
        B = if (method %in% c("cosi-m", "cosi-m+")) cosi$B else
          if (method == "HulC") benchmark$bH else NA_integer_,
        n_pilot = cosi$n_pilot,
        b = if (method %in% c("Subsampling", "Subsampling_est"))
          benchmark$b else NA_integer_,
        cutoff_m = cosi$cutoff_m, cutoff_mp = cosi$cutoff_mp,
        covered = .iso_band_covered(band, truth),
        avg_width = band$avg_width,
        width_ratio = band$avg_width / reference_width,
        hulc_width = reference_width, grid_avg_width = band$grid_avg_width,
        empty = band$empty,
        gamma_hat = if (method == "Subsampling_est") benchmark$gammahat else NA_real_,
        rate_status = if (method == "Subsampling_est") benchmark$rate_status else "",
        gamma_hat_min = if(method == "Subsampling_est") min(benchmark$gammahat_vector) else NA_real_,
        gamma_hat_median = if(method == "Subsampling_est") median(benchmark$gammahat_vector) else NA_real_,
        gamma_hat_max = if(method == "Subsampling_est") max(benchmark$gammahat_vector) else NA_real_,
        pilot_block_n_center = cosi$pilot_block_n[center],
        center_split_above = mean(cosi$W[, center] > truth[center]),
        center_split_below = mean(cosi$W[, center] < truth[center]),
        center_split_equal = mean(cosi$W[, center] == truth[center]),
        center_pilot_error = cosi$pilot[center] - truth[center],
        center_full_error = full$estimate[center] - truth[center],
        pilot_clamp_fraction = mean(cosi$pilot_clamped),
        inference_clamp_fraction = mean(cosi$inference_clamped),
        seed = rep_seed, stringsAsFactors = FALSE)
    }))

    differences <- sweep(cosi$W, 2L, truth, "-")
    grid_counts[, "split_below"] <- grid_counts[, "split_below"] + colSums(differences < 0)
    grid_counts[, "split_above"] <- grid_counts[, "split_above"] + colSums(differences > 0)
    grid_counts[, "split_equal"] <- grid_counts[, "split_equal"] + colSums(differences == 0)
    grid_counts[, "pilot_below"] <- grid_counts[, "pilot_below"] + (cosi$pilot < truth)
    grid_counts[, "pilot_above"] <- grid_counts[, "pilot_above"] + (cosi$pilot > truth)
    grid_counts[, "pilot_equal"] <- grid_counts[, "pilot_equal"] + (cosi$pilot == truth)
    grid_counts[, "full_below"] <- grid_counts[, "full_below"] + (full$estimate < truth)
    grid_counts[, "full_above"] <- grid_counts[, "full_above"] + (full$estimate > truth)
    grid_counts[, "full_equal"] <- grid_counts[, "full_equal"] + (full$estimate == truth)
    grid_counts[, "pilot_clamped"] <- grid_counts[, "pilot_clamped"] + cosi$pilot_clamped
    grid_counts[, "inference_clamped"] <- grid_counts[, "inference_clamped"] +
      colSums(cosi$inference_clamped)
    pilot_counts[replicate, ] <- cosi$pilot_block_n

    # The first replication is selected in advance, without searching for
    # an attractive or successfully covering example.
    if (replicate == 1L) example <- list(scenario = scenario, n = n,
      seed = rep_seed, x = dat$x, y = dat$y, grid = grid,
      truth = truth, full = full, fits = fits, cosi = cosi,
      benchmark = benchmark)
  }

  diagnostics <- data.frame(scenario = scenario,
    scenario_label = ISO_SCENARIO_LABELS[scenario], n = n,
    t = grid, truth = truth, B = cosi$B,
    split_below = grid_counts[, "split_below"] / (settings$R * cosi$B),
    split_above = grid_counts[, "split_above"] / (settings$R * cosi$B),
    split_equal = grid_counts[, "split_equal"] / (settings$R * cosi$B),
    pilot_below = grid_counts[, "pilot_below"] / settings$R,
    pilot_above = grid_counts[, "pilot_above"] / settings$R,
    pilot_equal = grid_counts[, "pilot_equal"] / settings$R,
    full_below = grid_counts[, "full_below"] / settings$R,
    full_above = grid_counts[, "full_above"] / settings$R,
    full_equal = grid_counts[, "full_equal"] / settings$R,
    median_pilot_block_n = apply(pilot_counts, 2L, stats::median),
    median_gamma_hat = apply(rate_estimates, 2L, stats::median),
    q10_gamma_hat = apply(rate_estimates, 2L, stats::quantile, probs=.1),
    q90_gamma_hat = apply(rate_estimates, 2L, stats::quantile, probs=.9),
    rate_fallback_fraction = colMeans(rate_fallbacks),
    pilot_clamp_fraction = grid_counts[, "pilot_clamped"] / settings$R,
    inference_clamp_fraction = grid_counts[, "inference_clamped"] /
      (settings$R * cosi$B), stringsAsFactors = FALSE)
  result <- list(signature = signature,
    records = do.call(rbind, records), diagnostics = diagnostics, example = example)
  saveRDS(result, checkpoint)
  message("Completed scenario ", scenario, ", n=", n, ": ", settings$R, " replications.")
  result
}

simulate_isotonic_study <- function(R = 500L,
    n_grid = c(100L, 250L, 500L, 1000L), scenarios = 1:12,
    alpha = .05, delta = .30, randomize = TRUE,
    R_boot = 500L, N_sub = 500L, N_rate = 200L,
    subsample_exp = 2/3, beta = c(.9, .75, .6, .45), seed = 100L,
    workers = 6L, output_dir = "isotonic_results", resume = TRUE,
    source_files = c("core_isotonic.R", "cosi_univariate.R",
      "simulation_isotonic_example.R", "isotonic_fast.c")) {
  stopifnot(length(R) == 1L, R >= 1L, R == floor(R),
    all(n_grid >= 20L & n_grid == floor(n_grid)),
    all(scenarios %in% seq_along(ISO_SCENARIO_LABELS)),
    length(seed) == 1L, is.finite(seed), seed >= 0, seed == floor(seed),
    length(workers) == 1L, workers >= 1L, workers == floor(workers))
  if (seed + 1000000 * max(scenarios) + 1000 * max(n_grid) + R +
      100000000 > .Machine$integer.max)
    stop("The requested settings exceed the integer seed range.")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  dir.create(file.path(output_dir, "checkpoints"), recursive = TRUE,
    showWarnings = FALSE)
  settings <- list(R = as.integer(R), n_grid = as.integer(n_grid),
    scenarios = as.integer(scenarios), alpha = alpha, delta = delta,
    randomize = randomize, R_boot = as.integer(R_boot), N_sub = as.integer(N_sub),
    N_rate = as.integer(N_rate), subsample_exp = subsample_exp, beta = beta,
    seed = seed, rng_kind = RNGkind(), source_md5 = .iso_source_hashes(source_files),
    R_version = R.version.string,
    ggplot2_version = as.character(utils::packageVersion("ggplot2")),
    native_acceleration = is.loaded("iso_estimates_many_native"),
    grid_trim = .1, grid_rule = "floor(sqrt(n)); equally spaced on [.1,.9]; include .5",
    width_reference = "HulC", width_domain = c(.1, .9))
  cells <- expand.grid(scenario = as.integer(scenarios), n = as.integer(n_grid))
  task <- function(i) .iso_cell(cells$scenario[i], cells$n[i], settings,
    output_dir = output_dir, resume = resume)
  if (workers > 1L && .Platform$OS.type != "windows") {
    fitted <- parallel::mclapply(seq_len(nrow(cells)), task,
      mc.cores = min(workers, nrow(cells)), mc.set.seed = FALSE)
  } else {
    fitted <- lapply(seq_len(nrow(cells)), task)
  }
  failed <- vapply(fitted, inherits, logical(1), "try-error")
  if (any(failed)) stop("Simulation cell failed: ",
    paste(as.character(fitted[failed]), collapse = "\n"))
  records <- do.call(rbind, lapply(fitted, `[[`, "records"))
  diagnostics <- do.call(rbind, lapply(fitted, `[[`, "diagnostics"))
  rownames(records) <- rownames(diagnostics) <- NULL
  list(records = records, diagnostics = diagnostics,
    examples = lapply(fitted, `[[`, "example"), settings = settings)
}

.iso_rate_summary <- function(records) {
  d <- records[records$method == "Subsampling_est", ]
  groups <- split(d, interaction(d$scenario, d$n, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) data.frame(
    scenario = z$scenario[1L], scenario_label = z$scenario_label[1L], n = z$n[1L],
    median_rate = stats::median(z$gamma_hat),
    q10_rate = unname(stats::quantile(z$gamma_hat, .1)),
    q90_rate = unname(stats::quantile(z$gamma_hat, .9)),
    minimum_rate = min(z$gamma_hat), maximum_rate = max(z$gamma_hat),
    fallback_fraction = mean(z$rate_status != "ok"), stringsAsFactors = FALSE)))
}

.iso_illustration_data <- function(examples) {
  bands <- vector("list", length(examples))
  observations <- vector("list", length(examples))
  for (i in seq_along(examples)) {
    example <- examples[[i]]
    dense <- sort(unique(c(seq(.1, .9, length.out = 501L), example$grid,
      example$grid[-length(example$grid)] + 1e-10,
      .5 - 1e-10, .5, .5 + 1e-10)))
    truth <- iso_truth(dense, example$scenario)
    estimate <- iso_fit(example$x, example$y, dense)$estimate
    bands[[i]] <- do.call(rbind, lapply(ISO_METHODS, function(method) {
      bound <- iso_band_at(example$fits[[method]], dense)
      data.frame(scenario = example$scenario, n = example$n,
        seed = example$seed, method = method, x = dense, truth = truth,
        estimate = estimate, lower = bound$lower, upper = bound$upper,
        stringsAsFactors = FALSE)
    }))
    observations[[i]] <- data.frame(scenario = example$scenario,
      n = example$n, seed = example$seed, x = example$x, y = example$y,
      truth = iso_truth(example$x, example$scenario))
  }
  list(bands = do.call(rbind, bands), observations = do.call(rbind, observations))
}

# All main figures are drawn from the saved CSV. The 12 scenarios are split
# into three groups of four, so each main figure has a legible 2-by-2 layout.
export_isotonic_plots <- function(csv_file,
    output_dir = dirname(csv_file), illustrations_file = NULL,
    observations_file = NULL, png = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  if (!capabilities("cairo")) stop("PDF export requires Cairo support.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  records <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
  groups <- list(1:4, 5:8, 9:12)
  for (g in seq_along(groups)) {
    scenarios <- intersect(groups[[g]], unique(records$scenario))
    if (!length(scenarios)) next
    plots <- list(
      empirical_coverage = plot_isotonic_coverage(records, scenarios,
        alpha = unique(records$alpha)),
      band_widths = plot_isotonic_width(records, scenarios, relative = FALSE),
      width_ratios = plot_isotonic_width(records, scenarios, relative = TRUE))
    for (name in names(plots)) {
      stem <- paste0(name, "_group_", g)
      ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), plots[[name]],
        device = grDevices::cairo_pdf, width = 7.2, height = 5.6, units = "in", bg = "white")
      if (png) ggplot2::ggsave(file.path(output_dir, paste0(stem, ".png")), plots[[name]],
        width = 7.2, height = 5.6, units = "in", dpi = 200, bg = "white")
    }
  }
  if (!is.null(illustrations_file)) {
    bands <- utils::read.csv(illustrations_file, stringsAsFactors = FALSE)
    observations <- if (is.null(observations_file)) NULL else
      utils::read.csv(observations_file, stringsAsFactors = FALSE)
    if (!is.null(observations)) observations <- observations[
      observations$x >= .1 & observations$x <= .9, ]
    scenarios <- sort(unique(bands$scenario))
    methods <- c("cosi-m+", "cosi-m")
    stems <- c("cosi_mplus_band_illustrations", "cosi_m_rectangle_band_illustrations")
    for (m in seq_along(methods)) {
      grDevices::cairo_pdf(file.path(output_dir, paste0(stems[m], ".pdf")),
        width = 7.2, height = 5.6, onefile = TRUE)
      tryCatch({
        for (scenario in scenarios) print(plot_isotonic_illustration(
          bands, scenario, observations, method = methods[m]))
      }, finally = grDevices::dev.off())
    }
    if (png) for (scenario in intersect(c(1L, 4L, 7L, 10L), scenarios)) {
      plot <- plot_isotonic_illustration(bands, scenario, observations, method = "cosi-m+")
      ggplot2::ggsave(file.path(output_dir,
        sprintf("band_illustration_scenario_%02d.png", scenario)), plot,
        width = 7.2, height = 5.6, units = "in", dpi = 200, bg = "white")
    }
  }
  invisible(records)
}

run_isotonic_study <- function(core_dir = ".", output_dir = "isotonic_results",
    R = 500L, n_grid = c(100L, 250L, 500L, 1000L), scenarios = 1:12,
    alpha = .05, delta = .30, randomize = TRUE, R_boot = 500L,
    N_sub = 500L, N_rate = 200L, subsample_exp = 2/3,
    beta = c(.9, .75, .6, .45), seed = 100L, workers = 6L,
    resume = TRUE, make_plots = TRUE,
    source_files = file.path(core_dir, c("core_isotonic.R", "cosi_univariate.R",
      "simulation_isotonic_example.R", "isotonic_fast.c"))) {
  if (!exists("cosi_isotonic", mode = "function"))
    source(file.path(core_dir, "core_isotonic.R"))
  load_isotonic_calibration(file.path(core_dir, "cosi_univariate.R"))
  iso_load_accelerator(core_dir, compile = TRUE, quiet = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  study <- simulate_isotonic_study(R = R, n_grid = n_grid, scenarios = scenarios,
    alpha = alpha, delta = delta, randomize = randomize, R_boot = R_boot,
    N_sub = N_sub, N_rate = N_rate, subsample_exp = subsample_exp,
    beta = beta, seed = seed, workers = workers, output_dir = output_dir,
    resume = resume, source_files = source_files)
  records <- study$records
  coverage <- stats::aggregate(covered ~ scenario + scenario_label + n + method,
    records, mean)
  names(coverage)[names(coverage) == "covered"] <- "coverage"
  widths <- stats::aggregate(cbind(avg_width, width_ratio, grid_avg_width) ~
    scenario + scenario_label + n + method, records, stats::median)
  status <- as.data.frame(table(scenario = records$scenario,
    n = records$n, method = records$method, empty = records$empty))
  status <- status[status$Freq > 0, ]
  rate_summary <- .iso_rate_summary(records)
  illustration <- .iso_illustration_data(study$examples)
  csv_file <- file.path(output_dir, "isotonic_simulation_records.csv")
  write <- function(value, name) utils::write.csv(value,
    file.path(output_dir, name), row.names = FALSE)
  write(records, "isotonic_simulation_records.csv")
  write(coverage, "empirical_coverage.csv")
  write(widths, "median_band_widths.csv")
  write(study$diagnostics, "grid_diagnostics.csv")
  write(rate_summary, "subsampling_rate_summary.csv")
  write(status, "empty_band_counts.csv")
  write(isotonic_scenario_table(), "scenario_settings.csv")
  write(illustration$bands, "illustration_bands.csv")
  write(illustration$observations, "illustration_observations.csv")
  saveRDS(study, file.path(output_dir, "simulation.rds"))
  saveRDS(study$settings, file.path(output_dir, "settings.rds"))
  dput(study$settings, file = file.path(output_dir, "settings.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  write(data.frame(file = names(study$settings$source_md5),
    md5 = unname(study$settings$source_md5)), "source_md5.csv")
  if (make_plots) export_isotonic_plots(csv_file, output_dir = output_dir,
    illustrations_file = file.path(output_dir, "illustration_bands.csv"),
    observations_file = file.path(output_dir, "illustration_observations.csv"))
  invisible(list(csv_file = csv_file, coverage = coverage, widths = widths,
    diagnostics = study$diagnostics, rate_summary = rate_summary,
    settings = study$settings))
}
