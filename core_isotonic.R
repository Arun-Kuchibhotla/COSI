# Isotonic COSI bands. Keep this file beside cosi_univariate.R.
# The supplied univariate core is read unchanged; its source-time installs
# are skipped because only two base-R calibration functions are needed.
load_isotonic_calibration <- function(core_file = "cosi_univariate.R") {
  needed <- c("smallest_B", "randomization_quantile_dist")
  if (!file.exists(core_file)) stop("Cannot find ", core_file)
  for (expr in parse(core_file)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% needed) eval(expr, envir = .GlobalEnv)
  }
  invisible(core_file)
}

# K is floor(sqrt(n)). Replace the nearest anchor by 1/2 so the kink,
# jump and zero-derivative scenarios are always evaluated at their center.
isotonic_grid <- function(n, trim = .1) {
  K <- floor(sqrt(n))
  if (K < 3L) stop("Use at least nine observations.")
  grid <- seq(trim, 1-trim, length.out = K)
  grid[which.min(abs(grid-.5))] <- .5
  sort(grid)
}

# Univariate isotonic least squares, with correct aggregation of tied X.
# Predictions are right-continuous steps, clamped to the first/last fitted
# value outside the observed design range. block_n counts original rows,
# including repeated pairs in a bootstrap sample.
iso_fit <- function(x, y, grid) {
  n <- length(x)
  if (length(y) != n || !n || any(!is.finite(x)) || any(!is.finite(y)))
    stop("x and y must be finite vectors of the same positive length.")
  if (is.unsorted(x)) {
    o <- order(x)
    x <- x[o]; y <- y[o]
  }
  if (anyDuplicated(x)) {
    group <- cumsum(c(TRUE, diff(x) != 0))
    counts <- tabulate(group)
    means <- as.numeric(rowsum(y, group, reorder = FALSE))/counts
    y <- rep(means, counts)
  }
  fitted <- stats::isoreg(y)$yf
  starts <- c(1L, which(diff(fitted) != 0)+1L)
  counts <- diff(c(starts, n+1L))
  block_id <- rep(seq_along(starts), counts)
  index <- pmax(1L, pmin(n, findInterval(grid, x)))
  list(estimate = unname(fitted[index]),
       block_n = unname(counts[block_id[index]]),
       clamped = grid < x[1L] | grid > x[n])
}

# Monotone tightening at the anchors. Between grid[j] and grid[j+1],
# the lower endpoint is lower[j] and the upper endpoint is upper[j+1].
# Outside the anchor span one endpoint is infinite. The reported finite
# average width is restricted to [grid[1], grid[K]].
iso_monotone_band <- function(lower, upper, grid) {
  K <- length(grid)
  stopifnot(K >= 2L, length(lower) == K, length(upper) == K,
    !anyNA(lower), !anyNA(upper), all(diff(grid) > 0))
  lower <- cummax(lower)
  upper <- rev(cummin(rev(upper)))
  empty <- any(lower > upper)
  average <- if (empty) 0 else
    sum(diff(grid)*(upper[-1L]-lower[-K]))/diff(range(grid))
  list(grid = grid, lower = lower, upper = upper, empty = empty,
    avg_width = average, grid_avg_width = if(empty) 0 else mean(upper-lower),
    full_domain_area = if(empty) 0 else Inf)
}

iso_band_at <- function(band, x) {
  grid <- band$grid; K <- length(grid)
  i <- findInterval(x, grid)
  lower <- c(-Inf, band$lower)[i+1L]
  right <- pmin(K+1L, i+1L)
  exact <- match(x, grid, nomatch = 0L)
  right[exact > 0L] <- exact[exact > 0L]
  upper <- c(band$upper, Inf)[right]
  data.frame(x = x, lower = lower, upper = upper)
}

iso_order_columns <- function(W, index) {
  B <- nrow(W)
  if (index <= 0L) return(rep(-Inf, ncol(W)))
  if (index > B) return(rep(Inf, ncol(W)))
  apply(W, 2L, function(z) sort(z, partial = index)[index])
}

iso_cutoff <- function(B, alpha, randomize = TRUE, one_sided = FALSE) {
  cal <- randomization_quantile_dist(B, if(one_sided) 2*alpha else alpha)
  cal$q_B-as.integer(randomize && runif(1) > cal$p)
}

# Coordinatewise infima/suprema of the selected-coordinate cosi-m set.
# Positive-radius cones are nonempty; at radius zero, smallest-index tie
# breaking only permits the first cone. Boundary attainment is irrelevant
# because the returned intervals are closed enclosures.
iso_m_rectangle <- function(lower, upper, pilot, scale) {
  K <- length(pilot)
  radius <- pmax(abs(lower-pilot), abs(upper-pilot))/scale
  eligible <- radius > 0 | seq_len(K) == 1L
  lo <- lower; hi <- upper
  for (j in seq_len(K)) {
    other <- radius[eligible & seq_len(K) != j]
    if (length(other)) {
      lo[j] <- min(lo[j], pilot[j]-scale[j]*max(other))
      hi[j] <- max(hi[j], pilot[j]+scale[j]*max(other))
    }
  }
  list(lower = lo, upper = hi)
}

cosi_isotonic <- function(x, y, grid = isotonic_grid(length(x)), alpha = .05,
    delta = .30, randomize = TRUE, seed = 2026) {
  if (!exists("smallest_B", mode = "function") ||
      !exists("randomization_quantile_dist", mode = "function"))
    load_isotonic_calibration()
  n <- length(x); n_pilot <- floor(n/2); B <- smallest_B(alpha, delta)
  if (n-n_pilot < B) stop("The inference half is too small for B blocks.")
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if(had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  permutation <- sample.int(n)
  S0 <- permutation[seq_len(n_pilot)]
  inference <- permutation[-seq_len(n_pilot)]
  blocks <- split(inference, rep(seq_len(B), length.out = length(inference)))
  pilot <- iso_fit(x[S0], y[S0], grid)
  scale <- 1/sqrt(pilot$block_n)
  fits <- lapply(blocks, function(i) iso_fit(x[i], y[i], grid))
  W <- do.call(rbind, lapply(fits, `[[`, "estimate"))
  cm <- iso_cutoff(B, alpha, randomize, one_sided = FALSE)
  cp <- iso_cutoff(B, alpha, randomize, one_sided = TRUE)
  coordinate_lower <- iso_order_columns(W, B-cm)
  coordinate_upper <- iso_order_columns(W, cm+1L)
  rectangle <- iso_m_rectangle(coordinate_lower, coordinate_upper, pilot$estimate, scale)
  median_band <- iso_monotone_band(rectangle$lower, rectangle$upper, grid)
  distances <- sweep(abs(sweep(W, 2L, pilot$estimate, "-")), 2L, scale, "/")
  # Maximum OF coordinatewise order statistics, as in the mean example.
  radius <- max(iso_order_columns(distances, cp+1L))
  plus_band <- iso_monotone_band(pilot$estimate-radius*scale,
    pilot$estimate+radius*scale, grid)
  list(fits = list("cosi-m" = median_band, "cosi-m+" = plus_band),
    grid = grid, W = W, pilot = pilot$estimate, scale = scale,
    pilot_block_n = pilot$block_n, pilot_clamped = pilot$clamped,
    inference_clamped = do.call(rbind, lapply(fits, `[[`, "clamped")),
    B = B, n_pilot = n_pilot, blocks = blocks, S0 = S0,
    cutoff_m = cm, cutoff_mp = cp, radius_mp = radius,
    coordinate_lower = coordinate_lower, coordinate_upper = coordinate_upper)
}
