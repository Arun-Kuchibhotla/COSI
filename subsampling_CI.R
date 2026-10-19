subsampling_ci <- function(
    data, estimate = mean, alpha = .05,
    b = round(NROW(data)^(2/3)),
    N = 500, beta = (9:1)/10, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  estimate <- match.fun(estimate)
  n <- NROW(data)
  theta_hat <- estimate(data)
  sizes <- round(n^beta)
  
  if (b < 1 || b >= n || b != floor(b))
    stop("b must be an integer between 1 and n - 1.")
  
  if (any(sizes < 1 | sizes >= n) || length(unique(sizes)) < 2)
    stop("beta must give at least two distinct subsample sizes below n.")
  
  rows <- function(i) {
    if (is.null(dim(data))) data[i] else data[i, , drop = FALSE]
  }
  
  # Quantiles of the subsampling distribution.
  subsample_quantiles <- function(size, probabilities, gamma = 0) {
    tmp <- replicate(
      N,
      estimate(rows(sample.int(n, size, replace = FALSE)))
    )
    
    root <- size^gamma * (tmp - theta_hat)
    grid <- seq(min(root), max(root), length.out = 100)
    G <- vapply(grid, function(z) mean(root <= z), numeric(1))
    
    vapply(
      probabilities,
      function(p) min(grid[G >= p]),
      numeric(1)
    )
  }
  
  # Estimate the convergence exponent using interquantile spacings.
  probabilities <- c(.1, .6, .2, .7, .3, .8, .4, .9)
  
  y <- vapply(sizes, function(size) {
    q <- subsample_quantiles(size, probabilities)
    mean(log(q[c(2, 4, 6, 8)] - q[c(1, 3, 5, 7)]))
  }, numeric(1))
  
  log_sizes <- log(sizes)
  gamma_hat <- -sum(
    (y - mean(y)) * (log_sizes - mean(log_sizes))
  ) / sum((log_sizes - mean(log_sizes))^2)
  
  if (is.nan(gamma_hat)) gamma_hat <- 0  # Original fallback.
  
  q <- subsample_quantiles(
    b, c(alpha/2, 1 - alpha/2), gamma = gamma_hat
  )
  
  list(
    CI = unname(theta_hat - q[c(2, 1)] / n^gamma_hat),
    gamma_hat = unname(gamma_hat),
    b = b
  )
}

# Basic bootstrap: the repository's boot()/boot.ci(type = "basic")
# implementation, with the statistic changed to the sample mean.
basic_bootstrap_mean <- function(x, alpha = .05, R_boot = 500) {
  fit <- boot::boot(x, function(data, indices) mean(data[indices]),
                    R = R_boot)
  list(CI = unname(boot::boot.ci(fit, conf = 1 - alpha,
                                 type = "basic")$basic[4:5]))
}