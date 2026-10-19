source("cosi_univariate.R")

COSI_METHODS <- c("median", "sym", "normal", "unimodal")
.cosi_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method))
    stop("method must be a single non-missing character string.")

  original <- method
  method <- sub("^cosi[[:space:]_.-]+", "", tolower(trimws(method)))

  distances <- vapply(COSI_METHODS, function(candidate) {
    aliases <- if (candidate == "sym") c("sym", "symmetric") else candidate
    min(utils::adist(method, aliases))
  }, numeric(1))

  chosen <- unname(COSI_METHODS[which.min(distances)])

  if (min(distances) > 0)
    warning(sprintf("Unknown method '%s'; using closest match '%s'.",
                    original, chosen), call. = FALSE)

  chosen
}
cosi1D <- function(data, estimate, alpha = 0.05, delta = 1/3,
                   randomize = TRUE, seed = NULL,
                   h = function(x) x, method = "median") {

  method <- .cosi_method(method)
  estimate <- match.fun(estimate)

  n <- NROW(data)
  B <- smallest_B(alpha = alpha, delta = delta) +
       (method == "unimodal")

  if (B > n)
    stop(sprintf("Need %d blocks, but only %d observations.", B, n))

  if (!is.null(seed)) set.seed(seed)

  TMP <- split(sample.int(n), rep(seq_len(B), length.out = n))

  w <- vapply(TMP, function(idx) {
    block <- if (is.null(dim(data))) {
      data[idx]
    } else {
      data[idx, , drop = FALSE]
    }
    estimate(block)
  }, numeric(1))

  pilot <- NULL
  if (method == "unimodal") {
    pilot <- w[1]
    w <- w[-1]
  }
  m <- length(w)

  if (method %in% c("median", "unimodal")) {

    cal <- randomization_quantile_dist(
      m,
      alpha = if (method == "median") alpha else 2 * alpha
    )

    # The larger binomial cutoff is conservative.
    cutoff <- cal$q_B -
              as.integer(randomize && runif(1) > cal$p)

    ans <- if (method == "median") {
      cosi_median(w, c = cutoff, alpha = alpha)
    } else {
      cosi_unimodal(w, c = cutoff, alpha = alpha, pilot = pilot)
    }

  } else {

    h_vec <- if (method == "sym") {
      h(seq_len(m))
    } else {
      qnorm(1/2 + seq_len(m) / (2 * (m + 1)))
    }

    cal <- weightedBernoulli::qweighted_bernoulli(
      p = alpha, weights = h_vec, prob = 0.5,
      lower.tail = FALSE, randomized = TRUE,
      two.sided = TRUE, draw = FALSE
    )

    # The smaller weighted cutoff is conservative.
    cutoff <- if (randomize && runif(1) < cal$tau) {
      cal$c_plus
    } else {
      cal$c_minus
    }

    ans <- if (method == "normal") {
      cosi_normal(w, c = cutoff, alpha = alpha)
    } else {
      cosi_sym(w, h = h_vec, c = cutoff, alpha = alpha)
    }
  }

  ans$B <- B             # Total data blocks, including any pilot
  ans$B_ci <- m          # Estimates used in the interval
  ans$delta <- delta
  ans$method <- method
  ans$w <- unname(w)
  if (!is.null(pilot)) ans$pilot <- unname(pilot)

  ans
}