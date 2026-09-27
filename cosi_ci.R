# Main workflow. source("cosi_ci.R") loads the three neighboring modules.
# Calibration and interval formulas are kept out of this file.
.cosi_load_env <- environment()
.cosi_source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
}, character(1))
.cosi_source_files <- .cosi_source_files[nzchar(.cosi_source_files)]
.cosi_source_dir <- if (length(.cosi_source_files))
  dirname(normalizePath(tail(.cosi_source_files,1L),mustWork=TRUE)) else getwd()
for (.cosi_module in c("cosi_smallest_B.R","cosi_calibration.R","cosi_intervals.R")) {
  sys.source(file.path(.cosi_source_dir,.cosi_module),envir=.cosi_load_env)
}
rm(.cosi_load_env,.cosi_source_files,.cosi_source_dir,.cosi_module)

# Validate gap as a nonnegative integer number of omitted rows between batches.
.cosi_check_gap <- function(gap) {
  if (!is.numeric(gap) || length(gap) != 1 || !is.finite(gap) ||
      gap < 0 || gap != floor(gap)) stop("gap must be a nonnegative integer")
  invisible(NULL)
}

# Validate a shared random plan; ordered gap splits require an identity permutation.
.cosi_check_randomness <- function(randomness, nn, B, signed_rank, gap = 0) {
  if (is.null(randomness)) return(invisible(NULL))
  required <- c("permutation","cutoff_uniform","uniforms")
  if (!is.list(randomness) || !setequal(names(randomness),required) || length(randomness) != 3) {
    stop("randomness needs permutation, cutoff_uniform, uniforms")
  }
  p <- randomness$permutation
  if (!is.numeric(p) || length(p) != nn || any(!is.finite(p)) ||
      !all(sort(p) == seq_len(nn))) stop("randomness$permutation must permute 1,...,n")
  if (gap > 0 && !all(p == seq_len(nn))) stop("gap>0 requires the original row order in randomness$permutation")
  u <- randomness$cutoff_uniform
  if (!is.numeric(u) || length(u) != 1 || !is.finite(u) || u < 0 || u >= 1) {
    stop("randomness$cutoff_uniform must lie in [0,1)")
  }
  if (signed_rank) {
    if (is.null(randomness$uniforms)) stop("randomness$uniforms must supply B tie-breakers")
    .cosi_tie_breakers(B,randomness$uniforms)
  }
  invisible(NULL)
}

# Generate reusable data-independent random inputs. Positive gaps preserve row
# order; seed controls this call locally. Share the resulting plan across R/Python.
make_cosi_randomness <- function(n, B, seed = NULL, gap = 0) {
  restore <- .cosi_seed_scope(seed)
  on.exit(restore(),add=TRUE)
  .cosi_check_gap(gap)
  if (!is.numeric(n) || !is.numeric(B) || length(n) != 1 || length(B) != 1 ||
      any(!is.finite(c(n,B))) || n < 1 || B < 1 || n != floor(n) || B != floor(B)) {
    stop("n and B must be positive integers")
  }
  p <- if (gap > 0) seq_len(n) else sample.int(n)
  list(permutation=p,cutoff_uniform=runif(1),uniforms=.cosi_tie_breakers(B))
}

# Compute balanced one-based row groups and omitted gap indices. At gap=0,
# preserve the original shuffled-split convention for allocating remainder rows.
.cosi_split_rows <- function(nn, batches, gap, permutation) {
  n_used <- nn-(batches-1)*gap
  if (n_used < batches) stop(sprintf(
    "Not enough rows: need at least %.0f for %.0f nonempty batches with gap %.0f; have %d",
    batches+(batches-1)*gap,batches,gap,nn))
  sizes <- rep(floor(n_used/batches),batches)
  remainder <- n_used %% batches
  if (remainder > 0) sizes[seq_len(remainder)+1L] <- sizes[seq_len(remainder)+1L]+1
  names(sizes) <- as.character(0:(batches-1))
  rows <- vector("list",batches)
  omitted <- vector("list",max(0,batches-1))
  at <- 1
  for (j in seq_len(batches)) {
    rows[[j]] <- permutation[seq.int(at,length.out=sizes[j])]
    at <- at+sizes[j]
    if (j < batches && gap > 0) {
      omitted[[j]] <- permutation[seq.int(at,length.out=gap)]
      at <- at+gap
    }
  }
  names(rows) <- names(sizes)
  list(batch_rows=rows,batch_sizes=sizes,
       gap_rows=as.integer(unlist(omitted,use.names=FALSE)),n_used=n_used)
}

# Choose B and its cutoff, split data into B (or B+1) balanced batches with
# optional gaps, compute scalar estimates, and call exactly one constructor.
# The first seven arguments are the requested interface. Advanced options retain
# custom scores, conservative cutoff selection, resource limits, and shared draws.
cosi_ci <- function(data, alpha = 0.05, gamma = 1/3, estimate,
                    method = c("median","symmetry","normality","unimodal-star"),
                    gap = 0, seed = NULL, h = function(i) i, randomized = TRUE,
                    control = list(), randomness = NULL) {
  restore <- .cosi_seed_scope(seed)
  on.exit(restore(),add=TRUE)
  if (!is.function(estimate)) stop("estimate must be a function")
  .cosi_check_alpha(alpha)
  .cosi_check_gap(gap)
  if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) ||
      gamma <= 0 || gamma >= 0.5) stop("gamma must lie strictly between 0 and 1/2")
  if (!is.logical(randomized) || length(randomized) != 1 || is.na(randomized)) {
    stop("randomized must be TRUE or FALSE")
  }
  if (identical(method,"symmetric")) method <- "symmetry"
  if (identical(method,"unimodal")) method <- "unimodal-star"
  method <- match.arg(method,c("median","symmetry","normality","unimodal-star"))
  control <- .cosi_ci_control(control)
  data <- as.matrix(data)
  nn <- nrow(data)
  if (is.null(nn) || nn < 1 || ncol(data) < 1) stop("data must have rows and columns")

  # 1. Select B and calibrate the cutoff independently of the data values.
  split_calibration <- smallest_B_distribution(alpha,gamma)
  B <- split_calibration$B
  # K is the actual number of batches: B+1 for an anchor, otherwise B.
  K <- B+as.integer(method == "unimodal-star")
  batches <- K
  if (nn < batches+(batches-1)*gap) stop(sprintf(
    "Not enough rows: need at least %.0f for %.0f batches with gap %.0f; have %d",
    batches+(batches-1)*gap,batches,gap,nn))
  signed_rank <- method %in% c("symmetry","normality")
  .cosi_check_randomness(randomness,nn,B,signed_rank,gap)
  cutoff_uniform <- if (is.null(randomness)) NULL else randomness$cutoff_uniform
  h_values <- NULL
  sign_scores <- FALSE
  if (method == "median") {
    law <- split_calibration
  } else if (method == "unimodal-star") {
    law <- binomial_cutoff_distribution(B,level=alpha)
  } else {
    h_values <- .cosi_score_values(B,method,h)
    sign_scores <- all(h_values[-1L] == 1)
    if (!sign_scores && B*(B+1)/2 > control$max_walsh) stop("Number of Walsh pairs exceeds control$max_walsh")
    calibrated <- cosi_kappa(h_values[-1L],alpha,randomized=randomized,
                             control=control,u=cutoff_uniform)
    law <- calibrated$distribution
    cutoff <- calibrated$kappa
  }
  if (!signed_rank) cutoff <- .cosi_draw_law(law,randomized,"max",cutoff_uniform)

  # 2. Positive gaps separate consecutive batches in the ORIGINAL row order.
  permutation <- if (gap > 0) seq_len(nn) else if (is.null(randomness))
    sample.int(nn) else randomness$permutation
  layout <- .cosi_split_rows(nn,batches,gap,permutation)
  all_estimates <- vapply(seq_len(batches),function(j) {
    value <- tryCatch(estimate(data[layout$batch_rows[[j]],,drop=FALSE]),
      error=function(e) stop(sprintf("estimate failed on batch %d: %s",j,conditionMessage(e))))
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
      stop(sprintf("estimate must return one finite numeric value; failed on batch %d",j))
    }
    as.numeric(value)
  },numeric(1))
  anchor <- if (method == "unimodal-star") all_estimates[1L] else NULL
  W <- if (method == "unimodal-star") all_estimates[-1L] else all_estimates
  uniforms <- if (signed_rank && !sign_scores) .cosi_tie_breakers(B,
    if (is.null(randomness)) NULL else randomness$uniforms) else NULL

  # 3. The constructor only needs estimates, cutoff, and method-specific inputs.
  interval <- switch(method,
    median=cosi_median(W,cutoff),
    symmetry=cosi_symmetry(W,cutoff,h=h_values,U=uniforms,control=control),
    normality=cosi_normality(W,cutoff,U=uniforms,control=control),
    `unimodal-star`=cosi_unimodal_star(W,anchor,cutoff))
  result <- c(interval,list(B=B,K=K,batches=batches,ci_est=W,all_estimates=all_estimates,
    anchor=anchor,uniforms=uniforms,c=cutoff,cutoff_distribution=law,
    kappa=if (signed_rank) cutoff else NULL,kappa_distribution=if (signed_rank) law else NULL,
    h_values=h_values,N_h=if (is.null(h_values)) NULL else sum(h_values),
    split_calibration=split_calibration,batch_sizes=layout$batch_sizes,batch_rows=layout$batch_rows,
    gap=gap,gap_rows=layout$gap_rows,n_used=layout$n_used,
    split_order=if (gap > 0) "ordered" else "shuffled",seed=seed,
    alpha=alpha,gamma=gamma,method=method,randomized=randomized))
  class(result) <- "cosi_ci"
  result
}

# Display the main result, including calibration, gap, and endpoint inclusion.
print.cosi_ci <- function(x, ...) {
  cat(sprintf("COSI %s interval; B=%g, batches=%g, alpha=%g, gap=%g\n",
              x$method,x$B,x$batches,x$alpha,x$gap))
  if (identical(x$kappa_distribution$engine,"concentration")) cat(sprintf(
    "Conservative concentration calibration (%s); rejection bound <= %.6g\n",
    x$kappa_distribution$bound_used,x$kappa_distribution$rejection_probability_bound))
  if (x$empty) cat("empty set\n") else cat(if (x$closed[1L]) "[" else "(",
    format(x$CI[1L]),", ",format(x$CI[2L]),if (x$closed[2L]) "]\n" else ")\n",sep="")
  invisible(x)
}
