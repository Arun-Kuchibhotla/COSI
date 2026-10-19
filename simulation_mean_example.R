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

## The following code simulates the coverage of all intervals including 
## bootstrap and subsampling as delta changes. 
source("subsampling_CI.R")
simulate_coverage_delta <- function(
    R = 500, n = 100, alpha = .05,
    deltas = c(.10, .15, .20, .25, .30, 1/3),
    randomize = TRUE, seed = 100,
    R_boot = 500, N_sub = 500, subsample_exp = 2/3) {
  
  if (!is.null(seed)) set.seed(seed)
  
  B_grid <- vapply(
    deltas,
    function(d) smallest_B(alpha, delta = d),
    numeric(1)
  )
  
  if (any(B_grid > n))
    warning("Dropping delta values requiring more than n blocks.")
  
  keep <- B_grid <= n
  deltas <- deltas[keep]
  B_grid <- B_grid[keep]
  stopifnot(length(deltas) > 0)
  
  b <- round(n^subsample_exp)
  out <- vector("list", 4L * R * length(deltas))
  at <- 0L
  
  for (idx in 1:4) {
    for (r in seq_len(R)) {
      x <- gen_data_mean(n, idx)
      split_seed <- sample.int(.Machine$integer.max, 1)
      
      # Each baseline is calculated ONCE on these same data, then reused
      # at every delta. The rate is estimated, not fixed at 1/2.
      sub <- subsampling_ci(x, estimate = mean,
                         b = b, N = N_sub, alpha = alpha)
      if (any(!is.finite(c(sub$CI, sub$gamma_hat))))
        stop("Nonfinite subsampling result: scenario ", idx,
             ", replication ", r, ".")
      
      baseline <- list(
        Wald = list(
          CI = mean(x) + c(-1, 1) *
            qnorm(1 - alpha/2) * sd(x)/sqrt(n)
        ),
        Bootstrap = basic_bootstrap_mean(x, alpha, R_boot),
        Subsampling = sub
      )
      
      for (g in seq_along(deltas)) {
        fits <- lapply(COSI_METHODS, function(m)
          cosi1D(
            x, mean,
            alpha = alpha,
            delta = deltas[g],
            randomize = randomize,
            seed = split_seed,
            method = m
          )
        )
        
        fits <- c(baseline, setNames(fits, COSI_METHODS))
        records <- do.call(rbind, lapply(fits, .cosi_record))
        
        at <- at + 1L
        out[[at]] <- cbind(
          data.frame(
            scenario = COSI_SCENARIO_LABELS[idx],
            replicate = r,
            n = n,
            delta = deltas[g],
            method = names(fits),
            B = c(rep(NA_integer_, 3), rep(B_grid[g], 4)),
            gamma_hat = c(NA, NA, sub[3], rep(NA_real_, 4))
          ),
          records
        )
      }
    }
    
    message("Completed scenario ", idx, ": ", R, " replications.")
  }
  
  results <- do.call(rbind, out)
  rownames(results) <- NULL
  
  coverage <- aggregate(
    covered ~ scenario + delta + method,
    data = results,
    FUN = mean
  )
  names(coverage)[names(coverage) == "covered"] <- "coverage"
  
  list(
    results = results,
    coverage = coverage,
    design = data.frame(delta = deltas, B = B_grid),
    alpha = alpha,
    subsample_size = b
  )
}

plot_coverage_delta <- function(sim) {
  labels <- c(Wald = "Wald", Bootstrap = "Basic bootstrap",
              Subsampling = "Subsampling (estimated rate)", median = "COSI median",
              sym = "COSI symmetric", normal = "COSI normal scores",
              unimodal = "COSI unimodal*")
  colors <- c(Wald = "#444444", Bootstrap = "#E69F00",
              Subsampling = "#CC79A7", median = "#0072B2", sym = "#009E73",
              normal = "#D55E00", unimodal = "#7B2CBF")
  d <- sim$coverage
  d$method <- factor(d$method, levels = names(labels))
  d$scenario <- factor(d$scenario, levels = COSI_SCENARIO_LABELS)
  ticks <- sort(unique(d$delta))
  tick_labels <- vapply(ticks, function(z)
    if (abs(z - 1/3) < 1e-9) "1/3" else format(z, trim = TRUE), character(1))
  
  ggplot2::ggplot(d, ggplot2::aes(delta, coverage,
                                  color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 1 - sim$alpha,
                        linetype = "dashed", color = "grey55") +
    ggplot2::geom_line(linewidth = .75) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::scale_color_manual(values = colors, labels = labels) +
    ggplot2::scale_x_continuous(breaks = ticks, labels = tick_labels) +
    ggplot2::scale_y_continuous(labels = function(z) paste0(round(100*z), "%")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = expression(delta), y = "Empirical coverage", color = NULL,
                  title = paste0("Coverage as delta changes: n = ", unique(sim$results$n)),
                  caption = "Dashed line: nominal coverage. All methods use the same data in each replication.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
}

sim <- simulate_coverage_delta(n = 100, R = 500, seed = 100)
p <- plot_coverage_delta(sim)
print(p)
