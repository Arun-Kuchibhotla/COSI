# Run from the COSI repository root after copying this docs directory there.
# All figures use existing saved results; no simulations or core edits.
# Requires ggplot2. Mean-PDF conversion additionally requires pdftoppm (Poppler).
if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Install ggplot2 first.")
dir.create("docs/assets", recursive=TRUE, showWarnings=FALSE)

file.copy("highdim_mean_results/empirical_coverage.png",
  "docs/assets/highdim_coverage.png", overwrite=TRUE)
file.copy("dependent_regression_results/empirical_coverage.png",
  "docs/assets/dependent_coverage.png", overwrite=TRUE)
file.copy("isotonic_results/band_illustration_scenario_04.png",
  "docs/assets/isotonic_band.png", overwrite=TRUE)

if (nzchar(Sys.which("pdftoppm"))) {
  status <- system2("pdftoppm", c("-png", "-scale-to", "1800", "-singlefile",
    "Illustration-of-cosi-intervals-mean-example.pdf", "docs/assets/mean_intervals"))
  if (status != 0L) stop("Mean illustration conversion failed.")
} else if (!file.exists("docs/assets/mean_intervals.png")) {
  stop("Install Poppler to convert the mean PDF, or use the included PNG.")
}

# Load only plotting functions and their named constants. The regression
# script normally starts a full study when sourced; that is not needed here.
constants <- c("COSI_METHODS", "REGRESSION_SCENARIO_LABELS",
  "REGRESSION_METHOD_LABELS", "REGRESSION_COLORS")
for (expr in parse("simulation_regression_example.R")) {
  assignment <- is.call(expr) && is.symbol(expr[[1]]) &&
    as.character(expr[[1]]) %in% c("<-", "=")
  if (!assignment) next
  function_definition <- is.call(expr[[3]]) &&
    identical(expr[[3]][[1]], as.name("function"))
  if (function_definition || as.character(expr[[2]]) %in% constants)
    eval(expr, envir=.GlobalEnv)
}

coverage <- plot_coverage_n_regression("regression_sample_size_intervals.csv") +
  ggplot2::labs(title=NULL, subtitle=NULL, caption=NULL) +
  ggplot2::facet_wrap(~scenario, ncol=2,
    labeller=ggplot2::label_wrap_gen(width=34))
# Use precomputed boxplot summaries, so infinite widths are not silently
# discarded. This drawing code supports ggplot2 3.4 as well as newer versions.
records <- read_regression_results("regression_sample_size_intervals.csv")
summary <- regression_boxplot_data(records[records$method != "unimodal", ])
summary$n_label <- factor(summary$n, levels=sort(unique(summary$n)))
order <- c("Wald","Bootstrap","Subsampling","spacer1",
  "median","sym","normal","spacer2")
colors <- c(REGRESSION_COLORS,spacer1=NA_character_,spacer2=NA_character_)
labels <- c(REGRESSION_METHOD_LABELS,spacer1="",spacer2="")
coordinate <- if ("coord_transform" %in% getNamespaceExports("ggplot2"))
  ggplot2::coord_transform else ggplot2::coord_trans
widths <- ggplot2::ggplot(summary,
  ggplot2::aes(n_label,fill=method,ymin=ymin,lower=lower,middle=middle,
    upper=upper,ymax=ymax))+
  ggplot2::geom_hline(yintercept=1,linetype="dashed",color="#7F8993")+
  ggplot2::geom_boxplot(stat="identity",width=.8,outlier.shape=NA,
    position=ggplot2::position_dodge(width=.86),linewidth=.3)+
  ggplot2::facet_wrap(~scenario,ncol=2,
    labeller=ggplot2::label_wrap_gen(width=34))+
  ggplot2::scale_fill_manual(values=colors,limits=order,breaks=order,
    labels=labels[order],drop=FALSE)+
  ggplot2::scale_y_continuous(breaks=c(1/16,1/8,1/4,1/2,1,2,4),
    labels=c("1/16","1/8","1/4","1/2","1","2","4"))+
  coordinate(y="log10",ylim=c(1/16,4))+
  ggplot2::labs(x="Sample size",y="Width / Wald width",fill=NULL,
    title=NULL,subtitle=NULL,caption=NULL)+regression_plot_theme()+
  ggplot2::guides(fill=ggplot2::guide_legend(nrow=2,byrow=TRUE,
    override.aes=list(alpha=c(1,1,1,0,1,1,1,0),
      colour=c(rep("#35424E",3),NA,rep("#35424E",3),NA))))
ggplot2::ggsave("docs/assets/regression_coverage.png", coverage,
  width=9.2, height=6.6, dpi=200, bg="white")
ggplot2::ggsave("docs/assets/regression_width_ratios.png", widths,
  width=9.2, height=6.6, dpi=200, bg="white")
