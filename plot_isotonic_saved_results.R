# Regenerate all figures from the supplied CSV files, without simulations.
source("simulation_isotonic_example.R")
directory <- "isotonic_results"
csv_file <- file.path(directory,"isotonic_simulation_records.csv")
export_isotonic_plots(csv_file,output_dir=directory)

# Three-page booklets: one four-scenario figure on each page.
records <- read.csv(csv_file)
groups <- list(1:4,5:8,9:12)
plotters <- list(
  empirical_coverage=function(s) plot_isotonic_coverage(records,s,alpha=unique(records$alpha)),
  band_widths=function(s) plot_isotonic_width(records,s,relative=FALSE),
  width_ratios=function(s) plot_isotonic_width(records,s,relative=TRUE))
for (name in names(plotters)) {
  grDevices::cairo_pdf(file.path(directory,paste0(name,".pdf")),
    width=7.2,height=5.6,onefile=TRUE)
  tryCatch({
    for(s in groups) print(plotters[[name]](s))
  },finally=grDevices::dev.off())
}

# Display every band boundary. Extremely large raw responses may lie outside
# these display limits; all observations remain in the saved CSV.
bands <- read.csv(file.path(directory,"illustration_bands.csv"))
observations <- read.csv(file.path(directory,"illustration_observations.csv"))
observations <- observations[observations$x>=.1 & observations$x<=.9,]
illustration_limits <- function(scenario,method) {
  z <- bands[bands$scenario==scenario & bands$method %in% c(method,"HulC"),]
  limits <- range(c(z$lower,z$upper,z$truth),finite=TRUE)
  limits+c(-1,1)*.05*max(diff(limits),.1)
}
for (method in c("cosi-m+","cosi-m")) {
  name <- if(method=="cosi-m+") "cosi_mplus_band_illustrations" else
    "cosi_m_rectangle_band_illustrations"
  grDevices::cairo_pdf(file.path(directory,paste0(name,".pdf")),
    width=7.2,height=5.6,onefile=TRUE)
  tryCatch({
    for(s in 1:12) print(plot_isotonic_illustration(bands,s,observations,
      method=method,ylim=illustration_limits(s,method)))
  },finally=grDevices::dev.off())
}
for(s in c(1L,4L,7L,10L)) {
  p <- plot_isotonic_illustration(bands,s,observations,method="cosi-m+",
    ylim=illustration_limits(s,"cosi-m+"))
  ggplot2::ggsave(file.path(directory,sprintf("band_illustration_scenario_%02d.png",s)),
    p,width=7.2,height=5.6,units="in",dpi=200,bg="white")
}

# Empirical sign balance of the inference-block isotonic estimates.
# This is a finite-sample diagnostic, with no asymptotic conclusion imposed.
diagnostics <- read.csv(file.path(directory,"grid_diagnostics.csv"))
diagnostic_plot <- function(scenarios) {
  z <- diagnostics[diagnostics$scenario %in% scenarios,]
  z$scenario_label <- factor(z$scenario_label,levels=ISO_SCENARIO_LABELS[scenarios])
  z$n_label <- factor(z$n,levels=sort(unique(z$n)))
  ggplot2::ggplot(z,ggplot2::aes(t,split_above,color=n_label,group=n_label))+
    ggplot2::geom_hline(yintercept=.5,linetype="dashed",color="#7F8993")+
    ggplot2::geom_line(linewidth=.5)+
    ggplot2::geom_point(size=.7)+
    ggplot2::facet_wrap(~scenario_label,ncol=2)+
    ggplot2::scale_color_manual(values=c("#A8B1BB","#56B4E9","#0072B2","#009E73"))+
    ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+
    ggplot2::labs(x="Anchor location",y="Fraction of split estimates above truth",
      color="Sample size",title=NULL,subtitle=NULL,caption=NULL)+
    isotonic_plot_theme()+
    ggplot2::guides(color=ggplot2::guide_legend(nrow=1))
}
grDevices::cairo_pdf(file.path(directory,"median_balance.pdf"),
  width=7.2,height=5.6,onefile=TRUE)
tryCatch({
  for(g in seq_along(groups)) {
    p <- diagnostic_plot(groups[[g]])
    print(p)
    ggplot2::ggsave(file.path(directory,paste0("median_balance_group_",g,".pdf")),p,
      device=grDevices::cairo_pdf,width=7.2,height=5.6,units="in",bg="white")
  }
},finally=grDevices::dev.off())
