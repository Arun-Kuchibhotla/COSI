# Optional small numerical checks; no claim about statistical coverage.
source("core_isotonic.R")
source("simulation_isotonic_example.R")
load_isotonic_calibration("cosi_univariate.R")
iso_load_accelerator(".")

# Tied covariates must receive equal fitted values and their original weights.
fit <- iso_fit(c(0,0,1), c(3,0,1), c(0,.5,1))
stopifnot(max(abs(fit$estimate-4/3)) < 1e-12,
  identical(fit$block_n, c(3L,3L,3L)))

# Check an independently specified monotone-envelope example and its tails.
band <- iso_monotone_band(c(.4,.2,.8),c(.6,.9,1),c(.1,.5,.9))
stopifnot(identical(band$lower,c(.4,.4,.8)),
  identical(band$upper,c(.6,.9,1)), abs(band$avg_width-.55)<1e-12)
pred <- iso_band_at(band,c(0,.1,.3,.5,.7,.9,1))
stopifnot(is.infinite(pred$lower[1]),is.infinite(pred$upper[7]),
  identical(pred$upper[-7],c(.6,.6,.9,.9,1,1)))

# The degenerate selected-coordinate tie example excludes its empty cone.
rect <- iso_m_rectangle(c(1,0),c(1,0),c(0,0),c(1,1))
stopifnot(identical(rect$lower,c(1,-1)),identical(rect$upper,c(1,1)))

# Native and base R must agree for repeated pairs, original ties and clamping.
set.seed(123)
x <- sort(sample(seq(0,1,length.out=20),60,replace=TRUE))
y <- x+rnorm(60)
indices <- replicate(40,sort(sample.int(60,60,replace=TRUE)))
grid <- c(-.1,seq(0,1,length.out=15),1.1)
options(isotonic.native=FALSE)
slow <- iso_estimates_many(x,y,indices,grid)
options(isotonic.native=TRUE)
fast <- iso_estimates_many(x,y,indices,grid)
stopifnot(isTRUE(all.equal(slow,fast,tolerance=1e-12)))

# Small complete run, including all twelve DGPs and every public plot exporter.
small <- run_isotonic_study(R=2L,n_grid=c(100L,250L),scenarios=1:12,
  R_boot=30L,N_sub=30L,N_rate=20L,workers=2L,
  output_dir="validation_results",resume=FALSE,make_plots=TRUE)
d <- read.csv(small$csv_file)
stopifnot(nrow(d)==2L*2L*12L*length(ISO_METHODS), !anyNA(d$covered),
  all(d$width_ratio[d$method=="HulC"]==1),
  all(is.finite(d$avg_width)&d$avg_width>=0))
cat("PASS: ties, envelopes, native fallback agreement, all DGPs, CSV and plots.\n")
