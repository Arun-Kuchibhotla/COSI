DEPENDENT LINEAR REGRESSION: FOUR SETTINGS

Run simulation_regression_dependent_example.R from the folder containing
cosi_univariate.R and cosi1d.R. Both supplied core files are unchanged.
The existing cosi1D_dependent() wrapper performs all COSI calculations.
Only function definitions are loaded, so source-time package-install commands
in the original files are not executed.

Required packages: weightedBernoulli0.2.0 and ggplot2>=3.4.0.
The ZIP includes the weightedBernoulli source archive used for this run:

  install.packages("ggplot2")
  install.packages("weightedBernoulli_0.2.0.tar.gz",repos=NULL,type="source")
  source("simulation_regression_dependent_example.R")

Set RUN_SIMULATION=FALSE near the end to load definitions without rerunning.
All functions, settings, summaries and ggplot code are in this single R file.
WORKERS=1 is portable to Windows; increase it on Linux/macOS if desired.

DESIGN

500 replications at n100,200,300,500,750,1000; alpha=.05; delta=.30.
Five covariates and an intercept in every OLS fit. Target: first slope.
The temporal AR coefficients are phi_x=phi_e=.6. Cross-coordinate Gaussian
design correlation is rho^|j-k| with rho=.7. These are different parameters.
Gaussian paths are initialized in their stationary distributions.
Independent design/error latent paths keep the conditional error mean zero.

1. Correlated Gaussian design: Gaussian covariates and Gaussian errors,
   correct linear conditional mean. First-slope target1.
2. Skewed design and errors: centered lognormal covariates; centered
   exponential errors; heteroscedastic scale .5+.5|X1|. Target1.
3. Nonlinear mean, Student t(3) errors: Gaussian covariates; add
   .5X1^2+.25X1^3+.5X1X2 to the linear mean; unit-variance t3 errors.
   Same heteroscedastic scale. First-slope projection target1.75.
4. Nonlinear mean, skewed Pareto errors: same mean/covariates/scales as3;
   centered, variance-standardized Pareto errors of shape2.5. Target1.75.

Non-Gaussian errors are quantile transforms of a separate stationary Gaussian
AR path. All four settings are stationary and geometrically mixing. All have
finite squared loss. Shape2.5 replaces the earlier independent-regression
shape1.5 specifically to preserve the squared-error projection in this task.
The nonlinear settings have population intercept1.85, not1.

GAPPED COSI

B=9 for median, symmetric and normal-score intervals. Unimodal* uses10 blocks,
reserving its first for a pilot and using9 estimates in the interval.
Gap g_n=max(1,floor(log(n)/2)) gives2,2,2,3,3,3 omitted observations between
neighboring blocks. Positive gaps prevent permutation and preserve time order.
The exact retained sample counts and smallest block sizes are in block_counts.csv.
All three nonpilot COSI methods share block estimates within each replication.
Unimodal* has a separate balanced partition because it uses one additional block.
The original cutoff randomizations and independent-pilot construction are used.
No core confidence-interval function was changed.

BENCHMARKS

Wald (HAC): the Bartlett/Newey-West sandwich covariance of Newey and West(1987).
Kernel K(u)=max(1-|u|,0); maximum lag L_n=floor(4*(n/100)^(2/9)).
The actual weight at lag h is1-h/(L_n+1), h=1,...,L_n, and the kernel
bandwidth is L_n+1. The six sample sizes use maximum lags4,4,5,5,6,6.
The saved hac_max_lag field contains L_n, and hac_bandwidth contains L_n+1.
No prewhitening; covariance multiplied by n/(n-6); normal critical value.
This is a fixed sample-size rule, documented by statsmodels cov_hac, not an
automatically estimated bandwidth. The implementation was verified against
sandwich::NeweyWest(fit,lag=L_n,prewhite=FALSE,adjust=TRUE).
References: Newey and West(1987), doi:10.2307/1913610; Zeileis(2004),
doi:10.18637/jss.v011.i10. Detailed equations and additional implementation
references are in the write-up and BibTeX file.
Moving-block bootstrap: resample overlapping pairs blocks of length round(n^(1/3)),
truncate the final block to obtain exactly n rows, and form the basic interval.
500 resamples; type6 quantiles use the (R+1)*p order-statistic convention.
Estimated-rate subsampling: 500 draws from overlapping CONSECUTIVE windows at
b=round(n^(2/3)) and each feasible rate-estimation size round(n^beta), beta=(9:1)/10.
Require at least18 rows per window. Keep Larry Wasserman's100-point inverse-CDF
grid and interquantile log-size regression. No rate clipping is applied.
Prefix cross-products only accelerate OLS: they match direct row-based QR fits.
All three benchmarks use the full sample, including rows omitted by COSI.

REPRODUCIBILITY AND OUTPUTS

Each setting/sample-size cell has seed100+1000000*idx+n. RNGkind is fixed to
Mersenne-Twister/Inversion/Rejection. Serial and parallel results agree exactly.
Every method uses the same generated dataset in each replication. Split seeds
are saved; wrapper seed resets do not affect the next generated dataset.
Rank-deficient fits are kept as unbounded intervals, never silently redrawn.

The output folder dependent_regression_results contains:
- dependent_regression_intervals.csv: replication-level bounds, widths,
  coverage, gaps, retained sample counts, estimated rates, seeds and status.
- simulation.rds: full-precision records, settings and block design.
- empirical_coverage.csv, median_width_ratios.csv: summaries from SAVED CSV.
- unimodal_median_width_ratios.csv: omitted-width medians for figure captions.
- fit_status_counts.csv, infinite_interval_rates.csv, block_counts.csv.
- subsampling_rate_summary.csv: estimated-rate quartiles and medians.
- empirical_coverage.pdf/png and width_ratios.pdf/png: same ggplot objects.
- sessionInfo.txt and checkpoints/: software versions and completed cells.

Coverage plots have four panels and no error bars. Widths are paired ratios
relative to Wald (HAC) from the same data, stored as wald_hac_width.
Width plots omit unimodal*, compute boxplot statistics before clipping, omit
outlier points, and use a logarithmic displayed range[1/16,4]. Both legends have
exactly two rows, with all displayed COSI methods on the second row. There are
no subtitles or bottom captions. The unimodal* medians are saved separately.

INTERPRETING THE THEORY

Gap growth plus within-block CLTs is an asymptotic argument. The finite gaps in
this experiment do not certify exact95% finite-sample coverage. Nonlinear
misspecification and skewness can also cause finite-block median bias.
Conditional on the random design, Gaussian-error block OLS estimates are
jointly Gaussian. Unconditionally they are generally Gaussian mixtures.

At rho=1 in the quoted pure covariance kernel sigma_i*sigma_j*rho^|i-j|,
errors equal sigma_i*Z. Ordinary differencing leaves (sigma_i-sigma_(i-1))*Z;
it does not generally restore mixing. Constant scales remove the errors
entirely, and differencing also removes an intercept column. A nondegenerate
common-factor-plus-idiosyncratic-noise model is different from that pure kernel.
General ARIMA differencing needs a correctly specified transformed regression,
identification of its coefficients, and a separate short-memory verification.
