HIGH-DIMENSIONAL MEAN: TWO READABLE R FILES

1. cosi_highdim_mean.R
   Core estimation and confidence sets only. No simulations or plots run when
   sourced. It reuses smallest_B() and randomization_quantile_dist() from
   cosi_univariate.R, sourcing that file only if either function is missing.

2. simulation_highdim_mean_example.R
   Sources the core, defines Gaussian/empirical-bootstrap/multiplier benchmarks
   first, then the four generators, simulation loop, summaries, and ggplot code.
   This file runs the full study using its editable settings block at the end.

CORE USAGE

  # Keep cosi_univariate.R in the same working directory, or source it first.
  source("cosi_highdim_mean.R")

  # X has observations in rows and coordinates in columns.
  fit <- cosi_hd(X, alpha=.05, delta=.30, pilot="mom", seed=2026)

  fit$median$contains(theta)       # Membership in cosi-m
  fit$mp$contains(theta)           # Membership in cosi-m+
  fit$mp$lower                    # cosi-m+ coordinate lower bounds
  fit$mp$upper                    # cosi-m+ coordinate upper bounds
  fit$volume_root                 # Named vector: Vol(set)^(1/d) for BOTH sets
  fit$median$log_volume           # Exact log volume of cosi-m
  fit$mp$log_volume               # Exact log volume of cosi-m+
  fit$vartheta                    # Pilot mean estimate
  fit$variance_hat                # Raw pilot covariance diagonal
  fit$sigma_hat                   # Positive scales used by the sets
  fit$Sigma_hat(1:5,1:5)           # Requested covariance block
  fit$Sigma_hat(full=TRUE)         # Full matrix, ONLY when explicitly requested

  empirical <- cosi_hd(X, pilot="empirical", seed=2026)
  both <- cosi_hd(X, pilot="both", seed=2026)
  both$empirical$volume_root
  both$mom$volume_root

The same explicit seed gives the same split and cutoff draws for both pilot
choices. pilot="both" shares these directly. The default seed is2026.
Explicit-seed calls preserve the caller's RNG state. seed=NULL instead uses
and advances the caller's RNG stream, as used by the simulation driver.

A default call uses floor(n/2) rows for the independent pilot and the other
rows for B inference blocks. At n100, alpha.05, delta.30, this means50 pilot
rows and nine inference blocks of sizes5/6. Core inputs can be a numeric
matrix or a numeric data.frame. No data-generating distribution is needed.

PILOT CHOICES

pilot="empirical": vartheta is the pilot sample mean. Sigma_hat is exactly
the usual pilot sample covariance with divisor n_pilot-1, computed in requested
blocks. This includes off-diagonal entries, not just the marginal variances.

pilot="mom": vartheta is the coordinatewise median of K pilot block means.
K is the smallest odd integer with

  P{Bin(K,1/8)>K/2} <= alpha.

The group count depends on alpha alone, with NO dimension factor. K is3 for
EVERY dimension at alpha.05, since P{Bin(3,1/8)>=2}=11/256<.05. The50 pilot
observations form location groups of sizes16/17; the25 paired differences
form covariance groups of sizes8/9.

Pair the pilot rows and define D_r=(X_(2r)-X_(2r-1))/sqrt(2). Each entry of
Sigma_hat is the median across the SAME K groups of group averages of
D_ri*D_rj. Its diagonal is precisely the previous paired-square MOM estimator.
The mean and covariance estimators use the same alpha input and grouping.

This scalar binomial calibration is not a simultaneous95% pilot accuracy claim.
Equation(2.7) of Minsker, arXiv:1812.03523, allows dimension-independent K:
the dimension enters a uniform error bound through the expected empirical
supremum, while the probability bound is1-exp(-c*K). Its universal constants
are unspecified, so it does not give an exact numerical95% calibration for K3.
Conditional COSI coverage allows any independent location pilot and positive
pilot scales; it does not need a separate pilot accuracy event.
The covariance estimator uses the same alpha-only group count; covariance
concentration requires additional product/fourth-moment assumptions. No joint
95% covariance-estimation guarantee is claimed for the t2.5 setting.
The entrywise MOM matrix need not be positive semidefinite. No inversion is
used: only positive marginal scales enter the confidence sets.

Sigma_hat is a FUNCTION representing the entire estimate. Only its diagonal
is needed to construct the confidence sets. A10000-by-10000 matrix needs about
800MB before intermediate calculations, so off-diagonals are computed on demand.
For MOM requests, both axes are chunked to avoid allocating K*d*d values.

variance_hat and Sigma_hat return RAW covariance estimates. sigma_hat is their
diagonal square root with a tiny relative floor for numerical nonpositivity.
scale_regularized records the number of affected coordinates. No coordinates
needed this floor in the completed simulation.

CONFIDENCE-SET DEFINITIONS

For a signed index k>d use coordinate k-d, and use the smallest signed index
in a tie. cosi-m is the selected-coordinate set in the question. It is NOT
the Cartesian product of its coordinate intervals; use median$contains(theta).
The coordinate_lower/coordinate_upper fields are descriptive, not box bounds.

cosi-m+ is the rectangle from the question. Its radius is the maximum OF the
coordinatewise order statistics, not the order statistic of vectorwise maxima.
Its membership can be tested directly using its returned lower/upper vectors.

Both use independent pilot information. One cutoff realization per set is
shared over all coordinates and, when requested together, both pilot choices.
The binomial calibration calls use the original smallest_B() and
randomization_quantile_dist() functions directly. If either is missing, the
core sources cosi_univariate.R in its current evaluation environment.
There are no duplicate calibration definitions, and cosi_univariate.R is
provided unchanged in the bundle.

The exact volume of cosi-m is computed by partitioning space into its signed
coordinate cones and integrating their cross-sections. Computation uses stable
logarithms even at d10000. Vol^(1/d) is the geometric mean coordinate width for
rectangles; the same volume-equivalent size is used for the nonrectangular set.
A bounding-box approximation is never substituted for its volume or membership.

Conditional coverage requires median-unbiased inference block means. The pilot
need not satisfy a separate confidence event for that conditional argument.
Finite covariance alone does NOT imply median unbiasedness. The first three
settings below have symmetric marginals; the fourth deliberately challenges
this assumption.

SIMULATIONS AND PLOTS

Keep the two high-dimensional R files beside cosi_univariate.R, and set that
folder as R's working directory:

  source("simulation_highdim_mean_example.R")

The editable block uses n100, d=(10,30,100,300,1000,3000,10000), alpha.05,
delta.30,500 replications, and500 draws for each bootstrap method.
workers=1 is portable to Windows; it can be increased on Linux/macOS.
Set RUN_SIMULATION=FALSE to load only the simulation and plotting functions.
Only ggplot2>=4.0.0 is required; matrixStats is optional.

The benchmark definitions are unchanged: known-covariance Gaussian reference,
full-data empirical max-t bootstrap, and Gaussian multiplier max-t bootstrap.
The original sample SDs are fixed across resamples. The Gaussian reference is
an oracle diagnostic, not a feasible unknown-covariance procedure.

FOUR SETTINGS: all means0, all marginal variances1, IID rows.
1. Correlated Gaussian: X_ik=(Z_ik+F_i)/sqrt(2), equicorrelation.5.
2. Symmetric heavy tails: independent t2.5/sqrt(5), covarianceI.
3. Rare symmetric shocks: X_ik=(Z_ik+sqrt(n*d)*R_ik*A_ik)/sqrt(2),
   A_ik~Bernoulli(1/(n*d)), R_ik~Rademacher; covarianceI.
4. Rare asymmetric shocks: X_ik=(Z_ik+30*(A_ik-.01))/sqrt(9.91),
   A_ik~Bernoulli(.01); covarianceI and substantial small-block median bias.

In setting3, exactly one shock has probability approaching exp(-1). On that
event, one coordinate of sqrt(n)*barX is approximately +/-sqrt(d/2), whereas
its Gaussian reference maximum is O(sqrt(log d)). This is a genuine failure
of the Gaussian approximation in a finite-covariance triangular array.
Studentization may change bootstrap behavior, so Gaussian-reference failure
and bootstrap coverage are reported separately.

All seven methods share the same dataset in each replication. Both COSI pilot
choices share the same pilot/inference split, block estimates, and cutoffs.
Each setting/dimension has a fixed seed, independent of parallel scheduling.
The simulation fixes Mersenne-Twister/Inversion/Rejection RNG algorithms and
saves the settings and software versions.
The data, empirical-pilot sets, and three benchmark streams are preserved
exactly. Changing the deterministic MOM group count consumes no random numbers.

OUTPUTS

highdim_mean_intervals.csv:98,000 rows =4 settings*7 dimensions*500*7 methods.
empirical_coverage.csv: simultaneous coverage by setting/dimension/method.
median_width_ratios.csv: paired size ratios relative to the oracle Gaussian box.
gaussian_approximation_diagnostics.csv: oracle coverage and raw-max CDF errors.
simulation_records.rds: full-precision replication records.
example_sets.rds: seven sets from the first replication in every cell.
settings.rds, sessionInfo.txt, and checkpoints/: settings and completed cells.
empirical_coverage.pdf/png and width_ratios.pdf/png: two2-by-2 figures.

The legends use cosi-m and cosi-m+, each with '(empirical)' or '(MOM)'. All four
COSI choices are on the second legend row; the three benchmarks are on the first.
There are no coverage error bars, subtitles, or bottom captions. Boxplot
statistics are computed before clipping; outlier points are omitted. The size
axis is logarithmic, with lower limit1/16 and data-based upper limit.

highdim_mean_writeup.tex: copy-paste manuscript write-up with the revised
dimension-independent tuning, numerical results, and figure captions.
highdim_mean_writeup.pdf: readable version with both updated figures.
highdim_mean_references.bib: BibTeX entry for the uniform MOM reference.
The write-up refers to highdim_mean_results/, the simulation's default folder
and the folder name used for results inside the ZIP bundle.

Primary bootstrap reference: Chernozhukov, Chetverikov, Kato and Koike,
High-dimensional Data Bootstrap, arXiv:2205.09691, Sections3.1-3.2.
https://arxiv.org/html/2205.09691v1

MOM reference: Minsker, Uniform bounds for robust mean estimators,
arXiv:1812.03523, Eq.(2.7).
https://arxiv.org/pdf/1812.03523
