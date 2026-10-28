MONOTONE REGRESSION BANDS: REPRODUCIBLE SIMULATION STUDY

Run from the directory containing this README and the R files.
The only additional R package needed is ggplot2 (version 3.4 or newer).

  install.packages("ggplot2")
  source("simulation_isotonic_example.R")
  study <- run_isotonic_study(core_dir=".", output_dir="isotonic_results")

The default study uses 500 Monte Carlo replications for each of twelve
scenarios at n=100,250,500,1000; alpha=.05 and delta=.30. All six methods
use the same data in a replication. COSI uses half for an isotonic pilot
and the other half for B=9 independently fitted blocks. The pilot scales
are exactly 1/sqrt(PAVA block count), without a noise-scale multiplier.

Methods:
  cosi-m+                  rectangular grid construction
  cosi-m (rectangle)       coordinatewise bounding rectangle of cosi-m
  HulC (Bonferroni)        deterministic B_H=ceil(log2(2K/alpha)) splits
  Pairs bootstrap         symmetric standardized supremum-root band
  Subsampling (1/3)        fixed cube-root rate
  Subsampling (estimated)  a separate estimated rate at each grid point

The bootstrap uses 500 resamples. Both subsampling methods share 500
construction subsets of size round(n^(2/3)); rate estimation uses 200
subsets at each feasible size round(n^beta), beta=.9,.75,.6,.45, min size8.
Quantile ranges use the 100-point inverse-CDF grid from the HulC code.
Finite estimated rates are not truncated; degenerate coordinates receive
rate0 with a recorded fallback flag.

WHAT WIDTH MEANS
The grid has K=floor(sqrt(n)) anchors in [.1,.9]. Its nearest anchor to
.5 is replaced by .5. Monotonicity extends the grid rectangle to a band.
Outside the anchor span the band has an infinite one-sided tail. Every
reported finite width is the average width integrated over [.1,.9].
width_ratio is paired with the SAME replication's HulC average width.
The log-scale ratio plots show [1/16,8]; boxplot statistics are calculated
before display clipping, and outlier points are omitted.

REPRODUCIBILITY AND SPEED
Each replication uses seed=100+1000000*scenario+1000*n+replicate.
COSI and benchmark streams have offsets50000000 and100000000.
RNG algorithms are explicitly Mersenne-Twister/Inversion/Rejection.
Serial and forked runs agree under the same software versions.
Use workers=1 on Windows or when running serially:

  study <- run_isotonic_study(workers=1)

The C source isotonic_fast.c is an optional batch PAVA accelerator. R
attempts to compile it automatically with R CMD SHLIB. A compiler or
Rtools is useful but not required: the base-R fallback remains available.
To force that fallback, run options(isotonic.native=FALSE) before the study.
The native version draws no random numbers and was checked against the
base-R predictions and resampling outputs. Platform-specific binaries
are intentionally omitted from this bundle.

Checkpoint RDS files are written after each completed setting/sample-size
cell. Existing checkpoints are reused only if the full settings and source
hashes match. Use a fresh output_dir or resume=FALSE for changed settings.
The bundle includes the final simulation.rds, rather than duplicate cell
checkpoints. Rerunning the full study regenerates those checkpoints.

PLOTS FROM SAVED CSV FILES
No simulation is needed to redraw the figures:

  source("simulation_isotonic_example.R")
  export_isotonic_plots(
    "isotonic_results/isotonic_simulation_records.csv",
    illustrations_file="isotonic_results/illustration_bands.csv",
    observations_file="isotonic_results/illustration_observations.csv")

To recreate the combined coverage/width PDF booklets, run:

  source("plot_isotonic_saved_results.R")

Main figures use groups1:4,5:8,9:12, each in a2-by-2 layout. Legends use
two rows, with both COSI methods on the second row. The PDFs and PNGs
are exported from the same ggplot objects. Coverage panels use separate
lower limits, rounded down with at least two percentage points of margin,
and a common upper limit of 101%; all empirical values remain visible.
Use plot_isotonic_coverage(..., ylim=c(0,1.01)) for a common full scale.
The illustration booklets have
one page per scenario, each showing the four sample sizes. They use the
first simulation replication, selected in advance without searching.

FILES
  core_isotonic.R                  fitting, pilot, COSI grid sets and bands
  simulation_isotonic_example.R    DGPs, comparisons, simulation and ggplot
  cosi_univariate.R                supplied calibration core, unchanged
  isotonic_fast.c                  optional native accelerator source
  plot_isotonic_saved_results.R    figures and combined PDF booklets
  validate_isotonic_example.R      small numerical/reproducibility checks
  isotonic_regression_report.pdf   methodology, results and illustrations
  isotonic_regression_report.tex   editable report source
  isotonic_references.bib          bibliography with specific locations
  isotonic_results/                raw CSV, summaries, RDS, PDFs and PNGs

IMPORTANT INTERPRETATION
The cosi-m output is the bounding rectangle used to make a band, not the
original nonrectangular selected-coordinate set or its volume. The two
COSI constructions use different randomized cutoffs and need not be nested.
The jump scenario fixes f0(.5)=.8; its behavior is sensitive to the explicit
last-design-point step convention. The t1.5 scenario has a conditional mean
but infinite error variance. Coverage in this report is empirical; the
report contains no theoretical coverage guarantees.

SOURCE PROVENANCE
The supplied cosi_univariate.R is copied unchanged from COSI commit
3d4bcaa3f2099ae0348c481e1ba203315a78af8c. Its package-installation commands
are not executed: only smallest_B and randomization_quantile_dist are read.
The report gives pinpoint references for the block-count scale and
monotone grid extension. Settings, package versions, core hashes, and
sessionInfo are saved beside the results.

PRESENTATION REVISION
The coverage plotting function and report presentation were revised after
the simulation run. Simulation algorithms, core functions and saved CSV
records are unchanged. source_md5.csv records the sources used for the
simulation; SHA256SUMS.txt records the current distribution files.
The report places one-replication illustrations immediately after the
settings, removes forced section page breaks, and discusses the bootstrap
comparison and the sources of the HulC/COSI width differences.
