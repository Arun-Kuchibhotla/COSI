# COSI confidence intervals in R and Python

COSI computes a scalar estimate on each of several disjoint data batches and
constructs a confidence interval from those estimates. This implementation
provides median, symmetry, normality, and unimodal-star methods, standalone
interval constructors, and cutoff-calibration helpers.

R uses base R. Python requires Python 3.10 or later, NumPy, and SciPy. Keep the
four source modules for each language together.

```sh
python -m pip install -r requirements.txt
```

## 1. Files and functions

| Module (available as `.R` and `.py`) | Purpose |
|---|---|
| `cosi_smallest_B` | Select B and compute binomial cutoff distributions. |
| `cosi_calibration` | Compute signed-rank cutoff distributions and draw a cutoff using `cosi_kappa`. |
| `cosi_intervals` | Construct intervals from supplied batch estimates and cutoffs. |
| `cosi_ci` | Manage seeds, select B, split data, compute estimates, and call the relevant constructor. |

| Public function | Goal |
|---|---|
| `cosi_ci` | Compute a complete interval directly from data and an estimator. |
| `cosi_median` | Construct the median interval from B estimates and a binomial cutoff. |
| `cosi_symmetry` | Construct the weighted Walsh interval from B estimates, scores, and a signed-rank cutoff. |
| `cosi_normality` | Construct the Walsh interval with fixed normal scores. |
| `cosi_unimodal_star` | Construct the expanded-anchor interval from B estimates and an independent anchor. |
| `cosi_kappa` | Calibrate and draw the signed-rank cutoff from a score vector, alpha, and a seed. |
| `smallest_B_distribution` | Return the smallest admissible B and its binomial cutoff law without random draws. |
| `binomial_cutoff_distribution` | Return the binomial cutoff law at a specified tail level. |
| `largest_support` | Return the largest support value of the binomial cutoff at alpha/2. |
| `smallest_B` | Return the selected B, its cutoff law, and a randomized binomial cutoff. |
| `make_cosi_randomness` | Generate shared random inputs for reproducibility and R/Python comparisons. |
| `contains_cosi_ci` | Check membership while respecting endpoint inclusion and empty intervals. |

Each named function has a purpose comment before its definition. Internal
helpers start with `.cosi_` in R or `_` in Python.

Load all R functions with `source("cosi_ci.R")`. This sources the three
neighboring modules automatically. To load calibration alone, source
`cosi_smallest_B.R` followed by `cosi_calibration.R`. Python functions can be
imported from their respective modules.

## 2. Main function

```r
cosi_ci(data, alpha=0.05, gamma=1/3, estimate,
        method=c("median", "symmetry", "normality", "unimodal-star"),
        gap=0, seed=NULL,
        h=function(i) i, randomized=TRUE, control=list(), randomness=NULL)
```

```python
cosi_ci(data, alpha=0.05, gamma=1/3, estimate=None,
        method="median", gap=0, seed=None, *,
        h=None, randomized=True, control=None, randomness=None)
```

`estimate` is required; pass it by name. `data` is converted to a matrix,
with a vector becoming a one-column matrix. Every call to `estimate` receives
a two-dimensional row subset, including when that subset has one row.
The estimator must return one finite real scalar. It can handle missing data
internally, but a nonfinite returned estimate causes an error.

`alpha` must lie in (0,1), and `gamma` in (0,1/2). `gap` must be a nonnegative
integer. The method aliases `"symmetric"` and `"unimodal"` are accepted.
The advanced arguments `h`, `randomized`, `control`, and `randomness` are
explained below. Only the symmetry method uses a supplied `h`.

```r
source("cosi_ci.R")
data <- cbind(x=sin(1:240), y=cos(1:240))
fit <- cosi_ci(data, alpha=.05, gamma=1/3,
               estimate=function(d) mean(d[,1]),
               method="normality", gap=2, seed=42)
print(fit)
fit$CI
fit$B
fit$K
fit$c
fit$batch_rows
fit$gap_rows
```

```python
import numpy as np
from cosi_ci import cosi_ci

data = np.column_stack((np.sin(np.arange(1,241)), np.cos(np.arange(1,241))))
fit = cosi_ci(data, alpha=.05, gamma=1/3,
              estimate=lambda d: np.mean(d[:,0]),
              method="normality", gap=2, seed=42)
print(fit)
fit["CI"]
fit["B"]
fit["K"]
fit["c"]
fit["batch_rows"]
fit["gap_rows"]
```

## 3. B, K, and splitting with gaps

B is the number of estimates used in the final order-statistic or Walsh
construction. **K is the actual number of data batches:**

$$
K=\begin{cases}
B+1,&\text{method = unimodal-star},\\
B,&\text{method = median, symmetry, or normality}.
\end{cases}
$$

For unimodal-star, the first batch supplies the pilot/anchor estimate a;
the remaining B batches supply W. The return fields `K` and `batches` are
aliases for the same batch count. `ci_est` always has length B;
`all_estimates` has length K.

B is the smallest positive integer satisfying

$$
\left|\frac{B-q_B}{B}-\frac{1}{2}\right|\le\gamma,
\qquad q_B<B,
\qquad q_B=\max\operatorname{supp}(c^*_{B,\alpha/2}).
$$

Here the independently randomized binomial cutoff obeys
$\mathbb{P}\{\operatorname{Bin}(B,1/2)\le c^*_{B,\alpha/2}\}=1-\alpha/2$.
The search uses finite-sample quantile bounds to remove impossible ranges,
then checks the remaining integers in increasing order. It does not assume
that feasibility is monotone in B. B selection uses neither data values nor
randomness. At alpha=.05 and gamma=1/3, B=6; at gamma=.1, B=90.

With n rows and gap g:

- **g=0:** shuffle rows once and split the resulting order into K balanced
  batches. All rows are used.
- **g>0:** do not shuffle. Use K contiguous batches in the original row order,
  leaving exactly g unused rows between adjacent batches. There are no
  leading or trailing gaps.

The number of used rows is n-(K-1)g. The sample-size requirement is

$$
n\ge K+(K-1)g,
$$

so every batch is nonempty. All available nongap rows are used, and batch
sizes differ by at most one. The gap does not change the definition of B.
For dependent observations, positive gaps preserve their supplied order;
separation alone does not establish the independence assumptions needed for
exact coverage.

Batch labels are `"0"`,...,`"K-1"`. Let r be the remainder of n-(K-1)g modulo K.
The r additional rows are assigned to labels 1,...,r. At gap zero this is the
R grouping `split(1:n, sort((1:n) %% K))` applied to the shuffled row positions.

For example, n=17, K=6, and g=1 produce batches `1:2`, `4:5`, `7:8`, `10:11`,
`13:14`, `16:17`; omitted indices are `3,6,9,12,15`.

`batch_rows` records original row indices for every batch; `gap_rows` records
omitted indices. **Both use one-based indexing in R and Python.** Subtract
one before using these indices to index a NumPy array.

## 4. Standalone interval constructors

These functions take already computed estimates. They neither choose B nor
split data nor calibrate the cutoff. B is inferred from the length of W.

| Function | Meaning of c | Construction |
|---|---|---|
| `cosi_median(W, c)` | A draw from the binomial cutoff at alpha/2 | Closed batch-estimate order-statistic interval. |
| `cosi_symmetry(W, c, h, U)` | Signed-rank cutoff kappa for the supplied h | Weighted augmented Walsh interval; unit scores use the closed median convention. |
| `cosi_normality(W, c, U)` | Signed-rank cutoff kappa for normal scores | Weighted augmented Walsh interval with fixed normal scores. |
| `cosi_unimodal_star(W, a, c)` | A draw from the binomial cutoff at alpha | Expanded-anchor interval. |

Each returns `CI` (lower and upper endpoints), `closed` (endpoint inclusion),
and `empty`. Unimodal-star also returns the vectors `L` and `U`.
General signed-rank constructors also return `augmented_endpoints` and
`mass_tolerance`; the unit-score median reduction returns the median result
without augmented endpoints.

### Median and unimodal-star formulas

Writing $W_{(1)}\le\cdots\le W_{(B)}$, the median interval is

$$
[W_{(B-c)},W_{(c+1)}].
$$

For unimodal-star define

$$
L_j=\min(W_j,a)-(1+\sqrt{2})|W_j-a|,\qquad
U_j=\max(W_j,a)+(1+\sqrt{2})|W_j-a|.
$$

Its interval is $[L_{(B-c)},U_{(c+1)}]$. B selection uses alpha/2,
but this method's actual cutoff is calibrated at alpha. Extended order
statistics use -infinity below index 1 and +infinity above index B. Some
large-alpha randomized outcomes can produce an empty interval; it is
represented by `empty=TRUE`/`True`.

### Symmetry and normality formulas

For symmetry, h defaults to h(i)=i. It can be supplied as a scalar function
or as the B+1 values h(0),...,h(B). These values must be finite, nonnegative,
and nondecreasing, with **h(0)=0**. The zero-origin condition ensures

$$
\sum_{i\le j}\{h(j-i+1)-h(j-i)\}=\sum_{i=1}^B h(i)=N_h.
$$

Normality fixes

$$
h(i)=\Phi^{-1}\left(\frac{1}{2}+\frac{i}{2(B+1)}\right),\quad i=0,\ldots,B.
$$

For general scores, U consists of B independent nonzero U(-1,1) tie-breakers.
A supplied U is reused; otherwise the constructor draws U. Sort $(W_i,U_i)$
lexicographically, construct augmented Walsh averages

$$
A_{ij}=\left((W_{[i]}+W_{[j]})/2,(U_{[i]}+U_{[j]})/2\right),
\quad \omega_{ij}=h(j-i+1)-h(j-i),\quad i\le j,
$$

and sort them lexicographically. With $Q_h$ the weighted quantile, the interval
represents

$$
\{t: Q_h(c+)\prec(t,0)\prec Q_h(N_h-c)\}.
$$

The lower quantile uses the right-hand limit, computed by a strict cumulative
weight comparison. It is not computed by adding one to a real-valued cutoff.
If the lower augmented endpoint is (l,u), l is included precisely when u<0.
If the upper endpoint is (r,v), r is included precisely when v>0. Equal real
endpoints can represent either a singleton or an empty set.

Use `contains_cosi_ci(result,t)` for membership; reading only `CI` loses this
endpoint information. Infinite endpoints are always open.

## 5. Unit scores: exact agreement with the median interval

Set h(0)=0 and h(i)=1 for every i>0. Then only diagonal Walsh averages have
positive weight, $N_h=B$, and the score-sum law is Binomial(B,1/2). The code
uses the relation

$$
\kappa=B-c^*_{B,\alpha/2}-1.
$$

For an integer kappa, the numerical Walsh endpoints are
$W_{(\kappa+1)}$ and $W_{(B-\kappa)}$. The unit-score branch returns these
as the **closed median interval**, including when estimates are tied. U is
unnecessary, and no Walsh-pair allocation or tie-breaker draw is performed.
The supplied cutoff argument remains kappa in `cosi_symmetry`; it is not the
binomial cutoff itself.

Cutoff randomization is coupled through the same binomial draw. Thus median
and unit-score symmetry calls have identical `CI`, `closed`, `empty`, and
batch estimates when given the same data, deterministic score function,
estimator, options, and seed within one language, or the same shared random
plan with a deterministic estimator. Equality also holds with
`randomized=FALSE`. Independent calls without matching random inputs need
not realize the same interval.

```r
source("cosi_ci.R")
x <- rep(0:4, 30)
est <- function(d) mean(d[,1])
a <- cosi_ci(x, estimate=est, method="median", seed=42)
b <- cosi_ci(x, estimate=est, method="symmetry",
             h=function(i) as.numeric(i>0), seed=42)
stopifnot(identical(a$CI,b$CI), identical(a$closed,b$closed),
          identical(a$empty,b$empty), b$c == b$B-a$c-1)

# Direct construction from an already chosen binomial cutoff:
W <- c(5,1,4,2,6,3)
c_binomial <- 4
stopifnot(identical(
  cosi_median(W,c_binomial),
  cosi_symmetry(W,length(W)-c_binomial-1,h=c(0,rep(1,length(W))))))
```

### Central-cutoff distinction

The strict augmented Walsh set and the closed median interval are different
endpoint conventions. Also, when B is even and the upper kappa support
point equals B/2, its two tails overlap. The binomial randomization needed
for median agreement can then be conservative relative to the canonical
exact two-sided kappa law.

For example, B=2 and alpha=.9 give binomial c values 0 and 1 with probabilities
.4 and .6. Median-compatible kappa therefore takes 0 and 1 with probabilities
.6 and .4. Its two-sided statistic rejection probability is .7. The canonical
exact-level kappa mixture has probabilities .2 and .8 and rejection
probability .9. These two randomizations cannot simultaneously be identical.

The default unit-score path retains the median-compatible law and closed
endpoints as specified above. Its law records `engine="sign-binomial"`,
`median_compatible=TRUE`, `nominal_level_exact`, and
`rejection_probability_bound`. The reported bound also covers the closed
interval, which can have additional conservatism at ties.

## 6. Public kappa-calibration helper

```r
cosi_kappa(h, alpha=0.05, seed=NULL, randomized=TRUE,
           control=list(), u=NULL, median_compatible=TRUE)
```

```python
cosi_kappa(h, alpha=0.05, seed=None, *, randomized=True,
           control=None, u=None, median_compatible=True)
```

**Here h is the vector h(1),...,h(B), of length B; omit h(0).** This differs
from the B+1 score vector accepted by `cosi_symmetry`, whose leading zero is
explicit. Scores must be finite, nonnegative, and nondecreasing. The helper
returns `kappa` and `distribution`; the latter includes support values,
probabilities, and the calibration engine. It does not access data or draw U.
Both symmetry and normality in `cosi_ci` call this helper.

For general scores, it calibrates

$$
G=\sum_{i=1}^B h(i)\xi_i,\qquad \xi_i\overset{\mathrm{iid}}{\sim}\operatorname{Ber}(1/2),
$$

using the two-sided probability
$P_h(c)=\mathbb{P}\{|G-N_h/2|\ge N_h/2-c\}$, with $P_h(-1)=0$.
When exact calibration is feasible, adjacent critical values bracket alpha
and are mixed with probability
$(\alpha-P_h(c_-))/(P_h(c_+)-P_h(c_-))$ on $c_+$.

Unit scores default to the median-compatible law in Section 5.
Set `median_compatible=FALSE`/`False` to compute the canonical $P_h$ mixture
instead, including at a central cutoff. This switch concerns calibration;
`cosi_symmetry` itself uses closed median endpoints for unit scores.

```r
source("cosi_ci.R")
W <- c(-2,-1,0,1,3,4)
B <- length(W)
scores <- qnorm(.5+(1:B)/(2*(B+1)))
cal <- cosi_kappa(scores, alpha=.05, seed=42)
cal$kappa
cal$distribution$support
cal$distribution$probabilities
U <- c(-.8,.2,-.4,.7,-.5,.9)
ci <- cosi_normality(W,c=cal$kappa,U=U)
```

```python
import numpy as np
from scipy.special import ndtri
from cosi_calibration import cosi_kappa
from cosi_intervals import cosi_normality
W = np.array([-2,-1,0,1,3,4])
B = len(W)
scores = ndtri(.5+np.arange(1,B+1)/(2*(B+1)))
cal = cosi_kappa(scores, alpha=.05, seed=42)
U = [-.8,.2,-.4,.7,-.5,.9]
ci = cosi_normality(W,c=cal["kappa"],U=U)
```

For cutoff selection, `randomized=FALSE` chooses the smallest kappa support
value. Median and unimodal-star instead choose the largest binomial c.
Concentration calibration has a singleton law, so both settings coincide.
The optional `u` supplies one cutoff uniform in [0,1) for coupling/testing.

## 7. Exact calibration and concentration bounds

The calibration engine tries closed-form endpoint cases, Wilcoxon formulas
for linear scores, and lattice dynamic programming when feasible. Small
general score problems use enumeration; larger general problems below the
configured threshold use meet-in-the-middle subset-sum counting.

For remaining general problems with **B >= 40**, concentration calibration
computes a conservative real cutoff. Exact Wilcoxon and manageable lattice
cases remain available above 40. State/counting/iteration limits can also
trigger concentration calibration below 40.

For $N_h>0$, set $w_i=h(i)/N_h$, $s^2=\sum_i w_i^2$, and
$L=\log(2/\alpha)$. For the normalized Rademacher sum
$T=\sum_i w_i\varepsilon_i$, the implementation takes the smallest computed
radius from four bounds on $\mathbb{P}(|T|\ge r)$:

- Hoeffding: $r=s\sqrt{2L}$.
- Chernoff: numerically minimize
  $[L+\sum_i\log\cosh(\lambda w_i)]/\lambda$ over finite positive lambda.
- Gaussian comparison:
  $r=s\Phi^{-1}(1-\alpha/(2C_*))$, with
  $C_*=[4\{1-\Phi(\sqrt{2})\}]^{-1}$.
- Chebyshev: $r=s/\sqrt{\alpha}$.

The Gaussian comparison is a finite-sample Rademacher inequality, not a
normal approximation. These are four implemented bounds; no universal
optimality over all possible concentration methods is asserted. The
Gaussian comparison is from [Bentkus and Dzindzalieta, Theorem 1.1](https://arxiv.org/abs/1307.3451).

The cutoff is $\kappa=(N_h/2)(1-r)$, decreased by a numerical cushion.
A negative result is replaced by -1, yielding the whole real line. This
calibration is deterministic and conservative. All computations use floating
point, including exact-distribution engines; they are not arbitrary-precision
certificates. A finite Chernoff parameter remains valid even if numerical
optimization reaches its iteration limit.

For a concentration result, `distribution` records `exact=FALSE`,
`bound_used`, `candidate_radii`, `radius`, `radius_fraction`, `lambda`,
`numerical_cushion`, `fallback_reason`, and a checked
`rejection_probability_bound` with its logarithm. Actual rejection
probabilities are not computed and are returned as NA/NaN. On exact-distribution
paths, `exact=TRUE` describes how the law was computed; it does not imply
exact interval coverage under ties or median-compatible central randomization.

The `control` list/dictionary accepts:

| Option | Default | Purpose |
|---|---:|---|
| `max_states` | $2^{20}$ | Limit lattice states and half-subset sums. |
| `max_dp_work` | 50,000,000 | Limit lattice work. |
| `max_walsh` | 2,000,000 | Limit Walsh pairs for general-score interval construction. |
| `max_iterations` | 256 | Limit pair-sum selection and Chernoff optimization. |
| `concentration_B` | 40 | Switch remaining general calibration problems to concentration. |
| `mass_tol` | 64 times machine epsilon | Tolerance for real-score cumulative Walsh weights. |

The unit-score median branch does not construct Walsh pairs and is exempt
from `max_walsh`. For exactly representable integer weights, cumulative-weight
tolerance is zero. Score distinctions near floating-point precision may
require higher-precision calculations.

## 8. Seeds and cross-language comparisons

In `cosi_ci`, seed is NULL/None or an integer in [0,2147483647]. A supplied
seed controls splitting, cutoff selection, tie-breakers, and ordinary RNG
calls inside `estimate` and h. The caller's RNG state is restored even if
validation or estimation fails. The same local-seed behavior applies to
`cosi_kappa`.

R uses a temporary `set.seed`. Python uses a local NumPy Generator for COSI
and temporarily seeds `numpy.random` and Python's `random` for callbacks.
A callback's independently managed generators, external libraries, and
parallel workers require their own seeds. With no seed, R uses its current
RNG; Python creates a fresh local generator and leaves callback RNGs alone.
The Python `cosi_kappa` helper also accepts a NumPy Generator as seed to
consume an existing stream; `cosi_ci` uses this internally.

Equal numeric seeds do not yield equal draws between R and Python. For
matching realizations, supply identical data, a deterministic estimator, and
the same `randomness` plan. Generate it with
`make_cosi_randomness(n,B,seed,gap)`. The plan has a one-based `permutation`,
a `cutoff_uniform`, and B nonzero tie-breakers in `uniforms`. With positive
gaps its permutation must be the identity; a nonidentity permutation is
rejected. Plans should be generated independently of the data.

A plan may be combined with a seed: the plan supplies COSI's random inputs,
and the seed controls callback randomness. Matching randomized callbacks
across languages requires supplying their random inputs as well.

## 9. Returned values

The main result is an S3 list in R and a `CosiResult` dictionary in Python.
Named R vectors/lists become Python dictionaries, unnamed numeric vectors
become NumPy arrays, and NULL becomes None.

| Field | Meaning |
|---|---|
| `CI`, `closed`, `empty` | Real endpoints, endpoint inclusion, and emptiness. |
| `B` | Number of final estimates W. |
| `K`, `batches` | Actual number of data batches: B or B+1 according to the method. |
| `ci_est`, `all_estimates` | B final estimates and all K estimates, respectively. |
| `anchor` | First batch estimate for unimodal-star; otherwise NULL/None. |
| `c`, `cutoff_distribution` | Actual constructor cutoff and its law for every method. |
| `kappa`, `kappa_distribution` | Aliases for signed-rank c and its law; otherwise NULL/None. |
| `split_calibration` | B-selection law at alpha/2. |
| `h_values`, `N_h` | h(0),...,h(B) and the total score for signed-rank methods. |
| `uniforms` | Used tie-breakers; NULL/None when unnecessary, including unit scores. |
| `batch_rows`, `batch_sizes` | Original one-based indices and sizes for all K batches. |
| `gap`, `gap_rows`, `n_used` | Gap size, omitted indices, and number of rows used. |
| `split_order`, `seed` | `"ordered"` or `"shuffled"`, and the supplied seed. |
| `L`, `U` | Expanded per-estimate endpoints for unimodal-star. |
| `augmented_endpoints`, `mass_tolerance` | General-score Walsh endpoints and numerical tolerance. |
| `alpha`, `gamma`, `method`, `randomized` | Analysis settings. |

Cutoff laws always contain `support` and corresponding `probabilities`;
zero-probability points are omitted. Exact signed-rank laws also expose
`c_minus`, `c_plus`, `P_minus`, `P_plus`, and `tau`. For unit scores with
median-compatible calibration, tau is the probability needed to match the
binomial draw; the central-cutoff exception is explained in Section 5.

Statistical coverage requires the relevant centering, independence,
symmetry, or unimodality assumptions on batch estimates. The software does
not infer those assumptions from data.

## 10. Tests

Run the test suite from the directory containing these files:

```sh
python -m unittest -v test_cosi_ci
```

The suite includes independent statistical checks and live comparisons of
all returned R/Python fields. It tests minimality of B, binomial mixtures,
exact calibration against enumeration, signed-rank inversion over sign
configurations, and 225 concentration-bound comparisons with enumerated laws.
It also checks unit-score/median equality with ties and central cutoffs,
native seeds and shared plans, public kappa calibration, local RNG
restoration, K and gap sizes, row-order preservation, and insufficient samples.

The runner uses native Rscript when available. Otherwise install Node.js and
WebR, or set `COSI_WEBR_MODULE` to an absolute `webr.mjs` path:

```sh
npm install webr@0.6.0
COSI_REQUIRE_R=1 python -m unittest -v test_cosi_ci
```

Without either R backend, the R comparison class is skipped.
`COSI_REQUIRE_R=1` makes an unavailable R backend an error.
`test_results.txt` records the delivered test run with live R comparisons.
Integer decisions and labels are compared exactly. Floating fields allow
relative tolerance 3e-10 and absolute tolerance 3e-11; the diagnostic Chernoff
parameter allows absolute tolerance 5e-7 because the root solvers differ.
