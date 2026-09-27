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
\qquad q_B \le B-1,
\qquad q_B=\max\text{supp}(c^*_{B,\alpha/2}).
$$

Here the independently randomized binomial cutoff obeys
$\mathbb{P}\{\text{Bin}(B,1/2)\le c^*_{B,\alpha/2}\}=1-\alpha/2$.
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
G=\sum_{i=1}^B h(i)\xi_i,\qquad \xi_i\overset{\mathrm{iid}}{\sim}\text{Ber}(1/2),
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
computes a deterministic conservative real cutoff. Exact Wilcoxon and
manageable lattice cases remain available above 40. State/counting/iteration
limits can also trigger concentration calibration below 40.

For $N_h>0$, write

$$
G=\sum_{i=1}^B h(i)\xi_i,\qquad
S=G-\frac{N_h}{2}=\sum_{i=1}^B b_i\varepsilon_i,\qquad
b_i=\frac{h(i)}2,\qquad v=\sum_i b_i^2,
$$

where the independent signs take values $-1$ and $1$ with equal probability.
Let $m$ be the number of positive scores, and let $L=\log(2/\alpha)$.
A radius $r$ satisfying $\mathbb{P}(|S|\ge r)\le\alpha$ gives the cutoff

$$
\kappa=\frac{N_h}{2}-r.
$$

The implementation computes the following candidates and uses the smallest:

| `bound_used` | Candidate and purpose |
|---|---|
| `chernoff` | Exact MGF evaluated at an analytical starting parameter and three fixed refinements. |
| `binomial_p2` | Matching-variance binomial comparison, inverted by explicit formulas on finitely many lattice cells. |
| `gaussian` | Sharp constant-factor Gaussian comparison for Rademacher sums. |
| `pinelis` | Explicit conservative inversion of an asymptotically Gaussian tail bound. |
| `hoeffding` | Sub-Gaussian radius $\sqrt{2vL}$. |
| `chebyshev` | Variance radius $\sqrt{v/\alpha}$. |
| `tomaszewski` | A radius strictly greater than $\sqrt v$, available when $\alpha\ge1/2$. |

These are finite-sample inequalities, not normal approximations. Section 11
gives their derivations, endpoint conditions, comparisons, and references.
There is no continuous root search or unrestricted optimization in the
concentration path. This collection is not claimed to be optimal among all
possible concentration inequalities.

To avoid overflow for very large scores and underflow for very small scores,
the code works with $T=S/(N_h/2)=\sum_i w_i\varepsilon_i$, where
$w_i=h(i)/N_h$. It converts candidate radii back to the units of $S$ when
returning `radius` and `candidate_radii`. The field `radius_fraction` and the
MGF parameters `lambda` and `lambda_candidates` refer to $T$. Thus a returned
MGF parameter divided by $N_h/2$ is the corresponding parameter for $S$.

The selected cutoff is decreased by

$$
\left(\mathtt{mass\_tol}+256\varepsilon_{\mathrm{mach}}\right)
\max(1,N_h),
$$

where $\varepsilon_{\mathrm{mach}}$ is machine epsilon. A negative cutoff is
replaced by -1, yielding the whole real line. The final rejection bound is
checked after accounting for the Walsh cumulative-weight tolerance. The
binomial calculation uses log probabilities and a conditional-variance
recurrence, and inflates its second moment by the factor
$\exp\{256\varepsilon_{\mathrm{mach}}(m+1)\}$ before inversion and checking.
These safeguards use floating-point arithmetic; they are not
arbitrary-precision certificates. The absolute part of the cutoff cushion
can lead to the whole real line when all scores are extremely small.

For a concentration result, `distribution` contains:

| Field | Meaning |
|---|---|
| `exact` | FALSE: calibration uses inequalities rather than a computed exact score law. |
| `bound_used` | Candidate giving the smallest computed radius. |
| `candidate_radii` | All seven candidate radii in the units of $S$; infinity marks an unavailable candidate. |
| `radius`, `radius_fraction` | Chosen radius in the units of $S$ and $T$, respectively. |
| `lambda`, `lambda_candidates` | Final Chernoff parameter and all evaluated parameters, in the normalization for $T$. |
| `chernoff_refinements` | Number of completed fixed refinements, at most three. |
| `chernoff_gap_bound` | Upper bound on the difference from the optimal Chernoff radius, in the units of $S$, with a numerical allowance; infinity when a finite optimizer does not exist. |
| `rejection_probability_bound`, `log_rejection_probability_bound` | Checked upper bound for the implemented cutoff, and its logarithm. |
| `numerical_cushion`, `fallback_reason` | Cutoff cushion and reason concentration calibration was selected. |

The actual rejection probability is not computed and is returned as NA/NaN.
The checked rejection bound can come from a different candidate than
`bound_used`, because all available bounds are evaluated at the final cutoff.
On exact-distribution paths, `exact=TRUE` describes how the law was computed;
it does not imply exact interval coverage under ties or median-compatible
central randomization.

The `control` list/dictionary accepts:

| Option | Default | Purpose |
|---|---:|---|
| `max_states` | $2^{20}$ | Limit lattice states and half-subset sums. |
| `max_dp_work` | 50,000,000 | Limit lattice work. |
| `max_walsh` | 2,000,000 | Limit Walsh pairs for general-score interval construction. |
| `max_iterations` | 256 | Limit exact pair-sum selection; concentration calibration always uses at most three refinements. |
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
python -m unittest -v test_cosi_ci test_calibration_analytic
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
COSI_REQUIRE_R=1 python -m unittest -v test_cosi_ci test_calibration_analytic
```

Without either R backend, the R comparison class is skipped.
`COSI_REQUIRE_R=1` makes an unavailable R backend an error.
`test_results.txt` records the delivered test run with live R comparisons.
Integer decisions and labels are compared exactly. Floating fields allow
relative tolerance 3e-10 and absolute tolerance 3e-11. The analytical
calibration tests also compare the finite checks with a numerical optimizer,
verify the error bound, check the binomial inversion independently, and
exercise very small probabilities and nearly equal scores in both languages.


## 11. Mathematical details and references

The formulas in this section use $S$, $b_i$, $v$, $m$, and $L$ from Section 7.
Zero scores do not affect the law. The all-zero case is handled exactly
before concentration calibration.

### 11.1. Chernoff radii and the finite optimizer

Independence gives the exact log MGF

$$
\psi(\lambda)=\log\mathbb{E}e^{\lambda S}
=\sum_i\log\cosh(\lambda b_i).
$$

Symmetry and the exponential Markov inequality imply, for every finite
$\lambda>0$,

$$
\mathbb{P}(|S|\ge r)
\le 2\exp\{\psi(\lambda)-\lambda r\},\qquad
r(\lambda)=\frac{L+\psi(\lambda)}{\lambda}.
$$

Define

$$
H(\lambda)=\lambda\psi'(\lambda)-\psi(\lambda),\qquad
\psi'(\lambda)=\sum_i b_i\tanh(\lambda b_i),\qquad
\psi''(\lambda)=\sum_i b_i^2\text{sech}^2(\lambda b_i).
$$

Then $r'(\lambda)=\{H(\lambda)-L\}/\lambda^2$ and
$H'(\lambda)=\lambda\psi''(\lambda)>0$. Also $H(0)=0$ and
$H(\infty)=m\log2$. Hence, when $L<m\log2$, there is a unique finite
minimizer $\lambda_*$ characterized by $H(\lambda_*)=L$.

For arbitrary real scores this characterization does not provide a general
closed-form inverse. For equal positive scores $b_i=b$, optimizing the
tail bound at a fixed $x\in(0,mb)$ gives

$$
\lambda_x=\frac1b\text{arctanh}\left(\frac{x}{mb}\right),\qquad
\inf_{\lambda>0}e^{\psi(\lambda)-\lambda x}
=\exp\left\{-m I\left(\frac{x}{mb}\right)\right\},
$$

$$
I(u)=\frac{(1+u)\log(1+u)+(1-u)\log(1-u)}2.
$$

Obtaining the corresponding quantile still requires inverting $I$. Equal
scores are handled through exact-distribution routines when feasible.

If $L\ge m\log2$, the infimum of $r(\lambda)$ is the support endpoint
$\sum_i b_i$ and has no finite minimizer. A limiting radius must not be
substituted into an inclusive tail event without checking its atom:

$$
\mathbb{P}\left(|S|\ge\sum_i b_i\right)=2^{1-m}.
$$

For the concentration helper this case retains a finite Hoeffding parameter
and performs no refinements. The exact calibration dispatcher separately
handles the randomized mixture of the -1 and 0 cutoff endpoints.

### 11.2. Analytical starting parameter

The Hoeffding parameter is $\lambda_H=\sqrt{2L/v}$. A sharper explicit
starting value, for $L<m\log2$, is

$$
t=\tan\left(\frac{L}{m\sqrt3}\right),\qquad
\lambda_0=\sqrt{\frac{2\sqrt3\,m\,t}{v(1-\sqrt3\,t)}}.
$$

The denominator is positive because $L/m<\log2<\pi/(2\sqrt3)$. In exact
arithmetic,

$$
\lambda_H\le\lambda_0\le\lambda_*,\qquad
r(\lambda_0)\le r(\lambda_H).
$$

To prove this, let $g(z)=z\tanh z-\log\cosh z$. Since

$$
g'(z)=\frac{z}{\cosh^2z},\qquad
\cosh^2z\ge1+z^2+\frac{z^4}{3},
$$

integration gives $g(z)\le f(z^2)$ for $z\ge0$, where

$$
f(u)=\int_0^u\frac{ds}{2(1+s+s^2/3)}
=\sqrt3\left[\arctan\left(\frac{2u+3}{\sqrt3}\right)-\frac\pi3\right].
$$

The function $f$ is increasing and concave and satisfies $f(u)\le u/2$.
Jensen's inequality therefore gives

$$
H(\lambda)=\sum_i g(\lambda b_i)
\le m f\left(\frac{\lambda^2v}{m}\right).
$$

Solving $m f(\lambda_0^2v/m)=L$ gives the stated parameter and ordering.
This is a direct analytical derivation used by this implementation; it is
not attributed to a theorem in the cited papers.

### 11.3. Three fixed refinements and an error bound

Starting at $\lambda_0$, the code performs three explicit updates:

$$
\lambda_{j+1}
=\sqrt{\lambda_j^2+
\frac{2\{L-H(\lambda_j)\}}{\psi''(\lambda_j)}},\qquad j=0,1,2.
$$

These are Newton steps in the variable $\lambda^2$, performed a fixed
number of times rather than until a numerical convergence criterion holds.
The calculation takes $O(B)$ time for each of the four MGF evaluations.

Because $\psi''$ decreases on $[0,\infty)$,

$$
H(\lambda_{j+1})-H(\lambda_j)
=\int_{\lambda_j}^{\lambda_{j+1}}u\psi''(u)\,du
\le\frac{\psi''(\lambda_j)}2(\lambda_{j+1}^2-\lambda_j^2)
=L-H(\lambda_j).
$$

It follows inductively that

$$
\lambda_j\le\lambda_{j+1}\le\lambda_*,\qquad
r(\lambda_{j+1})\le r(\lambda_j).
$$

At the optimum, $r(\lambda_*)=\psi'(\lambda_*)$. Monotonicity of $\psi'$
also yields the computable bound

$$
0\le r(\lambda_j)-r(\lambda_*)
\le r(\lambda_j)-\psi'(\lambda_j)
=\frac{L-H(\lambda_j)}{\lambda_j}.
$$

The implementation clips a numerically negative residual to zero, retains
the last finite parameter if a refinement cannot be evaluated, and adds a
floating-point allowance to the reported error bound. The exponential
Markov bound remains valid at any finite parameter; its validity never
depends on convergence of the refinements. The error estimate concerns the
optimal Chernoff radius, not the unknown exact quantile.

### 11.4. Binomial positive-part second-moment comparison

Bentkus and Juskevicius [3, Theorem 3 and Section 2], using the comparison
developed by Bentkus [4], imply the following. Set

$$
y=\max_i b_i,\qquad
p=\frac{v}{my^2+v},\qquad d=y+\frac{v}{my},\qquad
Y\sim\text{Bin}(m,p),\qquad T=d(Y-mp).
$$

Then $\text{Var}(T)=v$ and, for all real $u$,

$$
\mathbb{E}(S-u)_+^2\le\mathbb{E}(T-u)_+^2.
$$

Thus

$$
\mathbb{P}(S\ge x)\le P_2(T;x)
:=\inf_{u<x}\frac{\mathbb{E}(T-u)_+^2}{(x-u)^2}.
$$

Our summands have skewness zero and kurtosis one. After division by $y$,
their known variances are at most one, so the skewness and kurtosis bounds
in [3] do not further reduce these variance inputs. The useful additional
candidate is the positive-part moment comparison itself. Its explicit
binomial form is also studied in [5].

For $\beta=\alpha/2>p^m$, define

$$
q_\beta=\min_{u<m}\left\{u+
\sqrt{\frac{\mathbb{E}(Y-u)_+^2}{\beta}}\right\},\qquad
r_{\mathrm{bin}}=d(q_\beta-mp).
$$

The minimizing shift gives $P_2(T;r_{\mathrm{bin}})\le\beta$.
Symmetry of $S$ therefore gives
$\mathbb{P}(|S|\ge r_{\mathrm{bin}})\le\alpha$.

There is no continuous optimization in the implementation. On each cell
$u\in[k-1,k]$, define

$$
P_k=\mathbb{P}(Y\ge k),\qquad
\mu_k=\mathbb{E}(Y\mid Y\ge k),\qquad
s_k^2=\text{Var}(Y\mid Y\ge k).
$$

The first cell has lower endpoint $-\infty$. The objective is

$$
u+\sqrt{\frac{P_k}{\beta}\{(u-\mu_k)^2+s_k^2\}}.
$$

If $P_k>\beta$ and $s_k^2>0$, its unconstrained minimizer and minimum are

$$
u_k=\mu_k-\sqrt{\frac{\beta s_k^2}{P_k-\beta}},\qquad
q_k=\mu_k+\sqrt{s_k^2\left(\frac{P_k}{\beta}-1\right)}.
$$

Clip $u_k$ to its cell and evaluate the objective at the clipped value.
If $P_k\le\beta$, use the left endpoint. Taking the minimum over all cells
gives $q_\beta$. At $\beta=p^m$, the mathematical endpoint is $q_\beta=m$;
below it, an inclusive tail guarantee needs a threshold strictly above $m$.
The implementation omits this binomial candidate when the target is at or
below the guarded top-atom probability; the other candidates remain valid.

For reference, the raw tail moments can be expressed through binomial tails:

$$
\mathbb{E}[Y\mathbf{1}\{Y\ge k\}]
=mp\,\mathbb{P}\{\text{Bin}(m-1,p)\ge k-1\},
$$

$$
\mathbb{E}[Y(Y-1)\mathbf{1}\{Y\ge k\}]
=m(m-1)p^2\,\mathbb{P}\{\text{Bin}(m-2,p)\ge k-2\}.
$$

The second identity is used only for $m\ge2$. To avoid cancellation, the
code instead scans downward from $k=m$. Let $\delta_{k+1}=\mu_{k+1}-(k+1)$,
$a_k=\mathbb{P}(Y=k)$, and $t_k=P_{k+1}/P_k$. Then

$$
P_k=a_k+P_{k+1},\qquad
\delta_k=t_k(\delta_{k+1}+1),\qquad
s_k^2=t_ks_{k+1}^2+t_k(1-t_k)(\delta_{k+1}+1)^2.
$$

All variance contributions are nonnegative. Probabilities are accumulated
in log space, including the target probability, so even an underflowing
floating-point value of $\alpha/2$ is not needed. The retained shift and
guarded moment are used to check the tail bound at the final cutoff.
The computation takes $O(m)$ time and $O(m)$ storage for binomial log PMFs.

### 11.5. Gaussian comparisons and other explicit candidates

Bentkus and Dzindzalieta [1, Theorem 1.1], also presented in the thesis [2],
give the sharp constant-factor comparison

$$
\mathbb{P}(S/\sqrt v\ge x)\le C_*\overline\Phi(x),\qquad
C_*=[4\overline\Phi(\sqrt2)]^{-1}\approx3.1786556.
$$

Consequently,

$$
r_G=\sqrt v\,\Phi^{-1}\left(1-\frac{\alpha}{2C_*}\right).
$$

The constant is optimal uniformly over coefficient vectors and thresholds;
two equal nonzero coefficients at $x=\sqrt2$ attain it. This does not imply
that its quantile bound is best for each particular score vector.

Pinelis [6, Theorem 1.1] proves, for $x>0$,

$$
\mathbb{P}(S/\sqrt v\ge x)\le Q(x)
:=\overline\Phi(x)+C_P\frac{\phi(x)}{9+x^2},\qquad
C_P=5\sqrt{2\pi e}\{2\Phi(1)-1\}\approx14.1068613.
$$

An explicit conservative inversion follows from the Mills-ratio bound
$\phi(x)/\overline\Phi(x)\le x+x^{-1}$. Let

$$
A(x)=1+C_P\frac{x+x^{-1}}{9+x^2},\qquad
z=\Phi^{-1}(1-\alpha/2),\qquad
r_P=\sqrt v\,\Phi^{-1}\left(1-\frac{\alpha}{2A(z)}\right).
$$

Indeed, $Q(x)\le A(x)\overline\Phi(x)$ and

$$
A'(x)=-C_P\frac{(x^2-3)^2}{x^2(9+x^2)^2}\le0.
$$

Since $r_P/\sqrt v\ge z$, it follows that $2Q(r_P/\sqrt v)\le\alpha$.
For numerical stability when $\alpha>0.9$, the code substitutes the lower
bound $z_0=\sqrt{\pi/2}(1-\alpha)\le z$ in $A$. Monotonicity makes the
result conservative. Normal quantiles otherwise use log tail probabilities.

Keller and Klein [7] prove Tomaszewski's conjecture:

$$
\mathbb{P}(|S|\le\sqrt v)\ge\frac12.
$$

For $\alpha\ge1/2$, any radius strictly greater than $\sqrt v$ is therefore
valid for the inclusive event. The code uses
$r_T=(1+512\varepsilon_{\mathrm{mach}})\sqrt v$ and checks that the final
radius is strictly larger than $\sqrt v$. Replacing this by an inclusive
tail bound at exactly $\sqrt v$ is invalid, as a single random sign shows.

The remaining candidates are Hoeffding's bound [11] and Chebyshev's bound:

$$
r_H=\sqrt{2v\log(2/\alpha)},\qquad r_V=\sqrt{v/\alpha}.
$$

### 11.6. Why the Lambert-W bound is not an additional candidate

Proposition 3.1 of Pinelis's preprint [8] analytically optimizes a
Pinelis-Utev MGF bound. For our symmetric sum, put

$$
\rho=\frac{\sum_i b_i^3}{2yv},\qquad 0<\rho\le\frac12.
$$

That log-MGF bound and its optimizing parameter at a fixed threshold $x$ are

$$
\psi_{\mathrm{PU}}(\lambda)
=\frac{(1-\rho)v\lambda^2}{2}
+\frac{\rho v}{y^2}(e^{\lambda y}-1-\lambda y),
$$

$$
A_x=\frac{\rho+xy/v}{1-\rho},\qquad
\lambda_x=\frac1y\left[
A_x-W_0\left(\frac{\rho}{1-\rho}e^{A_x}\right)\right],
$$

where $W_0$ is the principal Lambert-W function. However,

$$
\psi(\lambda)\le\frac{v\lambda^2}{2}
\le\psi_{\mathrm{PU}}(\lambda),\qquad \lambda\ge0,
$$

because $e^u-1-u\ge u^2/2$ for $u\ge0$. After optimization this yields

$$
\inf_{\lambda>0}e^{\psi(\lambda)-\lambda x}
\le e^{-x^2/(2v)}
\le\inf_{\lambda>0}e^{\psi_{\mathrm{PU}}(\lambda)-\lambda x}.
$$

It cannot improve the Hoeffding candidate for this problem. The formula
also removes only the MGF-parameter optimization: the resulting tail bound
still has to be inverted to obtain a quantile.

### 11.7. Additional recent results and scope

Jakimiuk [9] strengthens Gaussian moment comparisons by a deficit depending
on $\sum_i(b_i/\sqrt v)^4$. Such results can supply moment-based Markov
bounds, but do not automatically dominate the implemented tail comparisons.
They are not additional calibration candidates in this implementation.

Hollom and Portier [10] give sharp lower bounds on Rademacher tail
probabilities, including the standardized two-sided value $7/32$ at one
standard deviation. Their anti-concentration results do not directly give
the upper-tail guarantee required for a conservative cutoff.

### 11.8. Numerical comparison

The following radii are divided by $\sqrt v$, with $B=40$. The first column
of radii evaluates the exact MGF at $\lambda_H$; it is not the larger
Hoeffding radius itself. Normal scores are those in Section 4. Nearly equal
scores mean $h(i)=0.9+0.1\sqrt{i/40}$.

| Scores | Level | Check at $\lambda_H$ | Analytical start plus three refinements | Optimized Chernoff benchmark | Binomial second-moment bound |
|---|---:|---:|---:|---:|---:|
| Normal | 0.05 | 2.624462 | 2.608949 | 2.608949 | 2.666355 |
| Normal | 0.000001 | 4.865173 | 4.543946 | 4.543946 | 6.081026 |
| Nearly equal | 0.000001 | 5.112439 | 5.027974 | 5.027974 | 4.900019 |

The benchmark optimizer is used only for testing. The table illustrates
accuracy for these inputs, not a universal error bound after three steps;
the computable bound in Section 11.3 controls that error. At level 0.05 the
Gaussian radius is about 2.415122 standard deviations and improves all the
normal-score candidates in the first row. The binomial comparison is useful
for the nearly equal scores in the last row. There is no uniform winner.

### 11.9. References

1. Bentkus, V. K., and Dzindzalieta, D. (2015). *A tight Gaussian bound for
   weighted sums of Rademacher random variables*. Bernoulli 21(2), 1231-1237.
   [Author preprint](https://arxiv.org/abs/1307.3451).
2. Dzindzalieta, D. (2014). *Tight Bernoulli Tail Probability Bounds*.
   Doctoral dissertation, Vilnius University. In particular, Theorems 1,
   11, and 12.
   [Thesis](https://www.mii.lt/files/doc/lt/doktorantura/apgintos_disertacijos/mii_dis_2014_dzindzalieta.pdf).
3. Bentkus, V., and Juskevicius, T. (2008). *Bounds for tail probabilities of
   martingales using skewness and kurtosis*. Lithuanian Mathematical Journal
   48(1), 30-37. Theorem 3 and Section 2.
   [DOI](https://doi.org/10.1007/s10986-008-0003-8).
4. Bentkus, V. (2004). *On Hoeffding's inequalities*. Annals of Probability
   32(2), 1650-1673.
5. Bentkus, V., Kalosha, N., and van Zuijlen, M. (2006). *On domination of
   tail probabilities of (super)martingales: explicit bounds*. Lithuanian
   Mathematical Journal 46(1), 3-54.
6. Pinelis, I. (2012). *An asymptotically Gaussian bound on the Rademacher
   tails*. Electronic Journal of Probability 17, paper 35, 1-22. Theorem 1.1.
   [DOI](https://doi.org/10.1214/EJP.v17-2026),
   [author preprint](https://arxiv.org/abs/1007.2137).
7. Keller, N., and Klein, O. (2022). *Proof of Tomaszewski's conjecture on
   randomly signed sums*. Advances in Mathematics 407, 108558.
   [DOI](https://doi.org/10.1016/j.aim.2022.108558),
   [author preprint](https://arxiv.org/abs/2006.16834).
8. Pinelis, I. (2014). *On the Bennett-Hoeffding inequality*. Annales de
   l'Institut Henri Poincare, Probabilites et Statistiques 50(1), 15-27.
   The Lambert-W formula cited here is Proposition 3.1 of the
   [2009 preprint](https://arxiv.org/abs/0902.4058).
9. Jakimiuk, J. (2025). *Stability of Khintchine inequalities with optimal
   constants between the second and the p-th moment for p >= 3*.
   [Preprint](https://arxiv.org/abs/2503.07001).
10. Hollom, L., and Portier, J. (2025). *Tight Anti-Concentration of
    Rademacher Sums*. Random Structures & Algorithms.
    [DOI](https://doi.org/10.1002/rsa.70024),
    [author preprint](https://arxiv.org/abs/2306.07811).
11. Hoeffding, W. (1963). *Probability inequalities for sums of bounded
    random variables*. Journal of the American Statistical Association
    58, 13-30.

