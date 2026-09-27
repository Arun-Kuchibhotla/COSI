# COSI confidence intervals in R

The implementation consists of two base-R source files:

- `cosi_smallest_B.R`: deterministic and randomized binomial calibration.
- `cosi_ci.R`: data splitting and the four confidence-interval constructions.

```r
source("cosi_smallest_B.R")
source("cosi_ci.R")

set.seed(2026)
x <- rnorm(120)
fit <- cosi_ci(x, estimate = function(d) mean(d[, 1]))
fit                 # prints the interval with correct endpoint brackets
fit$CI              # named numeric vector: lower, upper
fit$B               # smallest B selected
fit$c               # binomial cutoff actually used
fit$ci_est          # the B batch estimates
```

## Interface

```r
cosi_ci(
  data,
  estimate,
  alpha = 0.05,
  gamma = 1/3,
  method = c("median", "symmetry", "normality", "unimodal-star"),
  randomized = TRUE,
  h = function(x) x,
  control = list()
)
```

`data` is converted with `as.matrix(data)`. `estimate` is called on a row
subset, always with `drop=FALSE`, and must return one finite real numeric
value. A vector is therefore passed to the estimator as a one-column matrix.
Missing values in the input can be handled by the estimator, for example
`mean(d[,1], na.rm=TRUE)`; a nonfinite estimator result is rejected with its
batch number.

`alpha` must be in `(0,1)` and `gamma` in `(0,1/2)`. The last interval is
strict for this function, as requested. The separate smallest-B helper still
supports `gamma=1/2`.

The accepted method names are those above. The exact aliases `"symmetric"`
and `"unimodal"` are also accepted.

## Batch construction

The function first calls the RNG-free `smallest_B_distribution(alpha,gamma)`.
The selected B is the smallest satisfying both

\[
q_B/B-1/2\le\gamma,\qquad q_B<B,
\quad q_B=\max\operatorname{supp}(c^*_{B,\alpha/2}).
\]

It then permutes the rows once and uses the requested balanced grouping:

```r
permutation <- sample.int(nn)
shuffled <- data[permutation, , drop = FALSE]
TMP <- split(seq_len(nn), sort(seq_len(nn) %% batches))
```

For the first three methods, `batches=B`. For unimodal-star, `batches=B+1`;
the first batch estimate is the anchor a, and the remaining B are W. The
unimodal method directly performs this B+1 split, without computing an unused
preliminary collection of B estimates. Every row is used exactly once and
batch sizes differ by at most one.

The sample-size check uses the number of batches actually needed. Its error
message reports that number, the sample size, alpha, and gamma; it does not
refer to an undefined `Delta`.

## Methods and calibration levels

| Method | Construction | Cutoff law |
|---|---|---|
| `median` | Batch-estimate order statistics | c* at level alpha/2 |
| `symmetry` | Weighted augmented Walsh quantiles with supplied h | Exact randomized kappa, or a conservative concentration cutoff |
| `normality` | Same Walsh construction with normal scores | Exact randomized kappa, or a conservative concentration cutoff |
| `unimodal-star` | Order statistics of expanded intervals around an independent batch anchor | c* at level alpha |

For median, the returned endpoints are

\[
[\widehat\theta_{(B-c)},\widehat\theta_{(c+1)}].
\]

For unimodal-star, the code uses the final formula in the request, including
its **alpha rather than alpha/2** cutoff level:

\[
L_j=\min(W_j,a)-(1+\sqrt2)|W_j-a|,\qquad
U_j=\max(W_j,a)+(1+\sqrt2)|W_j-a|,
\]

\[
[L_{(B-c)},U_{(c+1)}],\qquad c\sim c^*_{B,\alpha}.
\]

B is still selected using the alpha/2 constraint from the original
smallest-B problem. The actual unimodal cutoff is recalibrated at alpha; the
binomial distribution helper accepts any level in (0,1), including alpha>1/2.

At large alpha, the randomized unimodal cutoff can equal -1. The code uses
the standard extended order-statistic conventions x_(0)=-Inf and
x_(B+1)=Inf, so this outcome produces the empty interval `[Inf,-Inf]` with
`empty=TRUE`. It is not passed to R as an invalid array index.

For normality, h is fixed to

\[
h(i)=\Phi^{-1}\{1/2+i/[2(B+1)]\},\quad i=0,\ldots,B.
\]

Thus `normality` is the normal-score signed-rank construction requested here;
it does not call a t-test. A supplied `h` is used only by `symmetry`.

## Requirements on h

The implementation checks h at 0,1,...,B; h need not be vectorized. Its values
must be finite, nonnegative, and nondecreasing. In addition, **h(0)=0 is
required by the supplied formulas**. Indeed,

\[
\sum_{1\le i\le j\le B}\{h(j-i+1)-h(j-i)\}
=\sum_{r=1}^B h(r)-B h(0).
\]

Without h(0)=0, the total Walsh weight is not the stated N_h and the proposed
quantile calibration is inconsistent. The code reports this rather than
silently changing h or N_h.

The default `h=function(i)i` gives ordinary Wilcoxon scores. Flat portions
and zero weights are allowed. If every score is zero, the stated definitions
give the whole real line with probability 1-alpha and the empty set with
probability alpha. With `randomized=FALSE`, that case returns the whole line.

## Randomization

`randomized` controls critical-value randomization:

- For median and unimodal-star, TRUE draws from the appropriate binomial
  cutoff law; FALSE chooses its **largest** support value.
- With exact calibration for symmetry and normality, TRUE draws kappa with the specified tau; FALSE
  chooses its **smallest positive-probability support value**. Increasing
  kappa increases rejection, so the smaller endpoint is conservative. When
  the law degenerates at c_plus, c_plus itself is used.
- With concentration calibration, both settings use the same deterministic
  conservative kappa. No extra random draw is made for this cutoff.

Rows are shuffled and the independent U(-1,1) tie-breakers are generated even
with `randomized=FALSE`. Use `set.seed(...)` before the call to reproduce the
whole construction. The deterministic smallest-B helper consumes no random
numbers, and unused binomial cutoffs are not drawn for signed-rank methods.

## Walsh quantiles and endpoint inclusion

The code sorts (W_i,U_i) lexicographically, constructs all i<=j augmented
Walsh averages and their lag weights, and sorts those averages
lexicographically. It computes Q_h(gamma+) using **strictly greater cumulative
weight**, rather than adding one to gamma. The latter shortcut is invalid
for general real-valued scores, including the normal scores.

If the augmented lower endpoint is (l,u_l), its real endpoint l is included
exactly when u_l<0. If the augmented upper endpoint is (r,u_r), r is included
exactly when u_r>0. Consequently the result includes:

- `CI`: real lower and upper endpoints;
- `closed`: two logical values specifying endpoint inclusion;
- `empty`: whether the represented set is empty;
- `augmented_endpoints`: the two augmented endpoints for signed-rank methods.

An interval with equal real endpoints can be a singleton or empty, depending
on the tie-breakers. `print(fit)` and `contains_cosi_ci(fit,t)` respect this.
Replacing every result by a closed interval using only `fit$CI` would change
the specified construction at atoms.

```r
contains_cosi_ci(fit, c(-1,0,1))
fit$closed
fit$empty
```

## Other returned fields

| Field | Meaning |
|---|---|
| `B`, `batches` | Selected B and actual number of batches |
| `ci_est` | B estimates used in the final construction |
| `all_estimates` | Includes the anchor estimate for unimodal-star |
| `anchor` | First estimate for unimodal-star; otherwise NULL |
| `c`, `cutoff_distribution` | Binomial cutoff and its law, where used |
| `kappa`, `kappa_distribution` | Signed-rank cutoff and its law, where used |
| `split_calibration` | B-selection cutoff law at alpha/2 for every method |
| `h_values`, `N_h`, `uniforms` | Scores, total score, and tie-breakers |
| `batch_sizes`, `batch_rows` | Batch sizes and original row indices |
| `L`, `U` | The expanded per-batch endpoints for unimodal-star |

With exact calibration, `kappa_distribution` also records c_minus, c_plus,
their rejection probabilities, tau, `exact=TRUE`, and the algorithm used.
With concentration calibration, it instead records:

- `engine="concentration"` and `exact=FALSE`;
- singleton `support` and `probabilities=1`;
- `bound_used`: the inequality giving the smallest computed radius;
- `candidate_radii`: the radii from all four inequalities;
- `radius` and `radius_fraction`: the selected radius in score units and
  after division by N_h/2, before the numerical cushion;
- `rejection_probability_bound` and its logarithm, evaluated for the final
  cutoff with the Walsh mass tolerance accounted for;
- `lambda`: the finite Chernoff parameter for the normalized statistic;
- `numerical_cushion` and `fallback_reason`.

The exact rejection probability is not computed on this path, so
`rejection_probabilities=NA`. The returned bound must not be interpreted as
the actual rejection probability. Printing the result explicitly identifies
the conservative calibration and its bound.

## Exact calibration and concentration fallback

The automatic calibration first tries the following fast exact routes:

1. A closed form handles the smallest critical value, and the all-zero case.
2. Linear Wilcoxon scores use base R's signed-rank distribution routines,
   with neighboring CDF checks to enforce the quantile definition.
3. Lattice scores use a dynamic program when its state count and work fit.
4. Small general problems enumerate subset sums.
5. Remaining general problems with B < `concentration_B` use
   meet-in-the-middle: only the two half-distributions are stored and pair-sum
   counts select the two critical support values.
6. Remaining general problems with B >= `concentration_B` use the
   concentration method below. The default threshold is **40**, inclusive.
   A state/counting limit or a pair-sum iteration limit also triggers this
   fallback, even below 40. A discrete tail underflow uses log-scale
   concentration calibration.

These algorithms are exact distributional calculations **in double precision**,
not arbitrary-precision certificates. The calculations retain the numerical
limitations of R's probability routines and floating-point subset sums.

The general real-score state count can grow exponentially, but reaching that
limit now returns a conservative interval rather than an error. For example,
alpha=0.05 and gamma=0.1 select B=90; the normal-score method now succeeds
using concentration calibration. Ordinary Wilcoxon scores and manageable
lattice scores can still use exact calibration above 40. Small problems with
at most 14 positive scores can also use enumeration above this B threshold.

### Four concentration bounds, with the tightest computed radius selected

Let N=sum h(i)>0, mu=N/2, w_i=h(i)/N, and

\[
T=(G_{B,h}-\mu)/\mu=\sum_i w_i\varepsilon_i,
\qquad s^2=\sum_i w_i^2,\qquad L=\log(2/\alpha).
\]

Each of the following radii r satisfies P(|T|>=r)<=alpha:

1. **Hoeffding:** r_H=s sqrt(2L).
2. **Exact-MGF Chernoff:** for every finite lambda>0,
   r_C(lambda)=[L+sum log cosh(lambda w_i)]/lambda.
3. **Sharp Gaussian comparison:**
   r_G=s Phi^{-1}(1-alpha/(2 C_*)), where
   C_*={4[1-Phi(sqrt(2))]}^{-1}, approximately 3.178.
   This is a rigorous Rademacher-tail inequality, not a normal approximation.
4. **Chebyshev:** r_V=s/sqrt(alpha). This can help at large alpha.

The code numerically minimizes the Chernoff expression over finite lambda.
It brackets the derivative's zero and solves a one-dimensional equation,
starting with the Hoeffding-optimal lambda=sqrt(2L)/s. Each objective or
derivative evaluation takes O(B) operations. If root finding does not
converge, the best evaluated finite lambda remains a valid choice. The
starting lambda already gives a radius no larger than Hoeffding's in exact
arithmetic.

The chosen radius is the minimum of the four computed radii. No multiple-test
correction is required: all bounds concern the same fixed score law and are
chosen independently of the data. This is the tightest among these implemented
bounds, not a claim of universal optimality among all concentration methods.
The constant C_* is optimal specifically for a universal multiplicative
comparison with the standard Gaussian tail.

Set kappa=mu(1-r), decreased by a small numerical cushion. When the resulting
cutoff is negative, use the sentinel -1 and return the whole real line. The
usual Q_h(kappa+) and Q_h(N-kappa) endpoints are then computed. There is no
need to locate kappa in the discrete subset-sum support: the weighted Walsh
quantiles accept real arguments. With the assumptions of the signed-rank
construction, coverage is **at least 1-alpha**, rather than necessarily equal
to 1-alpha. There is no simulation-based or uncorrected normal approximation.

The implementation normalizes scores to avoid overflow in sum h(i)^2, uses
log probabilities for very small alpha, evaluates log cosh stably, and never
substitutes an infinite Chernoff parameter at a discrete support boundary.
The numerical cushion covers the configured inward Walsh mass tolerance;
the reported bound is checked after this adjustment. The calculations use
ordinary double precision, not certified interval arithmetic. Very small
score scales relative to the absolute component of the mass tolerance can
produce the whole line conservatively.

References:
[Hoeffding's Rademacher inequality](https://ocw.mit.edu/courses/18-465-topics-in-statistics-statistical-learning-theory-spring-2007/8fc2f7e9e92e6ed03edf22dd35f37eb5_lecture07.pdf);
[Bentkus and Dzindzalieta (2015), Theorem 1.1](https://arxiv.org/abs/1307.3451).

### Remaining resource controls

The optional `control` list contains:

| Option | Default | Meaning |
|---|---:|---|
| `max_states` | 2^20 | Maximum lattice states or half-subset sums |
| `max_walsh` | 2,000,000 | Maximum number of Walsh pairs |
| `max_dp_work` | 50,000,000 | Limit on scores times lattice-state count |
| `max_iterations` | 256 | Pair-sum selection limit; also limits Chernoff bracketing/root solving |
| `concentration_B` | 40 | Use concentration for remaining general problems at or above this B |
| `mass_tol` | 64 times machine epsilon | Relative tolerance for real-score cumulative weights |

The separate Walsh-pair memory limit still applies. Increasing
`concentration_B` permits trying exact general calibration at larger B, but
exceeding its state/counting/iteration limits still produces a conservative
fallback. Invalid inputs and insufficient rows still raise informative errors.

For noninteger scores, a small tolerance is used only when comparing
cumulative Walsh weights with a score threshold, to absorb telescoping-sum
roundoff. It does not merge estimate values or alter their lexicographic
ordering. Integer cumulative weights use zero tolerance while exactly
representable. Set `control=list(mass_tol=0)` for literal floating-point mass
comparisons; score sums that should agree mathematically can then differ by
roundoff. Custom scores with distinctions near machine precision require
particular care or higher-precision arithmetic.

## Additional examples

```r
set.seed(9)
d <- cbind(x = rnorm(240), y = rnorm(240))
est <- function(z) mean(z[, 1])

cosi_ci(d, est, method = "median", randomized = FALSE)
cosi_ci(d, est, method = "symmetry", h = function(i) i)
cosi_ci(d, est, method = "symmetry", h = function(i) i^2)
cosi_ci(d, est, method = "normality")
cosi_ci(d, est, method = "unimodal-star")

# Previously blocked by general real-score calibration: now returns a CI.
large <- cosi_ci(d, est, alpha = .05, gamma = .1, method = "normality")
large$B                                      # 90
large$kappa_distribution$engine               # "concentration"
large$kappa_distribution$bound_used
large$kappa_distribution$rejection_probability_bound  # <= .05
large$kappa_distribution$candidate_radii

# Multivariate rows, but one scalar estimate per batch:
slope <- function(z) unname(coef(lm(z[, 2] ~ z[, 1]))[2])
cosi_ci(d, slope, method = "median")

# Distribution-only helper, with no RNG use:
smallest_B_distribution(alpha = 0.05, gamma = 1/3)

# The previous randomized helper remains available with the same API:
smallest_B(alpha = 0.05, gamma = 1/3, seed = 42)
```

The code implements the requested constructions. Coverage claims require the
appropriate centering, independence, symmetry, or unimodality assumptions on
the batch estimates; these assumptions are not inferred from the data.

## Validation

The R source files were executed in R 4.6.0 through WebR. Checks covered all
four end-to-end methods, matrix-preserving balanced splits, reproducibility,
RNG-free calibration, general binomial levels including support at -1,
conservative cutoff choices, and input/resource-limit errors. Exhaustive
conditional sign configurations checked both signed-rank rejection
probabilities and the resulting intervals, including tied estimates, zero
estimates, flat score functions, and degenerate scores. Meet-in-the-middle
critical values were compared with full subset enumeration.

The concentration update was additionally checked against 225 fully enumerated
Bernoulli score laws and 54 exhaustive conditional-sign interval checks, including
ties and zero estimates. Checks covered the inclusive B=40 switch, the previously
failing B=90 normal-score call, retained exact Wilcoxon/lattice paths, resource
fallbacks, deterministic cutoff behavior, scale robustness, and extreme alpha.

API references: [R signed-rank distribution](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/SignRank.html)
and [R binomial distribution](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html).
