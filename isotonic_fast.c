/* Optional batch isotonic predictor. Sampling and all randomization stay in R. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

SEXP iso_estimates_many_native(SEXP x_, SEXP y_, SEXP indices_, SEXP grid_) {
    if (TYPEOF(x_) != REALSXP || TYPEOF(y_) != REALSXP ||
        TYPEOF(grid_) != REALSXP || TYPEOF(indices_) != INTSXP)
        error("Native isotonic inputs have unexpected storage types.");

    const int n = LENGTH(x_);
    const int K = LENGTH(grid_);
    SEXP dimensions = getAttrib(indices_, R_DimSymbol);
    if (LENGTH(dimensions) != 2)
        error("indices must be an integer matrix with one sample per column.");
    const int m = INTEGER(dimensions)[0];
    const int R = INTEGER(dimensions)[1];
    if (n < 1 || LENGTH(y_) != n || K < 1 || m < 1 || R < 1)
        error("Native isotonic inputs must have positive compatible dimensions.");

    const double *x = REAL(x_);
    const double *y = REAL(y_);
    const double *grid = REAL(grid_);
    const int *indices = INTEGER(indices_);
    for (int i = 0; i < n; ++i) {
        if (!R_FINITE(x[i]) || !R_FINITE(y[i]) ||
            (i > 0 && x[i] < x[i - 1]))
            error("x must be sorted; x and y must be finite.");
    }
    for (int k = 0; k < K; ++k)
        if (!R_FINITE(grid[k])) error("grid must be finite.");

    SEXP answer_ = PROTECT(allocMatrix(REALSXP, K, R));
    double *answer = REAL(answer_);
    double *level = (double *) R_alloc(m, sizeof(double));
    double *left = (double *) R_alloc(m, sizeof(double));
    int *weight = (int *) R_alloc(m, sizeof(int));

    for (int r = 0; r < R; ++r) {
        if (r % 64 == 0) R_CheckUserInterrupt();
        const int *draw = indices + (R_xlen_t) m * r;
        for (int i = 0; i < m; ++i) {
            if (draw[i] < 1 || draw[i] > n || (i > 0 && draw[i] < draw[i - 1]))
                error("Each sampled index column must be sorted and lie in 1:n.");
        }

        int blocks = 0;
        int i = 0;
        while (i < m) {
            /* Aggregate all observations with the same covariate, including
               repeated pairs from the bootstrap and ties in the original X. */
            const double location = x[draw[i] - 1];
            double total = 0.0;
            int count = 0;
            do {
                total += y[draw[i] - 1];
                ++count;
                ++i;
            } while (i < m && x[draw[i] - 1] == location);

            level[blocks] = total / count;
            left[blocks] = location;
            weight[blocks] = count;
            ++blocks;

            /* Weighted pool-adjacent-violators algorithm. */
            while (blocks > 1 && level[blocks - 2] > level[blocks - 1]) {
                const int a = blocks - 2;
                const int b = blocks - 1;
                const int combined = weight[a] + weight[b];
                level[a] = (weight[a] * level[a] + weight[b] * level[b]) / combined;
                weight[a] = combined;
                --blocks;
            }
        }

        for (int k = 0; k < K; ++k) {
            /* The last fitted block starting at or before grid[k] gives
               the same right-continuous, endpoint-clamped prediction as R. */
            int lo = 0;
            int hi = blocks;
            while (lo < hi) {
                const int mid = lo + (hi - lo) / 2;
                if (left[mid] <= grid[k]) lo = mid + 1;
                else hi = mid;
            }
            const int block = lo > 0 ? lo - 1 : 0;
            answer[k + (R_xlen_t) K * r] = level[block];
        }
    }

    UNPROTECT(1);
    return answer_;
}
