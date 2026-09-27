"""COSI orchestration: calibration, batch splitting, estimation, and dispatch.

The four standalone constructors live in cosi_intervals.py. Returned row
indices are one-based so that results can be compared directly with R.
"""
import math

import numpy as np

from cosi_smallest_B import smallest_B_distribution, binomial_cutoff_distribution
from cosi_calibration import _control, _real_scalar, _draw_law, _seed_scope, cosi_kappa
from cosi_intervals import (
    cosi_median, cosi_symmetry, cosi_normality, cosi_unimodal_star,
    contains_cosi_ci, _score_values, _tie_breakers,
)


# Validate gap as a nonnegative integer number of omitted rows between batches.
def _check_gap(gap):
    if not _real_scalar(gap) or not math.isfinite(gap) or gap < 0 or gap != int(gap):
        raise ValueError("gap must be a nonnegative integer")


# Validate shared random inputs; ordered gap splits require an identity permutation.
def _check_randomness(randomness, n, B, signed_rank, gap=0):
    if randomness is None:
        return
    if not isinstance(randomness, dict) or set(randomness) != {
            "permutation", "cutoff_uniform", "uniforms"}:
        raise ValueError("randomness needs permutation, cutoff_uniform, uniforms")
    p = np.asarray(randomness["permutation"])
    if p.ndim != 1 or p.size != n or not np.array_equal(np.sort(p), np.arange(1,n+1)):
        raise ValueError("randomness['permutation'] must permute 1,...,n")
    if gap > 0 and not np.array_equal(p, np.arange(1,n+1)):
        raise ValueError("gap>0 requires the original row order in randomness['permutation']")
    u = randomness["cutoff_uniform"]
    if not _real_scalar(u) or not math.isfinite(u) or not 0 <= u < 1:
        raise ValueError("randomness['cutoff_uniform'] must lie in [0,1)")
    if signed_rank:
        if randomness["uniforms"] is None:
            raise ValueError("randomness['uniforms'] must supply B tie-breakers")
        _tie_breakers(B,randomness["uniforms"])


# Generate a reusable random plan for R/Python comparisons; positive gaps use
# the original row order. Row identifiers are one-based in both languages.
def make_cosi_randomness(n, B, seed=None, gap=0):
    _check_gap(gap)
    for x in (n,B):
        if not _real_scalar(x) or not math.isfinite(x) or x < 1 or x != int(x):
            raise ValueError("n and B must be positive integers")
    n, B = int(n), int(B)
    with _seed_scope(seed) as rng:
        permutation = np.arange(1,n+1) if gap > 0 else rng.permutation(n)+1
        return dict(permutation=permutation, cutoff_uniform=float(rng.random()),
                    uniforms=_tie_breakers(B,rng=rng))


# Compute balanced row groups and omitted indices, preserving the original
# shuffled-split convention for allocating remainder rows when gap is zero.
def _split_rows(n, batches, gap, permutation):
    n_used = n-(batches-1)*gap
    if n_used < batches:
        raise ValueError(f"Not enough rows: need at least {batches+(batches-1)*gap} "
                         f"for {batches} batches with gap {gap}; have {n}")
    sizes = np.full(batches,n_used//batches,dtype=int)
    sizes[1:1+n_used % batches] += 1
    rows, omitted, at = {}, [], 0
    for j,size in enumerate(sizes):
        rows[str(j)] = permutation[at:at+size].copy()
        at += size
        if j < batches-1 and gap > 0:
            omitted.extend(permutation[at:at+gap])
            at += gap
    return dict(batch_rows=rows,batch_sizes={str(j):int(v) for j,v in enumerate(sizes)},
                gap_rows=np.asarray(omitted,dtype=int),n_used=n_used)


# Store the full result as a dictionary with a concise interval representation.
class CosiResult(dict):
    # Format the interval, including endpoint inclusion and conservative calibration.
    def __repr__(self):
        header = (f"COSI {self['method']} interval; B={self['B']}, "
                  f"batches={self['batches']}, alpha={self['alpha']:g}, gap={self['gap']}")
        law = self.get("kappa_distribution")
        if law is not None and law["engine"] == "concentration":
            header += (f"\nConservative concentration calibration ({law['bound_used']}); "
                       f"rejection bound <= {law['rejection_probability_bound']:.6g}")
        if self["empty"]:
            return header+"\nempty set"
        left = "[" if self["closed"]["lower"] else "("
        right = "]" if self["closed"]["upper"] else ")"
        return header+f"\n{left}{self['CI']['lower']:g}, {self['CI']['upper']:g}{right}"


# Choose B and its cutoff, split data into B (or B+1) balanced batches with
# optional gaps, compute scalar estimates, and call the appropriate constructor.
# Extra keyword-only options retain custom scores, conservative cutoffs,
# resource limits, and shared draws for cross-language comparisons.
def cosi_ci(data, alpha=0.05, gamma=1/3, estimate=None, method="median", gap=0,
            seed=None, *, h=None, randomized=True, control=None, randomness=None):
    with _seed_scope(seed) as rng:
        if not callable(estimate):
            raise TypeError("estimate must be a function")
        if not _real_scalar(alpha) or not math.isfinite(alpha) or not 0 < alpha < 1:
            raise ValueError("alpha must lie strictly between 0 and 1")
        if not _real_scalar(gamma) or not math.isfinite(gamma) or not 0 < gamma < .5:
            raise ValueError("gamma must lie strictly between 0 and 1/2")
        _check_gap(gap)
        gap = int(gap)
        if not isinstance(randomized,(bool,np.bool_)):
            raise ValueError("randomized must be True or False")
        method = {"symmetric":"symmetry","unimodal":"unimodal-star"}.get(method,method)
        if method not in ("median","symmetry","normality","unimodal-star"):
            raise ValueError("Unknown method")
        control = _control(control)
        data = np.asarray(data)
        if data.ndim == 1:
            data = data[:,None]
        if data.ndim != 2 or min(data.shape) == 0:
            raise ValueError("data must be a nonempty vector or matrix")
        n = data.shape[0]

        # 1. Select B and calibrate a cutoff independently of the data values.
        split_calibration = smallest_B_distribution(alpha,gamma)
        B = split_calibration["B"]
        # K is the actual number of batches: B+1 for an anchor, otherwise B.
        K = B+int(method == "unimodal-star")
        batches = K
        if n < batches+(batches-1)*gap:
            raise ValueError(f"Not enough rows: need at least {batches+(batches-1)*gap} "
                             f"for {batches} batches with gap {gap}; have {n}")
        signed_rank = method in ("symmetry","normality")
        _check_randomness(randomness,n,B,signed_rank,gap)
        cutoff_uniform = None if randomness is None else randomness["cutoff_uniform"]
        h_values = None
        sign_scores = False
        if method == "median":
            law = split_calibration
        elif method == "unimodal-star":
            law = binomial_cutoff_distribution(B,alpha)
        else:
            h_values = _score_values(B,method,h)
            sign_scores = bool(np.all(h_values[1:] == 1))
            if not sign_scores and B*(B+1)//2 > control["max_walsh"]:
                raise ValueError("Number of Walsh pairs exceeds max_walsh")
            calibrated = cosi_kappa(h_values[1:],alpha,seed=rng,randomized=randomized,
                                    control=control,u=cutoff_uniform)
            law,cutoff = calibrated["distribution"],calibrated["kappa"]
        if not signed_rank:
            cutoff = _draw_law(law,randomized,"max",rng,cutoff_uniform)

        # 2. Positive gaps separate consecutive batches in the ORIGINAL row order.
        permutation = (np.arange(1,n+1) if gap > 0 else rng.permutation(n)+1
                       if randomness is None else np.asarray(randomness["permutation"],dtype=int))
        layout = _split_rows(n,batches,gap,permutation)
        all_estimates = np.empty(batches)
        for j,rows in enumerate(layout["batch_rows"].values()):
            try:
                value = estimate(data[rows-1,:])
            except Exception as error:
                raise ValueError(f"estimate failed on batch {j+1}: {error}") from error
            if not _real_scalar(value) or not math.isfinite(value):
                raise ValueError(f"estimate must return a finite real scalar; failed on batch {j+1}")
            all_estimates[j] = value
        anchor = float(all_estimates[0]) if method == "unimodal-star" else None
        W = all_estimates[1:] if method == "unimodal-star" else all_estimates
        uniforms = (_tie_breakers(B,None if randomness is None else randomness["uniforms"],rng)
                    if signed_rank and not sign_scores else None)

        # 3. Constructors only receive estimates, cutoff, and method-specific inputs.
        if method == "median":
            interval = cosi_median(W,cutoff)
        elif method == "unimodal-star":
            interval = cosi_unimodal_star(W,anchor,cutoff)
        elif method == "normality":
            interval = cosi_normality(W,cutoff,U=uniforms,control=control)
        else:
            interval = cosi_symmetry(W,cutoff,h=h_values,U=uniforms,control=control)
        interval.update(B=B,K=K,batches=batches,ci_est=W,all_estimates=all_estimates,
                        anchor=anchor,uniforms=uniforms,c=cutoff,cutoff_distribution=law,
                        kappa=cutoff if signed_rank else None,
                        kappa_distribution=law if signed_rank else None,h_values=h_values,
                        N_h=None if h_values is None else math.fsum(h_values),
                        split_calibration=split_calibration,**layout,gap=gap,
                        split_order="ordered" if gap > 0 else "shuffled",
                        seed=None if seed is None else int(seed),alpha=float(alpha),
                        gamma=float(gamma),method=method,randomized=bool(randomized))
        return CosiResult(interval)
