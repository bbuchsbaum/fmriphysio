# Fast Design-Free PHYCAA+-Style Denoising: Implementation Notes

Date: 2026-02-14  
Scope: Notes for implementing a fast R package denoiser using `neuroim2`
for IO/data structures, with optional dynamic extraction via `dipca`
(`dicca`/`dipca`).

## 1) Goal

Implement a PHYCAA+-class physiological denoiser that is:

- design-free by default,

- optionally design-aware as a guardrail,

- > 10x faster than a literal PHYCAA+ implementation for long runs
  > (large `T`),

- faithful to the key prior: nuisance components are temporally
  structured and non-neuronal (NN)-linked in space.

## 2) Core Simplifications to Preserve Behavior and Cut Runtime

1.  Remove expensive `K`-sweep over PCA sizes + repeated CAA.
2.  Replace per-voxel FFT scoring with time-domain high-frequency
    surrogate.
3.  Use one low-rank decomposition (`K0 << T`) + one-shot lagged
    decomposition.
4.  Compute component spatial scoring via batched correlations (no
    voxelwise regression loops).
5.  Limit outer denoising passes to 1-2.

## 3) Original PHYCAA+ Hotspots (What to Eliminate)

- Step 1 FFT high-frequency power per voxel.
- Step 2 repeated CAA/CCA across many `K` values.
- Repeated per-component variance-map computations inside `K` scan
  loops.
- Repeated full/near-full SVD and iterative loops until convergence.

Primary speed lever: replace repeated search loops with one-shot linear
algebra in a fixed low-rank space.

## 4) Mathematical Shortcut: One-Shot CAA in Whitened PC Space

Let data matrix be `X` (voxels x time), centered over time.  
Low-rank approximation: `X ≈ U Σ V^T`, rank `K0`.

- Score matrix: `Q = Σ V^T`.
- Whitened scores: `Z = Σ^{-1} Q = V^T` (shape `K0 x T`).

Define lagged score blocks:

- `Z0 = Z[, 1:(T-1)]`
- `Z1 = Z[, 2:T]`
- `M = Z0 %*% t(Z1) / (T-1)`

Then CAA/lag-1 CCA reduces to:

- `M = A diag(rho) B^T` (SVD),
- `rho` are canonical lag correlations,
- component candidate scores:
  - `s0_k = t(A[,k]) %*% Z0`
  - `s1_k = t(B[,k]) %*% Z1`
  - align and average `s0_k`, `s1_k` into a length-`T` timecourse `t_k`.

This replaces `K`-sweep CAA with one `K0 x K0` SVD.

## 5) Fast Step 1 (`wNN`) Without FFT or Segmentation Dependency

### 5.1 Time-domain high-frequency proxy

For each voxel `n`:

`e(n) = mean(diff(x_n)^2)`

Rationale: first difference is high-pass weighted, O(`N*T`), fully
vectorizable.

### 5.2 Build deviation-from-linearity `delta`

1.  Sort `e`.
2.  Fit robust line on lower 70-80% ranks.
3.  `delta(rank) = e_sorted - line(rank)`; map back to voxel order.

### 5.3 Thresholds (`deltaNT`, `deltaNN`)

Default design-free thresholds:

- `deltaNT = P74(delta)` -\> likely neuronal tissue (`wNN = 1`)
- `deltaNN = P95(delta)` -\> likely non-neuronal (`wNN = 0`)

Piecewise linear weights (PHYCAA+ style):

- `wNN = 1` if `delta < deltaNT`
- `wNN = 0` if `delta > deltaNN`
- linear ramp in between

Optional smarter variant:

- fit 2-component mixture on `log(e)` or `delta`,
- set `wNN = 1 - P(NN | delta)` as soft posterior weight.

## 6) Step 2: Fast Dynamic Component Extraction (Design-Free)

### 6.1 Weighted low-rank projection

Use `Xw = diag(wNN) %*% X` (or voxelwise multiply rows by `wNN`).

Compute truncated SVD (prefer randomized/partial):

- `K0 = min(100, floor(T/3))` default,
- optionally `K0 = min(200, floor(T/2))` for long runs.

### 6.2 Extract candidate dynamic regressors

Three pluggable methods:

1.  `method = "caa"` (one-shot lagged SVD above).
2.  `method = "dicca"` using `dipca::dicca(Z_t, s, l, inner=...)`.
3.  `method = "dipca"` using `dipca::dipca(Z_t, s, l, algorithm=...)`.

Important orientation:

- `dipca` API expects `X` as `time x variables`.
- If `Z` is `K0 x T`, pass `t(Z)` to `dicca`/`dipca`.

Useful `dipca` knobs observed in local API:

- `dicca(..., inner = c("classic","ar","arma","arima"))`
- `dipca(..., algorithm = c("I","II"), inner_power, inner_tol)`
- `scores(fit)` gives component scores (time x components)
- `fit$R2` available for component predictability ranking/filtering.

### 6.3 NN-link spatial selection (fast Eq.8-style criterion)

For candidate component `t_k` (length `T`), compute voxelwise univariate
`R^2`:

`R2(n) = corr(Xw[n,], t_k)^2`

Vectorized form:

- `c = Xw %*% t_k` (`N x 1`)
- `R2 = c^2 / (||x_n||^2 * ||t_k||^2)`

Selection statistic:

`RatioNN_NT = median(R2[wNN < 0.5]) / median(R2[wNN > 0.5])`

Keep component if:

- `RatioNN_NT > tau` (`tau = 1.0` or `1.1` default),
- and optional predictability filter (`rho_k` or DiCCA `R2_k` above
  threshold).

### 6.4 Design-free stability safeguard

Split time into 2 halves (or odd/even):

- recompute candidate components per split,
- keep only components reproducible across splits:
  - high temporal correlation after sign alignment, or
  - high correlation between split-specific `R2` spatial maps.

This is the design-free substitute for task-regression protection.

## 7) One-Shot Regression / Projection

Stack selected nuisance components in `Tphy` (`T x m`), QR-factorize:

- `Tphy = Q R` with orthonormal `Q`.

Project from all voxels at once:

- if `X` is `N x T`: `Xclean = X - (X %*% Q) %*% t(Q)`
- if `X` is `T x N`: `Xclean = X - Q %*% (t(Q) %*% X)`

Avoid per-voxel regressions; use BLAS GEMM.

## 8) Optional Design-Aware Guardrail (Not Required)

If design matrix `D` exists:

- either residualize data before extraction: `Xr = (I - P_D) X`,
- or orthogonalize nuisance regressors: `Tphy <- (I - P_D) Tphy`.

Use design only as protection, not as dependency.

## 9) Complexity and Why \>10x Is Realistic

Literal PHYCAA+-style cost drivers:

- repeated full/near-full decomposition scales like O(`N*T^2`) in
  practice,
- repeated `K` scan multiplies CAA and spatial scoring cost.

Fast variant:

- truncated randomized SVD: O(`N*T*K0`), `K0 << T`,
- one-shot lagged decomposition: O(`K0^3 + K0^2*T`),
- component scoring: O(`N*T*L`) with small `L`,
- one-shot regression: O(`N*T*m`) with small `m`.

Example (`N=100k`, `T=1200`, `K0=100`):

- `N*T^2` term: `~1.44e11`
- `N*T*K0` term: `~1.2e10`

Already about 12x on decomposition alone before removing `K`-sweep
overhead.

## 10) `neuroim2` Integration Plan

Local API checkpoints (verified from local package files):

- IO: `read_vec()`, `write_vec()`, `read_vol()`, `write_vol()`
- data containers: `NeuroVec`, `DenseNeuroVec`, `SparseNeuroVec`,
  `LogicalNeuroVol`
- memory modes in `read_vec`: `"normal"`, `"mmap"`, `"bigvec"`,
  `"filebacked"`
- matrix conversion available via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) methods for
  sparse/dense vectors.

Recommended execution shape:

1.  Read fMRI 4D as `NeuroVec` via `read_vec()`.
2.  Apply mask as `LogicalNeuroVol` if available to reduce `N`.
3.  Convert to matrix for linear algebra kernels (`N x T` or `T x N`,
    but keep consistent).
4.  Keep metadata/space in `NeuroVec` object for write-back.
5.  Write denoised output with `write_vec()`.

Memory notes:

- prefer `"mmap"` or `"filebacked"` when full dense matrix is too large,
- for very large `N`, compute norms/dot-products in blocks.

## 11) `dipca`/`dicca` Optional Integration Plan

Extractor plugin API:

- `extractor = "caa"` (internal one-shot lagged SVD),
- `extractor = "dicca"` (use
  [`dipca::dicca`](https://rdrr.io/pkg/dipca/man/dicca.html)),
- `extractor = "dipca"` (use
  [`dipca::dipca`](https://rdrr.io/pkg/dipca/man/dipca.html)).

Suggested defaults:

- `dicca`: `s = 1`, `l = 20`, `inner = "classic"`
- `dipca`: `s = 1`, `l = 20`, `algorithm = "II"`

Selection can combine:

- NN-link ratio (`RatioNN_NT`),
- temporal predictability (`rho` for CAA, `fit$R2` for DiCCA/DiPCA),
- split-half stability.

## 12) Proposed Package API (R)

``` r
fast_phy_denoise <- function(
  x,
  tr,
  design = NULL,
  mask = NULL,
  wnn = NULL,
  wnn_method = c("diff_energy", "fft"),
  wnn_threshold = c("percentile", "mixture"),
  delta_nt = 0.74,
  delta_nn = 0.95,
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  pca_rank = 100L,
  n_candidates = 30L,
  ratio_thresh = 1.0,
  pred_thresh = NULL,
  stability = c("none", "half"),
  max_passes = 1L,
  return_diagnostics = TRUE
) {}
```

Suggested return object fields:

- `x_clean`
- `wNN`
- `components_selected`
- `component_scores`
- `ratio_nn_nt`
- `predictability`
- `stability_scores`
- `diagnostics` (timings, ranks, thresholds, pass summaries).

## 13) Implementation Blueprint (Phased)

Phase 0: Baseline scaffolding

- data ingest/output wrappers for `NeuroVec` + matrix path,
- strict orientation helpers (`to_nt()`, `to_tn()`).

Phase 1: Fast `wNN`

- diff-energy score + robust linear deviation,
- percentile and mixture thresholds,
- diagnostics plots (sorted `e`, `delta`, threshold marks).

Phase 2: Extractor engines

- internal `caa_oneshot()`,
- plugin wrappers `extract_dicca()`, `extract_dipca()`,
- unified component table with method-specific metrics.

Phase 3: Selection + regression

- vectorized `R2` map computation,
- `RatioNN_NT` filter,
- optional split-half stability,
- one-shot projection regression.

Phase 4: Performance and quality harness

- timing benchmark suite over synthetic and real runs,
- no-design QC report generation,
- optional task-aware report when design supplied.

## 14) Numerical/Engineering Guardrails

- Always center each voxel time series before dynamic extraction.
- Handle zero-variance voxels safely in `R2` denominator.
- Sign-align split components before temporal correlation.
- Cap selected nuisance count (`m_max`) to avoid over-regression.
- Use deterministic seeds for randomized SVD and initialization.
- Save per-step walltime and FLOP-like counters for optimization
  tracking.

## 15) Validation Metrics (Design-Free First)

Primary no-design QC:

- NT tSNR increase,
- NN variance reduction,
- reduced NN-\>NT leakage (e.g., NN seed correlation in NT),
- split-half reproducibility of FC/ICA features,
- stability of selected nuisance components.

Secondary (when design exists, but not used to drive denoising):

- preserve task contrast sensitivity/specificity metrics.

## 16) Ultra-Fast Fallback Variant

If max speed is required:

1.  Build `wNN`.
2.  Use only NN-heavy voxels (`wNN < 0.5`).
3.  PCA on NN subset -\> top `r` components.
4.  Optional lag-1 autocorrelation filter.
5.  Regress selected components from full data.

This gives CompCor-like behavior with PHYCAA+-style NN prior.

## Details: CompCor Reduction Mode

The package also supports an explicit reduction to CompCor-style
denoising via [`compcor_denoise()`](reference/compcor_denoise.md):

- `mode = "acompcorr"`: select nuisance voxels from `nuisance_mask`
  (preferred), or from `wNN < wnn_thresh` (auto-computed `wNN` if not
  supplied), then extract top PCA nuisance regressors and project them
  out.
- `mode = "tcompcorr"`: select nuisance voxels using highest temporal
  variance (top fraction or fixed count), then extract top PCA nuisance
  regressors and project them out.

This is intended as a clean baseline that is directly comparable to
standard aCompCor/tCompCor behavior while reusing the same data path and
diagnostics.

Compatibility context:

- This CompCor mode is intentionally simplified and “CompCor-style”.
- It should not be interpreted as a strict reproduction of fMRIPrep
  confound generation without additional compatibility logic (for
  example, matching preprocessing, mask construction, and
  component-retention rules).

## 17) Open Decisions to Lock Early

1.  Matrix orientation convention package-wide (`N x T` internally
    recommended).
2.  Initial default extractor (`caa` vs `dicca`).
3.  Default thresholds (`tau`, predictability, stability cutoffs).
4.  Whether mixture-based `wNN` is default from day 1 or v2.
5.  Whether pass count default is fixed at 1 or adaptive up to 2.

## 18) Minimal Pseudocode Skeleton

``` r
# X: N x T
X <- center_rows(X)
wNN <- compute_wNN_diff_energy(X, nt_q = 0.74, nn_q = 0.95)
Xw <- X * wNN

sv <- truncated_svd(Xw, rank = K0)         # U, d, V
Z  <- t(sv$v)                               # K0 x T

cand <- switch(extractor,
  caa   = extract_caa_oneshot(Z),
  dicca = extract_dicca(t(Z), s = lag_order, l = n_candidates),
  dipca = extract_dipca(t(Z), s = lag_order, l = n_candidates)
)

stats <- score_components_nn_nt(Xw, cand$Tcomp, wNN)
keep  <- select_components(stats, ratio_thresh, pred_thresh)
keep  <- enforce_stability_if_needed(Xw, wNN, keep, extractor)

Tphy <- cand$Tcomp[, keep, drop = FALSE]   # T x m
Xc   <- regress_out_projection(X, Tphy)    # N x T
```

------------------------------------------------------------------------

These notes are intended as an implementation contract for a first
high-performance version, then iterative refinement via QC + benchmarks.
