# Changelog

## fmriphysio 0.0.0.9000

### Renamed

- The package is renamed from `phynd` to `fmriphysio`
  (<https://github.com/bbuchsbaum/fmriphysio>).

### neuroim2 integration

- neuroim2 is now imported.
  [`fast_phy_denoise()`](../reference/fast_phy_denoise.md) and
  [`compcor_denoise()`](../reference/compcor_denoise.md) accept
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  images and return the denoised data as the same kind of image in the
  same `NeuroSpace`: a `DenseNeuroVec`, or a `SparseNeuroVec` with its
  original mask. Only voxels inside a sparse image’s mask are analysed,
  and
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)
  masks are accepted wherever a mask is.

### Changed defaults in `fast_phy_denoise()`

- New argument `extraction_weight`, default `"non_neuronal"`: the data
  are weighted by `1 - wNN` before the low-rank decomposition. The
  previous behaviour, which is also the PHYCAA+ weighting, is
  `"neuronal"`; `"none"` uses unweighted data. Scoring always uses
  unweighted data.
- `pca_rank` now defaults to 20 (was 100), and `ratio_thresh` to 2 (was
  1). In simulations with known sources, run with the corrected
  algorithm, these settings removed more physiological signal and kept
  more neural signal on average than neuronal weighting with rank 100
  and ratio 1, though not in every run. Use
  `extraction_weight = "neuronal", pca_rank = 100, ratio_thresh = 1` to
  approximate the previous settings.

### API changes

- [`fast_phy_denoise()`](../reference/fast_phy_denoise.md): `regressors`
  replaces `components_selected`. The `component_table` column
  `component` is now `candidate`, and a `stability` column is added.
  `wNN` now has one value per voxel of the input, with `NA` for excluded
  voxels; previously it covered only the masked voxels. The result also
  has `active` (voxels analysed) and `params$rank_used`, and a
  [`print()`](https://rdrr.io/r/base/print.html) method.
- [`fast_phy_denoise()`](../reference/fast_phy_denoise.md): `tr` is
  optional and stored as metadata only.
- [`compcor_denoise()`](../reference/compcor_denoise.md): modes are
  `"acompcor"` and `"tcompcor"` (previously `"acompcorr"` and
  `"tcompcorr"`). `center_rows_first` is removed: voxel means are always
  removed before estimation and restored in the output. The result gains
  `nuisance_voxels` (replacing `diagnostics$selector`), `active`, `wNN`,
  `params` (previously `diagnostics$params`; now including
  `n_comp_used`), and a [`print()`](https://rdrr.io/r/base/print.html)
  method. `diagnostics$n_selected_voxels` is now
  `diagnostics$n_nuisance_voxels`.
- [`benchmark_fast_phy_denoise()`](../reference/benchmark_fast_phy_denoise.md):
  `engines` replaces `svd_engines`. New arguments `tr` and `extractor`;
  new result columns `dataset_seed`, `physio_kept_nt` and
  `neural_kept_nt`.
- [`simulate_phy_data()`](../reference/simulate_phy_data.md) is
  exported. It returns a list of class `"fmriphysio_sim"` with the data
  (`x`), the true non-neuronal voxels (`nn`), the physiological and
  neural source time courses and loadings, and `tr`, instead of a bare
  matrix. Cardiac and respiratory sources are defined in Hz and sampled
  at `tr`, so they alias as real signals do.
- [`fast_phy_denoise()`](../reference/fast_phy_denoise.md),
  [`compcor_denoise()`](../reference/compcor_denoise.md) and
  [`compute_design_free_qc()`](../reference/compute_design_free_qc.md)
  gain an `orientation` argument. A matrix with fewer rows than columns
  is still read as time × voxels by default, now with a message.

### Correctness fixes

- Voxel means are restored in `x_clean`. Previously the denoised data
  were returned mean-centred.
- Voxels that are constant, contain non-finite values, or lie outside
  `mask` are excluded from the analysis and returned unchanged, so
  `x_clean` keeps the dimensions of matrix and array input. Previously
  `mask` subset the output, and zero-variance voxels entered the R²
  medians used for selection.
- The subspace rank is capped at `floor(T / 3)`, so the lag-1 analysis
  can’t fit a subspace almost as large as the number of time points.
- With `max_passes > 1`, the union of all selected regressors is removed
  from the original data in a single orthogonal projection. Previously
  the passes were applied one after another, which is not an orthogonal
  projection when regressors from different passes are correlated.
- [`extract_caa_oneshot()`](../reference/extract_caa_oneshot.md)
  computes an exact lag-1 canonical correlation analysis, with both
  lagged blocks whitened, so `predictability` is a true canonical
  correlation. Component time courses come from the left canonical
  projection of the full series, so components with negative lag-1
  autocorrelation (for example aliased cardiac noise near the Nyquist
  frequency) are no longer cancelled. The result gains
  `autocorrelation`.
- Design guardrail: candidates are orthogonalised against the design
  plus an intercept, and predictability values stay aligned with their
  candidates when the guardrail drops a component. A design with
  non-finite values is an error.
- Split-half stability (`stability = "half"`) is redesigned. Candidates
  are re-extracted independently in each half with the same voxel
  weights, and a candidate is kept only if its spatial R² map is
  reproduced in both halves and the two halves agree with each other.
  The score is reported in `component_table$stability`.
- SVD: previously `svd_engine = "auto"` used randomized SVD whenever the
  requested rank was below the smaller data dimension, so results could
  vary between runs unless `seed` was set. `"auto"` now uses an exact
  eigendecomposition of the Gram matrix when the smaller data dimension
  (usually the number of time points) is at most 2000, and randomized
  SVD otherwise. Singular-vector signs are fixed, so results don’t flip
  sign between BLAS/LAPACK builds. Expect numerical results to change.
- Seeds are scoped. `seed` arguments no longer reset the global random
  number stream, and the caller’s `.Random.seed` is restored.
- Input handling: 3D/4D arrays with time last are accepted and returned
  with their original dimensions, as are array-based images such as
  [`neuroim2::DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.html)
  (returned as a base array). Other S4 image objects, such as
  [`neuroim2::SparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.html),
  are coerced with `methods::as(x, "matrix")` and returned as a voxels ×
  time matrix. Masks may be logical vectors, voxel indices, spatial
  arrays, or mask volumes.
- [`compcor_denoise()`](../reference/compcor_denoise.md): nuisance time
  series are detrended (mean and linear trend, plus the cosine basis
  with `pre_highpass = "dct"`) and scaled to unit variance before the
  PCA, as in Behzadi et al. (2007). tCompCor ranks voxels by their
  variance after detrending. Components that the nuisance data can’t
  support are dropped with a warning. The DCT high-pass cutoff is
  checked against the Nyquist frequency. `nuisance_mask` and `wnn` are
  indexed over all voxels, independently of `mask`.
- [`compute_design_free_qc()`](../reference/compute_design_free_qc.md):
  tSNR after denoising uses the raw voxel mean, so it’s meaningful for
  mean-centred output. Voxels with zero variance are skipped instead of
  producing extreme tSNR values. New fields `nt_variance_reduction_frac`
  and `n_unscored`.
- [`compute_wNN_diff_energy()`](../reference/compute_wNN_diff_energy.md):
  the mixture threshold is fitted to `log(delta)` and no longer depends
  on data units. Quantile arguments are validated.
- [`write_qc_artifact()`](../reference/write_qc_artifact.md): the CSV
  includes nested entries such as `ratio_summary` (flattened), and JSON
  writes missing values as `null` instead of the string `"NA"`.
- [`benchmark_fast_phy_denoise()`](../reference/benchmark_fast_phy_denoise.md):
  each simulated dataset is generated once and reused for every SVD
  engine, so engines are compared on identical data, and runs are scored
  against the simulation’s ground truth.
