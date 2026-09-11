# CompCor denoising (aCompCor / tCompCor)

Removes the principal temporal components of a nuisance voxel set from
every analysed voxel, following the component-based noise correction
(CompCor) method of Behzadi et al. (2007). It serves as a fixed-rank
baseline against which [`fast_phy_denoise()`](fast_phy_denoise.md) can
be compared.

## Usage

``` r
compcor_denoise(
  x,
  tr,
  mode = c("acompcor", "tcompcor"),
  n_comp = 5L,
  mask = NULL,
  nuisance_mask = NULL,
  wnn = NULL,
  wnn_thresh = 0.5,
  top_var_frac = 0.02,
  top_var_n = NULL,
  design = NULL,
  pre_highpass = c("none", "dct"),
  highpass_hz = 1/128,
  orientation = c("auto", "voxels_by_time", "time_by_voxels"),
  svd_engine = c("auto", "svd", "rsvd"),
  seed = NULL,
  return_diagnostics = TRUE
)
```

## Arguments

- x:

  Numeric data: a voxels x time matrix, a time x voxels matrix (see
  `orientation`), a 4D array with time last, or a
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  image (only voxels inside a `SparseNeuroVec`'s mask are analysed).

- tr:

  Repetition time in seconds.

- mode:

  `"acompcor"` (nuisance voxels from a mask or wNN map) or `"tcompcor"`
  (highest-variance voxels).

- n_comp:

  Number of principal components to remove (at least 1).

- mask:

  Optional analysis mask (logical vector, voxel indices, 3D array, or
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)).
  Only masked voxels are cleaned.

- nuisance_mask:

  Optional nuisance region for `"acompcor"`, in the same forms as
  `mask`, over all voxels of `x`.

- wnn:

  Optional wNN weights for `"acompcor"` (1 = neuronal, 0 =
  non-neuronal), over all voxels of `x` or over the masked voxels.

- wnn_thresh:

  Voxels with wNN below this value (in \[0, 1\]) are nuisance voxels in
  `"acompcor"`.

- top_var_frac:

  Fraction (in (0, 1\]) of analysed voxels used by `"tcompcor"` when
  `top_var_n` is `NULL`.

- top_var_n:

  Optional number of highest-variance voxels for `"tcompcor"`; clamped
  to between `n_comp` and the number of analysed voxels.

- design:

  Optional design matrix (time x predictors). See the design guardrail
  in Details.

- pre_highpass:

  `"none"` or `"dct"`: remove a discrete cosine basis with cutoff
  `highpass_hz` from the analysed voxels before estimating and removing
  components.

- highpass_hz:

  High-pass cutoff in Hz for `pre_highpass = "dct"` (default 1/128 Hz);
  must be below the Nyquist frequency `1 / (2 * tr)`.

- orientation:

  `"auto"` (voxels x time unless there are fewer rows than columns),
  `"voxels_by_time"`, or `"time_by_voxels"`.

- svd_engine:

  `"auto"`, `"svd"`, or `"rsvd"`; see
  [`fast_phy_denoise()`](fast_phy_denoise.md).

- seed:

  Optional integer seed for the randomized SVD. The caller's
  random-number stream is left untouched.

- return_diagnostics:

  Whether to include the `diagnostics` element.

## Value

An object of class `"compcor_denoise_result"`: a list with

- `x_clean`:

  Denoised data of the same type, shape, and orientation as `x` (a
  `NeuroVec` input gives a `DenseNeuroVec` or `SparseNeuroVec` in the
  same space), voxel means retained.

- `regressors`:

  Matrix, time x components: the removed components, each scaled to unit
  variance (zero columns if none remain).

- `wNN`:

  wNN weights over all voxels (`NA` where unavailable) when a wNN map
  was used to select nuisance voxels, otherwise `NULL`.

- `nuisance_voxels`:

  Logical vector over all voxels: the voxels whose time series formed
  the PCA.

- `active`:

  Logical vector over all voxels: the voxels that were cleaned.

- `mode`, `tr`:

  As supplied.

- `params`:

  List of settings, including `n_comp` (requested) and `n_comp_used`.

- `diagnostics`:

  When `return_diagnostics = TRUE`: input information, voxel counts, the
  source of the nuisance selection, singular values, the high-pass
  basis, the design columns dropped by the guardrail, and the wNN
  computation (for automatic wNN).

## Details

**Choosing nuisance voxels.**

- `mode = "acompcor"` (anatomical CompCor) uses the voxels in
  `nuisance_mask` when it is supplied; otherwise voxels whose wNN weight
  is below `wnn_thresh`, using `wnn` when supplied and otherwise a map
  computed with
  [`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md) from the
  analysed voxels.

- `mode = "tcompcor"` (temporal CompCor) uses the `top_var_n` analysed
  voxels with the highest temporal variance after detrending (by default
  the larger of `5 * n_comp` and `top_var_frac` of the analysed voxels).

`nuisance_mask` and `wnn` are indexed over *all* voxels of `x` and are
independent of `mask`, so a white-matter/CSF region outside the analysis
mask can supply the regressors. A `wnn` of length `sum(mask)` is also
accepted and read as covering the masked voxels; `NA` weights never
select a voxel.

**Estimating components.** Each nuisance voxel time series has its mean
and linear trend removed (and the cosine high-pass basis, when
`pre_highpass = "dct"`) and is scaled to unit variance before the PCA,
as in Behzadi et al. (2007). Components whose singular value is below
`1e-6` times the largest are discarded, so no more components are
removed than the nuisance data support; a warning is given when fewer
than `n_comp` remain.

**Cleaning.** The components are projected out of every analysed voxel
in one orthogonal projection, and each voxel's mean is added back.
Voxels outside `mask`, constant voxels, and voxels with non-finite
values are returned unchanged. With `pre_highpass = "dct"` the output is
also high-pass filtered (the cosine basis is removed from the analysed
voxels); the basis is returned in `diagnostics$highpass$basis`.

**Design guardrail.** When `design` is supplied, the regressors are made
orthogonal to the design (plus an intercept) before removal.
Least-squares task estimates computed from `x_clean` are then identical
to those from `x`: the guardrail protects task signal but also leaves
any task-correlated physiological noise in place. With
`pre_highpass = "dct"` the regressors are also kept orthogonal to the
cosine basis, so the output stays high-passed, and the guarantee holds
for the high-passed design.

**Relation to fMRIPrep.** This is a CompCor-style baseline, not a
reimplementation of the fMRIPrep/NiWorkflows confounds pipeline. Mask
construction and erosion, censoring, and the component-retention rule (a
fixed `n_comp` here, a variance-explained threshold in fMRIPrep) differ.
Use fMRIPrep's confounds directly when exact parity is required.

## References

Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A component
based noise correction method (CompCor) for BOLD and perfusion based
fMRI. *NeuroImage*, 37(1), 90–101.
[doi:10.1016/j.neuroimage.2007.04.042](https://doi.org/10.1016/j.neuroimage.2007.04.042)

## See also

[`fast_phy_denoise()`](fast_phy_denoise.md),
[`compute_design_free_qc()`](compute_design_free_qc.md)

## Examples

``` r
sim <- simulate_phy_data(n_vox = 400, n_time = 120, tr = 2, seed = 1)

# aCompCor with a known nuisance region
fit <- compcor_denoise(sim$x, tr = 2, mode = "acompcor",
                       nuisance_mask = sim$nn, n_comp = 3)
fit
#> <compcor_denoise_result> mode: acompcor
#>   voxels cleaned: 400 of 400; nuisance voxels: 80
#>   components removed: 3 (requested 3)

# tCompCor with a high-pass filter
fit_t <- compcor_denoise(sim$x, tr = 2, mode = "tcompcor", n_comp = 3,
                         pre_highpass = "dct")
dim(fit_t$regressors)
#> [1] 120   3
```
