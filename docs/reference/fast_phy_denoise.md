# Fast design-free physiological denoising (PHYCAA+-style)

Removes physiological noise from fMRI time series without a task design,
following the logic of PHYCAA+ with one-shot linear algebra in place of
its repeated searches:

## Usage

``` r
fast_phy_denoise(
  x,
  tr = NULL,
  design = NULL,
  mask = NULL,
  orientation = c("auto", "voxels_by_time", "time_by_voxels"),
  wnn = NULL,
  wnn_method = c("diff_energy"),
  wnn_threshold = c("percentile", "mixture"),
  delta_nt = 0.74,
  delta_nn = 0.95,
  extraction_weight = c("non_neuronal", "neuronal", "none"),
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  pca_rank = 20L,
  n_candidates = 30L,
  ratio_thresh = 2,
  pred_thresh = NULL,
  stability = c("none", "half"),
  stability_thresh = 0.3,
  max_passes = 1L,
  return_diagnostics = TRUE,
  svd_engine = c("auto", "svd", "rsvd"),
  use_cpp = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  Numeric data: a matrix (voxels x time by default; see `orientation`),
  a 3D/4D array with time as the last dimension, or a
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  image (for example a `DenseNeuroVec` or `SparseNeuroVec`).

- tr:

  Optional repetition time in seconds. Not used by the algorithm; stored
  in the result as metadata.

- design:

  Optional task design matrix (time points x predictors) used as a
  guardrail; see Details.

- mask:

  Optional voxel mask: a logical vector, whole-number voxel indices, a
  spatial array, or a
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html).
  Voxels outside the mask are returned unchanged.

- orientation:

  Layout of a matrix `x`: `"voxels_by_time"`, `"time_by_voxels"`, or
  `"auto"` (default), which treats a matrix with fewer rows than columns
  as time x voxels and says so with a message.

- wnn:

  Optional precomputed neuronal-tissue weights in \[0, 1\] (1 =
  neuronal), one per voxel or one per masked voxel. Overrides the
  internal computation.

- wnn_method:

  Method for computing `wNN`; currently only `"diff_energy"`.

- wnn_threshold:

  `"percentile"` or `"mixture"`; passed to
  [`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md) as
  `threshold`.

- delta_nt, delta_nn:

  Quantile probabilities in (0, 1) passed to
  [`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md) as `nt_q`
  and `nn_q`. Used only when `wnn` is not supplied; `delta_nn` is
  ignored by the mixture threshold.

- extraction_weight:

  Per-voxel weighting applied before the low-rank decomposition:
  `"non_neuronal"` (default; weight `1 - wNN`), `"neuronal"` (weight
  `wNN`, as in PHYCAA+), or `"none"`. Scoring always uses the unweighted
  data.

- extractor:

  Component extractor: `"caa"` (built in), or `"dicca"` / `"dipca"`,
  which require the 'dipca' package from
  <https://bbuchsbaum.r-universe.dev>.

- lag_order:

  Lag order for DiCCA/DiPCA (whole number \>= 1). The CAA extractor
  always uses lag 1.

- pca_rank:

  Maximum rank of the low-rank subspace (default 20). The rank used is
  `min(pca_rank, floor(T / 3), number of analysed voxels)`.

- n_candidates:

  Maximum number of candidate components per pass.

- ratio_thresh:

  Minimum ratio of median \\R^2\\ in non-neuronal voxels to median
  \\R^2\\ in neuronal voxels for a candidate to be selected (strictly
  greater). The default 2 selects components expressed more than twice
  as strongly in non-neuronal as in neuronal tissue.

- pred_thresh:

  Optional minimum predictability in \[0, 1\]: the canonical lag-1
  correlation for `"caa"`, or the DiCCA/DiPCA \\R^2\\. `NULL` (default)
  disables the filter.

- stability:

  `"none"` (default) or `"half"`. With `"half"`, candidates are
  re-extracted independently in each half of the run and a candidate is
  kept only if its voxelwise \\R^2\\ map is reproduced in both halves
  (correlation of maps \>= `stability_thresh`).

- stability_thresh:

  Minimum correlation between spatial \\R^2\\ maps for split-half
  stability, in \[0, 1\].

- max_passes:

  Number of extraction/removal passes (whole number \>= 1; 1 or 2
  recommended).

- return_diagnostics:

  Whether to include the `diagnostics` element.

- svd_engine:

  Low-rank decomposition: `"auto"` (exact Gram-matrix eigendecomposition
  when the smaller data dimension – usually the number of time points –
  is at most 2000, randomized SVD otherwise), `"svd"` (exact LAPACK
  SVD), or `"rsvd"` (randomized).

- use_cpp:

  Use the compiled scoring kernel (identical results; faster for large
  data).

- seed:

  Optional seed for the randomized SVD and the DiCCA/DiPCA extractors.
  The caller's random number stream is restored afterwards.

## Value

An object of class `fast_phy_denoise_result`, a list with

- `x_clean`:

  Denoised data of the same type, shape, and orientation as `x` (a
  `NeuroVec` input gives a `DenseNeuroVec` or `SparseNeuroVec` in the
  same space); voxel means are preserved and excluded voxels are
  unchanged.

- `wNN`:

  Neuronal-tissue weight per voxel (1 = neuronal), `NA` for voxels
  excluded from the analysis.

- `regressors`:

  Time points x selected components matrix of the nuisance regressors
  that were projected out (zero columns if none).

- `component_table`:

  Data frame with one row per scored candidate: `pass`, `candidate`,
  `ratio_nn_nt`, `med_nn`, `med_nt`, `predictability`, `stability`
  (split-half score, `NA` unless `stability = "half"`), and `selected`.

- `active`:

  Logical vector, one per voxel: `TRUE` for voxels that were analysed
  and denoised.

- `tr`:

  The supplied `tr`, or `NULL`.

- `params`:

  The main settings used, including the rank actually used.

- `diagnostics`:

  If requested: input and voxel-exclusion summaries, the full `wNN`
  computation, a per-pass log, and timings.

## Details

1.  **Tissue weights.** A neuronal-tissue weight `wNN` is computed per
    voxel from temporal-difference energy
    ([`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md)); 1 marks
    likely neuronal voxels, 0 likely non-neuronal voxels (vasculature,
    cerebrospinal fluid).

2.  **Low-rank dynamics.** The voxel-centred data are weighted per voxel
    – by default by `1 - wNN`, emphasising non-neuronal tissue, where
    physiological noise is strongest (see `extraction_weight`) – and
    reduced to their leading `min(pca_rank, floor(T / 3))` temporal
    singular vectors. Temporally predictable components are extracted
    from that subspace
    ([`extract_caa_oneshot()`](extract_caa_oneshot.md), or DiCCA/DiPCA
    via the 'dipca' package).

3.  **Selection.** Each candidate is scored by the ratio of its median
    voxelwise \\R^2\\ in non-neuronal voxels to that in neuronal voxels,
    computed on the unweighted data. Candidates with ratio above
    `ratio_thresh` (and optionally predictability above `pred_thresh`
    and split-half stability) are nuisance regressors.

4.  **Removal.** The selected regressors are projected out of every
    analysed voxel in one orthogonal projection, and voxel means are
    restored. With `max_passes > 1`, extraction and selection repeat on
    the cleaned data and the union of all selected regressors is
    projected out of the original data.

Voxels outside `mask`, voxels containing non-finite values, and voxels
that are constant over time are excluded from the analysis and returned
unchanged. `x_clean` has the type and shape of `x`: matrices and arrays
keep their dimensions and orientation, and a
[`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
is returned as a `DenseNeuroVec` – or, for sparse input, a
`SparseNeuroVec` with the same mask – in the same `NeuroSpace`. Only
voxels inside a sparse image's mask are analysed.

**Design guardrail.** When `design` is supplied, candidate regressors
are orthogonalised against the design (plus an intercept) before scoring
and removal. Denoising therefore cannot change task-effect estimates
from a GLM with that design; the trade-off is that task-correlated
physiological noise is left in place.

## Departures from PHYCAA+ and limitations

PHYCAA+ estimates components from data weighted toward neuronal tissue
(`extraction_weight = "neuronal"`) and searches over many subspace
ranks. In simulations with known physiological and neural sources,
weighting toward non-neuronal tissue with a small fixed rank
(`pca_rank = 20`) and `ratio_thresh = 2` removed substantially more
physiological signal while preserving more neural signal than the
PHYCAA+ weighting or larger ranks, so these are the defaults. Larger
ranks let the lag-1 canonical analysis fit spurious components.

The defaults are better on average, not in every case. When
physiological rhythms alias to nearby frequencies – to each other, or
into the low-frequency band occupied by neural signal – they can be
mixed within a candidate whose ratio falls below `ratio_thresh`, and
part of that noise is left in the data.

The selection criterion assumes the data contain physiological noise
concentrated in non-neuronal tissue. Without such noise, the
difference-energy weights no longer reflect tissue type (they follow
incidental differences in voxel noise or signal amplitude), and
components carrying neural signal can then be selected and removed.
Inspect `component_table` and the removed `regressors`, and do not apply
the method to data that have already been physiologically corrected.

## References

Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized,
adaptive procedure for measuring and controlling physiological noise in
BOLD fMRI. *NeuroImage*, 82, 306–325.
[doi:10.1016/j.neuroimage.2013.05.102](https://doi.org/10.1016/j.neuroimage.2013.05.102)

Churchill, N. W., Yourganov, G., Spring, R., Rasmussen, P. M., Lee, W.,
Ween, J. E., & Strother, S. C. (2012). PHYCAA: Data-driven measurement
and removal of physiological noise in BOLD fMRI. *NeuroImage*, 59(2),
1299–1314.
[doi:10.1016/j.neuroimage.2011.08.021](https://doi.org/10.1016/j.neuroimage.2011.08.021)

## See also

[`compcor_denoise()`](compcor_denoise.md) for a CompCor baseline and
[`compute_design_free_qc()`](compute_design_free_qc.md) for
quality-control summaries.

## Examples

``` r
sim <- simulate_phy_data(n_vox = 400, n_time = 150, seed = 1)
fit <- fast_phy_denoise(sim$x, tr = sim$tr)
fit
#> <fast_phy_denoise_result>
#>   data:       400 x 150 (400 of 400 voxels analysed, 150 time points)
#>   wNN:        67 non-neuronal (< 0.5), 333 neuronal (> 0.5) voxels
#>   extractor:  caa, rank 20
#>   removed:    11 of 20 candidate components
head(fit$component_table)
#>   pass candidate ratio_nn_nt     med_nn      med_nt predictability stability
#> 1    1         1   0.3745497 0.02737985 0.073100733      0.7561841        NA
#> 2    1         2   1.7049928 0.09577284 0.056171994      0.7018514        NA
#> 3    1         3   5.3797370 0.15563530 0.028929908      0.6596105        NA
#> 4    1         4   4.0807447 0.09306398 0.022805636      0.6236286        NA
#> 5    1         5   0.3199423 0.01223561 0.038243192      0.5745125        NA
#> 6    1         6   7.5529580 0.02629394 0.003481277      0.5396847        NA
#>   selected
#> 1    FALSE
#> 2    FALSE
#> 3     TRUE
#> 4     TRUE
#> 5    FALSE
#> 6     TRUE
```
