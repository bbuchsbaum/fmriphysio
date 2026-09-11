# Compute design-free QC metrics

Summarises how denoising changed the data without reference to any task
design, separately for neuronal (wNN \> 0.5) and non-neuronal (wNN \<
0.5) voxels.

## Usage

``` r
compute_design_free_qc(
  x_raw,
  x_clean,
  wNN,
  component_table = NULL,
  orientation = c("auto", "voxels_by_time", "time_by_voxels")
)
```

## Arguments

- x_raw:

  Raw data, in any form accepted by
  [`fast_phy_denoise()`](fast_phy_denoise.md).

- x_clean:

  Denoised data with the same shape as `x_raw`.

- wNN:

  Voxel weights (1 = neuronal, 0 = non-neuronal), one per voxel, e.g.
  `fit$wNN` from [`fast_phy_denoise()`](fast_phy_denoise.md). Voxels
  with `NA` weight, or with non-finite data, are not scored.

- component_table:

  Optional component table from
  [`fast_phy_denoise()`](fast_phy_denoise.md), summarised in
  `selected_components` and `ratio_summary`.

- orientation:

  Orientation of `x_raw`: `"auto"`, `"voxels_by_time"`, or
  `"time_by_voxels"`. `x_clean` is read in the same orientation.

## Value

A named list: `n_vox`, `n_time`, `n_nt`, `n_nn` (voxels scored in each
class), `n_unscored`; `tsnr_nt_before`, `tsnr_nt_after`,
`tsnr_nt_delta`, `tsnr_nn_before`, `tsnr_nn_after` (median tSNR);
`nt_variance_reduction_frac`, `nn_variance_before`, `nn_variance_after`,
`nn_variance_reduction_frac` (median variance and the fractional
reduction of the median); `selected_components` (`NA` without a
component table); and `ratio_summary`, a list with the `min`, `median`,
`mean` and `max` of the finite NN/NT ratios (`NA` when there are none),
or `NULL` without a component table.

## Details

Temporal SNR is `|mean| / sd` per voxel. Because denoisers may alter
voxel means, the "after" tSNR uses the *raw* voxel mean over the cleaned
standard deviation, so it measures noise reduction on a fixed signal
level. Voxels whose standard deviation is zero (relative to their mean)
have undefined tSNR and are skipped; a metric with no defined voxels is
`NA`.

## See also

[`write_qc_artifact()`](write_qc_artifact.md)

## Examples

``` r
sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
qc <- compute_design_free_qc(sim$x, fit$x_clean, wNN = ifelse(sim$nn, 0, 1))
qc$tsnr_nt_before
#> [1] 56.71821
qc$tsnr_nt_after
#> [1] 58.95709
```
