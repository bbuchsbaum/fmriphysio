# Benchmark fast_phy_denoise on simulated data

Times [`fast_phy_denoise()`](fast_phy_denoise.md) over a grid of data
sizes and SVD engines, and scores each run against the ground truth of
[`simulate_phy_data()`](simulate_phy_data.md). Each simulated dataset is
generated once and reused for every engine, so engines are compared on
identical data.

## Usage

``` r
benchmark_fast_phy_denoise(
  grid_n_vox = 4000L,
  grid_n_time = 400L,
  engines = c("svd", "rsvd"),
  reps = 1L,
  tr = 2,
  extractor = c("caa", "dicca", "dipca"),
  out_file = NULL,
  seed = 1L,
  ...
)
```

## Arguments

- grid_n_vox:

  Integer vector of voxel counts.

- grid_n_time:

  Integer vector of time-series lengths.

- engines:

  Character vector of SVD engines to compare (any of `"auto"`, `"svd"`,
  `"rsvd"`).

- reps:

  Number of simulated datasets per size (at least 1).

- tr:

  Repetition time in seconds used to simulate and denoise.

- extractor:

  Component extractor passed to
  [`fast_phy_denoise()`](fast_phy_denoise.md).

- out_file:

  Optional path; results are also written there as CSV.

- seed:

  Base seed. Dataset `i` is simulated with seed `seed + i`, and the same
  seed is passed to [`fast_phy_denoise()`](fast_phy_denoise.md). The
  caller's random-number stream is left untouched. Use `NULL` for
  unseeded runs.

- ...:

  Further arguments to [`fast_phy_denoise()`](fast_phy_denoise.md). `x`,
  `orientation`, `svd_engine`, and `return_diagnostics` are set by the
  benchmark and may not be given.

## Value

A data frame with one row per dataset and engine and columns `n_vox`,
`n_time`, `rep`, `dataset_seed`, `svd_engine`, `status` (`"ok"` or
`"error"`), `total_seconds`, `wnn_seconds`, `pass_total_seconds`,
`selected_components`, `physio_kept_nt` and `neural_kept_nt` (median
over neuronal voxels of the fraction of physiological / neural source
energy remaining after denoising), and `error_message`.

## See also

[`simulate_phy_data()`](simulate_phy_data.md),
[`fast_phy_denoise()`](fast_phy_denoise.md)

## Examples

``` r
res <- benchmark_fast_phy_denoise(
  grid_n_vox = 300, grid_n_time = 100, engines = "svd", reps = 1
)
res[, c("svd_engine", "status", "selected_components", "neural_kept_nt")]
#>   svd_engine status selected_components neural_kept_nt
#> 1        svd     ok                   4      0.9847498
```
