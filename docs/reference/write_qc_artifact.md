# Write a QC summary to disk

Saves the list returned by
[`compute_design_free_qc()`](compute_design_free_qc.md) as an R object
(`.rds`), a two-column `metric,value` table (`.csv`, nested entries such
as `ratio_summary` flattened to `ratio_summary.min` etc.), or JSON
(`.json`, requires the 'jsonlite' package; missing and non-finite values
are written as `null`).

## Usage

``` r
write_qc_artifact(qc, file)
```

## Arguments

- qc:

  QC list from [`compute_design_free_qc()`](compute_design_free_qc.md).

- file:

  Output path; the format is chosen from its extension (`.rds`, `.csv`,
  or `.json`, case-insensitive).

## Value

`file`, invisibly.

## Examples

``` r
sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
qc <- compute_design_free_qc(sim$x, fit$x_clean, wNN = ifelse(sim$nn, 0, 1))
path <- tempfile(fileext = ".csv")
write_qc_artifact(qc, path)
head(read.csv(path))
#>           metric     value
#> 1          n_vox 300.00000
#> 2         n_time 120.00000
#> 3           n_nt 240.00000
#> 4           n_nn  60.00000
#> 5     n_unscored   0.00000
#> 6 tsnr_nt_before  56.71821
```
