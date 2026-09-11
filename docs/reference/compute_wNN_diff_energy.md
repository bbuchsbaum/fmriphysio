# Neuronal-tissue weights from temporal-difference energy

Builds the voxel weighting map used by
[`fast_phy_denoise()`](fast_phy_denoise.md) without an FFT or tissue
segmentation. Each voxel's energy is the mean squared first difference
of its time series, a cheap high-pass power estimate. Energies are
sorted and a line is fitted to the lower 80% of ranks; a voxel's
*excess* energy (`delta`) is how far it lies above that line. Voxels
with large excess high-frequency energy – vasculature, cerebrospinal
fluid, and tissue edges – are treated as non-neuronal.

## Usage

``` r
compute_wNN_diff_energy(
  X_nt,
  nt_q = 0.74,
  nn_q = 0.95,
  threshold = c("percentile", "mixture")
)
```

## Arguments

- X_nt:

  Numeric matrix, voxels (rows) by time points (columns). Rows should be
  finite; constant rows get zero energy.

- nt_q:

  Quantile probability in (0, 1) of `delta` at and below which a voxel
  is fully neuronal (`wNN = 1`). Default 0.74.

- nn_q:

  Quantile probability in (0, 1), greater than `nt_q`, of `delta` at and
  above which a voxel is fully non-neuronal (`wNN = 0`). Default 0.95.

- threshold:

  `"percentile"` (default) or `"mixture"`; see Details.

## Value

A list with elements

- `wNN`:

  Per-voxel weights in \[0, 1\]; 1 = neuronal.

- `delta`:

  Per-voxel excess difference energy (\>= 0).

- `energy`:

  Per-voxel mean squared first difference.

- `posterior_nn`:

  Per-voxel non-neuronal weight, `1 - wNN`.

- `delta_nt`, `delta_nn`:

  The `delta` values at quantiles `nt_q` and `nn_q`.

- `threshold_method`:

  The thresholding method used.

## Details

The returned weight `wNN` is **1 for likely neuronal voxels and 0 for
likely non-neuronal voxels** (the polarity used by PHYCAA+). With
`threshold = "percentile"` the weight ramps linearly from 1 at the
`nt_q` quantile of `delta` to 0 at the `nn_q` quantile. With
`threshold = "mixture"` a two-component Gaussian mixture is fitted to
`log(delta)` (scale-free, so the result does not depend on data units)
and `wNN` is one minus the posterior probability of the high-energy
component; voxels with no excess energy get `wNN = 1`.

## References

Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized,
adaptive procedure for measuring and controlling physiological noise in
BOLD fMRI. *NeuroImage*, 82, 306–325.
[doi:10.1016/j.neuroimage.2013.05.102](https://doi.org/10.1016/j.neuroimage.2013.05.102)

## Examples

``` r
sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
w <- compute_wNN_diff_energy(sim$x)
# Non-neuronal voxels get lower weights:
tapply(w$wNN, sim$nn, mean)
#>     FALSE      TRUE 
#> 0.9983603 0.3569858 
```
