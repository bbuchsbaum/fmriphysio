# One-shot canonical autocorrelation analysis (CAA)

Finds the linear combinations of a set of reduced time series that are
maximally correlated with their own lag-1 values, by an exact lag-1
canonical correlation analysis between \\Z\_{1:T-1}\\ and \\Z\_{2:T}\\.
This is the one-shot replacement for the repeated CAA search in PHYCAA+:
one \\K \times K\\ decomposition in a fixed low-rank space.

## Usage

``` r
extract_caa_oneshot(Z_kt, n_candidates = 30L)
```

## Arguments

- Z_kt:

  Numeric matrix of reduced time series, components (rows) by time
  points (columns); typically the leading right singular vectors of the
  voxel data, transposed. Must have at least 4 time points.

- n_candidates:

  Maximum number of components to return; at most `nrow(Z_kt)` are
  returned.

## Value

A list with elements

- `scores`:

  Time points x components matrix of component time courses, each
  centred and scaled to unit variance.

- `predictability`:

  Canonical lag-1 correlations in \[0, 1\], in decreasing order; the
  magnitude of each component's lag-1 predictability.

- `autocorrelation`:

  Signed lag-1 autocorrelation of each returned time course.

- `method`:

  The string `"caa"`.

## Details

Each component's time course is the left canonical projection applied to
the full (centred) series, so components with negative lag-1
autocorrelation (for example aliased cardiac noise near the Nyquist
frequency) are recovered as well as positively autocorrelated ones.

## References

Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized,
adaptive procedure for measuring and controlling physiological noise in
BOLD fMRI. *NeuroImage*, 82, 306–325.
[doi:10.1016/j.neuroimage.2013.05.102](https://doi.org/10.1016/j.neuroimage.2013.05.102)

## Examples

``` r
set.seed(1)
n <- 200
ar <- as.numeric(stats::arima.sim(list(ar = 0.9), n))
Z <- rbind(ar, matrix(rnorm(4 * n), 4, n))
caa <- extract_caa_oneshot(Z, n_candidates = 3)
round(caa$predictability, 2)
#> [1] 0.91 0.23 0.13
abs(cor(caa$scores[, 1], ar))
#> [1] 0.9993381
```
