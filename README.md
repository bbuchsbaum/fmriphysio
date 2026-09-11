# fmriphysio

[Getting started](vignettes/fmriphysio.Rmd) · [Changelog](NEWS.md) · [Issues](https://github.com/bbuchsbaum/fmriphysio/issues)

fmriphysio is an R package for removing physiological noise from BOLD fMRI
without physiological recordings or a task design. It extracts candidate
components with a lag-1 canonical autocorrelation analysis, keeps those
expressed more than twice as strongly in non-neuronal tissue (vessels, CSF)
as in neuronal tissue, and projects them out of every voxel. It follows the
logic of PHYCAA+ (Churchill & Strother, 2013), with one-shot linear algebra
in place of its repeated searches. Use it on resting-state or task data that
lack usable cardiac and respiratory recordings.

> **Status:** early development (0.0.0.9000). The API may change. The
> evidence so far comes from simulation (see [Limitations](#limitations)).

## Installation

From r-universe:

```r
install.packages("fmriphysio",
  repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))
```

From GitHub (needs a C++ toolchain):

```r
# install.packages("pak")
pak::pak("bbuchsbaum/fmriphysio")
# or: remotes::install_github("bbuchsbaum/fmriphysio")
```

## Quick start

`simulate_phy_data()` generates data with known cardiac and respiratory
sources, so the result can be checked:

```r
library(fmriphysio)

sim <- simulate_phy_data(n_vox = 1500, n_time = 300, tr = 2, seed = 1)
fit <- fast_phy_denoise(sim$x, tr = sim$tr)
fit
#> <fast_phy_denoise_result>
#>   data:       1500 x 300 (1500 of 1500 voxels analysed, 300 time points)
#>   wNN:        187 non-neuronal (< 0.5), 1313 neuronal (> 0.5) voxels
#>   extractor:  caa, rank 20
#>   removed:    2 of 20 candidate components

# Each removed regressor matches one true physiological source
round(abs(cor(fit$regressors, sim$physio)), 2)
#>            cardiac respiratory
#> pass1_caa1    0.02        1.00
#> pass1_caa2    1.00        0.03
```

`fit$x_clean` is the denoised data, with the dimensions of the input matrix
and voxel means preserved. `fit$component_table` shows why each candidate was
kept or not.

## What's in the box

- `fast_phy_denoise()`: design-free denoising. It accepts voxels × time or
  time × voxels matrices, 4D arrays, and `neuroim2::NeuroVec` images (returned
  as the same kind of image in the same `NeuroSpace`), and takes an optional
  mask. An optional design guardrail leaves task estimates unchanged.
- `compcor_denoise()`: aCompCor and tCompCor baselines for comparison.
- `compute_wNN_diff_energy()` and `extract_caa_oneshot()`: the
  neuronal-tissue map and the lag-1 canonical autocorrelation extractor, as
  standalone functions. DiCCA/DiPCA extractors are available through the
  suggested 'dipca' package.
- `compute_design_free_qc()` and `write_qc_artifact()`: tSNR and variance
  summaries by tissue class, saved as `.rds`, `.csv` or `.json`.
- `simulate_phy_data()` and `benchmark_fast_phy_denoise()`: simulated data
  with ground truth, and timing and accuracy benchmarks.

## Limitations

- The defaults (weighting toward non-neuronal tissue, rank 20, ratio
  threshold 2) depart from PHYCAA+. They were chosen from simulations, not
  validated on acquired data.
- Selection assumes the data contain physiological noise. Without it, neural
  components can be selected and removed. Don't apply the method to data
  that have already been physiologically corrected, and inspect
  `component_table` and `regressors`.
- A physiological source can be only partly removed.
- The design guardrail leaves task-correlated noise in place.
- `compcor_denoise()` is a baseline, not a reproduction of fMRIPrep's
  confounds.

The [getting-started vignette](vignettes/fmriphysio.Rmd) works through the
first four on simulated data.

## References

Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized, adaptive
procedure for measuring and controlling physiological noise in BOLD fMRI.
*NeuroImage*, 82, 306–325. <https://doi.org/10.1016/j.neuroimage.2013.05.102>

Churchill, N. W., Yourganov, G., Spring, R., Rasmussen, P. M., Lee, W.,
Ween, J. E., & Strother, S. C. (2012). PHYCAA: Data-driven measurement and
removal of physiological noise in BOLD fMRI. *NeuroImage*, 59(2), 1299–1314.
<https://doi.org/10.1016/j.neuroimage.2011.08.021>

Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A component based
noise correction method (CompCor) for BOLD and perfusion based fMRI.
*NeuroImage*, 37(1), 90–101. <https://doi.org/10.1016/j.neuroimage.2007.04.042>

Dong, Y., & Qin, S. J. (2018). A novel dynamic PCA algorithm for dynamic data
modeling and process monitoring. *Journal of Process Control*, 67, 1–11.
<https://doi.org/10.1016/j.jprocont.2017.05.002>

## License

MIT © Bradley Buchsbaum
