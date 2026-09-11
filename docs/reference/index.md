# Package index

## Denoising

Remove physiological noise from fMRI time series.

- [`fast_phy_denoise()`](fast_phy_denoise.md) : Fast design-free
  physiological denoising (PHYCAA+-style)
- [`compcor_denoise()`](compcor_denoise.md) : CompCor denoising
  (aCompCor / tCompCor)

## Tissue weights and components

The building blocks of
[`fast_phy_denoise()`](../reference/fast_phy_denoise.md):
neuronal-tissue weights and lag-1 canonical autocorrelation components.

- [`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md) :
  Neuronal-tissue weights from temporal-difference energy
- [`extract_caa_oneshot()`](extract_caa_oneshot.md) : One-shot canonical
  autocorrelation analysis (CAA)

## Quality control and benchmarking

Summarise what denoising changed, save the summary, and time the
pipeline against simulated ground truth.

- [`compute_design_free_qc()`](compute_design_free_qc.md) : Compute
  design-free QC metrics
- [`write_qc_artifact()`](write_qc_artifact.md) : Write a QC summary to
  disk
- [`benchmark_fast_phy_denoise()`](benchmark_fast_phy_denoise.md) :
  Benchmark fast_phy_denoise on simulated data

## Simulation

Simulated BOLD data with known physiological and neural sources.

- [`simulate_phy_data()`](simulate_phy_data.md) : Simulate BOLD fMRI
  data with known physiological noise
