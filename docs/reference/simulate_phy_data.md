# Simulate BOLD fMRI data with known physiological noise

Generates a voxels x time matrix in which a known subset of voxels is
non-neuronal (vasculature / CSF-like) and carries strong physiological
noise, while the remaining (neuronal) voxels carry slow neural signal
plus weaker physiological noise. Because every source and loading is
returned, the data can be used to check how much physiological signal a
denoiser removes and how much neural signal it keeps.

## Usage

``` r
simulate_phy_data(
  n_vox,
  n_time,
  tr = 2,
  nn_frac = 0.2,
  baseline = 100,
  physio_amp_nn = 1.5,
  physio_amp_nt = 0.3,
  neural_amp_nt = 0.8,
  neural_amp_nn = 0.2,
  n_neural = 3L,
  noise_sd = 1,
  seed = NULL
)
```

## Arguments

- n_vox:

  Number of voxels (at least 20).

- n_time:

  Number of time points (at least 20).

- tr:

  Repetition time in seconds.

- nn_frac:

  Fraction of voxels that are non-neuronal, in (0, 1).

- baseline:

  Mean signal level. Each voxel's mean is `baseline` scaled by a random
  factor near 1; use `0` for zero-mean data.

- physio_amp_nn, physio_amp_nt:

  Typical absolute loading of the physiological sources in non-neuronal
  and neuronal voxels, in units of the noise standard deviation.

- neural_amp_nt, neural_amp_nn:

  Typical absolute loading of the neural sources in neuronal and
  non-neuronal voxels.

- n_neural:

  Number of neural sources.

- noise_sd:

  Median white-noise standard deviation.

- seed:

  Optional integer seed. The caller's random-number stream is left
  untouched when a seed is given.

## Value

A list of class `"fmriphysio_sim"` with elements:

- `x`:

  Numeric matrix, `n_vox` x `n_time`: the simulated data.

- `nn`:

  Logical vector of length `n_vox`; `TRUE` for non-neuronal voxels.

- `physio`:

  Matrix, `n_time` x 2: standardized cardiac and respiratory source time
  courses.

- `neural`:

  Matrix, `n_time` x `n_neural`: standardized neural source time
  courses.

- `physio_loadings`, `neural_loadings`:

  Matrices of voxel loadings (`n_vox` x 2 and `n_vox` x `n_neural`).

- `noise_sd`:

  Per-voxel white-noise standard deviation.

- `voxel_mean`:

  Per-voxel baseline.

- `physio_hz`:

  Named numeric vector of the base cardiac and respiratory rates in Hz.

- `tr`:

  The repetition time.

The data satisfy
`x = voxel_mean + physio_loadings %*% t(physio) + neural_loadings %*% t(neural) + noise`.

## Details

Two physiological sources are simulated in physical units and then
sampled at the repetition time `tr`, so they alias exactly as real
signals would:

- a cardiac source whose rate drifts slowly around a base rate drawn
  from 0.9–1.3 Hz (about 55–80 beats per minute), and

- a respiratory source whose rate drifts around a base rate drawn from
  0.22–0.32 Hz, with slow amplitude modulation.

Neural sources are smooth, low-pass (below about 0.08 Hz) random
processes. Each voxel receives white noise whose standard deviation
varies mildly across voxels, and a positive baseline.

## Examples

``` r
sim <- simulate_phy_data(n_vox = 300, n_time = 120, tr = 2, seed = 1)
dim(sim$x)
#> [1] 300 120
table(sim$nn)
#> 
#> FALSE  TRUE 
#>   240    60 
```
