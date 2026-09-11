# Simulated BOLD data with known physiological and neural ground truth.

#' Simulate BOLD fMRI data with known physiological noise
#'
#' Generates a voxels x time matrix in which a known subset of voxels is
#' non-neuronal (vasculature / CSF-like) and carries strong physiological
#' noise, while the remaining (neuronal) voxels carry slow neural signal plus
#' weaker physiological noise. Because every source and loading is returned,
#' the data can be used to check how much physiological signal a denoiser
#' removes and how much neural signal it keeps.
#'
#' Two physiological sources are simulated in physical units and then sampled
#' at the repetition time `tr`, so they alias exactly as real signals would:
#'
#' - a cardiac source whose rate drifts slowly around a base rate drawn from
#'   0.9--1.3 Hz (about 55--80 beats per minute), and
#' - a respiratory source whose rate drifts around a base rate drawn from
#'   0.22--0.32 Hz, with slow amplitude modulation.
#'
#' Neural sources are smooth, low-pass (below about 0.08 Hz) random
#' processes. Each voxel receives white noise whose standard deviation varies
#' mildly across voxels, and a positive baseline.
#'
#' @param n_vox Number of voxels (at least 20).
#' @param n_time Number of time points (at least 20).
#' @param tr Repetition time in seconds.
#' @param nn_frac Fraction of voxels that are non-neuronal, in (0, 1).
#' @param baseline Mean signal level. Each voxel's mean is `baseline` scaled
#'   by a random factor near 1; use `0` for zero-mean data.
#' @param physio_amp_nn,physio_amp_nt Typical absolute loading of the
#'   physiological sources in non-neuronal and neuronal voxels, in units of
#'   `noise_sd`, so source-to-noise ratios do not change with `noise_sd`.
#' @param neural_amp_nt,neural_amp_nn Typical absolute loading of the neural
#'   sources in neuronal and non-neuronal voxels, in units of `noise_sd`.
#' @param n_neural Number of neural sources.
#' @param noise_sd Median white-noise standard deviation.
#' @param seed Optional integer seed. The caller's random-number stream is
#'   left untouched when a seed is given.
#' @return A list of class `"fmriphysio_sim"` with elements:
#' \describe{
#'   \item{`x`}{Numeric matrix, `n_vox` x `n_time`: the simulated data.}
#'   \item{`nn`}{Logical vector of length `n_vox`; `TRUE` for non-neuronal
#'     voxels.}
#'   \item{`physio`}{Matrix, `n_time` x 2: standardized cardiac and
#'     respiratory source time courses.}
#'   \item{`neural`}{Matrix, `n_time` x `n_neural`: standardized neural source
#'     time courses.}
#'   \item{`physio_loadings`, `neural_loadings`}{Matrices of voxel loadings
#'     (`n_vox` x 2 and `n_vox` x `n_neural`).}
#'   \item{`noise_sd`}{Per-voxel white-noise standard deviation.}
#'   \item{`voxel_mean`}{Per-voxel baseline.}
#'   \item{`physio_hz`}{Named numeric vector of the base cardiac and
#'     respiratory rates in Hz.}
#'   \item{`tr`}{The repetition time.}
#' }
#' The data satisfy `x = voxel_mean + physio_loadings %*% t(physio) +
#' neural_loadings %*% t(neural) + noise`.
#' @examples
#' sim <- simulate_phy_data(n_vox = 300, n_time = 120, tr = 2, seed = 1)
#' dim(sim$x)
#' table(sim$nn)
#' @export
simulate_phy_data <- function(
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
) {
  n_vox <- check_count(n_vox, "n_vox", min = 20L)
  n_time <- check_count(n_time, "n_time", min = 20L)
  tr <- check_positive(tr, "tr")
  nn_frac <- check_number(nn_frac, "nn_frac", lower = 0, upper = 1)
  if (nn_frac <= 0 || nn_frac >= 1) {
    stop("`nn_frac` must be strictly between 0 and 1.", call. = FALSE)
  }
  baseline <- check_number(baseline, "baseline")
  physio_amp_nn <- check_number(physio_amp_nn, "physio_amp_nn", lower = 0)
  physio_amp_nt <- check_number(physio_amp_nt, "physio_amp_nt", lower = 0)
  neural_amp_nt <- check_number(neural_amp_nt, "neural_amp_nt", lower = 0)
  neural_amp_nn <- check_number(neural_amp_nn, "neural_amp_nn", lower = 0)
  n_neural <- check_count(n_neural, "n_neural", min = 1L)
  noise_sd <- check_positive(noise_sd, "noise_sd")
  seed <- check_number(seed, "seed", null_ok = TRUE)

  with_seed(seed, {
    t_sec <- (seq_len(n_time) - 1) * tr
    n_nn <- max(1L, min(n_vox - 1L, round(nn_frac * n_vox)))
    nn <- rep(FALSE, n_vox)
    nn[sample.int(n_vox, n_nn)] <- TRUE

    # Physiological sources: instantaneous rate (Hz) drifts slowly; the phase
    # is accumulated in continuous time and sampled every `tr` seconds, so
    # the sampled series aliases exactly as a real physiological signal does.
    slow_drift <- function(scale) {
      k <- 1:3
      amp <- stats::rnorm(3) / k
      ph <- stats::runif(3, 0, 2 * pi)
      span <- max(t_sec[n_time], tr)
      d <- colSums(amp * sin(outer(k, t_sec) * 2 * pi / (2 * span) + ph))
      scale * d / max(abs(d), 1e-12)
    }
    cardiac_hz <- stats::runif(1, 0.9, 1.3)
    resp_hz <- stats::runif(1, 0.22, 0.32)
    cardiac_rate <- cardiac_hz * (1 + slow_drift(0.05))
    resp_rate <- resp_hz * (1 + slow_drift(0.08))
    cardiac <- sin(2 * pi * cumsum(cardiac_rate) * tr + stats::runif(1, 0, 2 * pi))
    resp <- (1 + slow_drift(0.3)) *
      sin(2 * pi * cumsum(resp_rate) * tr + stats::runif(1, 0, 2 * pi))
    physio <- cbind(cardiac = cardiac, respiratory = resp)
    physio <- apply(physio, 2, safe_scale_vec)

    # Neural sources: white noise low-pass filtered below ~0.08 Hz with a
    # smooth cosine-taper in the frequency domain.
    cutoff_hz <- min(0.08, 0.4 / tr)
    freqs <- (seq_len(n_time) - 1) / (n_time * tr)
    freqs <- pmin(freqs, 1 / tr - freqs)
    taper <- ifelse(freqs <= cutoff_hz, 1, exp(-((freqs - cutoff_hz) / (0.25 * cutoff_hz))^2))
    neural <- vapply(seq_len(n_neural), function(j) {
      z <- Re(stats::fft(stats::fft(stats::rnorm(n_time)) * taper, inverse = TRUE))
      safe_scale_vec(z)
    }, numeric(n_time))
    neural <- matrix(neural, nrow = n_time, ncol = n_neural)
    colnames(neural) <- paste0("neural_", seq_len(n_neural))

    # Loadings: random sign, magnitude around the requested amplitude, which
    # is expressed in units of the noise standard deviation.
    loading <- function(n, amp, k) {
      matrix(sample(c(-1, 1), n * k, replace = TRUE) * amp * noise_sd * stats::runif(n * k, 0.6, 1.4), n, k)
    }
    Wp <- matrix(0, n_vox, 2)
    Wn <- matrix(0, n_vox, n_neural)
    Wp[nn, ] <- loading(n_nn, physio_amp_nn, 2)
    Wp[!nn, ] <- loading(n_vox - n_nn, physio_amp_nt, 2)
    Wn[!nn, ] <- loading(n_vox - n_nn, neural_amp_nt, n_neural)
    Wn[nn, ] <- loading(n_nn, neural_amp_nn, n_neural)
    colnames(Wp) <- colnames(physio)
    colnames(Wn) <- colnames(neural)

    sd_vox <- noise_sd * exp(stats::rnorm(n_vox, sd = 0.1))
    voxel_mean <- baseline * exp(stats::rnorm(n_vox, sd = 0.05))
    noise <- matrix(stats::rnorm(n_vox * n_time), n_vox, n_time) * sd_vox

    x <- voxel_mean + Wp %*% t(physio) + Wn %*% t(neural) + noise

    structure(
      list(
        x = x,
        nn = nn,
        physio = physio,
        neural = neural,
        physio_loadings = Wp,
        neural_loadings = Wn,
        noise_sd = sd_vox,
        voxel_mean = voxel_mean,
        physio_hz = c(cardiac = cardiac_hz, respiratory = resp_hz),
        tr = tr
      ),
      class = "fmriphysio_sim"
    )
  })
}

#' @method print fmriphysio_sim
#' @export
print.fmriphysio_sim <- function(x, ...) {
  cat(sprintf(
    "<fmriphysio_sim> %d voxels x %d time points (TR = %g s)\n",
    nrow(x$x), ncol(x$x), x$tr
  ))
  cat(sprintf(
    "  non-neuronal voxels: %d (%.0f%%)\n",
    sum(x$nn), 100 * mean(x$nn)
  ))
  cat(sprintf(
    "  physiological sources: cardiac %.2f Hz, respiratory %.2f Hz\n",
    x$physio_hz[["cardiac"]], x$physio_hz[["respiratory"]]
  ))
  cat(sprintf("  neural sources: %d\n", ncol(x$neural)))
  invisible(x)
}
