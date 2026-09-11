# Fast neuronal-tissue weighting map (wNN) construction.

#' Neuronal-tissue weights from temporal-difference energy
#'
#' Builds the voxel weighting map used by [fast_phy_denoise()] without an FFT
#' or tissue segmentation. Each voxel's energy is the mean squared first
#' difference of its time series, a cheap high-pass power estimate. Energies
#' are sorted and a line is fitted to the lower 80% of ranks; a voxel's
#' *excess* energy (`delta`) is how far it lies above that line. Voxels with
#' large excess high-frequency energy -- vasculature, cerebrospinal fluid, and
#' tissue edges -- are treated as non-neuronal.
#'
#' The returned weight `wNN` is **1 for likely neuronal voxels and 0 for
#' likely non-neuronal voxels** (the polarity used by PHYCAA+). With
#' `threshold = "percentile"` the weight ramps linearly from 1 at the
#' `nt_q` quantile of `delta` to 0 at the `nn_q` quantile. With
#' `threshold = "mixture"` a two-component Gaussian mixture is fitted to
#' `log(delta)` (scale-free, so the result does not depend on data units) and
#' `wNN` is one minus the posterior probability of the high-energy component;
#' voxels with no excess energy get `wNN = 1`.
#'
#' @param X_nt Numeric matrix, voxels (rows) by time points (columns). Rows
#'   should be finite; constant rows get zero energy.
#' @param nt_q Quantile probability in (0, 1) of `delta` at and below which a
#'   voxel is fully neuronal (`wNN = 1`). Default 0.74.
#' @param nn_q Quantile probability in (0, 1), greater than `nt_q`, of `delta`
#'   at and above which a voxel is fully non-neuronal (`wNN = 0`). Default
#'   0.95.
#' @param threshold `"percentile"` (default) or `"mixture"`; see Details.
#' @return A list with elements
#'   \describe{
#'     \item{`wNN`}{Per-voxel weights in \[0, 1\]; 1 = neuronal.}
#'     \item{`delta`}{Per-voxel excess difference energy (>= 0).}
#'     \item{`energy`}{Per-voxel mean squared first difference.}
#'     \item{`posterior_nn`}{Per-voxel non-neuronal weight, `1 - wNN`.}
#'     \item{`delta_nt`, `delta_nn`}{The `delta` values at quantiles `nt_q`
#'       and `nn_q`.}
#'     \item{`threshold_method`}{The thresholding method used.}
#'   }
#' @references
#' Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized, adaptive
#' procedure for measuring and controlling physiological noise in BOLD fMRI.
#' *NeuroImage*, 82, 306--325. \doi{10.1016/j.neuroimage.2013.05.102}
#' @examples
#' sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
#' w <- compute_wNN_diff_energy(sim$x)
#' # Non-neuronal voxels get lower weights:
#' tapply(w$wNN, sim$nn, mean)
#' @export
compute_wNN_diff_energy <- function(
  X_nt,
  nt_q = 0.74,
  nn_q = 0.95,
  threshold = c("percentile", "mixture")
) {
  threshold <- match.arg(threshold)
  X_nt <- as.matrix(X_nt)
  if (!is.numeric(X_nt) || any(!is.finite(X_nt))) {
    stop("X_nt must be a finite numeric matrix.", call. = FALSE)
  }
  nt_q <- check_number(nt_q, "nt_q", lower = 0, upper = 1)
  nn_q <- check_number(nn_q, "nn_q", lower = 0, upper = 1)
  if (nt_q >= nn_q) {
    stop("`nt_q` must be smaller than `nn_q`.", call. = FALSE)
  }
  if (ncol(X_nt) < 3L) {
    stop("Need at least 3 time points for differenced-energy scoring.", call. = FALSE)
  }
  n <- nrow(X_nt)
  if (n < 3L) {
    stop("Need at least 3 voxels for differenced-energy scoring.", call. = FALSE)
  }

  dX <- X_nt[, -1, drop = FALSE] - X_nt[, -ncol(X_nt), drop = FALSE]
  energy <- rowMeans(dX * dX)

  ord <- order(energy)
  e_sorted <- energy[ord]

  k_fit <- min(n, max(3L, floor(0.8 * n)))
  rank_fit <- seq_len(k_fit)
  fit <- stats::lm.fit(cbind(1, rank_fit), e_sorted[rank_fit])
  line_fit <- cbind(1, seq_len(n)) %*% fit$coefficients
  delta_sorted <- pmax(e_sorted - as.numeric(line_fit), 0)

  delta <- numeric(n)
  delta[ord] <- delta_sorted

  delta_nt <- as.numeric(stats::quantile(delta, nt_q, names = FALSE, type = 8))
  delta_nn <- as.numeric(stats::quantile(delta, nn_q, names = FALSE, type = 8))

  if (threshold == "percentile") {
    if (!is.finite(delta_nt) || !is.finite(delta_nn) || delta_nn <= delta_nt) {
      stop(
        "wNN percentile thresholds coincide: too few voxels have excess ",
        "high-frequency energy. Try threshold = \"mixture\" or supply `wnn`.",
        call. = FALSE
      )
    }
    wNN <- build_wnn_ramp(delta, delta_nt, delta_nn)
  } else {
    pos <- delta > 0
    post_nn <- numeric(n)
    if (sum(pos) >= 20L) {
      mix <- fit_gaussian_mixture_1d(log(delta[pos] / stats::median(delta[pos])))
      post_nn[pos] <- mix$post_high
    } else {
      stop("Too few voxels with excess energy for mixture thresholding.", call. = FALSE)
    }
    wNN <- 1 - post_nn
  }

  wNN <- pmin(pmax(wNN, 0), 1)
  list(
    wNN = wNN,
    delta = delta,
    energy = energy,
    posterior_nn = 1 - wNN,
    delta_nt = delta_nt,
    delta_nn = delta_nn,
    threshold_method = threshold
  )
}

build_wnn_ramp <- function(delta, delta_nt, delta_nn) {
  w <- numeric(length(delta))
  w[delta <= delta_nt] <- 1
  w[delta >= delta_nn] <- 0
  mid <- delta > delta_nt & delta < delta_nn
  if (any(mid)) {
    w[mid] <- 1 - (delta[mid] - delta_nt) / (delta_nn - delta_nt)
  }
  w
}

fit_gaussian_mixture_1d <- function(x, max_iter = 100L, tol = 1e-6) {
  x <- as.numeric(x)
  keep <- is.finite(x)
  x_fit <- x[keep]
  if (length(x_fit) < 20L) {
    stop("Too few finite values for mixture thresholding.", call. = FALSE)
  }

  q <- stats::quantile(x_fit, probs = c(0.25, 0.75), names = FALSE, type = 8)
  mu <- c(q[1], q[2])
  sdv <- rep(stats::sd(x_fit), 2)
  sdv[sdv < 1e-6] <- 1
  pi_k <- c(0.5, 0.5)

  ll_prev <- -Inf
  for (iter in seq_len(max_iter)) {
    dens1 <- pi_k[1] * stats::dnorm(x_fit, mu[1], sdv[1])
    dens2 <- pi_k[2] * stats::dnorm(x_fit, mu[2], sdv[2])
    denom <- pmax(dens1 + dens2, 1e-300)

    g1 <- dens1 / denom
    g2 <- dens2 / denom
    n1 <- sum(g1)
    n2 <- sum(g2)

    mu[1] <- sum(g1 * x_fit) / n1
    mu[2] <- sum(g2 * x_fit) / n2
    sdv[1] <- sqrt(sum(g1 * (x_fit - mu[1])^2) / n1)
    sdv[2] <- sqrt(sum(g2 * (x_fit - mu[2])^2) / n2)
    sdv <- pmax(sdv, 1e-6)
    pi_k <- c(n1, n2) / length(x_fit)

    ll <- sum(log(denom))
    if (abs(ll - ll_prev) < tol) break
    ll_prev <- ll
  }

  high <- which.max(mu)
  d1 <- pi_k[1] * stats::dnorm(x_fit, mu[1], sdv[1])
  d2 <- pi_k[2] * stats::dnorm(x_fit, mu[2], sdv[2])
  post2 <- d2 / pmax(d1 + d2, 1e-300)
  post_high_fit <- if (high == 2L) post2 else (1 - post2)

  post_high <- rep(0, length(x))
  post_high[keep] <- post_high_fit

  list(post_high = post_high, mu = mu, sd = sdv, pi = pi_k)
}
