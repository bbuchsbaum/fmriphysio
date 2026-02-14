# Fast non-neuronal weighting map (wNN) construction.

#' Compute fast wNN map using differenced energy
#'
#' @param X_nt Numeric matrix (voxels x time).
#' @param nt_q Quantile for neuronal cutoff (default 0.74).
#' @param nn_q Quantile for non-neuronal cutoff (default 0.95).
#' @param threshold One of "percentile" or "mixture".
#' @return A list containing wNN, delta, energy, and threshold diagnostics.
#' @export
compute_wNN_diff_energy <- function(
  X_nt,
  nt_q = 0.74,
  nn_q = 0.95,
  threshold = c("percentile", "mixture")
) {
  threshold <- match.arg(threshold)
  X_nt <- as.matrix(X_nt)

  if (ncol(X_nt) < 3L) {
    stop("Need at least 3 time points for differenced-energy scoring.")
  }

  dX <- X_nt[, -1, drop = FALSE] - X_nt[, -ncol(X_nt), drop = FALSE]
  energy <- rowMeans(dX * dX)

  n <- length(energy)
  ord <- order(energy)
  e_sorted <- energy[ord]

  k_fit <- max(10L, floor(0.8 * n))
  rank_all <- seq_len(n)
  rank_fit <- seq_len(k_fit)

  fit <- stats::lm.fit(cbind(1, rank_fit), e_sorted[rank_fit])
  line_fit <- cbind(1, rank_all) %*% fit$coefficients
  delta_sorted <- pmax(e_sorted - as.numeric(line_fit), 0)

  delta <- numeric(n)
  delta[ord] <- delta_sorted

  if (threshold == "percentile") {
    delta_nt <- as.numeric(stats::quantile(delta, nt_q, na.rm = TRUE, names = FALSE, type = 8))
    delta_nn <- as.numeric(stats::quantile(delta, nn_q, na.rm = TRUE, names = FALSE, type = 8))

    if (!is.finite(delta_nt) || !is.finite(delta_nn) || delta_nn <= delta_nt) {
      stop("Invalid percentile thresholds for wNN.")
    }

    wNN <- build_wnn_ramp(delta, delta_nt, delta_nn)
    post_nn <- 1 - wNN
  } else {
    mix <- fit_gaussian_mixture_1d(log1p(delta))
    post_nn <- mix$post_high
    wNN <- 1 - post_nn

    delta_nt <- as.numeric(stats::quantile(delta, nt_q, na.rm = TRUE, names = FALSE, type = 8))
    delta_nn <- as.numeric(stats::quantile(delta, nn_q, na.rm = TRUE, names = FALSE, type = 8))
  }

  list(
    wNN = pmin(pmax(wNN, 0), 1),
    delta = delta,
    energy = energy,
    posterior_nn = post_nn,
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
    stop("Too few finite values for mixture thresholding.")
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
