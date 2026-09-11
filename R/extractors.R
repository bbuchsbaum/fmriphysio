# Low-rank and dynamic component extraction helpers.

cpp_symbol_available <- function(symbol) {
  isTRUE(is.loaded(symbol))
}

# Singular values below this fraction of the largest are treated as null
# directions. The Gram-eigen SVD resolves null singular values only to about
# sqrt(machine epsilon) times the largest (~1e-8), so the cut must sit above
# that.
null_sv_tol <- 1e-6

# Truncated SVD returning the leading `rank` singular triplets.
#
# "auto" uses an exact eigendecomposition of the Gram matrix when the smaller
# dimension is at most `gram_max` (fast and deterministic for fMRI, where the
# number of time points is modest), and randomized SVD otherwise. Column signs
# are fixed so that the largest-magnitude entry of each right singular vector
# is positive, making results independent of the LAPACK/BLAS build.
truncated_svd <- function(
  X,
  rank,
  svd_engine = c("auto", "svd", "rsvd"),
  seed = NULL,
  return_u = TRUE,
  gram_max = 2000L
) {
  svd_engine <- match.arg(svd_engine)
  k <- min(as.integer(rank), min(dim(X)))
  if (k < 1L) stop("rank must be >= 1", call. = FALSE)

  engine <- svd_engine
  if (engine == "auto") {
    engine <- if (min(dim(X)) <= gram_max) "gram" else "rsvd"
  }

  out <- switch(
    engine,
    gram = gram_svd(X, k, return_u = return_u),
    svd = {
      s <- svd(X, nu = if (return_u) k else 0L, nv = k)
      list(u = if (return_u) s$u else NULL, d = s$d[seq_len(k)], v = s$v)
    },
    rsvd = {
      s <- with_seed(seed, rsvd::rsvd(X, k = k, nu = if (return_u) k else 0L, nv = k, q = 4L))
      list(u = if (return_u) s$u else NULL, d = s$d[seq_len(k)], v = s$v[, seq_len(k), drop = FALSE])
    }
  )

  flip <- apply(out$v, 2, function(col) {
    s <- sign(col[which.max(abs(col))])
    if (s == 0) 1 else s
  })
  out$v <- sweep(out$v, 2, flip, `*`)
  if (!is.null(out$u)) out$u <- sweep(out$u, 2, flip, `*`)
  out$engine <- engine
  out
}

gram_svd <- function(X, k, return_u = TRUE) {
  scale_d <- function(M, d) sweep(M, 2, ifelse(d > 0, d, 1), `/`)
  if (nrow(X) >= ncol(X)) {
    e <- eigen(crossprod(X), symmetric = TRUE)
    d <- sqrt(pmax(e$values[seq_len(k)], 0))
    v <- e$vectors[, seq_len(k), drop = FALSE]
    u <- if (return_u) scale_d(X %*% v, d) else NULL
  } else {
    e <- eigen(tcrossprod(X), symmetric = TRUE)
    d <- sqrt(pmax(e$values[seq_len(k)], 0))
    u <- e$vectors[, seq_len(k), drop = FALSE]
    v <- scale_d(crossprod(X, u), d)
    if (!return_u) u <- NULL
  }
  list(u = u, d = d, v = v)
}

# Inverse symmetric square root of a covariance matrix, with a small ridge
# relative to its average variance for numerical stability.
inv_sqrt_sym <- function(C, ridge = 1e-8) {
  C <- (C + t(C)) / 2
  lam_ridge <- ridge * max(mean(diag(C)), .Machine$double.eps)
  e <- eigen(C + diag(lam_ridge, nrow(C)), symmetric = TRUE)
  vals <- pmax(e$values, lam_ridge)
  e$vectors %*% (t(e$vectors) / sqrt(vals))
}

#' One-shot canonical autocorrelation analysis (CAA)
#'
#' Finds the linear combinations of a set of reduced time series that are
#' maximally correlated with their own lag-1 values, by an exact lag-1
#' canonical correlation analysis between \eqn{Z_{1:T-1}} and \eqn{Z_{2:T}}.
#' This is the one-shot replacement for the repeated CAA search in PHYCAA+:
#' one \eqn{K \times K} decomposition in a fixed low-rank space.
#'
#' Each component's time course is the left canonical projection applied to
#' the full (centred) series, so components with negative lag-1
#' autocorrelation (for example aliased cardiac noise near the Nyquist
#' frequency) are recovered as well as positively autocorrelated ones.
#'
#' @param Z_kt Numeric matrix of reduced time series, components (rows) by
#'   time points (columns); typically the leading right singular vectors of
#'   the voxel data, transposed. Must have at least 4 time points.
#' @param n_candidates Maximum number of components to return; at most
#'   `nrow(Z_kt)` are returned.
#' @return A list with elements
#'   \describe{
#'     \item{`scores`}{Time points x components matrix of component time
#'       courses, each centred and scaled to unit variance.}
#'     \item{`predictability`}{Canonical lag-1 correlations in \[0, 1\], in
#'       decreasing order; the magnitude of each component's lag-1
#'       predictability.}
#'     \item{`autocorrelation`}{Signed lag-1 autocorrelation of each returned
#'       time course.}
#'     \item{`method`}{The string `"caa"`.}
#'   }
#' @references
#' Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized, adaptive
#' procedure for measuring and controlling physiological noise in BOLD fMRI.
#' *NeuroImage*, 82, 306--325. \doi{10.1016/j.neuroimage.2013.05.102}
#' @examples
#' set.seed(1)
#' n <- 200
#' ar <- as.numeric(stats::arima.sim(list(ar = 0.9), n))
#' Z <- rbind(ar, matrix(rnorm(4 * n), 4, n))
#' caa <- extract_caa_oneshot(Z, n_candidates = 3)
#' round(caa$predictability, 2)
#' abs(cor(caa$scores[, 1], ar))
#' @export
extract_caa_oneshot <- function(Z_kt, n_candidates = 30L) {
  Z_kt <- as.matrix(Z_kt)
  if (!is.numeric(Z_kt) || any(!is.finite(Z_kt))) {
    stop("Z_kt must be a finite numeric matrix.", call. = FALSE)
  }
  n_candidates <- check_count(n_candidates, "n_candidates")
  K <- nrow(Z_kt)
  T <- ncol(Z_kt)
  if (K < 1L) stop("Z_kt must have at least one row.", call. = FALSE)
  if (T < 4L) stop("Need at least 4 time points for lagged extraction.", call. = FALSE)
  if (K > T - 2L) {
    stop(sprintf(
      "Z_kt has %d series but only %d time points; lag-1 canonical analysis needs at most T - 2 series.",
      K, T
    ), call. = FALSE)
  }

  Zc <- Z_kt - rowMeans(Z_kt)
  Z0 <- Zc[, 1:(T - 1), drop = FALSE]
  Z1 <- Zc[, 2:T, drop = FALSE]
  Z0 <- Z0 - rowMeans(Z0)
  Z1 <- Z1 - rowMeans(Z1)

  W0 <- inv_sqrt_sym(tcrossprod(Z0))
  W1 <- inv_sqrt_sym(tcrossprod(Z1))
  decomp <- svd(W0 %*% tcrossprod(Z0, Z1) %*% W1)

  L <- min(n_candidates, K)
  A <- W0 %*% decomp$u[, seq_len(L), drop = FALSE]
  comps <- crossprod(Zc, A)
  comps <- apply(comps, 2, safe_scale_vec)
  comps <- matrix(comps, nrow = T, ncol = L)
  colnames(comps) <- paste0("caa_", seq_len(L))

  autocor <- apply(comps, 2, function(s) sum(s[-1] * s[-T]) / max(sum(s * s), .Machine$double.eps))

  list(
    scores = comps,
    predictability = pmin(decomp$d[seq_len(L)], 1),
    autocorrelation = as.numeric(autocor),
    method = "caa"
  )
}

extract_dynamic_components <- function(
  Z_kt,
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  n_candidates = 30L,
  seed = NULL
) {
  extractor <- match.arg(extractor)
  if (extractor == "caa") {
    return(extract_caa_oneshot(Z_kt, n_candidates = n_candidates))
  }

  if (!requireNamespace("dipca", quietly = TRUE)) {
    stop("extractor = '", extractor, "' requires package 'dipca' (https://bbuchsbaum.r-universe.dev).", call. = FALSE)
  }

  Z_tk <- t(Z_kt)
  l <- min(as.integer(n_candidates), ncol(Z_tk))
  # dipca's fits draw random numbers; scope them so `seed` makes results
  # reproducible and the caller's RNG stream is left untouched.
  fit <- with_seed(seed, if (extractor == "dicca") {
    dipca::dicca(Z_tk, s = as.integer(lag_order), l = l)
  } else {
    dipca::dipca(Z_tk, s = as.integer(lag_order), l = l, algorithm = "II")
  })

  scores <- as.matrix(fit$s)
  pred <- fit$R2
  if (is.null(pred)) pred <- rep(NA_real_, ncol(scores))

  list(scores = scores, predictability = as.numeric(pred), method = extractor, fit = fit)
}

# Squared Pearson correlation between every voxel time series (rows of X_nt)
# and every component time course (columns of T_tl). NA where either series
# has (numerically) zero variance, judged relative to the largest voxel and
# component variance so the result does not depend on data units.
r2_maps <- function(X_nt, T_tl, eps = 1e-12) {
  Xc <- X_nt - rowMeans(X_nt)
  Tc <- sweep(T_tl, 2, colMeans(T_tl), `-`)
  xnorm2 <- rowSums(Xc * Xc)
  tnorm2 <- colSums(Tc * Tc)
  C <- Xc %*% Tc
  R2 <- (C * C) / outer(xnorm2, tnorm2)
  R2[xnorm2 <= eps * max(xnorm2, 0), ] <- NA_real_
  R2[, tnorm2 <= eps * max(tnorm2, 0)] <- NA_real_
  R2
}

score_components_nn_nt <- function(X_nt, T_tl, wNN, eps = 1e-12, use_cpp = TRUE) {
  X_nt <- as.matrix(X_nt)
  T_tl <- as.matrix(T_tl)
  wNN <- as.numeric(wNN)
  L <- ncol(T_tl)

  if (L == 0L) {
    return(list(ratio = numeric(0), med_nn = numeric(0), med_nt = numeric(0), r2_mean = numeric(0)))
  }
  if (!any(wNN < 0.5) || !any(wNN > 0.5)) {
    stop("Need both non-neuronal (wNN < 0.5) and neuronal (wNN > 0.5) voxels for ratio scoring. Check wNN thresholds.", call. = FALSE)
  }

  if (isTRUE(use_cpp) &&
      exists("score_components_nn_nt_cpp", mode = "function", inherits = TRUE) &&
      cpp_symbol_available("_fmriphysio_score_components_nn_nt_cpp")) {
    Tc <- sweep(T_tl, 2, colMeans(T_tl), `-`)
    return(score_components_nn_nt_cpp(X_nt = X_nt - rowMeans(X_nt), T_tl = Tc, wNN = wNN, eps = eps))
  }

  R2 <- r2_maps(X_nt, T_tl, eps = eps)
  nn <- wNN < 0.5
  nt <- wNN > 0.5
  med <- function(v) if (any(is.finite(v))) stats::median(v, na.rm = TRUE) else NA_real_
  med_nn <- apply(R2[nn, , drop = FALSE], 2, med)
  med_nt <- apply(R2[nt, , drop = FALSE], 2, med)
  r2_mean <- apply(R2, 2, function(v) if (any(is.finite(v))) mean(v, na.rm = TRUE) else NA_real_)
  ratio <- ifelse(is.na(med_nn) | is.na(med_nt), NA_real_, med_nn / pmax(med_nt, eps))

  list(ratio = as.numeric(ratio), med_nn = as.numeric(med_nn), med_nt = as.numeric(med_nt), r2_mean = as.numeric(r2_mean))
}

# Split-half reproducibility of candidate components (design notes, 6.4).
#
# Candidates are re-extracted independently in each half of the run, using
# the same per-voxel extraction `weights` as the full-data pass. A
# full-data candidate is stable when its voxelwise R^2 map matches a
# half-1 component and a half-2 component, and those two half components
# match each other -- the half-1/half-2 match is computed on disjoint data,
# so components that do not reproduce across time fail.
component_stability <- function(
  X_nt,
  weights,
  T_tl,
  pca_rank,
  extractor,
  lag_order,
  n_candidates,
  thresh = 0.30,
  svd_engine = "auto",
  seed = NULL
) {
  L <- ncol(T_tl)
  if (L == 0L) {
    return(list(stable = logical(0), score = numeric(0)))
  }
  T <- ncol(X_nt)
  halves <- list(seq_len(floor(T / 2)), seq.int(floor(T / 2) + 1L, T))

  half_maps <- lapply(halves, function(idx) {
    Xh <- center_rows(X_nt[, idx, drop = FALSE])
    k <- min(as.integer(pca_rank), floor(length(idx) / 3), sum(weights > 0))
    if (k < 1L) return(NULL)
    sv <- truncated_svd(Xh * weights, rank = k, svd_engine = svd_engine, seed = seed, return_u = FALSE)
    ok <- sv$d > null_sv_tol * max(sv$d[1], .Machine$double.eps)
    if (!any(ok)) return(NULL)
    ext <- extract_dynamic_components(
      Z_kt = t(sv$v[, ok, drop = FALSE]),
      extractor = extractor,
      lag_order = lag_order,
      n_candidates = n_candidates,
      seed = seed
    )
    r2_maps(Xh, as.matrix(ext$scores))
  })

  if (any(vapply(half_maps, is.null, logical(1)))) {
    return(list(stable = rep(FALSE, L), score = rep(NA_real_, L)))
  }

  map_cor <- function(A, B) {
    A[!is.finite(A)] <- 0
    B[!is.finite(B)] <- 0
    out <- suppressWarnings(stats::cor(A, B))
    out[!is.finite(out)] <- 0
    out
  }

  R_full <- r2_maps(X_nt, T_tl)
  C1 <- map_cor(R_full, half_maps[[1]])
  C2 <- map_cor(R_full, half_maps[[2]])
  C12 <- map_cor(half_maps[[1]], half_maps[[2]])

  score <- vapply(seq_len(L), function(k) {
    j1 <- which.max(C1[k, ])
    j2 <- which.max(C2[k, ])
    min(C1[k, j1], C2[k, j2], C12[j1, j2])
  }, numeric(1))

  list(stable = score >= thresh, score = score)
}
