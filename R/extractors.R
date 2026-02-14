# Low-rank and dynamic component extraction helpers.

cpp_symbol_available <- function(symbol) {
  isTRUE(is.loaded(symbol))
}

truncated_svd <- function(X, rank, svd_engine = c("auto", "rsvd", "svd"), seed = NULL) {
  svd_engine <- match.arg(svd_engine)
  k <- min(as.integer(rank), min(dim(X)))
  if (k < 1L) stop("rank must be >= 1")

  if (!is.null(seed)) set.seed(seed)

  use_rsvd <- identical(svd_engine, "rsvd") ||
    (identical(svd_engine, "auto") && k < min(dim(X)))

  if (use_rsvd) {
    out <- rsvd::rsvd(X, k = k, nu = k, nv = k)
    return(list(u = out$u, d = out$d, v = out$v, engine = "rsvd"))
  }

  out <- svd(X, nu = k, nv = k)
  list(u = out$u, d = out$d[seq_len(k)], v = out$v, engine = "svd")
}

#' One-shot CAA extraction in whitened PC subspace
#'
#' @param Z_kt Numeric matrix with shape K x T (whitened reduced scores).
#' @param n_candidates Number of candidates to return.
#' @return A list with `scores` (T x L) and `predictability` (length L).
#' @export
extract_caa_oneshot <- function(Z_kt, n_candidates = 30L) {
  Z_kt <- as.matrix(Z_kt)
  K <- nrow(Z_kt)
  T <- ncol(Z_kt)
  if (T < 4L) stop("Need at least 4 time points for lagged extraction.")

  Z0 <- Z_kt[, 1:(T - 1), drop = FALSE]
  Z1 <- Z_kt[, 2:T, drop = FALSE]
  M <- Z0 %*% t(Z1) / (T - 1)

  decomp <- svd(M)
  L <- min(as.integer(n_candidates), K, length(decomp$d))
  comps <- matrix(0, nrow = T, ncol = L)
  rho <- decomp$d[seq_len(L)]

  for (k in seq_len(L)) {
    s0 <- as.numeric(crossprod(decomp$u[, k], Z0))
    s1 <- as.numeric(crossprod(decomp$v[, k], Z1))

    tvec <- numeric(T)
    cnt <- numeric(T)
    tvec[1:(T - 1)] <- tvec[1:(T - 1)] + s0
    cnt[1:(T - 1)] <- cnt[1:(T - 1)] + 1
    tvec[2:T] <- tvec[2:T] + s1
    cnt[2:T] <- cnt[2:T] + 1
    tvec <- tvec / pmax(cnt, 1)
    comps[, k] <- safe_scale_vec(tvec)
  }

  colnames(comps) <- paste0("caa_", seq_len(L))
  list(scores = comps, predictability = rho, method = "caa")
}

extract_dynamic_components <- function(
  Z_kt,
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  n_candidates = 30L
) {
  extractor <- match.arg(extractor)
  if (extractor == "caa") {
    return(extract_caa_oneshot(Z_kt, n_candidates = n_candidates))
  }

  if (!requireNamespace("dipca", quietly = TRUE)) {
    stop("extractor='dicca' or 'dipca' requires package 'dipca'.")
  }

  Z_tk <- t(Z_kt)
  if (extractor == "dicca") {
    fit <- dipca::dicca(Z_tk, s = as.integer(lag_order), l = as.integer(n_candidates))
  } else {
    fit <- dipca::dipca(
      Z_tk,
      s = as.integer(lag_order),
      l = as.integer(n_candidates),
      algorithm = "II"
    )
  }

  scores <- as.matrix(fit$s)
  pred <- fit$R2
  if (is.null(pred)) pred <- rep(NA_real_, ncol(scores))

  list(scores = scores, predictability = as.numeric(pred), method = extractor, fit = fit)
}

score_components_nn_nt <- function(X_nt, T_tl, wNN, eps = 1e-12, use_cpp = TRUE) {
  X_nt <- as.matrix(X_nt)
  T_tl <- as.matrix(T_tl)

  if (isTRUE(use_cpp) &&
      exists("score_components_nn_nt_cpp", mode = "function", inherits = TRUE) &&
      cpp_symbol_available("_phynd_score_components_nn_nt_cpp")) {
    return(score_components_nn_nt_cpp(X_nt = X_nt, T_tl = T_tl, wNN = as.numeric(wNN), eps = eps))
  }

  if (ncol(T_tl) == 0L) {
    return(list(ratio = numeric(0), med_nn = numeric(0), med_nt = numeric(0), r2_mean = numeric(0)))
  }

  xnorm2 <- rowSums(X_nt * X_nt)
  nn_idx <- which(wNN < 0.5)
  nt_idx <- which(wNN > 0.5)
  if (length(nn_idx) == 0L || length(nt_idx) == 0L) {
    stop("Need both NN and NT voxels for ratio scoring. Check wNN thresholds.")
  }

  L <- ncol(T_tl)
  ratio <- numeric(L)
  med_nn <- numeric(L)
  med_nt <- numeric(L)
  r2_mean <- numeric(L)

  for (k in seq_len(L)) {
    tk <- T_tl[, k]
    tnorm2 <- sum(tk * tk)
    if (!is.finite(tnorm2) || tnorm2 < eps) {
      ratio[k] <- 0
      med_nn[k] <- 0
      med_nt[k] <- 0
      r2_mean[k] <- 0
      next
    }

    cvec <- as.numeric(X_nt %*% tk)
    r2 <- (cvec * cvec) / pmax(xnorm2 * tnorm2, eps)
    mn <- stats::median(r2[nn_idx], na.rm = TRUE)
    mt <- stats::median(r2[nt_idx], na.rm = TRUE)

    ratio[k] <- mn / pmax(mt, eps)
    med_nn[k] <- mn
    med_nt[k] <- mt
    r2_mean[k] <- mean(r2, na.rm = TRUE)
  }

  list(ratio = ratio, med_nn = med_nn, med_nt = med_nt, r2_mean = r2_mean)
}

compute_split_stability <- function(
  X_nt,
  pca_rank,
  extractor,
  lag_order,
  n_candidates,
  svd_engine = "auto",
  seed = NULL
) {
  T <- ncol(X_nt)
  i1 <- seq_len(floor(T / 2))
  i2 <- seq(from = floor(T / 2) + 1L, to = T)

  stab_half <- function(X_half, n_ref) {
    sv <- truncated_svd(X_half, rank = min(pca_rank, min(dim(X_half))), svd_engine = svd_engine, seed = seed)
    Z <- t(sv$v)
    ext <- extract_dynamic_components(
      Z_kt = Z,
      extractor = extractor,
      lag_order = lag_order,
      n_candidates = n_ref
    )
    ext$scores
  }

  S1 <- stab_half(X_nt[, i1, drop = FALSE], n_ref = n_candidates)
  S2 <- stab_half(X_nt[, i2, drop = FALSE], n_ref = n_candidates)

  list(scores_h1 = S1, scores_h2 = S2, idx_h1 = i1, idx_h2 = i2)
}

stability_filter <- function(T_tl, split_obj, thresh = 0.30, use_cpp = TRUE) {
  if (ncol(T_tl) == 0L) return(logical(0))
  S1 <- split_obj$scores_h1
  S2 <- split_obj$scores_h2
  i1 <- split_obj$idx_h1
  i2 <- split_obj$idx_h2

  if (isTRUE(use_cpp) &&
      exists("stability_filter_cpp", mode = "function", inherits = TRUE) &&
      cpp_symbol_available("_phynd_stability_filter_cpp")) {
    return(as.logical(stability_filter_cpp(
      T_tl = as.matrix(T_tl),
      S1 = as.matrix(S1),
      S2 = as.matrix(S2),
      i1 = as.integer(i1),
      i2 = as.integer(i2),
      thresh = thresh
    )))
  }

  keep <- rep(FALSE, ncol(T_tl))
  for (k in seq_len(ncol(T_tl))) {
    t1 <- safe_scale_vec(T_tl[i1, k])
    t2 <- safe_scale_vec(T_tl[i2, k])

    c1 <- apply(S1, 2, function(v) abs(stats::cor(t1, safe_scale_vec(v))))
    c2 <- apply(S2, 2, function(v) abs(stats::cor(t2, safe_scale_vec(v))))

    m1 <- if (length(c1)) max(c1, na.rm = TRUE) else 0
    m2 <- if (length(c2)) max(c2, na.rm = TRUE) else 0
    keep[k] <- is.finite(m1) && is.finite(m2) && min(m1, m2) >= thresh
  }

  keep
}
