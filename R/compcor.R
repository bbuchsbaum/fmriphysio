# Reduced CompCor-style denoising baselines.

#' CompCor-style denoising (aCompCor / tCompCor)
#'
#' @param x Matrix-like data (voxels x time or time x voxels).
#' @param tr Repetition time in seconds.
#' @param mode One of "acompcorr" or "tcompcorr".
#' @param n_comp Number of principal regressors to remove.
#' @param mask Optional voxel mask (logical or indices) applied before processing.
#' @param nuisance_mask Optional voxel mask for aCompCor ROI (logical or indices).
#' @param wnn Optional wNN vector; used for aCompCor when nuisance_mask is not provided.
#' @param wnn_thresh wNN threshold for aCompCor voxel selection.
#' @param top_var_frac Fraction of highest-variance voxels for tCompCor.
#' @param top_var_n Optional explicit number of voxels for tCompCor.
#' @param design Optional design matrix (time x predictors) for guardrail orthogonalization.
#' @param center_rows_first Whether to mean-center voxel time series first.
#' @param pre_highpass Optional pre-PCA high-pass mode ("none" or "dct").
#' @param highpass_hz High-pass cutoff in Hz for `pre_highpass = "dct"`.
#' @param svd_engine SVD backend for PCA ("auto", "rsvd", "svd").
#' @param seed Optional random seed.
#' @param return_diagnostics Return detailed diagnostics.
#' @details
#' This function provides an explicit reduced baseline that maps onto CompCor:
#'
#' - `mode = "acompcorr"` reduces to an anatomical/noise-ROI CompCor variant:
#'   nuisance voxels are selected from `nuisance_mask` when provided; otherwise
#'   from `wnn < wnn_thresh` (or auto-derived `wNN`).
#' - `mode = "tcompcorr"` reduces to temporal CompCor:
#'   nuisance voxels are selected by highest temporal variance.
#'
#' In both modes, principal nuisance regressors are extracted from the selected
#' voxel set and projected out from all voxels in one linear projection step.
#'
#' This is a CompCor-style baseline, not a strict reimplementation of the
#' fMRIPrep/NiWorkflows confounds pipeline. Important differences may include:
#'
#' - preprocessing before CompCor (for example, high-pass/censoring strategy),
#' - exact anatomical mask construction and erosion rules,
#' - component retention policy (fixed `n_comp` here vs variance-threshold
#'   policies often used in fMRIPrep outputs).
#'
#' Use this function for in-package baseline comparisons. If exact confound
#' parity with fMRIPrep is required, use fMRIPrep confounds directly or add a
#' dedicated compatibility mode.
#' @return A list with cleaned data, regressors, and diagnostics.
#' @export
compcor_denoise <- function(
  x,
  tr,
  mode = c("acompcorr", "tcompcorr"),
  n_comp = 5L,
  mask = NULL,
  nuisance_mask = NULL,
  wnn = NULL,
  wnn_thresh = 0.5,
  top_var_frac = 0.02,
  top_var_n = NULL,
  design = NULL,
  center_rows_first = TRUE,
  pre_highpass = c("none", "dct"),
  highpass_hz = 1 / 128,
  svd_engine = c("auto", "rsvd", "svd"),
  seed = NULL,
  return_diagnostics = TRUE
) {
  if (missing(tr) || !is.numeric(tr) || length(tr) != 1L || tr <= 0) {
    stop("tr must be a positive scalar.")
  }

  mode <- match.arg(mode)
  pre_highpass <- match.arg(pre_highpass)
  svd_engine <- match.arg(svd_engine)
  n_comp <- as.integer(n_comp)
  if (n_comp < 1L) stop("n_comp must be >= 1")
  if (!is.numeric(highpass_hz) || length(highpass_hz) != 1L || highpass_hz <= 0) {
    stop("highpass_hz must be a positive scalar.")
  }

  parsed <- as_nt_matrix(x, mask = mask)
  X0 <- parsed$X
  N <- nrow(X0)
  T <- ncol(X0)
  if (T < 8L) stop("Need at least 8 time points.")
  if (N < 20L) stop("Need at least 20 voxels.")

  if (isTRUE(center_rows_first)) {
    X0 <- center_rows(X0)
  }

  hp_info <- list(mode = pre_highpass, highpass_hz = highpass_hz, n_regressors = 0L)
  if (pre_highpass == "dct") {
    hp <- dct_highpass_residualize(X0, tr = tr, highpass_hz = highpass_hz)
    X0 <- hp$X
    hp_info$n_regressors <- hp$n_regressors
  }

  selector <- logical(N)
  selector_source <- "none"
  wres <- NULL

  if (mode == "acompcorr") {
    if (!is.null(nuisance_mask)) {
      selector <- normalize_selector(nuisance_mask, N)
      selector_source <- "nuisance_mask"
    } else if (!is.null(wnn)) {
      if (length(wnn) != N) stop("wnn length must match voxel dimension.")
      selector <- as.numeric(wnn) < wnn_thresh
      selector_source <- "wnn_given"
    } else {
      wres <- compute_wNN_diff_energy(X0, threshold = "percentile")
      selector <- wres$wNN < wnn_thresh
      selector_source <- "wnn_auto"
    }
  } else {
    sdev <- apply(X0, 1, stats::sd)
    if (is.null(top_var_n)) {
      top_var_n <- max(n_comp * 5L, ceiling(top_var_frac * N))
    }
    top_var_n <- min(max(as.integer(top_var_n), n_comp), N)
    idx <- order(sdev, decreasing = TRUE)[seq_len(top_var_n)]
    selector[idx] <- TRUE
    selector_source <- "top_variance"
  }

  if (!any(selector)) {
    stop("No voxels selected for CompCor regressors.")
  }

  Xsel <- X0[selector, , drop = FALSE]
  rank_use <- min(n_comp, min(dim(Xsel)))
  sv <- truncated_svd(Xsel, rank = rank_use, svd_engine = svd_engine, seed = seed)

  # Right singular vectors are temporal PCs (T x k).
  Tcomp <- sv$v[, seq_len(rank_use), drop = FALSE]
  Tcomp <- apply(Tcomp, 2, safe_scale_vec)
  Tcomp <- as.matrix(Tcomp)

  if (!is.null(design)) {
    Tcomp <- orthogonalize_to_design(Tcomp, design)
  }

  Xclean <- project_out_components(X0, Tcomp)
  Xout <- restore_orientation(Xclean, parsed$info)

  out <- list(
    x_clean = Xout,
    regressors = Tcomp,
    mode = mode,
    tr = tr
  )

  if (return_diagnostics) {
    out$diagnostics <- list(
      input = parsed$info,
      n_selected_voxels = sum(selector),
      selected_fraction = mean(selector),
      selector_source = selector_source,
      selector = selector,
      wnn = wres,
      params = list(
        n_comp = n_comp,
        wnn_thresh = wnn_thresh,
        top_var_frac = top_var_frac,
        top_var_n = top_var_n,
        pre_highpass = pre_highpass,
        highpass_hz = highpass_hz,
        highpass_n_regressors = hp_info$n_regressors,
        svd_engine = svd_engine
      )
    )
  }

  class(out) <- c("compcor_denoise_result", class(out))
  out
}

normalize_selector <- function(sel, n) {
  if (is.logical(sel)) {
    if (length(sel) != n) stop("Logical selector length must match voxel dimension.")
    return(sel)
  }
  if (is.numeric(sel)) {
    out <- rep(FALSE, n)
    idx <- as.integer(sel)
    if (length(idx) == 0L) return(out)
    if (any(idx < 1L | idx > n)) stop("Numeric selector contains out-of-range indices.")
    out[idx] <- TRUE
    return(out)
  }
  stop("Selector must be logical or numeric indices.")
}

dct_highpass_residualize <- function(X_nt, tr, highpass_hz) {
  X_nt <- as.matrix(X_nt)
  T <- ncol(X_nt)
  if (T < 4L) {
    return(list(X = X_nt, n_regressors = 0L))
  }

  kmax <- floor(2 * T * tr * highpass_hz)
  if (!is.finite(kmax) || kmax <= 0) {
    return(list(X = X_nt, n_regressors = 0L))
  }

  t_idx <- seq_len(T)
  D <- matrix(0, nrow = T, ncol = kmax + 1L)
  D[, 1] <- 1
  for (k in seq_len(kmax)) {
    D[, k + 1L] <- cos(pi * (2 * t_idx - 1) * k / (2 * T))
  }

  qd <- qr(D)
  Q <- qr.Q(qd)[, seq_len(qd$rank), drop = FALSE]
  X_hp <- X_nt - (X_nt %*% Q) %*% t(Q)
  list(X = X_hp, n_regressors = qd$rank)
}
