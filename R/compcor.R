# CompCor denoising baselines (aCompCor / tCompCor).

#' CompCor denoising (aCompCor / tCompCor)
#'
#' Removes the principal temporal components of a nuisance voxel set from
#' every analysed voxel, following the component-based noise correction
#' (CompCor) method of Behzadi et al. (2007). It serves as a fixed-rank
#' baseline against which [fast_phy_denoise()] can be compared.
#'
#' @details
#' **Choosing nuisance voxels.**
#' - `mode = "acompcor"` (anatomical CompCor) uses the voxels in
#'   `nuisance_mask` when it is supplied; otherwise voxels whose wNN weight is
#'   below `wnn_thresh`, using `wnn` when supplied and otherwise a map computed
#'   with [compute_wNN_diff_energy()] from the analysed voxels.
#' - `mode = "tcompcor"` (temporal CompCor) uses the `top_var_n` analysed
#'   voxels with the highest temporal variance after detrending (by default
#'   the larger of `5 * n_comp` and `top_var_frac` of the analysed voxels).
#'
#' `nuisance_mask` and `wnn` are indexed over *all* voxels of `x` and are
#' independent of `mask`, so a white-matter/CSF region outside the analysis
#' mask can supply the regressors. A `wnn` of length `sum(mask)` is also
#' accepted and read as covering the masked voxels; `NA` weights never select
#' a voxel.
#'
#' **Estimating components.** Each nuisance voxel time series has its mean and
#' linear trend removed (and the cosine high-pass basis, when
#' `pre_highpass = "dct"`) and is scaled to unit variance before the PCA, as
#' in Behzadi et al. (2007). Components whose singular value is below
#' `1e-6` times the largest are discarded, so no more components are removed
#' than the nuisance data support; a warning is given when fewer than
#' `n_comp` remain.
#'
#' **Cleaning.** The components are projected out of every analysed voxel in
#' one orthogonal projection, and each voxel's mean is added back. Voxels
#' outside `mask`, constant voxels, and voxels with non-finite values are
#' returned unchanged. With `pre_highpass = "dct"` the output is also high-pass
#' filtered (the cosine basis is removed from the analysed voxels); the basis
#' is returned in `diagnostics$highpass$basis`.
#'
#' **Design guardrail.** When `design` is supplied, the regressors are made
#' orthogonal to the design (plus an intercept) before removal. Least-squares
#' task estimates computed from `x_clean` are then identical to those from
#' `x`: the guardrail protects task signal but also leaves any
#' task-correlated physiological noise in place. With `pre_highpass = "dct"`
#' the regressors are also kept orthogonal to the cosine basis, so the output
#' stays high-passed, and the guarantee holds for the high-passed design.
#'
#' **Relation to fMRIPrep.** This is a CompCor-style baseline, not a
#' reimplementation of the fMRIPrep/NiWorkflows confounds pipeline. Mask
#' construction and erosion, censoring, and the component-retention rule
#' (a fixed `n_comp` here, a variance-explained threshold in fMRIPrep) differ.
#' Use fMRIPrep's confounds directly when exact parity is required.
#'
#' @param x Numeric data: a voxels x time matrix, a time x voxels matrix (see
#'   `orientation`), a 4D array with time last, or a `neuroim2::NeuroVec`
#'   image (only voxels inside a `SparseNeuroVec`'s mask are analysed).
#' @param tr Repetition time in seconds.
#' @param mode `"acompcor"` (nuisance voxels from a mask or wNN map) or
#'   `"tcompcor"` (highest-variance voxels).
#' @param n_comp Number of principal components to remove (at least 1).
#' @param mask Optional analysis mask (logical vector, voxel indices, 3D
#'   array, or `neuroim2::LogicalNeuroVol`). Only masked voxels are cleaned.
#' @param nuisance_mask Optional nuisance region for `"acompcor"`, in the same
#'   forms as `mask`, over all voxels of `x`.
#' @param wnn Optional wNN weights for `"acompcor"` (1 = neuronal, 0 =
#'   non-neuronal), over all voxels of `x` or over the masked voxels.
#' @param wnn_thresh Voxels with wNN below this value (in \[0, 1\]) are
#'   nuisance voxels in `"acompcor"`.
#' @param top_var_frac Fraction (in (0, 1\]) of analysed voxels used by
#'   `"tcompcor"` when `top_var_n` is `NULL`.
#' @param top_var_n Optional number of highest-variance voxels for
#'   `"tcompcor"`; clamped to between `n_comp` and the number of analysed
#'   voxels.
#' @param design Optional design matrix (time x predictors). See the design
#'   guardrail in Details.
#' @param pre_highpass `"none"` or `"dct"`: remove a discrete cosine basis
#'   with cutoff `highpass_hz` from the analysed voxels before estimating and
#'   removing components.
#' @param highpass_hz High-pass cutoff in Hz for `pre_highpass = "dct"`
#'   (default 1/128 Hz); must be below the Nyquist frequency `1 / (2 * tr)`.
#' @param orientation `"auto"` (voxels x time unless there are fewer rows than
#'   columns), `"voxels_by_time"`, or `"time_by_voxels"`.
#' @param svd_engine `"auto"`, `"svd"`, or `"rsvd"`; see
#'   [fast_phy_denoise()].
#' @param seed Optional integer seed for the randomized SVD. The caller's
#'   random-number stream is left untouched.
#' @param return_diagnostics Whether to include the `diagnostics` element.
#' @return An object of class `"compcor_denoise_result"`: a list with
#' \describe{
#'   \item{`x_clean`}{Denoised data of the same type, shape, and orientation
#'     as `x` (a `NeuroVec` input gives a `DenseNeuroVec` or
#'     `SparseNeuroVec` in the same space), voxel means retained.}
#'   \item{`regressors`}{Matrix, time x components: the removed components,
#'     each scaled to unit variance (zero columns if none remain).}
#'   \item{`wNN`}{wNN weights over all voxels (`NA` where unavailable) when a
#'     wNN map was used to select nuisance voxels, otherwise `NULL`.}
#'   \item{`nuisance_voxels`}{Logical vector over all voxels: the voxels whose
#'     time series formed the PCA.}
#'   \item{`active`}{Logical vector over all voxels: the voxels that were
#'     cleaned.}
#'   \item{`mode`, `tr`}{As supplied.}
#'   \item{`params`}{List of settings, including `n_comp` (requested) and
#'     `n_comp_used`.}
#'   \item{`diagnostics`}{When `return_diagnostics = TRUE`: input
#'     information, voxel counts, the source of the nuisance selection,
#'     singular values, the high-pass basis, the design columns dropped by the
#'     guardrail, and the wNN computation (for automatic wNN).}
#' }
#' @references
#' Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A component based
#' noise correction method (CompCor) for BOLD and perfusion based fMRI.
#' *NeuroImage*, 37(1), 90--101. \doi{10.1016/j.neuroimage.2007.04.042}
#' @seealso [fast_phy_denoise()], [compute_design_free_qc()]
#' @examples
#' sim <- simulate_phy_data(n_vox = 400, n_time = 120, tr = 2, seed = 1)
#'
#' # aCompCor with a known nuisance region
#' fit <- compcor_denoise(sim$x, tr = 2, mode = "acompcor",
#'                        nuisance_mask = sim$nn, n_comp = 3)
#' fit
#'
#' # tCompCor with a high-pass filter
#' fit_t <- compcor_denoise(sim$x, tr = 2, mode = "tcompcor", n_comp = 3,
#'                          pre_highpass = "dct")
#' dim(fit_t$regressors)
#' @export
compcor_denoise <- function(
  x,
  tr,
  mode = c("acompcor", "tcompcor"),
  n_comp = 5L,
  mask = NULL,
  nuisance_mask = NULL,
  wnn = NULL,
  wnn_thresh = 0.5,
  top_var_frac = 0.02,
  top_var_n = NULL,
  design = NULL,
  pre_highpass = c("none", "dct"),
  highpass_hz = 1 / 128,
  orientation = c("auto", "voxels_by_time", "time_by_voxels"),
  svd_engine = c("auto", "svd", "rsvd"),
  seed = NULL,
  return_diagnostics = TRUE
) {
  if (missing(tr)) {
    stop("`tr` (repetition time in seconds) is required.", call. = FALSE)
  }
  tr <- check_positive(tr, "tr")
  mode <- match.arg(mode)
  pre_highpass <- match.arg(pre_highpass)
  orientation <- match.arg(orientation)
  svd_engine <- match.arg(svd_engine)
  n_comp <- check_count(n_comp, "n_comp")
  wnn_thresh <- check_number(wnn_thresh, "wnn_thresh", lower = 0, upper = 1)
  top_var_frac <- check_number(top_var_frac, "top_var_frac", lower = 0, upper = 1)
  if (top_var_frac <= 0) {
    stop("`top_var_frac` must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(top_var_n)) {
    top_var_n <- check_count(top_var_n, "top_var_n")
  }
  highpass_hz <- check_positive(highpass_hz, "highpass_hz")
  if (pre_highpass == "dct" && highpass_hz >= 1 / (2 * tr)) {
    stop(sprintf(
      "`highpass_hz` (%g Hz) must be below the Nyquist frequency 1/(2*tr) = %g Hz.",
      highpass_hz, 1 / (2 * tr)
    ), call. = FALSE)
  }
  seed <- check_number(seed, "seed", null_ok = TRUE)
  return_diagnostics <- check_flag(return_diagnostics, "return_diagnostics")
  if (mode == "tcompcor" && (!is.null(nuisance_mask) || !is.null(wnn))) {
    stop("`nuisance_mask` and `wnn` apply only to mode = \"acompcor\".", call. = FALSE)
  }
  if (!is.null(nuisance_mask) && !is.null(wnn)) {
    stop("Supply either `nuisance_mask` or `wnn`, not both.", call. = FALSE)
  }

  parsed <- as_nt_matrix(x, mask = mask, orientation = orientation)
  X <- parsed$X
  N <- nrow(X)
  T <- ncol(X)
  if (T < 8L) {
    stop("Need at least 8 time points.", call. = FALSE)
  }

  # Candidate nuisance voxels over the full voxel set.
  wnn_full <- NULL
  wres <- NULL
  source_keep <- parsed$keep
  if (mode == "tcompcor") {
    selector_source <- "top_variance"
  } else if (!is.null(nuisance_mask)) {
    source_keep <- tryCatch(
      normalize_mask(nuisance_mask, N),
      error = function(e) {
        stop(sprintf("Invalid `nuisance_mask`: %s", conditionMessage(e)), call. = FALSE)
      }
    )
    selector_source <- "nuisance_mask"
  } else if (!is.null(wnn)) {
    wnn_full <- expand_given_wnn(wnn, parsed$keep)
    source_keep <- !is.na(wnn_full) & wnn_full < wnn_thresh
    selector_source <- "wnn_given"
  } else {
    selector_source <- "wnn_auto"
  }

  vox <- analysis_voxels(X, parsed$keep | source_keep)
  analyzable <- vox$active
  active <- parsed$keep & analyzable
  if (sum(active) < 20L) {
    stop("Need at least 20 analysable voxels (inside the mask, finite, and non-constant).", call. = FALSE)
  }

  # Analysed voxels: centred, optionally high-passed.
  hp_basis <- NULL
  if (pre_highpass == "dct") {
    hp_basis <- dct_basis(T, tr, highpass_hz)
  }
  prep <- function(rows) {
    Z <- X[rows, , drop = FALSE]
    Z <- Z - rowMeans(Z)
    if (!is.null(hp_basis)) {
      Z <- Z - (Z %*% hp_basis) %*% t(hp_basis)
    }
    Z
  }
  Xa <- X[active, , drop = FALSE]
  mu <- rowMeans(Xa)
  Xc <- prep(active)

  # Constant and linear trend (plus the cosine basis when high-passing) are
  # removed in one joint projection, so the nuisance data stay orthogonal to
  # every one of them.
  trend <- cbind(1, seq_len(T) - (T + 1) / 2)
  if (!is.null(hp_basis)) {
    trend <- cbind(hp_basis, trend[, 2])
  }
  trend_qr <- qr(trend)
  trend_q <- qr.Q(trend_qr)[, seq_len(trend_qr$rank), drop = FALSE]
  detrend <- function(Z) Z - (Z %*% trend_q) %*% t(trend_q)

  if (selector_source == "wnn_auto") {
    wres <- compute_wNN_diff_energy(Xa - mu, threshold = "percentile")
    wnn_full <- expand_voxel_vector(wres$wNN, active)
    source <- active & !is.na(wnn_full) & wnn_full < wnn_thresh
  } else if (selector_source == "top_variance") {
    v <- rowSums(detrend(Xc)^2)
    n_active <- sum(active)
    if (is.null(top_var_n)) {
      top_var_n <- max(n_comp * 5L, ceiling(top_var_frac * n_active))
    }
    top_var_n <- as.integer(min(max(top_var_n, n_comp), n_active))
    source <- rep(FALSE, N)
    source[which(active)[order(v, decreasing = TRUE)[seq_len(top_var_n)]]] <- TRUE
  } else {
    source <- source_keep & analyzable
  }

  if (!any(source)) {
    stop("No usable voxels selected for CompCor regressors.", call. = FALSE)
  }

  # Nuisance matrix: detrended and variance-normalized (Behzadi et al., 2007).
  # A series that is (numerically) all trend or low frequency -- judged
  # against its own variance before filtering -- has nothing left to
  # estimate from; scaling its roundoff residue to unit variance would
  # manufacture spurious components.
  Xsrc <- X[source, , drop = FALSE]
  s_raw <- sqrt(rowSums((Xsrc - rowMeans(Xsrc))^2) / (T - 1))
  Xs <- detrend(prep(source))
  s <- sqrt(rowSums(Xs * Xs) / (T - 1))
  usable <- s > 1e-8 * s_raw
  Xs <- Xs[usable, , drop = FALSE] / s[usable]

  dof <- T - 2L - (if (is.null(hp_basis)) 0L else ncol(hp_basis) - 1L)
  rank_req <- min(n_comp, nrow(Xs), max(dof, 1L))
  d <- numeric(0)
  svd_engine_used <- NA_character_
  Tcomp <- matrix(0, nrow = T, ncol = 0)
  if (rank_req >= 1L) {
    sv <- truncated_svd(Xs, rank = rank_req, svd_engine = svd_engine, seed = seed)
    d <- sv$d[seq_len(min(rank_req, length(sv$d)))]
    svd_engine_used <- sv$engine
    keep_comp <- which(d > 1e-6 * max(d, 0))
    if (length(keep_comp)) {
      Tcomp <- sv$v[, keep_comp, drop = FALSE]
      Tcomp <- matrix(apply(Tcomp, 2, safe_scale_vec), nrow = T, ncol = length(keep_comp))
    }
  }
  n_comp_used <- ncol(Tcomp)
  if (n_comp_used < n_comp) {
    warning(sprintf(
      "Only %d of the requested %d CompCor components could be estimated from %d nuisance voxel(s); using %d.",
      n_comp_used, n_comp, sum(source), n_comp_used
    ), call. = FALSE)
  }
  n_before_design <- ncol(Tcomp)
  design_dropped <- integer(0)
  if (!is.null(design)) {
    # With high-passing, guard against the design and the cosine basis
    # together; guarding against the design alone reintroduces low
    # frequencies into the output.
    guard <- as.matrix(design)
    if (!is.null(hp_basis)) {
      guard <- cbind(guard, hp_basis)
    }
    Tcomp <- orthogonalize_to_design(Tcomp, guard)
    design_dropped <- setdiff(seq_len(n_before_design), attr(Tcomp, "kept"))
    attr(Tcomp, "kept") <- NULL
  }
  if (ncol(Tcomp)) {
    colnames(Tcomp) <- paste0("compcor_", seq_len(ncol(Tcomp)))
  }

  Xclean <- project_out_components(Xc, Tcomp) + mu
  Xout <- restore_output(Xclean, parsed, active)

  out <- list(
    x_clean = Xout,
    regressors = Tcomp,
    wNN = wnn_full,
    nuisance_voxels = source,
    active = active,
    mode = mode,
    tr = tr,
    params = list(
      mode = mode,
      n_comp = n_comp,
      n_comp_used = ncol(Tcomp),
      wnn_thresh = wnn_thresh,
      top_var_frac = top_var_frac,
      top_var_n = if (mode == "tcompcor") top_var_n else NULL,
      design = !is.null(design),
      pre_highpass = pre_highpass,
      highpass_hz = highpass_hz,
      svd_engine = svd_engine
    )
  )

  if (return_diagnostics) {
    out$diagnostics <- list(
      input = parsed$info,
      selector_source = selector_source,
      n_nuisance_voxels = sum(source),
      n_active = sum(active),
      n_excluded_nonfinite = vox$n_nonfinite,
      n_excluded_constant = vox$n_constant,
      singular_values = d,
      svd_engine_used = svd_engine_used,
      highpass = list(
        basis = hp_basis,
        n_regressors = if (is.null(hp_basis)) 0L else ncol(hp_basis)
      ),
      design_dropped = design_dropped,
      wnn = wres
    )
  }

  class(out) <- c("compcor_denoise_result", "list")
  out
}

#' @method print compcor_denoise_result
#' @export
print.compcor_denoise_result <- function(x, ...) {
  p <- x$params
  cat(sprintf("<compcor_denoise_result> mode: %s\n", x$mode))
  cat(sprintf(
    "  voxels cleaned: %d of %d; nuisance voxels: %d\n",
    sum(x$active), length(x$active), sum(x$nuisance_voxels)
  ))
  cat(sprintf(
    "  components removed: %d (requested %d)\n",
    ncol(x$regressors), p$n_comp
  ))
  if (identical(p$pre_highpass, "dct")) {
    cat(sprintf("  high-pass: DCT, cutoff %g Hz\n", p$highpass_hz))
  }
  if (isTRUE(p$design)) {
    cat("  design guardrail: regressors orthogonal to design\n")
  }
  invisible(x)
}

# Validate a user-supplied wNN vector and return it over all voxels.
expand_given_wnn <- function(wnn, keep) {
  if (!is.numeric(wnn)) {
    stop("`wnn` must be a numeric vector.", call. = FALSE)
  }
  wnn <- as.numeric(wnn)
  n <- length(keep)
  if (length(wnn) == n) {
    out <- wnn
  } else if (length(wnn) == sum(keep)) {
    out <- rep(NA_real_, n)
    out[keep] <- wnn
  } else {
    stop(sprintf(
      "`wnn` has length %d; expected %d (all voxels) or %d (masked voxels).",
      length(wnn), n, sum(keep)
    ), call. = FALSE)
  }
  fin <- out[!is.na(out)]
  if (!length(fin)) {
    stop("`wnn` contains no non-missing values.", call. = FALSE)
  }
  if (any(!is.finite(fin)) || any(fin < 0 | fin > 1)) {
    stop("`wnn` values must lie in [0, 1] (NA marks voxels without a weight).", call. = FALSE)
  }
  out
}

# Orthonormal DCT-II high-pass basis (constant plus cosines below the cutoff),
# matching SPM's spm_filter and nipype's cosine drift regressors.
dct_basis <- function(n_time, tr, highpass_hz) {
  kmax <- floor(2 * n_time * tr * highpass_hz)
  if (!is.finite(kmax) || kmax < 0) {
    kmax <- 0
  }
  if (kmax + 1 > n_time - 3) {
    stop(sprintf(
      "High-pass cutoff %g Hz would remove %d of %d temporal degrees of freedom.",
      highpass_hz, as.integer(kmax + 1), n_time
    ), call. = FALSE)
  }
  t_idx <- seq_len(n_time)
  B <- matrix(1 / sqrt(n_time), nrow = n_time, ncol = kmax + 1L)
  for (k in seq_len(kmax)) {
    B[, k + 1L] <- sqrt(2 / n_time) * cos(pi * (2 * t_idx - 1) * k / (2 * n_time))
  }
  B
}

dct_highpass_residualize <- function(X_nt, tr, highpass_hz) {
  X_nt <- as.matrix(X_nt)
  B <- dct_basis(ncol(X_nt), tr, highpass_hz)
  list(X = X_nt - (X_nt %*% B) %*% t(B), n_regressors = ncol(B), basis = B)
}
