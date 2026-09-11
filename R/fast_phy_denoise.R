#' Fast design-free physiological denoising (PHYCAA+-style)
#'
#' Removes physiological noise from fMRI time series without a task design,
#' following the logic of PHYCAA+ with one-shot linear algebra in place of its
#' repeated searches:
#'
#' 1. **Tissue weights.** A neuronal-tissue weight `wNN` is computed per voxel
#'    from temporal-difference energy ([compute_wNN_diff_energy()]); 1 marks
#'    likely neuronal voxels, 0 likely non-neuronal voxels (vasculature,
#'    cerebrospinal fluid).
#' 2. **Low-rank dynamics.** The voxel-centred data are weighted per voxel --
#'    by default by `1 - wNN`, emphasising non-neuronal tissue, where
#'    physiological noise is strongest (see `extraction_weight`) -- and reduced
#'    to their leading `min(pca_rank, floor(T / 3))` temporal singular vectors.
#'    Temporally predictable components are extracted from that subspace
#'    ([extract_caa_oneshot()], or DiCCA/DiPCA via the 'dipca' package).
#' 3. **Selection.** Each candidate is scored by the ratio of its median
#'    voxelwise \eqn{R^2} in non-neuronal voxels to that in neuronal voxels,
#'    computed on the unweighted data. Candidates with ratio above
#'    `ratio_thresh` (and optionally predictability above `pred_thresh` and
#'    split-half stability) are nuisance regressors.
#' 4. **Removal.** The selected regressors are projected out of every
#'    analysed voxel in one orthogonal projection, and voxel means are
#'    restored. With `max_passes > 1`, extraction and selection repeat on the
#'    cleaned data and the union of all selected regressors is projected out
#'    of the original data.
#'
#' Voxels outside `mask`, voxels containing non-finite values, and voxels
#' that are constant over time are excluded from the analysis and returned
#' unchanged. `x_clean` has the type and shape of `x`: matrices and arrays
#' keep their dimensions and orientation, and a `neuroim2::NeuroVec` is
#' returned as a `DenseNeuroVec` -- or, for sparse input, a `SparseNeuroVec`
#' with the same mask -- in the same `NeuroSpace`. Only voxels inside a
#' sparse image's mask are analysed.
#'
#' **Design guardrail.** When `design` is supplied, candidate regressors are
#' orthogonalised against the design (plus an intercept) before scoring and
#' removal. Denoising therefore cannot change task-effect estimates from a
#' GLM with that design; the trade-off is that task-correlated physiological
#' noise is left in place.
#'
#' @section Departures from PHYCAA+ and limitations:
#' PHYCAA+ estimates components from data weighted toward neuronal tissue
#' (`extraction_weight = "neuronal"`) and searches over many subspace ranks.
#' In simulations with known physiological and neural sources, weighting
#' toward non-neuronal tissue with a small fixed rank (`pca_rank = 20`) and
#' `ratio_thresh = 2` removed substantially more physiological signal while
#' preserving more neural signal than the PHYCAA+ weighting or larger ranks,
#' so these are the defaults. Larger ranks let the lag-1 canonical analysis
#' fit spurious components.
#'
#' The defaults are better on average, not in every case. When physiological
#' rhythms alias to nearby frequencies -- to each other, or into the
#' low-frequency band occupied by neural signal -- they can be mixed within a
#' candidate whose ratio falls below `ratio_thresh`, and part of that noise is
#' left in the data.
#'
#' The selection criterion assumes the data contain physiological noise
#' concentrated in non-neuronal tissue. Without such noise, the
#' difference-energy weights no longer reflect tissue type (they follow
#' incidental differences in voxel noise or signal amplitude), and components
#' carrying neural signal can then be selected and removed. Inspect
#' `component_table` and the removed `regressors`, and do not apply the method
#' to data that have already been physiologically corrected.
#'
#' @param x Numeric data: a matrix (voxels x time by default; see
#'   `orientation`), a 3D/4D array with time as the last dimension, or a
#'   `neuroim2::NeuroVec` image (for example a `DenseNeuroVec` or
#'   `SparseNeuroVec`).
#' @param tr Optional repetition time in seconds. Not used by the algorithm;
#'   stored in the result as metadata.
#' @param design Optional task design matrix (time points x predictors) used
#'   as a guardrail; see Details.
#' @param mask Optional voxel mask: a logical vector, whole-number voxel
#'   indices, a spatial array, or a `neuroim2::LogicalNeuroVol`. Voxels
#'   outside the mask are returned unchanged.
#' @param orientation Layout of a matrix `x`: `"voxels_by_time"`,
#'   `"time_by_voxels"`, or `"auto"` (default), which treats a matrix with
#'   fewer rows than columns as time x voxels and says so with a message.
#' @param wnn Optional precomputed neuronal-tissue weights in \[0, 1\]
#'   (1 = neuronal), one per voxel or one per masked voxel. Overrides the
#'   internal computation.
#' @param wnn_method Method for computing `wNN`; currently only
#'   `"diff_energy"`.
#' @param wnn_threshold `"percentile"` or `"mixture"`; passed to
#'   [compute_wNN_diff_energy()] as `threshold`.
#' @param delta_nt,delta_nn Quantile probabilities in (0, 1) passed to
#'   [compute_wNN_diff_energy()] as `nt_q` and `nn_q`. Used only when `wnn`
#'   is not supplied; `delta_nn` is ignored by the mixture threshold.
#' @param extraction_weight Per-voxel weighting applied before the low-rank
#'   decomposition: `"non_neuronal"` (default; weight `1 - wNN`),
#'   `"neuronal"` (weight `wNN`, as in PHYCAA+), or `"none"`. Scoring always
#'   uses the unweighted data.
#' @param extractor Component extractor: `"caa"` (built in), or `"dicca"` /
#'   `"dipca"`, which require the 'dipca' package from
#'   <https://bbuchsbaum.r-universe.dev>.
#' @param lag_order Lag order for DiCCA/DiPCA (whole number >= 1). The CAA
#'   extractor always uses lag 1.
#' @param pca_rank Maximum rank of the low-rank subspace (default 20). The
#'   rank used is `min(pca_rank, floor(T / 3), number of analysed voxels)`.
#' @param n_candidates Maximum number of candidate components per pass.
#' @param ratio_thresh Minimum ratio of median \eqn{R^2} in non-neuronal
#'   voxels to median \eqn{R^2} in neuronal voxels for a candidate to be
#'   selected (strictly greater). The default 2 selects components expressed
#'   more than twice as strongly in non-neuronal as in neuronal tissue.
#' @param pred_thresh Optional minimum predictability in \[0, 1\]: the
#'   canonical lag-1 correlation for `"caa"`, or the DiCCA/DiPCA \eqn{R^2}.
#'   `NULL` (default) disables the filter.
#' @param stability `"none"` (default) or `"half"`. With `"half"`, candidates
#'   are re-extracted independently in each half of the run and a candidate
#'   is kept only if its voxelwise \eqn{R^2} map is reproduced in both halves
#'   (correlation of maps >= `stability_thresh`).
#' @param stability_thresh Minimum correlation between spatial \eqn{R^2} maps
#'   for split-half stability, in \[0, 1\].
#' @param max_passes Number of extraction/removal passes (whole number >= 1;
#'   1 or 2 recommended).
#' @param return_diagnostics Whether to include the `diagnostics` element.
#' @param svd_engine Low-rank decomposition: `"auto"` (exact Gram-matrix
#'   eigendecomposition when the smaller data dimension -- usually the number
#'   of time points -- is at most 2000, randomized SVD otherwise), `"svd"`
#'   (exact LAPACK SVD), or `"rsvd"` (randomized).
#' @param use_cpp Use the compiled scoring kernel (identical results; faster
#'   for large data).
#' @param seed Optional seed for the randomized SVD and the DiCCA/DiPCA
#'   extractors. The caller's random number stream is restored afterwards.
#' @return An object of class `fast_phy_denoise_result`, a list with
#'   \describe{
#'     \item{`x_clean`}{Denoised data of the same type, shape, and
#'       orientation as `x` (a `NeuroVec` input gives a `DenseNeuroVec` or
#'       `SparseNeuroVec` in the same space); voxel means are preserved and
#'       excluded voxels are unchanged.}
#'     \item{`wNN`}{Neuronal-tissue weight per voxel (1 = neuronal), `NA` for
#'       voxels excluded from the analysis.}
#'     \item{`regressors`}{Time points x selected components matrix of the
#'       nuisance regressors that were projected out (zero columns if none).}
#'     \item{`component_table`}{Data frame with one row per scored candidate:
#'       `pass`, `candidate`, `ratio_nn_nt`, `med_nn`, `med_nt`,
#'       `predictability`, `stability` (split-half score, `NA` unless
#'       `stability = "half"`), and `selected`.}
#'     \item{`active`}{Logical vector, one per voxel: `TRUE` for voxels that
#'       were analysed and denoised.}
#'     \item{`tr`}{The supplied `tr`, or `NULL`.}
#'     \item{`params`}{The main settings used, including the rank actually
#'       used.}
#'     \item{`diagnostics`}{If requested: input and voxel-exclusion summaries,
#'       the full `wNN` computation, a per-pass log, and timings.}
#'   }
#' @references
#' Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized, adaptive
#' procedure for measuring and controlling physiological noise in BOLD fMRI.
#' *NeuroImage*, 82, 306--325. \doi{10.1016/j.neuroimage.2013.05.102}
#'
#' Churchill, N. W., Yourganov, G., Spring, R., Rasmussen, P. M., Lee, W.,
#' Ween, J. E., & Strother, S. C. (2012). PHYCAA: Data-driven measurement and
#' removal of physiological noise in BOLD fMRI. *NeuroImage*, 59(2),
#' 1299--1314. \doi{10.1016/j.neuroimage.2011.08.021}
#' @seealso [compcor_denoise()] for a CompCor baseline and
#'   [compute_design_free_qc()] for quality-control summaries.
#' @examples
#' sim <- simulate_phy_data(n_vox = 400, n_time = 150, seed = 1)
#' fit <- fast_phy_denoise(sim$x, tr = sim$tr)
#' fit
#' head(fit$component_table)
#' @export
fast_phy_denoise <- function(
  x,
  tr = NULL,
  design = NULL,
  mask = NULL,
  orientation = c("auto", "voxels_by_time", "time_by_voxels"),
  wnn = NULL,
  wnn_method = c("diff_energy"),
  wnn_threshold = c("percentile", "mixture"),
  delta_nt = 0.74,
  delta_nn = 0.95,
  extraction_weight = c("non_neuronal", "neuronal", "none"),
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  pca_rank = 20L,
  n_candidates = 30L,
  ratio_thresh = 2,
  pred_thresh = NULL,
  stability = c("none", "half"),
  stability_thresh = 0.30,
  max_passes = 1L,
  return_diagnostics = TRUE,
  svd_engine = c("auto", "svd", "rsvd"),
  use_cpp = TRUE,
  seed = NULL
) {
  t0_total <- proc.time()[["elapsed"]]

  tr <- check_number(tr, "tr", lower = .Machine$double.eps, null_ok = TRUE)
  orientation <- match.arg(orientation)
  wnn_method <- match.arg(wnn_method)
  wnn_threshold <- match.arg(wnn_threshold)
  extraction_weight <- match.arg(extraction_weight)
  extractor <- match.arg(extractor)
  stability <- match.arg(stability)
  svd_engine <- match.arg(svd_engine)
  delta_nt <- check_number(delta_nt, "delta_nt", lower = 0, upper = 1)
  delta_nn <- check_number(delta_nn, "delta_nn", lower = 0, upper = 1)
  lag_order <- check_count(lag_order, "lag_order")
  pca_rank <- check_count(pca_rank, "pca_rank")
  n_candidates <- check_count(n_candidates, "n_candidates")
  ratio_thresh <- check_number(ratio_thresh, "ratio_thresh", lower = 0)
  pred_thresh <- check_number(pred_thresh, "pred_thresh", lower = 0, upper = 1, null_ok = TRUE)
  stability_thresh <- check_number(stability_thresh, "stability_thresh", lower = 0, upper = 1)
  max_passes <- check_count(max_passes, "max_passes")
  return_diagnostics <- check_flag(return_diagnostics, "return_diagnostics")
  use_cpp <- check_flag(use_cpp, "use_cpp")

  parsed <- as_nt_matrix(x, mask = mask, orientation = orientation)
  vox <- analysis_voxels(parsed$X, parsed$keep)
  active <- vox$active
  X_act <- parsed$X[active, , drop = FALSE]
  N <- nrow(X_act)
  T <- ncol(X_act)

  if (T < 8L) stop("Need at least 8 time points.", call. = FALSE)
  if (N < 20L) {
    stop(sprintf("Need at least 20 analysed voxels; found %d after masking and excluding constant/non-finite voxels.", N), call. = FALSE)
  }

  mu <- rowMeans(X_act)
  X0 <- X_act - mu

  t0_wnn <- proc.time()[["elapsed"]]
  if (is.null(wnn)) {
    wres <- compute_wNN_diff_energy(
      X_nt = X0,
      nt_q = delta_nt,
      nn_q = delta_nn,
      threshold = wnn_threshold
    )
    wNN <- wres$wNN
  } else {
    wNN <- wnn_for_active_voxels(wnn, parsed, active)
    wres <- list(wNN = wNN, threshold_method = "given")
  }
  wnn_seconds <- proc.time()[["elapsed"]] - t0_wnn

  if (!any(wNN < 0.5) || !any(wNN > 0.5)) {
    stop("wNN must classify some voxels as non-neuronal (< 0.5) and some as neuronal (> 0.5).", call. = FALSE)
  }
  w_extract <- switch(
    extraction_weight,
    non_neuronal = 1 - wNN,
    neuronal = wNN,
    none = rep(1, N)
  )
  # Voxels with zero extraction weight contribute nothing to the subspace.
  rank_use <- min(pca_rank, floor(T / 3), sum(w_extract > 0))

  pass_log <- list()
  pass_timing <- list()
  all_selected <- list()
  all_component_tables <- list()
  Xcur <- X0

  for (pass in seq_len(max_passes)) {
    t_pass_start <- proc.time()[["elapsed"]]
    timing <- c(svd_seconds = 0, extract_seconds = 0, score_seconds = 0, stability_seconds = 0, regress_seconds = 0)
    Xw <- Xcur * w_extract

    t0 <- proc.time()[["elapsed"]]
    sv <- truncated_svd(Xw, rank = rank_use, svd_engine = svd_engine, seed = seed, return_u = FALSE)
    nonnull <- sv$d > null_sv_tol * max(sv$d[1], .Machine$double.eps)
    timing[["svd_seconds"]] <- proc.time()[["elapsed"]] - t0

    reason <- NULL
    Tcomp <- matrix(0, nrow = T, ncol = 0)
    if (any(nonnull)) {
      t0 <- proc.time()[["elapsed"]]
      ext <- extract_dynamic_components(
        Z_kt = t(sv$v[, nonnull, drop = FALSE]),
        extractor = extractor,
        lag_order = lag_order,
        n_candidates = n_candidates,
        seed = seed
      )
      Tcomp <- as.matrix(ext$scores)
      candidate <- seq_len(ncol(Tcomp))
      pred <- rep_len(c(ext$predictability, rep(NA_real_, ncol(Tcomp))), ncol(Tcomp))
      timing[["extract_seconds"]] <- proc.time()[["elapsed"]] - t0

      if (!is.null(design)) {
        Tcomp <- orthogonalize_to_design(Tcomp, design)
        kept <- attr(Tcomp, "kept")
        candidate <- candidate[kept]
        pred <- pred[kept]
        if (ncol(Tcomp) == 0L) reason <- "no_components_after_guardrail"
      }
    } else {
      reason <- "no_variance_left"
    }

    if (!is.null(reason)) {
      pass_log[[pass]] <- list(pass = pass, selected = integer(0), reason = reason)
      pass_timing[[pass]] <- data.frame(pass = pass, as.list(timing), pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start)
      break
    }

    t0 <- proc.time()[["elapsed"]]
    # Score enrichment on unweighted data: weighting would suppress the
    # non-neuronal variance that the ratio is meant to detect.
    ss <- score_components_nn_nt(X_nt = Xcur, T_tl = Tcomp, wNN = wNN, use_cpp = use_cpp)
    timing[["score_seconds"]] <- proc.time()[["elapsed"]] - t0

    keep <- is.finite(ss$ratio) & ss$ratio > ratio_thresh
    if (!is.null(pred_thresh)) {
      keep <- keep & !is.na(pred) & pred >= pred_thresh
    }

    stab_score <- rep(NA_real_, ncol(Tcomp))
    if (stability == "half" && any(keep)) {
      t0 <- proc.time()[["elapsed"]]
      stab <- component_stability(
        X_nt = Xcur,
        weights = w_extract,
        T_tl = Tcomp,
        pca_rank = pca_rank,
        extractor = extractor,
        lag_order = lag_order,
        n_candidates = n_candidates,
        thresh = stability_thresh,
        svd_engine = svd_engine,
        seed = seed
      )
      stab_score <- stab$score
      keep <- keep & stab$stable
      timing[["stability_seconds"]] <- proc.time()[["elapsed"]] - t0
    }

    all_component_tables[[pass]] <- data.frame(
      pass = pass,
      candidate = candidate,
      ratio_nn_nt = ss$ratio,
      med_nn = ss$med_nn,
      med_nt = ss$med_nt,
      predictability = pred,
      stability = stab_score,
      selected = keep,
      stringsAsFactors = FALSE
    )

    idx <- which(keep)
    if (length(idx) == 0L) {
      pass_log[[pass]] <- list(pass = pass, selected = integer(0), reason = "no_components_selected")
      pass_timing[[pass]] <- data.frame(pass = pass, as.list(timing), pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start)
      break
    }

    Tsel <- Tcomp[, idx, drop = FALSE]
    colnames(Tsel) <- sprintf("pass%d_%s%d", pass, extractor, candidate[idx])
    all_selected[[pass]] <- Tsel

    t0 <- proc.time()[["elapsed"]]
    # Project the union of all regressors out of the original data, so the
    # result is one orthogonal projection however many passes ran.
    Xcur <- project_out_components(X0, do.call(cbind, all_selected))
    timing[["regress_seconds"]] <- proc.time()[["elapsed"]] - t0

    pass_log[[pass]] <- list(pass = pass, selected = candidate[idx], n_selected = length(idx))
    pass_timing[[pass]] <- data.frame(pass = pass, as.list(timing), pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start)
  }

  regressors <- if (length(all_selected)) do.call(cbind, all_selected) else matrix(0, nrow = T, ncol = 0)
  component_table <- if (length(all_component_tables)) {
    do.call(rbind, all_component_tables)
  } else {
    data.frame(
      pass = integer(0), candidate = integer(0), ratio_nn_nt = numeric(0),
      med_nn = numeric(0), med_nt = numeric(0), predictability = numeric(0),
      stability = numeric(0), selected = logical(0)
    )
  }

  out <- list(
    x_clean = restore_output(Xcur + mu, parsed, active),
    wNN = expand_voxel_vector(wNN, active),
    regressors = regressors,
    component_table = component_table,
    active = active,
    tr = tr,
    params = list(
      extraction_weight = extraction_weight,
      extractor = extractor,
      lag_order = lag_order,
      pca_rank = pca_rank,
      rank_used = rank_use,
      n_candidates = n_candidates,
      ratio_thresh = ratio_thresh,
      pred_thresh = pred_thresh,
      stability = stability,
      stability_thresh = stability_thresh,
      max_passes = max_passes,
      wnn_threshold = if (is.null(wnn)) wnn_threshold else "given",
      svd_engine = svd_engine,
      use_cpp = use_cpp
    )
  )

  if (return_diagnostics) {
    out$diagnostics <- list(
      input = parsed$info[c("source_type", "transposed_input", "array_input", "original_dim", "n_voxels", "n_time")],
      voxels = list(
        n_voxels = length(active),
        n_masked_out = sum(!parsed$keep),
        n_constant = vox$n_constant,
        n_nonfinite = vox$n_nonfinite,
        n_analysed = N
      ),
      wnn = wres,
      pass_log = pass_log,
      n_passes_executed = length(pass_log),
      timings = list(
        wnn_seconds = wnn_seconds,
        pass_table = if (length(pass_timing)) do.call(rbind, pass_timing) else data.frame(),
        total_seconds = proc.time()[["elapsed"]] - t0_total
      )
    )
  }

  class(out) <- c("fast_phy_denoise_result", "list")
  out
}

# Map a user-supplied wNN (one value per voxel, or per masked voxel) onto the
# analysed voxels, with validation.
wnn_for_active_voxels <- function(wnn, parsed, active) {
  if (!is.numeric(wnn)) stop("`wnn` must be numeric.", call. = FALSE)
  wnn <- as.numeric(wnn)
  n_full <- length(parsed$keep)
  n_keep <- sum(parsed$keep)
  w <- if (length(wnn) == n_full) {
    wnn[active]
  } else if (length(wnn) == n_keep) {
    wnn[active[parsed$keep]]
  } else {
    stop(sprintf(
      "`wnn` must have one value per voxel (%d) or per masked voxel (%d); got %d.",
      n_full, n_keep, length(wnn)
    ), call. = FALSE)
  }
  if (any(!is.finite(w)) || any(w < 0 | w > 1)) {
    stop("`wnn` values for analysed voxels must be finite and in [0, 1].", call. = FALSE)
  }
  w
}

#' Print a fast_phy_denoise result
#'
#' @param x A `fast_phy_denoise_result` from [fast_phy_denoise()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @method print fast_phy_denoise_result
#' @export
print.fast_phy_denoise_result <- function(x, ...) {
  w <- x$wNN[x$active]
  n_cand <- nrow(x$component_table)
  cat("<fast_phy_denoise_result>\n")
  cat(sprintf(
    "  data:       %s (%d of %d voxels analysed, %d time points)\n",
    paste(dim(x$x_clean), collapse = " x "), sum(x$active), length(x$active), nrow(x$regressors)
  ))
  cat(sprintf(
    "  wNN:        %d non-neuronal (< 0.5), %d neuronal (> 0.5) voxels\n",
    sum(w < 0.5), sum(w > 0.5)
  ))
  cat(sprintf(
    "  extractor:  %s, rank %d\n",
    x$params$extractor, x$params$rank_used
  ))
  cat(sprintf(
    "  removed:    %d of %d candidate components\n",
    ncol(x$regressors), n_cand
  ))
  invisible(x)
}
