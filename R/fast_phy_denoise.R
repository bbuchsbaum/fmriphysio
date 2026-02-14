#' Fast design-free PHYCAA+-style denoising
#'
#' @param x Matrix-like data. Default internal orientation is voxels x time.
#' @param tr Repetition time in seconds.
#' @param design Optional design matrix (time x predictors) used as guardrail.
#' @param mask Optional voxel mask (logical or indices).
#' @param wnn Optional precomputed wNN vector.
#' @param wnn_method Currently supports "diff_energy".
#' @param wnn_threshold "percentile" or "mixture".
#' @param delta_nt Percentile for neuronal cutoff.
#' @param delta_nn Percentile for non-neuronal cutoff.
#' @param extractor One of "caa", "dicca", or "dipca".
#' @param lag_order Dynamic lag order.
#' @param pca_rank Low-rank dimension.
#' @param n_candidates Max candidate components.
#' @param ratio_thresh NN-vs-NT median R2 ratio threshold.
#' @param pred_thresh Optional predictability threshold.
#' @param stability "none" or "half" split-half filter.
#' @param stability_thresh Correlation threshold for split-half stability.
#' @param max_passes Number of denoising passes (1-2 recommended).
#' @param return_diagnostics Whether to return diagnostics.
#' @param svd_engine "auto", "rsvd", or "svd".
#' @param use_cpp Use compiled Rcpp kernels when available.
#' @param seed Optional random seed.
#' @return A list with cleaned data, selected regressors, and diagnostics.
#' @export
fast_phy_denoise <- function(
  x,
  tr,
  design = NULL,
  mask = NULL,
  wnn = NULL,
  wnn_method = c("diff_energy"),
  wnn_threshold = c("percentile", "mixture"),
  delta_nt = 0.74,
  delta_nn = 0.95,
  extractor = c("caa", "dicca", "dipca"),
  lag_order = 1L,
  pca_rank = 100L,
  n_candidates = 30L,
  ratio_thresh = 1.0,
  pred_thresh = NULL,
  stability = c("none", "half"),
  stability_thresh = 0.30,
  max_passes = 1L,
  return_diagnostics = TRUE,
  svd_engine = c("auto", "rsvd", "svd"),
  use_cpp = TRUE,
  seed = NULL
) {
  t0_total <- proc.time()[["elapsed"]]

  if (missing(tr) || !is.numeric(tr) || length(tr) != 1L || tr <= 0) {
    stop("tr must be a positive scalar.")
  }

  wnn_method <- match.arg(wnn_method)
  wnn_threshold <- match.arg(wnn_threshold)
  extractor <- match.arg(extractor)
  stability <- match.arg(stability)
  svd_engine <- match.arg(svd_engine)

  parsed <- as_nt_matrix(x, mask = mask)
  X0 <- parsed$X
  N <- nrow(X0)
  T <- ncol(X0)

  if (T < 8L) stop("Need at least 8 time points.")
  if (N < 20L) stop("Need at least 20 voxels.")

  X0 <- center_rows(X0)
  Xcur <- X0

  t0_wnn <- proc.time()[["elapsed"]]
  if (is.null(wnn)) {
    if (wnn_method != "diff_energy") {
      stop("Only wnn_method='diff_energy' is currently implemented.")
    }
    wres <- compute_wNN_diff_energy(
      X_nt = Xcur,
      nt_q = delta_nt,
      nn_q = delta_nn,
      threshold = wnn_threshold
    )
    wNN <- wres$wNN
  } else {
    wNN <- as.numeric(wnn)
    if (length(wNN) != N) {
      stop("Length of wnn must match voxel dimension.")
    }
    wres <- list(
      wNN = wNN,
      delta = NA_real_,
      energy = NA_real_,
      delta_nt = NA_real_,
      delta_nn = NA_real_,
      threshold_method = "given"
    )
  }
  wnn_seconds <- proc.time()[["elapsed"]] - t0_wnn

  pass_log <- vector("list", max_passes)
  pass_timing <- vector("list", max_passes)
  all_selected <- list()
  all_component_tables <- list()

  for (pass in seq_len(max_passes)) {
    t_pass_start <- proc.time()[["elapsed"]]
    Xw <- Xcur * wNN
    rank_use <- min(as.integer(pca_rank), min(dim(Xw)))

    t0 <- proc.time()[["elapsed"]]
    sv <- truncated_svd(Xw, rank = rank_use, svd_engine = svd_engine, seed = seed)
    Z <- t(sv$v)
    svd_seconds <- proc.time()[["elapsed"]] - t0

    t0 <- proc.time()[["elapsed"]]
    ext <- extract_dynamic_components(
      Z_kt = Z,
      extractor = extractor,
      lag_order = lag_order,
      n_candidates = n_candidates
    )
    Tcomp <- as.matrix(ext$scores)
    extract_seconds <- proc.time()[["elapsed"]] - t0

    # Optional design guardrail: remove design-explainable subspace from nuisances.
    if (!is.null(design)) {
      Tcomp <- orthogonalize_to_design(Tcomp, design)
    }

    if (ncol(Tcomp) == 0L) {
      pass_log[[pass]] <- list(pass = pass, selected = integer(0), reason = "no_components_after_guardrail")
      pass_timing[[pass]] <- data.frame(
        pass = pass,
        svd_seconds = svd_seconds,
        extract_seconds = extract_seconds,
        score_seconds = 0,
        stability_seconds = 0,
        regress_seconds = 0,
        pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start
      )
      break
    }

    t0 <- proc.time()[["elapsed"]]
    # Score NN-vs-NT enrichment on unweighted data to avoid suppressing NN variance by construction.
    ss <- score_components_nn_nt(X_nt = Xcur, T_tl = Tcomp, wNN = wNN, use_cpp = use_cpp)
    score_seconds <- proc.time()[["elapsed"]] - t0
    pred <- ext$predictability
    if (length(pred) < ncol(Tcomp)) {
      pred <- c(pred, rep(NA_real_, ncol(Tcomp) - length(pred)))
    }
    pred <- pred[seq_len(ncol(Tcomp))]

    keep <- ss$ratio > ratio_thresh
    if (!is.null(pred_thresh)) {
      keep <- keep & (pred >= pred_thresh)
    }

    split_scores <- NULL
    stability_seconds <- 0
    if (stability == "half" && any(keep)) {
      t0 <- proc.time()[["elapsed"]]
      split_scores <- compute_split_stability(
        X_nt = Xw,
        pca_rank = rank_use,
        extractor = extractor,
        lag_order = lag_order,
        n_candidates = n_candidates,
        svd_engine = svd_engine,
        seed = seed
      )
      stable <- stability_filter(Tcomp, split_scores, thresh = stability_thresh, use_cpp = use_cpp)
      keep <- keep & stable
      stability_seconds <- proc.time()[["elapsed"]] - t0
    }

    idx <- which(keep)
    comp_table <- data.frame(
      pass = pass,
      component = seq_len(ncol(Tcomp)),
      ratio_nn_nt = ss$ratio,
      med_nn = ss$med_nn,
      med_nt = ss$med_nt,
      predictability = pred,
      selected = keep,
      stringsAsFactors = FALSE
    )
    all_component_tables[[pass]] <- comp_table

    if (length(idx) == 0L) {
      pass_log[[pass]] <- list(pass = pass, selected = integer(0), reason = "no_components_selected")
      pass_timing[[pass]] <- data.frame(
        pass = pass,
        svd_seconds = svd_seconds,
        extract_seconds = extract_seconds,
        score_seconds = score_seconds,
        stability_seconds = stability_seconds,
        regress_seconds = 0,
        pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start
      )
      break
    }

    Tsel <- Tcomp[, idx, drop = FALSE]
    t0 <- proc.time()[["elapsed"]]
    Xnext <- project_out_components(Xcur, Tsel)
    regress_seconds <- proc.time()[["elapsed"]] - t0

    pass_log[[pass]] <- list(
      pass = pass,
      selected = idx,
      n_selected = length(idx),
      total_ratio_kept = sum(ss$ratio[idx], na.rm = TRUE)
    )
    pass_timing[[pass]] <- data.frame(
      pass = pass,
      svd_seconds = svd_seconds,
      extract_seconds = extract_seconds,
      score_seconds = score_seconds,
      stability_seconds = stability_seconds,
      regress_seconds = regress_seconds,
      pass_total_seconds = proc.time()[["elapsed"]] - t_pass_start
    )

    all_selected[[pass]] <- Tsel
    Xcur <- Xnext
  }

  used_passes <- which(vapply(pass_log, Negate(is.null), logical(1)))
  if (length(used_passes) == 0L) {
    used_passes <- integer(0)
  }

  components_selected <- if (length(all_selected)) do.call(cbind, all_selected) else matrix(0, nrow = T, ncol = 0)
  component_table <- if (length(all_component_tables)) do.call(rbind, all_component_tables) else data.frame()
  Xout <- restore_orientation(Xcur, parsed$info)

  out <- list(
    x_clean = Xout,
    wNN = wNN,
    components_selected = components_selected,
    component_table = component_table,
    tr = tr,
    params = list(
      extractor = extractor,
      lag_order = lag_order,
      pca_rank = pca_rank,
      ratio_thresh = ratio_thresh,
      pred_thresh = pred_thresh,
      stability = stability,
      max_passes = max_passes,
      svd_engine = svd_engine,
      use_cpp = use_cpp
    )
  )

  if (return_diagnostics) {
    used_timing <- pass_timing[used_passes]
    timing_table <- if (length(used_timing)) do.call(rbind, used_timing) else data.frame()
    total_seconds <- proc.time()[["elapsed"]] - t0_total
    out$diagnostics <- list(
      input = parsed$info,
      wnn = wres,
      pass_log = pass_log[used_passes],
      n_passes_executed = length(used_passes),
      timings = list(
        wnn_seconds = wnn_seconds,
        pass_table = timing_table,
        total_seconds = total_seconds
      )
    )
  }

  class(out) <- c("fast_phy_denoise_result", class(out))
  out
}
