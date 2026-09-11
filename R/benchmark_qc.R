# Benchmark and QC utilities.

#' Benchmark fast_phy_denoise on simulated data
#'
#' Times [fast_phy_denoise()] over a grid of data sizes and SVD engines, and
#' scores each run against the ground truth of [simulate_phy_data()]. Each
#' simulated dataset is generated once and reused for every engine, so engines
#' are compared on identical data.
#'
#' @param grid_n_vox Integer vector of voxel counts.
#' @param grid_n_time Integer vector of time-series lengths.
#' @param engines Character vector of SVD engines to compare (any of
#'   `"auto"`, `"svd"`, `"rsvd"`).
#' @param reps Number of simulated datasets per size (at least 1).
#' @param tr Repetition time in seconds used to simulate and denoise.
#' @param extractor Component extractor passed to [fast_phy_denoise()].
#' @param out_file Optional path; results are also written there as CSV.
#' @param seed Base seed. Dataset `i` is simulated with seed `seed + i`, and
#'   the same seed is passed to [fast_phy_denoise()]. The caller's
#'   random-number stream is left untouched. Use `NULL` for unseeded runs.
#' @param ... Further arguments to [fast_phy_denoise()]. `x`, `orientation`,
#'   `svd_engine`, and `return_diagnostics` are set by the benchmark and may
#'   not be given.
#' @return A data frame with one row per dataset and engine and columns
#'   `n_vox`, `n_time`, `rep`, `dataset_seed`, `svd_engine`, `status`
#'   (`"ok"` or `"error"`), `total_seconds`, `wnn_seconds`,
#'   `pass_total_seconds`, `selected_components`, `physio_kept_nt` and
#'   `neural_kept_nt` (median over neuronal voxels of the fraction of
#'   physiological / neural source energy remaining after denoising), and
#'   `error_message`.
#' @seealso [simulate_phy_data()], [fast_phy_denoise()]
#' @examples
#' res <- benchmark_fast_phy_denoise(
#'   grid_n_vox = 300, grid_n_time = 100, engines = "svd", reps = 1
#' )
#' res[, c("svd_engine", "status", "selected_components", "neural_kept_nt")]
#' @export
benchmark_fast_phy_denoise <- function(
  grid_n_vox = 4000L,
  grid_n_time = 400L,
  engines = c("svd", "rsvd"),
  reps = 1L,
  tr = 2,
  extractor = c("caa", "dicca", "dipca"),
  out_file = NULL,
  seed = 1L,
  ...
) {
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All arguments passed through `...` must be named.", call. = FALSE)
  }
  reserved <- intersect(names(dots), c("x", "orientation", "svd_engine", "return_diagnostics"))
  if (length(reserved)) {
    stop(sprintf(
      "`%s` is set by the benchmark and cannot be passed through `...`.",
      reserved[1]
    ), call. = FALSE)
  }
  grid_n_vox <- vapply(grid_n_vox, check_count, integer(1), name = "grid_n_vox", min = 20L)
  grid_n_time <- vapply(grid_n_time, check_count, integer(1), name = "grid_n_time", min = 20L)
  if (!is.character(engines) || !length(engines) ||
      !all(engines %in% c("auto", "svd", "rsvd"))) {
    stop("`engines` must contain only \"auto\", \"svd\", or \"rsvd\".", call. = FALSE)
  }
  reps <- check_count(reps, "reps")
  tr <- check_positive(tr, "tr")
  extractor <- match.arg(extractor)
  seed <- check_number(seed, "seed", null_ok = TRUE)

  rows <- list()
  dataset <- 0L
  for (n_vox in grid_n_vox) {
    for (n_time in grid_n_time) {
      for (rep in seq_len(reps)) {
        dataset <- dataset + 1L
        ds_seed <- if (is.null(seed)) NULL else seed + dataset
        sim <- simulate_phy_data(n_vox, n_time, tr = tr, seed = ds_seed)
        nt <- !sim$nn

        for (eng in engines) {
          fit <- tryCatch(
            do.call(fast_phy_denoise, c(
              list(
                x = sim$x,
                tr = tr,
                orientation = "voxels_by_time",
                extractor = extractor,
                svd_engine = eng,
                seed = ds_seed,
                return_diagnostics = TRUE
              ),
              dots
            )),
            error = function(e) e
          )

          row <- data.frame(
            n_vox = n_vox,
            n_time = n_time,
            rep = rep,
            dataset_seed = if (is.null(ds_seed)) NA_real_ else ds_seed,
            svd_engine = eng,
            status = "ok",
            total_seconds = NA_real_,
            wnn_seconds = NA_real_,
            pass_total_seconds = NA_real_,
            selected_components = NA_integer_,
            physio_kept_nt = NA_real_,
            neural_kept_nt = NA_real_,
            error_message = NA_character_,
            stringsAsFactors = FALSE
          )
          if (inherits(fit, "error")) {
            row$status <- "error"
            row$error_message <- conditionMessage(fit)
          } else {
            timing <- fit$diagnostics$timings
            ptab <- timing$pass_table
            row$total_seconds <- as.numeric(timing$total_seconds)
            row$wnn_seconds <- as.numeric(timing$wnn_seconds)
            row$pass_total_seconds <- if (is.data.frame(ptab) && nrow(ptab)) sum(ptab$pass_total_seconds) else 0
            row$selected_components <- ncol(fit$regressors)
            row$physio_kept_nt <- stats::median(
              source_energy_kept(sim$x[nt, , drop = FALSE], fit$x_clean[nt, , drop = FALSE], sim$physio)
            )
            row$neural_kept_nt <- stats::median(
              source_energy_kept(sim$x[nt, , drop = FALSE], fit$x_clean[nt, , drop = FALSE], sim$neural)
            )
          }
          rows[[length(rows) + 1L]] <- row
        }
      }
    }
  }

  res <- do.call(rbind, rows)
  if (!is.null(out_file)) {
    utils::write.csv(res, out_file, row.names = FALSE)
  }
  res
}

# Per-voxel fraction of the energy in the span of the source time courses
# `S` (time x k) that survives from `X_raw` to `X_clean` (both voxels x time).
source_energy_kept <- function(X_raw, X_clean, S) {
  S <- scale(as.matrix(S), center = TRUE, scale = FALSE)
  Q <- qr.Q(qr(S))
  before <- rowSums((center_rows(X_raw) %*% Q)^2)
  after <- rowSums((center_rows(X_clean) %*% Q)^2)
  ifelse(before > 0, after / before, NA_real_)
}

#' Compute design-free QC metrics
#'
#' Summarises how denoising changed the data without reference to any task
#' design, separately for neuronal (wNN > 0.5) and non-neuronal (wNN < 0.5)
#' voxels.
#'
#' Temporal SNR is `|mean| / sd` per voxel. Because denoisers may alter
#' voxel means, the "after" tSNR uses the *raw* voxel mean over the cleaned
#' standard deviation, so it measures noise reduction on a fixed signal level.
#' Voxels whose standard deviation is zero (relative to their mean) have
#' undefined tSNR and are skipped; a metric with no defined voxels is `NA`.
#'
#' @param x_raw Raw data, in any form accepted by [fast_phy_denoise()].
#' @param x_clean Denoised data with the same shape as `x_raw`.
#' @param wNN Voxel weights (1 = neuronal, 0 = non-neuronal), one per voxel,
#'   e.g. `fit$wNN` from [fast_phy_denoise()]. Voxels with `NA` weight, or
#'   with non-finite data, are not scored.
#' @param component_table Optional component table from
#'   [fast_phy_denoise()], summarised in `selected_components` and
#'   `ratio_summary`.
#' @param orientation Orientation of `x_raw`: `"auto"`, `"voxels_by_time"`,
#'   or `"time_by_voxels"`. `x_clean` is read in the same orientation.
#' @return A named list: `n_vox`, `n_time`, `n_nt`, `n_nn` (voxels scored in
#'   each class), `n_unscored`; `tsnr_nt_before`, `tsnr_nt_after`,
#'   `tsnr_nt_delta`, `tsnr_nn_before`, `tsnr_nn_after` (median tSNR);
#'   `nt_variance_reduction_frac`, `nn_variance_before`, `nn_variance_after`,
#'   `nn_variance_reduction_frac` (median variance and the fractional
#'   reduction of the median); `selected_components` (`NA` without a
#'   component table); and `ratio_summary`, a list with the `min`, `median`,
#'   `mean` and `max` of the finite NN/NT ratios (`NA` when there are none), or
#'   `NULL` without a component table.
#' @seealso [write_qc_artifact()]
#' @examples
#' sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
#' fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
#' qc <- compute_design_free_qc(sim$x, fit$x_clean, wNN = ifelse(sim$nn, 0, 1))
#' qc$tsnr_nt_before
#' qc$tsnr_nt_after
#' @export
compute_design_free_qc <- function(
  x_raw,
  x_clean,
  wNN,
  component_table = NULL,
  orientation = c("auto", "voxels_by_time", "time_by_voxels")
) {
  orientation <- match.arg(orientation)
  raw_p <- as_nt_matrix(x_raw, orientation = orientation)
  clean_orientation <- if (isTRUE(raw_p$info$transposed_input)) "time_by_voxels" else "voxels_by_time"
  cln_p <- as_nt_matrix(x_clean, orientation = clean_orientation)
  Xr <- raw_p$X
  Xc <- cln_p$X

  if (!identical(dim(Xr), dim(Xc))) {
    stop("x_raw and x_clean must have matching dimensions.", call. = FALSE)
  }
  if (!is.numeric(wNN) || length(wNN) != nrow(Xr)) {
    stop(sprintf("`wNN` must be numeric with one value per voxel (%d).", nrow(Xr)), call. = FALSE)
  }
  wNN <- as.numeric(wNN)

  finite_rows <- is.finite(rowSums(Xr)) & is.finite(rowSums(Xc))
  scored <- finite_rows & !is.na(wNN)
  nn_idx <- which(scored & wNN < 0.5)
  nt_idx <- which(scored & wNN > 0.5)
  if (length(nn_idx) == 0L || length(nt_idx) == 0L) {
    stop("Need both neuronal (wNN > 0.5) and non-neuronal (wNN < 0.5) voxels for QC.", call. = FALSE)
  }

  T <- ncol(Xr)
  row_sd <- function(X) sqrt(rowSums((X - rowMeans(X))^2) / (T - 1))
  mu <- rowMeans(Xr)
  sd_r <- row_sd(Xr)
  sd_c <- row_sd(Xc)
  floor_sd <- 1e-10 * pmax(abs(mu), 1)
  tsnr_r <- ifelse(sd_r > floor_sd, abs(mu) / sd_r, NA_real_)
  tsnr_c <- ifelse(sd_c > floor_sd, abs(mu) / sd_c, NA_real_)
  var_r <- sd_r^2
  var_c <- sd_c^2

  med <- function(v) {
    v <- v[is.finite(v)]
    if (length(v)) stats::median(v) else NA_real_
  }
  reduction <- function(idx) {
    before <- med(var_r[idx])
    if (!is.finite(before) || before <= 0) {
      return(NA_real_)
    }
    1 - med(var_c[idx]) / before
  }

  selected_n <- NA_integer_
  ratio_summary <- NULL
  if (!is.null(component_table) && nrow(component_table) > 0L) {
    if ("selected" %in% names(component_table)) {
      selected_n <- sum(component_table$selected %in% TRUE)
    }
    r <- component_table$ratio_nn_nt
    r <- r[is.finite(r)]
    ratio_summary <- if (length(r)) {
      list(min = min(r), median = stats::median(r), mean = mean(r), max = max(r))
    } else {
      list(min = NA_real_, median = NA_real_, mean = NA_real_, max = NA_real_)
    }
  }

  tsnr_nt_before <- med(tsnr_r[nt_idx])
  tsnr_nt_after <- med(tsnr_c[nt_idx])
  list(
    n_vox = nrow(Xr),
    n_time = T,
    n_nt = length(nt_idx),
    n_nn = length(nn_idx),
    n_unscored = nrow(Xr) - length(nt_idx) - length(nn_idx),
    tsnr_nt_before = tsnr_nt_before,
    tsnr_nt_after = tsnr_nt_after,
    tsnr_nt_delta = tsnr_nt_after - tsnr_nt_before,
    tsnr_nn_before = med(tsnr_r[nn_idx]),
    tsnr_nn_after = med(tsnr_c[nn_idx]),
    nt_variance_reduction_frac = reduction(nt_idx),
    nn_variance_before = med(var_r[nn_idx]),
    nn_variance_after = med(var_c[nn_idx]),
    nn_variance_reduction_frac = reduction(nn_idx),
    selected_components = selected_n,
    ratio_summary = ratio_summary
  )
}

#' Write a QC summary to disk
#'
#' Saves the list returned by [compute_design_free_qc()] as an R object
#' (`.rds`), a two-column `metric,value` table (`.csv`, nested entries such
#' as `ratio_summary` flattened to `ratio_summary.min` etc.), or JSON
#' (`.json`, requires the 'jsonlite' package; missing and non-finite values
#' are written as `null`).
#'
#' @param qc QC list from [compute_design_free_qc()].
#' @param file Output path; the format is chosen from its extension (`.rds`,
#'   `.csv`, or `.json`, case-insensitive).
#' @return `file`, invisibly.
#' @examples
#' sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 1)
#' fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
#' qc <- compute_design_free_qc(sim$x, fit$x_clean, wNN = ifelse(sim$nn, 0, 1))
#' path <- tempfile(fileext = ".csv")
#' write_qc_artifact(qc, path)
#' head(read.csv(path))
#' @export
write_qc_artifact <- function(qc, file) {
  if (!is.list(qc)) {
    stop("`qc` must be a list, as returned by compute_design_free_qc().", call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a single file path.", call. = FALSE)
  }

  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    saveRDS(qc, file = file)
    return(invisible(file))
  }

  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    flat <- unlist(qc, use.names = TRUE)
    df <- data.frame(
      metric = names(flat),
      value = format(flat, digits = 15, trim = TRUE, scientific = FALSE),
      stringsAsFactors = FALSE
    )
    df$value[is.na(flat)] <- NA_character_
    utils::write.csv(df, file, row.names = FALSE)
    return(invisible(file))
  }

  if (grepl("\\.json$", file, ignore.case = TRUE)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Writing a JSON QC artifact requires the 'jsonlite' package.", call. = FALSE)
    }
    txt <- jsonlite::toJSON(
      qc,
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null",
      na = "null",
      digits = NA
    )
    writeLines(txt, con = file)
    return(invisible(file))
  }

  stop("Unsupported artifact extension. Use .rds, .csv, or .json.", call. = FALSE)
}
