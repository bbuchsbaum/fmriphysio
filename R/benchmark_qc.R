# Benchmark and QC utilities for fast PHY denoising.

simulate_phy_data <- function(n_vox, n_time, nn_frac = 0.20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_vox <- as.integer(n_vox)
  n_time <- as.integer(n_time)
  if (n_vox < 20L || n_time < 20L) stop("simulate_phy_data requires at least 20 voxels and 20 time points.")

  idx <- seq_len(n_vox)
  n_nn <- max(1L, floor(nn_frac * n_vox))
  nn_idx <- idx[seq_len(n_nn)]

  t <- seq(0, 1, length.out = n_time)
  phys1 <- sin(2 * pi * 9 * t)
  phys2 <- cos(2 * pi * 13 * t + 0.2)
  neural1 <- sin(2 * pi * 2 * t + 0.3)
  neural2 <- cos(2 * pi * 3 * t - 0.1)

  Wp <- matrix(rnorm(n_vox * 2, sd = 0.6), n_vox, 2)
  Wn <- matrix(rnorm(n_vox * 2, sd = 0.3), n_vox, 2)
  Wp[-nn_idx, ] <- Wp[-nn_idx, ] * 0.25
  Wn[nn_idx, ] <- Wn[nn_idx, ] * 0.30

  X <- Wp %*% rbind(phys1, phys2) + Wn %*% rbind(neural1, neural2)
  X <- X + matrix(rnorm(n_vox * n_time, sd = 0.8), n_vox, n_time)
  X
}

#' Benchmark fast_phy_denoise across dimensions and SVD engines
#'
#' @param grid_n_vox Integer vector of voxel counts.
#' @param grid_n_time Integer vector of time lengths.
#' @param svd_engines Character vector of SVD engines to compare.
#' @param reps Number of repetitions per configuration.
#' @param out_file Optional CSV path for machine-readable results.
#' @param seed Seed base.
#' @param ... Additional arguments passed to fast_phy_denoise.
#' @return Data frame of benchmark results.
#' @export
benchmark_fast_phy_denoise <- function(
  grid_n_vox = c(4000L),
  grid_n_time = c(400L),
  svd_engines = c("svd", "rsvd"),
  reps = 1L,
  out_file = NULL,
  seed = 1L,
  ...
) {
  reps <- as.integer(reps)
  run <- 0L
  rows <- vector("list", length(grid_n_vox) * length(grid_n_time) * length(svd_engines) * reps)

  for (n_vox in grid_n_vox) {
    for (n_time in grid_n_time) {
      for (eng in svd_engines) {
        for (rep in seq_len(reps)) {
          run <- run + 1L
          set.seed(as.integer(seed) + run)
          X <- simulate_phy_data(as.integer(n_vox), as.integer(n_time))

          fit <- tryCatch(
            fast_phy_denoise(
              x = X,
              tr = 0.8,
              extractor = "caa",
              svd_engine = eng,
              return_diagnostics = TRUE,
              ...
            ),
            error = function(e) e
          )

          if (inherits(fit, "error")) {
            rows[[run]] <- data.frame(
              n_vox = as.integer(n_vox),
              n_time = as.integer(n_time),
              svd_engine = eng,
              rep = rep,
              status = "error",
              total_seconds = NA_real_,
              wnn_seconds = NA_real_,
              pass_total_seconds = NA_real_,
              selected_components = NA_integer_,
              error_message = conditionMessage(fit),
              stringsAsFactors = FALSE
            )
          } else {
            timing <- fit$diagnostics$timings
            ptab <- timing$pass_table
            rows[[run]] <- data.frame(
              n_vox = as.integer(n_vox),
              n_time = as.integer(n_time),
              svd_engine = eng,
              rep = rep,
              status = "ok",
              total_seconds = as.numeric(timing$total_seconds),
              wnn_seconds = as.numeric(timing$wnn_seconds),
              pass_total_seconds = if (nrow(ptab)) sum(ptab$pass_total_seconds) else 0,
              selected_components = ncol(fit$components_selected),
              error_message = NA_character_,
              stringsAsFactors = FALSE
            )
          }
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

#' Compute design-free QC summary
#'
#' @param x_raw Raw data (matrix-like).
#' @param x_clean Denoised data (matrix-like).
#' @param wNN Non-neuronal weighting vector.
#' @param component_table Optional component summary table from fast_phy_denoise.
#' @return Named list of QC metrics.
#' @export
compute_design_free_qc <- function(x_raw, x_clean, wNN, component_table = NULL) {
  raw_p <- as_nt_matrix(x_raw)
  cln_p <- as_nt_matrix(x_clean)
  Xr <- raw_p$X
  Xc <- cln_p$X

  if (!all(dim(Xr) == dim(Xc))) {
    stop("x_raw and x_clean must have matching dimensions after orientation parsing.")
  }
  if (length(wNN) != nrow(Xr)) {
    stop("wNN length must match voxel dimension.")
  }

  nn_idx <- which(wNN < 0.5)
  nt_idx <- which(wNN > 0.5)
  if (length(nn_idx) == 0L || length(nt_idx) == 0L) {
    stop("Need both NN and NT voxels for QC.")
  }

  tsnr <- function(X) {
    mu <- rowMeans(X)
    s <- apply(X, 1, stats::sd)
    mu / pmax(s, 1e-8)
  }
  row_var <- function(X) apply(X, 1, stats::var)

  tsnr_r <- tsnr(Xr)
  tsnr_c <- tsnr(Xc)
  var_r <- row_var(Xr)
  var_c <- row_var(Xc)

  selected_n <- NA_integer_
  ratio_summary <- NULL
  if (!is.null(component_table) && nrow(component_table) > 0L) {
    selected_n <- sum(component_table$selected %in% TRUE)
    ratio_summary <- list(
      min = min(component_table$ratio_nn_nt, na.rm = TRUE),
      median = stats::median(component_table$ratio_nn_nt, na.rm = TRUE),
      mean = mean(component_table$ratio_nn_nt, na.rm = TRUE),
      max = max(component_table$ratio_nn_nt, na.rm = TRUE)
    )
  }

  list(
    n_vox = nrow(Xr),
    n_time = ncol(Xr),
    n_nt = length(nt_idx),
    n_nn = length(nn_idx),
    tsnr_nt_before = stats::median(tsnr_r[nt_idx], na.rm = TRUE),
    tsnr_nt_after = stats::median(tsnr_c[nt_idx], na.rm = TRUE),
    tsnr_nn_before = stats::median(tsnr_r[nn_idx], na.rm = TRUE),
    tsnr_nn_after = stats::median(tsnr_c[nn_idx], na.rm = TRUE),
    tsnr_nt_delta = stats::median(tsnr_c[nt_idx], na.rm = TRUE) - stats::median(tsnr_r[nt_idx], na.rm = TRUE),
    nn_variance_before = stats::median(var_r[nn_idx], na.rm = TRUE),
    nn_variance_after = stats::median(var_c[nn_idx], na.rm = TRUE),
    nn_variance_reduction_frac = 1 - stats::median(var_c[nn_idx], na.rm = TRUE) / pmax(stats::median(var_r[nn_idx], na.rm = TRUE), 1e-8),
    selected_components = selected_n,
    ratio_summary = ratio_summary
  )
}

#' Write QC artifact to file
#'
#' @param qc QC list from compute_design_free_qc().
#' @param file Path ending in .rds, .csv, or .json.
#' @return Invisibly returns file path.
#' @export
write_qc_artifact <- function(qc, file) {
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    saveRDS(qc, file = file)
    return(invisible(file))
  }

  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    flat <- unlist(qc[!vapply(qc, is.list, logical(1))], use.names = TRUE)
    df <- data.frame(metric = names(flat), value = as.character(flat), stringsAsFactors = FALSE)
    utils::write.csv(df, file, row.names = FALSE)
    return(invisible(file))
  }

  if (grepl("\\.json$", file, ignore.case = TRUE)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Writing JSON QC artifact requires package 'jsonlite'.")
    }
    txt <- jsonlite::toJSON(qc, auto_unbox = TRUE, pretty = TRUE, null = "null")
    writeLines(txt, con = file)
    return(invisible(file))
  }

  stop("Unsupported artifact extension. Use .rds, .csv, or .json.")
}
