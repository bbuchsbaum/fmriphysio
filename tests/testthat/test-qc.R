qc_fixture <- function(seed = 30) {
  sim <- simulate_phy_data(n_vox = 400, n_time = 120, tr = 2, baseline = 100, seed = seed)
  fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
  list(sim = sim, fit = fit, w = ifelse(sim$nn, 0, 1))
}

test_that("tSNR after denoising uses the raw mean, so noise removal raises it", {
  f <- qc_fixture()
  qc <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w)
  expect_gt(qc$tsnr_nt_before, 20)
  expect_gt(qc$tsnr_nt_after, qc$tsnr_nt_before)
  expect_equal(qc$tsnr_nt_delta, qc$tsnr_nt_after - qc$tsnr_nt_before)
  expect_gt(qc$nn_variance_reduction_frac, 0.5)
  expect_gt(qc$nt_variance_reduction_frac, 0)
  expect_lt(qc$nt_variance_reduction_frac, qc$nn_variance_reduction_frac)

  # A denoiser that returns mean-removed data gets the same tSNR, not ~0.
  centered <- f$fit$x_clean - rowMeans(f$fit$x_clean)
  qc_c <- compute_design_free_qc(f$sim$x, centered, f$w)
  expect_equal(qc_c$tsnr_nt_after, qc$tsnr_nt_after, tolerance = 1e-10)
})

test_that("identity denoising reports no change", {
  f <- qc_fixture()
  qc <- compute_design_free_qc(f$sim$x, f$sim$x, f$w)
  expect_equal(qc$tsnr_nt_delta, 0)
  expect_equal(qc$nn_variance_reduction_frac, 0)
  expect_equal(qc$nt_variance_reduction_frac, 0)
  expect_identical(qc$n_nt + qc$n_nn, 400L)
  expect_identical(qc$n_unscored, 0L)
})

test_that("constant voxels give NA metrics instead of huge tSNR or full reduction", {
  set.seed(31)
  X <- matrix(100, 60, 50)
  X[11:60, ] <- X[11:60, ] + matrix(rnorm(50 * 50), 50)
  w <- c(rep(0, 10), rep(1, 50))
  expect_no_warning(qc <- compute_design_free_qc(X, X, w, orientation = "voxels_by_time"))
  expect_true(is.na(qc$tsnr_nn_before))
  expect_true(is.na(qc$tsnr_nn_after))
  expect_true(is.na(qc$nn_variance_reduction_frac))
  expect_true(is.finite(qc$tsnr_nt_before))
})

test_that("voxels with NA weight or non-finite data are not scored", {
  f <- qc_fixture()
  w <- f$w
  w[1:5] <- NA
  xr <- f$sim$x
  xr[6, 3] <- NA
  qc <- compute_design_free_qc(xr, f$fit$x_clean, w)
  expect_identical(qc$n_unscored, 6L)
  expect_true(is.finite(qc$tsnr_nt_after))
})

test_that("time x voxels input gives the same metrics", {
  f <- qc_fixture()
  qc <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w)
  qc_t <- compute_design_free_qc(t(f$sim$x), t(f$fit$x_clean), f$w, orientation = "time_by_voxels")
  expect_equal(qc_t, qc)
})

test_that("component summaries handle selections and all-missing ratios", {
  f <- qc_fixture()
  tab <- data.frame(ratio_nn_nt = c(0.5, 2, 4), selected = c(FALSE, TRUE, TRUE))
  qc <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w, component_table = tab)
  expect_identical(qc$selected_components, 2L)
  expect_equal(qc$ratio_summary, list(min = 0.5, median = 2, mean = 6.5 / 3, max = 4))

  tab_na <- data.frame(ratio_nn_nt = c(NA_real_, NA_real_), selected = c(FALSE, FALSE))
  expect_no_warning(
    qc_na <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w, component_table = tab_na)
  )
  expect_identical(qc_na$selected_components, 0L)
  expect_true(all(is.na(unlist(qc_na$ratio_summary))))

  qc_none <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w)
  expect_true(is.na(qc_none$selected_components))
  expect_null(qc_none$ratio_summary)
})

test_that("compute_design_free_qc validates its inputs", {
  f <- qc_fixture()
  expect_error(compute_design_free_qc(f$sim$x, f$sim$x[-1, ], f$w), "matching dimensions")
  expect_error(compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w[-1]), "one value per voxel")
  expect_error(compute_design_free_qc(f$sim$x, f$fit$x_clean, rep(1, 400)), "Need both")
})

test_that("QC artifacts round-trip through RDS and CSV", {
  f <- qc_fixture()
  tab <- data.frame(ratio_nn_nt = c(0.5, 2), selected = c(FALSE, TRUE))
  qc <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w, component_table = tab)

  rds <- tempfile(fileext = ".rds")
  expect_identical(write_qc_artifact(qc, rds), rds)
  expect_identical(readRDS(rds), qc)

  csv <- tempfile(fileext = ".CSV")
  write_qc_artifact(qc, csv)
  df <- utils::read.csv(csv, stringsAsFactors = FALSE)
  expect_true(all(c("ratio_summary.min", "ratio_summary.max", "tsnr_nt_after") %in% df$metric))
  val <- function(m) as.numeric(df$value[df$metric == m])
  expect_equal(val("tsnr_nt_after"), qc$tsnr_nt_after, tolerance = 1e-12)
  expect_equal(val("ratio_summary.median"), qc$ratio_summary$median, tolerance = 1e-12)
})

test_that("JSON artifacts write missing and non-finite values as null", {
  skip_if_not_installed("jsonlite")
  f <- qc_fixture()
  tab_na <- data.frame(ratio_nn_nt = NA_real_, selected = FALSE)
  qc <- compute_design_free_qc(f$sim$x, f$fit$x_clean, f$w, component_table = tab_na)
  qc$inf_metric <- Inf
  json <- tempfile(fileext = ".json")
  write_qc_artifact(qc, json)
  txt <- paste(readLines(json), collapse = "\n")
  expect_false(grepl("\"NA\"|\"Inf\"|\"-Inf\"|\"NaN\"", txt))
  expect_match(txt, "\"min\": null")
  expect_match(txt, "\"inf_metric\": null")
  back <- jsonlite::fromJSON(json)
  expect_equal(back$tsnr_nt_after, qc$tsnr_nt_after, tolerance = 1e-12)
  expect_equal(back$n_vox, 400L)
})

test_that("write_qc_artifact rejects bad inputs", {
  expect_error(write_qc_artifact(list(a = 1), tempfile(fileext = ".txt")), "Unsupported")
  expect_error(write_qc_artifact(1:3, tempfile(fileext = ".rds")), "must be a list")
  expect_error(write_qc_artifact(list(a = 1), c("a.rds", "b.rds")), "single file path")
})
