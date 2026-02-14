test_that("compute_wNN_diff_energy returns bounded weights", {
  set.seed(1)
  X <- matrix(rnorm(200 * 120), nrow = 200, ncol = 120)
  out <- compute_wNN_diff_energy(X, threshold = "percentile")

  expect_length(out$wNN, nrow(X))
  expect_true(all(out$wNN >= 0 & out$wNN <= 1))
  expect_true(is.finite(out$delta_nt))
  expect_true(is.finite(out$delta_nn))
})

test_that("fast_phy_denoise runs and preserves matrix shape", {
  set.seed(2)
  n_vox <- 250
  n_time <- 160

  s1 <- sin(seq(0, 8 * pi, length.out = n_time))
  s2 <- cos(seq(0, 6 * pi, length.out = n_time))
  S <- cbind(s1, s2)
  W <- matrix(rnorm(n_vox * 2), n_vox, 2)
  X <- W %*% t(S) + matrix(rnorm(n_vox * n_time, sd = 0.8), n_vox, n_time)

  fit <- fast_phy_denoise(
    x = X,
    tr = 0.8,
    extractor = "caa",
    pca_rank = 30,
    n_candidates = 8,
    max_passes = 1,
    stability = "none"
  )

  expect_equal(dim(fit$x_clean), dim(X))
  expect_length(fit$wNN, nrow(X))
  expect_true(is.data.frame(fit$component_table))
})

test_that("fast_phy_denoise handles T x N input and restores orientation", {
  set.seed(3)
  X_tn <- matrix(rnorm(140 * 80), nrow = 140, ncol = 80) # T x N
  fit <- fast_phy_denoise(X_tn, tr = 1.0, extractor = "caa", pca_rank = 20, n_candidates = 5)
  expect_equal(dim(fit$x_clean), dim(X_tn))
})

test_that("fast_phy_denoise selects nonzero components on simulated phy data", {
  set.seed(42)
  X <- simulate_phy_data(n_vox = 500, n_time = 200, nn_frac = 0.20)
  fit <- fast_phy_denoise(
    x = X,
    tr = 2.0,
    extractor = "caa",
    stability = "none",
    max_passes = 1
  )

  expect_true(ncol(fit$components_selected) > 0)
})

test_that("benchmark_fast_phy_denoise returns machine-readable table", {
  set.seed(4)
  out_csv <- tempfile(fileext = ".csv")
  bench <- benchmark_fast_phy_denoise(
    grid_n_vox = c(120L),
    grid_n_time = c(80L),
    svd_engines = c("svd"),
    reps = 1L,
    out_file = out_csv,
    pca_rank = 10L,
    n_candidates = 4L,
    max_passes = 1L,
    stability = "none"
  )

  expect_true(is.data.frame(bench))
  expect_true(nrow(bench) >= 1L)
  expect_true(file.exists(out_csv))
  expect_true(all(c("status", "total_seconds", "svd_engine") %in% names(bench)))
})

test_that("compute_design_free_qc and artifact writer work", {
  set.seed(5)
  X <- matrix(rnorm(180 * 100), nrow = 180, ncol = 100)
  fit <- fast_phy_denoise(X, tr = 0.8, extractor = "caa", pca_rank = 20, n_candidates = 5, stability = "none")

  qc <- compute_design_free_qc(
    x_raw = X,
    x_clean = fit$x_clean,
    wNN = fit$wNN,
    component_table = fit$component_table
  )

  expect_true(is.list(qc))
  expect_true(all(c("tsnr_nt_before", "tsnr_nt_after", "nn_variance_reduction_frac") %in% names(qc)))

  out_rds <- tempfile(fileext = ".rds")
  write_qc_artifact(qc, out_rds)
  expect_true(file.exists(out_rds))

  out_csv <- tempfile(fileext = ".csv")
  write_qc_artifact(qc, out_csv)
  expect_true(file.exists(out_csv))
})

test_that("compcor_denoise runs in aCompCor mode with auto wNN", {
  set.seed(6)
  X <- matrix(rnorm(220 * 120), nrow = 220, ncol = 120)
  fit <- compcor_denoise(
    x = X,
    tr = 0.8,
    mode = "acompcorr",
    n_comp = 5L,
    wnn_thresh = 0.5
  )

  expect_equal(dim(fit$x_clean), dim(X))
  expect_true(is.matrix(fit$regressors))
  expect_true(ncol(fit$regressors) <= 5L)
  expect_true(fit$diagnostics$n_selected_voxels > 0)
})

test_that("compcor_denoise runs in tCompCor mode with top variance selection", {
  set.seed(7)
  X <- matrix(rnorm(300 * 140), nrow = 300, ncol = 140)
  fit <- compcor_denoise(
    x = X,
    tr = 1.0,
    mode = "tcompcorr",
    n_comp = 6L,
    top_var_frac = 0.05
  )

  expect_equal(dim(fit$x_clean), dim(X))
  expect_true(is.matrix(fit$regressors))
  expect_true(ncol(fit$regressors) <= 6L)
  expect_true(fit$diagnostics$selector_source == "top_variance")
})

test_that("compcor_denoise supports DCT pre-highpass", {
  set.seed(8)
  X <- matrix(rnorm(260 * 160), nrow = 260, ncol = 160)
  fit <- compcor_denoise(
    x = X,
    tr = 2.0,
    mode = "tcompcorr",
    n_comp = 5L,
    pre_highpass = "dct",
    highpass_hz = 1 / 128
  )

  expect_equal(dim(fit$x_clean), dim(X))
  expect_true(fit$diagnostics$params$pre_highpass == "dct")
  expect_true(fit$diagnostics$params$highpass_n_regressors >= 1L)
})
