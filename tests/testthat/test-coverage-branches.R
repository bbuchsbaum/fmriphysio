test_that("truncated_svd covers engine branches and rank checks", {
  set.seed(11)
  X <- matrix(rnorm(60 * 20), 60, 20)

  out_svd <- truncated_svd(X, rank = 5, svd_engine = "svd")
  expect_true(out_svd$engine == "svd")
  expect_equal(length(out_svd$d), 5)

  out_rsvd <- truncated_svd(X, rank = 5, svd_engine = "rsvd")
  expect_true(out_rsvd$engine == "rsvd")
  expect_equal(length(out_rsvd$d), 5)

  expect_error(truncated_svd(X, rank = 0), "rank must be >= 1")
})

test_that("extractor helpers cover error and alternate paths", {
  set.seed(12)
  Z <- matrix(rnorm(8 * 50), 8, 50)

  out <- extract_caa_oneshot(Z, n_candidates = 6)
  expect_equal(ncol(out$scores), 6)
  expect_equal(length(out$predictability), 6)

  expect_error(extract_caa_oneshot(matrix(rnorm(8 * 3), 8, 3)), "Need at least 4 time points")

  if (requireNamespace("dipca", quietly = TRUE)) {
    d1 <- extract_dynamic_components(Z, extractor = "dicca", lag_order = 1L, n_candidates = 3L)
    expect_true(d1$method == "dicca")
    expect_true(is.matrix(d1$scores))

    d2 <- extract_dynamic_components(Z, extractor = "dipca", lag_order = 1L, n_candidates = 3L)
    expect_true(d2$method == "dipca")
    expect_true(is.matrix(d2$scores))
  }
})

test_that("score_components_nn_nt covers edge branches", {
  set.seed(13)
  X <- matrix(rnorm(40 * 60), 40, 60)
  T0 <- matrix(0, 60, 0)
  w <- c(rep(0.2, 10), rep(0.8, 30))

  out0 <- score_components_nn_nt(X, T0, w, use_cpp = FALSE)
  expect_length(out0$ratio, 0)

  expect_error(score_components_nn_nt(X, matrix(rnorm(60 * 2), 60, 2), rep(1, 40), use_cpp = FALSE), "Need both NN and NT")

  Tz <- matrix(0, 60, 2)
  outz <- score_components_nn_nt(X, Tz, w, use_cpp = FALSE)
  expect_true(all(outz$ratio == 0))
})

test_that("stability_filter covers R branch and empty input", {
  set.seed(14)
  T_tl <- matrix(rnorm(100 * 4), 100, 4)
  split_obj <- list(
    scores_h1 = matrix(rnorm(50 * 4), 50, 4),
    scores_h2 = matrix(rnorm(50 * 4), 50, 4),
    idx_h1 = 1:50,
    idx_h2 = 51:100
  )

  keep <- stability_filter(T_tl, split_obj, thresh = 0.1, use_cpp = FALSE)
  expect_equal(length(keep), 4)
  expect_true(is.logical(keep))

  expect_equal(length(stability_filter(matrix(0, 100, 0), split_obj, use_cpp = FALSE)), 0)
})

test_that("utils_matrix helpers cover mask and error branches", {
  set.seed(15)
  X <- matrix(rnorm(30 * 10), 30, 10)

  p1 <- as_nt_matrix(X, mask = 1:10)
  expect_equal(dim(p1$X), c(10, 10))

  p2 <- as_nt_matrix(t(X)) # triggers transpose to N x T
  expect_equal(dim(restore_orientation(p2$X, p2$info)), dim(t(X)))

  expect_error(as_nt_matrix(X, mask = rep(TRUE, 5)), "Logical mask length")
  expect_error(as_nt_matrix(matrix(letters[1:20], 4, 5)), "must be numeric")

  z <- safe_scale_vec(rep(1, 10))
  expect_true(all(z == 0))

  Tcomp <- matrix(rnorm(20 * 3), 20, 3)
  expect_error(orthogonalize_to_design(Tcomp, matrix(rnorm(10 * 2), 10, 2)), "design rows")

  Dz <- matrix(0, 20, 2) # rank zero
  Tproj <- orthogonalize_to_design(Tcomp, Dz)
  expect_equal(dim(Tproj), dim(Tcomp))

  Xproj <- project_out_components(X, matrix(0, nrow(X), 0))
  expect_equal(Xproj, X)
})

test_that("wnn helpers cover mixture and error branches", {
  set.seed(16)
  X <- matrix(rnorm(80 * 60), 80, 60)

  out_mix <- compute_wNN_diff_energy(X, threshold = "mixture")
  expect_equal(length(out_mix$wNN), 80)
  expect_true(all(is.finite(out_mix$posterior_nn)))

  expect_error(compute_wNN_diff_energy(matrix(rnorm(20 * 2), 20, 2)), "at least 3 time points")
  expect_error(compute_wNN_diff_energy(X, nt_q = 0.95, nn_q = 0.74), "Invalid percentile")

  r <- build_wnn_ramp(c(0, 0.5, 1), delta_nt = 0.25, delta_nn = 0.75)
  expect_true(r[1] == 1 && r[3] == 0 && r[2] < 1 && r[2] > 0)

  expect_error(fit_gaussian_mixture_1d(rnorm(10)), "Too few finite values")
})

test_that("compcor_denoise covers nuisance mask and wnn branches", {
  set.seed(17)
  X <- matrix(rnorm(120 * 90), 120, 90)

  fit_mask <- compcor_denoise(
    x = X, tr = 1.0, mode = "acompcorr", n_comp = 4,
    nuisance_mask = c(1:15, 30:40)
  )
  expect_true(fit_mask$diagnostics$selector_source == "nuisance_mask")

  wnn <- c(rep(0.1, 20), rep(0.9, 100))
  fit_wnn <- compcor_denoise(
    x = X, tr = 1.0, mode = "acompcorr", n_comp = 4,
    wnn = wnn, wnn_thresh = 0.5
  )
  expect_true(fit_wnn$diagnostics$selector_source == "wnn_given")

  expect_error(
    compcor_denoise(x = X, tr = 1.0, mode = "acompcorr", wnn = rep(1, 120), wnn_thresh = 0.0),
    "No voxels selected"
  )
})

test_that("fast_phy_denoise covers wnn-provided, pred_thresh and half-stability paths", {
  set.seed(18)
  X <- simulate_phy_data(400, 180, nn_frac = 0.2)
  wnn <- compute_wNN_diff_energy(X)$wNN

  fit_half <- fast_phy_denoise(
    x = X, tr = 1.5, wnn = wnn, extractor = "caa",
    stability = "half", ratio_thresh = 0.0, max_passes = 1,
    use_cpp = FALSE
  )
  expect_true(is.list(fit_half$diagnostics$pass_log))
  expect_equal(dim(fit_half$x_clean), dim(X))

  fit_pred <- fast_phy_denoise(
    x = X, tr = 1.5, wnn = wnn, extractor = "caa",
    pred_thresh = 1.0, ratio_thresh = 0.0, stability = "none",
    max_passes = 1, use_cpp = FALSE
  )
  expect_equal(dim(fit_pred$x_clean), dim(X))
})

test_that("QC and artifact utilities cover additional branches", {
  set.seed(19)
  X <- matrix(rnorm(100 * 80), 100, 80)
  fit <- fast_phy_denoise(X, tr = 1.0, extractor = "caa", pca_rank = 12, n_candidates = 4, stability = "none")
  qc <- compute_design_free_qc(X, fit$x_clean, fit$wNN, fit$component_table)

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    fjson <- tempfile(fileext = ".json")
    write_qc_artifact(qc, fjson)
    expect_true(file.exists(fjson))
  }

  expect_error(write_qc_artifact(qc, tempfile(fileext = ".txt")), "Unsupported artifact")
  expect_error(compute_design_free_qc(X, X[, -1], fit$wNN), "matching dimensions")
  expect_error(compute_design_free_qc(X, X, fit$wNN[-1]), "wNN length")
})

test_that("Rcpp kernels can be called directly when symbols are loaded", {
  skip_if_not(is.loaded("_phynd_score_components_nn_nt_cpp"))
  skip_if_not(exists("score_components_nn_nt_cpp", mode = "function"))
  skip_if_not(is.loaded("_phynd_stability_filter_cpp"))
  skip_if_not(exists("stability_filter_cpp", mode = "function"))

  set.seed(20)
  X <- matrix(rnorm(30 * 40), 30, 40)
  T <- matrix(rnorm(40 * 3), 40, 3)
  w <- c(rep(0.2, 10), rep(0.8, 20))

  out <- score_components_nn_nt_cpp(X, T, w, 1e-12)
  expect_equal(length(out$ratio), 3)

  out0 <- score_components_nn_nt_cpp(X, matrix(0, 40, 2), w, 1e-12)
  expect_true(all(out0$ratio == 0))

  keep0 <- stability_filter_cpp(
    matrix(0, 40, 0),
    matrix(0, 20, 2),
    matrix(0, 20, 2),
    1:20, 21:40, 0.3
  )
  expect_equal(length(keep0), 0)

  keep <- stability_filter_cpp(
    matrix(rnorm(40 * 2), 40, 2),
    matrix(rnorm(20 * 3), 20, 3),
    matrix(rnorm(20 * 3), 20, 3),
    1:20, 21:40, 0.0
  )
  expect_equal(length(keep), 2)
})
