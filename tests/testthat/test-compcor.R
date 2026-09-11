# Principal angles (degrees) between the column spans of A and B.
principal_angles <- function(A, B) {
  qa <- qr.Q(qr(scale(A, scale = FALSE)))
  qb <- qr.Q(qr(scale(B, scale = FALSE)))
  acos(pmin(1, svd(crossprod(qa, qb))$d)) * 180 / pi
}

# A nuisance region driven by three known sources; the remaining voxels carry
# the same sources weakly. Optional high-variance outlier voxels in the region
# and voxel-specific linear drift.
make_roi_data <- function(seed, n_time = 200, n_roi = 100, n_rest = 200,
                          outlier_sd = 0, drift = 0) {
  set.seed(seed)
  S <- matrix(rnorm(n_time * 3), n_time, 3)
  roi <- matrix(rnorm(n_roi * 3), n_roi, 3) %*% t(S) +
    matrix(rnorm(n_roi * n_time, sd = 0.5), n_roi)
  if (outlier_sd > 0) {
    roi[1:6, ] <- roi[1:6, ] + matrix(rnorm(6 * n_time, sd = outlier_sd), 6)
  }
  rest <- matrix(rnorm(n_rest * 3, sd = 0.3), n_rest, 3) %*% t(S) +
    matrix(rnorm(n_rest * n_time), n_rest)
  X <- rbind(roi, rest) + 50
  if (drift > 0) {
    X <- X + outer(runif(n_roi + n_rest, -1, 1) * drift, seq(-1, 1, length.out = n_time))
  }
  list(x = X, roi = c(rep(TRUE, n_roi), rep(FALSE, n_rest)), sources = S)
}

test_that("aCompCor recovers the nuisance sources despite uneven voxel variance", {
  # Six ROI voxels carry large independent noise. Without per-voxel variance
  # normalization they dominate the PCA (angles of ~80 degrees).
  d <- make_roi_data(10, outlier_sd = 25)
  fit <- compcor_denoise(d$x, tr = 2, nuisance_mask = d$roi, n_comp = 3)
  expect_lt(max(principal_angles(fit$regressors, d$sources)), 15)
})

test_that("aCompCor recovers all sources in the presence of linear drift", {
  # Without detrending, the drift takes one component and a source is lost.
  d <- make_roi_data(11, drift = 8)
  fit <- compcor_denoise(d$x, tr = 2, nuisance_mask = d$roi, n_comp = 3)
  expect_lt(max(principal_angles(fit$regressors, d$sources)), 15)
  # The regressors carry no linear trend.
  tt <- seq_len(ncol(d$x))
  expect_lt(max(abs(stats::cor(fit$regressors, tt))), 1e-8)
})

test_that("tCompCor recovers sources under drift from the highest-variance voxels", {
  d <- make_roi_data(12, drift = 8)
  fit <- compcor_denoise(d$x, tr = 2, mode = "tcompcor", n_comp = 3, top_var_n = 100)
  expect_lt(max(principal_angles(fit$regressors, d$sources)), 15)
  expect_equal(sum(fit$nuisance_voxels), 100L)
  expect_identical(fit$params$top_var_n, 100L)
})

test_that("aCompCor matches an independent textbook implementation", {
  d <- make_roi_data(13)
  k <- 3
  X <- d$x
  tt <- seq_len(ncol(X))
  Y <- t(stats::residuals(stats::lm(t(X[d$roi, ]) ~ tt)))
  Y <- Y / apply(Y, 1, stats::sd)
  V <- svd(Y, nu = 0, nv = k)$v
  Xc <- X - rowMeans(X)
  expected <- Xc - (Xc %*% V) %*% t(V) + rowMeans(X)

  for (engine in c("auto", "svd")) {
    fit <- compcor_denoise(X, tr = 2, nuisance_mask = d$roi, n_comp = k, svd_engine = engine)
    expect_lt(max(principal_angles(fit$regressors, V)), 1e-4)
    expect_equal(fit$x_clean, expected, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("voxel means are retained in both modes", {
  d <- make_roi_data(14)
  for (mode in c("acompcor", "tcompcor")) {
    fit <- if (mode == "acompcor") {
      compcor_denoise(d$x, tr = 2, nuisance_mask = d$roi, n_comp = 2)
    } else {
      compcor_denoise(d$x, tr = 2, mode = mode, n_comp = 2)
    }
    expect_equal(rowMeans(fit$x_clean), rowMeans(d$x), tolerance = 1e-10)
  }
})

test_that("a rank-deficient nuisance region removes only its true rank, whatever the seed", {
  set.seed(15)
  n_time <- 200
  B <- matrix(rnorm(2 * n_time), n_time, 2)
  roi <- matrix(rnorm(10 * 2), 10, 2) %*% t(B)
  X <- rbind(roi, matrix(rnorm(100 * n_time), 100)) + 10
  nm <- c(rep(TRUE, 10), rep(FALSE, 100))

  expect_warning(
    fit1 <- compcor_denoise(X, tr = 2, nuisance_mask = nm, n_comp = 5,
                            orientation = "voxels_by_time", svd_engine = "rsvd", seed = 1),
    "Only 2 of the requested 5"
  )
  fit2 <- suppressWarnings(compcor_denoise(X, tr = 2, nuisance_mask = nm, n_comp = 5,
                                           orientation = "voxels_by_time", svd_engine = "rsvd", seed = 2))
  expect_equal(ncol(fit1$regressors), 2L)
  expect_identical(fit1$params$n_comp_used, 2L)
  expect_identical(fit1$params$n_comp, 5L)
  expect_equal(fit1$x_clean, fit2$x_clean, tolerance = 1e-8)
  # Exactly a rank-2 subspace was removed from the other voxels.
  removed <- (X - fit1$x_clean)[!nm, ]
  sv <- svd(removed)$d
  expect_lt(sv[3] / sv[1], 1e-8)
})

test_that("masked-out, constant, and non-finite voxels are returned unchanged", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 16)
  X <- sim$x
  X[5, 10] <- NA
  X[6, ] <- 7
  mask <- rep(TRUE, 300)
  mask[1:3] <- FALSE
  expect_warning(
    fit <- compcor_denoise(X, tr = 2, mask = mask, nuisance_mask = sim$nn, n_comp = 2),
    "non-finite"
  )
  expect_equal(dim(fit$x_clean), dim(X))
  for (i in c(1:3, 5, 6)) {
    expect_identical(fit$x_clean[i, ], X[i, ])
  }
  expect_false(any(fit$active[c(1:3, 5, 6)]))
  expect_true(all(fit$active[-c(1:3, 5, 6)]))
  expect_false(any(fit$nuisance_voxels[c(5, 6)]))
  expect_identical(fit$diagnostics$n_excluded_nonfinite, 1L)
  expect_identical(fit$diagnostics$n_excluded_constant, 1L)
})

test_that("the nuisance region may lie outside the analysis mask", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 17)
  full <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
  gm_only <- compcor_denoise(sim$x, tr = 2, mask = !sim$nn, nuisance_mask = sim$nn, n_comp = 2)
  expect_identical(gm_only$nuisance_voxels, sim$nn)
  expect_identical(gm_only$active, !sim$nn)
  expect_lt(max(principal_angles(gm_only$regressors, full$regressors)), 1e-6)
  expect_identical(gm_only$x_clean[sim$nn, ], sim$x[sim$nn, ])
  expect_equal(gm_only$x_clean[!sim$nn, ], full$x_clean[!sim$nn, ], tolerance = 1e-8)
})

test_that("output keeps the input's shape and orientation", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 18)
  X <- sim$x
  ref <- compcor_denoise(X, tr = 2, nuisance_mask = sim$nn, n_comp = 2, orientation = "voxels_by_time")

  tv <- compcor_denoise(t(X), tr = 2, nuisance_mask = sim$nn, n_comp = 2, orientation = "time_by_voxels")
  expect_equal(tv$x_clean, t(ref$x_clean), tolerance = 1e-10)

  expect_message(
    auto <- compcor_denoise(t(X), tr = 2, nuisance_mask = sim$nn, n_comp = 2),
    "time x voxels"
  )
  expect_equal(auto$x_clean, t(ref$x_clean), tolerance = 1e-10)

  arr <- array(X, dim = c(10, 10, 3, 120))
  fa <- compcor_denoise(arr, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
  expect_equal(dim(fa$x_clean), c(10L, 10L, 3L, 120L))
  expect_equal(fa$x_clean, array(ref$x_clean, dim = dim(arr)), tolerance = 1e-10)
})

test_that("the design guardrail leaves task estimates unchanged", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 19)
  D <- rep(rep(c(0, 1), each = 10), length.out = 120)
  fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 3, design = D)
  G <- cbind(1, D)
  expect_lt(max(abs(crossprod(G, fit$regressors))), 1e-8)
  beta_raw <- qr.coef(qr(G), t(sim$x))
  beta_clean <- qr.coef(qr(G), t(fit$x_clean))
  expect_equal(beta_clean, beta_raw, tolerance = 1e-8)
  expect_true(fit$params$design)

  # A nuisance region whose only signal lies in the design's span yields no
  # regressors. (The design includes the linear trend because the nuisance
  # data are detrended before the PCA.)
  Dc <- D - mean(D)
  X <- rbind(outer(rnorm(20, sd = 2), Dc) + 30, matrix(rnorm(100 * 120), 100))
  fit2 <- compcor_denoise(X, tr = 2, nuisance_mask = c(rep(TRUE, 20), rep(FALSE, 100)),
                          n_comp = 1, design = cbind(D, seq_len(120)),
                          orientation = "voxels_by_time")
  expect_equal(ncol(fit2$regressors), 0L)
  expect_identical(fit2$diagnostics$design_dropped, 1L)
  expect_equal(fit2$x_clean, X, tolerance = 1e-10)
})

test_that("the DCT high-pass matches SPM's basis size, keeps the constant, and respects Nyquist", {
  B <- dct_basis(120, tr = 2, highpass_hz = 1 / 128)
  expect_equal(ncol(B), floor(2 * 120 * 2 / 128) + 1)
  expect_equal(crossprod(B), diag(ncol(B)), tolerance = 1e-12)

  # Short runs still remove the constant (kmax = 0).
  expect_equal(ncol(dct_basis(20, tr = 1, highpass_hz = 1 / 128)), 1L)

  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 20)
  fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2, pre_highpass = "dct")
  basis <- fit$diagnostics$highpass$basis
  expect_equal(ncol(basis), 4L)
  expect_identical(fit$diagnostics$highpass$n_regressors, 4L)
  Xc <- fit$x_clean - rowMeans(fit$x_clean)
  expect_lt(max(abs(Xc %*% basis[, -1])), 1e-8)
  expect_equal(rowMeans(fit$x_clean), rowMeans(sim$x), tolerance = 1e-10)

  expect_error(
    compcor_denoise(sim$x, tr = 2, pre_highpass = "dct", highpass_hz = 0.25),
    "Nyquist"
  )
  expect_error(dct_basis(20, tr = 1, highpass_hz = 0.49), "degrees of freedom")
})

test_that("automatic wNN selects non-neuronal voxels on simulated data", {
  sim <- simulate_phy_data(n_vox = 600, n_time = 150, seed = 21)
  fit <- compcor_denoise(sim$x, tr = 2, n_comp = 2)
  expect_identical(fit$diagnostics$selector_source, "wnn_auto")
  expect_length(fit$wNN, 600L)
  expect_gt(mean(sim$nn[fit$nuisance_voxels]), 0.9)
  expect_lt(min(principal_angles(fit$regressors, sim$physio)), 10)
})

test_that("a supplied wnn may contain NA and cover all or only masked voxels", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 22)
  w <- ifelse(sim$nn, 0, 1)
  w[1:10] <- NA
  fit <- compcor_denoise(sim$x, tr = 2, wnn = w, n_comp = 2)
  expect_identical(fit$diagnostics$selector_source, "wnn_given")
  expect_false(any(fit$nuisance_voxels[1:10]))
  expect_identical(fit$nuisance_voxels[-(1:10)], sim$nn[-(1:10)])

  mask <- rep(c(TRUE, FALSE), 150)
  fit_m <- compcor_denoise(sim$x, tr = 2, mask = mask, wnn = ifelse(sim$nn, 0, 1)[mask], n_comp = 2)
  expect_identical(fit_m$nuisance_voxels, mask & sim$nn)
})

test_that("compcor_denoise validates inputs with clear messages", {
  sim <- simulate_phy_data(n_vox = 100, n_time = 60, seed = 23)
  x <- sim$x
  nm <- sim$nn
  expect_error(compcor_denoise(x), "`tr`")
  expect_error(compcor_denoise(x, tr = 0), "`tr`")
  expect_error(compcor_denoise(x, tr = 2, n_comp = 0), "n_comp")
  expect_error(compcor_denoise(x, tr = 2, n_comp = NA), "n_comp")
  expect_error(compcor_denoise(x, tr = 2, mode = "acompcorr"), "should be one of")
  nm_na <- nm
  nm_na[1] <- NA
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = nm_na), "Invalid `nuisance_mask`.*NA")
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = c(2.7, 5)), "Invalid `nuisance_mask`.*whole-number")
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = c(0, 5)), "Invalid `nuisance_mask`")
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = nm[-1]), "Invalid `nuisance_mask`.*length")
  expect_error(compcor_denoise(x, tr = 2, wnn = rep(NA_real_, 100)), "no non-missing")
  expect_error(compcor_denoise(x, tr = 2, wnn = rep(2, 100)), "\\[0, 1\\]")
  expect_error(compcor_denoise(x, tr = 2, wnn = rep(0.5, 7)), "length 7")
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = nm, wnn = rep(1, 100)), "not both")
  expect_error(compcor_denoise(x, tr = 2, mode = "tcompcor", nuisance_mask = nm), "acompcor")
  expect_error(compcor_denoise(x, tr = 2, wnn_thresh = 1.5), "wnn_thresh")
  expect_error(compcor_denoise(x, tr = 2, top_var_frac = 0), "top_var_frac")
  expect_error(compcor_denoise(x, tr = 2, nuisance_mask = !logical(100) & FALSE), "Invalid `nuisance_mask`")
})

test_that("a seeded randomized SVD leaves the caller's RNG untouched", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 120, seed = 24)
  set.seed(99)
  before <- .Random.seed
  a <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2, svd_engine = "rsvd", seed = 3)
  expect_identical(.Random.seed, before)
  b <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2, svd_engine = "rsvd", seed = 3)
  expect_identical(a$x_clean, b$x_clean)
})

test_that("results have a consistent structure and print compactly", {
  sim <- simulate_phy_data(n_vox = 200, n_time = 80, seed = 25)
  fit <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2, return_diagnostics = FALSE)
  expect_s3_class(fit, "compcor_denoise_result")
  expect_named(fit, c("x_clean", "regressors", "wNN", "nuisance_voxels", "active", "mode", "tr", "params"),
               ignore.order = TRUE)
  expect_null(fit$diagnostics)
  expect_equal(colnames(fit$regressors), c("compcor_1", "compcor_2"))
  expect_equal(apply(fit$regressors, 2, stats::sd), c(compcor_1 = 1, compcor_2 = 1), tolerance = 1e-10)
  out <- capture.output(print(fit))
  expect_match(out[1], "compcor_denoise_result")
  expect_lt(length(out), 10L)
})

test_that("the design guardrail and the DCT high-pass hold together", {
  sim <- simulate_phy_data(n_vox = 400, n_time = 200, tr = 2, seed = 30)
  task <- sim$neural[, 1]
  fit <- compcor_denoise(
    sim$x, tr = 2, mode = "tcompcor", n_comp = 5,
    design = cbind(task), pre_highpass = "dct"
  )
  B <- fit$diagnostics$highpass$basis
  Xc <- fit$x_clean - rowMeans(fit$x_clean)
  # Guarding against the design alone put low frequencies back (max 0.59).
  expect_lt(max(abs(Xc %*% B)), 1e-8)

  task_hp <- task - B %*% crossprod(B, task)
  D <- cbind(1, task_hp)
  b_raw <- qr.coef(qr(D), t(sim$x))
  b_cln <- qr.coef(qr(D), t(fit$x_clean))
  expect_equal(b_cln[2, ], b_raw[2, ], tolerance = 1e-8)
})

test_that("compcor_denoise returns neuroim2 images in the same space", {
  sim <- simulate_phy_data(n_vox = 384, n_time = 80, seed = 31)
  sp <- neuroim2::NeuroSpace(c(8, 8, 6, 80), spacing = c(2, 2, 3))
  nv <- neuroim2::DenseNeuroVec(array(sim$x, c(8, 8, 6, 80)), sp)
  fit <- compcor_denoise(nv, tr = 2, nuisance_mask = sim$nn, n_comp = 2)
  ref <- compcor_denoise(sim$x, tr = 2, nuisance_mask = sim$nn, n_comp = 2)

  expect_s4_class(fit$x_clean, "DenseNeuroVec")
  expect_identical(neuroim2::space(fit$x_clean), sp)
  expect_equal(unname(methods::as(fit$x_clean, "matrix")), unname(ref$x_clean), tolerance = 1e-8)
})

test_that("nuisance voxels that are pure trend yield no regressors, not spurious ones", {
  set.seed(1)
  n_t <- 120
  X <- matrix(rnorm(200 * n_t), 200, n_t) + 50
  X[1:20, ] <- outer(rnorm(20, 5), rep(1, n_t)) + outer(rnorm(20), seq_len(n_t))
  nm <- rep(c(TRUE, FALSE), c(20, 180))
  # The roundoff left after detrending was scaled to unit variance and removed
  # as 3 spurious components, changing other voxels by up to 1.7.
  expect_warning(
    fit <- compcor_denoise(X, tr = 2, nuisance_mask = nm, n_comp = 3),
    "Only 0 of the requested 3"
  )
  expect_equal(ncol(fit$regressors), 0L)
  expect_equal(fit$x_clean, X, tolerance = 1e-10)
  expect_length(fit$diagnostics$singular_values, 0L)

  # Likewise for nuisance series lying entirely in the DCT high-pass basis.
  B <- dct_basis(n_t, 2, 1 / 128)
  X2 <- X
  X2[1:20, ] <- 50 + matrix(rnorm(20 * (ncol(B) - 1)), 20) %*% t(B[, -1, drop = FALSE]) * 3
  expect_warning(
    fit2 <- compcor_denoise(X2, tr = 2, nuisance_mask = nm, n_comp = 3, pre_highpass = "dct"),
    "Only 0 of the requested 3"
  )
  expect_equal(ncol(fit2$regressors), 0L)
})
