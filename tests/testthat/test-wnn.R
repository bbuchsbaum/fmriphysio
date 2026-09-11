test_that("wNN is bounded, one per voxel, and 1 means neuronal", {
  set.seed(1)
  X <- matrix(rnorm(200 * 120), 200, 120)
  # Voxels 1:30 get extra high-frequency variance (non-neuronal-like).
  X[1:30, ] <- X[1:30, ] + 3 * matrix(rnorm(30 * 120), 30, 120)
  out <- compute_wNN_diff_energy(X)

  expect_length(out$wNN, nrow(X))
  expect_true(all(out$wNN >= 0 & out$wNN <= 1))
  expect_equal(out$posterior_nn, 1 - out$wNN)
  expect_lt(mean(out$wNN[1:30]), 0.2)
  expect_gt(mean(out$wNN[31:200]), 0.8)
})

test_that("difference energy and excess energy are computed exactly", {
  set.seed(2)
  X <- matrix(rnorm(50 * 40), 50, 40)
  out <- compute_wNN_diff_energy(X)
  expect_equal(out$energy, apply(X, 1, function(v) mean(diff(v)^2)))
  expect_true(all(out$delta >= 0))
})

test_that("percentile ramp maps the quantile thresholds to 1 and 0", {
  set.seed(3)
  X <- matrix(rnorm(300 * 60), 300, 60) * rep(stats::runif(300, 0.5, 3), 60)
  out <- compute_wNN_diff_energy(X, nt_q = 0.5, nn_q = 0.9)
  expect_true(all(out$wNN[out$delta <= out$delta_nt] == 1))
  expect_true(all(out$wNN[out$delta >= out$delta_nn] == 0))
  mid <- out$delta > out$delta_nt & out$delta < out$delta_nn
  expect_equal(out$wNN[mid], 1 - (out$delta[mid] - out$delta_nt) / (out$delta_nn - out$delta_nt))
})

test_that("wNN does not depend on the units of the data", {
  set.seed(4)
  X <- matrix(rnorm(400 * 80), 400, 80)
  X[1:60, ] <- X[1:60, ] * 4
  for (thr in c("percentile", "mixture")) {
    ref <- compute_wNN_diff_energy(X, threshold = thr)$wNN
    for (s in c(1e-3, 1e3)) {
      # The mixture fit used log1p(delta): 0% vs 59% flagged across scales.
      expect_equal(compute_wNN_diff_energy(X * s, threshold = thr)$wNN, ref, tolerance = 1e-6)
    }
  }
})

test_that("mixture thresholding separates a high-energy subpopulation", {
  set.seed(5)
  X <- matrix(rnorm(400 * 100), 400, 100)
  X[1:80, ] <- X[1:80, ] * 3
  out <- compute_wNN_diff_energy(X, threshold = "mixture")
  expect_lt(mean(out$wNN[1:80]), 0.2)
  expect_gt(mean(out$wNN[81:400]), 0.9)
})

test_that("wNN works for small numbers of voxels", {
  set.seed(6)
  # k_fit = max(10, ...) made the line fit fail for fewer than 10 voxels.
  out <- compute_wNN_diff_energy(matrix(rnorm(5 * 30), 5, 30), nt_q = 0.2, nn_q = 0.8)
  expect_length(out$wNN, 5)
  expect_true(all(is.finite(out$wNN)))
})

test_that("compute_wNN_diff_energy validates its arguments", {
  X <- matrix(rnorm(80 * 60), 80, 60)
  expect_error(compute_wNN_diff_energy(matrix(rnorm(20 * 2), 20, 2)), "at least 3 time points")
  expect_error(compute_wNN_diff_energy(matrix(rnorm(2 * 20), 2, 20)), "at least 3 voxels")
  expect_error(compute_wNN_diff_energy(X, nt_q = 0.95, nn_q = 0.74), "smaller than")
  expect_error(compute_wNN_diff_energy(X, nt_q = 74, nn_q = 95), "nt_q")
  X[3, 4] <- NA
  expect_error(compute_wNN_diff_energy(X), "finite numeric")
})

test_that("coincident percentile thresholds give an actionable error", {
  X <- matrix(rep(sin(1:60), each = 50), 50, 60)
  expect_error(compute_wNN_diff_energy(X), "thresholds coincide")
})

test_that("build_wnn_ramp is piecewise linear", {
  r <- build_wnn_ramp(c(0, 0.25, 0.5, 0.75, 1), delta_nt = 0.25, delta_nn = 0.75)
  expect_equal(r, c(1, 1, 0.5, 0, 0))
})

test_that("fit_gaussian_mixture_1d finds two separated clusters", {
  set.seed(7)
  x <- c(rnorm(200, 0, 1), rnorm(100, 6, 1))
  fit <- fit_gaussian_mixture_1d(x)
  expect_equal(sort(fit$mu), c(0, 6), tolerance = 0.3)
  expect_gt(mean(fit$post_high[201:300]), 0.95)
  expect_lt(mean(fit$post_high[1:200]), 0.05)
  expect_error(fit_gaussian_mixture_1d(rnorm(10)), "Too few finite values")
})
