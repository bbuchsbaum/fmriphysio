test_that("CAA recovers a positively autocorrelated source and reports its lag-1 correlation", {
  set.seed(4)
  n <- 300
  a <- as.numeric(stats::arima.sim(list(ar = 0.88), n))
  Z <- rbind(a, matrix(rnorm(5 * n), 5, n))
  caa <- extract_caa_oneshot(Z, n_candidates = 3)

  expect_gt(abs(cor(caa$scores[, 1], a)), 0.95)
  # Predictability is a canonical correlation, not rho / (T - 1) (old bug gave ~0.002).
  expect_equal(caa$predictability[1], abs(stats::acf(a, plot = FALSE)$acf[2]), tolerance = 0.1)
  expect_gt(caa$autocorrelation[1], 0.7)
})

test_that("CAA recovers negatively autocorrelated (near-Nyquist) sources", {
  set.seed(3)
  n <- 300
  a <- as.numeric(stats::arima.sim(list(ar = -0.8), n))
  Z <- rbind(a, matrix(rnorm(5 * n), 5, n))
  caa <- extract_caa_oneshot(Z, n_candidates = 3)

  # Averaging the two lagged projections cancelled these components (|cor| 0.04).
  expect_gt(abs(cor(caa$scores[, 1], a)), 0.95)
  expect_lt(caa$autocorrelation[1], -0.6)
  expect_gt(caa$predictability[1], 0.6)
})

test_that("CAA recovers an aliased cardiac-like oscillation", {
  set.seed(9)
  n <- 240
  card <- sin(2 * pi * cumsum(0.42 + 0.01 * rnorm(n)))
  Z <- rbind(card + 0.3 * rnorm(n), matrix(rnorm(6 * n), 6, n))
  caa <- extract_caa_oneshot(Z, n_candidates = 4)
  best <- max(abs(cor(caa$scores, card)))
  expect_gt(best, 0.9)
})

test_that("CAA output has the documented structure", {
  set.seed(12)
  Z <- matrix(rnorm(8 * 50), 8, 50)
  out <- extract_caa_oneshot(Z, n_candidates = 6)

  expect_equal(dim(out$scores), c(50L, 6L))
  expect_equal(colnames(out$scores), paste0("caa_", 1:6))
  expect_length(out$predictability, 6)
  expect_length(out$autocorrelation, 6)
  expect_identical(out$method, "caa")
  expect_true(all(out$predictability >= 0 & out$predictability <= 1))
  expect_false(is.unsorted(rev(out$predictability)))
  expect_equal(unname(colMeans(out$scores)), rep(0, 6), tolerance = 1e-10)
  expect_equal(unname(apply(out$scores, 2, stats::sd)), rep(1, 6), tolerance = 1e-10)
})

test_that("CAA caps n_candidates at the number of input series", {
  set.seed(13)
  out <- extract_caa_oneshot(matrix(rnorm(3 * 40), 3, 40), n_candidates = 10)
  expect_equal(ncol(out$scores), 3L)
})

test_that("CAA validates its input", {
  expect_error(extract_caa_oneshot(matrix(rnorm(8 * 3), 8, 3)), "at least 4 time points")
  expect_error(extract_caa_oneshot(matrix(rnorm(20), 4, 5), n_candidates = 0), "n_candidates")
  # 10 series over 6 time points gave spurious canonical correlations of 1.
  expect_error(extract_caa_oneshot(matrix(rnorm(60), 10, 6)), "at most T - 2 series")
  Z <- matrix(rnorm(40), 4, 10)
  Z[2, 3] <- NA
  expect_error(extract_caa_oneshot(Z), "finite numeric")
})

test_that("dynamic extractors from dipca return scores and predictability", {
  skip_if_not_installed("dipca")
  set.seed(12)
  Z <- matrix(rnorm(8 * 60), 8, 60)

  d1 <- extract_dynamic_components(Z, extractor = "dicca", lag_order = 1L, n_candidates = 3L)
  expect_identical(d1$method, "dicca")
  expect_true(is.matrix(d1$scores))
  expect_equal(nrow(d1$scores), 60L)
  expect_length(d1$predictability, ncol(d1$scores))

  d2 <- extract_dynamic_components(Z, extractor = "dipca", lag_order = 1L, n_candidates = 3L)
  expect_identical(d2$method, "dipca")
  expect_true(is.matrix(d2$scores))
})

test_that("dynamic extractors cap n_candidates at the reduced rank", {
  skip_if_not_installed("dipca")
  set.seed(14)
  Z <- matrix(rnorm(4 * 60), 4, 60)
  expect_no_error(extract_dynamic_components(Z, extractor = "dicca", n_candidates = 30L))
})
