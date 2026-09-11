auc_lower_is_positive <- function(score, positive) {
  # AUC for "lower score indicates the positive class".
  r <- rank(-score)
  n1 <- sum(positive)
  n2 <- sum(!positive)
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

test_that("simulate_phy_data returns the documented structure", {
  sim <- simulate_phy_data(n_vox = 200, n_time = 60, tr = 1.5, nn_frac = 0.25, n_neural = 2, seed = 1)
  expect_s3_class(sim, "fmriphysio_sim")
  expect_equal(dim(sim$x), c(200L, 60L))
  expect_type(sim$nn, "logical")
  expect_length(sim$nn, 200L)
  expect_equal(sum(sim$nn), 50L)
  expect_equal(dim(sim$physio), c(60L, 2L))
  expect_equal(colnames(sim$physio), c("cardiac", "respiratory"))
  expect_equal(dim(sim$neural), c(60L, 2L))
  expect_equal(dim(sim$physio_loadings), c(200L, 2L))
  expect_equal(dim(sim$neural_loadings), c(200L, 2L))
  expect_length(sim$noise_sd, 200L)
  expect_length(sim$voxel_mean, 200L)
  expect_identical(sim$tr, 1.5)
  expect_output(print(sim), "fmriphysio_sim")
})

test_that("simulated data follow the documented generative model", {
  sim <- simulate_phy_data(n_vox = 300, n_time = 200, seed = 2)
  resid <- sim$x - sim$voxel_mean -
    sim$physio_loadings %*% t(sim$physio) -
    sim$neural_loadings %*% t(sim$neural)
  # What remains is white noise with the stated per-voxel sd.
  ratio <- apply(resid, 1, stats::sd) / sim$noise_sd
  expect_gt(stats::median(ratio), 0.93)
  expect_lt(stats::median(ratio), 1.07)
  expect_lt(max(abs(rowMeans(resid))), 5 / sqrt(200) * max(sim$noise_sd))
  # Baseline sets the voxel means; physiological loadings are larger in NN voxels.
  expect_equal(mean(sim$voxel_mean), 100, tolerance = 0.02)
  expect_gt(
    stats::median(abs(sim$physio_loadings[sim$nn, ])),
    3 * stats::median(abs(sim$physio_loadings[!sim$nn, ]))
  )
  expect_gt(
    stats::median(abs(sim$neural_loadings[!sim$nn, ])),
    3 * stats::median(abs(sim$neural_loadings[sim$nn, ]))
  )
})

test_that("seeded simulation is reproducible and leaves the caller's RNG untouched", {
  set.seed(42)
  before <- .Random.seed
  a <- simulate_phy_data(100, 40, seed = 7)
  expect_identical(.Random.seed, before)
  b <- simulate_phy_data(100, 40, seed = 7)
  expect_identical(a, b)
  c <- simulate_phy_data(100, 40, seed = 8)
  expect_false(isTRUE(all.equal(a$x, c$x)))
  # Without a seed the global stream is used (and advanced).
  simulate_phy_data(100, 40)
  expect_false(identical(.Random.seed, before))
})

test_that("wNN separates simulated non-neuronal voxels at short and long TR", {
  for (tr in c(0.8, 2)) {
    for (s in 1:3) {
      sim <- simulate_phy_data(n_vox = 500, n_time = 150, tr = tr, seed = s)
      w <- compute_wNN_diff_energy(sim$x)$wNN
      expect_gte(auc_lower_is_positive(w, sim$nn), 0.85)
    }
  }
})

test_that("physiological sources are sampled at tr and alias accordingly", {
  alias <- function(f, fs) abs(f - round(f / fs) * fs)
  for (tr in c(0.8, 2)) {
    sim <- simulate_phy_data(n_vox = 50, n_time = 400, tr = tr, seed = 3)
    resp <- sim$physio[, "respiratory"]
    spec <- Mod(stats::fft(resp - mean(resp)))^2
    freqs <- (seq_along(resp) - 1) / (length(resp) * tr)
    half <- seq_len(floor(length(resp) / 2))
    peak <- freqs[half][which.max(spec[half])]
    expected <- alias(sim$physio_hz[["respiratory"]], 1 / tr)
    expect_lt(abs(peak - expected), 0.035 + 1 / (length(resp) * tr))
  }
})

test_that("neural sources are low-pass", {
  sim <- simulate_phy_data(n_vox = 50, n_time = 300, tr = 1, seed = 4)
  freqs <- (seq_len(300) - 1) / (300 * 1)
  freqs <- pmin(freqs, 1 - freqs)
  for (j in seq_len(ncol(sim$neural))) {
    spec <- Mod(stats::fft(sim$neural[, j]))^2
    expect_gt(sum(spec[freqs <= 0.12]) / sum(spec), 0.95)
  }
})

test_that("simulate_phy_data validates its arguments", {
  expect_error(simulate_phy_data(10, 50), "n_vox")
  expect_error(simulate_phy_data(50, 10), "n_time")
  expect_error(simulate_phy_data(50, 50, tr = -1), "tr")
  expect_error(simulate_phy_data(50, 50, nn_frac = 0), "nn_frac")
  expect_error(simulate_phy_data(50, 50, nn_frac = 1), "nn_frac")
  expect_error(simulate_phy_data(50, 50, noise_sd = 0), "noise_sd")
})
