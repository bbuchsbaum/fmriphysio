fake_fit <- function(x) {
  list(
    x_clean = x,
    regressors = matrix(0, ncol(x), 0),
    diagnostics = list(timings = list(total_seconds = 0.5, wnn_seconds = 0.1, pass_table = data.frame()))
  )
}

test_that("every engine is benchmarked on the same simulated dataset", {
  seen <- list()
  local_mocked_bindings(
    fast_phy_denoise = function(x, svd_engine, seed, ...) {
      seen[[length(seen) + 1L]] <<- list(x = x, engine = svd_engine, seed = seed)
      fake_fit(x)
    }
  )
  res <- benchmark_fast_phy_denoise(
    grid_n_vox = 100, grid_n_time = 40, engines = c("svd", "rsvd", "auto"), reps = 2, seed = 5
  )
  expect_equal(nrow(res), 6L)
  expect_length(seen, 6L)
  # Within a dataset, all engines receive identical data and seed ...
  expect_identical(seen[[1]]$x, seen[[2]]$x)
  expect_identical(seen[[1]]$x, seen[[3]]$x)
  expect_identical(seen[[4]]$x, seen[[6]]$x)
  expect_identical(vapply(seen, `[[`, character(1), "engine"), rep(c("svd", "rsvd", "auto"), 2))
  # ... and repetitions use different data.
  expect_false(identical(seen[[1]]$x, seen[[4]]$x))
  expect_identical(res$dataset_seed, rep(c(6, 7), each = 3))
  # The data are exactly what simulate_phy_data produces for that seed.
  expect_identical(seen[[1]]$x, simulate_phy_data(100, 40, tr = 2, seed = 6)$x)
  expect_true(all(res$status == "ok"))
  expect_equal(res$physio_kept_nt, rep(1, 6))
})

test_that("benchmark rejects reserved and unnamed arguments and validates inputs", {
  expect_error(benchmark_fast_phy_denoise(x = 1), "`x` is set by the benchmark")
  expect_error(benchmark_fast_phy_denoise(svd_engine = "svd"), "`svd_engine` is set")
  expect_error(benchmark_fast_phy_denoise(return_diagnostics = FALSE), "`return_diagnostics` is set")
  # An unnamed value reaches `...` only after every formal is filled.
  expect_error(
    benchmark_fast_phy_denoise(100, 40, "svd", 1, 2, "caa", NULL, 1, 3),
    "must be named"
  )
  expect_error(benchmark_fast_phy_denoise(reps = 0), "reps")
  expect_error(benchmark_fast_phy_denoise(engines = "qr"), "engines")
  expect_error(benchmark_fast_phy_denoise(grid_n_vox = 5), "grid_n_vox")
})

test_that("errors inside the denoiser are reported per row", {
  local_mocked_bindings(fast_phy_denoise = function(...) stop("boom"))
  res <- benchmark_fast_phy_denoise(grid_n_vox = 60, grid_n_time = 30, engines = "svd")
  expect_identical(res$status, "error")
  expect_identical(res$error_message, "boom")
  expect_true(is.na(res$total_seconds))
})

test_that("an end-to-end benchmark runs, scores against ground truth, and leaves the RNG untouched", {
  out <- tempfile(fileext = ".csv")
  set.seed(123)
  before <- .Random.seed
  res <- benchmark_fast_phy_denoise(
    grid_n_vox = 300, grid_n_time = 100, engines = "svd", reps = 1, seed = 1, out_file = out
  )
  expect_identical(.Random.seed, before)
  expect_identical(res$status, "ok")
  expect_true(all(c("physio_kept_nt", "neural_kept_nt", "selected_components", "total_seconds") %in% names(res)))
  expect_true(res$physio_kept_nt >= 0 && res$physio_kept_nt <= 1)
  expect_true(res$neural_kept_nt >= 0 && res$neural_kept_nt <= 1)
  expect_gte(res$total_seconds, 0)
  expect_equal(nrow(utils::read.csv(out)), 1L)
})

test_that("source_energy_kept measures the surviving source energy", {
  set.seed(32)
  S <- matrix(rnorm(100 * 2), 100, 2)
  X <- matrix(rnorm(5 * 2), 5, 2) %*% t(S) + 10
  expect_equal(source_energy_kept(X, X, S), rep(1, 5))
  Q <- qr.Q(qr(scale(S, scale = FALSE)))
  Xc <- X - rowMeans(X)
  cleaned <- Xc - (Xc %*% Q) %*% t(Q)
  expect_equal(source_energy_kept(X, cleaned, S), rep(0, 5), tolerance = 1e-10)
})

test_that("the benchmark reserves orientation", {
  # Passing it through `...` used to turn every row into an error.
  expect_error(
    benchmark_fast_phy_denoise(
      grid_n_vox = 100, grid_n_time = 60, engines = "svd", orientation = "voxels_by_time"
    ),
    "`orientation` is set by the benchmark"
  )
})
