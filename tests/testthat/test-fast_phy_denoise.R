test_that("defaults remove physiology, keep neural signal, and beat PHYCAA+ settings", {
  kept <- function(fit, sim, src) {
    subspace_energy(fit$x_clean, sim[[src]], !sim$nn) / subspace_energy(sim$x, sim[[src]], !sim$nn)
  }
  runs <- expand.grid(tr = c(0.8, 2), seed = 1:12)
  res <- t(mapply(function(tr, seed) {
    sim <- simulate_phy_data(n_vox = 1000, n_time = 300, tr = tr, seed = seed)
    d <- fast_phy_denoise(sim$x, tr = tr, return_diagnostics = FALSE)
    p <- fast_phy_denoise(
      sim$x, tr = tr, extraction_weight = "neuronal", pca_rank = 100, ratio_thresh = 1,
      return_diagnostics = FALSE
    )
    c(
      d_phys = kept(d, sim, "physio"), d_neur = kept(d, sim, "neural"),
      p_phys = kept(p, sim, "physio"), p_neur = kept(p, sim, "neural")
    )
  }, runs$tr, runs$seed))
  med <- apply(res, 2, stats::median)

  # Energy fractions in neuronal voxels relative to the raw data, over 24
  # simulated runs. Medians, because a rhythm that aliases into the neural
  # band is occasionally missed (see ?fast_phy_denoise).
  expect_lt(med[["d_phys"]], 0.10)
  expect_gt(med[["d_neur"]], 0.92)
  expect_gt(min(res[, "d_neur"]), 0.75)
  # The PHYCAA+-weighted, rank-100, ratio > 1 settings (the previous defaults)
  # removed less physiology and more neural signal.
  expect_lt(med[["d_phys"]], med[["p_phys"]])
  expect_gt(med[["d_neur"]], med[["p_neur"]])
  expect_gte(mean(res[, "d_phys"] < res[, "p_phys"]), 0.75)
})

test_that("output keeps the input's shape, orientation, and voxel means", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 2, baseline = 1000)
  fit <- fast_phy_denoise(fx$x, tr = 2)

  expect_equal(dim(fit$x_clean), dim(fx$x))
  # Means used to be removed and never restored (row means ~1e-13).
  expect_equal(rowMeans(fit$x_clean), rowMeans(fx$x), tolerance = 1e-8)

  fit_t <- fast_phy_denoise(t(fx$x), orientation = "time_by_voxels")
  expect_equal(fit_t$x_clean, t(fit$x_clean), tolerance = 1e-8)

  arr <- array(fx$x, c(10, 6, 5, 120))
  fit_a <- fast_phy_denoise(arr)
  expect_equal(dim(fit_a$x_clean), dim(arr))
  expect_equal(matrix(fit_a$x_clean, 300, 120), fit$x_clean, tolerance = 1e-8)
})

test_that("voxels x time input with fewer voxels than time points is not transposed", {
  x <- core_fixture(n_vox = 60, n_time = 200, seed = 3)$x
  fit <- fast_phy_denoise(x, orientation = "voxels_by_time", ratio_thresh = 1)
  expect_length(fit$wNN, 60)
  expect_equal(nrow(fit$regressors), 200L)
  expect_message(fast_phy_denoise(x, ratio_thresh = 1), "treating it as time x voxels")
})

test_that("constant background voxels do not change the result for brain voxels", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 4)
  ref <- fast_phy_denoise(fx$x)
  xz <- rbind(fx$x, matrix(0, 300, 120), matrix(500, 50, 120))
  fit <- fast_phy_denoise(xz)

  # Zero rows used to make every candidate look like nuisance (30/30 selected).
  expect_equal(fit$x_clean[1:300, ], ref$x_clean, tolerance = 1e-8)
  expect_equal(fit$regressors, ref$regressors, tolerance = 1e-8)
  expect_identical(fit$x_clean[301:650, ], xz[301:650, ])
  expect_equal(sum(fit$active), 300L)
  expect_true(all(is.na(fit$wNN[301:650])))
  expect_equal(fit$diagnostics$voxels$n_constant, 350L)
})

test_that("masked-out and non-finite voxels are returned unchanged", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 5)
  x <- fx$x
  x[100, 7] <- NA
  mask <- rep(TRUE, 300)
  mask[201:220] <- FALSE

  expect_warning(fit <- fast_phy_denoise(x, mask = mask), "non-finite")
  expect_equal(dim(fit$x_clean), dim(x))
  expect_identical(fit$x_clean[201:220, ], x[201:220, ])
  expect_identical(fit$x_clean[100, ], x[100, ])
  expect_false(fit$active[100])

  fit_idx <- suppressWarnings(fast_phy_denoise(x, mask = which(mask)))
  expect_equal(fit_idx$x_clean, fit$x_clean)
})

test_that("wnn may be given per voxel or per masked voxel, and is validated", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 6)
  x <- fx$x
  x[4, ] <- 7 # constant voxel inside the mask: excluded from the analysis
  mask <- rep(c(TRUE, TRUE, FALSE), 100)
  w <- ifelse(fx$nn, 0, 1)

  f1 <- fast_phy_denoise(x, mask = mask, wnn = w)
  f2 <- fast_phy_denoise(x, mask = mask, wnn = w[mask])
  expect_false(f1$active[4])
  expect_equal(f1$x_clean, f2$x_clean)
  expect_equal(f1$wNN[f1$active], w[f1$active])
  expect_equal(f2$wNN[f2$active], w[f2$active])

  expect_error(fast_phy_denoise(fx$x, wnn = w[-1]), "one value per voxel")
  expect_error(fast_phy_denoise(fx$x, wnn = w * 2), "in \\[0, 1\\]")
  expect_error(fast_phy_denoise(fx$x, wnn = rep(1, 300)), "must classify some voxels")
})

test_that("the subspace rank is capped at floor(T / 3)", {
  fx <- core_fixture(n_vox = 300, n_time = 60, seed = 7)
  fit <- fast_phy_denoise(fx$x, pca_rank = 100)
  expect_equal(fit$params$rank_used, 20L)
  expect_equal(fast_phy_denoise(fx$x, pca_rank = 5)$params$rank_used, 5L)
})

test_that("removed regressors are orthogonal to the output over multiple passes", {
  fx <- core_fixture(n_vox = 400, n_time = 160, seed = 8)
  fit <- fast_phy_denoise(fx$x, max_passes = 2, ratio_thresh = 0.5)
  expect_true(any(grepl("^pass2_", colnames(fit$regressors))))

  Xc <- fit$x_clean - rowMeans(fit$x_clean)
  # Invariant: the output is orthogonal to every removed regressor.
  expect_lt(max(abs(Xc %*% qr.Q(qr(fit$regressors)))), 1e-8)
  expect_equal(fit$diagnostics$n_passes_executed, 2L)
})

test_that("the design guardrail leaves task effects unchanged", {
  fx <- core_fixture(n_vox = 400, n_time = 160, seed = 9)
  task <- fx$neur[, 1]
  fit <- fast_phy_denoise(fx$x, design = cbind(task), ratio_thresh = 1)
  expect_gt(ncol(fit$regressors), 0)
  expect_lt(max(abs(stats::cor(fit$regressors, task))), 1e-8)

  D <- cbind(1, task)
  b_raw <- qr.coef(qr(D), t(fx$x))
  b_cln <- qr.coef(qr(D), t(fit$x_clean))
  expect_equal(b_cln[2, ], b_raw[2, ], tolerance = 1e-8)
})

test_that("component table rows stay aligned with candidates when the guardrail drops one", {
  fx <- core_fixture(n_vox = 400, n_time = 160, seed = 10)
  X0 <- fx$x - rowMeans(fx$x)
  w <- compute_wNN_diff_energy(X0)$wNN
  sv <- truncated_svd(X0 * (1 - w), 20, return_u = FALSE)
  ext <- extract_caa_oneshot(t(sv$v), 30)

  fit <- fast_phy_denoise(fx$x, design = ext$scores[, 2])
  tb <- fit$component_table
  expect_false(2 %in% tb$candidate)
  # Predictability used to shift onto the wrong rows after a drop.
  expect_equal(tb$predictability, ext$predictability[tb$candidate], tolerance = 1e-8)
})

test_that("split-half stability removes nothing from pure noise", {
  set.seed(11)
  X <- matrix(rnorm(500 * 200), 500, 200)
  fit <- fast_phy_denoise(X, wnn = rep(c(0, 1), c(100, 400)), ratio_thresh = 0, stability = "half")
  expect_equal(ncol(fit$regressors), 0L)
  expect_equal(fit$x_clean, X, tolerance = 1e-10)
  expect_true(all(is.finite(fit$component_table$stability)))
})

test_that("pred_thresh filters on canonical lag-1 correlation", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 12)
  fit <- fast_phy_denoise(fx$x, ratio_thresh = 0, pred_thresh = 0.9)
  tb <- fit$component_table
  expect_true(all(tb$predictability >= 0 & tb$predictability <= 1))
  expect_true(all(tb$predictability[tb$selected] >= 0.9))
  expect_true(any(!tb$selected))
})

test_that("extraction_weight sets the weighting of the decomposition", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 13)
  X0 <- fx$x - rowMeans(fx$x)
  w <- compute_wNN_diff_energy(X0)$wNN
  weights <- list(non_neuronal = 1 - w, neuronal = w, none = rep(1, length(w)))
  for (ew in names(weights)) {
    fit <- fast_phy_denoise(fx$x, extraction_weight = ew, ratio_thresh = 0)
    k <- min(20, sum(weights[[ew]] > 0))
    ext <- extract_caa_oneshot(t(truncated_svd(X0 * weights[[ew]], k, return_u = FALSE)$v), 30)
    expect_identical(fit$params$extraction_weight, ew)
    expect_equal(fit$component_table$predictability, ext$predictability, tolerance = 1e-8)
  }
})

test_that("compiled and R paths give identical results", {
  skip_if_not(cpp_kernel_loaded())
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 14)
  a <- fast_phy_denoise(fx$x, ratio_thresh = 1, use_cpp = TRUE, return_diagnostics = FALSE)
  b <- fast_phy_denoise(fx$x, ratio_thresh = 1, use_cpp = FALSE, return_diagnostics = FALSE)
  a$params$use_cpp <- b$params$use_cpp <- NULL
  expect_equal(a, b, tolerance = 1e-10)
})

test_that("results are deterministic and never touch the caller's RNG", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 15)
  set.seed(99)
  before <- .Random.seed

  a <- fast_phy_denoise(fx$x, return_diagnostics = FALSE)
  b <- fast_phy_denoise(fx$x, return_diagnostics = FALSE)
  expect_identical(a$x_clean, b$x_clean)

  r1 <- fast_phy_denoise(fx$x, svd_engine = "rsvd", seed = 3, return_diagnostics = FALSE)
  r2 <- fast_phy_denoise(fx$x, svd_engine = "rsvd", seed = 3, return_diagnostics = FALSE)
  expect_identical(r1$x_clean, r2$x_clean)
  # set.seed() inside the function used to reset the caller's stream.
  expect_identical(.Random.seed, before)
})

test_that("an empty selection returns the input unchanged", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 16, baseline = 50)
  fit <- fast_phy_denoise(fx$x, ratio_thresh = 1e6)
  expect_equal(ncol(fit$regressors), 0L)
  expect_equal(fit$x_clean, fx$x, tolerance = 1e-10)
  expect_equal(fit$diagnostics$pass_log[[1]]$reason, "no_components_selected")
})

test_that("neuroim2 images come back as neuroim2 images in the same space", {
  fx <- core_fixture(n_vox = 384, n_time = 80, seed = 17)
  arr <- array(fx$x, c(8, 8, 6, 80))
  sp <- neuroim2::NeuroSpace(c(8, 8, 6, 80), spacing = c(2, 2, 3), origin = c(-8, -8, -9))
  nv <- neuroim2::DenseNeuroVec(arr, sp)

  fit_nv <- fast_phy_denoise(nv)
  fit_arr <- fast_phy_denoise(arr)
  # NeuroVec input used to come back as a plain array without its NeuroSpace.
  expect_s4_class(fit_nv$x_clean, "DenseNeuroVec")
  expect_identical(neuroim2::space(fit_nv$x_clean), sp)
  expect_equal(as.vector(fit_nv$x_clean@.Data), as.vector(fit_arr$x_clean), tolerance = 1e-8)
  expect_equal(
    compute_design_free_qc(nv, fit_nv$x_clean, fit_nv$wNN),
    compute_design_free_qc(arr, fit_arr$x_clean, fit_arr$wNN)
  )
  expect_error(fast_phy_denoise(nv, orientation = "time_by_voxels"), "voxels x time")
})

test_that("mask volumes and sparse images restrict the analysis and are preserved", {
  fx <- core_fixture(n_vox = 384, n_time = 80, seed = 18)
  arr <- array(fx$x, c(8, 8, 6, 80))
  sp <- neuroim2::NeuroSpace(c(8, 8, 6, 80), spacing = c(2, 2, 3))
  nv <- neuroim2::DenseNeuroVec(arr, sp)
  mvol <- neuroim2::LogicalNeuroVol(
    array(rep(c(TRUE, FALSE), c(300, 84)), c(8, 8, 6)),
    neuroim2::NeuroSpace(c(8, 8, 6), spacing = c(2, 2, 3))
  )

  fit_m <- fast_phy_denoise(nv, mask = mvol)
  expect_equal(sum(fit_m$active), 300L)
  inp <- methods::as(nv, "matrix")
  out_m <- methods::as(fit_m$x_clean, "matrix")
  expect_identical(out_m[301:384, ], inp[301:384, ])

  sv <- neuroim2::SparseNeuroVec(t(matrix(arr, 384, 80)[1:300, ]), sp, mvol)
  fit_sv <- fast_phy_denoise(sv)
  expect_s4_class(fit_sv$x_clean, "SparseNeuroVec")
  expect_identical(neuroim2::space(fit_sv$x_clean), sp)
  expect_identical(as.vector(neuroim2::mask(fit_sv$x_clean)), as.vector(mvol))
  expect_equal(sum(fit_sv$active), 300L)
  # A sparse image is analysed exactly like the dense image under its mask.
  expect_equal(methods::as(fit_sv$x_clean, "matrix")[1:300, ], out_m[1:300, ], tolerance = 1e-8)
})

test_that("small or sparsely weighted data do not produce spurious null candidates", {
  fx <- core_fixture(n_vox = 70, n_time = 66, seed = 20)
  a <- fast_phy_denoise(fx$x, ratio_thresh = 1)
  b <- fast_phy_denoise(fx$x, ratio_thresh = 1, svd_engine = "svd")
  # Gram-eigen null directions (~1e-8 * d1) used to pass as extra candidates.
  expect_equal(nrow(a$component_table), nrow(b$component_table))
  expect_equal(a$component_table$selected, b$component_table$selected)
  expect_equal(a$x_clean, b$x_clean, tolerance = 1e-6)
  expect_lte(a$params$rank_used, sum(a$wNN < 1, na.rm = TRUE))
})

test_that("the guardrail can remove every candidate", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 21)
  X0 <- fx$x - rowMeans(fx$x)
  w <- compute_wNN_diff_energy(X0)$wNN
  ext <- extract_caa_oneshot(t(truncated_svd(X0 * (1 - w), 20, return_u = FALSE)$v), 30)

  fit <- fast_phy_denoise(fx$x, design = ext$scores)
  expect_equal(ncol(fit$regressors), 0L)
  expect_equal(nrow(fit$component_table), 0L)
  expect_equal(fit$diagnostics$pass_log[[1]]$reason, "no_components_after_guardrail")
  expect_equal(fit$x_clean, fx$x, tolerance = 1e-10)
})

test_that("results do not depend on the units of the data", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 22)
  ref <- fast_phy_denoise(fx$x, return_diagnostics = FALSE)
  for (s in c(1e-8, 1e6)) {
    fit <- fast_phy_denoise(fx$x * s, return_diagnostics = FALSE)
    # An absolute variance cut in scoring made every ratio NA at small scales.
    expect_equal(fit$component_table$selected, ref$component_table$selected)
    expect_equal(fit$x_clean / s, ref$x_clean, tolerance = 1e-6)
  }
})

test_that("DiCCA/DiPCA extraction runs through the pipeline reproducibly", {
  skip_if_not_installed("dipca")
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 23)
  set.seed(5)
  before <- .Random.seed
  a <- fast_phy_denoise(fx$x, extractor = "dicca", seed = 1, return_diagnostics = FALSE)
  b <- fast_phy_denoise(fx$x, extractor = "dicca", seed = 1, return_diagnostics = FALSE)
  d <- fast_phy_denoise(fx$x, extractor = "dipca", seed = 1, return_diagnostics = FALSE)
  # dipca's random draws used to ignore `seed` and advance the caller's RNG.
  expect_identical(a$x_clean, b$x_clean)
  expect_identical(.Random.seed, before)
  expect_equal(dim(d$x_clean), dim(fx$x))
  expect_true(all(is.finite(d$component_table$predictability)))
})

test_that("invalid arguments give clear errors", {
  x <- core_fixture(n_vox = 100, n_time = 60, seed = 18)$x
  expect_error(fast_phy_denoise(x, tr = NA), "`tr`")
  expect_error(fast_phy_denoise(x, tr = -1), "`tr`")
  expect_error(fast_phy_denoise(x, n_candidates = 0), "`n_candidates`")
  expect_error(fast_phy_denoise(x, max_passes = 0), "`max_passes`")
  expect_error(fast_phy_denoise(x, ratio_thresh = "a"), "`ratio_thresh`")
  expect_error(fast_phy_denoise(x, delta_nt = 74), "`delta_nt`")
  expect_error(fast_phy_denoise(x, pred_thresh = 2), "`pred_thresh`")
  expect_error(fast_phy_denoise(x, use_cpp = NA), "`use_cpp`")
  expect_error(fast_phy_denoise(x, mask = c(0, -1)), "voxel indices")
  expect_error(fast_phy_denoise(x, mask = c(TRUE, NA, rep(TRUE, 98))), "must not contain NA")
  expect_error(fast_phy_denoise(x[, 1:6]), "at least 8 time points")
  expect_error(fast_phy_denoise(x[1:10, ], orientation = "voxels_by_time"), "at least 20 analysed voxels")
  expect_error(fast_phy_denoise(x, design = matrix(1, 10, 1)), "design has 10 rows")
  expect_error(fast_phy_denoise(x, extraction_weight = "grey"), "should be one of")
})

test_that("print method summarises the fit", {
  fx <- core_fixture(n_vox = 300, n_time = 120, seed = 19)
  fit <- fast_phy_denoise(fx$x)
  expect_output(print(fit), "<fast_phy_denoise_result>")
  expect_output(print(fit), "removed:")
  expect_output(res <- withVisible(print(fit)))
  expect_false(res$visible)
  expect_identical(res$value, fit)
})
