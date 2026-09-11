test_that("as_nt_matrix only transposes when told to, or says so", {
  X <- matrix(rnorm(30 * 10), 30, 10)

  p <- as_nt_matrix(X)
  expect_false(p$info$transposed_input)
  expect_equal(dim(p$X), c(30L, 10L))

  expect_message(p2 <- as_nt_matrix(t(X)), "time x voxels")
  expect_true(p2$info$transposed_input)
  expect_equal(p2$X, X)

  # A 10-voxel x 30-TR matrix used to be transposed silently.
  expect_no_message(p3 <- as_nt_matrix(t(X), orientation = "voxels_by_time"))
  expect_equal(dim(p3$X), c(10L, 30L))

  p4 <- as_nt_matrix(t(X), orientation = "time_by_voxels")
  expect_equal(p4$X, X)
})

test_that("as_nt_matrix accepts arrays, data frames, and S4 image objects", {
  a <- array(rnorm(2 * 3 * 4 * 5), c(2, 3, 4, 5))
  pa <- as_nt_matrix(a)
  expect_true(pa$info$array_input)
  expect_equal(pa$X, matrix(a, 24, 5))
  expect_error(as_nt_matrix(a, orientation = "time_by_voxels"), "last dimension")

  X <- matrix(rnorm(40), 10, 4)
  expect_equal(unname(as_nt_matrix(as.data.frame(X))$X), X)

  expect_error(as_nt_matrix(matrix(letters[1:20], 4, 5)), "must be numeric")
  expect_error(as_nt_matrix(list(1, 2)), "Unsupported input type")

  arr <- array(rnorm(4 * 4 * 3 * 20), c(4, 4, 3, 20))
  sp <- neuroim2::NeuroSpace(c(4, 4, 3, 20), spacing = c(2, 2, 2))
  nv <- neuroim2::DenseNeuroVec(arr, sp)
  # base::as.matrix() gave a 20 x 1 matrix for NeuroVec input.
  pn <- as_nt_matrix(nv)
  expect_equal(dim(pn$X), c(48L, 20L))
  expect_equal(unname(pn$X), matrix(arr, 48, 20))
  expect_false(pn$info$neuro$sparse)
})

test_that("restore_output rebuilds dense and sparse neuroim2 images", {
  arr <- array(rnorm(4 * 4 * 3 * 20), c(4, 4, 3, 20))
  sp <- neuroim2::NeuroSpace(c(4, 4, 3, 20), spacing = c(2, 2, 2))
  nv <- neuroim2::DenseNeuroVec(arr, sp)
  pn <- as_nt_matrix(nv)
  back <- restore_output(pn$X, pn, rep(TRUE, 48))
  expect_s4_class(back, "DenseNeuroVec")
  expect_identical(neuroim2::space(back), sp)
  expect_equal(as.vector(back@.Data), as.vector(arr))

  mk <- neuroim2::LogicalNeuroVol(
    array(rep(c(TRUE, FALSE), c(30, 18)), c(4, 4, 3)),
    neuroim2::NeuroSpace(c(4, 4, 3), spacing = c(2, 2, 2))
  )
  sv <- neuroim2::SparseNeuroVec(t(matrix(arr, 48, 20)[1:30, ]), sp, mk)
  ps <- as_nt_matrix(sv)
  expect_true(ps$info$neuro$sparse)
  # Voxels outside a sparse image's mask carry no data and are not analysed.
  expect_equal(ps$keep, rep(c(TRUE, FALSE), c(30, 18)))
  back_s <- restore_output(ps$X[ps$keep, ], ps, ps$keep)
  expect_s4_class(back_s, "SparseNeuroVec")
  expect_identical(neuroim2::space(back_s), sp)
  expect_equal(back_s@data, sv@data)
})

test_that("normalize_mask accepts the documented forms and rejects bad masks", {
  expect_equal(normalize_mask(NULL, 3), rep(TRUE, 3))
  expect_equal(normalize_mask(c(1, 3), 4), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(normalize_mask(c(TRUE, FALSE, TRUE), 3), c(TRUE, FALSE, TRUE))
  expect_equal(normalize_mask(array(c(1, 0, 1, 1), c(2, 2)), 4), c(TRUE, FALSE, TRUE, TRUE))

  expect_error(normalize_mask(c(TRUE, NA, TRUE), 3), "must not contain NA")
  expect_error(normalize_mask(rep(TRUE, 5), 3), "length 5 but x has 3")
  expect_error(normalize_mask(rep(FALSE, 3), 3), "selects no voxels")
  # Non-positive indices used to drop voxels silently.
  expect_error(normalize_mask(c(0, -1), 4), "voxel indices in 1..4")
  expect_error(normalize_mask(c(1, 9), 4), "voxel indices")
  expect_error(normalize_mask(2.5, 4), "whole-number")
  expect_error(normalize_mask(c(1, NA), 4), "voxel indices")
  expect_error(normalize_mask("a", 4), "must be NULL")
})

test_that("analysis_voxels excludes masked, constant, and non-finite voxels", {
  X <- rbind(rnorm(10), rep(0, 10), rep(1000, 10), c(NA, rnorm(9)), c(Inf, rnorm(9)), rnorm(10))
  keep <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  expect_warning(v <- analysis_voxels(X, keep), "Excluded 2 voxel")
  expect_equal(v$active, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(v$n_constant, 2L)
  expect_equal(v$n_nonfinite, 2L)
})

test_that("restore_output writes rows back and restores shape and orientation", {
  X <- matrix(rnorm(12 * 5), 12, 5)
  p <- suppressMessages(as_nt_matrix(t(X)))
  active <- rep(c(TRUE, FALSE), 6)
  out <- restore_output(p$X[active, ] * 0, p, active)
  expect_equal(dim(out), dim(t(X)))
  expect_equal(out[, !active], t(X)[, !active])
  expect_true(all(out[, active] == 0))

  a <- array(rnorm(2 * 3 * 2 * 5), c(2, 3, 2, 5), dimnames = list(NULL, letters[1:3], NULL, NULL))
  pa <- as_nt_matrix(a)
  expect_equal(restore_output(pa$X, pa, rep(TRUE, 12)), a)

  expect_equal(expand_voxel_vector(c(5, 6), c(FALSE, TRUE, TRUE)), c(NA, 5, 6))
})

test_that("orthogonalize_to_design removes the design span and reports kept columns", {
  set.seed(1)
  D <- matrix(rnorm(50 * 2), 50, 2)
  Tc <- cbind(matrix(rnorm(50 * 2), 50, 2), 2 * D[, 1] + 3)
  out <- orthogonalize_to_design(Tc, D)

  expect_equal(attr(out, "kept"), 1:2)
  expect_lt(max(abs(crossprod(cbind(1, D), out))), 1e-10)
  # Components come back centred because an intercept is always included.
  expect_equal(unname(colMeans(out)), c(0, 0), tolerance = 1e-12)

  expect_equal(attr(orthogonalize_to_design(Tc, NULL), "kept"), 1:3)
  expect_error(orthogonalize_to_design(Tc, D[1:10, ]), "design has 10 rows but the data have 50")
  D[1, 1] <- NA
  expect_error(orthogonalize_to_design(Tc, D), "missing or infinite")
})

test_that("project_out_components removes exactly the span of the components", {
  set.seed(2)
  a <- rnorm(40)
  Y <- matrix(rnorm(10 * 40), 10, 40)
  p1 <- project_out_components(Y, cbind(a))
  # Duplicate columns used to remove a rank-2 (or rank-3) subspace.
  expect_equal(project_out_components(Y, cbind(a, a, 2 * a)), p1)
  expect_lt(max(abs(p1 %*% a)), 1e-10)
  expect_equal(project_out_components(p1, cbind(a)), p1)
  expect_equal(project_out_components(Y, matrix(0, 40, 0)), Y)
  expect_equal(project_out_components(Y, matrix(0, 40, 2)), Y)
})

test_that("with_seed is reproducible and restores the caller's RNG state", {
  set.seed(1)
  before <- .Random.seed
  x1 <- with_seed(5, stats::runif(3))
  expect_identical(.Random.seed, before)
  expect_identical(with_seed(5, stats::runif(3)), x1)
  expect_identical(with_seed(NULL, 42), 42)

  saved <- .Random.seed
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  rm(".Random.seed", envir = globalenv())
  with_seed(1, stats::runif(1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("truncated_svd engines agree and signs are deterministic", {
  set.seed(3)
  L <- matrix(rnorm(200 * 5), 200) %*% diag(c(20, 15, 10, 8, 6)) %*% matrix(rnorm(5 * 40), 5)
  X <- L + matrix(rnorm(200 * 40), 200)

  g <- truncated_svd(X, 5)
  s <- truncated_svd(X, 5, svd_engine = "svd")
  r <- truncated_svd(X, 5, svd_engine = "rsvd", seed = 1)
  expect_identical(g$engine, "gram")
  expect_equal(g$d, s$d, tolerance = 1e-8)
  expect_equal(g$v, s$v, tolerance = 1e-6)
  expect_equal(g$u, s$u, tolerance = 1e-6)
  expect_gt(min_principal_cos(g$v, r$v), 0.999)
  expect_true(all(apply(g$v, 2, function(v) v[which.max(abs(v))] > 0)))

  wide <- truncated_svd(t(X), 5)
  expect_equal(wide$d, s$d, tolerance = 1e-8)
  expect_equal(abs(wide$u), abs(s$v), tolerance = 1e-6)
  expect_equal(abs(wide$v), abs(s$u), tolerance = 1e-6)
  expect_equal(crossprod(wide$v), diag(5), tolerance = 1e-8)

  expect_null(truncated_svd(X, 5, return_u = FALSE)$u)
  expect_length(truncated_svd(X, 100)$d, 40)
  expect_error(truncated_svd(X, 0), "rank must be >= 1")
})

test_that("randomized SVD with a seed does not touch the caller's RNG", {
  set.seed(4)
  before <- .Random.seed
  X <- matrix(rnorm(300 * 30), 300)
  set.seed(4)
  before <- .Random.seed
  a <- truncated_svd(X, 5, svd_engine = "rsvd", seed = 9)
  b <- truncated_svd(X, 5, svd_engine = "rsvd", seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(a$v, b$v)
})

test_that("argument validators give clear messages", {
  expect_identical(check_count(3, "n"), 3L)
  expect_error(check_count(0, "n"), "`n` must be a whole number >= 1")
  expect_error(check_count(2.5, "n"), "whole number")
  expect_error(check_count(NA, "n"), "whole number")
  expect_equal(check_number(0.5, "p", 0, 1), 0.5)
  expect_null(check_number(NULL, "p", null_ok = TRUE))
  expect_error(check_number(2, "p", 0, 1), "in \\[0, 1\\]")
  expect_error(check_number("a", "p"), "`p` must be a finite number")
  expect_error(check_number(-1, "p", lower = 0), ">= 0")
  expect_error(check_flag(NA, "f"), "TRUE or FALSE")
})
