test_that("ratio scoring matches a direct computation of median R^2", {
  set.seed(20)
  X <- matrix(rnorm(60 * 50), 60, 50)
  Tc <- matrix(rnorm(50 * 3), 50, 3)
  w <- c(rep(0.1, 20), rep(0.9, 40))
  out <- score_components_nn_nt(X, Tc, w, use_cpp = FALSE)

  r2 <- stats::cor(t(X), Tc)^2
  expect_equal(out$med_nn, apply(r2[1:20, ], 2, stats::median))
  expect_equal(out$med_nt, apply(r2[21:60, ], 2, stats::median))
  expect_equal(out$ratio, out$med_nn / out$med_nt)
  expect_equal(out$r2_mean, colMeans(r2))
})

test_that("constant voxels are ignored rather than counted as R^2 = 0", {
  set.seed(21)
  X <- matrix(rnorm(60 * 50), 60, 50)
  Tc <- matrix(rnorm(50 * 2), 50, 2)
  w <- c(rep(0.1, 20), rep(0.9, 40))
  ref <- score_components_nn_nt(X, Tc, w, use_cpp = FALSE)

  # 60 zero rows labelled neuronal drove med_nt to 0 and the ratio to ~1e9.
  Xz <- rbind(X, matrix(0, 60, 50))
  out <- score_components_nn_nt(Xz, Tc, c(w, rep(1, 60)), use_cpp = FALSE)
  expect_equal(out$ratio, ref$ratio)
})

test_that("scoring does not depend on data or component units", {
  set.seed(25)
  X <- matrix(rnorm(60 * 50), 60, 50)
  Tc <- matrix(rnorm(50 * 3), 50, 3)
  w <- c(rep(0.1, 20), rep(0.9, 40))
  ref <- score_components_nn_nt(X, Tc, w, use_cpp = FALSE)
  for (s in c(1e-8, 1e8)) {
    expect_equal(score_components_nn_nt(X * s, Tc, w, use_cpp = FALSE), ref)
    expect_equal(score_components_nn_nt(X, Tc * s, w, use_cpp = FALSE), ref)
    if (cpp_kernel_loaded()) {
      expect_equal(score_components_nn_nt(X * s, Tc, w, use_cpp = TRUE), ref, tolerance = 1e-10)
    }
  }
})

test_that("degenerate components and missing classes are handled", {
  set.seed(22)
  X <- matrix(rnorm(40 * 60), 40, 60)
  w <- c(rep(0.2, 10), rep(0.8, 30))

  expect_length(score_components_nn_nt(X, matrix(0, 60, 0), w, use_cpp = FALSE)$ratio, 0)
  zero <- score_components_nn_nt(X, matrix(0, 60, 2), w, use_cpp = FALSE)
  expect_true(all(is.na(zero$ratio)))
  expect_error(
    score_components_nn_nt(X, matrix(rnorm(120), 60, 2), rep(1, 40), use_cpp = FALSE),
    "Need both non-neuronal"
  )
})

test_that("compiled and R scoring agree, including edge cases", {
  skip_if_not(cpp_kernel_loaded())
  set.seed(23)
  for (i in 1:20) {
    n <- sample(30:80, 1)
    tt <- sample(20:60, 1)
    X <- matrix(rnorm(n * tt), n, tt)
    X[sample(n, 3), ] <- 0
    Tc <- matrix(rnorm(tt * 4), tt, 4)
    Tc[, 4] <- 0
    w <- stats::runif(n)
    w[1:2] <- c(0.1, 0.9)
    r_out <- score_components_nn_nt(X, Tc, w, use_cpp = FALSE)
    c_out <- score_components_nn_nt(X, Tc, w, use_cpp = TRUE)
    expect_equal(c_out, r_out, tolerance = 1e-10)
  }
})

test_that("the compiled kernel returns plain vectors", {
  skip_if_not(cpp_kernel_loaded())
  set.seed(24)
  out <- score_components_nn_nt_cpp(
    matrix(rnorm(30 * 40), 30, 40), matrix(rnorm(40 * 3), 40, 3),
    c(rep(0.2, 10), rep(0.8, 20)), 1e-12
  )
  expect_null(dim(out$ratio))
  expect_length(out$ratio, 3)
})

test_that("split-half stability rejects noise and keeps reproducible structure", {
  set.seed(5)
  Xn <- matrix(rnorm(400 * 200), 400, 200)
  w <- rep(c(0, 1), c(80, 320))
  sv <- truncated_svd(Xn * w, 20, return_u = FALSE)
  cand <- extract_caa_oneshot(t(sv$v), 10)$scores
  st <- component_stability(Xn, w, cand, 20, "caa", 1L, 10L)
  # The old filter compared each half with itself and passed 6/30 noise components.
  expect_false(any(st$stable))

  fx <- core_fixture(n_vox = 400, n_time = 200, seed = 2)
  X0 <- fx$x - rowMeans(fx$x)
  w <- ifelse(fx$nn, 0, 1)
  sv <- truncated_svd(X0 * (1 - w), 10, return_u = FALSE)
  cand <- extract_caa_oneshot(t(sv$v), 6)$scores
  st <- component_stability(X0, 1 - w, cand, 10, "caa", 1L, 6L)
  physio_like <- apply(cand, 2, function(v) summary(stats::lm(v ~ fx$phys))$r.squared) > 0.8
  expect_true(any(physio_like))
  expect_true(all(st$stable[physio_like]))
})
