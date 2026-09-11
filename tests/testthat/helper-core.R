# Shared fixtures for core-algorithm tests.

# Voxels x time data with known structure. The first `nn_frac` of voxels are
# non-neuronal: they carry strong aliased cardiac (near Nyquist) and
# respiratory signals. The remaining voxels are neuronal: they carry smooth
# neural signals and weak physiology. Everything is returned so tests can
# measure what was removed.
core_fixture <- function(n_vox = 600, n_time = 200, seed = 1, nn_frac = 0.2,
                         nn_amp = 2, nt_phys_amp = 0.3, baseline = 0) {
  set.seed(seed)
  nn <- seq_len(n_vox) <= nn_frac * n_vox
  card <- sin(2 * pi * cumsum(0.42 + 0.01 * rnorm(n_time)))
  resp <- sin(2 * pi * cumsum(0.18 + 0.005 * rnorm(n_time)))
  phys <- cbind(card = card, resp = resp)
  neur <- apply(matrix(rnorm(n_time * 3), n_time), 2, function(e) {
    as.numeric(stats::filter(e, rep(1 / 8, 8), circular = TRUE))
  })
  neur <- scale(neur)
  Wp <- matrix(rnorm(n_vox * 2), n_vox) * ifelse(nn, nn_amp, nt_phys_amp)
  Wn <- matrix(rnorm(n_vox * 3), n_vox) * ifelse(nn, 0.1, 0.8)
  x <- Wp %*% t(phys) + Wn %*% t(neur) + matrix(rnorm(n_vox * n_time), n_vox) + baseline
  list(x = x, nn = nn, phys = phys, neur = unclass(neur))
}

# Energy of the (row-centred) data in `rows` that lies in the span of `S`.
subspace_energy <- function(X, S, rows = seq_len(nrow(X))) {
  Xc <- X - rowMeans(X)
  q <- qr.Q(qr(S))
  sum((Xc[rows, , drop = FALSE] %*% q)^2)
}

# Smallest cosine of the principal angles between two column spaces
# (1 = identical subspaces).
min_principal_cos <- function(A, B) {
  qa <- qr.Q(qr(A))
  qb <- qr.Q(qr(B))
  min(svd(crossprod(qa, qb))$d)
}

cpp_kernel_loaded <- function() {
  is.loaded("_fmriphysio_score_components_nn_nt_cpp")
}
