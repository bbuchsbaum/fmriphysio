# Utilities for matrix orientation, centering, and projection.

as_nt_matrix <- function(x, mask = NULL) {
  original <- x
  source_type <- class(x)[1]

  if (is.matrix(x)) {
    X <- x
  } else if (is.numeric(x) && is.array(x) && length(dim(x)) == 2L) {
    X <- unclass(x)
  } else {
    if (!is.null(getS3method("as.matrix", class(x)[1], optional = TRUE)) ||
        methods::hasMethod("as.matrix", class(x)[1])) {
      X <- as.matrix(x)
    } else {
      stop("Unsupported input type for x. Provide a matrix-like object.")
    }
  }

  if (!is.numeric(X)) {
    stop("Input data must be numeric.")
  }

  transposed <- FALSE
  # fMRI convention is usually voxels >> time; if not, treat as T x N.
  if (nrow(X) < ncol(X)) {
    X <- t(X)
    transposed <- TRUE
  }

  if (!is.null(mask)) {
    if (is.logical(mask)) {
      if (length(mask) != nrow(X)) {
        stop("Logical mask length must match the voxel dimension of x.")
      }
      X <- X[mask, , drop = FALSE]
    } else if (is.numeric(mask)) {
      X <- X[mask, , drop = FALSE]
    } else {
      stop("mask must be NULL, logical, or numeric indices.")
    }
  }

  list(
    X = X,
    info = list(
      source_type = source_type,
      transposed_input = transposed,
      original_dim = dim(original),
      matrix_dim = dim(X)
    )
  )
}

restore_orientation <- function(X_nt, info) {
  if (isTRUE(info$transposed_input)) {
    return(t(X_nt))
  }
  X_nt
}

center_rows <- function(X) {
  X - rowMeans(X)
}

safe_scale_vec <- function(x, eps = 1e-12) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < eps) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x)) / s)
}

orthogonalize_to_design <- function(Tcomp, design, eps = 1e-10) {
  if (is.null(design) || ncol(Tcomp) == 0L) {
    return(Tcomp)
  }

  D <- as.matrix(design)
  if (!is.numeric(D)) stop("design must be numeric")
  if (nrow(D) != nrow(Tcomp)) {
    stop("design rows must match number of time points")
  }

  qd <- qr(D)
  rank_d <- qd$rank
  if (rank_d <= 0L) {
    return(Tcomp)
  }

  Q <- qr.Q(qd)[, seq_len(rank_d), drop = FALSE]
  Tproj <- Tcomp - Q %*% (crossprod(Q, Tcomp))

  keep <- apply(Tproj, 2, function(v) stats::sd(v) > eps)
  Tproj[, keep, drop = FALSE]
}

project_out_components <- function(X_nt, Tcomp_tl) {
  if (ncol(Tcomp_tl) == 0L) {
    return(X_nt)
  }

  q <- qr.Q(qr(Tcomp_tl))
  X_nt - (X_nt %*% q) %*% t(q)
}
