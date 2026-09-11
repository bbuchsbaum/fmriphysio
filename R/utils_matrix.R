# Utilities for input parsing, voxel selection, orientation, and projection.

# Parse matrix-like input into a full voxels x time matrix plus a logical
# voxel-inclusion vector. The matrix is not subset, so denoised rows can be
# written back and the output keeps the input's shape.
as_nt_matrix <- function(
  x,
  mask = NULL,
  orientation = c("auto", "voxels_by_time", "time_by_voxels")
) {
  orientation <- match.arg(orientation)
  source_type <- class(x)[1]
  original_dim <- dim(x)
  original_dimnames <- dimnames(x)
  array_input <- FALSE
  neuro <- NULL

  if (methods::is(x, "NeuroVec")) {
    # neuroim2 images: voxels x time over the full grid (zeros outside the
    # mask of a sparse image). The space, label, and mask are kept so the
    # denoised data can be returned as the same kind of image.
    if (orientation == "time_by_voxels") {
      stop("neuroim2 images are voxels x time; `orientation` must be \"auto\" or \"voxels_by_time\".", call. = FALSE)
    }
    orientation <- "voxels_by_time"
    X <- methods::as(x, "matrix")
    sparse <- methods::is(x, "AbstractSparseNeuroVec")
    neuro <- list(
      sparse = sparse,
      space = neuroim2::space(x),
      mask = if (sparse) neuroim2::mask(x) else NULL,
      label = if ("label" %in% methods::slotNames(x)) x@label else NULL
    )
  } else if (is.matrix(x)) {
    X <- x
  } else if (is.array(x) && length(dim(x)) >= 3L) {
    # Spatial dimensions first, time last (e.g. a 4D volume series).
    d <- dim(x)
    X <- matrix(as.vector(x), nrow = prod(d[-length(d)]), ncol = d[length(d)])
    array_input <- TRUE
    if (orientation == "time_by_voxels") {
      stop("Array input must have time as its last dimension.", call. = FALSE)
    }
    orientation <- "voxels_by_time"
  } else if (isS4(x)) {
    # S4 containers such as neuroim2::NeuroVec define a coercion to a
    # voxels x time matrix; base::as.matrix() would not dispatch to it.
    X <- tryCatch(methods::as(x, "matrix"), error = function(e) NULL)
    if (is.null(X)) {
      stop(sprintf("Cannot convert an object of class '%s' to a matrix.", source_type), call. = FALSE)
    }
  } else if (is.data.frame(x) || length(dim(x)) == 2L) {
    X <- as.matrix(x)
  } else {
    stop("Unsupported input type for x. Provide a matrix, 3D/4D array, or matrix-like object.", call. = FALSE)
  }

  if (!is.numeric(X)) {
    stop("Input data must be numeric.", call. = FALSE)
  }

  transposed <- FALSE
  if (orientation == "time_by_voxels") {
    X <- t(X)
    transposed <- TRUE
  } else if (orientation == "auto" && nrow(X) < ncol(X)) {
    message(
      "x has fewer rows than columns; treating it as time x voxels. ",
      "Set `orientation` explicitly to silence this message."
    )
    X <- t(X)
    transposed <- TRUE
  }

  keep <- normalize_mask(mask, nrow(X))
  if (!is.null(neuro$mask)) {
    # Only voxels stored in a sparse image carry data.
    keep <- keep & as.logical(as.vector(neuro$mask))
  }

  list(
    X = X,
    keep = keep,
    info = list(
      source_type = source_type,
      transposed_input = transposed,
      array_input = array_input,
      original_dim = original_dim,
      original_dimnames = original_dimnames,
      n_voxels = nrow(X),
      n_time = ncol(X),
      neuro = neuro
    )
  )
}

# Convert a mask (logical vector, voxel indices, 3D array, or an S4 volume such
# as neuroim2::LogicalNeuroVol) into a logical vector of length `n`.
normalize_mask <- function(mask, n) {
  if (is.null(mask)) {
    return(rep(TRUE, n))
  }

  spatial <- isS4(mask) || is.array(mask)
  if (spatial) {
    mask <- as.vector(mask)
    if (is.numeric(mask) && length(mask) == n && all(mask %in% c(0, 1))) {
      mask <- as.logical(mask)
    }
  }

  if (is.logical(mask)) {
    if (length(mask) != n) {
      stop(sprintf("Logical mask has length %d but x has %d voxels.", length(mask), n), call. = FALSE)
    }
    if (anyNA(mask)) {
      stop("mask must not contain NA.", call. = FALSE)
    }
    if (!any(mask)) {
      stop("mask selects no voxels.", call. = FALSE)
    }
    return(mask)
  }

  if (is.numeric(mask)) {
    if (!length(mask) || anyNA(mask) || any(mask != round(mask)) || any(mask < 1) || any(mask > n)) {
      stop(sprintf("Numeric mask must contain whole-number voxel indices in 1..%d.", n), call. = FALSE)
    }
    keep <- rep(FALSE, n)
    keep[mask] <- TRUE
    return(keep)
  }

  stop("mask must be NULL, a logical vector, or numeric voxel indices.", call. = FALSE)
}

# Voxels that enter the analysis: inside the mask, all values finite, and not
# constant over time. Excluded voxels are returned unchanged by the denoisers.
analysis_voxels <- function(X, keep) {
  mu <- rowMeans(X)
  ss_dev <- rowSums((X - mu)^2)
  finite <- is.finite(ss_dev)
  constant <- finite & (ss_dev <= 1e-12 * rowSums(X * X))

  n_nonfinite <- sum(keep & !finite)
  if (n_nonfinite > 0L) {
    warning(sprintf(
      "Excluded %d voxel(s) containing non-finite values; they are returned unchanged.",
      n_nonfinite
    ), call. = FALSE)
  }

  active <- keep & finite & !constant
  list(
    active = active,
    n_nonfinite = n_nonfinite,
    n_constant = sum(keep & constant)
  )
}

# Write denoised rows back into the full matrix and restore the input's
# orientation, array shape, or neuroim2 image class.
restore_output <- function(X_active, parsed, active) {
  X <- parsed$X
  X[active, ] <- X_active
  info <- parsed$info
  if (!is.null(info$neuro)) {
    return(rebuild_neurovec(X, info$neuro))
  }
  if (isTRUE(info$transposed_input)) {
    X <- t(X)
  }
  if (isTRUE(info$array_input)) {
    X <- array(X, dim = info$original_dim, dimnames = info$original_dimnames)
  }
  X
}

# Rebuild a neuroim2 image from a full-grid voxels x time matrix: a
# SparseNeuroVec with the original mask for sparse input, otherwise a
# DenseNeuroVec, always in the original NeuroSpace.
rebuild_neurovec <- function(X, neuro) {
  sp <- neuro$space
  out <- if (isTRUE(neuro$sparse)) {
    idx <- which(as.logical(as.vector(neuro$mask)))
    neuroim2::SparseNeuroVec(t(X[idx, , drop = FALSE]), sp, neuro$mask)
  } else {
    neuroim2::DenseNeuroVec(array(X, dim = dim(sp)), sp)
  }
  if (!is.null(neuro$label) && "label" %in% methods::slotNames(out)) {
    out@label <- neuro$label
  }
  out
}

# Expand a per-active-voxel vector to full voxel length, NA elsewhere.
expand_voxel_vector <- function(v, active) {
  out <- rep(NA_real_, length(active))
  out[active] <- v
  out
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

# Remove the span of the design (plus an intercept) from candidate nuisance
# time courses. Columns that are (numerically) fully explained by the design
# are dropped; the indices of retained columns are returned as attr "kept".
orthogonalize_to_design <- function(Tcomp, design, tol = 1e-6) {
  L <- ncol(Tcomp)
  if (is.null(design) || L == 0L) {
    return(structure(Tcomp, kept = seq_len(L)))
  }

  D <- as.matrix(design)
  if (!is.numeric(D)) stop("design must be numeric.", call. = FALSE)
  if (nrow(D) != nrow(Tcomp)) {
    stop(sprintf(
      "design has %d rows but the data have %d time points.",
      nrow(D), nrow(Tcomp)
    ), call. = FALSE)
  }
  if (any(!is.finite(D))) stop("design must not contain missing or infinite values.", call. = FALSE)

  qd <- qr(cbind(1, D))
  Q <- qr.Q(qd)[, seq_len(qd$rank), drop = FALSE]
  Tproj <- Tcomp - Q %*% crossprod(Q, Tcomp)

  norm_in <- sqrt(colSums(Tcomp * Tcomp))
  norm_out <- sqrt(colSums(Tproj * Tproj))
  kept <- which(norm_out > tol * pmax(norm_in, .Machine$double.eps))
  out <- Tproj[, kept, drop = FALSE]
  attr(out, "kept") <- kept
  out
}

# Orthogonal projection of each voxel time series onto the complement of the
# column space of `Tcomp_tl` (time x components). Rank-deficient component
# sets remove only their actual span.
project_out_components <- function(X_nt, Tcomp_tl, tol = 1e-7) {
  if (ncol(Tcomp_tl) == 0L) {
    return(X_nt)
  }
  qd <- qr(Tcomp_tl, tol = tol)
  if (qd$rank == 0L) {
    return(X_nt)
  }
  q <- qr.Q(qd)[, seq_len(qd$rank), drop = FALSE]
  X_nt - (X_nt %*% q) %*% t(q)
}

# Evaluate `code` with a fixed RNG seed, then restore the caller's RNG state.
with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  code
}

# Argument validators ---------------------------------------------------------

check_count <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) || x < min) {
    stop(sprintf("`%s` must be a whole number >= %d.", name, min), call. = FALSE)
  }
  as.integer(x)
}

check_number <- function(x, name, lower = -Inf, upper = Inf, null_ok = FALSE) {
  if (is.null(x) && null_ok) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    range_txt <- if (is.finite(lower) && is.finite(upper)) {
      sprintf(" in [%g, %g]", lower, upper)
    } else if (is.finite(lower)) {
      sprintf(" >= %g", lower)
    } else if (is.finite(upper)) {
      sprintf(" <= %g", upper)
    } else {
      ""
    }
    stop(sprintf("`%s` must be a finite number%s.", name, range_txt), call. = FALSE)
  }
  as.numeric(x)
}

check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a positive finite number.", name), call. = FALSE)
  }
  as.numeric(x)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}
