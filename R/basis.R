#' Polynomial basis
#'
#' @param x Numeric vector or matrix.
#' @param degree Maximum marginal polynomial degree.
#' @param intercept Include an intercept column.
#' @return Numeric design matrix.
#' @export
poly_basis <- function(x, degree = 3L, intercept = TRUE) {
  x <- .as_matrix(x)
  degree <- as.integer(degree)
  if (degree < 1L) stop("degree must be at least 1.", call. = FALSE)
  out <- if (isTRUE(intercept)) matrix(1, nrow(x), 1L) else NULL
  for (j in seq_len(ncol(x))) {
    for (k in seq_len(degree)) out <- cbind(out, x[, j]^k)
  }
  colnames(out) <- NULL
  out
}

.fit_basis_spec <- function(x, type = c("bs", "poly"), degree = 3L, df = 6L) {
  type <- match.arg(type)
  x <- .as_matrix(x)
  if (anyNA(x)) stop("basis fitting data cannot contain missing values.", call. = FALSE)
  if (type == "poly") {
    return(list(type = type, degree = as.integer(degree), p = ncol(x)))
  }
  specs <- vector("list", ncol(x))
  for (j in seq_len(ncol(x))) {
    b <- splines::bs(x[, j], df = df, degree = degree, intercept = FALSE)
    specs[[j]] <- list(
      knots = attr(b, "knots"),
      Boundary.knots = attr(b, "Boundary.knots"),
      degree = attr(b, "degree")
    )
  }
  list(type = type, specs = specs, p = ncol(x))
}

.eval_basis_spec <- function(x, spec) {
  x <- .as_matrix(x)
  if (ncol(x) != spec$p) stop("new data have the wrong number of columns.", call. = FALSE)
  if (spec$type == "poly") return(poly_basis(x, degree = spec$degree, intercept = TRUE))
  out <- matrix(1, nrow(x), 1L)
  for (j in seq_len(ncol(x))) {
    sj <- spec$specs[[j]]
    bj <- splines::bs(
      x[, j], knots = sj$knots, Boundary.knots = sj$Boundary.knots,
      degree = sj$degree, intercept = FALSE
    )
    out <- cbind(out, bj)
  }
  out
}

#' Gaussian radial-basis kernel
#'
#' @param x Numeric vector or matrix.
#' @param centers Numeric vector or matrix of kernel centers.
#' @param bandwidth Positive Gaussian bandwidth.
#' @return Kernel design matrix.
#' @export
rbf_kernel <- function(x, centers, bandwidth) {
  x <- .as_matrix(x)
  centers <- .as_matrix(centers)
  if (ncol(x) != ncol(centers)) stop("x and centers must have the same number of columns.", call. = FALSE)
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("bandwidth must be positive.", call. = FALSE)
  x2 <- rowSums(x^2)
  c2 <- rowSums(centers^2)
  d2 <- outer(x2, c2, "+") - 2 * tcrossprod(x, centers)
  d2[d2 < 0] <- 0
  exp(-d2 / (2 * bandwidth^2))
}
