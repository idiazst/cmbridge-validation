.fit_sieve_md <- function(y, d, x, z, control) {
  ctrl <- .merge_control(list(
    target_basis = "bs", instrument_basis = "bs",
    target_degree = 3L, instrument_degree = 3L,
    target_df = 6L, instrument_df = 9L,
    lambda = 1e-8, weight_ridge = 1e-8
  ), control)

  active <- abs(d) > 0
  if (!any(active)) stop("diagonal is zero for every observation.", call. = FALSE)
  if (anyNA(x[active, , drop = FALSE])) stop("target may be missing only where diagonal is zero.", call. = FALSE)
  if (anyNA(z)) stop("instrument cannot be missing for sieve minimum distance.", call. = FALSE)

  xspec <- .fit_basis_spec(x[active, , drop = FALSE], ctrl$target_basis,
                           ctrl$target_degree, ctrl$target_df)
  zspec <- .fit_basis_spec(z, ctrl$instrument_basis,
                           ctrl$instrument_degree, ctrl$instrument_df)
  H <- matrix(0, nrow(x), ncol(.eval_basis_spec(x[active, , drop = FALSE], xspec)))
  H[active, ] <- .eval_basis_spec(x[active, , drop = FALSE], xspec)
  Q <- .eval_basis_spec(z, zspec)

  n <- length(y)
  A <- crossprod(Q, d * H) / n
  cvec <- crossprod(Q, y) / n
  W <- .safe_inverse(crossprod(Q) / n + ctrl$weight_ridge * diag(ncol(Q)), ridge = ctrl$weight_ridge)
  P <- diag(ncol(H)); P[1L, 1L] <- 0
  lhs <- crossprod(A, W %*% A) + ctrl$lambda * P
  rhs <- crossprod(A, W %*% cvec)
  coef <- as.numeric(.safe_solve(lhs, rhs))

  pred_train <- as.numeric(H %*% coef)
  pred_train[!active] <- NA_real_
  residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  moment <- as.numeric(crossprod(Q, residual) / n)

  list(
    coefficients = coef,
    target_spec = xspec,
    instrument_spec = zspec,
    target_design_cols = ncol(H),
    fitted = pred_train,
    residual = residual,
    moment_loss = as.numeric(crossprod(moment, W %*% moment)),
    tuning = ctrl,
    predict_fun = function(newx) {
      newx <- .as_matrix(newx)
      as.numeric(.eval_basis_spec(newx, xspec) %*% coef)
    },
    moment_weight = W,
    instrument_features = function(newz) .eval_basis_spec(newz, zspec)
  )
}
