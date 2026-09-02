## Large-sample validation for cmbridge.
## Run after installing the package:
##   Rscript inst/simulations/validate_large_sample.R

library(cmbridge)

out_dir <- Sys.getenv("CMBRIDGE_VALIDATION_DIR", unset = "")
if (!nzchar(out_dir)) out_dir <- file.path(tempdir(), "cmbridge_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

poly_truth <- function(v) 1.7 + 0.25 * v + 0.10 * v^2 - 0.08 * v^3
rbf_truth_factory <- function(centers, bandwidth, coef, intercept = 1.6) {
  force(centers); force(bandwidth); force(coef); force(intercept)
  function(v) intercept + as.numeric(rbf_kernel(v, centers, bandwidth) %*% coef)
}

generate_bridge <- function(n, truth, seed, strength = 0.97) {
  set.seed(seed)
  B <- runif(n, -1, 1)
  Vfull <- strength * B + (1 - strength) * runif(n, -1, 1)
  beta <- truth(Vfull)
  stopifnot(min(beta) > 1)
  M <- rbinom(n, 1, 1 / beta)
  Vobs <- Vfull
  Vobs[M == 0] <- NA_real_
  list(B = B, Vfull = Vfull, Vobs = Vobs, M = M, beta = beta)
}

score_fit <- function(fit, truth, grid, method, seed) {
  est <- predict(fit, grid)
  tru <- truth(grid)
  data.frame(
    method = method, seed = seed,
    rmse = sqrt(mean((est - tru)^2)),
    max_abs = max(abs(est - tru)),
    moment_loss = fit$moment_loss
  )
}

seeds <- 1:5
grid <- seq(-1, 1, length.out = 201)
results <- list()
curves <- list()

## 1) Sieve minimum distance. The true bridge is cubic, hence exactly in the
## degree-3 polynomial sieve used for this diagnostic.
for (seed in seeds) {
  dat <- generate_bridge(300000, poly_truth, seed)
  fit <- fit_bridge(
    dat$B, dat$Vobs, dat$M, method = "sieve_md",
    control = list(
      target_basis = "poly", instrument_basis = "poly",
      target_degree = 3L, instrument_degree = 7L,
      lambda = 1e-10, weight_ridge = 1e-8
    )
  )
  results[[length(results) + 1L]] <- score_fit(fit, poly_truth, grid, "sieve_md", seed)
  if (seed == 1) curves[["sieve_md"]] <- cbind(truth = poly_truth(grid), estimate = predict(fit, grid))
}

## 2) Landweber regularization on the same exact polynomial class.
for (seed in seeds) {
  dat <- generate_bridge(300000, poly_truth, seed)
  fit <- fit_bridge(
    dat$B, dat$Vobs, dat$M, method = "landweber",
    control = list(
      target_basis = "poly", instrument_basis = "poly",
      target_degree = 3L, instrument_degree = 7L,
      n_iter = 1500L, step_fraction = 0.95, weight_ridge = 1e-8,
      tol = 0
    )
  )
  results[[length(results) + 1L]] <- score_fit(fit, poly_truth, grid, "landweber", seed)
  if (seed == 1) curves[["landweber"]] <- cbind(truth = poly_truth(grid), estimate = predict(fit, grid))
}

## 3) PMMR. This diagnostic targets the actual PMMR objective of Mastouri et al.
## The true bridge is a finite Gaussian-RKHS expansion with the same target
## kernel centers and bandwidth used by the learner. The implementation uses
## Nyström approximations to both target and instrument kernels, so the
## validation truth lies exactly in the fitted target subspace.
centers <- matrix(seq(-0.8, 0.8, length.out = 5), ncol = 1)
bw <- 0.70
coef <- c(0.55, 0.50, 0.60, 0.50, 0.55)
rbf_truth <- rbf_truth_factory(centers, bw, coef, intercept = 0)
critic_centers <- matrix(seq(-1, 1, length.out = 81), ncol = 1)
for (seed in seeds) {
  dat <- generate_bridge(300000, rbf_truth, seed)
  fit <- fit_bridge(
    dat$B, dat$Vobs, dat$M, method = "pmmr",
    control = list(
      scale = FALSE,
      target_centers_matrix = centers,
      instrument_centers_matrix = critic_centers,
      target_bandwidth = bw,
      instrument_bandwidth = 0.40,
      lambda = 1e-5,
      nystrom_ridge = 1e-10,
      nystrom_tol = 1e-10
    )
  )
  results[[length(results) + 1L]] <- score_fit(fit, rbf_truth, grid, "pmmr", seed)
  if (seed == 1) curves[["pmmr"]] <- cbind(truth = rbf_truth(grid), estimate = predict(fit, grid))
}

res <- do.call(rbind, results)
print(res)
print(aggregate(cbind(rmse, max_abs) ~ method, res, function(x) c(mean = mean(x), max = max(x))))

## These tolerances are intentionally wider than the observed Monte Carlo error
## and are meant to catch implementation failures, not benchmark efficiency.
stopifnot(all(res$rmse < 0.03))
stopifnot(all(res$max_abs < 0.07))

utils::write.csv(res, file.path(out_dir, "validation_results.csv"), row.names = FALSE)
pdf(file.path(out_dir, "truth_vs_estimate.pdf"), width = 6, height = 6)
for (nm in names(curves)) {
  cc <- curves[[nm]]
  plot(cc[, "truth"], cc[, "estimate"], pch = 19, cex = 0.45,
       xlab = "True bridge", ylab = "Estimated bridge", main = nm)
  abline(0, 1, lty = 2)
}
dev.off()

cat("Validation passed. Outputs written to:", out_dir, "\n")
