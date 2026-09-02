## Large-sample validation for the adjoint learners in cmbridge.
## Run after installing the package:
##   Rscript inst/simulations/validate_adjoint_large_sample.R

library(cmbridge)

out_dir <- Sys.getenv("CMBRIDGE_ADJOINT_VALIDATION_DIR", unset = "")
if (!nzchar(out_dir)) out_dir <- file.path(tempdir(), "cmbridge_adjoint_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## We test the adjoint equation
##
##   E[M {lambda(B) - phi(V)} | V] = 0.
##
## The design deliberately makes M depend on B even conditional on V.
##
## Let S in {-1,+1}, B = V + delta*S, and
##   P(M=1 | V,S=+1) = p_plus,
##   P(M=1 | V,S=-1) = p_minus.
##
## Since P(S=+1|V)=P(S=-1|V)=1/2,
##
##   P(S=+1 | V,M=1) = p_plus/(p_plus+p_minus).
##
## Therefore defining
##
##   phi(V) = E[lambda(B) | V,M=1]
##
## makes the adjoint moment restriction hold exactly.

delta <- 0.15
p_plus <- 0.80
p_minus <- 0.40
w_plus <- p_plus / (p_plus + p_minus)
w_minus <- p_minus / (p_plus + p_minus)

generate_adjoint <- function(n, lambda_truth, seed) {
  set.seed(seed)
  V <- runif(n, -0.75, 0.75)
  S <- sample(c(-1, 1), n, replace = TRUE)
  B <- V + delta * S

  pM <- ifelse(S == 1, p_plus, p_minus)
  M <- rbinom(n, 1, pM)

  phi <- w_plus * lambda_truth(V + delta) +
    w_minus * lambda_truth(V - delta)

  list(B = B, V = V, M = M, phi = phi, S = S)
}

score_fit <- function(fit, truth, grid, method, seed) {
  est <- predict(fit, grid)
  tru <- truth(grid)

  data.frame(
    method = method,
    seed = seed,
    rmse = sqrt(mean((est - tru)^2)),
    max_abs = max(abs(est - tru)),
    moment_loss = fit$moment_loss
  )
}

summarize_results <- function(res) {
  do.call(rbind, lapply(split(res, res$method), function(x) {
    data.frame(
      method = x$method[[1L]],
      mean_rmse = mean(x$rmse),
      max_rmse = max(x$rmse),
      mean_max_abs = mean(x$max_abs),
      max_abs = max(x$max_abs),
      mean_moment_loss = mean(x$moment_loss)
    )
  }))
}

seeds <- 1:5
grid <- seq(-0.90, 0.90, length.out = 201)
results <- list()
curves <- list()

## 1) Sieve minimum distance.
## The true adjoint is cubic, so it lies exactly in the degree-3 target sieve.
lambda_poly <- function(b) {
  1.2 + 0.40 * b + 0.15 * b^2 - 0.08 * b^3
}

for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_poly, seed)

  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi,
    method = "sieve_md",
    control = list(
      target_basis = "poly",
      instrument_basis = "poly",
      target_degree = 3L,
      instrument_degree = 7L,
      lambda = 1e-10,
      weight_ridge = 1e-8
    )
  )

  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_poly, grid, "sieve_md", seed)

  if (seed == 1L) {
    curves[["sieve_md"]] <- data.frame(
      truth = lambda_poly(grid),
      estimate = predict(fit, grid)
    )
  }
}

## 2) Landweber regularization.
## Same exact polynomial target class.
for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_poly, seed)

  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi,
    method = "landweber",
    control = list(
      target_basis = "poly",
      instrument_basis = "poly",
      target_degree = 3L,
      instrument_degree = 7L,
      n_iter = 1500L,
      step_fraction = 0.95,
      weight_ridge = 1e-8,
      tol = 0
    )
  )

  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_poly, grid, "landweber", seed)

  if (seed == 1L) {
    curves[["landweber"]] <- data.frame(
      truth = lambda_poly(grid),
      estimate = predict(fit, grid)
    )
  }
}

## 3) PMMR.
## The true adjoint is exactly a finite Gaussian-RKHS expansion using the same
## target kernel centers and bandwidth supplied to the learner.

rbf_truth_factory <- function(centers, bandwidth, coef) {
  force(centers)
  force(bandwidth)
  force(coef)

  function(b) {
    as.numeric(rbf_kernel(b, centers, bandwidth) %*% coef)
  }
}

target_centers <- matrix(seq(-0.75, 0.75, length.out = 5), ncol = 1)
target_bw <- 0.65
target_coef <- c(0.50, 0.62, 0.55, 0.60, 0.48)
lambda_rbf <- rbf_truth_factory(target_centers, target_bw, target_coef)

critic_centers <- matrix(seq(-0.75, 0.75, length.out = 81), ncol = 1)

for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_rbf, seed)

  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi,
    method = "pmmr",
    control = list(
      scale = FALSE,
      target_centers_matrix = target_centers,
      instrument_centers_matrix = critic_centers,
      target_bandwidth = target_bw,
      instrument_bandwidth = 0.35,
      lambda = 1e-5,
      nystrom_ridge = 1e-10,
      nystrom_tol = 1e-10
    )
  )

  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_rbf, grid, "pmmr", seed)

  if (seed == 1L) {
    curves[["pmmr"]] <- data.frame(
      truth = lambda_rbf(grid),
      estimate = predict(fit, grid)
    )
  }
}

res <- do.call(rbind, results)
summary <- summarize_results(res)

print(res)
print(summary)

## Verify empirically that missingness really depends on B beyond V through S.
mechanism_dat <- generate_adjoint(300000, lambda_poly, 991)
mechanism_check <- data.frame(
  group = c("S=-1", "S=+1"),
  empirical_P_M1 = c(
    mean(mechanism_dat$M[mechanism_dat$S == -1]),
    mean(mechanism_dat$M[mechanism_dat$S == 1])
  ),
  target_P_M1 = c(p_minus, p_plus)
)
print(mechanism_check)

utils::write.csv(
  res,
  file.path(out_dir, "adjoint_validation_results.csv"),
  row.names = FALSE
)

utils::write.csv(
  summary,
  file.path(out_dir, "adjoint_validation_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  mechanism_check,
  file.path(out_dir, "observation_mechanism_check.csv"),
  row.names = FALSE
)

curve_df <- do.call(
  rbind,
  lapply(names(curves), function(nm) {
    data.frame(method = nm, curves[[nm]])
  })
)
utils::write.csv(
  curve_df,
  file.path(out_dir, "adjoint_seed1_curves.csv"),
  row.names = FALSE
)

pdf(file.path(out_dir, "adjoint_truth_vs_estimate.pdf"), width = 6, height = 6)
for (nm in names(curves)) {
  cc <- curves[[nm]]
  plot(
    cc$truth,
    cc$estimate,
    pch = 19,
    cex = 0.45,
    xlab = "True adjoint",
    ylab = "Estimated adjoint",
    main = nm
  )
  abline(0, 1, lty = 2)
}
dev.off()

## Tolerances are intentionally much wider than expected Monte Carlo error.
## They are implementation-failure diagnostics, not efficiency benchmarks.
stopifnot(all(res$rmse < 0.03))
stopifnot(all(res$max_abs < 0.07))

cat("Adjoint validation passed. Outputs written to:", out_dir, "\n")
