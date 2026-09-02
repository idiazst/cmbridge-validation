test_that("adjoint wrapper uses the complete-case conditional moment", {
  set.seed(91)
  n <- 30000
  V <- runif(n, -0.8, 0.8)
  eps <- runif(n, -0.05, 0.05)
  B <- V + eps
  lambda_true <- function(b) 1.3 + 0.4*b + 0.1*b^2
  phi <- 1.3 + 0.4*V + 0.1*(V^2 + 0.05^2/3)
  p <- plogis(0.3 + 0.4*V)
  M <- rbinom(n, 1, p)
  f <- fit_adjoint(B, V, M, phi, "sieve_md",
                   list(target_basis="poly", instrument_basis="poly",
                        target_degree=2L, instrument_degree=5L, lambda=1e-8))
  g <- seq(-0.8, 0.8, length.out=51)
  expect_lt(sqrt(mean((predict(f,g)-lambda_true(g))^2)), 0.05)
})
