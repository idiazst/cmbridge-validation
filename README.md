# cmbridge

`cmbridge` standardizes three learners for conditional moment equations of the form

`E[Y - D h(X) | Z] = 0`.

The initial learner library contains:

1. `sieve_md`: sieve minimum distance with polynomial or additive B-spline bases;
2. `landweber`: Landweber iterative regularization of the same empirical conditional-moment operator;
3. `pmmr`: scalable Gaussian-kernel proxy/maximum-moment restriction using low-rank target and critic features.

The package exposes a common generic interface and two manuscript-oriented wrappers:

```r
# Generic conditional moment
fit <- fit_cm(response = Y, diagonal = D,
              target = X, instrument = Z,
              method = "sieve_md")

# Inverse-probability bridge: E[M beta(V)-1 | B] = 0
fit_beta <- fit_bridge(B = B, V = V, M = M, method = "pmmr")

# Adjoint: E[M{lambda(B)-phi(V)} | V] = 0
fit_lambda <- fit_adjoint(B = B, V = V, M = M, phi = phi,
                          method = "landweber")

predict(fit_beta, newdata)
moment_loss(fit_beta)
```

`V` is allowed to be missing when `M=0` in `fit_bridge()`. The adjoint wrapper uses complete cases, which is equivalent because its conditional moment is multiplied by `M`.

See `REVIEW.md` for the implementation review and method-selection rationale, `REFERENCES.bib` for citations, and `inst/simulations/validate_large_sample.R` for the large-sample recovery study.
