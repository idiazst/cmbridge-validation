# R implementation review and method selection

## Statistical target

The package standardizes estimators for the conditional moment equation

`E[Y - D h(X) | Z] = 0`.

For the bridge equation, use `Y=1`, `D=M`, `X=V`, and `Z=B`. For the adjoint,
restrict to complete cases (`M=1`) and use `Y=phi(V)`, `D=1`, `X=B`, and `Z=V`.
This avoids forcing every learner into a special-purpose Fredholm API.

## Search and selection

I reviewed the conditional-moment methods cited in the manuscript and searched
for maintained R implementations. The strongest R paths I found were:

1. **Sieve minimum distance / sieve NPIV.** Ai and Chen (2003) formulate
   conditional-moment models with unknown functions as sieve minimum-distance
   problems. The modern CRAN package `npiv` (Racine and Christensen, version
   0.1.3) implements B-spline sieve NPIV with data-driven sieve dimension based
   on Chen, Christensen, and Kankanala (2025), and has a documented API,
   binaries, CRAN checks, and a reverse importer. This is the strongest
   production R implementation I found for the sieve family.

2. **Landweber--Fridman iterative regularization.** The CRAN package `np`
   provides `npregiv(..., method="Landweber-Fridman")`, explicitly implementing
   the Horowitz (2011) approach, with bandwidth selection, an optimal stopping
   rule, reusable tuning objects, and evaluation-data support. The package is
   mature and extensively documented, although the `npregiv` help page itself
   still labels that function as beta-test software.

3. **Kernel maximum moment restriction / PMMR.** Muandet, Jitkrittum, and
   Kuebler (2020) show that RKHS maximum-moment restriction characterizes a
   conditional moment restriction. Mastouri et al. (2021) develop a kernel
   maximum-moment estimator for proximal bridge functions. I did not find a
   maintained CRAN implementation of PMMR, adversarial GMM, VMM, or DeepGMM.
   I did verify an R replication implementation of PMMR in the public PCSE
   replication repository for Park, Stensrud, and Tchetgen Tchetgen. That code
   has the particularly relevant interface `Y`, `Target`, `Perturb`, and
   `Diagonal`, and evaluates the held-out kernel moment loss. Among the modern
   nonlinear conditional-moment learners, PMMR therefore had the clearest R
   implementation path even though it is research replication code rather than
   a package.

### Methods not selected

- `gmm` is exceptionally mature and well documented, and is useful when a
  finite set of test functions and a finite-dimensional sieve have already been
  chosen. It is a strong implementation option, but I treated it as the engine
  underlying a finite-series approximation rather than as a distinct
  nonparametric learner for this three-method library.
- Tikhonov NPIV is also implemented in `np::npregiv` and `crs::crsiv`. I chose
  Landweber rather than Tikhonov because the R documentation explicitly notes
  the `n x n` memory burden of Tikhonov and recommends Landweber for larger
  samples.
- The official/reference implementations I found for adversarial GMM, VMM,
  DeepGMM, and neural maximum-moment approaches are Python-based. The CRAN
  package named `deepgmm` is a Gaussian-mixture-model package and is unrelated
  to DeepGMM for conditional moments.

## Reimplementation decision

I did **not** simply call `npiv::npiv` or `np::npregiv` inside `fit_bridge`.
Those APIs solve the standard NPIV residual `Y-h(X)`. The bridge residual is
`1-M h(V)`. One can reduce the bridge to standard NPIV on complete cases by
estimating `P(M=1|B)` and using its inverse as a pseudo-response, but that adds
an unnecessary nuisance estimator and partly defeats the direct conditional-
moment formulation of the manuscript.

Instead, the package implements the three regularization ideas directly for
`Y-D h(X)` while giving them a common interface and output class:

- `sieve_md`: weighted sieve minimum distance using polynomial or additive
  B-spline bases;
- `landweber`: Landweber iteration applied to the same empirical conditional-
  moment operator, with iteration count acting as the regularizer;
- `pmmr`: a scalable low-rank Gaussian-kernel maximum-moment learner with
  target and critic centers. This is a package-quality reimplementation of the
  PMMR criterion, not copied source code from the replication repository.

The direct implementation is intentional: it preserves the exact bridge
moment, permits `V` to be missing when `M=0`, and makes all three learners
interchangeable in stacking code.

## References

See `REFERENCES.bib` for complete BibTeX entries.

## PMMR implementation correction

The PMMR implementation was audited against Mastouri et al. (2021) and the authors' public code. PMMR minimizes a kernel maximum-moment V-statistic plus the RKHS norm of the target function. The first package draft instead used inverse-covariance weighting of finite critic features and an ordinary Euclidean penalty on raw target RBF coefficients; that is closer to a kernel-feature sieve minimum-distance estimator than to PMMR. The corrected implementation uses Nyström-normalized features for both RKHSs, the unweighted kernel maximum-moment V-statistic in the Nyström feature space, and the corresponding RKHS norm penalty. This is a scalable low-rank approximation of the PMMR criterion.
