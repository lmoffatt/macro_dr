# Implementation Note: Posterior Information Distortion

## Purpose

This note describes how the Posterior Information Distortion framework
(defined in `supplement_posterior_main.tex` and companion supplements in
this folder) is realized in code, what infrastructure already exists, and
what remains to be wired up.

## Conceptual contrast with the Likelihood folder

The Likelihood Information Distortion framework, documented in
`theory/macroir/docs/Likelihood_Information_Distortion/`, requires
retained-eigenspace projection (`implementation_subspace_information_distortion.md`)
because `H_lik` can be singular or near-singular. The whole subspace
machinery exists to give the diagnostics a meaning when likelihood-only
inversion is ill-posed.

The Posterior Information Distortion framework removes this problem
structurally. With a proper Gaussian prior, the posterior Hessian
`H_post = H_lik + H_prior` is positive definite under a checkable
dominance condition (Section 2.1 of the main supplement), so
`H_post^{-1/2}` is always well-defined. No subspace projection is
required when the posterior framework is in use; the existing low-level
spectral utilities (`compute_psd_decomp`, `apply_normalized_congruence`,
etc.) still apply, but they will succeed on `H_post` for any well-chosen
proper prior.

The choice between the Likelihood path and the Posterior path is then a
choice about what diagnostic the user wants to read, not about whether
the math is well-posed. See the planned `supplement_safety_categorization.tex`
for the bootstrap-level policy that decides when to report the
Likelihood IDM (strict-or-NaN per replicate) and when to fall back to
the Posterior IDM (always defined).

## Prior infrastructure (already in place)

The prior type and CSV loader were built for the evidence-computation
pipeline (parallel tempering, thermodynamic integration) and are
production-ready:

- `var::Parameters_Normal_Distribution` in
  [`legacy/parameters_distribution.h`](../../../legacy/parameters_distribution.h):
  multivariate Gaussian with diagonal covariance in the transformed
  parameter space (log10 by convention). Methods:
  - `logPrior(prior, params)` — log p(theta)
  - `prior.dlogP(params)` — score grad log pi(theta)
  - `prior.hessian(params)` — Hessian of log pi (= -Sigma_prior^{-1})
  - `prior.FIM(params)` — precision matrix Sigma_prior^{-1}
    (= the H_prior of the Posterior framework)
- `var::load_Prior(filename, sep, model_name, names)` for CSV-based
  prior loading (legacy).
- `var::create_prior(ModelName, ParamNames, prior_info)` and the cmd
  wrapper `macrodr::cmd::create_prior(model, prior_info)`, added in the
  `cmd: add create_prior DSL function` commit, allow constructing a
  prior in-line from a `.macroir` script. Each tuple is
  `(name, transformation, value, transformed_variance)`; the prior
  mean in transformed space is derived as `tr(value)`.

## Numerical Fisher Information (already in place)

`calculate_mnumerical_fisher_information` in
[`src/core/likelihood.cpp`](../../../src/core/likelihood.cpp) computes
`H_lik` via finite differences of the score, returning a
`parameter_spd_payload` wrapping a `SymPosDefMatrix<double>`. The
returned matrix is symmetrized but not validated for positive
semidefiniteness — finite-difference bias can produce indefinite
matrices in weakly identified directions, which is the case the
Posterior framework handles via the dominance condition.

## What is not yet wired up

These pieces are documented as future commits, not present in code as of
this writing:

- `calc_numerical_fisher_information_posterior`: combines
  `calculate_mnumerical_fisher_information` (giving `H_lik`) with the
  prior contribution `prior.FIM()` (giving `H_prior`) into
  `H_post = H_lik + H_prior`. The bootstrap pipeline currently
  averages `F_per_recording[i]` over a resampled index set to build
  `F_b`; the posterior counterpart adds `H_prior` once to that average.
- `Likelihood_Information_Distortion` / `Posterior_Information_Distortion`
  type pair in `legacy/distributions.h`: the existing
  `Likelihood_Information_Distortion` (renamed from the original
  `Information_Distortion_Matrix` in the `code: rename diagnostic classes
  with Likelihood_ prefix` commit) needs a Posterior counterpart whose
  construction uses `H_post` and `JT_post = JT_lik + H_prior`.
- The same Posterior counterparts for sibling diagnostics:
  `Posterior_Sample_Distortion`, `Posterior_Correlation_Distortion`,
  `Posterior_Distortion_Corrected_Covariance`,
  `Posterior_Distortion_Induced_Bias`,
  `Posterior_Information_Distortion_Reconstituted`,
  `Posterior_Fisher_Covariance`,
  `Posterior_Gaussian_Fisher_Distortion`,
  `Posterior_Numerical_Fisher_Information`.
- Binary writer registration (mirror of the Likelihood add_dense / add_scalar
  calls in `write_bootstrap_samples_binary` at
  [`src/core/likelihood.cpp`](../../../src/core/likelihood.cpp)) so that
  Posterior diagnostics appear as separate CSV columns alongside the
  Likelihood ones.
- The bootstrap safety categorization (SAFE / MARGINAL / UNRELIABLE /
  BIASED) discussed in
  `program/notes/diagnostics_golden_test_plan.md` and to be formalized
  in `supplement_safety_categorization.tex`.

## Numerical-stability notes

- The whitening transform in Section 4.5 of the main supplement uses the
  spectral decomposition `H_post = U_post Lambda_post U_post^T` and the
  symmetric inverse square root `U_post Lambda_post^{-1/2} U_post^T`,
  which is what `apply_normalized_congruence` already computes via
  eigendecomposition (it does NOT use Cholesky despite the name suggesting
  otherwise; see the comment chain in `legacy/lapack_headers.h` around
  `compute_psd_decomp`). Per-parameter diagonals `[C_post]_{ii}` are
  therefore basis-invariant under this implementation, consistent with
  the supplement's Section 4.5 claim.
- The dominance condition `H_prior - H_lik^- > 0` in Section 2.1 of the
  main supplement is the necessary-and-sufficient guarantee that `H_post`
  is positive definite when `H_lik` is finite-difference indefinite.
  A simple sufficient bound
  `lambda_min(H_prior) > max(0, -lambda_min(H_lik))`
  can be evaluated from two scalar eigenvalue computations and used as
  a fast monitoring statistic; a rejection by the sufficient bound
  triggers the matrix-level check, which costs one eigendecomposition
  of `H_lik` plus an eigenvalue test of `H_prior - H_lik^-`.

## Scope of this note

This note covers Sections 1-4 of the main supplement (definitions of
posterior log-posterior, score, Hessian decomposition, score fluctuation
matrices, and the central PID diagnostic). Sections 5+ of the main
supplement, the evidence-correction supplement, the information-gain
supplement, and the safety-categorization supplement will be written in
subsequent commits, and this implementation note will be extended to
match.
