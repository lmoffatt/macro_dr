# Implementation Note: Subspace Information Distortion

## Purpose

This note describes how singular `H`, `J_sample`, and `J_total` are handled in the codebase for information distortion diagnostics.

The scientific document uses an identifiable-subspace formulation. The implementation mirrors that directly: diagonalize the reference PSD matrix, retain only eigenmodes above a cutoff, and evaluate inverse powers only on that retained eigenspace.

## Why the Previous Path Failed

The original helpers:

- `idm_matrix(H, J)`
- `dcc_matrix(H, J)`
- `c_h_r_c_h_matrix(C, R)`

used Cholesky-based LAPACK routines. That requires the reference matrix to be strictly positive definite.

This fails when:

- `H` is singular or nearly singular
- `J_sample` is singular
- `C_sample` becomes singular after subspace projection

In the diagnostics path in [`src/core/likelihood.cpp`](/home/lmoffatt/Code/macro_dr/macro_dr/src/core/likelihood.cpp), those failures previously produced empty/default matrices through `.value_or(SymPosDefMatrix<double>{})`.

## Current Subspace Helpers

The singular-capable helpers are:

- `idm_matrix_subspace(H, J_total, rtol, atol)`
- `dcc_matrix_subspace(H, J_total, rtol, atol)`
- `sample_distortion_matrix_subspace(H, J_sample, rtol, atol)`
- `correlation_distortion_matrix_subspace(J_sample, J_total, rtol, atol)`
- `c_h_r_c_h_matrix_subspace(C_sample, R, rtol, atol)`
- `distortion_induced_bias_subspace(H, g, rtol, atol)`

The matrix helpers return `Maybe_error<SymPosDefMatrix<double>>`.

The DIB helper returns `Maybe_error<Matrix<double>>`.

These are semantic wrappers. Internally they use two numerical kernels:

- retained-eigenspace normalized congruence:
  `A^{-1/2} B A^{-1/2}` on the retained eigenspace of `A`
- retained-eigenspace congruence:
  `A^{1/2} B A^{1/2}` on the retained eigenspace of `A`
- retained-eigenspace inverse application:
  `U_A \Lambda_A^{-1} U_A^T g` on the retained eigenspace of `A`

## Retained-Eigenspace Algorithm

Given a PSD reference matrix `A`:

1. Compute a symmetric eigendecomposition of `A`.
2. Let `lambda_max` be the largest eigenvalue.
3. Retain eigenvalues satisfying
   `lambda_i > max(atol, rtol * lambda_max)`.
4. Build the retained basis `U` from the corresponding eigenvectors.
5. Form reduced-space matrices with `U^T B U`.
6. Apply inverse square roots or square roots only to the retained eigenvalues.
7. Reconstruct the ambient-coordinate matrix as `U * reduced_result * U^T`.

The implementation keeps the full eigendecomposition before truncation, so the retained-space work is derived from a decomposition object that still reconstructs the original PSD operator and its retained approximation.

Default tolerances:

- `rtol = 1e-10`
- `atol = 0.0`

## Failure Behavior

The subspace helpers still return `Maybe_error` for legitimate runtime failures:

- incompatible dimensions
- non-finite entries
- PSD contract violations beyond tolerance
- LAPACK eigensolver failures
- empty retained subspace

The diagnostics code keeps the existing caller behavior:

- compute the matrix helper
- call `.value_or(SymPosDefMatrix<double>{})`
- continue producing the rest of the diagnostics

This preserves current robustness while allowing singular but meaningful cases to succeed.

## Mapping to Likelihood Diagnostics

In `calculate_Likelihood_diagnostics_evolution_f`:

- total information distortion uses `idm_matrix_subspace(H, J_total)`
- sample distortion uses `sample_distortion_matrix_subspace(H, J_sample)`
- correlation distortion uses `correlation_distortion_matrix_subspace(J_sample, J_total)`
- reconstituted distortion uses `c_h_r_c_h_matrix_subspace(C_sample, R)`
- corrected covariance uses `dcc_matrix_subspace(H, J_total)`
- distortion-induced bias uses `distortion_induced_bias_subspace(H, g)`

This matches the supplement notation:

- `H`: Gaussian Fisher curvature
- `J_total`: covariance of the total score
- `J_sample`: sum of within-interval score covariances
- `C_sample`: sample distortion on the retained `H`-subspace
- `R`: correlation distortion on the retained `J_sample`-subspace
- `Sigma_DCC`: corrected covariance on the retained `H`-subspace
- `g`: expected score vector
- `b_DIB`: distortion-induced bias on the retained `H`-subspace

## Interpretation

Ambient-coordinate outputs are returned as full-size symmetric matrices, but they should be interpreted as living on the retained informative subspace.

Zero directions outside the retained subspace do not indicate finite certainty. They indicate directions that were not retained as numerically informative for the corresponding distortion object.
