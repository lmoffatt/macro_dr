# ADR-003 — Eigenvector Derivatives, Gauge Ambiguity, and Testing Strategy

Status: Accepted

## Context

MacroDR uses eigendecompositions of nonsymmetric real matrices (e.g., CTMC generators `Q`) in several
pipelines. While eigenvalues are uniquely determined (for simple spectra), eigenvectors are not:

- For simple eigenvalues, each right eigenvector `v_j` and the corresponding left eigenvector `u_j`
  are defined only up to a nonzero scalar factor. Any rescaling `v_j ← α v_j`, `u_j ← u_j / α`
  preserves the identity `A = V Λ U^T` and the biorthogonality constraint `u_j^T v_j = 1`.
- For repeated eigenvalues, any basis of the eigenspace is valid. Only the invariant subspace is
  intrinsic; individual columns may rotate under arbitrarily small perturbations.

Because of this gauge freedom, the derivative of eigenvectors is not intrinsic: it depends on the
normalization (gauge) chosen as a function of the parameters. Different gauges lead to different
(valid) derivatives of `V` and `U`.

In practice, LAPACK (DGEEVX) returns right/left eigenvectors normalized to unit 2‑norm with an
additional phase choice (the component of largest magnitude is made real). Our analytic derivative
formulas, instead, enforce biorthogonality `U^T V = I` and typically fix the gauge by setting the
diagonal of the coupling matrix `Ω` to zero (`diag(Ω) = 0`). These are different gauges. Finite‑
difference estimates approximate the derivative of LAPACK’s particular gauge, while the analytic
formula computes the derivative under the biorthogonal/`diag(Ω)=0` gauge. Therefore the two can
disagree even with well‑separated eigenvalues.

Additional practical sources of mismatch:

- Sign/phase flips: LAPACK’s phase rule can change discontinuously as parameters vary (e.g., when the
  largest component changes index or crosses zero). This introduces jumps between `x+h` and `x−h` that
  finite differences will capture but analytic derivatives (which assume a smooth gauge) will not.
- Sorting/tie behavior: we reorder eigenpairs by eigenvalue (more positive first). Within any
  numerically equal (or near‑equal) block, column order is not stable unless additionally aligned to a
  reference. Tiny perturbations may swap columns across `±h`.
- Stabilization steps: any post‑processing that clamps, masks, or reprojects eigendata is typically
  non‑smooth; it breaks derivative tests unless disabled.

By contrast, gauge‑invariant functions of the eigensystem—such as the exponential `Qdt = V e^{Λ dt} U^T`
and spectral projectors `P_j = v_j u_j^T`—have unique derivatives (away from eigenvalue crossings). We
observe that derivative tests for `Qdt` pass when all stabilizations are disabled.

## Decision

1) Do not assert derivative tests on raw eigenvectors `V`/`U` in production. Prefer testing
   gauge‑invariant quantities built from the eigensystem (e.g., `Qdt`, spectral projectors, full
   reconstruction `A ≈ V Λ U^T`).

2) For derivative testing of eigendecomposition‑based flows, disable non‑smooth stabilizers
   (`StabilizerPolicyDisabled`) so that the mapping is as smooth as possible.

3) If and when direct testing of `V`/`U` is necessary, align the numerical eigenpairs to a reference
   gauge before differencing:
   - permute columns by maximum correlation with the reference columns; and
   - flip column signs to ensure positive correlation with the reference.

4) Keep analytic formulas for eigenpair derivatives in biorthogonal gauge with `diag(Ω)=0`:
   - `G = U^T (dA) V`
   - for `i ≠ j`, `Ω_{ij} = G_{ij} / (λ_j − λ_i)`; set `Ω_{ii} = 0`
   - `dV = V Ω`, `dU = − U Ω^T`

   These are standard and numerically sound away from repeated eigenvalues.

## Rationale

Comparing derivatives of gauge‑dependent objects to finite differences of a potentially different,
discontinuous gauge (LAPACK’s) is ill‑posed and fragile. Basing tests on gauge‑invariant quantities
eliminates this ambiguity and yields meaningful, stable checks of the codepaths that matter (e.g., the
transition kernel over a small time step).

## Consequences

- Unit/integration tests should target `Qdt` (or downstream observables) rather than raw `V`/`U`.
- Where a numeric stability policy applies non‑smooth transforms, tests must explicitly disable it for
  derivative checks.
- The derivative of `V`/`U` remains available for internal use but should be interpreted as gauge‑
  dependent; discrepancies with finite differences are not necessarily bugs.

## Notes on Implementation

- The left‑eigenvector derivative uses `dU = − U Ω^T`. Ensure `Ω` is constructed off‑diagonal from
  `G = U^T (dA) V` and apply the transpose in the left update.
- If at some point we need to compare `V`/`U` across parameters, optional alignment helpers (permute
  by correlation; sign‑flip by inner product with reference) can be used to bring `±h` to a common
  gauge prior to differencing.

## References

- docs/math/eigensystem.md — computation, conventions, and reconstruction identities.
- docs/testing.md — derivative test guidance and practical procedures.

