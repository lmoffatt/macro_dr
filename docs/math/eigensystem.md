# Eigensystem Decomposition (Left/Right Vectors, Sorted Eigenvalues)

This note specifies how eigenpairs are computed and returned across the LAPACK wrappers, the high‑level `eigs` API, and the model‑layer `calc_eigen` helper.

## Summary

- We compute both right (`v_i`) and left (`u_i`) eigenvectors and scale them so that `u_i^T v_i = 1` (bi‑orthonormality for the nonsymmetric case; orthonormality for the symmetric case).
- Eigenvalues are ordered with the more positive first; the corresponding columns of left/right eigenvectors are permuted consistently.
- Public return convention for eigen‑systems is the triple `(L, VR, VL)`:
  - `L` — diagonal matrix with eigenvalues `λ_i` (non‑increasing: more positive first)
  - `VR` — columns are right eigenvectors `v_i`
  - `VL` — columns are left eigenvectors `u_i`

With this convention, the matrix reconstructs as `A = VR · L · VL^T`.

## LAPACK vs C++ Orientation and Naming

LAPACK arrays are column‑major; our `Matrix<double>` uses row‑major storage. To avoid copying the input matrix, we pass the buffer “as is” to LAPACK. LAPACK therefore operates on `A^T` (the transpose), which swaps the roles of left/right eigenvectors. We then transpose LAPACK’s outputs to recover C++ column semantics for the original `A`:

- `VL_lapack`, `VR_lapack` — left/right eigenvector arrays reported by LAPACK (column‑wise, for `A^T`).
- `VR_cpp = tr(VL_lapack)` — right eigenvectors of the original `A` (C++ layout).
- `VL_cpp = tr(VR_lapack)` — left  eigenvectors of the original `A` (C++ layout).

After this mapping we rescale `VL_cpp` and `VR_cpp` to enforce `u_i^T v_i = 1`.

## LAPACK Eigensystem

- For general (nonsymmetric) matrices, we call the LAPACK nonsymmetric eigensolver and return `(L, VR, VL)` with the properties above:
  - Both left and right eigenvectors are computed.
  - Left/right eigenvectors are scaled so `u_i^T v_i = 1` for each eigenpair.
  - Eigenvalues are sorted with more positive first; `VR`/`VL` columns are permuted accordingly.
- For symmetric matrices, we call the symmetric routine and return the same triple `(L, VR, VL)` where eigenvectors are orthonormal and `VL = VR` up to transposition conventions.

## `eigs` Function

- Primitive input `Matrix<double>`:
  - Returns `(L, VR, VL)` with left/right vectors and eigenvalues as above.
  - Example usage: `auto [L, VR, VL] = eigs(A).value();` then `A ≈ VR · L · VL^T`.
- Derivative input `Derivative<Matrix<double>, Parameters_transformed>`:
  - Returns `(dL, dVR, dVL)` to match the primitive order `(L, VR, VL)`.
  - First‑order identities used internally include `dλ = diag(tr(VL) · dA · VR)` and standard off‑diagonal resolvent corrections for the eigenvectors.

## Real‑Only Spectrum

- Scope: the implementation supports real matrices whose spectra are strictly real. Complex
  eigenvalues (nonzero imaginary parts) are not supported in MacroDR at this time (no complex
  arithmetic throughout the stack).
- Behavior: the LAPACK nonsymmetric path (`dgeevx`) returns `WR` and `WI` (real/imaginary parts). If
  the imaginary parts are non‑negligible, the eigensystem builder returns an error instead of
  attempting to construct complex eigenvectors.
  - Practical check used in code: if `norm_1(WI) > sqrt(eps)` → error
- Test guidance: choose matrices with strictly real spectra (e.g., symmetric matrices, CTMC
  generators `Q`, upper‑triangular with real diagonal). Do not feed matrices with complex eigenpairs
  into `eigs`.

## `calc_eigen` (Model Layer)

`calc_eigen` adapts the eigensystem into the form used by downstream model code:

- Starting from `(L, VR, VL)` returned by `eigs`:
  - Set `V = VR` (right eigenvectors as columns).
  - Set `W = tr(VL)` (left eigenvectors transposed).
- Because of the enforced bi‑orthonormality (`u_i^T v_i = 1`), we have `W = V^{-1}`. Thus, downstream code can use `A = V · L · W` interchangeably with `A = V · L · V^{-1}`.

## Practical Notes

- Shapes: `VL` and `VR` are `n×n`, `L` is a diagonal `n×n`.
- Degeneracies: within a degenerate eigenvalue cluster, the basis is not unique; only the invariant subspace is well‑defined. The sort order is deterministic (by eigenvalue), but individual columns may vary under perturbations.
- For CTMC generators `Q`, a dedicated post‑processing may reorder the zero‑mode and apply a consistent gauge; the core contracts above remain unchanged.
