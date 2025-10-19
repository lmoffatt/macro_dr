# LAPACK Bridging With Row‑Major Matrices (No‑Copy Pattern)

This note documents how we call Fortran LAPACK from our row‑major `Matrix<T>` type without copying
large buffers. It explains the orientation conventions used across products, factorizations,
inverses, and eigensystems, and what identities tests should verify.

## Key Idea

- Fortran/LAPACK reads column‑major arrays; our `Matrix<T>` is row‑major.
- The same contiguous memory interpreted in the other order corresponds to the transpose.
- We either:
  - Use transpose‑in/transpose‑out for routines where it’s clearer (e.g., QR), or
  - Choose transposition flags and/or swap operands so LAPACK computes on transposed operands and we
    reinterpret results back in row‑major without extra copies.

## Implemented Patterns

- Full product `Z = X · Y` (dgemm)
  - We compute `Z^T = Y^T · X^T` with LAPACK, choosing `TRANSA/TRANSB` accordingly.
  - Code: `legacy/lapack_headers.h:156` (`Lapack_Full_Product`).

- QR factorization (dgeqrf/dorgqr)
  - We explicitly transpose into a temporary column‑major buffer, run the factorization, then
    transpose outputs back. This keeps semantics simple for `Q` and `R`.
  - Code: `legacy/lapack_headers.h:39–153` (`Lapack_QR`).

- Inverse (dgetrf/dgetri)
  - We pass our row‑major buffer directly. LAPACK computes `(A^T)^{-1}` in column‑major; when we read
    this buffer as row‑major, we obtain `((A^T)^{-1})^T = A^{-1}`.
  - Code: `legacy/lapack_headers.h:529–585` (`Lapack_Full_inv`).

- Nonsymmetric eigensystem (dgeevx)
  - LAPACK sees `A^T`. Right eigenvectors of `A^T` are left eigenvectors of `A`, and vice‑versa.
  - Mapping back to eigenvectors of `A`:
    - `VR_cpp = tr(VL_lapack)` → right eigenvectors of `A` (columns)
    - `VL_cpp = tr(VR_lapack)` → left  eigenvectors of `A` (columns)
  - We rescale columns so `u_i^T v_i = 1`, then sort eigenvalues (more positive first) and permute
    VR/VL consistently.
  - Code: `legacy/lapack_headers.h:1554–1561` (mapping), `legacy/lapack_headers.h:1143–1175` (sorting).

- Symmetric eigensystem (dsyevx)
  - Since `A^T = A`, the orientation trick is moot; we still sort eigenvalues and return `(L, VR, VL)`
    consistently with the nonsymmetric path.
  - Code: `legacy/lapack_headers.h:1717+`.

## Reconstruction Identities (What Tests Should Check)

Let the public return order be `(L, VR, VL)` where columns of `VR`/`VL` are right/left eigenvectors
of `A`, and `L` holds the eigenvalues on the diagonal (sorted, more positive first).

- Right‑only reconstruction (robust):
  - `A = VR · L · inv(VR)`
- Left‑only reconstruction (robust):
  - `A = inv(tr(VL)) · L · tr(VL)`
- Mixed (bi‑orthonormal case):
  - If `tr(VL) · VR = I`, then `A = VR · L · tr(VL)`
- Mixed (general case):
  - Define Gram `G = tr(VL) · VR`. Then `A = VR · L · G^{-1} · tr(VL)`

Notes:
- Degenerate clusters (repeated eigenvalues) make individual columns non‑unique; prefer the right‑only
  or left‑only reconstructions as primary checks. The Gram‑corrected mixed form is stable too.
- For symmetric `A`, `VR = VL` (orthonormal) and `A = VR · L · tr(VR)`.

## Sorting and Permutations

- We sort eigenvalues in non‑increasing order (more positive first) and apply the same column
  permutation to `VR` and `VL`.
- Code: `legacy/lapack_headers.h:1143–1175`.

## Practical Guidance

- When writing tests, prefer:
  - `A ≈ VR · L · inv(VR)` and/or `A ≈ inv(tr(VL)) · L · tr(VL)`
  - For mixed reconstruction, use the Gram correction if you don’t enforce full bi‑orthonormality.
- Keep in mind the mapping between LAPACK outputs (`VL_lapack`, `VR_lapack`) and C++ (`VR_cpp`,
  `VL_cpp`). We always publish the C++‑layout `(L, VR, VL)`.

## Spectrum Constraints (Real‑Only)

- MacroDR does not use complex arithmetic; eigensystem code assumes a strictly real spectrum and
  returns an error when LAPACK reports non‑zero imaginary parts for the eigenvalues.
- In the nonsymmetric path, after `dgeevx` returns `WR/WI`, we reject the decomposition if
  `norm_1(WI) > sqrt(eps)`.
- Tests and examples should therefore use matrices with real eigenvalues (e.g., symmetric matrices,
  CTMC generators `Q`, or upper‑triangular with a real diagonal).
