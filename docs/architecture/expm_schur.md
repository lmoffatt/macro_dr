# Schur-Based Matrix Exponential for Generators Q

Status: draft / exploration

This note captures a practical, production-ready path to compute the matrix exponential `exp(t Q)` using the real Schur decomposition. It is intended for conservative CTMC generators `Q` (row sums zero, spectrum with Re(λ) ≤ 0 and at least one zero eigenvalue) but remains applicable to general real matrices.

Goals
- Numerically robust on non-normal matrices (typical for Q).
- Avoid explicit eigenvectors; cope well with clustered/complex eigenvalues.
- Provide both full `exp(tQ)` and action-only `exp(tQ)·x` variants.
- Leave room to enforce CTMC invariants (nonnegativity, row-sum 1) post hoc.

Overview
1) Compute the real Schur decomposition
   - `Q = U T U^T`, with `U` orthogonal and `T` quasi-upper-triangular (1×1 and 2×2 blocks).
   - LAPACK: `dgees` (JOBVS='V', SORT='N') for real Schur + Schur vectors.
2) Compute `F = exp(t T)` using a method specialized for (block) triangular matrices.
   - Preferred: Scaling-and-squaring with Padé approximants on the Schur form (Higham’s algorithm).
   - Alternative: Block Parlett recurrence on `T` (solves Sylvester equations for off-diagonal blocks).
3) Recover `exp(tQ) = U F U^T`.

Why Schur?
- Backward-stable orthogonal similarity transformations minimize roundoff amplification.
- Handles complex conjugate eigenpairs as 2×2 real blocks in `T` without complex arithmetic.
- Avoids forming/inverting eigenvector matrices, which can be ill-conditioned for non-normal `Q`.

Algorithm Details

Option A — Scaling-and-Squaring + Padé on Schur form (recommended)
1) Choose Padé degree `m` ∈ {3, 5, 7, 9, 13} and scaling `s` via 1-norm bounds (Higham’s θ_m).
   - Compute `A = t T` and `α = ||A||_1`. If `α > θ_m`, set `s = ceil(log2(α/θ_m))`; scale `A ← A / 2^s`.
2) Build Padé [m/m] approximant
   - `R_m(A) = (I + ∑ c_odd A^odd) · (I − ∑ c_even A^even)^{-1}`.
   - Since `A` is (block) upper triangular, powers `A^k` remain triangular; denominator solve reduces to triangular back/forward substitutions (efficient and stable).
3) Square `s` times: `F ← R_m(A)` then repeat `F ← F·F` for `k = 1..s`.
4) Return `exp(tQ) = U F U^T`.

Pros: High robustness across spectra; widely accepted as state-of-the-art. Cons: Slightly more code (Padé coefficients, θ_m thresholds, triangular solves).

Option B — Block Parlett Recurrence
1) Partition `T` into diagonal 1×1 or 2×2 blocks; compute block exponentials on the diagonal:
   - 1×1: `F_ii = exp(A_ii)`.
   - 2×2: Stable real formula using trace/shift: let `τ = tr(A_ii)/2`, `B = A_ii − τI`, then
     `exp(A_ii) = e^τ [cosh(δ) I + (sinh(δ)/δ) B]`, with `δ^2 = det(B)` (use robust real handling).
2) For off-diagonal blocks `i < j`, solve Sylvester-type equations
   - `A_ii F_ij − F_ij A_jj = RHS`,
   - `RHS = Σ_{k=i+1}^{j−1} (F_ik A_kj − A_ik F_kj)`
   - Proceed by increasing block distance `(j−i)` (forward substitution on blocks).
3) Compose `exp(tQ) = U F U^T`.

Pros: Elegant, no scaling/squaring; good on moderate sizes. Cons: More delicate near repeated/clustered eigenvalues; requires careful 2×2 handling and Sylvester solvers.

CTMC-Specific Notes (Q is a generator)
- Invariants: `Q·1 = 0` ⇒ there is a zero eigenvalue; `exp(tQ)` should be row-stochastic (row sums = 1, nonnegative entries).
- Schur-based `exp` preserves these only approximately due to rounding; you may optionally
  - Clamp tiny negative entries to 0, and renormalize rows to sum to 1 (tolerance-based).
  - Provide this step as an opt-in post-processor for workflows requiring exact stochasticity.
- Zero mode: handled naturally by Schur; no special code beyond standard numerics is required.

Uniformization (for contrast)
- For sparse/large CTMCs, the action `exp(tQ)·v` is often computed via uniformization:
  - Choose `μ ≥ max_i(−Q_ii)`, define `P = I + Q/μ` (stochastic), then
  - `exp(tQ) = e^{−μ t} Σ_{k=0}^∞ (μ t)^k/k! · P^k`.
  - Truncate using Poisson tail bounds. Guarantees nonnegativity/row sums; great for action-only.
- This is complementary to Schur; we can support both.

Proposed APIs

Low-level (LAPACK wrappers)
- `Maybe_error<std::pair<Matrix<double>, Matrix<double>>> Lapack_Schur(const Matrix<double>& Q);`
  - Returns `(U, T)` such that `Q = U T U^T`.

Triangular exponential on Schur form
- `Matrix<double> expm_schur_triangular(const Matrix<double>& A)`
  - Computes `exp(A)` where `A` is quasi-upper-triangular.
  - Implementation: Padé 13 + scaling/squaring with triangular solves (preferred) or block Parlett.

User-facing
- `Maybe_error<Matrix<double>> expm_schur(const Matrix<double>& Q, double t = 1.0)`
  - Steps: `(U,T) = Schur(Q)`, `F = expm_schur_triangular(t·T)`, return `U F U^T`.
- `Maybe_error<Matrix<double>> expm_multiply_schur(const Matrix<double>& Q, const Matrix<double>& x, double t = 1.0)`
  - Action-only: `y = U^T x`, `y ← expm_schur_triangular(t·T)·y` (triangular algorithm), return `U y`.

CTMC helpers (optional, opt-in)
- `Matrix<double> project_row_stochastic(const Matrix<double>& P, double eps = 1e-14)`
  - Clamp tiny negatives to 0, renormalize rows to sum to 1 (if deviation ≤ tolerance).

Numerical Details to Get Right
- Real Schur via `dgees` with Schur vectors (JOBVS='V').
- Padé coefficients and `θ_m` thresholds per Higham; choose degree adaptively from `||A||_1`.
- Use triangular solves for `(I − Σ c_even A^even)^{-1}·(I + Σ c_odd A^odd)` to avoid general solves.
- Stable 2×2 block exponential formula; avoid cancellation.
- Squaring uses triangular multiplication (still O(n^3) but cache friendly).

Complexity
- Full dense `exp(Q)`: O(n^3), dominated by Schur and triangular multiplications.
- Action-only: O(n^2) per multiply once Schur is available; cheaper if reusing `(U,T)`.

Testing Strategy
- Algebraic identities: small to moderate matrices
  - Compare `expm_schur(Q)·1` against `1` (row sums) within tolerance.
  - Check `expm_schur(tQ) ≈ expm_schur(0.5 t Q)^2` (scaling/squaring consistency).
  - Compare against uniformization on small Q.
- Edge cases
  - Near-repeated eigenvalues.
  - Mixed real/complex 2×2 blocks.
  - Very small/large `t` (stress scaling choice).

Concrete constants (Higham 2005)
- θ_m 1‑norm bounds (double precision):
  - m=3: 1.495585217958292e-002
  - m=5: 2.539398330063230e-001
  - m=7: 9.504178996162932e-001
  - m=9: 2.097847961257068e+000
  - m=13: 5.371920351148152e+000
- Padé [m/m] coefficients (descending powers), e.g. for m=13:
  - b0=64764752532480000, b2=32382376266240000, b4=7771770303897600,
    b6=1187353796428800, b8=129060195264000, b10=10559470521600,
    b12=670442572800, b14=33522128640, b16=1323241920, b18=40840800,
    b20=960960, b22=16380, b24=182, b26=1
  - Standard construction uses even/odd separation into numerator/denominator polynomials.

LAPACK Schur call details
- Use `dgees` (real Schur):
  - JOBVS = 'V' (compute Schur vectors U)
  - SORT = 'N' (no reordering here)
  - Select = nullptr (unused when SORT='N')
  - Workspace query with LWORK = -1, then allocate and call.
  - Outputs: T in-place (quasi-upper triangular), WR/WI (eigenvalues), and VS (U).

2×2 block exponential (real, robust)
- For 2×2 A: let τ = tr(A)/2, B = A − τ I, δ2 = det(B) = (a11−τ)(a22−τ) − a12 a21.
- Define series‑safe s = sinh(δ)/δ and c = cosh(δ), using stable evaluations for small |δ|:
  - For |δ| < 1e-6, use Taylor: s ≈ 1 + δ²/6 + …, c ≈ 1 + δ²/2 + …
- Then exp(A) = e^τ ( c I + s B ).

Triangular Padé pseudocode (Schur form)
```
function expm_schur_triangular(A):  # A quasi-upper triangular
  m, theta = choose_degree_and_theta(norm1(A))
  s = max(0, ceil(log2(norm1(A)/theta)))
  A = A / 2^s
  # Build powers
  A2 = A*A; A4 = A2*A2; A6 = A4*A2; ... as needed by m
  # Form P (odd) and Q (even) using Padé coefficients b_k
  P = b1*A + b3*A^3 + ...
  Q = b0*I + b2*A^2 + b4*A^4 + ...
  # R = (Q − P)^{-1} (Q + P); solve triangular system, do not invert
  Y = Q + P
  X = Q − P
  R = solve_triangular(X, Y)  # two triangular solves (left/right) as appropriate
  F = R
  repeat s times: F = F*F  # triangular mult retains structure
  return F
```

Acceptance and sanity checks
- Row-sum preservation (CTMC): max_i |(exp(tQ)·1)_i − 1| ≤ 1e−12 (tunable).
- Squaring consistency: exp(tQ) ≈ exp(0.5 t Q)^2 within tolerance tied to ‖Q‖ and machine eps.
- Cross-check against uniformization for small dense Q.

Integration Plan (incremental)
1) Add `Lapack_Schur(Q)` wrapper (dgees).
2) Implement `expm_schur_triangular(A)` using Padé 13 + scaling/squaring with triangular solves.
3) Expose `expm_schur(Q,t)` and (optionally) `expm_multiply_schur(Q,x,t)`.
4) Optional CTMC post-processor `project_row_stochastic`.
5) Tests: unit tests for triangular exp and full exp; CTMC sanity tests vs uniformization.

Risks / Trade-offs
- Dense algorithm is O(n^3); on large/sparse Q, use uniformization or Krylov methods for action-only.
- Schur-based full exp does not enforce row-stochasticity exactly; projection step is optional.
- Careful 2×2 handling and Padé selection are critical for robustness.

Open Questions
- Do we need a sparse/action-only path (uniformization/Krylov) in the near term?
- Should we expose a toggle to project to row-stochastic matrices post exp?
- Reuse of Schur factors across multiple `t` (cache `(U,T)` per Q) — API ergonomics?

References
- N. J. Higham, “The Scaling and Squaring Method for the Matrix Exponential Revisited,” SIAM J. Matrix Anal. Appl., 2005.
- N. J. Higham, Functions of Matrices: Theory and Computation, SIAM, 2008.
