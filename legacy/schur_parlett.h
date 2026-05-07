#ifndef SCHUR_PARLETT_H
#define SCHUR_PARLETT_H

// Schur+Parlett block evaluation of a matrix function f(Q).
//
// Algorithm (Higham, Functions of Matrices, ch. 9):
//   1) Schur decompose Q = U · T · U^T (T upper quasi-triangular, U orthogonal).
//   2) Detect the natural block structure of T: 1×1 blocks for real
//      eigenvalues, 2×2 blocks for complex-conjugate pairs.
//   3) Optionally re-cluster: merge adjacent blocks whose eigenvalues are
//      close so the off-diagonal Sylvester solves stay well-conditioned. The
//      naive Parlett recurrence has a (T_jj − T_ii) divisor that blows up at
//      clustered eigenvalues; grouping puts the troublesome eigenvalues
//      inside one block where the function is evaluated by a different,
//      stable method (Pade or Taylor on the block itself).
//   4) Evaluate f on each diagonal block via a user-supplied callable
//      (eval_block). Off-diagonal blocks F_ij are obtained by solving the
//      Sylvester equation
//          T_ii · F_ij − F_ij · T_jj
//                  = T_ij · F_jj − F_ii · T_ij + Σ_{i<k<j} (F_ik · T_kj − T_ik · F_kj)
//      reading the right-hand side bottom-up so each F_ij depends only on
//      already-computed entries.
//   5) Transform back: f(Q) = U · F · U^T.
//
// This header exposes the building blocks. The next layer (expm_schur,
// frechet_integral_schur, calc_Qdt_schur) will compose them with specific
// eval_block choices for each matrix function.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "exponential_matrix.h"
#include "lapack_schur.h"
#include "matrix.h"
#include "maybe_error.h"
// inside_out / outside_in for Derivative<Matrix<double>, Parameters_transformed>
// — used by the templated Van Loan augmented-matrix builders below to keep
// the same code path for plain double and Derivative inputs.
#include "parameters_derivative.h"

namespace lapack {

// Identify the Schur block structure of an upper quasi-triangular T. A 2×2
// block lives at (i, i+1) when the (i+1, i) sub-diagonal entry is nonzero.
// Returns the starting row/col index of each block, plus a sentinel = T.nrows().
//
// Threshold: anything below `tol · max|T|` is treated as zero. With dgees output
// the sub-diagonals of 1×1 blocks are exactly zero (set by the algorithm), so
// the threshold is generous and only matters if tiny numerical fuzz creeps in.
inline std::vector<std::size_t> detect_schur_block_starts(Matrix<double> const& T,
                                                           double tol = 1e-12) {
    const std::size_t N = T.nrows();
    std::vector<std::size_t> starts;
    if (N == 0) {
        starts.push_back(0);
        return starts;
    }

    double max_T = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < T.ncols(); ++j) max_T = std::max(max_T, std::abs(T(i, j)));
    const double zero_tol = tol * std::max(max_T, 1.0);

    std::size_t i = 0;
    while (i < N) {
        starts.push_back(i);
        if (i + 1 < N && std::abs(T(i + 1, i)) > zero_tol) {
            i += 2;  // 2×2 block (complex-conjugate eigenvalue pair)
        } else {
            i += 1;  // 1×1 block (real eigenvalue)
        }
    }
    starts.push_back(N);  // sentinel
    return starts;
}

// Extract a square block T[r0..r1, r0..r1] from T into a freshly allocated
// matrix.
inline Matrix<double> extract_block(Matrix<double> const& T, std::size_t r0, std::size_t r1) {
    const std::size_t n = r1 - r0;
    Matrix<double> B(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) B(i, j) = T(r0 + i, r0 + j);
    return B;
}

// Extract a rectangular off-diagonal block T[r0..r1, c0..c1].
inline Matrix<double> extract_offdiag(Matrix<double> const& T,
                                       std::size_t r0, std::size_t r1,
                                       std::size_t c0, std::size_t c1) {
    const std::size_t nr = r1 - r0;
    const std::size_t nc = c1 - c0;
    Matrix<double> B(nr, nc, 0.0);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) B(i, j) = T(r0 + i, c0 + j);
    return B;
}

// Write block X back into M[r0..r0+X.nrows, c0..c0+X.ncols].
inline void write_block(Matrix<double>& M, std::size_t r0, std::size_t c0,
                         Matrix<double> const& X) {
    for (std::size_t i = 0; i < X.nrows(); ++i)
        for (std::size_t j = 0; j < X.ncols(); ++j) M(r0 + i, c0 + j) = X(i, j);
}

// Solve A·X − X·B = C for X, with A (n_a×n_a) and B (n_b×n_b) small and
// well-separated in spectrum (so the solve is well-conditioned). Used for the
// Parlett off-diagonal updates. n_a, n_b ∈ {1, 2}: the largest system is
// 2×2 X giving 4 unknowns. Solve via direct vec(X) = (I⊗A − B^T⊗I)^{-1} vec(C).
inline Maybe_error<Matrix<double>> solve_sylvester_small(
    Matrix<double> const& A, Matrix<double> const& B, Matrix<double> const& C) {
    const std::size_t n_a = A.nrows();
    const std::size_t n_b = B.nrows();
    if (n_a == 0 || n_b == 0) return Matrix<double>(n_a, n_b, 0.0);

    // Assemble the n_a·n_b × n_a·n_b linear system M · vec(X) = vec(C),
    // where M = I_{n_b} ⊗ A − B^T ⊗ I_{n_a}, and vec is column-major.
    const std::size_t n = n_a * n_b;
    Matrix<double> M(n, n, 0.0);
    Matrix<double> rhs(n, 1, 0.0);
    auto idx = [n_a](std::size_t i, std::size_t j) { return i + j * n_a; };
    // Build M:
    for (std::size_t i = 0; i < n_a; ++i)
        for (std::size_t j = 0; j < n_b; ++j) {
            // Equation for X[i,j]: Σ_k A[i,k]·X[k,j] − Σ_l X[i,l]·B[l,j] = C[i,j]
            std::size_t row = idx(i, j);
            for (std::size_t k = 0; k < n_a; ++k) M(row, idx(k, j)) += A(i, k);
            for (std::size_t l = 0; l < n_b; ++l) M(row, idx(i, l)) -= B(l, j);
            rhs(row, 0) = C(i, j);
        }

    // Solve via partial-pivot Gaussian elimination on the small system.
    // (Up to 4×4, hand-rolled is fine and avoids LAPACK plumbing.)
    Matrix<double> aug(n, n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) aug(i, j) = M(i, j);
        aug(i, n) = rhs(i, 0);
    }
    for (std::size_t k = 0; k < n; ++k) {
        // Pivot
        std::size_t piv = k;
        double piv_v = std::abs(aug(k, k));
        for (std::size_t r = k + 1; r < n; ++r) {
            if (std::abs(aug(r, k)) > piv_v) {
                piv = r;
                piv_v = std::abs(aug(r, k));
            }
        }
        if (piv_v < 1e-30) {
            return error_message(std::string{"solve_sylvester_small: singular system"
                                              " (clustered eigenvalues across blocks?)"});
        }
        if (piv != k) {
            for (std::size_t j = 0; j <= n; ++j) std::swap(aug(k, j), aug(piv, j));
        }
        // Eliminate
        for (std::size_t r = k + 1; r < n; ++r) {
            double f = aug(r, k) / aug(k, k);
            for (std::size_t j = k; j <= n; ++j) aug(r, j) -= f * aug(k, j);
        }
    }
    // Back-substitute
    Matrix<double> x(n, 1, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double s = aug(i, n);
        for (std::size_t j = i + 1; j < n; ++j) s -= aug(i, j) * x(j, 0);
        x(i, 0) = s / aug(i, i);
    }
    // Repack to n_a × n_b
    Matrix<double> X(n_a, n_b, 0.0);
    for (std::size_t i = 0; i < n_a; ++i)
        for (std::size_t j = 0; j < n_b; ++j) X(i, j) = x(idx(i, j), 0);
    return X;
}

// Block-Parlett evaluation of f(T) for upper quasi-triangular T. eval_block
// is a callable Block -> Block returning f(T_kk). Caller pre-detects block
// structure (or uses detect_schur_block_starts).
//
// Returns F = f(T) as an upper-quasi-triangular matrix (sub-diagonal entries
// from the original 2×2 Schur blocks are preserved through the same block
// evaluation, so F has the same block-zero structure as T below the first
// sub-diagonal).
template <class EvalBlock>
Maybe_error<Matrix<double>> parlett_block_eval(
    Matrix<double> const& T,
    std::vector<std::size_t> const& starts,  // includes sentinel = T.nrows()
    EvalBlock&& eval_block) {
    const std::size_t N = T.nrows();
    Matrix<double> F(N, N, 0.0);
    if (N == 0) return F;
    if (starts.size() < 2 || starts.back() != N) {
        return error_message(std::string{"parlett_block_eval: starts[] must end with N"});
    }
    const std::size_t nb = starts.size() - 1;  // number of blocks

    // Diagonal blocks.
    std::vector<Matrix<double>> F_diag;
    F_diag.reserve(nb);
    std::vector<Matrix<double>> T_diag;
    T_diag.reserve(nb);
    for (std::size_t b = 0; b < nb; ++b) {
        auto T_bb = extract_block(T, starts[b], starts[b + 1]);
        T_diag.push_back(T_bb);
        auto F_bb = eval_block(T_bb);
        F_diag.push_back(F_bb);
        write_block(F, starts[b], starts[b], F_bb);
    }

    // Off-diagonal blocks, computed by super-diagonal index (j − i = 1, 2, ...).
    // F_ij depends on F_kk, T_ij, and (F_ik, T_kj, T_ik, F_kj) for i < k < j.
    for (std::size_t d = 1; d < nb; ++d) {
        for (std::size_t i = 0; i + d < nb; ++i) {
            const std::size_t j = i + d;
            const std::size_t r0 = starts[i], r1 = starts[i + 1];
            const std::size_t c0 = starts[j], c1 = starts[j + 1];

            // From F·T = T·F (commutativity of f(T) with T) we get the
            // Sylvester relation T_ii · F_ij − F_ij · T_jj
            //     = F_ii · T_ij − T_ij · F_jj + Σ_{i<k<j} (F_ik · T_kj − T_ik · F_kj).
            auto T_ij = extract_offdiag(T, r0, r1, c0, c1);
            Matrix<double> rhs = F_diag[i] * T_ij - T_ij * F_diag[j];
            for (std::size_t k = i + 1; k < j; ++k) {
                const std::size_t mr0 = starts[k], mr1 = starts[k + 1];
                auto F_ik = extract_offdiag(F, r0, r1, mr0, mr1);
                auto T_kj = extract_offdiag(T, mr0, mr1, c0, c1);
                auto T_ik = extract_offdiag(T, r0, r1, mr0, mr1);
                auto F_kj = extract_offdiag(F, mr0, mr1, c0, c1);
                rhs = rhs + F_ik * T_kj - T_ik * F_kj;
            }

            auto Maybe_X = solve_sylvester_small(T_diag[i], T_diag[j], rhs);
            if (!Maybe_X) return Maybe_X.error();
            write_block(F, r0, c0, Maybe_X.value());
        }
    }

    return F;
}

// Eigenvalues of a 1×1 or 2×2 diagonal block. For larger merged blocks the
// Bavely-Stewart criterion only needs the *spectrum*, so we walk the 1×1/2×2
// natural sub-blocks of the merged piece and pool their eigenvalues.
inline std::vector<std::complex<double>> block_eigenvalues_natural(Matrix<double> const& B) {
    std::vector<std::complex<double>> out;
    if (B.nrows() == 1) {
        out.emplace_back(B(std::size_t{0}, std::size_t{0}), 0.0);
        return out;
    }
    if (B.nrows() == 2) {
        double a = B(std::size_t{0}, std::size_t{0});
        double b = B(std::size_t{0}, std::size_t{1});
        double c = B(std::size_t{1}, std::size_t{0});
        double d = B(std::size_t{1}, std::size_t{1});
        double tr = a + d;
        double det = a * d - b * c;
        double disc = (tr * tr) / 4.0 - det;
        if (disc >= 0) {
            double r = std::sqrt(disc);
            out.emplace_back(tr / 2.0 + r, 0.0);
            out.emplace_back(tr / 2.0 - r, 0.0);
        } else {
            double r = std::sqrt(-disc);
            out.emplace_back(tr / 2.0, r);
            out.emplace_back(tr / 2.0, -r);
        }
        return out;
    }
    // For larger merged blocks: not invoked directly during clustering since
    // we cluster from the natural blocks. If called anyway, return empty so
    // callers can detect non-natural input.
    return out;
}

// Bavely-Stewart clustering: starting from the natural 1×1/2×2 Schur blocks,
// merge two adjacent blocks if any eigenvalue of one is within
// `delta · max(|λ₁|, |λ₂|, 1)` of an eigenvalue of the other. Iterate until
// no further merges occur. Returns merged block starts (always includes the
// sentinel = T.nrows() at the end).
inline std::vector<std::size_t> cluster_adjacent_blocks(
    Matrix<double> const& T,
    std::vector<std::size_t> const& natural_starts,
    double delta = 0.1) {
    if (natural_starts.size() < 3) return natural_starts;  // <2 blocks, nothing to merge

    std::vector<std::size_t> starts = natural_starts;
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t b = 0; b + 2 < starts.size(); ++b) {
            // Block b is starts[b]..starts[b+1], block b+1 is starts[b+1]..starts[b+2].
            auto B1 = extract_block(T, starts[b], starts[b + 1]);
            auto B2 = extract_block(T, starts[b + 1], starts[b + 2]);
            // Compute eigenvalues of each piece by walking the natural 1×1/2×2
            // sub-blocks (every merged block is a contiguous slice of the
            // natural decomposition).
            auto walk_natural = [&](std::size_t s_lo, std::size_t s_hi) {
                std::vector<std::complex<double>> evs;
                std::size_t k = s_lo;
                while (k < s_hi) {
                    std::size_t step = 1;
                    if (k + 1 < s_hi && std::abs(T(k + 1, k)) > 1e-12) step = 2;
                    auto sub = extract_block(T, k, k + step);
                    auto e = block_eigenvalues_natural(sub);
                    evs.insert(evs.end(), e.begin(), e.end());
                    k += step;
                }
                return evs;
            };
            auto evs1 = walk_natural(starts[b], starts[b + 1]);
            auto evs2 = walk_natural(starts[b + 1], starts[b + 2]);

            bool merge = false;
            for (auto const& l1 : evs1) {
                for (auto const& l2 : evs2) {
                    double scale = std::max({std::abs(l1), std::abs(l2), 1.0});
                    if (std::abs(l1 - l2) < delta * scale) {
                        merge = true;
                        break;
                    }
                }
                if (merge) break;
            }
            if (merge) {
                starts.erase(starts.begin() + (b + 1));
                changed = true;
                break;  // restart the scan since indices shifted
            }
        }
    }
    return starts;
}

// Raw matrix exponential by Pade + scaling-and-squaring, no
// to_transition_Probability projection (which is only valid for stochastic
// matrices). Templated on CMatrix so both plain Matrix<double> and
// Derivative<Matrix<double>, Parameters_transformed> work uniformly — Pade
// is built from matrix arithmetic + one inverse, all of which propagate
// derivatives through the codebase's overloaded operators. Used as the
// per-block evaluator inside expm_schur_parlett (plain only there, since
// Lapack_Schur is plain), and as the unconditionally-stable expm for the
// Frechet integrals' augmented matrices (works in both modes).
template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline Maybe_error<CMatrix> expm_pade_scaling_squaring(CMatrix const& M) {
    auto const& M_prim = var::primitive(M);
    double max_abs = 0.0;
    for (std::size_t i = 0; i < M_prim.nrows(); ++i)
        for (std::size_t j = 0; j < M_prim.ncols(); ++j) {
            max_abs = std::max(max_abs, std::abs(M_prim(i, j)));
        }
    int n = 0;
    if (max_abs > 0.5) {
        n = static_cast<int>(std::ceil(std::log2(max_abs / 0.5)));
    }
    if (n < 0) n = 0;
    double scale = std::pow(2.0, -n);
    auto scaled = M * scale;
    auto Maybe_E = expm_pade(scaled);
    if (!Maybe_E) return Maybe_E.error();
    auto E = std::move(Maybe_E.value());
    for (int i = 0; i < n; ++i) E = E * E;
    return E;
}

// Public API: matrix exponential of Q via Schur+Parlett. For 1×1 diagonal
// blocks we use scalar exp; for 2×2 / merged blocks we fall through to
// expm_pade_scaling_squaring. Clustering threshold δ=0.1 is the
// Bavely-Stewart default — close eigenvalues are absorbed into super-blocks
// so the inter-block Sylvester solves stay well-conditioned.
inline Maybe_error<Matrix<double>> expm_schur_parlett(Matrix<double> const& Q,
                                                       double cluster_delta = 0.1) {
    auto eval_block = [](Matrix<double> const& B) -> Matrix<double> {
        if (B.nrows() == 1) {
            return Matrix<double>(1, 1, std::exp(B(std::size_t{0}, std::size_t{0})));
        }
        return expm_pade_scaling_squaring(B).value();
    };
    auto Maybe_schur = Lapack_Schur(Q);
    if (!Maybe_schur) return Maybe_schur.error();
    auto& [T, U] = Maybe_schur.value();

    auto natural = detect_schur_block_starts(T);
    auto starts = cluster_adjacent_blocks(T, natural, cluster_delta);
    auto Maybe_F = parlett_block_eval(T, starts, eval_block);
    if (!Maybe_F) return Maybe_F.error();

    return U * Maybe_F.value() * tr(U);
}

// Build the 2N×2N Van Loan augmented matrix [[Q·dt, G·dt], [0, Q·dt]],
// templated on the matrix type. Plain Matrix<double> goes through directly;
// Derivative<Matrix<double>, Parameters_transformed> uses inside_out (to
// expose Matrix<Derivative<double>>), element-fills the augmented layout,
// then outside_in to repack to Derivative<Matrix>. The roundtrip preserves
// parameter sensitivities cleanly without poking at the d_d_ internals.
template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline auto build_van_loan_2x2_block(CMatrix const& Q, CMatrix const& G, double dt) {
    auto Q_dt = Q * dt;
    auto G_dt = G * dt;
    auto Q_inner = var::inside_out(Q_dt);  // Matrix<scalar_t>
    auto G_inner = var::inside_out(G_dt);
    using scalar_t = std::decay_t<decltype(Q_inner(std::size_t{0}, std::size_t{0}))>;
    const std::size_t N = Q_inner.nrows();
    scalar_t zero = Q_inner(std::size_t{0}, std::size_t{0}) -
                     Q_inner(std::size_t{0}, std::size_t{0});
    Matrix<scalar_t> M_aug_inner(2 * N, 2 * N, zero);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            M_aug_inner(i, j) = Q_inner(i, j);
            M_aug_inner(i, j + N) = G_inner(i, j);
            M_aug_inner(i + N, j + N) = Q_inner(i, j);
        }
    }
    if constexpr (var::is_derivative_v<CMatrix>) {
        return var::outside_in(M_aug_inner, Q.dx());
    } else {
        return M_aug_inner;
    }
}

// Build the 3N×3N Van Loan augmented matrix
//   [[Q·dt, G·dt,    0   ],
//    [   0, Q·dt,  G·dt ],
//    [   0,    0,  Q·dt ]]
// with the same templated construction strategy.
template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline auto build_van_loan_3x3_block(CMatrix const& Q, CMatrix const& G, double dt) {
    auto Q_dt = Q * dt;
    auto G_dt = G * dt;
    auto Q_inner = var::inside_out(Q_dt);
    auto G_inner = var::inside_out(G_dt);
    using scalar_t = std::decay_t<decltype(Q_inner(std::size_t{0}, std::size_t{0}))>;
    const std::size_t N = Q_inner.nrows();
    scalar_t zero = Q_inner(std::size_t{0}, std::size_t{0}) -
                     Q_inner(std::size_t{0}, std::size_t{0});
    Matrix<scalar_t> M_aug_inner(3 * N, 3 * N, zero);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            M_aug_inner(i, j) = Q_inner(i, j);
            M_aug_inner(i + N, j + N) = Q_inner(i, j);
            M_aug_inner(i + 2 * N, j + 2 * N) = Q_inner(i, j);
            M_aug_inner(i, j + N) = G_inner(i, j);
            M_aug_inner(i + N, j + 2 * N) = G_inner(i, j);
        }
    }
    if constexpr (var::is_derivative_v<CMatrix>) {
        return var::outside_in(M_aug_inner, Q.dx());
    } else {
        return M_aug_inner;
    }
}

// Frechet derivative of expm in direction G, equivalent to the time integral
//
//   A(dt) = ∫₀^dt expm(Q·s) · G · expm(Q·(dt − s)) ds.
//
// Computed via Van Loan's block trick: with M_aug = [[Q·dt, G·dt], [0, Q·dt]]
// (size 2N×2N), expm(M_aug) = [[expm(Q·dt), A(dt)], [0, expm(Q·dt)]] — the
// upper-right block IS A(dt). Templated on the matrix type so plain double
// and Derivative<Matrix<double>> both go through the same Pade
// scaling-and-squaring (Schur+Parlett would require eigenvalue reordering
// for the duplicated-eigenvalue augmented structure).
template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline Maybe_error<CMatrix> frechet_integral_schur(CMatrix const& Q, CMatrix const& G, double dt,
                                                    double cluster_delta = 0.1) {
    auto const& Q_prim = var::primitive(Q);
    auto const& G_prim = var::primitive(G);
    const std::size_t N = Q_prim.nrows();
    if (N == 0 || N != Q_prim.ncols() || N != G_prim.nrows() || N != G_prim.ncols()) {
        return error_message(std::string{"frechet_integral_schur: Q and G must be N×N "
                                          "with the same N"});
    }

    auto M_aug = build_van_loan_2x2_block(Q, G, dt);
    auto Maybe_E = expm_pade_scaling_squaring(M_aug);
    (void)cluster_delta;
    if (!Maybe_E) return Maybe_E.error();
    auto& E = Maybe_E.value();

    // Extract upper-right N×N block. Use inside_out → element-extract →
    // outside_in to keep the templated path uniform across plain / Derivative.
    auto E_inner = var::inside_out(E);
    using scalar_t = std::decay_t<decltype(E_inner(std::size_t{0}, std::size_t{0}))>;
    Matrix<scalar_t> A_inner(N, N,
                              E_inner(std::size_t{0}, std::size_t{0}) -
                                  E_inner(std::size_t{0}, std::size_t{0}));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) A_inner(i, j) = E_inner(i, j + N);

    if constexpr (var::is_derivative_v<CMatrix>) {
        return var::outside_in(A_inner, Q.dx());
    } else {
        return A_inner;
    }
}

// Second Frechet integral / iterated time integral
//
//   B(dt) = ∫₀^dt ∫₀^{t1} expm(Q·s)·G·expm(Q·(t1−s))·G·expm(Q·(dt−t1)) ds dt1.
//
// Van Loan's block trick again: the 3N×3N augmented matrix
//
//   M_aug = [[Q·dt, G·dt,    0   ],
//            [   0, Q·dt,  G·dt ],
//            [   0,    0,  Q·dt ]]
//
// has expm with upper-right N×N block equal to B(dt). Pade
// scaling-and-squaring on the augmented matrix is unconditionally stable
// (eigenvalues of M_aug are Q·dt's eigenvalues each tripled — Schur+Parlett
// would require reordering to handle the duplicates, which we skip).
// Consolidated Frechet pipeline: one Pade scaling-and-squaring on the 3N×3N
// Van Loan augmented matrix
//   M_aug = [[Q·dt, G·dt,    0   ],
//            [   0, Q·dt,  G·dt ],
//            [   0,    0,  Q·dt ]]
// yields P, A and B simultaneously as the (0,0), (0,1), (0,2) blocks of
// expm(M_aug). Replaces three separate expm calls (N×N for P, 2N×2N for A,
// 3N×3N for B) with one — saves the N×N and 2N×2N work entirely.
template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline Maybe_error<std::tuple<CMatrix, CMatrix, CMatrix>>
frechet_p_a_b_combined(CMatrix const& Q, CMatrix const& G, double dt) {
    auto const& Q_prim = var::primitive(Q);
    auto const& G_prim = var::primitive(G);
    const std::size_t N = Q_prim.nrows();
    if (N == 0 || N != Q_prim.ncols() || N != G_prim.nrows() || N != G_prim.ncols()) {
        return error_message(std::string{"frechet_p_a_b_combined: Q and G must be N×N "
                                          "with the same N"});
    }

    auto M_aug = build_van_loan_3x3_block(Q, G, dt);
    auto Maybe_E = expm_pade_scaling_squaring(M_aug);
    if (!Maybe_E) return Maybe_E.error();
    auto& E = Maybe_E.value();

    auto E_inner = var::inside_out(E);
    using scalar_t = std::decay_t<decltype(E_inner(std::size_t{0}, std::size_t{0}))>;
    scalar_t zero = E_inner(std::size_t{0}, std::size_t{0}) -
                     E_inner(std::size_t{0}, std::size_t{0});
    Matrix<scalar_t> P_inner(N, N, zero);
    Matrix<scalar_t> A_inner(N, N, zero);
    Matrix<scalar_t> B_inner(N, N, zero);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            P_inner(i, j) = E_inner(i, j);                   // top-left  = P
            A_inner(i, j) = E_inner(i, j + N);               // (0,1)     = A
            B_inner(i, j) = E_inner(i, j + 2 * N);           // (0,2)     = B
        }
    }

    if constexpr (var::is_derivative_v<CMatrix>) {
        return std::tuple{var::outside_in(P_inner, Q.dx()),
                           var::outside_in(A_inner, Q.dx()),
                           var::outside_in(B_inner, Q.dx())};
    } else {
        return std::tuple{P_inner, A_inner, B_inner};
    }
}

template <class CMatrix>
    requires(var::U<CMatrix, Matrix<double>>)
inline Maybe_error<CMatrix> frechet_double_integral_schur(CMatrix const& Q, CMatrix const& G,
                                                            double dt) {
    auto const& Q_prim = var::primitive(Q);
    auto const& G_prim = var::primitive(G);
    const std::size_t N = Q_prim.nrows();
    if (N == 0 || N != Q_prim.ncols() || N != G_prim.nrows() || N != G_prim.ncols()) {
        return error_message(std::string{"frechet_double_integral_schur: Q and G must "
                                          "be N×N with the same N"});
    }

    auto M_aug = build_van_loan_3x3_block(Q, G, dt);
    auto Maybe_E = expm_pade_scaling_squaring(M_aug);
    if (!Maybe_E) return Maybe_E.error();
    auto& E = Maybe_E.value();

    // Extract top-row, right-column N×N block (rows 0..N, cols 2N..3N).
    auto E_inner = var::inside_out(E);
    using scalar_t = std::decay_t<decltype(E_inner(std::size_t{0}, std::size_t{0}))>;
    Matrix<scalar_t> B_inner(N, N,
                              E_inner(std::size_t{0}, std::size_t{0}) -
                                  E_inner(std::size_t{0}, std::size_t{0}));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) B_inner(i, j) = E_inner(i, j + 2 * N);

    if constexpr (var::is_derivative_v<CMatrix>) {
        return var::outside_in(B_inner, Q.dx());
    } else {
        return B_inner;
    }
}

// Full Schur + Parlett pipeline with optional clustering. Q → f(Q) given a
// per-block evaluator. The evaluator must handle blocks of any size that
// clustering can produce — pass an evaluator that uses Pade or
// scaling-squaring on larger merged blocks.
template <class EvalBlock>
Maybe_error<Matrix<double>> schur_parlett(Matrix<double> const& Q, EvalBlock&& eval_block,
                                            double cluster_delta = 0.1) {
    auto Maybe_schur = Lapack_Schur(Q);
    if (!Maybe_schur) return Maybe_schur.error();
    auto& [T, U] = Maybe_schur.value();

    auto natural = detect_schur_block_starts(T);
    auto starts = cluster_adjacent_blocks(T, natural, cluster_delta);
    auto Maybe_F = parlett_block_eval(T, starts, std::forward<EvalBlock>(eval_block));
    if (!Maybe_F) return Maybe_F.error();

    return U * Maybe_F.value() * tr(U);
}

}  // namespace lapack

#endif  // SCHUR_PARLETT_H
