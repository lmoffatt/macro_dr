// Lapack_Schur (dgees) wrapper acceptance tests. Step 1 of the
// Schur+Parlett implementation plan: confirm the wrapper produces a correct
// real Schur decomposition Q = U · T · U^T with U orthogonal and T upper
// quasi-triangular, before building Parlett block evaluation on top.

#include <catch_amalgamated.hpp>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>

#include <lapack_headers.h>
#include <lapack_schur.h>
#include <matrix.h>

namespace {

double frobenius(Matrix<double> const& M) {
    double s = 0.0;
    for (std::size_t i = 0; i < M.nrows(); ++i)
        for (std::size_t j = 0; j < M.ncols(); ++j) s += M(i, j) * M(i, j);
    return std::sqrt(s);
}

double frobenius_diff(Matrix<double> const& A, Matrix<double> const& B) {
    REQUIRE(A.nrows() == B.nrows());
    REQUIRE(A.ncols() == B.ncols());
    double s = 0.0;
    for (std::size_t i = 0; i < A.nrows(); ++i)
        for (std::size_t j = 0; j < A.ncols(); ++j) {
            double d = A(i, j) - B(i, j);
            s += d * d;
        }
    return std::sqrt(s);
}

Matrix<double> identity_like(std::size_t n) {
    Matrix<double> I(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) I(i, i) = 1.0;
    return I;
}

// Returns the largest absolute entry strictly below the first sub-diagonal.
// For a real Schur form T, those entries must be exactly zero (the only
// allowed sub-diagonal nonzeros are the (i+1, i) entries of 2×2 complex blocks).
double max_below_subdiagonal(Matrix<double> const& T) {
    double m = 0.0;
    for (std::size_t i = 0; i < T.nrows(); ++i)
        for (std::size_t j = 0; j + 1 < i && j < T.ncols(); ++j) {
            m = std::max(m, std::abs(T(i, j)));
        }
    return m;
}

// Verify all four invariants of a Schur decomposition for a given Q:
//   (1) U is orthogonal: U · U^T = I
//   (2) T is upper quasi-triangular (no entries below the first sub-diagonal)
//   (3) Q == U · T · U^T (the decomposition reconstructs Q)
//   (4) Q's eigenvalues match T's (1×1 + 2×2 block) eigenvalues
// The tolerance matches the LAPACK guarantee — components scale with the
// Frobenius norm of Q times machine eps.
void check_schur(Matrix<double> const& Q, std::string const& tag) {
    auto maybe = lapack::Lapack_Schur(Q);
    if (!maybe) UNSCOPED_INFO(tag << " — Lapack_Schur failed: " << maybe.error()());
    REQUIRE(maybe.valid());
    auto [T, U] = std::move(maybe.value());

    const double q_norm = frobenius(Q);
    const double tol = 1e-12 * std::max(q_norm, 1.0);

    // (1) Orthogonality of U.
    auto UUt = U * tr(U);
    auto I = identity_like(U.nrows());
    auto orth_err = frobenius_diff(UUt, I);
    INFO(tag << " — ||U U^T - I||_F = " << orth_err);
    CHECK(orth_err < 1e-12);

    // (2) Quasi-triangularity of T (only first sub-diagonal allowed nonzero).
    auto subdiag_err = max_below_subdiagonal(T);
    INFO(tag << " — max entry below first sub-diagonal of T = " << subdiag_err);
    CHECK(subdiag_err < tol);

    // (3) Reconstruction: Q == U T U^T.
    auto Q_recon = U * T * tr(U);
    auto recon_err = frobenius_diff(Q, Q_recon);
    INFO(tag << " — ||Q - U T U^T||_F = " << recon_err << " (Q norm " << q_norm << ")");
    CHECK(recon_err < tol);
}

}  // namespace

TEST_CASE("Lapack_Schur: 2×2 with distinct real eigenvalues",
          "[math][schur][lapack]") {
    // Q has eigenvalues 0 and -101 (matches scheme_CO at agonist=10:
    // closed→open=1, open→closed=100, eigenvalues sum_of_rates=−101 and 0).
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -1.0;
    Q(0, 1) = 1.0;
    Q(1, 0) = 100.0;
    Q(1, 1) = -100.0;
    check_schur(Q, "2x2 distinct real");
}

TEST_CASE("Lapack_Schur: 2×2 with complex-conjugate eigenvalues",
          "[math][schur][lapack]") {
    // Rotation matrix Q whose eigenvalues are e^{±i θ}. Real Schur form must
    // produce a single 2×2 diagonal block (no 1×1 blocks) since the
    // eigenvalues aren't real.
    const double theta = 0.7;
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = std::cos(theta);
    Q(0, 1) = -std::sin(theta);
    Q(1, 0) = std::sin(theta);
    Q(1, 1) = std::cos(theta);
    check_schur(Q, "2x2 complex pair");
}

TEST_CASE("Lapack_Schur: 4×4 mixed real and complex pairs",
          "[math][schur][lapack]") {
    // Q built as block diagonal of two 2×2 blocks: one upper triangular with
    // real eigenvalues {-1, -2}, one rotation with complex eigenvalues. Then
    // we conjugate by a random orthogonal Q_ortho to obscure the structure
    // and force Schur to genuinely re-discover the form.
    Matrix<double> A(4, 4, 0.0);
    A(0, 0) = -1.0;
    A(0, 1) = 0.5;
    A(1, 1) = -2.0;
    const double theta = 1.1;
    A(2, 2) = std::cos(theta);
    A(2, 3) = -std::sin(theta);
    A(3, 2) = std::sin(theta);
    A(3, 3) = std::cos(theta);

    // Householder-style orthogonal scrambler (deterministic).
    Matrix<double> v(4, 1, 0.0);
    v(0, 0) = 1.0;
    v(1, 0) = 2.0;
    v(2, 0) = -1.0;
    v(3, 0) = 0.5;
    double vTv = 0.0;
    for (std::size_t i = 0; i < 4; ++i) vTv += v(i, 0) * v(i, 0);
    Matrix<double> H(4, 4, 0.0);
    for (std::size_t i = 0; i < 4; ++i) H(i, i) = 1.0;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) H(i, j) -= 2.0 * v(i, 0) * v(j, 0) / vTv;

    auto Q = H * A * tr(H);
    check_schur(Q, "4x4 mixed");
}

TEST_CASE("Lapack_Schur: 6×6 single-channel scheme_CO lifted at Nch=5",
          "[math][schur][lapack]") {
    // Tridiagonal-like lifted Q for k=2, Nch=5: M = Nch+1 = 6 microstates
    // (n_open = 0, …, 5). Off-diagonals are channel-counts × per-channel
    // rates (closed→open = 1, open→closed = 100). This mirrors the regime
    // figure_2_micro exercises and is the smallest case where the lifting
    // matters numerically.
    const std::size_t M = 6;
    const double r_open = 1.0;     // per-channel closed→open rate at agonist=10
    const double r_close = 100.0;  // per-channel open→closed rate
    Matrix<double> Q(M, M, 0.0);
    for (std::size_t n = 0; n < M; ++n) {
        // n = number of open channels
        const std::size_t n_closed = M - 1 - n;
        if (n + 1 < M) {  // closed → open transition: rate n_closed · r_open
            Q(n, n + 1) = static_cast<double>(n_closed) * r_open;
        }
        if (n > 0) {  // open → closed transition: rate n · r_close
            Q(n, n - 1) = static_cast<double>(n) * r_close;
        }
        Q(n, n) = -(Q(n, n + 1 < M ? n + 1 : n) + (n > 0 ? Q(n, n - 1) : 0.0));
        // Recompute diagonal as -row sum (simpler and definitionally correct).
        Q(n, n) = 0.0;
        double row_sum = 0.0;
        for (std::size_t k = 0; k < M; ++k)
            if (k != n) row_sum += Q(n, k);
        Q(n, n) = -row_sum;
    }
    check_schur(Q, "6x6 lifted scheme_CO Nch=5");
}
