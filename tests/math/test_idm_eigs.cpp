#include "catch_amalgamated.hpp"

#include <cmath>
#include <limits>

#include "lapack_headers.h"
#include "matrix.h"

using ::Matrix;
using ::SymPosDefMatrix;
using ::SymmetricMatrix;

namespace {

double rel_residual_1(const Matrix<double>& Aref, const Matrix<double>& A) {
    using ::norm_1;
    auto N = norm_1(Aref);
    if (N == 0.0)
        return norm_1(A);
    return norm_1(A - Aref) / N;
}

Matrix<double> identity(std::size_t n) {
    Matrix<double> I(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) I(i, i) = 1.0;
    return I;
}

SymPosDefMatrix<double> make_diag_spd(std::initializer_list<double> diag_vals) {
    const std::size_t n = diag_vals.size();
    SymPosDefMatrix<double> out(n, n, false);
    std::size_t i = 0;
    for (auto v : diag_vals) {
        out.set(i, i, v);
        ++i;
    }
    return out;
}

double max_abs_diff(const SymPosDefMatrix<double>& A, const SymPosDefMatrix<double>& B) {
    double m = 0.0;
    for (std::size_t i = 0; i < A.nrows(); ++i)
        for (std::size_t j = 0; j < A.ncols(); ++j)
            m = std::max(m, std::abs(A(i, j) - B(i, j)));
    return m;
}

}  // namespace

TEST_CASE("eigs_idm identity sanity", "[idm][eigs]") {
    auto H = make_diag_spd({1.0, 1.0, 1.0});
    auto J = make_diag_spd({1.0, 1.0, 1.0});
    auto H0 = H;
    auto J0 = J;

    auto me = eigs_idm(H, J);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    const double tol = 1e-12;
    for (std::size_t i = 0; i < L.size(); ++i) REQUIRE(std::abs(L[i] - 1.0) < tol);

    auto I = identity(3);
    REQUIRE(rel_residual_1(I, tr(VR) * VR) < 1e-10);
    REQUIRE(rel_residual_1(I, VR * L * tr(VL)) < 1e-10);

    REQUIRE(max_abs_diff(H, H0) < std::numeric_limits<double>::epsilon());
    REQUIRE(max_abs_diff(J, J0) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("eigs_idm diagonal closed form", "[idm][eigs]") {
    auto H = make_diag_spd({1.0, 2.0, 4.0});
    auto J = make_diag_spd({4.0, 6.0, 8.0});

    auto me = eigs_idm(H, J);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    REQUIRE(std::abs(L[0] - 4.0) < 1e-12);
    REQUIRE(std::abs(L[1] - 3.0) < 1e-12);
    REQUIRE(std::abs(L[2] - 2.0) < 1e-12);

    Matrix<double> Cref(3, 3, 0.0);
    Cref(0, 0) = 4.0;
    Cref(1, 1) = 3.0;
    Cref(2, 2) = 2.0;

    REQUIRE(rel_residual_1(identity(3), tr(VR) * VR) < 1e-10);
    REQUIRE(rel_residual_1(Cref, VR * L * tr(VL)) < 1e-10);
}

TEST_CASE("eigs_idm matches Cholesky-whitened IDM", "[idm][eigs]") {
    SymPosDefMatrix<double> H(3, 3, false);
    H.set(0, 0, 4.0);
    H.set(0, 1, 1.0);
    H.set(0, 2, 0.5);
    H.set(1, 1, 3.0);
    H.set(1, 2, 0.2);
    H.set(2, 2, 2.0);

    SymPosDefMatrix<double> J(3, 3, false);
    J.set(0, 0, 5.0);
    J.set(0, 1, 0.3);
    J.set(0, 2, 0.1);
    J.set(1, 1, 2.5);
    J.set(1, 2, 0.4);
    J.set(2, 2, 1.8);

    auto me = eigs_idm(H, J);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    for (std::size_t i = 0; i < L.size(); ++i) REQUIRE(L[i] > 0.0);
    REQUIRE(rel_residual_1(identity(3), tr(VR) * VR) < 1e-10);

    auto H_dense = to_dense(H);
    int N = static_cast<int>(H.nrows());
    int LDA = N;
    int INFO = 0;
    char UPLO = kSymmetricUplo;
    lapack::dpotrf_(&UPLO, &N, &H_dense[0], &LDA, &INFO);
    REQUIRE(INFO == 0);

    DownTrianMatrix<double> L_factor(H.nrows(), true);
    if (UPLO == 'U') {
        for (std::size_t i = 0; i < H.nrows(); ++i)
            for (std::size_t j = 0; j <= i; ++j) L_factor.set(i, j, H_dense(i, j));
    } else if (UPLO == 'L') {
        for (std::size_t i = 0; i < H.nrows(); ++i)
            for (std::size_t j = 0; j <= i; ++j) L_factor.set(i, j, H_dense(j, i));
    } else {
        REQUIRE(false);
    }

    auto Maybe_L_inv = inv(L_factor);
    REQUIRE(Maybe_L_inv.valid());

    auto J_dense = to_dense(J);
    auto Cref = Maybe_L_inv.value() * J_dense * tr(Maybe_L_inv.value());
    auto Crec = VR * L * tr(VL);
    REQUIRE(rel_residual_1(Cref, Crec) < 1e-10);
}

TEST_CASE("eigs_idm fails when H is not SPD", "[idm][eigs]") {
    SymmetricMatrix<double> bad_H(2);
    bad_H.set(0, 0, 1.0);
    bad_H.set(0, 1, 0.0);
    bad_H.set(1, 1, -1.0);
    auto H = SymPosDefMatrix<double>::I_sware_it_is_possitive(std::move(bad_H));
    auto J = make_diag_spd({1.0, 2.0});

    auto me = eigs_idm(H, J);
    REQUIRE_FALSE(me.valid());
}
