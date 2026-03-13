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

TEST_CASE("idm_matrix identity sanity", "[idm][matrix]") {
    auto H = make_diag_spd({1.0, 1.0, 1.0});
    auto J = make_diag_spd({1.0, 1.0, 1.0});
    auto H0 = H;
    auto J0 = J;

    auto me = idm_matrix(H, J);
    REQUIRE(me.valid());
    auto C = me.value();

    REQUIRE(rel_residual_1(identity(3), to_dense(C)) < 1e-12);
    REQUIRE(max_abs_diff(H, H0) < std::numeric_limits<double>::epsilon());
    REQUIRE(max_abs_diff(J, J0) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("idm_matrix diagonal closed form", "[idm][matrix]") {
    auto H = make_diag_spd({1.0, 2.0, 4.0});
    auto J = make_diag_spd({4.0, 6.0, 8.0});

    auto me = idm_matrix(H, J);
    REQUIRE(me.valid());
    auto C = me.value();

    Matrix<double> Cref(3, 3, 0.0);
    Cref(0, 0) = 4.0;
    Cref(1, 1) = 3.0;
    Cref(2, 2) = 2.0;

    REQUIRE(rel_residual_1(Cref, to_dense(C)) < 1e-12);
}

TEST_CASE("idm_matrix matches eigs_idm reconstruction", "[idm][matrix]") {
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

    auto Cmaybe = idm_matrix(H, J);
    REQUIRE(Cmaybe.valid());
    auto C = Cmaybe.value();

    auto E = eigs_idm(H, J);
    REQUIRE(E.valid());
    auto [L, VR, VL] = E.value();

    auto Crec = VR * L * tr(VL);
    REQUIRE(rel_residual_1(Crec, to_dense(C)) < 1e-10);
}

TEST_CASE("idm_matrix fails when H is not SPD", "[idm][matrix]") {
    SymmetricMatrix<double> bad_H(2);
    bad_H.set(0, 0, 1.0);
    bad_H.set(0, 1, 0.0);
    bad_H.set(1, 1, -1.0);
    auto H = SymPosDefMatrix<double>::I_sware_it_is_possitive(std::move(bad_H));
    auto J = make_diag_spd({1.0, 2.0});

    auto me = idm_matrix(H, J);
    REQUIRE_FALSE(me.valid());
}

TEST_CASE("idm_matrix fails when J is not SPD", "[idm][matrix]") {
    auto H = make_diag_spd({1.0, 2.0});
    SymmetricMatrix<double> bad_J(2);
    bad_J.set(0, 0, 1.0);
    bad_J.set(0, 1, 0.0);
    bad_J.set(1, 1, -1.0);
    auto J = SymPosDefMatrix<double>::I_sware_it_is_possitive(std::move(bad_J));

    auto me = idm_matrix(H, J);
    REQUIRE_FALSE(me.valid());
}
