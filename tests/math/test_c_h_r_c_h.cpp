#include "catch_amalgamated.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "lapack_headers.h"
#include "matrix.h"

using ::Matrix;
using ::SymPosDefMatrix;

namespace {

double rel_residual_1(const Matrix<double>& Aref, const Matrix<double>& A) {
    using ::norm_1;
    auto N = norm_1(Aref);
    if (N == 0.0)
        return norm_1(A);
    return norm_1(A - Aref) / N;
}

double symmetry_defect_1(const Matrix<double>& A) {
    using ::norm_1;
    const auto N = std::max(1.0, norm_1(A));
    return norm_1(A - tr(A)) / N;
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

TEST_CASE("c_h_r_c_h identity sanity", "[c_h_r_c_h][matrix]") {
    auto C = make_diag_spd({1.0, 1.0, 1.0});
    auto R = make_diag_spd({1.0, 1.0, 1.0});
    auto C0 = C;
    auto R0 = R;

    auto out = lapack::Lapack_C_h_R_C_h(C, R);

    REQUIRE(rel_residual_1(identity(3), to_dense(out)) < 1e-12);
    REQUIRE(max_abs_diff(C, C0) < std::numeric_limits<double>::epsilon());
    REQUIRE(max_abs_diff(R, R0) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("c_h_r_c_h diagonal closed form", "[c_h_r_c_h][matrix]") {
    auto C = make_diag_spd({1.0, 4.0, 9.0});
    auto R = make_diag_spd({4.0, 6.0, 8.0});

    auto out = lapack::Lapack_C_h_R_C_h(C, R);

    Matrix<double> Cref(3, 3, 0.0);
    Cref(0, 0) = 4.0;
    Cref(1, 1) = 24.0;
    Cref(2, 2) = 72.0;

    REQUIRE(rel_residual_1(Cref, to_dense(out)) < 1e-12);
}

TEST_CASE("c_h_r_c_h matches Cholesky reference", "[c_h_r_c_h][matrix]") {
    SymPosDefMatrix<double> C(3, 3, false);
    C.set(0, 0, 4.0);
    C.set(0, 1, 1.0);
    C.set(0, 2, 0.5);
    C.set(1, 1, 3.0);
    C.set(1, 2, 0.3);
    C.set(2, 2, 2.5);

    SymPosDefMatrix<double> R(3, 3, false);
    R.set(0, 0, 5.0);
    R.set(0, 1, 0.2);
    R.set(0, 2, 0.1);
    R.set(1, 1, 2.2);
    R.set(1, 2, 0.4);
    R.set(2, 2, 1.9);

    auto out = lapack::Lapack_C_h_R_C_h(C, R);

    auto Maybe_L = cholesky(C);
    REQUIRE(Maybe_L.valid());
    Matrix<double> L = static_cast<Matrix<double> const&>(Maybe_L.value());
    auto Cref = L * to_dense(R) * tr(L);

    auto out_dense = to_dense(out);
    REQUIRE(rel_residual_1(Cref, out_dense) < 1e-10);
    REQUIRE(symmetry_defect_1(out_dense) < 1e-12);
}
