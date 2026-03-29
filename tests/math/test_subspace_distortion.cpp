#include "catch_amalgamated.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "lapack_headers.h"
#include "matrix.h"

using ::Matrix;
using ::SymmetricMatrix;
using ::SymPosDefMatrix;

namespace {

double rel_residual_1(const Matrix<double>& Aref, const Matrix<double>& A) {
    using ::norm_1;
    auto N = norm_1(Aref);
    if (N == 0.0)
        return norm_1(A);
    return norm_1(A - Aref) / N;
}

SymPosDefMatrix<double> make_diag_psd(std::initializer_list<double> diag_vals) {
    const std::size_t n = diag_vals.size();
    SymPosDefMatrix<double> out(n, n, false);
    std::size_t i = 0;
    for (auto v : diag_vals) {
        out.set(i, i, v);
        ++i;
    }
    return out;
}

Matrix<double> make_diag_dense(std::initializer_list<double> diag_vals) {
    const std::size_t n = diag_vals.size();
    Matrix<double> out(n, n, 0.0);
    std::size_t i = 0;
    for (auto v : diag_vals) {
        out(i, i) = v;
        ++i;
    }
    return out;
}

Matrix<double> make_col_vector(std::initializer_list<double> vals) {
    Matrix<double> out(vals.size(), 1, false);
    std::size_t i = 0;
    for (auto v : vals) {
        out(i, 0) = v;
        ++i;
    }
    return out;
}

}  // namespace

TEST_CASE("idm_matrix_subspace matches idm_matrix on SPD inputs", "[idm][subspace]") {
    auto H = make_diag_psd({1.0, 2.0, 4.0});
    auto J = make_diag_psd({4.0, 6.0, 8.0});

    auto maybe_old = idm_matrix(H, J);
    auto maybe_new = idm_matrix_subspace(H, J);

    REQUIRE(maybe_old.valid());
    REQUIRE(maybe_new.valid());
    REQUIRE(rel_residual_1(to_dense(maybe_old.value()), to_dense(maybe_new.value())) < 1e-12);
}

TEST_CASE("dcc_matrix_subspace matches dcc_matrix on SPD inputs", "[dcc][subspace]") {
    auto H = make_diag_psd({1.0, 2.0, 4.0});
    auto J = make_diag_psd({4.0, 6.0, 8.0});

    auto maybe_old = dcc_matrix(H, J);
    auto maybe_new = dcc_matrix_subspace(H, J);

    REQUIRE(maybe_old.valid());
    REQUIRE(maybe_new.valid());
    REQUIRE(rel_residual_1(to_dense(maybe_old.value()), to_dense(maybe_new.value())) < 1e-12);
}

TEST_CASE("subspace IDM and DCC succeed when H is singular", "[idm][dcc][subspace]") {
    auto H = make_diag_psd({4.0, 1.0, 0.0});
    auto J = make_diag_psd({8.0, 3.0, 5.0});

    auto maybe_idm = idm_matrix_subspace(H, J);
    auto maybe_dcc = dcc_matrix_subspace(H, J);

    REQUIRE(maybe_idm.valid());
    REQUIRE(maybe_dcc.valid());

    auto Cref = make_diag_dense({2.0, 3.0, 0.0});
    auto Dref = make_diag_dense({0.5, 3.0, 0.0});

    REQUIRE(rel_residual_1(Cref, to_dense(maybe_idm.value())) < 1e-12);
    REQUIRE(rel_residual_1(Dref, to_dense(maybe_dcc.value())) < 1e-12);
}

TEST_CASE("correlation_distortion_matrix_subspace succeeds when J_sample is singular",
          "[correlation][subspace]") {
    auto J_sample = make_diag_psd({9.0, 0.0, 4.0});
    auto J_total = make_diag_psd({18.0, 5.0, 12.0});

    auto maybe_corr = correlation_distortion_matrix_subspace(J_sample, J_total);
    REQUIRE(maybe_corr.valid());

    auto Rref = make_diag_dense({2.0, 0.0, 3.0});
    REQUIRE(rel_residual_1(Rref, to_dense(maybe_corr.value())) < 1e-12);
}

TEST_CASE("sample/correlation recomposition is consistent on retained subspace",
          "[recompose][subspace]") {
    auto C_sample = make_diag_psd({4.0, 0.0, 9.0});
    auto R = make_diag_psd({2.0, 0.0, 3.0});

    auto maybe_recomposed = c_h_r_c_h_matrix_subspace(C_sample, R);
    REQUIRE(maybe_recomposed.valid());

    auto Crecomp_ref = make_diag_dense({8.0, 0.0, 27.0});
    REQUIRE(rel_residual_1(Crecomp_ref, to_dense(maybe_recomposed.value())) < 1e-12);
}

TEST_CASE("distortion_induced_bias_subspace matches dense inverse on SPD inputs",
          "[dib][subspace]") {
    auto H = make_diag_psd({2.0, 4.0});
    auto g = make_col_vector({6.0, 8.0});

    auto maybe_dib = distortion_induced_bias_subspace(H, g);
    REQUIRE(maybe_dib.valid());

    auto bref = make_col_vector({3.0, 2.0});
    REQUIRE(rel_residual_1(bref, maybe_dib.value()) < 1e-12);
}

TEST_CASE("distortion_induced_bias_subspace projects onto retained H subspace",
          "[dib][subspace]") {
    auto H = make_diag_psd({4.0, 1.0, 0.0});
    auto g = make_col_vector({8.0, 3.0, 5.0});

    auto maybe_dib = distortion_induced_bias_subspace(H, g);
    REQUIRE(maybe_dib.valid());

    auto bref = make_col_vector({2.0, 3.0, 0.0});
    REQUIRE(rel_residual_1(bref, maybe_dib.value()) < 1e-12);
}

TEST_CASE("empty informative subspace returns Maybe_error", "[subspace][error]") {
    auto zero3 = make_diag_psd({0.0, 0.0, 0.0});
    auto J = make_diag_psd({1.0, 2.0, 3.0});
    auto g = make_col_vector({1.0, 2.0, 3.0});

    REQUIRE_FALSE(idm_matrix_subspace(zero3, J).valid());
    REQUIRE_FALSE(dcc_matrix_subspace(zero3, J).valid());
    REQUIRE_FALSE(sample_distortion_matrix_subspace(zero3, J).valid());
    REQUIRE_FALSE(correlation_distortion_matrix_subspace(zero3, J).valid());
    REQUIRE_FALSE(distortion_induced_bias_subspace(zero3, g).valid());
}

TEST_CASE("dimension mismatch returns Maybe_error", "[subspace][error]") {
    auto H2 = make_diag_psd({1.0, 2.0});
    auto J3 = make_diag_psd({1.0, 2.0, 3.0});
    auto g3 = make_col_vector({1.0, 2.0, 3.0});

    REQUIRE_FALSE(idm_matrix_subspace(H2, J3).valid());
    REQUIRE_FALSE(dcc_matrix_subspace(H2, J3).valid());
    REQUIRE_FALSE(sample_distortion_matrix_subspace(H2, J3).valid());
    REQUIRE_FALSE(correlation_distortion_matrix_subspace(H2, J3).valid());
    REQUIRE_FALSE(c_h_r_c_h_matrix_subspace(H2, J3).valid());
    REQUIRE_FALSE(distortion_induced_bias_subspace(H2, g3).valid());
}

TEST_CASE("non-finite distortion-induced bias inputs return Maybe_error", "[dib][subspace][error]") {
    auto H = make_diag_psd({2.0, 4.0});
    auto g = make_col_vector({1.0, std::numeric_limits<double>::quiet_NaN()});

    REQUIRE_FALSE(distortion_induced_bias_subspace(H, g).valid());
}

TEST_CASE("sqr_X false keeps symmetric matrices in matrix form", "[matrix][symmetric]") {
    SymPosDefMatrix<double> x(2, 2, false);
    x.set(0, 0, 2.0);
    x.set(0, 1, 3.0);
    x.set(1, 1, 5.0);

    auto squared = sqr_X<false>(x);

    REQUIRE(squared.nrows() == 2);
    REQUIRE(squared.ncols() == 2);
    CHECK(squared(0, 0) == Catch::Approx(4.0));
    CHECK(squared(0, 1) == Catch::Approx(9.0));
    CHECK(squared(1, 0) == Catch::Approx(9.0));
    CHECK(squared(1, 1) == Catch::Approx(25.0));
}

TEST_CASE("to_vector_half packs the upper triangle in order", "[matrix][symmetric]") {
    SymmetricMatrix<double> x(3);
    x.set(0, 0, 1.0);
    x.set(0, 1, 2.0);
    x.set(0, 2, 3.0);
    x.set(1, 1, 4.0);
    x.set(1, 2, 5.0);
    x.set(2, 2, 6.0);

    auto packed = to_vector_half(x);

    REQUIRE(packed.nrows() == 6);
    REQUIRE(packed.ncols() == 1);
    CHECK(packed[0] == Catch::Approx(1.0));
    CHECK(packed[1] == Catch::Approx(2.0));
    CHECK(packed[2] == Catch::Approx(3.0));
    CHECK(packed[3] == Catch::Approx(4.0));
    CHECK(packed[4] == Catch::Approx(5.0));
    CHECK(packed[5] == Catch::Approx(6.0));
}
