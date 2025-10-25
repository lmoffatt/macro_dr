#include "catch_amalgamated.hpp"
#include "matrix.h"
#include "lapack_headers.h"

using ::DiagonalMatrix;
using ::Matrix;

namespace {

double rel_residual_1(const Matrix<double>& Aref, const Matrix<double>& A) {
    using ::norm_1;
    auto N = norm_1(Aref);
    if (N == 0.0)
        return norm_1(A);
    return norm_1(A - Aref) / N;
}

}  // namespace

TEST_CASE("eigs reconstructs nonsymmetric matrix via right eigenvalues", "[eigs]") {
    // A fixed, nonsymmetric 3x3 with distinct eigenvalues
    Matrix<double> A(3, 3, 0.0);
    A(0, 0) = -30;
    A(0, 1) = 10;
    A(0, 2) = 20;
    A(1, 0) = 50;
    A(1, 1) = -70;
    A(1, 2) = 20;
    A(2, 0) = 40;
    A(2, 1) = 60;
    A(2, 2) = -100;

    auto me = eigs(A);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    auto A2 = VR * L * inv(VR).value();

    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * 1e3;
    REQUIRE(rel_residual_1(A, A2) < tol);
}

TEST_CASE("eigs reconstructs via left eigenvectors", "[eigs][left]") {
    // Same nonsymmetric matrix
    Matrix<double> A(3, 3, 0.0);
    A(0, 0) = -30;
    A(0, 1) = 10;
    A(0, 2) = 20;
    A(1, 0) = 50;
    A(1, 1) = -70;
    A(1, 2) = 20;
    A(2, 0) = 40;
    A(2, 1) = 60;
    A(2, 2) = -100;

    auto me = eigs(A);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    auto Arec = inv(tr(VL)).value() * L * tr(VL);
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * 1e3;
    REQUIRE(rel_residual_1(A, Arec) < tol);
}

TEST_CASE("eigs reconstructs nonsymmetric matrix", "[eigs]") {
    // A fixed, nonsymmetric 3x3 with distinct eigenvalues
    Matrix<double> A(3, 3, 0.0);
    A(0, 0) = -30;
    A(0, 1) = 10;
    A(0, 2) = 20;
    A(1, 0) = 50;
    A(1, 1) = -70;
    A(1, 2) = 20;
    A(2, 0) = 40;
    A(2, 1) = 60;
    A(2, 2) = -100;

    auto me = eigs(A);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    auto A2 = VR * L * tr(VL);

    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * 1e3;
    REQUIRE(rel_residual_1(A, A2) < tol);
}

TEST_CASE("eigs detects zero-mode for generator Q", "[eigs][Q]") {
    // Small CTMC generator with distinct eigenvalues
    Matrix<double> Q(3, 3, 0.0);
    // Off-diagonals
    Q(0, 1) = 0.7;
    Q(0, 2) = 0.3;
    Q(1, 0) = 0.5;
    Q(1, 2) = 1.0;
    Q(2, 0) = 0.25;
    Q(2, 1) = 0.75;
    // Diagonals to enforce row sums zero
    for (std::size_t i = 0; i < 3; ++i) {
        double row = 0.0;
        for (std::size_t j = 0; j < 3; ++j)
            if (j != i)
                row += Q(i, j);
        Q(i, i) = -row;
    }

    auto me = eigs(Q);
    REQUIRE(me.valid());
    auto [L, VR, VL] = me.value();

    // Find eigenvalue closest to zero
    std::size_t k0 = 0;
    double amin = std::abs(L[0]);
    for (std::size_t i = 1; i < L.size(); ++i) {
        if (std::abs(L[i]) < amin) {
            amin = std::abs(L[i]);
            k0 = i;
        }
    }
    // Zero eigenvalue within tolerance
    REQUIRE(amin < 1e-12);

    // Right eigenvector should align with ones vector
    auto u0 = VR(":", k0);
    Matrix<double> ones(u0.nrows(), u0.ncols(), 1.0);
    auto dot = getvalue(tr(ones) * u0);
    auto nu = std::sqrt(getvalue(tr(u0) * u0));
    auto n1 = std::sqrt(getvalue(tr(ones) * ones));
    double corr = std::abs(dot) / (nu * n1);
    REQUIRE(corr > 1.0 - 1e-10);
}
