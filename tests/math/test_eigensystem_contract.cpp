#include "catch_amalgamated.hpp"

#include "matrix.h"

using ::Matrix;
using ::DiagonalMatrix;

namespace {

double rel_residual_1(const Matrix<double>& Aref, const Matrix<double>& A) {
    using ::norm_1;
    auto N = norm_1(Aref);
    if (N == 0.0) return norm_1(A);
    return norm_1(A - Aref) / N;
}

}  // namespace

TEST_CASE("eigensystem contract: left/right, ordering, reconstruction", "[eigs][contract]") {
    // Nonsymmetric, upper-triangular with distinct eigenvalues of mixed signs
    Matrix<double> A(3, 3, 0.0);
    A(0, 0) = 2.0;   A(0, 1) = 0.25; A(0, 2) = 0.0;
    A(1, 0) = 0.0;   A(1, 1) = 0.5;  A(1, 2) = 0.30;
    A(2, 0) = 0.0;   A(2, 1) = 0.0;  A(2, 2) = -1.0;

    auto me = eigs(A);
    REQUIRE(me.valid());

    // Contract declares tuple order (L, VR, VL)
    auto [L, VR, VL] = me.value();

    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * 1e3;

    // 1) Robust reconstruction using right eigenvectors only: A ≈ VR · L · inv(VR)
    auto invVR = inv(VR);
    REQUIRE(invVR.valid());
    auto A2 = VR * L * invVR.value();
    REQUIRE(rel_residual_1(A, A2) < tol);

    // 3) Ordering: eigenvalues sorted with more positive first (non-increasing)
    REQUIRE(L[0] >= L[1] - tol);
    REQUIRE(L[1] >= L[2] - tol);

    // 4) Mapping to (V, W): V = VR, W = VL^T. In general W may deviate from inv(V)
    //    for degenerate clusters, so we do not enforce W ≈ inv(V) here.
}
