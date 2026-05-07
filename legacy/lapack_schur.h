#ifndef LAPACK_SCHUR_H
#define LAPACK_SCHUR_H

// Real Schur decomposition wrapper, Q = U · T · U^T with U orthogonal and T
// upper quasi-triangular (1×1 and 2×2 diagonal blocks). LAPACK routine: dgees.
//
// Used by the planned Schur+Parlett path for calc_Qdt; kept in its own header
// so it can be included narrowly without dragging in the rest of LAPACK_HEADERS.

#include <string>
#include <vector>

#include "matrix.h"
#include "maybe_error.h"

namespace lapack {

extern "C" void dgees_(char* JOBVS, char* SORT,
                       int (*SELECT)(double*, double*),  // unused when SORT='N'
                       int* N, double* A, int* LDA, int* SDIM,
                       double* WR, double* WI,
                       double* VS, int* LDVS,
                       double* WORK, int* LWORK,
                       int* BWORK, int* INFO);

// Compute the real Schur decomposition Q = U · T · U^T.
//
// Returns std::tuple<T, U> with:
//   - T: row-major, upper quasi-triangular (1×1 blocks for real eigenvalues,
//        2×2 blocks for complex-conjugate pairs).
//   - U: row-major, orthogonal (U · U^T = I).
//
// Identity satisfied in row-major matmul: Q == U · T · U^T to ~ machine eps.
//
// Implementation note: dgees works on the matrix as LAPACK sees it (column-
// major). Our row-major Q reads as Q^T column-major, so we pre-transpose Q
// before passing to LAPACK and post-transpose the outputs. The element-wise
// quasi-triangular structure of T is preserved through the transpose.
inline Maybe_error<std::tuple<Matrix<double>, Matrix<double>>>
Lapack_Schur(Matrix<double> const& Q) {
    int N = static_cast<int>(Q.nrows());
    if (N <= 0 || N != static_cast<int>(Q.ncols())) {
        return error_message(std::string{"Lapack_Schur: input must be square and non-empty (got "}
                              + std::to_string(Q.nrows()) + "x" + std::to_string(Q.ncols()) + ")");
    }

    // Pre-transpose: row-major view of tr(Q) is column-major view of Q.
    Matrix<double> A = tr(Q);

    Matrix<double> U(static_cast<std::size_t>(N), static_cast<std::size_t>(N), 0.0);
    std::vector<double> WR(static_cast<std::size_t>(N), 0.0);
    std::vector<double> WI(static_cast<std::size_t>(N), 0.0);

    char JOBVS = 'V';  // compute Schur vectors
    char SORT = 'N';   // no sorting
    int SDIM = 0;
    int LDA = N;
    int LDVS = N;
    int INFO = 0;

    // Workspace query: LWORK = -1 makes dgees report the optimal size in WORK[0].
    int LWORK = -1;
    double WORK_QUERY = 0.0;
    dgees_(&JOBVS, &SORT, /*SELECT=*/nullptr, &N, &A(0, 0), &LDA, &SDIM,
           WR.data(), WI.data(), &U(0, 0), &LDVS,
           &WORK_QUERY, &LWORK, /*BWORK=*/nullptr, &INFO);
    if (INFO != 0) {
        return error_message(std::string{"Lapack_Schur: workspace query failed, INFO="} +
                              std::to_string(INFO));
    }
    LWORK = static_cast<int>(WORK_QUERY);
    if (LWORK < 3 * N) LWORK = 3 * N;  // dgees minimum requirement
    std::vector<double> WORK(static_cast<std::size_t>(LWORK), 0.0);

    dgees_(&JOBVS, &SORT, /*SELECT=*/nullptr, &N, &A(0, 0), &LDA, &SDIM,
           WR.data(), WI.data(), &U(0, 0), &LDVS,
           WORK.data(), &LWORK, /*BWORK=*/nullptr, &INFO);

    if (INFO < 0) {
        return error_message(std::string{"Lapack_Schur: argument "} +
                              std::to_string(-INFO) + " had an illegal value");
    }
    if (INFO > 0) {
        return error_message(std::string{"Lapack_Schur: QR algorithm failed to converge, INFO="} +
                              std::to_string(INFO));
    }

    // Un-transpose: A's memory now holds T_F column-major; row-major view is T_F^T.
    // tr() flips axes so T_cpp is element-wise equal to T_F (upper quasi-triangular
    // in row-major indexing).
    auto T_cpp = tr(A);
    auto U_cpp = tr(U);

    return std::tuple<Matrix<double>, Matrix<double>>{std::move(T_cpp), std::move(U_cpp)};
}

}  // namespace lapack

#endif  // LAPACK_SCHUR_H
