#ifndef LAPACK_HEADERS_H
#define LAPACK_HEADERS_H

#include "matrix.h"
#include "maybe_error.h"
#include <algorithm>
#include <iostream>
#include <sstream>

namespace lapack {

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
                       double *ALPHA, double *A, int *LDA, double *B, int *LDB,
                       double *BETA, double *C, int *LDC);

extern "C" void dgeqrf_(int *M, int *N, double *A, int *LDA, double *TAU,
                        double *WORK, int *LWORK, int *INFO);

extern "C" void dorgqr_(int *M, int *N, int *K, double *A, int *LDA,
                        double *TAU, double *WORK, int *LWORK, int *INFO);

inline auto Lapack_UT(const Matrix<double> &x) {
  Matrix<double> out(x.nrows(), x.ncols(), false);
  for (std::size_t i = 0; i < out.ncols(); ++i)
    for (std::size_t j = 0; j < out.nrows(); ++j)
      out(j, i) = j >= i ? x(j, i) : 0;
  return out;
}

inline auto Lapack_LT(const Matrix<double> &x) {
  Matrix<double> out(x.nrows(), x.ncols(), false);
  for (std::size_t i = 0; i < x.ncols(); ++i)
    for (std::size_t j = 0; j < x.nrows(); ++j)
      out(j, i) = i >= j ? x(j, i) : 0;
  return out;
}

inline std::pair<Matrix<double>, Matrix<double>> Lapack_QR(const Matrix<double> &x) {
  int M = x.nrows();
  /*    [in]	M

                M is INTEGER
                    The number of rows of the matrix A.  M >= 0.
*/

  int N = x.ncols();
  /*    [in]	N

              N is INTEGER
              The number of columns of the matrix A.  N >= 0.

*/
  auto a = tr(x);
  double &A = a[0];

  /*
[in,out]	A

          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit, the elements on and above the diagonal of the array
          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          upper triangular if m >= n); the elements below the diagonal,
          with the array TAU, represent the orthogonal matrix Q as a
          product of min(m,n) elementary reflectors (see Further
          Details).
*/

  int LDA = M;
  /*

[in]	LDA

          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

*/
  Matrix<double> tau(std::min(M, N), 1, false);
  double &TAU = tau[0];

  /*
[out]	TAU

          TAU is DOUBLE PRECISION array, dimension (min(M,N))
          The scalar factors of the elementary reflectors (see Further
          Details).

*/

  Matrix<double> work(N * M, 1, false);
  double &WORK = work[0];
  int LWORK = N * M;
  /*
[out]	WORK

          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

*/
  /*
[in]	LWORK

          LWORK is INTEGER
          The dimension of the array WORK.
          LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
          For optimum performance LWORK >= N*NB, where NB is
          the optimal blocksize.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
*/

  int INFO;
  /*

[out]	INFO

          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value

*/

  dgeqrf_(&M, &N, &A, &LDA, &TAU, &WORK, &LWORK, &INFO);

  /*

Purpose:

   DGEQRF computes a QR factorization of a real M-by-N matrix A:

      A = Q * ( R ),
              ( 0 )

   where:

      Q is a M-by-M orthogonal matrix;
      R is an upper-triangular N-by-N matrix;
      0 is a (M-N)-by-N zero matrix, if M > N.

*/

  auto r = tr(Lapack_UT(a));

  int K = std::min(N, M);

  dorgqr_(&M, &N, &K, &A, &LDA, &TAU, &WORK, &LWORK, &INFO);

  auto q = tr(a);
  return std::pair(std::move(q), std::move(r));
}

inline Matrix<double> &Lapack_Full_Product(const Matrix<double> &x,
                                    const Matrix<double> &y, Matrix<double> &z,
                                    bool transpose_x, bool transpose_y,
                                    double alpha = 1.0, double beta = 0.0) {

  std::size_t cols_i, rows_e, cols_e;
  if (transpose_x) {
    rows_e = x.ncols();
    cols_i = x.nrows();
  } else {
    rows_e = x.nrows();
    cols_i = x.ncols();
  }

  if (transpose_y) {
    //  rows_i=y.ncols();
    cols_e = y.nrows();
  } else {
    //  rows_i=y.nrows();
    cols_e = y.ncols();
  }

  // assert(rows_i==cols_i);
  assert(((transpose_y ? y.ncols() : y.nrows()) == cols_i));
  // First it has to find out if the last dimension of x matches the
  // first of y
  // now we build the M_Matrix result
  if (z.nrows() != rows_e || z.ncols() != cols_e)
    z = Matrix<double>(rows_e, cols_e);

  /***  as fortran uses the reverse order for matrices and we want to
    avoid a copying operation, we calculate
        Transpose(Z)=Transpose(y)*Transpose(x)

        Transpose(matrix)=just plain matrix in C++ format


    */
  char TRANSA;
  char TRANSB;

  if (transpose_y)
    TRANSA = 'T';
  else
    TRANSA = 'N';

  if (transpose_x)
    TRANSB = 'T';
  else
    TRANSB = 'N';

  int M = cols_e;
  int N = rows_e;
  int K = cols_i;

  double ALPHA = alpha;
  double *A = const_cast<double *>(&y[0]);
  int LDA;
  if (transpose_y)
    LDA = K;
  else
    LDA = M;

  double *B = const_cast<double *>(&x[0]);

  int LDB;
  if (transpose_x)
    LDB = N;
  else
    LDB = K;

  double BETA = beta;

  double *C = &z[0];

  int LDC = M;

  try {
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C,
           &LDC);
  } catch (...) {
    assert(false);
  }
  return z;
}

inline Matrix<double> Lapack_Full_Product(const Matrix<double> &x,
                                   const Matrix<double> &y, bool transpose_x,
                                   bool transpose_y, double alpha,
                                   double beta) {
  using lapack::dgemm_;

  std::size_t cols_i, rows_e, cols_e;
  if (transpose_x) {
    rows_e = x.ncols();
    cols_i = x.nrows();
  } else {
    rows_e = x.nrows();
    cols_i = x.ncols();
  }

  if (transpose_y) {
    //  rows_i=y.ncols();
    cols_e = y.nrows();
  } else {
    //  rows_i=y.nrows();
    cols_e = y.ncols();
  }

  // assert(rows_i==cols_i);
  assert(((transpose_y ? y.ncols() : y.nrows()) == cols_i));
  // First it has to find out if the last dimension of x matches the
  // first of y
  // now we build the M_Matrix result
  Matrix<double> z(rows_e, cols_e);

  /***  as fortran uses the reverse order for matrices and we want to
    avoid a copying operation, we calculate
        Transpose(Z)=Transpose(y)*Transpose(x)

        Transpose(matrix)=just plain matrix in C++ format


    */
  char TRANSA;
  char TRANSB;

  if (transpose_y)
    TRANSA = 'T';
  else
    TRANSA = 'N';

  if (transpose_x)
    TRANSB = 'T';
  else
    TRANSB = 'N';

  int M = cols_e;
  int N = rows_e;
  int K = cols_i;

  double ALPHA = alpha;
  double *A = const_cast<double *>(&y[0]);
  int LDA;
  if (transpose_y)
    LDA = K;
  else
    LDA = M;

  double *B = const_cast<double *>(&x[0]);

  int LDB;
  if (transpose_x)
    LDB = N;
  else
    LDB = K;

  double BETA = beta;

  double *C = &z[0];

  int LDC = M;

  try {
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C,
           &LDC);
  } catch (...) {
    assert(false);
  }
  return z;
}

extern "C" void dsymm_(char *SIDE, char *UPLO, int *M, int *N, double *ALPHA,
                       double *A, int *LDA, double *B, int *LDB, double *BETA,
                       double *C, int *LDC);

inline Matrix<double> Lapack_Sym_Product(const SymmetricMatrix<double> &x,
                                  const Matrix<double> &y,
                                  bool first_symmetric_c, double alpha,
                                  double beta) {

  /*
   DSYMM  performs one of the matrix-matrix operations

      C := alpha*A*B + beta*C,

   or

      C := alpha*B*A + beta*C,

   where alpha and beta are scalars,  A is a symmetric matrix and  B and
   C are  m by n matrices.

Parameters
*/

  auto out = first_symmetric_c ? Matrix<double>(x.nrows(), y.ncols())
                               : Matrix<double>(y.nrows(), x.ncols());

  char SIDE = first_symmetric_c ? 'R' : 'L';

  /*  [in]	SIDE

            SIDE is CHARACTER*1
             On entry,  SIDE  specifies whether  the  symmetric matrix  A
             appears on the  left or right  in the  operation as follows:

                SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,

                SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*/

  char UPLO = 'L';
  /*
  [in]	UPLO

            UPLO is CHARACTER*1
             On  entry,   UPLO  specifies  whether  the  upper  or  lower
             triangular  part  of  the  symmetric  matrix   A  is  to  be
             referenced as follows:

                UPLO = 'U' or 'u'   Only the upper triangular part of the
                                    symmetric matrix is to be referenced.

                UPLO = 'L' or 'l'   Only the lower triangular part of the
                                    symmetric matrix is to be referenced.


*/
  int M = first_symmetric_c ? y.ncols() : x.ncols();
  /*

  [in]	M

            M is INTEGER
             On entry,  M  specifies the number of rows of the matrix  C.
             M  must be at least zero.

*/

  int N = first_symmetric_c ? x.nrows() : y.nrows();

  /*

  [in]	N

            N is INTEGER
             On entry, N specifies the number of columns of the matrix C.
             N  must be at least zero.

*/
  double ALPHA = alpha;
  /*

  [in]	ALPHA

            ALPHA is DOUBLE PRECISION.
             On entry, ALPHA specifies the scalar alpha.

*/
  
  
  double &A = const_cast<double &>(x[0]);

  /*

  [in]	A

            A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
             m  when  SIDE = 'L' or 'l'  and is  n otherwise.
             Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
             the array  A  must contain the  symmetric matrix,  such that
             when  UPLO = 'U' or 'u', the leading m by m upper triangular
             part of the array  A  must contain the upper triangular part
             of the  symmetric matrix and the  strictly  lower triangular
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
             the leading  m by m  lower triangular part  of the  array  A
             must  contain  the  lower triangular part  of the  symmetric
             matrix and the  strictly upper triangular part of  A  is not
             referenced.
             Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
             the array  A  must contain the  symmetric matrix,  such that
             when  UPLO = 'U' or 'u', the leading n by n upper triangular
             part of the array  A  must contain the upper triangular part
             of the  symmetric matrix and the  strictly  lower triangular
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
             the leading  n by n  lower triangular part  of the  array  A
             must  contain  the  lower triangular part  of the  symmetric
             matrix and the  strictly upper triangular part of  A  is not
             referenced.
*/

  int LDA = x.nrows();

  /*
  [in]	LDA

            LDA is INTEGER
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
             LDA must be at least  max( 1, m ), otherwise  LDA must be at
             least  max( 1, n ).

*/
  double &B = const_cast<double &>(y[0]);
  /*
  [in]	B

            B is DOUBLE PRECISION array, dimension ( LDB, N )
             Before entry, the leading  m by n part of the array  B  must
             contain the matrix B.
*/

  int LDB = M;

  /*

  [in]	LDB

            LDB is INTEGER
             On entry, LDB specifies the first dimension of B as declared
             in  the  calling  (sub)  program.   LDB  must  be  at  least
             max( 1, m ).
*/
  double BETA = beta;

  /*
  [in]	BETA

            BETA is DOUBLE PRECISION.
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is
             supplied as zero then C need not be set on input.

*/

  double &C = out[0];

  /*

  [in,out]	C

            C is DOUBLE PRECISION array, dimension ( LDC, N )
             Before entry, the leading  m by n  part of the array  C must
             contain the matrix  C,  except when  beta  is zero, in which
             case C need not be set on entry.
             On exit, the array  C  is overwritten by the  m by n updated
             matrix.

*/

  int LDC = M;
  /*


  [in]	LDC

            LDC is INTEGER
             On entry, LDC specifies the first dimension of C as declared
             in  the  calling  (sub)  program.   LDC  must  be  at  least
             max( 1, m ).

*/

  dsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, &A, &LDA, &B, &LDB, &BETA, &C, &LDC);
  return out;
}

extern "C" void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                        double *RCOND, double *WORK, int *IWORK, int *INFO);
extern "C" void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV,
                        int *INFO);
extern "C" void dgetri_(int *n, double *B, int *dla, int *ipiv, double *work1,
                        int *lwork, int *info);
extern "C" double dlange_(char *NORM, int *M, int *N, double *A, int *LDA,
                          double *WORK);

inline Maybe_error<Matrix<double>> Lapack_Full_inv(const Matrix<double> &a);
} // namespace lapack
template <> constexpr std::string function_name<&lapack::Lapack_Full_inv>() {
  return "Lapack_Full_inv";
}

namespace lapack {

inline Maybe_error<Matrix<double>> Lapack_Full_inv(const Matrix<double> &a)

{
  return_error<Matrix<double>, Lapack_Full_inv> Error;
  const double min_inv_condition_number = 1e-12;
  // using Matrix_Unary_Transformations::Transpose;
  using lapack::dgecon_;
  using lapack::dgetrf_;
  using lapack::dgetri_;
  using lapack::dlange_;
  if (a.size() == 0)
    return Error("EMPTY MATRIX");
  else {
    assert(a.ncols() == a.nrows());

    char NORM = '1';
    int N = a.ncols();
    int M = N;

    int INFO = 0;
    //  char msg[101];
    auto IPIV = std::make_unique<int[]>(N);
    int LWORK;
    //        M_Matrix<double> B=Transpose(a);
    Matrix<double> B = a;
    int LDA = N;
    // A=new double[n*n];
    double *A = &B[0]; // more efficient code
    auto WORK_lange = std::make_unique<double[]>(N);

    double ANORM = dlange_(&NORM, &M, &N, A, &LDA, WORK_lange.get());

    dgetrf_(&N, &M, A, &LDA, IPIV.get(), &INFO);

    double RCOND;
    auto WORK_cond = std::make_unique<double[]>(N * 4);
    auto IWORK = std::make_unique<int[]>(N);
    int INFO_con;

    dgecon_(&NORM, &N, A, &LDA, &ANORM, &RCOND, WORK_cond.get(), IWORK.get(),
            &INFO_con);

    LWORK = N * N;
    Matrix<double> W(N, N);
    double *WORK = &W[0];

    dgetri_(&N, A, &LDA, IPIV.get(), WORK, &LWORK, &INFO);

    if (RCOND < min_inv_condition_number)
      return Error("bad condition number RCOND=" + std::to_string(RCOND));
    if (INFO == 0)
      return B;
    //  return Op({B,RCOND});
    else
      return Error("Singular Matrix on i=" + std::to_string(INFO));
    ;
  }
}

extern "C" void dsycon_(char *UPLO, int *N,
                        double * /* precision, dimension( lda, * ) */ A,
                        int *LDA, int * /*, dimension( * ) */ IPIV,
                        double *ANORM, double *RCOND,
                        double * /* dimension( * )  */ WORK,
                        int * /* dimension( * ) */ IWORK, int *INFO);

extern "C" void dsytrf_(char *UPLO, int *N, double *A, int *LDA, int *IPIV,
                        double *WORK, int *LWORK, int *INFO);

extern "C" void dsytri_(char *UPLO, int *N, double * /*dimension( lda, * )*/ A,
                        int *LDA, int * /* dimension( * ) */ IPIV,
                        double * /*dimension( * )*/ WORK, int *INFO);

Maybe_error<SymmetricMatrix<double>>
Lapack_Symm_inv(const SymmetricMatrix<double> &a);
Maybe_error<SymPosDefMatrix<double>>
Lapack_SymmPosDef_inv(const SymPosDefMatrix<double> &a);

Maybe_error<DownTrianMatrix<double>>
Lapack_chol(const SymPosDefMatrix<double> &x);

inline Maybe_error<
    std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>>
Lapack_EigenSystem(const Matrix<double> &x, bool does_permutations,
                   bool does_diagonal_scaling,
                   bool computes_eigenvalues_condition_numbers,
                   bool computes_eigenvectors_condition_numbers);

} // namespace lapack
template <> constexpr std::string function_name<&lapack::Lapack_Symm_inv>() {
  return "Lapack_Symm_inv";
}

template <> constexpr std::string function_name<&lapack::Lapack_chol>() {
  return "Lapack_chol";
}

template <>
constexpr std::string function_name<&lapack::Lapack_SymmPosDef_inv>() {
  return "Lapack_SymmPosDef_inv";
}

template <>
constexpr std::string
function_name<static_cast<Maybe_error<SymPosDefMatrix<double>> (*)(
    const DownTrianMatrix<double> &)>(lapack::Lapack_LT_Cholesky_inv)>() {
  return "Lapack_LT_Cholesky_inv";
}

template <>
constexpr std::string
function_name<static_cast<Maybe_error<SymPosDefMatrix<double>> (*)(
    const UpTrianMatrix<double> &)>(lapack::Lapack_UT_Cholesky_inv)>() {
  return "Lapack_UT_Cholesky_inv";
}
template <> constexpr std::string function_name<&lapack::Lapack_LT_inv>() {
  return "Lapack_LT_inv";
}

template <> constexpr std::string function_name<&lapack::Lapack_UT_inv>() {
  return "Lapack_UT_inv";
}
template <> constexpr std::string function_name<&lapack::Lapack_EigenSystem>() {
  return "Lapack_EigenSystem";
}

namespace lapack {

inline Maybe_error<SymmetricMatrix<double>>
Lapack_Symm_inv(const SymmetricMatrix<double> &a) {
  return_error<SymmetricMatrix<double>, Lapack_Full_inv> Error;

  if (a.size() == 0)
    return Error("EMPTY MATRIX");
  else {
    assert(a.nrows() == a.ncols());

    /**
Purpose:

DSYTRF computes the factorization of a real symmetric matrix A using
the Bunch-Kaufman diagonal pivoting method.  The form of the
factorization is

  A = U*D*U**T  or  A = L*D*L**T

where U (or L) is a product of permutation and unit upper (lower)
triangular matrices, and D is symmetric and block diagonal with
1-by-1 and 2-by-2 diagonal blocks.

This is the blocked version of the algorithm, calling Level 3 BLAS.

Parameters
[in]	UPLO

        UPLO is CHARACTER*1
        = 'U':  Upper triangle of A is stored;
        = 'L':  Lower triangle of A is stored.

[in]	N

        N is INTEGER
        The order of the matrix A.  N >= 0.

[in,out]	A

        A is DOUBLE PRECISION array, dimension (LDA,N)
        On entry, the symmetric matrix A.  If UPLO = 'U', the leading
        N-by-N upper triangular part of A contains the upper
        triangular part of the matrix A, and the strictly lower
        triangular part of A is not referenced.  If UPLO = 'L', the
        leading N-by-N lower triangular part of A contains the lower
        triangular part of the matrix A, and the strictly upper
        triangular part of A is not referenced.

        On exit, the block diagonal matrix D and the multipliers used
        to obtain the factor U or L (see below for further details).

[in]	LDA

        LDA is INTEGER
        The leading dimension of the array A.  LDA >= max(1,N).

[out]	IPIV

        IPIV is INTEGER array, dimension (N)
        Details of the interchanges and the block structure of D.
        If IPIV(k) > 0, then rows and columns k and IPIV(k) were
        interchanged and D(k,k) is a 1-by-1 diagonal block.
        If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
        columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
        is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
        IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
        interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.

[out]	WORK

        WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
        On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

[in]	LWORK

        LWORK is INTEGER
        The length of WORK.  LWORK >=1.  For best performance
        LWORK >= N*NB, where NB is the block size returned by ILAENV.

        If LWORK = -1, then a workspace query is assumed; the routine
        only calculates the optimal size of the WORK array, returns
        this value as the first entry of the WORK array, and no error
        message related to LWORK is issued by XERBLA.

[out]	INFO

        INFO is INTEGER
        = 0:  successful exit
        < 0:  if INFO = -i, the i-th argument had an illegal value
        > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
              has been completed, but the block diagonal matrix D is
              exactly singular, and division by zero will occur if it
              is used to solve a system of equations.

*/

    char UPLO = 'L';

    int INFO = 0;
    int N = a.ncols();
    auto IPIV = std::make_unique<int[]>(N);
    int LWORK = N * N; //
    SymmetricMatrix<double> B(a);

    int LDA = N;
    double *A = &B[0]; // more efficient code
    Matrix<double> W(N, N);
    double *WORK = &W[0];
    dsytrf_(&UPLO, &N, A, &LDA, IPIV.get(), WORK, &LWORK, &INFO);
    double RCOND;
    auto WORK_cond = std::make_unique<double[]>(N * 4);
    auto IWORK = std::make_unique<int[]>(N);
    int INFO_con;

    /**
DSYCON

Download DSYCON + dependencies [TGZ] [ZIP] [TXT]

Purpose:

DSYCON estimates the reciprocal of the condition number (in the
1-norm) of a real symmetric matrix A using the factorization
A = U*D*U**T or A = L*D*L**T computed by DSYTRF.

An estimate is obtained for norm(inv(A)), and the reciprocal of the
condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

Parameters
[in]	UPLO

  UPLO is CHARACTER*1
  Specifies whether the details of the factorization are stored
  as an upper or lower triangular matrix.
  = 'U':  Upper triangular, form is A = U*D*U**T;
  = 'L':  Lower triangular, form is A = L*D*L**T.

[in]	N

  N is INTEGER
  The order of the matrix A.  N >= 0.

[in]	A

  A is DOUBLE PRECISION array, dimension (LDA,N)
  The block diagonal matrix D and the multipliers used to
  obtain the factor U or L as computed by DSYTRF.

[in]	LDA

  LDA is INTEGER
  The leading dimension of the array A.  LDA >= max(1,N).

[in]	IPIV

  IPIV is INTEGER array, dimension (N)
  Details of the interchanges and the block structure of D
  as determined by DSYTRF.

[in]	ANORM

  ANORM is DOUBLE PRECISION
  The 1-norm of the original matrix A.

[out]	RCOND

  RCOND is DOUBLE PRECISION
  The reciprocal of the condition number of the matrix A,
  computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
  estimate of the 1-norm of inv(A) computed in this routine.

[out]	WORK

  WORK is DOUBLE PRECISION array, dimension (2*N)

[out]	IWORK

  IWORK is INTEGER array, dimension (N)

[out]	INFO

  INFO is INTEGER
  = 0:  successful exit
  < 0:  if INFO = -i, the i-th argument had an illegal value*/

    char NORM = '1';
    int M = N;

    auto WORK_lange = std::make_unique<double[]>(N);

    double ANORM = dlange_(&NORM, &M, &N, A, &LDA, WORK_lange.get());

    //        dsycon_( 	char*  	UPLO, int *  	N,double */* precision,
    //        dimension( lda, * ) */  	A,
    //                                    int *  	LDA,double * ANORM,
    //                                    double *   	RCOND, double */*
    //                                    dimension( * )  */	WORK, int */*
    //                                    dimension( * ) */ 	IWORK, int *
    //                                    INFO);

    dsycon_(&UPLO, &N, A, &LDA, IPIV.get(), &ANORM, &RCOND, WORK_cond.get(),
            IWORK.get(), &INFO_con);

    if (INFO < 0) {
      std::string argNames[] = {"UPLO",  "N",     "A",         "LDA",   "IPIV",
                                "ANORM", "RCOND", "WORK_cond", "IWORK", "INFO"};
      return Error("INVALID ARGUMENT " + std::to_string(INFO) +
                   argNames[-INFO]);
    } else if (INFO > 0) {
      return Error("SINGULAR MATRIX ON " + std::to_string(INFO));
    } else {
      /**
 * dsytri()
subroutine dsytri 	( 	character  	UPLO,
      integer  	N,
      double precision, dimension( lda, * )  	A,
      integer  	LDA,
      integer, dimension( * )  	IPIV,
      double precision, dimension( * )  	WORK,
      integer  	INFO
)

DSYTRI

Download DSYTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

DSYTRI computes the inverse of a real symmetric indefinite matrix
A using the factorization A = U*D*U**T or A = L*D*L**T computed by
DSYTRF.

Parameters
[in]	UPLO

    UPLO is CHARACTER*1
    Specifies whether the details of the factorization are stored
    as an upper or lower triangular matrix.
    = 'U':  Upper triangular, form is A = U*D*U**T;
    = 'L':  Lower triangular, form is A = L*D*L**T.

[in]	N

    N is INTEGER
    The order of the matrix A.  N >= 0.

[in,out]	A

    A is DOUBLE PRECISION array, dimension (LDA,N)
    On entry, the block diagonal matrix D and the multipliers
    used to obtain the factor U or L as computed by DSYTRF.

    On exit, if INFO = 0, the (symmetric) inverse of the original
    matrix.  If UPLO = 'U', the upper triangular part of the
    inverse is formed and the part of A below the diagonal is not
    referenced; if UPLO = 'L' the lower triangular part of the
    inverse is formed and the part of A above the diagonal is
    not referenced.

[in]	LDA

    LDA is INTEGER
    The leading dimension of the array A.  LDA >= max(1,N).

[in]	IPIV

    IPIV is INTEGER array, dimension (N)
    Details of the interchanges and the block structure of D
    as determined by DSYTRF.

[out]	WORK

    WORK is DOUBLE PRECISION array, dimension (N)

[out]	INFO

    INFO is INTEGER
    = 0: successful exit
    < 0: if INFO = -i, the i-th argument had an illegal value
    > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
         inverse could not be computed.
 */
      dsytri_(&UPLO, &N, A, &LDA, IPIV.get(), WORK, &INFO);

      if (INFO != 0) {
        return Error("cannot invert a singular matrix " + std::to_string(INFO));
      } else {
        SymmetricMatrix<double> out(B.nrows());
        for (std::size_t i = 0; i < B.nrows(); ++i)
          for (std::size_t j = i; j < B.ncols(); ++j)
                out.set(i, j, B(j, i));

        /*    auto aa=a;
M_Matrix<double> test(a.nrows(),a.ncols(),Matrix_TYPE::FULL,a);

auto invtest=inv(test).first;

auto test_it=test*out;
auto test_a=aa*out;
auto invv=invtest*test;
*/
        copy_UT_to_LT(out);
        return out;
      }
    }
  }
}

extern "C" void dpocon_(char *UPLO, int *N,
                        double * /* precision, dimension( lda, * ) */ A,
                        int *LDA, double *ANORM, double *RCOND,
                        double * /* dimension( * )  */ WORK,
                        int * /* dimension( * ) */ IWORK, int *INFO);

extern "C" void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
extern "C" void dpotri_(char *UPLO, int *N, double *A, int *LDA, int *INFO);

inline Maybe_error<SymPosDefMatrix<double>>
Lapack_SymmPosDef_inv(const SymPosDefMatrix<double> &x) {
  return_error<SymPosDefMatrix<double>, Lapack_SymmPosDef_inv> Error;

  /**
 *
dpotrf()
subroutine dpotrf 	( 	character  	UPLO,
              integer  	N,
              double precision, dimension( lda, * )  	A,
              integer  	LDA,
              integer  	INFO
      )

DPOTRF

DPOTRF VARIANT: top-looking block version of the algorithm, calling Level 3
BLAS.

Download DPOTRF + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DPOTRF computes the Cholesky factorization of a real symmetric
   positive definite matrix A.

   The factorization has the form
      A = U**T * U,  if UPLO = 'U', or
      A = L  * L**T,  if UPLO = 'L',
   where U is an upper triangular matrix and L is lower triangular.

   This is the block version of the algorithm, calling Level 3 BLAS.

Parameters
*/
  char UPLO = 'L';

  /*
 * [in]	UPLO

            UPLO is CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.
*/

  int N = x.nrows();

  /*
  [in]	N

            N is INTEGER
            The order of the matrix A.  N >= 0.
*/
  auto a = x;
  double &A = a[0];
  /*
  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.

            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T.
*/
  int LDA = N;
  /*
 *
   [in]	LDA

             LDA is INTEGER
             The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;
  /*
  [out]	INFO

            INFO is INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

*/
  try {
    dpotrf_(&UPLO, &N, &A, &LDA, &INFO);
  } catch (...) {
    assert(false);
  }
  if (INFO < 0) {
    return Error(std::to_string(INFO) + "th argument for cholesky failed");
  } else if (INFO > 0) {
    return Error("the leading minor of order" + std::to_string(INFO) +
                 " is not positive definite, and the factorization "
                 "could not be completed.");

  } else {
    dpotri_(&UPLO, &N, &A, &LDA, &INFO);
    /*
DPOTRI

Download DPOTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DPOTRI computes the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF.

Parameters
*/
    /*
 * [in]	UPLO

            UPLO is CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.
*/
    /*
 *
  [in]	N

            N is INTEGER
            The order of the matrix A.  N >= 0.
*/
    /*
  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            DPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.

  [in]	LDA

            LDA is INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

  [out]	INFO

            INFO is INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the (i,i) element of the factor U or L is
                  zero, and the inverse could not be computed.*
*/

    if (INFO < 0) {
      return Error("the " + std::to_string(INFO) +
                   "th argument for inversehad an illegal value");
    } else if (INFO > 0) {
      return Error(" the (" + std::to_string(INFO) + "," +
                   std::to_string(INFO) +
                   ") element of the factor U or L is "
                   "zero, and the inverse could not be computed");
    } else {
      copy_UT_to_LT(a);
      return a;
    }
  }
}

extern "C" void
dgeevx_(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE, int *N,
        double * /* precision, dimension( lda, * ) */ A, int *LDA,
        double * /* precision, dimension( * ) */ WR,
        double * /* precision, dimension( * ) */ WI,
        double * /* precision, dimension( ldvl, * ) */ VL, int *LDVL,
        double * /* precision, dimension( ldvr, * ) */ VR, int *LDVR, int *ILO,
        int *IHI, double * /* precision, dimension( * ) */ SCALE, double *ABNRM,
        double * /* precision, dimension( * ) */ RCONDE,
        double * /* precision, dimension( * ) */ RCONDV,
        double * /* precision, dimension( * ) */ WORK, int *LWORK,
        int * /*dimension( * ) */ IWORK, int *INFO);

inline auto sort_by_eigenvalue(
    std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>> const
        &e) {
  auto &[VL, L, VR] = e;
  std::vector<std::pair<double, std::size_t>> la(L.size());
  for (std::size_t i = 0; i < L.size(); ++i)
    la[i] = std::pair(L[i], i);
  std::sort(la.begin(), la.end());
  DiagonalMatrix<double> Ls(L.nrows(), L.ncols());

  Matrix<double> VRs(VR.nrows(), VR.ncols());
  Matrix<double> VLs(VL.nrows(), VL.ncols());
  for (std::size_t i = 0; i < la.size(); ++i) {
    auto is = la[i].second;
    Ls[i] = L[is];
    for (std::size_t j = 0; j < la.size(); ++j) {
      VRs(j, i) = VR(j, is);
      VLs(i, j) = VL(is, j);
    }
  }
  
  return std::tuple(VRs, Ls, VLs);
}

inline Matrix<double> &
Right_Eigenvector_Nelson_Normalization(Matrix<double> &X) {
  assert(X.nrows() == X.ncols());
  auto n = X.nrows();
  for (std::size_t k = 0; k < n; ++k) {
    double max = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (std::abs(X(i, k)) >= std::abs(max)) {
        max = X(i, k);
      }
    for (std::size_t i = 0; i < n; ++i)
      X(i, k) = X(i, k) / max;
  }
  return X;
}
inline bool Nelson_Normalization(Matrix<double> &VR, Matrix<double> &VL) {
  Right_Eigenvector_Nelson_Normalization(VR);
  for (std::size_t i = 0; i < VR.nrows(); ++i) {
    double sum = 0;
    for (std::size_t j = 0; j < VR.ncols(); ++j)
      sum += VL(i, j) * VR(j, i);
    for (std::size_t j = 0; j < VR.ncols(); ++j)
      VL(i, j) = VL(i, j) / sum;
    return sum!=0;
  }
  //   assert((are_Equal<true, Matrix<double>>().test_prod(
  //       VR * VL, Matrix_Generators::eye<double>(VR.ncols()), std::cerr)));
}

inline Maybe_error<
    std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>>
Lapack_EigenSystem(const Matrix<double> &x, bool does_permutations,
                   bool does_diagonal_scaling,
                   bool computes_eigenvalues_condition_numbers,
                   bool computes_eigenvectors_condition_numbers) {

  return_error<
      std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>,
      Lapack_EigenSystem>
      Error;

  using lapack::dgeevx_;

  /**
  DGEEVX computes the eigenvalues and, optionally, the left and/or right
 eigenvectors for GE matrices

 Download DGEEVX + dependencies [TGZ] [ZIP] [TXT]

 Purpose:

      DGEEVX computes for an N-by-N real nonsymmetric matrix A, the
      eigenvalues and, optionally, the left and/or right eigenvectors.

      Optionally also, it computes a balancing transformation to improve
      the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
      SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
      (RCONDE), and reciprocal condition numbers for the right
      eigenvectors (RCONDV).

      The right eigenvector v(j) of A satisfies
                       A * v(j) = lambda(j) * v(j)
      where lambda(j) is its eigenvalue.
      The left eigenvector u(j) of A satisfies
                    u(j)**H * A = lambda(j) * u(j)**H
      where u(j)**H denotes the conjugate-transpose of u(j).

      The computed eigenvectors are normalized to have Euclidean norm
      equal to 1 and largest component real.

      Balancing a matrix means permuting the rows and columns to make it
      more nearly upper triangular, and applying a diagonal similarity
      transformation D * A * D**(-1), where D is a diagonal matrix, to
      make its rows and columns closer in norm and the condition numbers
      of its eigenvalues and eigenvectors smaller.  The computed
      reciprocal condition numbers correspond to the balanced matrix.
      Permuting rows and columns will not change the condition numbers
      (in exact arithmetic) but diagonal scaling will.  For further
      explanation of balancing, see section 4.10.2 of the LAPACK
      Users' Guide.
 **/

  char BALANC = does_diagonal_scaling ? (does_permutations ? 'B' : 'S')
                                      : (does_permutations ? 'P' : 'N');
  /**
  Parameters
      [in]	BALANC

                BALANC is CHARACTER*1
                Indicates how the input matrix should be diagonally scaled
                and/or permuted to improve the conditioning of its
                eigenvalues.
                = 'N': Do not diagonally scale or permute;
                = 'P': Perform permutations to make the matrix more nearly
                       upper triangular. Do not diagonally scale;
                = 'S': Diagonally scale the matrix, i.e. replace A by
                       D*A*D**(-1), where D is a diagonal matrix chosen
                       to make the rows and columns of A more equal in
                       norm. Do not permute;
                = 'B': Both diagonally scale and permute A.

                Computed reciprocal condition numbers will be for the matrix
                after balancing and/or permuting. Permuting does not change
                condition numbers (in exact arithmetic), but balancing does.
  **/

  bool calculates_right_eigen = true;
  char JOBVL = calculates_right_eigen ? 'V' : 'N';
  /**
      [in]	JOBVL

                JOBVL is CHARACTER*1
                = 'N': left eigenvectors of A are not computed;
                = 'V': left eigenvectors of A are computed.
                If SENSE = 'E' or 'B', JOBVL must = 'V'.
  **/
  bool calculates_left_eigen = true;
  char JOBVR = calculates_left_eigen ? 'V' : 'N';

  /**

  [in]	JOBVR

            JOBVR is CHARACTER*1
            = 'N': right eigenvectors of A are not computed;
            = 'V': right eigenvectors of A are computed.
            If SENSE = 'E' or 'B', JOBVR must = 'V'.

**/

  char SENSE = computes_eigenvalues_condition_numbers
                   ? (computes_eigenvectors_condition_numbers ? 'B' : 'E')
                   : (computes_eigenvectors_condition_numbers ? 'V' : 'N');

  /**
  [in]	SENSE

            SENSE is CHARACTER*1
            Determines which reciprocal condition numbers are computed.
            = 'N': None are computed;
            = 'E': Computed for eigenvalues only;
            = 'V': Computed for right eigenvectors only;
            = 'B': Computed for eigenvalues and right eigenvectors.

            If SENSE = 'E' or 'B', both left and right eigenvectors
            must also be computed (JOBVL = 'V' and JOBVR = 'V').
**/

  int N = x.nrows();
  /**
  [in]	N

            N is INTEGER
            The order of the matrix A. N >= 0.
**/

  Matrix<double> A(x);

  /**
[in]	LDA

          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).

*/
  int LDA = N;

  /**
[out]	WR

          WR is DOUBLE PRECISION array, dimension (N)

*/
  DiagonalMatrix<double> WR(x.nrows(), x.nrows(), false);

  /**
[out]	WI

          WI is DOUBLE PRECISION array, dimension (N)
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues.  Complex
          conjugate pairs of eigenvalues will appear consecutively
          with the eigenvalue having the positive imaginary part
          first.

*/
  DiagonalMatrix<double> WI(x.nrows(), x.nrows(), false);

  /**
[out]	VL

          VL is DOUBLE PRECISION array, dimension (LDVL,N)
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues.
          If JOBVL = 'N', VL is not referenced.
          If the j-th eigenvalue is real, then u(j) = VL(:,j),
          the j-th column of VL.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
          u(j+1) = VL(:,j) - i*VL(:,j+1).

*/
  Matrix<double> /* precision; dimension( ldvl; * ) */ VL_lapack(x.nrows(),
                                                                 x.nrows());

  /**
[in]	LDVL

          LDVL is INTEGER
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N.

*/
  int LDVL = N;

  /**
[out]	VR

          VR is DOUBLE PRECISION array, dimension (LDVR,N)
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues.
          If JOBVR = 'N', VR is not referenced.
          If the j-th eigenvalue is real, then v(j) = VR(:,j),
          the j-th column of VR.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
          v(j+1) = VR(:,j) - i*VR(:,j+1).

*/
  Matrix<double> /* precision; dimension( ldvl; * ) */ VR_lapack(x.nrows(),
                                                                 x.nrows());

  /**
[in]	LDVR

          LDVR is INTEGER
          The leading dimension of the array VR.  LDVR >= 1, and if
          JOBVR = 'V', LDVR >= N.
*/
  int LDVR = N;

  /**

[out]	ILO

          ILO is INTEGER

*/
  int ILO;

  /**
[out]	IHI

          IHI is INTEGER
          ILO and IHI are integer values determined when A was
          balanced.  The balanced A(i,j) = 0 if I > J and
          J = 1,...,ILO-1 or I = IHI+1,...,N.

*/
  int IHI;

  /**
[out]	SCALE

          SCALE is DOUBLE PRECISION array, dimension (N)
          Details of the permutations and scaling factors applied
          when balancing A.  If P(j) is the index of the row and column
          interchanged with row and column j, and D(j) is the scaling
          factor applied to row and column j, then
          SCALE(J) = P(J),    for J = 1,...,ILO-1
                   = D(J),    for J = ILO,...,IHI
                   = P(J)     for J = IHI+1,...,N.
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1.
*/
  DiagonalMatrix<double> /* precision; dimension( * ) */ SCALE(x.nrows(),
                                                               x.nrows());

  /**

[out]	ABNRM

          ABNRM is DOUBLE PRECISION
          The one-norm of the balanced matrix (the maximum
          of the sum of absolute values of elements of any column).

*/
  double ABNRM;

  /**
[out]	RCONDE

          RCONDE is DOUBLE PRECISION array, dimension (N)
          RCONDE(j) is the reciprocal condition number of the j-th
          eigenvalue.

*/
  DiagonalMatrix<double> /* precision; dimension( * ) */ RCONDE(x.nrows(),
                                                                x.nrows());

  /**
[out]	RCONDV

          RCONDV is DOUBLE PRECISION array, dimension (N)
          RCONDV(j) is the reciprocal condition number of the j-th
          right eigenvector.
*/
  DiagonalMatrix<double> /* precision; dimension( * ) */ RCONDV(x.nrows(),
                                                                x.nrows());

  /**

[out]	WORK

          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

*/
  Matrix<double> /* precision; dimension( * ) */ WORK(1ul, 1ul, 0.0);

  /**
[in]	LWORK

          LWORK is INTEGER
          The dimension of the array WORK.   If SENSE = 'N' or 'E',
          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',
          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).
          For good performance, LWORK must generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.

*/
  int LWORK = -1;

  /**
[out]	IWORK

          IWORK is INTEGER array, dimension (2*N-2)
          If SENSE = 'N' or 'E', not referenced.

*/
  Matrix<int> /*dimension( * ) */ IWORK(2, x.nrows() - 2);

  /**
[out]	INFO

          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors or condition numbers
                have been computed; elements 1:ILO-1 and i+1:N of WR
                and WI contain eigenvalues which have converged.

     */

  int INFO;
  dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0), &LDA, &WR[0], &WI[0],
          &VL_lapack[0], &LDVL, &VR_lapack[0], &LDVR, &ILO, &IHI, &SCALE[0],
          &ABNRM, &RCONDE[0], &RCONDV[0], &WORK[0], &LWORK, &IWORK[0], &INFO);

  Matrix<double> WORK_OPT(1, WORK[0]);
  LWORK = WORK[0];
  dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0), &LDA, &WR[0], &WI[0],
          &VL_lapack[0], &LDVL, &VR_lapack[0], &LDVR, &ILO, &IHI, &SCALE[0],
          &ABNRM, &RCONDE[0], &RCONDV[0], &WORK_OPT[0], &LWORK, &IWORK[0],
          &INFO);

  if (INFO != 0) {
    if (INFO > 0)
      return Error(
          std::string("the QR algorithm failed to compute all the  "
                      "eigenvalues, and no eigenvectors or condition numbers "
                      "have been computed; elements 1:ILO-1 and i+1:N of WR "
                      "and WI contain eigenvalues which have converged ILO=") +
          std::to_string(ILO) + " i=" + std::to_string(INFO));
    else {
      auto args =
          std::make_tuple(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0), &LDA,
                          &WR[0], &WI[0], &VL_lapack[0], &LDVL, &VR_lapack[0],
                          &LDVR, &ILO, &IHI, &SCALE[0], &ABNRM, &RCONDE[0],
                          &RCONDV[0], &WORK_OPT[0], &LWORK, &IWORK[0], &INFO);
      std::string argumentsNames[] = {
          "BALANC", "JOBVL", "JOBVR", "SENSE", "N",     "A",
          "LDA",    "WR",    "WI",    "VL",    "LDVL",  "VR",
          "LDVR",   "ILO",   "IHI",   "SCALE", "ABNRM", "RCONDE",
          "RCONDV", "WORK",  "LWORK", "IWORK", "INFO"};

      std::stringstream ss;
      std::size_t i = -INFO;
      write_tuple_i(ss, args, i);
      return Error("the" + std::to_string(i) + "-th argument " +
                   argumentsNames[i] + " had the illegal value =" + ss.str());
    }

  } else {
    auto VL_cpp = tr(VR_lapack);
    auto VR_cpp = tr(VL_lapack);
  // Nelson_Normalization(VR_cpp, VL_cpp);
    
    
    
//    std::cerr << "\ntr(VL_cpp)*VR_cpp\n" << tr(VL_cpp)*VR_cpp ;
//    std::cerr << "\nVL_cpp*VR_cpp\n" << VL_cpp*VR_cpp  ;
    
//    std::cerr << "\nx\n" << x;
//    std::cerr << "\nVR_cpp\n" << VR_cpp;
//    std::cerr << "\nWR\n" << WR;
//    std::cerr << "\nVL_cpp\n" << VL_cpp;
    
//    std::cerr << "\nVR_cpp*VL_cpp\n" << VR_cpp * VL_cpp;
//    std::cerr << "\nVR_cpp * WR * VL_cpp-x\n" << VR_cpp * WR * VL_cpp - x;
//    std::cerr << "\nVR_cpp * WR * inv(VR_cpp).value()\n" << VR_cpp * WR * inv(VR_cpp).value();
//    std::cerr << "\nVR_cpp * WR * inv(VR_cpp).value()-x\n" << VR_cpp * WR * inv(VR_cpp).value()-x;
//    std::cerr << "\ninv(tr(VR_cpp)).value() * WR * tr(VR_cpp)\n" << inv(tr(VR_cpp)).value() * WR * tr(VR_cpp);
//    std::cerr << "\ninv(tr(VR_cpp)).value() * WR * tr(VR_cpp)-x \n" << inv(tr(VR_cpp)).value() * WR * tr(VR_cpp);
    
    
//    std::cerr << "\ninv(VL_cpp).value() * WR * VL_cpp\n" << inv(VL_cpp).value() * WR * VL_cpp;
//    std::cerr << "\nVL_cpp * WR * inv(VL_cpp).value()\n" << VL_cpp * WR * inv(VL_cpp).value();
    
//    std::cerr << "\ninv(VL_cpp).value() * WR * VL_cpp\n" << inv(tr(VL_cpp)).value() * WR * tr(VL_cpp);
//    std::cerr << "\ntr(VL_cpp) * WR * inv(tr(VL_cpp)).value()\n" << tr(VL_cpp) * WR * inv(tr(VL_cpp)).value();
    
//    std::cerr << "\nnorm_1(VR_cpp * WR * VL_cpp-x)\n"
//              << norm_1(VR_cpp * WR * VL_cpp - x);
//    std::cerr << "\neps\n" << eps << "\n";
//    std::cerr << "\nstd::sqrt(eps)\n" << std::sqrt(eps) << "\n";
    
//    std::abort();
    
    auto inv_VR=inv(VR_cpp);
    if (!inv_VR)
      return error_message("eigenvector singular: "+inv_VR.error()());
    
    if (norm_1(VR_cpp * WR * inv(VR_cpp).value() - x)/norm_1(x) >  std::sqrt(eps)) {
      // std::cerr << "\nx\n" << x;
      // std::cerr << "\nVR_cpp\n" << VR_cpp;
      // std::cerr << "\nWR\n" << WR;
      // std::cerr << "\ninv(VR_cpp)\n" << inv(VR_cpp);
      // std::cerr << "\nVL_cpp\n" << VL_cpp;
      
      // std::cerr << "\nVR_cpp * WR * inv(VR_cpp).value()\n" << VR_cpp * WR * inv(VR_cpp).value() ;
      // std::cerr << "\nVR_cpp * WR * inv(VR_cpp).value() -x\n" << VR_cpp * WR * inv(VR_cpp).value()  - x;
      // std::cerr << "\neps\n" << eps << "\n";
      // std::cerr << "\nstd::sqrt(eps)\n" << std::sqrt(eps) << "\n";
      return error_message("bad eigenvalue decomposition");
    }
    
    assert(
        (norm_1(VR_cpp * WR * inv(VR_cpp).value() - x)/norm_1(x) <  std::sqrt(eps)*1000) &&
           " fails in eigendecomposition");
//    return sort_by_eigenvalue(std::make_tuple(
    return std::make_tuple(
        VR_cpp, WR,
        VL_cpp); // in reality VL, L, VR because of transposition
  }
}


extern "C" void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, int *N,
                        double * /*precision, dimension(lda, *) */ A, int *LDA,
                        double *VL, double *VU, int *IL, int *IU,
                        double *ABSTOL, int *M,
                        double * /*precision, dimension(*) */ W,
                        double * /*
                         dimension(ldz, *) */
                            Z,
                        int *LDZ, double * /*precision, dimension(*) */ WORK,
                        
                        
                        int *LWORK, int * /*, dimension(*) */ IWORK,
                        int * /*, dimension(*) */ IFAIL, int *INFO);



extern "C" void ddisna_(char *JOB, int * M, int *N,double * D, double *SEP,int * INFO );


inline  Maybe_error<std::tuple<Matrix<double>,DiagonalMatrix<double>,Matrix<double>>> Lapack_Symm_EigenSystem(const SymmetricMatrix<double> &x, std::string kind) {
    
    using lapack::dsyevx_;
    
    /** DSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices

Download DSYEVX + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DSYEVX computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of indices
     for the desired eigenvalues.
*/
    /**
Parameters
    [in]	JOBZ

              JOBZ is CHARACTER*1
              = 'N':  Compute eigenvalues only;
              = 'V':  Compute eigenvalues and eigenvectors.
*/
    
    char JOBZ = 'V';
    
    /**
    [in]	RANGE

              RANGE is CHARACTER*1
              = 'A': all eigenvalues will be found.
              = 'V': all eigenvalues in the half-open interval (VL,VU]
                     will be found.
              = 'I': the IL-th through IU-th eigenvalues will be found.


*/
    
    char RANGE = 'A';
    
    /**
    [in]	UPLO

              UPLO is CHARACTER*1
              = 'U':  Upper triangle of A is stored;
              = 'L':  Lower triangle of A is stored.


*/
    char UPLO;
    if (kind == "lower")
        UPLO = 'U';
    else
        UPLO = 'L';
    
    /*
    [in]	N

              N is INTEGER
              The order of the matrix A.  N >= 0.

*/
    
    int N = x.nrows();
    
    /**



    [in,out]	A

              A is DOUBLE PRECISION array, dimension (LDA, N)
              On entry, the symmetric matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of A contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A.
              On exit, the lower triangle (if UPLO='L') or the upper
              triangle (if UPLO='U') of A, including the diagonal, is
              destroyed.

*/
    double &A = const_cast<double &>(x[0]);
    
    /*   M_Matrix<double> A(x.nrows(), x.ncols());
  if (kind != "lower") {
    for (std::size_t i = 0; i < x.nrows(); ++i)
      for (std::size_t j = i; j < x.ncols(); ++j)
        A(i, j) = x(i, j);
  } else {
    for (std::size_t i = 0; i < x.nrows(); ++i)
      for (std::size_t j = 0; j < i + 1; ++j)
        A(i, j) = x(i, j);
  }
*/
    
    
    /**

    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

*/
    
    int LDA = N;
    
    
    /**


    [in]	VL

              VL is DOUBLE PRECISION
              If RANGE='V', the lower bound of the interval to
              be searched for eigenvalues. VL < VU.
              Not referenced if RANGE = 'A' or 'I'.
*/
    
    double VL;
    
    
    /**


    [in]	VU

              VU is DOUBLE PRECISION
              If RANGE='V', the upper bound of the interval to
              be searched for eigenvalues. VL < VU.
              Not referenced if RANGE = 'A' or 'I'.
*/
    
    double VU;
    
    /**




    [in]	IL

              IL is INTEGER
              If RANGE='I', the index of the
              smallest eigenvalue to be returned.
              1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
              Not referenced if RANGE = 'A' or 'V'.


*/
    
    int IL;
    
    
    /**
    [in]	IU

              IU is INTEGER
              If RANGE='I', the index of the
              largest eigenvalue to be returned.
              1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
              Not referenced if RANGE = 'A' or 'V'.
*/
    int IU;
    
    
    /**

    [in]	ABSTOL

              ABSTOL is DOUBLE PRECISION
              The absolute error tolerance for the eigenvalues.
              An approximate eigenvalue is accepted as converged
              when it is determined to lie in an interval [a,b]
              of width less than or equal to

                      ABSTOL + EPS *   max( |a|,|b| ) ,

              where EPS is the machine precision.  If ABSTOL is less than
              or equal to zero, then  EPS*|T|  will be used in its place,
              where |T| is the 1-norm of the tridiagonal matrix obtained
              by reducing A to tridiagonal form.

              Eigenvalues will be computed most accurately when ABSTOL is
              set to twice the underflow threshold 2*DLAMCH('S'), not zero.
              If this routine returns with INFO>0, indicating that some
              eigenvectors did not converge, try setting ABSTOL to
              2*DLAMCH('S').

              See "Computing Small Singular Values of Bidiagonal Matrices
              with Guaranteed High Relative Accuracy," by Demmel and
              Kahan, LAPACK Working Note #3.
*/
    
    double ABSTOL = std ::numeric_limits<double>::epsilon()*100;
    
    /**
    [out]	M

              M is INTEGER
              The total number of eigenvalues found.  0 <= M <= N.
              If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*/
    
    int M;
    
    /**
    [out]	W

              W is DOUBLE PRECISION array, dimension (N)
              On normal exit, the first M elements contain the selected
              eigenvalues in ascending order.
*/
    
    DiagonalMatrix<double> W(x.nrows(), x.ncols());
    
    
    /**
    [out]	Z

              Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
              If JOBZ = 'V', then if INFO = 0, the first M columns of Z
              contain the orthonormal eigenvectors of the matrix A
              corresponding to the selected eigenvalues, with the i-th
              column of Z holding the eigenvector associated with W(i).
              If an eigenvector fails to converge, then that column of Z
              contains the latest approximation to the eigenvector, and the
              index of the eigenvector is returned in IFAIL.
              If JOBZ = 'N', then Z is not referenced.
              Note: the user must ensure that at least max(1,M) columns are
              supplied in the array Z; if RANGE = 'V', the exact value of M
              is not known in advance and an upper bound must be used.



    [in]	LDZ

              LDZ is INTEGER
              The leading dimension of the array Z.  LDZ >= 1, and if
              JOBZ = 'V', LDZ >= max(1,N).



*/
    
    int LDZ = N;
    Matrix<double> Z(x.nrows(), x.nrows());
    /*

    [out]	WORK

              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

*/
    
    std ::vector<double> WORK(1);
    
    /**
    [in]	LWORK

              LWORK is INTEGER
              The length of the array WORK.  LWORK >= 1, when N <= 1;
              otherwise 8*N.
              For optimal efficiency, LWORK >= (NB+3)*N,
              where NB is the max of the blocksize for DSYTRD and DORMTR
              returned by ILAENV.

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA.
*/
    
    int LWORK = -1;
    
    /**
    [out]	IWORK

              IWORK is INTEGER array, dimension (5*N)

*/
    
    std ::vector<int> IWORK(5 * x.nrows());
    
    
    /**
    [out]	IFAIL

              IFAIL is INTEGER array, dimension (N)
              If JOBZ = 'V', then if INFO = 0, the first M elements of
              IFAIL are zero.  If INFO > 0, then IFAIL contains the
              indices of the eigenvectors that failed to converge.
              If JOBZ = 'N', then IFAIL is not referenced.
*/
    std::vector<int> IFAIL(x.nrows());
    
    /**
    [out]	INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, then i eigenvectors failed to converge.
                    Their indices are stored in array IFAIL.
*/
    
    int INFO;
    
    
    dsyevx_(&JOBZ, &RANGE, &UPLO, &N, &A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0],&INFO );
    
    LWORK = WORK[0];
    WORK.resize(LWORK);
    dsyevx_(&JOBZ, &RANGE, &UPLO, &N, &A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0], &INFO);
    
    
    
    
    if (INFO != 0) {
        if (INFO > 0)
        {
            std::stringstream ss;
            ss<<INFO<<" eigenvectors failed to converge. Their indices are: ";
            std::for_each(IFAIL.begin(),IFAIL.end(),[&ss](auto const& x){ss<<x;});
            return error_message(ss.str());
        }
        else {
            auto args = std::make_tuple(&JOBZ, &RANGE, &UPLO, &N, &A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W, &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0],&INFO);
            std::string argumentsNames[] = {
                                            "JOBZ", "RANGE", "UPLO", "N", "A[0]", "LDA", "VL", "VU", "IL", "IU", "ABSTOL", "M", "W", "Z[0]", "LDZ", "WORK[0]", "LWORK",  "IWORK[0]","IFAIL[0]", "INFO"    };
            
            std::stringstream ss;
            std::size_t i = -INFO;
            write_tuple_i(ss, args, i);
            return error_message("the" + std::to_string(i) + "-th argument " +
                                 argumentsNames[i] +
                                 " had the illegal value =" + ss.str());
        }
        
    } else {
        /**
 * extern "C" void ddisna_(char *JOB, int * M, int *N,double * D, double *SEP,int * INFO );
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), SEP( * )
*     ..
*
*  Purpose
*  =======
*
*  DDISNA computes the reciprocal condition numbers for the eigenvectors
*  of a real symmetric or complex Hermitian matrix or for the left or
*  right singular vectors of a general m-by-n matrix. The reciprocal
*  condition number is the 'gap' between the corresponding eigenvalue or
*  singular value and the nearest other one.
*
*  The bound on the error, measured by angle in radians, in the I-th
*  computed vector is given by
*
*         DLAMCH( 'E' ) * ( ANORM / SEP( I ) )
*
*  where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed
*  to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of
*  the error bound.
*
*  DDISNA may also be used to compute error bounds for eigenvectors of
*  the generalized symmetric definite eigenproblem.
*/
        
        /**
*  JOB     (input) CHARACTER*1
*          Specifies for which problem the reciprocal condition numbers
*          should be computed:
*          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
*          = 'L':  the left singular vectors of a general matrix;
*          = 'R':  the right singular vectors of a general matrix.
*/
        char JOB='E';
        /**  M       (input) INTEGER
*          The number of rows of the matrix. M >= 0.
*
*/
        
        /**
*  N       (input) INTEGER
*          If JOB = 'L' or 'R', the number of columns of the matrix,
*          in which case N >= 0. Ignored if JOB = 'E'.
*/
        
        /**  D       (input) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
*                              dimension (min(M,N)) if JOB = 'L' or 'R'
*          The eigenvalues (if JOB = 'E') or singular values (if JOB =
*          'L' or 'R') of the matrix, in either increasing or decreasing
*          order. If singular values, they must be non-negative.
*
*/
        
        /**  SEP     (output) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
*                               dimension (min(M,N)) if JOB = 'L' or 'R'
*          The reciprocal condition numbers of the vectors.
*/
        DiagonalMatrix<double> RCOND(x.nrows(), x.ncols());
        
        /**  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*/
        lapack::ddisna_(&JOB,&M,&N,&W[0],&RCOND[0] ,&INFO );
        auto VL_cpp = Z;
        auto VR_cpp = tr(Z);
        if (!Nelson_Normalization(VR_cpp, VL_cpp))
            return error_message("singular normalization");
        
        return std::tuple(
            VR_cpp, W,
            VL_cpp);
    }
}





template <class T>
Maybe_error<SymPosDefMatrix<T>>
Lapack_UT_Cholesky_inv(const UpTrianMatrix<T> &x) {
  return_error<SymPosDefMatrix<double>, Lapack_SymmPosDef_inv> Error;

  /*
DPOTRI

Download DPOTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DPOTRI computes the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF.

Parameters
*/
  char UPLO = 'L';
  /*
 * [in]	UPLO

            UPLO is CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.
*/

  int N = x.nrows();

  /*
 *
  [in]	N

            N is INTEGER
            The order of the matrix A.  N >= 0.
*/
  auto a = x;
  double &A = const_cast<double &>(a[0]);


  /*
  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            DPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.
*/
  int LDA = N;

  /*
    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;

  /*
    [out]	INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, the (i,i) element of the factor U or L is
                    zero, and the inverse could not be computed.*
*/
  dpotri_(&UPLO, &N, &A, &LDA, &INFO);

  if (INFO < 0) {
    return Error("the " + std::to_string(INFO) +
                 "th argument for inversehad an illegal value");
  } else if (INFO > 0) {
    return Error(" the (" + std::to_string(INFO) + "," + std::to_string(INFO) +
                 ") element of the factor U or L is "
                 "zero, and the inverse could not be computed");
  } else {
    copy_UT_to_LT(a);
    return SymPosDefMatrix<double>(std::move(a));
  }
}

template <class T>
Maybe_error<SymPosDefMatrix<T>>
Lapack_LT_Cholesky_inv(const DownTrianMatrix<T> &x) {
  return_error<SymPosDefMatrix<double>> Error{"Lapack_LT_Cholesky_inv :"};

  /*
DPOTRI

Download DPOTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DPOTRI computes the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF.

Parameters
*/
  char UPLO = 'L';
  /*
 * [in]	UPLO

            UPLO is CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.
*/

  int N = x.nrows();

  /*
 *
  [in]	N

            N is INTEGER
            The order of the matrix A.  N >= 0.
*/
  auto a = x;
  double &A = const_cast<double &>(a[0]);


  /*
  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            DPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.
*/
  int LDA = N;

  /*
    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;

  /*
    [out]	INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, the (i,i) element of the factor U or L is
                    zero, and the inverse could not be computed.*
*/
  dpotri_(&UPLO, &N, &A, &LDA, &INFO);

  if (INFO < 0) {
    return Error("the " + std::to_string(INFO) +
                 "th argument for inversehad an illegal value");
  } else if (INFO > 0) {
    return Error(" the (" + std::to_string(INFO) + "," + std::to_string(INFO) +
                 ") element of the factor U or L is "
                 "zero, and the inverse could not be computed");
  } else {
    copy_UT_to_LT(a);
    return SymPosDefMatrix<double>(std::move(a));
  }
}

extern "C" void dsyrk_(char *UPLO, char *TRANS, int *N, int *K, double *ALPHA,
                       double *A, int *LDA, double *BETA, double *C, int *LDC);

inline auto &Lapack_Product_Self_Transpose_mod(const Matrix<double> &a,
                                        SymPosDefMatrix<double> &c,
                                        bool first_transposed_in_c,
                                        char UPLO_in_c = 'U', double alpha = 1,
                                        double beta = 0) {

  char UPLO = (UPLO_in_c == 'U') ? 'L' : 'U';
  /*
    [in]	UPLO

              UPLO is CHARACTER*1
               On  entry,   UPLO  specifies  whether  the  upper  or  lower
               triangular  part  of the  array  C  is to be  referenced  as
               follows:

                  UPLO = 'U' or 'u'   Only the  upper triangular part of  C
                                      is to be referenced.

                  UPLO = 'L' or 'l'   Only the  lower triangular part of  C
                                      is to be referenced.
*/
  char TRANS = first_transposed_in_c ? 'N' : 'T';
  /*
  [in]	TRANS

            TRANS is CHARACTER*1
             On entry,  TRANS  specifies the operation to be performed as
             follows:

                TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.

                TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.

                TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.

*/

  int N = c.nrows();
  /*      [in]	N

              N is INTEGER
                  On entry,  N specifies the order of the matrix C.  N must be
          at least zero.
  */
  int K = first_transposed_in_c ? a.nrows() : a.ncols();

  /*
                [in]	K

                        K is INTEGER
                            On entry with  TRANS = 'N' or 'n',  K  specifies
   the number of  columns   of  the   matrix   A,   and  on   entry   with
        TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
        of rows of the matrix  A.  K must be at least zero.

*/
  double ALPHA = alpha;

  /*
*     [in]	ALPHA

            ALPHA is DOUBLE PRECISION.
             On entry, ALPHA specifies the scalar alpha.

*/
  double &A = const_cast<double &>(a[0]);

  /*

  [in]	A

            A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
             part of the array  A  must contain the matrix  A,  otherwise
             the leading  k by n  part of the array  A  must contain  the
             matrix A.

* */

  int LDA = first_transposed_in_c ? N : K;

  /*
 *     [in]	LDA

            LDA is INTEGER
             On entry, LDA specifies the first dimension of A as declared
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
             then  LDA must be at least  max( 1, n ), otherwise  LDA must
             be at least  max( 1, k ).

 * */

  double BETA = beta;
  /*
 *     [in]	BETA

            BETA is DOUBLE PRECISION.
             On entry, BETA specifies the scalar beta.

*/

  double &C = c[0];

  /*
  [in,out]	C

            C is DOUBLE PRECISION array, dimension ( LDC, N )
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
             upper triangular part of the array C must contain the upper
             triangular part  of the  symmetric matrix  and the strictly
             lower triangular part of C is not referenced.  On exit, the
             upper triangular part of the array  C is overwritten by the
             upper triangular part of the updated matrix.
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
             lower triangular part of the array C must contain the lower
             triangular part  of the  symmetric matrix  and the strictly
             upper triangular part of C is not referenced.  On exit, the
             lower triangular part of the array  C is overwritten by the
             lower triangular part of the updated matrix.
*/

  int LDC = N;

  /*
  [in]	LDC

            LDC is INTEGER
             On entry, LDC specifies the first dimension of C as declared
             in  the  calling  (sub)  program.   LDC  must  be  at  least
             max( 1, n ).


 * */

  dsyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, &A, &LDA, &BETA, &C, &LDC);

  copy_UT_to_LT(c);
  return c;
}

/*
 * dsyrk 	( 	character  	UPLO,
                character  	TRANS,
                integer  	N,
                integer  	K,
                double precision  	ALPHA,
                double precision, dimension(lda,*)  	A,
                integer  	LDA,
                double precision  	BETA,
                double precision, dimension(ldc,*)  	C,
                integer  	LDC
        )

DSYRK

Purpose:

     DSYRK  performs one of the symmetric rank k operations

        C := alpha*A*A**T + beta*C,

     or

        C := alpha*A**T*A + beta*C,

     where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     in the second case.

Parameters




 * */

inline SymPosDefMatrix<double>
Lapack_Product_Self_Transpose(const Matrix<double> &a,
                              bool first_transposed_in_c, char UPLO_in_c,
                              double alpha, double beta) {
  std::size_t n = first_transposed_in_c ? a.ncols() : a.nrows();
  SymPosDefMatrix<double> c(n, n,false);
  c = Lapack_Product_Self_Transpose_mod(a, c, first_transposed_in_c, UPLO_in_c,
                                        alpha, beta);
  return c;
};

inline Maybe_error<DownTrianMatrix<double>>
Lapack_chol(const SymPosDefMatrix<double> &x) {
  return_error<DownTrianMatrix<double>, Lapack_chol> Error;
  assert(x.nrows() == x.ncols());
  auto a = x;
  /*
 *      DPOTRF computes the Cholesky factorization of a real symmetric
   positive definite matrix A.

   The factorization has the form
      A = U**T * U,  if UPLO = 'U', or
      A = L  * L**T,  if UPLO = 'L',
   where U is an upper triangular matrix and L is lower triangular.

   This is the block version of the algorithm, calling Level 3 BLAS.



Parameters
*/
  char UPLO = 'U';
  /*
 *
 * [in]	UPLO

            UPLO is CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.
*/
  int N = x.nrows();

  /*
  [in]	N

            N is INTEGER
            The order of the matrix A.  N >= 0.
*/

  double &A = a[0];

  /*



  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.

            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T.
*/
  int LDA = N;

  /*
  [in]	LDA

            LDA is INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;
  /*
  [out]	INFO

            INFO is INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

*/

  if (x.size() == 0)
    return Error(" ZERO MATRIX");
  try {
    lapack::dpotrf_(&UPLO, &N, &A, &LDA, &INFO);
  } catch (...) {
    assert(false);
  }

  if (INFO != 0) {
    if (INFO < 0)
      return Error("Cholesky fails, the" + std::to_string(-INFO) +
                   "-th argument had an illegal value");
    else
      return Error("Cholesky fails, zero diagonal at" + std::to_string(INFO));
  } else {

    return fill_UT_zeros(std::move(a));
  }
}

extern "C" void dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int *M,
                       int *N, double *ALPHA, double *A, int *LDA, double *B,
                       int *LDB);

inline Matrix<double>
Lapack_Triang_Product(const Matrix<double> &a, const Matrix<double> &b,
                      bool up_triangular_in_c, bool triangular_first_in_c,
                      bool transpose_A_in_c, bool ones_in_diag, double alpha) {
  /*
dtrmm()
subroutine dtrmm 	( 	character  	SIDE,
              character  	UPLO,
              character  	TRANSA,
              character  	DIAG,
              integer  	M,
              integer  	N,
              double precision  	ALPHA,
              double precision, dimension(lda,*)  	A,
              integer  	LDA,
              double precision, dimension(ldb,*)  	B,
              integer  	LDB
      )

DTRMM

Purpose:

   DTRMM  performs one of the matrix-matrix operations

      B := alpha*op( A )*B,   or   B := alpha*B*op( A ),

   where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
   non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

      op( A ) = A   or   op( A ) = A**T.

Parameters
*/
  char SIDE = triangular_first_in_c ? 'R' : 'L';

  /*
  [in]	SIDE

            SIDE is CHARACTER*1
             On entry,  SIDE specifies whether  op( A ) multiplies B from
             the left or right as follows:

                SIDE = 'L' or 'l'   B := alpha*op( A )*B.

                SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*/
  char UPLO = up_triangular_in_c ? 'L' : 'U';
  /*
  [in]	UPLO

            UPLO is CHARACTER*1
             On entry, UPLO specifies whether the matrix A is an upper or
             lower triangular matrix as follows:

                UPLO = 'U' or 'u'   A is an upper triangular matrix.

                UPLO = 'L' or 'l'   A is a lower triangular matrix.
*/
  char TRANSA = transpose_A_in_c ? 'T' : 'N';
  /*



  [in]	TRANSA

            TRANSA is CHARACTER*1
             On entry, TRANSA specifies the form of op( A ) to be used in
             the matrix multiplication as follows:

                TRANSA = 'N' or 'n'   op( A ) = A.

                TRANSA = 'T' or 't'   op( A ) = A**T.

                TRANSA = 'C' or 'c'   op( A ) = A**T.
*/
  char DIAG = ones_in_diag ? 'U' : 'N';
  /*
  [in]	DIAG

            DIAG is CHARACTER*1
             On entry, DIAG specifies whether or not A is unit triangular
             as follows:

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.

                DIAG = 'N' or 'n'   A is not assumed to be unit
                                    triangular.

*/
  int M = b.ncols();
  /*
  [in]	M

            M is INTEGER
             On entry, M specifies the number of rows of B. M must be at
             least zero.
*/
  int N = b.nrows();
  /*
  [in]	N

            N is INTEGER
             On entry, N specifies the number of columns of B.  N must be
             at least zero.
*/
  double ALPHA = alpha;
  /*
  [in]	ALPHA

            ALPHA is DOUBLE PRECISION.
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is
             zero then  A is not referenced and  B need not be set before
             entry.
*/
  double &A = const_cast<double &>(a[0]);


  /*
  [in]	A

             A is DOUBLE PRECISION array, dimension ( LDA, k ), where k is m
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
             upper triangular part of the array  A must contain the upper
             triangular matrix  and the strictly lower triangular part of
             A is not referenced.
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
             lower triangular part of the array  A must contain the lower
             triangular matrix  and the strictly upper triangular part of
             A is not referenced.
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of
             A  are not referenced either,  but are assumed to be  unity.
*/

  int LDA = SIDE == 'L' ? M : N;

  /*
  [in]	LDA

            LDA is INTEGER
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
             LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
             then LDA must be at least max( 1, n ).
*/

  auto out = b;

  double &B = out(0, 0);

  /*
  [in,out]	B

            B is DOUBLE PRECISION array, dimension ( LDB, N )
             Before entry,  the leading  m by n part of the array  B must
             contain the matrix  B,  and  on exit  is overwritten  by the
             transformed matrix.
*/
  int LDB = M;
  /*
  [in]	LDB

            LDB is INTEGER
             On entry, LDB specifies the first dimension of B as declared
             in  the  calling  (sub)  program.   LDB  must  be  at  least
             max( 1, m ). *
 * */

  try {
    dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, &A, &LDA, &B, &LDB);

  } catch (...) {
    assert(false);
  }
  return out;
}

extern "C" void dtrtri_(char *UPLO, char *DIAG, int *N, double *A, int *LDA,
                        int *INFO);

inline Maybe_error<DownTrianMatrix<double>>
Lapack_LT_inv(const DownTrianMatrix<double> &x, bool ones_in_diag) {

  return_error<DownTrianMatrix<double>, Lapack_LT_inv> Error;

  /*
subroutine dtrtri 	( 	character  	UPLO,
                character  	DIAG,
                integer  	N,
                double precision, dimension( lda, * )  	A,
                integer  	LDA,
                integer  	INFO
        )

DTRTRI

Download DTRTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DTRTRI computes the inverse of a real upper or lower triangular
     matrix A.

     This is the Level 3 BLAS version of the algorithm.

Parameters

*/
  char UPLO = 'U';
  /*
    [in]	UPLO

              UPLO is CHARACTER*1
              = 'U':  A is upper triangular;
              = 'L':  A is lower triangular.
*/

  char DIAG = ones_in_diag ? 'U' : 'N';
  /*
    [in]	DIAG

              DIAG is CHARACTER*1
              = 'N':  A is non-unit triangular;
              = 'U':  A is unit triangular.
*/
  int N = x.nrows();
  /*
    [in]	N

              N is INTEGER
              The order of the matrix A.  N >= 0.
*/
  auto a = x;
  double &A = const_cast<double &>(a[0]);


  /*
    [in,out]	A

              A is DOUBLE PRECISION array, dimension (LDA,N)
              On entry, the triangular matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of the array A contains
              the upper triangular matrix, and the strictly lower
              triangular part of A is not referenced.  If UPLO = 'L', the
              leading N-by-N lower triangular part of the array A contains
              the lower triangular matrix, and the strictly upper
              triangular part of A is not referenced.  If DIAG = 'U', the
              diagonal elements of A are also not referenced and are
              assumed to be 1.
              On exit, the (triangular) inverse of the original matrix, in
              the same storage format.
*/
  int LDA = N;
  /*
    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;
  /*

    [out]	INFO

              INFO is INTEGER
              = 0: successful exit
              < 0: if INFO = -i, the i-th argument had an illegal value
              > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                   matrix is singular and its inverse can not be computed.



*/

  dtrtri_(&UPLO, &DIAG, &N, &A, &LDA, &INFO);

  if (INFO < 0)
    return Error(std::to_string(INFO) + " argument had an illegal value");
  else if (INFO > 0)
    return Error("A(" + std::to_string(INFO) + "," + std::to_string(INFO) +
                 ") is exactly zero.  The triangular matrix is singular and "
                 "its inverse can not be computed");
  else
    return a;
}

inline Maybe_error<UpTrianMatrix<double>> Lapack_UT_inv(const UpTrianMatrix<double> &x,
                                                 bool ones_in_diag) {

  return_error<UpTrianMatrix<double>, Lapack_UT_inv> Error;

  /*
subroutine dtrtri 	( 	character  	UPLO,
                character  	DIAG,
                integer  	N,
                double precision, dimension( lda, * )  	A,
                integer  	LDA,
                integer  	INFO
        )

DTRTRI

Download DTRTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DTRTRI computes the inverse of a real upper or lower triangular
     matrix A.

     This is the Level 3 BLAS version of the algorithm.

Parameters

*/
  char UPLO = 'L';
  /*
    [in]	UPLO

              UPLO is CHARACTER*1
              = 'U':  A is upper triangular;
              = 'L':  A is lower triangular.
*/

  char DIAG = ones_in_diag ? 'U' : 'N';
  /*
    [in]	DIAG

              DIAG is CHARACTER*1
              = 'N':  A is non-unit triangular;
              = 'U':  A is unit triangular.
*/
  int N = x.nrows();
  /*
    [in]	N

              N is INTEGER
              The order of the matrix A.  N >= 0.
*/
  auto a = x;
  double &A = const_cast<double &>(a[0]);


  /*
    [in,out]	A

              A is DOUBLE PRECISION array, dimension (LDA,N)
              On entry, the triangular matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of the array A contains
              the upper triangular matrix, and the strictly lower
              triangular part of A is not referenced.  If UPLO = 'L', the
              leading N-by-N lower triangular part of the array A contains
              the lower triangular matrix, and the strictly upper
              triangular part of A is not referenced.  If DIAG = 'U', the
              diagonal elements of A are also not referenced and are
              assumed to be 1.
              On exit, the (triangular) inverse of the original matrix, in
              the same storage format.
*/
  int LDA = N;
  /*
    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).
*/
  int INFO;
  /*

    [out]	INFO

              INFO is INTEGER
              = 0: successful exit
              < 0: if INFO = -i, the i-th argument had an illegal value
              > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                   matrix is singular and its inverse can not be computed.



*/

  dtrtri_(&UPLO, &DIAG, &N, &A, &LDA, &INFO);

  if (INFO < 0)
    return Error(std::to_string(INFO) + " argument had an illegal value");
  else if (INFO > 0)
    return Error("A(" + std::to_string(INFO) + "," + std::to_string(INFO) +
                 ") is exactly zero.  The triangular matrix is singular and "
                 "its inverse can not be computed");
  else
    return a;
}




} // namespace lapack

#endif // LAPACK_HEADERS_H
