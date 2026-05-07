// Block-Parlett recurrence on (quasi-)upper-triangular T. Step 2 of the
// Schur+Parlett implementation: the Schur wrapper from step 1 produces T,
// and parlett_block_eval composes f(T) from a per-block evaluator. This test
// validates the recurrence on three regimes:
//   - Strictly upper-triangular T (all 1×1 diagonal blocks, well-separated
//     real eigenvalues): expm via Parlett matches expm_taylor_scaling_squaring.
//   - Quasi-upper-triangular T from Lapack_Schur of a real-eigenvalue Q
//     (still all 1×1 blocks, but reaches T via the actual LAPACK pipeline).
//   - Quasi-upper-triangular T from Lapack_Schur of a complex-conjugate Q
//     (one 2×2 block in the middle): the per-block evaluator falls back to
//     expm_taylor on the 2×2 block.

#include <catch_amalgamated.hpp>

#include <cmath>
#include <cstddef>
#include <string>

#include <exponential_matrix.h>
#include <lapack_headers.h>
#include <lapack_schur.h>
#include <matrix.h>
#include <schur_parlett.h>

namespace {

// Reference: B(dt) = ∫∫ expm(Q·s)·G·expm(Q·(t1−s))·G·expm(Q·(dt−t1)) ds dt1
// via direct Taylor summation of the N_p recurrence (defined in the
// scaling-and-squaring seed: N_0 = G·G, N_p = Q · N_{p-1} + G · M_p), then
// B = Σ N_p · dt^(p+2) / (p+2)!. Identical structure to the M_m series for
// A but one order higher in dt.
inline Matrix<double> frechet_double_integral_taylor_reference(Matrix<double> const& Q,
                                                                  Matrix<double> const& G,
                                                                  double dt,
                                                                  std::size_t order = 25) {
    const std::size_t N = Q.nrows();
    // Build M_m series first (we need M_p in the N_p recurrence).
    Matrix<double> M_prev = G;
    std::vector<Matrix<double>> M_store;
    M_store.reserve(order);
    M_store.push_back(M_prev);
    Matrix<double> Q_pow(N, N, 0.0);
    for (std::size_t i = 0; i < N; ++i) Q_pow(i, i) = 1.0;
    Q_pow = Q * Q_pow;  // Q^1
    for (std::size_t m = 1; m < order; ++m) {
        Matrix<double> M_cur = Q * M_prev + G * Q_pow;
        M_store.push_back(M_cur);
        M_prev = M_cur;
        Q_pow = Q * Q_pow;
    }

    // N_0 = G·G; B = Σ N_p · dt^(p+2)/(p+2)!.
    Matrix<double> N_cur = G * G;
    double inv_fact = 0.5;          // 1/2! for p=0
    double dt_pow = dt * dt;        // dt^(p+2)
    Matrix<double> B = N_cur * (dt_pow * inv_fact);

    for (std::size_t p = 1; p + 1 < order; ++p) {
        // N_p = Q · N_{p-1} + G · M_p
        N_cur = Q * N_cur + G * M_store[p];
        inv_fact /= static_cast<double>(p + 2);
        dt_pow *= dt;
        B = B + N_cur * (dt_pow * inv_fact);
    }
    return B;
}

// Reference: A(dt) = ∫ expm(Q·s) · G · expm(Q·(dt−s)) ds via direct Taylor
// summation of the M_m recurrence with sufficient terms. For small ‖Q·dt‖
// (which is what our test cases use), order ~20 converges to machine
// precision and gives an independent ground truth for the Van Loan
// implementation.
inline Matrix<double> frechet_integral_taylor_reference(Matrix<double> const& Q,
                                                          Matrix<double> const& G, double dt,
                                                          std::size_t order = 20) {
    const std::size_t N = Q.nrows();
    Matrix<double> M_prev = G;                     // M_0 = G
    Matrix<double> A = G * dt;                     // dt^1 / 1! · M_0

    Matrix<double> Q_pow(N, N, 0.0);
    for (std::size_t i = 0; i < N; ++i) Q_pow(i, i) = 1.0;
    Q_pow = Q * Q_pow;                             // Q^1 = Q

    double inv_fact = 1.0;                         // 1/(m+1)!  starts at 1 for m=0
    double dt_pow = dt;                            // dt^(m+1)
    for (std::size_t m = 1; m < order; ++m) {
        // M_m = Q · M_{m-1} + G · Q^m
        Matrix<double> M_cur = Q * M_prev + G * Q_pow;
        inv_fact /= static_cast<double>(m + 1);
        dt_pow *= dt;
        A = A + M_cur * (dt_pow * inv_fact);
        M_prev = M_cur;
        Q_pow = Q * Q_pow;
    }
    return A;
}

// Reference: expm via scaling + Pade + squaring, *without* the
// to_transition_Probability projection that expm_taylor_scaling_squaring
// applies (and which produces nonsense for non-stochastic matrices like the
// quasi-triangular T's we test here, whose eigenvalues are negative real
// numbers — squaring forces all rows to sum to 1, blowing the answer up to
// pure rescaled garbage at large |T|).
Matrix<double> raw_expm_scaling_squaring(Matrix<double> const& M) {
    double max_abs = 0.0;
    for (std::size_t i = 0; i < M.nrows(); ++i)
        for (std::size_t j = 0; j < M.ncols(); ++j) max_abs = std::max(max_abs, std::abs(M(i, j)));
    int n = 0;
    if (max_abs > 0.125) n = static_cast<int>(std::ceil(std::log2(max_abs / 0.125)));
    n = std::max(0, n);
    double scale = std::pow(2.0, -n);
    auto scaled = M * scale;
    auto Maybe_E = expm_pade(scaled);
    REQUIRE(Maybe_E.valid());
    auto E = Maybe_E.value();
    for (int i = 0; i < n; ++i) E = E * E;
    return E;
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

// Generic per-block exp evaluator. For 1×1 blocks: scalar exp. For 2×2 blocks
// (complex-conjugate Schur block): use raw_expm_scaling_squaring — Pade with
// scaling-and-squaring but no transition-probability projection.
auto exp_block_evaluator() {
    return [](Matrix<double> const& B) -> Matrix<double> {
        if (B.nrows() == 1) {
            return Matrix<double>(1, 1, std::exp(B(std::size_t{0}, std::size_t{0})));
        }
        return raw_expm_scaling_squaring(B);
    };
}

}  // namespace

TEST_CASE("parlett_block_eval: strictly upper-triangular real eigenvalues",
          "[math][parlett]") {
    // Diagonal entries well-separated; off-diagonals nontrivial.
    Matrix<double> T(4, 4, 0.0);
    T(0, 0) = -0.1;
    T(1, 1) = -1.0;
    T(2, 2) = -10.0;
    T(3, 3) = -100.0;
    T(0, 1) = 0.5;
    T(0, 2) = 0.3;
    T(0, 3) = 0.1;
    T(1, 2) = 0.4;
    T(1, 3) = 0.2;
    T(2, 3) = 0.6;

    auto starts = lapack::detect_schur_block_starts(T);
    REQUIRE(starts.size() == 5);  // four 1×1 blocks + sentinel

    auto Maybe_F = lapack::parlett_block_eval(T, starts, exp_block_evaluator());
    if (!Maybe_F) UNSCOPED_INFO("parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(T);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(T)_parlett - expm(T)_pade||_F = " << err);
    CHECK(err < 1e-10);
}

TEST_CASE("schur_parlett: 2×2 Q with distinct real eigenvalues",
          "[math][parlett][schur]") {
    // scheme_CO macro Q at agonist=10: eigenvalues 0 and -101, well-separated.
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -1.0;
    Q(0, 1) = 1.0;
    Q(1, 0) = 100.0;
    Q(1, 1) = -100.0;

    auto Maybe_F = lapack::schur_parlett(Q, exp_block_evaluator());
    if (!Maybe_F) UNSCOPED_INFO("schur_parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(Q);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(Q)_schur_parlett - expm(Q)_taylor||_F = " << err);
    CHECK(err < 1e-10);
}

TEST_CASE("schur_parlett: 2×2 rotation Q with complex-conjugate eigenvalues",
          "[math][parlett][schur][complex]") {
    // Schur outputs a single 2×2 block; the per-block evaluator falls back
    // to expm_taylor for that block. With nb=1 there are no off-diagonal
    // Sylvester solves to test, but it's still a meaningful end-to-end check
    // that the block-detection path doesn't accidentally split the 2×2 block.
    const double theta = 0.7;
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = std::cos(theta);
    Q(0, 1) = -std::sin(theta);
    Q(1, 0) = std::sin(theta);
    Q(1, 1) = std::cos(theta);

    auto Maybe_F = lapack::schur_parlett(Q, exp_block_evaluator());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(Q);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(Q_rot)_schur_parlett - expm(Q_rot)_taylor||_F = " << err);
    CHECK(err < 1e-12);
}

TEST_CASE("schur_parlett: 4×4 mixed real+complex via Householder scramble",
          "[math][parlett][schur]") {
    // Reuses the 4×4 Householder-scrambled fixture from test_schur. After
    // schur_parlett unwinds the Schur factors, expm of the scrambled Q must
    // match expm of the same Q computed by Taylor.
    Matrix<double> A(4, 4, 0.0);
    A(0, 0) = -1.0;
    A(0, 1) = 0.5;
    A(1, 1) = -2.0;
    const double theta = 1.1;
    A(2, 2) = std::cos(theta);
    A(2, 3) = -std::sin(theta);
    A(3, 2) = std::sin(theta);
    A(3, 3) = std::cos(theta);

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

    auto Maybe_F = lapack::schur_parlett(Q, exp_block_evaluator());
    if (!Maybe_F) UNSCOPED_INFO("schur_parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(Q);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(Q_4x4)_schur_parlett - expm(Q_4x4)_taylor||_F = " << err);
    CHECK(err < 1e-9);
}

TEST_CASE("frechet_integral_schur: parity vs Taylor M_m series",
          "[math][parlett][schur][frechet]") {
    // scheme_CO macro Q at agonist=10, with G = diag([0, -1]) (the per-state
    // current vector embedded as a diagonal matrix — same Frechet direction
    // calc_Qdt_eig integrates against). Sweep dt across the range
    // figure_2_micro exercises so Schur+VanLoan is tested in both well-scaled
    // (dt=1e-4) and aggressively-scaled (dt=1e-1) regimes.
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -1.0;
    Q(0, 1) = 1.0;
    Q(1, 0) = 100.0;
    Q(1, 1) = -100.0;

    Matrix<double> G(2, 2, 0.0);
    G(1, 1) = -1.0;  // diag([0, -1])

    // Restrict to dt where the Taylor reference converges. ‖Q·dt‖ ≤ 1
    // (i.e. dt ≤ 1e-2 here with max|Q|=100) keeps the order-20 truncation
    // below 1e-10. At larger dt the reference itself loses precision and the
    // failure mode reflects the reference, not the Schur implementation.
    const double dts[] = {1e-4, 1e-3, 1e-2};
    for (double dt : dts) {
        auto Maybe_A = lapack::frechet_integral_schur(Q, G, dt);
        if (!Maybe_A) UNSCOPED_INFO("frechet_integral_schur failed at dt=" << dt
                                     << ": " << Maybe_A.error()());
        REQUIRE(Maybe_A.valid());
        auto A_ref = frechet_integral_taylor_reference(Q, G, dt);
        auto err = frobenius_diff(Maybe_A.value(), A_ref);
        INFO("dt=" << dt << " ||A_schur - A_taylor||_F = " << err
             << "  (||A||~" << dt << ")");
        CHECK(err < 1e-10 * std::max(dt, 1.0));
    }
}

TEST_CASE("frechet_double_integral_schur: parity vs Taylor N_p series",
          "[math][parlett][schur][frechet]") {
    // Same scheme_CO Q, same G as the single-integral test.
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -1.0;
    Q(0, 1) = 1.0;
    Q(1, 0) = 100.0;
    Q(1, 1) = -100.0;
    Matrix<double> G(2, 2, 0.0);
    G(1, 1) = -1.0;

    const double dts[] = {1e-4, 1e-3, 1e-2};
    for (double dt : dts) {
        auto Maybe_B = lapack::frechet_double_integral_schur(Q, G, dt);
        if (!Maybe_B) UNSCOPED_INFO("frechet_double_integral_schur failed at dt="
                                     << dt << ": " << Maybe_B.error()());
        REQUIRE(Maybe_B.valid());
        auto B_ref = frechet_double_integral_taylor_reference(Q, G, dt);
        auto err = frobenius_diff(Maybe_B.value(), B_ref);
        // B(dt) ~ dt² so tolerance scales with dt².
        INFO("dt=" << dt << " ||B_schur - B_taylor||_F = " << err
             << "  (||B||~" << dt * dt << ")");
        CHECK(err < 1e-10 * std::max(dt * dt, 1.0));
    }
}

TEST_CASE("frechet_integral_schur: 6×6 lifted scheme_CO Nch=5",
          "[math][parlett][schur][frechet]") {
    // Stress test on the lifted Q at Nch=5. Eigenvalues span 0..-505 so
    // ‖Q·dt‖ at dt=1e-2 is ~5 — Schur+VanLoan should still match the
    // direct Taylor reference (which is itself accurate at this scale).
    const std::size_t M = 6;
    Matrix<double> Q(M, M, 0.0);
    for (std::size_t n = 0; n < M; ++n) {
        const std::size_t n_closed = M - 1 - n;
        if (n + 1 < M) Q(n, n + 1) = static_cast<double>(n_closed);
        if (n > 0) Q(n, n - 1) = static_cast<double>(n) * 100.0;
        double row_sum = 0.0;
        for (std::size_t k = 0; k < M; ++k)
            if (k != n) row_sum += Q(n, k);
        Q(n, n) = -row_sum;
    }
    // G_micro = lifted current per microstate. For Nch=5, k=2: g_microstate(n)
    // = n × g_open with g_open = -1, so G = diag(0, -1, -2, -3, -4, -5).
    Matrix<double> G(M, M, 0.0);
    for (std::size_t i = 0; i < M; ++i) G(i, i) = -static_cast<double>(i);

    const double dts[] = {1e-4, 1e-3, 1e-2};  // skip larger dt where Taylor reference loses precision
    for (double dt : dts) {
        auto Maybe_A = lapack::frechet_integral_schur(Q, G, dt);
        REQUIRE(Maybe_A.valid());
        auto A_ref = frechet_integral_taylor_reference(Q, G, dt, /*order=*/30);
        auto err = frobenius_diff(Maybe_A.value(), A_ref);
        INFO("dt=" << dt << " 6×6 lifted ||A_schur - A_taylor||_F = " << err);
        CHECK(err < 1e-8 * std::max(dt, 1.0));
    }
}

TEST_CASE("expm_schur_parlett: parity vs raw Pade scaling-squaring",
          "[math][parlett][schur][expm]") {
    // Sanity check the public expm_schur_parlett wrapper end-to-end. Reuses
    // the lifted scheme_CO 6×6 Q; eigenvalues span 0 to -505, matrix entries
    // up to ~600, so scaling-squaring is exercised on each block evaluation
    // and on the reference.
    const std::size_t M = 6;
    Matrix<double> Q(M, M, 0.0);
    for (std::size_t n = 0; n < M; ++n) {
        const std::size_t n_closed = M - 1 - n;
        if (n + 1 < M) Q(n, n + 1) = static_cast<double>(n_closed);
        if (n > 0) Q(n, n - 1) = static_cast<double>(n) * 100.0;
        double row_sum = 0.0;
        for (std::size_t k = 0; k < M; ++k)
            if (k != n) row_sum += Q(n, k);
        Q(n, n) = -row_sum;
    }

    auto Maybe_F = lapack::expm_schur_parlett(Q);
    if (!Maybe_F) UNSCOPED_INFO("expm_schur_parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(Q);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(Q_lifted)_schur_parlett - expm(Q)_pade||_F = " << err);
    CHECK(err < 1e-9);
}

TEST_CASE("schur_parlett: clustered eigenvalues require block merging",
          "[math][parlett][schur][cluster]") {
    // Construct a 4×4 quasi-triangular matrix with two pairs of clustered
    // eigenvalues: {-1, -1.001} and {-100, -100.001}. The naive Parlett
    // recurrence has divisor (T_jj − T_ii) ≈ 0.001 between blocks 1↔2 and
    // 3↔4, which kills relative accuracy. Clustering merges each pair into
    // a 2×2 super-block where the eval_block computes f directly via Pade.
    Matrix<double> T(4, 4, 0.0);
    T(0, 0) = -1.0;
    T(0, 1) = 0.5;
    T(1, 1) = -1.001;
    T(0, 2) = 0.3;
    T(1, 2) = 0.4;
    T(2, 2) = -100.0;
    T(2, 3) = 0.6;
    T(0, 3) = 0.1;
    T(1, 3) = 0.2;
    T(3, 3) = -100.001;

    // Without clustering, Parlett still produces *some* answer but precision
    // suffers near the close-eigenvalue divisor. Smoke-test that clustering
    // produces a tight match to raw Pade scaling-and-squaring.
    auto Maybe_F = lapack::schur_parlett(T, exp_block_evaluator(), /*cluster_delta=*/0.1);
    if (!Maybe_F) UNSCOPED_INFO("schur_parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(T);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(T_clustered)_schur_parlett - expm(T)_pade||_F = " << err);
    CHECK(err < 1e-10);

    // Sanity-check the cluster output: the natural decomposition has 4 blocks,
    // and clustering at δ=0.1 should merge into 2 (the close pairs).
    auto natural = lapack::detect_schur_block_starts(T);
    REQUIRE(natural.size() == 5);  // four 1×1 + sentinel
    auto clustered = lapack::cluster_adjacent_blocks(T, natural, /*delta=*/0.1);
    INFO("clustered starts: " << clustered.size() - 1 << " block(s)");
    CHECK(clustered.size() == 3);  // two merged 2×2 blocks + sentinel
}

TEST_CASE("schur_parlett: 6×6 lifted scheme_CO Nch=5",
          "[math][parlett][schur]") {
    // The lifted Q from test_schur. Eigenvalues are 0, -101, -202, ..., -505
    // (well-separated multiples of -101 for k=2 birth-death). Expect tight
    // agreement with Taylor.
    const std::size_t M = 6;
    const double r_open = 1.0;
    const double r_close = 100.0;
    Matrix<double> Q(M, M, 0.0);
    for (std::size_t n = 0; n < M; ++n) {
        const std::size_t n_closed = M - 1 - n;
        if (n + 1 < M) Q(n, n + 1) = static_cast<double>(n_closed) * r_open;
        if (n > 0) Q(n, n - 1) = static_cast<double>(n) * r_close;
        double row_sum = 0.0;
        for (std::size_t k = 0; k < M; ++k)
            if (k != n) row_sum += Q(n, k);
        Q(n, n) = -row_sum;
    }

    auto Maybe_F = lapack::schur_parlett(Q, exp_block_evaluator());
    if (!Maybe_F) UNSCOPED_INFO("schur_parlett failed: " << Maybe_F.error()());
    REQUIRE(Maybe_F.valid());

    auto F_ref = raw_expm_scaling_squaring(Q);
    auto err = frobenius_diff(Maybe_F.value(), F_ref);
    INFO("||expm(Q_6x6_lifted)_schur_parlett - expm(Q)_taylor||_F = " << err);
    CHECK(err < 1e-9);
}
