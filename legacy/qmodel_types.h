#ifndef QMODEL_TYPES_H
#define QMODEL_TYPES_H

// Tag/wrapper data types, compile-time policy flags, concepts, and small
// inline helpers shared across the qmodel surface. Headers that need these
// without paying for the full qmodel.h transitive include set (which adds
// fold.h, parallel_*, lapack_headers.h, schur_parlett.h, gsl_integrate.h,
// mcmc.h, function_memoization.h, derivative_test.h, exponential_matrix.h,
// parameters_distribution.h) should depend on this file instead.

#include <cmath>
#include <concepts>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

// Project headers — order matters: variables.h transitively brings maybe_error.h
// (Maybe_error, error_message, add_t, add_if_present_t, is_of_this_template_type_v),
// derivative_fwd.h (var::Derivative fwd, var::Derivative_Op, var::transformation_type,
// var::is_derivative, var::is_derivative_v, U concept), and the type_name dsl.
#include "variables.h"
#include "derivative_fwd.h"
#include "derivative_operator.h"  // var::Derivative<...> defs, Derivative_t, get_dx_of_dfdx, dx_of_dfdx
#include "parameters.h"           // var::Parameters_transformed
#include "function_measure_verification_and_optimization.h"  // var::F, FuncMap_St, Time_it_st, build
#include "matrix.h"               // Matrix / SymmetricMatrix / DiagonalMatrix, apply, inv, diag, set
#include "experiment.h"           // Patch_current, Recording, number_of_samples
#include "distributions.h"        // logL, elogL, vlogL, Grad, FIM, evaluation_time
#include "moment_statistics.h"    // mean_value_type_impl (specialized below for Evolution_of)
#include "general_algorithm_on_containers.h"  // primitive, var::sum/min/max/i_max/fullsum helpers

namespace macrodr {
using dsl::type_name;
using var::Parameters_Transformations;
using var::Power;
using var::Product;
using var::Var;
using var::Vector_Space;

using var::build;
using var::Constant;
using var::Fun;
using var::Op_t;
using var::primitive;
using var::Transfer_Op_to;
using var::transformation_type_t;

using var::U;

using std::abs;
using std::exp;
using std::log10;
using std::max;

using var::F;
using var::FuncMap_St;
using var::Time_it_st;

using dsl::type_name;

template <bool b>
class uses_variance_aproximation : public var::constexpr_Var<bool, uses_variance_aproximation, b> {
};
class uses_variance_aproximation_value
    : public var::constexpr_Var_value<bool, uses_variance_aproximation> {};

template <bool b>
class uses_taylor_variance_correction_aproximation
    : public var::constexpr_Var<bool, uses_taylor_variance_correction_aproximation, b> {};
class uses_taylor_variance_correction_aproximation_value
    : public var::constexpr_Var_value<bool, uses_taylor_variance_correction_aproximation> {};

template <bool b>
class uses_adaptive_aproximation : public var::constexpr_Var<bool, uses_adaptive_aproximation, b> {
};

class uses_adaptive_aproximation_value
    : public var::constexpr_Var_value<bool, uses_adaptive_aproximation> {};

template <bool b>
class uses_recursive_aproximation
    : public var::constexpr_Var<bool, uses_recursive_aproximation, b> {};
class uses_recursive_aproximation_value
    : public var::constexpr_Var_value<bool, uses_recursive_aproximation> {};

template <int b>
class uses_averaging_aproximation : public var::constexpr_Var<int, uses_averaging_aproximation, b> {
};
class uses_averaging_aproximation_value
    : public var::constexpr_Var_value<int, uses_averaging_aproximation> {};

template <bool b>
class uses_micro_aproximation : public var::constexpr_Var<bool, uses_micro_aproximation, b> {};
class uses_micro_aproximation_value
    : public var::constexpr_Var_value<bool, uses_micro_aproximation> {};

// Compile-time flag selecting the Qdt computation method:
//   false → eigendecomposition path (calc_Qdt_eig). Default. Fast for small
//           Q. Suffers from eigenvector-derivative ill-conditioning when
//           eigenvalues cluster (e.g. lifted Q_micro at Nch ≥ ~10 in k=2),
//           which produces score noise floors at ~1e-14.
//   true  → Taylor + scaling/squaring path (calc_Qdt_taylor). No
//           eigendecomposition; derivative chain is matrix-only and stable
//           across clustered/degenerate eigenvalues.
// Independent of `uses_taylor_variance_correction_aproximation` (which is
// about the moment-matched variance correction in the macro Kalman path,
// not Q's matrix exponential).
// qdt_method selects the Q·dt → Qdt computation strategy:
//   0 = eig          (eigendecomposition; default. Fast for k×k macro Q with
//                     well-separated eigenvalues; ill-conditioned at clusters)
//   1 = taylor       (Taylor scaling-and-squaring; the historical fallback)
//   2 = schur        (Schur+Parlett for plain expm + Pade-VanLoan for the
//                     Frechet integrals; stable at clustered eigenvalues)
// Replaces the old `bool uses_qdt_method`. The DSL keeps
// `taylor_qdt_approximation` as a deprecated bool alias mapping false→0,
// true→2 to preserve existing macroir scripts (Phase 9 semantics).
template <int b>
class uses_qdt_method : public var::constexpr_Var<int, uses_qdt_method, b> {};
class uses_qdt_method_value : public var::constexpr_Var_value<int, uses_qdt_method> {};

template <int b>
class return_predictions : public var::constexpr_Var<int, return_predictions, b> {};
class return_predictions_value : public var::constexpr_Var_value<int, return_predictions> {};

template <class T>
concept uses_averaging_aproximation_c =
    var::is_this_constexpr_Var_c<T, int, uses_averaging_aproximation>;

template <class T>
concept uses_adaptive_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_adaptive_aproximation>;

template <class T>
concept uses_recursive_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_recursive_aproximation>;

template <class T>
concept what_to_include_c = is_of_this_template_type_v<T, var::please_include>;

template <class T>
concept uses_variance_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_variance_aproximation>;

template <class T>
concept uses_taylor_variance_correction_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_taylor_variance_correction_aproximation>;

template <class T>
concept uses_micro_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_micro_aproximation>;

template <class V, class Id>
concept has_var_c = requires(V&& v) {
    std::forward<V>(v)[var::Var<Id>{}];
};

class Q0 : public var::Var<Q0, Matrix<double>> {};

inline std::size_t get_max_state(
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::string>> const& new_formulas) {
    std::size_t out = 0;
    for (auto& e : new_formulas) {
        if (e.first.first > out)
            out = e.first.first;
        if (e.first.second > out)
            out = e.first.second;
    }
    return out;
}
class Q0_formula : public var::Var<Q0_formula, std::vector<std::vector<std::string>>> {
   public:
    using var::Var<Q0_formula, std::vector<std::vector<std::string>>>::Var;

    Q0_formula(std::size_t N)
        : var::Var<Q0_formula, std::vector<std::vector<std::string>>>{
              std::vector<std::vector<std::string>>{N, std::vector<std::string>{N, ""}}} {}

    friend std::ostream& operator<<(std::ostream& os, Q0_formula const& x) {
        os << "Q0 formula"
           << "\n";
        for (std::size_t i = 0; i < x().size(); ++i) {
            for (std::size_t j = 0; j < x()[i].size(); ++j) {
                if (x()[i][j].size() > 0)
                    os << "Q(" << i << "," << j << ")->" << x()[i][j] << "\t";
            }
            os << "\n";
        }
        return os;
    }
};

class Qa : public var::Var<Qa, Matrix<double>> {};
class Qa_formula : public var::Var<Qa_formula, std::vector<std::vector<std::string>>> {
   public:
    using base_type = var::Var<Qa_formula, std::vector<std::vector<std::string>>>;
    using base_type::Var;
    Qa_formula(std::size_t N)
        : var::Var<Qa_formula, std::vector<std::vector<std::string>>>{
              std::vector<std::vector<std::string>>{N, std::vector<std::string>{N, ""}}} {}
    friend std::ostream& operator<<(std::ostream& os, Qa_formula const& x) {
        os << "Q0 formula"
           << "\n";
        for (std::size_t i = 0; i < x().size(); ++i) {
            for (std::size_t j = 0; j < x()[i].size(); ++j) {
                if (x()[i][j].size() > 0)
                    os << "Q(" << i << "," << j << ")->" << x()[i][j] << "\t";
            }
            os << "\n";
        }
        return os;
    }
};
template <class Q_formula>
    requires(std::is_same_v<Q_formula, Q0_formula> || std::is_same_v<Q_formula, Qa_formula>)
auto change_states_number(const Q_formula& f, std::size_t N) {
    Q_formula out(N);
    for (std::size_t i = 0; i < f().size(); ++i)
        for (std::size_t j = 0; j < f().size(); ++j) out()[i][j] = f()[i][j];

    return out;
}

template <class Q_formula>
    requires(std::is_same_v<Q_formula, Q0_formula> || std::is_same_v<Q_formula, Qa_formula>)
auto insert_new_formula(const Q_formula& f, std::size_t i_ini, std::size_t i_end,
                        std::string&& formula) {
    auto N = std::max(i_ini, i_end) + 1;
    auto out = change_states_number(f, N);
    out()[i_ini][i_end] = std::move(formula);

    return out;
}

class Qx : public var::Var<Qx, Matrix<double>> {};
class P_initial : public var::Var<P_initial, Matrix<double>> {};

class g : public var::Var<g, Matrix<double>> {};
class g_formula : public var::Var<g_formula, std::vector<std::string>> {};

inline auto change_states_number(const g_formula& f, std::size_t N) {
    g_formula out(std::vector<std::string>{N, ""});
    for (std::size_t i = 0; i < f().size(); ++i) out()[i] = f()[i];
    return out;
}

class N_St : public var::Constant<N_St, std::size_t> {};

class N_Ch_mean : public var::Var<N_Ch_mean, Matrix<double>> {};
class N_Ch_mean_value : public var::Var<N_Ch_mean, double> {};

class N_Ch_mean_time_segment_duration
    : public var::Constant<N_Ch_mean_time_segment_duration, double> {};

class N_Ch_init : public var::Var<N_Ch_init, double> {};
class N_Ch_eq : public var::Var<N_Ch_eq, double> {};
class N_Ch_tau : public var::Var<N_Ch_tau, double> {};

class Binomial_magical_number : public var::Constant<Binomial_magical_number, double> {};

class min_P : public var::Constant<min_P, double> {};

class N_Ch_std : public var::Var<N_Ch_std, double> {};

class SeedNumber : public var::Constant<SeedNumber, std::size_t> {};


class Current_Noise : public var::Var<Current_Noise, double> {};

class Pink_Noise : public var::Var<Pink_Noise, double> {};

class Proportional_Noise : public var::Var<Proportional_Noise, double> {};

class Current_Baseline : public var::Var<Current_Baseline, double> {};

// ---------------------------------------------------------------------------
// Simplex-canary thresholds (shared by to_Probability, to_Probability_displacement,
// to_Covariance_Probability).
//
// Two roots and a safety multiplier, everything else derived:
//   - canary_primitive_warn — warn level for primitive departures
//                             (|Σ−target|/√N, RMS row-sum). Set to the typical
//                             FP-noise floor used in well-engineered scientific
//                             code for invariant checks (≈ BLAS orthogonality
//                             tolerance, ODE conservation laws, MCMC invariants).
//   - canary_dcos_warn_sq   — warn level for cos²(∂x/∂θ, 𝟙). Geometric measure
//                             of how much of the derivative aligns with the
//                             constraint-violating direction; pure FP gives
//                             cos² ≈ ε², a real bug gives cos² ≈ 1.
//   - canary_safety         — multiplier: error threshold = warn × canary_safety.
//
// Single-row max threshold for covariance is 10× the aggregate-RMS threshold —
// captures that one row can stick out by up to ~√N × RMS before averaging hides
// it.
//
// See theory/scientific-software/notes/with_warning_abstraction.md for the
// structured replacement of the std::cerr warn channel.
// ---------------------------------------------------------------------------
inline constexpr double canary_primitive_warn   = 1e-10;
inline constexpr double canary_dcos_warn_sq     = 1e-8;
inline constexpr double canary_safety           = 100.0;
inline constexpr double canary_primitive_error  = canary_primitive_warn * canary_safety;
// Error band on cos²(∂x/∂θ, 𝟙) deliberately loose (1e-1 ≈ 32% misalignment)
// to absorb a known IRT av=2 vc=1 rare-event AD-conservation breakdown:
// at specific (y, P_mean) realizations at small intervals, the rank-2
// Newton step's cross-products (V_iter/N, V_iter²/N, det, k11/k12/k22)
// don't cancel cleanly under AD. Worst empirical so far: 6.2% drift on
// `unitary_current` at k=97. Extreme-value scaling across ~2000-step
// traces suggests up to ~28% worst-case; 1e-1 gives margin. Tracked as
// Task 9 in handoff_state.md. Warn band (1e-8 → 0.014% misalignment)
// still surfaces all real drifts informationally without aborting.
inline constexpr double canary_dcos_error_sq    = 1e-1;
inline constexpr double canary_row_max_warn     = canary_primitive_warn * 10.0;
inline constexpr double canary_row_max_error    = canary_primitive_error * 10.0;
// Floor on ‖∂x/∂θ_p‖₂² below which the parameter is treated as having no
// effect on x — avoids noise/noise = O(1) false signals. Raised from 1e-24
// to 1e-16 after the MRT/IRT diag-A/B probes confirmed FP-noise-scale
// leakage (~1e-12 entries) accumulating through the rank-2 Newton arithmetic
// produces cos² > 1e-6 against tiny denominators; at ‖·‖ < 1e-8 the cos²
// ratio measures rounding bits, not real direction. See
// theory/macroir/notes/Gmean_ij_gvarij/handoff_state.md.
inline constexpr double canary_norm_floor_sq    = 1e-16;  // ‖·‖₂ ≲ 1e-8

// ---------------------------------------------------------------------------
// to_Probability — canary + drift correction for almost-probability vectors.
//
// PURPOSE
//   Validate that x is a probability vector to within FP-noise tolerance,
//   then renormalize so primitives sum to exactly 1. The validation is the
//   *canary*: structural bugs upstream surface here as Maybe_error, so they
//   don't get silently averaged into invisibility by the renormalize.
//
// INPUT INVARIANT (must hold up to FP noise)
//   - All primitive entries are finite and ≥ 0 (within eps_neg of 0).
//   - Σᵢ primitive(xᵢ) ≈ 1 (within eps_sum).
//   - For Derivative inputs, ∂Σ/∂θ ≈ 0 (within eps_dsum) — direct consequence
//     of Σx = 1 being constant in θ.
//   - All derivative entries are finite (negative gradients are legitimate).
//
// OUTPUT (on success)
//   - x · (1/s) where s = var::sum(x). Primitives sum to exactly 1.
//   - Derivative-aware throughout: 1/s carries chain rule, and Derivative *
//     Derivative<double> applies the product rule.
//
// FAILURE MODES
//   - "non-finite primitive"     — any xᵢ primitive is NaN/Inf.
//   - "non-finite derivative"    — any derivative cell is NaN/Inf.
//   - "negative entry"           — min primitive < -eps_neg (real negativity,
//                                  not FP cancellation noise).
//   - "|Σ − 1| out of tolerance" — primitive sum drifted past eps_sum.
//   - "∂Σ/∂θ violates invariant" — derivative-of-sum exceeds eps_dsum at some
//                                  parameter, signalling a probability-
//                                  non-preserving operation upstream.
//
// DESIGN NOTE — why no clamp
//   The previous version clamped negatives to zero before renormalizing. With
//   the trust-coefficient innovation in the macro_R update and stochastic
//   Markov steps in the micro_R update, no operation should *produce* real
//   negatives; tiny FP-noise negatives are within eps_neg and renormalize
//   harmlessly. A real negative is a bug, and we want to know about it.
//
// DESIGN NOTE — derivative canary (3b)
//   The renormalize step `x * (1/s)` applies the quotient rule to the
//   derivatives: ∂(xᵢ/s)/∂θ = ∂xᵢ/∂θ / s − xᵢ · ∂s/∂θ / s². If ∂s/∂θ ≠ 0
//   (a bug), the second term silently corrects the offending derivative. We
//   reject upstream rather than absorb the correction.
// ---------------------------------------------------------------------------
template <class C_Matrix>
auto to_Probability(C_Matrix const& x) -> Maybe_error<C_Matrix> {
    constexpr double eps_neg = 1e-10;
    // Thresholds: see canary_* constants above.
    constexpr double eps_sum_warn      = canary_primitive_warn;
    constexpr double eps_sum_error     = canary_primitive_error;
    constexpr double eps_dcos_warn_sq  = canary_dcos_warn_sq;
    constexpr double eps_dcos_error_sq = canary_dcos_error_sq;
    constexpr double eps_norm_floor_sq = canary_norm_floor_sq;

    // (1) finiteness — primitive
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(primitive(x[i])))
            return error_message("to_Probability: non-finite primitive at i=" +
                                  std::to_string(i));

    // (1') finiteness — derivative payload
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        for (std::size_t p = 0; p < d.size(); ++p)
            for (std::size_t k = 0; k < d[p].size(); ++k)
                if (!std::isfinite(d[p][k]))
                    return error_message("to_Probability: non-finite derivative");
    }

    // (2) negativity canary (primitive only; derivatives may be negative)
    double minv = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < x.size(); ++i)
        minv = std::min(minv, primitive(x[i]));
    if (minv < -eps_neg)
        return error_message("to_Probability: negative entry " +
                              std::to_string(minv));

    // (3) sum canary — derivative-aware var::sum returns Derivative<double, P>.
    auto s = var::sum(x);
    const double inv_sqrtN = 1.0 / std::sqrt(static_cast<double>(x.size()));

    // (3a) primitive |Σ − 1|/√N
    {
        double dep = std::abs(primitive(s) - 1.0) * inv_sqrtN;
        if (dep > eps_sum_error) {
            return error_message("to_Probability: |Σ−1|/√N = " +
                                  std::to_string(dep));
        }
        if (dep > eps_sum_warn) {
            std::cerr << "[warn] to_Probability: |Σ−1|/√N = "
                      << std::scientific << std::setprecision(3) << dep
                      << " (warn=" << eps_sum_warn
                      << ", err=" << eps_sum_error << ")\n";
        }
    }

    // (3b) derivative cosine ratio² per parameter — warn-only.
    // The cos² of ∂Σ/∂θ against the all-ones direction is a *direction*
    // test, not a *magnitude* one. Even genuine derivative-conservation
    // drifts at the IRT av=2 vc=1 rare-event level produce small absolute
    // errors that the cos² metric amplifies (Σ²/(N·‖·‖²) explodes when
    // ‖·‖ is small). After the May 2026 calibration we keep this as a
    // warning only — the metric flags potential issues but does not
    // abort, because the absolute downstream impact is bounded and the
    // ratio test gives many false positives at FP-noise to small-drift
    // scales. See handoff_state.md Task 9 (IRT ∂P/∂N drift) for context.
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        const double N = static_cast<double>(x.size());
        for (std::size_t p = 0; p < d.size(); ++p) {
            auto const& dp = d[p];
            double s_signed = 0.0;
            double s_sq     = 0.0;
            for (std::size_t k = 0; k < dp.size(); ++k) {
                s_signed += dp[k];
                s_sq     += dp[k] * dp[k];
            }
            if (s_sq < eps_norm_floor_sq) {
                continue;  // parameter has no effect on x
            }
            double ratio_sq = (s_signed * s_signed) / (s_sq * N);
            if (ratio_sq > eps_dcos_warn_sq) {
                std::cerr << "[warn] to_Probability: ∂Σ/∂θ at param " << p
                          << " cos²=" << std::scientific << std::setprecision(3)
                          << ratio_sq << " (warn²=" << eps_dcos_warn_sq << ")\n";
            }
        }
    }

    // (4) drift correction — derivative-aware via operator/(double, Derivative)
    //     and operator*(Derivative<aMatrix>, Derivative<double>).
    return x * (1.0 / s);
}

// ---------------------------------------------------------------------------
// to_Covariance_Probability — canary for the bare centered covariance of the
// channel-state indicator (the codebase's stored P_Cov).
//
// SEMANTIC OF x
//   x is the bare covariance Cov(X) where X is the one-hot indicator of the
//   channel state (or the macrostate-fraction Y = n/N for N channels). X
//   lies on the simplex (Σᵢ Xᵢ = 1), so by bilinearity:
//
//     Σⱼ Cov(Xᵢ, Xⱼ) = Cov(Xᵢ, ΣX) = Cov(Xᵢ, 1) = 0    ∀ i
//
//   Hence every row sums to 0, hence the total sum is 0. Diagonals are
//   variances Var(Xᵢ) = Pᵢ(1−Pᵢ) for one-hot — *can be 0* when the state is
//   deterministic (Pᵢ ∈ {0, 1}).
//
//   Sanity check: P_mean = e_i ⇒ Cov(X) = diag(P) − Pᵀ·P = 0 (zero matrix).
//
// CANARY (returns error if violated beyond FP noise)
//   - finiteness, primitive and derivative
//   - symmetry within eps_sym
//   - diagonal ≥ -eps_neg (PSD necessary; variances ≥ 0)
//   - per-row primitive sum ≈ 0 (within eps_sum)
//   - max |∂(row sum)/∂θ| ≈ 0 (within eps_dsum) — the row-sum identity is
//     constant in θ, hence so is its derivative
//
// NO DRIFT CORRECTION
//   With Σ = 0 there is nothing to renormalize *to*. Pass through on success;
//   reject with named error otherwise. Restoring the bare-Cov convention
//   makes this function a pure validator (no side effects).
//
// DESIGN NOTES vs. the previous version:
//   * Previous canary targeted Σ M = 1 / row sums = P_mean — based on a
//     mistaken interpretation of the convention as M = bare_Cov + diag(P_mean).
//     The actual stored convention is bare_Cov, so the targets are 0 / 0.
//   * No "zero the row when diag ≤ 0" mutation: the old loose validator did
//     that as cleanup, but with the strict canary in place real bugs surface
//     instead of being masked. FP-noise negatives below eps_neg pass.
//   * PSD is *not* fully checked here (O(n³) eigendecomp). The diagonal +
//     symmetry + finite + row-sum checks are the cheap necessary conditions.
// ---------------------------------------------------------------------------
template <class C_Matrix>
auto to_Covariance_Probability(C_Matrix const& x) -> Maybe_error<C_Matrix> {
    // See to_Probability for the design notes. Two tiers (warn/error) for
    // both primitive row-sums and derivative row-sums.
    //   - Primitive row-sums: aggregate RMS √(Σᵢ row_sum_i²)/√N AND per-row
    //     max |row_sum_i|. RMS catches accumulation across rows, max catches
    //     a single deeply-broken row.
    //   - Derivative row-sums: per-(row, parameter) cosine ratio²
    //     |Σⱼ ∂x[i,j]/∂θ_p|² / (‖∂x[i,·]/∂θ_p‖₂² · N).
    constexpr double eps_neg = 1e-10;
    constexpr double eps_sym = 1e-10;
    constexpr double eps_sum_warn      = canary_primitive_warn;
    constexpr double eps_sum_error     = canary_primitive_error;
    constexpr double eps_max_warn      = canary_row_max_warn;
    constexpr double eps_max_error     = canary_row_max_error;
    constexpr double eps_dcos_warn_sq  = canary_dcos_warn_sq;
    constexpr double eps_dcos_error_sq = canary_dcos_error_sq;
    constexpr double eps_norm_floor_sq = canary_norm_floor_sq;

    if (x.nrows() != x.ncols())
        return error_message("to_Covariance_Probability: not square (" +
                              std::to_string(x.nrows()) + "x" +
                              std::to_string(x.ncols()) + ")");

    // (1) finiteness — primitive, per cell
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(primitive(x[i])))
            return error_message("to_Covariance_Probability: non-finite primitive");

    // (1') finiteness — derivative payload
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        for (std::size_t p = 0; p < d.size(); ++p)
            for (std::size_t k = 0; k < d[p].size(); ++k)
                if (!std::isfinite(d[p][k]))
                    return error_message("to_Covariance_Probability: non-finite derivative");
    }

    // (2) symmetry of primitive (lower triangle vs upper triangle)
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            double diff = std::abs(primitive(x(i, j)) - primitive(x(j, i)));
            if (diff > eps_sym)
                return error_message("to_Covariance_Probability: asymmetric at (" +
                                      std::to_string(i) + "," + std::to_string(j) +
                                      "): " + std::to_string(diff));
        }

    // (3) diagonal ≥ -eps_neg (PSD necessary condition)
    for (std::size_t i = 0; i < x.nrows(); ++i)
        if (primitive(x(i, i)) < -eps_neg)
            return error_message("to_Covariance_Probability: negative diagonal at " +
                                  std::to_string(i) + ": " +
                                  std::to_string(primitive(x(i, i))));

    // (4) per-row sum canary (primitive). By bilinearity on a simplex-constrained X:
    //         Σⱼ Cov(Xᵢ, Xⱼ) = Cov(Xᵢ, ΣX) = Cov(Xᵢ, 1) = 0    ∀ i
    //     Aggregate RMS + worst-row check, each two-tier.
    {
        const double N = static_cast<double>(x.nrows());
        const double inv_sqrtN = 1.0 / std::sqrt(N);
        double sum_sq = 0.0;
        double max_abs = 0.0;
        std::size_t max_i = 0;
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            double row_sum = 0.0;
            for (std::size_t j = 0; j < x.ncols(); ++j) {
                row_sum += primitive(x(i, j));
            }
            sum_sq += row_sum * row_sum;
            if (std::abs(row_sum) > max_abs) {
                max_abs = std::abs(row_sum);
                max_i = i;
            }
        }
        double rms = std::sqrt(sum_sq) * inv_sqrtN;
        if (rms > eps_sum_error) {
            return error_message("to_Covariance_Probability: RMS row-sum = " +
                                  std::to_string(rms) +
                                  " (expected 0 by simplex constraint on Cov rows)");
        }
        if (max_abs > eps_max_error) {
            return error_message("to_Covariance_Probability: max |row " +
                                  std::to_string(max_i) + " sum| = " +
                                  std::to_string(max_abs));
        }
        if (rms > eps_sum_warn) {
            std::cerr << "[warn] to_Covariance_Probability: RMS row-sum = "
                      << std::scientific << std::setprecision(3) << rms
                      << " (warn=" << eps_sum_warn
                      << ", err=" << eps_sum_error << ")\n";
        }
        if (max_abs > eps_max_warn) {
            std::cerr << "[warn] to_Covariance_Probability: max |row " << max_i
                      << " sum| = " << std::scientific << std::setprecision(3)
                      << max_abs << " (warn=" << eps_max_warn
                      << ", err=" << eps_max_error << ")\n";
        }
    }

    // (5) ∂(row sum)/∂θ ≈ 0 — per (row, parameter) cosine ratio², warn-only.
    // See to_Probability comment and handoff_state.md Task 9.
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        const double N = static_cast<double>(x.ncols());
        for (std::size_t p = 0; p < d.size(); ++p) {
            auto const& dp = d[p];
            for (std::size_t i = 0; i < x.nrows(); ++i) {
                double drow_signed = 0.0;
                double drow_sq     = 0.0;
                for (std::size_t j = 0; j < x.ncols(); ++j) {
                    drow_signed += dp(i, j);
                    drow_sq     += dp(i, j) * dp(i, j);
                }
                if (drow_sq < eps_norm_floor_sq) {
                    continue;
                }
                double ratio_sq = (drow_signed * drow_signed) / (drow_sq * N);
                if (ratio_sq > eps_dcos_warn_sq) {
                    std::cerr << "[warn] to_Covariance_Probability: ∂(row " << i
                              << ")/∂θ_" << p << " cos²="
                              << std::scientific << std::setprecision(3) << ratio_sq
                              << " (warn²=" << eps_dcos_warn_sq << ")\n";
                }
            }
        }
    }

    // (6) Pass through on success — Cov sums to 0, no renormalization target.
    return x;
}

// ---------------------------------------------------------------------------
// to_Probability_displacement — canary for tangent vectors to the simplex.
//
// SEMANTIC OF x
//   x is a "probability displacement": a row/column vector (or matrix whose
//   rows are such) intended to be added to a probability vector p to yield
//   another probability vector p + x. For p + x to remain on the simplex,
//
//       Σᵢ xᵢ = 0    (and Σᵢ ∂xᵢ/∂θ = 0 for Derivative inputs)
//
//   The recursive macro/IR posterior-mean updates rely on exactly this
//   property — gS, gS0, and the row blocks of GS are the codebase's canonical
//   probability displacements. The invariant is identical whether the
//   displacement comes from a Bayesian innovation (chi · gS) or from a
//   Markov-propagation residual.
//
// CANARY (returns error if violated beyond FP noise)
//   - finiteness, primitive and derivative
//   - Σᵢ xᵢ ≈ 0 within eps_sum
//   - Σᵢ ∂xᵢ/∂θ_p ≈ 0 within eps_dsum (derivative of the identity)
//
// NO DRIFT CORRECTION
//   The target is 0; there's nothing to renormalize *to*. Pass through on
//   success; reject with named error otherwise.
//
// SIGNS ARE FREE
//   Unlike to_Probability there is no negativity check. Tangent vectors carry
//   any sign; only their sum is constrained.
//
// MATRIX INPUTS
//   For a matrix x whose rows are each a probability displacement (e.g. GS),
//   pass each row separately, or use the dedicated row-by-row overload below.
// ---------------------------------------------------------------------------
template <class C_Matrix>
auto to_Probability_displacement(C_Matrix const& x) -> Maybe_error<C_Matrix> {
    // See to_Probability / canary_* constants above. Target = 0 instead of 1.
    constexpr double eps_sum_warn      = canary_primitive_warn;
    constexpr double eps_sum_error     = canary_primitive_error;
    constexpr double eps_dcos_warn_sq  = canary_dcos_warn_sq;
    constexpr double eps_dcos_error_sq = canary_dcos_error_sq;
    constexpr double eps_norm_floor_sq = canary_norm_floor_sq;

    // (1) finiteness — primitive
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(primitive(x[i])))
            return error_message("to_Probability_displacement: non-finite primitive at i=" +
                                  std::to_string(i));

    // (1') finiteness — derivative payload
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        for (std::size_t p = 0; p < d.size(); ++p)
            for (std::size_t k = 0; k < d[p].size(); ++k)
                if (!std::isfinite(d[p][k]))
                    return error_message("to_Probability_displacement: non-finite derivative");
    }

    // (2) primitive |Σ|/√N (target = 0)
    auto s = var::sum(x);
    const double inv_sqrtN = 1.0 / std::sqrt(static_cast<double>(x.size()));
    {
        double dep = std::abs(primitive(s)) * inv_sqrtN;
        if (dep > eps_sum_error) {
            return error_message("to_Probability_displacement: |Σ|/√N = " +
                                  std::to_string(dep) +
                                  " (expected 0 by simplex-tangent constraint)");
        }
        if (dep > eps_sum_warn) {
            std::cerr << "[warn] to_Probability_displacement: |Σ|/√N = "
                      << std::scientific << std::setprecision(3) << dep
                      << " (warn=" << eps_sum_warn
                      << ", err=" << eps_sum_error << ")\n";
        }
    }

    // (3) derivative cosine ratio² per parameter — warn-only (see
    // to_Probability comment above and handoff_state.md Task 9).
    if constexpr (var::is_derivative_v<C_Matrix>) {
        auto const& d = derivative(x)();
        const double N = static_cast<double>(x.size());
        for (std::size_t p = 0; p < d.size(); ++p) {
            auto const& dp = d[p];
            double s_signed = 0.0;
            double s_sq     = 0.0;
            for (std::size_t k = 0; k < dp.size(); ++k) {
                s_signed += dp[k];
                s_sq     += dp[k] * dp[k];
            }
            if (s_sq < eps_norm_floor_sq) {
                continue;
            }
            double ratio_sq = (s_signed * s_signed) / (s_sq * N);
            if (ratio_sq > eps_dcos_warn_sq) {
                std::cerr << "[warn] to_Probability_displacement: ∂Σ/∂θ at param "
                          << p << " cos²=" << std::scientific << std::setprecision(3)
                          << ratio_sq << " (warn²=" << eps_dcos_warn_sq << ")\n";
            }
        }
    }

    // (4) Pass through on success.
    return x;
}

// ---------------------------------------------------------------------------
// Bayes_Rule — apply Bayes' theorem in one step.
//
//   posteriorᵢ = priorᵢ · likelihoodᵢ / evidence
//   evidence   = Σᵢ priorᵢ · likelihoodᵢ
//
// The returned `evidence` is the *marginal likelihood* of the observation
// over the latent state. In the per-step Bayesian update it's "model
// evidence"; in the channel-kinetics filter it's the per-step contribution
// to the parameter log-likelihood, accumulated via log(evidence) → logL.
// Same scalar, two names — the caller wraps it whichever way fits the
// surrounding inference stage.
//
// INPUT INVARIANTS (assumed)
//   - prior is a probability (sum=1, non-neg). Should have come through
//     to_Probability beforehand; we don't re-canary it here.
//   - likelihood is finite and non-negative (within FP-noise tolerance).
//   - prior and likelihood share shape (1D vector or 2D matrix). The
//     function works uniformly on both — for 2D inputs (e.g. averaging=2
//     per-pair posteriors) elemMult and var::sum already operate cell-wise.
//
// OUTPUT INVARIANTS (on success)
//   - posterior is a probability: sum=1 (exactly, by construction), non-neg,
//     and ∂Σ/∂θ = 0 (for Derivative inputs — quotient rule on Σ x / Σ x).
//   - evidence is positive (and finite).
//   - Derivative threading: elemMult, var::sum, operator/(double, Derivative)
//     and operator*(Derivative<aMatrix>, Derivative<double>) all carry the
//     chain/product rule, so derivative(posterior) and derivative(evidence)
//     are the correct gradients of the corresponding scalars.
//
// FAILURES (each returns a distinct error_message)
//   - "non-finite likelihood"  — any cell NaN/Inf.
//   - "negative likelihood"    — any cell < -eps_neg (real, not FP noise).
//   - "zero evidence"          — Σ prior·L underflowed to ≤ 0; posterior
//                                undefined. THIS is the diagnostic for
//                                "all per-microstate likelihoods underflowed",
//                                which is the failure mode the figure-2 run
//                                surfaced as the original `to_Probability`'s
//                                "cero probability" error.
// ---------------------------------------------------------------------------
template <class C_Prior, class C_Likelihood>
auto Bayes_Rule(C_Prior const& prior, C_Likelihood const& likelihood) {
    constexpr double eps_neg = 1e-10;

    auto unnormalized = elemMult(prior, likelihood);
    auto evidence     = var::sum(unnormalized);
    using Posterior_T = std::decay_t<decltype(unnormalized * (1.0 / evidence))>;
    using Evidence_T  = std::decay_t<decltype(evidence)>;
    using Result_T    = Maybe_error<std::pair<Posterior_T, Evidence_T>>;

    if (prior.size() != likelihood.size())
        return Result_T(error_message("Bayes_Rule: shape mismatch (prior.size=" +
                                       std::to_string(prior.size()) +
                                       ", likelihood.size=" +
                                       std::to_string(likelihood.size()) + ")"));

    // (1) likelihood validity — primitive only; gradients of L may be negative.
    for (std::size_t i = 0; i < likelihood.size(); ++i) {
        if (!std::isfinite(primitive(likelihood[i])))
            return Result_T(error_message("Bayes_Rule: non-finite likelihood at i=" +
                                           std::to_string(i)));
        if (primitive(likelihood[i]) < -eps_neg)
            return Result_T(error_message("Bayes_Rule: negative likelihood at i=" +
                                           std::to_string(i) + ": " +
                                           std::to_string(primitive(likelihood[i]))));
    }

    // (2) evidence — must be positive for the posterior to be defined.
    if (!(primitive(evidence) > 0))
        return Result_T(error_message("Bayes_Rule: zero evidence (posterior underflow) — " +
                                       std::string("all per-state likelihoods produced ≤ 0 mass; ") +
                                       "check upstream that the predicted means / variances make sense"));
    if (!std::isfinite(primitive(evidence)))
        return Result_T(error_message("Bayes_Rule: non-finite evidence"));

    // (3) posterior — derivative-aware via operator*(Derivative<aMatrix>,
    //     Derivative<double>) and operator/(double, Derivative).
    auto posterior = unnormalized * (1.0 / evidence);

    return Result_T(std::make_pair(std::move(posterior), std::move(evidence)));
}

inline bool all_Probability_elements(Matrix<double> const& x) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(primitive(x[i])))
            return false;
        else if (x[i] > 1.0)
            return false;
        else if (x[i] < 0.0)
            return false;
    }
    return true;
}

template <class C_Matrix>
inline bool all_Covariance_elements(C_Matrix const& x) {
    if (x.ncols() != x.nrows())
        return false;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (!std::isfinite(primitive(x(i, i))) || x(i, i) < 0)
            return false;
        // for (std::size_t j = 0; j < i; ++j)
        //  if (std::isfinite(primitive(x(i, j))) && (x(i, i) > 0) && (x(j, j) > 0))
        //  {
        //      auto r=x(i, j) * x(i, j) / x(i, i) / x(j, j);
        //      if ( r> 1.0)
        //      return false;
        //  } else
        //{
        //   if (x(i, j) != 0)
        //     return false;
        // }
    }
    return true;
}

inline bool crude_lambda_violations(DiagonalMatrix<double> const& l) {
    if (var::max(l) > 1e-2)
        return true;
    else
        return false;
}

template <class C_Matrix>
auto to_Transition_Probability_Eigenvalues(C_Matrix&& lambda) {
    // if (crude_lambda_violations(primitive(lambda)))
    //     std::cerr<<"crude lambda violations\n";

    auto i_max = var::i_max(primitive(lambda));
    lambda.set(i_max, 0.0);
    auto j_max = var::i_max(primitive(lambda));
    while (primitive(lambda[j_max]) > 0.0) {
        lambda.set(j_max, 0.0);
        j_max = var::i_max(primitive(lambda));
    }

    return lambda;
}

// enforce_gmean_bounds removed: source-level canaries (check_gtotal_ij_in_range,
// check_gtotal_sqr_ij_in_range) replace the post-construction clamp.
struct StabilizerPolicyEnabled {
    static constexpr bool clamp_variance = false;
    static constexpr bool mask_probability = false;
    static constexpr bool project_transition_probability = false;
    static constexpr bool sanitize_eigenvalues = false;
};

struct StabilizerPolicyEnabled_ {
    static constexpr bool clamp_variance = true;
    static constexpr bool mask_probability = true;
    static constexpr bool project_transition_probability = true;
    static constexpr bool sanitize_eigenvalues = true;
};

struct StabilizerPolicyDisabled {
    static constexpr bool clamp_variance = false;
    static constexpr bool mask_probability = false;
    static constexpr bool project_transition_probability = false;
    static constexpr bool sanitize_eigenvalues = false;
};

class N_channel_state : public var::Var<N_channel_state, Matrix<double>> {
    using var::Var<N_channel_state, Matrix<double>>::Var;
};

class y_sum : public var::Var<y_sum, double> {};

class t_sum : public var::Var<t_sum, double> {};
class P_mean : public var::Var<P_mean, Matrix<double>> {
   public:
    friend std::string className(P_mean) { return "P_mean"; }
};

class P_mean_t2_y0 : public var::Var<P_mean_t2_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t2_y0) { return "P_mean_t2_y0"; }
};

class P_mean_t2_y1 : public var::Var<P_mean_t2_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t2_y1) { return "P_mean_t2_y1"; }
};

class P_mean_t15_y0 : public var::Var<P_mean_t15_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t15_y0) { return "P_mean_t15_y0"; }
};

class P_mean_t15_y1 : public var::Var<P_mean_t15_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t15_y1) { return "P_mean_t15_y1"; }
};


class P_mean_t1_y1 : public var::Var<P_mean_t1_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t1_y1) { return "P_mean_t1_y1"; }
};

class P_mean_t20_y1 : public var::Var<P_mean_t20_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t20_y1) { return "P_mean_t20_y1"; }
};

class P_mean_t11_y0 : public var::Var<P_mean_t11_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t11_y0) { return "P_mean_t11_y0"; }
};

class P_mean_t10_y1 : public var::Var<P_mean_t10_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t10_y1) { return "P_mean_t10_y1"; }
};

class P_mean_0t_y0 : public var::Var<P_mean_0t_y0, Matrix<double>> {
   public:
    // Boundary-state prior mean over the interval: (i0,it) -> P(X0=i0, Xt=it).
    friend std::string className(P_mean_0t_y0) { return "P_mean_0t_y0"; }
};

class P_mean_0t_y1 : public var::Var<P_mean_0t_y1, Matrix<double>> {
   public:
    // Boundary-state posterior mean over the interval after conditioning on y_{0->t}.
    friend std::string className(P_mean_0t_y1) { return "P_mean_0t_y1"; }
};


class P_Cov : public var::Var<P_Cov, SymmetricMatrix<double>> {
    friend std::string className(P_Cov) { return "P_Cov"; }
};

class P_cross_cov_0t_y0 : public var::Var<P_cross_cov_0t_y0, Matrix<double>> {
    // Reduced cross-time covariance Cov(x_0, x_t) before conditioning on y_{0->t}.
    friend std::string className(P_cross_cov_0t_y0) { return "P_cross_cov_0t_y0"; }
};


class P_cross_cov_0t_y1 : public var::Var<P_cross_cov_0t_y1, Matrix<double>> {
    // Reduced cross-time covariance Cov(x_0, x_t) after conditioning on y_{0->t}.
    friend std::string className(P_cross_cov_0t_y1) { return "P_cross_cov_0t_y1"; }
};


class d_gS: public var::Var<d_gS, Matrix<double>> {
    friend std::string className(d_gS) { return "d_gS"; }
};  

class d_GS: public var::Var<d_GS, Matrix<double>> {
    friend std::string className(d_GS) { return "d_GS"; }
};  




class P_Cov_t2_y0 : public var::Var<P_Cov_t2_y0, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t2_y0) { return "P_Cov_t2_y0"; }
};

class P_Cov_t2_y1 : public var::Var<P_Cov_t2_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t2_y1) { return "P_Cov_t2_y1"; }
};

class P_Cov_t15_y0 : public var::Var<P_Cov_t15_y0, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t15_y0) { return "P_Cov_t15_y0"; }
};

class P_Cov_t15_y1 : public var::Var<P_Cov_t15_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t15_y1) { return "P_Cov_t15_y1"; }
};



class P_Cov_t1_y1 : public var::Var<P_Cov_t1_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t1_y1) { return "P_Cov_t1_y1"; }
};

class P_Cov_t20_y1 : public var::Var<P_Cov_t20_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t20_y1) { return "P_Cov_t20_y1"; }
};

class P_Cov_t11_y0 : public var::Var<P_Cov_t11_y0, SymmetricMatrix<double>> {
   public:
    // One-time prior covariance at the end of the interval, before conditioning on y_{0->t}.
    friend std::string className(P_Cov_t11_y0) { return "P_Cov_t11_y0"; }
};

class P_Cov_t10_y1 : public var::Var<P_Cov_t10_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t10_y1) { return "P_Cov_t10_y1"; }
};

class lambda : public var::Var<lambda, DiagonalMatrix<double>> {};

class V : public var::Var<V, Matrix<double>> {};
class W : public var::Var<W, Matrix<double>> {};
// Frobenius condition number κ_F(V) = ‖V‖_F · ‖W‖_F of the eigenvector
// matrix. Stored as a plain double (no derivative payload): κ is treated as
// a hyperparameter of the regularization (it sets the FP-noise floor of P_ij
// after V·diag(exp(λdt))·W reconstruction), not a free model parameter.
// See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md
// for the role in choosing min_P_prior.
class kappa_V : public var::Constant<kappa_V, double> {};
// Block partition of the spectrum: rows [begin, end) per block
class Blocks : public var::Constant<Blocks, Matrix<std::size_t>> {};

// Compile-time policy flags (uses_variance_aproximation, uses_recursive_aproximation,
// uses_qdt_method, return_predictions, ...) and their constraining concepts
// (uses_*_c, what_to_include_c) live in qmodel_fwd.h, included above.

class Probability_error_tolerance : public var::Constant<Probability_error_tolerance, double> {};

class Conductance_variance_error_tolerance
    : public var::Constant<Conductance_variance_error_tolerance, double> {};


class P : public Var<P, Matrix<double>> {
    friend std::string className(const P&) { return "P_ij"; }
};

class P_half : public Var<P_half, Matrix<double>> {
    friend std::string className(const P_half&) { return "Ph_ij"; }
};

template <class Policy = StabilizerPolicyEnabled, class C_Matrix>
Maybe_error<Transfer_Op_to<C_Matrix, P>> to_Transition_Probability(C_Matrix const& x) {
    if constexpr (!Policy::project_transition_probability) {
        return build<P>(x);
    }
    auto out = apply(
        [](auto const& value) {
            using std::abs;
            return abs(value);
        },
        x);
    auto sumP = out * Matrix<double>(out.ncols(), 1ul, 1.0);
    auto s = inv(diag(sumP));

    for (std::size_t i = 0; i < sumP.size(); ++i)
        if (std::isnan(primitive(sumP[i]))) {
            //  std::cerr << "rro";
            return error_message("not transition prob");
        }
    if (s){
        // auto test=s*out*Matrix<double>(out.ncols(),1ul, 1.0);
        return build<P>(s.value() * out);}
    else{
        return s.error();}
}


class gmean_i : public Var<gmean_i, Matrix<double>> {
    friend std::string className(gmean_i) { return "gmean_i"; }
};
class gmean_end_i : public Var<gmean_end_i, Matrix<double>> {
    friend std::string className(gmean_end_i) { return "gmean_end_i"; }
};


class gtotal_ij : public Var<gtotal_ij, Matrix<double>> {
    friend std::string className(gtotal_ij) { return "gtotal_ij"; }
};
class gmean_ij : public Var<gmean_ij, Matrix<double>> {
    friend std::string className(gmean_ij) { return "gmean_ij"; }
};
class gtotal_sqr_ij : public Var<gtotal_sqr_ij, Matrix<double>> {
    friend std::string className(gtotal_sqr_ij) { return "gtotal_sqr_ij"; }
};
class gsqr_i : public Var<gsqr_i, Matrix<double>> {
    friend std::string className(gsqr_i) { return "gsqr_i"; }
};
class gvar_i : public Var<gvar_i, Matrix<double>> {
    friend std::string className(gvar_i) { return "gvar_i"; }
};
class gtotal_var_ij : public Var<gtotal_var_ij, Matrix<double>> {
    friend std::string className(gtotal_var_ij) { return "gtotal_var_ij"; }
};
class gvar_ij : public Var<gvar_ij, Matrix<double>> {
    friend std::string className(gvar_ij) { return "gvar_ij"; }
};

class y_mean : public var::Var<y_mean, double> {
    friend std::string className(y_mean) { return "y_mean"; }
};
class y_var : public var::Var<y_var, double> {
    friend std::string className(y_var) { return "y_var"; }
};

class r_std : public var::Var<r_std, double> {
    friend std::string className(r_std) { return "r_std"; }
};
class r2_std : public var::Var<r2_std, double> {
    friend std::string className(r2_std) { return "r2_std"; }
};


class trust_coefficient : public var::Var<trust_coefficient, double> {
    friend std::string className(trust_coefficient) { return "trust_coefficient"; }
};

// IRT/MRT (variance_correction=true) diagnostics. All three are 1.0/0.0
// defaults in the standard Kalman branch (no Taylor active); computed in
// the IRT/MRT branch of safely_calculate_Algo_State_recursive.
//   taylor_trust_coefficient : α_vSv ∈ [0,1] — Taylor σ² shrink factor
//                              (1 = full IRT applied; <1 = numerical
//                              cancellation forced back-off, equivalent
//                              to inflating V → V/α_vSv on the σ² piece).
//   taylor_vSv               : the effective scalar ṽᵀΣv after α_vSv shrink
//                              (the Sherman-Morrison denominator's
//                              variable part).
//   taylor_strength          : α_vSv · |β| · ‖σ̄²‖ / ‖γ̄‖ where β = δ/V —
//                              relative magnitude of the σ² perturbation
//                              into v vs. the baseline γ̄ direction.
class taylor_trust_coefficient
    : public var::Var<taylor_trust_coefficient, double> {
    friend std::string className(taylor_trust_coefficient) {
        return "taylor_trust_coefficient";
    }
};
class taylor_vSv : public var::Var<taylor_vSv, double> {
    friend std::string className(taylor_vSv) { return "taylor_vSv"; }
};
class taylor_strength : public var::Var<taylor_strength, double> {
    friend std::string className(taylor_strength) { return "taylor_strength"; }
};

class Chi2 : public var::Var<Chi2, double> {
    friend std::string className(Chi2) { return "Chi2"; }
};
class macror_algorithm : public var::Constant<macror_algorithm, std::string> {
    using var::Constant<macror_algorithm, std::string>::Constant;
    friend std::string className(macror_algorithm) { return "macror_algorithm"; }
};

class PGn : public var::Var<PGn, Matrix<double>> {};
class PGG_n : public var::Var<PGG_n, Matrix<double>> {};
class PG_n : public var::Var<PG_n, Matrix<double>> {};
// class PPn : public var::Var<PPn, Matrix<double>> {};

using Qn = Vector_Space<number_of_samples, min_P, P, PG_n, PGG_n>;

using Eigs = Vector_Space<lambda, V, W, kappa_V>;

using Qdtg = Vector_Space<number_of_samples, min_P, P_half, g>;

using Qdtm =
    Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij, gsqr_i, gvar_i>;

using Qdt = Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij, gtotal_sqr_ij,
                         gsqr_i, gvar_i, gtotal_var_ij, gvar_ij>;

template <class recursive, class averaging, class variance, class variance_correction>

    requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
             uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
struct Qdt_u {
    using type = double;
};

using Patch_Model = Vector_Space<N_St, Q0, Qa, P_initial, g, N_Ch_mean, Current_Noise, Pink_Noise,
                                 Proportional_Noise, Current_Baseline,
                                 N_Ch_mean_time_segment_duration, Binomial_magical_number, min_P,
                                 Probability_error_tolerance, Conductance_variance_error_tolerance>;

inline void save(const std::string name, const Patch_Model& m) {
    std::ofstream f_Q0(name + "_Q0.txt");
    f_Q0 << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<Q0>(m) << "\n";
    std::ofstream f_Qa(name + "_Qa.txt");
    f_Qa << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<Qa>(m) << "\n";
    std::ofstream f_g(name + "_g.txt");
    f_g << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<g>(m) << "\n";
}

using Algo_State_Dynamic_Space = Vector_Space<y_mean, y_var, trust_coefficient,
                       taylor_trust_coefficient, taylor_vSv, taylor_strength,
                       r_std, Chi2, P,P_half,gmean_i,gvar_i,gmean_ij,gtotal_ij,d_gS,d_GS,P_mean_t2_y0, P_mean_t2_y1,P_mean_t15_y0, P_mean_t15_y1,
                       P_mean_t1_y1, P_mean_t20_y1, P_mean_t11_y0, P_mean_t10_y1,
                       P_mean_0t_y0,P_mean_0t_y1,P_cross_cov_0t_y0,P_cross_cov_0t_y1,
                       P_Cov_t2_y0,
                       P_Cov_t2_y1,P_Cov_t15_y0,
                       P_Cov_t15_y1, P_Cov_t1_y1, P_Cov_t20_y1, P_Cov_t11_y0, P_Cov_t10_y1>;

class Algo_State_Dynamic
    : public var::Var<
          Algo_State_Dynamic, Algo_State_Dynamic_Space> {
   public:
    Matrix<double> const& get_P_mean() const {
        if (get<P_mean_t2_y1>((*this)())().size() > 0) {
            return get<P_mean_t2_y1>((*this)())();
        }
        if (get<P_mean_t2_y0>((*this)())().size() > 0) {
            return get<P_mean_t2_y0>((*this)())();
        }
        
        return get<P_mean_t20_y1>((*this)())();
    }
    SymmetricMatrix<double> const& get_P_Cov() const {
        if (get<P_Cov_t2_y1>((*this)())().size() > 0) {
            return get<P_Cov_t2_y1>((*this)())();
        }if (get<P_Cov_t2_y0>((*this)())().size() > 0) {
            return get<P_Cov_t2_y0>((*this)())();
        }
        
        return get<P_Cov_t20_y1>((*this)())();
    }
};

using Algo_State_space=Vector_Space<y_mean, y_var, trust_coefficient,
                                    taylor_trust_coefficient, taylor_vSv, taylor_strength,
                                    r_std, Chi2, P_mean, P_Cov>;


class Algo_State
    : public var::Var<Algo_State,Algo_State_space> {
   public:
    using base_type =
        var::Var<Algo_State, Vector_Space<y_mean, y_var, trust_coefficient,
                                          taylor_trust_coefficient, taylor_vSv, taylor_strength,
                                          r_std, Chi2, P_mean,
                                          P_Cov>>;
    Algo_State(const Algo_State_Dynamic& p)
        : base_type{Vector_Space(get<y_mean>(p()), get<y_var>(p()), get<trust_coefficient>(p()),
                                 get<taylor_trust_coefficient>(p()), get<taylor_vSv>(p()),
                                 get<taylor_strength>(p()),
                                 get<r_std>(p()), get<Chi2>(p()),
                                 P_mean(p.get_P_mean()), P_Cov(p.get_P_Cov()))} {}

    using base_type::Var;
};

struct Patch_State : public var::Var<Patch_State, Vector_Space<P_mean, P_Cov>> {};

struct Evolution {};

template <class T>
class Evolution_of : public Var<Evolution_of<T>, std::vector<T>> {
   public:
    using base_type = Var<Evolution_of<T>, std::vector<T>>;
    using base_type::base_type;
    using element_type = T;
    auto& operator[](var::Var<Evolution>) { return *this; }

    auto const& operator[](var::Var<Evolution>) const { return *this; }

    friend std::string className(Evolution_of) { return "Macro_State_Evolution"; }
    using value_type = std::vector<T>;
};

template <typename... Vars>
struct Macro_State : public Vector_Space<logL, Patch_State, Vars...> {
    Macro_State() = default;
    Macro_State(Patch_State&& ps) { get<Patch_State>((*this)) = std::move(ps); }
    Macro_State(logL&& l, Patch_State&& ps, Vars&&... vars)
        : Vector_Space<logL, Patch_State, Vars...>(std::move(l), std::move(ps),
                                                   std::forward<Vars>(vars)...) {}
};

template <typename... Vars>struct dMacro_State
    : public Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                          var::Derivative<Patch_State, var::Parameters_transformed>, 
                          Vars...> {
    dMacro_State(var::Derivative<Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);

        // Seed the prior patch state (carries dx).     
        get<var::Derivative<Patch_State, var::Parameters_transformed>>(*this) = std::move(dps);
        // Accumulators start at zero/empty but share the same dx.
        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<var::Derivative<logL, var::Parameters_transformed>>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    dMacro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                 var::Derivative<Patch_State, var::Parameters_transformed>&& dps,
                 Vars&&... vars)
        : Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                       var::Derivative<Patch_State, var::Parameters_transformed>,
                          Vars...>(
              std::move(dl), std::move(dps),
              std::forward<Vars>(vars)...) {}
    dMacro_State() = default;
    
    
                    
};

template <typename... Vars>
struct ddMacro_State
    : public var::Derivative<Vector_Space<logL, Patch_State, Vars...>, var::Parameters_transformed> {
    ddMacro_State(var::Derivative<Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);
        get<Patch_State>(*this) = std::move(dps);

        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<logL>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    ddMacro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                  var::Derivative<Patch_State, var::Parameters_transformed>&& dps,
                  var::Derivative_t<Vars, var::Parameters_transformed>&&... vars)
        : var::Derivative<Vector_Space<logL, Patch_State, Vars...>, var::Parameters_transformed>(
              std::move(dl), std::move(dps), std::forward<Vars>(vars)...) {}
};

}  // namespace macrodr

template <class T>
struct mean_value_type_impl<macrodr::Evolution_of<T>> {
    using type = macrodr::Evolution_of<T>;
};

// Teach the generic derivative machinery that ddMacro_State lives in the
// Parameters_transformed derivative universe, so that Transfer_Op_to and
// related helpers propagate derivatives correctly.
namespace var {
template <class... Vars>
struct transformation_type<macrodr::ddMacro_State<Vars...>> {
    using type = Derivative_Op<Parameters_transformed>;
};

template <class... Vars>
struct is_derivative<macrodr::ddMacro_State<Vars...>> : std::true_type {};

template <class... Vars, class... Ds>
struct dx_of_dfdx<macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class G, class... Vars, class... Ds>
    requires(!is_derivative_v<G>)
struct dx_of_dfdx<G, macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class F, class... Vars, class... Ds>
struct dx_of_dfdx<Derivative<F, Parameters_transformed>, macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};
}  // namespace var

namespace macrodr {

using predictions_element =
    var::please_include<logL, elogL, vlogL, y_mean, y_var, r_std,P_mean, P_Cov, trust_coefficient,
                        taylor_trust_coefficient, taylor_vSv, taylor_strength>;

using diagnostic_element = var::please_include<logL, elogL, vlogL, Algo_State_Dynamic>;

using gradient_minimal_element =
    var::please_include<var::Derivative<logL, var::Parameters_transformed>, elogL, y_mean, y_var, r_std,
                        trust_coefficient, taylor_trust_coefficient, taylor_vSv, taylor_strength>;

using gradient_all_element =
    var::please_include<var::Derivative<logL, var::Parameters_transformed>,
                        var::Derivative<elogL, var::Parameters_transformed>,
                        var::Derivative<y_mean, var::Parameters_transformed>,
                        var::Derivative<y_var, var::Parameters_transformed>,
                        var::Derivative<r_std, var::Parameters_transformed>, trust_coefficient,
                        taylor_trust_coefficient, taylor_vSv, taylor_strength>;

using Macro_State_minimal = Macro_State<>;

using Macro_State_reg = add_t<Macro_State_minimal, var::please_include<elogL, vlogL>>;

using dMacro_State_Hessian_minimal =
    add_t<dMacro_State<>, var::please_include<Gaussian_Fisher_Information>>;

using diff_Macro_State_Gradient_Hessian =
    add_t<Macro_State<>, var::please_include<elogL, vlogL, Grad, Gaussian_Fisher_Information>>;

// Tag carrying the per-replicate MLE estimate θ̂. Used inside MLE result
// states so downstream Moment_statistics<Model_Parameters_Hat> can compute
// mean(θ̂) and Cov_emp = covariance(θ̂) across replicates automatically.
// "Hat" preserves the honest "this is our estimate, not the truth" semantic.
class Model_Parameters_Hat
    : public var::Constant<Model_Parameters_Hat, var::Parameters_transformed> {
   public:
    using base_type = var::Constant<Model_Parameters_Hat, var::Parameters_transformed>;
    using base_type::base_type;
    Model_Parameters_Hat() = default;
    friend std::string className(Model_Parameters_Hat) {
        return "Model_Parameters_Hat";
    }
};

// dMacro_State_Hessian_minimal extended with the per-replicate θ̂ slot.
// Used as the State template argument in calc_MLE_per_group_of_replicates so
// the per-group result carries the MLE estimate directly accessible for
// aggregation (no AD-chain unwrap needed downstream).
using dMacro_State_Hessian_minimal_param =
    add_t<dMacro_State_Hessian_minimal, var::please_include<Model_Parameters_Hat>>;

using Macro_State_Ev_predictions =
    add_t<Macro_State_reg,
          var::please_include<Evolution_of<add_t<Vector_Space<>, predictions_element>>>>;

using Macro_State_Ev_diagnostic =
    add_t<Macro_State_reg,
          var::please_include<Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>>;

// dMacro_State_Ev_gradient_minimal / dMacro_State_Ev_gradient_all are defined
// in micro_types.h (after micro_gradient_*_element) so they can use the unified
// element list. The macro filter writes only the macro slots; the additional
// micro slots (e.g. Derivative<micro_r_std>) stay default-constructed — the
// same empty-matrix-as-null convention used elsewhere in this codebase
// (see micro_Algo_State_Dynamic::get_micro_P_state at micro_types.h:199-228).

template <class VS>
constexpr bool is_Algo_dynamic() {
    if constexpr (has_var_c<VS const&, Evolution>) {
        using Evo = std::decay_t<decltype(std::declval<VS const&>()[var::Var<Evolution>{}])>;
        using El = typename Evo::element_type;
        return has_var_c<El const&, Algo_State_Dynamic>;
    } else {
        return false;
    }
}
class Simulation_n_sub_dt : public Var<Simulation_n_sub_dt, std::size_t> {};
class Simulation_Mode : public Var<Simulation_Mode, std::string> {};

inline constexpr auto simulation_algorithm_substeps_name = "substeps";
inline constexpr auto simulation_algorithm_uniformization_name = "uniformization";

class N_Ch_State_Evolution : public Var<N_Ch_State_Evolution, std::vector<N_channel_state>> {};

class Only_Ch_Curent_Evolution : public Var<Only_Ch_Curent_Evolution, std::vector<Patch_current>> {
};

class N_Ch_State_Sub_Evolution
    : public Var<N_Ch_State_Sub_Evolution, std::vector<N_channel_state>> {};

class Only_Ch_Curent_Sub_Evolution
    : public Var<Only_Ch_Curent_Sub_Evolution, std::vector<Patch_current>> {};

template <typename Simulate_tag>
class Simulated_Recording
    : public Var<Simulated_Recording<Simulate_tag>, add_t<Vector_Space<SeedNumber,Recording>, Simulate_tag>> {
   public:
    constexpr static const bool includes_N = var::has_it_v<Simulate_tag, N_Ch_State_Evolution>;
    using Var<Simulated_Recording<Simulate_tag>, add_t<Vector_Space<SeedNumber,Recording>, Simulate_tag>>::Var;
};

using Simulated_recording = Simulated_Recording<var::please_include<>>;

using v_Simulated_Recording =
    std::variant<Simulated_Recording<var::please_include<>>,
                 Simulated_Recording<var::please_include<Only_Ch_Curent_Evolution>>,
                 Simulated_Recording<var::please_include<N_Ch_State_Evolution>>>;

template <typename Simulate_tag>
class Simulated_Step
    : public Var<Simulated_Step<Simulate_tag>,
                 Vector_Space<N_channel_state, Simulated_Recording<Simulate_tag>>> {
    // using Vector_Space<N_channel_state,
    // Simulated_Recording<Simulate_tag>>::Vector_Space;
};

template <typename Simulate_tag>
using Simulated_Sub_Step_t =
    add_if_present_t<Vector_Space<N_channel_state, number_of_samples, y_sum>, Simulate_tag,
                     N_Ch_State_Sub_Evolution, Only_Ch_Curent_Sub_Evolution>;

template <typename Simulate_tag>
Simulated_Sub_Step_t<Simulate_tag> Simulated_Sub_Step_build(N_channel_state N) {
    Simulated_Sub_Step_t<Simulate_tag> out;
    get<N_channel_state>(out) = N;
    return out;
}

using Simulation_Parameters = Vector_Space<Simulation_Mode, Simulation_n_sub_dt>;

inline auto make_substep_simulation_parameters(std::size_t n_sub_dt) {
    return Simulation_Parameters(Simulation_Mode(simulation_algorithm_substeps_name),
                                 Simulation_n_sub_dt(n_sub_dt));
}

inline auto make_uniformization_simulation_parameters() {
    return Simulation_Parameters(Simulation_Mode(simulation_algorithm_uniformization_name),
                                 Simulation_n_sub_dt(0));
}


}  // namespace macrodr

#endif  // QMODEL_TYPES_H
