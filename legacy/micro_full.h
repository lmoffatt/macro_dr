#ifndef MICRO_FULL_H
#define MICRO_FULL_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

#include "matrix.h"
#include "maybe_error.h"
#include "multivariate_normal_distribution.h"
#include "qmodel.h"
#include "variables.h"

namespace macrodr {

// Per-microstate count of channels in each single-channel state.
// Shape: (num_full_states, k_states). Each row is an occupation vector summing to N_channels.
class Micro_state_Num_ch : public var::Constant<Micro_state_Num_ch, Matrix<std::size_t>> {
   public:
    using var::Constant<Micro_state_Num_ch, Matrix<std::size_t>>::Constant;
    friend std::string className(Micro_state_Num_ch) { return "Micro_state_Num_ch"; }
};

// Probability distribution over microstates.
// Shape: (1, num_full_states). One probability per row of Micro_state_Num_ch.
// Parallel to P_mean (which at the macro level has shape (1, k_states)).
class Micro_P_mean : public var::Var<Micro_P_mean, Matrix<double>> {
   public:
    using var::Var<Micro_P_mean, Matrix<double>>::Var;
    friend std::string className(Micro_P_mean) { return "Micro_P_mean"; }
};

using Patch_Model = Vector_Space<N_St, Q0, Qa, P_initial, g, N_Ch_mean, Current_Noise, Pink_Noise,
                                 Proportional_Noise, Current_Baseline,
                                 N_Ch_mean_time_segment_duration, Binomial_magical_number, min_P,
                                 Probability_error_tolerance, Conductance_variance_error_tolerance>;


inline std::size_t num_full_states_of(std::size_t N_channels, std::size_t k_states) {
    if (k_states == 0) {
        return 0;
    }
    std::size_t r = 1;
    for (std::size_t i = 1; i < k_states; ++i) {
        r = (r * (N_channels + i)) / i;
    }
    return r;
}

namespace micro_full_detail {
inline void fill_rows(Matrix<std::size_t>& table, std::size_t& row, std::vector<std::size_t>& buf,
                      std::size_t pos, std::size_t k_states, std::size_t remaining) {
    if (pos + 1 == k_states) {
        buf[pos] = remaining;
        for (std::size_t j = 0; j < k_states; ++j) {
            table(row, j) = buf[j];
        }
        ++row;
        return;
    }
    for (std::size_t v = 0; v <= remaining; ++v) {
        buf[pos] = v;
        fill_rows(table, row, buf, pos + 1, k_states, remaining - v);
    }
}
}  // namespace micro_full_detail

inline Micro_state_Num_ch create_Micro_state_Num_ch(std::size_t N_channels, std::size_t k_states) {
    if (k_states == 0) {
        return Micro_state_Num_ch(Matrix<std::size_t>{});
    }
    Matrix<std::size_t> table(num_full_states_of(N_channels, k_states), k_states);
    std::vector<std::size_t> buf(k_states, 0);
    std::size_t row = 0;
    micro_full_detail::fill_rows(table, row, buf, 0, k_states, N_channels);
    return Micro_state_Num_ch(std::move(table));
}

inline std::size_t N_channels_of(Micro_state_Num_ch const& micro) {
    auto const& m = micro();
    if (m.nrows() == 0) {
        return 0;
    }
    std::size_t s = 0;
    for (std::size_t j = 0; j < m.ncols(); ++j) {
        s += m(std::size_t{0}, j);
    }
    return s;
}

// Convert integer counts (num_full_states, k_states) to the probability-vector representation
// (num_full_states, k_states) of doubles, where each row is n / N_channels.
inline Matrix<double> as_probability_rows(Micro_state_Num_ch const& micro, double N_channels) {
    auto const& Ns = micro();
    Matrix<double> out(Ns.nrows(), Ns.ncols());
    double inv_N = 1.0 / N_channels;
    for (std::size_t i = 0; i < Ns.nrows(); ++i) {
        for (std::size_t j = 0; j < Ns.ncols(); ++j) {
            out(i, j) = static_cast<double>(Ns(i, j)) * inv_N;
        }
    }
    return out;
}

// State of the micro filter for a single observation step.
// Wraps the probability distribution over microstates as a Micro_P_mean.
// Parallel to Patch_State (which holds P_mean (1,k) + P_Cov (k,k)).
struct Micro_Patch_State : public var::Var<Micro_Patch_State, Vector_Space<Micro_P_mean>> {};

// Full state of the Micro_full fold, parallel to Macro_State<Vars...>.
// The micro probability distribution lives inside Micro_Patch_State; additional slots
// (y_mean, y_var, elogL, vlogL, ...) are carried as Vars... just like in Macro_State.
template <typename... Vars>
struct Micro_State : public Vector_Space<logL, Micro_Patch_State, Vars...> {
    Micro_State() = default;
    Micro_State(Micro_Patch_State&& mps) {
        get<Micro_Patch_State>((*this)) = std::move(mps);
    }
    Micro_State(logL&& l, Micro_Patch_State&& mps, Vars&&... vars)
        : Vector_Space<logL, Micro_Patch_State, Vars...>(std::move(l), std::move(mps),
                                                         std::forward<Vars>(vars)...) {}
};

template <class C_Micro_Patch_State>
auto project_Micro_to_Macro(Micro_state_Num_ch const& micro,
                            C_Micro_Patch_State const& micro_ps) {
    // Templated per-argument: when micro_ps is a Derivative-typed Micro_Patch_State,
    // the returned Patch_State carries derivatives through the matrix-level ops below.
    //
    // Shapes: prob is (1, size), P_rows is (size, k). We build the marginal mean
    // P_mean = prob · P_rows  (1, k) and the marginal covariance
    //     P_Cov = P_rowsᵀ · [diag(prob) - probᵀ · prob] · P_rows   (k, k)
    // entirely at the monolithic Derivative<Matrix> level (qmodel.h convention,
    // i.e. Derivative<Matrix<double>> rather than Matrix<Derivative<double>>).
    auto const& prob = get<Micro_P_mean>(micro_ps());
    double N = static_cast<double>(N_channels_of(micro));

    Matrix<double> P_rows = as_probability_rows(micro, N);

    auto mean = prob() * P_rows;

    // Multinomial-covariance form for one trial: M = diag(prob) − probᵀ·prob.
    // diagpos returns DiagPosDetMatrix (or Derivative<DiagPosDetMatrix>); XTX
    // returns SymPosDefMatrix (or Derivative<SymPosDefMatrix>); the subtraction
    // resolves to SymPosDefMatrix via the existing matrix-type arithmetic.
    auto M = diagpos(prob()) - XTX(prob());

    // cov_asym = P_rowsᵀ · M · P_rows — (k, k). This is mathematically symmetric
    // but typed as plain Matrix (or Derivative<Matrix>). We avoid AT_B_A here
    // because its return type (SymPosDefMatrix) doesn't line up cleanly with
    // P_Cov's underlying SymmetricMatrix<double> for the derivative path:
    // Matrix<SymPosDefMatrix> is not constructible from Matrix<SymmetricMatrix>
    // inside d_d_, so the conversion ctor chain fails. Going through a plain
    // Matrix intermediate plus X_plus_XT avoids the issue — X_plus_XT returns
    // the right type (SymmetricMatrix / Derivative<SymmetricMatrix>) in both
    // modes, and multiplying by 0.5 undoes the X + Xᵀ doubling.
    auto cov_asym = TranspMult(P_rows, M * P_rows);
    auto cov = X_plus_XT(cov_asym) * 0.5;

    Transfer_Op_to<C_Micro_Patch_State, Patch_State> out;
    get<P_mean>(out())() = std::move(mean);
    get<P_Cov>(out())() = std::move(cov);
    return out;
}

template <class C_Patch_State>
auto lift_Macro_to_Micro(Micro_state_Num_ch const& micro, C_Patch_State const& patch)
    -> Maybe_error<Transfer_Op_to<C_Patch_State, Micro_Patch_State>> {
    using std::exp;
    using std::log;
    // At init, Macro state is exactly a multinomial(N, p): P_Cov = diag(p) − ppᵀ.
    // Use the multinomial log-PMF directly — this is exact for init and handles
    // degenerate p (pure states, zero diagonals in P_Cov) without Cholesky.
    // Templated per-argument so p_mean can be Derivative<Matrix<double>, Parameters>;
    // derivatives propagate through log/exp/multiply/normalize.
    auto const& p_mean = get<P_mean>(patch());

    auto const& Ns = micro();
    std::size_t k = Ns.ncols();
    std::size_t size = Ns.nrows();
    double N = static_cast<double>(N_channels_of(micro));

    using elem_t = std::decay_t<decltype(p_mean()[std::size_t{0}])>;
    auto type_zero = p_mean()[std::size_t{0}] - p_mean()[std::size_t{0}];

    std::vector<elem_t> log_p;
    log_p.reserve(k);
    std::vector<bool> p_is_zero(k, false);
    for (std::size_t j = 0; j < k; ++j) {
        auto p_j = p_mean()[j];
        if (var::primitive(p_j) > 0.0) {
            log_p.push_back(log(p_j));
        } else {
            p_is_zero[j] = true;
            log_p.push_back(type_zero);  // sentinel; only used when n_ij == 0
        }
    }

    double log_N_fact = std::lgamma(N + 1.0);

    std::vector<elem_t> logs;
    logs.reserve(size);
    double max_log_primitive = -std::numeric_limits<double>::infinity();
    std::vector<bool> impossible(size, false);

    for (std::size_t i = 0; i < size; ++i) {
        auto lp = type_zero + log_N_fact;
        bool bad = false;
        for (std::size_t j = 0; j < k; ++j) {
            std::size_t n_ij = Ns(i, j);
            lp = lp - std::lgamma(static_cast<double>(n_ij) + 1.0);
            if (n_ij > std::size_t{0}) {
                if (p_is_zero[j]) {
                    bad = true;
                    break;
                }
                lp = lp + static_cast<double>(n_ij) * log_p[j];
            }
        }
        impossible[i] = bad;
        if (!bad) {
            double lp_prim = var::primitive(lp);
            if (lp_prim > max_log_primitive) max_log_primitive = lp_prim;
        }
        logs.push_back(std::move(lp));
    }

    Matrix<elem_t> result(1, size);
    auto sum = type_zero;
    for (std::size_t i = 0; i < size; ++i) {
        if (impossible[i]) {
            result[i] = type_zero;
        } else {
            result[i] = exp(logs[i] - max_log_primitive);
            sum = sum + result[i];
        }
    }
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = result[i] / sum;
    }

    Transfer_Op_to<C_Patch_State, Micro_Patch_State> mps;
    // result has per-element type elem_t: double for plain, Derivative<double> for derivative.
    // The Micro_P_mean field underneath is monolithic (Matrix<double> or Derivative<Matrix<double>>),
    // so in derivative mode we need to re-pack via outside_in; in plain mode result is already
    // the right type and can be moved in directly.
    if constexpr (var::is_derivative_v<elem_t>) {
        auto const& dx = var::get_dx_of_dfdx_container(result);
        get<Micro_P_mean>(mps())() = var::outside_in(result, dx);
    } else {
        get<Micro_P_mean>(mps())() = std::move(result);
    }
    return mps;
}

// -----------------------------------------------------------------------------
// k=2 averaging=2 helpers: conditional Nij moments given (n_s, n_e).
//
// For k=2 the contingency table Nij has a single free parameter N00. Given
// n_s = (j, N-j) and n_e = (m, N-m), the conditional distribution of N00
// follows Fisher's non-central hypergeometric:
//
//     P(N00 = x | n_s, n_e) ∝ C(j, x) * C(N-j, m-x) * ω^x,
//     ω = P00·P11 / (P01·P10).
//
// The three returned quantities are:
//     log_P_trans : log P(n_e | n_s) including the single-channel probabilities
//     E_N00       : conditional mean of N00 given (n_s, n_e)
//     Var_N00     : conditional variance of N00 given (n_s, n_e)
//
// All the other conditional moments E[Nij | n_s, n_e] and Cov[Nij, Nkl | n_s, n_e]
// are linear functions of (E_N00, Var_N00) because the other Nij are linear
// combinations of N00 (via the fixed row/column sums).
// -----------------------------------------------------------------------------
template <class T>
struct Fisher_moments_k2 {
    T log_P_trans;
    T E_N00;
    T Var_N00;
};

// Templated on the matrix type of P_single so the deduced scalar T flows through:
// T = double when P_single is Matrix<double>,
// T = Derivative<double, ...> when P_single is Matrix<Derivative<double, ...>>.
// Counts (j, m, N, x) and their lgamma values stay in plain double — they have zero derivative.
template <class C_P_single>
auto fisher_hypergeom_moments_k2(std::size_t j, std::size_t m, std::size_t N,
                                 C_P_single const& P_single) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(P_single(std::size_t{0}, std::size_t{0}))>;

    // Valid range for N00: [max(0, j+m-N), min(j, m)]
    std::size_t x_lo = (j + m > N) ? (j + m - N) : std::size_t{0};
    std::size_t x_hi = std::min(j, m);

    if (x_lo > x_hi) {
        // Degenerate pair: P(n_e | n_s) = 0. Build T-typed sentinels via algebraic
        // zero on an input element so derivative context is preserved.
        auto zero = P_single(std::size_t{0}, std::size_t{0}) -
                    P_single(std::size_t{0}, std::size_t{0});
        auto neg_inf = zero + (-std::numeric_limits<double>::infinity());
        return Fisher_moments_k2<T>{std::move(neg_inf), zero, zero};
    }

    // log_omega is well-defined only if all four P_single entries are > 0.
    // If any entry is exactly zero, log() gives -inf primitive and inf in
    // the derivative payload, which later becomes 0*inf = NaN through exp.
    // In that degenerate case the Fisher pmf collapses to a single x value
    // (the one where the zero entry is not raised to a positive power), so
    // we return a sentinel and let the caller treat the pair as impossible.
    auto P00 = P_single(std::size_t{0}, std::size_t{0});
    auto P01 = P_single(std::size_t{0}, std::size_t{1});
    auto P10 = P_single(std::size_t{1}, std::size_t{0});
    auto P11 = P_single(std::size_t{1}, std::size_t{1});
    if (var::primitive(P00) <= 0.0 || var::primitive(P01) <= 0.0 ||
        var::primitive(P10) <= 0.0 || var::primitive(P11) <= 0.0) {
        auto zero = P_single(std::size_t{0}, std::size_t{0}) -
                    P_single(std::size_t{0}, std::size_t{0});
        auto neg_inf = zero + (-std::numeric_limits<double>::infinity());
        return Fisher_moments_k2<T>{std::move(neg_inf), zero, zero};
    }
    auto log_P00 = log(P00);
    auto log_P01 = log(P01);
    auto log_P10 = log(P10);
    auto log_P11 = log(P11);
    auto log_omega = log_P00 + log_P11 - log_P01 - log_P10;

    double log_gamma_j_plus_1 = std::lgamma(static_cast<double>(j) + 1.0);
    double log_gamma_Nj_plus_1 = std::lgamma(static_cast<double>(N - j) + 1.0);

    std::size_t num_x = x_hi - x_lo + 1;
    std::vector<T> log_weights;
    log_weights.reserve(num_x);
    double max_log_primitive = -std::numeric_limits<double>::infinity();

    auto zero_T_lw = P_single(std::size_t{0}, std::size_t{0}) -
                     P_single(std::size_t{0}, std::size_t{0});
    for (std::size_t i = 0; i < num_x; ++i) {
        double x = static_cast<double>(x_lo + i);
        double log_C_j_x = log_gamma_j_plus_1 - std::lgamma(x + 1.0) -
                           std::lgamma(static_cast<double>(j) - x + 1.0);
        double log_C_Nj_mx = log_gamma_Nj_plus_1 -
                             std::lgamma(static_cast<double>(m) - x + 1.0) -
                             std::lgamma(static_cast<double>(N - j - m) + x + 1.0);
        T lw = zero_T_lw + (log_C_j_x + log_C_Nj_mx);
        if (x > 0.0) lw = lw + x * log_omega;
        double lw_prim = var::primitive(lw);
        if (lw_prim > max_log_primitive) {
            max_log_primitive = lw_prim;
        }
        log_weights.push_back(std::move(lw));
    }

    // Rescale with a primitive max for stability (AD treats it as a constant).
    auto zero = P_single(std::size_t{0}, std::size_t{0}) -
                P_single(std::size_t{0}, std::size_t{0});
    T sum_w = zero;
    T sum_x_w = zero;
    T sum_x2_w = zero;
    for (std::size_t i = 0; i < num_x; ++i) {
        double x = static_cast<double>(x_lo + i);
        T w = exp(log_weights[i] - max_log_primitive);
        sum_w = sum_w + w;
        sum_x_w = sum_x_w + x * w;
        sum_x2_w = sum_x2_w + (x * x) * w;
    }

    T E_N00 = sum_x_w / sum_w;
    T E_N00_sq = sum_x2_w / sum_w;
    T Var_N00 = E_N00_sq - E_N00 * E_N00;

    // log P_trans = log_Z + j·log(P01) + m·log(P10) + (N-j-m)·log(P11).
    // Guard against 0·log(0) = NaN by skipping terms whose coefficient is zero
    // (the multinomial/hypergeom convention: 0^0 = 1, so log contribution is 0).
    auto log_Z = max_log_primitive + log(sum_w);
    auto zero_T = P_single(std::size_t{0}, std::size_t{0}) -
                  P_single(std::size_t{0}, std::size_t{0});
    T const_term = zero_T;
    if (j > std::size_t{0}) const_term = const_term + static_cast<double>(j) * log_P01;
    if (m > std::size_t{0}) const_term = const_term + static_cast<double>(m) * log_P10;
    double rem = static_cast<double>(N) - static_cast<double>(j) - static_cast<double>(m);
    if (rem > 0.0) const_term = const_term + rem * log_P11;
    T log_P_trans = log_Z + const_term;

    return Fisher_moments_k2<T>{std::move(log_P_trans), std::move(E_N00), std::move(Var_N00)};
}

// -----------------------------------------------------------------------------
// Generic Micro_full step at averaging=2, for any k_states.
//
// Iterates over pairs (n_s, Nᵢⱼ) where Nᵢⱼ is a k×k transition count matrix
// with row sums matching n_s. This factorization is k-generic because each
// row of Nᵢⱼ is an independent multinomial given the row sum n_s[i] and the
// transition probabilities P_single[i, :]. The emission per (n_s, Nᵢⱼ) is
// computed via matrix operations on Nᵢⱼ and the boundary-conditional moments
// gmean_ij, gvar_ij.
//
// Posterior accumulation uses a running-max rescaling for log-sum-exp stability
// without storing per-(n_s, Nᵢⱼ) log-joints.
// -----------------------------------------------------------------------------
// Result type of a single Micro_full step. Templated on the matrix and scalar
// types that flow from the inputs — so for plain double inputs this is
// Micro_full_step_result<Matrix<double>, double>, and for derivative-tracked
// inputs it becomes Micro_full_step_result<Matrix<Derivative<double, ...>>,
// Derivative<double, ...>>.
template <class T_Matrix, class T_scalar>
struct Micro_full_step_result {
    T_Matrix new_probs;          // posterior distribution (1, num_full_states)
    T_scalar logL_contribution;  // log likelihood of this observation
};

namespace micro_full_detail {

// Log multinomial probability: log P(row | n, p).
// row is a 1×k vector of counts summing to n; log_p is 1×k with per-category log probabilities.
inline double multinomial_log_prob(Matrix<std::size_t> const& row, std::size_t n,
                                   Matrix<double> const& log_p) {
    double r = std::lgamma(static_cast<double>(n) + 1.0);
    for (std::size_t j = 0; j < row.ncols(); ++j) {
        double a = static_cast<double>(row(std::size_t{0}, j));
        r -= std::lgamma(a + 1.0);
        r += a * log_p(std::size_t{0}, j);
    }
    return r;
}

// Lex-order index of microstate n in the enumeration produced by
// create_Micro_state_Num_ch(N, k). Cost O(N · k).
inline std::size_t index_of_microstate(Matrix<std::size_t> const& n, std::size_t N,
                                       std::size_t k) {
    std::size_t idx = 0;
    std::size_t remaining = N;
    for (std::size_t i = 0; i + 1 < k; ++i) {
        std::size_t ni = n(std::size_t{0}, i);
        for (std::size_t j = 0; j < ni; ++j) {
            idx += num_full_states_of(remaining - j, k - 1 - i);
        }
        remaining -= ni;
    }
    return idx;
}

// Enumerate all k × k Nᵢⱼ matrices whose row i sums to n_s[i].
// For each valid Nᵢⱼ, invokes cb(Nᵢⱼ). Uses row_allocs[i] which enumerates
// all k-tuples summing to n_s[i].
template <class Callback>
void enumerate_Nij_rows(std::vector<Micro_state_Num_ch> const& row_allocs,
                        Matrix<std::size_t>& Nij, std::size_t row, Callback const& cb) {
    std::size_t k = Nij.ncols();
    if (row == Nij.nrows()) {
        cb(Nij);
        return;
    }
    auto const& allocs = row_allocs[row]();
    std::size_t num_allocs = allocs.nrows();
    for (std::size_t a = 0; a < num_allocs; ++a) {
        for (std::size_t j = 0; j < k; ++j) {
            Nij(row, j) = allocs(a, j);
        }
        enumerate_Nij_rows(row_allocs, Nij, row + 1, cb);
    }
}

}  // namespace micro_full_detail

template <class C_Prior_probs, class C_P_single, class C_gmean_ij, class C_gvar_ij, class C_sigma2,
          class C_y>
auto micro_full_step_avg2(C_Prior_probs const& prior_prob, C_P_single const& P_single,
                          C_gmean_ij const& gmean_ij_matrix, C_gvar_ij const& gvar_ij_matrix,
                          C_sigma2 const& sigma2_obs, Micro_state_Num_ch const& full_states,
                          C_y const& y_obs) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(log(P_single(std::size_t{0}, std::size_t{0})))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    auto const& Ns = full_states();
    std::size_t size = prior_prob.size();
    std::size_t k = Ns.ncols();
    std::size_t N = N_channels_of(full_states);

    if (P_single.nrows() != k || P_single.ncols() != k) {
        return Maybe_error<Result>(error_message("P_single shape mismatches k_states"));
    }
    if (gmean_ij_matrix.nrows() != k || gmean_ij_matrix.ncols() != k) {
        return Maybe_error<Result>(error_message("gmean_ij shape mismatches k_states"));
    }
    if (gvar_ij_matrix.nrows() != k || gvar_ij_matrix.ncols() != k) {
        return Maybe_error<Result>(error_message("gvar_ij shape mismatches k_states"));
    }

    // Zero sentinel that preserves the (possibly Derivative) type context.
    auto type_zero = P_single(std::size_t{0}, std::size_t{0}) -
                     P_single(std::size_t{0}, std::size_t{0});

    // Precompute log of single-channel transition for multinomial probs.
    // Gate log(P) by P>0: log(0) gives -inf primitive AND inf in the
    // derivative chain, which later turns into NaN through 0*inf at exp().
    Matrix<T> log_P_single(k, k);
    std::vector<bool> p_zero(k * k, false);
    auto p_zero_at = [&](std::size_t i, std::size_t j) -> std::vector<bool>::reference {
        return p_zero[i * k + j];
    };
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) {
            auto P_ij = P_single(i, j);
            if (var::primitive(P_ij) > 0.0) {
                log_P_single(i, j) = log(P_ij);
            } else {
                p_zero_at(i, j) = true;
                log_P_single(i, j) = type_zero;
            }
        }
    }

    // Accumulators with running-max rescaling for numerical stability.
    Matrix<T> post(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        post[i] = type_zero;
    }
    double max_log_primitive = -std::numeric_limits<double>::infinity();
    T total_mass = type_zero;

    Matrix<std::size_t> Nij(k, k);
    Matrix<std::size_t> n_e(1, k);
    Matrix<double> Nij_double(k, k);

    for (std::size_t s_idx = 0; s_idx < size; ++s_idx) {
        auto const& p_prior = prior_prob[s_idx];
        if (var::primitive(p_prior) <= 0.0) {
            continue;
        }
        auto log_prior = log(p_prior);

        std::vector<Micro_state_Num_ch> row_allocs;
        row_allocs.reserve(k);
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t n_s_i = Ns(s_idx, i);
            row_allocs.push_back(create_Micro_state_Num_ch(n_s_i, k));
        }

        auto callback = [&](Matrix<std::size_t> const& table) {
            T log_p_trans = type_zero;
            bool impossible = false;
            for (std::size_t i = 0; i < k && !impossible; ++i) {
                double lgam_row_double = std::lgamma(static_cast<double>(Ns(s_idx, i)) + 1.0);
                T lgam_row = type_zero + lgam_row_double;
                for (std::size_t j = 0; j < k; ++j) {
                    double a = static_cast<double>(table(i, j));
                    lgam_row = lgam_row - std::lgamma(a + 1.0);
                    if (a > 0.0) {
                        if (p_zero_at(i, j)) {
                            impossible = true;
                            break;
                        }
                        lgam_row = lgam_row + a * log_P_single(i, j);
                    }
                }
                if (!impossible) log_p_trans = log_p_trans + lgam_row;
            }
            if (impossible) return;

            for (std::size_t i = 0; i < k; ++i) {
                for (std::size_t j = 0; j < k; ++j) {
                    Nij_double(i, j) = static_cast<double>(table(i, j));
                }
            }
            auto y_mean = var::sum(elemMult(gmean_ij_matrix, Nij_double));
            auto y_var = sigma2_obs + var::sum(elemMult(gvar_ij_matrix, Nij_double));

            if (var::primitive(y_var) <= 0.0 || !std::isfinite(var::primitive(y_var))) {
                return;
            }

            auto residual = y_obs - y_mean;
            auto log_L = -0.5 * log(2.0 * std::numbers::pi * y_var) -
                         0.5 * residual * residual / y_var;
            auto log_joint = log_prior + log_p_trans + log_L;

            double log_joint_prim = var::primitive(log_joint);
            if (!std::isfinite(log_joint_prim)) {
                return;
            }

            for (std::size_t j = 0; j < k; ++j) {
                std::size_t s = 0;
                for (std::size_t i = 0; i < k; ++i) {
                    s += table(i, j);
                }
                n_e(std::size_t{0}, j) = s;
            }
            std::size_t e_idx = micro_full_detail::index_of_microstate(n_e, N, k);

            if (log_joint_prim > max_log_primitive) {
                if (std::isfinite(max_log_primitive)) {
                    double scale = std::exp(max_log_primitive - log_joint_prim);
                    for (std::size_t i = 0; i < size; ++i) {
                        post[i] = post[i] * scale;
                    }
                    total_mass = total_mass * scale;
                }
                max_log_primitive = log_joint_prim;
            }
            auto w = exp(log_joint - max_log_primitive);
            post[e_idx] = post[e_idx] + w;
            total_mass = total_mass + w;
        };

        micro_full_detail::enumerate_Nij_rows(row_allocs, Nij, 0, callback);
    }

    if (!std::isfinite(max_log_primitive) || var::primitive(total_mass) <= 0.0) {
        return Maybe_error<Result>(
            error_message("micro_full_step_avg2: no valid (n_s, Nij) with finite likelihood"));
    }

    for (std::size_t i = 0; i < size; ++i) {
        post[i] = post[i] / total_mass;
    }

    T logL_step = max_log_primitive + log(total_mass);

    return Maybe_error<Result>(Result{std::move(post), std::move(logL_step)});
}

// -----------------------------------------------------------------------------
// k=2 specialization of the Micro_full step at averaging=2.
//
// Iterates pairs (n_s, n_e) directly instead of (n_s, Nᵢⱼ) as the generic
// version does. Per pair, Fisher non-central hypergeometric moments over the
// single free parameter N00 give E[N00 | n_s, n_e], Var[N00 | n_s, n_e], and
// log P(n_e | n_s) in closed form (a single summation over min(j, m)+1 terms).
//
// Emission moments use the E[Nᵢⱼ | n_s, n_e] matrix (2×2 affine function of
// E[N00]) with var::sum(elemMult(·, ·)) operations, matching the idiom of
// Micror_stochastic's likelihood.
// -----------------------------------------------------------------------------
template <class C_Prior_probs, class C_P_single, class C_gmean_ij, class C_gvar_ij, class C_sigma2,
          class C_y>
auto micro_full_step_avg2_k2(C_Prior_probs const& prior_prob, C_P_single const& P_single,
                             C_gmean_ij const& gmean_ij_matrix, C_gvar_ij const& gvar_ij_matrix,
                             C_sigma2 const& sigma2_obs, Micro_state_Num_ch const& full_states,
                             C_y const& y_obs) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(log(P_single(std::size_t{0}, std::size_t{0})))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    auto const& Ns = full_states();
    std::size_t size = prior_prob.size();
    std::size_t k = Ns.ncols();
    std::size_t N = N_channels_of(full_states);

    if (k != 2) {
        return Maybe_error<Result>(error_message("micro_full_step_avg2_k2 requires k_states == 2"));
    }

    auto type_zero = P_single(std::size_t{0}, std::size_t{0}) -
                     P_single(std::size_t{0}, std::size_t{0});

    // Δg for the variance spread term, via matrix op over the sign pattern (+,-,-,+).
    Matrix<double> sign_matrix(std::size_t{2}, std::size_t{2});
    sign_matrix(std::size_t{0}, std::size_t{0}) = 1.0;
    sign_matrix(std::size_t{0}, std::size_t{1}) = -1.0;
    sign_matrix(std::size_t{1}, std::size_t{0}) = -1.0;
    sign_matrix(std::size_t{1}, std::size_t{1}) = 1.0;
    auto delta_g = var::sum(elemMult(gmean_ij_matrix, sign_matrix));

    Matrix<T> post(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        post[i] = type_zero;
    }
    double max_log_primitive = -std::numeric_limits<double>::infinity();
    T total_mass = type_zero;

    Matrix<T> E_Nij(std::size_t{2}, std::size_t{2});
    double N_d = static_cast<double>(N);

    for (std::size_t s_idx = 0; s_idx < size; ++s_idx) {
        auto const& p_prior = prior_prob[s_idx];
        if (var::primitive(p_prior) <= 0.0) {
            continue;
        }
        auto log_prior = log(p_prior);
        std::size_t j = Ns(s_idx, std::size_t{0});
        double j_d = static_cast<double>(j);

        for (std::size_t e_idx = 0; e_idx < size; ++e_idx) {
            std::size_t m = Ns(e_idx, std::size_t{0});
            double m_d = static_cast<double>(m);

            auto fm = fisher_hypergeom_moments_k2(j, m, N, P_single);
            if (!std::isfinite(var::primitive(fm.log_P_trans))) {
                continue;
            }

            E_Nij(std::size_t{0}, std::size_t{0}) = fm.E_N00;
            E_Nij(std::size_t{0}, std::size_t{1}) = j_d - fm.E_N00;
            E_Nij(std::size_t{1}, std::size_t{0}) = m_d - fm.E_N00;
            E_Nij(std::size_t{1}, std::size_t{1}) = N_d - j_d - m_d + fm.E_N00;

            auto y_mean = var::sum(elemMult(gmean_ij_matrix, E_Nij));
            auto y_var = sigma2_obs + var::sum(elemMult(gvar_ij_matrix, E_Nij)) +
                         delta_g * delta_g * fm.Var_N00;

            if (var::primitive(y_var) <= 0.0 || !std::isfinite(var::primitive(y_var))) {
                continue;
            }

            auto residual = y_obs - y_mean;
            auto log_L = -0.5 * log(2.0 * std::numbers::pi * y_var) -
                         0.5 * residual * residual / y_var;
            auto log_joint = log_prior + fm.log_P_trans + log_L;

            double log_joint_prim = var::primitive(log_joint);
            if (!std::isfinite(log_joint_prim)) {
                continue;
            }

            if (log_joint_prim > max_log_primitive) {
                if (std::isfinite(max_log_primitive)) {
                    double scale = std::exp(max_log_primitive - log_joint_prim);
                    for (std::size_t i = 0; i < size; ++i) {
                        post[i] = post[i] * scale;
                    }
                    total_mass = total_mass * scale;
                }
                max_log_primitive = log_joint_prim;
            }
            auto w = exp(log_joint - max_log_primitive);
            post[e_idx] = post[e_idx] + w;
            total_mass = total_mass + w;
        }
    }

    if (!std::isfinite(max_log_primitive) || var::primitive(total_mass) <= 0.0) {
        return Maybe_error<Result>(
            error_message("micro_full_step_avg2_k2: no valid pairs with finite likelihood"));
    }

    for (std::size_t i = 0; i < size; ++i) {
        post[i] = post[i] / total_mass;
    }

    T logL_step = max_log_primitive + log(total_mass);

    return Maybe_error<Result>(Result{std::move(post), std::move(logL_step)});
}

// -----------------------------------------------------------------------------
// Multinomial propagation kernel: applies a single-channel transition matrix
// P_trans (k×k) to a distribution P(n_s) over microstates, producing P(n_e)
// where n_e is distributed according to the multinomial convolution of per-row
// transitions. k-generic via enumerate_Nij_rows.
//
// Used by avg=0 (applied twice with P_half) and avg=1 (applied once with P).
// -----------------------------------------------------------------------------
template <class C_P_trans, class C_Prior_probs>
auto propagate_multinomial(Micro_state_Num_ch const& full_states, C_P_trans const& P_trans,
                           C_Prior_probs const& prior_prob) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(log(P_trans(std::size_t{0}, std::size_t{0})))>;

    auto const& Ns = full_states();
    std::size_t size = prior_prob.size();
    std::size_t k = Ns.ncols();
    std::size_t N = N_channels_of(full_states);

    auto type_zero = P_trans(std::size_t{0}, std::size_t{0}) -
                     P_trans(std::size_t{0}, std::size_t{0});

    // Mirror lift_Macro_to_Micro's pattern: gate log(P) by P>0 to keep the
    // derivative payload finite. log(0) gives -inf primitive AND inf in the
    // derivative chain (1/0); a later exp() would underflow the primitive to
    // 0 while the derivative becomes 0*inf = NaN. Tables that would multiply
    // zero-P entries by a>0 are physically impossible (probability zero), so
    // we skip them entirely.
    Matrix<T> log_P_trans(k, k);
    std::vector<bool> p_zero(k * k, false);
    auto p_zero_at = [&](std::size_t i, std::size_t j) -> std::vector<bool>::reference {
        return p_zero[i * k + j];
    };
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) {
            auto P_ij = P_trans(i, j);
            if (var::primitive(P_ij) > 0.0) {
                log_P_trans(i, j) = log(P_ij);
            } else {
                p_zero_at(i, j) = true;
                log_P_trans(i, j) = type_zero;  // sentinel; only consumed when a==0
            }
        }
    }

    Matrix<T> out(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        out[i] = type_zero;
    }

    Matrix<std::size_t> Nij(k, k);
    Matrix<std::size_t> n_e(1, k);

    for (std::size_t s_idx = 0; s_idx < size; ++s_idx) {
        auto const& p_s = prior_prob[s_idx];
        if (var::primitive(p_s) <= 0.0) {
            continue;
        }

        std::vector<Micro_state_Num_ch> row_allocs;
        row_allocs.reserve(k);
        for (std::size_t i = 0; i < k; ++i) {
            row_allocs.push_back(create_Micro_state_Num_ch(Ns(s_idx, i), k));
        }

        auto callback = [&](Matrix<std::size_t> const& table) {
            T log_p_trans = type_zero;
            bool impossible = false;
            for (std::size_t i = 0; i < k && !impossible; ++i) {
                double lgam_row_double = std::lgamma(static_cast<double>(Ns(s_idx, i)) + 1.0);
                T lgam_row = type_zero + lgam_row_double;
                for (std::size_t j = 0; j < k; ++j) {
                    double a = static_cast<double>(table(i, j));
                    lgam_row = lgam_row - std::lgamma(a + 1.0);
                    if (a > 0.0) {
                        if (p_zero_at(i, j)) {
                            impossible = true;
                            break;
                        }
                        lgam_row = lgam_row + a * log_P_trans(i, j);
                    }
                }
                if (!impossible) log_p_trans = log_p_trans + lgam_row;
            }
            if (impossible) return;

            for (std::size_t j = 0; j < k; ++j) {
                std::size_t s = 0;
                for (std::size_t i = 0; i < k; ++i) {
                    s += table(i, j);
                }
                n_e(std::size_t{0}, j) = s;
            }
            std::size_t e_idx = micro_full_detail::index_of_microstate(n_e, N, k);
            out[e_idx] = out[e_idx] + p_s * exp(log_p_trans);
        };

        micro_full_detail::enumerate_Nij_rows(row_allocs, Nij, 0, callback);
    }

    return out;
}

// -----------------------------------------------------------------------------
// Micro_full step at averaging=0 (pointwise current, "g-style" — parallel to Qdtg).
// Sampling is at the midpoint of the Δt interval, so the structure is:
//   half-propagate → emit → update → half-propagate.
// -----------------------------------------------------------------------------
template <class C_Prior_probs, class C_P_half, class C_g, class C_sigma2, class C_y>
auto micro_full_step_avg0(C_Prior_probs const& prior_prob, C_P_half const& P_half,
                          C_g const& g_vector, C_sigma2 const& sigma2_obs,
                          Micro_state_Num_ch const& full_states, C_y const& y_obs) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(log(P_half(std::size_t{0}, std::size_t{0})))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    auto const& Ns = full_states();
    std::size_t size = prior_prob.size();
    std::size_t k = Ns.ncols();

    if (P_half.nrows() != k || P_half.ncols() != k) {
        return Maybe_error<Result>(error_message("P_half shape mismatches k_states"));
    }
    if (g_vector.size() != k) {
        return Maybe_error<Result>(error_message("g_vector size mismatches k_states"));
    }
    if (var::primitive(sigma2_obs) <= 0.0) {
        return Maybe_error<Result>(error_message("sigma2_obs must be positive"));
    }

    auto type_zero = P_half(std::size_t{0}, std::size_t{0}) -
                     P_half(std::size_t{0}, std::size_t{0});

    // Step 1: half-propagate prior to the midpoint.
    auto mid = propagate_multinomial(full_states, P_half, prior_prob);

    // Step 2: per-microstate pointwise emission L[s] = Normal(y; g·n_s, σ²).
    double max_log_primitive = -std::numeric_limits<double>::infinity();
    T total_mass = type_zero;
    Matrix<T> mid_post(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        mid_post[i] = type_zero;
    }

    for (std::size_t s_idx = 0; s_idx < size; ++s_idx) {
        auto const& p_mid = mid[s_idx];
        if (var::primitive(p_mid) <= 0.0) {
            continue;
        }
        T y_mean = type_zero;
        for (std::size_t j = 0; j < k; ++j) {
            y_mean = y_mean + static_cast<double>(Ns(s_idx, j)) * g_vector[j];
        }
        auto residual = y_obs - y_mean;
        auto log_L = -0.5 * log(2.0 * std::numbers::pi * sigma2_obs) -
                     0.5 * residual * residual / sigma2_obs;
        auto log_joint = log(p_mid) + log_L;

        double log_joint_prim = var::primitive(log_joint);
        if (!std::isfinite(log_joint_prim)) {
            continue;
        }

        if (log_joint_prim > max_log_primitive) {
            if (std::isfinite(max_log_primitive)) {
                double scale = std::exp(max_log_primitive - log_joint_prim);
                for (std::size_t i = 0; i < size; ++i) {
                    mid_post[i] = mid_post[i] * scale;
                }
                total_mass = total_mass * scale;
            }
            max_log_primitive = log_joint_prim;
        }
        auto w = exp(log_joint - max_log_primitive);
        mid_post[s_idx] = mid_post[s_idx] + w;
        total_mass = total_mass + w;
    }

    if (!std::isfinite(max_log_primitive) || var::primitive(total_mass) <= 0.0) {
        return Maybe_error<Result>(
            error_message("micro_full_step_avg0: no valid microstate with finite likelihood"));
    }

    for (std::size_t i = 0; i < size; ++i) {
        mid_post[i] = mid_post[i] / total_mass;
    }

    T logL_step = max_log_primitive + log(total_mass);

    // Step 3: half-propagate the post-update distribution to the end of the interval.
    auto post_end = propagate_multinomial(full_states, P_half, mid_post);

    return Maybe_error<Result>(Result{std::move(post_end), std::move(logL_step)});
}

// -----------------------------------------------------------------------------
// Micro_full step at averaging=1 (start-conditional integrated current, "m-style" — parallel to Qdtm).
// Emission uses gmean_i, gvar_i; update on start state, then propagate once with P.
// -----------------------------------------------------------------------------
template <class C_Prior_probs, class C_P_single, class C_gmean_i, class C_gvar_i, class C_sigma2,
          class C_y>
auto micro_full_step_avg1(C_Prior_probs const& prior_prob, C_P_single const& P_single,
                          C_gmean_i const& gmean_i_vector, C_gvar_i const& gvar_i_vector,
                          C_sigma2 const& sigma2_obs, Micro_state_Num_ch const& full_states,
                          C_y const& y_obs) {
    using std::exp;
    using std::log;
    using T = std::decay_t<decltype(log(P_single(std::size_t{0}, std::size_t{0})))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    auto const& Ns = full_states();
    std::size_t size = prior_prob.size();
    std::size_t k = Ns.ncols();

    if (P_single.nrows() != k || P_single.ncols() != k) {
        return Maybe_error<Result>(error_message("P_single shape mismatches k_states"));
    }
    if (gmean_i_vector.size() != k || gvar_i_vector.size() != k) {
        return Maybe_error<Result>(error_message("gmean_i/gvar_i size mismatches k_states"));
    }

    auto type_zero = P_single(std::size_t{0}, std::size_t{0}) -
                     P_single(std::size_t{0}, std::size_t{0});

    // Step 1: per-microstate start-conditional emission.
    double max_log_primitive = -std::numeric_limits<double>::infinity();
    T total_mass = type_zero;
    Matrix<T> start_post(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        start_post[i] = type_zero;
    }

    for (std::size_t s_idx = 0; s_idx < size; ++s_idx) {
        auto const& p_s = prior_prob[s_idx];
        if (var::primitive(p_s) <= 0.0) {
            continue;
        }
        T y_mean = type_zero;
        auto y_var = sigma2_obs;  // type deduces from sigma2_obs (could be T or scalar)
        for (std::size_t j = 0; j < k; ++j) {
            double n_sj = static_cast<double>(Ns(s_idx, j));
            y_mean = y_mean + n_sj * gmean_i_vector[j];
            y_var = y_var + n_sj * gvar_i_vector[j];
        }

        if (var::primitive(y_var) <= 0.0 || !std::isfinite(var::primitive(y_var))) {
            continue;
        }

        auto residual = y_obs - y_mean;
        auto log_L = -0.5 * log(2.0 * std::numbers::pi * y_var) -
                     0.5 * residual * residual / y_var;
        auto log_joint = log(p_s) + log_L;

        double log_joint_prim = var::primitive(log_joint);
        if (!std::isfinite(log_joint_prim)) {
            continue;
        }

        if (log_joint_prim > max_log_primitive) {
            if (std::isfinite(max_log_primitive)) {
                double scale = std::exp(max_log_primitive - log_joint_prim);
                for (std::size_t i = 0; i < size; ++i) {
                    start_post[i] = start_post[i] * scale;
                }
                total_mass = total_mass * scale;
            }
            max_log_primitive = log_joint_prim;
        }
        auto w = exp(log_joint - max_log_primitive);
        start_post[s_idx] = start_post[s_idx] + w;
        total_mass = total_mass + w;
    }

    if (!std::isfinite(max_log_primitive) || var::primitive(total_mass) <= 0.0) {
        return Maybe_error<Result>(
            error_message("micro_full_step_avg1: no valid microstate with finite likelihood"));
    }

    for (std::size_t i = 0; i < size; ++i) {
        start_post[i] = start_post[i] / total_mass;
    }

    T logL_step = max_log_primitive + log(total_mass);

    // Step 2: propagate the updated start distribution to end of interval.
    auto post_end = propagate_multinomial(full_states, P_single, start_post);

    return Maybe_error<Result>(Result{std::move(post_end), std::move(logL_step)});
}

// -----------------------------------------------------------------------------
// Tag struct paralelo a MacroR2 (legacy/qmodel.h:6383-6414).
// Despacha por averaging::value al step function correspondiente.
//
// El caller (log_Likelihood_micro) es responsable de pasarle los argumentos
// correctos para cada averaging:
//   averaging=0 → (prior, P_half, g_vector, sigma2, full_states, y_obs)
//   averaging=1 → (prior, P_single, gmean_i, gvar_i, sigma2, full_states, y_obs)
//   averaging=2 → (prior, P_single, gmean_ij, gvar_ij, sigma2, full_states, y_obs)
//
// Nota: la forma exacta del concept constraint y los template parameters
// se alinea con MacroR2 cuando agreguemos las flags `uses_*_aproximation`
// compartidas. Por ahora, requerir solo lo que está definido.
// -----------------------------------------------------------------------------
template <class recursive, class averaging, class variance>
struct MicroR2 {
    friend std::string ToString(MicroR2) {
        std::string out = "MicroR";
        if (recursive::value) {
            out += "_R";
        } else {
            out += "_NR";
        }
        if (averaging::value == 2) {
            out += "_2";
        } else if (averaging::value == 1) {
            out += "_1";
        } else {
            out += "_0";
        }
        if (variance::value) {
            out += "_V";
        } else {
            out += "_M";
        }
        return out;
    }

    template <class... Ts>
    auto operator()(Ts&&... args) {
        if constexpr (averaging::value == 0) {
            return micro_full_step_avg0(std::forward<Ts>(args)...);
        } else if constexpr (averaging::value == 1) {
            return micro_full_step_avg1(std::forward<Ts>(args)...);
        } else {
            return micro_full_step_avg2(std::forward<Ts>(args)...);
        }
    }
};

// =============================================================================
// Lifted micro path
//
// Strategy: build a Patch_Model at the microstate granularity (M = number of
// microstates, single-channel per "patch") and reuse the macro Calc_Qdt_step
// memoization on it. The per-step update is plain Bayes on the joint
// (s_start, s_end) microstate space — direct multiplication of prior,
// transition probability, and Gaussian emission likelihood. No multinomial
// enumeration, no log of any transition probability, no NaN injection in the
// derivative chain.
// =============================================================================

// build_Q_micro: assemble the (M × M) microstate generator from the single-
// channel rate matrix. M = num_full_states. For each microstate s with count
// vector n = (n_1, ..., n_k), every off-diagonal entry corresponds to a one-
// channel transition (i → j) with rate n_i · Q_single[i, j] taking n →
// n - e_i + e_j. Diagonal entries are the negative of the row sum.
template <class C_Q_single>
auto build_Q_micro(C_Q_single const& Q_single, Micro_state_Num_ch const& full_states) {
    using elem_t =
        std::decay_t<decltype(Q_single(std::size_t{0}, std::size_t{0}))>;
    auto const& Ns = full_states();
    std::size_t M = Ns.nrows();
    std::size_t k = Ns.ncols();
    std::size_t N = N_channels_of(full_states);

    auto type_zero = Q_single(std::size_t{0}, std::size_t{0}) -
                     Q_single(std::size_t{0}, std::size_t{0});

    Matrix<elem_t> Q_micro(M, M);
    for (std::size_t s = 0; s < M; ++s)
        for (std::size_t s2 = 0; s2 < M; ++s2) Q_micro(s, s2) = type_zero;

    Matrix<std::size_t> n_prime(1, k);
    for (std::size_t s = 0; s < M; ++s) {
        elem_t row_sum = type_zero;
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t n_s_i = Ns(s, i);
            if (n_s_i == 0) continue;
            for (std::size_t j = 0; j < k; ++j) {
                if (i == j) continue;
                for (std::size_t a = 0; a < k; ++a) n_prime(std::size_t{0}, a) = Ns(s, a);
                n_prime(std::size_t{0}, i) -= 1;
                n_prime(std::size_t{0}, j) += 1;
                std::size_t s_prime =
                    micro_full_detail::index_of_microstate(n_prime, N, k);
                auto rate = static_cast<double>(n_s_i) * Q_single(i, j);
                Q_micro(s, s_prime) = Q_micro(s, s_prime) + rate;
                row_sum = row_sum + rate;
            }
        }
        Q_micro(s, s) = type_zero - row_sum;
    }

    if constexpr (var::is_derivative_v<elem_t>) {
        // Per-element Matrix<Derivative<double>> → monolithic
        // Derivative<Matrix<double>> for downstream matrix-level ops.
        auto const& dx = var::get_dx_of_dfdx_container(Q_micro);
        return var::outside_in(Q_micro, dx);
    } else {
        return Q_micro;
    }
}

// lift_Patch_Model_to_Micro — assemble Patch_Model at microstate granularity.
//   Q0_micro / Qa_micro: lifted via build_Q_micro from the macro single-channel
//                        rate matrices.
//   g_micro:             full_states · g_macro (ensemble integrated current
//                        per microstate).
//   N_Ch_mean = 1:       the patch IS one microstate-channel.
//   P_initial:           uniform 1/M placeholder. Macro_DMR::init(m_micro) is
//                        never called; the micro fold computes the initial
//                        distribution separately via lift_Macro_to_Micro on the
//                        macro init Patch_State.
//   All other fields:    carried over from m_macro unchanged.
template <class C_Patch_Model>
auto lift_Patch_Model_to_Micro(C_Patch_Model const& m_macro,
                                Micro_state_Num_ch const& full_states) {
    auto Q0_micro_v = build_Q_micro(get<Q0>(m_macro)(), full_states);
    auto Qa_micro_v = build_Q_micro(get<Qa>(m_macro)(), full_states);

    std::size_t M = full_states().nrows();
    std::size_t k = full_states().ncols();
    Matrix<double> Ns_double(M, k);
    for (std::size_t s = 0; s < M; ++s)
        for (std::size_t j = 0; j < k; ++j)
            Ns_double(s, j) = static_cast<double>(full_states()(s, j));
    auto g_micro_v = Ns_double * get<g>(m_macro)();

    Matrix<double> P_initial_micro_v(1, M, 1.0 / static_cast<double>(M));
    auto N_ch_macro_v = get<N_Ch_mean>(m_macro)();
    Matrix<double> N_Ch_mean_micro_v(N_ch_macro_v.nrows(), N_ch_macro_v.ncols(), 1.0);

    return build<Patch_Model>(
        N_St(M),
        build<Q0>(std::move(Q0_micro_v)),
        build<Qa>(std::move(Qa_micro_v)),
        build<P_initial>(std::move(P_initial_micro_v)),
        build<g>(std::move(g_micro_v)),
        build<N_Ch_mean>(std::move(N_Ch_mean_micro_v)),
        get<Current_Noise>(m_macro),
        get<Pink_Noise>(m_macro),
        get<Proportional_Noise>(m_macro),
        get<Current_Baseline>(m_macro),
        get<N_Ch_mean_time_segment_duration>(m_macro),
        get<Binomial_magical_number>(m_macro),
        get<min_P>(m_macro),
        get<Probability_error_tolerance>(m_macro),
        get<Conductance_variance_error_tolerance>(m_macro));
}

// -----------------------------------------------------------------------------
// Per-step micro-fold inner functions for the lifted path.
//
// Each takes the micro-level matrices/vectors directly (from macro
// Calc_Qdt_step on the lifted Patch_Model_micro) and updates the posterior
// over microstates by direct probability multiplication. One log() per step
// at the end (logL accumulator); no per-cell logarithms.
//
// Inputs/outputs are matrix types (Matrix<double> or
// Derivative<Matrix<double>>); scalar indexing on monolithic Derivative<Matrix>
// uses the (i, j) operator which returns Derivative<double>.
// -----------------------------------------------------------------------------

// avg=2 (micro_IR): plain Bayes on (s_start, s_end) joint space.
//   joint(s, s') = π(s) · P_micro[s, s'] · N(y_obs; gmean_ij[s, s'], gvar_ij[s, s'] + σ²_obs)
//   π'(s')       = (Σ_s joint(s, s')) / Σ_{s,s'} joint
//   logL_step    = log(Σ_{s,s'} joint)
//
// At high N_ch (M = N+1 large) and small σ² the Gaussian densities can underflow
// for every (s, s') pair if `y_obs` is more than a few σ from any predicted mean.
// Anchor: find max log-density across pairs (primitive only), then evaluate
// exp(log_L − max_log_L) — that ratio is in [0, 1] for the dominant pair, so
// no underflow on the leading contributor. π and P_micro stay in linear space:
// no log of any transition or prior probability, no log(0) NaN hazard.
//
// Internal representation: per-element Matrix<T> where T is the deduced scalar
// type — Matrix<double> in plain mode, Matrix<Derivative<double>> in derivative
// mode. Caller (fold body) re-packs to monolithic Derivative<Matrix<double>>
// via outside_in for the next iteration's prior_prob, identity-move otherwise.
template <class C_Prior_probs, class C_P, class C_gmean_ij, class C_gvar_ij, class C_sigma2,
          class C_y>
auto micro_full_step_avg2_lifted(C_Prior_probs const& prior_prob, C_P const& P_micro,
                                  C_gmean_ij const& gmean_ij_m, C_gvar_ij const& gvar_ij_m,
                                  C_sigma2 const& sigma2_obs, C_y const& y_obs) {
    using std::exp;
    using std::log;
    using std::sqrt;
    using T = std::decay_t<decltype(sigma2_obs *
                                     P_micro(std::size_t{0}, std::size_t{0}))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    if (var::primitive(sigma2_obs) <= 0.0)
        return Maybe_error<Result>(error_message("sigma2_obs must be positive"));

    auto type_zero = sigma2_obs - sigma2_obs;
    std::size_t M = prior_prob.ncols();
    auto two_pi = 2.0 * std::numbers::pi;
    double sigma2_p = var::primitive(sigma2_obs);

    // Pass 1: max log-density (primitive only) across (s, s') pairs with positive
    // prior · P_micro. Anchors the rescaling in pass 2.
    double y_obs_p = var::primitive(y_obs);
    double max_log_L = -std::numeric_limits<double>::infinity();
    for (std::size_t s = 0; s < M; ++s) {
        if (var::primitive(prior_prob[s]) <= 0.0) continue;
        for (std::size_t s2 = 0; s2 < M; ++s2) {
            if (var::primitive(P_micro(s, s2)) <= 0.0) continue;
            double y_var_p = sigma2_p + var::primitive(gvar_ij_m)(s, s2);
            if (y_var_p <= 0.0 || !std::isfinite(y_var_p)) continue;
            double residual_p = y_obs_p - var::primitive(gmean_ij_m)(s, s2);
            double log_L = -0.5 * std::log(two_pi * y_var_p) -
                           0.5 * residual_p * residual_p / y_var_p;
            if (log_L > max_log_L) max_log_L = log_L;
        }
    }
    if (!std::isfinite(max_log_L))
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg2_lifted: no valid (s, s') pair with finite likelihood"));

    Matrix<T> post(1, M);
    for (std::size_t s = 0; s < M; ++s) post[s] = type_zero;
    T total_mass = type_zero;

    for (std::size_t s = 0; s < M; ++s) {
        if (var::primitive(prior_prob[s]) <= 0.0) continue;
        for (std::size_t s2 = 0; s2 < M; ++s2) {
            if (var::primitive(P_micro(s, s2)) <= 0.0) continue;
            auto y_var_pair = sigma2_obs + gvar_ij_m(s, s2);
            if (var::primitive(y_var_pair) <= 0.0 ||
                !std::isfinite(var::primitive(y_var_pair))) {
                continue;
            }
            auto residual = y_obs - gmean_ij_m(s, s2);
            auto log_L = -0.5 * log(two_pi * y_var_pair) -
                         0.5 * residual * residual / y_var_pair;
            auto rescaled = exp(log_L - max_log_L);
            auto w = prior_prob[s] * P_micro(s, s2) * rescaled;
            post[s2] = post[s2] + w;
            total_mass = total_mass + w;
        }
    }

    if (var::primitive(total_mass) <= 0.0)
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg2_lifted: total emission mass is zero"));

    for (std::size_t s = 0; s < M; ++s) post[s] = post[s] / total_mass;
    auto logL_step = max_log_L + log(total_mass);
    return Maybe_error<Result>(Result{std::move(post), std::move(logL_step)});
}

// avg=1 (micro_MR): emission depends only on the start state. Same
// max-anchored Gaussian rescaling as avg=2 to avoid underflow at high N.
template <class C_Prior_probs, class C_P, class C_gmean_i, class C_gvar_i, class C_sigma2,
          class C_y>
auto micro_full_step_avg1_lifted(C_Prior_probs const& prior_prob, C_P const& P_micro,
                                  C_gmean_i const& gmean_i_m, C_gvar_i const& gvar_i_m,
                                  C_sigma2 const& sigma2_obs, C_y const& y_obs) {
    using std::exp;
    using std::log;
    using std::sqrt;
    using T = std::decay_t<decltype(sigma2_obs *
                                     P_micro(std::size_t{0}, std::size_t{0}))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    if (var::primitive(sigma2_obs) <= 0.0)
        return Maybe_error<Result>(error_message("sigma2_obs must be positive"));

    auto type_zero = sigma2_obs - sigma2_obs;
    std::size_t M = prior_prob.ncols();
    auto two_pi = 2.0 * std::numbers::pi;
    double sigma2_p = var::primitive(sigma2_obs);

    double y_obs_p = var::primitive(y_obs);
    double max_log_L = -std::numeric_limits<double>::infinity();
    for (std::size_t s = 0; s < M; ++s) {
        if (var::primitive(prior_prob[s]) <= 0.0) continue;
        double y_var_p = sigma2_p + var::primitive(gvar_i_m(s, std::size_t{0}));
        if (y_var_p <= 0.0 || !std::isfinite(y_var_p)) continue;
        double residual_p = y_obs_p - var::primitive(gmean_i_m(s, std::size_t{0}));
        double log_L = -0.5 * std::log(two_pi * y_var_p) -
                       0.5 * residual_p * residual_p / y_var_p;
        if (log_L > max_log_L) max_log_L = log_L;
    }
    if (!std::isfinite(max_log_L))
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg1_lifted: no valid microstate with finite likelihood"));

    Matrix<T> post(1, M);
    for (std::size_t s = 0; s < M; ++s) post[s] = type_zero;
    T total_mass = type_zero;

    for (std::size_t s = 0; s < M; ++s) {
        if (var::primitive(prior_prob[s]) <= 0.0) continue;
        auto y_var_s = sigma2_obs + gvar_i_m(s, std::size_t{0});
        if (var::primitive(y_var_s) <= 0.0 || !std::isfinite(var::primitive(y_var_s))) {
            continue;
        }
        auto residual = y_obs - gmean_i_m(s, std::size_t{0});
        auto log_L_s = -0.5 * log(two_pi * y_var_s) -
                       0.5 * residual * residual / y_var_s;
        auto L_s = exp(log_L_s - max_log_L);
        auto v_s = prior_prob[s] * L_s;
        for (std::size_t s2 = 0; s2 < M; ++s2) {
            if (var::primitive(P_micro(s, s2)) <= 0.0) continue;
            auto w = v_s * P_micro(s, s2);
            post[s2] = post[s2] + w;
            total_mass = total_mass + w;
        }
    }

    if (var::primitive(total_mass) <= 0.0)
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg1_lifted: total emission mass is zero"));

    for (std::size_t s = 0; s < M; ++s) post[s] = post[s] / total_mass;
    auto logL_step = max_log_L + log(total_mass);
    return Maybe_error<Result>(Result{std::move(post), std::move(logL_step)});
}

// avg=0 (micro_R): midpoint sample. Two half-steps; emission depends on the
//   midpoint microstate's instantaneous current g_micro[s_mid]. Same
//   max-anchored Gaussian rescaling as avg=1 / avg=2.
template <class C_Prior_probs, class C_P_half, class C_g_micro, class C_sigma2, class C_y>
auto micro_full_step_avg0_lifted(C_Prior_probs const& prior_prob, C_P_half const& P_half_m,
                                  C_g_micro const& g_micro_m, C_sigma2 const& sigma2_obs,
                                  C_y const& y_obs) {
    using std::exp;
    using std::log;
    using std::sqrt;
    using T = std::decay_t<decltype(sigma2_obs *
                                     P_half_m(std::size_t{0}, std::size_t{0}))>;
    using Result = Micro_full_step_result<Matrix<T>, T>;

    if (var::primitive(sigma2_obs) <= 0.0)
        return Maybe_error<Result>(error_message("sigma2_obs must be positive"));

    auto type_zero = sigma2_obs - sigma2_obs;
    std::size_t M = prior_prob.ncols();
    auto two_pi = 2.0 * std::numbers::pi;
    double sigma2_p = var::primitive(sigma2_obs);

    // First half-step: prior → mid_prob (per-element matrix-vector product).
    Matrix<T> mid_prob(1, M);
    for (std::size_t s_mid = 0; s_mid < M; ++s_mid) mid_prob[s_mid] = type_zero;
    for (std::size_t s_start = 0; s_start < M; ++s_start) {
        if (var::primitive(prior_prob[s_start]) <= 0.0) continue;
        for (std::size_t s_mid = 0; s_mid < M; ++s_mid) {
            if (var::primitive(P_half_m(s_start, s_mid)) <= 0.0) continue;
            mid_prob[s_mid] =
                mid_prob[s_mid] + prior_prob[s_start] * P_half_m(s_start, s_mid);
        }
    }

    // Pass 1: max log-density (primitive only) across midpoint microstates.
    double y_obs_p = var::primitive(y_obs);
    double max_log_L = -std::numeric_limits<double>::infinity();
    for (std::size_t s_mid = 0; s_mid < M; ++s_mid) {
        if (var::primitive(mid_prob[s_mid]) <= 0.0) continue;
        double residual_p = y_obs_p - var::primitive(g_micro_m(s_mid, std::size_t{0}));
        double log_L = -0.5 * std::log(two_pi * sigma2_p) -
                       0.5 * residual_p * residual_p / sigma2_p;
        if (log_L > max_log_L) max_log_L = log_L;
    }
    if (!std::isfinite(max_log_L))
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg0_lifted: no valid microstate with finite likelihood"));

    // Pass 2: emission at the midpoint microstate, with rescaled Gaussian.
    Matrix<T> mid_post(1, M);
    for (std::size_t s = 0; s < M; ++s) mid_post[s] = type_zero;
    T total_mass = type_zero;
    for (std::size_t s_mid = 0; s_mid < M; ++s_mid) {
        if (var::primitive(mid_prob[s_mid]) <= 0.0) continue;
        auto residual = y_obs - g_micro_m(s_mid, std::size_t{0});
        auto log_L_s = -0.5 * log(two_pi * sigma2_obs) -
                       0.5 * residual * residual / sigma2_obs;
        auto L_s = exp(log_L_s - max_log_L);
        auto w = mid_prob[s_mid] * L_s;
        mid_post[s_mid] = w;
        total_mass = total_mass + w;
    }

    if (var::primitive(total_mass) <= 0.0)
        return Maybe_error<Result>(error_message(
            "micro_full_step_avg0_lifted: total emission mass is zero"));

    for (std::size_t s = 0; s < M; ++s) mid_post[s] = mid_post[s] / total_mass;

    // Second half-step: mid_post → post (per-element).
    Matrix<T> post(1, M);
    for (std::size_t s_end = 0; s_end < M; ++s_end) post[s_end] = type_zero;
    for (std::size_t s_mid = 0; s_mid < M; ++s_mid) {
        if (var::primitive(mid_post[s_mid]) <= 0.0) continue;
        for (std::size_t s_end = 0; s_end < M; ++s_end) {
            if (var::primitive(P_half_m(s_mid, s_end)) <= 0.0) continue;
            post[s_end] = post[s_end] + mid_post[s_mid] * P_half_m(s_mid, s_end);
        }
    }

    auto logL_step = max_log_L + log(total_mass);
    return Maybe_error<Result>(Result{std::move(post), std::move(logL_step)});
}


#if 0  // legacy cell-weighted FIM helpers; superseded by the moment-match
       // formula on the (mixture-corrected) y_var emitted from run_step.
template <class C_Prior_probs, class C_P, class C_gmean_ij, class C_gvar_ij, class C_sigma2>
auto compute_cell_weighted_gfi_avg2(C_Prior_probs const& prior_prob, C_P const& P_micro,
                                     C_gmean_ij const& gmean_ij_m, C_gvar_ij const& gvar_ij_m,
                                     C_sigma2 const& sigma2_obs) -> SymPosDefMatrix<double> {
    // Tier 2 complete-data FIM with prior weights. Cells are (s, s'); cell
    // weight w_c = π(s)·P(s,s'); per-cell score:
    //   score_c = ∂log π(s) + ∂log P(s,s') + ∂log 𝒩(y; μ_c, σ²_c)
    // Cross terms with ∂log 𝒩 average to zero under E_{y|s,s'}, and the
    // (∂log π)·(∂log P)^⊤ cross term collapses via Σ_{s'} P(s,s') = 1, so:
    //   G_step = Σ_s π(s)·(∂log π[s])(∂log π[s])^⊤
    //          + Σ_{s,s'} π(s)·P(s,s') · [(∂log P[s,s'])(...)^⊤ + I_𝒩(s,s')]
    std::size_t M = prior_prob.ncols();
    double sigma2_p = var::primitive(sigma2_obs);
    auto sample_d = derivative(gmean_ij_m(std::size_t{0}, std::size_t{0}))();
    std::size_t npar = sample_d.size();
    SymPosDefMatrix<double> H_step(npar, npar, 0.0);

    for (std::size_t s = 0; s < M; ++s) {
        double prior_s = var::primitive(prior_prob[s]);
        if (prior_s <= 0.0) continue;
        // π(s)·(∂log π[s])^⊗²: contribution from cell-prior score, summed
        // once per s.
        auto d_log_pi = derivative(prior_prob[s])() * (1.0 / prior_s);
        auto pi_outer = sqr_X<true>(d_log_pi);
        for (std::size_t i = 0; i < npar; ++i)
            for (std::size_t j = 0; j <= i; ++j)
                H_step.set(i, j, H_step(i, j) + prior_s * pi_outer(i, j));

        for (std::size_t s2 = 0; s2 < M; ++s2) {
            double Pss2 = var::primitive(P_micro(s, s2));
            if (Pss2 <= 0.0) continue;
            double weight = prior_s * Pss2;
            double sigma2_pair = sigma2_p + var::primitive(gvar_ij_m(s, s2));
            if (sigma2_pair <= 0.0 || !std::isfinite(sigma2_pair)) continue;

            // (∂log P[s,s'])^⊗² : transition score outer product, weighted π·P.
            auto d_log_P = derivative(P_micro(s, s2))() * (1.0 / Pss2);
            auto p_outer = sqr_X<true>(d_log_P);

            auto d_mu = derivative(gmean_ij_m(s, s2))();
            auto d_var = derivative(gvar_ij_m(s, s2))();
            auto first = sqr_X<true>(d_mu) * (1.0 / sigma2_pair);
            auto second = sqr_X<true>(d_var) * (1.0 / (2.0 * sigma2_pair * sigma2_pair));

            for (std::size_t i = 0; i < npar; ++i)
                for (std::size_t j = 0; j <= i; ++j)
                    H_step.set(i, j, H_step(i, j) +
                                          weight * (p_outer(i, j) + first(i, j) + second(i, j)));
        }
    }
    if (!lapack::matrix_has_only_finite(H_step))
        return SymPosDefMatrix<double>(npar, npar, 0.0);
    return H_step;
}

template <class C_Prior_probs, class C_P, class C_gmean_i, class C_gvar_i, class C_sigma2>
auto compute_cell_weighted_gfi_avg1(C_Prior_probs const& prior_prob, C_P const& P_micro,
                                     C_gmean_i const& gmean_i_m, C_gvar_i const& gvar_i_m,
                                     C_sigma2 const& sigma2_obs) -> SymPosDefMatrix<double> {
    // Tier 2 complete-data FIM with prior weights. Cells are (s, s'); weight
    // π(s)·P(s,s'); emission depends only on s, so I_𝒩(s) collapses out of
    // the s'-sum to π(s)·I_𝒩(s). Cross term (∂log π)·(∂log P)^⊤ vanishes
    // via Σ_{s'} P(s,s') = 1.
    //   G_step = Σ_s π(s) · [(∂log π[s])(...)^⊤ + I_𝒩(s)]
    //          + Σ_{s,s'} π(s)·P(s,s') · (∂log P[s,s'])(...)^⊤
    std::size_t M = prior_prob.ncols();
    double sigma2_p = var::primitive(sigma2_obs);
    auto sample_d = derivative(gmean_i_m(std::size_t{0}, std::size_t{0}))();
    std::size_t npar = sample_d.size();
    SymPosDefMatrix<double> H_step(npar, npar, 0.0);

    for (std::size_t s = 0; s < M; ++s) {
        double prior_s = var::primitive(prior_prob[s]);
        if (prior_s <= 0.0) continue;
        double sigma2_eff = sigma2_p + var::primitive(gvar_i_m(s, std::size_t{0}));
        if (sigma2_eff <= 0.0 || !std::isfinite(sigma2_eff)) continue;

        // π(s) · [(∂log π[s])^⊗² + I_𝒩(s)]
        auto d_log_pi = derivative(prior_prob[s])() * (1.0 / prior_s);
        auto pi_outer = sqr_X<true>(d_log_pi);
        auto d_mu = derivative(gmean_i_m(s, std::size_t{0}))();
        auto d_var = derivative(gvar_i_m(s, std::size_t{0}))();
        auto first = sqr_X<true>(d_mu) * (1.0 / sigma2_eff);
        auto second = sqr_X<true>(d_var) * (1.0 / (2.0 * sigma2_eff * sigma2_eff));

        for (std::size_t i = 0; i < npar; ++i)
            for (std::size_t j = 0; j <= i; ++j)
                H_step.set(i, j, H_step(i, j) +
                                      prior_s * (pi_outer(i, j) + first(i, j) + second(i, j)));

        // Σ_{s'} π(s)·P(s,s') · (∂log P[s,s'])^⊗²
        for (std::size_t s2 = 0; s2 < M; ++s2) {
            double Pss2 = var::primitive(P_micro(s, s2));
            if (Pss2 <= 0.0) continue;
            auto d_log_P = derivative(P_micro(s, s2))() * (1.0 / Pss2);
            auto p_outer = sqr_X<true>(d_log_P);
            double w = prior_s * Pss2;
            for (std::size_t i = 0; i < npar; ++i)
                for (std::size_t j = 0; j <= i; ++j)
                    H_step.set(i, j, H_step(i, j) + w * p_outer(i, j));
        }
    }
    if (!lapack::matrix_has_only_finite(H_step))
        return SymPosDefMatrix<double>(npar, npar, 0.0);
    return H_step;
}

template <class C_Mid_probs, class C_g_micro, class C_sigma2>
auto compute_cell_weighted_gfi_avg0(C_Mid_probs const& mid_prob, C_g_micro const& g_micro_m,
                                     C_sigma2 const& sigma2_obs) -> SymPosDefMatrix<double> {
    // Tier 2 complete-data FIM. Cells are midpoint microstates s_mid; weight
    // mid_prob[s_mid] = (prior @ P_half)[s_mid]; emission 𝒩(y; g_micro[s_mid],
    // σ²_obs) (σ²_obs constant per cell, so I_𝒩 has only the ∂μ term):
    //   G_step = Σ_{s_mid} mid_prob[s_mid] · [(∂log mid_prob[s_mid])(...)^⊤
    //                                       + (∂g_micro[s_mid])(...)^⊤/σ²_obs]
    std::size_t M = mid_prob.ncols();
    double sigma2_p = var::primitive(sigma2_obs);
    if (sigma2_p <= 0.0)
        return SymPosDefMatrix<double>{};
    auto sample_d = derivative(g_micro_m(std::size_t{0}, std::size_t{0}))();
    std::size_t npar = sample_d.size();
    SymPosDefMatrix<double> H_step(npar, npar, 0.0);

    for (std::size_t s_mid = 0; s_mid < M; ++s_mid) {
        double weight = var::primitive(mid_prob[s_mid]);
        if (weight <= 0.0) continue;

        auto d_log_mid = derivative(mid_prob[s_mid])() * (1.0 / weight);
        auto cell_weight_outer = sqr_X<true>(d_log_mid);

        auto d_g = derivative(g_micro_m(s_mid, std::size_t{0}))();
        auto emission_fim = sqr_X<true>(d_g) * (1.0 / sigma2_p);

        for (std::size_t i = 0; i < npar; ++i)
            for (std::size_t j = 0; j <= i; ++j)
                H_step.set(i, j,
                           H_step(i, j) +
                               weight * (cell_weight_outer(i, j) + emission_fim(i, j)));
    }
    if (!lapack::matrix_has_only_finite(H_step))
        return SymPosDefMatrix<double>(npar, npar, 0.0);
    return H_step;
}
#endif  // legacy cell-weighted FIM helpers


// -----------------------------------------------------------------------------
// log_Likelihood_micro — scalar likelihood via the lifted micro path.
//
// Lifts the macro Patch_Model to a microstate-granular Patch_Model_micro once,
// then drives the existing macro Calc_Qdt_step pipeline on it. Per-step update
// is plain Bayes on (s_start, s_end) microstate pairs via the *_lifted step
// functions defined above. Outer structure matches dlog_Likelihood_micro;
// derivative variant lives below.
// -----------------------------------------------------------------------------
template <class recursive, class averaging, class variance, class variance_correction,
          class qdt_method, class MacroState, class FuncTable, class C_Parameters, class Model>
auto log_Likelihood_micro(FuncTable& f, Model const& model, C_Parameters const& par,
                          Recording const& y, Experiment const& e) -> Maybe_error<MacroState> {
    auto Maybe_m = model(par);
    if (!is_valid(Maybe_m)) return get_error(Maybe_m);
    auto m_macro = std::move(get_value(Maybe_m));
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto f_local = f.create("_lik_micro");

    Macro_DMR macro_dmr{};

    auto ini = macro_dmr.init(m_macro);
    if (!ini) return ini.error();

    std::size_t k_states = get<N_St>(m_macro)();
    auto Nchs = get<N_Ch_mean>(m_macro)();
    std::size_t N_ch = static_cast<std::size_t>(std::round(var::primitive(Nchs[0])));
    auto full_states = create_Micro_state_Num_ch(N_ch, k_states);

    // Initial microstate distribution from the multinomial lift of macro init.
    auto Maybe_mps = lift_Macro_to_Micro(full_states, ini.value());
    if (!Maybe_mps) return Maybe_mps.error();

    // Lift the macro Patch_Model once — feeds calc_Qdt at the microstate level.
    auto m_micro = lift_Patch_Model_to_Micro(m_macro, full_states);

    auto initial_probs = get<Micro_P_mean>(Maybe_mps.value()())();
    using T_probs = std::decay_t<decltype(initial_probs)>;
    auto zero_logL = initial_probs[std::size_t{0}] - initial_probs[std::size_t{0}];
    using T_logL = std::decay_t<decltype(zero_logL)>;

    struct Micro_fold_state {
        T_probs probs;
        T_logL logL_acc;
    };

    Micro_fold_state init_state{std::move(initial_probs), std::move(zero_logL)};

    auto Maybe_run = fold(
        0ul, y().size(), std::move(init_state),
        [&](Micro_fold_state&& prior, std::size_t i_step) -> Maybe_error<Micro_fold_state> {
            Agonist_evolution const& t_step =
                get<Agonist_evolution>(get<Recording_conditions>(e)()[i_step]);
            if (t_step().size() != 1)
                return error_message(
                    "micro lifted path requires single-substep Agonist_evolution");
            Agonist_step const& sub_step = t_step()[0];
            // y_obs is the raw recorded current; the macro pipeline's gmean_ij /
            // gmean_i / g are *baseline-free* per-microstate currents. Subtract
            // baseline once so the step function's residuals are y_obs − μ_pair
            // exactly as the macro path computes them via y_mean = N·... + baseline.
            auto y_obs_centered =
                y()[i_step].value() - get<Current_Baseline>(m_micro).value();

            // Helper: select Qdtg / Qdt computation method based on the
            // qdt_method flag (int): 0=eig, 1=taylor, 2=schur.
            // For Nch ≥ ~10 with k=2 the lifted Q_micro has clustered
            // eigenvalues, so 0 (eig) is ill-conditioned and 2 (schur) is
            // the recommended path. 1 (taylor) is kept for benchmarking.
            if constexpr (averaging::value == 0) {
                auto Maybe_Qdtg = [&] {
                    if constexpr (qdt_method::value == 1) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdtg_taylor(m_micro, t_Qx,
                                                           get<number_of_samples>(sub_step),
                                                           dt_sub);
                    } else if constexpr (qdt_method::value == 2) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdtg_schur(m_micro, t_Qx,
                                                          get<number_of_samples>(sub_step),
                                                          dt_sub);
                    } else {
                        return macro_dmr.calc_Qdtg(f_local, m_micro, sub_step, fs);
                    }
                }();
                if (!Maybe_Qdtg) return Maybe_Qdtg.error();
                auto const& t_Qdtg = Maybe_Qdtg.value();
                auto const& P_half_m = get<P_half>(t_Qdtg)();
                auto const& g_micro_v = get<g>(m_micro)();
                double n_samples = static_cast<double>(get<number_of_samples>(t_Qdtg).value());
                auto sigma2_obs = get<Current_Noise>(m_micro).value() * fs / n_samples;

                auto Maybe_step = micro_full_step_avg0_lifted(prior.probs, P_half_m, g_micro_v,
                                                               sigma2_obs, y_obs_centered);
                if (!Maybe_step) return Maybe_step.error();
                auto step = std::move(Maybe_step.value());
                return Micro_fold_state{std::move(step.new_probs),
                                        prior.logL_acc + step.logL_contribution};
            } else {
                auto Maybe_Qdt = [&] {
                    if constexpr (qdt_method::value == 1) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdt_taylor(m_micro, t_Qx,
                                                         get<number_of_samples>(sub_step),
                                                         dt_sub);
                    } else if constexpr (qdt_method::value == 2) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        // Plain double uses Schur+Parlett + Pade-VanLoan;
                        // Derivative uses Pade scaling-and-squaring through
                        // the templated Frechet integrals (native Phase 1
                        // path — no Taylor delegation).
                        return macro_dmr.calc_Qdt_schur(m_micro, t_Qx,
                                                         get<number_of_samples>(sub_step),
                                                         dt_sub);
                    } else {
                        return macro_dmr.calc_Qdt(f_local, m_micro, sub_step, fs);
                    }
                }();
                if (!Maybe_Qdt) return Maybe_Qdt.error();
                auto const& t_Qdt = Maybe_Qdt.value();
                auto const& P_micro_v = get<P>(t_Qdt)();
                double n_samples = static_cast<double>(get<number_of_samples>(t_Qdt).value());
                auto sigma2_obs = get<Current_Noise>(m_micro).value() * fs / n_samples;

                if constexpr (averaging::value == 1) {
                    auto const& gmean_i_v = get<gmean_i>(t_Qdt)();
                    auto const& gvar_i_v = get<gvar_i>(t_Qdt)();
                    auto Maybe_step =
                        micro_full_step_avg1_lifted(prior.probs, P_micro_v, gmean_i_v, gvar_i_v,
                                                     sigma2_obs, y_obs_centered);
                    if (!Maybe_step) return Maybe_step.error();
                    auto step = std::move(Maybe_step.value());
                    return Micro_fold_state{std::move(step.new_probs),
                                            prior.logL_acc + step.logL_contribution};
                } else {  // averaging::value == 2
                    auto const& gmean_ij_v = get<gmean_ij>(t_Qdt)();
                    auto const& gvar_ij_v = get<gvar_ij>(t_Qdt)();
                    auto Maybe_step =
                        micro_full_step_avg2_lifted(prior.probs, P_micro_v, gmean_ij_v, gvar_ij_v,
                                                     sigma2_obs, y_obs_centered);
                    if (!Maybe_step) return Maybe_step.error();
                    auto step = std::move(Maybe_step.value());
                    return Micro_fold_state{std::move(step.new_probs),
                                            prior.logL_acc + step.logL_contribution};
                }
            }
        });

    f += f_local;
    if (!Maybe_run) return Maybe_run.error();
    auto final_state = std::move(Maybe_run.value());

    Micro_Patch_State final_mps;
    get<Micro_P_mean>(final_mps()) = Micro_P_mean(std::move(final_state.probs));
    auto final_patch_state = project_Micro_to_Macro(full_states, final_mps);

    MacroState out;
    get<logL>(out) = logL(final_state.logL_acc);
    get<Patch_State>(out) = std::move(final_patch_state);
    return Maybe_error<MacroState>(std::move(out));
}


// -----------------------------------------------------------------------------
// dlog_Likelihood_micro — derivative-aware outer function paralleling
// log_Likelihood_micro. Lifts m_macro → m_micro once, drives the macro
// Calc_Qdt_step / calc_Qdtg pipeline on the lifted model, runs plain-Bayes
// per-step updates via micro_full_step_avg{0,1,2}_lifted. Per-step
// bookkeeping (Evolution_of) uses the per-microstate current vector V
// directly (V = g_micro for avg=0, V = gmean_i_micro for avg≥1).
// -----------------------------------------------------------------------------
template <class recursive, class averaging, class variance, class variance_correction,
          class qdt_method, class MacroState, class FuncTable, class Model>
auto dlog_Likelihood_micro(FuncTable& f, Model const& model,
                           var::Parameters_transformed const& par,
                           Recording const& y, Experiment const& e) -> Maybe_error<MacroState> {
    using std::sqrt;
    using std::abs;
    auto dp = var::selfDerivative(par);
    auto dpp = dp.to_value();

    auto Maybe_m = model(dpp);
    if (!is_valid(Maybe_m)) return get_error(Maybe_m);
    auto m_macro = std::move(get_value(Maybe_m));
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto f_local = f.create("_dlik_micro");

    Macro_DMR macro_dmr{};

    auto ini = macro_dmr.init(m_macro);
    if (!ini) return ini.error();

    std::size_t k_states = get<N_St>(m_macro)();
    auto Nchs = get<N_Ch_mean>(m_macro)();
    std::size_t N_ch = static_cast<std::size_t>(std::round(var::primitive(Nchs[std::size_t{0}])));
    auto full_states = create_Micro_state_Num_ch(N_ch, k_states);

    // Initial microstate distribution from the multinomial lift of macro init.
    auto Maybe_mps = lift_Macro_to_Micro(full_states, ini.value());
    if (!Maybe_mps) return Maybe_mps.error();

    // Lift the macro Patch_Model once.
    auto m_micro = lift_Patch_Model_to_Micro(m_macro, full_states);

    auto init_macro = MacroState(std::move(ini.value()));

    auto initial_probs = get<Micro_P_mean>(Maybe_mps.value()())();
    using T_probs = std::decay_t<decltype(initial_probs)>;

    struct Fold_state {
        MacroState macro_state;
        T_probs probs;
    };

    Fold_state init_state{std::move(init_macro), std::move(initial_probs)};

    auto Maybe_run = fold(
        0ul, y().size(), std::move(init_state),
        [&](Fold_state&& prior, std::size_t i_step) -> Maybe_error<Fold_state> {
            Agonist_evolution const& t_step =
                get<Agonist_evolution>(get<Recording_conditions>(e)()[i_step]);
            if (t_step().size() != 1)
                return error_message(
                    "micro lifted path requires single-substep Agonist_evolution");
            Agonist_step const& sub_step = t_step()[0];

            double y_obs = y()[i_step].value();
            bool y_is_nan = std::isnan(y_obs);
            auto y_baseline = get<Current_Baseline>(m_macro);
            // gmean_ij / gmean_i / g_micro from the macro pipeline are baseline-free
            // per-microstate currents. Subtract baseline once so the step function's
            // residuals match the macro y_mean = N·... + baseline construction.
            auto y_obs_centered = y_obs - y_baseline.value();

            // Per-step Evolution_of bookkeeping. V is the per-microstate
            // current vector (M-vector); `mean_v` and `var_v` are the mixture
            // moments of V under prior.probs. `within_v` is the prior-
            // weighted mean of the per-cell within-emission variance:
            //   avg=0 (micro_R):  within_v = 0       (point Gaussians per cell)
            //   avg=1 (micro_MR): within_v = Σ_s π(s)·gvar_i[s]
            //   avg=2 (micro_IR): within_v = Σ_{s,s'} π(s)·P(s,s')·gvar_ij[s,s']
            // y_var = e_noise + var_v + within_v matches the marginal mixture
            // variance Var[y] = E_c[Var[y|c]] + Var_c[E[y|c]], which makes the
            // diagnostic moment-match formula well-defined as the FIM of the
            // single Gaussian whose mean and variance match the mixture.
            // The per-step Gaussian-formula FIM consumed by GFD is computed
            // downstream by the diagnostic preset directly from y_mean/y_var
            // derivatives — no per-step FIM slot is needed.
            auto run_step = [&](auto const& V, auto const& sigma2_obs,
                                auto const& within_v,
                                auto&& step) -> Maybe_error<Fold_state> {
                auto V_sq = elemMult(V, V);
                auto mean_v = getvalue(prior.probs * V);
                auto mean_v_sq = getvalue(prior.probs * V_sq);
                auto var_v = mean_v_sq - mean_v * mean_v;

                auto r_y_mean = build<y_mean>(mean_v + y_baseline());
                auto dy = y_obs - r_y_mean();
                auto e_noise = sigma2_obs +
                               get<Pink_Noise>(m_macro).value() +
                               get<Proportional_Noise>(m_macro).value() * abs(dy);
                auto r_y_var = build<y_var>(e_noise + var_v + within_v);
                auto r_r_std = build<r_std>(dy / sqrt(r_y_var()));
                auto chi_t = dy / r_y_var();
                auto r_chi2 = build<Chi2>(dy * chi_t);
                auto alfa = trust_coefficient(1.0);
                auto d_logL = build<logL>(step.logL_contribution);
                auto d_elogL = macro_dmr.calculate_elogL(y_is_nan, r_y_var, m_macro);

                Transfer_Op_to<T_probs, Algo_State> r_algo;
                get<y_mean>(r_algo()) = r_y_mean;
                get<y_var>(r_algo()) = r_y_var;
                get<r_std>(r_algo()) = r_r_std;
                get<Chi2>(r_algo()) = r_chi2;
                get<trust_coefficient>(r_algo()) = alfa;

                auto Maybe_upd = macro_dmr.update(std::move(prior.macro_state),
                                                  std::move(r_algo), d_logL, d_elogL,
                                                  r_y_mean, r_y_var);
                if (!Maybe_upd) return Maybe_upd.error();
                // _lifted step functions return per-element Matrix<Derivative<double>>
                // in derivative mode; T_probs is monolithic Derivative<Matrix<double>>.
                // outside_in re-packs; identity in plain mode.
                if constexpr (var::is_derivative_v<T_probs>) {
                    return Fold_state{std::move(Maybe_upd.value()),
                                      var::outside_in(step.new_probs, prior.probs.dx())};
                } else {
                    return Fold_state{std::move(Maybe_upd.value()),
                                      std::move(step.new_probs)};
                }
            };

            // Qdtg / Qdt method: 0=eig (eigen dispatcher), 1=taylor (Taylor
            // scaling+squaring), 2=schur (Schur+Parlett+Pade-VanLoan). 2 is
            // the recommended path at high Nch where Q_micro's eigenvalues
            // cluster and the eigenvector-derivative inverse becomes
            // ill-conditioned.
            if constexpr (averaging::value == 0) {
                auto Maybe_Qdtg = [&] {
                    if constexpr (qdt_method::value == 1) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdtg_taylor(m_micro, t_Qx,
                                                           get<number_of_samples>(sub_step),
                                                           dt_sub);
                    } else if constexpr (qdt_method::value == 2) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdtg_schur(m_micro, t_Qx,
                                                          get<number_of_samples>(sub_step),
                                                          dt_sub);
                    } else {
                        return macro_dmr.calc_Qdtg(f_local, m_micro, sub_step, fs);
                    }
                }();
                if (!Maybe_Qdtg) return Maybe_Qdtg.error();
                auto const& t_Qdtg = Maybe_Qdtg.value();
                auto const& P_half_m = get<P_half>(t_Qdtg)();
                auto const& g_micro_v = get<g>(m_micro)();
                double n_samples =
                    static_cast<double>(get<number_of_samples>(t_Qdtg).value());
                auto sigma2_obs = get<Current_Noise>(m_micro).value() * fs / n_samples;

                auto Maybe_step = micro_full_step_avg0_lifted(prior.probs, P_half_m, g_micro_v,
                                                               sigma2_obs, y_obs_centered);
                if (!Maybe_step) return Maybe_step.error();
                // micro_R: per-cell emission is a point-Gaussian at g_micro[s_mid]
                // with variance σ²_obs only — no within-cell variance term.
                auto within_v = sigma2_obs - sigma2_obs;
                return run_step(g_micro_v, sigma2_obs, within_v,
                                std::move(Maybe_step.value()));
            } else {
                auto Maybe_Qdt = [&] {
                    if constexpr (qdt_method::value == 1) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdt_taylor(m_micro, t_Qx,
                                                         get<number_of_samples>(sub_step),
                                                         dt_sub);
                    } else if constexpr (qdt_method::value == 2) {
                        auto t_Qx = build<Qx>(macro_dmr.calc_Qx(
                            m_micro, get<Agonist_concentration>(sub_step)));
                        double dt_sub = get<number_of_samples>(sub_step)() / fs;
                        return macro_dmr.calc_Qdt_schur(m_micro, t_Qx,
                                                         get<number_of_samples>(sub_step),
                                                         dt_sub);
                    } else {
                        return macro_dmr.calc_Qdt(f_local, m_micro, sub_step, fs);
                    }
                }();
                if (!Maybe_Qdt) return Maybe_Qdt.error();
                auto const& t_Qdt = Maybe_Qdt.value();
                auto const& P_micro_v = get<P>(t_Qdt)();
                auto const& t_gmean_i = get<gmean_i>(t_Qdt)();
                double n_samples =
                    static_cast<double>(get<number_of_samples>(t_Qdt).value());
                auto sigma2_obs = get<Current_Noise>(m_micro).value() * fs / n_samples;

                if constexpr (averaging::value == 1) {
                    auto const& gvar_i_v = get<gvar_i>(t_Qdt)();
                    auto Maybe_step = micro_full_step_avg1_lifted(
                        prior.probs, P_micro_v, t_gmean_i, gvar_i_v, sigma2_obs, y_obs_centered);
                    if (!Maybe_step) return Maybe_step.error();
                    // micro_MR: within_v = Σ_s π(s) · gvar_i[s].
                    auto within_v = getvalue(prior.probs * gvar_i_v);
                    return run_step(t_gmean_i, sigma2_obs, within_v,
                                    std::move(Maybe_step.value()));
                } else {  // averaging::value == 2
                    auto const& gmean_ij_v = get<gmean_ij>(t_Qdt)();
                    auto const& gvar_ij_v = get<gvar_ij>(t_Qdt)();
                    auto Maybe_step = micro_full_step_avg2_lifted(
                        prior.probs, P_micro_v, gmean_ij_v, gvar_ij_v, sigma2_obs, y_obs_centered);
                    if (!Maybe_step) return Maybe_step.error();
                    // micro_IR: within_v = Σ_{s, s'} π(s)·P(s, s')·gvar_ij[s, s']
                    //          = sum(prior · elemMult(P, gvar_ij))
                    // Single matrix expression — consolidates the auto-diff
                    // chain into elemMult / matmul / sum (three ops) instead
                    // of M² scalar accumulations, which at Nch=50 is ~2600
                    // chained Derivative additions and amplifies derivative
                    // round-off.
                    auto within_v =
                        var::sum(prior.probs * elemMult(P_micro_v, gvar_ij_v));
                    return run_step(t_gmean_i, sigma2_obs, within_v,
                                    std::move(Maybe_step.value()));
                }
            }
        });

    f += f_local;
    if (!Maybe_run) {
        return Maybe_run.error();
    }
    auto final_state = std::move(Maybe_run.value());

    // Single project at the end of the fold: convert the final Micro
    // distribution to a (P_mean, P_Cov) Patch_State for the returned MacroState.
    Transfer_Op_to<T_probs, Micro_Patch_State> final_mps;
    get<Micro_P_mean>(final_mps())() = std::move(final_state.probs);
    auto final_patch_state = project_Micro_to_Macro(full_states, final_mps);
    get<var::Derivative<Patch_State, var::Parameters_transformed>>(final_state.macro_state) =
        std::move(final_patch_state);

    return Maybe_error<MacroState>(std::move(final_state.macro_state));
}


}  // namespace macrodr

#endif  // MICRO_FULL_H
