#pragma once

#include <derivative_fwd.h>
#include <matrix.h>
#include <micro_types.h>
#include <parameters_derivative.h>
#include <qmodel.h>
#include <qmodel_types.h>
#include <cstddef>
#include <utility>
#include "abstract_algebra.h"

    namespace macrodr {

// Bring var::sum into scope so unqualified `sum(elemMult(..., ...))` calls in
// the safely_calculate_micro_Algo_State_* bodies resolve. ADL via Matrix<double>
// (declared in global namespace) cannot find var::sum on its own.
using var::sum;

inline std::size_t num_full_states_of_new(std::size_t N_channels, std::size_t k_states) {
    if (k_states == 0) {
        return 0;
    }
    std::size_t r = 1;
    for (std::size_t i = 1; i < k_states; ++i) {
        r = (r * (N_channels + i)) / i;
    }
    return r;
}

// Inverse of micro_full_detail::index_of_microstate: given a lex-order index
// into the enumeration produced by create_Micro_state_Num_ch(N, k), return the
// occupation vector of length k whose entries sum to N.
inline std::vector<std::size_t> index_to_microstate(std::size_t idx, std::size_t N,
                                                    std::size_t k) {
    if (k == 0) {
        return {};
    }
    std::vector<std::size_t> out(k, 0);
    std::size_t remaining = N;
    for (std::size_t i = 0; i + 1 < k; ++i) {
        std::size_t ni = 0;
        while (idx >= num_full_states_of_new(remaining - ni, k - 1 - i)) {
            idx -= num_full_states_of_new(remaining - ni, k - 1 - i);
            ++ni;
        }
        out[i] = ni;
        remaining -= ni;
    }
    out[k - 1] = remaining;
    return out;
}

inline std::size_t microstate_to_index(std::vector<std::size_t> const& n,
                                        std::size_t N,
                                       std::size_t k) {
    std::size_t idx = 0;
    std::size_t remaining = N;
    for (std::size_t i = 0; i + 1 < k; ++i) {
        std::size_t ni = n[i];
        for (std::size_t j = 0; j < ni; ++j) {
            idx += num_full_states_of_new(remaining - j, k - 1 - i);
        }
        remaining -= ni;
    }
    return idx;

}

// Multinomial coefficient  N! / (n_0! · n_1! · … · n_{k-1}!)  for an unordered
// channel-count vector n with total = N.  Computed multiplicatively to avoid
// overflowing factorials for moderate N.
//
// Used in P_to_micro_P (and downstream) to undo the multinomial over-counting
// that arises when projecting an ordered N-channel transition matrix to the
// unordered count space: the codebase's convolution
//     P_out(i,j) += P0(i1,j1) · P1(i2,j2)        (line 195)
// sums over BOTH start orderings (the M(i) decompositions of i_states) and
// end orderings (the M(j) decompositions of j_states). Summing over end
// orderings is correct (gives the unordered transition probability — the
// multinomial weighting of the end is required); summing over start orderings
// over-counts each row by M(i). The final out.first thus has row i sum = M(i).
// Dividing each row by M(i) yields the proper row-stochastic ensemble
// transition matrix that  unordered_probs · P  consumers expect.
//
// Apply ONLY at the final level (after power_semigroup completes). Dividing
// at intermediate levels would feed row-stochastic matrices into the next
// squaring, where the convolution would then yield row sums equal to the
// number of unordered decompositions of the next-level state — strictly
// less than M(parent) when the parent has nested multinomial structure.
// One division at the end recovers the correct M(i_at_Nchannels) compensation.
inline double multinomial_count(std::vector<std::size_t> const& n) {
    std::size_t total = 0;
    for (auto x : n) total += x;
    double result = 1.0;
    std::size_t remaining = total;
    for (auto x : n) {
        for (std::size_t k = 0; k < x; ++k)
            result *= static_cast<double>(remaining - k) / static_cast<double>(k + 1);
        remaining -= x;
    }
    return result;
}


template <class>
inline constexpr bool always_false_v = false;

class Micro_DMR {

   public:
    
     
     template <class C_P_mean>
     requires (U<C_P_mean, P_mean>||U<C_P_mean, P_initial>)
     static auto P_mean_to_micro_P_state(const C_P_mean& t_P_mean, std::size_t Nchannels) {
            auto k= t_P_mean().ncols();
            using Matr=Matrix<Transfer_Op_to<C_P_mean, double>>;

             auto out= power_semigroup([k]( std::pair<Matr, std::size_t>  PN0, std::pair<Matr, std::size_t> PN1){
                auto N0= PN0.second;
                auto N1= PN1.second;
                auto& P0= PN0.first;
                auto& P1= PN1.first;

                auto num_states = num_full_states_of_new(N0+N1, k);
                Matr P_out(1, num_states);
                for (std::size_t i = 0; i < P0.ncols(); ++i) {
                 auto i_states= index_to_microstate(i,  N0,  k)  ;
                 for (std::size_t j = 0; j < P1.ncols(); ++j) {
                    auto j_states= index_to_microstate(j,  N1,  k);
                    auto ij_states = i_states + j_states;
                    auto ij_index= microstate_to_index(ij_states, N0+N1, k);
                    P_out[ij_index] += P0[i] * P1[j];
                }
            }
                return std::make_pair(P_out, N0+N1);
            },
              std::make_pair(var::inside_out(t_P_mean()),std::size_t{1}), Nchannels);
              if constexpr (U<C_P_mean, P_mean>) {
                return build<micro_P>(var::outside_in(std::move(out.first), t_P_mean));
              } else {
                static_assert(always_false_v<C_P_mean>, "unhandled P_mean type");
              }
     }

template <class C_g>
     requires (U<C_g, g>|| U<C_g,gmean_i> ||  U<C_g,gvar_i> )
     static auto g_to_micro_g(const C_g& t_g, std::size_t Nchannels) {
            auto k = t_g().nrows();
            using Matr = Matrix<Transfer_Op_to<C_g, double>>;

            auto out = power_semigroup(
                [k](std::pair<Matr, std::size_t> gN0,
                    std::pair<Matr, std::size_t> gN1) {
                    auto N0 = gN0.second;
                    auto N1 = gN1.second;
                    auto& g0 = gN0.first;
                    auto& g1 = gN1.first;

                    auto num_states = num_full_states_of_new(N0 + N1, k);
                    Matr g_out(1, num_states);
                    Matrix<double> state_count(num_states, 1UL, 0.0);
                    // g/gmean_i/gvar_i seeds are column vectors (k,1); after
                    // the first squaring g0/g1 become row vectors
                    // (1, num_states_so_far). Iterate by .size() so the
                    // microstate enumeration is independent of layout.
                    for (std::size_t i = 0; i < g0.size(); ++i) {
                        auto i_states = index_to_microstate(i, N0, k);
                        for (std::size_t j = 0; j < g1.size(); ++j) {
                            auto j_states = index_to_microstate(j, N1, k);
                            auto ij_states = i_states + j_states;
                            auto ij_index =
                                microstate_to_index(ij_states, N0 + N1, k);
                            // Sum over all decomposition paths (i, j) that
                            // map to ij_index — each contributes g(i)+g(j).
                            g_out[ij_index] += g0[i] + g1[j];
                            state_count(ij_index, 0) += 1;
                        }
                    }
                    // Per-merge: divide each cell by the number of (i, j) decompositions
                    // that landed there. For 1D additive (g, gmean_i, gvar_i) the cell
                    // already sums g_total over orderings; dividing by ordering count
                    // recovers the per-state total conductance. The mixed-cell elemDiv
                    // overload (matrix.h) handles g_out: Matrix<Derivative<double>>
                    // vs state_count: Matrix<double> in the derivative path.
                    g_out = elemDiv(g_out, state_count);
                    return std::make_pair(g_out, N0 + N1);
                },
                std::make_pair(var::inside_out(t_g()), std::size_t{1}), Nchannels);

            if constexpr (U<C_g, g>) {
                return build<micro_g>(var::outside_in(std::move(out.first), t_g));
            } else if constexpr (U<C_g, gmean_i>) {
                return build<micro_gmean_i>(var::outside_in(std::move(out.first), t_g));
            } else if constexpr (U<C_g, gvar_i>) {
                return build<micro_gvar_i>  (var::outside_in(std::move(out.first), t_g));
            } else {
                static_assert(always_false_v<C_g>, "unhandled g type");
            }
     }

     template <class C_P>
     requires (U<C_P, P>||U<C_P, P_half>)
     static auto P_to_micro_P(const C_P& t_P, std::size_t Nchannels) {
            auto k = t_P().ncols();
            using Matr = Matrix<Transfer_Op_to<C_P, double>>;

            auto out = power_semigroup(
                [k](std::pair<Matr, std::size_t> PN0,
                    std::pair<Matr, std::size_t> PN1) {
                    auto N0 = PN0.second;
                    auto N1 = PN1.second;
                    auto& P0 = PN0.first;
                    auto& P1 = PN1.first;

                    auto num_states = num_full_states_of_new(N0 + N1, k);
                    Matr P_out(num_states, num_states);
                    Matrix<double> state_count(num_states, 1UL, 0.0);
                    // i1, j1 index N0-microstates (rows/cols of P0).
                    // i2, j2 index N1-microstates (rows/cols of P1).
                    for (std::size_t i1 = 0; i1 < P0.nrows(); ++i1) {
                        auto i1_states = index_to_microstate(i1, N0, k);
                        for (std::size_t i2 = 0; i2 < P1.nrows(); ++i2) {
                                auto i2_states =
                                    index_to_microstate(i2, N1, k);
                                    auto i_states = i1_states + i2_states;
                                    auto i_index = microstate_to_index(
                                        i_states, N0 + N1, k);
                                    state_count(i_index, 0) += 1;  
                            for (std::size_t j1 = 0; j1 < P0.ncols(); ++j1) {
                            auto j1_states = index_to_microstate(j1, N0, k);
                                for (std::size_t j2 = 0; j2 < P1.ncols();
                                     ++j2) {
                                    auto j2_states =
                                        index_to_microstate(j2, N1, k);
                                    auto j_states = j1_states + j2_states;
                                    auto j_index = microstate_to_index(
                                        j_states, N0 + N1, k);
                                    // Joint transition probability convolves
                                    // over decompositions: += is correct.
                                    P_out(i_index, j_index) +=
                                        P0(i1, j1) * P1(i2, j2);
                                }
                            }
                        }
                    }
                    // state_count(i, 0) counts the number of (i1, i2) decompositions
                    // of i_states. Divide row i of P_out by state_count(i, 0) so each
                    // row sums to 1 (row-stochastic ensemble transition).
                    for (std::size_t i = 0; i < P_out.nrows(); ++i) {
                        double n = state_count(i, 0);
                        if (n > 1.0) {
                            double inv_n = 1.0 / n;
                            for (std::size_t j = 0; j < P_out.ncols(); ++j) {
                                P_out(i, j) = P_out(i, j) * inv_n;
                            }
                        }
                    }
                    return std::make_pair(P_out, N0 + N1);
                },
                std::make_pair(var::inside_out(t_P()), std::size_t{1}), Nchannels);

            if constexpr (U<C_P, P>) {
                return build<micro_P>(var::outside_in(std::move(out.first), t_P));
            } else if constexpr (U<C_P, P_half>) {
                return build<micro_P_half>(var::outside_in(std::move(out.first), t_P));
            } else {
                static_assert(always_false_v<C_P>, "unhandled P type");
            }
     }
    
    static auto calc_micro_P_to_P_mean(std::size_t k, std::size_t Nchannels) {
            auto num_states = num_full_states_of_new(Nchannels, k);
            Matrix<double> P_to_Pmean(num_states, k);
            for (std::size_t s = 0; s < num_states; ++s) {
                auto n = index_to_microstate(s, Nchannels, k);
                for (std::size_t i = 0; i < k; ++i) {
                    P_to_Pmean(s, i) = static_cast<double>(n[i]) / Nchannels;
                }
            }
            return build<micro_P_state_to_P_mean>   (std::move(P_to_Pmean));
        }

     template <class C_gij>
     requires (U<C_gij, gmean_ij>|| U<C_gij, gvar_ij> )
     static auto gij_to_micro_gij(const C_gij t_gij, std::size_t Nchannels) {
            auto k = t_gij().ncols();
            using Matr = Matrix<Transfer_Op_to<C_gij, double>>;

            auto out = power_semigroup(
                [k](std::pair<Matr, std::size_t> g_ij_N0,
                    std::pair<Matr, std::size_t> g_ij_N1) {
                    auto N0 = g_ij_N0.second;
                    auto N1 = g_ij_N1.second;
                    auto& gij_0 = g_ij_N0.first;
                    auto& gij_1 = g_ij_N1.first;

                    auto num_states = num_full_states_of_new(N0 + N1, k);
                    Matr gij_out(num_states, num_states);
                    Matrix<double> state_count(num_states, num_states, 0.0);
                    for (std::size_t i1 = 0; i1 < gij_0.nrows(); ++i1) {
                        auto i1_states = index_to_microstate(i1, N0, k);
                        for (std::size_t j1 = 0; j1 < gij_0.ncols(); ++j1) {
                            auto j1_states = index_to_microstate(j1, N0, k);
                            for (std::size_t i2 = 0; i2 < gij_1.nrows();
                                 ++i2) {
                                auto i2_states =
                                    index_to_microstate(i2, N1, k);
                                for (std::size_t j2 = 0; j2 < gij_1.ncols();
                                     ++j2) {
                                    auto j2_states =
                                        index_to_microstate(j2, N1, k);
                                    auto i_states = i1_states + i2_states;
                                    auto j_states = j1_states + j2_states;
                                    auto i_index = microstate_to_index(
                                        i_states, N0 + N1, k);
                                    auto j_index = microstate_to_index(
                                        j_states, N0 + N1, k);
                                    // Sum over all decomposition paths
                                    // (i1, j1, i2, j2) that reach
                                    // (i_index, j_index) — each contributes
                                    // gij_0(i1,j1) + gij_1(i2,j2).
                                    gij_out(i_index, j_index) +=
                                        gij_0(i1, j1) + gij_1(i2, j2);
                                    state_count(i_index, j_index) += 1;    
                                }
                            }
                        }
                    }
                    // Per-merge: state_count(i, j) counts (i1, j1, i2, j2) tuples that
                    // landed in (i_index, j_index). For 2D additive, divide each cell by
                    // its tuple count to recover the per-(start, end) joint conductance.
                    gij_out = elemDiv(gij_out, state_count);
                    return std::make_pair(gij_out, N0 + N1);
                },
                std::make_pair(var::inside_out(t_gij()), std::size_t{1}), Nchannels);


            if constexpr (U<C_gij, gmean_ij>) {
                return build<micro_gmean_ij>(var::outside_in(std::move(out.first), t_gij));
            } else if constexpr (U<C_gij, gvar_ij>) {
                return build<micro_gvar_ij>(var::outside_in(std::move(out.first), t_gij));
            } else {
                static_assert(always_false_v<C_gij>, "unhandled g_ij type");
            }
     }
    



     template <class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtg(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs, std::size_t Num_channels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtg>> {
        if constexpr (std::is_same_v<Nothing, decltype(f[Calc_micro_Qdtg_step{}])>)
            return calc_micro_Qdtg_agonist_step(f, m, t_step, fs, Num_channels);
        else
            return f.f(Calc_micro_Qdtg_step{}, m, t_step, fs, Num_channels );
    }


   template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_micro_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtg_agonist_step(FunctionTable& f, const C_micro_Patch_Model& m,
                                const Agonist_step& t_step, double fs, std::size_t Num_chanels)
        -> Maybe_error<Transfer_Op_to<C_micro_Patch_Model, micro_Qdtg>> {
        
        auto Maybe_Qdt= Macro_DMR{}.calc_Qdtg(f, m, t_step, fs);
        if (!Maybe_Qdt) return Maybe_Qdt.error();
        auto const& t_Qdt = Maybe_Qdt.value();
        auto  P_half_m = P_to_micro_P(get<P_half>(t_Qdt), Num_chanels);
        auto  t_P_to_Pmean = calc_micro_P_to_P_mean(get<P_half>(t_Qdt)().ncols(), Num_chanels);
        auto  g_micro_v = g_to_micro_g(get<g>(m), Num_chanels);
        return build<micro_Qdtg>(get<number_of_samples>(t_Qdt), get<min_P>(t_Qdt), std::move(P_half_m),
                                 std::move(t_P_to_Pmean), std::move(g_micro_v));
    }



    
    template <class Policy = StabilizerPolicyEnabled, class C_Q0, class C_Qa>
        requires U<C_Q0, Q0>
    Maybe_error<Transfer_Op_to<C_Q0, P_initial>> calc_micro_Pinitial(const C_Q0& t_Q0, const C_Qa& t_Qa,
                                                               Agonist_concentration x,
                                                               N_St nstates, std::size_t Nchannels) {

        auto Maybe_Pinitial = Macro_DMR{}.calc_Pinitial<Policy>(t_Q0, t_Qa, x, nstates);
        if (!Maybe_Pinitial) return Maybe_Pinitial.error();
        auto const& t_Pinitial = Maybe_Pinitial.value();
        return P_mean_to_micro_P_state(t_Pinitial, Nchannels);
    }



    

    
    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtm_agonist_step(FunctionTable& f, const C_Patch_Model& m,
                                const Agonist_step& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtm>> {
        auto maybe_Qdt = Macro_DMR{}.calc_Qdtm_agonist_step(f, m, t_step, fs);
        if (!maybe_Qdt)
            return maybe_Qdt.error();
        auto& Qdt = maybe_Qdt.value();
        auto t_micro_P = P_to_micro_P(get<P>(Qdt), Nchannels);
        auto t_P_to_Pmean = calc_micro_P_to_P_mean(get<P>(Qdt)().ncols(), Nchannels);
        auto t_gmean_i   = g_to_micro_g(get<gmean_i>(Qdt), Nchannels);
        auto t_gvar_i    = g_to_micro_g(get<gvar_i>(Qdt), Nchannels);
        return build<micro_Qdtm>(get<number_of_samples>(Qdt), get<min_P>(Qdt), std::move(t_micro_P), std::move(t_P_to_Pmean),
                                 std::move(t_gmean_i), std::move(t_gvar_i));
    }

    // Per-Agonist_step micro_Qdt builder for the averaging=2 (joint) algorithms.
    // Mirrors calc_micro_Qdtm_agonist_step but consumes the joint-moment Qdt
    // (gmean_ij / gvar_ij) and produces the joint-moment micro_Qdt.
    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    auto calc_micro_Qdt_agonist_step(FunctionTable& f, const C_Patch_Model& m,
                                      const Agonist_step& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdt>> {
        auto maybe_Qdt = Macro_DMR{}.calc_Qdt_agonist_step(f, m, t_step, fs);
        if (!maybe_Qdt)
            return maybe_Qdt.error();
        auto& Qdt = maybe_Qdt.value();
        auto t_micro_P    = P_to_micro_P(get<P>(Qdt), Nchannels);
        auto t_P_to_Pmean = calc_micro_P_to_P_mean(get<P>(Qdt)().ncols(), Nchannels);
        auto t_gmean_ij   = gij_to_micro_gij(get<gmean_ij>(Qdt), Nchannels);
        auto t_gvar_ij    = gij_to_micro_gij(get<gvar_ij>(Qdt), Nchannels);
        return build<micro_Qdt>(get<number_of_samples>(Qdt), get<min_P>(Qdt),
                                std::move(t_micro_P), std::move(t_P_to_Pmean),
                                std::move(t_gmean_ij), std::move(t_gvar_ij));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdt(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdt>> {
        if constexpr (std::is_same_v<Policy, StabilizerPolicyEnabled>) {
            if constexpr (std::is_same_v<Nothing, decltype(f[Calc_micro_Qdt_step{}])>)
                return calc_micro_Qdt_agonist_step<Policy>(f, m, t_step, fs,Nchannels);
            else
                return f.f(Calc_micro_Qdt_step{}, m, t_step, fs,Nchannels);
        } else {
            return calc_micro_Qdt_agonist_step<Policy>(f, m, t_step, fs,Nchannels);
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtm(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtm>> {
        if constexpr (std::is_same_v<Policy, StabilizerPolicyEnabled>) {
            if constexpr (std::is_same_v<Nothing, decltype(f[Calc_micro_Qdtm_step{}])>)
                return calc_micro_Qdtm_agonist_step<Policy>(f, m, t_step, fs, Nchannels);
            else
                return f.f(Calc_micro_Qdtm_step{}, m, t_step, fs, Nchannels);
        } else {
            return calc_micro_Qdtm_agonist_step<Policy>(f, m, t_step, fs, Nchannels);
        }
    }

   
    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_micro_Qdt(FunctionTable& f, const C_Patch_Model& m, const std::vector<Agonist_step>& t_step,
                  double fs, std::size_t Nchannels) -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdt>> {
        // Single-step short-circuit hits the memoizable single-Agonist_step overload
        // (mirrors Macro_DMR::calc_Qdt at qmodel.h:2890-2891). Multi-step recordings
        // still take the recompute path below; handling those with caching requires
        // a per-Nchannels micro-Qdt composition that doesn't exist yet.
        if (t_step.size() == 1)
            return calc_micro_Qdt<Policy>(f, m, t_step[0], fs, Nchannels);

        auto maybe_Qdt = Macro_DMR{}.calc_Qdt(f, m, t_step, fs);
        if (!maybe_Qdt)
            return maybe_Qdt.error();
        auto& Qdt = maybe_Qdt.value();
        auto t_micro_P = P_to_micro_P(get<P>(Qdt), Nchannels);
        auto t_P_to_Pmean = calc_micro_P_to_P_mean(get<P>(Qdt)().ncols(), Nchannels);
        auto t_gmean_ij = gij_to_micro_gij(get<gmean_ij>(Qdt), Nchannels);
        auto t_gvar_ij = gij_to_micro_gij(get<gvar_ij>(Qdt), Nchannels);
        return build<micro_Qdt>(get<number_of_samples>(Qdt), get<min_P>(Qdt), std::move(t_micro_P), std::move(t_P_to_Pmean),std::move(t_gmean_ij), std::move(t_gvar_ij));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_micro_Qdtg(FunctionTable& f, const C_Patch_Model& m,
                   const std::vector<Agonist_step>& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtg>> {
        // Single-step short-circuit; see notes on calc_micro_Qdt above.
        if (t_step.size() == 1)
            return calc_micro_Qdtg(f, m, t_step[0], fs, Nchannels);

        auto maybe_Qdt = Macro_DMR{}.calc_Qdtg(f, m, t_step, fs);
        if (!maybe_Qdt)
            return maybe_Qdt.error();
        auto& Qdt = maybe_Qdt.value();
        auto t_P = P_to_micro_P(get<P_half>(Qdt), Nchannels);
        auto t_P_to_Pmean = calc_micro_P_to_P_mean(get<P_half>(Qdt)().ncols(), Nchannels);
        auto t_g = g_to_micro_g(get<g>(Qdt), Nchannels);
        return build<micro_Qdtg>(get<number_of_samples>(Qdt), get<min_P>(Qdt), std::move(t_P), std::move(t_P_to_Pmean), std::move(t_g));
    }
   
    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_micro_Qdtm(FunctionTable& f, const C_Patch_Model& m,
                   const std::vector<Agonist_step>& t_step, double fs, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtm>> {
        // Single-step short-circuit; see notes on calc_micro_Qdt above.
        if (t_step.size() == 1)
            return calc_micro_Qdtm<Policy>(f, m, t_step[0], fs, Nchannels);

         auto maybe_Qdt = Macro_DMR{}.calc_Qdtm(f, m, t_step, fs);
        if (!maybe_Qdt) 
            return maybe_Qdt.error();
        auto& Qdt = maybe_Qdt.value();
        auto t_micro_P = P_to_micro_P(get<P>(Qdt), Nchannels);
        auto t_P_to_Pmean = calc_micro_P_to_P_mean(get<P>(Qdt)().ncols(), Nchannels);
        auto t_gmean_i = g_to_micro_g(get<gmean_i>(Qdt), Nchannels);
        auto t_gvar_i = g_to_micro_g(get<gvar_i>(Qdt), Nchannels);
        return build<micro_Qdtm>(get<number_of_samples>(Qdt), get<min_P>(Qdt), std::move(t_micro_P), std::move(t_P_to_Pmean), std::move(t_gmean_i), std::move(t_gvar_i));
        
}


    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdt(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                  double fs, std::size_t Nchannels) -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdt>> {
        return calc_micro_Qdt<Policy>(f, m, t_step(), fs, Nchannels );
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtm(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                   double fs, std::size_t Nchannels) -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtm>> {
        return calc_micro_Qdtm<Policy>(f, m, t_step(), fs, Nchannels);
    }

 template <class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_micro_Qdtg(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                   double fs, std::size_t Nchannels) -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Qdtg>> {
        return calc_micro_Qdtg(f, m, t_step(), fs, Nchannels);
    }




    template <class recursive, class averaging, class variance, class variance_correction,
              class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 (U<C_Patch_State, Algo_State> || U<C_Patch_State, Algo_State_Dynamic>))
    Maybe_error<C_Patch_State> safely_calculate_y_mean_yvar_Pmean_PCov(
        C_Patch_State const& t_prior, C_Qdt const& t_Qdt, C_Patch_Model const& m, C_double const& N,
        const Patch_current& p_y, double fs) const {
        constexpr bool PoissonDif = true;
        using Transf = transformation_type_t<C_Qdt>;

        auto& p_P_mean = get<P_mean>(t_prior());
        auto SmD = get<P_Cov>(t_prior())() - diag(p_P_mean());
        auto& y = p_y.value();
        bool is_y_nan = std::isnan(y);
        auto y_baseline = get<Current_Baseline>(m);
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value();
        Matrix<double> u(p_P_mean().size(), 1, 1.0);
        if constexpr (averaging::value == 0) {
            auto& t_g = get<g>(m);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_g()) + y_baseline());
            auto gSg = getvalue(TranspMult(t_g(), SmD) * t_g());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);

            auto r_y_var = build<y_var>(e);
            if (primitive(gSg) > 0)
                r_y_var() = r_y_var() + N * gSg;

            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        } else if constexpr (averaging::value == 2) {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
            auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
            if constexpr (variance_correction::value) {
                auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
                auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                auto& t_gvar_i = get<gvar_i>(t_Qdt);
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

                auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
                auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

                auto delta_emu = sqr(ms + e / N) - 2.0 / N * sSs;
                if (!std::isfinite(primitive(delta_emu)) || primitive(delta_emu) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

                auto e_mu = e + N * ms0;
                auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) -
                                              N * 0.5 / e_mu * sSg + y_baseline());
                auto zeta = N / (2 * sqr(e_mu) + N * sSs);

                if (!std::isfinite(primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg))) ||
                    primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg)) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto r_y_var = build<y_var>(e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
                auto& t_P = get<P>(t_Qdt);
                auto dy = y - r_y_mean();
                auto chi = dy / r_y_var();
                auto chi2 = dy * chi;

                auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
                auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();

                auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS -
                                                     (chi * zeta * sSg + 0.5 / e_mu) * sS);

                if (!Maybe_r_P_mean)
                    return Maybe_r_P_mean.error();
                auto r_P_mean = std::move(Maybe_r_P_mean.value());

                auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD) + diag(r_P_mean * t_P()) -
                                            (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
                                            (2.0 * N / r_y_var() * zeta * sSg) *
                                                X_plus_XT(TranspMult(sS, gS)) -
                                            (N / r_y_var()) * XTX(gS));

                if (!all_Probability_elements(primitive(r_P_mean)) ||
                    !all_Covariance_elements(primitive(r_P_cov())))
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                else {
                    auto r_macro_algo = macror_algorithm(
                        ToString(MacroR2<recursive, averaging, variance, variance_correction>{}));
                    return std::tuple(std::move(r_y_mean), std::move(r_y_var), std::move(r_P_mean),
                                      std::move(r_P_cov), std::move(chi2), std::move(r_macro_algo));
                }

            } else {
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
                if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, uses_averaging_aproximation<1>, variance, variance_correction>(
                        t_prior, t_Qdt, m, N, p_y, fs);
                auto r_y_mean =
                    build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
                if constexpr (PoissonDif)
                    e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
                else
                    e = e + get<Proportional_Noise>(m).value() * abs(y);
                auto r_y_var = build<y_var>(e + N * gSg);
                if constexpr (variance::value) {
                    auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                    if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                        r_y_var() = r_y_var() + N * ms;
                    } else {
                        return safely_calculate_Algo_Pmean_Pcov<recursive, averaging,
                                                                uses_variance_aproximation<false>,
                                                                variance_correction>(
                            is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                            p_y, fs, SmD);
                    }
                }
                return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                        variance_correction>(
                    is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y,
                    fs, SmD);
            }
        } else /* if constexpr (averaging::value == 1) */ {
            auto& t_gmean_i = [&t_Qdt](){ if constexpr(averaging::value>0) {return get<gmean_i>(t_Qdt);} else {return get<g>(t_Qdt);} }();    
            auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i());
            if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                return safely_calculate_y_mean_yvar_Pmean_PCov<
                    recursive, uses_averaging_aproximation<0>, uses_variance_aproximation<false>,
                    uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N, p_y,
                                                                         fs);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);
            auto r_y_var = build<y_var>(e + N * gSg);
            if constexpr (variance::value) {
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                if (std::isfinite(primitive(ms)) && primitive(ms) > 0) {
                    r_y_var() = r_y_var() + N * ms;
                } else {
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, uses_variance_aproximation<false>,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                }
            }
            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        }
    }

    template <class recursive, class averaging, class variance, class variance_correction,
              class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction>)
    Maybe_error<std::tuple<Transfer_Op_to<C_Qdt, y_mean>, Transfer_Op_to<C_Qdt, y_var>,
                           Transfer_Op_to<C_Qdt, P_mean>, Transfer_Op_to<C_Qdt, P_Cov>,
                           Transfer_Op_to<C_Qdt, double>, macror_algorithm>>
        safely_calculate_y_mean_yvar_Pmean_PCov(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                                C_Patch_Model const& m, C_double const& N,
                                                const Patch_current& p_y, double fs) const {
        constexpr bool PoissonDif = true;
        using Transf = transformation_type_t<C_Qdt>;

        auto& p_P_mean = get<P_mean>(t_prior);
        auto SmD = get<P_Cov>(t_prior)() - diag(p_P_mean());
        auto& y = p_y.value();
        bool is_y_nan = std::isnan(y);
        auto y_baseline = get<Current_Baseline>(m);
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value();
        Matrix<double> u(p_P_mean().size(), 1, 1.0);
        if constexpr (averaging::value == 0) {
            auto& t_g = get<g>(m);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_g()) + y_baseline());
            auto gSg = getvalue(TranspMult(t_g(), SmD) * t_g());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);

            auto r_y_var = build<y_var>(e);
            if (primitive(gSg) > 0)
                r_y_var() = r_y_var() + N * gSg;

            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        } else if constexpr (averaging::value == 2) {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
            auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
            if constexpr (variance_correction::value) {
                auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
                auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                auto& t_gvar_i = get<gvar_i>(t_Qdt);
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

                auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
                auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

                auto delta_emu = sqr(ms + e / N) - 2.0 / N * sSs;
                if (!std::isfinite(primitive(delta_emu)) || primitive(delta_emu) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

                auto e_mu = e + N * ms0;
                auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) -
                                              N * 0.5 / e_mu * sSg + y_baseline());
                auto zeta = N / (2 * sqr(e_mu) + N * sSs);

                if (!std::isfinite(primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg))) ||
                    primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg)) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto r_y_var = build<y_var>(e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
                auto& t_P = get<P>(t_Qdt);
                auto dy = y - r_y_mean();
                auto chi = dy / r_y_var();
                auto chi2 = dy * chi;

                auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
                auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();

                auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS -
                                                     (chi * zeta * sSg + 0.5 / e_mu) * sS);

                if (!Maybe_r_P_mean)
                    return Maybe_r_P_mean.error();
                auto r_P_mean = std::move(Maybe_r_P_mean.value());

                auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD) + diag(r_P_mean * t_P()) -
                                            (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
                                            (2.0 * N / r_y_var() * zeta * sSg) *
                                                X_plus_XT(TranspMult(sS, gS)) -
                                            (N / r_y_var()) * XTX(gS));

                if (!all_Probability_elements(primitive(r_P_mean)) ||
                    !all_Covariance_elements(primitive(r_P_cov())))
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                else {
                    auto r_macro_algo = macror_algorithm(
                        ToString(MacroR2<recursive, averaging, variance, variance_correction>{}));
                    return std::tuple(std::move(r_y_mean), std::move(r_y_var), std::move(r_P_mean),
                                      std::move(r_P_cov), std::move(chi2), std::move(r_macro_algo));
                }

            } else {
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
                if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, uses_averaging_aproximation<1>, variance, variance_correction>(
                        t_prior, t_Qdt, m, N, p_y, fs);
                auto r_y_mean =
                    build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
                if constexpr (PoissonDif)
                    e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
                else
                    e = e + get<Proportional_Noise>(m).value() * abs(y);
                auto r_y_var = build<y_var>(e + N * gSg);
                if constexpr (variance::value) {
                    auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                    if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                        r_y_var() = r_y_var() + N * ms;
                    } else {
                        return safely_calculate_Algo_Pmean_Pcov<recursive, averaging,
                                                                uses_variance_aproximation<false>,
                                                                variance_correction>(
                            is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                            p_y, fs, SmD);
                    }
                }
                return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                        variance_correction>(
                    is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y,
                    fs, SmD);
            }
        } else /* if constexpr (averaging::value == 1) */ {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i());
            if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                return safely_calculate_y_mean_yvar_Pmean_PCov<
                    recursive, uses_averaging_aproximation<0>, uses_variance_aproximation<false>,
                    uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N, p_y,
                                                                         fs);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);
            auto r_y_var = build<y_var>(e + N * gSg);
            if constexpr (variance::value) {
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                if (std::isfinite(primitive(ms)) && primitive(ms) > 0) {
                    r_y_var() = r_y_var() + N * ms;
                } else {
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, uses_variance_aproximation<false>,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                }
            }
            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        }
    }

    template <class C_y_var>
    auto calculate_logL(bool y_is_nan, C_y_var const& r_y_var, auto const& chi2,
                         auto& m) const -> Transfer_Op_to<C_y_var, logL> {
        using DX = var::dx_of_dfdx_t<C_y_var>;
        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<logL>(std::move(base));
        }
        if (get<Proportional_Noise>(m).value() == 0) {
            return build<logL>(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
                               0.5 * chi2());
        } else {
            return build<logL>(0.5 * chi2() - log(var::Poisson_noise_normalization(
                                                  primitive(r_y_var()),
                                                  primitive(get<Proportional_Noise>(m).value()))));
        }
    }

	    template <class C_y_var>
	    auto calculate_elogL(bool y_is_nan, C_y_var const& r_y_var, auto& m) const {
	        using DX = var::dx_of_dfdx_t<C_y_var>;
	        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<elogL>(std::move(base));
        }
	        if (get<Proportional_Noise>(m).value() == 0) {
	            return build<elogL>(-0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5);
	        } else {
	            const auto pn_value = get<Proportional_Noise>(m).value();
	            if constexpr (var::is_derivative_v<std::decay_t<decltype(pn_value)>>) {
	                return build<elogL>(var::Poisson_noise_expected_logL(r_y_var() ,
	                                                                    pn_value));
	            } else {
	                auto pn = var::init_with_dx<DX>(pn_value, dx);
	                return build<elogL>(
	                    var::Poisson_noise_expected_logL(r_y_var() , pn));
	            }
	        }
	    }

    auto calculate_trust_coefficient(Matrix<double> const& t_pmean, Matrix<double> const& d,
                                     double factor) const {
        auto alfa = 1.0;
        for (std::size_t i = 0; i < d.size(); ++i) {
            double d_i = d[i];
            double p_i = t_pmean[i];
            if (d_i > 0) {
                double alfa_i = (1.0 - p_i) / d_i;
                alfa = std::min(alfa_i, alfa);
            } else if (d_i < 0) {
                double alfa_i = -(p_i) / d_i;
                alfa = std::min(alfa_i, alfa);
            }
        }
        if (alfa < 1)
            alfa = alfa * factor;
        return trust_coefficient(alfa);
    }


    template <bool dynamic, class averaging, class C_micro_Patch_State, class C_micro_Qdt,
              class C_micro_Patch_Model, class C_double>
        requires(uses_averaging_aproximation_c<averaging> && (U<C_micro_Patch_State, micro_Patch_State>))
    auto safely_calculate_micro_Algo_State_non_recursive(C_micro_Patch_State const& t_micro_prior, C_micro_Qdt const& t_micro_Qdt,
                                                   C_micro_Patch_Model const& m,
                                                   const Patch_current& p_y, double fs, C_double const& N) const
        -> Maybe_error<Transfer_Op_to<
            C_micro_Patch_State, std::conditional_t<dynamic, micro_Algo_State_Dynamic, micro_Algo_State>>> {
        constexpr bool PoissonDif = true;
        auto& y = p_y.value();

        auto const& t_micro_P = [&t_micro_Qdt]() { 
            if constexpr(averaging::value > 0) 
               {return get<micro_P>(t_micro_Qdt);} 
            else {return get<micro_P_half>(t_micro_Qdt); }}
            ();
        auto  p_micro_P_state = get<micro_P_state>(t_micro_prior());
        if constexpr (averaging::value == 0) {
            p_micro_P_state() = p_micro_P_state() * t_micro_P();
        }


        auto y_baseline = get<Current_Baseline>(m);
        
        auto& t_g_st = [&t_micro_Qdt, &m]() -> decltype(auto) {
            if constexpr (averaging::value == 0)
                return get<micro_g>(t_micro_Qdt);
            else if constexpr (averaging::value == 1)
                return get<micro_gmean_i>(t_micro_Qdt);
            else if constexpr (averaging::value == 2)
                return get<micro_gmean_ij>(t_micro_Qdt);
            else static_assert("invalid averaging");    
        }();
       // avg=0 has no per-microstate variance — build a zero placeholder of the same
       // shape as t_g_st. Derive the zero from t_g_st itself so that under
       // derivative parameters it remains Derivative<micro_g, P> (with both
       // primitive and derivative parts zeroed), keeping var::zip's
       // all-Derivatives constraint satisfied.
       auto t_gvar_zero_for_avg0 = [&t_g_st]() {
            if constexpr (averaging::value == 0)
                return build<micro_g>(t_g_st() * 0.0);
            else
                return 0;
        }();
       auto& t_gvar_st = [&t_micro_Qdt, &t_gvar_zero_for_avg0]() -> decltype(auto) {
            if constexpr (averaging::value == 0)
                return (t_gvar_zero_for_avg0);
            else if constexpr (averaging::value == 1)
                return get<micro_gvar_i>(t_micro_Qdt);
            else if constexpr (averaging::value == 2)
                return get<micro_gvar_ij>(t_micro_Qdt);
            else static_assert("invalid averaging");
        }();
       // For avg<2, return a reference to the existing micro_P_state Var.
       // For avg=2, build a new micro_P_state by-value from diag(prior)·micro_P.
       // Use auto&& to bind to both the lvalue ref (avg<2) and the rvalue (avg=2),
       // while keeping the consumer call shape `t_P_state()` uniform.
       // avg=0 evaluates evidence at the interval midpoint, so use the half-stepped
       // p_micro_P_state (start→middle). avg=1 uses the start-of-interval prior with
       // gmean_i/gvar_i (full-interval averages). avg=2 builds the joint (start,end).
       auto&& t_P_state = [&t_micro_Qdt, &t_micro_prior, &p_micro_P_state]() -> decltype(auto) {
            if constexpr (averaging::value == 0)
                return (p_micro_P_state);
            else if constexpr (averaging::value == 1)
                return get<micro_P_state>(t_micro_prior());
            else if constexpr (averaging::value == 2)
                return build<micro_P_state>(diag(get<micro_P_state>(t_micro_prior())()) * get<micro_P>(t_micro_Qdt)());
            else static_assert("invalid averaging");
        }();

        // --- Mixture-of-Gaussians per-microstate predictive moments ---
        // V_i = per-microstate predicted current (baseline-free).  σ²_i = gvar_i + e.
        // For averaging<2 these are M-shaped; for averaging=2 they are M×M (per (start,end) pair).
        // Sum-via-elemMult-then-sum works uniformly across both shapes.
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_micro_Qdt).value() +
                 get<Pink_Noise>(m).value();

        
        auto V_sq = elemMult(t_g_st(), t_g_st());
        auto mean_v    = sum(elemMult(t_P_state(), t_g_st()));    // Σ p_i V_i
        auto mean_v_sq = sum(elemMult(t_P_state(), V_sq));        // Σ p_i V_i²
        auto var_v     = mean_v_sq - mean_v * mean_v;             // Var_p(V)
        auto within_v  = sum(elemMult(t_P_state(), t_gvar_st())); // Σ p_i gvar_i

        auto r_y_mean = build<y_mean>(mean_v + y_baseline());
        auto r_y_var  = build<y_var>(e + within_v + var_v);

        // --- Macro-slot diagnostics (moment-matched, mirrors macro path) ---
        using std::sqrt;
        auto dy       = y - r_y_mean();
        auto r_r_std  = build<r_std>(dy / sqrt(r_y_var()));
        auto chi2     = build<Chi2>(dy * dy / r_y_var());

        // --- Micro-slot diagnostics (Rosenblatt PIT residual) ---
        // F_y = Σ_i p_i · Φ((y − μ_i)/σ_i)  with  μ_i = V_i + baseline,  σ²_i = gvar_i + e.
        // Under correct specification micro_r_std is i.i.d. N(0,1), so micro_Chi2 ~ χ²₁.
        auto y_corr = y - y_baseline();
        auto cdf_per_micro = zip(
            [&y_corr, &e](auto const& V_i, auto const& gvar_i) {
                using std::sqrt;
                return normal_cdf((y_corr - V_i) / sqrt(gvar_i + e));
            },
            t_g_st(), t_gvar_st());
        auto F_y = sum(elemMult(t_P_state(), cdf_per_micro));
        auto r_micro_r_std = build<micro_r_std>(normal_quantile(F_y));
        auto r_micro_chi2  = build<micro_Chi2>(r_micro_r_std() * r_micro_r_std());

        // === Per-step logL / elogL — mixture forms stored in the algo state ===
        //   logL          : log Σ_n π(n) · L_n          (Bayesian sum of contribs)
        //   elogL         : -½ - ½ E_π[log(2π σ²_n)]    (per-state Gaussian
        //                   entropy averaged by prior — mixture analog of macro
        //                   `-½log(2π σ²) - ½`)
        //   gaussian_logL : moment-match form via calculate_logL(y_var, Chi2, m).
        //                   Diagnostic comparator for the negentropy/propagation
        //                   decomposition.
        // y may be NaN here (gap observations); mirror calculate_logL's NaN guard.
        bool y_is_nan = std::isnan(y);
        auto r_logL = [&]() {
            if (y_is_nan) {
                using DX = var::dx_of_dfdx_t<decltype(r_y_var)>;
                auto const& dx = var::get_dx_of_dfdx(r_y_var);
                return build<logL>(var::init_with_dx<DX>(0.0, dx));
            }
            // log Σ_n π(n) · (1/√(2π σ²_n)) · exp(-½ (y − μ_n)²/σ²_n).
            // Per-microstate setup mirrors the Rosenblatt PIT block above —
            // same y_corr, e, t_g_st, t_gvar_st, t_P_state.
            auto t_posterior_y = zip(
                [&y_corr, &e](auto const& prior, auto const& g_mean, auto const& g_var) {
                    using std::sqrt;
                    using std::exp;
                    return prior * exp(-0.5 * (y_corr - g_mean) * (y_corr - g_mean) / (g_var + e)) /
                           sqrt(2 * std::numbers::pi * (g_var + e));
                },
                t_P_state(), t_g_st(), t_gvar_st());
            using std::log;
            return build<logL>(log(sum(t_posterior_y)));
        }();
        auto r_elogL = [&]() {
            if (y_is_nan) {
                using DX = var::dx_of_dfdx_t<decltype(r_y_var)>;
                auto const& dx = var::get_dx_of_dfdx(r_y_var);
                return build<elogL>(var::init_with_dx<DX>(0.0, dx));
            }
            auto log_2pi_sigma2_per_micro = zip(
                [&e](auto const& g_var) {
                    using std::log;
                    return log(2 * std::numbers::pi * (g_var + e));
                },
                t_gvar_st());
            return build<elogL>(
                -0.5 - 0.5 * sum(elemMult(t_P_state(), log_2pi_sigma2_per_micro)));
        }();
        auto r_gaussian_logL = calculate_logL(y_is_nan, r_y_var, chi2, m);

        auto alfa = trust_coefficient(1.0);

        // --- Project to macro-shape P_mean / P_Cov ---
        // Mirrors the macro non-recursive sister at qmodel.h:3735-3789: apply one
        // Markov step to the microstate distribution, then take all macro-shape
        // moments from the post-step distribution.
        //   r_micro_P_state = p_micro_P_state · t_micro_P
        //     (for avg=0, p_micro_P_state was half-advanced at line 893 and t_micro_P is
        //      the second half-step → full step from prior; for avg>0, full step.)
        // Π = micro_P_state_to_P_mean : (M_micro, k_states); each row is the
        // macrostate-fraction vector for that microstate.
        // Random vector X = macrostate-fraction-vector (random because microstate is random).
        //   E[X]   = probs · Π                              →  P_mean   (1, k_states)
        //   E[XXᵀ] = Πᵀ · diag(probs) · Π   = AT_D_A(Π, DiagPosDet(probs))
        //   Cov(X) = E[XXᵀ] − E[X]ᵀE[X]     = AT_D_A(...) − XTX(P_mean)
        // P_Cov is stored as bare centered covariance (rows sum to 0), matching the macro
        // convention. Micro uses AT_D_A (exact second moment from the full microstate
        // distribution) where macro uses AT_B_A on a moment-closed SmD (approximation) — same
        // structural shape, exact in the micro universe.  Validation goes through
        // to_Covariance_Probability directly; no "+ diag(P_mean)" correction is needed because
        // AT_D_A − XTX already produces bare cov.
        auto Maybe_r_micro_P_state = to_Probability(p_micro_P_state() * t_micro_P());
        if (!Maybe_r_micro_P_state.valid())
            return Maybe_r_micro_P_state.error();
        auto r_micro_P_state = build<micro_P_state>(std::move(Maybe_r_micro_P_state.value()));

        auto& Pi = get<micro_P_state_to_P_mean>(t_micro_Qdt)();
        auto& probs = r_micro_P_state();

        auto Maybe_r_P_mean = to_Probability(probs * Pi);
        if (!Maybe_r_P_mean.valid())
            return Maybe_r_P_mean.error();
        auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

        // AT_D_A(Π, diag(probs)) = Πᵀ·diag(probs)·Π = E[YYᵀ] (M-form, rows sum to μ).
        // Subtracting XTX(μ) = μμᵀ gives bare centered covariance directly (rows sum to 0).
        // No `+ diag(μ)` — that's the macro-path correction needed when r_P_cov was constructed
        // *without* the diag(μ) term; here the construction already produces full bare_cov.
        auto r_P_cov = build<P_Cov>(
            AT_D_A(Pi, diagpos(probs)) - XTX(r_P_mean()));
        auto Maybe_r_P_cov = to_Covariance_Probability(r_P_cov());
        if (!Maybe_r_P_cov.valid())
            return Maybe_r_P_cov.error();
        r_P_cov() = std::move(Maybe_r_P_cov.value());

        if constexpr (!dynamic) {
            Transfer_Op_to<C_micro_Patch_State, micro_Algo_State> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out())  = std::move(r_y_var);
            get<r_std>(out())  = std::move(r_r_std);
            get<Chi2>(out())   = std::move(chi2);
            get<trust_coefficient>(out()) = alfa;

            get<micro_r_std>(out()) = std::move(r_micro_r_std);
            get<micro_Chi2>(out())  = std::move(r_micro_chi2);
            get<logL>(out()) = std::move(r_logL);
            get<elogL>(out()) = std::move(r_elogL);
            get<gaussian_logL>(out())() = std::move(r_gaussian_logL)();
            get<micro_P_state>(out()) = std::move(r_micro_P_state);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())()  = std::move(r_P_cov());
            return out;
        } else {
            Transfer_Op_to<C_micro_Patch_State, micro_Algo_State_Dynamic> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out())  = std::move(r_y_var);
            get<r_std>(out())  = std::move(r_r_std);
            get<Chi2>(out())   = std::move(chi2);
            get<trust_coefficient>(out()) = alfa;

            get<micro_r_std>(out()) = std::move(r_micro_r_std);
            get<micro_Chi2>(out())  = std::move(r_micro_chi2);
            get<logL>(out()) = std::move(r_logL);
            get<elogL>(out()) = std::move(r_elogL);
            get<gaussian_logL>(out())() = std::move(r_gaussian_logL)();
            // Predictive (pre-observation) slots — mirrors macro non-recursive at
            // qmodel.h:3815-3816 which fills only the _t2_y0 pair for this path.
            get<P_mean_t2_y0>(out())() = std::move(r_P_mean());
            get<P_Cov_t2_y0>(out())()  = std::move(r_P_cov());
            get<micro_P_state_t2_y0>(out())() = std::move(r_micro_P_state());
            return out;
        }
    }



template <bool dynamic, class averaging,  class C_micro_Patch_State, class C_micro_Qdt,
              class C_micro_Patch_Model, class C_double>

        requires(uses_averaging_aproximation_c<averaging>  && U<C_micro_Patch_State, micro_Patch_State>)
    auto safely_calculate_micro_Algo_State_recursive(C_micro_Patch_State const& t_micro_prior, 
        C_micro_Qdt const& t_micro_Qdt,
                                               C_micro_Patch_Model const& m, 
                                               C_double const& N,
                                               const Patch_current& p_y, double fs) const
        -> Maybe_error<Transfer_Op_to<
            C_micro_Patch_State, std::conditional_t<dynamic, micro_Algo_State_Dynamic, micro_Algo_State>>> {
        auto& y = p_y.value();
        if (std::isnan(y)) {
            return safely_calculate_micro_Algo_State_non_recursive<dynamic, averaging>(
                t_micro_prior, t_micro_Qdt, m, p_y, fs, N);
        }

        // === Setup (mirrors non-recursive sister) ===
        auto const& t_micro_P = [&t_micro_Qdt]() -> decltype(auto) {
            if constexpr (averaging::value == 0) return get<micro_P_half>(t_micro_Qdt);
            else                                  return get<micro_P>(t_micro_Qdt);
        }();
        auto p_micro_P_state = get<micro_P_state>(t_micro_prior());
        if constexpr (averaging::value == 0) {
            p_micro_P_state() = p_micro_P_state() * t_micro_P();
        }

        auto y_baseline = get<Current_Baseline>(m);

        auto& t_g_st = [&t_micro_Qdt]() -> decltype(auto) {
            if constexpr (averaging::value == 0)      return get<micro_g>(t_micro_Qdt);
            else if constexpr (averaging::value == 1) return get<micro_gmean_i>(t_micro_Qdt);
            else if constexpr (averaging::value == 2) return get<micro_gmean_ij>(t_micro_Qdt);
            else static_assert("invalid averaging");
        }();
        // avg=0 has no per-microstate variance — build a zero placeholder of the same
        // shape as t_g_st so the dispatch can return a uniform reference. We
        // derive the zero from t_g_st itself (rather than constructing a plain
        // Matrix<double>) so that under derivative parameters the result is
        // Derivative<micro_g, P> with primitive=zeros and derivative=zeros, not
        // a plain Matrix that would break the all-Derivatives constraint on
        // var::zip.
        auto t_gvar_zero_for_avg0 = [&t_g_st]() {
            if constexpr (averaging::value == 0)
                return build<micro_g>(t_g_st() * 0.0);
            else
                return 0;
        }();
        auto& t_gvar_st = [&t_micro_Qdt, &t_gvar_zero_for_avg0]() -> decltype(auto) {
            if constexpr (averaging::value == 0)      return (t_gvar_zero_for_avg0);
            else if constexpr (averaging::value == 1) return get<micro_gvar_i>(t_micro_Qdt);
            else if constexpr (averaging::value == 2) return get<micro_gvar_ij>(t_micro_Qdt);
            else static_assert("invalid averaging");
        }();
        // avg=0 evaluates evidence at the interval midpoint, so use the half-stepped
        // p_micro_P_state (start→middle). avg=1 uses the start-of-interval prior with
        // gmean_i/gvar_i (full-interval averages). avg=2 builds the joint (start,end).
        auto&& t_P_state = [&t_micro_Qdt, &t_micro_prior, &p_micro_P_state]() -> decltype(auto) {
            if constexpr (averaging::value == 0)
                return (p_micro_P_state);
            else if constexpr (averaging::value == 1)
                return get<micro_P_state>(t_micro_prior());
            else
                return build<micro_P_state>(diag(get<micro_P_state>(t_micro_prior())()) * get<micro_P>(t_micro_Qdt)());
        }();

        // === Per-microstate predictive moments (mirrors non-recursive) ===
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_micro_Qdt).value() +
                 get<Pink_Noise>(m).value();

        auto V_sq = elemMult(t_g_st(), t_g_st());
        auto mean_v    = sum(elemMult(t_P_state(), t_g_st()));
        auto mean_v_sq = sum(elemMult(t_P_state(), V_sq));
        auto var_v     = mean_v_sq - mean_v * mean_v;
        auto within_v  = sum(elemMult(t_P_state(), t_gvar_st()));

        auto r_y_mean = build<y_mean>(mean_v + y_baseline());
        auto r_y_var  = build<y_var>(e + within_v + var_v);

        // === Macro-slot diagnostics (moment-matched, mirrors macro path) ===
        using std::sqrt;
        auto dy       = y - r_y_mean();
        auto r_r_std  = build<r_std>(dy / sqrt(r_y_var()));
        auto chi2     = build<Chi2>(dy * dy / r_y_var());

        // === Micro-slot diagnostics (Rosenblatt PIT residual) ===
        // F_y = Σᵢ pᵢ · Φ((y − μᵢ)/σᵢ) — same as non-recursive.
        auto y_corr = y - y_baseline();
        auto cdf_per_micro = zip(
            [&y_corr, &e](auto const& V_i, auto const& gvar_i) {
                using std::sqrt;
                return normal_cdf((y_corr - V_i) / sqrt(gvar_i + e));
            },
            t_g_st(), t_gvar_st());
        auto F_y = sum(elemMult(t_P_state(), cdf_per_micro));
        auto r_micro_r_std = build<micro_r_std>(normal_quantile(F_y));
        auto r_micro_chi2  = build<micro_Chi2>(r_micro_r_std() * r_micro_r_std());

        // === Bayesian update: per-microstate likelihood-weighted posterior ===
        // Build the pure Gaussian likelihood (no prior factor); Bayes_Rule
        // handles the prior multiplication, normalization and evidence.
        //   avg<2: 1D — per-microstate i.
        //   avg=2: 2D — per-pair (start, end).
        auto t_likelihood = zip(
            [&y_corr, &e](auto const& g_mean, auto const& g_var) {
                using std::sqrt;
                using std::exp;
                return exp(-0.5 * (y_corr - g_mean) * (y_corr - g_mean) / (g_var + e)) /
                       sqrt(2 * std::numbers::pi * (g_var + e));
            },
            t_g_st(), t_gvar_st());

        // Bayes_Rule returns (posterior, evidence) where evidence is the
        // marginal likelihood Σ priorᵢ · Lᵢ. In channel-kinetics terms,
        // log(evidence) is the per-step contribution to the parameter
        // log-likelihood — that's where the Bayesian → channel-kinetics
        // renaming happens (the `build<logL>(log(evidence))` line below).
        // Variant tag for breadcrumbs: encodes the compile-time choice of this
        // recursive-Bayes specialization. Single string reused across the
        // Maybe-returning sites in this function.
        std::string const variant_tag =
            std::string("[recursive,avg=") + std::to_string(averaging::value) + "]";

        auto Maybe_bayes = Bayes_Rule(t_P_state(), t_likelihood);
        if (!Maybe_bayes)
            return error_message(variant_tag + " | Bayes_Rule | " + Maybe_bayes.error()());
        auto bayes_pair  = std::move(Maybe_bayes.value());  // own posterior + evidence
        auto& bayes_post = bayes_pair.first;
        auto& evidence   = bayes_pair.second;

        // === Per-step logL / elogL — mixture form, stored in the algo state ===
        //   logL  : log evidence  =  log Σ_n π(n) · L_n
        //   elogL : Σ_n π(n) · (-½ - ½ log(2π σ²_n))  =  -½ - ½ E_π[log(2π σ²_n)]
        //           (per-state Gaussian-entropy averaged by prior — the natural
        //           mixture analog of the macro `-½log(2π σ²) - ½`).
        //   gaussian_logL : moment-match form via calculate_logL(y_var, Chi2, m).
        //                   Diagnostic comparator for the negentropy / propagation
        //                   decomposition.
        // y is non-NaN here (recursive function returns to non_recursive on NaN at
        // line 1052); no NaN guard needed in this branch.
        using std::log;
        auto r_logL = build<logL>(log(evidence));
        auto log_2pi_sigma2_per_micro = zip(
            [&e](auto const& g_var) {
                using std::log;
                return log(2 * std::numbers::pi * (g_var + e));
            },
            t_gvar_st());
        auto r_elogL = build<elogL>(
            -0.5 - 0.5 * sum(elemMult(t_P_state(), log_2pi_sigma2_per_micro)));
        auto r_gaussian_logL = calculate_logL(false, r_y_var, chi2, m);

        // Posterior microstate distribution at the END of the interval:
        //   avg<2: posterior · t_micro_P                 (Bayes then Markov)
        //   avg=2: uᵀ · posterior                        (marginalize joint over start)
        // Both inputs to to_Probability now genuinely sum to 1 by construction
        // (posterior is normalized, t_micro_P is row-stochastic, Π is row-
        // stochastic), so to_Probability acts as a pure canary here.
        auto Maybe_r_micro_P_state = [&bayes_post, &t_micro_P]() {
            if constexpr (averaging::value == 2) {
                Matrix<double> uT(1ul, bayes_post.ncols(), 1.0);
                return to_Probability(uT * bayes_post);
            } else {
                return to_Probability(bayes_post * t_micro_P());
            }
        }();
        if (!Maybe_r_micro_P_state) {
            // Two distinct call sites collapse here; label which.
            std::string const site = (averaging::value == 2)
                ? "marginalize-joint"
                : "Markov-step-on-posterior";
            return error_message(variant_tag + " | " + site + " | " +
                                  Maybe_r_micro_P_state.error()());
        }
        auto r_micro_P_state = build<micro_P_state>(std::move(Maybe_r_micro_P_state.value()));

        if (!all_Probability_elements(primitive(r_micro_P_state())))
            return error_message(variant_tag + " | micro_P_state-post-update-elements-check | "
                                  "non-Probability cell after Bayes+Markov update");

        // === Project to macro shape (same pattern as non-recursive, but using the
        // post-observation r_micro_P_state) ===
        auto& Pi = get<micro_P_state_to_P_mean>(t_micro_Qdt)();
        auto& probs = r_micro_P_state();

        auto Maybe_r_P_mean = to_Probability(probs * Pi);
        if (!Maybe_r_P_mean.valid())
            return error_message(variant_tag + " | macro-shape:P_mean projection (probs·Π) | " +
                                  Maybe_r_P_mean.error()());
        auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

        // AT_D_A − XTX = bare centered covariance directly (rows sum to 0). No `+ diag(μ)`.
        auto r_P_cov = build<P_Cov>(
            AT_D_A(Pi, diagpos(probs)) - XTX(r_P_mean()));
        auto Maybe_r_P_cov = to_Covariance_Probability(r_P_cov());
        if (!Maybe_r_P_cov.valid())
            return error_message(variant_tag + " | macro-shape:P_Cov(bare) | " +
                                  Maybe_r_P_cov.error()());
        r_P_cov() = std::move(Maybe_r_P_cov.value());

        auto alfa = trust_coefficient(1.0);

        if constexpr (!dynamic) {
            Transfer_Op_to<C_micro_Patch_State, micro_Algo_State> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out())  = std::move(r_y_var);
            get<r_std>(out())  = std::move(r_r_std);
            get<Chi2>(out())   = std::move(chi2);
            get<trust_coefficient>(out()) = alfa;

            get<micro_r_std>(out()) = std::move(r_micro_r_std);
            get<micro_Chi2>(out())  = std::move(r_micro_chi2);
            get<logL>(out()) = std::move(r_logL);
            get<elogL>(out()) = std::move(r_elogL);
            get<gaussian_logL>(out())() = std::move(r_gaussian_logL)();
            get<micro_P_state>(out()) = std::move(r_micro_P_state);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())()  = std::move(r_P_cov());
            return out;
        } else {
            Transfer_Op_to<C_micro_Patch_State, micro_Algo_State_Dynamic> out;

            // === Common scalar diagnostics ===
            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out())  = std::move(r_y_var);
            get<r_std>(out())  = std::move(r_r_std);
            get<Chi2>(out())   = std::move(chi2);
            get<trust_coefficient>(out()) = alfa;
            get<micro_r_std>(out()) = std::move(r_micro_r_std);
            get<micro_Chi2>(out())  = std::move(r_micro_chi2);
            get<logL>(out()) = std::move(r_logL);
            get<elogL>(out()) = std::move(r_elogL);
            get<gaussian_logL>(out())() = std::move(r_gaussian_logL)();

            // === Helpers: project a microstate distribution / joint to macro shape ===
            // Single distribution q (1, M_micro) → (P_mean (1, k), P_Cov (k, k) in stored form).
            auto project_to_macro =
                [&Pi, &variant_tag](Matrix<double> const& q)
                -> Maybe_error<std::pair<Matrix<double>, SymmetricMatrix<double>>> {
                auto Maybe_pm = to_Probability(q * Pi);
                if (!Maybe_pm.valid())
                    return error_message(variant_tag + " | project_to_macro:P_mean (q·Π) | " +
                                          Maybe_pm.error()());
                auto pm = std::move(Maybe_pm.value());
                // AT_D_A − XTX = bare centered covariance directly (rows sum to 0).
                auto raw = AT_D_A(Pi, diagpos(q)) - XTX(pm);
                auto Maybe_pc = to_Covariance_Probability(raw);
                if (!Maybe_pc.valid())
                    return error_message(variant_tag + " | project_to_macro:P_Cov(bare) | " +
                                          Maybe_pc.error()());
                return std::make_pair(std::move(pm), std::move(Maybe_pc.value()));
            };
            // Joint (M, M) → cross-second-moment Πᵀ J Π (k, k).
            auto project_joint_to_macro =
                [&Pi](Matrix<double> const& J) -> Matrix<double> {
                return TranspMult(Pi, J * Pi);
            };

            // === Predictive (no Bayes) microstate distribution at end of interval ===
            auto Maybe_pred = to_Probability(p_micro_P_state() * t_micro_P());
            if (!Maybe_pred)
                return error_message(variant_tag + " | dynamic:predictive Markov | " +
                                      Maybe_pred.error()());
            auto pred_micro_P_state = std::move(Maybe_pred.value());
            auto Maybe_pred_macro = project_to_macro(pred_micro_P_state);
            if (!Maybe_pred_macro)
                return error_message("dynamic:predictive macro-projection | " +
                                      Maybe_pred_macro.error()());
            auto& [pred_P_mean, pred_P_cov] = Maybe_pred_macro.value();

            if constexpr (averaging::value == 2) {
                // === IR mode: integrated current over interval [0, t] ===
                // Joint prior  P(start=i, end=j) = p_micro_P_state(i)·P_micro(i,j)  = t_P_state.
                // Joint posterior P(start, end | y) = bayes_post                     (from Bayes_Rule).
                // _t11_y0: full-step Markov-stepped predictive (= end-state marginal of joint prior)
                // _t10_y1: start-state posterior marginal (instant Bayes, no Markov)
                // _t20_y1: end-state posterior marginal = r_micro_P_state (already projected)
                // _0t_y0/y1: joint micro distributions; macro projection = Πᵀ J Π
                Matrix<double> uT(1ul, p_micro_P_state().size(), 1.0);

                Matrix<double> joint_prior = t_P_state();   // copy: used for assertions and projection
                auto& joint_posterior = bayes_post;          // already normalized

                auto Maybe_start_post = to_Probability(MultTransp(uT, joint_posterior));
                if (!Maybe_start_post)
                    return error_message(variant_tag + " | dynamic:start-state-marginal-of-joint | " +
                                          Maybe_start_post.error()());
                auto start_micro_P_state_post = std::move(Maybe_start_post.value());
                auto Maybe_start_macro = project_to_macro(start_micro_P_state_post);
                if (!Maybe_start_macro)
                    return error_message("dynamic:start-state macro-projection | " +
                                          Maybe_start_macro.error()());

                // _t11_y0 (end-state Markov-stepped predictive)
                get<micro_P_state_t11_y0>(out())() = pred_micro_P_state;
                get<P_mean_t11_y0>(out())()        = pred_P_mean;
                get<P_Cov_t11_y0>(out())()         = pred_P_cov;

                // _t10_y1 (start-state instant Bayes posterior marginal)
                get<micro_P_state_t10_y1>(out())() = std::move(start_micro_P_state_post);
                get<P_mean_t10_y1>(out())()        = std::move(Maybe_start_macro.value().first);
                get<P_Cov_t10_y1>(out())()         = std::move(Maybe_start_macro.value().second);

                // _t20_y1 (end-state Bayes+integrate posterior marginal)
                get<micro_P_state_t20_y1>(out())() = r_micro_P_state();
                get<P_mean_t20_y1>(out())()        = r_P_mean();
                get<P_Cov_t20_y1>(out())()         = r_P_cov();

                // _0t_y0 / _0t_y1 — joint distributions and macro cross-second-moment projections.
                get<micro_P_state_0t_y0>(out())() = joint_prior;
                get<micro_P_state_0t_y1>(out())() = joint_posterior;
                get<P_mean_0t_y0>(out())()        = project_joint_to_macro(joint_prior);
                get<P_mean_0t_y1>(out())()        = project_joint_to_macro(joint_posterior);
                // Cross-cov in second-moment-form (mirrors macro storage = P_Cov_stored · t_P).
                get<P_cross_cov_0t_y0>(out())()   = project_joint_to_macro(joint_prior);
                get<P_cross_cov_0t_y1>(out())()   = project_joint_to_macro(joint_posterior);

                // Carry through Qdt-derived quantities for downstream consumers (mirrors
                // qmodel.h:4039-4040).
                get<micro_gtotal_ij>(out()) = get<micro_gtotal_ij>(t_micro_Qdt);
                get<micro_gmean_ij>(out())  = get<micro_gmean_ij>(t_micro_Qdt);

            } else if constexpr (averaging::value == 1) {
                // === MR mode ===
                // _t2_y0: full-step Markov-stepped predictive
                // _t2_y1: full-step Bayes+Markov posterior  (= r_micro_P_state, already projected)
                // _t1_y1: instant Bayes-only posterior (no Markov) = bayes_post
                // _0t_y0/y1: joint (start, end) distributions in micro / macro space.
                auto inst_micro_P_state_post = bayes_post;   // copy — used downstream by joint-builders
                auto Maybe_inst_macro = project_to_macro(inst_micro_P_state_post);
                if (!Maybe_inst_macro)
                    return error_message("dynamic:instant-Bayes macro-projection | " +
                                          Maybe_inst_macro.error()());

                // _t2_y0
                get<micro_P_state_t2_y0>(out())() = pred_micro_P_state;
                get<P_mean_t2_y0>(out())()        = pred_P_mean;
                get<P_Cov_t2_y0>(out())()         = pred_P_cov;

                // _t2_y1
                get<micro_P_state_t2_y1>(out())() = r_micro_P_state();
                get<P_mean_t2_y1>(out())()        = r_P_mean();
                get<P_Cov_t2_y1>(out())()         = r_P_cov();

                // _t1_y1 (instant Bayes-only at start)
                get<micro_P_state_t1_y1>(out())() = inst_micro_P_state_post;
                get<P_mean_t1_y1>(out())()        = std::move(Maybe_inst_macro.value().first);
                get<P_Cov_t1_y1>(out())()         = std::move(Maybe_inst_macro.value().second);

                // _0t_y0 / _0t_y1 — joint = diag(start_dist) · t_micro_P (per-row scaling).
                Matrix<double> joint_prior = diag(p_micro_P_state())     * t_micro_P();
                Matrix<double> joint_posterior = diag(inst_micro_P_state_post) * t_micro_P();
                get<micro_P_state_0t_y0>(out())() = joint_prior;
                get<micro_P_state_0t_y1>(out())() = joint_posterior;
                get<P_mean_0t_y0>(out())()        = project_joint_to_macro(joint_prior);
                get<P_mean_0t_y1>(out())()        = project_joint_to_macro(joint_posterior);
                get<P_cross_cov_0t_y0>(out())()   = project_joint_to_macro(joint_prior);
                get<P_cross_cov_0t_y1>(out())()   = project_joint_to_macro(joint_posterior);

            } else {  // averaging::value == 0
                // === R mode: single instant at midpoint ===
                // _t15_y0: midpoint prior (= p_micro_P_state, already half-advanced at line 893)
                // _t15_y1: midpoint Bayes-only posterior = bayes_post
                // _t2_y0:  full-step Markov-stepped predictive
                // _t2_y1:  full-step Bayes+Markov posterior (= r_micro_P_state)
                auto Maybe_mid_pri_macro = project_to_macro(p_micro_P_state());
                if (!Maybe_mid_pri_macro) {
                    return error_message("dynamic:midpoint-prior macro-projection | " +
                                          Maybe_mid_pri_macro.error()());
                }
                auto mid_post_micro = bayes_post;            // copy — used by storage and projection below
                auto Maybe_mid_post_macro = project_to_macro(mid_post_micro);
                if (!Maybe_mid_post_macro) {
                    return error_message("dynamic:midpoint-posterior macro-projection | " +
                                          Maybe_mid_post_macro.error()());
                }

                // _t15_y0 (midpoint prior)
                get<micro_P_state_t15_y0>(out())() = p_micro_P_state();
                get<P_mean_t15_y0>(out())()        = std::move(Maybe_mid_pri_macro.value().first);
                get<P_Cov_t15_y0>(out())()         = std::move(Maybe_mid_pri_macro.value().second);

                // _t15_y1 (midpoint Bayes-only)
                get<micro_P_state_t15_y1>(out())() = std::move(mid_post_micro);
                get<P_mean_t15_y1>(out())()        = std::move(Maybe_mid_post_macro.value().first);
                get<P_Cov_t15_y1>(out())()         = std::move(Maybe_mid_post_macro.value().second);

                // _t2_y0 (full-step predictive)
                get<micro_P_state_t2_y0>(out())() = pred_micro_P_state;
                get<P_mean_t2_y0>(out())()        = pred_P_mean;
                get<P_Cov_t2_y0>(out())()         = pred_P_cov;

                // _t2_y1 (full-step post-observation)
                get<micro_P_state_t2_y1>(out())() = r_micro_P_state();
                get<P_mean_t2_y1>(out())()        = r_P_mean();
                get<P_Cov_t2_y1>(out())()         = r_P_cov();
            }

            return out;
        }
    }


    template <bool dynamic, class recursive, class averaging, class variance, class C_micro_Patch_State,
              class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> && (U<C_micro_Patch_State, micro_Patch_State>))
    auto safely_calculate_micro_Algo_State(C_micro_Patch_State const& t_micro_prior, C_Qdt const& t_Qdt,
                                     C_Patch_Model const& m, C_double const& N,
                                     const Patch_current& p_y, double fs) const {
        if constexpr (!recursive::value) {
            return safely_calculate_micro_Algo_State_non_recursive<dynamic, averaging>(
                t_micro_prior, t_Qdt, m, p_y, fs, N);
        }
        return safely_calculate_micro_Algo_State_recursive<dynamic, averaging>(
            t_micro_prior, t_Qdt, m, N, p_y, fs);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<Macro_State<vVars...>> update_macro_state(Macro_State<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, logL const& t_logL, elogL const& t_elogL,... ) const 
                                              {
        // Update patch state for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            } else if constexpr (requires { algo.get_P_mean(); algo.get_P_Cov(); }) {
                get<P_mean>(ps())() = algo.get_P_mean();
                get<P_Cov>(ps())() = algo.get_P_Cov();
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        // Optional accumulators when both state and Algo_State provide them.
        if constexpr (has_var_c<Macro_State<vVars...>&, elogL> ) {
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        }
        if constexpr (has_var_c<Macro_State<vVars...>&, vlogL> ) {
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;
        }

        if constexpr (has_var_c<Macro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            // Avoid uninitialized Constant defaults.
            if constexpr (has_var_c<Element&, elogL>) {
                get<elogL>(el)() = t_elogL();   
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                get<vlogL>(el)() = 0.5;
            }

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            // Common predicted fields.
            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            copy_component(std::type_identity<Chi2>{});

            if constexpr (has_var_c<Element&, Algo_State_Dynamic>) {
                if constexpr (requires { get<Algo_State_Dynamic>(el) = algo; }) {
                    get<Algo_State_Dynamic>(el) = algo;
                }
            }

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    
    // MMicro_State overload — mirrors the Vector_Space one but selected explicitly
    // because MMicro_State's derived-to-base relation isn't considered by template
    // argument deduction.  Body intentionally minimal: accumulate logL/elogL into
    // the running state and (if the type carries Evolution) push the per-step element.
    template <class... vVars, class C_Algo_State, class C_logL, class C_elogL>
        requires((U<C_Algo_State, micro_Algo_State> || U<C_Algo_State, micro_Algo_State_Dynamic>) ||
                 U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<MMicro_State<vVars...>> update_macro_state(
        MMicro_State<vVars...>&& t_prior_all, C_Algo_State&& algo, C_logL const& t_logL,
        C_elogL const& t_elogL, ...) const {
        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();
        if constexpr (var::has_it_v<MMicro_State<vVars...>, elogL>)
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL();
        if constexpr (var::has_it_v<MMicro_State<vVars...>, vlogL>)
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;
        if constexpr (var::has_it_v<MMicro_State<vVars...>, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all);
            evo.emplace_back(algo);
        }
        // Carry the post-step microstate distribution from algo into the running state.
        if constexpr (has_var_c<decltype(algo()) const&, micro_P_state>) {
            get<micro_Patch_State>(t_prior_all)() = get<micro_P_state>(algo());
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_logL, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<Vector_Space<vVars...>> update_macro_state(Vector_Space<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, C_logL const& t_logL,C_elogL const& t_elogL,...) const {
        // If this Vector_Space carries a Patch_State, keep recursion semantics consistent with
        // Macro_State/dMacro_State/ddMacro_State updates.
        if constexpr (has_var_c<Vector_Space<vVars...>&, Patch_State>) {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            } else if constexpr (requires { algo.get_P_mean(); algo.get_P_Cov(); }) {
                get<P_mean>(ps())() = algo.get_P_mean();
                get<P_Cov>(ps())() = algo.get_P_Cov();
            }
        }

        get<logL>(t_prior_all) ()= get<logL>(t_prior_all)() + t_logL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, elogL>)
            get<elogL>(t_prior_all) ()= get<elogL>(t_prior_all)() + t_elogL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, vlogL>)
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (var::has_it_v<Vector_Space<vVars...>, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all);

            evo.emplace_back(algo);
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_elogL >
    Maybe_error<var::Derivative<Vector_Space<vVars...>, var::Parameters_transformed>>
    update(var::Derivative<Vector_Space<vVars...>, var::Parameters_transformed>&& t_prior_all,
           C_Algo_State&& algo,
           var::Derivative<logL, var::Parameters_transformed> const& t_logL, C_elogL const& t_elogL,...) const {
        get<logL>(t_prior_all) ()= get<logL>(t_prior_all)() + t_logL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, elogL> )
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        if constexpr (var::has_it_v<Macro_State<vVars...>, vlogL> )
            get<vlogL>(t_prior_all) ()= get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (var::has_it_v<Vector_Space<vVars...>, Evolution> &&
                      var::has_it_v<std::decay_t<C_Algo_State>, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all);

            evo.emplace_back(algo);
        }
        return std::move(t_prior_all);
    }


    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<dMacro_State<vVars...>> update_macro_state(
        dMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL ,
        C_elogL const& t_elogL, 
        var::Derivative<y_mean, var::Parameters_transformed> const& t_ymean, 
        var::Derivative<y_var, var::Parameters_transformed> const& t_yvar) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }
        auto d_y_mean= t_ymean.derivative()();
        auto d_y_var= t_yvar.derivative()();

        auto r_y_var=t_yvar.primitive()();
        
        
        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();
       if constexpr (var::has_it_v<dMacro_State<vVars...>, covariance<Grad>> ){
        auto t_CovGradient = covariance<Grad>(
            parameter_spd_payload(XXT(t_logL.derivative()()), var::get_dx_of_dfdx(t_logL)));
        get<covariance<Grad>>(t_prior_all)() = get<covariance<Grad>>(t_prior_all)() + t_CovGradient();
       }

       if constexpr (var::has_it_v<dMacro_State<vVars...>, Hessian> ){
        auto t_Hessian = parameter_spd_payload(
            XXT(d_y_mean) / r_y_var + XXT(d_y_var) / (2 * r_y_var * r_y_var),
            var::get_dx_of_dfdx(t_logL));
         get<Hessian>(t_prior_all)() = get<Hessian>(t_prior_all)() + t_Hessian;
       }
       if constexpr (var::has_it_v<dMacro_State<vVars...>, elogL> )
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        if constexpr (var::has_it_v<dMacro_State<vVars...>, vlogL> )
            get<vlogL>(t_prior_all) ()= get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (has_var_c<dMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            if constexpr (has_var_c<Element&, elogL>) {
                get<elogL>(el) = t_elogL;
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                if constexpr (requires { get<vlogL>(el)() = 0.5; }) {
                    get<vlogL>(el)() = 0.5;
                } else {
                    get<vlogL>(el) = 0.5;
                }
            }



            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<r_std>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    // Derivative-aware micro analog. Mirrors the dMacro_State version above
    // but carries micro_Patch_State (just micro_P_state). The Evolution element
    // is the unified one (micro_gradient_*_element), so both r_std (from the
    // moment-matched residual) and micro_r_std (from Rosenblatt PIT) get
    // populated when the algo state holds them.
    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, micro_Algo_State> || U<C_Algo_State, micro_Algo_State_Dynamic>)
    Maybe_error<dMicro_State<vVars...>> update_macro_state(
        dMicro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,
        C_elogL const& t_elogL,
        var::Derivative<y_mean, var::Parameters_transformed> const& t_ymean,
        var::Derivative<y_var, var::Parameters_transformed> const& t_yvar) const {
        // Update micro patch state (the propagating microstate distribution) for recursion.
        {
            auto& ps = get<var::Derivative<micro_Patch_State, var::Parameters_transformed>>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, micro_P_state>) {
                get<micro_P_state>(ps()) = get<micro_P_state>(algo());
            }
        }
        auto d_y_mean = t_ymean.derivative()();
        auto d_y_var = t_yvar.derivative()();
        auto r_y_var = t_yvar.primitive()();

        get<var::Derivative<logL, var::Parameters_transformed>>(t_prior_all)() =
            get<var::Derivative<logL, var::Parameters_transformed>>(t_prior_all)() + t_logL();

        if constexpr (var::has_it_v<dMicro_State<vVars...>, covariance<Grad>>) {
            auto t_CovGradient = covariance<Grad>(parameter_spd_payload(
                XXT(t_logL.derivative()()), var::get_dx_of_dfdx(t_logL)));
            get<covariance<Grad>>(t_prior_all)() =
                get<covariance<Grad>>(t_prior_all)() + t_CovGradient();
        }
        if constexpr (var::has_it_v<dMicro_State<vVars...>, Hessian>) {
            auto t_Hessian = parameter_spd_payload(
                XXT(d_y_mean) / r_y_var + XXT(d_y_var) / (2 * r_y_var * r_y_var),
                var::get_dx_of_dfdx(t_logL));
            get<Hessian>(t_prior_all)() = get<Hessian>(t_prior_all)() + t_Hessian;
        }
        if constexpr (var::has_it_v<dMicro_State<vVars...>, elogL>)
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL();
        if constexpr (var::has_it_v<dMicro_State<vVars...>, vlogL>)
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (has_var_c<dMicro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) get<logL>(el) = t_logL;
            if constexpr (has_var_c<Element&, elogL>) get<elogL>(el) = t_elogL;
            if constexpr (has_var_c<Element&, vlogL>) {
                if constexpr (requires { get<vlogL>(el)() = 0.5; }) {
                    get<vlogL>(el)() = 0.5;
                } else {
                    get<vlogL>(el) = 0.5;
                }
            }

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<r_std>{});
            copy_component(std::type_identity<micro_r_std>{});

            copy_component(std::type_identity<micro_P_state>{});
            copy_component(std::type_identity<Chi2>{});
            copy_component(std::type_identity<micro_Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update_macro_state(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,...) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        if constexpr (has_var_c<ddMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            auto seed_zero = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id>) {
                    using Comp = std::decay_t<decltype(get<Id>(el))>;
                    if constexpr (var::is_derivative_v<Comp>) {
                        if constexpr (std::constructible_from<Comp, decltype(t_logL.dx()) const&>) {
                            get<Id>(el) = Comp(t_logL.dx());
                        }
                    } else if constexpr (requires { get<Id>(el)() = 0.0; }) {
                        get<Id>(el)() = 0.0;
                    }
                }
            };

            seed_zero(std::type_identity<elogL>{});
            seed_zero(std::type_identity<vlogL>{});

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update_macro_state(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,
        C_elogL const& t_elogL,...) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        if constexpr (var::has_it_v<ddMacro_State<vVars...>, elogL>) {
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL();
        }
        if constexpr (var::has_it_v<ddMacro_State<vVars...>, vlogL>) {
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;
        }

        if constexpr (has_var_c<ddMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            if constexpr (has_var_c<Element&, elogL>) {
                if constexpr (requires { get<elogL>(el) = t_elogL; }) {
                    get<elogL>(el) = t_elogL;
                } else if constexpr (requires { get<elogL>(el) = var::primitive(t_elogL); }) {
                    get<elogL>(el) = var::primitive(t_elogL);
                } else if constexpr (requires { get<elogL>(el)() = var::primitive(t_elogL)(); }) {
                    get<elogL>(el)() = var::primitive(t_elogL)();
                }
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                if constexpr (requires { get<vlogL>(el)() = 0.5; }) {
                    get<vlogL>(el)() = 0.5;
                }
            }

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class recursive, class averaging, class variance, class variance_correction,
              class FunctionTable, class C_Macro_State, class C_micro_Qdt, class C_Patch_Model,
              class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 /*(U<std::decay_t<C_Patch_State>,
                                                         Patch_State>||U<std::decay_t<C_Patch_State>,
                                                         Patch_State_and_Evolution>
               )&& U<C_Patch_Model, Patch_Model>
                                                   &&*/
                 U<C_double, double> && (U<C_micro_Qdt, micro_Qdt> || U<C_micro_Qdt, micro_Qdtm> || U<C_micro_Qdt, micro_Qdtg>))

        Maybe_error<C_Macro_State> Micror(FunctionTable&, C_Macro_State&& t_micro_prior_all,
                                        C_micro_Qdt const& t_micro_Qdt, C_Patch_Model const& m,
                                        C_double const& Nch, const Patch_current& p_y,
                                        double fs) const {
        using Transf = transformation_type_t<C_micro_Qdt>;

        auto& t_prior = get<micro_Patch_State>(t_micro_prior_all);

        auto Maybe_Algo =
            safely_calculate_micro_Algo_State<is_Algo_dynamic<C_Macro_State>(), recursive, averaging,
                                        variance>(t_prior, t_micro_Qdt, m, Nch, p_y, fs);
        if (!Maybe_Algo)
            return Maybe_Algo.error();

        auto r_Algo_state = std::move(Maybe_Algo.value());
        auto r_y_mean = get<y_mean>(r_Algo_state());
        auto r_y_var = get<y_var>(r_Algo_state());
        auto y = p_y.value();
        bool y_is_nan = std::isnan(y);

        // Per-step logL / elogL — both in mixture form, computed and stored inside
        // safely_calculate_micro_Algo_State_*.  Just retrieve here.  The diagnostic
        // moment-match form is stored alongside in gaussian_logL.
        auto r_logL = get<logL>(r_Algo_state());
        auto r_elogL = get<elogL>(r_Algo_state());

        auto r_prior_all =
            update_macro_state(std::move(t_micro_prior_all), std::move(r_Algo_state),
                               std::move(r_logL), std::move(r_elogL), r_y_mean, r_y_var);
        return {std::move(r_prior_all)};
    }

    // template <class recursive, class averaging, class variance, class variance_correction,
    //           class FunctionTable, class C_Patch_State, class C_Qdt, class C_Patch_Model,
    //           class C_double>

    //     requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
    //              uses_variance_aproximation_c<variance> &&
    //              uses_taylor_variance_correction_aproximation_c<variance_correction> &&
    //              /*(U<std::decay_t<C_Patch_State>,
    //                                                          Patch_State>||U<std::decay_t<C_Patch_State>,
    //                                                          Patch_State_and_Evolution>
    //                )&& U<C_Patch_Model, Patch_Model>
    //                                                    &&*/
    //              U<C_double, double> && (U<C_Qdt, Qdt> || U<C_Qdt, Qdtm> || U<C_Qdt, Qdtg>))

    // Maybe_error<C_Patch_State> Macror_old(FunctionTable&, C_Patch_State&& t_prior, C_Qdt const& t_Qdt,
    //                                       C_Patch_Model const& m, C_double const& Nch,
    //                                       const Patch_current& p_y, double fs) const {
    //     get<macror_algorithm>(t_prior)() =
    //         ToString(MacroR2<recursive, averaging, variance, variance_correction>{});
    //     using Transf = transformation_type_t<C_Qdt>;

    //     auto& p_P_mean = get<P_mean>(t_prior);
    //     auto SmD = get<P_Cov>(t_prior)() - diag(p_P_mean());
    //     auto& y = p_y.value();
    //     auto& t_tolerance = get<Probability_error_tolerance>(m);
    //     auto& t_min_P = get<min_P>(m);
    //     auto y_baseline = get<Current_Baseline>(m);
    //     auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
    //              get<Pink_Noise>(m).value();

    //     auto N = Nch;
    //     Matrix<double> u(p_P_mean().size(), 1, 1.0);

    //     auto N_states = p_P_mean().ncols();

    //     Op_t<Transf, double> ms = 0;
    //     if constexpr (variance::value)
    //         ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

    //     auto& t_gmean_i = get<gmean_i>(t_Qdt);
    //     auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    //     auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
    //     auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
    //                getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    //     Op_t<Transf, y_mean> r_y_mean;
    //     Op_t<Transf, y_var> r_y_var;

    //     Op_t<Transf, double> sSg;
    //     auto t_P = get<P>(t_Qdt);

    //     r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());

    //     if (std::isnan(y)) {
    //         get<macror_algorithm>(t_prior)() =
    //             ToString(MacroR2<uses_recursive_aproximation<false>, averaging, variance,
    //                              variance_correction>{});

    //         auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();

    //         auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

    //         auto Maybe_P_cov = to_Covariance_Probability(r_P_cov() + diag(r_P_mean()));
    //         if (!Maybe_P_cov)
    //             return Maybe_P_cov.error();
    //         r_P_cov() = std::move(Maybe_P_cov.value());
    //         if constexpr (U<C_Patch_State, Patch_State_and_Evolution>) {
    //             auto& ev = get<Macro_State_Evolution>(t_prior);
    //             ev().push_back(build<Patch_State>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), r_P_mean, r_P_cov, r_y_mean, r_y_var,
    //                 plogL(NaN), eplogL(NaN), vplogL(NaN), get<macror_algorithm>(t_prior)));
    //             return build<Patch_State_and_Evolution>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), std::move(ev));
    //         } else if constexpr (U<C_Patch_State, Patch_State_and_y_Evolution>) {
    //             auto& yev = get<ymean_Evolution>(t_prior);
    //             yev().push_back(r_y_mean);
    //             auto& yvev = get<yvar_Evolution>(t_prior);
    //             yvev().push_back(r_y_var);

    //             return build<Patch_State_and_y_Evolution>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), std::move(yev), std::move(yvev));
    //         } else if constexpr (var::is_derivative_v<C_Patch_State>) {
    //             return build<Patch_State_and_Hessian>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), get<FIM>(t_prior));

    //         } else
    //             return Patch_State(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior));
    //     }

    //     constexpr bool PoissonDif = true;
    //     using std::abs;
    //     if constexpr (PoissonDif)
    //         e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
    //     else
    //         e = e + get<Proportional_Noise>(m).value() * abs(y);

    //     auto r_y_mean_max = max_possible_value_of_ymean(N_Ch_mean_value(primitive(Nch)),
    //                                                     primitive(get<g>(m)), primitive(y_baseline));

    //     auto r_y_mean_min = min_possible_value_of_ymean(N_Ch_mean_value(primitive(Nch)),
    //                                                     primitive(get<g>(m)), primitive(y_baseline));

    //     // if ((primitive(r_y_mean()) - r_y_mean_max()) >
    //     //     std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_max())) * 1e-3)
    //     //     std::cerr << "\n max violation" << r_y_mean() << "  vs  max: " << r_y_mean_max();
    //     // if ((r_y_mean_min() - primitive(r_y_mean())) >
    //     //     std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_min())) *
    //     //         1e-1)
    //     // std::cerr << "\n min violation\n"
    //     //          << r_y_mean() << "  vs  min: " << r_y_mean_min();

    //     if (std::isfinite(primitive(gSg)) && primitive(gSg) > 0) {
    //         if (isfinite(primitive(ms)) && primitive(ms) > 0) {
    //             r_y_var = build<y_var>(e + N * gSg + N * ms);
    //         } else {
    //             r_y_var = build<y_var>(e + N * gSg);
    //         }
    //     } else {
    //         if (isfinite(primitive(ms)) && primitive(ms) > 0)
    //             r_y_var = build<y_var>(e + N * ms);
    //         else
    //             r_y_var = build<y_var>(e);
    //     }

    //     auto dy = y - r_y_mean();
    //     auto chi = dy / r_y_var();
    //     Op_t<Transf, P_mean> r_P_mean;
    //     Op_t<Transf, P_Cov> r_P_cov;

    //     if constexpr (!recursive::value) {
    //         r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));

    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean.valid())
    //             return Maybe_r_P_mean.error();

    //         r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

    //         auto Maybe_r_P_cov = to_Covariance_Probability(r_P_cov() + diag(r_P_mean()));

    //         if (!Maybe_r_P_cov.valid())
    //             return Maybe_r_P_cov.error();
    //         r_P_cov() = std::move(Maybe_r_P_cov.value());

    //     } else if constexpr (!variance_correction::value) {
    //         auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();

    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS);
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();
    //         r_P_mean() = std::move(Maybe_r_P_mean.value());

    //         auto Maybe_r_P_cov = to_Covariance_Probability(
    //             AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) - (N / r_y_var()) * XTX(gS));
    //         if (!Maybe_r_P_cov)
    //             return Maybe_r_P_cov.error();

    //         r_P_cov() = std::move(Maybe_r_P_cov.value());
    //     } else {
    //         auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
    //         auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
    //         auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    //         auto& t_gvar_i = get<gvar_i>(t_Qdt);
    //         auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    //         auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
    //         auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));

    //         auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    //         auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    //         auto e_mu = e + N * ms0;
    //         r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) - N * 0.5 / e_mu * sSg + y_baseline();
    //         auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    //         r_y_var() = var::max(e, e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
    //         auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
    //         auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();
    //         r_P_mean() =
    //             to_Probability(p_P_mean() * t_P() + chi * gS - (chi * zeta * sSg + 0.5 / e_mu) * sS);

    //         r_P_cov() = AT_B_A(t_P(), SmD) + diag(r_P_mean() * t_P()) -
    //                     (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
    //                     (2.0 * N / r_y_var() * zeta * sSg) * X_plus_XT(TranspMult(sS, gS)) -
    //                     (N / r_y_var()) * XTX(gS);
    //     }

    //     if (!all_Probability_elements(primitive(r_P_mean())) ||
    //         !all_Covariance_elements(primitive(r_P_cov()))) {
    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();

    //         r_P_mean() = Maybe_r_P_mean.value();
    //         auto Maybe_r_P_cov =
    //             to_Covariance_Probability(AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));

    //         if (!Maybe_r_P_cov)
    //             return Maybe_r_P_cov.error();

    //         r_P_cov() = std::move(Maybe_r_P_cov.value());
    //         get<macror_algorithm>(t_prior)() =
    //             ToString(MacroR2<uses_recursive_aproximation<false>, averaging, variance,
    //                              variance_correction>{});
    //     }

    //     auto chi2 = dy * chi;

    //     Op_t<Transf, plogL> r_plogL;
    //     Op_t<Transf, eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5);
    //     if (primitive(r_y_var()) > 0.0) {
    //         if (get<Proportional_Noise>(m).value() == 0) {
    //             r_plogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2;
    //             if constexpr (var::is_derivative_v<std::decay_t<decltype(r_plogL)>>)
    //                 if (std::isnan(var::derivative(r_plogL())()[0]))
    //                     std::cerr << "nan der\n";
    //             r_eplogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5;
    //         } else {
    //             r_plogL() = -log(var::Poisson_noise_normalization(
    //                             primitive(r_y_var()), primitive(get<Proportional_Noise>(m).value()))) -
    //                         0.5 * chi2;
    //             r_eplogL() = var::Poisson_noise_expected_logL(
    //                 primitive(r_y_var()), primitive(get<Proportional_Noise>(m).value()));
    //         }
    //     } else {
    //         std::stringstream ss;
    //         ss << "Negative variance!!\n";
    //         ss << "\nr_y_var=\t" << r_y_var;
    //         ss << "\ngSg=\t" << gSg;
    //         return error_message(ss.str());
    //     }

    //     vplogL r_vlogL(0.5);
    //     if (std::isnan(primitive(r_plogL()))) {
    //         std::stringstream ss;
    //         ss << "likelihood is nan \n patch current=";
    //         print(ss, p_y) << "Qdt";
    //         print(ss, primitive(t_Qdt)) << "tprior";
    //         print(ss, primitive(t_prior));
    //         return error_message(ss.str());
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_Evolution>) {
    //         auto& ev = get<Macro_State_Evolution>(t_prior);
    //         ev().push_back(build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
    //                                           build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //                                           build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), r_P_mean,
    //                                           r_P_cov, r_y_mean, r_y_var, r_plogL, r_eplogL, r_vlogL,
    //                                           get<macror_algorithm>(t_prior)));
    //         return build<Patch_State_and_Evolution>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), std::move(ev));
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_y_Evolution>) {
    //         auto& yev = get<ymean_Evolution>(t_prior);
    //         yev().push_back(r_y_mean);
    //         auto& yvev = get<yvar_Evolution>(t_prior);
    //         yvev().push_back(r_y_var);
    //         return build<Patch_State_and_y_Evolution>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), std::move(yev), std::move(yvev));
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_Hessian>) {
    //         auto r_J = derivative(r_y_mean)();
    //         auto r_JS = derivative(r_y_var)();

    //         auto r_FIM =
    //             XXT(r_J) / primitive(r_y_var()) + XXT(r_JS) / (2.0 * sqr(primitive(r_y_var())));

    //         return build<Patch_State_and_Hessian>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), FIM(get<FIM>(t_prior)() + r_FIM));
    //     } else {
    //         return build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
    //                                   build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //                                   build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
    //                                   std::move(r_P_mean), std::move(r_P_cov), std::move(r_y_mean),
    //                                   std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //                                   get<macror_algorithm>(t_prior));
    //     }
    // }

    // === Initial micro_Patch_State ===
    // Builds the initial microstate distribution from the model's P_initial,
    // expanded to the N-channel multinomial via P_mean_to_micro_P_state.
    // Replaces the macro `init` (which returns Patch_State with P_mean / P_Cov);
    // micro carries only the microstate distribution in micro_Patch_State.
    template <class C_Patch_Model>
    auto init_micro(const C_Patch_Model& m, std::size_t Nchannels)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, micro_Patch_State>> {
        auto t_P_mean = build<P_mean>(get<P_initial>(m)());
        auto r_micro_P_state = P_mean_to_micro_P_state(t_P_mean, Nchannels);
        Transfer_Op_to<C_Patch_Model, micro_Patch_State> out;
        get<micro_P_state>(out())() = std::move(r_micro_P_state());
        return out;
    }

    // === Per-recording micro log-likelihood ===
    // Mirrors Macro_DMR::log_Likelihood at qmodel.h:4919 in shape, but:
    //   - init -> init_micro(m, Nchannels)        (microstate seed)
    //   - calc_Qdt/Qdtm/Qdtg -> calc_micro_*      (per-step Qdt with per-state moments)
    //   - MacroR2<>{}(...) -> Micror<>(...)       (member fn instead of struct callable)
    //   - the `adaptive` template param is dropped (micro uses the exact mixture, no
    //     binomial/Gaussian-approx switching)
    //   - macro-specific test_Qx / test_eigen / test_macroir_derivative blocks are
    //     stripped (they referenced macro Kalman shapes; if you want derivative
    //     diagnostics here, Macro_DMR{}.calc_Qx(...) is still callable from inside).
    template <class recursive, class averaging, class variance, class variance_correction,
              class MicroState, class FuncTable, class C_Parameters, class Model>
        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 is_of_this_template_type_v<FuncTable, FuncMap_St>)
    auto log_Likelihood(FuncTable& f, const Model& model, const C_Parameters& par,
                        const Recording& y, const Experiment& e) -> Maybe_error<MicroState> {
        auto Maybe_m = model(par);
        if (!is_valid(Maybe_m)) {
            return get_error(Maybe_m);
        }
        auto m = std::move(get_value(Maybe_m));
        auto fs = get<Frequency_of_Sampling>(e).value();

        // Microstate enumeration size is fixed at the start (combinatorial in N_channels).
        // Per-step Nch is still interpolated per segment for the algorithm; only the
        // *enumeration* uses a representative integer Nchannels.
        auto Nchs = get<N_Ch_mean>(m)();
        std::size_t Nchannels =
        
            static_cast<std::size_t>(std::round(primitive(Nchs)[0]));

        auto ini = init_micro(m, Nchannels);
        if (!ini) return ini.error();
        auto f_local = f.create("_lik");

        using DX = var::dx_of_dfdx_t<C_Parameters>;
        if constexpr (!std::is_same_v<DX, var::NoDerivative>) {
            MACRODR_DX_ASSERT(var::has_dx(par) &&
                              "log_Likelihood: derivative parameters missing dx()");
            MACRODR_DX_ASSERT(var::has_dx(ini.value()) &&
                              "log_Likelihood: init_micro result missing dx in derivative mode");
        }

        auto t_micro = MicroState(std::move(ini.value()));

        if constexpr (!std::is_same_v<DX, var::NoDerivative>) {
            MACRODR_DX_ASSERT(var::has_dx(t_micro) &&
                              "log_Likelihood: MicroState constructed without dx in derivative mode");
        }

        auto Maybe_run = fold(
            0ul, y().size(), std::move(t_micro),
            [this, &f_local, &m, fs, &e, &y, Nchannels](MicroState&& t_prior, std::size_t i_step) {
                Agonist_evolution const& t_step =
                    get<Agonist_evolution>(get<Recording_conditions>(e)()[i_step]);

                auto time = get<Time>(get<Recording_conditions>(e)()[i_step])();
                auto time_segment = get<N_Ch_mean_time_segment_duration>(m)();
                auto Nchs = get<N_Ch_mean>(m)();
                std::size_t i_segment =
                    std::min(Nchs.size() - 1.0, std::floor(time / time_segment));
                auto j_segment = std::min(Nchs.size() - 1, i_segment + 1);
                auto r = std::max(1.0, time / time_segment - i_segment);
                auto Nch = Nchs[i_segment] * (1 - r) + r * Nchs[j_segment];

                // Compute the per-step micro Qdt according to averaging mode:
                //   avg=0 (R)  -> micro_Qdtg : half-step micro_P_half + per-microstate g
                //   avg=1 (MR) -> micro_Qdtm : full-step P + per-microstate gmean_i / gvar_i
                //   avg=2 (IR) -> micro_Qdt  : full-step P + per-pair gmean_ij / gvar_ij
                auto Maybe_t_Qdt =
                    [this, &f_local, &m, &t_step, fs, Nchannels]() {
                        if constexpr (averaging::value == 0)
                            return calc_micro_Qdtg(f_local, m, t_step, fs, Nchannels);
                        else if constexpr (averaging::value == 1)
                            return calc_micro_Qdtm(f_local, m, t_step, fs, Nchannels);
                        else
                            return calc_micro_Qdt(f_local, m, t_step, fs, Nchannels);
                    }();
                if (!Maybe_t_Qdt)
                    return Maybe_error<MicroState>(error_message(
                        "k=" + std::to_string(i_step) + " | calc_micro_Qdt | " +
                        Maybe_t_Qdt.error()()));
                auto t_Qdt = std::move(Maybe_t_Qdt.value());

                auto Maybe_micror = Micror<recursive, averaging, variance, variance_correction>(
                    f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                if (!Maybe_micror)
                    return Maybe_error<MicroState>(error_message(
                        "k=" + std::to_string(i_step) + " | Micror | " +
                        Maybe_micror.error()()));
                return Maybe_error<MicroState>(std::move(Maybe_micror.value()));
            });
        f += f_local;
        if (!Maybe_run)
            return Maybe_run.error();
        return std::move(Maybe_run.value());
    }


};

// Definitions for the function-table lambda wrappers forward-declared in
// CLI_function_table.h. We can't include micro_monoid.h there (circular with
// qmodel.h), so the wrappers live here, after Micro_DMR is fully defined, and
// are instantiated only when the function-table memoizer lambdas actually fire.
template <class FT, class M>
auto micro_Qdtg_step_call(FT& f, const M& m, const Agonist_step& t_step, double fs,
                          std::size_t Nchannels) {
    return Micro_DMR{}.calc_micro_Qdtg_agonist_step(f, m, t_step, fs, Nchannels);
}
template <class FT, class M>
auto micro_Qdtm_step_call(FT& f, const M& m, const Agonist_step& t_step, double fs,
                          std::size_t Nchannels) {
    return Micro_DMR{}.calc_micro_Qdtm_agonist_step(f, m, t_step, fs, Nchannels);
}
template <class FT, class M>
auto micro_Qdt_step_call(FT& f, const M& m, const Agonist_step& t_step, double fs,
                         std::size_t Nchannels) {
    return Micro_DMR{}.calc_micro_Qdt_agonist_step(f, m, t_step, fs, Nchannels);
}



 













}  // namespace macrodr
