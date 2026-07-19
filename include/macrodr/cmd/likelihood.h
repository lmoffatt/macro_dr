#pragma once

#include <derivative_fwd.h>
#include <distributions.h>
#include <macrodr/interface/IModel.h>
#include <macrodr/cmd/detail/write_csv_common.h>
#include <gauss_newton.h>
#include <probit_samples.h>

#include <type_traits>
#include <utility>

#include "patch_model.h"
#include "qmodel.h"
// micro_types.h defines the unified per-step Evolution element shape (a
// superset of the macro one) and re-exports dMacro_State_Ev_gradient_* with
// that wider element list. Must come AFTER qmodel.h because micro_types.h
// transitively pulls in micro_full.h, which references Macro_DMR (defined in
// qmodel.h).
#include <micro_types.h>

namespace macrodr::cmd {

// Core likelihood-algorithm builder. The compute family is an int selector —
// family_macro / family_micro / family_nonlinearsqr — merging FOUR
// Likelihood_Model_regular branches (macro, macro-taylor, micro, and the lean
// nonlinear-least-squares fold). build_likelihood_function below is a thin bool
// wrapper over this for existing scripts and call sites.
inline auto build_likelihood_function_with_family(
    const ModelPtr& model0, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation, int family_approximation,
    bool taylor_qdt_approximation = false) {
    // Map deprecated bool taylor_qdt_approximation to the new int qdt_method:
    //   false → 0 (eig), true → 2 (schur). Phase 9 wired the previous "true"
    //   semantics through calc_Qdt_schur, so this preserves user behavior.
    int qdt_method_int = taylor_qdt_approximation ? 2 : 0;
    auto nsub = Simulation_n_sub_dt(100);
    const interface::IModel<var::Parameters_values>& model_ref = *model0;

    return merge_Maybe_variant(
        merge_Maybe_variant(
            merge_Maybe_variant(
                // Macro non-taylor-variance-correction branch: taylor_qdt fixed
                // to false (eigen path is well-conditioned for k×k macro Q).
                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              false>,
                    var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                    decltype(model_ref),
                    var::constexpr_Var_domain<int, uses_qdt_method, 0>>(
                    model_ref, nsub,
                    uses_adaptive_aproximation_value(adaptive_approximation),
                    uses_recursive_aproximation_value(recursive_approximation),
                    uses_averaging_aproximation_value(averaging_approximation),
                    uses_variance_aproximation_value(variance_approximation),
                    uses_taylor_variance_correction_aproximation_value(
                        taylor_variance_correction_approximation),
                    uses_family_aproximation_value(family_approximation),
                    uses_qdt_method_value(0))
                    .get_variant(),
                // Macro taylor-variance-correction branch: taylor_qdt fixed false.
                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 1, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              true>,
                    var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                    decltype(model_ref),
                    var::constexpr_Var_domain<int, uses_qdt_method, 0>>(
                    model_ref, nsub,
                    uses_adaptive_aproximation_value(adaptive_approximation),
                    uses_recursive_aproximation_value(recursive_approximation),
                    uses_averaging_aproximation_value(averaging_approximation),
                    uses_variance_aproximation_value(variance_approximation),
                    uses_taylor_variance_correction_aproximation_value(
                        taylor_variance_correction_approximation),
                    uses_family_aproximation_value(family_approximation),
                    uses_qdt_method_value(0))
                    .get_variant()),
            // Micro branch: recursive=true, taylor_variance_correction=false,
            // averaging ∈ {0, 1, 2}, taylor_qdt ∈ {false, true} — the eigen vs
            // Taylor-expm flag for Q's matrix exponential. Only the micro branch
            // exposes both taylor_qdt values; macro paths fix it false to keep
            // the variant count manageable.
            Likelihood_Model_regular<
                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
                var::constexpr_Var_domain<int, uses_family_aproximation, family_micro>,
                decltype(model_ref),
                var::constexpr_Var_domain<int, uses_qdt_method, 0, 1, 2>>(
                model_ref, nsub,
                uses_adaptive_aproximation_value(adaptive_approximation),
                uses_recursive_aproximation_value(recursive_approximation),
                uses_averaging_aproximation_value(averaging_approximation),
                uses_variance_aproximation_value(variance_approximation),
                uses_taylor_variance_correction_aproximation_value(
                    taylor_variance_correction_approximation),
                uses_family_aproximation_value(family_approximation),
                uses_qdt_method_value(qdt_method_int))
                .get_variant()),
        // Nonlinear-least-squares (LSE) branch: the lean Moffatt & Hume 2007 JGP
        // fold — mean only, no Kalman covariance. adaptive/recursive/taylor_vc
        // fixed false, qdt fixed to eig (0), averaging ∈ {0, 1}. variance accepts
        // BOTH values but is ignored in compute: the live .macroir scripts
        // hardcode variance=true and cannot inject it, so a {false}-only domain
        // would leave merge_Maybe_variant with no runtime match.
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, false>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            var::constexpr_Var_domain<int, uses_family_aproximation, family_nonlinearsqr>,
            decltype(model_ref),
            var::constexpr_Var_domain<int, uses_qdt_method, 0>>(
            model_ref, nsub,
            uses_adaptive_aproximation_value(adaptive_approximation),
            uses_recursive_aproximation_value(recursive_approximation),
            uses_averaging_aproximation_value(averaging_approximation),
            uses_variance_aproximation_value(variance_approximation),
            uses_taylor_variance_correction_aproximation_value(
                taylor_variance_correction_approximation),
            uses_family_aproximation_value(family_approximation),
            uses_qdt_method_value(0))
            .get_variant());
}

// Thin bool wrapper preserved so existing scripts and C++ call sites keep
// working: maps the deprecated bool micro_approximation onto the int family
// selector (true → family_micro, false → family_macro).
inline auto build_likelihood_function(const ModelPtr& model0,
                                      bool adaptive_approximation, bool recursive_approximation,
                                      int averaging_approximation, bool variance_approximation,
                                      bool taylor_variance_correction_approximation,
                                      bool micro_approximation,
                                      bool taylor_qdt_approximation = false) {
    return build_likelihood_function_with_family(
        model0, adaptive_approximation, recursive_approximation, averaging_approximation,
        variance_approximation, taylor_variance_correction_approximation,
        micro_approximation ? family_micro : family_macro, taylor_qdt_approximation);
}

using likelihood_algorithm_type =
    var::untransformed_type_t<decltype(build_likelihood_function_with_family(
        std::declval<const ModelPtr&>(), false, false, 2, true, false, family_macro, false))>;

auto calculate_mlikelihood(const likelihood_algorithm_type& likelihood_algorithm,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r) -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

inline auto calculate_simulation_mlikelihood(
    const likelihood_algorithm_type& likelihood_algorithm,
    const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& r)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    return calculate_mlikelihood(likelihood_algorithm, par, e, get<Recording>(r()));
}

auto calculate_mdlikelihood(const likelihood_algorithm_type& likelihood_algorithm,
                            const var::Parameters_transformed& par, const Experiment& e,
                            const Recording& r) -> Maybe_error<dMacro_State_Hessian_minimal>;

auto calculate_mdiff_likelihood(const likelihood_algorithm_type& likelihood_algorithm,
                                const var::Parameters_transformed& par, const Experiment& e,
                                const Recording& r, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian>;

auto calculate_mlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_predictions>;

inline auto calculate_simulation_mlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Simulated_Recording<var::please_include<>>& r)
    -> Maybe_error<Macro_State_Ev_predictions>
    {
        return calculate_mlikelihood_predictions   (likelihood_algorithm, par, e,
                                             get<Recording>(r()));

    }



auto calculate_mlikelihood_diagnostics(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_diagnostic>;

inline auto calculate_simulated_mlikelihood_diagnostics(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par,const Experiment& e,
                                       const Simulated_Recording<var::please_include<>>& simulation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_mlikelihood_diagnostics(likelihood_algorithm, par, e,
                                             get<Recording>(simulation()));
}

// Memoization is on by default (per-call cache shared across intervals at the
// same θ) but can be disabled at runtime by setting MACRODR_MEMOIZE=0 — useful
// for correctness cross-checks against the cached path.
auto calculate_mdlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                        const var::Parameters_transformed& par, const Experiment& e,
                                        const Recording& r)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_mdlikelihood_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_mdlikelihood_predictions(likelihood_algorithm, par, e,
                                              get<Recording>(simulation()));
}

// Multi-simulation analytic dlogL evaluation. All n simulations share the same
// θ; the implementation in the cpp constructs a single FuncMap_St and reuses
// it across the loop, so the macro and micro Qdt caches stay populated across
// simulations. Set MACRODR_MEMOIZE=0 in the environment to disable memoization
// for correctness cross-checks.
auto calculate_n_simulation_mdlikelihood_predictions_impl(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<dMacro_State_Ev_gradient_all>>;

inline auto calculate_n_simulation_mdlikelihood_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<dMacro_State_Ev_gradient_all>> {
    return calculate_n_simulation_mdlikelihood_predictions_impl(likelihood_algorithm, par, e,
                                                                 simulation);
}

// Numerical Fisher information via the likelihood_algorithm_type variant —
// matches the script-level usage pattern of calc_dlikelihood_predictions.
auto calculate_mnumerical_fisher_information(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r,
    double h_rel) -> Maybe_error<parameter_spd_payload>;

inline auto calculate_simulation_mnumerical_fisher_information(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    double h_rel) -> Maybe_error<parameter_spd_payload> {
    return calculate_mnumerical_fisher_information(likelihood_algorithm, par, e,
                                                    get<Recording>(simulation()), h_rel);
}

// Per-replica per-sample numerical Fisher Information. Decomposes the global
// F per replica into one F_t contribution per timestep, via per-step FD on
// the cumulative AD score. By linearity of d/dθ on the additive logL, the
// sum over timesteps of F_t equals the global F returned by
// calculate_n_simulation_mnumerical_fisher_information at the same (par, h_rel).
//
// Output: one dMacro_State_Ev_per_sample_F per replica, whose Evolution_of<>
// slot is populated with the per-step F at each timestep. Caller writes the
// result with the existing write_csv overloads for vector<State> with
// Evolution (axes columns propagate via Indexed-aware variants).
//
// Cost: 2·n_params calls to calculate_n_simulation_mdlikelihood_predictions
// (each batched across all replicas — the FuncMap_St cache is reused). Same
// per-h_rel cost as calc_numerical_fisher_information but produces per-sample
// data instead of the aggregate.
auto calculate_per_sample_n_simulation_mnumerical_fisher_information(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation,
    double h_rel) -> Maybe_error<std::vector<dMacro_State_Ev_per_sample_F>>;

// Detailed per-sample diagnostic. Runs the dlikelihood ONCE at θ and dumps the
// SIGNIFICANT recursion variables (logL, y_mean, y_var, P_mean, P_Cov,
// trust_coefficient) per sample, each as a Derivative carrying the value X AND
// the regular AD gradient ∂X/∂θ over ALL parameters (no per-parameter
// perturbation — the Jacobian already covers every parameter). The perturbation
// for an FD-instability hunt is built OUTSIDE, by passing `par` as an
// Indexed<Parameters_transformed> over the perturbation axes (see
// by_parameter_coordinate + apply_relative_perturbation); the DSL lifts this
// command over them.
//
// Each replica's Evolution is windowed to the inclusive ABSOLUTE sample range
// [sample_min, sample_max] BEFORE returning — this caps memory when the command
// is lifted over many perturbations (2·p legs × replicas), instead of holding
// full Evolutions until write time. The matching write_csv overload takes the
// same sample_min and maps the windowed elements back to absolute samples.
auto calculate_n_simulation_mdetailed_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation,
    std::size_t sample_min,
    std::size_t sample_max) -> Maybe_error<std::vector<dMacro_State_Ev_detailed>>;

inline auto calculate_n_simulation_mnumerical_fisher_information(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation,
    double h_rel,
    std::size_t decimate) -> Maybe_error<std::vector<parameter_spd_payload>> {
    std::size_t step = decimate == 0 ? 1 : decimate;
    // Decimated simulation indices to evaluate F at.
    std::vector<std::size_t> sel;
    for (std::size_t i = 0; i < simulation.size(); i += step) sel.push_back(i);
    // calculate_mnumerical_fisher_information builds its own function table per
    // call (no shared state — see calculate_mdlikelihood_predictions), so this
    // loop parallelizes with no fork needed. Active OpenMP level only when the
    // combo loop is serial (MACRODR_AXIS_SERIAL=1); else it nests and serializes.
    std::vector<std::optional<parameter_spd_payload>> slots(sel.size());
    std::vector<std::string> errs(sel.size());
#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < sel.size(); ++k) {
        auto res = calculate_mnumerical_fisher_information(
            likelihood_algorithm, par, e, get<Recording>(simulation[sel[k]]()), h_rel);
        if (res)
            slots[k].emplace(std::move(res.value()));
        else
            errs[k] = res.error()();
    }
    std::vector<parameter_spd_payload> results;
    results.reserve(sel.size());
    for (std::size_t k = 0; k < sel.size(); ++k) {
        if (!errs[k].empty()) return error_message(errs[k]);
        results.push_back(std::move(*slots[k]));
    }
    return results;
}


auto calculate_likelihood(const ModelPtr& model0,
                          const var::Parameters_transformed& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation,
                          bool micro_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

inline auto calculate_simulation_likelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    return calculate_likelihood(model0, par, e, get<Recording>(simulation()),
                                adaptive_approximation, recursive_approximation,
                                averaging_approximation, variance_approximation,
                                taylor_variance_correction_approximation, micro_approximation);
}

auto calculate_dlikelihood(const ModelPtr& model0,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation,
                           bool micro_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal>;

inline auto calculate_simulation_dlikelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal> {
    return calculate_dlikelihood(model0, par, e, get<Recording>(simulation()),
                                 adaptive_approximation, recursive_approximation,
                                 averaging_approximation, variance_approximation,
                                 taylor_variance_correction_approximation, micro_approximation);
}
auto calculate_diff_likelihood(const ModelPtr& model0,
                               const var::Parameters_transformed& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation,
                               bool micro_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian>;

inline auto calculate_simulation_diff_likelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    return calculate_diff_likelihood(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation, delta_param);
}

auto calculate_likelihood_predictions(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation,
                                      bool micro_approximation)
    -> Maybe_error<Macro_State_Ev_predictions>;

inline auto calculate_simulation_likelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<Macro_State_Ev_predictions> {
    return calculate_likelihood_predictions(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation);
}

auto calculate_likelihood_diagnostics(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation,
                                      bool micro_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic>;

inline auto calculate_simulation_likelihood_diagnostics(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_likelihood_diagnostics(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation);
}

auto calculate_dlikelihood_predictions(const ModelPtr& model0,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation,
                                       bool micro_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation);
}

inline auto calculate_simulation_sub_dlikelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e,
    const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation);
}

auto calculate_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation,
    bool micro_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& simulation, bool adaptive_approximation,
    bool recursive_approximation, int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation,
    bool micro_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions_model(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation);
}

// Numerical (observed) Fisher information for one recording. Returns the
// positive-definite F = -∂²(logL)/∂θ² as a parameter_spd_payload. Cost:
// 2·n_params calls to the dlikelihood path at θ ± h_i·e_i. See the
// implementation comment in src/core/likelihood.cpp for the central-difference
// derivation. h_rel defaults to 1e-5 (textbook double-precision). Decimation,
// if desired, is the caller's job — pass a subsetted Recording.
auto calculate_numerical_fisher_information(
    const ModelPtr& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation, bool micro_approximation,
    double h_rel) -> Maybe_error<parameter_spd_payload>;

inline auto calculate_simulation_numerical_fisher_information(
    const ModelPtr& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& simulation, bool adaptive_approximation,
    bool recursive_approximation, int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation, bool micro_approximation,
    double h_rel) -> Maybe_error<parameter_spd_payload> {
    return calculate_numerical_fisher_information(
        model0, par, e, get<Recording>(simulation()), adaptive_approximation,
        recursive_approximation, averaging_approximation, variance_approximation,
        taylor_variance_correction_approximation, micro_approximation, h_rel);
}

// Multi-recording form: one F matrix per simulation, parallel to the dy vector
// produced by calculate_n_simulation_mdlikelihood_predictions. Bootstrap
// consumes this independently from dy (independent index draws). `decimate`
// selects every K-th simulation for the F computation — the bootstrap then
// resamples a vector of length ⌈n_sims / decimate⌉. decimate=1 → all sims.
inline auto calculate_n_simulation_numerical_fisher_information(
    const ModelPtr& model0, const var::Parameters_transformed& par, const Experiment& e,
    const std::vector<Simulated_Recording<var::please_include<>>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation,
    bool micro_approximation, double h_rel,
    std::size_t decimate) -> Maybe_error<std::vector<parameter_spd_payload>> {
    std::size_t step = decimate == 0 ? 1 : decimate;
    std::vector<parameter_spd_payload> out;
    out.reserve((simulation.size() + step - 1) / step);
    for (std::size_t i = 0; i < simulation.size(); i += step) {
        auto res = calculate_simulation_numerical_fisher_information(
            model0, par, e, simulation[i], adaptive_approximation, recursive_approximation,
            averaging_approximation, variance_approximation,
            taylor_variance_correction_approximation, micro_approximation, h_rel);
        if (!res) return res.error();
        out.push_back(std::move(res.value()));
    }
    return out;
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path);

template <class... Vs>
inline Maybe_error<std::string> write_csv(var::Vector_Space<Vs...> const& lik, std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

template <class T>
inline Maybe_error<std::string> write_csv(var::Indexed<T> const& indexed, std::string path) {
    return detail::write_summary_csv(indexed, std::move(path), "summary");
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
                                   std::vector<TMacro_State<vVars...>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    std::vector<TMacro_State<vVars...>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path);

// Windowed sibling of the (Experiment, Indexed-sims, Indexed-states) overload:
// emits only the per-step Evolution rows for the inclusive absolute sample
// window [sample_min, sample_max]. Used by the detailed FD-instability
// diagnostic to keep the CSV small while keeping absolute sample indices /
// times / agonist / patch_current correct.
template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::size_t sample_min,
    std::size_t sample_max, std::string path);

template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path);

// Indexed state-batch with NO Experiment/Simulation. Emits axis columns +
// simulation_index per replica + the state's Evolution (sample_index per
// timestep). Used for the per-sample numerical Fisher CSV: passing the
// Indexed states directly (no plain-sim arg) lets the DSL match an Indexed
// param without lifting, so the axis columns survive into the CSV.
template <template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path);

// Per-replica summary CSV: emits ALL slots of each state EXCEPT Evolution_of
// (skipped via CsvContext::skip_evolution=true so the Evolution branch in
// emit_any returns early), and tags the numerical Fisher per replica as
// Likelihood_Numerical_Fisher_Information. Designed for bug-hunting on
// numeric_Fisher_information variability across replicates.
template <template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    std::vector<TMacro_State<vVars...>> const& dlikelihood_predictions,
    std::vector<parameter_spd_payload> const& numerical_fisher_information,
    std::string path);

// Indexed (multi-axis) overload of the above. Iterates the axis space and
// populates the CSV's axis_values columns per coord so the output matches the
// schema of the analysis CSV (algorithm, noise_in_conductance_tau, Num_ch,
// interval_in_tau, simulation_algorithm, axis_h_fim, ...). Used when the DSL
// hands dlikelihood_predictions / numerical_fisher_information as Indexed
// values from a multi-axis context.
template <template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& dlikelihood_predictions,
    var::Indexed<std::vector<parameter_spd_payload>> const& numerical_fisher_information,
    std::string path);


auto calculate_boot_Likelihood_diagnostics(const std::vector<dMacro_State_Ev_gradient_all>& dy,
               const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<Macro_State_Ev_diagnostic>>;


using Analisis_derivative_diagnostic_base = var::Vector_Space
<
        Probit_statistics<Moment_statistics<Sum<logL>, false>>,
        Probit_statistics<Moment_statistics<Sum<elogL>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::r_std>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::r2_std>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::trust_coefficient>, false>>,
        // TODO: Probit aggregations for taylor_trust_coefficient / taylor_vSv /
        // taylor_strength (added to Algo_State_Dynamic per Layer 1+2). Requires
        // mirroring the addition into the internal sum_moments / evol_moments
        // type_identity Vector_Spaces in
        // src/core/likelihood.cpp:calculate_Likelihood_diagnostics_preset_f
        // (~line 2521+) plus matching extractor lambdas. Per-interval values are
        // already exposed via the dynamic-state framework; aggregations can be
        // wired when figure_2 results call for per-trace aggregations of α_vSv.
        Probit_statistics<Moment_statistics<Sum<dlogL>, true>>,
        Probit_statistics<Moment_statistics<Sum<Gaussian_Fisher_Information>, false>>,
        Probit_statistics<Sum<Moment_statistics<macrodr::r_std, false>>>,
        Probit_statistics<Sum<Moment_statistics<dlogL, true>>>,
        Probit_statistics<Sum<Moment_statistics<Gaussian_Fisher_Information, false>>>,
        Probit_statistics<Moment_statistics<evaluation_time, false>>,
        Probit_statistics<Likelihood_Numerical_Fisher_Information>,
        Probit_statistics<Likelihood_Information_Distortion>,
        Probit_statistics<log_Det<Likelihood_Information_Distortion>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Information_Distortion>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Information_Distortion>>,
        Probit_statistics<Max_Eigenvalue<Likelihood_Information_Distortion>>,
        Probit_statistics<Min_Eigenvalue<Likelihood_Information_Distortion>>,
        Probit_statistics<Mean_Log_Eigenvalue<Likelihood_Information_Distortion>>,
        Probit_statistics<Log_Eigenvalue_Variance<Likelihood_Information_Distortion>>,
        Probit_statistics<Affine_Invariant_Distance<Likelihood_Information_Distortion>>,
        Probit_statistics<Symmetrized_KL_Distortion<Likelihood_Information_Distortion>>,
        Probit_statistics<Likelihood_Gaussian_Information_Distortion>,
        Probit_statistics<log_Det<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Max_Eigenvalue<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Min_Eigenvalue<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Mean_Log_Eigenvalue<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Log_Eigenvalue_Variance<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Affine_Invariant_Distance<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Symmetrized_KL_Distortion<Likelihood_Gaussian_Information_Distortion>>,
        Probit_statistics<Likelihood_Gaussian_Fisher_Distortion>,
        Probit_statistics<log_Det<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Max_Eigenvalue<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Min_Eigenvalue<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Mean_Log_Eigenvalue<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Log_Eigenvalue_Variance<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Affine_Invariant_Distance<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Symmetrized_KL_Distortion<Likelihood_Gaussian_Fisher_Distortion>>,
        Probit_statistics<Likelihood_Information_Distortion_Reconstituted>,
        Probit_statistics<Likelihood_Sample_Distortion>,
        Probit_statistics<log_Det<Likelihood_Sample_Distortion>>,
        Probit_statistics<Gaussian_Sample_Distortion>,
        Probit_statistics<log_Det<Gaussian_Sample_Distortion>>,
        Probit_statistics<Likelihood_Correlation_Distortion>,
        Probit_statistics<log_Det<Likelihood_Correlation_Distortion>>,
        Probit_statistics<Likelihood_Fisher_Covariance>,
        Probit_statistics<log_Det<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Correlation_Of<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Effective_Rank<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Null_Space_Projector<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Likelihood_Fisher_Covariance>>,
        Probit_statistics<Gaussian_Fisher_Covariance>,
        Probit_statistics<log_Det<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Correlation_Of<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Effective_Rank<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Null_Space_Projector<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Gaussian_Fisher_Covariance>>,
        Probit_statistics<Likelihood_Distortion_Corrected_Covariance>,
        Probit_statistics<log_Det<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Correlation_Of<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Effective_Rank<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Null_Space_Projector<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Gaussian_Distortion_Corrected_Covariance>,
        Probit_statistics<log_Det<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Correlation_Of<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Effective_Rank<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Null_Space_Projector<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Gaussian_Distortion_Corrected_Covariance>>,
        Probit_statistics<Likelihood_Distortion_Induced_Bias>,
        Probit_statistics<Gaussian_Distortion_Induced_Bias>>;

// Per-sample derived diagnostics: F- and G-anchored Sample_Distortion and
// Distortion_Induced_Bias evaluated at each sample. Included by presets that
// ask for trace-local detail (series_var, series_cov, series_kernel_full).
// Field order MUST match build_per_sample_derived's derived_vs in
// src/core/likelihood.cpp.
using Per_sample_derived_diagnostics =
    Probit_statistics<macrodr::Evolution_of<var::Vector_Space<
        Likelihood_Sample_Distortion, Likelihood_Distortion_Induced_Bias,
        Gaussian_Sample_Distortion, Gaussian_Distortion_Induced_Bias>>>;

// Preset 1: basic — one integral correlation lag per headline observable.
// Used for lattice-scale runs; keeps output at megabyte scale.
using Analisis_derivative_diagnostic_basic = var::concatenate_t<
    Analisis_derivative_diagnostic_base,
    var::Vector_Space<
        Probit_statistics<Report_integral<logL>>,
        Probit_statistics<Report_integral<macrodr::y_mean>>,
        Probit_statistics<Report_integral<macrodr::y_var>>,
        Probit_statistics<Report_integral<macrodr::r_std>>,
        Probit_statistics<Report_integral<dlogL>>>>;

// Preset 2: series_var — per-sample moments (diagonal variance only), plus
// forward and integral lags, for 7 leaf observables + Per_sample_derived.
using Analisis_derivative_diagnostic_series_var = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_local_var<logL>>,
            Probit_statistics<Report_local_var<elogL>>,
            Probit_statistics<Report_local_var<macrodr::y_mean>>,
            Probit_statistics<Report_local_var<macrodr::y_var>>,
            Probit_statistics<Report_local_var<macrodr::r_std>>,
            Probit_statistics<Report_local_var<macrodr::trust_coefficient>>,
            Probit_statistics<Report_local_var<dlogL>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;

// Preset 3: series_cov — same as series_var but with full per-sample
// covariance blocks for vector V's; adds GFI at local_var level
// (matrix-of-matrix covariance is impractical for GFI).
using Analisis_derivative_diagnostic_series_cov = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_local_cov<logL>>,
            Probit_statistics<Report_local_cov<elogL>>,
            Probit_statistics<Report_local_cov<macrodr::y_mean>>,
            Probit_statistics<Report_local_cov<macrodr::y_var>>,
            Probit_statistics<Report_local_cov<macrodr::r_std>>,
            Probit_statistics<Report_local_cov<macrodr::trust_coefficient>>,
            Probit_statistics<Report_local_cov<dlogL>>,
            Probit_statistics<Report_local_var<Gaussian_Fisher_Information>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;

// Preset 4: series_kernel — full cross-correlation kernel for the 5 core
// observables. No auxiliaries (elogL, trust, GFI), no Per_sample_derived.
// Focused single-condition mechanism view.
using Analisis_derivative_diagnostic_series_kernel = var::concatenate_t<
    Analisis_derivative_diagnostic_base,
    var::Vector_Space<
        Probit_statistics<Report_cross<logL>>,
        Probit_statistics<Report_cross<macrodr::y_mean>>,
        Probit_statistics<Report_cross<macrodr::y_var>>,
        Probit_statistics<Report_cross<macrodr::r_std>>,
        Probit_statistics<Report_cross<dlogL>>>>;

// Preset 5: series_kernel_full — exhaustive single-condition: kernel for 7
// leaf observables, GFI at local_var, and Per_sample_derived included.
using Analisis_derivative_diagnostic_series_kernel_full = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_cross<logL>>,
            Probit_statistics<Report_cross<elogL>>,
            Probit_statistics<Report_cross<macrodr::y_mean>>,
            Probit_statistics<Report_cross<macrodr::y_var>>,
            Probit_statistics<Report_cross<macrodr::r_std>>,
            Probit_statistics<Report_cross<macrodr::trust_coefficient>>,
            Probit_statistics<Report_cross<dlogL>>,
            Probit_statistics<Report_local_var<Gaussian_Fisher_Information>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;


// =============================================================================
// MLE_Run — the COMPLETE per-group Gauss-Newton refit record (one Vector_Space
// per refit group of n replicas). Reuses the GN result_value dlogPs slots at θ̂
// (logL / Grad / Gaussian_Fisher_Information — same meaning, same types) plus
// θ̂ itself and the run diagnostics. Built as a Vector_Space so the existing
// write_csv machinery emits its columns; the gradient Grad = ∇logL(θ̂) is ≈0 at
// a clean MLE, a direct convergence probe.
using MLE_Run = var::Vector_Space<
    Model_Parameters_Hat,          // θ̂_group (p×1, parameter-named)
    logL,                          // maximum logL at θ̂  (result.max_value)
    Grad,                          // ∇logL(θ̂) (p×1, named) — the gradient
    Gaussian_Fisher_Information,   // G_group at θ̂
    Newton_Iterations,             // GN iterations to convergence
    Convergence_Status_Code,       // status → numeric code (see the tag)
    Recovery_Distance>;            // ‖θ̂ − θ_reference‖₂

// Cloud slot holding every group's MLE_Run. Marked is_write_csv_companion_slot
// so write_csv emits it to the <path>_runs.csv companion (one row per group ×
// component, group keyed by sample_index) and the flat summary writer SKIPS it
// (a vector of Vector_Space records would otherwise collide the group index with
// the per-parameter value_row). See detail::write_runs_csv / emit_any.
struct MLE_Run_Records : public var::Constant<MLE_Run_Records, std::vector<MLE_Run>> {
    using base_type = var::Constant<MLE_Run_Records, std::vector<MLE_Run>>;
    using base_type::base_type;
    static constexpr bool is_write_csv_companion_slot = true;
    MLE_Run_Records() = default;
    friend std::string className(MLE_Run_Records) { return "MLE_Run_Records"; }
};

// =============================================================================
// MLE per-group-of-replicates output (Path A: minimal State) — SLIM "cloud".
//
// Composite output of calc_MLE_per_group_of_replicates at one group_size n.
// This is the SLIMMED output: only the θ̂ cloud (mean + empirical covariance),
// the recovery Wald T² of (θ̄ − θ_reference) in the aggregated G_lik metric, a
// BARE point estimate of θ̄ (the handoff for get_parameters_mean), and the
// (currently empty) representative-group sample stub. The heavy diagnostic
// battery (numerical Fisher F̄_b, FC family, GFD, IDM, DCC, the Empirical-
// Covariance distortions, eigenvector slots) is NOT computed here anymore — it
// is recomputed downstream by composing the existing figure_2 commands at θ̄
// (calc_numerical_fisher_information / calc_dlikelihood_predictions /
// likelihood_derivative_basic_diagnostics_paired). The distortion / eigenvector
// TAG types are intentionally retained (distributions.h / moment_statistics.h)
// for that downstream reuse; only THIS command stops computing them.
//
// Templated on State so the saved probit samples carry the State the command
// was instantiated with — for Path A this is dMacro_State_Hessian_minimal_param.
template <class State>
using MLE_Group_Cloud = var::Vector_Space<
    Group_Size,

    // θ̂ cloud over the bootstrapped groups: θ̄ (mean) + Cov_emp (covariance):
    Probit_statistics<Moment_statistics<Model_Parameters_Hat, true>>,

    // Recovery Wald T² = Δθᵀ·Ḡ_lik·Δθ, Δθ = θ̄ − theta_reference, Ḡ_lik = group-
    // mean of the per-group Gaussian-formula FIM (the precision-like curvature):
    Probit_statistics<Wald_T2<Gaussian_Fisher_Information>>,

    // BARE point estimate of θ̄ (full-sample mean over ALL valid groups, NOT
    // bootstrapped) — the handoff consumed by get_parameters_mean:
    Model_Parameters_Hat,

    // BARE raw full-sample empirical covariance Cov_emp of the θ̂ cloud over ALL
    // valid groups (the covariance of the SAME all-group Moment_statistics whose
    // mean is the θ̄ point above). Handoff consumed by calc_empirical_distortion:
    Empirical_Parameter_Covariance,

    // Complete per-group GN refit log (θ̂ + gradient + Fisher + status + n_iter +
    // recovery distance), one MLE_Run per group. Companion-only: emitted to
    // <path>_runs.csv, skipped in the summary:
    MLE_Run_Records,

    // K representative groups per ranking variable at each probit height —
    // full per-group State preserved for downstream inspection (empty stub):
    Probit_Samples_at_Group_Size<State>>;

// Run per-group MLE optimisation at a single group_size, then bootstrap-
// aggregate the per-group F + state into the figure_2 battery applied at
// θ̂_group (NOT at fixed θ_sim — F is per-group). Save K representative
// groups per ranking variable at the requested probit heights for later
// inspection (full state preserved at each saved group).
//
// group_size = 1 recovers per-replicate analysis. Larger group_size pools n
// recordings into one MLE optimisation via combined-likelihood. Independent
// runs at multiple group_sizes are script-level (call this command once per n).
//
// Bootstrap is on GROUPS, not on recordings — N_groups = N_total / group_size.
// If N_groups < min_groups_for_bootstrap, bootstrap-derived probit slots are
// NaN-filled; point estimates and probit samples are still produced.
//
// theta_warmstart = the GN initial guess; theta_reference = the value θ̂ is
// tested against in the Wald T² (θ₀, typically the simulation truth θ_sim).
// They are decoupled so the optimizer may start anywhere while the Wald tests
// recovery of the true generating parameter.
//
// State template parameter (used for the saved probit samples only; the GN
// inner loop always uses dMacro_State_Hessian_minimal for speed):
//   - dMacro_State_Hessian_minimal_param : Path A, minimal memory footprint
//   - dMacro_State_Ev_gradient_all_param : Path B, full Evolution_of preserved
template <class State>
auto calc_MLE_per_group_of_replicates(
    const likelihood_algorithm_type& likelihood_algorithm,
    const var::Parameters_transformed& theta_warmstart,
    const var::Parameters_transformed& theta_reference,
    const Experiment& experiment,
    const std::vector<Recording>& recordings,
    std::size_t group_size,
    std::size_t n_bootstrap_samples,
    std::size_t min_groups_for_bootstrap,
    const std::set<double>& probit_cis,
    const std::set<double>& probit_sample_heights,
    const std::vector<std::string>& ranking_variables,
    std::size_t seed,
    const macrodr::optimization::gauss_newton_options& gn_opts,
    double F_h_relative = 1e-5)
    -> Maybe_error<MLE_Group_Cloud<State>>;

// Rebuild θ̄ (the bare full-sample mean of θ̂) from the cloud's Model_Parameters_Hat
// slot as a Parameters_transformed bound to the same parameter metadata. The
// handoff to the downstream figure_2 battery (calc_*_fisher / calc_dlikelihood
// at θ̄). The slot's parameter_vector_payload carries the Matrix (value()) and
// the Parameters_transformed* (parameters_ptr()); create() rebinds the names.
template <class State>
auto get_parameters_mean(const MLE_Group_Cloud<State>& cloud)
    -> Maybe_error<var::Parameters_transformed>;

// =============================================================================
// Empirical-vs-theoretical distortion capstone (figure_3 Fase 2).
//
// Composes the MLE cloud's bare empirical covariance Cov_emp with the figure_2
// outputs (numerical Fisher at θ_sim and at θ̄, score states at θ̄) into three
// distortion matrices, each with the same scalar suite the diagnostic battery
// emits:
//   - Empirical_Covariance_Fisher_Distortion  = F̄_bar^{1/2} · Cov_emp · F̄_bar^{1/2}
//   - Empirical_Covariance_Corrected_Distortion ~ DCC^{-1/2} · Cov_emp · DCC^{-1/2},
//                                                  DCC = F̄_bar^{-1} J F̄_bar^{-1}
//     computed WITHOUT forming/decomposing the doubly-ill-conditioned DCC, via the
//     identity eig(DCC^{-1/2}·Cov_emp·DCC^{-1/2}) = eig(IDM^{-1}·ECD_Fisher) as
//     IDM^{-1/2}·ECD_Fisher·IDM^{-1/2}, IDM = F̄_bar^{-1/2} J F̄_bar^{-1/2} (≈ I ⟺
//     Cov_emp ≈ DCC; the test that bites when IDM ≠ I, i.e. under misspecification)
//   - Optimum_Fisher_Distortion = F̄_bar^{-1/2} · F̄_sim · F̄_bar^{-1/2}
// where F̄_bar / F̄_sim are the means of the per-recording numerical Fisher
// vectors and J is the HAC score covariance at θ̄ (the same J the battery
// computes via get<covariance<Sum<dlogL>>>). Point estimates only — NO
// bootstrap in v1 (the cloud already carries the bootstrap of the θ̂ cloud; this
// command is a deterministic comparison of aggregate matrices).
using Empirical_Distortion_Analysis = var::Vector_Space<
    Empirical_Covariance_Fisher_Distortion,
    Affine_Invariant_Distance<Empirical_Covariance_Fisher_Distortion>,
    log_Det<Empirical_Covariance_Fisher_Distortion>,
    Eigenvalue_Spectrum<Empirical_Covariance_Fisher_Distortion>,
    Spectrum_Condition_Number<Empirical_Covariance_Fisher_Distortion>,
    Empirical_Covariance_Corrected_Distortion,
    Affine_Invariant_Distance<Empirical_Covariance_Corrected_Distortion>,
    log_Det<Empirical_Covariance_Corrected_Distortion>,
    Eigenvalue_Spectrum<Empirical_Covariance_Corrected_Distortion>,
    Spectrum_Condition_Number<Empirical_Covariance_Corrected_Distortion>,
    Optimum_Fisher_Distortion,
    Affine_Invariant_Distance<Optimum_Fisher_Distortion>,
    log_Det<Optimum_Fisher_Distortion>,
    Eigenvalue_Spectrum<Optimum_Fisher_Distortion>,
    Min_Eigenvalue<Optimum_Fisher_Distortion>,  // the key "FIM_sim indefinite?" readout
    Spectrum_Condition_Number<Optimum_Fisher_Distortion>>;

// Bootstrap output: the point estimate (over ALL groups) CONCATENATED with the
// per-component bootstrap CIs (every distortion matrix entry, every sorted
// eigenvalue, every scalar wrapped in Probit_statistics — apply_Probit_statistics
// produces exactly var::Vector_Space<Probit_statistics<Vs>...>). The bare point
// and the Probit_statistics<bare> are DISTINCT types, so the flat concatenation
// has no component-path collision (mirrors MLE_Group_Cloud's bare ++ probit).
template <class VS>
struct probit_wrapped;
template <class... Vs>
struct probit_wrapped<var::Vector_Space<Vs...>> {
    using type = var::Vector_Space<Probit_statistics<Vs>...>;
};
template <class VS>
using probit_wrapped_t = typename probit_wrapped<std::remove_cvref_t<VS>>::type;

using Empirical_Distortion_Bootstrap =
    decltype(concatenate(std::declval<Empirical_Distortion_Analysis>(),
                         std::declval<probit_wrapped_t<Empirical_Distortion_Analysis>>()));

// Gaussian-Fisher-anchored twin of Empirical_Distortion_Analysis (same 16-field
// 5/5/6 shape; Min_Eigenvalue only on the Gaussian Optimum). Anchored on the
// Gaussian Fisher Ḡ built from the dlikelihood, so it needs NO numerical Fisher.
using Empirical_Distortion_Gaussian_Analysis = var::Vector_Space<
    Empirical_Covariance_Gaussian_Distortion,
    Affine_Invariant_Distance<Empirical_Covariance_Gaussian_Distortion>,
    log_Det<Empirical_Covariance_Gaussian_Distortion>,
    Eigenvalue_Spectrum<Empirical_Covariance_Gaussian_Distortion>,
    Spectrum_Condition_Number<Empirical_Covariance_Gaussian_Distortion>,
    Empirical_Covariance_Gaussian_Corrected_Distortion,
    Affine_Invariant_Distance<Empirical_Covariance_Gaussian_Corrected_Distortion>,
    log_Det<Empirical_Covariance_Gaussian_Corrected_Distortion>,
    Eigenvalue_Spectrum<Empirical_Covariance_Gaussian_Corrected_Distortion>,
    Spectrum_Condition_Number<Empirical_Covariance_Gaussian_Corrected_Distortion>,
    Optimum_Gaussian_Distortion,
    Affine_Invariant_Distance<Optimum_Gaussian_Distortion>,
    log_Det<Optimum_Gaussian_Distortion>,
    Eigenvalue_Spectrum<Optimum_Gaussian_Distortion>,
    Min_Eigenvalue<Optimum_Gaussian_Distortion>,  // sim-vs-pool frame gap (NOT indefiniteness)
    Spectrum_Condition_Number<Optimum_Gaussian_Distortion>>;

using Empirical_Distortion_Gaussian_Bootstrap =
    decltype(concatenate(
        std::declval<Empirical_Distortion_Gaussian_Analysis>(),
        std::declval<probit_wrapped_t<Empirical_Distortion_Gaussian_Analysis>>()));

// Empirical-vs-theoretical capstone command (figure_3 Fase 2). State templates
// only the cloud (the figure_2 vectors are State-agnostic).
//   fim_sim  : numerical Fisher at θ_sim, per recording (decimate=1)
//   fim_bar  : numerical Fisher at θ_pool (the anchor), per recording (decimate=1)
//   dlik_bar : score states at the anchor, per recording (the J / HAC Ω source)
// Bootstrap is over GROUPS (the cloud's MLE_Run_Records θ̂ are the iid units);
// each resampled group is EXPANDED to its group_size recordings (contiguous,
// reconstructed from Group_Size) for the per-recording F̄ / J. Returns the point
// estimate ++ per-component probit CIs.
template <class State>
auto calc_empirical_distortion(
    const MLE_Group_Cloud<State>& cloud,
    const std::vector<parameter_spd_payload>& fim_sim,
    const std::vector<parameter_spd_payload>& fim_bar,
    const std::vector<dMacro_State_Ev_gradient_all>& dlik_bar,
    std::size_t n_bootstrap, std::size_t seed, const std::set<double>& probit_cis,
    double rtol = 1e-10, double atol = 0.0)
    -> Maybe_error<Empirical_Distortion_Bootstrap>;

// Gaussian-Fisher-anchored twin of calc_empirical_distortion — NO numerical
// Fisher. The anchor Ḡ and the Optimum numerator are the Gaussian Fisher
// Σ_t GFI_t (× gsize) built from the dlikelihood:
//   dlik_sim : score states at θ_sim   (Optimum-Gaussian numerator, Ḡ_sim)
//   dlik_bar : score states at θ̄       (anchor Ḡ_bar + the J / score-cov source)
// Same GROUP-scale + group bootstrap as the F version. Returns point ++ probit CIs.
template <class State>
auto calc_empirical_distortion_gaussian(
    const MLE_Group_Cloud<State>& cloud,
    const std::vector<dMacro_State_Ev_gradient_all>& dlik_sim,
    const std::vector<dMacro_State_Ev_gradient_all>& dlik_bar,
    std::size_t n_bootstrap, std::size_t seed, const std::set<double>& probit_cis,
    double rtol = 1e-10, double atol = 0.0)
    -> Maybe_error<Empirical_Distortion_Gaussian_Bootstrap>;

// =============================================================================

auto calculate_Likelihood_derivative_basic_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_basic;

// Paired-bootstrap variant: requires F_per_recording.size() == dy.size()
// (i.e. caller must use decimate=1 when building F_per_recording). Uses one
// shared index vector per bootstrap replicate, removing the independent-draw
// variance inflation that makes IDM = F⁻¹ᐟ²·J·F⁻¹ᐟ² blow up in the worst
// eigendirection. Yields tighter ratio-of-quadratic-forms diagnostics under
// the information identity J ≈ F at the recording level.
//
// samples_dir (default ""): when non-empty, the raw bootstrap samples are
// also written as flat binary doubles inside this directory (one .bin per
// observable: IDM, GFD, SDM, CDM, FC, DCC, IDR, F_b matrices; their
// eigenvalue spectra; the new IDM/GFD scalar suite; log-dets; condition
// numbers; effective ranks). Sidecar files in the same directory:
// `manifest.csv` (observable shapes), `cells.csv` (one row per cell call),
// `reader.R` (auto-generated R reader). User joins binary cells to the
// analysis CSV by row order in R.
auto calculate_Likelihood_derivative_basic_diagnostics_paired(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag,
    std::string samples_dir = "")
    -> Analisis_derivative_diagnostic_basic;

auto calculate_Likelihood_derivative_series_var_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_var;

auto calculate_Likelihood_derivative_series_cov_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_cov;

auto calculate_Likelihood_derivative_series_kernel_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel;

auto calculate_Likelihood_derivative_series_kernel_full_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy,
    const std::vector<parameter_spd_payload>& F_per_recording, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel_full;


inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_base const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_basic const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_var const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_cov const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_kernel const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(
    Analisis_derivative_diagnostic_series_kernel_full const& lik, std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

// Dedicated overload for the per-group MLE cloud. MLE_Group_Cloud<State>
// is a var::Vector_Space alias, so it would otherwise be ambiguous between the
// generic write_csv(Vector_Space<Vs...>) and write_csv(TMacro_State<vVars...>)
// templates; templating on State (not on the full pack) makes this overload the
// more specialized — and therefore unambiguous — match.
template <class State>
inline Maybe_error<std::string> write_csv(MLE_Group_Cloud<State> const& lik,
                                          std::string path) {
    // The cloud's only CSV output is the raw per-group run log (<path>_runs.csv):
    // one row per (group × component), the group keyed by sample_index. The reduced
    // view (θ̄ / Cov_emp / bootstrap CIs / Wald T²) is fully derivable from it in R,
    // so no separate reduced summary is written (raw in C++, derived in R).
    auto param_names = detail::get_param_names_if_any(lik);
    return detail::write_runs_csv(get<MLE_Run_Records>(lik)(), param_names, path + "_runs");
}

// write_csv for the INDEXED per-group MLE cloud. figure_3 broadcasts the cloud
// over the algorithm / interval axes, so the DSL hands write_csv an
// Indexed<MLE_Group_Cloud<State>>. The cloud's only CSV output is the per-group
// run log (<path>_runs.csv) with the cell's axis columns, iterating each indexed
// cell's MLE_Run_Records (one row per cell × group × component). The reduced view
// is derived in R from it (raw in C++, derived in R), so no reduced summary is
// written. Registered explicitly in command_manager.cpp for the Indexed cloud case.
template <class State>
inline Maybe_error<std::string> write_csv_indexed_cloud(
    var::Indexed<MLE_Group_Cloud<State>> const& indexed, std::string path) {
    // Parameter names from the first cell's cloud (shared θ̂ metadata frame).
    std::vector<std::string> param_names;
    auto maybe_first = indexed.begin();
    if (maybe_first) {
        auto first_cell = indexed.at(maybe_first.value());
        if (first_cell) param_names = detail::get_param_names_if_any(first_cell.value().get());
    }

    return detail::write_indexed_rows_csv(
        indexed.index_space(), param_names, path + "_runs",
        [&indexed](auto& writer, detail::CsvContext base_ctx,
                   auto const& coord) -> Maybe_error<bool> {
            auto cell = indexed.at(coord);
            if (!cell) return cell.error();
            auto const& recs = get<MLE_Run_Records>(cell.value().get())();
            for (std::size_t g = 0; g < recs.size(); ++g) {
                detail::CsvContext ctx = base_ctx;
                ctx.scope = "mle_run";
                ctx.sample_index = g;
                auto ok = detail::emit_any(writer, std::move(ctx), recs[g]);
                if (!ok || !ok.value()) return ok;
            }
            return Maybe_error<bool>(true);
        });
}

// write_csv for the figure_3 Fase-2 empirical-distortion capstone output.
inline Maybe_error<std::string> write_csv(Empirical_Distortion_Analysis const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

// write_csv for the BOOTSTRAPPED capstone (point estimate ++ per-component probit
// CIs). It is a flat var::Vector_Space, so the generic summary writer emits the
// bare point components and the Probit_statistics<·> components (with the
// probit/quantile columns) side by side, distinguished by component_path.
inline Maybe_error<std::string> write_csv(Empirical_Distortion_Bootstrap const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

// write_csv for the Gaussian-Fisher-anchored capstone (mle_G variant).
inline Maybe_error<std::string> write_csv(Empirical_Distortion_Gaussian_Analysis const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

inline Maybe_error<std::string> write_csv(Empirical_Distortion_Gaussian_Bootstrap const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}


}  // namespace macrodr::cmd
