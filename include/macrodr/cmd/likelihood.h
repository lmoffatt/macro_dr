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

inline auto build_likelihood_function(const ModelPtr& model0,
                                      bool adaptive_approximation, bool recursive_approximation,
                                      int averaging_approximation, bool variance_approximation,
                                      bool taylor_variance_correction_approximation,
                                      bool micro_approximation,
                                      bool taylor_qdt_approximation = false) {
    // Map deprecated bool taylor_qdt_approximation to the new int qdt_method:
    //   false → 0 (eig), true → 2 (schur). Phase 9 wired the previous "true"
    //   semantics through calc_Qdt_schur, so this preserves user behavior.
    //   New code should call build_likelihood_function_with_method() directly
    //   with int qdt_method ∈ {0, 1, 2} (overload below).
    int qdt_method_int = taylor_qdt_approximation ? 2 : 0;
    auto nsub = Simulation_n_sub_dt(100);
    const interface::IModel<var::Parameters_values>& model_ref = *model0;

    return merge_Maybe_variant(
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
                var::constexpr_Var_domain<bool, uses_micro_aproximation, false>,
                decltype(model_ref),
                var::constexpr_Var_domain<int, uses_qdt_method, 0>>(
                model_ref, nsub,
                uses_adaptive_aproximation_value(adaptive_approximation),
                uses_recursive_aproximation_value(recursive_approximation),
                uses_averaging_aproximation_value(averaging_approximation),
                uses_variance_aproximation_value(variance_approximation),
                uses_taylor_variance_correction_aproximation_value(
                    taylor_variance_correction_approximation),
                uses_micro_aproximation_value(micro_approximation),
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
                var::constexpr_Var_domain<bool, uses_micro_aproximation, false>,
                decltype(model_ref),
                var::constexpr_Var_domain<int, uses_qdt_method, 0>>(
                model_ref, nsub,
                uses_adaptive_aproximation_value(adaptive_approximation),
                uses_recursive_aproximation_value(recursive_approximation),
                uses_averaging_aproximation_value(averaging_approximation),
                uses_variance_aproximation_value(variance_approximation),
                uses_taylor_variance_correction_aproximation_value(
                    taylor_variance_correction_approximation),
                uses_micro_aproximation_value(micro_approximation),
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
            var::constexpr_Var_domain<bool, uses_micro_aproximation, true>,
            decltype(model_ref),
            var::constexpr_Var_domain<int, uses_qdt_method, 0, 1, 2>>(
            model_ref, nsub,
            uses_adaptive_aproximation_value(adaptive_approximation),
            uses_recursive_aproximation_value(recursive_approximation),
            uses_averaging_aproximation_value(averaging_approximation),
            uses_variance_aproximation_value(variance_approximation),
            uses_taylor_variance_correction_aproximation_value(
                taylor_variance_correction_approximation),
            uses_micro_aproximation_value(micro_approximation),
            uses_qdt_method_value(qdt_method_int))
            .get_variant());
}

using likelihood_algorithm_type = var::untransformed_type_t<decltype(build_likelihood_function(
    std::declval<const ModelPtr&>(), false, false, 2, true, false, false, false))>;

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

template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path);


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
        Probit_statistics<Likelihood_Distortion_Corrected_Covariance>,
        Probit_statistics<log_Det<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Correlation_Of<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Effective_Rank<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Null_Space_Projector<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Likelihood_Distortion_Corrected_Covariance>>,
        Probit_statistics<Likelihood_Distortion_Induced_Bias>>;

// Per-sample derived diagnostics: Likelihood_Sample_Distortion and
// Likelihood_Distortion_Induced_Bias evaluated at each sample. Included by presets that
// ask for trace-local detail (series_var, series_cov, series_kernel_full).
using Per_sample_derived_diagnostics =
    Probit_statistics<macrodr::Evolution_of<var::Vector_Space<
        Likelihood_Sample_Distortion, Likelihood_Distortion_Induced_Bias>>>;

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
// MLE per-group-of-replicates output (Path A: minimal State).
//
// Composite output of calc_MLE_per_group_of_replicates at one group_size n.
// Only recording-level diagnostics (no per-step Evolution_of) — the state
// stored per group is dMacro_State_Hessian_minimal which lacks the per-step
// data feeding the Moment_statistics<Sum<logL>>, Moment_statistics<dlogL>,
// etc. slots in Analisis_derivative_diagnostic_base. A Path B variant for
// group_size=1 with the full state will be added later, separate type.
//
// Templated on State so the saved probit samples carry the State the command
// was instantiated with — for Path A this is dMacro_State_Hessian_minimal_param.
template <class State>
using MLE_Group_Analysis = var::Vector_Space<
    Group_Size,

    // MLE-specific aggregates over per-group θ̂:
    Probit_statistics<Moment_statistics<Model_Parameters_Hat, true>>,

    // Recording-level Fisher (per group, aggregated to F_b_n at this n):
    Probit_statistics<Likelihood_Numerical_Fisher_Information>,
    Probit_statistics<Likelihood_Fisher_Covariance>,
    Probit_statistics<log_Det<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Eigenvalue_Spectrum<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Correlation_Of<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Spectrum_Condition_Number<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Effective_Rank<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Min_Eigenvalue<Likelihood_Fisher_Covariance>>,

    // Sandwich-corrected covariance (DCC, F⁻¹·J·F⁻¹) — needs the score
    // covariance J across groups; computed from per-group total scores:
    Probit_statistics<Likelihood_Distortion_Corrected_Covariance>,
    Probit_statistics<log_Det<Likelihood_Distortion_Corrected_Covariance>>,
    Probit_statistics<Eigenvalue_Spectrum<Likelihood_Distortion_Corrected_Covariance>>,
    Probit_statistics<Spectrum_Condition_Number<Likelihood_Distortion_Corrected_Covariance>>,

    // Wald T² tests against two metrics:
    Probit_statistics<Wald_T2<Likelihood_Fisher_Covariance>>,
    Probit_statistics<Wald_T2<Likelihood_Distortion_Corrected_Covariance>>,

    // Empirical-covariance vs F⁻¹ distortion + derived (only valid if
    // N_groups ≥ p+1, else NaN-filled by the analysis function):
    Probit_statistics<Empirical_Covariance_Distortion>,
    Probit_statistics<Affine_Invariant_Distance<Empirical_Covariance_Distortion>>,
    Probit_statistics<log_Det<Empirical_Covariance_Distortion>>,
    Probit_statistics<Eigenvalue_Spectrum<Empirical_Covariance_Distortion>>,
    Probit_statistics<Spectrum_Condition_Number<Empirical_Covariance_Distortion>>,

    // K representative groups per ranking variable at each probit height —
    // full per-group State preserved for downstream inspection:
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
// State template parameter (used for the saved probit samples only; the GN
// inner loop always uses dMacro_State_Hessian_minimal for speed):
//   - dMacro_State_Hessian_minimal_param : Path A, minimal memory footprint
//   - dMacro_State_Ev_gradient_all_param : Path B, full Evolution_of preserved
template <class State>
auto calc_MLE_per_group_of_replicates(
    const likelihood_algorithm_type& likelihood_algorithm,
    const var::Parameters_transformed& theta_warmstart,
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
    -> Maybe_error<MLE_Group_Analysis<State>>;

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


}  // namespace macrodr::cmd
