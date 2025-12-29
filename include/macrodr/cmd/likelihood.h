#pragma once

#include <derivative_fwd.h>
#include <distributions.h>
#include <macrodr/interface/IModel.h>

#include <type_traits>

#include "patch_model.h"
#include "qmodel.h"

namespace macrodr::cmd {

inline auto build_likelihood_function(const ModelPtr& model0,
                                      bool adaptive_approximation, bool recursive_approximation,
                                      int averaging_approximation, bool variance_approximation,
                                      bool taylor_variance_correction_approximation) {
    auto nsub = Simulation_n_sub_dt(100);
    const interface::IModel<var::Parameters_values>& model_ref = *model0;

    return Likelihood_Model_regular<
               var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
               var::constexpr_Var_domain<bool, uses_recursive_aproximation, false,true>,
               var::constexpr_Var_domain<int, uses_averaging_aproximation,0,1, 2>,
               var::constexpr_Var_domain<bool, uses_variance_aproximation, false,true>,
               var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
               decltype(model_ref)>(model_ref, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
        .get_variant();
}

using likelihood_algorithm_type = var::untransformed_type_t<decltype(build_likelihood_function(
    std::declval<const ModelPtr&>(), false, false, 2, true, false))>;

auto calculate_mlikelihood(const likelihood_algorithm_type& likelihood_algorithm,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r) -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

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

inline auto calculate_n_simulation_mdlikelihood_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<dMacro_State_Ev_gradient_all>> {
    
    std::vector<dMacro_State_Ev_gradient_all> results;
    results.reserve(simulation.size());
    for (const  auto& sim : simulation){
        auto res = calculate_mdlikelihood_predictions(likelihood_algorithm, par, e,
                                              get<Recording>(sim()));
        if (!res){
            return res.error();
        }
        results.push_back(std::move(res.value()));        
    }
    return results;
}




auto calculate_likelihood(const interface::IModel<var::Parameters_values>& model0,
                          const var::Parameters_transformed& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

inline auto calculate_simulation_likelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    return calculate_likelihood(model0, par, e, get<Recording>(simulation()),
                                adaptive_approximation, recursive_approximation,
                                averaging_approximation, variance_approximation,
                                taylor_variance_correction_approximation);
}

auto calculate_dlikelihood(const interface::IModel<var::Parameters_values>& model0,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal>;

inline auto calculate_simulation_dlikelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal> {
    return calculate_dlikelihood(model0, par, e, get<Recording>(simulation()),
                                 adaptive_approximation, recursive_approximation,
                                 averaging_approximation, variance_approximation,
                                 taylor_variance_correction_approximation);
}
auto calculate_diff_likelihood(const interface::IModel<var::Parameters_values>& model0,
                               const var::Parameters_transformed& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian>;

inline auto calculate_simulation_diff_likelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    return calculate_diff_likelihood(model0, par, e, get<Recording>(simulation()),
                                     adaptive_approximation, recursive_approximation,
                                     averaging_approximation, variance_approximation,
                                     taylor_variance_correction_approximation, delta_param);
}

auto calculate_likelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions>;

inline auto calculate_simulation_likelihood_predictions(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions> {
    return calculate_likelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                            adaptive_approximation, recursive_approximation,
                                            averaging_approximation, variance_approximation,
                                            taylor_variance_correction_approximation);
}

auto calculate_likelihood_diagnostics(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic>;

inline auto calculate_simulation_likelihood_diagnostics(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_likelihood_diagnostics(model0, par, e, get<Recording>(simulation()),
                                            adaptive_approximation, recursive_approximation,
                                            averaging_approximation, variance_approximation,
                                            taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}

inline auto calculate_simulation_sub_dlikelihood_predictions(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_transformed& par,
    const Experiment& e,
    const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& simulation, bool adaptive_approximation,
    bool recursive_approximation, int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions_model(model0, par, e, get<Recording>(simulation()),
                                                   adaptive_approximation, recursive_approximation,
                                                   averaging_approximation, variance_approximation,
                                                   taylor_variance_correction_approximation);
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);


template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
                                   std::vector<TMacro_State<vVars...>> const& liks, std::string path);

template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path);

}  // namespace macrodr::cmd
