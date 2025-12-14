#pragma once

#include <distributions.h>
#include <macrodr/interface/IModel.h>

#include "qmodel.h"

namespace macrodr::cmd {
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
                           bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Hessian_minimal>;

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
    const Experiment& e,const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_likelihood_diagnostics(model0, par, e,get<Recording>(simulation()),
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
    const Experiment& e, const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}



auto calculate_dlikelihood_predictions_model(const std::string& model0,
                                             const var::Parameters_transformed& par, const Experiment& e,
                                             const Recording& r, bool adaptive_approximation,
                                             bool recursive_approximation,
                                             int averaging_approximation,
                                             bool variance_approximation,
                                             bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions_model(model0, par, e, get<Recording>(simulation()),
                                                   adaptive_approximation, recursive_approximation,
                                                   averaging_approximation, variance_approximation,
                                                   taylor_variance_correction_approximation);
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);

template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path);

}  // namespace macrodr::cmd
