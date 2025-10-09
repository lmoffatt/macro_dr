#pragma once

#include <distributions.h>
#include <macrodr/interface/IModel.h>

#include "qmodel.h"

namespace macrodr::cmd {
auto calculate_likelihood(const interface::IModel<var::Parameters_values>& model0,
                          const var::Parameters_values& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

inline auto calculate_simulation_likelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_values& par,
    const Experiment& e, const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    return calculate_likelihood(model0, par, e, get<Recording>(simulation()),
                                adaptive_approximation, recursive_approximation,
                                averaging_approximation, variance_approximation,
                                taylor_variance_correction_approximation);
}

auto calculate_dlikelihood(const interface::IModel<var::Parameters_values>& model0,
                           const var::Parameters_values& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation) -> Maybe_error<dlogLs>;

inline auto calculate_simulation_dlikelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_values& par,
    const Experiment& e, const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dlogLs> {
    return calculate_dlikelihood(model0, par, e, get<Recording>(simulation()),
                                 adaptive_approximation, recursive_approximation,
                                 averaging_approximation, variance_approximation,
                                 taylor_variance_correction_approximation);
}
auto calculate_diff_likelihood(const interface::IModel<var::Parameters_values>& model0,
                               const var::Parameters_values& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<dlogLs>;

inline auto calculate_simulation_diff_likelihood(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_values& par,
    const Experiment& e, const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<dlogLs> {
    return calculate_diff_likelihood(model0, par, e, get<Recording>(simulation()),
                                     adaptive_approximation, recursive_approximation,
                                     averaging_approximation, variance_approximation,
                                     taylor_variance_correction_approximation, delta_param);
}

auto calculate_likelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_values& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Patch_State_Evolution>;

inline auto calculate_simulation_likelihood_predictions(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_values& par,
    const Experiment& e, const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Patch_State_Evolution> {
    return calculate_likelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                            adaptive_approximation, recursive_approximation,
                                            averaging_approximation, variance_approximation,
                                            taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                       const var::Parameters_values& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>>;

inline auto calculate_simulation_dlikelihood_predictions(
    const interface::IModel<var::Parameters_values>& model0, const var::Parameters_values& par,
    const Experiment& e, const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions_model(const std::string& model0,
                                             const var::Parameters_values& par, const Experiment& e,
                                             const Recording& r, bool adaptive_approximation,
                                             bool recursive_approximation,
                                             int averaging_approximation,
                                             bool variance_approximation,
                                             bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>>;

inline auto calculate_simulation_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_values& par, const Experiment& e,
    const Simulated_Recording<includes_N_state_evolution<false>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
    return calculate_dlikelihood_predictions_model(model0, par, e, get<Recording>(simulation()),
                                                   adaptive_approximation, recursive_approximation,
                                                   averaging_approximation, variance_approximation,
                                                   taylor_variance_correction_approximation);
}

}  // namespace macrodr::cmd
