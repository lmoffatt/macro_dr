#pragma once

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

}  // namespace macrodr::cmd
