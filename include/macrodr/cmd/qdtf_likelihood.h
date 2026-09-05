#pragma once

// cmd surface for the Bessel member (qdtf): direct likelihood evaluation.
// n_poles = 0 → box configuration (the regression anchor against the
// existing av=2 member); n_poles = 4 or 8 → Bessel filter with the given
// −3 dB cutoff in Hz. The filter parameters are rig metadata: fixed inputs,
// never fitted. See legacy/qdtf_member.h for the member itself.

#include <macrodr/cmd/likelihood.h>

#include "qdtf_member.h"

namespace macrodr::cmd {

inline auto calculate_qdtf_likelihood(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, std::size_t n_poles, double cutoff_hz)
    -> Maybe_error<logLs> {
    const interface::IModel<var::Parameters_values>& model_ref = *model0;
    auto par_values = par.to_value();
    return qdtf_log_Likelihood(model_ref, par_values, r, e, int(n_poles), cutoff_hz);
}

inline auto calculate_simulation_qdtf_likelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& sim, std::size_t n_poles, double cutoff_hz)
    -> Maybe_error<logLs> {
    return calculate_qdtf_likelihood(model0, par, e, get<Recording>(sim()), n_poles, cutoff_hz);
}

}  // namespace macrodr::cmd
