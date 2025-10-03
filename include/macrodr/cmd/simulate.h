#ifndef SIMULATE_H
#define SIMULATE_H

#include <cstddef>
#include <string>

#include "CLI_macro_dr_base.h"
//#include "cuevi.h"
#include "function_memoization.h"
#include "maybe_error.h"
#include "models_MoffattHume_allosteric.h"
#include "parallel_tempering.h"
#include "parameters.h"
//#include "parameters_derivative.h"
#include <macrodr/interface/IModel.h>

#include "qmodel.h"

namespace macrodr::cmd {

Maybe_error<Simulated_Recording<includes_N_state_evolution<false>>> run_simulations(
    const interface::IModel<var::Parameters_values>& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<std::string> runsimulation(std::string filename_prefix, recording_type recording_file,
                                       experiment_type experiment, std::size_t myseed,
                                       const std::string& modelName,
                                       parameters_value_type parameter_files,
                                       simulation_algo_type sim_algo_type);

}  // namespace macrodr::cmd

#endif  // SIMULATE_H
