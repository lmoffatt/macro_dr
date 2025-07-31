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
#include "qmodel.h"
namespace macrodr {
namespace cmd {

std::string run_simulation(std::string filename_prefix, recording_type recording_file,
                           experiment_type experiment, std::size_t myseed, std::string modelName,
                           parameters_value_type parameter_files,
                           simulation_algo_type sim_algo_type);
}
}  // namespace macrodr

#endif  // SIMULATE_H
