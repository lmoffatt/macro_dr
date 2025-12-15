#ifndef SIMULATE_H
#define SIMULATE_H

#include <cstddef>
#include <string>

#include "CLI_macro_dr_base.h"
#include "maybe_error.h"
#include "parameters.h"
//#include "parameters_derivative.h"
#include <macrodr/interface/IModel.h>

#include "qmodel.h"

namespace macrodr::cmd {

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const interface::IModel<var::Parameters_values>& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const interface::IModel<var::Parameters_values>& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> run_simulations_with_sub_intervals(
    const interface::IModel<var::Parameters_values>& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> run_simulations_with_sub_intervals(
    const interface::IModel<var::Parameters_values>& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

inline Simulated_Recording<var::please_include<>> remove_intervals(
    const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation){
    return Simulated_Recording<var::please_include<>>(
        get<Recording>(simulation()));         
    }

Maybe_error<std::string> runsimulation(std::string filename_prefix, recording_type recording_file,
                                       experiment_type experiment, std::size_t myseed,
                                       const std::string& modelName,
                                       parameters_value_type parameter_files,
                                       simulation_algo_type sim_algo_type);

Maybe_error<std::string> write_csv(Experiment const& e,
    Simulated_Recording<var::please_include<>> const& simulation, std::string  path);

Maybe_error<std::string> write_csv(Experiment const& e,
    Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const& simulation, std::size_t n_sub, std::string  path);

}  // namespace macrodr::cmd

#endif  // SIMULATE_H
