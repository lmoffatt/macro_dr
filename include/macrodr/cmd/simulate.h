#ifndef SIMULATE_H
#define SIMULATE_H

#include <cstddef>
#include <string>
#include <vector>

#include "CLI_macro_dr_base.h"
#include "maybe_error.h"
#include "parameters.h"
#include "patch_model.h"
//#include "parameters_derivative.h"

#include "qmodel.h"

namespace macrodr::cmd {

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::string simulation_algorithm, std::size_t number_of_substeps,
    std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::string simulation_algorithm, std::size_t number_of_substeps,
    std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> run_simulations_with_sub_intervals(
    const ModelPtr& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> run_simulations_with_sub_intervals(
    const ModelPtr& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
run_simulations_with_sub_intervals(const ModelPtr& model, const var::Parameters_values& par,
                                   const Experiment& e, const Recording& r,
                                   std::string simulation_algorithm,
                                   std::size_t myseed);

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
run_simulations_with_sub_intervals(const ModelPtr& model,
                                   const var::Parameters_transformed& par, const Experiment& e,
                                   const Recording& r, std::string simulation_algorithm,
                                   std::size_t myseed);

inline Simulated_Recording<var::please_include<>> remove_intervals(
    const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation){

    return Simulated_Recording<var::please_include<> >{{get<SeedNumber>(simulation()), get<Recording>(simulation())}};
}

// Load simulations from a CSV file written by write_csv(experiment, simulations,
// path). Returns one Simulated_Recording per requested simulation_index, with
// patch_current values populated from the "value" column. If
// replica_indices is empty, loads ALL simulations found in the file; otherwise
// only loads those listed. SeedNumber is set to 0 in the loaded records (the
// original seed is not persisted in the CSV; downstream uses only need the
// recorded patch_current values).
// Index-aware sibling of simulate(): reads a simulation CSV and reconstructs
// the SAME var::Indexed<vector<Simulated_Recording>> type simulate() produces
// via DSL lifting, carrying the same axes (reconstructed from the CSV's axis
// columns — any header column outside the fixed schema). Keeps the axis chain
// alive so downstream stays Indexed and axis columns propagate.
//
// Args BY VALUE (not const-ref) so the DSL accepts literal arguments inline
// (load_simulations(filename="path", replica_indices=[0,44])). Errors if any
// reconstructed axis coordinate is missing data.
Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>>
load_simulations(std::string filename,
                 std::vector<std::size_t> replica_indices);

// Overload that loads all simulations (no index filter).
Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>>
load_simulations(std::string filename);


Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_values& par,std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::size_t number_of_substeps, std::size_t myseed);

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_values& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::string simulation_algorithm,std::size_t number_of_substeps, 
    std::size_t myseed);

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::string simulation_algorithm,std::size_t number_of_substeps, 
    std::size_t myseed);

Maybe_error<std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>> run_n_simulations_with_sub_intervals(
    const ModelPtr& model, const var::Parameters_values& par,  std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

Maybe_error<std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>> run_simulations_with_sub_intervals(
    const ModelPtr& model, const var::Parameters_transformed& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed);

inline std::vector<Simulated_Recording<var::please_include<>>> remove_n_intervals(
    const std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>& simulation){

std::vector<Simulated_Recording<var::please_include<> >> result;
result.reserve(simulation.size());
for (const auto& sim : simulation){
    result.push_back(Simulated_Recording<var::please_include<> >{{get<SeedNumber>(sim()), get<Recording>(sim())}});
}
    return result;
}

Maybe_error<std::string> runsimulation(std::string filename_prefix, recording_type recording_file,
                                       experiment_type experiment, std::size_t myseed,
                                       const std::string& modelName,
                                       parameters_value_type parameter_files,
                                       simulation_algo_type sim_algo_type);



Maybe_error<std::string> write_csv(Experiment const& e,
    std::vector<Simulated_Recording<var::please_include<>>> const& simulation, std::string  path);

Maybe_error<std::string> write_csv(
    Experiment const& e,
    var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>> const& simulation,
    std::string path);

Maybe_error<std::string> write_csv(Experiment const& e, 
    std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const& simulation, std::size_t n_sub, std::string  path);

Maybe_error<std::string> write_csv(Experiment const& e,
    Simulated_Recording<var::please_include<> > const& simulation, std::string  path);

Maybe_error<std::string> write_csv(
    Experiment const& e,
    var::Indexed<Simulated_Recording<var::please_include<>>> const& simulation,
    std::string path);

Maybe_error<std::string> write_csv(Experiment const& e, 
    Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const& simulation, std::size_t n_sub, std::string  path);

Maybe_error<std::string> write_csv(
    Experiment const& e,
    var::Indexed<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const&
        simulation,
    std::size_t n_sub, std::string path);

Maybe_error<std::string> write_csv(Experiment const& e, 
    std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const& simulation, std::size_t n_sub, const std::string& path);

Maybe_error<std::string> write_csv(
    Experiment const& e,
    var::Indexed<
        std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>> const&
        simulation,
    std::size_t n_sub, std::string path);


}  // namespace macrodr::cmd

#endif  // SIMULATE_H
