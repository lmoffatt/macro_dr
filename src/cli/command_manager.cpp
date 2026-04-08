#include <experiment.h>
#include <macrodr/cli/command_manager.h>
#include <macrodr/cmd/cli_meta.h>
#include <macrodr/cmd/indexed_construction.h>
#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/patch_model.h>
#include <macrodr/cmd/simulate.h>
#include <maybe_error.h>
#include <parameters.h>
#include <qmodel.h>
#include <variables.h>

#include <cstddef>
#include <set>

// Legacy registry builders used until migrated
#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
#include "function_builder.h"
#include "lapack_headers.h"

namespace macrodr::cli {
using std::size_t;
using ModelPtr = macrodr::cmd::ModelPtr;

namespace {

template<bool include_evolution=true>
auto calculate_likelihood_derivative_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed) {
    return cmd::calculate_Likelihood_derivative_diagnostics<include_evolution>(
        dlikelihood_predictions, n_boostrap_samples, probits, seed);
}

}

inline macrodr::dsl::Compiler<macrodr::dsl::Lexer> make_simulations_compiler() {
    dsl::Compiler<dsl::Lexer> cm;
    cm.push_function("get_num_parameters",
                     dsl::to_typed_function<std::string>(&get_num_parameters, "model"));

    cm.push_function("simulate",
                     dsl::to_typed_function<std::string, cmd::recording_type, cmd::experiment_type,
                                            std::size_t, const std::string&,
                                            cmd::parameters_value_type, cmd::simulation_algo_type>(
                         &cmd::runsimulation, "filename_prefix", "recording", "experiment",
                         "init_seed", "modelName", "parameter_values", "simulation_algorithm"));

    {
        using SimFromValues = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromValues>(&cmd::run_simulations), "model", "parameter_values",
                "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimUniformFromValues = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_values&, const Experiment&, const Recording&,
            std::string, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_values&,
                                   const Experiment&, const Recording&, std::string,std::size_t, std::size_t>
                                   (
                static_cast<SimUniformFromValues>(&cmd::run_simulations), "model",
                "parameter_values", "experiment", "observations", "simulation_algorithm","number_of_substeps",
                "seed"));
    }

    {
        using SimFromTransformed = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromTransformed>(&cmd::run_simulations), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimUniformFromTransformed = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&, const Experiment&,
            const Recording&, std::string, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                                   const Experiment&, const Recording&, std::string,std::size_t, std::size_t>(
                static_cast<SimUniformFromTransformed>(&cmd::run_simulations), "model",
                "parameter_transformed", "experiment", "observations",
                "simulation_algorithm", "number_of_substeps", "seed"));
    }
    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&,
            const var::Parameters_transformed&, std::size_t,const Experiment&,
            const Recording&, std::size_t, std::size_t>(
                static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                    const ModelPtr&, const var::Parameters_transformed&, std::size_t,
                    const Experiment&, const Recording&, std::size_t, std::size_t)>(
                    &cmd::run_n_simulations),
                "model",
                "parameter_transformed", "n_simulations", "experiment", "observations", "number_of_substeps", "seed"));

    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&, const var::Parameters_transformed&, std::size_t,
            const Experiment&, const Recording&, std::string, std::size_t, std::size_t>(
            static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                const ModelPtr&, const var::Parameters_transformed&, std::size_t,
                const Experiment&, const Recording&, std::string, std::size_t, std::size_t)>(
                &cmd::run_n_simulations),
            "model", "parameter_transformed", "n_simulations", "experiment", "observations",
            "simulation_algorithm", "number_of_substeps", "seed"));
    
    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&,
            const var::Parameters_values&, std::size_t,const Experiment&,
            const Recording&, std::size_t, std::size_t>(
                static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                    const ModelPtr&, const var::Parameters_values&, std::size_t,
                    const Experiment&, const Recording&, std::size_t, std::size_t)>(
                    &cmd::run_n_simulations),
                "model",
                "parameter_values", "n_simulations", "experiment", "observations", "number_of_substeps", "seed"));

    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&, const var::Parameters_values&, std::size_t,
            const Experiment&, const Recording&, std::string, std::size_t, std::size_t>(
            static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                const ModelPtr&, const var::Parameters_values&, std::size_t,
                const Experiment&, const Recording&, std::string, std::size_t, std::size_t)>(
                &cmd::run_n_simulations),
            "model", "parameter_values", "n_simulations", "experiment", "observations",
            "simulation_algorithm", "number_of_substeps", "seed"));
    

    {
        using SimSubFromValues = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromValues>(&cmd::run_simulations_with_sub_intervals), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimSubUniformFromValues = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_values&, const Experiment&, const Recording&,
            std::string, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_values&,
                                   const Experiment&, const Recording&, std::string,
                                   std::size_t>(
                static_cast<SimSubUniformFromValues>(&cmd::run_simulations_with_sub_intervals),
                "model", "parameter_values", "experiment", "observations",
                "simulation_algorithm", "seed"));
    }

    {
        using SimSubFromTransformed = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromTransformed>(&cmd::run_simulations_with_sub_intervals),
                "model", "parameters", "experiment", "observations", "number_of_substeps",
                "seed"));
    }

    {
        using SimSubUniformFromTransformed = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&, const Experiment&,
            const Recording&, std::string, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                                   const Experiment&, const Recording&, std::string,
                                   std::size_t>(
                static_cast<SimSubUniformFromTransformed>(
                    &cmd::run_simulations_with_sub_intervals),
                "model", "parameter_transformed", "experiment", "observations",
                "simulation_algorithm", "seed"));
    }

    cm.push_function(
        "remove_intervals",
        dsl::to_typed_function<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&>(
            &cmd::remove_intervals, "simulation_with_sub_intervals"));  
  
    return cm;
}

dsl::Compiler<dsl::Lexer> make_compiler_new() {
    dsl::Compiler<dsl::Lexer> cm;
    // Meta commands (help/version) available everywhere
    cm.merge(macrodr::cmd::make_cli_meta_compiler());
    cm.merge(macrodr::cmd::make_utilities_compiler());
    cm.merge(macrodr::cmd::make_io_compiler());
    cm.merge(macrodr::cmd::make_experiment_compiler());
    cm.merge(macrodr::cmd::make_model_compiler());
    cm.merge(make_simulations_compiler());
    cm.merge(macrodr::cmd::make_likelihood_compiler());
    cm.merge(macrodr::cmd::make_dts_compiler());
    cm.push_function("load_model",
                     dsl::to_typed_function<std::string>(&cmd::load_model, "model_name"));
    cm.push_function("load_experiment", dsl::to_typed_function<std::string, double, double>(
                                            &macrodr::cmd::load_experiment, "filename",
                                            "frequency_of_sampling", "initial_agonist"));

    cm.push_function(
        "create_experiment",
        dsl::to_typed_function<std::vector<std::tuple<std::size_t, std::size_t, double>>, double,
                               double, double>(&macrodr::cmd::create_experiment,
                                               "experiment_structure", "frequency_of_sampling",
                                               "initial_agonist", "initial_time"));

    cm.push_function("load_observations", dsl::to_typed_function<std::string>(
                                              &macrodr::cmd::load_recording, "filename"));

    
    cm.push_function("set_observations", dsl::to_typed_function<std::vector<double>>(
                                              &macrodr::cmd::define_recording, "values_set"));

    cm.push_function("axis",
                     dsl::to_typed_function<std::string, std::vector<std::string>>(
                         &cmd::axis, "name", "labels"));
    cm.push_function("indexed_bool_by",
                     dsl::to_typed_function<var::Axis, std::vector<bool>>(
                         &cmd::indexed_by<bool>, "axis", "values"));
    cm.push_function("indexed_int_by",
                     dsl::to_typed_function<var::Axis, std::vector<int>>(
                         &cmd::indexed_by<int>, "axis", "values"));
    cm.push_function("indexed_size_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                         &cmd::indexed_by<std::size_t>, "axis", "values"));
    cm.push_function("indexed_double_by",
                     dsl::to_typed_function<var::Axis, std::vector<double>>(
                         &cmd::indexed_by<double>, "axis", "values"));
    cm.push_function("indexed_string_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::string>>(
                         &cmd::indexed_by<std::string>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                         &cmd::indexed_by<std::size_t>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<int>>(
                         &cmd::indexed_by<int>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<double>>(
                         &cmd::indexed_by<double>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::string>>(
                         &cmd::indexed_by<std::string>, "axis", "values"));

    using SimulationBatch = std::vector<Simulated_recording>;
    using IndexedSimulation = var::Indexed<Simulated_recording>;
    using IndexedSimulationBatch = var::Indexed<SimulationBatch>;
    using SubSimulation =
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>;
    using IndexedSubSimulation = var::Indexed<SubSimulation>;
    using GradientBatch = std::vector<dMacro_State_Ev_gradient_all>;
    using IndexedPredictions = var::Indexed<Macro_State_Ev_predictions>;
    using IndexedDiagnostics = var::Indexed<Macro_State_Ev_diagnostic>;
    using IndexedGradients = var::Indexed<dMacro_State_Ev_gradient_all>;
    using IndexedGradientBatch = var::Indexed<GradientBatch>;
    using IndexedDerivativeDiagnosticsEvolution =
        var::Indexed<macrodr::cmd::Analisis_derivative_diagnostic<true>>;
     using IndexedDerivativeDiagnostics =
        var::Indexed<macrodr::cmd::Analisis_derivative_diagnostic<false>>;

    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,Simulated_recording const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","simulation", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "path"));
    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,std::vector<Simulated_recording> const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","simulations", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulations", "path"));
    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,Recording const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","observations", "path"));

    cm.push_function("write_csv", 
        dsl::to_typed_return_function<
        Maybe_error<std::string>,
        Experiment const&,
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&, 
        std::size_t,
        std::string >(
            &macrodr::cmd::write_csv,
            "experiment",
            "simulation", 
            "number_of_substates",
            "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&, std::size_t, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "number_of_substates",
            "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_predictions const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      Macro_State_Ev_predictions const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedPredictions const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedPredictions const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_diagnostic const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      Macro_State_Ev_diagnostic const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedDiagnostics const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedDiagnostics const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, GradientBatch const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, GradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, IndexedGradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, IndexedGradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));


    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, Experiment const&,
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
            dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SubSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_evolution const& ,
                                          std::string >(
            &macrodr::cmd::write_csv, "analysis",  "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_base const& ,
                                          std::string >(
            &macrodr::cmd::write_csv, "analysis",  "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, IndexedDerivativeDiagnostics const&,
                                      std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, IndexedDerivativeDiagnosticsEvolution const&,
                                      std::string>(&macrodr::cmd::write_csv, "analysis", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<std::size_t>&,
                                      std::string>(&macrodr::cmd::write_csv<std::size_t>,
                                                   "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<int>&,
                                      std::string>(&macrodr::cmd::write_csv<int>, "analysis",
                                                   "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<double>&,
                                      std::string>(&macrodr::cmd::write_csv<double>, "analysis",
                                                   "path"));

            

                                              // Expose low-level model/qmodel helpers
    // Load parameters file (filename, sep) using model’s schema
    cm.push_function("load_parameters",
                     dsl::to_typed_function<ModelPtr, std::string>(&macrodr::cmd::load_parameters,
                                                                   "model", "parameter_file"));

    cm.push_function("create_parameters",
                     dsl::to_typed_function<ModelPtr,  std::vector<std::tuple<std::string, std::string, double>>>
                     (&macrodr::cmd::create_parameters,"model", "parameters_info"));

    
                                                                   // load_parameter_values returns transformations (safer to persist)
    cm.push_function("get_standard_parameter_values",
                     dsl::to_typed_function<var::Parameters_Transformations const&>(
                         &macrodr::cmd::get_standard_parameter_values, "parameters"));

    cm.push_function("get_standard_parameter_transformed_values",
                     dsl::to_typed_function<var::Parameters_Transformations const&>(
                         &macrodr::cmd::get_standard_parameter_transformed_values, "parameters"));

    
    cm.push_function(
        "patch_model",
        dsl::to_typed_return_function<Maybe_error<macrodr::cmd::PatchModel>, ModelPtr const&,
                                      std::pair<std::string, std::string> const&>(
            &macrodr::cmd::patch_model, "model", "parameter_values"));
    // Overload: patch from in-memory values
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromValues = Return (*)(const ModelPtr&, const var::Parameters_values&);
        cm.push_function(
            "patch_model",
            dsl::to_typed_function<ModelPtr, var::Parameters_values>(
                static_cast<FnPatchFromValues>(&macrodr::cmd::patch_model), "model", "values"));
    }
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromTr = Return (*)(const ModelPtr&, const var::Parameters_Transformations&);
        cm.push_function("patch_model",
                         dsl::to_typed_function<ModelPtr, var::Parameters_Transformations>(
                             static_cast<FnPatchFromTr>(&macrodr::cmd::patch_model), "model",
                             "parameter_values"));
    }
    // Back-compat: accept 'parameter_values' as the argument name for in-memory values as well
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromValues = Return (*)(const ModelPtr&, const var::Parameters_values&);
        cm.push_function("patch_model",
                         dsl::to_typed_function<ModelPtr, var::Parameters_values>(
                             static_cast<FnPatchFromValues>(&macrodr::cmd::patch_model), "model",
                             "parameter_values"));
    }

    using PatchModel = macrodr::Transfer_Op_to<var::Parameters_values, macrodr::Patch_Model>;
    cm.push_function("calc_eigen", dsl::to_typed_function<PatchModel const&, double>(
                                       &macrodr::cmd::calc_eigen, "patch_model", "agonist"));

    using QxEig = macrodr::Transfer_Op_to<PatchModel, macrodr::Eigs>;
    cm.push_function("calc_peq", dsl::to_typed_function<QxEig const&, PatchModel const&>(
                                     &macrodr::cmd::calc_peq, "qx_eig", "patch_model"));

    // Patch state initial condition helpers
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromPM = RetPS (*)(const macrodr::cmd::PatchModel&, double);
        cm.push_function("path_state", dsl::to_typed_function<PatchModel const&, double>(
                                           static_cast<FnPSfromPM>(&macrodr::cmd::path_state),
                                           "patch_model", "initial_agonist"));
    }
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromVals = RetPS (*)(const ModelPtr&, const var::Parameters_values&, double);
        cm.push_function("path_state",
                         dsl::to_typed_function<ModelPtr, var::Parameters_values const&, double>(
                             static_cast<FnPSfromVals>(&macrodr::cmd::path_state), "model",
                             "values", "initial_agonist"));
    }

    cm.push_function("build_likelihood_function",
                     dsl::to_typed_function<
                         ModelPtr const&, bool, bool, int, bool, bool>(
                         &macrodr::cmd::build_likelihood_function, "model",
                                "adaptive_approximation", "recursive_approximation",
                                "averaging_approximation", "variance_approximation",
                                "taylor_variance_correction")
                            );

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood, "likelihood_algorithm", "parameters", "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mdlikelihood, "likelihood_algorithm", "parameters", "experiment",
            "data"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&, const var::Parameters_transformed&,
                               const Experiment&, const Recording&, double>(
            &cmd::calculate_mdiff_likelihood, "likelihood_algorithm", "parameters", "experiment",
            "data", "delta_param"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data"));

        cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Simulated_Recording<var::please_include<>>&>(
            &cmd::calculate_simulation_mlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data")); 


    cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood_diagnostics, "likelihood_algorithm", "parameters",
            "experiment", "data"));
    cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Simulated_Recording<var::please_include<>>&&>(
            &cmd::calculate_simulated_mlikelihood_diagnostics, "likelihood_algorithm", "parameters",
            "experiment", "simulation"));


    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mdlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&>(
            &cmd::calculate_simulation_mdlikelihood_predictions, "likelihood_algorithm",
            "parameters", "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<> >>&>(
            &cmd::calculate_n_simulation_mdlikelihood_predictions, "likelihood_algorithm",
            "parameters", "experiment", "data_series"));


            cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, double>(
            &cmd::calculate_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "delta_param"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, double>(
            &cmd::calculate_simulation_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "delta_param"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));
    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_dlikelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood_diagnostics, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, 
                               const  Simulated_Recording<var::please_include<>>&&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood_diagnostics, "model", "parameters", "experiment", "simulation",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));


    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_dlikelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_sub_dlikelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

   
            cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const std::string&, const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_dlikelihood_predictions_model, "model", "parameters",
            "experiment", "data", "adaptive_approximation", "recursive_approximation",
            "averaging_approximation", "variance_approximation", "taylor_variance_correction"));

      
            cm.push_function(
        "likelihood_derivative_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>& , 
    std::size_t, std::set<double>,  std::size_t>(
            &calculate_likelihood_derivative_diagnostics_dsl<false>, "dlikelihood_predictions", "n_boostrap_samples",
            "probits", "seed"));
      cm.push_function(
        "likelihood_derivative_evolution_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>& , 
    std::size_t, std::set<double>,  std::size_t>(
            &calculate_likelihood_derivative_diagnostics_dsl<true>, "dlikelihood_predictions", "n_boostrap_samples",
            "probits", "seed"));


    return cm;
}

}  // namespace macrodr::cli
