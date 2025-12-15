#include <experiment.h>
#include <macrodr/cli/command_manager.h>
#include <macrodr/cmd/cli_meta.h>
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

// Legacy registry builders used until migrated
#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
#include "function_builder.h"
#include "lapack_headers.h"

namespace macrodr::cli {
using std::size_t;

inline macrodr::dsl::Compiler make_simulations_compiler() {
    dsl::Compiler cm;
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
            interface::IModel<var::Parameters_values> const&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<interface::IModel<var::Parameters_values> const&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromValues>(&cmd::run_simulations), "model", "parameter_values",
                "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimFromTransformed = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            interface::IModel<var::Parameters_values> const&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<interface::IModel<var::Parameters_values> const&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromTransformed>(&cmd::run_simulations), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimSubFromValues = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            interface::IModel<var::Parameters_values> const&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<interface::IModel<var::Parameters_values> const&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromValues>(&cmd::run_simulations_with_sub_intervals), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimSubFromTransformed = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            interface::IModel<var::Parameters_values> const&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<interface::IModel<var::Parameters_values> const&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromTransformed>(&cmd::run_simulations_with_sub_intervals),
                "model", "parameter_values", "experiment", "observations", "number_of_substeps",
                "seed"));
    }

    cm.push_function(
        "remove_intervals",
        dsl::to_typed_function<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&>(
            &cmd::remove_intervals, "simulation_with_sub_intervals"));  
  
    return cm;
}

dsl::Compiler make_compiler_new() {
    dsl::Compiler cm;
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

    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,Simulated_recording const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","simulation", "path"));
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
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_predictions const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_diagnostic const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, Experiment const&,
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
            dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));

            

                                              // Expose low-level model/qmodel helpers
    using ModelPtr = std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>;

    // Load parameters file (filename, sep) using model’s schema
    cm.push_function("load_parameters",
                     dsl::to_typed_function<ModelPtr, std::string>(&macrodr::cmd::load_parameters,
                                                                   "model", "parameter_file"));

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

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, double>(
            &cmd::calculate_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "delta_param"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, double>(
            &cmd::calculate_simulation_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "delta_param"));

    cm.push_function(
        "cal_likelihood_predictions",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));
    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_dlikelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_likelihood_diagnostics, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&, 
                               const  Simulated_Recording<var::please_include<>>&&,
                               bool, bool, int, bool, bool>(
            &cmd::calculate_simulation_likelihood_diagnostics, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));


    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool>(
            &cmd::calculate_simulation_dlikelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const interface::IModel<var::Parameters_values>&,
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

    return cm;
}

}  // namespace macrodr::cli
