#include <macrodr/cli/command_manager.h>
#include <macrodr/cmd/cli_meta.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/patch_model.h>
#include <macrodr/cmd/simulate.h>

// Legacy registry builders used until migrated
#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
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

    cm.push_function("simulate",
                     dsl::to_typed_function<interface::IModel<var::Parameters_values> const&,
                                            const var::Parameters_values&, const Experiment&,
                                            const Simulation_Parameters&, const Recording&, size_t>(
                         &cmd::run_simulations, "model", "parameter_values", "experiment",
                         "Simulation_Parameters", "Recording", "seed"));

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
                                            "frequency_of_sampling", "initial_ATP"));

    // Expose low-level model/qmodel helpers
    using ModelPtr = std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>;

    // Load parameters file (filename, sep) using model’s schema
    cm.push_function("load_parameters",
                     dsl::to_typed_function<ModelPtr, std::pair<std::string, std::string>>(
                         &macrodr::cmd::load_parameters, "model", "parameter_values"));

    // load_parameter_values returns transformations (safer to persist)
    cm.push_function("load_parameter_values",
                     dsl::to_typed_function<ModelPtr, std::pair<std::string, std::string>>(
                         &macrodr::cmd::load_parameter_values, "model", "parameter_values"));

    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromFile =
            Return (*)(const ModelPtr&, const std::pair<std::string, std::string>&);
        cm.push_function("patch_model",
                         dsl::to_typed_function<ModelPtr, std::pair<std::string, std::string>>(
                             static_cast<FnPatchFromFile>(&macrodr::cmd::patch_model), "model",
                             "parameter_values"));
    }
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
    cm.push_function("calc_eigen", dsl::to_typed_function<PatchModel, double>(
                                       &macrodr::cmd::calc_eigen, "patch_model", "atp"));

    using QxEig = macrodr::Transfer_Op_to<PatchModel, macrodr::Qx_eig>;
    cm.push_function("calc_peq", dsl::to_typed_function<QxEig, PatchModel>(
                                     &macrodr::cmd::calc_peq, "qx_eig", "patch_model"));

    // Patch state initial condition helpers
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromPM = RetPS (*)(const macrodr::cmd::PatchModel&, double);
        cm.push_function("path_state", dsl::to_typed_function<PatchModel, double>(
                                           static_cast<FnPSfromPM>(&macrodr::cmd::path_state),
                                           "patch_model", "initial_atp"));
    }
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromVals = RetPS (*)(const ModelPtr&, const var::Parameters_values&, double);
        cm.push_function("path_state",
                         dsl::to_typed_function<ModelPtr, var::Parameters_values, double>(
                             static_cast<FnPSfromVals>(&macrodr::cmd::path_state), "model",
                             "values", "initial_atp"));
    }

    return cm;
}

}  // namespace macrodr::cli
