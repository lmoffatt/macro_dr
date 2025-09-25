#include <macrodr/cli/command_manager.h>
#include <macrodr/cmd/cli_meta.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/patch_model.h>

// Legacy registry builders used until migrated
#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
#include "lapack_headers.h"

namespace macrodr::cli {

dsl::Compiler make_compiler_new() {
    dsl::Compiler cm;
    // Meta commands (help/version) available everywhere
    cm.merge(macrodr::cmd::make_cli_meta_compiler());
    cm.merge(macrodr::cmd::make_utilities_compiler());
    cm.merge(macrodr::cmd::make_io_compiler());
    cm.merge(macrodr::cmd::make_experiment_compiler());
    cm.merge(macrodr::cmd::make_model_compiler());
    cm.merge(macrodr::cmd::make_simulation_compiler());
    cm.merge(macrodr::cmd::make_likelihood_compiler());
    cm.merge(macrodr::cmd::make_dts_compiler());
    cm.push_function("load_model",
                     dsl::to_typed_function<std::string>(&cmd::load_model, "model_name"));
    cm.push_function("load_experiment", dsl::to_typed_function<std::string, double, double>(
                                            &macrodr::cmd::load_experiment, "filename",
                                            "frequency_of_sampling", "initial_ATP"));

    // Expose low-level model/qmodel helpers
    using ModelPtr = std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>;
    cm.push_function("patch_model",
                     dsl::to_typed_function<ModelPtr, std::pair<std::string, std::string>>(
                         &macrodr::cmd::patch_model, "model", "parameter_values"));

    using PatchModel = macrodr::Transfer_Op_to<var::Parameters_values, macrodr::Patch_Model>;
    cm.push_function("calc_eigen",
                     dsl::to_typed_function<PatchModel, double>(
                         &macrodr::cmd::calc_eigen, "patch_model", "atp"));

    using QxEig = macrodr::Transfer_Op_to<PatchModel, macrodr::Qx_eig>;
    cm.push_function("calc_peq",
                     dsl::to_typed_function<QxEig, PatchModel>(
                         &macrodr::cmd::calc_peq, "qx_eig", "patch_model"));

    return cm;
}

}  // namespace macrodr::cli
