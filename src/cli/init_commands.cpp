// Aggregates all CMD builders into a single compiler
#include <macrodr/cmd/init_commands.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/dsl/lexer_typed.h>

// Bring inline LAPACK-backed functions into this TU

// Legacy registry builders used until migrated

#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
#include "lapack_headers.h"

namespace macrodr {

dsl::Compiler make_compiler_new() {
    dsl::Compiler cm;
    cm.merge(macrodr::cmd::make_utilities_compiler());
    cm.merge(macrodr::cmd::make_io_compiler());
    cm.merge(macrodr::cmd::make_experiment_compiler());
    cm.merge(macrodr::cmd::make_model_compiler());
    cm.merge(macrodr::cmd::make_simulation_compiler());
    cm.merge(macrodr::cmd::make_likelihood_compiler());
    cm.merge(macrodr::cmd::make_dts_compiler());
    cm.push_function("load_experiment",
                     dsl::to_typed_function<std::string, double, double>(
                         &macrodr::cmd::load_experiment, "filename",
                         "frequency_of_sampling", "initial_ATP"));
    return cm;
}

}  // namespace macrodr
