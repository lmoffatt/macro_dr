#ifndef COMMAND_MANAGER_H
#define COMMAND_MANAGER_H

#include <macrodr/cmd/init_commands.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/dsl/grammar_Identifier.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/lexer_untyped.h>

#include "CLI_macro_dr.h"
namespace macrodr {

namespace cli {

inline dsl::Compiler make_compiler_new() {
    dsl::Compiler cm;
    cm.push_function("load_experiment", dsl::to_typed_function<std::string, double, double>(
                                            &cmd::load_experiment, "filename",
                                            "frequency_of_sampling", "initial_ATP"));

    return cm;
}

}  // namespace cli
}  // namespace macrodr

#endif  // COMMAND_MANAGER_H
