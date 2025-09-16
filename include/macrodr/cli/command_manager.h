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


namespace macrodr::cli {

inline dsl::Compiler make_compiler_new() {
    dsl::Compiler cm;
    // Intentionally empty: command registrations live under macrodr::cmd builders
    return cm;
}

} // namespace macrodr::cli


#endif  // COMMAND_MANAGER_H
