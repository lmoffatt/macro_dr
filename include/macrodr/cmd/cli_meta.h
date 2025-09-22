// Lightweight CLI metadata commands exposed to the DSL
// - help(): prints brief usage
// - version(): prints and returns version string

#pragma once

#include <string>
#include <macrodr/dsl/lexer_typed.h>

namespace macrodr::cmd {

// Prints a brief usage overview to stdout.
void cli_help();

// Prints version to stdout and also returns it.
std::string cli_version();

// Registry builder that registers help() and version() as DSL functions.
// This returns a dsl::Compiler containing only these commands; aggregate it
// via merge(...) in the CLI command manager.
dsl::Compiler make_cli_meta_compiler();

}  // namespace macrodr::cmd

