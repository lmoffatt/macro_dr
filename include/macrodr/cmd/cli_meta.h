// Lightweight CLI metadata commands exposed to the DSL
// - help(): prints brief usage
// - version(): prints and returns version string

#pragma once

#include <string>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/Lexer.h>
#include <macrodr/dsl/Compiler.h>

namespace macrodr::cmd {

// Prints a brief usage overview to stdout.
void cli_help();

// Prints version to stdout and also returns it.
std::string cli_version();

// Prints ONLY the build's git commit hash (no prefix) to stdout and returns it.
// The bare form is deliberate: dispatchers capture it as `commit=$(bin --commit)`
// to name the per-commit output directory (see ops/.../dispatch_figure_3.sh).
std::string cli_commit();

// Registry builder that registers help() and version() as DSL functions.
// This returns a dsl::Compiler<dsl::Lexer> containing only these commands; aggregate it
// via merge(...) in the CLI command manager.
dsl::Compiler<dsl::Lexer> make_cli_meta_compiler();

}  // namespace macrodr::cmd

