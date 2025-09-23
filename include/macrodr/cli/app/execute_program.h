#pragma once

#include <string>

#include <macrodr/cli/cli_parser.h>
#include <macrodr/dsl/lexer_typed.h>

namespace macrodr::cli::app {

int execute_program(const std::string& script,
                    const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler& compiler);

}  // namespace macrodr::cli::app

