#pragma once

#include <string>
#include <filesystem>

#include <macrodr/cli/cli_parser.h>
#include <macrodr/dsl/lexer_typed.h>

namespace macrodr::cli::app {

int execute_program(const std::string& script,
                    const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler& compiler,
                    const std::filesystem::path& run_dir);

}  // namespace macrodr::cli::app
