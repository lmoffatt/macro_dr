#pragma once

#include <string>
#include <filesystem>

#include <macrodr/cli/cli_parser.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/lexer_untyped.h>
#include<macrodr/dsl/Compiler.h>
#include<macrodr/dsl/Lexer.h>

namespace macrodr::cli::app {

int execute_program(const std::string& script,
                    const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler<macrodr::dsl::Lexer>& compiler,
                    const std::filesystem::path& run_dir);

}  // namespace macrodr::cli::app
