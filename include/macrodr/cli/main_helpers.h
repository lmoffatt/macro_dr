#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <macrodr/cli/cli_parser.h>

namespace macrodr::dsl {
class Compiler;
}

namespace macrodr::cli::app {

std::vector<std::string> collect_arguments(int argc, char** argv);
void emit_cli_diagnostics(const macrodr::cli::CliOptions& opts);
std::optional<int> handle_immediate_flags(const macrodr::cli::CliOptions& opts);
bool apply_working_directory(const macrodr::cli::CliOptions& opts);
std::filesystem::path persist_run_workspace(const std::string& script,
                                            const macrodr::cli::CliOptions& opts,
                                            const std::vector<std::string>& arguments);
int execute_program(const std::string& script,
                    const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler& compiler);

}  // namespace macrodr::cli::app

