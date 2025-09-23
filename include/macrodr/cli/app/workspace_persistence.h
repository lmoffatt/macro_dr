#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <macrodr/cli/cli_parser.h>

namespace macrodr::cli::app {

std::filesystem::path persist_run_workspace(const std::string& script,
                                            const macrodr::cli::CliOptions& opts,
                                            const std::vector<std::string>& arguments);

}  // namespace macrodr::cli::app

