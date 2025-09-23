#pragma once

#include <optional>

#include <macrodr/cli/cli_parser.h>

namespace macrodr::cli::app {

std::optional<int> handle_immediate_flags(const macrodr::cli::CliOptions& opts);

}  // namespace macrodr::cli::app

