#pragma once

#include <variant>

#include <macrodr/cli/cli_parser.h>

#include "maybe_error.h"

namespace macrodr::cli::app {

Maybe_error<std::monostate> apply_working_directory(const macrodr::cli::CliOptions& opts);

}  // namespace macrodr::cli::app

