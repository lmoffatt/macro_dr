#pragma once

#include <macrodr/cli/app/arguments.h>
#include <macrodr/cli/cli_parser.h>

#include "maybe_error.h"

namespace macrodr::cli::app {

Maybe_error<macrodr::cli::CliOptions> parse_cli_options(const CliArguments& arguments);

}  // namespace macrodr::cli::app

