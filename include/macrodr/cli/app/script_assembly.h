#pragma once

#include <string>

#include <macrodr/cli/cli_parser.h>
#include <macrodr/io/path_resolver.h>

#include "maybe_error.h"

namespace macrodr::cli::app {

Maybe_error<std::string> assemble_program(const macrodr::cli::CliOptions& opts,
                                          const macrodr::io::PathResolver& resolver);

void log_assembled_program(const std::string& script, const macrodr::cli::CliOptions& opts);

}  // namespace macrodr::cli::app

