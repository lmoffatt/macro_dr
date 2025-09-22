// Assemble a MacroIR script from an ordered sequence of CLI events
// (files and inline eval lines), resolving file paths with PathResolver.

#pragma once

#include <string>
#include <vector>

#include <macrodr/cli/cli_parser.h>
#include <macrodr/io/path_resolver.h>
#include "maybe_error.h"

namespace macrodr::cli {

// Reads the given file using the resolver and appends its contents with trailing newlines.
Maybe_error<std::string> append_file(std::string&& acc, const std::string& path,
                                     const macrodr::io::PathResolver& resolver);

// Returns the assembled script or an error message referring to the failing file.
Maybe_error<std::string> assemble_script(const std::vector<CliEvent>& sequence,
                                         const macrodr::io::PathResolver& resolver);

}  // namespace macrodr::cli

