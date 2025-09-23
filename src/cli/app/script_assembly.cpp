#include <macrodr/cli/app/script_assembly.h>

#include <iostream>

#include <macrodr/cli/script_loader.h>

namespace macrodr::cli::app {

Maybe_error<std::string> assemble_program(const macrodr::cli::CliOptions& opts,
                                          const macrodr::io::PathResolver& resolver) {
    auto maybe_script = macrodr::cli::assemble_script(opts.sequence, resolver);
    if (!maybe_script) {
        return maybe_script.error();
    }
    return maybe_script.value();
}

void log_assembled_program(const std::string& script, const macrodr::cli::CliOptions& opts) {
    if (opts.verbosity > 0) {
        std::cout << "[macro_dr] assembled script:\n" << script << '\n';
    }
}

}  // namespace macrodr::cli::app

