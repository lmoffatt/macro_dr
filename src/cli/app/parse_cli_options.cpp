#include <macrodr/cli/app/parse_cli_options.h>

#include <iostream>

namespace macrodr::cli::app {
namespace {

void emit_cli_diagnostics(const macrodr::cli::CliOptions& opts) {
    for (const auto& diag : opts.diagnostics) {
        std::cerr << "[cli] " << diag << '\n';
    }
}

}  // namespace

Maybe_error<macrodr::cli::CliOptions> parse_cli_options(const CliArguments& arguments) {
    auto opts = macrodr::cli::parse_cli(arguments);
    emit_cli_diagnostics(opts);
    if (!opts.ok) {
        return error_message("[cli] failed to parse command line");
    }
    return opts;
}

}  // namespace macrodr::cli::app

