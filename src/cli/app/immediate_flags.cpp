#include <macrodr/cli/app/immediate_flags.h>

#include <macrodr/cmd/cli_meta.h>

namespace macrodr::cli::app {

std::optional<int> handle_immediate_flags(const macrodr::cli::CliOptions& opts) {
    if (opts.help) {
        macrodr::cmd::cli_help();
        return 0;
    }
    if (opts.version) {
        macrodr::cmd::cli_version();
        return 0;
    }
    return std::nullopt;
}

}  // namespace macrodr::cli::app

