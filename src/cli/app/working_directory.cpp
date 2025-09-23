#include <macrodr/cli/app/working_directory.h>

#include <filesystem>
#include <iostream>

namespace macrodr::cli::app {

Maybe_error<std::monostate> apply_working_directory(const macrodr::cli::CliOptions& opts) {
    if (!opts.has_chdir) {
        return std::monostate{};
    }

    std::error_code ec;
    std::filesystem::current_path(opts.chdir, ec);
    if (ec) {
        return error_message("Error: cannot chdir to '" + opts.chdir + "'");
    }
    return std::monostate{};
}

}  // namespace macrodr::cli::app

