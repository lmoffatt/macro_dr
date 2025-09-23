#include <macrodr/cli/app/arguments.h>

#include <span>

namespace macrodr::cli::app {

CliArguments collect_arguments(int argc, char** argv) {
    CliArguments out;
    out.reserve(static_cast<std::size_t>(argc));
    for (const char* arg : std::span(argv, static_cast<std::size_t>(argc))) {
        out.emplace_back(arg ? arg : "");
    }
    return out;
}

}  // namespace macrodr::cli::app

