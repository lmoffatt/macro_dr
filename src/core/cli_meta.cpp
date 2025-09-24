#include <macrodr/cmd/cli_meta.h>

#include <iostream>
#include <string>

// Factory helpers for DSL function registration
#include "CLI_macro_dr.h"

namespace {

constexpr const char* kUsage =
    "Usage:\n"
    "  macro_dr [options] <script1> [script2 ...] [-e \"<dsl>\"]...\n"
    "\nOptions (process-level):\n"
    "  -e, --eval <dsl>      Append a DSL line (repeatable, order-preserving)\n"
    "  -n, --check-syntax    Parse/compile only; do not run\n"
    "  -v, --verbose         Increase logging verbosity (repeatable)\n"
    "  -C, --chdir <dir>     Change working directory before running\n"
    "      --path <dir>      Add a search path for relative assets (repeatable)\n"
    "      --env-save <mode> Save environment: off|end|step (default: off)\n"
    "      --env-save-path <dir>\n"
    "                         Override directory for environment snapshots\n"
    "      --env-load <file>  Load environment JSON before execution\n"
    "      --env-load-mode <mode>\n"
    "                         append (default) or replace when loading\n"
    "      --version         Print version and exit\n"
    "  -h, --help            Print this help and exit\n"
    "  --                    End of options; remaining args are script files\n"
    "\nNotes:\n"
    "  - Scripts and --eval lines are processed in order.\n"
    "  - This help/version are also available as DSL commands: help(), version().\n"
    "  - Step mode writes snapshots to runs/<run-id>/snapshots by default.\n";

#ifndef MACRODR_VERSION
#define MACRODR_VERSION "0.0.0-dev"
#endif

}  // namespace

namespace macrodr::cmd {

void cli_help() {
    std::cout << kUsage;
}

std::string cli_version() {
    std::string v = std::string{"macro_dr "} + MACRODR_VERSION;
    std::cout << v << '\n';
    return v;
}

dsl::Compiler make_cli_meta_compiler() {
    dsl::Compiler cm;
    cm.push_function("help", dsl::to_typed_function<>(&cli_help));
    cm.push_function("version", dsl::to_typed_function<>(&cli_version));
    return cm;
}

}  // namespace macrodr::cmd
