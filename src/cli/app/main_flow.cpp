#include <macrodr/cli/app/main_flow.h>

#include <iostream>
#include <utility>

#include <macrodr/cli/app/arguments.h>
#include <macrodr/cli/app/execute_program.h>
#include <macrodr/cli/app/immediate_flags.h>
#include <macrodr/cli/app/parse_cli_options.h>
#include <macrodr/cli/app/script_assembly.h>
#include <macrodr/cli/app/workspace_persistence.h>
#include <macrodr/cli/app/working_directory.h>
#include <macrodr/cli/command_manager.h>
#include <macrodr/io/path_resolver.h>

#include "models_used.h"

namespace macrodr::cli::app {

int main_flow(int argc, char** argv) {
    print_model_Priors(2.0);

    auto arguments = collect_arguments(argc, argv);

    auto options = parse_cli_options(arguments);
    if (!options) {
        std::cerr << options.error()();
        std::cerr << '\n';
        return 1;
    }

    auto opts = std::move(options.value());

    if (auto exit_code = handle_immediate_flags(opts); exit_code) {
        return *exit_code;
    }

    auto chdir_result = apply_working_directory(opts);
    if (!chdir_result) {
        std::cerr << chdir_result.error()();
        std::cerr << '\n';
        return 2;
    }

    macrodr::io::PathResolver resolver(opts.search_paths);
    auto script = assemble_program(opts, resolver);
    if (!script) {
        std::cerr << "Error: \n" << script.error()();
        return 1;
    }

    auto script_value = std::move(script.value());

    log_assembled_program(script_value, opts);

    auto run_dir = persist_run_workspace(script_value, opts, arguments);

    auto compiler = macrodr::cli::make_compiler_new();
    return execute_program(script_value, opts, compiler, run_dir);
}

}  // namespace macrodr::cli::app
