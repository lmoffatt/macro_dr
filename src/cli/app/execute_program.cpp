#include <macrodr/cli/app/execute_program.h>

#include <iostream>

#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_untyped.h>

#include "CLI_macro_dr_base.h"  // for cmd::write_text

namespace macrodr::cli::app {

int execute_program(const std::string& script,
                    const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler& compiler) {
    cmd::write_text(cmd::temp_script_file, script);

    auto parsed = dsl::extract_program(script);
    std::cerr << parsed;
    if (!parsed) {
        std::cerr << parsed.error()();
        return 1;
    }

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);
    auto compiled = dsl::compile_program(env, parsed.value());
    if (!compiled) {
        std::cerr << "\n --------------Error--------\n"
                  << compiled.error()() << "\n --------------Error--------\n";
        return 1;
    }

    if (opts.check_syntax) {
        if (opts.verbosity > 0) {
            std::cout << "[macro_dr] syntax OK (compile-only)\n";
        }
        return 0;
    }

    auto exec = compiled.value().run(env);
    (void)exec;

    auto ss = parsed.value().str();
    std::cerr << ss;
    return 0;
}

}  // namespace macrodr::cli::app

