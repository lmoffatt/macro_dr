#include <macrodr/cli/app/execute_program.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_untyped.h>

#include <iostream>

#include "CLI_macro_dr_base.h"  // for cmd::write_text

#include <macrodr/io/json/environment_io.h>

namespace macrodr::cli::app {

int execute_program(const std::string& script, const macrodr::cli::CliOptions& opts,
                    macrodr::dsl::Compiler& compiler,
                    const std::filesystem::path& run_dir) {
    cmd::write_text(cmd::temp_script_file, script);

    auto parsed = dsl::extract_program(script);
    std::cerr << parsed;
    if (!parsed) {
        std::cerr << parsed.error()();
        return 1;
    }

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);

    // Optional: load environment before executing
    if (opts.has_env_load) {
        std::string mode = (opts.env_load_mode == macrodr::cli::CliOptions::EnvLoadMode::Append)
                               ? std::string("append")
                               : std::string("replace");
        auto loaded = macrodr::io::json::envio::load_environment_json<dsl::Lexer, dsl::Compiler>(
            opts.env_load, env, mode);
        if (!loaded) {
            std::cerr << "[env] load failed: " << loaded.error()() << "\n";
            return 1;
        }
        if (opts.verbosity > 0) {
            std::cout << "[env] " << loaded.value() << '\n';
        }
    }

    // Compile-only path
    {
        auto compiled = dsl::compile_program(env, parsed.value());
        if (!compiled) {
            std::cerr << "\n --------------Error--------\n" << compiled.error()()
                      << "\n --------------Error--------\n";
            return 1;
        }

        if (opts.check_syntax) {
            if (opts.verbosity > 0) {
                std::cout << "[macro_dr] syntax OK (compile-only)\n";
            }
            return 0;
        }
    }
    // Execute: step mode or full-program
    if (opts.env_save == macrodr::cli::CliOptions::EnvSaveMode::Step) {
        // Create snapshots dir
        std::filesystem::path snap_dir = run_dir / "snapshots";
        std::error_code ec;
        std::filesystem::create_directories(snap_dir, ec);
        (void)ec;

        int step = 0;
        for (const auto& stmt : parsed.value().statements()) {
            auto maybe_stmt = stmt->compile_statement(env);
            if (!maybe_stmt) {
                std::cerr << maybe_stmt.error()();
                return 1;
            }
            auto run_ok = maybe_stmt.value()->run_statement(env);
            if (!run_ok) {
                std::cerr << run_ok.error()();
                return 1;
            }
            // Save snapshot after each step
            auto snap_path = snap_dir / (std::string("step_") + std::to_string(step++) + ".json");
            auto saved = macrodr::io::json::envio::save_environment_json<dsl::Lexer, dsl::Compiler>(
                snap_path, env, compiler, /*functions*/ false, /*identifiers*/ true, "step");
            if (!saved && opts.verbosity > 0) {
                std::cerr << "[env] snapshot failed: " << saved.error()() << "\n";
            }
        }
        // Final save as well
        auto final_path = (opts.has_env_save_path ? std::filesystem::path(opts.env_save_path)
                                                  : run_dir) /
                          "environment.json";
        auto saved = macrodr::io::json::envio::save_environment_json<dsl::Lexer, dsl::Compiler>(
            final_path, env, compiler, /*functions*/ true, /*identifiers*/ true, "end");
        if (!saved && opts.verbosity > 0) {
            std::cerr << "[env] save failed: " << saved.error()() << "\n";
        }
    } else {
        // Full run via typed_program
        auto compiled2 = dsl::compile_program(env, parsed.value());
        if (!compiled2) {
            std::cerr << compiled2.error()();
            return 1;
        }
        auto exec = compiled2.value().run(env);
        if (!exec) {
            std::cerr << exec.error()();
            return 1;
        }
        if (opts.env_save == macrodr::cli::CliOptions::EnvSaveMode::End) {
            auto out_dir = opts.has_env_save_path ? std::filesystem::path(opts.env_save_path)
                                                  : run_dir;
            auto final_path = out_dir / "environment.json";
            auto saved = macrodr::io::json::envio::save_environment_json<dsl::Lexer, dsl::Compiler>(
                final_path, env, compiler, /*functions*/ true, /*identifiers*/ true, "end");
            if (!saved && opts.verbosity > 0) {
                std::cerr << "[env] save failed: " << saved.error()() << "\n";
            }
        }
    }

    auto ss = parsed.value().str();
    std::cerr << ss;
    return 0;
}

}  // namespace macrodr::cli::app
