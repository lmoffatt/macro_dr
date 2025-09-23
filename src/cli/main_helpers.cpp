#include <macrodr/cli/main_helpers.h>

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>

#include <macrodr/cli/script_loader.h>
#include <macrodr/cmd/cli_meta.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_untyped.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/io/path_resolver.h>

#include "CLI_macro_dr_base.h"  // for cmd::write_text

namespace macrodr::cli::app {

std::vector<std::string> collect_arguments(int argc, char** argv) {
    std::vector<std::string> out;
    out.reserve(static_cast<std::size_t>(argc));
    for (const char* a : std::span(argv, static_cast<std::size_t>(argc))) {
        out.emplace_back(a ? a : "");
    }
    return out;
}

void emit_cli_diagnostics(const macrodr::cli::CliOptions& opts) {
    for (const auto& diag : opts.diagnostics) {
        std::cerr << "[cli] " << diag << '\n';
    }
}

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

bool apply_working_directory(const macrodr::cli::CliOptions& opts) {
    if (!opts.has_chdir) {
        return true;
    }
    std::error_code ec;
    std::filesystem::current_path(opts.chdir, ec);
    if (ec) {
        std::cerr << "Error: cannot chdir to '" << opts.chdir << "'\n";
        return false;
    }
    return true;
}

std::filesystem::path persist_run_workspace(const std::string& script,
                                            const macrodr::cli::CliOptions& opts,
                                            const std::vector<std::string>& arguments) {
    const auto now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    tm = *std::localtime(&now);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "run-%Y%m%d-%H%M%S");
    const auto run_dir = std::filesystem::path("runs") / oss.str();

    std::error_code ec;
    std::filesystem::create_directories(run_dir, ec);
    if (ec) {
        std::cerr << "Warning: could not create run directory '" << run_dir.string()
                  << "' (" << ec.message() << ")\n";
        return {};
    }

    if (opts.verbosity > 0) {
        std::cout << "[macro_dr] run dir: " << run_dir.string() << '\n';
    }

    {
        std::ofstream script_out(run_dir / "script.macroir");
        script_out << script;
    }

    {
        std::ofstream meta(run_dir / "meta.json");
        meta << "{\n";
        meta << "  \"cwd\": \"" << std::filesystem::current_path().string() << "\",\n";
        meta << "  \"timestamp\": " << static_cast<long long>(now) << ",\n";
        meta << "  \"args\": [";
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            meta << (i == 0 ? "" : ", ") << "\"" << arguments[i] << "\"";
        }
        meta << "]\n";
        meta << "}\n";
    }

    return run_dir;
}

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

