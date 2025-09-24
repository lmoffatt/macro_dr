// Minimal CLI parser: converts argv into ordered events and process options.
// Keeps domain semantics in the DSL by only providing file/eval events.

#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace macrodr::cli {

enum class CliEventKind { File, Eval };

struct CliEvent {
    CliEventKind kind{};
    std::string payload;  // path for File, DSL line for Eval
};

struct CliOptions {
    // Process-level flags
    bool help = false;
    bool version = false;
    bool check_syntax = false;  // parse/compile only
    int verbosity = 0;          // -v increments

    // Working-dir and search paths
    bool has_chdir = false;
    std::string chdir;                  // if has_chdir
    std::vector<std::string> search_paths;  // --path entries in order

    // Ordered sequence of events to assemble the script
    std::vector<CliEvent> sequence;

    // Diagnostics collected while parsing (warnings/errors)
    std::vector<std::string> diagnostics;
    bool ok = true;  // set false on hard parse errors

    // Environment persistence (JSON)
    enum class EnvSaveMode { Off, End, Step };
    EnvSaveMode env_save = EnvSaveMode::Off;  // --env-save off|end|step
    bool has_env_save_path = false;
    std::string env_save_path;               // --env-save-path <dir>
    bool has_env_load = false;
    std::string env_load;                    // --env-load <file>
    enum class EnvLoadMode { Append, Replace };
    EnvLoadMode env_load_mode = EnvLoadMode::Append;  // --env-load-mode append|replace
};

// Parse argv→CliOptions. Order-preserving for files and -e/--eval.
CliOptions parse_cli(int argc, const char* const argv[]);

// Convenience overload for vector<string>.
CliOptions parse_cli(const std::vector<std::string>& args);

}  // namespace macrodr::cli
