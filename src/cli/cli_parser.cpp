#include <macrodr/cli/cli_parser.h>

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace macrodr::cli {

namespace {

inline bool starts_with(std::string_view s, std::string_view p) {
    return s.substr(0, p.size()) == p;
}

inline bool is_long_option(std::string_view s) { return starts_with(s, "--"); }
inline bool is_short_option(std::string_view s) {
    return s.size() >= 2 && s[0] == '-' && s[1] != '-';
}

}  // namespace

CliOptions parse_cli(int argc, const char* const argv[]) {
    CliOptions out;
    if (argc <= 0) return out;

    bool end_of_options = false;
    for (int i = 1; i < argc; ++i) {
        std::string_view tok(argv[i] ? argv[i] : "");
        if (!end_of_options && tok == "--") {
            end_of_options = true;
            continue;
        }

        auto need_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                out.ok = false;
                out.diagnostics.emplace_back(std::string(name) + " requires a value");
                return {};
            }
            ++i;
            return std::string(argv[i] ? argv[i] : "");
        };

        if (!end_of_options && (tok == "-h" || tok == "--help")) {
            out.help = true;
            continue;
        }
        if (!end_of_options && tok == "--version") {
            out.version = true;
            continue;
        }
        if (!end_of_options && (tok == "-n" || tok == "--check-syntax")) {
            out.check_syntax = true;
            continue;
        }
        if (!end_of_options && (tok == "-v" || tok == "--verbose")) {
            out.verbosity += 1;
            continue;
        }
        if (!end_of_options && (tok == "-C" || tok == "--chdir")) {
            auto v = need_value("--chdir");
            if (v.empty()) continue;
            out.has_chdir = true;
            out.chdir = std::move(v);
            continue;
        }
        if (!end_of_options && tok == "--path") {
            auto v = need_value("--path");
            if (v.empty()) continue;
            out.search_paths.emplace_back(std::move(v));
            continue;
        }
        if (!end_of_options && (tok == "-e" || tok == "--eval")) {
            auto line = need_value("--eval");
            if (line.empty() && !out.ok) continue;
            out.sequence.push_back(CliEvent{CliEventKind::Eval, std::move(line)});
            continue;
        }
        if (!end_of_options && (tok == "-f" || tok == "--file")) {
            auto path = need_value("--file");
            if (path.empty() && !out.ok) continue;
            out.sequence.push_back(CliEvent{CliEventKind::File, std::move(path)});
            continue;
        }

        // Compatibility shim: unknown long options that start with `--` are
        // treated as inline DSL lines (strip the `--`), with a warning.
        if (!end_of_options && is_long_option(tok)) {
            std::string dsl_line = std::string(tok.substr(2));
            out.diagnostics.emplace_back(
                std::string("treating legacy inline DSL as --eval: ") + std::string(tok));
            out.sequence.push_back(CliEvent{CliEventKind::Eval, std::move(dsl_line)});
            continue;
        }

        // Positional or after `--`: treat as script file
        out.sequence.push_back(CliEvent{CliEventKind::File, std::string(tok)});
    }

    return out;
}

CliOptions parse_cli(const std::vector<std::string>& args) {
    CliOptions out;
    if (args.empty()) return out;
    std::vector<const char*> cargs;
    cargs.reserve(args.size());
    for (auto const& s : args) cargs.push_back(s.c_str());
    return parse_cli(static_cast<int>(cargs.size()), cargs.data());
}

}  // namespace macrodr::cli

