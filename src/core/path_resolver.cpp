#include <macrodr/io/path_resolver.h>

#include <cstdlib>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

std::string expand_home(const std::string& path) {
    if (path.empty() || path[0] != '~') return path;
    const char* home = std::getenv("HOME");
    if (!home) return path;  // fallback: leave as-is
    if (path.size() == 1) return std::string(home);
    if (path[1] == '/') return std::string(home) + path.substr(1);
    return path;  // ~user not supported; leave as-is
}

std::vector<std::string> split_paths(const char* env) {
    std::vector<std::string> out;
    if (!env) return out;
    std::string s(env);
    char sep1 = ':';  // POSIX
    char sep2 = ';';  // Windows compatibility
    std::string cur;
    for (char c : s) {
        if (c == sep1 || c == sep2) {
            if (!cur.empty()) out.emplace_back(cur), cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.emplace_back(cur);
    return out;
}

}  // namespace

namespace macrodr::io {

PathResolver::PathResolver() : m_search_paths{} {}

PathResolver::PathResolver(std::vector<std::string> search_paths)
    : m_search_paths(std::move(search_paths)) {}

std::vector<std::string> PathResolver::env_paths() {
    return split_paths(std::getenv("MACRODR_PATH"));
}

Maybe_error<std::string> PathResolver::resolve_existing(const std::string& path) const {
    if (path.empty()) {
        return error_message("empty path");
    }

    auto try_path = [](const fs::path& p) -> Maybe_error<std::string> {
        std::error_code ec;
        if (fs::exists(p, ec)) {
            fs::path abs = fs::absolute(p, ec);
            if (!ec) return abs.string();
        }
        return error_message("");
    };

    const fs::path p_in = fs::path(expand_home(path));

    // Absolute path as-is
    if (p_in.is_absolute()) {
        auto r = try_path(p_in);
        if (r) return r.value();
        return error_message(path + " does not exist");
    }

    // CWD
    if (auto r = try_path(p_in); r) return r.value();

    // User-provided search paths
    for (const auto& base : m_search_paths) {
        fs::path candidate = fs::path(expand_home(base)) / p_in;
        if (auto r = try_path(candidate); r) return r.value();
    }

    // Environment search paths
    for (const auto& base : env_paths()) {
        fs::path candidate = fs::path(expand_home(base)) / p_in;
        if (auto r = try_path(candidate); r) return r.value();
    }

    return error_message(path + " not found in CWD or search paths");
}

}  // namespace macrodr::io

