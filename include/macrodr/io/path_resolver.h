// Path resolution utility for MacroDR file IO.
// Resolves relative paths against CWD, user-specified search paths, and MACRODR_PATH.

#pragma once

#include <string>
#include <vector>

#include "maybe_error.h"

namespace macrodr::io {

class PathResolver {
    std::vector<std::string> m_search_paths;  // ordered

   public:
    PathResolver();
    explicit PathResolver(std::vector<std::string> search_paths);

    // Returns MACRODR_PATH split by ':' (or ';') in order.
    static std::vector<std::string> env_paths();

    // Returns absolute path if the input exists after searching; otherwise error.
    // Search order: absolute as-is; CWD; each m_search_paths; each env_paths().
    Maybe_error<std::string> resolve_existing(const std::string& path) const;
};

}  // namespace macrodr::io

