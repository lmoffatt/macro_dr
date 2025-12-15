// Lightweight, portable(ish) type-name utilities for diagnostics
#pragma once

#include <cctype>
#include <string>
#include <string_view>
#include <typeinfo>

#if defined(__GNUG__)
#include <cxxabi.h>

#include <cstdlib>
#endif

namespace macrodr::dsl {

inline std::string demangle(const char* name) {
#if defined(__GNUG__)
    int status = 0;
    char* dem = abi::__cxa_demangle(name, nullptr, nullptr, &status);
    std::string out = (status == 0 && (dem != nullptr)) ? dem : name;
    std::free(dem);
    return out;
#else
    // Fallback – returns compiler-provided name
    return name;
#endif
}

inline std::string strip_namespaces(std::string_view name) {
    std::string out;
    out.reserve(name.size());

    constexpr std::string_view anon_ns = "(anonymous namespace)::";

    for (std::size_t i = 0; i < name.size();) {
        if (name.substr(i, anon_ns.size()) == anon_ns) {
            i += anon_ns.size();
            continue;
        }

        const unsigned char c0 = static_cast<unsigned char>(name[i]);
        if (std::isalpha(c0) || name[i] == '_') {
            std::size_t j = i + 1;
            while (j < name.size()) {
                const unsigned char cj = static_cast<unsigned char>(name[j]);
                if (!(std::isalnum(cj) || name[j] == '_')) {
                    break;
                }
                ++j;
            }
            if (j + 1 < name.size() && name[j] == ':' && name[j + 1] == ':') {
                i = j + 2;
                continue;
            }
        }

        out.push_back(name[i]);
        ++i;
    }

    return out;
}

inline std::string demangle_no_namespace(const char* name) {
    return strip_namespaces(demangle(name));
}

template <class T>
inline std::string type_name() {
    return demangle(typeid(T).name());
}

template <class T>
inline std::string type_name(const T& v) {
    return demangle(typeid(v).name());
}

template <class T>
inline std::string type_name_no_namespace() {
    return demangle_no_namespace(typeid(T).name());
}

template <class T>
inline std::string type_name_no_namespace(const T& v) {
    return demangle_no_namespace(typeid(v).name());
}

}  // namespace macrodr::dsl
