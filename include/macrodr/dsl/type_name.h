// Lightweight, portable(ish) type-name utilities for diagnostics
#pragma once

#include <string>
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

template <class T>
inline std::string type_name() {
    return demangle(typeid(T).name());
}

template <class T>
inline std::string type_name(const T& v) {
    return demangle(typeid(v).name());
}

}  // namespace macrodr::dsl
