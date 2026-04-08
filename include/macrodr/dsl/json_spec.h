#pragma once
#include <string>
#include <macrodr/io/json/convert.h>
#include "maybe_error.h"

namespace macrodr::dsl {

template <class L>
struct json_spec {
    using TagPolicy = macrodr::io::json::conv::TagPolicy;
    using Json = macrodr::io::json::Json;

    template <class T>
    static Json to_json(const T& v, TagPolicy policy) {
        return macrodr::io::json::conv::to_json(v, policy);
    }

    template <class T>
    static Maybe_error<void> from_json(
        const Json& j, T& out, const std::string& path, TagPolicy policy) {
        return macrodr::io::json::conv::from_json(j, out, path, policy);
    }
};

}