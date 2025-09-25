#pragma once

#include <string>
#include <type_traits>
#include <utility>

#include "macrodr/io/json/minijson.h"
#include "macrodr/io/json/convert.h"
#include "maybe_error.h"

namespace macrodr::dsl {

// Trait that indicates whether a type T can be constructed from a JSON literal.
// Defaults to false and is enabled when a matching conv::from_json overload exists.
namespace detail {

template <class T, class = void>
struct has_json_decoder : std::false_type {};

template <class T>
struct has_json_decoder<
    T, std::void_t<decltype(macrodr::io::json::conv::from_json(
        std::declval<const macrodr::io::json::Json&>(), std::declval<T&>(),
        std::declval<const std::string&>(), macrodr::io::json::conv::TagPolicy::None))>>
    : std::true_type {};

inline std::string trim_copy(const std::string& s) {
    constexpr const char* whitespace = " \t\r\n";
    const auto first = s.find_first_not_of(whitespace);
    if (first == std::string::npos) return std::string{};
    const auto last = s.find_last_not_of(whitespace);
    auto cropped = s.substr(first, last - first + 1);
    if (cropped.size() >= 2) {
        const char quote = cropped.front();
        if ((quote == '\'' || quote == '"') && cropped.back() == quote) {
            cropped = cropped.substr(1, cropped.size() - 2);
        }
    }
    return cropped;
}

}  // namespace detail

template <class T, class = void>
struct literal_decodable : std::false_type {};

template <class T>
struct literal_decodable<T, std::enable_if_t<detail::has_json_decoder<std::remove_cv_t<T>>::value>>
    : std::true_type {};

template <class T>
Maybe_error<std::remove_cv_t<T>> from_literal(const std::string& text)
    requires literal_decodable<T>::value {
    using Decayed = std::remove_cv_t<T>;
    using namespace macrodr::io::json;

    auto maybe_json = parse(text);
    if (!maybe_json) {
        if constexpr (std::is_same_v<Decayed, std::string>) {
            return detail::trim_copy(text);
        }
        return maybe_json.error();
    }

    Decayed value{};
    auto status = conv::from_json(maybe_json.value(), value, std::string{"$"},
                                  conv::TagPolicy::None);
    if (!status) {
        return status.error();
    }
    return value;
}

template <class T>
Maybe_error<std::remove_cv_t<T>> from_literal(const std::string&)
    requires(!literal_decodable<T>::value) = delete;

}  // namespace macrodr::dsl
