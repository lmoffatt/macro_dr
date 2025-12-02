#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "minijson.h"
#include "maybe_error.h"
#include "matrix.h"
#include "macrodr/dsl/type_name.h"
#include "variables.h"
#include "derivative_operator.h"
#include "matrix_derivative.h"
// IModel interface for model pointer serialization
#include <macrodr/interface/IModel.h>
// Bring legacy parameter types for optional JSON views
#include "parameters.h"

namespace macrodr::io::json::conv {

enum class TagPolicy { None, Minimal, Full };

template <class T, class = void>
struct has_json_codec : std::false_type {};

template <class T>
struct has_json_codec<
    T, std::void_t<decltype(to_json(std::declval<const T&>(), TagPolicy::None)),
                   decltype(from_json(std::declval<const Json&>(), std::declval<T&>(),
                                      std::declval<const std::string&>(), TagPolicy::None))>>
    : std::true_type {};

template <class T>
inline constexpr bool has_json_codec_v = has_json_codec<T>::value;

namespace detail {

inline Json number_or_null(double v) {
    if (std::isfinite(v)) {
        return Json::number(v);
    }
    return Json::null();
}

template <class T>
Json value_to_json(const T& value, TagPolicy policy) {
    if constexpr (std::is_same_v<std::decay_t<T>, double>) {
        return number_or_null(static_cast<double>(value));
    } else if constexpr (requires { to_json(value, policy); }) {
        return to_json(value, policy);
    } else if constexpr (requires { value(); to_json(value(), policy); }) {
        return to_json(value(), policy);
    } else {
        Json obj = Json::object();
        obj.obj.emplace("type", Json::string(macrodr::dsl::type_name<T>()));
        obj.obj.emplace("value", Json::null());
        return obj;
    }
}

inline std::string json_type_name(Json::Type t) {
    switch (t) {
        case Json::Type::Null: return "null";
        case Json::Type::Bool: return "bool";
        case Json::Type::Number: return "number";
        case Json::Type::String: return "string";
        case Json::Type::Array: return "array";
        case Json::Type::Object: return "object";
    }
    return "unknown";
}

inline Maybe_error<void> type_error(const std::string& path, const std::string& expected,
                                    Json::Type actual) {
    return error_message(path + ": expected " + expected + ", got " + json_type_name(actual));
}

template <class T, class = void>
struct has_class_name : std::false_type {};

template <class T>
struct has_class_name<T, std::void_t<decltype(className(T{}))>> : std::true_type {};

template <class T>
std::string view_class_name() {
    if constexpr (has_class_name<T>::value) {
        return className(T{});
    } else {
        return macrodr::dsl::type_name<T>();
    }
}

template <class T>
struct is_matrix : std::false_type {};

template <class V>
struct is_matrix<Matrix<V>> : std::true_type {};

template <class T>
constexpr bool is_matrix_v = is_matrix<std::decay_t<T>>::value;

template <class T>
struct is_vector_space : std::false_type {};

template <class... Vars>
struct is_vector_space<var::Vector_Space<Vars...>> : std::true_type {};

template <class T>
constexpr bool is_vector_space_v = is_vector_space<std::decay_t<T>>::value;

template <class T, class = void>
struct is_var : std::false_type {};

template <class T>
struct is_var<T, std::void_t<typename T::variable_type>> : std::true_type {};

template <class Id, class Value>
struct var_traits_impl {
    using id = Id;
    using value = Value;
};

template <class T>
struct var_traits;

template <class Id, class Value>
struct var_traits<var::Var<Id, Value>> : var_traits_impl<Id, Value> {};

template <class Id, class Value>
struct var_traits<var::Constant<Id, Value>> : var_traits_impl<Id, Value> {};

// Derivative wrapper: forward to the underlying untransformed type
template <class T, class X>
struct var_traits<var::Derivative<T, X>> : var_traits<var::untransformed_type_t<var::Derivative<T, X>>> {};

// Access the underlying stored value of a Var/Constant by explicitly casting to its base
// to avoid ambiguity with other inherited members (e.g., derivative wrappers).
template <class VarType>
auto& var_value(VarType& v) {
    using Base = typename VarType::variable_type;
    return static_cast<Base&>(v).value();
}
template <class VarType>
auto const& var_value(const VarType& v) {
    using Base = typename VarType::variable_type;
    return static_cast<const Base&>(v).value();
}

template <class T>
struct var_traits {
    using base = typename T::variable_type;
    using id = typename var_traits<base>::id;
    using value = typename var_traits<base>::value;
};

template <class T>
std::string var_id_name() {
    using Id = typename var_traits<T>::id;
    return view_class_name<Id>();
}

template <class T>
std::string value_type_name() {
    using Value = typename var_traits<T>::value;
    return macrodr::dsl::type_name<Value>();
}

template <class Tuple, std::size_t... Is>
void tuple_push_json(Json& arr, const Tuple& tu, TagPolicy policy, std::index_sequence<Is...>);

template <class Tuple, std::size_t... Is>
Maybe_error<void> tuple_from_json(const Json& arr, Tuple& tu, const std::string& path,
                                  TagPolicy policy, std::index_sequence<Is...>);

template <std::size_t I, class VS, class... Vars>
Maybe_error<void> assign_vector_space_entry(const Json& value, const std::string& id,
                                            VS& vs, std::array<bool, sizeof...(Vars)>& seen,
                                            const std::string& path, TagPolicy policy);

template <std::size_t I, class... Vars>
Maybe_error<void> ensure_all_seen(const std::array<bool, sizeof...(Vars)>& seen,
                                  const std::string& path);

}  // namespace detail

// Primitive overloads
Json to_json(bool value, TagPolicy);
Json to_json(double value, TagPolicy);
Json to_json(const std::string& value, TagPolicy);

Maybe_error<void> from_json(const Json& j, bool& out, const std::string& path, TagPolicy);
Maybe_error<void> from_json(const Json& j, double& out, const std::string& path, TagPolicy);
Maybe_error<void> from_json(const Json& j, std::string& out, const std::string& path, TagPolicy);

template <class Int>
Json to_json(Int value, TagPolicy policy)
    requires(std::is_integral_v<Int> && std::is_signed_v<Int> && !std::is_same_v<Int, bool>);

template <class Int>
Maybe_error<void> from_json(const Json& j, Int& out, const std::string& path, TagPolicy policy)
    requires(std::is_integral_v<Int> && std::is_signed_v<Int> && !std::is_same_v<Int, bool>);

template <class UInt>
Json to_json(UInt value, TagPolicy policy)
    requires(std::is_integral_v<UInt> && std::is_unsigned_v<UInt> &&
             !std::is_same_v<UInt, bool>);

template <class UInt>
Maybe_error<void> from_json(const Json& j, UInt& out, const std::string& path, TagPolicy policy)
    requires(std::is_integral_v<UInt> && std::is_unsigned_v<UInt> &&
             !std::is_same_v<UInt, bool>);

// Vector conversions
template <class T>
Json to_json(const std::vector<T>& values, TagPolicy policy = TagPolicy::None);

template <class T>
Maybe_error<void> from_json(const Json& j, std::vector<T>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

// Map conversions
template <class V>
Json to_json(const std::map<std::string, V>& m, TagPolicy policy = TagPolicy::None);

template <class V>
Maybe_error<void> from_json(const Json& j, std::map<std::string, V>& out,
                            const std::string& path, TagPolicy policy = TagPolicy::None);

template <class K, class V>
Json to_json(const std::map<K, V>& m, TagPolicy policy = TagPolicy::None)
    requires(!std::is_same_v<K, std::string>);

template <class K, class V>
Maybe_error<void> from_json(const Json& j, std::map<K, V>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None)
    requires(!std::is_same_v<K, std::string>);

// Matrix conversions
Json to_json(const Matrix<double>& m, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, Matrix<double>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

Json to_json(const DiagonalMatrix<double>& m, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, DiagonalMatrix<double>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

// Pair and tuple
template <class A, class B>
Json to_json(const std::pair<A, B>& value, TagPolicy policy = TagPolicy::None);

template <class A, class B>
Maybe_error<void> from_json(const Json& j, std::pair<A, B>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

template <class... Ts>
Json to_json(const std::tuple<Ts...>& value, TagPolicy policy = TagPolicy::None);

template <class... Ts>
Maybe_error<void> from_json(const Json& j, std::tuple<Ts...>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

// Unique pointer – serialize pointee if possible, else type-tag/null (no from_json)
template <class T>
Json to_json(const std::unique_ptr<T>& p, TagPolicy policy = TagPolicy::None);

// Parameters_values – reversible shape references schema by id
Json to_json(const var::Parameters_values& pv, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, var::Parameters_values& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

// Transformations and schema (reversible)
Json to_json(const var::transformations_vector& tr, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, var::transformations_vector& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

Json to_json(const var::Parameters_Transformations& pt, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, var::Parameters_Transformations& out,
                            const std::string& path, TagPolicy policy = TagPolicy::None);

// Legacy matrix specializations frequently found inside Vars
Json to_json(const SymPosDefMatrix<double>& m, TagPolicy policy = TagPolicy::None);
Json to_json(const SymmetricMatrix<double>& m, TagPolicy policy = TagPolicy::None);
Json to_json(const DiagPosDetMatrix<double>& m, TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, SymPosDefMatrix<double>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, SymmetricMatrix<double>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);
Maybe_error<void> from_json(const Json& j, DiagPosDetMatrix<double>& out, const std::string& path,
                            TagPolicy policy = TagPolicy::None);

// var::Var and Vector_Space
template <class VarType>
Json to_json(const VarType& value, TagPolicy policy = TagPolicy::None)
    requires(detail::is_var<VarType>::value && !detail::is_vector_space_v<VarType>);

template <class VarType>
Maybe_error<void> from_json(const Json& j, VarType& out, const std::string& path,
                            TagPolicy policy)
    requires(detail::is_var<VarType>::value && !detail::is_vector_space_v<VarType>);

template <class... Vars>
Json to_json(const var::Vector_Space<Vars...>& value, TagPolicy policy = TagPolicy::None);

template <class... Vars>
Maybe_error<void> from_json(const Json& j, var::Vector_Space<Vars...>& out,
                            const std::string& path, TagPolicy policy = TagPolicy::None);

// -----------------------------------------------------------------------------
// Implementations
// -----------------------------------------------------------------------------

inline Json to_json(bool value, TagPolicy) { return Json::boolean(value); }
inline Json to_json(double value, TagPolicy) { return Json::number(value); }
inline Json to_json(const std::string& value, TagPolicy) { return Json::string(value); }

template <class Int>
Json to_json(Int value, TagPolicy)
    requires(std::is_integral_v<Int> && std::is_signed_v<Int> && !std::is_same_v<Int, bool>) {
    return Json::number(static_cast<double>(value));
}

inline Maybe_error<void> from_json(const Json& j, bool& out, const std::string& path,
                                   TagPolicy) {
    if (j.type == Json::Type::Bool) {
        out = j.b;
        return {};
    }
    if (j.type == Json::Type::Number) {
        double rounded = std::round(j.num);
        if (std::abs(rounded - j.num) > std::numeric_limits<double>::epsilon()) {
            return error_message(path + ": expected boolean (0/1), got " + std::to_string(j.num));
        }
        if (rounded != 0.0 && rounded != 1.0) {
            return error_message(path + ": expected boolean (0/1), got " + std::to_string(j.num));
        }
        out = (rounded != 0.0);
        return {};
    }
    return detail::type_error(path, "bool", j.type);
}

inline Maybe_error<void> from_json(const Json& j, double& out, const std::string& path,
                                   TagPolicy) {
    if (j.type != Json::Type::Number) {
        return detail::type_error(path, "number", j.type);
    }
    out = j.num;
    return {};
}

template <class Int>
Maybe_error<void> from_json(const Json& j, Int& out, const std::string& path, TagPolicy)
    requires(std::is_integral_v<Int> && std::is_signed_v<Int> && !std::is_same_v<Int, bool>) {
    if (j.type != Json::Type::Number) {
        return detail::type_error(path, "number", j.type);
    }
    double v = j.num;
    double rounded = std::round(v);
    if (std::abs(rounded - v) > std::numeric_limits<double>::epsilon()) {
        return error_message(path + ": expected integer number, got " + std::to_string(v));
    }
    if (rounded < static_cast<double>(std::numeric_limits<Int>::lowest()) ||
        rounded > static_cast<double>(std::numeric_limits<Int>::max())) {
        return error_message(path + ": integer out of range");
    }
    out = static_cast<Int>(rounded);
    return {};
}

inline Maybe_error<void> from_json(const Json& j, std::string& out, const std::string& path,
                                   TagPolicy) {
    if (j.type != Json::Type::String) {
        return detail::type_error(path, "string", j.type);
    }
    out = j.str;
    return {};
}

template <class UInt>
Json to_json(UInt value, TagPolicy)
    requires(std::is_integral_v<UInt> && std::is_unsigned_v<UInt> &&
             !std::is_same_v<UInt, bool>) {
    return Json::number(static_cast<double>(value));
}

template <class UInt>
Maybe_error<void> from_json(const Json& j, UInt& out, const std::string& path, TagPolicy)
    requires(std::is_integral_v<UInt> && std::is_unsigned_v<UInt> &&
             !std::is_same_v<UInt, bool>) {
    if (j.type != Json::Type::Number) {
        return detail::type_error(path, "number", j.type);
    }
    double v = j.num;
    double rounded = std::round(v);
    if (std::abs(rounded - v) > std::numeric_limits<double>::epsilon()) {
        return error_message(path + ": expected integer number, got " + std::to_string(v));
    }
    if (rounded < 0.0) {
        return error_message(path + ": expected non-negative integer");
    }
    if (rounded > static_cast<double>(std::numeric_limits<UInt>::max())) {
        return error_message(path + ": integer out of range");
    }
    out = static_cast<UInt>(rounded);
    return {};
}

template <class T>
Json to_json(const std::vector<T>& values, TagPolicy policy) {
    Json arr = Json::array();
    for (const auto& v : values) {
        arr.arr.push_back(to_json(v, policy));
    }
    return arr;
}

template <class T>
Maybe_error<void> from_json(const Json& j, std::vector<T>& out, const std::string& path,
                            TagPolicy policy) {
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    out.clear();
    for (std::size_t i = 0; i < j.arr.size(); ++i) {
        T element{};
        auto status = from_json(j.arr[i], element, path + "[" + std::to_string(i) + "]", policy);
        if (!status) {
            return status;
        }
        out.push_back(std::move(element));
    }
    return {};
}

template <class V>
Json to_json(const std::map<std::string, V>& m, TagPolicy policy) {
    Json obj = Json::object();
    for (const auto& [key, value] : m) {
        obj.obj.emplace(key, to_json(value, policy));
    }
    return obj;
}

template <class V>
Maybe_error<void> from_json(const Json& j, std::map<std::string, V>& out,
                            const std::string& path, TagPolicy policy) {
    out.clear();
    if (j.type == Json::Type::Object) {
        for (const auto& [key, value] : j.obj) {
            V decoded{};
            auto status = from_json(value, decoded, path + "." + key, policy);
            if (!status) {
                return status;
            }
            out.emplace(key, std::move(decoded));
        }
        return {};
    }
    if (j.type == Json::Type::Array) {
        for (std::size_t i = 0; i < j.arr.size(); ++i) {
            const auto& entry = j.arr[i];
            if (entry.type != Json::Type::Object) {
                return detail::type_error(path + "[" + std::to_string(i) + "]", "object", entry.type);
            }
            const auto* key_json = entry.find("key");
            const auto* value_json = entry.find("value");
            if (!key_json || !value_json) {
                return error_message(path + "[" + std::to_string(i) + "] missing 'key' or 'value'");
            }
            std::string key;
            auto key_status = from_json(*key_json, key, path + "[" + std::to_string(i) + "].key", policy);
            if (!key_status) {
                return key_status;
            }
            V decoded{};
            auto value_status = from_json(*value_json, decoded,
                                          path + "[" + std::to_string(i) + "].value", policy);
            if (!value_status) {
                return value_status;
            }
            out.emplace(std::move(key), std::move(decoded));
        }
        return {};
    }
    return detail::type_error(path, "object", j.type);
}

template <class K, class V>
Json to_json(const std::map<K, V>& m, TagPolicy policy)
    requires(!std::is_same_v<K, std::string>) {
    Json arr = Json::array();
    arr.arr.reserve(m.size());
    for (const auto& [key, value] : m) {
        Json entry = Json::object();
        entry.obj.emplace("key", to_json(key, policy));
        entry.obj.emplace("value", to_json(value, policy));
        arr.arr.push_back(std::move(entry));
    }
    return arr;
}

template <class K, class V>
Maybe_error<void> from_json(const Json& j, std::map<K, V>& out, const std::string& path,
                            TagPolicy policy)
    requires(!std::is_same_v<K, std::string>) {
    out.clear();
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    for (std::size_t i = 0; i < j.arr.size(); ++i) {
        const auto& entry = j.arr[i];
        if (entry.type != Json::Type::Object) {
            return detail::type_error(path + "[" + std::to_string(i) + "]", "object", entry.type);
        }
        const auto* key_json = entry.find("key");
        const auto* value_json = entry.find("value");
        if (!key_json || !value_json) {
            return error_message(path + "[" + std::to_string(i) + "] missing 'key' or 'value'");
        }
        K key{};
        auto key_status = from_json(*key_json, key, path + "[" + std::to_string(i) + "].key", policy);
        if (!key_status) {
            return key_status;
        }
        V decoded{};
        auto value_status = from_json(*value_json, decoded,
                                      path + "[" + std::to_string(i) + "].value", policy);
        if (!value_status) {
            return value_status;
        }
        out.emplace(std::move(key), std::move(decoded));
    }
    return {};
}

inline Json to_json(const Matrix<double>& m, TagPolicy) {
    Json rows = Json::array();
    for (std::size_t i = 0; i < m.nrows(); ++i) {
        Json row = Json::array();
        for (std::size_t jcol = 0; jcol < m.ncols(); ++jcol) {
            row.arr.push_back(detail::number_or_null(m(i, jcol)));
        }
        rows.arr.push_back(std::move(row));
    }
    return rows;
}

// Minimal serializer for index matrices used by spectral Blocks (std::size_t)
inline Json to_json(const Matrix<std::size_t>& m, TagPolicy) {
    Json rows = Json::array();
    for (std::size_t i = 0; i < m.nrows(); ++i) {
        Json row = Json::array();
        for (std::size_t jcol = 0; jcol < m.ncols(); ++jcol) {
            row.arr.push_back(Json::number(static_cast<double>(m(i, jcol))));
        }
        rows.arr.push_back(std::move(row));
    }
    return rows;
}

inline Maybe_error<void> from_json(const Json& j, Matrix<double>& out, const std::string& path,
                                   TagPolicy) {
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    if (j.arr.empty()) {
        out = Matrix<double>();
        return {};
    }
    if (j.arr.front().type != Json::Type::Array) {
        return detail::type_error(path + "[0]", "array", j.arr.front().type);
    }
    std::size_t rows = j.arr.size();
    std::size_t cols = j.arr.front().arr.size();
    Matrix<double> m(rows, cols, false);
    for (std::size_t i = 0; i < rows; ++i) {
        const auto& row = j.arr[i];
        if (row.type != Json::Type::Array) {
            return detail::type_error(path + "[" + std::to_string(i) + "]", "array", row.type);
        }
        if (row.arr.size() != cols) {
            return error_message(path + ": inconsistent row lengths");
        }
        for (std::size_t jcol = 0; jcol < cols; ++jcol) {
            const auto& cell = row.arr[jcol];
            if (cell.type != Json::Type::Number) {
                return detail::type_error(path + "[" + std::to_string(i) + "][" +
                                             std::to_string(jcol) + "]",
                                         "number", cell.type);
            }
            m(i, jcol) = cell.num;
        }
    }
    out = std::move(m);
    return {};
}

inline Maybe_error<void> from_json(const Json& j, Matrix<std::size_t>& out, const std::string& path,
                                   TagPolicy) {
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    if (j.arr.empty()) {
        out = Matrix<std::size_t>();
        return {};
    }
    if (j.arr.front().type != Json::Type::Array) {
        return detail::type_error(path + "[0]", "array", j.arr.front().type);
    }
    std::size_t rows = j.arr.size();
    std::size_t cols = j.arr.front().arr.size();
    Matrix<std::size_t> m(rows, cols, false);
    for (std::size_t i = 0; i < rows; ++i) {
        const auto& row = j.arr[i];
        if (row.type != Json::Type::Array) {
            return detail::type_error(path + "[" + std::to_string(i) + "]", "array", row.type);
        }
        if (row.arr.size() != cols) {
            return error_message(path + ": inconsistent row lengths");
        }
        for (std::size_t jcol = 0; jcol < cols; ++jcol) {
            const auto& cell = row.arr[jcol];
            if (cell.type != Json::Type::Number) {
                return detail::type_error(path + "[" + std::to_string(i) + "][" +
                                             std::to_string(jcol) + "]",
                                         "number", cell.type);
            }
            double v = cell.num;
            if (v < 0) v = 0;  // clamp
            m(i, jcol) = static_cast<std::size_t>(v);
        }
    }
    out = std::move(m);
    return {};
}

inline Json to_json(const DiagonalMatrix<double>& m, TagPolicy) {
    Json obj = Json::object();
    obj.obj.emplace("rows", Json::number(static_cast<double>(m.nrows())));
    obj.obj.emplace("cols", Json::number(static_cast<double>(m.ncols())));
    Json diag = Json::array();
    for (std::size_t i = 0; i < m.size(); ++i)
        diag.arr.push_back(detail::number_or_null(m[i]));
    obj.obj.emplace("diag", std::move(diag));
    return obj;
}

inline Maybe_error<void> from_json(const Json& j, DiagonalMatrix<double>& out, const std::string& path,
                                   TagPolicy policy) {
    if (j.type != Json::Type::Object && j.type != Json::Type::Array) {
        return detail::type_error(path, "object", j.type);
    }

    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> diag_values;

    if (j.type == Json::Type::Array) {
        auto status = from_json(j, diag_values, path, policy);
        if (!status) {
            return status;
        }
        rows = cols = diag_values.size();
    } else {
        const auto* diag_json = j.find("diag");
        if (!diag_json || diag_json->type != Json::Type::Array) {
            return error_message(path + ": missing 'diag' array");
        }
        auto status = from_json(*diag_json, diag_values, path + ".diag", policy);
        if (!status) {
            return status;
        }
        rows = cols = diag_values.size();
        if (const auto* rows_json = j.find("rows")) {
            auto rows_status = from_json(*rows_json, rows, path + ".rows", policy);
            if (!rows_status) {
                return rows_status;
            }
        }
        if (const auto* cols_json = j.find("cols")) {
            auto cols_status = from_json(*cols_json, cols, path + ".cols", policy);
            if (!cols_status) {
                return cols_status;
            }
        }
    }

    if (diag_values.empty() && (rows == 0 || cols == 0)) {
        out = DiagonalMatrix<double>();
        return {};
    }

    if (rows == 0 || cols == 0) {
        return error_message(path + ": invalid matrix dimensions");
    }

    const std::size_t diag_size = std::min(rows, cols);
    if (diag_values.size() != diag_size) {
        return error_message(path + ": expected " + std::to_string(diag_size) +
                             " diagonal entries, got " + std::to_string(diag_values.size()));
    }

    DiagonalMatrix<double> m(rows, cols, true);
    for (std::size_t i = 0; i < diag_size; ++i) {
        m[i] = diag_values[i];
    }
    out = std::move(m);
    return {};
}

// -----------------------------------------------------------------------------
// Legacy matrix flavors → serialize densely like Matrix<double>, or diag-only
// -----------------------------------------------------------------------------
inline Json to_json(const SymPosDefMatrix<double>& m, TagPolicy policy) {
    (void)policy;
    Json rows = Json::array();
    // Avoid reserve on potentially bad sizes; push back progressively
    for (std::size_t i = 0; i < m.nrows(); ++i) {
        Json row = Json::array();
        for (std::size_t j = 0; j < m.ncols(); ++j) {
            row.arr.push_back(detail::number_or_null(m(i, j)));
        }
        rows.arr.push_back(std::move(row));
    }
    return rows;
}

inline Json to_json(const SymmetricMatrix<double>& m, TagPolicy policy) {
    (void)policy;
    Json rows = Json::array();
    for (std::size_t i = 0; i < m.nrows(); ++i) {
        Json row = Json::array();
        for (std::size_t j = 0; j < m.ncols(); ++j) {
            row.arr.push_back(detail::number_or_null(m(i, j)));
        }
        rows.arr.push_back(std::move(row));
    }
    return rows;
}

inline Json to_json(const DiagPosDetMatrix<double>& m, TagPolicy policy) {
    (void)policy;
    Json obj = Json::object();
    obj.obj.emplace("rows", Json::number(static_cast<double>(m.nrows())));
    obj.obj.emplace("cols", Json::number(static_cast<double>(m.ncols())));
    Json diag = Json::array();
    for (std::size_t i = 0; i < m.size(); ++i)
        diag.arr.push_back(detail::number_or_null(m[i]));
    obj.obj.emplace("diag", std::move(diag));
    return obj;
}

inline Maybe_error<void> from_json(const Json& j, SymPosDefMatrix<double>& out, const std::string& path,
                                   TagPolicy policy) {
    (void)j; (void)out; (void)policy;
    return error_message(path + ": deserialization for SymPosDefMatrix<double> not supported");
}
inline Maybe_error<void> from_json(const Json& j, SymmetricMatrix<double>& out, const std::string& path,
                                   TagPolicy policy) {
    (void)j; (void)out; (void)policy;
    return error_message(path + ": deserialization for SymmetricMatrix<double> not supported");
}
inline Maybe_error<void> from_json(const Json& j, DiagPosDetMatrix<double>& out, const std::string& path,
                                   TagPolicy policy) {
    (void)j; (void)out; (void)policy;
    return error_message(path + ": deserialization for DiagPosDetMatrix<double> not supported");
}

template <class Value, class DX>
Json to_json(const var::d_d_<Value, DX>& d, TagPolicy policy) {
    return detail::value_to_json(d(), policy);
}

template <class Value, class DX>
Maybe_error<void> from_json(const Json& j, var::d_d_<Value, DX>& out, const std::string& path,
                           TagPolicy policy) {
    (void)j;
    (void)out;
    (void)policy;
    return error_message(path + ": deserialization for derivative components not supported");
}

template <class Primitive, class Differential>
Json to_json(const var::Derivative<Primitive, Differential>& value, TagPolicy policy) {
    Json obj = Json::object();
    obj.obj.emplace("primitive", detail::value_to_json(value.primitive(), policy));
    obj.obj.emplace("derivative", detail::value_to_json(value.derivative(), policy));
    return obj;
}

// Specialization for derivatives of vectors: serialize elementwise derivatives
template <class T, class DX>
Json to_json(const var::Derivative<std::vector<T>, DX>& value, TagPolicy policy) {
    Json obj = Json::object();
    // Primitive is the underlying std::vector<T>
    obj.obj.emplace("primitive", detail::value_to_json(value.primitive(), policy));

    // Derivative: array with each element's derivative JSON
    Json darr = Json::array();
    for (const auto& elem : value) {
        // elem is Derivative<T, DX>
        darr.arr.push_back(detail::value_to_json(elem.derivative(), policy));
    }
    obj.obj.emplace("derivative", std::move(darr));
    return obj;
}

template <class Primitive, class Differential>
Maybe_error<void> from_json(const Json& j, var::Derivative<Primitive, Differential>& out,
                           const std::string& path, TagPolicy policy) {
    (void)j;
    (void)out;
    (void)policy;
    return error_message(path + ": deserialization for derivatives is not supported");
}

template <class A, class B>
Json to_json(const std::pair<A, B>& value, TagPolicy policy) {
    Json arr = Json::array();
    arr.arr.reserve(2);
    arr.arr.push_back(to_json(value.first, policy));
    arr.arr.push_back(to_json(value.second, policy));
    return arr;
}

template <class A, class B>
Maybe_error<void> from_json(const Json& j, std::pair<A, B>& out, const std::string& path,
                            TagPolicy policy) {
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    if (j.arr.size() != 2) {
        return error_message(path + ": expected array of length 2");
    }
    auto first_status = from_json(j.arr[0], out.first, path + "[0]", policy);
    if (!first_status) {
        return first_status;
    }
    return from_json(j.arr[1], out.second, path + "[1]", policy);
}

template <class... Ts>
Json to_json(const std::tuple<Ts...>& value, TagPolicy policy) {
    Json arr = Json::array();
    arr.arr.reserve(sizeof...(Ts));
    detail::tuple_push_json(arr, value, policy, std::index_sequence_for<Ts...>{});
    return arr;
}

template <class... Ts>
Maybe_error<void> from_json(const Json& j, std::tuple<Ts...>& out, const std::string& path,
                            TagPolicy policy) {
    if (j.type != Json::Type::Array) {
        return detail::type_error(path, "array", j.type);
    }
    if (j.arr.size() != sizeof...(Ts)) {
        return error_message(path + ": expected array of length " +
                             std::to_string(sizeof...(Ts)));
    }
    return detail::tuple_from_json(j, out, path, policy, std::index_sequence_for<Ts...>{});
}

// -----------------------------------------------------------------------------
// Unique ptr codec (pointee-first if available)
// -----------------------------------------------------------------------------
template <class T>
Json to_json(const std::unique_ptr<T>& p, TagPolicy policy) {
    if (!p) {
        return Json::null();
    }
    if constexpr (has_json_codec_v<T>) {
        return to_json(*p, policy);
    } else {
        Json obj = Json::object();
        obj.obj.emplace("type", Json::string(macrodr::dsl::type_name<T>()));
        obj.obj.emplace("value", Json::null());
        return obj;
    }
}

// Specialized model pointer serializer: save by id and binary
inline Json to_json(
    const std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>& m, TagPolicy) {
    Json obj = Json::object();
    obj.obj.emplace("model_id", Json::string(m ? m->model_name() : std::string{}));
    obj.obj.emplace("binary", Json::string(GIT_COMMIT_HASH));
    return obj;
}

inline Maybe_error<void> from_json(
    const Json& j,
    std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>& out,
    const std::string& path,
    TagPolicy) {
    (void)j;
    (void)out;
    return error_message(path + ": deserialization for model pointer requires compiler registry");
}

// no generic from_json for unique_ptr<T>: loader is unsupported; saver uses to_json only

// -----------------------------------------------------------------------------
// Parameters_values codec (serialize names + values)
// -----------------------------------------------------------------------------
inline Json to_json(const var::Parameters_values& pv, TagPolicy policy) {
    (void)policy;
    Json obj = Json::object();
    // Refer to schema by id; values by matrix
    // We compute id in terms of the underlying schema; must match schema codec
    obj.obj.emplace("schema_id", Json::string("")); // placeholder; writers of schemas should fill top-level
    // Since we cannot thread Environment here, store only values; environment_io will gather schemas
    obj.obj.emplace("values", to_json(pv(), TagPolicy::None));
    return obj;
}

inline Json to_json(const var::Parameters_transformed& pt, TagPolicy policy) {
    (void)policy;
    Json obj = Json::object();
    obj.obj.emplace("schema_id", Json::string(""));
    obj.obj.emplace("values", to_json(pt(), TagPolicy::None));
    return obj;
}

inline Maybe_error<void> from_json(const Json& j, var::Parameters_values& out,
                                   const std::string& path, TagPolicy policy) {
    (void)j;
    (void)out;
    (void)policy;
    return error_message(path + ": deserialization for var::Parameters_values requires model context");
}

inline Maybe_error<void> from_json(const Json& j, var::Parameters_transformed& out,
                                   const std::string& path, TagPolicy policy) {
    (void)j;
    (void)out;
    (void)policy;
    return error_message(path + ": deserialization for var::Parameters_transformed requires model context");
}

// Transformations vector: array of tags
inline Json to_json(const var::transformations_vector& tr, TagPolicy) {
    Json arr = Json::array();
    for (std::size_t i = 0; i < tr.size(); ++i) arr.arr.push_back(Json::string(tr[i]->to_string()));
    return arr;
}

inline Maybe_error<void> from_json(const Json& j, var::transformations_vector& out,
                                   const std::string& path, TagPolicy) {
    if (j.type != Json::Type::Array) return detail::type_error(path, "array", j.type);
    std::vector<std::unique_ptr<var::base_transformation>> vec;
    vec.reserve(j.arr.size());
    for (std::size_t i = 0; i < j.arr.size(); ++i) {
        const auto& e = j.arr[i];
        if (e.type != Json::Type::String)
            return detail::type_error(path + "[" + std::to_string(i) + "]", "string", e.type);
        auto mt = var::MyTranformations::from_string(e.str);
        if (!mt) return error_message(path + "[" + std::to_string(i) + "]: " + mt.error()());
        vec.emplace_back(std::move(mt.value()));
    }
    out = var::transformations_vector(std::move(vec));
    return {};
}

// Schema id helper
inline std::string schema_id_of(const var::Parameters_Transformations& pt) {
    uint64_t h = 1469598103934665603ull; // FNV offset basis
    auto mix = [&](const void* p, std::size_t n) {
        const unsigned char* d = static_cast<const unsigned char*>(p);
        const uint64_t prime = 1099511628211ull;
        for (std::size_t i = 0; i < n; ++i) {
            h ^= static_cast<uint64_t>(d[i]);
            h *= prime;
        }
    };
    auto mix_str = [&](const std::string& s) { mix(s.data(), s.size()); };
    mix_str(pt.IdName());
    for (auto& n : pt.names()) mix_str(n);
    for (std::size_t i = 0; i < pt.transf().size(); ++i) mix_str(pt.transf()[i]->to_string());
    const auto& m = pt.standard_values();
    for (std::size_t i = 0; i < m.size(); ++i) mix(&m[i], sizeof(double));
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return std::string(buf);
}

inline Json to_json(const var::Parameters_Transformations& pt, TagPolicy) {
    Json obj = Json::object();
    obj.obj.emplace("id", Json::string(schema_id_of(pt)));
    obj.obj.emplace("model_name", Json::string(pt.IdName()));
    obj.obj.emplace("names", to_json(pt.names(), TagPolicy::None));
    obj.obj.emplace("transformations", to_json(pt.transf(), TagPolicy::None));
    obj.obj.emplace("standard_values", to_json(pt.standard_values(), TagPolicy::None));
    return obj;
}

inline Maybe_error<void> from_json(const Json& j, var::Parameters_Transformations& out,
                                   const std::string& path, TagPolicy) {
    if (j.type != Json::Type::Object) return detail::type_error(path, "object", j.type);
    const Json* model = j.find("model_name");
    const Json* names = j.find("names");
    const Json* trs = j.find("transformations");
    const Json* vals = j.find("standard_values");
    if (!model || !names || !trs || !vals) return error_message(path + ": missing fields");
    std::string model_name;
    auto s1 = from_json(*model, model_name, path + ".model_name", TagPolicy::None);
    if (!s1) return s1;
    std::vector<std::string> nn;
    auto s2 = from_json(*names, nn, path + ".names", TagPolicy::None);
    if (!s2) return s2;
    var::transformations_vector tv;
    auto s3 = from_json(*trs, tv, path + ".transformations", TagPolicy::None);
    if (!s3) return s3;
    Matrix<double> mv;
    auto s4 = from_json(*vals, mv, path + ".standard_values", TagPolicy::None);
    if (!s4) return s4;
    out = var::Parameters_Transformations(model_name, nn, std::move(tv), mv);
    return {};
}

template <class VarType>
Json to_json(const VarType& value, TagPolicy policy)
    requires(detail::is_var<VarType>::value && !detail::is_vector_space_v<VarType>) {
    Json obj = Json::object();
    obj.obj.emplace("id", Json::string(detail::var_id_name<VarType>()));
    obj.obj.emplace("value", to_json(detail::var_value(value), policy));
    if (policy != TagPolicy::None) {
        obj.obj.emplace("type", Json::string(detail::value_type_name<VarType>()));
    }
    return obj;
}

template <class VarType>
Maybe_error<void> from_json(const Json& j, VarType& out, const std::string& path,
                            TagPolicy policy)
    requires(detail::is_var<VarType>::value && !detail::is_vector_space_v<VarType>) {
    if (j.type != Json::Type::Object) {
        return detail::type_error(path, "object", j.type);
    }
    const auto* id_json = j.find("id");
    const auto* value_json = j.find("value");
    if (!id_json || !value_json) {
        return error_message(path + ": missing 'id' or 'value'");
    }
    std::string id;
    auto id_status = from_json(*id_json, id, path + ".id", policy);
    if (!id_status) {
        return id_status;
    }
    auto expected = detail::var_id_name<VarType>();
    if (id != expected) {
        return error_message(path + ": unexpected id '" + id + "' (expected '" + expected + "')");
    }
    auto status = from_json(*value_json, detail::var_value(out), path + ".value", policy);
    if (!status) {
        return status;
    }
    return {};
}

template <class... Vars>
Json to_json(const var::Vector_Space<Vars...>& value, TagPolicy policy) {
    Json obj = Json::object();
    Json vars_array = Json::array();
    (vars_array.arr.push_back(to_json(static_cast<const Vars&>(value), policy)), ...);
    obj.obj.emplace("vars", std::move(vars_array));
    if (policy != TagPolicy::None) {
        obj.obj.emplace("type", Json::string(macrodr::dsl::type_name<var::Vector_Space<Vars...>>()));
    }
    return obj;
}

template <class... Vars>
Maybe_error<void> from_json(const Json& j, var::Vector_Space<Vars...>& out,
                            const std::string& path, TagPolicy policy) {
    if (j.type != Json::Type::Object) {
        return detail::type_error(path, "object", j.type);
    }
    const auto* vars_json = j.find("vars");
    if (!vars_json) {
        return error_message(path + ": missing 'vars'");
    }
    if (vars_json->type != Json::Type::Array) {
        return detail::type_error(path + ".vars", "array", vars_json->type);
    }
    std::array<bool, sizeof...(Vars)> seen{};
    for (std::size_t i = 0; i < vars_json->arr.size(); ++i) {
        const auto& entry = vars_json->arr[i];
        if (entry.type != Json::Type::Object) {
            return detail::type_error(path + ".vars[" + std::to_string(i) + "]", "object",
                                      entry.type);
        }
        const auto* id_json = entry.find("id");
        const auto* value_json = entry.find("value");
        if (!id_json || !value_json) {
            return error_message(path + ".vars[" + std::to_string(i) + "] missing 'id' or 'value'");
        }
        std::string id;
        auto id_status = from_json(*id_json, id, path + ".vars[" + std::to_string(i) + "].id", policy);
        if (!id_status) {
            return id_status;
        }
        auto status = detail::assign_vector_space_entry<0, var::Vector_Space<Vars...>, Vars...>(
            *value_json, id, out, seen, path + ".vars[" + std::to_string(i) + "]", policy);
        if (!status) {
            return status;
        }
    }
    return detail::ensure_all_seen<0, Vars...>(seen, path);
}

// Convenience overloads with default path/policy
template <class T>
Maybe_error<void> from_json(const Json& j, T& out, const std::string& path) {
    return from_json(j, out, path, TagPolicy::None);
}

// Use '$' to denote the root when no explicit path is supplied
template <class T>
Maybe_error<void> from_json(const Json& j, T& out) {
    return from_json(j, out, std::string("$"), TagPolicy::None);
}

// -----------------------------------------------------------------------------
// detail helper implementations
// -----------------------------------------------------------------------------

namespace detail {

template <class Tuple, std::size_t... Is>
void tuple_push_json(Json& arr, const Tuple& tu, TagPolicy policy,
                     std::index_sequence<Is...>) {
    (arr.arr.push_back(to_json(std::get<Is>(tu), policy)), ...);
}

template <std::size_t I, class Tuple>
Maybe_error<void> tuple_from_json_index(const Json& arr_json, Tuple& tu, const std::string& path,
                                        TagPolicy policy) {
    if constexpr (I == std::tuple_size_v<Tuple>) {
        (void)arr_json;
        (void)tu;
        (void)path;
        (void)policy;
        return {};
    } else {
        auto status = from_json(arr_json.arr[I], std::get<I>(tu),
                                 path + "[" + std::to_string(I) + "]", policy);
        if (!status) {
            return status;
        }
        return tuple_from_json_index<I + 1>(arr_json, tu, path, policy);
    }
}

template <class Tuple, std::size_t... Is>
Maybe_error<void> tuple_from_json(const Json& arr_json, Tuple& tu, const std::string& path,
                                  TagPolicy policy, std::index_sequence<Is...>) {
    (void)sizeof...(Is);
    return tuple_from_json_index<0>(arr_json, tu, path, policy);
}

template <std::size_t I, class VS, class... Vars>
Maybe_error<void> assign_vector_space_entry(const Json& value, const std::string& id,
                                            VS& vs, std::array<bool, sizeof...(Vars)>& seen,
                                            const std::string& path, TagPolicy policy) {
    if constexpr (I == sizeof...(Vars)) {
        return error_message(path + ": unknown id '" + id + "'");
    } else {
        using VarType = std::tuple_element_t<I, std::tuple<Vars...>>;
        auto expected = var_id_name<VarType>();
        if (id == expected) {
            if (seen[I]) {
                return error_message(path + ": duplicate id '" + id + "'");
            }
            seen[I] = true;
            auto& component = static_cast<VarType&>(vs);
            return from_json(value, component(), path + ".value", policy);
        }
        return assign_vector_space_entry<I + 1, VS, Vars...>(value, id, vs, seen, path, policy);
    }
}

template <std::size_t I, class... Vars>
Maybe_error<void> ensure_all_seen(const std::array<bool, sizeof...(Vars)>& seen,
                                  const std::string& path) {
    if constexpr (I == sizeof...(Vars)) {
        (void)seen;
        (void)path;
        return {};
    } else {
        if (!seen[I]) {
            using VarType = std::tuple_element_t<I, std::tuple<Vars...>>;
            return error_message(path + ": missing entry for '" + var_id_name<VarType>() + "'");
        }
        return ensure_all_seen<I + 1, Vars...>(seen, path);
    }
}

}  // namespace detail

}  // namespace macrodr::io::json::conv
