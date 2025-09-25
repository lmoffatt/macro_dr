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

namespace macrodr::io::json::conv {

enum class TagPolicy { None, Minimal, Full };

namespace detail {

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

// var::Var and Vector_Space
template <class VarType>
Json to_json(const VarType& value, TagPolicy policy = TagPolicy::None)
    requires detail::is_var<VarType>::value;

template <class VarType>
Maybe_error<void> from_json(const Json& j, VarType& out, const std::string& path,
                            TagPolicy policy)
    requires detail::is_var<VarType>::value;

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
    arr.arr.reserve(values.size());
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
    out.reserve(j.arr.size());
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
    rows.arr.reserve(m.nrows());
    for (std::size_t i = 0; i < m.nrows(); ++i) {
        Json row = Json::array();
        row.arr.reserve(m.ncols());
        for (std::size_t jcol = 0; jcol < m.ncols(); ++jcol) {
            row.arr.push_back(Json::number(m(i, jcol)));
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

inline Json to_json(const DiagonalMatrix<double>& m, TagPolicy) {
    Json obj = Json::object();
    obj.obj.emplace("rows", Json::number(static_cast<double>(m.nrows())));
    obj.obj.emplace("cols", Json::number(static_cast<double>(m.ncols())));
    Json diag = Json::array();
    diag.arr.reserve(m.size());
    for (std::size_t i = 0; i < m.size(); ++i) {
        diag.arr.push_back(Json::number(m[i]));
    }
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

template <class VarType>
Json to_json(const VarType& value, TagPolicy policy)
    requires detail::is_var<VarType>::value {
    Json obj = Json::object();
    obj.obj.emplace("id", Json::string(detail::var_id_name<VarType>()));
    obj.obj.emplace("value", to_json(value(), policy));
    if (policy != TagPolicy::None) {
        obj.obj.emplace("type", Json::string(detail::value_type_name<VarType>()));
    }
    return obj;
}

template <class VarType>
Maybe_error<void> from_json(const Json& j, VarType& out, const std::string& path,
                            TagPolicy policy)
    requires detail::is_var<VarType>::value {
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
    auto status = from_json(*value_json, out(), path + ".value", policy);
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
