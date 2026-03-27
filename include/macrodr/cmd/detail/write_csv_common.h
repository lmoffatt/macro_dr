#pragma once

#include <derivative_operator.h>
#include <macrodr/dsl/type_name.h>

#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "maybe_error.h"
#include "moment_statistics.h"
#include "parameter_indexed.h"
#include "parameters_derivative.h"
#include "qmodel.h"

namespace macrodr::cmd::detail {

template <class VS>
struct vector_space_types;

template <class... Cs>
struct vector_space_types<var::Vector_Space<Cs...>> {
    using types = std::tuple<std::type_identity<Cs>...>;
};

template <class T>
struct is_vector_space : std::false_type {};

template <class... Cs>
struct is_vector_space<var::Vector_Space<Cs...>> : std::true_type {};

template <class T>
inline constexpr bool is_vector_space_v = is_vector_space<std::remove_cvref_t<T>>::value;

template <class T>
struct is_moment_statistics : std::false_type {};

template <class Id, bool IncludeCovariance>
struct is_moment_statistics<Moment_statistics<Id, IncludeCovariance>> : std::true_type {};

template <class T>
inline constexpr bool is_moment_statistics_v =
    is_moment_statistics<std::remove_cvref_t<T>>::value;

template <class T>
struct moment_statistics_traits;

template <class Id, bool IncludeCovariance>
struct moment_statistics_traits<Moment_statistics<Id, IncludeCovariance>> {
    using id = Id;
    using spread_type = variance_t<Id, IncludeCovariance>;
    static constexpr bool include_covariance = IncludeCovariance;
    static constexpr bool componentwise = is_vector_space_v<Id>;
};

template <class T>
struct is_probit_statistics : std::false_type {};

template <class Id>
struct is_probit_statistics<Probit_statistics<Id>> : std::true_type {};

template <class T>
inline constexpr bool is_probit_statistics_v =
    is_probit_statistics<std::remove_cvref_t<T>>::value;

template <class T>
struct probit_statistics_traits;

template <class Id>
struct probit_statistics_traits<Probit_statistics<Id>> {
    using id = Id;
};

template <class T>
struct is_evolution_of : std::false_type {};

template <class T>
struct is_evolution_of<Evolution_of<T>> : std::true_type {};

template <class T>
inline constexpr bool is_evolution_of_v = is_evolution_of<std::remove_cvref_t<T>>::value;

template <class T>
std::string sanitize_component_segment(T&& raw) {
    std::string in(std::forward<T>(raw));
    std::string out;
    out.reserve(in.size());
    bool last_was_underscore = false;
    for (unsigned char ch : in) {
        if (std::isalnum(ch) != 0) {
            out.push_back(static_cast<char>(ch));
            last_was_underscore = false;
        } else if (!last_was_underscore) {
            out.push_back('_');
            last_was_underscore = true;
        }
    }
    while (!out.empty() && out.front() == '_') {
        out.erase(out.begin());
    }
    while (!out.empty() && out.back() == '_') {
        out.pop_back();
    }
    return out.empty() ? std::string("value") : out;
}

template <class T>
std::string component_label() {
    return sanitize_component_segment(
        macrodr::dsl::type_name_no_namespace<var::untransformed_type_t<std::remove_cvref_t<T>>>());
}

inline std::string append_component_path(const std::string& base, const std::string& segment) {
    return base.empty() ? segment : base + "." + segment;
}

struct CsvContext {
    std::string scope = "summary";
    std::optional<std::size_t> simulation_index;
    std::optional<std::size_t> sample_index;
    std::optional<std::size_t> segment_index;
    std::optional<std::size_t> sub_index;
    std::optional<double> n_step;
    std::optional<double> time_start;
    std::optional<double> time_end;
    std::optional<double> time_middle;
    std::optional<double> agonist;
    std::optional<double> patch_current;
    std::string component_path;
    std::optional<std::size_t> value_row;
    std::optional<std::size_t> value_col;
    std::string probit = "point";
    std::string statistic = "value";
    std::optional<double> quantile_level;
    std::optional<std::size_t> param_index;
    std::optional<std::size_t> param_col;
    std::optional<std::string> param_name;
};

class CsvWriter {
  public:
    CsvWriter(std::ostream& out, const std::vector<std::string>& param_names)
        : out_(out), param_names_(param_names) {
        out_ << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        write_header();
    }

    const std::vector<std::string>& param_names() const { return param_names_; }

    void write_row(const CsvContext& ctx, std::string_view calculus, double value) {
        write_string(ctx.scope);
        out_ << ",";
        write_optional(ctx.simulation_index);
        out_ << ",";
        write_optional(ctx.sample_index);
        out_ << ",";
        write_optional(ctx.segment_index);
        out_ << ",";
        write_optional(ctx.sub_index);
        out_ << ",";
        write_optional(ctx.n_step);
        out_ << ",";
        write_optional(ctx.time_start);
        out_ << ",";
        write_optional(ctx.time_end);
        out_ << ",";
        write_optional(ctx.time_middle);
        out_ << ",";
        write_optional(ctx.agonist);
        out_ << ",";
        write_optional(ctx.patch_current);
        out_ << ",";
        write_string(ctx.component_path);
        out_ << ",";
        write_optional(ctx.value_row);
        out_ << ",";
        write_optional(ctx.value_col);
        out_ << ",";
        write_string(ctx.probit);
        out_ << ",";
        write_string(std::string(calculus));
        out_ << ",";
        write_string(ctx.statistic);
        out_ << ",";
        write_optional(ctx.quantile_level);
        out_ << ",";
        write_optional(ctx.param_index);
        out_ << ",";
        write_optional(ctx.param_col);
        out_ << ",";
        if (ctx.param_name.has_value()) {
            write_string(*ctx.param_name);
        }
        out_ << "," << value << "\n";
    }

  private:
    template <class T>
    void write_optional(const std::optional<T>& v) {
        if (v.has_value()) {
            out_ << *v;
        }
    }

    void write_string(const std::string& v) { out_ << v; }

    void write_header() {
        out_ << "scope,simulation_index,sample_index,segment_index,sub_index,n_step,"
                "step_start,step_end,step_middle,agonist,patch_current,"
                "component_path,value_row,value_col,"
                "probit,calculus,statistic,quantile_level,"
                "param_index,param_col,param_name,"
                "value\n";
    }

    std::ostream& out_;
    const std::vector<std::string>& param_names_;
};

template <class T>
std::optional<std::vector<std::string>> find_param_names(const T& x);

inline void assign_vector_param_metadata(CsvContext& ctx, const std::vector<std::string>& names,
                                         std::size_t index) {
    ctx.param_index = index;
    ctx.param_col.reset();
    if (index < names.size()) {
        ctx.param_name = names[index];
    } else {
        ctx.param_name.reset();
    }
}

inline void assign_matrix_param_metadata(CsvContext& ctx, const std::vector<std::string>& names,
                                         std::size_t row, std::size_t col) {
    ctx.param_index = row;
    ctx.param_col = col;
    if (row < names.size() && col < names.size()) {
        ctx.param_name = names[row] + "_" + names[col];
    } else {
        ctx.param_name.reset();
    }
}

template <class Writer, class T>
Maybe_error<bool> emit_any(Writer& w, CsvContext ctx, const T& x);

template <class Writer, class T>
Maybe_error<bool> emit_calculus_value(Writer& w, CsvContext ctx, std::string_view calculus,
                                      const T& x);

template <class Writer, class T>
Maybe_error<bool> emit_named_component(Writer& w, CsvContext ctx, const std::string& label,
                                       const T& x) {
    ctx.component_path = append_component_path(ctx.component_path, sanitize_component_segment(label));
    return emit_any(w, std::move(ctx), x);
}

template <class Writer, class T>
Maybe_error<bool> emit_calculus_named_component(Writer& w, CsvContext ctx,
                                                std::string_view calculus,
                                                const std::string& label, const T& x) {
    ctx.component_path = append_component_path(ctx.component_path, sanitize_component_segment(label));
    return emit_calculus_value(w, std::move(ctx), calculus, x);
}

template <class T>
std::optional<std::vector<std::string>> find_param_names_in_vector_space(const T& x) {
    using Tuple = typename vector_space_types<std::remove_cvref_t<T>>::types;
    std::optional<std::vector<std::string>> found;
    []<std::size_t... Is>(const T& value, std::optional<std::vector<std::string>>& out,
                          std::index_sequence<Is...>) {
        (([&] {
            if (out.has_value()) {
                return;
            }
            using Comp = typename std::tuple_element_t<Is, Tuple>::type;
            out = find_param_names(get<Comp>(value));
        }()),
         ...);
    }(x, found, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
    return found;
}

template <class T>
std::optional<std::vector<std::string>> find_param_names(const T& x) {
    if constexpr (var::is_parameter_indexed_v<T>) {
        auto names = var::parameter_names(x);
        if (!names.empty()) {
            return names;
        }
        return std::nullopt;
    } else if constexpr (var::is_derivative_v<std::remove_cvref_t<T>> &&
                  requires { var::get_dx_of_dfdx(x).parameters().names(); }) {
        const auto& names = var::get_dx_of_dfdx(x).parameters().names();
        return std::vector<std::string>(names.begin(), names.end());
    } else if constexpr (is_probit_statistics_v<T>) {
        using Id = typename probit_statistics_traits<std::remove_cvref_t<T>>::id;
        if (auto found = find_param_names(get<mean<Id>>(x())); found.has_value()) {
            return found;
        }
        const auto& probits = get<Probits<Id>>(x())();
        for (const auto& [level, value] : probits) {
            (void)level;
            if (auto found = find_param_names(value); found.has_value()) {
                return found;
            }
        }
        return std::nullopt;
    } else if constexpr (is_moment_statistics_v<T>) {
        using Traits = moment_statistics_traits<std::remove_cvref_t<T>>;
        if constexpr (Traits::componentwise) {
            using Base = typename std::remove_cvref_t<T>::base_type;
            return find_param_names(static_cast<const Base&>(x));
        } else {
            using Id = typename Traits::id;
            if (auto found = find_param_names(get<mean<Id>>(x())); found.has_value()) {
                return found;
            }
            return find_param_names(get<typename Traits::spread_type>(x()));
        }
    } else if constexpr (is_evolution_of_v<T>) {
        const auto& values = x();
        if (!values.empty()) {
            return find_param_names(values.front());
        }
        return std::nullopt;
    } else if constexpr (is_vector_space_v<T>) {
        return find_param_names_in_vector_space(x);
    } else if constexpr (requires { x(); }) {
        return find_param_names(x());
    } else {
        return std::nullopt;
    }
}

template <class T>
std::vector<std::string> get_param_names_if_any(const T& x) {
    if (auto found = find_param_names(x); found.has_value()) {
        return *found;
    }
    return {};
}

template <class Writer, class T>
Maybe_error<bool> emit_vector_like(Writer& w, CsvContext ctx, const T& x) {
    const auto param_names = find_param_names(x);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto item_ctx = ctx;
        item_ctx.value_row = i;
        item_ctx.value_col.reset();
        if (param_names.has_value()) {
            assign_vector_param_metadata(item_ctx, *param_names, i);
        }
        auto ok = emit_any(w, std::move(item_ctx), x[i]);
        if (!ok || !ok.value()) {
            return ok;
        }
    }
    return true;
}

template <class Writer, class T>
Maybe_error<bool> emit_matrix_like(Writer& w, CsvContext ctx, const T& x) {
    const auto param_names = find_param_names(x);
    if (x.nrows() == 1 || x.ncols() == 1) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            auto item_ctx = ctx;
            item_ctx.value_row = i;
            item_ctx.value_col.reset();
            if (param_names.has_value()) {
                assign_vector_param_metadata(item_ctx, *param_names, i);
            }
            auto ok = emit_any(w, std::move(item_ctx), x[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
        }
        return true;
    }
    for (std::size_t r = 0; r < x.nrows(); ++r) {
        for (std::size_t c = 0; c < x.ncols(); ++c) {
            auto item_ctx = ctx;
            item_ctx.value_row = r;
            item_ctx.value_col = c;
            if (param_names.has_value()) {
                assign_matrix_param_metadata(item_ctx, *param_names, r, c);
            }
            auto ok = emit_any(w, std::move(item_ctx), x(r, c));
            if (!ok || !ok.value()) {
                return ok;
            }
        }
    }
    return true;
}

template <class Writer, class T>
Maybe_error<bool> emit_vector_space(Writer& w, CsvContext ctx, const T& x) {
    using Tuple = typename vector_space_types<std::remove_cvref_t<T>>::types;
    Maybe_error<bool> ok = true;
    []<std::size_t... Is>(Writer& writer, CsvContext base_ctx, const T& value,
                          Maybe_error<bool>& ok_ref, std::index_sequence<Is...>) {
        ((ok_ref =
              ok_ref
                  ? emit_named_component(
                        writer, base_ctx,
                        component_label<typename std::tuple_element_t<Is, Tuple>::type>(),
                        get<typename std::tuple_element_t<Is, Tuple>::type>(value))
                  : ok_ref),
         ...);
    }(w, std::move(ctx), x, ok, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
    return ok;
}

template <class Writer, class T>
Maybe_error<bool> emit_calculus_value(Writer& w, CsvContext ctx, std::string_view calculus,
                                      const T& x) {
    if constexpr (is_vector_space_v<T>) {
        using Tuple = typename vector_space_types<std::remove_cvref_t<T>>::types;
        Maybe_error<bool> ok = true;
        []<std::size_t... Is>(Writer& writer, CsvContext base_ctx, std::string_view calc,
                              const T& value, Maybe_error<bool>& ok_ref,
                              std::index_sequence<Is...>) {
            ((ok_ref =
                  ok_ref
                      ? emit_calculus_named_component(
                            writer, base_ctx, calc,
                            component_label<typename std::tuple_element_t<Is, Tuple>::type>(),
                            get<typename std::tuple_element_t<Is, Tuple>::type>(value))
                      : ok_ref),
             ...);
        }(w, std::move(ctx), calculus, x, ok, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
        return ok;
    } else if constexpr (requires { x(); }) {
        return emit_calculus_value(w, std::move(ctx), calculus, x());
    } else if constexpr (requires {
                             x.nrows();
                             x.ncols();
                             x(0ul, 0ul);
                         }) {
        if (x.nrows() == 1 || x.ncols() == 1) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                auto item_ctx = ctx;
                item_ctx.value_row = i;
                item_ctx.value_col.reset();
                auto ok = emit_calculus_value(w, std::move(item_ctx), calculus, x[i]);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        }
        for (std::size_t r = 0; r < x.nrows(); ++r) {
            for (std::size_t c = 0; c < x.ncols(); ++c) {
                auto item_ctx = ctx;
                item_ctx.value_row = r;
                item_ctx.value_col = c;
                auto ok = emit_calculus_value(w, std::move(item_ctx), calculus, x(r, c));
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
        }
        return true;
    } else if constexpr (requires {
                             x.size();
                             x[0];
                         }) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            auto item_ctx = ctx;
            item_ctx.value_row = i;
            item_ctx.value_col.reset();
            auto ok = emit_calculus_value(w, std::move(item_ctx), calculus, x[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
        }
        return true;
    } else if constexpr (std::convertible_to<T, double>) {
        w.write_row(ctx, calculus, static_cast<double>(x));
        return true;
    } else {
        return error_message("write_csv_rows: unsupported ", calculus, " value type ",
                             ctx.component_path, " (", macrodr::dsl::type_name<T>(), ")");
    }
}

template <class Writer, class Der>
Maybe_error<bool> emit_derivative_object(Writer& w, CsvContext ctx, const Der& d) {
    if constexpr (var::is_derivative_v<std::remove_cvref_t<Der>> &&
                  is_vector_space_v<
                      var::untransformed_type_t<std::remove_cvref_t<Der>>>) {
        using VS = var::untransformed_type_t<std::remove_cvref_t<Der>>;
        using Tuple = typename vector_space_types<VS>::types;
        Maybe_error<bool> ok = true;
        []<std::size_t... Is>(Writer& writer, CsvContext base_ctx, const Der& value,
                              Maybe_error<bool>& ok_ref, std::index_sequence<Is...>) {
            ((ok_ref =
                  ok_ref
                      ? emit_named_component(
                            writer, base_ctx,
                            component_label<typename std::tuple_element_t<Is, Tuple>::type>(),
                            get<typename std::tuple_element_t<Is, Tuple>::type>(value))
                      : ok_ref),
             ...);
        }(w, std::move(ctx), d, ok, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
        return ok;
    } else if constexpr (requires { var::inside_out(d); primitive(d); }) {
        auto primitive_ok = emit_any(w, ctx, primitive(d));
        if (!primitive_ok || !primitive_ok.value()) {
            return primitive_ok;
        }
        auto inside = var::inside_out(d);
        if constexpr (requires {
                          inside.nrows();
                          inside.ncols();
                          inside(0ul, 0ul);
                      }) {
            for (std::size_t r = 0; r < inside.nrows(); ++r) {
                for (std::size_t c = 0; c < inside.ncols(); ++c) {
                    auto item_ctx = ctx;
                    item_ctx.value_row = r;
                    item_ctx.value_col = c;
                    auto ok = emit_derivative_object(w, std::move(item_ctx), inside(r, c));
                    if (!ok || !ok.value()) {
                        return ok;
                    }
                }
            }
            return true;
        } else if constexpr (requires {
                                 inside.size();
                                 inside[0];
                             }) {
            for (std::size_t i = 0; i < inside.size(); ++i) {
                auto item_ctx = ctx;
                item_ctx.value_row = i;
                item_ctx.value_col.reset();
                auto ok = emit_derivative_object(w, std::move(item_ctx), inside[i]);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        } else {
            return error_message("write_csv_rows: unsupported inside_out container for ",
                                 ctx.component_path, " (",
                                 macrodr::dsl::type_name<decltype(inside)>(), ")");
        }
    } else if constexpr (requires { d.derivative()(); primitive(d); }) {
        auto primitive_ok = emit_any(w, ctx, primitive(d));
        if (!primitive_ok || !primitive_ok.value()) {
            return primitive_ok;
        }
        const auto& m = d.derivative()();
        for (std::size_t r = 0; r < m.nrows(); ++r) {
            for (std::size_t c = 0; c < m.ncols(); ++c) {
                auto row_ctx = ctx;
                row_ctx.param_index = r;
                row_ctx.param_col = c;
                if (r < w.param_names().size()) {
                    row_ctx.param_name = w.param_names()[r];
                } else {
                    row_ctx.param_name.reset();
                }
                auto ok = emit_calculus_value(w, std::move(row_ctx), "derivative", m(r, c));
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
        }
        return true;
    } else if constexpr (requires { d(); }) {
        return emit_derivative_object(w, std::move(ctx), d());
    } else {
        return error_message("write_csv_rows: unsupported derivative type ", ctx.component_path,
                             " (", macrodr::dsl::type_name<Der>(), ")");
    }
}

template <class Writer, class T>
Maybe_error<bool> emit_any(Writer& w, CsvContext ctx, const T& x) {
    if constexpr (is_probit_statistics_v<T>) {
        using Id = typename probit_statistics_traits<std::remove_cvref_t<T>>::id;
        auto mean_ctx = ctx;
        mean_ctx.probit = "mean";
        mean_ctx.quantile_level.reset();
        auto ok = emit_any(w, std::move(mean_ctx), get<mean<Id>>(x()));
        if (!ok || !ok.value()) {
            return ok;
        }
        const auto& probits = get<Probits<Id>>(x())();
        for (const auto& [quantile, value] : probits) {
            auto quantile_ctx = ctx;
            quantile_ctx.probit = "quantile";
            quantile_ctx.quantile_level = quantile;
            ok = emit_any(w, std::move(quantile_ctx), value);
            if (!ok || !ok.value()) {
                return ok;
            }
        }
        return true;
    } else if constexpr (is_moment_statistics_v<T>) {
        using Traits = moment_statistics_traits<std::remove_cvref_t<T>>;
        if constexpr (Traits::componentwise) {
            using Base = typename std::remove_cvref_t<T>::base_type;
            return emit_any(w, std::move(ctx), static_cast<const Base&>(x));
        } else {
            using Id = typename Traits::id;
            auto count_ctx = ctx;
            count_ctx.statistic = "count";
            auto ok = emit_any(w, std::move(count_ctx), get<count<Id>>(x()));
            if (!ok || !ok.value()) {
                return ok;
            }
            auto mean_ctx = ctx;
            mean_ctx.statistic = "mean";
            ok = emit_any(w, std::move(mean_ctx), get<mean<Id>>(x()));
            if (!ok || !ok.value()) {
                return ok;
            }
            auto spread_ctx = ctx;
            spread_ctx.statistic = Traits::include_covariance ? "covariance" : "variance";
            return emit_any(w, std::move(spread_ctx), get<typename Traits::spread_type>(x()));
        }
    } else if constexpr (var::is_derivative_v<std::remove_cvref_t<T>>) {
        return emit_derivative_object(w, std::move(ctx), x);
    } else if constexpr (is_evolution_of_v<T>) {
        const auto& values = x();
        for (std::size_t i = 0; i < values.size(); ++i) {
            auto item_ctx = ctx;
            item_ctx.scope = "evolution";
            item_ctx.sample_index = i;
            item_ctx.segment_index.reset();
            item_ctx.sub_index.reset();
            item_ctx.n_step.reset();
            item_ctx.time_start.reset();
            item_ctx.time_end.reset();
            item_ctx.time_middle.reset();
            item_ctx.agonist.reset();
            item_ctx.patch_current.reset();
            auto ok = emit_any(w, std::move(item_ctx), values[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
        }
        return true;
    } else if constexpr (is_vector_space_v<T>) {
        return emit_vector_space(w, std::move(ctx), x);
    } else if constexpr (requires { x(); }) {
        return emit_any(w, std::move(ctx), x());
    } else if constexpr (requires {
                             x.nrows();
                             x.ncols();
                             x(0ul, 0ul);
                         }) {
        return emit_matrix_like(w, std::move(ctx), x);
    } else if constexpr (requires {
                             x.size();
                             x[0];
                         }) {
        return emit_vector_like(w, std::move(ctx), x);
    } else if constexpr (std::convertible_to<T, double>) {
        w.write_row(ctx, "primitive", static_cast<double>(x));
        return true;
    } else {
        return error_message("write_csv_rows: unsupported value type ", ctx.component_path, " (",
                             macrodr::dsl::type_name<T>(), ")");
    }
}

template <class T>
Maybe_error<std::string> write_summary_csv(const T& x, std::string path,
                                           std::string scope = "summary") {
    const auto path_with_extension = path + ".csv";
    std::ofstream f(path_with_extension);
    if (!f.is_open()) {
        return error_message("cannot open ", path_with_extension);
    }

    const auto param_names = get_param_names_if_any(x);
    CsvWriter writer(f, param_names);
    CsvContext ctx;
    ctx.scope = std::move(scope);

    auto ok = emit_any(writer, std::move(ctx), x);
    if (!ok || !ok.value()) {
        return ok.error()();
    }
    return path_with_extension;
}

}  // namespace macrodr::cmd::detail
