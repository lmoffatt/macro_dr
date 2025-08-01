#pragma once
#include <macrodr/core/fp/map.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

#include "maybe_error.h"
namespace macrodr::core::io::csv {

using DateTime = std::chrono::sys_time<std::chrono::milliseconds>;

using CsvColumn = std::variant<std::vector<bool>, std::vector<int64_t>,
                               std::vector<double>,      // canonical scalar types
                               std::vector<DateTime>,    // timestamp
                               std::vector<std::string>  // fallback and categorical
                               >;

using CsvScalar = std::variant<std::monostate,         // null or missing
                               bool, int64_t, double,  // canonical scalar types
                               DateTime,               // timestamp
                               std::string             // fallback and categorical
                               >;
struct BoolParseConfig {
    std::unordered_set<std::string> true_values = {"true"};
    std::unordered_set<std::string> false_values = {"false"};
    std::unordered_set<std::string> null_values = {"", "nan", "null"};
};

struct ParseOptions {
    char separator = ',';
    BoolParseConfig bools;
    std::vector<std::string> datetime_formats = {"%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S",
                                                 "%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%d %H:%M:%S.%f"};
    std::unordered_map<std::string, std::string> overrides;  // col → "int"|"string"|...
};

struct DataFrame {
    std::vector<std::string> titles;
    std::vector<CsvColumn> columns;
    auto ncols() const {
        return titles.size();
    }
    std::size_t nrows() const {
        if (ncols() > 0)
            return std::visit([](auto& col) { return col.size(); }, columns[0]);
        else
            return 0ul;
    }
};

struct RawCSV {
    std::vector<std::string> titles;
    std::vector<std::vector<std::string>> lines;
};


template <class R>
std::optional<std::vector<R>> col_to_type_strict(const std::vector<std::string>& col,
                                                 const ParseOptions& opt);

template <class R>
std::pair<std::vector<R>, error_message> col_to_type_full(const std::vector<std::string>& col,
                                                          const ParseOptions& opt);

Maybe_error<RawCSV> read_lines(const std::string& filename, char sep = ',');
Maybe_error<std::vector<std::vector<std::string>>> transpose(
    std::vector<std::vector<std::string>> const& lines);
std::string to_lower(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

Maybe_error<CsvColumn> col_to_variant_full(const std::vector<std::string>& col,
                                           const std::string& colname, const ParseOptions& opt);

template <typename... Rs>
std::variant<std::vector<Rs>...> col_to_variant_strict(std::variant<std::monostate, Rs...>,
                                                       const std::vector<std::string>& col,
                                                       const ParseOptions& opt);

std::ostream& print(std::ostream& ss, const DataFrame& df, const ParseOptions& options = {}) {
    ss << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
    auto ncols = df.ncols();
    if (ncols == 0)
        return ss;
    auto nrows = df.nrows();
    for (std::size_t j = 0; j < ncols; ++j) ss << df.titles[j] << (j + 1 < ncols ? "," : "");
    ss << "\n";
    for (std::size_t i = 0; i < nrows; ++i) {
        for (std::size_t j = 0; j < ncols; ++j)
            std::visit(
                [i, j, ncols, &ss](auto&& col) {
                    auto&& e = col[i];
                    using T = std::decay_t<decltype(e)>;
                    if constexpr (std::is_same_v<T, std::string>) {
                        if (e.find(',') != std::string::npos)
                            ss << "\"" << e << "\"";
                        else
                            ss << e;
                    } else if constexpr (std::is_same_v<T, bool>) {
                        ss << (e ? "true" : "false");
                    } else if constexpr (std::is_same_v<T, DateTime>) {
                        ss << std::format("{:%Y-%m-%dT%H:%M:%S}", e);  // requires C++20
                    } else {
                        ss << e;
                    }
                    if (j + 1 < ncols)
                        ss << ",";
                },
                df.columns[j]);
        ss << "\n";
    }
    return ss;
}

std::string to_string(const DataFrame& df, const ParseOptions& options = {}) {
    std::stringstream ss;
    print(ss, df, options);
    return ss.str();
}
Maybe_error<std::string> write_dataframe(const std::string& filename, const DataFrame& df,
                                         const ParseOptions opt = {}) {
    std::ofstream of(filename);
    if (!of)
        return error_message(filename + "could not be opened");
    else
        print(of, df, opt);
    return std::to_string(df.nrows()) + " rows written to " + filename;
}

Maybe_error<DataFrame> load_dataframe(const std::string filename, const ParseOptions opt = {}) {
    auto maybe_Raw = read_lines(filename);
    if (!maybe_Raw)
        return maybe_Raw.error();
    auto raw = std::move(maybe_Raw.value());

    auto maybe_col = transpose(raw.lines);
    if (!maybe_col)
        return maybe_col.error();
    auto cols = std::move(maybe_col.value());

    auto colnames = raw.titles;

    std::vector<CsvColumn> datacols;
    auto ncols = raw.titles.size();
    for (std::size_t i = 0; i < ncols; ++i) {
        auto& colname = colnames[i];
        if (opt.overrides.contains(colname)) {
            auto coltype = opt.overrides.at(colname);
            auto Maybe_col = col_to_variant_full(cols[i], coltype, opt);
            if (!Maybe_col)
                return Maybe_col.error();
            else
                datacols[i] = std::move(Maybe_col.value());
        } else {
            datacols[i] = col_to_variant_strict(CsvScalar{}, cols[i], opt);
        }
    }

    return DataFrame{std::move(colnames), std::move(datacols)};
}

Maybe_error<CsvColumn> col_to_variant_full(const std::vector<std::string>& col,
                                           const std::string& coltype, const ParseOptions& opt) {
    if (to_lower(coltype) == "bool") {
        auto [col_t, error] = col_to_type_full<bool>(col, opt);
        if (error().size() > 0)
            return error;
        else
            return std::move(col_t);
    } else if (to_lower(coltype) == "datetime") {
        auto [col_t, error] = col_to_type_full<DateTime>(col, opt);
        if (error().size() > 0)
            return error;
        else
            return std::move(col_t);

    } else if (to_lower(coltype) == "int") {
        auto [col_t, error] = col_to_type_full<int64_t>(col, opt);
        if (error().size() > 0)
            return error;
        else
            return std::move(col_t);

    } else if (to_lower(coltype) == "double") {
        auto [col_t, error] = col_to_type_full<double>(col, opt);
        if (error().size() > 0)
            return error;
        else
            return std::move(col_t);
    } else {
        return col;
    }
}

std::string to_string(const CsvScalar& x) {
    return std::visit(
        [](auto&& e) -> std::string {
            using T = std::decay_t<decltype(e)>;
            if constexpr (std::is_same_v<T, std::monostate>) {
                return "";
            } else if constexpr (std::is_same_v<T, std::string>) {
                return e;
            } else if constexpr (std::is_same_v<T, DateTime>) {
                return std::format("{:%Y-%m-%dT%H:%M:%S}", e);  // requires C++20
            } else {
                return std::to_string(e);
            }
        },
        x);
}

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t");
    auto end = s.find_last_not_of(" \t");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

std::vector<std::string> parse_csv_line(const std::string& line, char sep = ',') {
    std::vector<std::string> cells;
    std::string cell;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];

        if (in_quotes) {
            if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') {
                cell.push_back('"');  // escaped quote
                ++i;
            } else if (c == '"') {
                in_quotes = false;
            } else {
                cell.push_back(c);
            }
        } else {
            if (c == '"') {
                in_quotes = true;
            } else if (c == sep) {
                cells.push_back(trim(cell));
                cell.clear();
            } else {
                cell.push_back(c);
            }
        }
    }
    cells.push_back(trim(cell));  // last cell
    return cells;
}

Maybe_error<RawCSV> read_lines(const std::string& filename, char sep) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        return error_message("Could not open file: " + filename);
    }

    std::string line;
    RawCSV result;

    // Read header
    if (!std::getline(f, line)) {
        return error_message("Empty file or failed to read header line");
    }

    std::istringstream header_stream(line);
    std::string cell;
    while (std::getline(header_stream, cell, sep)) {
        result.titles.push_back(trim(cell));
    }

    // Read data lines
    while (std::getline(f, line)) {
        auto row = parse_csv_line(line, sep);
        result.lines.push_back(std::move(row));
    }
    return result;
}

template <class R>
std::optional<std::vector<R>> col_to_type_strict(const std::vector<std::string>& col,
                                                 const ParseOptions& opt) {
    std::vector<R> out;
    out.reserve(col.size());

    for (std::size_t i = 0; i < col.size(); ++i) {
        const auto& s = col[i];
        auto lx = to_lower(s);

        if constexpr (std::is_same_v<R, std::string>) {
            out.push_back(s);
        }

        else if constexpr (std::is_same_v<R, int64_t>) {
            int64_t value = 0;
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            if (ec == std::errc() && ptr == s.data() + s.size())
                out.push_back(value);
            else {
                return std::nullopt;
            }
        }

        else if constexpr (std::is_same_v<R, double>) {
            if (opt.bools.null_values.contains(lx)) {
                out.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                double value = 0.0;
                auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
                if (ec == std::errc() && ptr == s.data() + s.size())
                    out.push_back(value);
                else {
                    return std::nullopt;
                    ;
                }
            }
        }

        else if constexpr (std::is_same_v<R, bool>) {
            if (opt.bools.false_values.contains(lx))
                out.push_back(false);
            else if (opt.bools.true_values.contains(lx))
                out.push_back(true);
            else {
                return std::nullopt;
            }
        }

        else if constexpr (std::is_same_v<R, DateTime>) {
            bool parsed = false;
            for (const std::string& fmt : opt.datetime_formats) {
                std::istringstream in(s);
                DateTime tp;
                in >> std::chrono::parse(fmt, tp);
                if (!in.fail()) {
                    out.push_back(tp);
                    parsed = true;
                    break;
                }
            }
            if (!parsed) {
                return std::nullopt;
            }
        }
    }
    return out;
}

template <class... Ts>
struct Classes {};

template <typename... Rs, typename T, typename... Ts>
std::variant<std::vector<Rs>...> col_to_variant_strict(Classes<std::string>,
                                                       const std::vector<std::string>& col,
                                                       const ParseOptions&) {
    return col;
}
template <typename... Rs, typename T, typename... Ts>
std::variant<std::vector<Rs>...> col_to_variant_strict(Classes<T, Ts...>,
                                                       const std::vector<std::string>& col,
                                                       const ParseOptions& opt) {
    auto o = col_to_type_strict<T>(col, opt);
    if (o)
        return std::move(*o);
    else
        return col_to_variant_strict<Rs...>(Classes<Ts...>{}, col, opt);
}

template <typename... Rs>
std::variant<std::vector<Rs>...> col_to_variant(std::variant<std::monostate, Rs...>,
                                                const std::vector<std::string>& col,
                                                const ParseOptions& opt) {
    return col_to_variant_strict<Rs...>(Classes<Rs...>{}, col, opt);
}

template <class R>
std::pair<std::vector<R>, error_message> col_to_type_full(const std::vector<std::string>& col,
                                                          const ParseOptions& opt) {
    std::vector<R> out;
    out.reserve(col.size());
    std::string error;

    for (std::size_t i = 0; i < col.size(); ++i) {
        const auto& s = col[i];
        auto lx = to_lower(s);

        if constexpr (std::is_same_v<R, std::string>) {
            out.push_back(s);
        }

        else if constexpr (std::is_same_v<R, int64_t>) {
            int64_t value = 0;
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            if (ec == std::errc() && ptr == s.data() + s.size())
                out.push_back(value);
            else {
                error += "row " + std::to_string(i) + ": '" + s + "' is not an int\n";
                out.push_back(0);
            }
        }

        else if constexpr (std::is_same_v<R, double>) {
            if (opt.bools.null_values.contains(lx)) {
                out.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                double value = 0.0;
                auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
                if (ec == std::errc() && ptr == s.data() + s.size())
                    out.push_back(value);
                else {
                    error += "row " + std::to_string(i) + ": '" + s + "' is not a double\n";
                    out.push_back(std::numeric_limits<double>::quiet_NaN());
                }
            }
        }

        else if constexpr (std::is_same_v<R, bool>) {
            if (opt.bools.false_values.contains(lx))
                out.push_back(false);
            else if (opt.bools.true_values.contains(lx))
                out.push_back(true);
            else {
                error += "row " + std::to_string(i) + ": '" + s + "' is not a boolean\n";
                out.push_back(false);
            }
        }

        else if constexpr (std::is_same_v<R, DateTime>) {
            bool parsed = false;
            for (const std::string& fmt : opt.datetime_formats) {
                std::istringstream in(s);
                DateTime tp;
                in >> std::chrono::parse(fmt, tp);
                if (!in.fail()) {
                    out.push_back(tp);
                    parsed = true;
                    break;
                }
            }
            if (!parsed) {
                error += "row " + std::to_string(i) + ": '" + s + "' is not a datetime\n";
                out.push_back(DateTime{});
            }
        }
    }
    return {out, error_message(error)};
}

Maybe_error<std::vector<std::vector<std::string>>> transpose(
    std::vector<std::vector<std::string>> const& lines) {
    if (lines.size() == 0)
        return lines;
    auto nrows = lines.size();
    auto ncols = lines.front().size();
    std::string err;
    for (std::size_t i = 0; i < lines.size(); ++i)
        if (lines[i].size() != ncols) {
            err += "at " + std::to_string(i) + " row, " + std::to_string(lines[i].size()) +
                   " columns \n";
            if (lines[i].size() > 0)
                err += lines[i][0];
            for (std::size_t j = 1; j < lines[i].size(); ++j) err += ", " + lines[i][j];
            err += "\n";
        }
    if (err.size() > 0)
        return error_message(err);
    std::vector<std::vector<std::string>> out(ncols, std::vector<std::string>(lines.size()));
    for (std::size_t j = 0; j < ncols; ++j)
        for (std::size_t i = 0; i < nrows; ++i) out[j][i] = lines[i][j];
    return out;
}
namespace zombie {
CsvScalar parse_cell(const std::string& s, const ParseOptions& opt) {
    auto lx = to_lower(s);
    if (opt.bools.null_values.contains(lx))
        return std::monostate{};
    else if (opt.bools.false_values.contains(lx))
        return false;
    else if (opt.bools.true_values.contains(lx))
        return true;
    {
        int64_t value = 0;
        const char* begin = s.data();
        const char* end = s.data() + s.size();
        auto [ptr, ec] = std::from_chars(begin, end, value);
        if (ec == std::errc() && ptr == end)
            return value;
    }
    {
        double value = 0.0;
        const char* begin = s.data();
        const char* end = s.data() + s.size();
        auto [ptr, ec] = std::from_chars(begin, end, value);
        if (ec == std::errc() && ptr == end)
            return value;
    }
    for (const std::string& fmt : opt.datetime_formats) {
        std::istringstream in(s);
        DateTime tp;
        in >> std::chrono::parse(fmt, tp);
        if (!in.fail())
            return tp;
    }
    return s;
}
struct TypedCSV {
    std::vector<std::string> titles;
    std::vector<std::vector<CsvScalar>> lines;
};

std::vector<CsvScalar> cell_typed(const std::vector<std::string>& raw,
                                  const ParseOptions& options) {
    return core::fp::map_to_vector(raw, [&options](auto& s) { return parse_cell(s, options); });
}

template <class R>
Maybe_error<std::vector<R>> col_to_type(const std::vector<CsvScalar>& col) {
    std::vector<R> out;
    out.reserve(col.size());
    std::string error;

    for (std::size_t i = 0; i < col.size(); ++i) {
        const auto& cell = col[i];
        bool match = std::visit(
            [&out, i](auto&& e) -> bool {
                using T = std::decay_t<decltype(e)>;
                if constexpr (std::is_same_v<T, R>) {
                    out[i] = e;
                    return true;
                } else if constexpr (std::is_same_v<R, double> && std::is_same_v<T, int64_t>) {
                    out[i] = static_cast<R>(e);
                    return true;
                } else if constexpr (std::is_same_v<R, double> &&
                                     std::is_same_v<T, std::monostate>) {
                    out[i] = std::numeric_limits<R>::quiet_NaN();
                    return true;
                } else {
                    return false;
                }
            },
            cell);
        if (!match)
            error += "at " + std::to_string(i) + " row, " + to_string(cell) + " is not int\n";
    }
    if (error.size() > 0)
        return error_message(error);
    return out;
}
}  // namespace zombie

}  // namespace macrodr::core::io::csv
