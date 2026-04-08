#pragma once
#include<string>
#include <string_view>
#include <maybe_error.h>


namespace macrodr::dsl {

class Lexer {
   public:
    constexpr static const std::string_view whitespace = " \t";
    constexpr static const std::string_view whitespaceline = " \t\n";

    constexpr static const std::string_view alfa =
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";
    constexpr static const std::string_view alfanum =
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";

    constexpr static const std::string_view numerals = "0123456789";
    constexpr static const std::string_view decimal_sep = ".";
    constexpr static const std::string_view numeric_exponential_allowed = "Ee";

    constexpr static const std::string_view numeric_plus = "+";
    constexpr static const std::string_view numeric_minus = "-";
    constexpr static const std::string_view numeric_plus_or_minus = "-+";

    constexpr static const std::string_view argument_list_start = "(";
    constexpr static const std::string_view argument_list_sep = ",";
    constexpr static const std::string_view argument_list_end = ")";

    constexpr static const std::string_view vector_list_start = "[";
    constexpr static const std::string_view vector_list_sep = ",";
    constexpr static const std::string_view vector_list_end = "]";

    constexpr static const std::string_view tuple_list_start = "{";
    constexpr static const std::string_view tuple_list_sep = ",";
    constexpr static const std::string_view tuple_list_end = "}";



    constexpr static const std::string_view statement_sep = "\n";

    constexpr static const std::string_view argument_sep = ",";

    constexpr static const std::string_view assignment_operator = "=";

    constexpr static const std::string_view sum_operator = "+";
    constexpr static const std::string_view product_operator = "*";

    constexpr static const std::string_view substraction_operator = "-";
    constexpr static const std::string_view division_operator = "/";

    constexpr static const std::string_view potentiation_operator = "^";

    constexpr static const std::string_view and_operator = "&";
    constexpr static const std::string_view or_operator = "|";

    constexpr static const std::string_view equality_operator = "==";

    constexpr static const std::string_view greater_operator = ">";

    constexpr static const std::string_view assignment_sep = " ";

    constexpr static const std::string_view label_delimiters = "'\"";

    constexpr static const std::string_view comment_start = "//";

    constexpr static const std::string_view comment_end = "/n";

    static auto skip_comment(const std::string& s, std::size_t pos) {
        if (s.substr(pos, comment_start.size()) == comment_start) {
            return s.find_first_of('\n', pos);
        }             return pos;
    }
    static bool is_end_of_statement(const std::string& s, std::size_t pos) {
        return (s.size() == pos + 1) || (s.substr(pos, statement_sep.size()) == statement_sep);
    }

    static bool is_end_of_argument(const std::string& s, std::size_t pos) {
        return ((s.substr(pos, Lexer::argument_list_sep.size()) == Lexer::argument_list_sep) ||
                (s.substr(pos, Lexer::argument_list_end.size()) == Lexer::argument_list_end));
    }

    static bool is_start_of_argument_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::argument_list_start.size()) == Lexer::argument_list_start);
    }

    static bool is_start_of_vector_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::vector_list_start.size()) == Lexer::vector_list_start);
    }

    static bool is_start_of_tuple_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::tuple_list_start.size()) == Lexer::tuple_list_start);
    }
    static bool is_start_of_comment(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::comment_start.size()) == Lexer::comment_start);
    }

    static bool is_argument_separator(const std::string& s, std::size_t pos) {
        return ((pos != std::string::npos) &&
                (s.substr(pos, Lexer::argument_list_sep.size()) == Lexer::argument_list_sep));
    }

    static bool is_end_of_argument_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::argument_list_end.size()) == Lexer::argument_list_end);
    }

    static bool is_end_of_vector_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::vector_list_end.size()) == Lexer::vector_list_end);
    }

    static bool is_end_of_tuple_list(const std::string& s, std::size_t pos) {
        return (s.substr(pos, Lexer::tuple_list_end.size()) == Lexer::tuple_list_end);
    }



    static auto skip_whitespace(const std::string& s, std::size_t pos) {
        return s.find_first_not_of(whitespace, pos);
    }
    static auto skip_whitespaceline(const std::string& s, std::size_t pos) {
        return s.find_first_not_of(whitespaceline, pos);
    }

    static auto skip_whiteline(const std::string& s, std::size_t pos) {
        return s.find_first_not_of('\n', pos);
    }

    static auto skip_plus_or_minus(const std::string& s, std::size_t pos) {
        return s.find_first_not_of(numeric_plus_or_minus, pos);
    }

    static auto to_end_of_line(const std::string& s, std::size_t pos) {
        return s.substr(pos, s.find('\n', pos) - pos);
    }

    static bool is_numeral(const std::string& s, std::size_t pos) {
        return (numerals.find(s[pos]) != std::string::npos);
    }
    static bool is_Exponent(const std::string& s, std::size_t pos) {
        return (numeric_exponential_allowed.find(s[pos]) != std::string::npos);
    }

    template <class T>
        requires(!is_of_this_template_type_v<T, Maybe_error>)
    static Maybe_error<T> get(const std::string& s) {
        std::stringstream ss(s);
        T x;
        if (ss >> x) {
            return x;
        }             return error_message("extraction fails");
    }

    template <class T>
        requires(is_of_this_template_type_v<T, Maybe_error>)
    static T get(const std::string& s) {
        return get<typename T::value_type>(s);
    }

    template <class T>
    static std::string put(const T& e) {
        std::string const s;
        std::stringstream ss(s);
        ss << e;
        return ss.str();
    }
};
}