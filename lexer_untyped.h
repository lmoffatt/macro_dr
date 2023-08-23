#ifndef LEXER_UNTYPED_H
#define LEXER_UNTYPED_H

#include "grammar_untyped.h"
#include <iostream>

namespace dcli {

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

  constexpr static const std::string_view statement_sep = "\n";

  constexpr static const std::string_view argument_sep = ",";

  constexpr static const std::string_view assignment_operator = "=";

  //  constexpr static const std::string_view sum_operator = "+";
  //  constexpr static const std::string_view product_operator = "*";

  //  constexpr static const std::string_view substraction_operator = "-";
  //  constexpr static const std::string_view division_operator = "/";

  //  constexpr static const std::string_view potentiation_operator = "^";

  //  constexpr static const std::string_view and_operator = "&";
  //  constexpr static const std::string_view or_operator = "|";

  //  constexpr static const std::string_view equality_operator = "==";

  //  constexpr static const std::string_view greater_operator = ">";

  constexpr static const std::string_view assignment_sep = " ";

  constexpr static const std::string_view label_delimiters = "'\"";

  constexpr static const std::string_view comment_start = "//";

  static auto skip_comment(const std::string &s, std::size_t pos) {
    if (s.substr(pos, comment_start.size()) == comment_start)
      return s.find_first_of("\n", pos);
    else
      return pos;
  }
  static bool is_end_of_statement(const std::string &s, std::size_t pos) {

    return (s.size() == pos + 1) ||
           (s.substr(pos, statement_sep.size()) == statement_sep);
  }

  static bool is_end_of_argument(const std::string &s, std::size_t pos) {

    return ((s.substr(pos, Lexer::argument_list_sep.size()) ==
             Lexer::argument_list_sep) ||
            (s.substr(pos, Lexer::argument_list_end.size()) ==
             Lexer::argument_list_end));
  }

  static bool is_start_of_argument_list(const std::string &s, std::size_t pos) {
    return (s.substr(pos, Lexer::argument_list_start.size()) ==
            Lexer::argument_list_start);
  }

  static bool is_argument_separator(const std::string &s, std::size_t pos) {
    return (s.substr(pos, Lexer::argument_list_sep.size()) ==
            Lexer::argument_list_sep);
  }

  static bool is_end_of_argument_list(const std::string &s, std::size_t pos) {
    return (s.substr(pos, Lexer::argument_list_end.size()) ==
            Lexer::argument_list_end);
  }
  static auto skip_whitespace(const std::string &s, std::size_t pos) {
    return s.find_first_not_of(whitespace, pos);
  }
  static auto skip_whitespaceline(const std::string &s, std::size_t pos) {
    return s.find_first_not_of(whitespaceline, pos);
  }

  static auto skip_whiteline(const std::string &s, std::size_t pos) {
    return s.find_first_not_of("\n", pos);
  }

  static auto skip_plus_or_minus(const std::string &s, std::size_t pos) {
    return s.find_first_not_of(numeric_plus_or_minus, pos);
  }

  static auto to_end_of_line(const std::string &s, std::size_t pos) {
    return s.substr(pos, s.find("\n", pos) - pos);
  }

  static bool is_numeral(const std::string &s, std::size_t pos) {
    return (numerals.find(s[pos]) != std::string::npos);
  }
  static bool is_Exponent(const std::string &s, std::size_t pos) {
    return (numeric_exponential_allowed.find(s[pos]) != std::string::npos);
  }
};

// template<class Lexer>

Maybe_error<std::pair<std::unique_ptr<untyped_identifier<Lexer>>, std::size_t>>
extract_identifier(const std::string &s, std::size_t pos) {

  auto pos0 = Lexer::alfa.find(s[pos]);
  if (pos0 == std::string::npos)
    return error_message("");
  else {
    auto posf = s.find_first_not_of(Lexer::alfanum, pos);
    auto id = to_Identifier<Lexer>(s.substr(pos, posf - pos));
    if (!id)
      return id.error();
    else
      return std::pair{std::unique_ptr<untyped_identifier<Lexer>>(
                           new untyped_identifier<Lexer>((*id))),
                       posf};
  }
}

Maybe_error<
    std::pair<std::unique_ptr<untyped_string_value<Lexer>>, std::size_t>>
extract_string(const std::string &s, std::size_t pos) {
  auto ch = std::string(Lexer::label_delimiters);
  auto delimiter_n = Lexer::label_delimiters.find(s[pos]);
  if (delimiter_n == std::string::npos)
    return error_message("");
  else {
    auto posf = s.find(Lexer::label_delimiters[delimiter_n], pos + 1);
    if (posf == std::string::npos) {
      return error_message(Lexer::to_end_of_line(s, pos) +
                           " does not end with an label delimiter");
    } else {
      auto label = s.substr(pos, posf - pos + 1);
      return std::pair(std::unique_ptr<untyped_string_value<Lexer>>(
                           new untyped_string_value<Lexer>(label)),
                       posf + 1);
    }
  }
}

bool is_equal_at_pos(const std::string &s, std::size_t pos,
                     std::string_view n) {
  return s.substr(pos, n.size()) == n;
}

Maybe_error<
    std::pair<std::unique_ptr<untyped_numeric_value<Lexer>>, std::size_t>>
extract_numeric(const std::string &s, std::size_t pos) {
  auto last_pos = Lexer::skip_plus_or_minus(s, pos);
  if (!Lexer::is_numeral(s, last_pos)) {
    return error_message("");
  } else {
    last_pos = s.find_first_not_of(Lexer::numerals, last_pos);

    // now remove the fraction
    if (is_equal_at_pos(s, last_pos, Lexer::decimal_sep)) {
      last_pos = s.find_first_not_of(Lexer::numerals, last_pos + 1);
    }

    // now remove the exponent
    if (Lexer::is_Exponent(s, last_pos)) {
      last_pos = Lexer::skip_plus_or_minus(s, last_pos + 1);
      if (!Lexer::is_numeral(s, last_pos)) {
        return error_message(Lexer::to_end_of_line(s, pos) +
                             " the exponent is not numeric");
      } else {
        last_pos = s.find_first_not_of(Lexer::numerals, last_pos);
      }
    }
    return std::pair{std::make_unique<untyped_numeric_value<Lexer>>(
                         s.substr(pos, last_pos - pos)),
                     last_pos};
  }
}

Maybe_error<std::pair<std::unique_ptr<untyped_expression<Lexer>>, std::size_t>>
extract_expression(const std::string &s, std::size_t pos);

Maybe_error<std::pair<std::unique_ptr<untyped_assignment<Lexer>>, std::size_t>>
extract_assignment(const std::unique_ptr<untyped_identifier<Lexer>> &id,
                   const std::string &s, std::size_t pos) {
  auto last_pos = Lexer::skip_whitespace(s, pos);

  auto may_a = s.substr(last_pos, Lexer::assignment_operator.size());
  if (s.substr(last_pos, Lexer::assignment_operator.size()) !=
      Lexer::assignment_operator)
    return error_message("");
  else {
    last_pos = last_pos + Lexer::assignment_operator.size();
  }
  auto maybe_expression = extract_expression(s, last_pos);
  if (!maybe_expression)
    return maybe_expression.error();
  else {
    return std::pair(std::make_unique<untyped_assignment<Lexer>>(
                         id->clone(), (*maybe_expression).first.release()),
                     (*maybe_expression).second);
  }
}

Maybe_error<
    std::pair<std::unique_ptr<untyped_argument_list<Lexer>>, std::size_t>>
extract_argument_list(const std::string &s, std::size_t pos) {
  auto last_pos = Lexer::skip_whitespace(s, pos);
  if (!Lexer::is_start_of_argument_list(s, last_pos)) {
    return error_message("");
  } else {
    last_pos = last_pos + Lexer::argument_list_start.size();
    last_pos = Lexer::skip_whitespaceline(s, last_pos);

    auto out = std::make_unique<untyped_argument_list<Lexer>>();

    bool is_end = Lexer::is_end_of_argument_list(s,last_pos);
    while (!is_end) {
      auto maybe_expression = extract_expression(s, last_pos);
      if (!maybe_expression) {
        return error_message(s.substr(pos, last_pos - pos) +
                             "\t error in the next argument:\n " +
                             maybe_expression.error()() + "\n");
      } else {
        out->push_back((*maybe_expression).first.release());
        last_pos = (*maybe_expression).second;
        last_pos = Lexer::skip_whitespaceline(s, last_pos);

        if (Lexer::is_argument_separator(s,last_pos)) {
          last_pos += Lexer::argument_list_sep.size();
        } else {
          is_end = Lexer::is_end_of_argument_list(s,last_pos);
          if (!is_end) {
            return error_message(
                s.substr(pos, last_pos - pos - 1) +
                " expected argument separator \"" +
                std::string(Lexer::argument_list_sep) +
                "\"or argument end: \""+
                std::string(Lexer::argument_list_end) +
                "\" found: \"" +
                s.substr(last_pos, Lexer::argument_list_end.size()) + "\"\n");
          }
        }
      }
    }
    return std::pair(std::move(out),
                     last_pos + Lexer::argument_list_end.size());
  }
}

Maybe_error<
    std::pair<std::unique_ptr<untyped_function_evaluation<Lexer>>, std::size_t>>
extract_function_evaluation(
    const std::unique_ptr<untyped_identifier<Lexer>> &id, const std::string &s,
    std::size_t pos) {
  auto maybe_arguments = extract_argument_list(s, pos);
  if (!maybe_arguments)
    return maybe_arguments.error();
  else
    return std::pair(std::unique_ptr<untyped_function_evaluation<Lexer>>(
                         new untyped_function_evaluation<Lexer>(
                             id->clone(), (*maybe_arguments).first.release())),
                     (*maybe_arguments).second);
}

Maybe_error<std::pair<std::unique_ptr<untyped_expression<Lexer>>, std::size_t>>
extract_expression(const std::string &s, std::size_t pos) {
  auto last_pos = Lexer::skip_whitespaceline(s, pos);
  auto maybe_identifier = extract_identifier(s, last_pos);

  if (maybe_identifier) {
    last_pos = (*maybe_identifier).second;
    last_pos = Lexer::skip_whitespace(s, last_pos);
    if (Lexer::is_end_of_statement(s, last_pos) ||
        Lexer::is_end_of_argument(s, last_pos))
      return std::pair{std::unique_ptr<untyped_expression<Lexer>>{
                           (*maybe_identifier).first.release()},
                       last_pos};
    else {
      auto maybe_assigment =
          extract_assignment((*maybe_identifier).first, s, last_pos);
      if (maybe_assigment) {
        return std::pair{std::unique_ptr<untyped_expression<Lexer>>{
                             (*maybe_assigment).first.release()},
                         (*maybe_assigment).second};
      } else if (!maybe_assigment.error()().empty()) {
        return maybe_assigment.error();
      } else {
        auto maybe_function_evaluation =
            extract_function_evaluation((*maybe_identifier).first, s, last_pos);
        if (maybe_function_evaluation) {
          return std::pair{std::unique_ptr<untyped_expression<Lexer>>{
                               (*maybe_function_evaluation).first.release()},
                           (*maybe_function_evaluation).second};
        } else if (!maybe_function_evaluation.error()().empty()) {
          return maybe_function_evaluation.error();
        } else {
          // auto maybe_operation_evaluation here should go
          return error_message(s.substr(pos, last_pos - pos) +
                               "\n is not a well formed statment\n");
        }
      }
    }
  }

  auto maybe_numeric = extract_numeric(s, pos);
  if (maybe_numeric) {
    // in the future check for operations
    return std::pair(std::unique_ptr<untyped_expression<Lexer>>(
                         (*maybe_numeric).first.release()),
                     (*maybe_numeric).second);
  } else if (!maybe_numeric.error()().empty())
    return maybe_numeric.error();
  auto maybe_string = extract_string(s, pos);
  if (maybe_string) {
    // in the future check for operations
    return std::pair(std::unique_ptr<untyped_expression<Lexer>>(
                         (*maybe_string).first.release()),
                     (*maybe_string).second);
  } else if (!maybe_string.error()().empty())
    return maybe_string.error();

  auto maybe_argument_list = extract_argument_list(s, pos);
  if (maybe_argument_list) {
    // in the future check for operations
    return std::pair(std::unique_ptr<untyped_expression<Lexer>>(
                         (*maybe_argument_list).first.release()),
                     (*maybe_argument_list).second);
  } else if (!maybe_argument_list.error()().empty())
    return maybe_argument_list.error();

  return error_message(s.substr(pos, s.find_first_of("\n", last_pos) - pos) +
                       "\n is not a well formed expression");
}

Maybe_error<untyped_program<Lexer>> extract_program(const std::string &s) {
  untyped_program<Lexer> out;
  auto last_pos = 0ul;
  last_pos = Lexer::skip_whitespaceline(s, last_pos);
  while (last_pos < s.size()) {
    auto maybe_expr = extract_expression(s, last_pos);
    if (maybe_expr) {
      out.push_back((*maybe_expr).first.release());
      last_pos = (*maybe_expr).second;
      last_pos = Lexer::skip_whitespaceline(s, last_pos);

    } else {
      return maybe_expr.error();
    }
  }
  return out;
}

} // namespace dcli

#endif // LEXER_UNTYPED_H
