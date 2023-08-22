#ifndef GRAMMAR_UNTYPED_H
#define GRAMMAR_UNTYPED_H

#include "grammar_Identifier.h"
#include <memory>
#include <string>
#include <vector>
namespace dcli {

template <class Abstract>
auto clone(const std::vector<std::unique_ptr<Abstract>>& x) {
  std::vector<std::unique_ptr<Abstract>> out;
  for (auto &ele : x)
      out.push_back(std::unique_ptr<Abstract>(ele->clone()));
  return out;
}

template <class Abstract>
auto clone(const std::unique_ptr<Abstract>& x) {
  return std::unique_ptr<Abstract>(x->clone());
}


template <class Lexer> class untyped_statement {
public:
  virtual ~untyped_statement(){};
  virtual std::string str() const = 0;
  virtual untyped_statement *clone() const = 0;
};

template <class Lexer>
class untyped_expression : public untyped_statement<Lexer> {
public:
  virtual ~untyped_expression(){};
  virtual untyped_expression *clone() const = 0;
};

template <class Lexer> class untyped_program {
  std::vector<std::unique_ptr<untyped_statement<Lexer>>> m_statements;

public:
  untyped_program(){}
  untyped_program(const untyped_program &other)
      : m_statements{clone(other.m_statements)} {}

  untyped_program(untyped_program &&other)
      : m_statements{std::move(other.m_statements)} {}

  auto &push_back(untyped_statement<Lexer> *t_expr) {
    m_statements.emplace_back(t_expr);
    return *this;
  }
  std::string str() const  {
    std::string out = "";
    for (auto &elem : m_statements)
        out += elem->str() +std::string( Lexer::statement_sep);
    return out;
  };



};

template <class Lexer>
class untyped_simple_expression : public untyped_expression<Lexer> {
  std::string m_expression;

public:
  virtual ~untyped_simple_expression(){};
  virtual untyped_simple_expression *clone() const override {
    return new untyped_simple_expression(*this);
  };

  untyped_simple_expression(std::string t_expression)
      : m_expression{std::move(t_expression)} {}

  virtual std::string str() const override final { return m_expression; };
  auto &operator()() const { return m_expression; }

};

template <class Lexer>
class untyped_numeric_value : public untyped_simple_expression<Lexer> {
public:
  using untyped_simple_expression<Lexer>::untyped_simple_expression;
  virtual untyped_numeric_value *clone() const override {
    return new untyped_numeric_value(*this);
  };
  virtual ~untyped_numeric_value(){};
};

template <class Lexer>
class untyped_string_value : public untyped_simple_expression<Lexer> {

public:
  using untyped_simple_expression<Lexer>::untyped_simple_expression;

  virtual untyped_string_value *clone() const override {
    return new untyped_string_value(*this);
  };
  virtual ~untyped_string_value(){};
};

template <class Lexer>
class untyped_identifier : public untyped_expression<Lexer> {
  Identifier<Lexer> m_id;

public:
  virtual ~untyped_identifier(){};
    untyped_identifier(Identifier<Lexer> t_id) : m_id{t_id} {}
  virtual untyped_identifier *clone() const override {
    return new untyped_identifier(*this);
  };

  virtual std::string str() const override final { return m_id();}
  auto& operator()() const { return m_id; }
};

template <class Lexer> class untyped_argument_list : public untyped_expression<Lexer> {
  std::vector<std::unique_ptr<untyped_expression<Lexer>>> m_list;

public:
  untyped_argument_list(
      const std::vector<std::unique_ptr<untyped_expression<Lexer>>> &t_list)
        : m_list{dcli::clone(t_list)} {}

  untyped_argument_list(
       std::vector<std::unique_ptr<untyped_expression<Lexer>>> &&t_list)
      : m_list{t_list} {}
  untyped_argument_list(
      const untyped_argument_list &other)
      : m_list{dcli::clone(other.m_list)} {}

  untyped_argument_list(
       untyped_argument_list &&other)
      : m_list{std::move(other.t_list)} {}

  untyped_argument_list(){}

  untyped_argument_list &push_back(untyped_expression<Lexer> *t_expr) {
    m_list.emplace_back(t_expr);
    return *this;
  }
  virtual untyped_argument_list *clone() const override {return new untyped_argument_list(*this);}

  virtual std::string str() const override final {
    std::string out = std::string(Lexer::argument_list_start);
    if (m_list.size() > 0) {
      out += m_list[0]->str();
      for (auto i = 1ul; i < m_list.size(); ++i) {
          out += std::string(Lexer::argument_sep) + m_list[i]->str();
      }
    }
    out += Lexer::argument_list_end;
    return out;
  }
};

template <class Lexer>
class untyped_assignment : public untyped_expression<Lexer> {
  std::unique_ptr<untyped_identifier<Lexer>> m_id;
  std::unique_ptr<untyped_expression<Lexer>> m_expr;

public:
  virtual ~untyped_assignment() override {}
  untyped_assignment(untyped_identifier<Lexer> *t_id,
                     untyped_expression<Lexer> *t_expr)
      : m_id{t_id}, m_expr{t_expr} {}

  untyped_assignment() {}
  untyped_assignment(untyped_assignment &&other)
      : m_id{std::move(other.m_id)}, m_expr{std::move(other.m_expr)} {}
  untyped_assignment(untyped_assignment const &other)
      : m_id{other.m_id->clone()}, m_expr{other.m_expr->clone()} {}

  virtual std::string str() const override final {
    return m_id->str() + std::string(Lexer::assignment_operator) + std::string(Lexer::assignment_sep) +
           m_expr->str();
  }
  virtual untyped_assignment *clone() const override {
    return new untyped_assignment(*this);
  };
};

template <class Lexer>
class untyped_function_evaluation : public untyped_expression<Lexer> {
  std::unique_ptr<untyped_identifier<Lexer>> m_fid;
  std::unique_ptr<untyped_argument_list<Lexer>> m_arg;

public:
  untyped_function_evaluation(untyped_identifier<Lexer> *t_fid,
                                untyped_argument_list<Lexer> *t_arg)
        : m_fid{t_fid}, m_arg{t_arg }{}
  virtual ~untyped_function_evaluation(){};
  untyped_function_evaluation(untyped_function_evaluation&& other):
      m_fid{std::move(other.m_fid)},m_arg{std::move(other.m_arg)}{}

  untyped_function_evaluation(untyped_function_evaluation const & other):
      m_fid{other.m_fid->clone()},m_arg{other.m_arg->clone()}{}

  virtual untyped_function_evaluation* clone() const override {return new untyped_function_evaluation(*this);}
  virtual std::string str() const override final {
      return m_fid->str() + m_arg->str();
  }


};

} // namespace dcli

#endif // GRAMMAR_UNTYPED_H
