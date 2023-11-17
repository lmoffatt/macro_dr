#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
#include "grammar_typed.h"
#include "grammar_untyped.h"
#include "lexer_untyped.h"
#include <iostream>
#include <map>

namespace dcli {

template <class Lexer, class Compiler> class base_Identifier_compiler {

public:
  virtual ~base_Identifier_compiler(){};
  virtual base_Identifier_compiler *clone() const = 0;

  virtual Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_Identifier(const Identifier<Lexer> &id) const = 0;
};

template <class Lexer, class Compiler, class T>
class Identifier_compiler : public base_Identifier_compiler<Lexer, Compiler> {
  std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

public:
  virtual ~Identifier_compiler(){};
  virtual Identifier_compiler *clone() const {
    return new Identifier_compiler(*this);
  };

  Identifier_compiler(typed_expression<Lexer, Compiler, T> &t_expr)
      : m_expr{t_expr} {}

  virtual Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_Identifier(const Identifier<Lexer> &id) const {
    return new typed_identifier<Lexer, Compiler, T>(id);
  }
};

template <class Lexer, class Compiler> class base_function_compiler {
public:
  virtual ~base_function_compiler() {}

  virtual Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_function_evaluation(
      Compiler const &cm,
      const typed_argument_list<Lexer, Compiler> *args) const = 0;
};

template <class Lexer, class Compiler, class F, class... Args>
class function_compiler : public base_function_compiler<Lexer, Compiler> {
  std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, Args>>...>
      m_args;
  F m_f;
  
  
  template <class It, class Map,
           class... Args0>
  Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_function_evaluation_impl(
      std::index_sequence<>, Compiler const &cm, It begin, It end,
      const Map &map,
      typed_expression<Lexer, Compiler, Args0>*...args)
      const {
    return new typed_function_evaluation<Lexer,Compiler,F,Args...>(m_f,std::move(args)...);
  }
  
  
  template <std::size_t I, std::size_t... Is, class It, class Map,
            class... Args0>
  Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_function_evaluation_impl(
      std::index_sequence<I, Is...>, Compiler const &cm, It begin, It end,
      const Map &map,
      typed_expression<Lexer, Compiler, Args0>*...args)
      const {
    if (begin != end) {
      auto Maybe_arg = std::get<I>(m_args)->compile_this_argument(*begin);
      if (!Maybe_arg)
        return Maybe_arg.error();
      else
        return compile_function_evaluation_impl(
            std::index_sequence<Is...>(), cm, ++begin, end, map, args...,
            Maybe_arg.value());
    }
    else
    {
      auto Maybe_arg = std::get<I>(m_args)->compile_this_argument_in_this_map(map);
      if (!Maybe_arg)
        return Maybe_arg.error();
      else
        return compile_function_evaluation_impl(
            std::index_sequence<Is...>(), cm, begin, end, map, args...,
            Maybe_arg.value());
      
    }
  }

public:
  function_compiler(F t_f, typed_expression<Lexer, Compiler, Args> *...t_args)
      : m_f{t_f}, m_args{t_args...} {}

  virtual Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_function_evaluation(
      Compiler const &cm,
      const typed_argument_list<Lexer, Compiler>*args) const {
    
    return compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm,args->arg_vector().begin(),args->arg_vector().end(),args->arg_map());
    
  }
};

class Compiler {
  std::map<Identifier<Lexer>,
           std::unique_ptr<base_function_compiler<Lexer, Compiler>>>
      m_func;
  std::map<Identifier<Lexer>,
           std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>>
      m_id;

public:
  Compiler() {}
  Compiler(
      std::map<Identifier<Lexer>,
               std::unique_ptr<base_function_compiler<Lexer, Compiler>>> &&func)
      : m_func{std::move(func)}, m_id{} {}

  Maybe_error<base_function_compiler<Lexer, Compiler> const *>
  get_function(const Identifier<Lexer> &id) const {
    auto it = m_func.find(id);
    if (it == m_func.end())
      return error_message(id() + " function is not defined");
    else {
      auto ptr = (*it).second.get();
      if (ptr == nullptr)
        return error_message(id() + " function is null");
      else
        return ptr;
    }
  }

  Maybe_error<base_typed_expression<Lexer, Compiler>  *>
  get_Identifier(const Identifier<Lexer> &id) const {
    auto it = m_id.find(id);
    if (it == m_id.end())
      return error_message(id() + " function is not defined");
    else {
      auto ptr = (*it).second.get();
      if (ptr == nullptr)
        return error_message(id() + " function is null");
      else
        return ptr->compile_Identifier(id);
    }
  }

  Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_Identifier(const Identifier<Lexer> &id) const {
    auto Maybe_id = get_Identifier(id);
    if (!Maybe_id)
      return Maybe_id.error();
    else
      return Maybe_id.value();
  }

  void push_back(Identifier<Lexer> id,
                 base_Identifier_compiler<Lexer, Compiler> *expr) {
    m_id.emplace(id, expr);
  }

  Maybe_error<bool>
  push_function(std::string id_candidate,
                base_function_compiler<Lexer, Compiler> *fun) {
    auto may_id = to_Identifier<Lexer>(id_candidate);
    if (!may_id)
      return may_id.error();
    m_func.emplace(may_id.value(), fun);
    return true;
  }
};

Maybe_error<typed_program<Lexer, Compiler>>
compile_program(Compiler &cm, const untyped_program<Lexer, Compiler> &s) {
  typed_program<Lexer, Compiler> out;
  for (auto &e : s.statements()) {
    auto Maybe_compiled_statement = e->compile_statement(cm);
    if (!Maybe_compiled_statement)
      return Maybe_compiled_statement.error();
    else
      out.push_back(Maybe_compiled_statement.value());
  }
  return out;
}

} // namespace dcli

#endif // LEXER_TYPED_H
