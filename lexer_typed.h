#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
//#include "grammar_typed.h"
#include "maybe_error.h"
#include <iostream>
#include <map>
#include <memory>
#include "grammar_Identifier.h"
namespace dcli {

template <class Lexer, class Compiler> class base_typed_expression;
template <class Lexer, class Compiler> class base_Identifier_compiler;
template <class Lexer, class Compiler, class T> class typed_expression;

template <class Lexer, class Compiler, class T>
class typed_identifier;
template <class Lexer, class Compiler>
class untyped_argument_list;
class Lexer;
template <class Lexer, class Compiler> class typed_argument_list;

template <class Lexer, class Compiler, class F, class... Args>
class typed_function_evaluation;

template <class Lexer, class Compiler> class untyped_program ;

template <class Lexer, class Compiler> class typed_program ;

template <class Lexer, class Compiler, class...T>
class typed_argument_typed_list; 


template <class Lexer, class Compiler> class base_Identifier_compiler {

public:
  virtual ~base_Identifier_compiler(){};
  virtual base_Identifier_compiler *clone() const = 0;
  
  virtual Maybe_error<std::unique_ptr<base_typed_expression<Lexer, Compiler>>>
  compile_Identifier(const Identifier<Lexer> &id) const = 0;
};

template <class Lexer, class Compiler, class T>
class Identifier_compiler : public base_Identifier_compiler<Lexer, Compiler> {
  std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

public:
  virtual ~Identifier_compiler(){};
  virtual Identifier_compiler *clone() const {
      return new Identifier_compiler(m_expr->clone());
  };

  Identifier_compiler(typed_expression<Lexer, Compiler, T> *t_expr)
      : m_expr{t_expr} {}
  
  virtual Maybe_error<std::unique_ptr<base_typed_expression<Lexer, Compiler>>>
  compile_Identifier(const Identifier<Lexer> &id) const {
    return new typed_identifier<Lexer, Compiler, T>(id);
  }
};





template <class Lexer, class Compiler> class base_function_compiler {
public:
  virtual ~base_function_compiler() {}
    
    virtual base_function_compiler* clone()const =0;
  
  virtual Maybe_unique<base_typed_expression<Lexer, Compiler> >
  compile_function_evaluation(
      Compiler const &cm,
      const untyped_argument_list<Lexer, Compiler> &args) const = 0;
  
  
};

template <class Lexer, class Compiler, class F, class... Args>
class function_compiler : public base_function_compiler<Lexer, Compiler> {
  std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, Args>>...>
      m_args;
  F m_f;
  
  
  template <std::size_t... Is>
Maybe_unique<base_typed_expression<Lexer, Compiler>>
  compile_function_evaluation_impl(
      std::index_sequence<Is...>, Compiler const &cm, 
      const untyped_argument_list<Lexer, Compiler>& args)
      const {
      auto Maybe_tuple=promote_Maybe_error(std::tuple(std::get<Is>(m_args)->compile_this_argument(cm,args.arg()[Is])...));
      if (!Maybe_tuple)
          return Maybe_tuple.error();
      else
      {
          return new typed_argument_typed_list<Lexer,Compiler,Args...>(std::move(Maybe_tuple.value()));
      }
  }
  
  
  
  
  
public:
  function_compiler(F t_f, typed_expression<Lexer, Compiler, Args> *...t_args)
      : m_f{t_f}, m_args{t_args...} {}
  
  function_compiler(function_compiler const & x)
      : m_f{x.m_f}, m_args{clone_tuple(x.m_args)} {}
  
  virtual Maybe_unique<base_typed_expression<Lexer, Compiler>>
  compile_function_evaluation(
      Compiler const &cm,
      const untyped_argument_list<Lexer, Compiler>& args) const {
      
     // using jeje=typename decltype(compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm,args.arg()))::kk;
      return compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm,args.arg());
      
  }
  
  // base_function_compiler interface
  public:
  virtual base_function_compiler<Lexer,Compiler> *clone() const override
  {
      return new function_compiler(*this);
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
  
  Compiler(const Compiler& cm)
      : m_func{clone_map(cm.m_func)}, m_id{clone_map(cm.m_id)} {}
  
  
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
  
  Maybe_unique<base_typed_expression<Lexer, Compiler>>
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
  
  Maybe_unique<base_typed_expression<Lexer, Compiler> >
  compile_Identifier(const Identifier<Lexer> &id) const {
    auto Maybe_id = get_Identifier(id);
    if (!Maybe_id)
      return Maybe_id.error();
    else
        return std::move(Maybe_id.value());
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





} // namespace dcli

#endif // LEXER_TYPED_H
