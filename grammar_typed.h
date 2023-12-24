#ifndef GRAMMAR_TYPED_H
#define GRAMMAR_TYPED_H
#include "grammar_Identifier.h"
//#include "grammar_untyped.h"
#include "lexer_typed.h"
#include "maybe_error.h"
#include "qmodel.h"
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

namespace dcli {


template <class Lexer, class Compiler> class Environment {
  std::map<Identifier<Lexer>,
           std::unique_ptr<base_typed_expression<Lexer, Compiler>>>
      m_var;

public:
  void insert(Identifier<Lexer> const &id,
              base_typed_expression<Lexer, Compiler> *expr) {
      m_var.insert_or_assign(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(expr));
  }

  Maybe_error<base_typed_expression<Lexer, Compiler> const *>
  get(Identifier<Lexer> const &id) const {
    auto it = m_var.find(id);
    if (it == m_var.end())
      return error_message("identifier not found");
    else
      return (*it).second.get();
  }
};

template <class Lexer, class Compiler> class typed_argument_list;

template <class Lexer, class Compiler> class base_typed_statement {
public:
  virtual ~base_typed_statement(){};
  virtual base_typed_statement *clone() const = 0;
  virtual Maybe_error<typed_argument_list<Lexer, Compiler>*>
  compile_argument_list(
      typed_argument_list<Lexer, Compiler> *t_growing_list) const = 0;
  virtual Maybe_error<bool> run_statement(Environment<Lexer,Compiler>& env)const=0;
  
  virtual base_Identifier_compiler<Lexer, Compiler>*
  compile_identifier()const=0;
  
};

template <class Lexer, class Compiler>
class base_typed_assigment : public base_typed_statement<Lexer, Compiler> {
public:
  virtual ~base_typed_assigment(){};
  virtual base_typed_assigment *clone() const = 0;
  virtual Identifier<Lexer> const &id() const = 0;
  virtual base_typed_expression<Lexer, Compiler> *expr() const = 0;

  virtual Maybe_error<typed_argument_list<Lexer, Compiler>*>
  compile_argument_list(
      typed_argument_list<Lexer, Compiler> *t_growing_list) const;
  
 
};

template <class Lexer, class Compiler>
class base_typed_expression : public base_typed_statement<Lexer, Compiler> {
public:
  virtual ~base_typed_expression(){};
  virtual base_typed_expression *clone() const = 0;
  
  virtual base_typed_expression<Lexer, Compiler> *
  compile_identifier(Identifier<Lexer> id) const= 0;

  virtual base_typed_assigment<Lexer, Compiler> *
  compile_assigment(Identifier<Lexer> id) = 0;
  
  virtual Maybe_error<typed_argument_list<Lexer, Compiler>*>
  compile_argument_list(
      typed_argument_list<Lexer, Compiler> *t_growing_list) const;
  
  virtual Maybe_error<base_typed_expression<Lexer,Compiler>*> run_expression(Environment<Lexer,Compiler>& env)const=0;
 
  virtual Maybe_error<bool> run_statement(Environment<Lexer,Compiler>& env)const
  {
    auto Maybe_exp=run_expression(env);
    if (!Maybe_exp) return Maybe_exp.error();
    else
     return true;
  }
  
  virtual base_Identifier_compiler<Lexer, Compiler>*
  compile_identifier() const=0; 
  
};

template <class Lexer, class Compiler>
class typed_argument_list : public base_typed_expression<Lexer, Compiler> {
  std::vector<std::unique_ptr<base_typed_expression<Lexer, Compiler>>> m_vec;
  std::map<Identifier<Lexer>,
           std::unique_ptr<base_typed_expression<Lexer, Compiler>>>
      m_map;
    
public:
  virtual ~typed_argument_list(){};
    virtual typed_argument_list *clone() const {
    return new typed_argument_list(*this);
    }
    typed_argument_list(const typed_argument_list& other):m_vec(clone_vector(other.m_vec)), m_map{clone_map(m_map)}{}
  auto& arg_vector()const {return m_vec;}
  auto& arg_map()const{return m_map;}
  void insert_assigment(base_typed_assigment<Lexer, Compiler> const &ass) {
    m_map[ass.id()] = std::unique_ptr<base_typed_expression<Lexer, Compiler>>(ass.expr()->clone());
  }

  Maybe_error<bool>
  push_back_expression(base_typed_expression<Lexer, Compiler> const &expr) {
    if (!m_map.empty())
      return error_message("assigment after expression in argument list");
    else
    {
      m_vec.emplace_back(expr.clone());
      return true;
    }
  }

  typed_argument_list() {}
  base_typed_assigment<Lexer, Compiler> *
  compile_assigment(Identifier<Lexer> id) {
 return nullptr;
  }
  
  virtual base_Identifier_compiler<Lexer, Compiler>*
  compile_identifier()const override {return nullptr;}
  
  virtual base_typed_expression<Lexer, Compiler> *
  compile_identifier(Identifier<Lexer> id) const{ return nullptr;};
  
  
  // base_typed_expression interface
  public:
  virtual Maybe_error<dcli::base_typed_expression<Lexer, Compiler> *> run_expression(dcli::Environment<Lexer, Compiler> &env) const override
  {
 return error_message("no reason to get here");
  }
};


template <class Lexer, class Compiler, class...T>
class typed_argument_typed_list : public base_typed_expression<Lexer, Compiler> {
    std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,T>>...> m_args;
    
    
public:
    virtual ~typed_argument_typed_list(){};
    virtual typed_argument_typed_list *clone() const {
        return new typed_argument_typed_list(clone_tuple(m_args));
    }
    typed_argument_typed_list(std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,T>>...> && t):m_args{std::move(t)}{}
    auto& args()const {return m_args;}
    
    
    virtual base_Identifier_compiler<Lexer, Compiler>*
    compile_identifier() const{return nullptr;};
    
    virtual base_typed_expression<Lexer, Compiler> *
    compile_identifier(Identifier<Lexer> id) const{return nullptr;};
    
    virtual base_typed_assigment<Lexer, Compiler> *
    compile_assigment(Identifier<Lexer> id){return nullptr;};
    
    
    // base_typed_statement interface
public:
    
    // base_typed_expression interface
public:
    virtual Maybe_error<dcli::base_typed_expression<Lexer, Compiler> *> run_expression(dcli::Environment<Lexer, Compiler> &env) const override
    {
        return error_message("still unresolved");
    }
};




template <class Lexer, class Compiler, class T> class typed_assigment;

template <class Lexer, class Compiler, class T>
class typed_expression : public base_typed_expression<Lexer, Compiler> {

public:
  virtual ~typed_expression(){};
  virtual typed_expression *clone() const = 0;

  virtual T run(Environment<Lexer, Compiler> const &) const = 0;
  
  
  
  virtual base_typed_assigment<Lexer, Compiler> *
  compile_assigment(Identifier<Lexer> id);

  virtual base_typed_expression<Lexer, Compiler> *
  compile_identifier(Identifier<Lexer> id) const override {
      return new typed_identifier<Lexer, Compiler, T>(id);}


   
  
  
  
  virtual base_Identifier_compiler<Lexer, Compiler>*
   compile_identifier() const override{
      return new Identifier_compiler<Lexer,Compiler,T>(this->clone());   
  }
   
   
   
   // base_typed_expression interface
   public:
  virtual Maybe_error<dcli::base_typed_expression<Lexer, Compiler> *> run_expression(dcli::Environment<Lexer, Compiler> &env) const override
  {
      if constexpr (std::is_void_v<T>)
      {
          run(env);
           return new typed_literal<Lexer,Compiler,T>();
      }
      else   
      return new typed_literal<Lexer,Compiler,T>(run(env));
  }
       base_typed_expression<Lexer, Compiler> *compile_identifier(Identifier<Lexer> id) {
           return new typed_identifier<Lexer, Compiler, T>(id);
       }
};

template <class Lexer, class Compiler, class T>
class typed_assigment : public base_typed_assigment<Lexer, Compiler> {
  Identifier<Lexer> m_id;
  std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

public:
  virtual ~typed_assigment(){};
  virtual typed_assigment *clone() const { return new typed_assigment(*this); };
  
  typed_assigment(const typed_assigment& other):m_id{other.m_id}, m_expr{other.m_expr->clone()}{}

  typed_assigment(Identifier<Lexer> const &t_id,
                  typed_expression<Lexer, Compiler, T> *t_expr)
      : m_id{t_id}, m_expr{t_expr} {}
  
  virtual Identifier<Lexer> const &id() const override { return m_id; }
  virtual typed_expression<Lexer, Compiler, T> *expr() const override {
    return m_expr.get();
  }

    
  virtual Maybe_error<bool> run_statement(Environment<Lexer,Compiler>& env)const
  {
    auto Maybe_exp=expr()->run_expression(env);
    if (!Maybe_exp) return Maybe_exp.error();
    else
    {
      env.insert(id(), Maybe_exp.value());
      return true;
    } 
  }
  
  virtual base_Identifier_compiler<Lexer, Compiler>*
  compile_identifier() const override {
      return new Identifier_compiler<Lexer, Compiler,T>(this->expr());
  }
  
  
};

template <class Lexer, class Compiler, class F, class... Args>
class typed_function_evaluation
    : public typed_expression<Lexer, Compiler,
                              std::invoke_result_t<F, Args...>> {
  std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, Args>>...>
      m_args;
  F m_f;

public:
  using T=std::invoke_result_t<F, Args...>;
  typed_function_evaluation(F t_f,
                            typed_expression<Lexer, Compiler, Args> *...t_args)
      : m_f{t_f}, m_args{t_args...} {}
  
  typed_function_evaluation(F t_f,
                            std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,Args>>...> && t):m_f{t_f},m_args{std::move(t)}{}
  
  // typed_expression interface
  public:
  typed_function_evaluation(const typed_function_evaluation& other):
          m_args{std::apply([](auto const&  ...args){
              return std::tuple(std::unique_ptr<typed_expression<Lexer, Compiler, Args>>(args->clone())...);},other.m_args)},m_f{other.m_f}{}
      
  virtual typed_function_evaluation *clone() const {
    return new typed_function_evaluation(*this);
  }
  
  virtual T run(const Environment<Lexer, Compiler> &env) const override
  {
   // return T{};
      return std::apply([this, &env](auto& ...args){return this->m_f(args->run(env)...);},m_args);
//   return std::apply([this, &env](auto& ...args){return std::invoke(this->m_f,args->run(env)...);},m_args);
  }
  
    
};

template <class Lexer, class Compiler, class T>
base_typed_assigment<Lexer, Compiler> *
typed_expression<Lexer, Compiler, T>::compile_assigment(Identifier<Lexer> id) {
  return new typed_assigment<Lexer, Compiler, T>(id, this->clone());
}

template <class Lexer, class Compiler, class T>
    requires (!std::is_void_v<T>)
class typed_literal<Lexer,Compiler,T> : public typed_expression<Lexer, Compiler, T> {
  T m_value;

public:
  typed_literal(T &&x) : m_value(std::move(x)) {}
  typed_literal(T const &x) : m_value(x) {}

  virtual ~typed_literal(){};

  virtual T run(Environment<Lexer, Compiler> const &) const override{
    return m_value;
  }
  
  // typed_expression interface
  public:
  virtual typed_literal *clone() const override
  {
    return new typed_literal(*this);
  }
  
};

template <class Lexer, class Compiler>
class typed_literal<Lexer, Compiler, void> : public typed_expression<Lexer, Compiler, void> {
    
public:
    
    virtual ~typed_literal(){};
    
    virtual void run(Environment<Lexer, Compiler> const &) const override{
    }
    
    // typed_expression interface
public:
    virtual typed_literal *clone() const override
    {
        return new typed_literal(*this);
    }
    
};


template <class Lexer, class Compiler, class T>
class typed_identifier : public typed_expression<Lexer, Compiler, T > {
  Identifier<Lexer> m_id;

public:
  typed_identifier(Identifier<Lexer> &&x) : m_id(std::move(x)) {}
  typed_identifier(Identifier<Lexer> const &x) : m_id(x) {}
  
  auto& id()const { return m_id;}
  virtual ~typed_identifier(){};
  virtual typed_identifier *clone() const override {
    return new typed_identifier(*this);
  }
  
  virtual T run(Environment<Lexer, Compiler> const &env) const {
    auto May_x = env.get(m_id);
    auto exp = dynamic_cast<typed_expression<Lexer, Compiler, T> const *>(
        May_x.value());
    return exp->run(env);
  }
};







 

template <class Lexer, class Compiler> class typed_program {
  std::vector<std::unique_ptr<base_typed_statement<Lexer, Compiler>>>
      m_statements;
  
  std::map<Identifier<Lexer>, base_typed_expression<Lexer, Compiler> *>
      m_identifier_table;
  
  public:
  typed_program() {}
  typed_program(const typed_program &other)
      : m_statements{clone(other.m_statements)} {}

  typed_program(typed_program &&other)
      : m_statements{std::move(other.m_statements)} {}

  auto &push_back(base_typed_statement<Lexer, Compiler> *t_expr) {
    m_statements.emplace_back(t_expr);
    return *this;
  }

  auto &insert(Identifier<Lexer> id,
               base_typed_expression<Lexer, Compiler> *t_expr) {
    m_identifier_table.insert_or_assign(id, t_expr);
    return *this;
  }

  Maybe_error<base_typed_expression<Lexer, Compiler> *>
  compile_identifier(Identifier<Lexer> const &id) const {
    auto it = m_identifier_table.find(id);
    if (it == m_identifier_table.end())
      return error_message("identifier " + id() + " is not defined");
    else
      return (*it)->compile_identifier(id);
  }
  
  Maybe_error<Environment<Lexer,Compiler>> run() {
    Environment<Lexer, Compiler> env;
    for (auto &e : m_statements)
    {
      auto Maybe_stat=e->run_statement(env);
      if (!Maybe_stat) return Maybe_stat.error();
    }
    return env;
  }
};

template <class Lexer, class Compiler>
Maybe_error<typed_argument_list<Lexer, Compiler>*>
base_typed_assigment<Lexer, Compiler>::compile_argument_list(
    typed_argument_list<Lexer, Compiler> *t_growing_list) const {
  t_growing_list->insert_assigment(*this);
  return t_growing_list;
}

template <class Lexer, class Compiler>
Maybe_error<typed_argument_list<Lexer, Compiler>*>
base_typed_expression<Lexer, Compiler>::compile_argument_list(
    typed_argument_list<Lexer, Compiler> *t_growing_list) const {
  auto May_be = t_growing_list->push_back_expression(*this);
  if (!May_be)
    return May_be.error();
  else
    return t_growing_list;
}






/*

template<class Lexer, class Compiler>
    class base_typed_comment: public base_typed_literal<Lexer, Compiler>{
    virtual ~base_typed_comment(){};
    virtual std::string str() const = 0;
    virtual base_typed_comment *clone() const = 0;

};


template<class Lexer, class Compiler>
class base_typed_assigment: public base_typed_expression<Lexer, Compiler>{
    virtual ~base_typed_assigment(){};
    virtual std::string str() const = 0;
    virtual base_typed_assigment *clone() const = 0;

};

template<class Lexer, class Compiler>
class base_typed_argument_list: public base_typed_expression<Lexer, Compiler>{
    virtual ~base_typed_argument_list(){};
    virtual std::string str() const = 0;
    virtual base_typed_argument_list *clone() const = 0;

};

*/

/*

template <class Lexer, class Compiler, class T>
class typed_assignment : public typed_expression<Lexer, Compiler,T> {
    std::unique_ptr<untyped_identifier<Lexer, Compiler>> m_id;
    std::unique_ptr<typed_expression<Lexer, Compiler,T>> m_expr;

public:
    virtual ~typed_assignment(){};

    // base_typed_statement interface
public:


    typed_assignment(untyped_identifier<Lexer, Compiler> *t_id,
                       typed_expression<Lexer, Compiler,T> *t_expr)
        : m_id{t_id}, m_expr{t_expr} {}

    typed_assignment() {}
    typed_assignment(typed_assignment &&other)
        : m_id{std::move(other.m_id)}, m_expr{std::move(other.m_expr)} {}
    typed_assignment(typed_assignment const &other)
        : m_id{other.m_id->clone()}, m_expr{other.m_expr->clone()} {}

    auto& id()const{ return *m_id;}
    auto& expression()const{ return *m_expr;}


    virtual std::string str() const override final {
        return m_id->str() + std::string(Lexer::assignment_operator) +
std::string(Lexer::assignment_sep) + m_expr->str();
    }
    virtual typed_assignment *clone() const override {
        return new typed_assignment(*this);
    };
    virtual Maybe_error<typed_expression<Lexer,Compiler,T>*>
    compile(const Compiler& t_cm, const untyped_expression<Lexer,Compiler>&
t_expr)const{ if
(t_expr.type_of_statement()!=untyped_assignment<Lexer,Compiler>::statement_type())
        {
            return error_message("");
        }else
        {
            auto& v_asigm=dynamic_cast<untyped_assignment<Lexer,Compiler>const
&>(t_expr); if (id().str()!=v_asigm.id().str())
            {
                return error_message("");
            }else{
                auto v_expr=expression().compile(t_cm,v_asigm.expression());
                if (!v_expr)
                    return v_expr.error();
                else
                {
                    new
typed_assignment<Lexer,Compiler,T>(id().clone(),v_expr.value());
                }
            }
        }
    }



};



template <class Lexer, class Compiler, class T>
class typed_identifier : public typed_expression<Lexer, Compiler,T> {
    Identifier<Lexer> m_id;

public:
    virtual ~typed_identifier(){};
    typed_identifier(Identifier<Lexer> t_id) : m_id{t_id} {}
    virtual typed_identifier *clone() const override {
        return new typed_identifier(*this);
    };

    virtual std::string str() const override final { return m_id();}
    virtual Maybe_error<typed_expression<Lexer,Compiler,T>*>
    compile(const Compiler& t_cm, const untyped_expression<Lexer,Compiler>&
t_expr)const{ return error_message("");
    }


    auto& operator()() const { return m_id; }
};








/**


template <class Lexer, class Compiler, class T>
    requires std::is_arithmetic_v<T>
class typed_numeric_value : public typed_simple_expression<Lexer,Compiler,T> {
    T m_number;
public:
    typed_numeric_value(T t_number): m_number{t_number}{}
    virtual typed_numeric_value *clone() const override {
        return new typed_numeric_value(*this);
    };


    virtual ~typed_numeric_value(){};
};

template < class Lexer, class Compiler>
class typed_string_value : public typed_simple_expression<Lexer,std::string> {
    std::string m_expression;

public:

    typed_string_value(std::string t_expression)
        : m_expression{std::move(t_expression)} {}

    virtual std::string str() const override  { return m_expression; };
    auto &operator()() const { return m_expression; }
};


template < class Lexer, class Compiler>
class typed_comment : public typed_statement<Lexer, Compiler> {
    std::string m_comment;

public:
    typed_comment(std::string t_comment)
        : m_comment{std::move(t_comment)} {}

    virtual typed_comment *clone() const override {
        return new typed_comment(*this);
    };

    virtual std::string str()const {
        return std::string(Lexer::comment_start)+m_comment;
    }
    virtual ~typed_comment(){};
};



template <class Lexer, class Compiler, class T>
class typed_identifier : public typed_expression<Lexer,Compiler,  T> {
    Identifier<Lexer> m_id;

public:
    virtual ~typed_identifier(){};
    typed_identifier(Identifier<Lexer> t_id) : m_id{t_id} {}
    virtual typed_identifier *clone() const override {
        return new typed_identifier(*this);
    };

    virtual std::string str() const override final { return m_id();}
    auto& operator()() const { return m_id; }
};




template <class Lexer, class Compiler, class...Ts> class typed_argument_list :
public base_typed_expression<Lexer, Compiler> {
    std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,Ts>>...> m_tuple;

public:
    typed_argument_list(
         std::unique_ptr<typed_expression<Lexer,Compiler,Ts>>const &... t_args)
        : m_tuple{dcli::clone(t_args)...} {}

    typed_argument_list(
        std::unique_ptr<typed_expression<Lexer,Compiler,Ts>>&&... t_args)
        : m_tuple{std::move(t_args)...} {}
    typed_argument_list(
        const typed_argument_list &other)
        : m_tuple{std::apply( [](auto const& ...arg){ return
std::tuple(dcli::clone(arg)...);},other)}{}

    typed_argument_list(
        typed_argument_list &&other)
        : m_tuple{std::move(other.m_tuple)} {}

    typed_argument_list(){}

    virtual typed_argument_list *clone() const override {return new
typed_argument_list(*this);}

    virtual std::string str() const override final {
        std::string out = std::string(Lexer::argument_list_start);
//        if (m_tuple.size() > 0) {
//            out += m_list[0]->str();
//            for (auto i = 1ul; i < m_list.size(); ++i) {
//                out += std::string(Lexer::argument_sep) + m_list[i]->str();
//            }
//        }
//        out += Lexer::argument_list_end;
        return out;
    }
};

template <class Lexer, class Compiler, class T>
class typed_assignment : public typed_expression<Lexer,Compiler,T> {
    std::unique_ptr<typed_identifier<Lexer,Compiler,T>> m_id;
    std::unique_ptr<typed_expression<Lexer,Compiler,T>> m_expr;

public:
    virtual ~typed_assignment() override {}
    typed_assignment(typed_identifier<Lexer,Compiler,T> *t_id,
                     typed_expression<Lexer,Compiler,T> *t_expr)
        : m_id{t_id}, m_expr{t_expr} {}

    typed_assignment() {}
    typed_assignment(typed_assignment &&other)
        : m_id{std::move(other.m_id)}, m_expr{std::move(other.m_expr)} {}
    typed_assignment(typed_assignment const &other)
        : m_id{other.m_id->clone()}, m_expr{other.m_expr->clone()} {}

    virtual std::string str() const override final {
        return m_id->str() + std::string(Lexer::assignment_operator) +
std::string(Lexer::assignment_sep) + m_expr->str();
    }
    virtual typed_assignment *clone() const override {
        return new typed_assignment(*this);
    };
};

template <class Lexer, class Compiler, class T, class F, class... Args>
    class typed_function_evaluation : public typed_expression<Lexer,Compiler,T>
{ std::unique_ptr<typed_identifier<Lexer,Compiler,F>> m_fid;
    std::unique_ptr<typed_argument_list<Lexer,Args...>> m_arg;

public:
    typed_function_evaluation(typed_identifier<Lexer,F> *t_fid,
                                typed_argument_list<Lexer,Args...> *t_arg)
        : m_fid{t_fid}, m_arg{t_arg }{}
    virtual ~typed_function_evaluation(){};
    typed_function_evaluation(typed_function_evaluation&& other):
        m_fid{std::move(other.m_fid)},m_arg{std::move(other.m_arg)}{}

    typed_function_evaluation(typed_function_evaluation const & other):
        m_fid{other.m_fid->clone()},m_arg{other.m_arg->clone()}{}

    virtual typed_function_evaluation* clone() const override {return new
typed_function_evaluation(*this);} virtual std::string str() const override
final { return m_fid->str() + m_arg->str();
    }


};
*/
} // namespace dcli

#endif // GRAMMAR_TYPED_H
