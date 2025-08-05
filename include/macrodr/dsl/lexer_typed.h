#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
//#include "grammar_typed.h"
#include <iostream>
#include <map>
#include <memory>
#include <type_traits>

#include "grammar_Identifier.h"
//#include "grammar_typed.h"
//#include "grammar_typed.h"
#include "maybe_error.h"
namespace macrodr::dsl {

template <class Lexer, class Compiler>
class Environment;
template <class Lexer, class Compiler>
class base_typed_expression;
template <class Lexer, class Compiler>
class base_Identifier_compiler;
template <class Lexer, class Compiler, class T>
class typed_expression;

template <class Lexer, class Compiler, class T>
class typed_identifier;

template <class Lexer, class Compiler, class T>
class typed_literal;
template <class Lexer, class Compiler>
class untyped_argument_list;
class Lexer;
template <class Lexer, class Compiler>
class typed_argument_list;

template <class Lexer, class Compiler, class T, class S>
    requires(std::convertible_to<S, T>)
class typed_conversion;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_function_evaluation;

template <class Lexer, class Compiler>
class untyped_program;

template <class Lexer, class Compiler>
class typed_program;

template <class Lexer, class Compiler, class... T>
class typed_argument_typed_list;

template <class Lexer, class Compiler, class T>
class typed_expression;
template <class Lexer, class Compiler>
class untyped_expression;

template <class Lexer, class Compiler>
class untyped_statement;

template <class Lexer, class Compiler>
class untyped_numeric_literal;

template <class Lexer, class Compiler>
class untyped_identifier;

template <class Lexer, class Compiler>
class untyped_assignment;

template <class Lexer, class Compiler>
class untyped_literal;

template <class Lexer, class Compiler>
class untyped_function_evaluation;

template <class Lexer, class Compiler, class P, class T>
    requires(std::is_same_v<Maybe_error<T>,std::invoke_result_t<P, T>>)
class typed_predicate_evaluation;

template <class Lexer, class Compiler>
class base_Identifier_compiler {
   public:
    virtual ~base_Identifier_compiler() {};
    virtual base_Identifier_compiler* clone() const = 0;

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const = 0;
};

template <class Lexer, class Compiler, class T>
class Identifier_compiler : public base_Identifier_compiler<Lexer, Compiler> {
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

   public:
    virtual ~Identifier_compiler() {};
    virtual Identifier_compiler* clone() const {
        return new Identifier_compiler(m_expr->clone());
    };

    Identifier_compiler(typed_expression<Lexer, Compiler, T>* t_expr) : m_expr{t_expr} {
    }

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const {
        return new typed_identifier<Lexer, Compiler, T>(id);
    }
};

template <class Lexer, class Compiler>
class base_function_compiler {
   public:
    virtual ~base_function_compiler() {
    }

    virtual base_function_compiler* clone() const = 0;

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_function_evaluation(
        Environment<Lexer,Compiler> const& cm, const untyped_argument_list<Lexer, Compiler>& args) const = 0;
};

template <class Lexer, class Compiler>
class base_predicate_compiler {
   public:
    virtual ~base_predicate_compiler() {
    }
    
    virtual base_predicate_compiler* clone() const = 0;
    
        
    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_predicate_evaluation(
            Environment<Lexer,Compiler> const& cm, const untyped_expression<Lexer, Compiler>& expr) const =0;
    

};


template <class T, class... Sources>
struct Convertible_To;
template <class T>
struct Convertible_To<T> {};

class Compiler;
template <class T>
Maybe_unique<typed_expression<Lexer, Compiler, T>> get_typed_expresion(
    std::unique_ptr<base_typed_expression<Lexer, Compiler>>& expr) {
    auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
    if (ptr != nullptr) {
        return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
    } else {
        return error_message("unexpected type");
    }
}

template <class T, class S, class... Ss>
Maybe_unique<typed_expression<Lexer, Compiler, T>> get_typed_expresion(
    std::unique_ptr<base_typed_expression<Lexer, Compiler>>& expr) {
    auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, S>*>(expr.get());
    if (ptr != nullptr) {
        return new typed_conversion<Lexer, Compiler, T, S>(
            dynamic_cast<typed_expression<Lexer, Compiler, S>*>(expr.release()));
    } else {
        return get_typed_expresion<T, Ss...>(expr);
    }
}

template <class Lexer, class Compiler, class T>
class field_compiler {
    Identifier<Lexer> m_id;

   public:
    field_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {
    }
    field_compiler(Identifier<Lexer> const& x) : m_id(x) {
    }

    auto& id() const {
        return m_id;
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto Maybe_id = cm.get_Identifier(t_arg());
        if (!Maybe_id)
            return Maybe_id.error();
        else {
            auto expr = std::move(Maybe_id.value());
            auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
            if (ptr != nullptr) {
                return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
            } else {
                return error_message("unexpected type");
            }
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const&, untyped_literal<Lexer, Compiler> const& t_arg) const {
        auto Maybe_T = Lexer::template get<T>(t_arg());
        if (!Maybe_T)
            return Maybe_T.error();
        else {
            return new typed_literal<Lexer, Compiler, T>(std::move(Maybe_T.value()));
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto environment_local_copy = cm;
        auto Maybe_expr = t_arg.compile_expression(environment_local_copy);

        if (!Maybe_expr)
            return Maybe_expr.error();
        else {
            auto expr = std::move(Maybe_expr.value());
            auto exp_ptr = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(expr.get());
            if (exp_ptr != nullptr) {
                return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
            } else
                return error_message("type mismatch");
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_expression<Lexer, Compiler> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg);
        if (ptr != nullptr)
            return compile_this_argument(cm, *ptr);
        else {
            auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
            if (li_ptr != nullptr)
                return compile_this_argument(cm, *li_ptr);
            else {
                auto f_ptr =
                    dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg);
                if (f_ptr != nullptr)
                    return compile_this_argument(cm, *f_ptr);
                else
                    return error_message("nj");
            }
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            bool is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent)
                return error_message("the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            else {
                return compile_this_argument(cm, ptr->expression());
            }
        } else {
            auto ptr2 = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get());

            if (ptr2 != nullptr) {
                return compile_this_argument(cm, *ptr2);
            } else
                return error_message(
                    "if it is  neither an assigment nor an expression, what the fuck is it?");
        }
    }
};

template <class Lexer, class Compiler, class T, class P>
    requires(std::is_same_v<std::invoke_result_t<P, T>, Maybe_error<T>>)
class field_compiler_precondition {
    Identifier<Lexer> m_id;
    P m_P;

   public:
    field_compiler_precondition(Identifier<Lexer>&& x, P&& p)
        : m_id(std::move(x)), m_P{std::move(p)} {
    }
    field_compiler_precondition(Identifier<Lexer> const& x, P const& p) : m_id(x), m_P{p} {
    }

    auto& id() const {
        return m_id;
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto Maybe_id = cm.get_Identifier(t_arg());
        if (!Maybe_id)
            return Maybe_id.error();
        else {
            auto expr = std::move(Maybe_id.value());
            auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
            if (ptr != nullptr) {
                return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
            } else {
                return error_message("unexpected type");
            }
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const&, untyped_literal<Lexer, Compiler> const& t_arg) const {
        auto Maybe_T = Lexer::template get<T>(t_arg());
        if (!Maybe_T)
            return Maybe_T.error();
        else {
            return new typed_literal<Lexer, Compiler, T>(std::move(Maybe_T.value()));
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto cmc = cm;
        auto Maybe_expr = t_arg.compile_expression(cmc);

        if (!Maybe_expr)
            return Maybe_expr.error();
        else {
            auto expr = std::move(Maybe_expr.value());
            auto exp_ptr = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(expr.get());
            if (exp_ptr != nullptr) {
                return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
            } else
                return error_message("type mismatch");
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm, untyped_expression<Lexer, Compiler> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg);
        if (ptr != nullptr)
            return compile_this_argument(cm, *ptr);
        else {
            auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
            if (li_ptr != nullptr)
                return compile_this_argument(cm, *li_ptr);
            else {
                auto f_ptr =
                    dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg);
                if (f_ptr != nullptr)
                    return compile_this_argument(cm, *f_ptr);
                else
                    return error_message("nj");
            }
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer,Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            bool is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent)
                return error_message("the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            else {
                return compile_this_argument(cm, ptr->expression());
            }
        } else {
            auto ptr2 = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get());

            if (ptr2 != nullptr) {
                return compile_this_argument(cm, *ptr2);
            } else
                return error_message(
                    "if it is  neither an assigment nor an expression, what the fuck is it?");
        }
    }
};









template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_void_v<std::invoke_result_t<F, Args...>> ||
             std::is_object_v<std::invoke_result_t<F, Args...>>)
class function_compiler : public base_function_compiler<Lexer, Compiler> {
    std::tuple<field_compiler<Lexer, Compiler, Args>...> m_args;
    F m_f;

    template <std::size_t... Is>
    Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_function_evaluation_impl(
        std::index_sequence<Is...>, Environment<Lexer,Compiler> const& cm,
        const untyped_argument_list<Lexer, Compiler>& args) const {
        auto Maybe_tuple = promote_Maybe_error(
            std::tuple(std::get<Is>(m_args).compile_this_argument(cm, args.arg()[Is])...));
        if (!Maybe_tuple)
            return Maybe_tuple.error();
        else {
            return new typed_function_evaluation<Lexer, Compiler, F, Args...>(
                m_f, std::move(Maybe_tuple.value()));
        }
    }

   public:
    function_compiler(F t_f, field_compiler<Lexer, Compiler, Args>&&... t_args)
        : m_f{t_f}, m_args{std::move(t_args)...} {
    }

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_function_evaluation(
        Environment<Lexer,Compiler> const& cm, const untyped_argument_list<Lexer, Compiler>& args) const override {
        return compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm, args.arg());
    }

    // base_function_compiler interface
   public:
    virtual base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }
};


template <class Lexer, class Compiler, class F, class T>
    requires(std::is_same_v<Maybe_error<T>,std::invoke_result_t<F, T>> )
class predicate_compiler : public base_function_compiler<Lexer, Compiler> {
    field_compiler<Lexer, Compiler, T> m_arg;
    F m_f;
    
   public:
    predicate_compiler(F t_f, field_compiler<Lexer, Compiler, T>&& t_arg)
        : m_f{t_f}, m_arg{std::move(t_arg)} {
    }
    
    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_predicate_evaluation(
        Environment<Lexer,Compiler> const& cm, const untyped_expression<Lexer, Compiler>& expr) const override {
        auto x=m_arg.compile_this_argument(cm,expr);
        if (!x)
            return x.error();
        else
            return new typed_predicate_evaluation<Lexer,Compiler,F,T>(m_f, x.value().release());
    }
    
    // base_function_compiler interface
   public:
    virtual base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }
};



class Compiler {
    std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>> m_func;
  
   public:
    Compiler() {
    }
    Compiler(std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>>&&
                 func)
        : m_func{std::move(func)}{
    }

    Compiler(const Compiler& cm) : m_func{clone_map(cm.m_func)} {
    }

    Maybe_error<base_function_compiler<Lexer, Compiler> const*> get_function(
        const Identifier<Lexer>& id) const {
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

    
    Maybe_error<bool> push_function(std::string id_candidate,
                                    base_function_compiler<Lexer, Compiler>* fun) {
        auto may_id = to_Identifier<Lexer>(id_candidate);
        if (!may_id)
            return may_id.error();
        m_func.emplace(may_id.value(), fun);
        return true;
    }
    void merge(const Compiler& other) {
        // Merge functions
        for (const auto& [name, func] : other.m_func) {
            m_func[name] = std::unique_ptr<base_function_compiler<Lexer, Compiler>>(func->clone());
        }
    }
    void merge(Compiler&& other) {
        // Move functions
        for (auto& [name, func] : other.m_func) {
            m_func[name] = std::move(func);  // Mueve el unique_ptr
        }
        other.m_func.clear();  // Opcional: limpia el mapa fuente

    }
};

}  // namespace macrodr::dsl

#endif  // LEXER_TYPED_H
