#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
//#include "grammar_typed.h"
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include "grammar_Identifier.h"
#include "literal_decode.h"
#include "type_name.h"
//#include "grammar_typed.h"
//#include "grammar_typed.h"
#include "maybe_error.h"
#include <macrodr/io/json/convert.h>
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

template <class Lexer, class Compiler, class T>
Maybe_error<void> load_literal_from_json_helper(
    const macrodr::io::json::Json& value, const std::string& path,
    macrodr::io::json::conv::TagPolicy policy, const Identifier<Lexer>& id,
    Environment<Lexer, Compiler>& env);

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
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<P, T>>)
class typed_predicate_evaluation;

template <class Lexer, class Compiler>
class base_Identifier_compiler {
   public:
    virtual ~base_Identifier_compiler() = default;
    [[nodiscard]] virtual base_Identifier_compiler* clone() const = 0;

    [[nodiscard]] virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const = 0;
};

template <class Lexer, class Compiler, class T>
class Identifier_compiler : public base_Identifier_compiler<Lexer, Compiler> {
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

   public:
    ~Identifier_compiler() override = default;
    [[nodiscard]] Identifier_compiler* clone() const override {
        return new Identifier_compiler(m_expr->clone());
    };

    Identifier_compiler(typed_expression<Lexer, Compiler, T>* t_expr) : m_expr{t_expr} {}

    [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const override {
        return new typed_identifier<Lexer, Compiler, T>(id);
    }
};

template <class Lexer, class Compiler>
class base_function_compiler {
   public:
    virtual ~base_function_compiler() = default;

    [[nodiscard]] virtual base_function_compiler* clone() const = 0;

    [[nodiscard]] virtual Maybe_unique<base_typed_expression<Lexer, Compiler>>
        compile_function_evaluation(Environment<Lexer, Compiler> const& cm,
                                    const untyped_argument_list<Lexer, Compiler>& args) const = 0;

    virtual void register_types(Compiler& registry) const = 0;
};

template <class Lexer, class Compiler>
class base_predicate_compiler {
   public:
    virtual ~base_predicate_compiler() = default;

    virtual base_predicate_compiler* clone() const = 0;

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_predicate_evaluation(
        Environment<Lexer, Compiler> const& cm,
        const untyped_expression<Lexer, Compiler>& expr) const = 0;
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
    }
    const auto expected = type_name<T>();
    const auto actual = expr ? type_name(*expr) : std::string{"<null>"};
    return error_message("\nunexpected type: expected ", expected, ", actual ", actual);
}

template <class T, class S, class... Ss>
Maybe_unique<typed_expression<Lexer, Compiler, T>> get_typed_expresion(
    std::unique_ptr<base_typed_expression<Lexer, Compiler>>& expr) {
    auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, S>*>(expr.get());
    if (ptr != nullptr) {
        return new typed_conversion<Lexer, Compiler, T, S>(
            dynamic_cast<typed_expression<Lexer, Compiler, S>*>(expr.release()));
    }
    return get_typed_expresion<T, Ss...>(expr);
}

template <class Lexer, class Compiler, class T>
class field_compiler {
    Identifier<Lexer> m_id;

   public:
    field_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    field_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto Maybe_id = cm.get_Identifier(t_arg());
        if (!Maybe_id) {
            return Maybe_id.error();
        }
        auto expr = std::move(Maybe_id.value());
        auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
        if (ptr != nullptr) {
            return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
        }
        const auto expected = type_name<T>();
        const auto actual = expr->type_name();
        return error_message(std::string("unexpected type: expected ") + expected + ", got " +
                             actual);
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value {
        auto maybe_value = from_literal<T>(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        return new typed_literal<Lexer, Compiler, T>(std::move(maybe_value.value()));
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value) {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto environment_local_copy = cm;
        auto Maybe_expr = t_arg.compile_expression(environment_local_copy);

        if (!Maybe_expr) {
            return Maybe_expr.error();
        }
        auto expr = std::move(Maybe_expr.value());
        auto exp_ptr = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(expr.get());
        if (exp_ptr != nullptr) {
            return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
        } else
            return error_message("type mismatch");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_expression<Lexer, Compiler> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg);
        if (ptr != nullptr) {
            return compile_this_argument(cm, *ptr);
        }
        auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
        if (li_ptr != nullptr)
            return compile_this_argument(cm, *li_ptr);
        else {
            auto f_ptr = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg);
            if (f_ptr != nullptr)
                return compile_this_argument(cm, *f_ptr);
            else
                return error_message("nj");
        }
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            bool const is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent) {
                return error_message("the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            }
            auto out = compile_this_argument(cm, ptr->expression());
            if (!out) {
                return error_message(std::string("\nargument ") + ptr->id()()() + ": \n" +
                                     out.error()());
            }
            return out;
        }
        auto ptr2 = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get());

        if (ptr2 != nullptr) {
            return compile_this_argument(cm, *ptr2);
        }
        return error_message(
            "if it is  neither an assigment nor an expression, what the fuck is it?");
    }
};

template <class Lexer, class Compiler, class T, class P>
    requires(std::is_same_v<std::invoke_result_t<P, T>, Maybe_error<T>>)
class field_compiler_precondition {
    Identifier<Lexer> m_id;
    P m_P;

   public:
    field_compiler_precondition(Identifier<Lexer>&& x, P&& p)
        : m_id(std::move(x)), m_P{std::move(p)} {}
    field_compiler_precondition(Identifier<Lexer> const& x, P const& p) : m_id(x), m_P{p} {}

    auto& id() const { return m_id; }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto Maybe_id = cm.get_Identifier(t_arg());
        if (!Maybe_id) {
            return Maybe_id.error();
        }
        auto expr = std::move(Maybe_id.value());
        auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
        if (ptr != nullptr) {
            return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
        }
        const auto expected = type_name<T>();
        const auto actual = expr ? type_name(*expr) : std::string{"<null>"};
        return error_message(std::string("fail to compile this argument type: "), expected,
                             ", got ", actual);
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const& /*unused*/, untyped_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value {
        auto maybe_value = from_literal<T>(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        return new typed_literal<Lexer, Compiler, T>(std::move(maybe_value.value()));
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const& /*unused*/, untyped_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value) {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto cmc = cm;
        auto Maybe_expr = t_arg.compile_expression(cmc);

        if (!Maybe_expr) {
            return Maybe_expr.error();
        }
        auto expr = std::move(Maybe_expr.value());
        auto exp_ptr = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(expr.get());
        if (exp_ptr != nullptr) {
            return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
        } else {
            const auto expected = type_name<T>();
            const auto actual = expr ? type_name(*expr) : std::string{"<null>"};
            return error_message("type mismatch: expected ", expected, ", got ", actual);
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_expression<Lexer, Compiler> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg);
        if (ptr != nullptr) {
            return compile_this_argument(cm, *ptr);
        }
        auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
        if (li_ptr != nullptr)
            return compile_this_argument(cm, *li_ptr);
        else {
            auto f_ptr = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg);
            if (f_ptr != nullptr)
                return compile_this_argument(cm, *f_ptr);
            else
                return error_message(
                    "unexpected node: expected identifier, literal, or function call");
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            bool is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent) {
                return error_message("the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            }
            return compile_this_argument(cm, ptr->expression());

        } else {
            auto ptr2 = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get());

            if (ptr2 != nullptr) {
                return compile_this_argument(cm, *ptr2);
            }
            return error_message("invalid argument: expected assignment or expression");
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
    [[nodiscard]] [[nodiscard]] [[nodiscard]] [[nodiscard]] [[nodiscard]] [[nodiscard]] [[nodiscard]] Maybe_unique<
        base_typed_expression<Lexer, Compiler>>
        compile_function_evaluation_impl(std::index_sequence<Is...> /*unused*/,
                                         Environment<Lexer, Compiler> const& cm,
                                         const untyped_argument_list<Lexer, Compiler>& args) const {
        if constexpr (sizeof...(Is) == 0) {
            // Zero-argument function: no argument compilation; pass empty tuple
            return new typed_function_evaluation<Lexer, Compiler, F, Args...>(m_f, std::tuple<>{});
        } else {
            auto Maybe_tuple = promote_Maybe_error(
                std::tuple(std::get<Is>(m_args).compile_this_argument(cm, args.arg()[Is])...));
            if (!Maybe_tuple) {
                return Maybe_tuple.error();
            }
            return new typed_function_evaluation<Lexer, Compiler, F, Args...>(
                m_f, std::move(Maybe_tuple.value()));
        }
    }

   public:
    function_compiler(F t_f, field_compiler<Lexer, Compiler, Args>&&... t_args)
        : m_args{std::move(t_args)...}, m_f{t_f} {}

    [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_function_evaluation(
        Environment<Lexer, Compiler> const& cm,
        const untyped_argument_list<Lexer, Compiler>& args) const override {
        return compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm, args.arg());
    }

    // base_function_compiler interface

    [[nodiscard]] base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }

    void register_types(Compiler& registry) const override {
        using Return = underlying_value_type_t<std::invoke_result_t<F, Args...>>;
        (registry.template ensure_type_registered<Args>(), ...);
        if constexpr (!std::is_void_v<Return>) {
            registry.template ensure_type_registered<Return>();
        }
    }
};

template <class Lexer, class Compiler, class F, class T>
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<F, T>>)
class predicate_compiler : public base_function_compiler<Lexer, Compiler> {
    field_compiler<Lexer, Compiler, T> m_arg;
    F m_f;

   public:
    predicate_compiler(F t_f, field_compiler<Lexer, Compiler, T>&& t_arg)
        : m_f{t_f}, m_arg{std::move(t_arg)} {}

    Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_predicate_evaluation(
        Environment<Lexer, Compiler> const& cm,
        const untyped_expression<Lexer, Compiler>& expr) const override {
        auto x = m_arg.compile_this_argument(cm, expr);
        if (!x) {
            return x.error();
        }
        return new typed_predicate_evaluation<Lexer, Compiler, F, T>(m_f, x.value().release());
    }

    // base_function_compiler interface

    base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }

    void register_types(Compiler& registry) const override {
        registry.template ensure_type_registered<T>();
    }
};

class Compiler {
    using Json = macrodr::io::json::Json;
    using TagPolicy = macrodr::io::json::conv::TagPolicy;

    struct TypeEntry {
        using LoadFn = Maybe_error<void> (*)(const Json&, const std::string&, TagPolicy,
                                             const Identifier<Lexer>&,
                                             Environment<Lexer, Compiler>&);
        LoadFn load = nullptr;
    };

    std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>> m_func;
    std::map<std::string, TypeEntry> m_type_registry;

   public:
    Compiler() { register_builtin_types(); }
    Compiler(std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>>&&
                 func)
        : m_func{std::move(func)} {
        register_builtin_types();
    }

    Compiler(const Compiler& cm)
        : m_func{clone_map(cm.m_func)}, m_type_registry{cm.m_type_registry} {}

    [[nodiscard]] Maybe_error<base_function_compiler<Lexer, Compiler> const*> get_function(
        const Identifier<Lexer>& id) const {
        auto it = m_func.find(id);
        if (it == m_func.end()) {
            return error_message(id() + " function is not defined");
        }
        auto ptr = (*it).second.get();
        if (ptr == nullptr)
            return error_message(id() + " function is null");
        else
            return ptr;
    }

    bool has_registered_type(const std::string& type_name) const {
        return m_type_registry.find(type_name) != m_type_registry.end();
    }

    Maybe_error<bool> push_function(std::string id_candidate,
                                    base_function_compiler<Lexer, Compiler>* fun) {
        auto may_id = to_Identifier<Lexer>(std::move(id_candidate));
        if (!may_id) {
            return may_id.error();
        }
        fun->register_types(*this);
        m_func.emplace(may_id.value(), fun);
        return true;
    }
    Maybe_error<bool> push_function(
        std::string id_candidate, std::unique_ptr<base_function_compiler<Lexer, Compiler>>&& fun) {
        auto may_id = to_Identifier<Lexer>(std::move(id_candidate));
        if (!may_id) {
            return may_id.error();
        }
        fun->register_types(*this);
        m_func.emplace(may_id.value(), std::move(fun));
        return true;
    }
    void merge(const Compiler& other) {
        // Merge functions
        for (const auto& [name, func] : other.m_func) {
            m_func[name] = std::unique_ptr<base_function_compiler<Lexer, Compiler>>(func->clone());
            if (auto it = m_func.find(name); it != m_func.end()) {
                it->second->register_types(*this);
            }
        }
        for (const auto& [key, entry] : other.m_type_registry) {
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, entry);
            }
        }
    }
    void merge(Compiler&& other) {
        // Move functions
        for (auto& [name, func] : std::move(other).m_func) {
            m_func[name] = std::move(func);  // Mueve el unique_ptr
            if (auto it = m_func.find(name); it != m_func.end() && it->second) {
                it->second->register_types(*this);
            }
        }
        //other.m_func.clear();  // Opcional: limpia el mapa fuente
        for (auto& [key, entry] : other.m_type_registry) {
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, entry);
            }
        }
    }

    // Introspection helper for environment persistence
    [[nodiscard]] std::vector<Identifier<Lexer>> list_functions() const {
        std::vector<Identifier<Lexer>> out;
        out.reserve(m_func.size());
        for (const auto& kv : m_func) out.push_back(kv.first);
        return out;
    }

  private:
    template <class Value>
    static TypeEntry make_entry() {
        return TypeEntry{
            &load_literal_from_json_helper<Lexer, Compiler, std::remove_cvref_t<Value>>};
    }

    void register_type_alias(const std::string& alias, const TypeEntry& entry) {
        if (!alias.empty() && m_type_registry.find(alias) == m_type_registry.end()) {
            m_type_registry.emplace(alias, entry);
        }
    }

    void register_builtin_types() {
        ensure_type_registered<double>();
        ensure_type_registered<int64_t>();
        ensure_type_registered<unsigned long>();
        ensure_type_registered<bool>();
        ensure_type_registered<std::string>();
        ensure_type_registered<Matrix<double>>();
        ensure_type_registered<DiagonalMatrix<double>>();

        const auto string_key = type_name<std::string>();
        auto string_it = m_type_registry.find(string_key);
        if (string_it != m_type_registry.end()) {
            register_type_alias("string", string_it->second);
        }

        const auto int_key = type_name<int64_t>();
        auto int_it = m_type_registry.find(int_key);
        if (int_it != m_type_registry.end()) {
            register_type_alias("int", int_it->second);
            register_type_alias("int64", int_it->second);
        }
    }

  public:
    template <class T>
    void ensure_type_registered() {
        using Value = std::remove_cvref_t<T>;
        if constexpr (std::is_void_v<Value>) {
            return;
        } else if constexpr (!macrodr::io::json::conv::has_json_codec_v<Value>) {
            return;
        } else {
            const std::string key = type_name<Value>();
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, make_entry<Value>());
            }
        }
    }

    Maybe_error<void> load_variable_from_json(const std::string& type_name, const Json& value,
                                              const std::string& path, TagPolicy policy,
                                              const Identifier<Lexer>& id,
                                              Environment<Lexer, Compiler>& env) const {
        auto it = m_type_registry.find(type_name);
        if (it == m_type_registry.end() || it->second.load == nullptr) {
            return error_message(path + ": unknown type '" + type_name + "'");
        }
        return it->second.load(value, path, policy, id, env);
    }
};

}  // namespace macrodr::dsl

#endif  // LEXER_TYPED_H
