#ifndef GRAMMAR_TYPED_H
#define GRAMMAR_TYPED_H
#include "grammar_Identifier.h"
//#include "grammar_untyped.h"
#include <concepts>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "deprecation.h"
#include "lexer_typed.h"
#include "maybe_error.h"
#include "type_name.h"
// Schemas
#include "parameters.h"

namespace macrodr::dsl {
// Forward declarations to break subtle include ordering issues
template <class Lexer, class Compiler>
class base_function_compiler;
template <class Lexer, class Compiler, class T>
class Identifier_compiler; // defined in lexer_typed.h
template <class Lexer, class Compiler, class T>
class field_compiler;
template <class Lexer, class Compiler>
class base_typed_statement;
template <class Lexer, class Compiler>
class base_typed_expression;
struct SerializedExpression {
    std::string type;
    macrodr::io::json::Json value;
};
template <class Lexer, class Compiler>
class base_Identifier_compiler;
template <class Lexer, class Compiler, class T>
class typed_expression;
// Forward declare typed_literal for earlier use in identifier_ref helpers
template <class Lexer, class Compiler, class T>
class typed_literal;

template <class Lexer, class Compiler>
class Environment {
    std::map<Identifier<Lexer>, std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>> m_id;
    std::map<Identifier<Lexer>, std::unique_ptr<base_typed_expression<Lexer, Compiler>>> m_var;
    Compiler const* cm_;
    // Identity map for parameter schemas (by id)
    std::map<std::string, std::shared_ptr<var::Parameters_Transformations>> m_param_schemas;

   public:
    Environment(Compiler const& cm) : cm_{&cm} {}

    Environment(const Environment& cm)
        : m_id{clone_map(cm.m_id)}, m_var{clone_map(cm.m_var)}, cm_{cm.cm_},
          m_param_schemas{cm.m_param_schemas} {}

    [[nodiscard]] Compiler const& compiler() const { return *cm_; }

    // Schema registry API
    void register_parameter_schema(const std::string& id,
                                   std::shared_ptr<var::Parameters_Transformations> schema) {
        if (!id.empty() && schema) m_param_schemas[id] = std::move(schema);
    }
    bool has_parameter_schema(const std::string& id) const {
        return m_param_schemas.find(id) != m_param_schemas.end();
    }
    std::shared_ptr<const var::Parameters_Transformations> get_parameter_schema(
        const std::string& id) const {
        auto it = m_param_schemas.find(id);
        if (it == m_param_schemas.end()) return {};
        return it->second;
    }
    std::vector<std::string> list_parameter_schema_ids() const {
        std::vector<std::string> ids;
        ids.reserve(m_param_schemas.size());
        for (auto& kv : m_param_schemas) ids.push_back(kv.first);
        return ids;
    }

    [[nodiscard]] Maybe_error<base_function_compiler<Lexer, Compiler> const*> get_function(
        const Identifier<Lexer>& id) const {
        return cm_->get_function(id);
    }
    [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>> get_Identifier(
        const Identifier<Lexer>& id) const {
        auto it = m_id.find(id);
        if (it == m_id.end()) {
            return error_message(std::string("\n") + id() + " function is not defined");
        }
        auto ptr = (*it).second.get();
        if (ptr == nullptr) {
            return error_message(std::string("\n") + id() + " function is null");
        }
        return ptr->compile_Identifier(id);
    }
    void push_back(Identifier<Lexer> id,
                   std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>> expr) {
        m_id.emplace(std::move(id), std::move(expr));
    }
    void push_back(Identifier<Lexer> id, base_Identifier_compiler<Lexer, Compiler>* expr) {
        push_back(std::move(id), std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>(expr));
    }
    [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const {
        auto Maybe_id = get_Identifier(id);
        if (!Maybe_id) {
            return Maybe_id.error();
        }
        return std::move(Maybe_id.value());
    }
    void insert(const Identifier<Lexer>& id,
                std::unique_ptr<base_typed_expression<Lexer, Compiler>> expr) {
        m_var.insert_or_assign(id, std::move(expr));
    }
    void insert(const Identifier<Lexer>& id, base_typed_expression<Lexer, Compiler>* expr) {
        insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(expr));
    }

    [[nodiscard]] Maybe_error<base_typed_expression<Lexer, Compiler> const*> get(
        Identifier<Lexer> const& id) const {
        auto it = m_var.find(id);
        if (it == m_var.end()) {
            return error_message(std::string("Identifier ") + id() + " not found");
        }
        return (*it).second.get();
    }

    // Introspection helpers for environment persistence
    [[nodiscard]] std::vector<Identifier<Lexer>> list_variables() const {
        std::vector<Identifier<Lexer>> out;
        out.reserve(m_var.size());
        for (const auto& kv : m_var) out.push_back(kv.first);
        return out;
    }

    [[nodiscard]] std::vector<Identifier<Lexer>> list_identifiers() const {
        std::vector<Identifier<Lexer>> out;
        out.reserve(m_id.size());
        for (const auto& kv : m_id) out.push_back(kv.first);
        return out;
    }

    void clear_variables() { m_var.clear(); }
    void clear_identifiers() { m_id.clear(); }
};

template <class Lexer, class Compiler>
class typed_argument_list;
template <class Lexer, class Compiler>
class base_Identifier_compiler;

template <class Lexer, class Compiler>
class base_typed_statement {
   public:
    virtual ~base_typed_statement() = default;
    [[nodiscard]] MACRODR_DEPRECATED(
        "Use clone_unique()") virtual base_typed_statement* clone() const = 0;
    [[nodiscard]] virtual std::unique_ptr<base_typed_statement> clone_unique() const {
        return std::unique_ptr<base_typed_statement>(this->clone());
    }
    virtual Maybe_error<typed_argument_list<Lexer, Compiler>*> compile_argument_list(
        typed_argument_list<Lexer, Compiler>* t_growing_list) const = 0;
    virtual Maybe_error<bool> run_statement(Environment<Lexer, Compiler>& env) const = 0;

    [[nodiscard]] virtual std::string type_name() const = 0;

    [[nodiscard]] MACRODR_DEPRECATED(
        "Use compile_identifier_unique()") virtual base_Identifier_compiler<Lexer,
                                                                            Compiler>* compile_identifier()
        const = 0;
    [[nodiscard]] virtual std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const {
        return std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>(
            this->compile_identifier());
    }
};

template <class Lexer, class Compiler>
class base_typed_assigment : public base_typed_statement<Lexer, Compiler> {
   public:
    ~base_typed_assigment() override = default;
    [[nodiscard]] MACRODR_DEPRECATED("Use clone_unique()")
        base_typed_assigment* clone() const override = 0;
    [[nodiscard]] virtual Identifier<Lexer> const& id() const = 0;
    [[nodiscard]] virtual base_typed_expression<Lexer, Compiler>* expr() const = 0;

    Maybe_error<typed_argument_list<Lexer, Compiler>*> compile_argument_list(
        typed_argument_list<Lexer, Compiler>* t_growing_list) const override;
};

template <class Lexer, class Compiler>
class base_typed_expression : public base_typed_statement<Lexer, Compiler> {
   public:
    ~base_typed_expression() override = default;
    [[nodiscard]] MACRODR_DEPRECATED("Use clone_unique()")
        base_typed_expression* clone() const override = 0;

    using base_typed_statement<Lexer, Compiler>::compile_identifier_unique;

    [[nodiscard]] MACRODR_DEPRECATED(
        "Use compile_identifier_unique()") virtual base_typed_expression<Lexer,
                                                                         Compiler>* compile_identifier(Identifier<Lexer>
                                                                                                           id)
        const = 0;

    [[nodiscard]] virtual std::unique_ptr<base_typed_expression<Lexer, Compiler>>
        compile_identifier_unique(Identifier<Lexer> id) const {
        return std::unique_ptr<base_typed_expression<Lexer, Compiler>>(
            this->compile_identifier(id));
    }

    MACRODR_DEPRECATED("Use compile_assigment_unique()")
    virtual base_typed_assigment<Lexer, Compiler>* compile_assigment(Identifier<Lexer> id) = 0;

    virtual std::unique_ptr<base_typed_assigment<Lexer, Compiler>> compile_assigment_unique(
        Identifier<Lexer> id) {
        return std::unique_ptr<base_typed_assigment<Lexer, Compiler>>(this->compile_assigment(id));
    }

    MACRODR_DEPRECATED("Use compile_argument_list_unique() where applicable")
    Maybe_error<typed_argument_list<Lexer, Compiler>*> compile_argument_list(
        typed_argument_list<Lexer, Compiler>* t_growing_list) const override;

    MACRODR_DEPRECATED("Use run_expression_unique()")
    virtual Maybe_error<base_typed_expression<Lexer, Compiler>*> run_expression(
        Environment<Lexer, Compiler>& env) const = 0;

    [[nodiscard]] virtual Maybe_unique<base_typed_expression<Lexer, Compiler>>
        run_expression_unique(Environment<Lexer, Compiler>& env) const {
        auto maybe_expr = run_expression(env);
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return std::unique_ptr<base_typed_expression<Lexer, Compiler>>(maybe_expr.value());
    }

    Maybe_error<bool> run_statement(Environment<Lexer, Compiler>& env) const override {
        auto Maybe_exp = run_expression(env);
        if (!Maybe_exp) {
            return Maybe_exp.error();
        }
        return true;
    }

    [[nodiscard]] base_Identifier_compiler<Lexer, Compiler>* compile_identifier() const override =
        0;

    virtual Maybe_error<SerializedExpression> serialize_json(
        const Environment<Lexer, Compiler>& env,
        typename json_spec<Lexer>::TagPolicy policy) const {
        return error_message(std::string{"serialization unavailable for type "} +
                             this->type_name());
    }
};

template <class Lexer, class Compiler>
class typed_argument_list : public base_typed_expression<Lexer, Compiler> {
    std::vector<std::unique_ptr<base_typed_expression<Lexer, Compiler>>> m_vec;
    std::map<Identifier<Lexer>, std::unique_ptr<base_typed_expression<Lexer, Compiler>>> m_map;

   public:
    [[nodiscard]] std::string type_name() const override { return "argument_list"; }
    ~typed_argument_list() override = default;
    [[nodiscard]] typed_argument_list* clone() const override {
        return new typed_argument_list(*this);
    }
    typed_argument_list(const typed_argument_list& other)
        : m_vec(clone_vector(other.m_vec)), m_map{clone_map(other.m_map)} {}
    typed_argument_list(typed_argument_list&&) noexcept = default;
    typed_argument_list& operator=(typed_argument_list&&) noexcept = default;
    typed_argument_list& operator=(const typed_argument_list& other) {
        if (this != &other) {
            m_vec = clone_vector(other.m_vec);
            m_map = clone_map(other.m_map);
        }
        return *this;
    }
    [[nodiscard]] auto& arg_vector() const { return m_vec; }
    [[nodiscard]] auto& arg_map() const { return m_map; }
    void insert_assigment(base_typed_assigment<Lexer, Compiler> const& ass) {
        m_map[ass.id()] =
            std::unique_ptr<base_typed_expression<Lexer, Compiler>>(ass.expr()->clone());
    }

    Maybe_error<bool> push_back_expression(base_typed_expression<Lexer, Compiler> const& expr) {
        if (!m_map.empty()) {
            return error_message("assigment after expression in argument list");
        }
        m_vec.emplace_back(expr.clone());
        return true;
    }

    typed_argument_list() = default;
    base_typed_assigment<Lexer, Compiler>* compile_assigment(Identifier<Lexer> /*id*/) override {
        return nullptr;
    }

    [[nodiscard]] base_Identifier_compiler<Lexer, Compiler>* compile_identifier() const override {
        return nullptr;
    }

    [[nodiscard]] base_typed_expression<Lexer, Compiler>* compile_identifier(
        Identifier<Lexer> /*id*/) const override {
        return nullptr;
    };

    // base_typed_expression interface

    Maybe_error<dsl::base_typed_expression<Lexer, Compiler>*> run_expression(
        dsl::Environment<Lexer, Compiler>& /*env*/) const override {
        return error_message("no reason to get here");
    }
};

template <class Lexer, class Compiler, class... T>
class typed_argument_typed_list : public base_typed_expression<Lexer, Compiler> {
    std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, T>>...> m_args;

   public:
    std::string type_name() const override { return "typed_argument_list"; }
    virtual ~typed_argument_typed_list() = default;
    virtual typed_argument_typed_list* clone() const override {
        return new typed_argument_typed_list(clone_tuple(m_args));
    }
    typed_argument_typed_list(
        std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, T>>...>&& t)
        : m_args{std::move(t)} {}
    auto& args() const { return m_args; }

    virtual base_Identifier_compiler<Lexer, Compiler>* compile_identifier() const override {
        return nullptr;
    };

    virtual base_typed_expression<Lexer, Compiler>* compile_identifier(
        Identifier<Lexer> /*id*/) const override {
        return nullptr;
    };

    virtual base_typed_assigment<Lexer, Compiler>* compile_assigment(
        Identifier<Lexer> /*id*/) override {
        return nullptr;
    };

    // base_typed_statement interface

    // base_typed_expression interface

    Maybe_error<dsl::base_typed_expression<Lexer, Compiler>*> run_expression(
        dsl::Environment<Lexer, Compiler>& /*env*/) const override {
        return error_message("still unresolved");
    }
};

template <class Lexer, class Compiler, class T>
class typed_assigment;
template <class Lexer, class Compiler, class T>
class typed_identifier : public typed_expression<Lexer, Compiler, T> {
    Identifier<Lexer> m_id;

   public:
    typed_identifier(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    typed_identifier(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }
    ~typed_identifier() override = default;
    [[nodiscard]] typed_identifier* clone() const override { return new typed_identifier(*this); }

    [[nodiscard]] Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const override {
        auto May_x = env.get(m_id);
        if (!May_x) {
            return May_x.error();
        }
        auto exp = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(May_x.value());
        if (!exp) {
            // TODO: we can make this message tell both types in the mismatch
            return error_message(std::string("type mismatch for identifier ") + m_id());
        }
        return exp->run(env);
    }
};

template <class Lexer, class Compiler, class T>
class typed_identifier_ref_const
    : public typed_expression<Lexer, Compiler, std::reference_wrapper<const T>> {
    Identifier<Lexer> m_id;

   public:
    explicit typed_identifier_ref_const(Identifier<Lexer> id) : m_id(std::move(id)) {}

    typed_identifier_ref_const* clone() const override {
        return new typed_identifier_ref_const(*this);
    }

    Maybe_error<std::reference_wrapper<const T>> run(
        const Environment<Lexer, Compiler>& env) const override {
        auto May_x = env.get(m_id);
        if (!May_x) {
            return May_x.error();
        }
        auto literal =
            dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(May_x.value());
        if (!literal) {
            auto actual = type_name(*May_x.value());
            return error_message(std::string("identifier '") + m_id() +
                                 "' is not bound to the expected reference type (actual: " +
                                 actual + ")");
        }
        return std::cref(literal->value_ref());
    }
};

template <class Lexer, class Compiler, class T>
class typed_identifier_ref
    : public typed_expression<Lexer, Compiler, std::reference_wrapper<T>> {
    Identifier<Lexer> m_id;

   public:
    explicit typed_identifier_ref(Identifier<Lexer> id) : m_id(std::move(id)) {}

    typed_identifier_ref* clone() const override { return new typed_identifier_ref(*this); }

    Maybe_error<std::reference_wrapper<T>> run(
        const Environment<Lexer, Compiler>& env) const override {
        auto May_x = env.get(m_id);
        if (!May_x) {
            return May_x.error();
        }
        auto literal_const =
            dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(May_x.value());
        if (!literal_const) {
            auto actual = type_name(*May_x.value());
            return error_message(std::string("identifier '") + m_id() +
                                 "' is not bound to the expected reference type (actual: " +
                                 actual + ")");
        }
        auto literal = const_cast<typed_literal<Lexer, Compiler, T>*>(literal_const);
        return std::ref(literal->value_ref());
    }
};

template <class Lexer, class Compiler, class T>
class typed_expression : public base_typed_expression<Lexer, Compiler> {
   public:
    [[nodiscard]] std::string type_name() const override { return macrodr::dsl::type_name<T>(); }

    ~typed_expression() override = default;
    // Complete special members to satisfy rule-of-five guidance for abstract base
    typed_expression() = default;
    typed_expression(const typed_expression&) = default;
    typed_expression& operator=(const typed_expression&) = default;
    typed_expression(typed_expression&&) noexcept = default;
    typed_expression& operator=(typed_expression&&) noexcept = default;
    [[nodiscard]] typed_expression* clone() const override = 0;

    [[nodiscard]] virtual Maybe_error<T> run(Environment<Lexer, Compiler> const&) const = 0;

    base_typed_assigment<Lexer, Compiler>* compile_assigment(Identifier<Lexer> id) override;

    [[nodiscard]] base_typed_expression<Lexer, Compiler>* compile_identifier(
        Identifier<Lexer> id) const override {
        return new typed_identifier<Lexer, Compiler, T>(id);
    }

    [[nodiscard]] base_Identifier_compiler<Lexer, Compiler>* compile_identifier() const override {
        return new Identifier_compiler<Lexer, Compiler, T>(this->clone());
    }

    // base_typed_expression interface

    Maybe_error<dsl::base_typed_expression<Lexer, Compiler>*> run_expression(
        dsl::Environment<Lexer, Compiler>& env) const override {
        if constexpr (std::is_void_v<T>) {
            run(env);
            return new typed_literal<Lexer, Compiler, T>();
        } else {
            auto Maybe_x = run(env);
            if (!Maybe_x) {
                return Maybe_x.error();
            }
            return new typed_literal<Lexer, Compiler, T>(std::move(Maybe_x.value()));
        }
    }

    Maybe_error<SerializedExpression> serialize_json(
        const Environment<Lexer, Compiler>& env,
        typename json_spec<Lexer>::TagPolicy policy) const override {
        using ValueType = std::remove_cvref_t<T>;
        if constexpr (std::is_void_v<ValueType>) {
            return error_message(std::string{"cannot serialize void expression of type "} +
                                 this->type_name());
        } else if constexpr (!macrodr::io::json::conv::has_json_codec_v<ValueType>) {
            return error_message(std::string{"no JSON codec registered for type "} +
                                 this->type_name());
        } else {
            auto maybe_value = this->run(env);
            if (!maybe_value) {
                return maybe_value.error();
            }
            SerializedExpression out;
            out.type = this->type_name();
            out.value = json_spec<Lexer>::to_json(maybe_value.value(), policy);
            return out;
        }
    }
};

template <class Lexer, class Compiler, class T>
class typed_assigment : public base_typed_assigment<Lexer, Compiler> {
    Identifier<Lexer> m_id;
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

   public:
    [[nodiscard]] std::string type_name() const override { return macrodr::dsl::type_name<T>(); }

    ~typed_assigment() override = default;
    [[nodiscard]] typed_assigment* clone() const override { return new typed_assigment(*this); };

    typed_assigment(const typed_assigment& other)
        : m_id{other.m_id}, m_expr{other.m_expr->clone()} {}

    typed_assigment(typed_assigment&&) noexcept = default;
    typed_assigment& operator=(typed_assigment&&) noexcept = default;
    typed_assigment& operator=(const typed_assigment& other) {
        if (this != &other) {
            m_id = other.m_id;
            m_expr.reset(other.m_expr ? other.m_expr->clone() : nullptr);
        }
        return *this;
    }

    typed_assigment(Identifier<Lexer> const& t_id, typed_expression<Lexer, Compiler, T>* t_expr)
        : m_id{t_id}, m_expr{t_expr} {}

    [[nodiscard]] Identifier<Lexer> const& id() const override { return m_id; }
    [[nodiscard]] typed_expression<Lexer, Compiler, T>* expr() const override {
        return m_expr.get();
    }

    Maybe_error<bool> run_statement(Environment<Lexer, Compiler>& env) const override {
        auto maybe_expr = expr()->run_expression_unique(env);
        if (!maybe_expr) {
            return error_message(std::string("\n in assignment to ") + id()() + " : " +
                                 maybe_expr.error()());
        }
        env.insert(id(), std::move(maybe_expr.value()));
        return true;
    }

    [[nodiscard]] base_Identifier_compiler<Lexer, Compiler>* compile_identifier() const override {
        return new Identifier_compiler<Lexer, Compiler, T>(this->expr());
    }
};

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_function_evaluation
    : public typed_expression<Lexer, Compiler,
                              underlying_value_type_t<std::invoke_result_t<F, Args...>>> {
    std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, Args>>...> m_args;
    F m_f;

   public:
    using T = underlying_value_type_t<std::invoke_result_t<F, Args...>>;
    typed_function_evaluation(F t_f, typed_expression<Lexer, Compiler, Args>*... t_args)
        : m_f{t_f}, m_args{t_args...} {}

    typed_function_evaluation(
        F t_f, std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, Args>>...>&& t)
        : m_args{std::move(t)}, m_f{t_f} {}

    // typed_expression interface

    typed_function_evaluation(const typed_function_evaluation& other)
        : m_args{std::apply(
              [](auto const&... args) {
                  return std::tuple(
                      std::unique_ptr<typed_expression<Lexer, Compiler, Args>>(args->clone())...);
              },
              other.m_args)},
          m_f{other.m_f} {}

    typed_function_evaluation(typed_function_evaluation&&) noexcept = default;
    typed_function_evaluation& operator=(typed_function_evaluation&&) noexcept = default;
    typed_function_evaluation& operator=(const typed_function_evaluation& other) {
        if (this != &other) {
            auto cloned = std::apply(
                [](auto const&... args) {
                    return std::tuple(
                        std::unique_ptr<typed_expression<Lexer, Compiler, Args>>(args->clone())...);
                },
                other.m_args);
            m_args = std::move(cloned);
            m_f = other.m_f;
        }
        return *this;
    }

    [[nodiscard]] typed_function_evaluation* clone() const override {
        return new typed_function_evaluation(*this);
    }

    [[nodiscard]] Maybe_error<T> run(const Environment<Lexer, Compiler>& env) const override {
        if constexpr (sizeof...(Args) == 0) {
            if constexpr (std::is_void_v<T>) {
                this->m_f();
                return Maybe_error<void>();
            } else {
                return Maybe_error<T>(this->m_f());
            }
        } else {
            auto Maybe_arg_tuple = std::apply(
                [this, &env](auto&... args) { return std::tuple(args->run(env)...); }, m_args);
            return std::apply(
                [this](auto&&... Maybe_args) -> Maybe_error<T> {
                    if ((Maybe_args.valid() && ...)) {
                        if constexpr (std::is_void_v<T>) {
                            this->m_f(std::move(Maybe_args.value())...);
                            return Maybe_error<void>();
                        } else {
                            return Maybe_error<T>(this->m_f(std::move(Maybe_args.value())...));
                        }
                    } else {
                        return error_message(((std::string{}) + ... + Maybe_args.error()()));
                    }
                },
                Maybe_arg_tuple);
        }
    }
};

template <class Lexer, class Compiler, class P, class T>
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<P, T>>)
class typed_predicate_evaluation : public typed_expression<Lexer, Compiler, T> {
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_arg;
    P m_P;

   public:
    typed_predicate_evaluation(P t_f, typed_expression<Lexer, Compiler, T>* t_arg)
        : m_P{t_f}, m_arg{t_arg} {}

    typed_predicate_evaluation(P t_f, std::unique_ptr<typed_expression<Lexer, Compiler, T>>&& t)
        : m_P{t_f}, m_arg{std::move(t)} {}

    // typed_expression interface

    typed_predicate_evaluation(const typed_predicate_evaluation& other)
        : m_arg{other.m_arg->clone()}, m_P{other.m_P} {}

    typed_predicate_evaluation(typed_predicate_evaluation&&) noexcept = default;
    typed_predicate_evaluation& operator=(typed_predicate_evaluation&&) noexcept = default;
    typed_predicate_evaluation& operator=(const typed_predicate_evaluation& other) {
        if (this != &other) {
            m_arg.reset(other.m_arg ? other.m_arg->clone() : nullptr);
            m_P = other.m_P;
        }
        return *this;
    }

    virtual typed_predicate_evaluation* clone() const override {
        return new typed_predicate_evaluation(*this);
    }

    Maybe_error<T> run(const Environment<Lexer, Compiler>& env) const override {
        // return T{};
        auto out = m_arg->run(env);
        if (!out) {
            return out.error();
        }
        return m_P(std::move(out.value()));
    }
};


template <class Lexer, class Compiler, class T>
class typed_vector_construction
    : public typed_expression<Lexer, Compiler, std::vector<T>> {

    // Internal holder type (same rule as for function/tuple arguments)
    template <class U>
    using storage_t = detail::function_argument_storage_t<U>;

    // Expressions producing storage_t<T> elements
    std::vector<std::unique_ptr<typed_expression<Lexer, Compiler, storage_t<T>>>> m_args;

    // Same adaptation logic as function_compiler / typed_tuple_construction
    template <class Param, class Storage>
    static decltype(auto) adapt(Storage& storage) {
        using Base        = std::remove_reference_t<Param>;
        using StorageType = std::remove_reference_t<Storage>;

        if constexpr (std::is_same_v<StorageType, std::reference_wrapper<Base>> ||
                      std::is_same_v<StorageType, std::reference_wrapper<std::remove_const_t<Base>>>) {
            return static_cast<Param>(storage.get());
        } else if constexpr (std::is_lvalue_reference_v<Param>) {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return static_cast<Param>(*storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return static_cast<Param>(*storage);
            } else {
                if constexpr (std::is_const_v<Base>) {
                    return static_cast<const Base&>(storage);
                } else {
                    return static_cast<Base&>(storage);
                }
            }
        } else if constexpr (std::is_rvalue_reference_v<Param>) {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return static_cast<Base&&>(*storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return static_cast<Base&&>(*storage);
            } else {
                return static_cast<Base&&>(storage);
            }
        } else {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return std::move(storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return storage;
            } else {
                return static_cast<Param>(std::move(storage));
            }
        }
    }

   public:
    using storage_args =
        std::vector<std::unique_ptr<typed_expression<Lexer, Compiler, storage_t<T>>>>;

    explicit typed_vector_construction(storage_args&& args)
        : m_args{std::move(args)} {}

    typed_vector_construction(storage_args const& args)
        : m_args{clone_vector(args)} {}

    ~typed_vector_construction() override = default;

    typed_vector_construction(const typed_vector_construction& other)
        : m_args{clone_vector(other.m_args)} {}

    typed_vector_construction(typed_vector_construction&&) noexcept = default;
    typed_vector_construction& operator=(typed_vector_construction&&) noexcept = default;

    typed_vector_construction& operator=(const typed_vector_construction& other) {
        if (this != &other) {
            m_args = clone_vector(other.m_args);
        }
        return *this;
    }

    [[nodiscard]] typed_vector_construction* clone() const override {
        return new typed_vector_construction(*this);
    }

    [[nodiscard]] Maybe_error<std::vector<T>>
    run(const Environment<Lexer, Compiler>& env) const override {
        std::vector<T> out;
        out.reserve(m_args.size());
        std::string err;

        for (std::size_t i = 0; i < m_args.size(); ++i) {
            auto maybe_elem = m_args[i]->run(env);  // Maybe_error<storage_t<T>>
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                out.emplace_back(adapt<T>(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
        }
        return out;
    }
};


template <class Lexer, class Compiler, class... Ts>
class typed_tuple_construction
    : public typed_expression<Lexer, Compiler, std::tuple<Ts...>> {

    // Internal holder type (same rule as function_compiler)
    template <class U>
    using storage_t = detail::function_argument_storage_t<U>;

    // Store expressions producing storage_t<Ts>...
    std::tuple<std::unique_ptr<typed_expression<Lexer, Compiler, storage_t<Ts>>>...> m_args;

    // Same adaptation logic as in function_compiler
    template <class Param, class Storage>
    static decltype(auto) adapt(Storage& storage) {
        using Base        = std::remove_reference_t<Param>;
        using StorageType = std::remove_reference_t<Storage>;

        if constexpr (std::is_same_v<StorageType, std::reference_wrapper<Base>> ||
                      std::is_same_v<StorageType, std::reference_wrapper<std::remove_const_t<Base>>>) {
            return static_cast<Param>(storage.get());
        } else if constexpr (std::is_lvalue_reference_v<Param>) {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return static_cast<Param>(*storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return static_cast<Param>(*storage);
            } else {
                if constexpr (std::is_const_v<Base>) {
                    return static_cast<const Base&>(storage);
                } else {
                    return static_cast<Base&>(storage);
                }
            }
        } else if constexpr (std::is_rvalue_reference_v<Param>) {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return static_cast<Base&&>(*storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return static_cast<Base&&>(*storage);
            } else {
                return static_cast<Base&&>(storage);
            }
        } else {
            if constexpr (detail::is_unique_ptr<StorageType>::value) {
                return std::move(storage);
            } else if constexpr (std::is_pointer_v<StorageType>) {
                return storage;
            } else {
                return static_cast<Param>(std::move(storage));
            }
        }
    }

   public:
    using tuple_type   = std::tuple<Ts...>;
    using storage_args = std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,storage_t<Ts>>>...>;

    explicit typed_tuple_construction(storage_args&& t)
        : m_args{std::move(t)} {}

    ~typed_tuple_construction() override = default;

    typed_tuple_construction(const typed_tuple_construction& other)
        : m_args{ clone_tuple(other.m_args) } {}  // <── use of clone_tuple

    [[nodiscard]] typed_tuple_construction* clone() const override {
        return new typed_tuple_construction(*this);
    }

    [[nodiscard]] Maybe_error<std::tuple<Ts...>>
    run(const Environment<Lexer, Compiler>& env) const override {

        // Evaluate each element → Maybe_error<storage_t<T>>
        auto maybe_storage =
            std::apply([&](auto&... exprs) {
                return std::tuple(exprs->run(env)...);
            }, m_args);

        bool ok = true;
        std::string msg;

        // Aggregate errors
        std::apply([&](auto&... me) {
            (([&]{
                if (!me.valid()) {
                    ok = false;
                    msg += me.error()();
                }
            }()), ...);
        }, maybe_storage);

        if (!ok) return error_message(msg);

        // Build the real tuple<Ts...>
        auto builder = [&](auto&... me) {
            return std::tuple<Ts...>( adapt<Ts>(me.value())... );
        };

        return std::apply(builder, maybe_storage);
    }
};



template <class Lexer, class Compiler, class T>
base_typed_assigment<Lexer, Compiler>* typed_expression<Lexer, Compiler, T>::compile_assigment(
    Identifier<Lexer> id) {
    return new typed_assigment<Lexer, Compiler, T>(id, this->clone());
}

template <class Lexer, class Compiler, class T, class S>
    requires(std::convertible_to<S, T>)
class typed_conversion : public typed_expression<Lexer, Compiler, T> {
    std::unique_ptr<typed_expression<Lexer, Compiler, S>> m_expr;

   public:
    typed_conversion(std::unique_ptr<typed_expression<Lexer, Compiler, S>>&& x)
        : m_expr(std::move(x)) {}
    typed_conversion(std::unique_ptr<typed_expression<Lexer, Compiler, S>> const& x)
        : m_expr(x->clone()) {}
    typed_conversion(typed_expression<Lexer, Compiler, S>* x) : m_expr(x) {}

    virtual ~typed_conversion() = default;

    Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const override {
        auto r = m_expr->run(env);
        if (!r) {
            return r.error();
        }
        return static_cast<T>(std::move(r.value()));
    }

    // typed_expression interface

    typed_conversion* clone() const override { return new typed_conversion(*this); }
};

template <class Lexer, class Compiler, class T>
    requires(!std::is_void_v<T> && std::is_copy_constructible_v<T>)
class typed_literal<Lexer, Compiler, T> : public typed_expression<Lexer, Compiler, T> {
    T m_value;

   public:
    typed_literal(T&& x) : m_value(std::move(x)) {}
    typed_literal(T const& x) : m_value(x) {}

    ~typed_literal() override = default;

    const T& value_ref() const { return m_value; }
    T& value_ref() { return m_value; }

    [[nodiscard]] Maybe_error<T> run(
        Environment<Lexer, Compiler> const& /*unused*/) const override {
        return m_value;
    }

    // typed_expression interface

    [[nodiscard]] typed_literal* clone() const override { return new typed_literal(*this); }
};

template <class Lexer, class Compiler, class T>
    requires requires(const T& value) {
        { value.clone() } -> std::convertible_to<std::unique_ptr<T>>;
    }
class typed_literal<Lexer, Compiler, std::unique_ptr<T>>
    : public typed_expression<Lexer, Compiler, std::unique_ptr<T>> {
    std::unique_ptr<T> m_value;

    static std::unique_ptr<T> clone_ptr(const std::unique_ptr<T>& ptr) {
        return ptr ? ptr->clone() : std::unique_ptr<T>{};
    }

   public:
    explicit typed_literal(std::unique_ptr<T>&& x) : m_value(std::move(x)) {}
    explicit typed_literal(const std::unique_ptr<T>& x) : m_value(clone_ptr(x)) {}

    const std::unique_ptr<T>& value_ref() const { return m_value; }
    std::unique_ptr<T>& value_ref() { return m_value; }

    typed_literal(const typed_literal& other) : m_value(clone_ptr(other.m_value)) {}
    typed_literal& operator=(const typed_literal& other) {
        if (this != &other)
            m_value = clone_ptr(other.m_value);
        return *this;
    }
    typed_literal(typed_literal&&) noexcept = default;
    typed_literal& operator=(typed_literal&&) noexcept = default;

    [[nodiscard]] Maybe_error<std::unique_ptr<T>> run(
        Environment<Lexer, Compiler> const&) const override {
        return clone_ptr(m_value);
    }

    [[nodiscard]] typed_literal* clone() const override { return new typed_literal(*this); }
};

template <class Lexer, class Compiler, class T>
    requires(!std::is_void_v<T> && !std::is_copy_constructible_v<T> &&
             std::is_move_constructible_v<T> &&
             requires(const T& value) {
                 { value.clone() } -> std::same_as<std::unique_ptr<T>>;
             })
class typed_literal<Lexer, Compiler, Maybe_error<std::unique_ptr<T>>>
    : public typed_expression<Lexer, Compiler, Maybe_error<std::unique_ptr<T>>> {
    using value_type = Maybe_error<std::unique_ptr<T>>;

    value_type m_value;

    static value_type clone_value(const value_type& source) {
        if (!source) {
            return source.error();
        }
        const auto& ptr = source.value();
        if (!ptr) {
            return std::unique_ptr<T>{};
        }
        return ptr->clone();
    }

   public:
    explicit typed_literal(value_type&& x) : m_value(std::move(x)) {}
    explicit typed_literal(const value_type& x) : m_value(clone_value(x)) {}

    const value_type& value_ref() const { return m_value; }
    value_type& value_ref() { return m_value; }

    typed_literal(const typed_literal& other) : m_value(clone_value(other.m_value)) {}
    typed_literal& operator=(const typed_literal& other) {
        if (this != &other) {
            m_value = clone_value(other.m_value);
        }
        return *this;
    }
    typed_literal(typed_literal&&) noexcept = default;
    typed_literal& operator=(typed_literal&&) noexcept = default;

    [[nodiscard]] Maybe_error<value_type> run(Environment<Lexer, Compiler> const&) const override {
        return clone_value(m_value);
    }

    [[nodiscard]] typed_literal* clone() const override { return new typed_literal(*this); }
};

template <class Lexer, class Compiler>
class typed_literal<Lexer, Compiler, void> : public typed_expression<Lexer, Compiler, void> {
   public:
    ~typed_literal() override = default;

    [[nodiscard]] Maybe_error<void> run(
        Environment<Lexer, Compiler> const& /*unused*/) const override {
        return {};
    }

    // typed_expression interface

    [[nodiscard]] typed_literal* clone() const override { return new typed_literal(*this); }
};

template <class Lexer, class Compiler, class T>
Maybe_error<void> load_literal_from_json_helper(
    const typename json_spec<Lexer>::Json& value, const std::string& path,
    typename json_spec<Lexer>::TagPolicy policy, const Identifier<Lexer>& id,
    Environment<Lexer, Compiler>& env) {
    using Value = std::remove_cvref_t<T>;
    Value decoded{};
    auto status = json_spec<Lexer>::from_json(value, decoded, path, policy);
    if (!status) {
        return status;
    }
    auto literal = std::make_unique<typed_literal<Lexer, Compiler, Value>>(std::move(decoded));
    env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(literal->clone()));
    env.push_back(id, new Identifier_compiler<Lexer, Compiler, Value>(literal.release()));
    return {};
}

// Overload for Parameters_values: resolve schema_id in Environment
template <class Lexer, class Compiler, class T,
          std::enable_if_t<std::is_same_v<std::remove_cvref_t<T>, var::Parameters_values>, int> = 0>
inline Maybe_error<void> load_literal_from_json_helper(
    const typename json_spec<Lexer>::Json& value, const std::string& path,
    typename json_spec<Lexer>::TagPolicy /*policy*/, const Identifier<Lexer>& id,
    Environment<Lexer, Compiler>& env) {
    using Json = typename json_spec<Lexer>::Json;
    if (value.type != Json::Type::Object) {
        return error_message(path + ": expected object for Parameters_values");
    }
    const Json* sid = value.find("schema_id");
    const Json* vals = value.find("values");
    if (!sid || sid->type != Json::Type::String || !vals) {
        return error_message(path + ": missing 'schema_id' or 'values'");
    }
    std::string schema_id = sid->str;
    auto schema = env.get_parameter_schema(schema_id);
    if (!schema) return error_message(path + ": unknown schema id '" + schema_id + "'");
    Matrix<double> mv;
    auto st = json_spec<Lexer>::from_json(*vals, mv, path + ".values",
                                          typename json_spec<Lexer>::TagPolicy{});
    if (!st) return st;
    var::Parameters_values pv(*schema, mv);
    auto literal = std::make_unique<typed_literal<Lexer, Compiler, var::Parameters_values>>(pv);
    env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(literal->clone()));
    env.push_back(id, new Identifier_compiler<Lexer, Compiler, var::Parameters_values>(literal.release()));
    return {};
}

template <class Lexer, class Compiler>
class typed_program {
    std::vector<std::unique_ptr<base_typed_statement<Lexer, Compiler>>> m_statements;

   public:
    typed_program() = default;
    typed_program(const typed_program& other) : m_statements{clone(other.m_statements)} {}
    typed_program(typed_program&& other) noexcept : m_statements{std::move(other.m_statements)} {}
    typed_program& operator=(const typed_program& other) {
        if (this != &other) {
            m_statements = clone(other.m_statements);
        }
        return *this;
    }
    typed_program& operator=(typed_program&&) noexcept = default;
    ~typed_program() = default;

    auto& push_back(std::unique_ptr<base_typed_statement<Lexer, Compiler>> stmt) {
        m_statements.emplace_back(std::move(stmt));
        return *this;
    }

    auto& push_back(base_typed_statement<Lexer, Compiler>* t_expr) {
        return push_back(std::unique_ptr<base_typed_statement<Lexer, Compiler>>(t_expr));
    }

    Maybe_error<Environment<Lexer, Compiler>> run(Environment<Lexer, Compiler>& env) {
        for (auto& e : m_statements) {
            auto Maybe_stat = e->run_statement(env);
            if (!Maybe_stat) {
                return Maybe_stat.error();
            }
        }
        return env;
    }
};

template <class Lexer, class Compiler>
Maybe_error<typed_argument_list<Lexer, Compiler>*>
    base_typed_assigment<Lexer, Compiler>::compile_argument_list(
        typed_argument_list<Lexer, Compiler>* t_growing_list) const {
    t_growing_list->insert_assigment(*this);
    return t_growing_list;
}

template <class Lexer, class Compiler>
Maybe_error<typed_argument_list<Lexer, Compiler>*>
    base_typed_expression<Lexer, Compiler>::compile_argument_list(
        typed_argument_list<Lexer, Compiler>* t_growing_list) const {
    auto May_be = t_growing_list->push_back_expression(*this);
    if (!May_be) {
        return May_be.error();
    }
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
        : m_tuple{dsl::clone(t_args)...} {}

    typed_argument_list(
        std::unique_ptr<typed_expression<Lexer,Compiler,Ts>>&&... t_args)
        : m_tuple{std::move(t_args)...} {}
    typed_argument_list(
        const typed_argument_list &other)
        : m_tuple{std::apply( [](auto const& ...arg){ return
std::tuple(dsl::clone(arg)...);},other)}{}

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
}  // namespace macrodr::dsl

#endif  // GRAMMAR_TYPED_H
