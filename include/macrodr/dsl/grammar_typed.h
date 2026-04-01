#ifndef GRAMMAR_TYPED_H
#define GRAMMAR_TYPED_H
#include "grammar_Identifier.h"
//#include "grammar_untyped.h"
#include <concepts>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "deprecation.h"
#include "lexer_typed.h"
#include "macrodr/dsl/lexer_untyped.h"
#include "maybe_error.h"
#include "type_name.h"
// Schemas
#include "parameters.h"
#include "indexed.h"

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

template <class... Abstract>
std::tuple<std::unique_ptr<Abstract>...> clone_tuple_unique(
    const std::tuple<std::unique_ptr<Abstract>...>& x) {
    return std::apply(
        [](auto const&... e) -> std::tuple<std::unique_ptr<Abstract>...> {
            return std::tuple<std::unique_ptr<Abstract>...>(clone_strict(e.get())...);
        },
        x);
}

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

    [[nodiscard]] Maybe_error<base_typed_expression<Lexer, Compiler> const*> get_RunValue(
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
class base_Identifier_compiler;

template <class Lexer, class Compiler>
class base_typed_statement {
   public:
    virtual ~base_typed_statement() = default;
    [[nodiscard]] virtual std::unique_ptr<base_typed_statement> clone_unique() const = 0;
    virtual Maybe_error<bool> run_statement(Environment<Lexer, Compiler>& env) const = 0;

    [[nodiscard]] virtual std::string type_name() const = 0;

    [[nodiscard]] virtual std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const = 0;
};

template <class Lexer, class Compiler>
class base_typed_assigment : public base_typed_statement<Lexer, Compiler> {
   public:
    ~base_typed_assigment() override = default;
    [[nodiscard]] virtual Identifier<Lexer> const& id() const = 0;
    [[nodiscard]] virtual base_typed_expression<Lexer, Compiler>* expr() const = 0;

};

template <class Lexer, class Compiler>
class base_typed_expression : public base_typed_statement<Lexer, Compiler> {
   public:
    ~base_typed_expression() override = default;

    using base_typed_statement<Lexer, Compiler>::compile_identifier_unique;

    
    virtual std::unique_ptr<base_typed_assigment<Lexer, Compiler>> compile_assigment_unique(
        Identifier<Lexer> id) = 0;

    [[nodiscard]] virtual Maybe_unique<base_typed_expression<Lexer, Compiler>>
        run_expression_unique(Environment<Lexer, Compiler>& env) const = 0;

    Maybe_error<bool> run_statement(Environment<Lexer, Compiler>& env) const override {
        auto Maybe_exp = run_expression_unique(env);
        if (!Maybe_exp) {
            return Maybe_exp.error();
        }
        return true;
    }

    [[nodiscard]] std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const override = 0;

    virtual Maybe_error<SerializedExpression> serialize_json(
        const Environment<Lexer, Compiler>& env,
        typename json_spec<Lexer>::TagPolicy policy) const {
        return error_message(std::string{"serialization unavailable for type "} +
                             this->type_name());
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
    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_identifier>(*this);
    }

    [[nodiscard]] Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const override {
        auto May_x = env.get_RunValue(m_id);
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

    std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique() const override {
        return std::make_unique<typed_identifier_ref_const>(*this);
    }

    Maybe_error<std::reference_wrapper<const T>> run(
        const Environment<Lexer, Compiler>& env) const override {
        auto May_x = env.get_RunValue(m_id);
        if (!May_x) {
            return May_x.error();
        }
        auto literal =
            dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(May_x.value());
        if (!literal) {
            auto expected = type_name<T>();
            auto actual = type_name(*May_x.value());
                  return error_message(std::string("\n\tIdentifier '") + m_id() +
                                 "' is not bound to the expected reference type \n\t\texpected: \n\t\t\t" +
                                 expected + ", \n\t\tactual: \n\t\t\t" + actual+ "\n" );
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

    std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique() const override {
        return std::make_unique<typed_identifier_ref>(*this);
    }

    Maybe_error<std::reference_wrapper<T>> run(
        const Environment<Lexer, Compiler>& env) const override {
        auto May_x = env.get_RunValue(m_id);
        if (!May_x) {
            return May_x.error();
        }
        auto literal_const =
            dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(May_x.value());
        if (!literal_const) {
            auto expected = type_name<T>();
            auto actual = type_name(*May_x.value());
            return error_message(std::string("\nIdentifier '") + m_id() +
                                 "' is not bound to the expected reference type \n\texpected: \n\t\t" +
                                 expected + ", \n\tactual: \n\t\t" + actual+ "\n" );
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
    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override = 0;

    [[nodiscard]] virtual bool has_IndexSpace() const { return false;}


    [[nodiscard]] virtual var::IndexSpace get_IndexSpace(Environment<Lexer, Compiler> const&) const {return {};};


    [[nodiscard]] virtual Maybe_error<T> run(Environment<Lexer, Compiler> const&, var::Coordinate const&) const = 0;


    [[nodiscard]] virtual Maybe_error<T> run(Environment<Lexer, Compiler> const&) const = 0;



    std::unique_ptr<base_typed_assigment<Lexer, Compiler>> compile_assigment_unique(
        Identifier<Lexer> id) override;

    
    [[nodiscard]] std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const override {
        return std::make_unique<Identifier_compiler<Lexer, Compiler, T>>(clone_strict(this));
    }

    // base_typed_expression interface

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> run_expression_unique(
        dsl::Environment<Lexer, Compiler>& env) const override {
        if constexpr (std::is_void_v<T>) {
            if(!has_IndexSpace(env)) {
                run(env);
                return std::make_unique<typed_literal<Lexer, Compiler, T>>();   
            }
            auto indexSpace = get_IndexSpace(env);
            auto coord = indexSpace.begin();
               while(true){
                auto maybe_run = run(env, coord);
                if (!maybe_run) {
                    return error_message(std::string("\nIn expression at coordinate ") + coord.str() + ": " +
                                         maybe_run.error());
                }
                if(coord == indexSpace.end()){
                    break;
                }
                coord.next();
            }
         return std::make_unique<typed_literal<Lexer, Compiler, T>>();
        }
        if(!has_IndexSpace(env)) {
            auto Maybe_x = run(env);
            if (!Maybe_x) {
                return Maybe_x.error();
            }
            return std::make_unique<typed_literal<Lexer, Compiler, T>>(std::move(Maybe_x.value()));
            }
        auto indexSpace = get_IndexSpace(env);
            std::vector<T> values;
            auto coord = indexSpace.begin();
            while(true){
                auto maybe_run = run(env, coord);
                if (!maybe_run) {
                    return error_message(std::string("\nIn expression at coordinate ") + coord.str() + ": " +
                                         maybe_run.error());
                }
                values.push_back(std::move(maybe_run.value()));
                if(coord == indexSpace.end()){
                    break;
                }
                coord.next();
            }
            return std::make_unique<typed_literal<Lexer, Compiler, var::Indexed<T>>>({std::move(indexSpace), std::move(values)});                                         
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
class typed_expression<Lexer, Compiler, var::Indexed<T>> : public base_typed_expression<Lexer, Compiler> {
   public:
    [[nodiscard]] std::string type_name() const override { return macrodr::dsl::type_name<var::Indexed<T>>(); }

    ~typed_expression() override = default;
    // Complete special members to satisfy rule-of-five guidance for abstract base
    typed_expression() = default;
    typed_expression(const typed_expression&) = default;
    typed_expression& operator=(const typed_expression&) = default;
    typed_expression(typed_expression&&) noexcept = default;
    typed_expression& operator=(typed_expression&&) noexcept = default;
    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override = 0;

    [[nodiscard]] virtual bool has_IndexSpace() const { return false;}


    [[nodiscard]] virtual var::IndexSpace get_IndexSpace(Environment<Lexer, Compiler> const&) const {
        return {};
    };

    [[nodiscard]] virtual var::IndexSpace get_myIndexSpace(Environment<Lexer, Compiler> const&) const = 0;

    [[nodiscard]] virtual Maybe_error<var::Indexed<T>> run(Environment<Lexer, Compiler> const& env, var::Coordinate const&/**/) const {
        return run(env);    
    };

    [[nodiscard]] virtual Maybe_error<T> run_by_coordinate(Environment<Lexer, Compiler> const&, var::Coordinate const&) const = 0;


    [[nodiscard]] virtual Maybe_error<var::Indexed<T>> run(Environment<Lexer, Compiler> const&) const = 0;



    std::unique_ptr<base_typed_assigment<Lexer, Compiler>> compile_assigment_unique(
        Identifier<Lexer> id) override;

    
    [[nodiscard]] std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const override {
        return std::make_unique<Identifier_compiler<Lexer, Compiler, var::Indexed<T>>>(clone_strict(this));
    }

    // base_typed_expression interface

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> run_expression_unique(
        dsl::Environment<Lexer, Compiler>& env) const override {
        if constexpr (std::is_void_v<T>) {
            if(!has_IndexSpace(env)) {
                run(env);
                return std::make_unique<typed_literal<Lexer, Compiler, T>>();   
            }
            auto indexSpace = get_IndexSpace(env);
            auto coord = indexSpace.begin();
               while(true){
                auto maybe_run = run(env, coord);
                if (!maybe_run) {
                    return error_message(std::string("\nIn expression at coordinate ") + coord.str() + ": " +
                                         maybe_run.error());
                }
                if(coord == indexSpace.end()){
                    break;
                }
                coord.next();
            }
         return std::make_unique<typed_literal<Lexer, Compiler, T>>();
        }
        if(!has_IndexSpace(env)) {
            auto Maybe_x = run(env);
            if (!Maybe_x) {
                return Maybe_x.error();
            }
            return std::make_unique<typed_literal<Lexer, Compiler, T>>(std::move(Maybe_x.value()));
            }
        auto indexSpace = get_IndexSpace(env);
            std::vector<T> values;
            auto coord = indexSpace.begin();
            while(true){
                auto maybe_run = run(env, coord);
                if (!maybe_run) {
                    return error_message(std::string("\nIn expression at coordinate ") + coord.str() + ": " +
                                         maybe_run.error());
                }
                values.push_back(std::move(maybe_run.value()));
                if(coord == indexSpace.end()){
                    break;
                }
                coord.next();
            }
            return std::make_unique<typed_literal<Lexer, Compiler, var::Indexed<T>>>({std::move(indexSpace), std::move(values)});                                         
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
class typed_Indexed_expression : public typed_expression<Lexer, Compiler, T> {
     std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> m_indexed_expr;
    public:
    typed_Indexed_expression(std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> expr) : m_indexed_expr(std::move(expr)) {}
    ~typed_Indexed_expression() override = default;

    ~typed_Indexed_expression() override = default;

      [[nodiscard]] virtual bool has_IndexSpace() const { return true;}

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_Indexed_expression>(*this);
    }

    [[nodiscard]] virtual var::IndexSpace get_IndexSpace(Environment<Lexer, Compiler> const& env) const {
        return m_indexed_expr->get_myIndexSpace(env);
    }

    [[nodiscard]] Maybe_error<T> run(Environment<Lexer, Compiler> const& env, var::Coordinate coord) const override {
        return m_indexed_expr->run_by_coordinate(env, coord);
    }
};


template <class Lexer, class Compiler, class T>
class typed_assigment : public base_typed_assigment<Lexer, Compiler> {
    Identifier<Lexer> m_id;
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

   public:
    [[nodiscard]] std::string type_name() const override { return macrodr::dsl::type_name<T>(); }

    ~typed_assigment() override = default;
    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_assigment>(*this);
    }

    typed_assigment(const typed_assigment& other)
        : m_id{other.m_id}, m_expr{clone_strict(other.m_expr.get())} {}

    typed_assigment(typed_assigment&&) noexcept = default;
    typed_assigment& operator=(typed_assigment&&) noexcept = default;
    typed_assigment& operator=(const typed_assigment& other) {
        if (this != &other) {
            m_id = other.m_id;
            m_expr = clone_strict(other.m_expr.get());
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
            return error_message(std::string("\nIn assignment to ") + id()() + ": " +
                                 maybe_expr.error()()+"\n");
        }
        env.insert(id(), std::move(maybe_expr.value()));
        return true;
    }

    [[nodiscard]] std::unique_ptr<base_Identifier_compiler<Lexer, Compiler>>
        compile_identifier_unique() const override {
        return std::make_unique<Identifier_compiler<Lexer, Compiler, T>>(clone_strict(this->expr()));
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
        : m_args{clone_tuple_unique(other.m_args)},
          m_f{other.m_f} {}

    typed_function_evaluation(typed_function_evaluation&&) noexcept = default;
    typed_function_evaluation& operator=(typed_function_evaluation&&) noexcept = default;
    typed_function_evaluation& operator=(const typed_function_evaluation& other) {
        if (this != &other) {
            m_args = clone_tuple_unique(other.m_args);
            m_f = other.m_f;
        }
        return *this;
    }

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_function_evaluation>(*this);
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

    [[nodiscard]] virtual bool has_IndexSpace(Environment<Lexer, Compiler> const&) const{
        return std::apply([](auto const&... args){ return ((args->has_IndexSpace()||...)); }, m_args);
    }


    [[nodiscard]] virtual var::IndexSpace get_IndexSpace(Environment<Lexer, Compiler> const& cm) const{
        return std::apply([&cm](auto const&... args){
           return get_IndexSpaces(args->get_IndexSpace(cm)...); }, m_args);
    }
    

    [[nodiscard]] virtual Maybe_error<T> run(Environment<Lexer, Compiler> const& env, var::Coordinate const& coord) const {
        if constexpr (sizeof...(Args) == 0) {
            if constexpr (std::is_void_v<T>) {
                this->m_f();
                return Maybe_error<void>();
            } else {
                return Maybe_error<T>(this->m_f());
            }
        } else {
            auto Maybe_arg_tuple = std::apply(
                [this, &env, &coord](auto&... args) { return std::tuple(args->run(env, coord )...); }, m_args);
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
        : m_arg{clone_strict(other.m_arg.get())}, m_P{other.m_P} {}

    typed_predicate_evaluation(typed_predicate_evaluation&&) noexcept = default;
    typed_predicate_evaluation& operator=(typed_predicate_evaluation&&) noexcept = default;
    typed_predicate_evaluation& operator=(const typed_predicate_evaluation& other) {
        if (this != &other) {
            m_arg = clone_strict(other.m_arg.get());
            m_P = other.m_P;
        }
        return *this;
    }

    std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique() const override {
        return std::make_unique<typed_predicate_evaluation>(*this);
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

namespace detail {

template <class Container>
struct homogeneous_container_ops;

template <class T, class Alloc>
struct homogeneous_container_ops<std::vector<T, Alloc>> {
    static void prepare(std::vector<T, Alloc>& out, std::size_t n) { out.reserve(n); }

    template <class U>
    static void insert(std::vector<T, Alloc>& out, U&& value) {
        out.emplace_back(std::forward<U>(value));
    }
};

template <class T, class Compare, class Alloc>
struct homogeneous_container_ops<std::set<T, Compare, Alloc>> {
    static void prepare(std::set<T, Compare, Alloc>& /*out*/, std::size_t /*n*/) {}

    template <class U>
    static void insert(std::set<T, Compare, Alloc>& out, U&& value) {
        out.insert(std::forward<U>(value));
    }
};

template <class Lexer, class Compiler, class Container, class expressedType,class T>
class typed_homogeneous_container_construction
    : public typed_expression<Lexer, Compiler, expressedType> {
   public:
    template <class U>
    using storage_t = function_argument_storage_t<U>;

    using storage_args =
        std::vector<std::unique_ptr<typed_expression<Lexer, Compiler, storage_t<T>>>>;

    explicit typed_homogeneous_container_construction(storage_args&& args)
        : m_args{std::move(args)} {}

    explicit typed_homogeneous_container_construction(const storage_args& args)
        : m_args{clone_vector(args)} {}

    ~typed_homogeneous_container_construction() override = default;

    typed_homogeneous_container_construction(
        const typed_homogeneous_container_construction& other)
        : m_args{clone_vector(other.m_args)} {}

    typed_homogeneous_container_construction(
        typed_homogeneous_container_construction&&) noexcept = default;
    typed_homogeneous_container_construction& operator=(
        typed_homogeneous_container_construction&&) noexcept = default;

    typed_homogeneous_container_construction& operator=(
        const typed_homogeneous_container_construction& other) {
        if (this != &other) {
            m_args = clone_vector(other.m_args);
        }
        return *this;
    }

    [[nodiscard]] Maybe_error<Container>
    run(const Environment<Lexer, Compiler>& env) const override {
        Container out;
        homogeneous_container_ops<Container>::prepare(out, m_args.size());
        std::string err;

        for (std::size_t i = 0; i < m_args.size(); ++i) {
            auto maybe_elem = m_args[i]->run(env);
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                homogeneous_container_ops<Container>::insert(
                    out, adapt_argument_like<T>(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
        }
        return out;
    }

   protected:
    storage_args m_args;
};

}  // namespace detail

template <class Lexer, class Compiler, class T>
class typed_vector_construction
    : public detail::typed_homogeneous_container_construction<Lexer, Compiler, std::vector<T>, std::vector<T>, T> {
    using base_type =
        detail::typed_homogeneous_container_construction<Lexer, Compiler, std::vector<T>,    std::vector<T>, T>;

   public:
    using storage_args = typename base_type::storage_args;
    using base_type::base_type;

    ~typed_vector_construction() override = default;

    typed_vector_construction(const typed_vector_construction&) = default;
    typed_vector_construction(typed_vector_construction&&) noexcept = default;
    typed_vector_construction& operator=(const typed_vector_construction&) = default;
    typed_vector_construction& operator=(typed_vector_construction&&) noexcept = default;

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_vector_construction>(*this);
    }
};

template <class Lexer, class Compiler, class T>
class typed_set_construction
    : public detail::typed_homogeneous_container_construction<Lexer, Compiler, std::set<T>,std::set<T>, T> {
    using base_type =
        detail::typed_homogeneous_container_construction<Lexer, Compiler, std::set<T>,    std::set<T>,T>;

   public:
    using storage_args = typename base_type::storage_args;
    using base_type::base_type;

    ~typed_set_construction() override = default;

    typed_set_construction(const typed_set_construction&) = default;
    typed_set_construction(typed_set_construction&&) noexcept = default;
    typed_set_construction& operator=(const typed_set_construction&) = default;
    typed_set_construction& operator=(typed_set_construction&&) noexcept = default;

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_set_construction>(*this);
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

   public:
    using tuple_type   = std::tuple<Ts...>;
    using storage_args = std::tuple<std::unique_ptr<typed_expression<Lexer,Compiler,storage_t<Ts>>>...>;

    explicit typed_tuple_construction(storage_args&& t)
        : m_args{std::move(t)} {}

    ~typed_tuple_construction() override = default;

    typed_tuple_construction(const typed_tuple_construction& other)
        : m_args{clone_tuple_unique(other.m_args)} {}

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_tuple_construction>(*this);
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
            return std::tuple<Ts...>( detail::adapt_argument_like<Ts>(me.value())... );
        };

        return std::apply(builder, maybe_storage);
    }
};




template <class Lexer, class Compiler, class T>
std::unique_ptr<base_typed_assigment<Lexer, Compiler>>
typed_expression<Lexer, Compiler, T>::compile_assigment_unique(Identifier<Lexer> id) {
    return std::make_unique<typed_assigment<Lexer, Compiler, T>>(id, clone_strict(this).release());
}

template <class Lexer, class Compiler, class T, class S>
    requires(std::convertible_to<S, T>)
class typed_conversion : public typed_expression<Lexer, Compiler, T> {
    std::unique_ptr<typed_expression<Lexer, Compiler, S>> m_expr;

   public:
    typed_conversion(std::unique_ptr<typed_expression<Lexer, Compiler, S>>&& x)
        : m_expr(std::move(x)) {}
    typed_conversion(std::unique_ptr<typed_expression<Lexer, Compiler, S>> const& x)
        : m_expr(clone_strict(x.get())) {}
    typed_conversion(typed_expression<Lexer, Compiler, S>* x) : m_expr(x) {}
    typed_conversion(const typed_conversion& other) : m_expr(clone_strict(other.m_expr.get())) {}
    typed_conversion& operator=(const typed_conversion& other) {
        if (this != &other) {
            m_expr = clone_strict(other.m_expr.get());
        }
        return *this;
    }
    typed_conversion(typed_conversion&&) noexcept = default;
    typed_conversion& operator=(typed_conversion&&) noexcept = default;

    virtual ~typed_conversion() = default;

    Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const override {
        auto r = m_expr->run(env);
        if (!r) {
            return r.error();
        }
        return static_cast<T>(std::move(r.value()));
    }

    // typed_expression interface

    std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique() const override {
        return std::make_unique<typed_conversion>(*this);
    }
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

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_literal>(*this);
    }
    
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

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_literal>(*this);
    }
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

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_literal>(*this);
    }
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

    [[nodiscard]] std::unique_ptr<base_typed_statement<Lexer, Compiler>> clone_unique()
        const override {
        return std::make_unique<typed_literal>(*this);
    }
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
    env.insert(id, clone_strict(literal.get()));
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
    env.insert(id, clone_strict(literal.get()));
    env.push_back(id, new Identifier_compiler<Lexer, Compiler, var::Parameters_values>(literal.release()));
    return {};
}

template <class Lexer, class Compiler>
class typed_program {
    std::vector<std::unique_ptr<base_typed_statement<Lexer, Compiler>>> m_statements;

   public:
    typed_program() = default;
    typed_program(const typed_program& other) : m_statements{clone_vector(other.m_statements)} {}
    typed_program(typed_program&& other) noexcept : m_statements{std::move(other.m_statements)} {}
    typed_program& operator=(const typed_program& other) {
        if (this != &other) {
            m_statements = clone_vector(other.m_statements);
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

}  // namespace macrodr::dsl

#endif  // GRAMMAR_TYPED_H
