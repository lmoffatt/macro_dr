#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
//#include "grammar_typed.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "grammar_Identifier.h"
#include "literal_decode.h"
#include "type_name.h"
#include <indexed.h>
#include <macrodr/interface/IModel.h>
#include <macrodr/io/json/convert.h>

#include "maybe_error.h"

#include "json_spec.h"
#include "dsl_argument_traits.h"
#include "dsl_forward.h"

namespace macrodr::dsl {
// JSON specification hook: by default, map to macrodr::io::json::conv



 template<class Lexer, class Compiler,class T>
 struct tuple_compiler_impl{
    using type=T;
 };

template<class Lexer, class Compiler,class... Ts>
 struct tuple_compiler_impl<Lexer, Compiler,std::tuple<Ts...>>{
    using type=tuple_compiler<Lexer,Compiler,Ts...>;
 };

 template<class Lexer,class Compiler, class T>
 using tuple_compiler_t=typename tuple_compiler_impl<Lexer,Compiler, T>::type;



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
        return new Identifier_compiler(clone_strict(m_expr.get()));
    };

    Identifier_compiler(std::unique_ptr<typed_expression<Lexer, Compiler, T>> t_expr)
        : m_expr{std::move(t_expr)} {}
    Identifier_compiler(typed_expression<Lexer, Compiler, T>* t_expr) : m_expr{t_expr} {}

    [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_Identifier(
        const Identifier<Lexer>& id) const override {
        return new typed_identifier<Lexer, Compiler, T>(id);
    }
};

template <class Lexer, class Compiler>
struct compiled_function_candidate {
    std::unique_ptr<base_typed_expression<Lexer, Compiler>> expr;
    std::size_t overload_rank = 0;
};

template <class Lexer, class Compiler>
class base_typed_argument {
   public:
    virtual ~base_typed_argument() = default;
};

template <class Lexer, class Compiler, class T>
class typed_argument : public base_typed_argument<Lexer, Compiler> {
   public:
    virtual ~typed_argument() override = default;
    [[nodiscard]] virtual std::unique_ptr<typed_argument> clone_unique() const = 0;
    [[nodiscard]] virtual bool has_index_space(
        Environment<Lexer, Compiler> const& env) const = 0;
    [[nodiscard]] virtual Maybe_error<var::IndexSpace> index_space(
        Environment<Lexer, Compiler> const& env) const = 0;
    [[nodiscard]] virtual Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const = 0;
    [[nodiscard]] virtual Maybe_error<T> run_at(Environment<Lexer, Compiler> const& env,
                                                var::Coordinate const& coord) const = 0;
};

template <class Lexer, class Compiler, class T>
class expression_argument : public typed_argument<Lexer, Compiler, T> {
    std::unique_ptr<typed_expression<Lexer, Compiler, T>> m_expr;

   public:
    explicit expression_argument(std::unique_ptr<typed_expression<Lexer, Compiler, T>> expr)
        : m_expr(std::move(expr)) {}
    explicit expression_argument(typed_expression<Lexer, Compiler, T>* expr) : m_expr(expr) {}

    expression_argument(const expression_argument& other)
        : m_expr(clone_strict(other.m_expr.get())) {}
    expression_argument(expression_argument&&) noexcept = default;
    expression_argument& operator=(expression_argument&&) noexcept = default;
    expression_argument& operator=(const expression_argument& other) {
        if (this != &other) {
            m_expr = clone_strict(other.m_expr.get());
        }
        return *this;
    }

    [[nodiscard]] typed_expression<Lexer, Compiler, T> const* expression() const {
        return m_expr.get();
    }

    [[nodiscard]] std::unique_ptr<typed_argument<Lexer, Compiler, T>> clone_unique()
        const override {
        return std::make_unique<expression_argument>(*this);
    }

    [[nodiscard]] bool has_index_space(Environment<Lexer, Compiler> const&) const override {
        return false;
    }

    [[nodiscard]] Maybe_error<var::IndexSpace> index_space(
        Environment<Lexer, Compiler> const&) const override {
        return error_message("scalar argument does not expose an index space");
    }

    [[nodiscard]] Maybe_error<T> run(Environment<Lexer, Compiler> const& env) const override {
        return m_expr->run(env);
    }

    [[nodiscard]] Maybe_error<T> run_at(Environment<Lexer, Compiler> const& env,
                                        var::Coordinate const&) const override {
        return m_expr->run(env);
    }
};

template <class Lexer, class Compiler, class T>
class pointwise_argument : public typed_argument<Lexer, Compiler, T> {
    std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> m_expr;

   public:
    explicit pointwise_argument(
        std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> expr)
        : m_expr(std::move(expr)) {}
    explicit pointwise_argument(typed_expression<Lexer, Compiler, var::Indexed<T>>* expr)
        : m_expr(expr) {}

    pointwise_argument(const pointwise_argument& other)
        : m_expr(clone_strict(other.m_expr.get())) {}
    pointwise_argument(pointwise_argument&&) noexcept = default;
    pointwise_argument& operator=(pointwise_argument&&) noexcept = default;
    pointwise_argument& operator=(const pointwise_argument& other) {
        if (this != &other) {
            m_expr = clone_strict(other.m_expr.get());
        }
        return *this;
    }

    [[nodiscard]] typed_expression<Lexer, Compiler, var::Indexed<T>> const* expression() const {
        return m_expr.get();
    }

    [[nodiscard]] std::unique_ptr<typed_argument<Lexer, Compiler, T>> clone_unique()
        const override {
        return std::make_unique<pointwise_argument>(*this);
    }

    [[nodiscard]] bool has_index_space(
        Environment<Lexer, Compiler> const&) const override {
        return true;
    }

    [[nodiscard]] Maybe_error<var::IndexSpace> index_space(
        Environment<Lexer, Compiler> const& env) const override {
        return m_expr->index_space(env);
    }

    [[nodiscard]] Maybe_error<T> run(Environment<Lexer, Compiler> const&) const override {
        return error_message("pointwise argument requires a coordinate");
    }

    [[nodiscard]] Maybe_error<T> run_at(Environment<Lexer, Compiler> const& env,
                                        var::Coordinate const& coord) const override {
        return m_expr->run_at(env, coord);
    }
};

template <class Lexer, class Compiler, class T>
class borrowed_expression_argument
    : public typed_argument<Lexer, Compiler, std::reference_wrapper<const T>> {
    std::unique_ptr<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>> m_expr;

   public:
    explicit borrowed_expression_argument(
        std::unique_ptr<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>> expr)
        : m_expr(std::move(expr)) {}
    explicit borrowed_expression_argument(
        typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>* expr)
        : m_expr(expr) {}

    borrowed_expression_argument(const borrowed_expression_argument& other)
        : m_expr(clone_strict(other.m_expr.get())) {}
    borrowed_expression_argument(borrowed_expression_argument&&) noexcept = default;
    borrowed_expression_argument& operator=(borrowed_expression_argument&&) noexcept = default;
    borrowed_expression_argument& operator=(const borrowed_expression_argument& other) {
        if (this != &other) {
            m_expr = clone_strict(other.m_expr.get());
        }
        return *this;
    }

    [[nodiscard]] std::unique_ptr<
        typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    clone_unique() const override {
        return std::make_unique<borrowed_expression_argument>(*this);
    }

    [[nodiscard]] bool has_index_space(Environment<Lexer, Compiler> const&) const override {
        return false;
    }

    [[nodiscard]] Maybe_error<var::IndexSpace> index_space(
        Environment<Lexer, Compiler> const&) const override {
        return error_message("borrowed scalar argument does not expose an index space");
    }

    [[nodiscard]] Maybe_error<std::reference_wrapper<const T>> run(
        Environment<Lexer, Compiler> const& env) const override {
        return m_expr->run(env);
    }

    [[nodiscard]] Maybe_error<std::reference_wrapper<const T>> run_at(
        Environment<Lexer, Compiler> const& env, var::Coordinate const&) const override {
        return m_expr->run(env);
    }
};

template <class Lexer, class Compiler, class T>
class pointwise_borrowed_argument
    : public typed_argument<Lexer, Compiler, std::reference_wrapper<const T>> {
    std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> m_expr;

   public:
    explicit pointwise_borrowed_argument(
        std::unique_ptr<typed_expression<Lexer, Compiler, var::Indexed<T>>> expr)
        : m_expr(std::move(expr)) {}
    explicit pointwise_borrowed_argument(typed_expression<Lexer, Compiler, var::Indexed<T>>* expr)
        : m_expr(expr) {}

    pointwise_borrowed_argument(const pointwise_borrowed_argument& other)
        : m_expr(clone_strict(other.m_expr.get())) {}
    pointwise_borrowed_argument(pointwise_borrowed_argument&&) noexcept = default;
    pointwise_borrowed_argument& operator=(pointwise_borrowed_argument&&) noexcept = default;
    pointwise_borrowed_argument& operator=(const pointwise_borrowed_argument& other) {
        if (this != &other) {
            m_expr = clone_strict(other.m_expr.get());
        }
        return *this;
    }

    [[nodiscard]] std::unique_ptr<
        typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    clone_unique() const override {
        return std::make_unique<pointwise_borrowed_argument>(*this);
    }

    [[nodiscard]] bool has_index_space(Environment<Lexer, Compiler> const&) const override {
        return true;
    }

    [[nodiscard]] Maybe_error<var::IndexSpace> index_space(
        Environment<Lexer, Compiler> const& env) const override {
        return m_expr->index_space(env);
    }

    [[nodiscard]] Maybe_error<std::reference_wrapper<const T>> run(
        Environment<Lexer, Compiler> const&) const override {
        return error_message("pointwise borrowed argument requires a coordinate");
    }

    [[nodiscard]] Maybe_error<std::reference_wrapper<const T>> run_at(
        Environment<Lexer, Compiler> const& env, var::Coordinate const& coord) const override {
        return m_expr->run_at_ref(env, coord);
    }
};

template <class Lexer, class Compiler>
class base_function_compiler {
   public:
    virtual ~base_function_compiler() = default;

    [[nodiscard]] virtual base_function_compiler* clone() const = 0;

    [[nodiscard]] virtual Maybe_error<compiled_function_candidate<Lexer, Compiler>>
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

template <class Lexer, class Compiler, class T>
Maybe_unique<typed_expression<Lexer, Compiler, T>> get_typed_expresion(
    std::unique_ptr<base_typed_expression<Lexer, Compiler>>& expr) {
    auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get());
    if (ptr != nullptr) {
        return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
    }
    const auto expected = type_name<T>();
    const auto actual = expr ? type_name(*expr) : std::string{"<null>"};
    return error_message("\nunexpected type while getting typed expression: expected ", expected, ", actual ", actual);
}

template <class Lexer, class Compiler, class T, class S, class... Ss>
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
class vector_compiler;
template <class Lexer, class Compiler, class T>
class set_compiler;

template <class Lexer, class Compiler, class T>
class element_compiler_new {
    public:
    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> adapt_compiled_expression(
        std::unique_ptr<base_typed_statement<Lexer, Compiler>> expr) const {
        if (auto ptr = dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.get())) {
            return new expression_argument<Lexer, Compiler, T>(
                dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release()));
        }
        if constexpr (!is_of_this_template_type_v<T, std::reference_wrapper> &&
                      !is_of_this_template_type_v<T, var::Indexed>) {
            if (auto indexed =
                    dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.get())) {
                return new pointwise_argument<Lexer, Compiler, T>(
                    dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(
                        expr.release()));
            }
        }
        const auto expected = type_name<T>();
        const auto actual = expr ? expr->type_name() : std::string{"<null>"};
        return error_message(std::string("unexpected type in compilation of parameter: expected ") + expected + ", got " +
                             actual);
    }

   public:
    element_compiler_new()=default;    
    
    
    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_id = cm.get_Identifier(t_arg());
        if (!maybe_id) {
            return maybe_id.error();
        }
        return adapt_compiled_expression(std::move(maybe_id.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value
    {
        auto maybe_value = from_literal<T>(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        return new expression_argument<Lexer, Compiler, T>(
            new typed_literal<Lexer, Compiler, T>(std::move(maybe_value.value())));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, std::string>) {
            return new expression_argument<Lexer, Compiler, T>(
                new typed_literal<Lexer, Compiler, T>(t_arg()));
        } else {
            return error_message(std::string{"unexpected type while compiling this argument: expected "} + type_name<T>() +
                                 ", got string");
        }
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto maybe_expr = t_arg.compile_expression(cm);
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return adapt_compiled_expression(std::move(maybe_expr.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_vector_construction<Lexer, Compiler> const& t_arg) const
        requires(is_of_this_template_type_v<T, std::vector>)
    {
        auto maybe_expr =
            vector_compiler<Lexer, Compiler, typename T::value_type>{}.compile_vector_construction(
                cm, t_arg.args());
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return adapt_compiled_expression(std::move(maybe_expr.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_vector_construction<Lexer, Compiler> const& t_arg) const
        requires(is_of_this_template_type_v<T, std::set>)
    {
        auto maybe_expr =
            set_compiler<Lexer, Compiler, typename T::value_type>{}.compile_set_construction(
                cm, t_arg.args());
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return adapt_compiled_expression(std::move(maybe_expr.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_tuple_construction<Lexer, Compiler> const& t_arg) const
        requires(is_of_this_template_type_v<T, std::tuple>)
    {
        auto maybe_expr = tuple_compiler_t<Lexer, Compiler, T>{}.compile_tuple_construction(
            cm, t_arg);
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return adapt_compiled_expression(std::move(maybe_expr.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto s_lit = dynamic_cast<untyped_string_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *s_lit);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        if constexpr (is_of_this_template_type_v<T, std::vector>) {
            if (auto vec =
                    dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_argument(cm, *vec);
            }
        }
        if constexpr (is_of_this_template_type_v<T, std::set>) {
            if (auto vec =
                    dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_argument(cm, *vec);
            }
        }
        if constexpr (is_of_this_template_type_v<T, std::tuple>) {
            if (auto tuple =
                    dynamic_cast<untyped_tuple_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_argument(cm, *tuple);
            }
        }
        return error_message(t_arg.str(),
                             " is an untyped expression that is not an identifier nor a literal "
                             "nor a function evaluation");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message(
            "\nif it is neither an assignment nor an expression, what the fuck is it?");
    }
};


template <class Lexer, class Compiler, class T>
class parameter_compiler : public element_compiler_new<Lexer, Compiler, T>
{
    Identifier<Lexer> m_id;


   public:
    using base_type=element_compiler_new<Lexer, Compiler, T>;
    using  base_type::adapt_compiled_expression;
    using  base_type::compile_this_argument;
    

    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }


    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            bool const is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent) {
                return error_message(" the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            }
            auto out = compile_this_argument(cm, ptr->expression());
            if (!out) {
                return error_message(std::string("\n\t\t\targument ") + ptr->id()()() +
                                     ": \n\t\t\t\t" + out.error()());
            }
            return out;
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message(
            "\nif it is neither an assignment nor an expression, what the fuck is it?");
    }
};


template <class Lexer, class Compiler, class T>
class element_compiler_new<Lexer, Compiler, var::Indexed<T>> {

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    adapt_compiled_expression(std::unique_ptr<base_typed_expression<Lexer, Compiler>> expr) const {
        if (auto ptr =
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.get())) {
            return new expression_argument<Lexer, Compiler, var::Indexed<T>>(
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.release()));
        }
        const auto expected = type_name<var::Indexed<T>>();
        const auto actual = expr ? expr->type_name() : std::string{"<null>"};
        return error_message(std::string("unexpected type while adapting compiled expression: expected ") + expected + ", got " +
                             actual);
    }

   public:

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_id = cm.get_Identifier(t_arg());
        if (!maybe_id) {
            return maybe_id.error();
        }
        return adapt_compiled_expression(std::move(maybe_id.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto maybe_expr = t_arg.compile_expression(cm);
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        return adapt_compiled_expression(std::move(maybe_expr.value()));
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        return error_message("indexed arguments must be identifiers or indexed-valued calls");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
class parameter_compiler<Lexer, Compiler, var::Indexed<T>> : public element_compiler_new<Lexer, Compiler, var::Indexed<T>> {
    Identifier<Lexer> m_id;

    

   public:
    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

    

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            bool const is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent) {
                return error_message(" the argument id expected: " + id()() +
                                     " found: " + ptr->id()()());
            }
            return compile_this_argument(cm, ptr->expression());
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
class element_compiler_new<Lexer, Compiler, std::reference_wrapper<const T>> {
  
   public:
  
  
    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_expr = cm.get_Identifier(t_arg());
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        if (dynamic_cast<typed_expression<Lexer, Compiler, T>*>(maybe_expr.value().get()) !=
            nullptr) {
            return new borrowed_expression_argument<Lexer, Compiler, T>(
                new typed_identifier_ref_const<Lexer, Compiler, T>(t_arg()));
        }
        if constexpr (!is_of_this_template_type_v<T, var::Indexed>) {
            auto expr = std::move(maybe_expr.value());
            auto ptr =
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.get());
            if (ptr != nullptr) {
                return new pointwise_borrowed_argument<Lexer, Compiler, T>(
                    dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(
                        expr.release()));
            }
        }
        return error_message(std::string{"argument '"} + 
                             "' expects a borrowable scalar or indexed identifier");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument "  
                             "' expects a reference to an environment variable; literals are "
                             "not allowed"});
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} +
                             "' expects a borrowable environment value; computed expressions are "
                             "not allowed");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for const reference parameter");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
class parameter_compiler<Lexer, Compiler, std::reference_wrapper<const T>> :
public element_compiler_new<Lexer, Compiler, std::reference_wrapper<const T>> {
    Identifier<Lexer> m_id;

   public:
    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_expr = cm.get_Identifier(t_arg());
        if (!maybe_expr) {
            return maybe_expr.error();
        }
        if (dynamic_cast<typed_expression<Lexer, Compiler, T>*>(maybe_expr.value().get()) !=
            nullptr) {
            return new borrowed_expression_argument<Lexer, Compiler, T>(
                new typed_identifier_ref_const<Lexer, Compiler, T>(t_arg()));
        }
        if constexpr (!is_of_this_template_type_v<T, var::Indexed>) {
            auto expr = std::move(maybe_expr.value());
            auto ptr =
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.get());
            if (ptr != nullptr) {
                return new pointwise_borrowed_argument<Lexer, Compiler, T>(
                    dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(
                        expr.release()));
            }
        }
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a borrowable scalar or indexed identifier");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a reference to an environment variable; literals are "
                             "not allowed");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a borrowable environment value; computed expressions are "
                             "not allowed");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for const reference parameter");
    }

    [[nodiscard]]
    Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<const T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto assignment = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            bool congruent = assignment->id()()() == this->id()();
            if (!congruent) {
                return error_message(" the argument id expected: " + id()() +
                                     " found: " + assignment->id()()());
            }
            return compile_this_argument(cm, assignment->expression());
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
class element_compiler_new<Lexer, Compiler, std::reference_wrapper<T>> {
    
   public:
    
    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_value = cm.get_RunValue(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        if (dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(maybe_value.value()) !=
            nullptr) {
            return new expression_argument<Lexer, Compiler, std::reference_wrapper<T>>(
                new typed_identifier_ref<Lexer, Compiler, T>(t_arg()));
        }
        return error_message(std::string{"argument '"} +
                             "' expects a mutable environment variable");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} +
                             "' expects a reference to an environment variable; literals are "
                             "not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} +
                             "' expects a mutable environment variable; function calls are not "
                             "allowed");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for reference parameter");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};


template <class Lexer, class Compiler, class T>
class parameter_compiler<Lexer, Compiler, std::reference_wrapper<T>> {
    Identifier<Lexer> m_id;

   public:
    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_value = cm.get_RunValue(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        if (dynamic_cast<const typed_literal<Lexer, Compiler, T>*>(maybe_value.value()) !=
            nullptr) {
            return new expression_argument<Lexer, Compiler, std::reference_wrapper<T>>(
                new typed_identifier_ref<Lexer, Compiler, T>(t_arg()));
        }
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a mutable environment variable");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a reference to an environment variable; literals are "
                             "not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& /*cm*/,
                          untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(std::string{"argument '"} + id()() +
                             "' expects a mutable environment variable; function calls are not "
                             "allowed");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_argument(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for reference parameter");
    }

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, std::reference_wrapper<T>>>
    compile_this_argument(Environment<Lexer, Compiler> const& cm,
                          std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto assignment = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            bool congruent = assignment->id()()() == this->id()();
            if (!congruent) {
                return error_message(" the argument id expected: " + id()() +
                                     " found: " + assignment->id()()());
            }
            return compile_this_argument(cm, assignment->expression());
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_argument(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
using field_compiler = parameter_compiler<Lexer, Compiler, T>;


template <class Lexer, class Compiler, class T>
struct element_compiler {
   
   public:
   
    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
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
        return error_message(std::string("unexpected type found in element compilation of an identifier: expected ") + expected + ", got " +
                             actual);
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value
    {
        auto maybe_value = from_literal<T>(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        return new typed_literal<Lexer, Compiler, T>(std::move(maybe_value.value()));
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, std::string>) {
            return new typed_literal<Lexer, Compiler, T>(t_arg());
        } else {
            return error_message(std::string{"unexpected type in element compilation of a string literal: expected "} + type_name<T>() +
                                 ", got string");
        }
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& /*unused*/,
        untyped_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        untyped_function_evaluation<Lexer, Compiler> const& t_arg) const {
        auto Maybe_expr = t_arg.compile_expression(cm);

        if (!Maybe_expr) {
            return Maybe_expr.error();
        }
        auto expr = std::move(Maybe_expr.value());
        auto exp_ptr = dynamic_cast<typed_expression<Lexer, Compiler, T> const*>(expr.get());
        if (exp_ptr != nullptr) {
            return dynamic_cast<typed_expression<Lexer, Compiler, T>*>(expr.release());
        }
        return error_message(t_arg.str(), "is a function evalution but not of type",
                             type_name<T>());
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        untyped_vector_construction<Lexer, Compiler> const& t_arg) const 
        requires(is_of_this_template_type_v<T,std::vector>)
        {
        return vector_compiler<Lexer,Compiler,typename T::value_type>{}.compile_vector_construction(
            cm, t_arg.args());
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        untyped_vector_construction<Lexer, Compiler> const& t_arg) const
        requires(is_of_this_template_type_v<T,std::set>)
        {
        return set_compiler<Lexer,Compiler,typename T::value_type>{}.compile_set_construction(
            cm, t_arg.args());
    }

[[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        untyped_tuple_construction<Lexer, Compiler> const& t_arg) const 
        requires(is_of_this_template_type_v<T,std::tuple>)
        {
        return tuple_compiler_t<Lexer,Compiler,T>{}.compile_tuple_construction(cm, t_arg);
    }


    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        untyped_expression<Lexer, Compiler> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg);
        if (ptr != nullptr) {
            return compile_this_element(cm, *ptr);
        }
        auto str_ptr = dynamic_cast<untyped_string_literal<Lexer, Compiler> const*>(&t_arg);
        if (str_ptr != nullptr) {
            return compile_this_element(cm, *str_ptr);
        }
        auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
        if (li_ptr != nullptr) {
            return compile_this_element(cm, *li_ptr);
        }
        if (auto f_ptr =
                dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *f_ptr);
        }
        if constexpr (is_of_this_template_type_v<T, std::vector>) {
            if (auto fn =
                    dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_element(cm, *fn);
            }
        }
        if constexpr (is_of_this_template_type_v<T, std::set>) {
            if (auto fn =
                    dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_element(cm, *fn);
            }
        }
        if constexpr (is_of_this_template_type_v<T, std::tuple>) {
            if (auto fn =
                    dynamic_cast<untyped_tuple_construction<Lexer, Compiler> const*>(&t_arg)) {
                return compile_this_element(cm, *fn);
            }
        }

        return error_message(t_arg.str(),
                             " is an untyped expression that is not an identifier nor a literal  "
                             "nor a function evaluation");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_element(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            auto out = compile_this_element(cm, ptr->expression());
            if (!out) {
                return error_message(std::string("\n element  ") + ptr->id()()() + ": \n" +
                                     out.error()());
            }
            return out;
        }
        auto ptr2 = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get());

        if (ptr2 != nullptr) {
            return compile_this_element(cm, *ptr2);
        }
        return error_message(
            "if it is  neither an assigment nor an expression, what the fuck is it?");
    }
};

template <class Lexer, class Compiler, class T>
class element_compiler<Lexer, Compiler, std::reference_wrapper<const T>> {
   
   public:
   
    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>>
        compile_this_element(Environment<Lexer, Compiler> const& cm,
                              untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_stored = cm.get_Identifier(t_arg());
        if (maybe_stored) {
            auto stored = maybe_stored.value();
            auto literal = dynamic_cast<const typed_expression<Lexer, Compiler, T>*>(stored);
            if (!literal) {
                const auto actual = stored ? stored->type_name() : std::string{"<null>"};
                return error_message(std::string{"identifier '"} + t_arg()() +
                                     "' cannot bind to expected'" + type_name<const T&>() +
                                     "'; current value has type " + actual);
            }
        }
        return new typed_identifier_ref_const<Lexer, Compiler, T>(t_arg());
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>>
        compile_this_element(Environment<Lexer, Compiler> const& /*cm*/,
                              untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(
            std::string{"argument '"}  +
            "' expects a reference to an environment variable; literals are not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>>
        compile_this_element(Environment<Lexer, Compiler> const& /*cm*/,
                              untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(
            std::string{"argument '"} + 
            "' expects a reference to an environment variable; function calls are not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>>
        compile_this_element(Environment<Lexer, Compiler> const& cm,
                              untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *fn);
        }
        if (auto fn = dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *fn);
        }
        if (auto fn = dynamic_cast<untyped_tuple_construction<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for reference parameter");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<const T>>>
        compile_this_element(
            Environment<Lexer, Compiler> const& cm,
            std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto assignment =
                dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_element(cm, assignment->expression());
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_element(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
    }
};

template <class Lexer, class Compiler, class T>
class element_compiler<Lexer, Compiler, std::reference_wrapper<T>> {
   
   public:
    
    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<T>>>
        compile_this_element(Environment<Lexer, Compiler> const& cm,
                              untyped_identifier<Lexer, Compiler> const& t_arg) const {
        auto maybe_stored = cm.get_Identifier(t_arg());
        if (maybe_stored) {
            auto stored = maybe_stored.value();
            auto literal = dynamic_cast<const typed_expression<Lexer, Compiler, T>*>(stored);
            if (!literal) {
                const auto actual = stored ? stored->type_name() : std::string{"<null>"};
                return error_message(std::string{"identifier '"} + t_arg()() +
                                     "' cannot bind to expected '" + type_name<T&>() +
                                     "'; current value has type " + actual);
            }
        }
        return new typed_identifier_ref<Lexer, Compiler, T>(t_arg());
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<T>>>
        compile_this_element(Environment<Lexer, Compiler> const& /*cm*/,
                              untyped_literal<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(
            std::string{"argument '"} + 
            "' expects a reference to an environment variable; literals are not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<T>>>
        compile_this_element(Environment<Lexer, Compiler> const& /*cm*/,
                              untyped_function_evaluation<Lexer, Compiler> const& /*t_arg*/) const {
        return error_message(
            std::string{"argument '"} +
            "' expects a reference to an environment variable; function calls are not allowed");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<T>>>
        compile_this_element(Environment<Lexer, Compiler> const& cm,
                              untyped_expression<Lexer, Compiler> const& t_arg) const {
        if (auto ptr = dynamic_cast<untyped_identifier<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *ptr);
        }
        if (auto lit = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *lit);
        }
        if (auto fn = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg)) {
            return compile_this_element(cm, *fn);
        }
        return error_message(
            "invalid argument: expected identifier expression for reference parameter");
    }

    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::reference_wrapper<T>>>
        compile_this_element(
            Environment<Lexer, Compiler> const& cm,
            std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        if (auto assignment =
                dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_element(cm, assignment->expression());
        }
        if (auto expr = dynamic_cast<untyped_expression<Lexer, Compiler> const*>(t_arg.get())) {
            return compile_this_element(cm, *expr);
        }
        return error_message("invalid argument: expected assignment or expression");
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
        requires literal_decodable<T>::value
    {
        auto maybe_value = from_literal<T>(t_arg());
        if (!maybe_value) {
            return maybe_value.error();
        }
        return new typed_literal<Lexer, Compiler, T>(std::move(maybe_value.value()));
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& t_arg) const
        requires literal_decodable<T>::value
    {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, std::string>) {
            return new typed_literal<Lexer, Compiler, T>(t_arg());
        } else {
            return error_message(std::string{"unexpected type while compiling string literal : expected "} + type_name<T>() +
                                 ", got string");
        }
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const& /*unused*/,
        untyped_string_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
        return error_message(std::string{"literal arguments are not supported for type "} +
                             type_name<T>() + "; use a variable or JSON helper function");
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Compiler const& /*unused*/, untyped_literal<Lexer, Compiler> const& /*t_arg*/) const
        requires(!literal_decodable<T>::value)
    {
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
        auto s_ptr = dynamic_cast<untyped_string_literal<Lexer, Compiler> const*>(&t_arg);
        if (s_ptr != nullptr) {
            return compile_this_argument(cm, *s_ptr);
        }
        auto li_ptr = dynamic_cast<untyped_literal<Lexer, Compiler> const*>(&t_arg);
        if (li_ptr != nullptr) {
            return compile_this_argument(cm, *li_ptr);
        }
        auto f_ptr = dynamic_cast<untyped_function_evaluation<Lexer, Compiler> const*>(&t_arg);
        if (f_ptr != nullptr) {
            return compile_this_argument(cm, *f_ptr);
        }
        auto v_ptr = dynamic_cast<untyped_vector_construction<Lexer, Compiler> const*>(&t_arg);
        if (v_ptr != nullptr) {
            return compile_this_argument(cm, *v_ptr);
        }
        auto t_ptr = dynamic_cast<untyped_tuple_construction<Lexer, Compiler> const*>(&t_arg);
        if (t_ptr != nullptr) {
            return compile_this_argument(cm, *t_ptr);
        }
        return error_message("unexpected node: expected identifier, literal, or function call");
    }

    Maybe_unique<typed_expression<Lexer, Compiler, T>> compile_this_argument(
        Environment<Lexer, Compiler> const& cm,
        std::unique_ptr<untyped_statement<Lexer, Compiler>> const& t_arg) const {
        auto ptr = dynamic_cast<untyped_assignment<Lexer, Compiler> const*>(t_arg.get());
        if (ptr != nullptr) {
            bool is_argument_id_congruent = ptr->id()()() == this->id()();
            if (!is_argument_id_congruent) {
                return error_message(" the argument id expected: " + id()() +
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
    template <class T>
    using storage_t = detail::function_argument_storage_t<T>;

    using argument_compilers_t =
        std::tuple<parameter_compiler<Lexer, Compiler, storage_t<Args>>...>;

    argument_compilers_t m_args;
    F m_f;

    [[nodiscard]] auto make_invoker() const {
        return [fn = m_f](storage_t<Args>... values) -> std::invoke_result_t<F, Args...> {
            return std::invoke(fn, detail::adapt_argument_like<Args>(values)...);
        };
    }

    template <class Tuple>
    [[nodiscard]] bool uses_pointwise_arguments(Environment<Lexer, Compiler> const& env,
                                                Tuple const& tuple) const {
        return std::apply(
            [&env](auto const&... args) { return ((args->has_index_space(env)) || ...); }, tuple);
    }

    template <std::size_t... Is>
    [[nodiscard]] Maybe_error<compiled_function_candidate<Lexer, Compiler>>
        compile_function_evaluation_impl(std::index_sequence<Is...> /*unused*/,
                                         Environment<Lexer, Compiler> const& cm,
                                         const untyped_argument_list<Lexer, Compiler>& args) const {
        auto invoker = make_invoker();
        using WrappedF = decltype(invoker);
        using Return = underlying_value_type_t<std::invoke_result_t<F, Args...>>;
        if constexpr (sizeof...(Is) == 0) {
            compiled_function_candidate<Lexer, Compiler> out;
            if constexpr (is_of_this_template_type_v<Return, var::Indexed>) {
                out.expr = std::make_unique<
                    typed_exact_indexed_function_evaluation<Lexer, Compiler, WrappedF>>(
                    std::move(invoker), std::tuple<>{});
            } else {
                out.expr =
                    std::make_unique<typed_scalar_function_evaluation<Lexer, Compiler, WrappedF>>(
                        std::move(invoker), std::tuple<>{});
            }
            out.overload_rank = 0;
            return out;
        } else {
            if (args.arg().size() != sizeof...(Is)) {
                return error_message("\n\t\tcount mismatch: \n\t\t\texpected ", sizeof...(Is), "\n\t\t\tgot ",
                                     args.arg().size(),"\n\t\t\texpected arguments ",(get<Is>(m_args).id()()+" ")..., "\n\t\t\tfound:  ", args.str());
            }
            auto Maybe_tuple = promote_Maybe_error(
                std::tuple(std::get<Is>(m_args).compile_this_argument(cm, args.arg()[Is])...));
            if (!Maybe_tuple) {
                return error_message("In", str(),": ",Maybe_tuple.error()());
            }
            auto compiled_args = std::move(Maybe_tuple.value());
            compiled_function_candidate<Lexer, Compiler> out;
            out.overload_rank = uses_pointwise_arguments(cm, compiled_args) ? 1U : 0U;
            if (out.overload_rank == 0) {
                if constexpr (is_of_this_template_type_v<Return, var::Indexed>) {
                    out.expr = std::make_unique<
                        typed_exact_indexed_function_evaluation<Lexer, Compiler, WrappedF,
                                                                storage_t<Args>...>>(
                        std::move(invoker), std::move(compiled_args));
                } else {
                    out.expr = std::make_unique<
                        typed_scalar_function_evaluation<Lexer, Compiler, WrappedF,
                                                         storage_t<Args>...>>(
                        std::move(invoker), std::move(compiled_args));
                }
            } else {
                if constexpr (std::is_void_v<Return>) {
                    return error_message(
                        "lifted indexed calls are not supported for void-returning functions");
                } else if constexpr (is_of_this_template_type_v<Return, var::Indexed>) {
                    return error_message(
                        "generic function_compiler does not support lifted native indexed "
                        "return types");
                } else {
                    out.expr = std::make_unique<
                        typed_lifted_function_evaluation<Lexer, Compiler, WrappedF,
                                                         storage_t<Args>...>>(
                        std::move(invoker), std::move(compiled_args));
                }
            }
            return out;
        }
    }

   public:
    function_compiler(F t_f, parameter_compiler<Lexer, Compiler, storage_t<Args>>&&... t_args)
        : m_args{std::move(t_args)...}, m_f{t_f} {}

    [[nodiscard]] Maybe_error<compiled_function_candidate<Lexer, Compiler>>
    compile_function_evaluation(
        Environment<Lexer, Compiler> const& cm,
        const untyped_argument_list<Lexer, Compiler>& args) const override {
        return compile_function_evaluation_impl(std::index_sequence_for<Args...>(), cm, args.arg());
    }

    // base_function_compiler interface

    [[nodiscard]] base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }

    auto str()  const -> std::string  {
        return std::apply(
            [this](auto&&... compiled_args) {
                return (std::string{"function with arguments: "} +
                        ... + (compiled_args.id()() + " "));
            },
            m_args);
    }  

    void register_types(Compiler& registry) const override {
        using Return = underlying_value_type_t<std::invoke_result_t<F, Args...>>;
        (registry.template ensure_type_registered<storage_t<Args>>(), ...);
        if constexpr (!std::is_void_v<Return>) {
            registry.template ensure_type_registered<Return>();
        }
    }
};


template <class Lexer, class Compiler, class... Ts>
class tuple_compiler {

    // Holder types used in evaluation (same rule as for function arguments)
    template <class U>
    using storage_t = detail::function_argument_storage_t<U>;

    using arg_compilers_t =
        std::tuple<element_compiler_new<Lexer, Compiler, storage_t<Ts>>...>;

    arg_compilers_t m_args;

template <class Tuple>
    [[nodiscard]] bool uses_pointwise_arguments(Environment<Lexer, Compiler> const& env,
                                                Tuple const& tuple) const {
        return std::apply(
            [&env](auto const&... args) { return ((args->has_index_space(env)) || ...); }, tuple);
    }



    template <std::size_t... Is>
    [[nodiscard]] Maybe_unique<
        base_typed_statement<Lexer, Compiler>>
    compile_impl(std::index_sequence<Is...>,
                 Environment<Lexer, Compiler> const& env,
                 const untyped_argument_list<Lexer, Compiler>& args) const
    {
        if (args.arg().size() != sizeof...(Is)) {
            return error_message("tuple arity mismatch: expected ",
                                 sizeof...(Is), ", got ", args.arg().size());
        }

        // Compile each element → unique_ptr<typed_expression<storage_t<T>>>
        auto maybe_args =
            promote_Maybe_error(std::tuple(
                std::get<Is>(m_args).compile_this_argument(env, args.arg()[Is])...
            ));

        if (!maybe_args) return maybe_args.error();


        

        bool is_pointwise = uses_pointwise_arguments(env, maybe_args.value());
        
        if (is_pointwise) {
            return std::make_unique<
                typed_lifted_tuple_construction<Lexer, Compiler, Ts...>>(
                std::move(maybe_args.value()));
        }
        
        return std::make_unique<
            typed_tuple_construction<Lexer, Compiler, Ts...>>(
                std::move(maybe_args.value()));
    }

public:

    [[nodiscard]] Maybe_unique<
        base_typed_statement<Lexer, Compiler>>
    compile_tuple_construction(Environment<Lexer, Compiler> const& env,
                               const untyped_tuple_construction<Lexer, Compiler>& t) const
    {
        return compile_impl(std::index_sequence_for<Ts...>{}, env, t.args());
    }
};



template <class Lexer, class Compiler, class T>
class vector_compiler {
    // Holder type for a single element
    template <class K>
    using storage_t = detail::function_argument_storage_t<K>;

    using element_compiler_t = element_compiler_new<Lexer, Compiler, storage_t<T>>;

    element_compiler_t m_elem;

   public:
    [[nodiscard]] bool uses_pointwise_arguments(Environment<Lexer, Compiler> const& env,
                                                std::vector<std::unique_ptr<typed_argument<Lexer, Compiler, storage_t<T>>>> const& compiled) const {
        for (const auto& arg : compiled) {
            if (arg->has_index_space(env)) {
                return true;
            }       
        }
         return false;
     }


   [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>>
    compile_vector_construction(Environment<Lexer, Compiler> const& cm,
                                const untyped_argument_list<Lexer, Compiler>& args) const {
        using expr_t = typed_argument<Lexer, Compiler, storage_t<T>>;

        std::vector<std::unique_ptr<expr_t>> compiled;
        compiled.reserve(args.arg().size());

        std::string err;

        auto n = args.arg().size();
        for (std::size_t i = 0; i < n; ++i) {
            auto maybe_elem = m_elem.compile_this_argument(cm, args.arg()[i]);
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                compiled.emplace_back(std::move(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
        }

        if (uses_pointwise_arguments(cm, compiled)) {
            return std::make_unique<typed_lifted_vector_construction<Lexer, Compiler, T>>(
                std::move(compiled));
        }
        return std::make_unique<typed_vector_construction<Lexer, Compiler, T>>(
            std::move(compiled));
    }
};

template <class Lexer, class Compiler, class T>
class set_compiler {
    template <class K>
    using storage_t = detail::function_argument_storage_t<K>;

    using element_compiler_t = element_compiler_new<Lexer, Compiler, storage_t<T>>;

    element_compiler_t m_elem;

   public:
    [[nodiscard]] bool uses_pointwise_arguments(Environment<Lexer, Compiler> const& env,
                                                std::vector<std::unique_ptr<typed_argument<Lexer, Compiler, storage_t<T>>>> const& compiled) const {
        for (const auto& arg : compiled) {
            if (arg->has_index_space(env)) {
                return true;
            }       
        }
         return false;
     }

     [[nodiscard]] Maybe_unique<base_typed_expression<Lexer, Compiler>>
    compile_set_construction(Environment<Lexer, Compiler> const& cm,
                             const untyped_argument_list<Lexer, Compiler>& args) const {
        using expr_t = typed_argument<Lexer, Compiler, storage_t<T>>;

        std::vector<std::unique_ptr<expr_t>> compiled;
        compiled.reserve(args.arg().size());

        std::string err;

        auto n = args.arg().size();
        for (std::size_t i = 0; i < n; ++i) {
            auto maybe_elem = m_elem.compile_this_argument(cm, args.arg()[i]);
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                compiled.emplace_back(std::move(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
        }

        if (uses_pointwise_arguments(cm, compiled)) {
            return std::make_unique<typed_lifted_set_construction<Lexer, Compiler, T>>(
                std::move(compiled));
        }
        return std::make_unique<typed_set_construction<Lexer, Compiler, T>>(
            std::move(compiled));
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
        return std::make_unique<typed_predicate_evaluation<Lexer, Compiler, F, T>>(
            m_f, std::move(x.value()));
    }

    // base_function_compiler interface

    base_function_compiler<Lexer, Compiler>* clone() const override {
        return new function_compiler(*this);
    }

    void register_types(Compiler& registry) const override {
        registry.template ensure_type_registered<T>();
    }
};


}  // namespace macrodr::dsl

#endif  // LEXER_TYPED_H
