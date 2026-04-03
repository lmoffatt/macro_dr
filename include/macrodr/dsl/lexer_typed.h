#ifndef LEXER_TYPED_H
#define LEXER_TYPED_H
//#include "grammar_typed.h"
#include <algorithm>
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
//#include "grammar_typed.h"
//#include "grammar_typed.h"
#include <indexed.h>
#include <macrodr/interface/IModel.h>
#include <macrodr/io/json/convert.h>

#include "maybe_error.h"
namespace macrodr::dsl {
// JSON specification hook: by default, map to macrodr::io::json::conv
template <class L>
struct json_spec {
    using TagPolicy = macrodr::io::json::conv::TagPolicy;
    using Json = macrodr::io::json::Json;
    template <class T>
    static Json to_json(const T& v, TagPolicy policy) {
        return macrodr::io::json::conv::to_json(v, policy);
    }
    template <class T>
    static Maybe_error<void> from_json(const Json& j, T& out, const std::string& path,
                                       TagPolicy policy) {
        return macrodr::io::json::conv::from_json(j, out, path, policy);
    }
};
namespace detail {

template <class T>
struct is_unique_ptr : std::false_type {};

template <class T, class Deleter>
struct is_unique_ptr<std::unique_ptr<T, Deleter>> : std::true_type {};

template <class Arg, class Enable = void>
struct function_argument_storage {
    using type = std::remove_cvref_t<Arg>;
};

template <class T>
struct function_argument_storage<const T&> {
    using type = std::reference_wrapper<const T>;
};

template <class T>
struct function_argument_storage<T&> {
    using type = std::reference_wrapper<T>;
};

template <class Arg>
using function_argument_storage_t = typename function_argument_storage<Arg>::type;

template <class Param, class Storage>
decltype(auto) adapt_argument_like(Storage& storage) {
    using Base = std::remove_reference_t<Param>;
    using StorageType = std::remove_reference_t<Storage>;

    if constexpr (std::is_same_v<StorageType, std::reference_wrapper<Base>> ||
                  std::is_same_v<StorageType,
                                 std::reference_wrapper<std::remove_const_t<Base>>>) {
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

}  // namespace detail

template <class Lexer, class Compiler>
class Environment;
template <class Lexer, class Compiler>
class base_typed_expression;
template <class Lexer, class Compiler>
class base_Identifier_compiler;
template <class Lexer, class Compiler, class T>
class typed_expression;

template <class Lexer, class Compiler>
class base_typed_argument;

template <class Lexer, class Compiler, class T>
class typed_argument;

template <class Lexer, class Compiler, class T>
class expression_argument;

template <class Lexer, class Compiler, class T>
class pointwise_argument;

template <class Lexer, class Compiler, class T>
class borrowed_expression_argument;

template <class Lexer, class Compiler, class T>
class pointwise_borrowed_argument;

template <class Lexer, class Compiler, class T>
class typed_vector_construction;
template <class Lexer, class Compiler, class T>
class typed_set_construction;

template <class Lexer, class Compiler, class... T>
class typed_tuple_construction;

template <class Lexer, class Compiler, class T>
class typed_identifier;

template <class Lexer, class Compiler, class T>
class typed_identifier_ref_const;

template <class Lexer, class Compiler, class T>
class typed_identifier_ref;

template <class Lexer, class Compiler, class T>
class typed_literal;

template <class Lexer, class Compiler, class T>
Maybe_error<void> load_literal_from_json_helper(const typename json_spec<Lexer>::Json& value,
                                                const std::string& path,
                                                typename json_spec<Lexer>::TagPolicy policy,
                                                const Identifier<Lexer>& id,
                                                Environment<Lexer, Compiler>& env);

template <class Lexer, class Compiler>
class untyped_argument_list;
class Lexer;

template <class Lexer, class Compiler, class T, class S>
    requires(std::convertible_to<S, T>)
class typed_conversion;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_scalar_function_evaluation;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_lifted_function_evaluation;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_exact_indexed_function_evaluation;

template <class Lexer, class Compiler>
class untyped_program;

template <class Lexer, class Compiler>
class typed_program;

template <class Lexer, class Compiler, class T>
class typed_expression;
template <class Lexer, class Compiler>
class untyped_expression;

template <class Lexer, class Compiler>
class untyped_statement;

template <class Lexer, class Compiler>
class untyped_numeric_literal;

template <class Lexer, class Compiler>
class untyped_string_literal;

template <class Lexer, class Compiler>
class untyped_identifier;

template <class Lexer, class Compiler>
class untyped_assignment;

template <class Lexer, class Compiler>
class untyped_literal;

template <class Lexer, class Compiler>
class untyped_function_evaluation;

template <class Lexer, class Compiler>
class untyped_vector_construction;

template <class Lexer, class Compiler>
class untyped_tuple_construction;

template <class Lexer, class Compiler, class P, class T>
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<P, T>>)
class typed_predicate_evaluation;


template <class Lexer, class Compiler,  class... Ts>
 class tuple_compiler;

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
class vector_compiler;
template <class Lexer, class Compiler, class T>
class set_compiler;

template <class Lexer, class Compiler, class T>
class parameter_compiler {
    Identifier<Lexer> m_id;

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, T>> adapt_compiled_expression(
        std::unique_ptr<base_typed_expression<Lexer, Compiler>> expr) const {
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
        return error_message(std::string("unexpected type: expected ") + expected + ", got " +
                             actual);
    }

   public:
    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

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
            return error_message(std::string{"unexpected type: expected "} + type_name<T>() +
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
        return new expression_argument<Lexer, Compiler, T>(maybe_expr.value().release());
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
        return new expression_argument<Lexer, Compiler, T>(maybe_expr.value().release());
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
        return new expression_argument<Lexer, Compiler, T>(maybe_expr.value().release());
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
class parameter_compiler<Lexer, Compiler, var::Indexed<T>> {
    Identifier<Lexer> m_id;

    [[nodiscard]] Maybe_unique<typed_argument<Lexer, Compiler, var::Indexed<T>>>
    adapt_compiled_expression(std::unique_ptr<base_typed_expression<Lexer, Compiler>> expr) const {
        if (auto ptr =
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.get())) {
            return new expression_argument<Lexer, Compiler, var::Indexed<T>>(
                dynamic_cast<typed_expression<Lexer, Compiler, var::Indexed<T>>*>(expr.release()));
        }
        const auto expected = type_name<var::Indexed<T>>();
        const auto actual = expr ? expr->type_name() : std::string{"<null>"};
        return error_message(std::string("unexpected type: expected ") + expected + ", got " +
                             actual);
    }

   public:
    parameter_compiler(Identifier<Lexer>&& x) : m_id(std::move(x)) {}
    parameter_compiler(Identifier<Lexer> const& x) : m_id(x) {}

    [[nodiscard]] auto& id() const { return m_id; }

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
class parameter_compiler<Lexer, Compiler, std::reference_wrapper<const T>> {
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
        return error_message(std::string("unexpected type: expected ") + expected + ", got " +
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
            return error_message(std::string{"unexpected type: expected "} + type_name<T>() +
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
            return error_message(std::string{"unexpected type: expected "} + type_name<T>() +
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
        std::tuple<element_compiler<Lexer, Compiler, storage_t<Ts>>...>;

    arg_compilers_t m_args;

    template <std::size_t... Is>
    [[nodiscard]] Maybe_unique<
        typed_expression<Lexer, Compiler, std::tuple<Ts...>>>
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
                std::get<Is>(m_args).compile_this_element(env, args.arg()[Is])...
            ));

        if (!maybe_args) return maybe_args.error();

        // Construct typed_tuple_construction<Ts...>
        return std::make_unique<
            typed_tuple_construction<Lexer, Compiler, Ts...>>(
                std::move(maybe_args.value()));
    }

public:

    [[nodiscard]] Maybe_unique<
        typed_expression<Lexer, Compiler, std::tuple<Ts...>>>
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

    using element_compiler_t = element_compiler<Lexer, Compiler, storage_t<T>>;

    element_compiler_t m_elem;

   public:
    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::vector<T>>>
    compile_vector_construction(Environment<Lexer, Compiler> const& cm,
                                const untyped_argument_list<Lexer, Compiler>& args) const {
        using expr_t = typed_expression<Lexer, Compiler, storage_t<T>>;

        std::vector<std::unique_ptr<expr_t>> compiled;
        compiled.reserve(args.arg().size());

        std::string err;

        auto n = args.arg().size();
        for (std::size_t i = 0; i < n; ++i) {
            auto maybe_elem = m_elem.compile_this_element(cm, args.arg()[i]);
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                compiled.emplace_back(std::move(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
        }

        return std::make_unique<typed_vector_construction<Lexer, Compiler, T>>(
            std::move(compiled));
    }
};

template <class Lexer, class Compiler, class T>
class set_compiler {
    template <class K>
    using storage_t = detail::function_argument_storage_t<K>;

    using element_compiler_t = element_compiler<Lexer, Compiler, storage_t<T>>;

    element_compiler_t m_elem;

   public:
    [[nodiscard]] Maybe_unique<typed_expression<Lexer, Compiler, std::set<T>>>
    compile_set_construction(Environment<Lexer, Compiler> const& cm,
                             const untyped_argument_list<Lexer, Compiler>& args) const {
        using expr_t = typed_expression<Lexer, Compiler, storage_t<T>>;

        std::vector<std::unique_ptr<expr_t>> compiled;
        compiled.reserve(args.arg().size());

        std::string err;

        auto n = args.arg().size();
        for (std::size_t i = 0; i < n; ++i) {
            auto maybe_elem = m_elem.compile_this_element(cm, args.arg()[i]);
            if (!maybe_elem) {
                err += std::to_string(i) + ": " + maybe_elem.error()();
            } else {
                compiled.emplace_back(std::move(maybe_elem.value()));
            }
        }

        if (!err.empty()) {
            return error_message(err);
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

class Compiler {
    using Json = macrodr::io::json::Json;
    using TagPolicy = macrodr::io::json::conv::TagPolicy;

    struct TypeEntry {
        using LoadFn = Maybe_error<void> (*)(const Json&, const std::string&, TagPolicy,
                                             const Identifier<Lexer>&,
                                             Environment<Lexer, Compiler>&);
        LoadFn load = nullptr;
    };

    // Support overloaded functions: multiple compilers per identifier
    std::map<Identifier<Lexer>,
             std::vector<std::unique_ptr<base_function_compiler<Lexer, Compiler>>>>
        m_func;
    std::map<std::string, TypeEntry> m_type_registry;

   public:
    Compiler() { register_builtin_types(); }
    // Compatibility ctor: accept a single-implementation map and wrap entries as one-overload vectors
    Compiler(std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>>&&
                 func) {
        for (auto& [name, fn] : func) {
            auto& vec = m_func[name];
            std::unique_ptr<base_function_compiler<Lexer, Compiler>> moved = std::move(fn);
            if (moved) {
                moved->register_types(*this);
            }
            vec.emplace_back(std::move(moved));
        }
        register_builtin_types();
    }

    Compiler(const Compiler& cm) : m_type_registry{cm.m_type_registry} {
        // Deep-clone the overload vectors per identifier
        for (const auto& [name, vec] : cm.m_func) {
            auto& dst_vec = m_func[name];
            dst_vec.reserve(vec.size());
            for (const auto& fn : vec) {
                dst_vec.emplace_back(fn ? fn->clone() : nullptr);
            }
        }
    }

    [[nodiscard]] Maybe_error<base_function_compiler<Lexer, Compiler> const*> get_function(
        const Identifier<Lexer>& id) const {
        auto it = m_func.find(id);
        if (it == m_func.end()) {
            return error_message(id() + " function is not defined");
        }
        if (it->second.empty() || it->second.front().get() == nullptr) {
            return error_message(id() + " function is null");
        }
        return it->second.front().get();
    }

    [[nodiscard]] Maybe_error<std::vector<base_function_compiler<Lexer, Compiler> const*>>
        get_functions(const Identifier<Lexer>& id) const {
        auto it = m_func.find(id);
        if (it == m_func.end()) {
            return error_message(id() + " function is not defined");
        }
        std::vector<base_function_compiler<Lexer, Compiler> const*> out;
        out.reserve(it->second.size());
        for (const auto& fn : it->second) {
            out.push_back(fn.get());
        }
        return out;
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
        m_func[may_id.value()].emplace_back(fun);
        return true;
    }
    Maybe_error<bool> push_function(
        std::string id_candidate, std::unique_ptr<base_function_compiler<Lexer, Compiler>>&& fun) {
        auto may_id = to_Identifier<Lexer>(std::move(id_candidate));
        if (!may_id) {
            return may_id.error();
        }
        fun->register_types(*this);
        m_func[may_id.value()].emplace_back(std::move(fun));
        return true;
    }
    void merge(const Compiler& other) {
        // Merge functions: clone and append overloads
        for (const auto& [name, vec] : other.m_func) {
            auto& dst_vec = m_func[name];
            for (const auto& fn : vec) {
                if (fn) {
                    auto cloned =
                        std::unique_ptr<base_function_compiler<Lexer, Compiler>>(fn->clone());
                    cloned->register_types(*this);
                    dst_vec.emplace_back(std::move(cloned));
                }
            }
        }
        for (const auto& [key, entry] : other.m_type_registry) {
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, entry);
            }
        }
    }
    void merge(Compiler&& other) {
        // Move functions: append overload vectors
        for (auto& [name, vec] : other.m_func) {
            auto& dst_vec = m_func[name];
            for (auto& fn : vec) {
                if (fn) {
                    fn->register_types(*this);
                }
                dst_vec.emplace_back(std::move(fn));
            }
        }
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
