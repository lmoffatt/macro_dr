#pragma once

#ifndef DYNAMIC_CLI_H
#define DYNAMIC_CLI_H

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "maybe_error.h"
namespace dsl {

using logic::error_message;
using logic::Maybe_error;

class expression {
   public:
    virtual ~expression() {};
    virtual std::string str() const = 0;
};

namespace det {}  // namespace det

class typed_expression : public expression {
   public:
    virtual ~typed_expression() {};
    virtual Identifier type() const = 0;
};

template <typename T>
struct T_s;

template <>
struct T_s<int> {
    Identifier name() const {
        return *to_Identifier("integer");
    }
};

template <class... Ts>
struct Cs {};

template <typename T>
class Value : public typed_expression {
    T m_x;

   public:
    Value(const T& x) : m_x{x} {
    }
    virtual ~Value() override {};
    virtual Identifier type() const {
        return T_s<T>::name();
    }

    auto& operator()() const {
        return m_x;
    }
    auto& operator()() {
        return m_x;
    }
};

template <typename T>
class typed_Identifier : public typed_expression {
    Identifier m_id;

   public:
    typed_Identifier(Identifier t_id) : m_id{std::move(t_id)} {
    }
    virtual ~typed_Identifier() override {};
    virtual Identifier type() const {
        return T_s<T>::name();
    }

    auto& operator()() const {
        return m_id;
    }
    auto& operator()() {
        return m_id;
    }
};

class Base_function {
   public:
    virtual ~Base_function() {};
};

template <class T>
auto Identifiers_to_string(const std::map<Identifier, T>& ids) {
    std::string out = "";
    for (auto& ele : ids) out += ele.first();
    return out;
}

template <typename T>
class Arg {
    Identifier m_name;

   public:
    Arg(Identifier t_name) : m_name{t_name} {
    }
    auto& name() const {
        return m_name;
    }

    Maybe_error<T> extract(std::map<Identifier, std::unique_ptr<base_value>> const& args) {
        auto it = args.find(name());
        if (it == args.end())
            return error_message(name() + " argument was not found in " +
                                 Identifiers_to_string(args));
        else if (T_s<T>::name() != it->second->type())
            return error_message(name() + " argument with unexpercted type " + it->second->name() +
                                 " instead of " + T_s<T>::name());
        else
            return reinterpret_cast<Value<T>&>(*(it->second))();
    }
};

template <typename T>
class Arg_v : public Arg<T> {
    T m_default;

   public:
    Arg_v(Identifier t_name, T t_default) : Arg<T>{t_name}, m_default{std::move(t_default)} {
    }

    auto& name() const {
        return Arg<T>::name();
    }
    auto& value() const {
        return m_default;
    }
    Maybe_error<T> extract(std::map<Identifier, std::unique_ptr<base_value>> const& args) {
        auto it = args.find(name());
        if (it == args.end())
            return value();
        else if (T_s<T>::name() != it->second->type())
            return error_message(name() + " argument with unexpercted type " + it->second->name() +
                                 " instead of " + T_s<T>::name());
        else
            return reinterpret_cast<Value<T>&>(*(it->second))();
    }
};

template <class...>
class function;

template <class F, typename... As, typename... A_ds>
class function<F, Cs<As...>, Cs<A_ds...>> : public Base_function {
    Identifier m_name;
    F m_f;
    std::tuple<Arg<As>...> m_args;
    std::tuple<Arg_v<A_ds>...> m_args_v;

    using R = std::invoke_result_t<F, As..., A_ds...>;

    Maybe_error<R> evaluate(Maybe_error<As>... t_arg, Maybe_error<A_ds>... t_arg_d) {
        if (((t_arg) && ... && true) && ((t_arg_d) && ... && true))
            return std::invoke(m_f, *(t_arg)..., *(t_arg_d)...);
        else
            return error_message((t_arg.error()() + ...) + ((t_arg_d.error()() + ...)));
    }

   public:
    function(Identifier fun_name, F t_f, Arg<As>... t_args, Arg_v<A_ds>... t_args_v)
        : m_name{fun_name}, m_f{t_f}, m_args{t_args...}, m_args_v{t_args_v...} {
    }

    auto evaluate(std::map<Identifier, std::unique_ptr<base_value>> const& args) const {
        return std::apply(
            [this, args](Arg<As> const&... t_arg) {
                return std::apply(
                    [this, &args, &t_arg...](Arg_v<A_ds>&... t_arv_v) {
                        return evaluate(t_arg.extract(args)..., t_arv_v.extract(args)...);
                    },
                    m_args_v);
            },
            m_args);
    }
};

template <class F, typename... As, typename... A_ds>
function(Identifier fun_name, F t_f, Arg<As>... t_args, Arg_v<A_ds>... t_args_v)
    -> function<F, As..., A_ds...>;

}  // namespace dsl

#endif  // DYNAMIC_CLI_H
