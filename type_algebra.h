#ifndef TYPE_ALGEBRA_H
#define TYPE_ALGEBRA_H

#include <functional>
#include <type_traits>
#include <utility>
#include <variant>

#include "maybe_error.h"
namespace in_progress {
template <class... T>
struct P {
    template <class... T2>
    friend auto operator*(P, P<T2...>) {
        return P<T..., T2...>{};
    }

    template <class... aT>
        requires(std::is_convertible_v<aT, T> && ...)
    auto operator()(aT&&... x) const {
        return std::tuple<T...>(std::forward<aT>(x)...);
    }
};

template <class... T>
struct S {
    template <class aT>
        requires(std::is_convertible_v<aT, T> || ...)
    auto operator()(aT&& x) const {
        return std::variant<T...>(std::forward<aT>(x));
    }

    template <class... T2>
    friend auto operator+(S, S<T2...>) {
        return S<T..., T2...>{};
    }
    template <class... T2>
    friend auto operator+(S, P<T2...>) {
        return S<T..., P<T2...>>{};
    }

    template <class... T2>
    friend auto operator+(S, std::variant<T2...>) {
        return S<T..., T2...>{};
    }

    template <class... T2>
    friend auto operator+(P<T2...>, S) {
        return S<P<T2...>, T...>{};
    }

    template <class T2>
        requires(!is_of_this_template_type_v<T2, P>)
    friend auto operator*(S<T2>, S) {
        return (S<>{} + ... + P<T2, T>{});
    }

    template <class... T2>
    friend auto operator*(S<P<T2...>>, S) {
        return (S<>{} + ... + P<T2..., T>{});
    }
    template <class... T2>
    friend auto operator*(S<T2...>, S) {
        return (S<>{} + ... + (S<T2>{} * S{}));
    }
};

template <class... T>
struct Set_of_Types {
    template <class S>
        requires(std::is_same_v<T, S> || ... || false)
    auto operator||(Set_of_Types<S>) const {
        return Set_of_Types{};
    }

    template <class S>
        requires(!(std::is_same_v<T, S> || ... || false))
    auto operator||(Set_of_Types<S>) const {
        return Set_of_Types<T..., S>{};
    }

    template <template <class...> class Co>
    using transport = Co<T...>;
};

template <class S, template <class...> class Co>
using transport_t = typename S::template transport<Co>;

template <class... T>
using set_of_types_t = decltype((Set_of_Types<>{} || ... || Set_of_Types<T>{}));

}  // namespace in_progress

template <auto x>
struct V {
    constexpr static decltype(x) value = x;
};
template <class T0, class T1, class... Ts>
class constexpr_P {};
template <class T0, class T1, class... Ts>
class constexpr_S {};
template <template <class T0, class T1, class... Ts> class constexpr_P,
          template <class T0, class T1, class... Ts> class constexpr_S>
struct algebra {
    template <class X, class Y, class Z, class... Us>
    constexpr auto P_constexpr(X, Y, Z, Us...) {
        return P_constexpr(P_constexpr(X{}, Y{}), Z{}, Us{}...);
    }

    template <class X, class Y, class Z, class... Us>
    constexpr auto S_constexpr(X, Y, Z, Us...) {
        return S_constexpr(S_constexpr(X{}, Y{}), Z{}, Us{}...);
    }

    template <class X, class Y>
    constexpr auto P_constexpr(X, Y) {
        return constexpr_P<X, Y>{};
    }

    template <class X, class Y0, class Y1, class... Ys>
    constexpr auto P_constexpr(X, constexpr_P<Y0, Y1, Ys...>) {
        return constexpr_P<X, Y0, Y1, Ys...>{};
    }
    template <class X0, class X1, class... Xs, class Y>
    constexpr auto P_constexpr(constexpr_P<X0, X1, Xs...>, Y) {
        return constexpr_P<X0, X1, Xs..., Y>{};
    }

    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    constexpr auto P_constexpr(constexpr_P<X0, X1, Xs...>, constexpr_P<Y0, Y1, Ys...>) {
        return constexpr_P<X0, X1, Xs..., Y0, Y1, Ys...>{};
    }

    template <class X, class Y>
    constexpr auto S_constexpr(X, Y) {
        return constexpr_S<X, Y>{};
    }

    template <class X, class Y0, class Y1, class... Ys>
    constexpr auto S_constexpr(X, constexpr_S<Y0, Y1, Ys...>) {
        return constexpr_S<X, Y0, Y1, Ys...>{};
    }
    template <class X0, class X1, class... Xs, class Y>
    constexpr auto S_constexpr(constexpr_S<X0, X1, Xs...>, Y) {
        return constexpr_S<X0, X1, Xs..., Y>{};
    }

    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    constexpr auto S_constexpr(constexpr_S<X0, X1, Xs...>, constexpr_S<Y0, Y1, Ys...>) {
        return constexpr_S<X0, X1, Xs..., Y0, Y1, Ys...>{};
    }

    template <class X0, class X1, class... Xs, class Y>
    constexpr auto P_constexpr(constexpr_S<X0, X1, Xs...>, Y) {
        return S_constexpr(P_constexpr(X0{}, Y{}), P_constexpr(X1{}, Y{}),
                           P_constexpr(Xs{}, Y{})...);
    }

    template <class X, class Y0, class Y1, class... Ys>
    constexpr auto P_constexpr(X, constexpr_S<Y0, Y1, Ys...>) {
        return S_constexpr(P_constexpr(X{}, Y0{}), P_constexpr(X{}, Y1{}),
                           P_constexpr(X{}, Ys{})...);
    }
};

template <template <class T0, class... Ts> class constexpr_P,
          template <class T0, class... Ts> class constexpr_S>
struct algebra_2 {
    template <class T>
    static constexpr bool is_regular =
        !is_of_this_template_type_v<T, constexpr_P> && !is_of_this_template_type_v<T, constexpr_S>;

    template <class T>
    static constexpr bool is_sum = is_of_this_template_type_v<T, constexpr_S>;

    template <class X, class Y, class... Us>
    struct P_constexpr_impl;

    template <class X, class Y, class... Us>
    using P_constexpr = typename P_constexpr_impl<X, Y, Us...>::type;

    template <class X, class... Us>
    struct S_constexpr_impl;
    template <class X, class... Us>
    using S_constexpr = typename S_constexpr_impl<X, Us...>::type;

    template <class X, class Y, class Z, class... Us>
    struct P_constexpr_impl<X, Y, Z, Us...> {
        using type = P_constexpr<P_constexpr<X, Y>, Z, Us...>;
    };

    template <class X, class Y, class Z, class... Us>
    struct S_constexpr_impl<X, Y, Z, Us...> {
        using type = S_constexpr<S_constexpr<X, Y>, Z, Us...>;
    };

    template <class X, class Y>
        requires is_regular<X> && is_regular<Y>
    struct P_constexpr_impl<X, Y> {
        using type = constexpr_P<X, Y>;
    };

    template <class X, class Y0, class Y1, class... Ys>
        requires is_regular<X>
    struct P_constexpr_impl<X, constexpr_P<Y0, Y1, Ys...>> {
        using type = constexpr_P<X, Y0, Y1, Ys...>;
    };

    template <class X0, class X1, class... Xs, class Y>
        requires is_regular<Y>
    struct P_constexpr_impl<constexpr_P<X0, X1, Xs...>, Y> {
        using type = constexpr_P<X0, X1, Xs..., Y>;
    };

    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct P_constexpr_impl<constexpr_P<X0, X1, Xs...>, constexpr_P<Y0, Y1, Ys...>> {
        using type = constexpr_P<X0, X1, Xs..., Y0, Y1, Ys...>;
    };

    template <class X, class Y>
        requires(!is_sum<X>) && (!is_sum<Y>)
    struct S_constexpr_impl<X, Y> {
        using type = constexpr_S<X, Y>;
    };

    template <class X, class Y0, class Y1, class... Ys>
        requires(!is_sum<X>)
    struct S_constexpr_impl<X, constexpr_S<Y0, Y1, Ys...>> {
        using type = constexpr_S<X, Y0, Y1, Ys...>;
    };
    template <class X0, class X1, class... Xs, class Y>
        requires(!is_sum<Y>)
    struct S_constexpr_impl<constexpr_S<X0, X1, Xs...>, Y> {
        using type = constexpr_S<X0, X1, Xs..., Y>;
    };

    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct S_constexpr_impl<constexpr_S<X0, X1, Xs...>, constexpr_S<Y0, Y1, Ys...>> {
        using type = constexpr_S<X0, X1, Xs..., Y0, Y1, Ys...>;
    };

    template <class X0, class... Xs>
    struct S_constexpr_impl<constexpr_S<X0, Xs...>> {
        using type = constexpr_S<X0, Xs...>;
    };
    template <class X>
        requires(!is_sum<X>)
    struct S_constexpr_impl<X> {
        using type = constexpr_S<X>;
    };

    template <class X, class Y0, class... Ys>
        requires(!is_sum<X>)
    struct P_constexpr_impl<X, constexpr_S<Y0, Ys...>> {
        using type = S_constexpr<P_constexpr<X, Y0>, P_constexpr<X, Ys>...>;
    };
    template <class X0, class... Xs, class Y>
    struct P_constexpr_impl<constexpr_S<X0, Xs...>, Y> {
        using type = S_constexpr<P_constexpr<X0, Y>, P_constexpr<Xs, Y>...>;
    };
};

using res = algebra_2<constexpr_P, constexpr_S>::P_constexpr<constexpr_S<int, double>,
                                                             constexpr_S<int, double>>;

template <class T, T x>
struct constexpr_value {
    constexpr static decltype(x) value = x;
};

template <typename T, T... allowed>
struct constexpr_domain {
    static constexpr std::array<T, sizeof...(allowed)> values = {allowed...};
};

template <template <class...> class>
struct Co {};

template <class... T0>
auto prod_variant_of_tuples_and_variant(std::variant<T0...> const& x) {
    return x;
}
template <class... T0>
auto sum_variant_type(std::variant<T0...> const& x) {
    return x;
}

template <class... T0, class... T1, class... Vs>
    requires(is_of_this_template_type_v<Vs, std::variant> && ...)
auto sum_variant_type(std::variant<T0...> const& x, std::variant<T1...> const& y, const Vs&... vs) {
    return sum_variant_type(std::variant<T0..., T1...>{}, vs...);
}

template <class... T0, class... T1>

auto prod_tuple_variant(std::tuple<T0...> const& x, const std::variant<T1...>& y) {
    // using m7=typename decltype(y)::dentro_prod_tuple_variant;
    return std::visit(
        [&x](auto& e) -> std::variant<std::tuple<T0..., T1>...> {
            return std::apply([&e](auto&... t) { return std::tuple(t..., e); }, x);
        },
        y);
}

template <class... Tu, class... T1, class... Vs>
    requires((is_of_this_template_type_v<Tu, std::tuple> && ...) &&
             (is_of_this_template_type_v<Vs, std::variant> && ...))
auto prod_variant_of_tuples_and_variant(std::variant<Tu...> const& x, const std::variant<T1...>& y,
                                        const Vs&... vs) {
    auto next_x = std::visit(
        [&y](auto& tu) {
            auto xy = prod_tuple_variant(tu, y);
            using variant_all =
                std::decay_t<decltype(sum_variant_type(prod_tuple_variant(Tu{}, y)...))>;
            return std::visit([](auto e) { return variant_all(e); }, xy);
        },
        x);
    return prod_variant_of_tuples_and_variant(next_x, vs...);
}

template <class... T1>
std::variant<std::tuple<T1>...> variant_to_variant_of_tuples(const std::variant<T1...>& y) {
    return std::visit([](auto& e) -> std::variant<std::tuple<T1>...> { return std::tuple(e); }, y);
}

template <class... T1, class... Vs>
    requires(is_of_this_template_type_v<Vs, std::variant> && ...)

auto prod_variant(const std::variant<T1...>& y, const Vs&... vs) {
    return prod_variant_of_tuples_and_variant(variant_to_variant_of_tuples(y), vs...);
}

template <class F, class... T>
auto Apply_tuples(F&& f, std::tuple<T...> const& x) {
    return std::apply([&f](auto&... ts) { return f(ts...); }, x);
}

template <class F, class... Tus>
    requires(is_of_this_template_type_v<Tus, std::tuple> && ...)
auto Apply_variant_of_tuples(F&& f, std::variant<Tus...> const& x) {
    using return_type =
        std::variant<std::decay_t<decltype(Apply_tuples(std::forward<F>(f), Tus{}))>...>;

    return std::visit([&f](auto& tu) { return return_type(Apply_tuples(f, tu)); }, x);
}

template <class F, class... Vs>
    requires(is_of_this_template_type_v<Vs, std::variant> && ...)
auto Apply_variant(F&& f, std::tuple<Vs...> const& x) {
    auto variant_of_tuples = std::apply([](auto&... e) { return prod_variant(e...); }, x);
    //using m5=typename decltype(variant_of_tuples)::despues_variant_of_tuples;

    return Apply_variant_of_tuples(std::forward<F>(f), variant_of_tuples);
}

namespace safe_visitor {

template <class F, class... Ts>
auto safe_visit(F&& f, std::variant<Ts...> const& x) {
    using return_type =
        in_progress::transport_t<in_progress::set_of_types_t<std::invoke_result_t<F, Ts const&>...>,
                                 std::variant>;
    return std::visit(
        [f](auto const& e) { return return_type(std::invoke(std::forward<F>(f), e)); }, x);
}

}  // namespace safe_visitor

namespace in_progress {
template <class F, class... Ts>
auto Apply(F, P<Ts...>) {
    return F{}(Ts{}...);
}

template <class F, class... Ps>
auto Map(F, S<Ps...>) {
    return (S<>{} + ... + Apply(F{}, Ps{}));
}

template <class... Ts>
auto Variant(S<Ts...>) {
    return std::variant<Ts...>{};
}

template <class S>
using Variant_t = decltype(Variant(S{}));

}  // namespace in_progress

#endif  // TYPE_ALGEBRA_H
