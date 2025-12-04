#ifndef DERIVATIVE_OPERATOR_H
#define DERIVATIVE_OPERATOR_H

#include <general_algorithm_on_containers.h>

#include <cassert>
#include <concepts>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>
#include <derivative_fwd.h>
#include "matrix.h"
#include "maybe_error.h"
#include "variables.h"

// Toggleable assertion for derivative dx-pointer checks in the core
// derivative machinery. If MACRODR_STRICT_DX_ASSERT is defined,
// dx-related assertions are enabled; otherwise they are compiled out.
#ifndef MACRODR_DX_ASSERT
  #ifndef MACRODR_STRICT_DX_ASSERT
    #define MACRODR_DX_ASSERT(cond) ((void)0)
  #else
    #define MACRODR_DX_ASSERT(cond) assert(cond)
  #endif
#endif
namespace var {

/**
 * the idea is that we have a derivative object of a variable that takes the place 
 * of the original variable object in any function that uses this variable as argument 
 * 
 * So, we will have two distinct types 
 * 1. with the derivative of one variable against the other
 * 2. with the value of a variable and the derivative against a set of variables. 
 *
 * My problem now is which name to apply to the first and second options....
 * 
 * what did I do in the past?
 * 
 * poorly.., the code is messy. 
 * 
 * 
 * well, I can use capital Der for the derivative plus the primitive, and der just for the derivative...
 *  
 * we can say Taylor, but it gives the wrong idea...No
 * 
 * 
 * Der<X,Y>
 * 
 * 
 * 
 * 
 */


template <>
class d_d_<double, double> {
    double m_dydx;
    double const* ptr_dx=nullptr;

   public:
    using value_type = double;
    constexpr d_d_(double dydx, const double& x) : m_dydx{dydx}, ptr_dx{&x} {}
    constexpr auto& operator()() { return m_dydx; }
    constexpr auto operator()() const { return m_dydx; }
    constexpr d_d_() = default;
    [[nodiscard]] constexpr auto& dx() const { return *ptr_dx; }
    [[nodiscard]] constexpr bool has_dx() const { return ptr_dx != nullptr; }
};

template <>
class d_d_<double, Matrix<double>> {
    Matrix<double> m_dydx;
    Matrix<double> const* ptr_dx=nullptr;

   public:
    using value_type = double;
    d_d_(Matrix<double> const& dydx, const Matrix<double>& x) : m_dydx{dydx}, ptr_dx{&x} {}
    d_d_(Matrix<double>&& dydx, const Matrix<double>& x) : m_dydx{std::move(dydx)}, ptr_dx{&x} {}
    auto& operator()() { return m_dydx; }
    auto operator()() const { return m_dydx; }
    d_d_() = default;
    [[nodiscard]] d_d_(const d_d_&) = default;
    d_d_(d_d_&&) = default;
    d_d_& operator=(const d_d_&) = default;
    d_d_& operator=(d_d_&&) = default;
    auto& dx() const { return *ptr_dx; }
    [[nodiscard]] bool has_dx() const { return ptr_dx != nullptr; }
};

constexpr auto self_derivative(const double& x) {
    return d_d_<double, double>(1.0, x);
}

template <class T>
class Primitive {
    T m_y;

   public:
    template <class aT>
        requires std::is_same_v<T, std::decay_t<aT>>
    constexpr Primitive(aT&& y) : m_y{std::forward<aT>(y)} {}
    constexpr auto& operator()() { return m_y; }
    constexpr auto& operator()() const { return m_y; }
};

template <class T>
class Primitive<Matrix<T>> : public Matrix<T> {
   public:
    template <class aT>
        requires std::is_same_v<Matrix<T>, std::decay_t<aT>>
    constexpr Primitive(aT&& y) : Matrix<T>{std::forward<aT>(y)} {}
    constexpr auto& operator()() { return static_cast<Matrix<T>&>(*this); }
    constexpr auto& operator()() const { return static_cast<Matrix<T> const&>(*this); }
};

template <class N, class D>
class Derivative;

template <class N, class D>
std::ostream& operator<<(std::ostream& os, const Derivative<N, D>& d) {
    os << primitive(d) << "\n derivative: \n" << derivative(d)();
    return os;
}

template <class N, class D>
std::ostream& print(std::ostream& os, const Derivative<N, D>& d) {
    print(os, primitive(d));
    os << "\n derivative: \n";
    print(os, derivative(d)());
    return os;
}

template <class... N, class D>
std::ostream& print(std::ostream& os, const Derivative<var::Vector_Space<N...>, D>& d) {
    ((os << "\n"
         << type_name<N>(),
      ": \n", print(os, primitive(get<N>(d))), os << "\n derivative: \n",
      print(os, derivative(get<N>(d))())),
     ..., 1);
    return os;
}

namespace impl {

template <class F, class X>
struct Derivative_impl {
    using type = std::conditional_t<F::is_constant, F, Derivative<F, X>>;
};

}  // namespace impl

template <class F, class X>
using Derivative_t = impl::Derivative_impl<F, X>::type;

template <>
class Derivative<double, double> {
    double m_x;
    d_d_<double, double> m_d;

   public:
    auto& primitive() const { return m_x; }
    auto& derivative() const { return m_d; }

    constexpr Derivative(double fx, double dfdx, const double& dx) : m_x{fx}, m_d{dfdx, dx} {}

    constexpr Derivative(double x, const double& dx) : m_x{x}, m_d{0.0, dx} {}
    Derivative() = default;
    auto& dx() const {
        MACRODR_DX_ASSERT(m_d.has_dx() &&
                          "Derivative<double,double>::dx() called without valid dx pointer");
        return m_d.dx();
    }
    [[nodiscard]] bool has_dx() const { return m_d.has_dx(); }

    friend auto exp(const Derivative& x) {
        auto f = exp(x.primitive());
        MACRODR_DX_ASSERT(x.derivative().has_dx() &&
                          "exp(Derivative<double,double>): missing dx pointer");
        return Derivative(f, f * x.derivative()(), x.dx());
    }

    friend auto max(const Derivative& x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }

    friend auto min(const Derivative& x, double y) {
        if (x.primitive() >= y) {
            return x;
        } else {
            MACRODR_DX_ASSERT(x.derivative().has_dx() &&
                              "min(Derivative<double,double>): missing dx pointer");
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
        }
    }

    friend auto log(const Derivative& x) {
        auto f = log(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / x.primitive()), x.dx());
    }
    friend auto log10(const Derivative& x) {
        auto f = log10(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / (x.primitive() * std::log(10))));
    }

    friend auto pow(double base, const Derivative& x) {
        using std::pow;
        auto f = pow(base, x.primitive());
        return Derivative(f, x.derivative()() * f * std::log(base), x.dx());
    }
    friend auto pow(const Derivative& base, const Derivative& x) {
        using std::pow;
        auto f = pow(base.primitive(), x.primitive());
        return Derivative(
            f,
            x.derivative()() * f * std::log(base.primitive()) +
                base.derivative()() * x.primitive() * pow(base.primitive(), x.primitive() - 1.0),
            x.dx());
    }

    friend auto abs(const Derivative& x) {
        auto f = std::abs(x.primitive());
        return Derivative(
            f,
            ((x.primitive() > 0.0) ? 1.0 : ((x.primitive() < 0) ? -1.0 : 0.0)) * x.derivative()(),
            x.dx());
    }
};

template <>
class Derivative<double, Matrix<double>> {
    double m_x;
    d_d_<double, Matrix<double>> m_d;

   public:
    auto& primitive() const { return m_x; }
    auto& derivative() const { return m_d; }

    Derivative(double fx, Matrix<double> const& dfdx, const Matrix<double>& dx)
        : m_x{fx}, m_d{dfdx, dx} {}
    Derivative(double fx, Matrix<double>&& dfdx, const Matrix<double>& dx)
        : m_x{fx}, m_d{std::move(dfdx), dx} {}

    Derivative(double x, const Matrix<double>& dx) : m_x{x}, m_d{dx - dx, dx} {}
    Derivative() = default;

    auto& dx() const {
        MACRODR_DX_ASSERT(m_d.has_dx() &&
                          "Derivative<double,Matrix<double>>::dx() called without valid dx pointer");
        return m_d.dx();
    }

    friend auto exp(const Derivative& x) {
        auto f = exp(x.primitive());
        MACRODR_DX_ASSERT(x.derivative().has_dx() &&
                          "exp(Derivative<double,Matrix<double>>): missing dx pointer");
        return Derivative(f, f * x.derivative()(), x.dx());
    }

    friend auto max(const Derivative& x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }
    friend auto max(const Derivative& x, const Derivative& y) {
        if (x.primitive() < y.primitive()) {
            return y;
        }
        return x;
    }

    friend auto min(const Derivative& x, const Derivative& y) {
        if (x.primitive() < y.primitive()) {
            return y;
        }
        return x;
    }

    friend auto min(const Derivative& x, double y) {
        if (x.primitive() >= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }

    friend auto log(const Derivative& x) {
        auto f = log(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / x.primitive()), x.dx());
    }
    friend auto log10(const Derivative& x) {
        auto f = log10(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / (x.primitive() * std::log(10))));
    }
    friend auto pow(double base, const Derivative& x) {
        using std::pow;
        auto f = pow(base, x.primitive());
        return Derivative(f, x.derivative()() * f * std::log(base), x.dx());
    }
    friend auto pow(const Derivative& base, const Derivative& x) {
        using std::pow;
        auto f = pow(base.primitive(), x.primitive());
        return Derivative(
            f,
            x.derivative()() * f * std::log(base.primitive()) +
                base.derivative()() * x.primitive() * pow(base.primitive(), x.primitive() - 1.0),
            x.dx());
    }

    friend auto abs(const Derivative& x) {
        auto f = std::abs(x.primitive());
        return Derivative(
            f,
            ((x.primitive() > 0.0) ? 1.0 : ((x.primitive() < 0) ? -1.0 : 0.0)) * x.derivative()(),
            x.dx());
    }
};


template <class N, class D>
Maybe_error<bool> compare_contents(Derivative<N, D> const& s0, Derivative<N, D> const& s1,
                                   double = 0, double = 0, std::size_t = 1) {
    return compare_contents(primitive(s0), primitive(s1));
}

template <class...>
struct dx_of_dfdx;

template <class F, class X>
struct dx_of_dfdx<Derivative<F, X>> {
    using type = X;
};

template <class F>
    requires(!is_derivative_v<F>)
struct dx_of_dfdx<F> {
    using type = NoDerivative;
};
template <>
struct dx_of_dfdx<> {
    using type = NoDerivative;
};

template <class F, class X, class G, class... Ds>
struct dx_of_dfdx<Derivative<F, X>, Derivative<G, X>, Ds...> {
    using type = typename dx_of_dfdx<Derivative<F, X>, Ds...>::type;
};

template <class F, class X, class G, class... Ds>
    requires(!is_derivative_v<G>)
struct dx_of_dfdx<Derivative<F, X>, G, Ds...> {
    using type = typename dx_of_dfdx<Derivative<F, X>, Ds...>::type;
};

template <class F, class X, class G, class... Ds>
    requires(!is_derivative_v<G>)
struct dx_of_dfdx<G, Derivative<F, X>, Ds...> {
    using type = typename dx_of_dfdx<Derivative<F, X>, Ds...>::type;
};

template <class F, class G, class... Ds>
    requires(!is_derivative_v<G> && !is_derivative_v<F>)
struct dx_of_dfdx<G, F, Ds...> {
    using type = typename dx_of_dfdx<Ds...>::type;
};

template <class X, class Y>
struct dx_of_dfdx<d_d_<X, Y>> {
    using type = Y;
};

template <class... T>
using dx_of_dfdx_t = typename std::decay_t<dx_of_dfdx<T...>>::type;

inline auto get_dx_of_dfdx() {
    return NoDerivative{};
}

template <class T, class... Ts>
    requires(!is_derivative_v<T>&&...&&!is_derivative_v<Ts>)
auto get_dx_of_dfdx(const T&, const Ts&... ) {
    return NoDerivative{};
}


template <class T,class T0, class... Ts>
    requires(!is_derivative_v<T>&&(is_derivative_v<T0>||...||is_derivative_v<Ts>))
auto const & get_dx_of_dfdx(const T&, const T0& y,const Ts&... xs) {
    return get_dx_of_dfdx(y,xs...);
}






template <class F, class X, class T,class... Ts>
    requires(is_derivative_v<T>||...||is_derivative_v<Ts>)
auto const& get_dx_of_dfdx(const Derivative<F, X>& x,  const T& y, const Ts&...z) {
    if (x.has_dx()) {
        return get_dx_of_dfdx(x);
    }
    return get_dx_of_dfdx(y,z...);
}
template <class F, class X, class T,class... Ts>
    requires(!is_derivative_v<T>&&...&&!is_derivative_v<Ts>)
auto const& get_dx_of_dfdx(const Derivative<F, X>& x,  const T& , const Ts&...) {
        return get_dx_of_dfdx(x);
}



//auto a=Derivative<double,double>(9.9,9.9);

//auto b=primitive(a);

template <class var, class T, class X>
    requires(var::is_variable && std::constructible_from<var, T>)
auto build(Derivative<T, X>&& x) {
    return Derivative<var, X>(std::move(x));
}

template <class var, class... T>
    requires(std::constructible_from<var, untransformed_type_t<T>...> &&
             (is_derivative_v<T> || ... || false) && (!is_Maybe_error<T> && ... && true))
auto build(T... x) {
    using X = dx_of_dfdx_t<T...>;
    return Derivative<var, X>(std::forward<T>(x)...);
}

template <class var, class... T>
    requires(std::constructible_from<var, untransformed_type_t<T>...> &&
             (is_derivative_v<T> || ... || false) && (is_Maybe_error<T> || ... || false))
auto build(T... x) {
    using X = dx_of_dfdx_t<T...>;
    if ((is_valid(x) && ... && true))
        return Maybe_error<Derivative<var, X>>(std::forward<T>(x)...);
    else
        return Maybe_error<Derivative<var, X>>(error_message((get_error(x)() + ... + "")));
}

template <class var, class... T>
    requires(std::constructible_from<var, untransformed_type_t<T>...> &&
             (!is_derivative_v<T> && ... && true) && (is_Maybe_error<T> || ... || false))
auto build(T... x) {
    if ((is_valid(x) && ... && true))
        return Maybe_error<var>(std::forward<T>(x)...);
    else
        return Maybe_error<var>(error_message((get_error(x)() + ... + "")));
}

template <template <class...> class Vector_Space, class... T>
    requires(std::constructible_from<Vector_Space<untransformed_type_t<T>...>,
                                     untransformed_type_t<T>...> &&
             (is_derivative_v<T> || ... || false) && (!is_Maybe_error<T> && ... && true))
auto build(T... x) {
    using X = dx_of_dfdx_t<T...>;
    return Derivative<Vector_Space<untransformed_type_t<T>...>, X>(std::forward<T>(x)...);
}

template <template <class...> class Vector_Space, class... T>
    requires(std::constructible_from<Vector_Space<untransformed_type_t<T>...>,
                                     untransformed_type_t<T>...> &&
             (is_derivative_v<T> || ... || false) && (is_Maybe_error<T> || ... || false))
auto build(T... x) {
    using X = dx_of_dfdx_t<T...>;
    if ((is_valid(x) && ... && true))
        return Maybe_error<Derivative<Vector_Space<untransformed_type_t<T>...>, X>>(
            std::forward<T>(x)...);
    else
        return error_message((get_error(x)() + ... + ""));
}

template <template <class...> class Vector_Space, class... T>
    requires(std::constructible_from<Vector_Space<untransformed_type_t<T>...>,
                                     untransformed_type_t<T>...> &&
             (!is_derivative_v<T> && ... && true) && (!is_Maybe_error<T> && ... && true))
auto build(T... x) {
    return Vector_Space<untransformed_type_t<T>...>(std::forward<T>(x)...);
}

template <template <class...> class Vector_Space, class... T>
    requires(std::constructible_from<Vector_Space<untransformed_type_t<T>...>,
                                     untransformed_type_t<T>...> &&
             (!is_derivative_v<T> && ... && true) && (is_Maybe_error<T> || ... || true))
auto build(T... x) {
    if ((is_valid(x) && ... && true))
        return Maybe_error<Vector_Space<untransformed_type_t<T>...>>(std::forward<T>(x)...);
    else
        return error_message((get_error(x)() + ... + ""));
}

template <template <class...> class varr, class Id, class F, class... T>
    requires(std::constructible_from<varr<Id, F, T...>, var::Var<Id>, F, T...> &&
             (std::is_same_v<varr<Id, T...>, Fun<Id, T...>>) &&
             (!is_derivative_v<T> && ... && true) && (!is_Maybe_error<T> && ... && true))
varr<Id, F, std::decay_t<T>...> build(Var<Id>, F&& t_f, T&&... t_x) {
    return varr(Var<Id>{}, std::forward<F>(t_f), std::forward<T>(t_x)...);
}

template <template <class...> class var, class Id, class F, class... T>
    requires(std::constructible_from<var<Id, F, untransformed_type_t<T>...>, Var<Id>, F,
                                     untransformed_type_t<T>...> &&
             (std::is_same_v<var<Id, F, untransformed_type_t<T>...>,
                             Fun<Id, F, untransformed_type_t<T>...>>) &&
             (is_derivative_v<T> || ... || false) && (!is_Maybe_error<T> && ... && true))
auto build(Var<Id>, F&& t_f, T&&... t_x) {
    using X = dx_of_dfdx_t<std::decay_t<T>...>;
    return Derivative<Fun<Id, F, untransformed_type_t<std::decay_t<T>>...>, X>(
        Var<Id>{}, std::forward<F>(t_f), std::forward<T>(t_x)...);
}

template <class T, class S>
    requires(is_derivative_v<T> || is_derivative_v<S>)
auto operator*(const T& x, const S& y) {
    using X = dx_of_dfdx_t<T, S>;
    using F = decltype(primitive(x) * primitive(y));

    return Derivative<F, X>(primitive(x) * primitive(y),
                            derivative(x)() * primitive(y) + primitive(x) * derivative(y)(),
                            get_dx_of_dfdx(x, y));
}

template <class T, class S>
//    requires(is_derivative_v<T>||is_derivative_v<S>)
auto max(const T& x, const S& y) {
    using std::max;
    return max(primitive(x), primitive(y));
}

template <class T, class S>
    requires(is_derivative_v<T> || is_derivative_v<S>)
auto operator+(const T& x, const S& y) {
    using X = dx_of_dfdx_t<T, S>;
    using F = decltype(primitive(x) + primitive(y));

    return Derivative<F, X>(primitive(x) + primitive(y), derivative(x)() + derivative(y)(),
                            get_dx_of_dfdx(x, y));
}

template <class T, class S>
    requires(is_derivative_v<T> && is_derivative_v<S>)
auto operator-(const T& x, const S& y) {
    using X = dx_of_dfdx_t<T, S>;
    using F = decltype(primitive(x) - primitive(y));

    return Derivative<F, X>(primitive(x) - primitive(y), derivative(x)() - derivative(y)(),
                            get_dx_of_dfdx(x, y));
}

template <class T, class S>
    requires(is_derivative_v<T> && (!is_derivative_v<S>))
auto operator-(const T& x, const S& y) {
    using X = dx_of_dfdx_t<T, S>;
    using F = decltype(primitive(x) - y);

    return Derivative<F, X>(primitive(x) - y, derivative(x)(), x.dx());
}

template <class T, class S>
    requires(!is_derivative_v<T> && is_derivative_v<S>)
auto operator-(const T& x, const S& y) {
    using X = dx_of_dfdx_t<T, S>;
    using F = decltype(x - primitive(y));

    return Derivative<F, X>(x - primitive(y), derivative(y)() * (-1.0), y.dx());
}

// Helper: initialize a value `x` either as a plain primitive (when `input`
// does not live in the derivative universe) or as a Derivative<...,Dir>
// sharing the same dx() pointer as `input`.
//
// This is the disciplined entry point for "always have dx by construction":
// - If `input` is not a var::Derivative, we just return `x`.
// - If `input` is a var::Derivative and `T` is primitive, we lift `x` into
//   Derivative<T, dx_of_dfdx_t<Input>> using `input.dx()`.
//
// More specialized lifting (e.g. for Vector_Space-valued outputs) can be
// added later by providing additional overloads or partial specializations.

// Primitive path: input is not a derivative → keep T as-is.
template <class Input, class T>
    requires(!is_derivative_v<std::decay_t<Input>>)
auto initialize(T&& x, const Input&) {
    return std::forward<T>(x);
}

// Derivative path: input is a derivative (or wraps one) and T is primitive.
template <class Input, class T>
    requires(is_derivative_v<std::decay_t<Input>> && !is_derivative_v<std::decay_t<T>>)
auto initialize(T&& x, const Input& in) {
    using Dir = dx_of_dfdx_t<std::decay_t<Input>>;
    static_assert(!std::is_same_v<Dir, NoDerivative>,
                  "initialize(Input,T): derivative input has no dx_of_dfdx_t direction");

    using F = std::decay_t<T>;
    return Derivative<F, Dir>(std::forward<T>(x), in.dx());
}

auto max(is_Container auto const& c) {
    auto out = c[0];
    for (std::size_t i = 1; i < c.size(); ++i)
        if (std::isnan(primitive(c[i])) || primitive(out) < primitive(c[i]))
            out = c[i];
    return out;
}

}  // namespace var

#endif  // DERIVATIVE_OPERATOR_H
