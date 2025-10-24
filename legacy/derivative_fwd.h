#pragma once

#include <concepts>
#include <vector>
#include <maybe_error.h>
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

template <class>
struct untransformed_type;

template <class>
struct transformation_type;

template <class T>
using untransformed_type_t = typename untransformed_type<T>::type;

template <class T>
using transformation_type_t = typename transformation_type<T>::type;

struct Identity_Op {
    template <class T>
    using type = T;
};

template <class T>
struct untransformed_type {
    using type = T;
};

template <class T>
struct transformation_type {
    using type = Identity_Op;
};

template <class Op, class F>
using Op_t = typename Op::template type<F>;

template <class Op_on_F, class G>
using Transfer_Op_to = Op_t<transformation_type_t<Op_on_F>, G>;

template <class Tr, class T>
concept U = std::same_as<untransformed_type_t<Tr>, T>;

template <class, class>
class d_d_;

template <class T>
constexpr d_d_<T, T> self_derivative(const T&);

template <class, class>
class Derivative;

namespace detail {

// Primary: non-container → Derivative<F, X>
template <class F, class X>
struct derivative_op_impl {
    using type = Derivative<F, X>;
};

// Specialization: std::vector<T, Alloc> → std::vector<Derivative<T, X>>
template <class T, class Alloc, class X>
struct derivative_op_impl<std::vector<T, Alloc>, X> {
    using type = std::vector<Derivative<T, X>>;
};
// Forwarding wrapper to strip cv/ref before dispatching to impl
template <class F, class X>
struct derivative_op : derivative_op_impl<std::remove_cvref_t<F>, X> {};

}  // namespace detail

template <class X>
struct Derivative_Op {
    template <class F>
    using type = typename detail::derivative_op<F, X>::type;
};

struct Maybe_error_Op {
    template <class X>
    using type = Maybe_error<X>;
};

template <class X>
struct Maybe_error_Derivative_Op {
    template <class F>
    using type = Maybe_error<Derivative<F, X>>;
};

template <class T, class X>
struct untransformed_type<Derivative<T, X>> {
    using type = T;
};

template <class T>
struct untransformed_type<Maybe_error<T>> {
    using type = T;
};

template <class T, class X>
struct untransformed_type<Maybe_error<Derivative<T, X>>> {
    using type = T;
};

template <class T, class X>
struct transformation_type<Derivative<T, X>> {
    using type = Derivative_Op<X>;
};


template <class T>
concept HasDx = requires(const T& t) {
    { t.dx() };  // optionally: -> /* some constraint */;
};

template <class T>
    requires HasDx<std::remove_reference_t<T>>
auto& get_dx(T const& x) {
    return x.dx();
}

template <class T>
    requires(!HasDx<std::remove_reference_t<T>>)
auto& get_dx(const T& x) {
    return x;
}
struct NoDerivative {
    template <class T>
    friend NoDerivative operator*(NoDerivative, T) {
        return NoDerivative{};
    }
    template <class T>
    friend NoDerivative operator*(T, NoDerivative) {
        return NoDerivative{};
    }

    template <class T>
    friend T operator+(T&& x, NoDerivative) {
        return std::forward<T>(x);
    }

    template <class T>
    friend T operator-(T&& x, NoDerivative) {
        return std::forward<T>(x);
    }

    template <class T>
    friend T operator-(NoDerivative, T&& x) {
        return -std::forward<T>(x);
    }

    template <class T>
    friend T operator+(NoDerivative, T&& x) {
        return std::forward<T>(x);
    }
    NoDerivative operator()() const { return {}; }
};

inline std::ostream& print(std::ostream& os, const NoDerivative& /*unused*/) {
    return os;
}

template <class>
struct is_derivative : public std::false_type {};

template <class X, class Y>
struct is_derivative<Derivative<X, Y>> : public std::true_type {};

template <class T>
constexpr bool is_derivative_v = is_derivative<std::decay_t<T>>::value;

template <class X, class Y>
decltype(auto) primitive(const Derivative<X, Y>& d) {
    return d.primitive();
}
template <class X, class Y>
decltype(auto) primitive(Derivative<X, Y>& d) {
    return d.primitive();
}

template <class X>
    requires(!is_derivative_v<X>)
decltype(auto) primitive(X&& x) {
    return std::forward<X>(x);
}

template <class X>
    requires(!is_derivative_v<X>)
auto derivative(X&&) {
    return NoDerivative{};
}

template <class T>
    requires(!is_derivative_v<T>)
bool has_dx(const T& /*x*/) {
    return false;
}

template <class T>
    requires(is_derivative_v<T>)
bool has_dx(const T& x) {
    return x.has_dx();
}

template <class X, class Y>
decltype(auto) derivative(const Derivative<X, Y>& d) {
    return d.derivative();
}


}  // namespace var
