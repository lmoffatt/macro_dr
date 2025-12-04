#ifndef VARIABLES_DERIVATIVE_H
#define VARIABLES_DERIVATIVE_H
#include "derivative_operator.h"
#include "variables.h"

namespace var {

template <class Id, class T, class X>
class Derivative<Var<Id, T>, X> {
    Derivative<T, X> m_x;

   public:
    static constexpr bool is_variable = true;
    template <class S>
        requires(std::constructible_from<Derivative<T, X>, S>)
    constexpr Derivative(S&& t_x) : m_x{std::forward<S>(t_x)} {}
    constexpr auto& operator()() { return m_x; }
    constexpr auto& operator()() const { return m_x; }
    //constexpr auto& operator[](Var<Id>) const { return *this; }

    template <class... Ts>
    constexpr auto operator()(const Ts&...) const {
        return Derivative<Id, X>(*this);
    }

    constexpr Derivative() = default;
    constexpr auto& value() const { return m_x; }
    decltype(auto) primitive() const { return m_x.primitive(); }
    
    auto derivative() const 
    requires requires (const Derivative<T, X>& d) { d.derivative(); }
    { return m_x.derivative(); }

    auto& dx() const { return m_x.dx(); }
    bool has_dx() const { return m_x.has_dx(); }

    friend auto& print(std::ostream& os, const Derivative& x) {
        os << type_name<Id>() << ": \n";
        print(os, x.value());
        os << "\n";
        return os;
    }

    // friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
    // friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}
};

template <class Id, class F, class... T, class X>
class Derivative<Fun<Id, F, T...>, X> {
    std::tuple<Derivative<T, X>...> m_x;
    F m_f;

   public:
    static constexpr bool is_variable = true;
    template <class... S>
        requires((std::constructible_from<Derivative<T, X>, S>) && ...)
    constexpr Derivative(Var<Id>, F const& t_f, S&&... t_x)
        : m_x{std::forward<S>(t_x)...}, m_f{t_f} {}
    constexpr auto& operator()() { return m_x; }
    constexpr auto& operator()() const { return m_x; }
    constexpr auto& operator[](Var<Id>) const { return *this; }

    template <class... Ts>
    constexpr auto operator()(const Ts&... ts) {
        return Derivative<Id, X>(std::apply(
            [this, &ts...](auto&... xs) { return std::invoke(m_f, xs..., ts...); }, m_x));
    }

    template <class... Ts>
    constexpr auto operator()(const Ts&... ts) const {
        return Derivative<Id, X>(std::apply(
            [this, &ts...](auto&... xs) { return std::invoke(m_f, xs..., ts...); }, m_x));
    }

    constexpr Derivative() = default;

    decltype(auto) dx() const {
        return std::apply([](auto&... ds) -> decltype(auto) { return get_dx_of_dfdx(ds...); }, m_x);
    }
    bool has_dx() const { return (std::get<decltype(T())>(m_x).has_dx() || ... || false); }

    friend auto& print(std::ostream& os, const Derivative& x) {
        os << type_name<Id>() << ": \n";
        print(os, x.m_x);
        os << "\t";
        return os;
    }

    // friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
    // friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}
};

template <class... Ids, class X>
    requires(Ids::is_variable && ...)
class Derivative<Vector_Space<Ids...>, X> : public Vector_Space<Derivative_t<Ids, X>...> {
    using base_type = Vector_Space<Derivative_t<Ids, X>...>;

   public:
    template <class Id>
    friend auto& get(Derivative const& x) {
        return static_cast<Derivative_t<Id, X> const&>(x);
    }
    template <class Id>
    friend auto& get(Derivative& x) {
        return static_cast<Derivative_t<Id, X>&>(x);
    }

    template <class... SelIds>
    friend auto select(Derivative&& x) {
        return Derivative<Vector_Space<SelIds...>, X>(std::move(get<SelIds>(x))...);
    }


    static constexpr bool is_vector_space = true;

    Derivative() = default;
    Derivative(Derivative_t<Ids, X>&&... t_vars) : base_type{std::move(t_vars)...} {}
    Derivative(Derivative_t<Ids, X> const&... t_vars) : base_type{t_vars...} {}
    // Vector_Space(std::decay_t <decltype(std::declval<Vars const&>().value())> ... t_vars): Vars{std::move(t_vars)}...{}

    friend auto& operator<<(std::ostream& os, const Derivative& x) {
        (os << ... << get<Ids>(x));
        return os;
    }
    // friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}

    Vector_Space<Ids...> primitive() const {
        return Vector_Space<Ids...>(var::primitive(get<Ids>(*this))...);
    }

    auto dx() const { return get_dx_of_dfdx(get<Ids>(*this)...); }
    bool has_dx() const { return (var::has_dx(get<Ids>(*this)) || ... || false); }
};

template <class... Ids, class X>
auto const& get_dx_of_dfdx(const Derivative<Vector_Space<Ids...>, X>& f) {
    return get_dx_of_dfdx(get<Ids>(f)...);
}

template <class... Ids, class X, class Parameters>
auto Taylor_first(const Derivative<Vector_Space<Ids...>, X>& f, const Parameters& x, double eps) {
    return Vector_Space<Ids...>(Ids(Taylor_first(get<Ids>(f), x, eps).value())...);
}

// initialize overload for Vector_Space values:
// when the input lives in a derivative universe, lift each component Id into
// a Derivative_t<Id,X> that shares the same dx() pointer. This complements
// the scalar initialize(Input,T) in derivative_operator.h.
template <class Input, class... Ids>
    requires(var::is_derivative_v<std::decay_t<Input>>)
auto initialize(Vector_Space<Ids...>&& v, const Input& in) {
    using Dir = dx_of_dfdx_t<std::decay_t<Input>>;
    static_assert(!std::is_same_v<Dir, NoDerivative>,
                  "initialize(Vector_Space,Input): derivative input has no dx_of_dfdx_t direction");

    return Derivative<Vector_Space<Ids...>, Dir>(initialize(get<Ids>(std::move(v)), in)...);
}

template <class Id, class X>
    requires(Id::is_variable)
class Derivative<Id, X> : public Derivative<typename Id::variable_type, X> {
   public:
    using base_type = Derivative<typename Id::variable_type, X>;

    template <class IdT>
        requires std::is_same_v<Id, std::decay_t<IdT>>
    Derivative(IdT&& m) : base_type{std::forward<IdT>(m)()} {}
    decltype(auto) primitive() const { return Id(base_type::primitive()); }

    constexpr auto const& operator[](Var<Id>) const { return *this; }
    constexpr auto & operator[](Var<Id>)  { return *this; }

    Derivative(base_type&& m) : base_type{std::move(m)} {}
    Derivative(base_type const& m) : base_type{m} {}
    Derivative() = default;

    auto& dx() const { return base_type::dx(); }
    bool has_dx() const { return base_type::has_dx(); }
};

template <class Id, class T, class X>
auto get_dx_of_dfdx(const Derivative<Var<Id, T>, X>& f) {
    return get_dx_of_dfdx(f());
}

template <class Id, class T, class X>
Var<Id, T> Taylor_first(const Derivative<Var<Id, T>, X>& f, const X& x, double eps) {
    return Taylor_first(f(), x, eps);
}

template <class Id, class T, class X>
Var<Id, T> Taylor_first(const Var<Id, T>& f, const X& x, double eps) {
    return f;
}

template <class Id, class T, class X>
Constant<Id, T> Taylor_first(const Constant<Id, T>& f, const X& x, double eps) {
    return f;
}

}  // namespace var

#endif  // VARIABLES_DERIVATIVE_H
