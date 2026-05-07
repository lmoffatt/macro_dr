#pragma once

#include <concepts>
#include <utility>

// Algebraic structures live at two layers:
//
//   1. Concepts (signature-only) — constrain templates at compile time.
//      Names: Magma, Semigroup, Monoid, Group, ...
//
//   2. Law-bundle structs — stateful predicates carrying an Eq projection
//      and an operation F, callable as bool over a sample (or pack of
//      samples). Composed via Satisfies<...>.
//      Names: Associative_law, Has_Identity_law, Monoid_law, ...
//
// The two layers don't collide because they use different naming.
// A user typically constrains template code with the concept and runs
// runtime law checks via the law-bundle.

// ─── identity customization point ────────────────────────────────────────
// id_of<F, T>::value() returns the identity element of F under T. The default
// match below picks up any F with a static F::id() returning T. For carriers
// that don't expose F::id() — e.g. std::plus<int>, std::multiplies<double>,
// function pointers, lambdas — provide a partial specialization:
//
//   template <class T> struct id_of<std::plus<T>, T> {
//       static constexpr T value() { return T{0}; }
//   };

template <class F, class T>
struct id_of;

template <class F, class T>
    requires requires { { F::id() } -> std::same_as<T>; }
struct id_of<F, T> {
    static constexpr T value() { return F::id(); }
};

// ─── signature concepts ──────────────────────────────────────────────────

template <class F, class T>
concept Magma = requires(F f, T a, T b) {
    { f(a, b) } -> std::same_as<T>;
};

template <class F, class T>
concept Has_Identity = requires { { id_of<F, T>::value() } -> std::same_as<T>; };

template <class F, class T>
concept Has_Shaped_Identity = requires(T a) { { F::id(a) } -> std::same_as<T>; };

template <class F, class T>
concept Has_Inverse = requires(T a) { { F::inv(a) } -> std::same_as<T>; };

// ─── composite axiom-claim concepts ──────────────────────────────────────
// Refinement is by claim, not by extra signature. Concepts cannot verify
// these claims — that's the law-bundle / test layer below.

template <class F, class T>
concept Semigroup = Magma<F, T>;

template <class F, class T>
concept Monoid = Semigroup<F, T> && Has_Identity<F, T>;

template <class F, class T>
concept Shaped_Monoid = Semigroup<F, T> && Has_Shaped_Identity<F, T>;

template <class F, class T>
concept Commutative_Monoid = Monoid<F, T>;

template <class F, class T>
concept Group = Monoid<F, T> && Has_Inverse<F, T>;

template <class F, class T>
concept Abelian_Group = Group<F, T>;

// ─── law-bundle composition: Satisfies<P...> ─────────────────────────────
// AND-composition of stateful predicates. Inherits from each P so that the
// per-predicate state (eq, f) is stored as a base subobject; operator()
// forwards the same argument pack to each base and ANDs the results.
//
// All composed predicates must accept the same argument pack arity. If the
// arities differ, build separate Satisfies bundles or call the predicates
// individually.

template <class... P>
struct Satisfies : P... {
    Satisfies(P... p) : P(std::move(p))... {}

    template <class... T>
    constexpr bool operator()(T&&... t) const {
        return (... && static_cast<P const&>(*this)(std::forward<T>(t)...));
    }
};

// ─── per-axiom law bundles ───────────────────────────────────────────────

template <class Eq, class F, class T>
class Associative_law {
    Eq m_eq;
    F m_f;

   public:
    Associative_law(Eq eq, F f) : m_eq(eq), m_f(f) {}
    template <std::same_as<T>... Ts> 
    constexpr bool operator()(T const& a, T const& b, T const& c,Ts const&... d) const 
    {
        return m_eq(m_f(m_f(a, b), c), m_f(a, m_f(b, c))) &&
        (*this)(b,c,d...);
    }
    constexpr bool operator()(T const& , T const& ) const {
        return true;
    }
};

template <class Eq, class F, class T>
class Has_Identity_law {
    Eq m_eq;
    F m_f;

   public:
    Has_Identity_law(Eq eq, F f) : m_eq(eq), m_f(f) {}

    template <std::same_as<T>... Ts>
    constexpr bool operator()(Ts... a) const {
        return (...
                && (m_eq(m_f(id_of<F, T>::value(), a), a)
                    && m_eq(m_f(a, id_of<F, T>::value()), a)));
    }
};

template <class Eq, class F, class T>
    requires Monoid<F, T>
struct Monoid_law
    : Satisfies<Associative_law<Eq, F, T>, Has_Identity_law<Eq, F, T>> {
    using base_type =
        Satisfies<Associative_law<Eq, F, T>, Has_Identity_law<Eq, F, T>>;

    Monoid_law(Eq eq, F f)
        : base_type{Associative_law<Eq, F, T>(eq, f),
                    Has_Identity_law<Eq, F, T>(eq, f)} {}

    using base_type::operator();
};



// Stepanov fast exponentiation (Elements of Programming §7.1):
// O(log n) operations exploiting associativity, with a leading-zero-bit
// skip in the first loop that saves one squaring vs. naive binary-exp.

template <class F, class T>
    requires Semigroup<F, T>
T power_semigroup(F f, T a, std::size_t n) {
    // Precondition: n >= 1.
    while ((n & 1U) == 0U) {
        a = f(a, a);
        n >>= 1;
    }
    T r = a;
    n >>= 1;
    while (n != 0U) {
        a = f(a, a);
        if ((n & 1U) != 0U) {
            r = f(r, a);
        }
        n >>= 1;
    }
    return r;
}

template <class F, class T>
    requires Monoid<F, T>
T monoid_power(F f, T a, std::size_t n) {
    if (n == 0U) {
        return id_of<F, T>::value();
    }
    return power_semigroup(f, std::move(a), n);
}

