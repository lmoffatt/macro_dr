#ifndef PARAMETERS_DERIVATIVE_H
#define PARAMETERS_DERIVATIVE_H

#include <cassert>
#include <cmath>

// Toggleable assertion for derivative dx-pointer checks.
// Define MACRODR_STRICT_DX_ASSERT at compile time to enable hard asserts.
#ifndef MACRODR_DX_ASSERT
  #ifndef MACRODR_STRICT_DX_ASSERT
    #define MACRODR_DX_ASSERT(cond) ((void)0)
  #else
    #define MACRODR_DX_ASSERT(cond) assert(cond)
  #endif
#endif

#include <cstddef>
#include <limits>
#include <vector>

#include "derivative_operator.h"
#include "matrix.h"
#include "matrix_der_traits.h"
#include "maybe_error.h"
#include "parameters.h"
//#include "general_algorithm_on_containers.h"
// #include "matrix_derivative.h"

namespace var {

template <class Matrix, class T>
    requires std::is_same_v<element_of_Matrix_t<Matrix>, T>
auto build_(std::size_t nrows, std::size_t ncols,
            std::initializer_list<std::pair<std::size_t, std::size_t>> indexes,
            std::initializer_list<T> values) {
    Matrix x(nrows, ncols, true);
    auto it_ind = indexes.begin();
    auto it_val = values.begin();
    for (std::size_t k = 0; k < values.size(); ++k) {
        auto i = it_ind->first;
        auto j = it_ind->second;
        auto value = *it_val;
        x(i, j) = value;
        ++it_ind;
        ++it_val;
    }
    return x;
}

template <class Matrix, class Der>
    requires(is_derivative_v<Der>)
auto build_(std::size_t nrows, std::size_t ncols,
            std::initializer_list<std::pair<std::size_t, std::size_t>> indexes,
            std::initializer_list<Der> values) {
    using X = dx_of_dfdx_t<Der>;
    decltype(auto) dx = get_dx_of_dfdx(*values.begin());
    using T = element_of_Matrix_t<Matrix>;

    auto df_dX0 = derivative(*values.begin());
    auto n = values.size();
    auto dfdx = applyMap_i(
        [&values, &indexes, n, nrows, ncols](auto e, std::size_t ii) {
            auto it_ind = indexes.begin();
            auto it_val = values.begin();
            Matrix ddx(nrows, ncols, true);
            for (std::size_t k = 0; k < n; ++k) {
                auto i = it_ind->first;
                auto j = it_ind->second;
                auto& Dx = *it_val;
                ddx(i, j) = derivative(Dx)()[ii];
                ++it_ind;
                ++it_val;
            }
            return ddx;
        },
        df_dX0());
    Matrix x(nrows, ncols, true);
    auto it_ind = indexes.begin();
    auto it_val = values.begin();

    for (std::size_t k = 0; k < values.size(); ++k) {
        auto i = it_ind->first;
        auto j = it_ind->second;
        auto& Dx = *it_val;
        x(i, j) = primitive(Dx);
        ++it_ind;
        ++it_val;
    }

    return Derivative<Matrix, X>(x, dfdx, dx);
}

template <>
class d_d_<double, Parameters_transformed> {
    Matrix<double> m_dydx;
    Parameters_transformed const* ptr_dx=nullptr;

   public:
    using value_type = Matrix<double>;
    template <class aMatrix>
        requires std::is_same_v<Matrix<double>, std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix&& dydx, Parameters_transformed const& x)
        : m_dydx{std::forward<aMatrix>(dydx)}, ptr_dx{&x} {}

    template <class aMatrix>
        requires std::is_same_v<Matrix<double>, std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix&& dydx, Parameters_transformed const* x)
        : m_dydx{std::forward<aMatrix>(dydx)}, ptr_dx{x} {}
    d_d_() = default;

    auto& operator()() { return m_dydx; }
    auto& operator()() const { return m_dydx; }

    auto& dx() const {
        MACRODR_DX_ASSERT(ptr_dx != nullptr &&
                          "Derivative d_d_<double, Parameters_transformed> missing dx pointer");
        return *ptr_dx;
    }
    bool has_dx() const { return ptr_dx != nullptr; }
    void set_dx(Parameters_transformed const& x) { ptr_dx = &x; }

    friend auto operator*(d_d_ const& df, const Parameters_transformed& x) {
        double out = 0;
        for (std::size_t i = 0; i < df().size(); ++i) out += df()[i] * x()[i];
        return out;
    }
};

template <template <class> class aMatrix>
    requires(aMatrix<double>::is_Matrix)
class d_d_<aMatrix<double>, Parameters_transformed> {
    Matrix<aMatrix<double>> m_dydx;
    Parameters_transformed const* ptr_par=nullptr;

   public:
    using value_type = Matrix<aMatrix<double>>;
    constexpr d_d_() = default;
    [[nodiscard]] d_d_(const d_d_&) = default;
    d_d_(d_d_&&) = default;
    d_d_& operator=(const d_d_&) = default;

    constexpr ~d_d_() = default;
    template <class aaMatrix>
        requires std::constructible_from<Matrix<aMatrix<double>>, std::decay_t<aaMatrix>>
    constexpr d_d_(aaMatrix&& dydx, const Parameters_transformed& x)
        : m_dydx{std::forward<aaMatrix>(dydx)}, ptr_par{&x} {}

    constexpr d_d_(aMatrix<double> const& y, const Parameters_transformed& x)
        : m_dydx{x().nrows(), x().ncols(), aMatrix<double>(y.nrows(), y.ncols(), 0.0)},
          ptr_par{&x} {}

    constexpr auto& operator()() { return m_dydx; }
    constexpr auto& operator()() const { return m_dydx; }

    auto& dx() const {
        MACRODR_DX_ASSERT(ptr_par != nullptr &&
                          "Derivative d_d_<Matrix, Parameters_transformed> missing dx pointer");
        return *ptr_par;
    }
    bool has_dx() const { return ptr_par != nullptr; }
    void set_dx(Parameters_transformed const& x) { ptr_par = &x; }

    template <class F>
    friend auto apply_par(F&& f, d_d_ const& a) {
        using S = std::decay_t<std::invoke_result_t<F, aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a()[i]);
        MACRODR_DX_ASSERT(a.has_dx() &&
                          "apply_par: d_d_<Matrix<double>,Parameters_transformed> missing dx");
        return d_d_<S, Parameters_transformed>(x, a.dx());
    }

    template <class F, template <class> class bMatrix>
    friend auto zip_par(F&& f, d_d_ const& a,
                        d_d_<bMatrix<double>, Parameters_transformed> const& b) {
        using S = std::decay_t<std::invoke_result_t<F, aMatrix<double>, bMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());

        for (std::size_t i = 0; i < std::min(x.size(), b().size()); ++i) x[i] = f(a()[i], b()[i]);

        if (a.has_dx() && a.dx().size() > 0) {
            return d_d_<S, Parameters_transformed>(x, a.dx());
        }
        return d_d_<S, Parameters_transformed>(x, b.dx());
    }

    template <class F>
    friend auto zip_par(F&& f, d_d_ const& a, d_d_<double, Parameters_transformed> const& b) {
        using S = std::decay_t<std::invoke_result_t<F, aMatrix<double>, double>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < std::min(x.size(), b().size()); ++i) x[i] = f(a()[i], b()[i]);
        if (a.has_dx() && a.dx().size() > 0) {
            return d_d_<S, Parameters_transformed>(x, a.dx());
        }
        return d_d_<S, Parameters_transformed>(x, b.dx());
    }

    template <class F>
    friend auto zip_par(F&& f, d_d_<double, Parameters_transformed> const& b, d_d_ const& a) {
        using S = std::decay_t<std::invoke_result_t<F, double, aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(b()[i], a()[i]);
        if (a.has_dx() && a.dx().size() > 0) {
            return d_d_<S, Parameters_transformed>(x, a.dx());
        }
        return d_d_<S, Parameters_transformed>(x, b.dx());
    }

    friend auto operator*(d_d_ const& df, const Parameters_transformed& x) {
        if (df().size() == 0)
            return aMatrix<double>{};
        else {
            auto nrows = df()[0].nrows();
            auto ncols = df()[0].ncols();

            auto out = aMatrix<double>(nrows, ncols);

            for (std::size_t i = 0; i < std::min(df().size(), x.size()); ++i)
                out = out + df()[i] * x()[i];
            return out;
        }
    }
};

template <template <class> class notSymmetricMatrix>
    requires((notSymmetricMatrix<double>::is_Matrix) &&
             (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
                std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto inside_out(const d_d_<notSymmetricMatrix<double>, Parameters_transformed>& x) {
    if (x().size() == 0)
        return notSymmetricMatrix<d_d_<double, Parameters_transformed>>{};
    else {
        notSymmetricMatrix<d_d_<double, Parameters_transformed>> out(x()[0].nrows(),
                                                                     x()[0].ncols());
        for (std::size_t i = 0; i < out.size(); ++i) {
            Matrix<double> d_d_par(x().nrows(), x().ncols());
            for (std::size_t j = 0; j < d_d_par.size(); ++j) d_d_par[j] = x()[j][i];
            out[i] = d_d_<double, Parameters_transformed>(std::move(d_d_par), x.dx());
        }
        return out;
    }
}

template <template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto inside_out(const d_d_<aSymmetricMatrix<double>, Parameters_transformed>& x) {
    aSymmetricMatrix<d_d_<double, Parameters_transformed>> out(x()[0].nrows(), x()[0].ncols());
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            Matrix<double> d_d_par(x().nrows(), x().ncols());
            for (std::size_t k = 0; k < d_d_par.size(); ++k) d_d_par[k] = x()[k](i, j);
            out.set(i, j, d_d_<double, Parameters_transformed>(std::move(d_d_par), x.dx()));
        }
    return out;
}

template <>
class d_d_<Parameters_transformed, Parameters_transformed> {
    d_d_<Matrix<double>, Parameters_transformed> m_dydx;

   public:
    using value_type = d_d_<Matrix<double>, Parameters_transformed>;
    d_d_() = default;
    [[nodiscard]] d_d_(const d_d_&) = default;
    d_d_(d_d_&&) noexcept = default;
    d_d_& operator=(const d_d_&) = default;
    auto& dx() const { return m_dydx.dx(); }
    bool has_dx() const { return m_dydx.has_dx(); }
    void set_dx(Parameters_transformed const& x) { m_dydx.set_dx(x); }
    template <class aMatrix>
        requires std::constructible_from<d_d_<Matrix<double>, Parameters_transformed>, aMatrix>
    constexpr d_d_(aMatrix&& dydx) : m_dydx{std::forward<aMatrix>(dydx)} {}

    constexpr auto& operator()() { return m_dydx; }
    constexpr auto& operator()() const { return m_dydx; }
};

template <>
class Derivative<double, Parameters_transformed>  //: public Primitive<double>, public
//: d_d_<double,Parameters_transformed>
{
    using primitive_type = double;
    using d_type = Parameters_transformed;
    using derivative_type = d_d_<primitive_type, d_type>;
    primitive_type m_x;
    derivative_type m_d;

   public:
    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D>)

    Derivative(P t_x, D&& t_d) : m_x{t_x}, m_d{std::forward<D>(t_d)} {}

    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D, Parameters_transformed const&>)

    Derivative(P t_x, D&& t_d, Parameters_transformed const& x)
        : m_x{t_x}, m_d{std::forward<D>(t_d), x} {}
    Derivative() = default;

    Derivative(double t_x, Parameters_transformed const& x)
        : m_x{t_x}, m_d{Matrix<double>(x().nrows(), x().ncols(), 0.0), x} {}

    Derivative& operator=(const Derivative<double, Matrix<double>>& x) {
        m_x = x.primitive();
        m_d() = x.derivative()();
        return *this;
    }

    // Value-only constructor: used in a few generic
    // algorithms (e.g. clamping helpers) where a
    // derivative shell is constructed from a scalar.
    // These shells deliberately have no dx(); any
    // attempt to use them in a semantic derivative
    // context will trip MACRODR_DX_ASSERT when
    // dx() or get_dx_of_dfdx(...) is called.
    Derivative(double t_x) : m_x{t_x} {}

    Derivative(const Derivative& x) = default;
    Derivative(Derivative&& x) = default;
    Derivative& operator=(const Derivative& x) = default;
    Derivative& operator=(Derivative&& x) = default;

    auto& primitive() { return m_x; }
    auto& primitive() const { return m_x; }
    auto& derivative() const { return m_d; }
    auto& derivative() { return m_d; }
    auto& dx() const { return m_d.dx(); }
    bool has_dx() const { return m_d.has_dx(); }

    Derivative<double, Matrix<double>> operator()() const {
        return Derivative<double, Matrix<double>>(primitive(), derivative()(), dx()());
    }

    friend auto operator/(const Derivative& x, const Derivative& y) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(
                fx / fy,
                zip([fx, fy](auto dx, auto dy) { return dx / fy - fx / fy / fy * dy; },
                    x.derivative()(), y.derivative()()),
                x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(fx / fy, -fx / fy / fy * y.derivative()(), y.dx());
        else
            return Derivative(fx / fy, x.derivative()() / fy, x.dx());
    }

    friend auto operator/(double x, const Derivative& y) {
        auto fy = y.primitive();
        return Derivative(x / fy, -x / fy / fy * y.derivative()(), y.dx());
    }

    friend auto elemDivSafe(const Derivative& x, const Derivative& y, double eps) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(elemDivSafe(fx, fy, eps),
                              zip(
                                  [fx, fy, eps](auto dx, auto dy) {
                                      return elemDivSafe(1.0, fy, eps) * dx -
                                             elemDivSafe(fx, fy * fy, eps) * dy;
                                  },
                                  x.derivative()(), y.derivative()()),
                              x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(elemDivSafe(fx, fy, eps),
                              -elemDivSafe(fx, fy * fy, eps) * y.derivative()(), x.dx());
        else
            return Derivative(elemDivSafe(fx, fy, eps), elemDivSafe(x.derivative()(), fy, eps),
                              x.dx());
    }

    friend auto elemDivSoftAbs(const Derivative& x, const Derivative& y, double eps) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        auto denom = std::sqrt(fy * fy + eps * eps);
        if (denom == 0.0)
            denom = std::numeric_limits<double>::epsilon();
        auto inv_denom = 1.0 / denom;
        auto inv_denom3 = inv_denom * inv_denom * inv_denom;
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(
                elemDivSoftAbs(fx, fy, eps),
                zip([inv_denom, inv_denom3, fx, fy](
                        auto dx, auto dy) { return dx * inv_denom - fx * fy * inv_denom3 * dy; },
                    x.derivative()(), y.derivative()()),
                x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(elemDivSoftAbs(fx, fy, eps), -fx * fy * inv_denom3 * y.derivative()(),
                              x.dx());
        else
            return Derivative(elemDivSoftAbs(fx, fy, eps), x.derivative()() * inv_denom, x.dx());
    }

    friend auto operator*(const Derivative& x, const Derivative& y) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(fx * fy,
                              zip([fx, fy](auto dx, auto dy) { return dx * fy + fx * dy; },
                                  x.derivative()(), y.derivative()()),
                              x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(fx * fy, fx * y.derivative()(), x.dx());
        else
            return Derivative(fx * fy, x.derivative()() * fy, x.dx());
    }

    friend auto exp(const Derivative& x) {
        auto f = exp(x.primitive());
        return Derivative(f, f * x.derivative()(), x.dx());
    }

    friend auto max(const Derivative& x, double y) {
        if (x.primitive() >= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }

    friend auto min(const Derivative& x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }

    friend bool operator>(const Derivative& x, const Derivative& y) {
        return x.primitive() > y.primitive();
    }

    friend bool operator<(const Derivative& x, const Derivative& y) {
        return x.primitive() < y.primitive();
    }

    friend auto max(const Derivative& x, const Derivative& y) {
        if (x.primitive() >= y.primitive())
            return x;
        else
            return y;
    }

    friend auto min(const Derivative& x, const Derivative& y) {
        if (x.primitive() <= y.primitive())
            return x;
        else
            return y;
    }

    friend auto log(const Derivative& x) {
        auto f = log(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / x.primitive()), x.dx());
    }
    friend auto log10(const Derivative& x) {
        auto f = log10(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / (x.primitive() * std::log(10))), x.dx());
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

    friend auto sqrt(const Derivative& x) {
        auto p = x.primitive();
        auto root = std::sqrt(p);
        if (root == 0.0)
            return Derivative(root, 0.0 * x.derivative()(), x.dx());
        auto scale = 0.5 / root;
        return Derivative(root, x.derivative()() * scale, x.dx());
    }

    friend bool operator==(Derivative const& one, double val) { return one.primitive() == val; }
};

template <template <class> class T_Matrix>
    requires(T_Matrix<double>::is_Matrix)
class Derivative<T_Matrix<double>,
                 Parameters_transformed> {  //: public T_Matrix<double>{
    using primitive_type = T_Matrix<double>;
    using d_type = Parameters_transformed;
    using derivative_type = d_d_<M_der_t<primitive_type>, d_type>;
    primitive_type m_x;
    derivative_type m_d;

   public:
    auto ncols() const { return m_x.ncols(); }  //{return primitive_type::ncols();}
    auto nrows() const { return m_x.nrows(); }  // {return primitive_type::nrows();}
    auto size() const { return m_x.size(); }    //{return primitive_type::size();}
    Derivative() = default;

    //   using gserg=typename primitive_type::sgr;
    //   using gserug=typename derivative_type::sgrdd;

    Derivative(std::size_t nrows, std::size_t ncols) : m_x(nrows, ncols) {}
    Derivative(std::size_t nrows, std::size_t ncols, double v) : m_x(nrows, ncols, v) {}
    auto& dx() const { return m_d.dx(); }
    [[nodiscard]] bool has_dx() const { return m_d.has_dx(); }
    void set_dx(Parameters_transformed const& x) { m_d.set_dx(x); }

    template <class F>
    friend auto apply(F&& f, const Derivative& x) {
        return outside_in(apply(std::forward<F>(f), inside_out(x)), x.dx());
    }

    template <class F>
    friend auto zip(F&& f, const Derivative& x, const Derivative& y) {
        if (x.has_dx() && x.dx().size() > 0) {
            return outside_in(zip(std::forward<F>(f), inside_out(x), inside_out(y)), x.dx());
        }
        return outside_in(zip(std::forward<F>(f), inside_out(x), inside_out(y)), y.dx());
    }

    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D>)
    Derivative(P&& t_x, D&& t_d) : m_x{std::forward<P>(t_x)}, m_d{std::forward<D>(t_d)} {}

    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D, const Parameters_transformed&>)
    Derivative(P&& t_x, D&& t_d, const Parameters_transformed& x)
        : m_x{std::forward<P>(t_x)}, m_d{std::forward<D>(t_d), x} {}

    template <class anotherCompatibleMatrix>
    Derivative(Derivative<anotherCompatibleMatrix, d_type> const& t_x)
        : m_x{t_x.primitive()}, m_d{t_x.derivative()(), t_x.dx()} {}

    template <class P>
        requires(std::constructible_from<primitive_type, P>)
    // Value-only constructor: provides an empty-shell derivative
    // with no dx() for container / structural uses. Semantic paths
    // should prefer constructors that also set dx(), and misuse is
    // guarded at runtime via MACRODR_DX_ASSERT in dx accessors.
    Derivative(P&& t_x) : m_x{std::forward<P>(t_x)}, m_d{} {}

    auto& primitive_non_const() { return m_x; }  // {return static_cast<primitive_type&>(*this);}
    auto& primitive() { return m_x; }            // {return static_cast<primitive_type&>(*this);}
    auto& primitive() const { return m_x; }  //{return static_cast<primitive_type const&>(*this);}
    auto& derivative() { return m_d; }
    auto& derivative_non_const() { return m_d; }
    auto& derivative() const { return m_d; }

    auto& dx() { return derivative().dx(); }

    friend auto operator*(const Derivative& x, double y) {
        return Derivative(x.primitive() * y, x.derivative()() * y, x.dx());
    }

    template <template <class> class anotherCompatibleMatrix>
    friend auto operator*(const Derivative& x,
                          anotherCompatibleMatrix<Derivative<double, d_type>> const& t_y) {
        if (x.has_dx() && x.dx().size() > 0) {
            auto Y = outside_in(t_y, x.dx());
            return x * Y;
        }
        auto Y = outside_in(t_y, get_dx_of_dfdx_container(t_y));
        return x * Y;
    }

    template <template <class> class anotherCompatibleMatrix>
    friend auto operator*(anotherCompatibleMatrix<Derivative<double, d_type>> const& t_y,
                          const Derivative& x) {
        if (x.has_dx() && x.dx().size() > 0) {
            auto Y = outside_in(t_y, x.dx());
            return Y * x;
        }
        auto Y = outside_in(t_y, get_dx_of_dfdx_container(t_y));
        return Y * x;
    }

    auto operator()(std::size_t i, std::size_t j) const {
        return Derivative<double, Parameters_transformed>(
            primitive()(i, j), applyMap([i, j](auto const& m) { return m(i, j); }, derivative()()),
            dx());
    }

    void set(std::size_t i, std::size_t j, double x) {
        primitive_non_const()(i, j) = x;
        for (std::size_t ij = 0; ij < derivative()().size(); ++ij)
            derivative_non_const()()[ij](i, j) = 0;
    }
    void set(std::size_t i, double x) {
        primitive_non_const()[i] = x;
        for (std::size_t ij = 0; ij < derivative()().size(); ++ij)
            derivative_non_const()()[ij][i] = 0;
    }

    void set(std::size_t i, Derivative<double, Parameters_transformed> const& x) {
        primitive_non_const()[i] = x.primitive();
        derivative_non_const().set_dx(x.dx());
        for (std::size_t ij = 0; ij < derivative()().size(); ++ij)
            derivative_non_const()()[ij][i] = x.derivative()()[ij];
    }
    void set(std::size_t i, std::size_t j, Derivative<double, Parameters_transformed> const& x) {
        primitive_non_const()(i, j) = x.primitive();
        derivative_non_const().set_dx(x.dx());
        for (std::size_t ij = 0; ij < derivative()().size(); ++ij)
            derivative_non_const()()[ij](i, j) = x.derivative()()[ij];
    }

    auto operator[](std::size_t i) const {
        return Derivative<double, Parameters_transformed>(
            primitive()[i], applyMap([i](auto const& m) { return m[i]; }, derivative()()), dx());
    }

    auto operator[](std::pair<std::size_t, std::size_t> ij) const {
        return Derivative<Matrix<double>, Parameters_transformed>(
            primitive()[ij], applyMap([ij](auto const& m) { return m[ij]; }, derivative()()), dx());
    }
};

template <template <class> class T_Matrix>
    requires(T_Matrix<double>::is_Matrix)
double fullsum(const Derivative<T_Matrix<double>, Parameters_transformed>& x) {
    return fullsum(x.primitive()) + fullsum(x.derivative()());
}

#if false
template <class T>
class Derivative<std::vector<T>, Parameters_transformed>
    : public std::vector<Derivative<T, Parameters_transformed>> {  //: public Container<double>{

   public:
    using base_type = std::vector<Derivative<T, Parameters_transformed>>;
    using base_type::size;
    using base_type::vector;
    using base_type::operator[];
    using base_type::push_back;

    auto primitive() const {
        std::vector<T> out(size());
        for (std::size_t i = 0; i < size(); ++i) out[i] = (*this)[i].primitive();
        return out;
    }

    auto derivative() const {
        d_d_<std::vector<T>, Parameters_transformed> out(Matrix<double>(size(), 1, 0.0),
                                                          (*this)[0].dx());
        for (std::size_t i = 0; i < size(); ++i) out[i] = (*this)[i].derivative();
        return out;
    }

    auto& dx() { return (*this)[0].derivative().dx(); }
};
#endif

// Minimal vector-of-derivatives wrapper to enable derivative containers of std::vector
// without introducing a full d_d_<std::vector<...>> algebra.
template <class T>
class Derivative<std::vector<T>, Parameters_transformed>
    : public std::vector<Derivative<T, Parameters_transformed>> {
   public:
    using base_type = std::vector<Derivative<T, Parameters_transformed>>;
    using base_type::operator[];
    using base_type::push_back;
    using base_type::size;
    using base_type::vector;

    Derivative() = default;
    Derivative(const Derivative&) = default;
    Derivative(Derivative&&) = default;
    Derivative& operator=(const Derivative&) = default;
    Derivative& operator=(Derivative&&) = default;

    auto primitive() const {
        std::vector<T> out(this->size());
        for (std::size_t i = 0; i < this->size(); ++i) out[i] = (*this)[i].primitive();
        return out;
    }

    // Minimal JSON-friendly derivative payload (empty matrix)
    auto derivative() const { return Matrix<double>(); }

    auto& dx() const { return (*this)[0].derivative().dx(); }
};

template <class F>
    requires requires(Derivative<F, Parameters_transformed>& f) {
        { f.derivative()().nrows() };
    }
auto& get_dx_of_dfdx(const Derivative<F, Parameters_transformed>& f) {
    MACRODR_DX_ASSERT(f.has_dx() &&
                      "get_dx_of_dfdx: Derivative<F, Parameters_transformed> missing dx pointer");
    return f.dx();
}

template <class F, template <class> class Container>
Parameters_transformed const& get_dx_of_dfdx_container(
    const Container<Derivative<F, Parameters_transformed>>& f) {
    assert(f.size() > 0);
    for (std::size_t i = 0; i < f.size(); ++i) {
        if (f[i].has_dx() && f[i].dx().size() > 0) {
            return f[i].dx();
        }
    }
    assert(false);
    return f[0].dx();
}

template <template <class> class notSymmetricMatrix>
    requires((notSymmetricMatrix<double>::is_Matrix) &&
             (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
                std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto inside_out(const Derivative<notSymmetricMatrix<double>, Parameters_transformed>& x) {
    notSymmetricMatrix<Derivative<double, Parameters_transformed>> out(x.nrows(), x.ncols());

    auto der = inside_out(x.derivative());
    if (der.size() > 0) {
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = Derivative<double, Parameters_transformed>(x.primitive()[i], der[i]);
        }
    } else {
        // No derivative payload available yet; still carry the dx pointer
        // so downstream code can access parameter context safely.
        auto zero = Matrix<double>(x.derivative()().nrows(), x.derivative()().ncols(), 0.0);
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = Derivative<double, Parameters_transformed>(
                x.primitive()[i], d_d_<double, Parameters_transformed>(zero, x.dx()));
        }
    }
    return out;
}

template <template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto inside_out(const Derivative<aSymmetricMatrix<double>, Parameters_transformed>& x) {
    aSymmetricMatrix<Derivative<double, Parameters_transformed>> out(x.nrows(), x.ncols());
    auto der = inside_out(x.derivative());
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j)
            out.set(i, j,
                    Derivative<double, Parameters_transformed>(x.primitive()(i, j), der(i, j)));
    return out;
}

template <template <class> class Matrix>
auto& outside_in(const Matrix<double>& x, ...) {
    return x;
}

template <template <class> class Matrix>
auto& inside_out(const Matrix<double>& x) {
    return x;
}

template <template <class> class notSymmetricMatrix>
    requires((notSymmetricMatrix<double>::is_Matrix) &&
             (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
                std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto outside_in(const notSymmetricMatrix<Derivative<double, Parameters_transformed>>& x,
                const Parameters_transformed& dx) {
    auto prim = notSymmetricMatrix<double>(x.nrows(), x.ncols());
    for (std::size_t i = 0; i < prim.size(); ++i) prim[i] = x[i].primitive();

    auto der =
        Matrix<notSymmetricMatrix<double>>(x[0].derivative()().nrows(), x[0].derivative()().ncols(),
                                           notSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t i = 0; i < der.size(); ++i)
        for (std::size_t j = 0; j < der[i].size(); ++j) der[i][j] = x[j].derivative()()[i];

    return Derivative<notSymmetricMatrix<double>, Parameters_transformed>(std::move(prim),
                                                                          std::move(der), dx);
}

template <template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto outside_in(const aSymmetricMatrix<Derivative<double, Parameters_transformed>>& x,
                const Parameters_transformed& dx) {
    auto prim = aSymmetricMatrix<double>(x.nrows(), x.ncols());
    for (std::size_t i = 0; i < prim.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j) prim.set(i, j, x(i, j).primitive());

    auto der =
        Matrix<aSymmetricMatrix<double>>(x[0].derivative()().nrows(), x[0].derivative()().ncols(),
                                         aSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t i = 0; i < der.size(); ++i)
        for (std::size_t j1 = 0; j1 < der[i].nrows(); ++j1)
            for (std::size_t j2 = 0; j2 <= j1; ++j2)
                der[i].set(j1, j2, x(j1, j2).derivative()()[i]);

    // Preserve the parameter context from the elements (use first element's dx)
    // to ensure the resulting symmetric matrix derivative carries a valid dx pointer.
    return Derivative<aSymmetricMatrix<double>, Parameters_transformed>(std::move(prim),
                                                                        std::move(der), dx);
}

template <template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
void set(Derivative<aSymmetricMatrix<double>, Parameters_transformed>& x, std::size_t i,
         std::size_t j, double value) {
    d_d_<double, Parameters_transformed> der(
        Matrix<double>(x.derivative()().nrows(), x.derivative()().ncols(), 0.0), x.dx());
    Derivative<double, Parameters_transformed> dv(value, der);
    set(x, i, j, dv);
}

template <template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
void set(Derivative<aSymmetricMatrix<double>, Parameters_transformed>& x, std::size_t i,
         std::size_t j, const Derivative<double, Parameters_transformed>& value) {
    x.primitive().set(i, j, value.primitive());
    if (x.derivative()().size() == 0)
        x.derivative()() = Matrix<aSymmetricMatrix<double>>(
            value.derivative()().nrows(), value.derivative()().ncols(),
            aSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t k = 0; k < x.derivative()().size(); ++k)
        x.derivative()()[k].set(i, j, value.derivative()()[k]);
}
template <template <class> class aNonSymmetricMatrix>
    requires(std::is_same_v<Matrix<double>, aNonSymmetricMatrix<double>>)
void set(Derivative<aNonSymmetricMatrix<double>, Parameters_transformed>& x, std::size_t i,
         std::size_t j, double value) {
    // Preserve parameter context when injecting a scalar into a matrix derivative.
    // Without passing x.dx(), the temporary derivative would carry a null dx pointer
    // and downstream accesses to Y.dx().parameters() would segfault during diagnostics.
    d_d_<double, Parameters_transformed> der(
        Matrix<double>(x.derivative()().nrows(), x.derivative()().ncols(), 0.0), x.dx());
    Derivative<double, Parameters_transformed> dv(value, der);
    set(x, i, j, dv);
}

template <template <class> class aNonSymmetricMatrix>
    requires(std::is_same_v<Matrix<double>, aNonSymmetricMatrix<double>>)
void set(Derivative<aNonSymmetricMatrix<double>, Parameters_transformed>& x, std::size_t i,
         std::size_t j, const Derivative<double, Parameters_transformed>& value) {
    x.primitive()(i, j) = value.primitive();
    if (x.derivative()().size() == 0)
        x.derivative() = d_d_<aNonSymmetricMatrix<double>, Parameters_transformed>(
            Matrix<aNonSymmetricMatrix<double>>(value.derivative()().nrows(),
                                                value.derivative()().ncols(),
                                                aNonSymmetricMatrix<double>(x.nrows(), x.ncols())),
            value.dx());
    for (std::size_t k = 0; k < x.derivative()().size(); ++k)
        if (value.derivative()().size() > 0)
            x.derivative()()[k](i, j) = value.derivative()()[k];
}

template <>
class Derivative<Parameters_values, Parameters_transformed> {
    Parameters_Transformations const* m_par;

    Derivative<Matrix<double>, Parameters_transformed> m_x;

   public:
    operator Matrix<double>&() { return m_x.primitive(); }

    // operator Id &() { return m_x.primitive(); }
    auto ncols() const { return m_x.ncols(); }
    auto nrows() const { return m_x.nrows(); }
    auto size() const { return m_x.size(); }

    Derivative() = default;

    template <class dParam>
        requires(std::constructible_from<Derivative<Matrix<double>, Parameters_transformed>,
                                         dParam>)
    Derivative(Parameters_Transformations const& tr, dParam&& x)
        : m_par{&tr}, m_x{std::forward<dParam>(x)} {}

    auto& primitive() const { return m_x.primitive(); }
    auto& derivative() const { return m_x.derivative(); }

    auto& dx() const { return m_x.dx(); }

    auto& operator()() const { return m_x; }
    auto& operator()() { return m_x; }

    auto operator[](std::size_t i) const { return (*this)()[i]; }

    auto operator[](std::pair<std::size_t, std::size_t> ij) const { return (*this)()[ij]; }
};

template <>
class Derivative<Parameters_transformed, Parameters_transformed> {
   private:
    Parameters_Transformations const* m_par;

    Derivative<Matrix<double>, Parameters_transformed> m_x;

    auto free_to_all(const Matrix<Derivative<double, Parameters_transformed>>& x) const {
        Parameters_Transformations const& p = m_x.dx().parameters();
        auto& m_fixed = p.transf().fixed_set();
        Matrix<Derivative<double, Parameters_transformed>> out;
        if (x.nrows() > x.ncols())
            out = Matrix<Derivative<double, Parameters_transformed>>(x.nrows() + m_fixed.size(), 1,
                                                                     0.0);
        else
            out = Matrix<Derivative<double, Parameters_transformed>>(1, x.ncols() + m_fixed.size(),
                                                                     0.0);
        assert(out.size() == p.names().size());
        std::size_t i_in = 0;
        std::size_t i_fi = 0;

        for (std::size_t i_out = 0; i_out < out.size(); ++i_out) {
            if ((m_fixed.size() > i_fi) && (i_out == m_fixed[i_fi])) {
                out[i_out] = Derivative<double, Parameters_transformed>(p.standard_values()[i_out],
                                                                        m_x.dx());
                ++i_fi;
            } else {
                out[i_out] = p.transf()[i_out]->inv(x[i_in]());
                // transf() works with derivatives with respect to matrices
                // the information regarding the parameter has to be added
                out[i_out].derivative().set_dx(m_x.dx());
                ++i_in;
            }
        }
        return out;
    }

   public:
    operator Matrix<double>&() { return m_x.primitive(); }

    // operator Id &() { return m_x.primitive(); }
    auto ncols() const { return m_x.ncols(); }
    auto nrows() const { return m_x.nrows(); }
    auto size() const { return m_x.size(); }

    auto& parameters() const { return *m_par; }

    Derivative() = default;
    Derivative(const Parameters_Transformations& par,
               Derivative<Matrix<double>, Parameters_transformed>&& der)
        : m_par{&par}, m_x{std::move(der)} {}

    auto& primitive() const { return m_x.primitive(); }
    auto& derivative() const { return m_x.derivative(); }

    auto& dx() const { return m_x.dx(); }

    auto& operator()() const { return m_x; }

    Derivative<Parameters_values, Parameters_transformed> to_value() {
        auto v = inside_out((*this)());
        auto out = free_to_all(v);
        return Derivative<Parameters_values, Parameters_transformed>(parameters(),
                                                                     outside_in(out, m_x.dx()));
    }
};

inline Derivative<Parameters_transformed, Parameters_transformed> selfDerivative(
    const Parameters_transformed& x) {
    auto out = Matrix<Matrix<double>>(x().nrows(), x().ncols(),
                                      Matrix<double>(x().nrows(), x().ncols(), 0.0));
    for (std::size_t i = 0; i < x().size(); i++) {
        out[i][i] = 1.0;
    }

    return Derivative<Parameters_transformed, Parameters_transformed>(
        x.parameters(), Derivative<Matrix<double>, Parameters_transformed>(
                            x(), d_d_<Matrix<double>, Parameters_transformed>(out, x)));
}

template <class T>
auto diag(const Derivative<Matrix<T>, Parameters_transformed>& a) {
    return Derivative<DiagonalMatrix<T>, Parameters_transformed>(
        diag(primitive(a)), apply_par([](auto const& d) { return diag(d); }, derivative(a)));
}

template <class T>
auto diagpos(const Derivative<Matrix<T>, Parameters_transformed>& a) {
    return Derivative<DiagPosDetMatrix<T>, Parameters_transformed>(
        diagpos(primitive(a)), apply_par([](auto const& d) { return diag(d); }, derivative(a)));
}

inline auto XTX(const Derivative<Matrix<double>, Parameters_transformed>& a) {
    auto& f = primitive(a);
    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        XTX(f), apply_par([&f](auto const& d) { return X_plus_XT(tr(d) * f); }, derivative(a)));
}

inline auto X_plus_XT(const Derivative<Matrix<double>, Parameters_transformed>& a) {
    auto& f = primitive(a);
    return Derivative<SymmetricMatrix<double>, Parameters_transformed>(
        X_plus_XT(f), apply_par([](auto const& d) { return X_plus_XT(d); }, derivative(a)));
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto tr(const Derivative<aMatrix<double>, Parameters_transformed>& a) {
    auto& f = primitive(a);
    return Derivative<aMatrix<double>, Parameters_transformed>(
        tr(f), apply_par([](auto const& d) { return tr(d); }, derivative(a)));
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDiv(const Derivative<aMatrix<double>, Parameters_transformed>& a,
             const Derivative<aMatrix<double>, Parameters_transformed>& b) {
    return zip([](auto& x, auto& y) { return x / y; }, a, b);
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDivSafe(const Derivative<aMatrix<double>, Parameters_transformed>& a,
                 const Derivative<aMatrix<double>, Parameters_transformed>& b, double eps) {
    return zip([eps](auto& x, auto& y) { return elemDivSafe(x, y, eps); }, a, b);
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDivSoftAbs(const Derivative<aMatrix<double>, Parameters_transformed>& a,
                    const Derivative<aMatrix<double>, Parameters_transformed>& b, double eps) {
    return zip([eps](auto& x, auto& y) { return elemDivSoftAbs(x, y, eps); }, a, b);
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemMult(const Derivative<aMatrix<double>, Parameters_transformed>& a,
              const Derivative<aMatrix<double>, Parameters_transformed>& b) {
    return zip([](auto& x, auto& y) { return x * y; }, a, b);
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>, Parameters_transformed>& a,
                const Derivative<bMatrix<double>, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto& fa = primitive(a);
    auto& fb = primitive(b);
    if (derivative(b)().size() == 0)
        return TranspMult(a, primitive(b));
    else if (derivative(a)().size() == 0)
        return TranspMult(primitive(a), b);
    else
        return Derivative<S, Parameters_transformed>(
            TranspMult(fa, fb),
            zip_par([&fa, &fb](auto const& da,
                               auto const& db) { return TranspMult(fa, db) + TranspMult(da, fb); },
                    derivative(a), derivative(b)));
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator+(const Derivative<aMatrix<double>, Parameters_transformed>& a,
               const Derivative<bMatrix<double>, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(aMatrix<double>{} + bMatrix<double>{})>;
    auto& fa = primitive(a);
    auto& fb = primitive(b);
    return Derivative<S, Parameters_transformed>(
        fa + fb, zip_par([](auto const& da, auto const& db) { return da + db; }, derivative(a),
                         derivative(b)));
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>, Parameters_transformed>& a,
                const bMatrix<double>& b) {
    using S = std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto& fa = primitive(a);
    auto& fb = b;

    return Derivative<S, Parameters_transformed>(
        TranspMult(fa, fb),
        apply_par([&fb](auto const& da) { return TranspMult(da, fb); }, derivative(a)));
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const aMatrix<double>& a,
                const Derivative<bMatrix<double>, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto& fa = a;
    auto& fb = primitive(b);

    return Derivative<S, Parameters_transformed>(
        TranspMult(fa, fb),
        apply_par([&fa](auto const& db) { return TranspMult(fa, db); }, derivative(b)));
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed>& a,
               const Derivative<bMatrix<double>, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto& fa = primitive(a);
    auto& fb = primitive(b);

    if ((derivative(a)().size() > 0) && (derivative(b)().size() > 0))

        return Derivative<S, Parameters_transformed>(
            fa * fb,
            zip_par([&fa, &fb](auto const& da, auto const& db) { return fa * db + da * fb; },
                    derivative(a), derivative(b)));
    else if ((derivative(a)().size() == 0) && (derivative(b)().size() == 0)) {
        // Both differentials are structurally empty. Preserve parameter context to avoid
        // constructing a Derivative<> without a dx pointer, which would later trip assertions
        // in apply_par/operator*.
        if (a.has_dx()) {
            d_d_<S, Parameters_transformed> zero_df(static_cast<const S&>(fa * fb), a.dx());
            return Derivative<S, Parameters_transformed>(fa * fb, zero_df);
        } else if (b.has_dx()) {
            d_d_<S, Parameters_transformed> zero_df(static_cast<const S&>(fa * fb), b.dx());
            return Derivative<S, Parameters_transformed>(fa * fb, zero_df);
        } else {
            // No parameter context available; return a derivative wrapper with no dx.
            // Downstream code should treat this as non-differentiable.
            return Derivative<S, Parameters_transformed>(fa * fb);
        }
    } else if (derivative(a)().size() == 0)
        return Derivative<S, Parameters_transformed>(
            fa * fb, apply_par([&fa](auto const& db) { return fa * db; }, derivative(b)));
    else
        return Derivative<S, Parameters_transformed>(
            fa * fb, apply_par([&fb](auto const& da) { return da * fb; }, derivative(a)));
}
template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed>& a, bMatrix<double>& b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto& fa = primitive(a);
    auto& fb = b;

    return Derivative<S, Parameters_transformed>(
        fa * fb, apply_par([&fa, &fb](auto const& da) { return da * fb; }, derivative(a)));
}

template <template <class> class aMatrix, template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const bMatrix<double>& a,
               const Derivative<aMatrix<double>, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto& fa = a;
    auto& fb = primitive(b);

    return Derivative<S, Parameters_transformed>(
        fa * fb, apply_par([&fa](auto const& db) { return fa * db; }, derivative(b)));
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed>& a,
               const Derivative<double, Parameters_transformed>& b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * double{})>;
    auto& fa = primitive(a);
    auto& fb = primitive(b);

    return Derivative<S, Parameters_transformed>(
        fa * fb, zip_par([&fa, &fb](auto const& da, auto const& db) { return fa * db + da * fb; },
                         derivative(a), derivative(b)));
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<double, Parameters_transformed>& b,
               const Derivative<aMatrix<double>, Parameters_transformed>& a) {
    using S = std::decay_t<decltype(double{} * std::declval<aMatrix<double>>())>;
    auto& fa = primitive(a);
    auto& fb = primitive(b);

    return Derivative<S, Parameters_transformed>(
        fb * fa, zip_par([&fa, &fb](auto const& da, auto const& db) { return db * fa + fb * da; },
                         derivative(a), derivative(b)));
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator/(const Derivative<aMatrix<double>, Parameters_transformed>& a,
               const Derivative<double, Parameters_transformed>& b) {
    return a * (1.0 / b);
}

inline Derivative<SymPosDefMatrix<double>, Parameters_transformed> AT_B_A(
    const Derivative<Matrix<double>, Parameters_transformed>& a,
    const Derivative<SymmetricMatrix<double>, Parameters_transformed>& b) {
    auto& fa = primitive(a);
    auto& fb = primitive(b);

    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        AT_B_A(fa, fb), zip_par(
                            [&fa, &fb](auto const& da, auto const& db) {
                                return AT_B_A(fa, db) + X_plus_XT(TranspMult(da, fb * fa));
                            },
                            derivative(a), derivative(b)));
}

inline Derivative<SymPosDefMatrix<double>, Parameters_transformed> AT_B_A(
    const Derivative<Matrix<double>, Parameters_transformed>& a, const SymmetricMatrix<double>& b) {
    auto& fa = primitive(a);
    auto& fb = b;

    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        AT_B_A(fa, fb),
        apply_par([&fa, &fb](auto const& da) { return X_plus_XT(TranspMult(da, fb * fa)); },
                  derivative(a)));
}

template <template <class> class Matrix>
    requires Matrix<double>::is_Matrix
auto getvalue(const Derivative<Matrix<double>, Parameters_transformed>& x) {
    assert(x.primitive().size() == 1);
    Matrix<double> der(x.derivative()().nrows(), x.derivative()().ncols());
    for (std::size_t i = 0; i < der.size(); ++i) der[i] = x.derivative()()[i][0];
    return Derivative<double, Parameters_transformed>(x.primitive()[0], std::move(der), x.dx());
}

template <template <class> class Matrix>
    requires Matrix<double>::is_Matrix
Maybe_error<Derivative<Matrix<double>, Parameters_transformed>> inv(
    const Derivative<Matrix<double>, Parameters_transformed>& x) {
    auto inv_x = inv(x.primitive());
    if (!inv_x)
        return inv_x.error();
    else {
        auto dinv = apply([&inv_x](auto const& dx) { return -inv_x.value() * dx * inv_x.value(); },
                          x.derivative()());
        return Derivative<Matrix<double>, Parameters_transformed>(inv_x.value(), dinv, x.dx());
    }
}

inline Maybe_error<std::tuple<Derivative<DiagonalMatrix<double>, Parameters_transformed>,
                              Derivative<Matrix<double>, Parameters_transformed>,
                              Derivative<Matrix<double>, Parameters_transformed>>>
    eigs(const Derivative<Matrix<double>, Parameters_transformed>& x,
         bool does_permutations = false, bool does_diagonal_scaling = false,
         bool computes_eigenvalues_condition_numbers = false,
         bool computes_eigenvectors_condition_numbers = false) {
    auto res =
        eigs(x.primitive(), does_permutations, does_diagonal_scaling,
             computes_eigenvalues_condition_numbers, computes_eigenvectors_condition_numbers);

    if (!res) {
        return res.error();
    }
    auto [lambda, VR, VL] = std::move(res).value();

    auto G = apply([&VR, &VL](auto const& x_k) { return tr(VL) * x_k * VR; }, x.derivative()());

    auto der_lambda = apply([](auto const& dx) -> DiagonalMatrix<double> { return diag(dx); }, G);

    auto der_VR = apply(
        [&VR, &lambda](Matrix<double> const& Gk) {
            auto omega_k = Gk - Gk;
            for (std::size_t j = 0; j < Gk.nrows(); ++j) {
                for (std::size_t k = 0; k < Gk.ncols(); ++k) {
                    if (std::abs(lambda[k] - lambda[j]) > 100 * sqrt(eps)) {
                        omega_k(j, k) = Gk(j, k) / (lambda[k] - lambda[j]);
                    }
                }
            }
            return VR * omega_k;
        },
        G);

    auto der_VL = apply(
        [&VL, &lambda](Matrix<double> const& Gk) {
            auto omega_k = Gk - Gk;
            for (std::size_t i = 0; i < Gk.nrows(); ++i) {
                for (std::size_t j = 0; j < Gk.ncols(); ++j) {
                    if (std::abs(lambda[i] - lambda[j]) > 100 * sqrt(eps)) {
                        omega_k(i, j) = Gk(i, j) / (lambda[j] - lambda[i]);
                    }
                }
            }
            return VL * tr(omega_k) * (-1.0);
        },
        G);

    auto dlambda = Derivative<DiagonalMatrix<double>, Parameters_transformed>(
        std::move(lambda), std::move(der_lambda), x.dx());

    auto dVR = Derivative<Matrix<double>, Parameters_transformed>(std::move(VR), std::move(der_VR),
                                                                  x.dx());
    auto dVL = Derivative<Matrix<double>, Parameters_transformed>(std::move(VL), std::move(der_VL),
                                                                  x.dx());

    // Match primitive order: (L, VR, VL)
    return std::tuple(std::move(dlambda), std::move(dVR), std::move(dVL));
}

// Enforce CTMC-generator conventions on eigendecomposition (derivative):
// - deterministically place the zero eigenvalue at index 0
// - apply the corresponding column permutation to VR and VL, and reorder lambda
// - zero-mode gauge: set dλ0 = 0 and du0 = 0 for all directions
inline Maybe_error<std::tuple<Derivative<Matrix<double>, Parameters_transformed>,
                              Derivative<DiagonalMatrix<double>, Parameters_transformed>,
                              Derivative<Matrix<double>, Parameters_transformed>>>
    eig_enforce_q_mode(const Matrix<double>& /*Q*/,
                       const std::tuple<Derivative<Matrix<double>, Parameters_transformed>,
                                        Derivative<DiagonalMatrix<double>, Parameters_transformed>,
                                        Derivative<Matrix<double>, Parameters_transformed>>& de,
                       double tol = 1e-12) {
    auto dVR = std::get<0>(de);
    auto dL = std::get<1>(de);
    auto dVL = std::get<2>(de);

    auto& VR = dVR.primitive();
    auto& L = dL.primitive();
    auto& VL = dVL.primitive();
    const std::size_t n = VR.ncols();
    if (n == 0)
        return de;

    // Determine indices belonging to the (near-)zero eigenvalue cluster
    std::vector<std::size_t> zero_idx;
    zero_idx.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (std::abs(L[i]) <= tol)
            zero_idx.push_back(i);

    // If none detected (defensive), fall back to picking the smallest-in-magnitude
    if (zero_idx.empty()) {
        std::size_t k0 = 0;
        double amin = std::abs(L[0]);
        for (std::size_t i = 1; i < n; ++i) {
            if (std::abs(L[i]) < amin) {
                amin = std::abs(L[i]);
                k0 = i;
            }
        }
        zero_idx.push_back(k0);
    }

    // Build permutation bringing the entire zero cluster to the front, preserving order
    Matrix<double> Per(n, n, 0.0);
    std::vector<std::size_t> order(n);
    std::size_t w = 0;
    // first, the zero-cluster indices in their original order
    for (auto k : zero_idx) order[w++] = k;
    // then all remaining indices
    for (std::size_t i = 0; i < n; ++i) {
        bool is_zero = std::find(zero_idx.begin(), zero_idx.end(), i) != zero_idx.end();
        if (!is_zero)
            order[w++] = i;
    }
    for (std::size_t i = 0; i < n; ++i) Per(i, order[i]) = 1.0;
    auto PerT = tr(Per);

    // Permute primitives
    auto VRp = VR * PerT;
    auto VLp = VL * PerT;
    DiagonalMatrix<double> Lp(L.nrows(), L.ncols());
    for (std::size_t i = 0; i < n; ++i) Lp[i] = L[order[i]];

    // Permute derivatives
    auto permute_cols = [&PerT](const Matrix<double>& M) { return M * PerT; };

    auto dVR_new =
        apply_par([&permute_cols](auto const& d) { return permute_cols(d); }, dVR.derivative());
    auto dVL_new =
        apply_par([&permute_cols](auto const& d) { return permute_cols(d); }, dVL.derivative());
    auto dL_new = apply_par(
        [&order, n](auto const& d) {
            DiagonalMatrix<double> out(n, n);
            for (std::size_t i = 0; i < n; ++i) out[i] = d[order[i]];
            return out;
        },
        dL.derivative());

    // Zero-mode gauge: set all zero-cluster eigenvalue derivatives to 0
    const std::size_t m0 = zero_idx.size();
    for (std::size_t k = 0; k < dL_new().size(); ++k) {
        if (dL_new()[k].nrows() > 0) {
            for (std::size_t i = 0; i < std::min(m0, dL_new()[k].size()); ++i) dL_new()[k][i] = 0.0;
        }
    }

    // Rebuild derivatives with permutation applied (no scaling yet)
    auto dVR_out = Derivative<Matrix<double>, Parameters_transformed>(VRp, dVR_new(), dVR.dx());
    auto dL_out = Derivative<DiagonalMatrix<double>, Parameters_transformed>(Lp, dL_new(), dL.dx());
    auto dVL_out = Derivative<Matrix<double>, Parameters_transformed>(VLp, dVL_new(), dVL.dx());

    // Biorthogonal rescaling on primitive left eigenvectors so v_i^T u_i = 1,
    // with consistent derivative adjustment for the scaling factor s = v_i^T u_i.
    for (std::size_t i = 0; i < n; ++i) {
        // Compute inner product s before scaling primitives
        auto ui = dVR_out.primitive()(":", i);
        auto vcol = dVL_out.primitive()(":", i);
        auto vTi_pre = tr(vcol);
        double s = getvalue(vTi_pre * ui);
        if (std::abs(s) <= tol)
            return error_message("eig_enforce_q_mode(d): singular (v_i^T u_i == 0)");

        // Scale left eigenvector primitive column by 1/s
        for (std::size_t r = 0; r < dVL_out.primitive().nrows(); ++r)
            dVL_out.primitive()(r, i) = dVL_out.primitive()(r, i) / s;

        // Adjust derivatives to account for scaling: v'_i = v_i / s
        // d v'_i = (1/s) d v_i - (ds/s) v'_i, where ds = dv_i^T u_i + v_i^T du_i (computed pre-scale)
        for (std::size_t k = 0; k < dVL_out.derivative()().size(); ++k) {
            // Skip empty derivative directions
            if (dVL_out.derivative()()[k].nrows() == 0)
                continue;

            // Column i of derivative matrices before scaling adjustment
            auto dvi_col = dVL_out.derivative()()[k](
                ":", i);  // currently unscaled derivative; we'll overwrite
            auto dui_col = dVR_out.derivative()()[k](":", i);

            double ds_k = getvalue(tr(dvi_col) * ui + vTi_pre * dui_col);

            // Update column entries in-place using v'_i now stored in primitive
            for (std::size_t r = 0; r < dVL_out.derivative()()[k].nrows(); ++r) {
                double dv_ir = dVL_out.derivative()()[k](
                    r, i);  // this is still unscaled value from previous line
                double vprime_r = dVL_out.primitive()(r, i);  // scaled primitive value (v_i / s)
                dVL_out.derivative()()[k](r, i) = dv_ir / s - (ds_k / s) * vprime_r;
            }
        }
    }

    return std::tuple(dVR_out, dL_out, dVL_out);
}

inline auto Taylor_first(Derivative<double, Parameters_transformed> const& f,
                         const Parameters_transformed& x, double eps) {
    auto dx = x;
    dx() = dx() * eps;
    return primitive(f) + f.derivative() * dx;
}

template <template <class> class aMatrix>
auto Taylor_first(Derivative<aMatrix<double>, Parameters_transformed> const& f,
                  const Parameters_transformed& x, double eps) {
    auto dx = x;
    dx() = dx() * eps;
    return primitive(f) + f.derivative() * dx;
}

inline auto Taylor_first(Derivative<Parameters_transformed, Parameters_transformed> const& f,
                         const Parameters_transformed& x, double eps) {
    return Parameters_transformed(x.parameters(), primitive(f) + f.derivative() * x * eps);
}

}  // namespace var

#endif  // PARAMETERS_DERIVATIVE_H
