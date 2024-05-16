#ifndef PARAMETERS_DERIVATIVE_H
#define PARAMETERS_DERIVATIVE_H

#include "derivative_operator.h"
#include "matrix.h"
#include "maybe_error.h"
#include "parameters.h"
#include <cstddef>
#include <vector>
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
            auto &Dx = *it_val;
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
        auto &Dx = *it_val;
        x(i, j) = primitive(Dx);
        ++it_ind;
        ++it_val;
    }
    
    return Derivative<Matrix, X>(x, dfdx, dx);
}

template <class M> struct M_der {
    using type = M;
};

template <> struct M_der<SymPosDefMatrix<double>> {
    using type = SymmetricMatrix<double>;
};

template <> struct M_der<DiagPosDetMatrix<double>> {
    using type = DiagonalMatrix<double>;
};

template <class M> using M_der_t = typename M_der<M>::type;

 template<> class d_d_<double, Parameters_transformed> {
    Matrix<double> m_dydx;
    Parameters_transformed const *ptr_dx;
    
public:
    using value_type = Matrix<double>;
    template <class aMatrix>
        requires std::is_same_v<Matrix<double>, std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix &&dydx, Parameters_transformed const &x)
        : m_dydx{std::forward<aMatrix>(dydx)}, ptr_dx{&x} {}
    
    template <class aMatrix>
        requires std::is_same_v<Matrix<double>, std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix &&dydx, Parameters_transformed const *x)
        : m_dydx{std::forward<aMatrix>(dydx)}, ptr_dx{x} {}
    d_d_() {}
    
    auto &operator()() { return m_dydx; }
    auto &operator()() const { return m_dydx; }
    
    auto &dx() const { return *ptr_dx; }
    bool has_dx()const {return ptr_dx!=nullptr;}
    void set_dx(Parameters_transformed const& x)  {
        ptr_dx=&x;
    }
    
    
    friend auto operator*(d_d_ const &df, const Parameters_transformed &x) {
        double out = 0;
        for (std::size_t i = 0; i < df().size(); ++i)
            out += df()[i] * x()[i];
        return out;
    }
};

template < template <class> class aMatrix>
    requires(aMatrix<double>::is_Matrix)
class d_d_<aMatrix<double>, Parameters_transformed> {
    Matrix<aMatrix<double>> m_dydx;
    Parameters_transformed const *ptr_par;
    
public:
    using value_type = Matrix<aMatrix<double>>;
    d_d_() {}
    constexpr ~d_d_()=default;
    template <class aaMatrix>
        requires std::constructible_from<Matrix<aMatrix<double>>,
                                         std::decay_t<aaMatrix>>
    constexpr d_d_(aaMatrix &&dydx, const Parameters_transformed &x)
        : m_dydx{std::forward<aaMatrix>(dydx)}, ptr_par{&x} {}
    
    constexpr d_d_(aMatrix<double> const &y, const Parameters_transformed &x)
        : m_dydx{x().nrows(), x().ncols(),aMatrix<double>(y.nrows(),y.ncols(),0.0)}, ptr_par{&x} {}
    
    
    
    
    constexpr auto &operator()() { return m_dydx; }
    constexpr auto &operator()() const { return m_dydx; }
    
    auto &dx() const { return *ptr_par; }
    bool has_dx()const {return ptr_par!=nullptr;}
    void set_dx(Parameters_transformed const& x)  {
        ptr_par=&x;
    }
    
    
    template <class F> friend auto apply_par(F &&f, d_d_ const &a) {
        using S = std::decay_t<std::invoke_result_t<F, aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(a()[i]);
        return d_d_<S, Parameters_transformed>(x, a.dx());
    }
    
    template <class F, template <class> class bMatrix>
    friend auto zip_par(F &&f, d_d_ const &a,
                        d_d_<bMatrix<double>, Parameters_transformed> const &b) {
        using S =
            std::decay_t<std::invoke_result_t<F, aMatrix<double>, bMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        
        for (std::size_t i = 0; i < std::min(x.size(), b().size()); ++i)
            x[i] = f(a()[i], b()[i]);
        return d_d_<S, Parameters_transformed>(x, b.dx());
    }
    
    template <class F>
    friend auto zip_par(F &&f, d_d_ const &a,
                        d_d_<double, Parameters_transformed> const &b) {
        using S = std::decay_t<std::invoke_result_t<F, aMatrix<double>, double>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < std::min(x.size(), b().size()); ++i)
            x[i] = f(a()[i], b()[i]);
        return d_d_<S, Parameters_transformed>(x, b.dx());
    }
    
    template <class F>
    friend auto zip_par(F &&f, d_d_<double, Parameters_transformed> const &b,
                        d_d_ const &a) {
        using S = std::decay_t<std::invoke_result_t<F, double, aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(b()[i], a()[i]);
        return d_d_<S, Parameters_transformed>(x, a.dx());
    }
    
    friend auto operator*(d_d_ const &df, const Parameters_transformed &x) {
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

template < template <class> class notSymmetricMatrix>
    requires(
        (notSymmetricMatrix<double>::is_Matrix) &&
        (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
           std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto inside_out(const d_d_<notSymmetricMatrix<double>, Parameters_transformed> &x) {
    if (x().size() == 0)
        return notSymmetricMatrix<d_d_<double, Parameters_transformed>>{};
    else {
        notSymmetricMatrix<d_d_<double, Parameters_transformed>> out(x()[0].nrows(),
                                                             x()[0].ncols());
        for (std::size_t i = 0; i < out.size(); ++i) {
            Matrix<double> d_d_par(x().nrows(), x().ncols());
            for (std::size_t j = 0; j < d_d_par.size(); ++j)
                d_d_par[j] = x()[j][i];
            out[i] = d_d_<double, Parameters_transformed>(std::move(d_d_par), x.dx());
        }
        return out;
    }
}

template < template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto inside_out(const d_d_<aSymmetricMatrix<double>, Parameters_transformed> &x) {
    aSymmetricMatrix<d_d_<double, Parameters_transformed>> out(x()[0].nrows(),
                                                       x()[0].ncols());
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            Matrix<double> d_d_par(x().nrows(), x().ncols());
            for (std::size_t k = 0; k < d_d_par.size(); ++k)
                d_d_par[k] = x()[k](i, j);
            out.set(i, j, d_d_<double, Parameters_transformed>(std::move(d_d_par), x.dx()));
        }
    return out;
}

template <> class d_d_<Parameters_transformed, Parameters_transformed> {
    d_d_<Matrix<double>, Parameters_transformed> m_dydx;
    
public:
    using value_type = d_d_<Matrix<double>, Parameters_transformed>;
    d_d_() {}
    auto &dx() const { return m_dydx.dx(); }
    template <class aMatrix>
        requires std::constructible_from<d_d_<Matrix<double>, Parameters_transformed>,
                                         aMatrix>
    constexpr d_d_(aMatrix &&dydx) : m_dydx{std::forward<aMatrix>(dydx)} {}
    
    constexpr auto &operator()() { return m_dydx; }
    constexpr auto &operator()() const { return m_dydx; }
};

template <> class Derivative<double, Parameters_transformed> //: public Primitive<double>, public
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
    
    Derivative(P t_x, D &&t_d) : m_x{t_x}, m_d{std::forward<D>(t_d)} {}
    
    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D,Parameters_transformed const &>)
    
    Derivative(P t_x, D &&t_d, Parameters_transformed const &x)
        : m_x{t_x}, m_d{std::forward<D>(t_d), x} {}
    Derivative() {}
    
    Derivative(double t_x,  Parameters_transformed const &x)
        : m_x{t_x}, m_d{Matrix<double>(x().nrows(),x().ncols(),0.0), x} {}
    
    Derivative& operator=(const Derivative<double,Matrix<double>>& x)
    {
        m_x=x.primitive();
        m_d()=x.derivative()();
        return *this;        
    }
    
    Derivative(double t_x) : m_x{t_x}, m_d{} {}
    
    Derivative(const Derivative& x)=default;
    Derivative(Derivative&& x)=default;
    Derivative& operator=(const Derivative& x)=default;
    Derivative& operator=(Derivative&& x)=default;
    
    auto &primitive() { return m_x; }
    auto &primitive() const { return m_x; }
    auto &derivative() const { return m_d; }
    auto &derivative()  { return m_d; }
    auto &dx() const { return m_d.dx(); }
    
    Derivative<double, Matrix<double>> operator()()const
    {
        return Derivative<double, Matrix<double>>(primitive(),derivative()(),dx()());
    }
    
    friend auto operator/(const Derivative &x, const Derivative &y) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(
                fx / fy,
                zip([fx, fy](auto dx,
                             auto dy) { return dx / fy - fx / fy / fy * dy; },
                    x.derivative()(), y.derivative()()),
                x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(fx / fy, -fx / fy / fy * y.derivative()(), x.dx());
        else
            return Derivative(fx / fy, x.derivative()() / fy, x.dx());
    }
    
    friend auto elemDivSafe(const Derivative &x, const Derivative &y,
                            double eps) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(elemDivSafe(fx, fy, eps),
                              zip(
                                  [fx, fy, eps](auto dx, auto dy) {
                                      return elemDivSafe(dx, fy, eps) -
                                             elemDivSafe(fx, fy * fy, eps) * dy;
                                  },
                                  x.derivative()(), y.derivative()()),
                              x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(elemDivSafe(fx, fy, eps),
                              -elemDivSafe(fx, fy * fy, eps) * y.derivative()(),
                              x.dx());
        else
            return Derivative(elemDivSafe(fx, fy, eps),
                              elemDivSafe(x.derivative()(), fy, eps), x.dx());
    }
    
    friend auto operator*(const Derivative &x, const Derivative &y) {
        auto fx = x.primitive();
        auto fy = y.primitive();
        if ((x.derivative()().size() > 0) && (y.derivative()().size() > 0))
            return Derivative(
                fx * fy,
                zip([fx, fy](auto dx, auto dy) { return dx * fy + fx * dy; },
                    x.derivative()(), y.derivative()()),
                x.dx());
        else if (x.derivative()().size() == 0)
            return Derivative(fx * fy, fx * y.derivative()(), x.dx());
        else
            return Derivative(fx * fy, x.derivative()() * fy, x.dx());
    }
    
    friend auto exp(const Derivative &x) {
        auto f = exp(x.primitive());
        return Derivative(f, f * x.derivative()(), x.dx());
    }
    
    friend auto max(const Derivative &x, double y) {
        if (x.primitive() >= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }
    
    friend auto min(const Derivative &x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }
    
    friend auto log(const Derivative &x) {
        auto f = log(x.primitive());
        return Derivative(f, x.derivative()() * (1.0 / x.primitive()), x.dx());
    }
    friend auto log10(const Derivative &x) {
        auto f = log10(x.primitive());
        return Derivative(f, x.derivative()() *
                                 (1.0 / (x.primitive() * std::log(10))), x.dx());
    }
    friend auto pow(double base, const Derivative &x) {
        using std::pow;
        auto f = pow(base, x.primitive());
        return Derivative(f, x.derivative()() * f * std::log(base), x.dx());
    }
    friend auto pow(const Derivative &base, const Derivative &x) {
        using std::pow;
        auto f = pow(base.primitive(), x.primitive());
        return Derivative(f,
                          x.derivative()() * f * std::log(base.primitive()) +
                              base.derivative()() * x.primitive() *
                                                                                  pow(base.primitive(), x.primitive() - 1.0),
                          x.dx());
    }
    
    friend auto abs(const Derivative &x) {
        auto f = std::abs(x.primitive());
        return Derivative(
            f, ((x.primitive() > 0.0) ? 1.0 : ((x.primitive() < 0) ? -1.0 : 0.0)) *
                   x.derivative()(),x.dx());
    }
    
    friend bool operator==(Derivative const &one, double val) {
        return one.primitive() == val;
    }
};

template < template <class> class T_Matrix>
    requires(T_Matrix<double>::is_Matrix)
class Derivative<T_Matrix<double>,
                 Parameters_transformed> { //: public T_Matrix<double>{
    using primitive_type = T_Matrix<double>;
    using d_type = Parameters_transformed;
    using derivative_type = d_d_<M_der_t<primitive_type>, d_type>;
    primitive_type m_x;
    derivative_type m_d;
    
public:
    
    auto ncols() const { return m_x.ncols(); } //{return primitive_type::ncols();}
    auto nrows() const {
        return m_x.nrows();
    }                                        // {return primitive_type::nrows();}
    auto size() const { return m_x.size(); } //{return primitive_type::size();}
    Derivative() {}
    
    //   using gserg=typename primitive_type::sgr;
    //   using gserug=typename derivative_type::sgrdd;
    
    Derivative(std::size_t nrows, std::size_t ncols) : m_x(nrows, ncols) {}
    Derivative(std::size_t nrows, std::size_t ncols, double v)
        : m_x(nrows, ncols, v) {}
    auto &dx() const { return m_d.dx(); }
    
    template <class F> friend auto apply(F &&f, const Derivative &x) {
        return outside_in(apply(std::forward<F>(f), inside_out(x)));
    }
    
    template <class F>
    friend auto zip(F &&f, const Derivative &x, const Derivative &y) {
        
        return outside_in(zip(std::forward<F>(f), inside_out(x), inside_out(y)));
    }
    
    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D>)
    Derivative(P &&t_x, D &&t_d)
        : m_x{std::forward<P>(t_x)}, m_d{std::forward<D>(t_d)} {}
    
    template <class P, class D>
        requires(std::constructible_from<primitive_type, P> &&
                 std::constructible_from<derivative_type, D,
                                                                                       const Parameters_transformed &>)
    Derivative(P &&t_x, D &&t_d, const Parameters_transformed &x)
        : m_x{std::forward<P>(t_x)}, m_d{std::forward<D>(t_d), x} {}
    
    template <class anotherCompatibleMatrix>
    Derivative(Derivative<anotherCompatibleMatrix, d_type> const &t_x)
        : m_x{t_x.primitive()}, m_d{t_x.derivative()(), t_x.dx()} {}
    
    template <class P>
        requires(std::constructible_from<primitive_type, P>)
    Derivative(P &&t_x) : m_x{std::forward<P>(t_x)}, m_d{} {}
    
    auto &primitive_non_const() {
        return m_x;
    } // {return static_cast<primitive_type&>(*this);}
    auto &primitive() {
        return m_x;
    } // {return static_cast<primitive_type&>(*this);}
    auto &primitive() const {
        return m_x;
    } //{return static_cast<primitive_type const&>(*this);}
    auto &derivative() { return m_d; }
    auto &derivative_non_const() { return m_d; }
    auto &derivative() const { return m_d; }
    
    auto& dx(){return derivative().dx();}
    
    friend auto operator*(const Derivative &x, double y) {
        return Derivative(x.primitive() * y, x.derivative()() * y, x.dx());
    }
    
    
    template <template <class> class anotherCompatibleMatrix>
    friend auto
    operator*(const Derivative &x,
              anotherCompatibleMatrix<Derivative<double, d_type>> const &t_y) {
        auto Y = outside_in(t_y);
        return x * Y;
    }
    
    template <template <class> class anotherCompatibleMatrix>
    friend auto
    operator*(anotherCompatibleMatrix<Derivative<double, d_type>> const &t_y,
              const Derivative &x) {
        auto Y = outside_in(t_y);
        return Y * x;
    }
    
    auto operator()(std::size_t i, std::size_t j) const {
        return Derivative<double, Parameters_transformed>(
            primitive()(i, j),
            applyMap([i, j](auto const &m) { return m(i, j); }, derivative()()),
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
    void set(std::size_t i, std::size_t j, Derivative<double, Parameters_transformed> const& x ) {
        primitive_non_const()(i, j) = x.primitive();
        derivative_non_const().set_dx(x.dx());
        for (std::size_t ij = 0; ij < derivative()().size(); ++ij)
            derivative_non_const()()[ij](i, j) = x.derivative()()[ij];
    }
    
    
    auto operator[](std::size_t i) const {
        return Derivative<double, Parameters_transformed>(
            primitive()[i],
            applyMap([i](auto const &m) { return m[i]; }, derivative()()), dx());
    }
   
    
    auto operator[](std::pair<std::size_t, std::size_t> ij) const {
        return Derivative<Matrix<double>, Parameters_transformed>(
            primitive()[ij],
            applyMap([ij](auto const &m) { return m[ij]; }, derivative()()), dx());
    }
};

template < template <class> class T_Matrix>
    requires(T_Matrix<double>::is_Matrix)
double fullsum(const Derivative<T_Matrix<double>,Parameters_transformed>& x)
{
    return fullsum(x.primitive())+ fullsum(x.derivative()());
}


template <  class T>
class Derivative<std::vector<T>,
                 Parameters_transformed>: public std::vector<Derivative<T,Parameters_transformed>>
{ //: public Container<double>{
    
    
public:
    using base_type=std::vector<Derivative<T,Parameters_transformed>>;
    using base_type::size;
    using base_type::vector;
    using base_type::operator[];
    using base_type::push_back;
     
    auto primitive()const {
        std::vector<T> out(size());
        for (std::size_t i=0; i<size(); ++i)
            out[i]=primitive((*this)[i]);
        return out; 
    } 
    auto& dx(){return (*this)[0].derivative().dx();}
    
};





template <class F>
    requires requires(Derivative<F, Parameters_transformed> &f) {
        { f.derivative()().nrows() };
    }
auto &get_dx_of_dfdx(const Derivative<F, Parameters_transformed> &f) {
    return f.dx();
}

template < template <class> class notSymmetricMatrix>
    requires(
        (notSymmetricMatrix<double>::is_Matrix) &&
        (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
           std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto inside_out(
    const Derivative<notSymmetricMatrix<double>, Parameters_transformed> &x) {
    notSymmetricMatrix<Derivative<double, Parameters_transformed>> out(x.nrows(),
                                                               x.ncols());
    
    auto der = inside_out(x.derivative());
    if (der.size() > 0)
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = Derivative<double, Parameters_transformed>(x.primitive()[i], der[i]);
    else
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = Derivative<double, Parameters_transformed>(x.primitive()[i]);
    return out;
}

template < template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto inside_out(const Derivative<aSymmetricMatrix<double>, Parameters_transformed> &x) {
    aSymmetricMatrix<Derivative<double, Parameters_transformed>> out(x.nrows(),
                                                             x.ncols());
    auto der = inside_out(x.derivative());
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j)
            out.set(
                i, j,
                Derivative<double, Parameters_transformed>(x.primitive()(i, j), der(i, j)));
    return out;
}

template <template <class> class Matrix>
auto &outside_in(const Matrix<double> &x) {
    return x;
}

template <template <class> class Matrix>
auto &inside_out(const Matrix<double> &x) {
    return x;
}

template < template <class> class notSymmetricMatrix>
    requires(
        (notSymmetricMatrix<double>::is_Matrix) &&
        (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>> ||
           std::is_same_v<SymPosDefMatrix<double>, notSymmetricMatrix<double>>)))

auto outside_in(
    const notSymmetricMatrix<Derivative<double, Parameters_transformed>> &x) {
    auto prim = notSymmetricMatrix<double>(x.nrows(), x.ncols());
    for (std::size_t i = 0; i < prim.size(); ++i)
        prim[i] = x[i].primitive();
    
    auto der = Matrix<notSymmetricMatrix<double>>(
        x[0].derivative()().nrows(), x[0].derivative()().ncols(),
        notSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t i = 0; i < der.size(); ++i)
        for (std::size_t j = 0; j < der[i].size(); ++j)
            der[i][j] = x[j].derivative()()[i];
    
    return Derivative<notSymmetricMatrix<double>, Parameters_transformed>(
        std::move(prim), std::move(der), x[0].dx());
}

template < template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
auto outside_in(const aSymmetricMatrix<Derivative<double, Parameters_transformed>> &x) {
    auto prim = aSymmetricMatrix<double>(x.nrows(), x.ncols());
    for (std::size_t i = 0; i < prim.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j)
            prim.set(i, j, x(i, j).primitive());
    
    auto der = Matrix<aSymmetricMatrix<double>>(
        x[0].derivative()().nrows(), x[0].derivative()().ncols(),
        aSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t i = 0; i < der.size(); ++i)
        for (std::size_t j1 = 0; j1 < der[i].nrows(); ++j1)
            for (std::size_t j2 = 0; j2 <= j1; ++j2)
                der[i].set(j1, j2, x(j1, j2).derivative()()[i]);
    
    return Derivative<aSymmetricMatrix<double>, Parameters_transformed>(std::move(prim),
                                                                std::move(der));
}

template < template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
void set(Derivative<aSymmetricMatrix<double>, Parameters_transformed> &x, std::size_t i,
         std::size_t j, double value) {
    d_d_<double, Parameters_transformed> der(
        Matrix<double>(x.derivative()().nrows(), x.derivative()().ncols(), 0.0),x.dx());
    Derivative<double, Parameters_transformed> dv(value, der);
    set(x, i, j, dv);
}

template < template <class> class aSymmetricMatrix>
    requires(std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>> ||
             std::is_same_v<SymPosDefMatrix<double>, aSymmetricMatrix<double>>)
void set(Derivative<aSymmetricMatrix<double>, Parameters_transformed> &x, std::size_t i,
         std::size_t j, const Derivative<double, Parameters_transformed> &value) {
    x.primitive().set(i, j, value.primitive());
    if (x.derivative()().size() == 0)
        x.derivative()() = Matrix<aSymmetricMatrix<double>>(
            value.derivative()().nrows(), value.derivative()().ncols(),
            aSymmetricMatrix<double>(x.nrows(), x.ncols()));
    for (std::size_t k = 0; k < x.derivative()().size(); ++k)
        x.derivative()()[k].set(i, j, value.derivative()()[k]);
}
template < template <class> class aNonSymmetricMatrix>
    requires(std::is_same_v<Matrix<double>, aNonSymmetricMatrix<double>>)
void set(Derivative<aNonSymmetricMatrix<double>, Parameters_transformed> &x,
         std::size_t i, std::size_t j, double value) {
    d_d_<double, Parameters_transformed> der(
        Matrix<double>(x.derivative()().nrows(), x.derivative()().ncols(), 0.0));
    Derivative<double, Parameters_transformed> dv(value, der);
    set(x, i, j, dv);
}

template < template <class> class aNonSymmetricMatrix>
    requires(std::is_same_v<Matrix<double>, aNonSymmetricMatrix<double>>)
void set(Derivative<aNonSymmetricMatrix<double>, Parameters_transformed> &x,
         std::size_t i, std::size_t j,
         const Derivative<double, Parameters_transformed> &value) {
    x.primitive()(i, j) = value.primitive();
    if (x.derivative()().size() == 0)
        x.derivative() = d_d_<aNonSymmetricMatrix<double> , Parameters_transformed>(
            Matrix<aNonSymmetricMatrix<double>>(value.derivative()().nrows(), value.derivative()().ncols(),
                                                aNonSymmetricMatrix<double>(x.nrows(), x.ncols())), value.dx());
    for (std::size_t k = 0; k < x.derivative()().size(); ++k)
        if (value.derivative()().size() > 0)
            x.derivative()()[k](i, j) = value.derivative()()[k];
}



template<> class Derivative<Parameters_values, Parameters_transformed> {
    Parameters_Transformations const *m_par;
    
    Derivative<Matrix<double>, Parameters_transformed> m_x;
    
public:
    operator Matrix<double> &() { return m_x.primitive(); }
    
   // operator Id &() { return m_x.primitive(); }
    auto ncols() const { return m_x.ncols(); }
    auto nrows() const { return m_x.nrows(); }
    auto size() const { return m_x.size(); }
    
    
    
    Derivative() {}
    
    template <class dParam>
        requires(std::constructible_from<
                    Derivative<Matrix<double>, Parameters_transformed>, dParam>)
    Derivative(Parameters_Transformations const& tr, dParam &&x)
        : m_par{&tr},m_x{std::forward<dParam>(x)} {}
    
 
    
    
       
    auto &primitive() const { return m_x.primitive(); }
    auto &derivative() const { return m_x.derivative(); }
    
    auto &dx() const { return m_x.dx(); }
    
    auto &operator()() const { return m_x; }
    auto &operator()()  { return m_x; }
    
    auto operator[](std::size_t i) const { return (*this)()[i]; }
    
    auto operator[](std::pair<std::size_t, std::size_t> ij) const {
       return (*this)()[ij];
    }
    
};




template<> class Derivative<Parameters_transformed, Parameters_transformed> {
private:
    Parameters_Transformations const *m_par;
    
    
    Derivative<Matrix<double>, Parameters_transformed> m_x;
    
    
    auto free_to_all(const Matrix<Derivative<double, Parameters_transformed>> &x) const {
        Parameters_Transformations const & p=m_x.dx().parameters();
        auto &m_fixed = p.transf().fixed_set();
        Matrix<Derivative<double, Parameters_transformed>> out;
        if (x.nrows() > x.ncols())
            out = Matrix<Derivative<double, Parameters_transformed>>(x.nrows() + m_fixed.size(), 1,0.0);
        else
            out = Matrix<Derivative<double, Parameters_transformed>>(1, x.ncols() + m_fixed.size(),0.0);
        assert(out.size() == p.names().size());
        std::size_t i_in = 0;
        std::size_t i_fi = 0;
        
        for (std::size_t i_out = 0; i_out < out.size(); ++i_out) {
            if ((m_fixed.size() > i_fi) && (i_out == m_fixed[i_fi])) {
                out[i_out] = Derivative<double, Parameters_transformed>(p.standard_values()[i_out],x[0].dx());
                ++i_fi;
            } else {
                out[i_out] = p.transf()[i_out]->inv(x[i_in]());
                // transf() works with derivatives with respect to matrices
                // the information regarding the parameter has to be added
                out[i_out].derivative().set_dx(x[i_in].dx());
                ++i_in;
            }
        }
        return out;
    }
    
    
    
    
    
    
    
public:
    operator Matrix<double> &() { return m_x.primitive(); }
    
   // operator Id &() { return m_x.primitive(); }
    auto ncols() const { return m_x.ncols(); }
    auto nrows() const { return m_x.nrows(); }
    auto size() const { return m_x.size(); }
    
    auto &parameters() const { return *m_par; }
    
    
    Derivative() {}
    Derivative(const Parameters_Transformations& par,Derivative<Matrix<double>, Parameters_transformed>&& der )
        : m_par{&par}, m_x{std::move(der)}    {}
    
    
    auto &primitive() const { return m_x.primitive(); }
    auto &derivative() const { return m_x.derivative(); }
    
    auto &dx() const { return m_x.dx(); }
    
    auto &operator()() const { return m_x; }
    
    
    
    
    Derivative<Parameters_values, Parameters_transformed>
    to_value()
    {
        auto v=inside_out((*this)());
        auto out=free_to_all(v);
        return Derivative<Parameters_values, Parameters_transformed>(parameters(),outside_in(out));
    }
};


Derivative<Parameters_transformed, Parameters_transformed>
selfDerivative(const Parameters_transformed &x) {
    auto out = Matrix<Matrix<double>>(
        x().nrows(), x().ncols(), Matrix<double>(x().nrows(), x().ncols(), 0.0));
    for (std::size_t i = 0; i < x().size(); i++) {
        out[i][i] = 1.0;
    }
    
    return Derivative<Parameters_transformed, Parameters_transformed>(
        x.parameters(),
        Derivative<Matrix<double>, Parameters_transformed>(
            x(), d_d_<Matrix<double>, Parameters_transformed>(out, x)));
}

template <class T>
auto diag(const Derivative<Matrix<T>, Parameters_transformed> &a) {
    
    return Derivative<DiagonalMatrix<T>, Parameters_transformed>(
        diag(primitive(a)),
        apply_par([](auto const &d) { return diag(d); }, derivative(a)));
}

template <class T>
auto diagpos(const Derivative<Matrix<T>, Parameters_transformed> &a) {
    
    return Derivative<DiagPosDetMatrix<T>, Parameters_transformed>(
        diagpos(primitive(a)),
        apply_par([](auto const &d) { return diag(d); }, derivative(a)));
}


auto XTX(const Derivative<Matrix<double>, Parameters_transformed> &a) {
    
    auto &f = primitive(a);
    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        XTX(f), apply_par([&f](auto const &d) { return X_plus_XT(tr(d) * f); },
                  derivative(a)));
}


auto X_plus_XT(const Derivative<Matrix<double>, Parameters_transformed> &a) {
    
    auto &f = primitive(a);
    return Derivative<SymmetricMatrix<double>, Parameters_transformed>(
        X_plus_XT(f),
        apply_par([](auto const &d) { return X_plus_XT(d); }, derivative(a)));
}

template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto tr(const Derivative<aMatrix<double>, Parameters_transformed> &a) {
    
    auto &f = primitive(a);
    return Derivative<aMatrix<double>, Parameters_transformed>(
        tr(f), apply_par([](auto const &d) { return tr(d); }, derivative(a)));
}

template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDiv(const Derivative<aMatrix<double>, Parameters_transformed> &a,
             const Derivative<aMatrix<double>, Parameters_transformed> &b) {
    
    return zip([](auto &x, auto &y) { return x / y; }, a, b);
}

template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDivSafe(const Derivative<aMatrix<double>, Parameters_transformed> &a,
                 const Derivative<aMatrix<double>, Parameters_transformed> &b,
                 double eps) {
    
    return zip([eps](auto &x, auto &y) { return elemDivSafe(x, y, eps); }, a, b);
}

template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemMult(const Derivative<aMatrix<double>, Parameters_transformed> &a,
              const Derivative<aMatrix<double>, Parameters_transformed> &b) {
    
    return zip([](auto &x, auto &y) { return x * y; }, a, b);
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>, Parameters_transformed> &a,
                const Derivative<bMatrix<double>, Parameters_transformed> &b) {
    using S =
        std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    if (derivative(b)().size() == 0)
        return TranspMult(a, primitive(b));
    else if (derivative(a)().size() == 0)
        return TranspMult(primitive(a), b);
    else
        return Derivative<S, Parameters_transformed>(
            TranspMult(fa, fb), zip_par(
                                    [&fa, &fb](auto const &da, auto const &db) {
                return TranspMult(fa, db) +
                       TranspMult(da, fb);
            },
                                        derivative(a), derivative(b)));
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator+(const Derivative<aMatrix<double>, Parameters_transformed> &a,
               const Derivative<bMatrix<double>, Parameters_transformed> &b) {
    using S = std::decay_t<decltype(aMatrix<double>{} + bMatrix<double>{})>;
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    return Derivative<S, Parameters_transformed>(
        fa + fb, zip_par([](auto const &da, auto const &db) { return da + db; },
                         derivative(a), derivative(b)));
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>, Parameters_transformed> &a,
                const bMatrix<double> &b) {
    using S =
        std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto &fa = primitive(a);
    auto &fb = b;
    
    return Derivative<S, Parameters_transformed>(
        TranspMult(fa, fb),
        apply_par([&fb](auto const &da) { return TranspMult(da, fb); },
                  derivative(a)));
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto TranspMult(const aMatrix<double> &a,
                const Derivative<bMatrix<double>, Parameters_transformed> &b) {
    using S =
        std::decay_t<decltype(TranspMult(aMatrix<double>{}, bMatrix<double>{}))>;
    auto &fa = a;
    auto &fb = primitive(b);
    
    return Derivative<S, Parameters_transformed>(
        TranspMult(fa, fb),
        apply_par([&fa](auto const &db) { return TranspMult(fa, db); },
                  derivative(b)));
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed> &a,
               const Derivative<bMatrix<double>, Parameters_transformed> &b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    
    if ((derivative(a)().size() > 0) && (derivative(b)().size() > 0))
        
        return Derivative<S, Parameters_transformed>(
            fa * fb,
            zip_par([&fa, &fb](auto const &da,
                               auto const &db) { return fa * db + da * fb; },
                    derivative(a), derivative(b)));
    else if ((derivative(a)().size() == 0) && (derivative(a)().size() == 0))
        return Derivative<S, Parameters_transformed>(fa * fb);
    else if (derivative(a)().size() == 0)
        return Derivative<S, Parameters_transformed>(
            fa * fb,
            apply_par([&fa](auto const &db) { return fa * db; }, derivative(b)));
    else
        return Derivative<S, Parameters_transformed>(
            fa * fb,
            apply_par([&fb](auto const &da) { return da * fb; }, derivative(a)));
}
template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed> &a,
               bMatrix<double> &b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto &fa = primitive(a);
    auto &fb = b;
    
    return Derivative<S, Parameters_transformed>(
        fa * fb,
        apply_par([&fa, &fb](auto const &da) { return da * fb; }, derivative(a)));
}

template < template <class> class aMatrix,
         template <class> class bMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const bMatrix<double> &b,
               const Derivative<aMatrix<double>, Parameters_transformed> &a) {
    using S = std::decay_t<decltype(aMatrix<double>{} * bMatrix<double>{})>;
    auto &fa = primitive(a);
    auto &fb = b;
    
    return Derivative<S, Parameters_transformed>(
        fb * fa,
        apply_par([&fa, &fb](auto const &da) { return fb * da; }, derivative(a)));
}

template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>, Parameters_transformed> &a,
               const Derivative<double, Parameters_transformed> &b) {
    using S = std::decay_t<decltype(aMatrix<double>{} * double{})>;
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    
    return Derivative<S, Parameters_transformed>(
        fa * fb, zip_par([&fa, &fb](auto const &da,
                                    auto const &db) { return fa * db + da * fb; },
                         derivative(a), derivative(b)));
}




template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<double, Parameters_transformed> &b,
               const Derivative<aMatrix<double>, Parameters_transformed> &a) {
    using S = std::decay_t<decltype(double{} * std::declval<aMatrix<double>>())>;
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    
    return Derivative<S, Parameters_transformed>(
        fb * fa, zip_par([&fa, &fb](auto const &da,
                                    auto const &db) { return db * fa + fb * da; },
                         derivative(a), derivative(b)));
}


template < template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator/(const Derivative<aMatrix<double>, Parameters_transformed> &a,
               const Derivative<double, Parameters_transformed> &b) {
    
    return a *(1.0/b);
}




Derivative<SymPosDefMatrix<double>, Parameters_transformed>
AT_B_A(const Derivative<Matrix<double>, Parameters_transformed> &a,
       const Derivative<SymmetricMatrix<double>, Parameters_transformed> &b) {
    auto &fa = primitive(a);
    auto &fb = primitive(b);
    
    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        AT_B_A(fa, fb), zip_par(
                            [&fa, &fb](auto const &da, auto const &db) {
            return AT_B_A(fa, db) +
                   X_plus_XT(TranspMult(da, fb * fa));
        },
                                derivative(a), derivative(b)));
}


Derivative<SymPosDefMatrix<double>, Parameters_transformed>
AT_B_A(const Derivative<Matrix<double>, Parameters_transformed> &a,
       const SymmetricMatrix<double> &b) {
    auto &fa = primitive(a);
    auto &fb = b;
    
    return Derivative<SymPosDefMatrix<double>, Parameters_transformed>(
        AT_B_A(fa, fb), apply_par(
                            [&fa, &fb](auto const &da) {
            return X_plus_XT(TranspMult(da, fb * fa));
        },
                                  derivative(a)));
}

template <template <class> class Matrix>
    requires Matrix<double>::is_Matrix
auto getvalue(const Derivative<Matrix<double>, Parameters_transformed> &x) {
    assert(x.primitive().size() == 1);
    Matrix<double> der(x.derivative()().nrows(), x.derivative()().ncols());
    for (std::size_t i = 0; i < der.size(); ++i)
        der[i] = x.derivative()()[i][0];
    return Derivative<double, Parameters_transformed>(x.primitive()[0], std::move(der),
                                              x.dx());
}

template < template <class> class Matrix>
    requires Matrix<double>::is_Matrix
Maybe_error<Derivative<Matrix<double>, Parameters_transformed>>
inv(const Derivative<Matrix<double>, Parameters_transformed> &x) {
    auto inv_x = inv(x.primitive());
    if (!inv_x)
        return inv_x.error();
    else {
        auto dinv = apply(
            [&inv_x](auto const &dx) {
                return -inv_x.value() * dx * inv_x.value();
            },
            x.derivative()());
        return Derivative<Matrix<double>, Parameters_transformed>(inv_x.value(), dinv,
                                                          x.dx());
    }
}


Maybe_error<std::tuple<Derivative<Matrix<double>, Parameters_transformed>,
                       Derivative<DiagonalMatrix<double>, Parameters_transformed>,
                       Derivative<Matrix<double>, Parameters_transformed>>>
eigs(const Derivative<Matrix<double>, Parameters_transformed> &x,
     bool does_permutations = false, bool does_diagonal_scaling = false,
     bool computes_eigenvalues_condition_numbers = false,
     bool computes_eigenvectors_condition_numbers = false) {
    
    bool do_Nelson_Normalization=true;
    auto res = eigs(x.primitive(), does_permutations, does_diagonal_scaling,
                    computes_eigenvalues_condition_numbers,
                    computes_eigenvectors_condition_numbers,do_Nelson_Normalization);
    
    if (!res)
        return res.error();
    else {
        auto [VR, lambda, VL] = std::move(res).value();
        
        auto Maybe_derlambda = apply_maybe(
            [VR, lambda, VL](auto const &dx) ->Maybe_error<DiagonalMatrix<double>>{
                auto out = DiagonalMatrix<double>(lambda.nrows(), lambda.ncols());
                for (std::size_t i = 0; i < lambda.size(); ++i) {
                    auto vT = tr(VL(":", i));
                    auto u = VR(":", i);
                    auto vT_u=getvalue(vT * u);
                    if (vT_u==0)
                        return error_message("nan value for derivative");
                    out[i] = getvalue(vT * dx * u) / getvalue(vT * u);
                }
                return out;
            },
            x.derivative()());
        if (!Maybe_derlambda)
            return Maybe_derlambda.error();
        //using test=typename decltype(Maybe_derlambda)::thet;
        
        auto dLambda = Derivative<DiagonalMatrix<double>, Parameters_transformed>(
            lambda, std::move(Maybe_derlambda.value()), x.dx()) ;
        
        auto VRRV = XTX(VR);
        
        auto Maybe_derVR = apply_maybe(
            [&VR, &lambda, &VL, &VRRV](auto const &dx)->Maybe_error<Matrix<double>>
            {
                Matrix<double> C(VR.nrows(), VR.ncols());
                for (std::size_t k = 0; k < VR.nrows(); ++k) {
                    auto uk = VR(":", k);
                    for (std::size_t j = 0; j < VR.ncols(); ++j) {
                        if (k != j) {
                            auto vTj = tr(VL(":", j));
                            auto uj = VR(":", j);
                            double dl = lambda[k] - lambda[j];
                            auto dl_vTJvuj=dl*getvalue(vTj * uj);
                            if (dl_vTJvuj==0)
                                return error_message("nan value for derivative");
                            C(k, j) = getvalue(vTj * dx * uk) / dl_vTJvuj;
                        }
                    }
                    C(k, k) = 0;
                    for (std::size_t j = 0; j < VR.ncols(); ++j) {
                        if (k != j)
                            C(k, k) -= VRRV(k, j) * C(k, j);
                    }
                }
                return tr(C) * VR;
            },
            x.derivative()());
        if (!Maybe_derVR)
            return Maybe_derVR.error();
        auto dVR = Derivative<Matrix<double>, Parameters_transformed>(VR, Maybe_derVR.value(), x.dx());
        auto dVL = inv(dVR);
        if (!dVL)
            return dVL.error();
        else
            return std::tuple(dVR, dLambda, dVL.value());
    }
}


auto Taylor_first(Derivative<double, Parameters_transformed> const &f,
                  const Parameters_transformed &x, double eps) {
    auto dx=x;
    dx()=dx()*eps;
    return primitive(f) + f.derivative() * dx;
}

template <template <class> class aMatrix>
auto Taylor_first(Derivative<aMatrix<double>, Parameters_transformed> const &f,
                  const Parameters_transformed &x, double eps) {
    auto dx=x;
    dx()=dx()*eps;
    return primitive(f) + f.derivative() * dx;
}


auto Taylor_first(Derivative<Parameters_transformed, Parameters_transformed> const &f,
                  const Parameters_transformed &x, double eps) {
    return Parameters_transformed(x.parameters(),
                          primitive(f) + f.derivative() * x * eps);
}

} // namespace var

#endif // PARAMETERS_DERIVATIVE_H
