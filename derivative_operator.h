#ifndef DERIVATIVE_OPERATOR_H
#define DERIVATIVE_OPERATOR_H

#include <concepts>
#include <utility>
#include <ostream>
#include "matrix.h"
#include "maybe_error.h"
#include "variables.h"
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

template<class>
struct untransformed_type;

template<class>
struct transformation_type;


template<class T>
using untransformed_type_t=typename untransformed_type<T>::type;

template<class T>
using transformation_type_t=typename transformation_type<T>::type;



struct Identity_Op{ template<class T> using type=T;};


template<class T>
struct untransformed_type{
    using type=T;
};


template<class T>
struct transformation_type{
    using type=Identity_Op;
};

template<class Op, class F>
using Op_t=typename Op::template type<F>;


template<class Op_on_F,class G>
using Transfer_Op_to=Op_t<transformation_type_t<Op_on_F>,G>;



template<class Tr, class T>
concept U=std::same_as<untransformed_type_t<Tr>,T>;



template<class, class>
class d_d_;

template<class T>
constexpr d_d_<T,T> self_derivative(const T&);

    

template<class, class>
class Derivative;

template<class X>
struct Derivative_Op{
    template<class F>
    using type=Derivative<F,X>;
};


struct Maybe_error_Op{
    template<class X>
    using type=Maybe_error<X>;
};

template<class X>
struct Maybe_error_Derivative_Op{
    template<class F>
    using type=Maybe_error<Derivative<F,X>>;
};



template<class T, class X>
struct untransformed_type<Derivative<T,X>>{
    using type=T;
};

template<class T>
struct untransformed_type<Maybe_error<T>>{
    using type=T;
};

template<class T, class X>
struct untransformed_type<Maybe_error<Derivative<T,X>>>{
    using type=T;
};



template<class T, class X>
struct transformation_type<Derivative<T,X>>{
    using type=Derivative_Op<X>;
};





template<>
class d_d_<double,double>
{
    double m_dydx;
    double const* ptr_dx;
public:
    using value_type=double;
    constexpr d_d_(double dydx, const double& x):m_dydx{dydx}, ptr_dx{&x}{}
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
    constexpr d_d_(){}
    constexpr auto& dx()const{return *ptr_dx;}
};



template<>
class d_d_<double,Matrix<double>>
{
    Matrix<double> m_dydx;
    Matrix<double> const* ptr_dx;
public:
    using value_type=double;
     d_d_(Matrix<double>const& dydx, const Matrix<double>& x):m_dydx{dydx}, ptr_dx{&x}{}
     d_d_(Matrix<double>&& dydx, const Matrix<double>& x):m_dydx{std::move(dydx)}, ptr_dx{&x}{}
     auto& operator()(){return m_dydx;}
     auto operator()()const{return m_dydx;}
     d_d_(){}
     auto& dx()const{return *ptr_dx;}
};





constexpr auto self_derivative(const double& x)
{
    return d_d_<double,double>(1.0,x);
}    

template<class T>
class Primitive
{
    T m_y;
public:
    template<class aT>
        requires std::is_same_v<T,std::decay_t<aT>>
    constexpr Primitive(aT&& y):m_y{std::forward<aT>(y)}{}
    constexpr auto& operator()(){return m_y;}
    constexpr auto& operator()()const{return m_y;}
};


template<class T>
class Primitive<Matrix<T>>: public Matrix<T>
{
public:
    
    template<class aT>
        requires std::is_same_v<Matrix<T>,std::decay_t<aT>>
    constexpr Primitive(aT&& y):Matrix<T>{std::forward<aT>(y)}{}
    constexpr auto& operator()(){return static_cast<Matrix<T>&>(*this);}
    constexpr auto& operator()()const {return static_cast<Matrix<T>const &>(*this);}
 };





template<class N,class D>
class Derivative;
 
 template<class N,class D>
std::ostream& operator<<(std::ostream& os, const Derivative<N,D>& d)
{
    os<<primitive(d)<<"\n derivative: \n"<<derivative(d)();
    return os;
}

template<class N,class D>
std::ostream& print(std::ostream& os, const Derivative<N,D>& d)
{
    print(os, primitive(d));
    os<<"\n derivative: \n";
    print(os,derivative(d)());
    return os;
}



namespace impl{

template<class F,class X>
struct Derivative_impl
{
    using type=std::conditional_t<F::is_constant,F,Derivative<F,X>>;  
};


}

template<class F,class X>
using Derivative_t=impl::Derivative_impl<F,X>::type;


template<>
class Derivative<double,double>{
    double m_x;
    d_d_<double,double> m_d;
public:
    
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}
    
    constexpr Derivative(double fx, double dfdx, const double& dx ):m_x{fx},
        m_d{dfdx,dx}{}
    
    constexpr Derivative(double x, const double& dx):m_x{x},
        m_d{0.0,dx}{}
    Derivative(){}
    auto& dx()const {return m_d.dx();}
    
    
    friend auto exp(const Derivative &x) {
        auto f = exp(x.primitive());
        return Derivative(f, f * x.derivative()(), x.dx());
    }
    
    friend auto max(const Derivative &x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }
    
    friend auto min(const Derivative &x, double y) {
        if (x.primitive() >= y)
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
                                 (1.0 / (x.primitive() * std::log(10))));
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
    
};




template<> class Derivative<double,Matrix<double>>{
    double m_x;
    d_d_<double,Matrix<double>> m_d;
public:
    
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}
    
     Derivative(double fx, Matrix<double>const & dfdx, const Matrix<double>& dx ):m_x{fx},
        m_d{dfdx,dx}{}
     Derivative(double fx, Matrix<double>&& dfdx, const Matrix<double>& dx ):m_x{fx},
         m_d{std::move(dfdx),dx}{}
     
     
     Derivative(double x, const Matrix<double>& dx):m_x{x},
        m_d{dx-dx,dx}{}
    Derivative(){}
    auto& dx()const {return m_d.dx();}
    
    
    friend auto exp(const Derivative &x) {
        auto f = exp(x.primitive());
        return Derivative(f, f * x.derivative()(), x.dx());
    }
    
    friend auto max(const Derivative &x, double y) {
        if (x.primitive() <= y)
            return x;
        else
            return Derivative(y, 0.0 * x.derivative()(), x.dx());
    }
    
    friend auto min(const Derivative &x, double y) {
        if (x.primitive() >= y)
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
                                 (1.0 / (x.primitive() * std::log(10))));
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
    
};






struct NoDerivative{
    template<class T>
    friend NoDerivative operator*(NoDerivative,T){return NoDerivative{};}
    template<class T>
    friend NoDerivative operator*(T,NoDerivative){return NoDerivative{};}
    
    template<class T>
    friend T operator+(T&& x,NoDerivative){return std::forward<T>(x);}
    
    template<class T>
    friend T operator-(T&& x,NoDerivative){return std::forward<T>(x);}
    
    template<class T>
    friend T operator-(NoDerivative,T&& x){return -std::forward<T>(x);}
    
    template<class T>
    friend T operator+(NoDerivative,T&& x){return std::forward<T>(x);}
    NoDerivative operator()()const {return {};}
    
};



template<class>
struct is_derivative: public std::false_type{};

template<class X,class Y>
struct is_derivative<Derivative<X,Y>>: public std::true_type{};

template<class T>
constexpr bool is_derivative_v=is_derivative<std::decay_t<T>>::value;


template<class X,class Y>
decltype(auto) primitive(const Derivative<X,Y>& d)
{
    return d.primitive();
}
template<class X,class Y>
decltype(auto) primitive(Derivative<X,Y>& d)
{
    return d.primitive();
}



template <class X>
    requires (!is_derivative_v<X>)
decltype(auto) primitive(X&& x) { return std::forward<X>(x);}

template <class X>
    requires (!is_derivative_v<X>)
auto derivative(X&& ) { return NoDerivative{};}


template<class X,class Y>
decltype(auto) derivative (const Derivative<X,Y>& d)
{
    return d.derivative();
}



template<class...>
struct dx_of_dfdx;

template<class F,class X>
struct dx_of_dfdx<Derivative<F,X>>
{
    using type=X;  
};

template<class F>
    requires(!is_derivative_v<F>)
struct dx_of_dfdx<F>
{
    using type=NoDerivative;  
};
template<>
struct dx_of_dfdx<>
{
    using type=NoDerivative;  
};


template<class F,class X, class G, class... Ds>
struct dx_of_dfdx<Derivative<F,X>,Derivative<G,X>,Ds...>
{
    using type=typename dx_of_dfdx<Derivative<F,X>,Ds...>::type;  
};

template<class F,class X, class G, class... Ds>
    requires (!is_derivative_v<G>)
struct dx_of_dfdx<Derivative<F,X>,G,Ds...>
{
    using type=typename dx_of_dfdx<Derivative<F,X>,Ds...>::type;  
};

template<class F,class X,class G, class... Ds>
    requires (!is_derivative_v<G>)
struct dx_of_dfdx<G,Derivative<F,X>, Ds...>
{
    using type=typename dx_of_dfdx<Derivative<F,X>,Ds...>::type;  
};

template<class F,class G, class... Ds>
    requires (!is_derivative_v<G>&&!is_derivative_v<F>)
struct dx_of_dfdx<G,F, Ds...>
{
    using type=typename dx_of_dfdx<Ds...>::type;  
};


template<class X,class Y>
struct dx_of_dfdx<d_d_<X,Y>>
{
    using type=Y;  
};




template<class ...T>
using    dx_of_dfdx_t=typename std::decay_t<dx_of_dfdx<T...>>::type;

inline auto get_dx_of_dfdx()
{
    return NoDerivative{};
}



template<class T,class ...Ts>
    requires (!is_derivative_v<T>)
decltype(auto) get_dx_of_dfdx(const T& , const Ts&...xs)
{
    return get_dx_of_dfdx(xs...);
}




template<class F,class X,class T,class ...Ts>
decltype(auto) get_dx_of_dfdx(const Derivative<F,X>& x, const T&,const Ts&...)
{
    return get_dx_of_dfdx(x);
}




 //auto a=Derivative<double,double>(9.9,9.9);


 //auto b=primitive(a);



template<class var,class T, class X >
    requires (var::is_variable&& std::constructible_from<var,T>)
auto build(Derivative<T,X>&& x){return Derivative<var,X>(std::move(x));}





template<class var, class... T>
    requires (std::constructible_from<var,untransformed_type_t<T>...>&&
             (is_derivative_v<T>||...||false)&&(!is_Maybe_error<T>&&...&&true))
auto build(T...x){
    using X=dx_of_dfdx_t<T...>;
    return Derivative<var,X>(std::forward<T>(x)...);}


template<class var, class... T>
    requires (std::constructible_from<var,untransformed_type_t<T>...>&&
             (is_derivative_v<T>||...||false)&&(is_Maybe_error<T>||...||false))
auto build(T...x){
    using X=dx_of_dfdx_t<T...>;
    if ((is_valid(x)&&...&&true))
        return Maybe_error<Derivative<var,X>>(std::forward<T>(x)...);
    else
        return Maybe_error<Derivative<var,X>>(error_message((get_error(x)()+...+"")));        
}

template<class var, class... T>
    requires (std::constructible_from<var,untransformed_type_t<T>...>&&
             (!is_derivative_v<T>&&...&&true)&&(is_Maybe_error<T>||...||false))
auto build(T...x){
    if ((is_valid(x)&&...&&true))
        return Maybe_error<var>(std::forward<T>(x)...);
    else
        return Maybe_error<var>(error_message((get_error(x)()+...+"")));        
}


template<template<class...>class Vector_Space, class... T>
    requires (std::constructible_from<Vector_Space<untransformed_type_t<T>...>,untransformed_type_t<T>...>&&
             (is_derivative_v<T>||...||false)&&(!is_Maybe_error<T>&&...&&true))
auto build(T...x){
    using X=dx_of_dfdx_t<T...>;
    return Derivative<Vector_Space<untransformed_type_t<T>...>,X>(std::forward<T>(x)...);
}

template<template<class...>class Vector_Space, class... T>
    requires (std::constructible_from<Vector_Space<untransformed_type_t<T>...>,untransformed_type_t<T>...>&&
             (is_derivative_v<T>||...||false)&&(is_Maybe_error<T>||...||false))
auto build(T...x){
    using X=dx_of_dfdx_t<T...>;
    if ((is_valid(x)&&...&&true))
        return Maybe_error<Derivative<Vector_Space<untransformed_type_t<T>...>,X>>(std::forward<T>(x)...);
    else
        return error_message((get_error(x)()+...+""));        
}



template<template<class...>class Vector_Space, class... T>
    requires (std::constructible_from<Vector_Space<untransformed_type_t<T>...>,untransformed_type_t<T>...>&&
             (!is_derivative_v<T>&&...&&true)&&(!is_Maybe_error<T>&&...&&true))
auto build(T...x){
    return Vector_Space<untransformed_type_t<T>...>(std::forward<T>(x)...);
}


template<template<class...>class Vector_Space, class... T>
    requires (std::constructible_from<Vector_Space<untransformed_type_t<T>...>,untransformed_type_t<T>...>&&
             (!is_derivative_v<T>&&...&&true)&&(is_Maybe_error<T>||...||true))
auto build(T...x){
    if ((is_valid(x)&&...&&true))
        return Maybe_error<Vector_Space<untransformed_type_t<T>...>>(std::forward<T>(x)...);
    else
        return error_message((get_error(x)()+...+""));        
}






template<template<class...> class varr, class Id, class F, class... T>
    requires (std::constructible_from<varr<Id,F,T...>,var::Var<Id>,F,T...>&&(std::is_same_v<varr<Id,T...>,Fun<Id,T...>>)&&
             (!is_derivative_v<T>&&...&&true)&&(!is_Maybe_error<T>&&...&&true))
varr<Id,F,std::decay_t<T>...> build(Var<Id>,F&& t_f,T&&...t_x){return varr(Var<Id>{},std::forward<F>(t_f),std::forward<T>(t_x)...);}


template<template<class...> class var, class Id, class F,class... T>
    requires (std::constructible_from<var<Id,F,untransformed_type_t<T>...>,Var<Id>,F,untransformed_type_t<T>...>&&
             (std::is_same_v<var<Id,F,untransformed_type_t<T>...>,Fun<Id,F,untransformed_type_t<T>...>>)&&
             (is_derivative_v<T>||...||false)&&(!is_Maybe_error<T>&&...&&true))
auto build(Var<Id>,F&& t_f,T&&...t_x){
    using X=dx_of_dfdx_t<std::decay_t<T>...>;
    return Derivative<Fun<Id,F,untransformed_type_t<std::decay_t<T>>...>,X>(
        Var<Id>{},std::forward<F>(t_f),std::forward<T>(t_x)...);
}




template<class T, class S>
    requires(is_derivative_v<T>||is_derivative_v<S>)
auto operator*(const T& x, const S& y)
{
    using X=dx_of_dfdx_t<T,S>;
    using F=decltype(primitive(x)*primitive(y));
    
    return Derivative<F,X>(primitive(x)*primitive(y),derivative(x)()*primitive(y)+primitive(x)*derivative(y)(), get_dx_of_dfdx(x,y));
}

template<class T, class S>
//    requires(is_derivative_v<T>||is_derivative_v<S>)
auto max(const T& x, const S& y)
{
    using std::max;
    return max(primitive(x),primitive(y));
}





template<class T, class S>
    requires(is_derivative_v<T>||is_derivative_v<S>)
auto operator+(const T& x, const S& y)
{
    using X=dx_of_dfdx_t<T,S>;
    using F=decltype(primitive(x)+primitive(y));
    
    return Derivative<F,X>(primitive(x)+primitive(y),derivative(x)()+derivative(y)(),get_dx_of_dfdx(x,y));
}

template<class T, class S>
    requires(is_derivative_v<T>&&is_derivative_v<S>)
auto operator-(const T& x, const S& y)
{
    using X=dx_of_dfdx_t<T,S>;
    using F=decltype(primitive(x)-primitive(y));
    
    return Derivative<F,X>(primitive(x)-primitive(y),derivative(x)()-derivative(y)(),get_dx_of_dfdx(x,y));
}

template<class T, class S>
    requires(is_derivative_v<T>&&(!is_derivative_v<S>))
auto operator-(const T& x, const S& y)
{
    using X=dx_of_dfdx_t<T,S>;
    using F=decltype(primitive(x)-y);
    
    return Derivative<F,X>(primitive(x)-y,derivative(x)(),x.dx());
}

template<class T, class S>
    requires(!is_derivative_v<T>&&is_derivative_v<S>)
auto operator-(const T& x, const S& y)
{
    using X=dx_of_dfdx_t<T,S>;
    using F=decltype(x-primitive(y));
    
    return Derivative<F,X>(x-primitive(y),derivative(y)()* (-1.0),y.dx());
}





auto max(is_Container auto const& c)
{
    auto out=c[0];
    for (std::size_t i=1; i<c.size(); ++i)
        if (std::isnan(primitive(c[i]))||primitive(out)<primitive(c[i]))
            out=c[i];
    return out; 
}




}



#endif // DERIVATIVE_OPERATOR_H
