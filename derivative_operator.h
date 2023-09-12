#ifndef DERIVATIVE_OPERATOR_H
#define DERIVATIVE_OPERATOR_H

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

template<class T>
using untransformed_type_t=typename untransformed_type<T>::type;


template<class T>
struct untransformed_type{
    using type=T;
};




template<class, class>
class d_d_;

template<class T>
constexpr d_d_<T,T> selfDerivative(const T&);

    

template<class, class>
class Derivative;

template<class T, class X>
struct untransformed_type<Derivative<T,X>>{
    using type=T;
};




template<>
class d_d_<double,double>
{
    double m_dydx;
public:
    using value_type=double;
    constexpr d_d_(double dydx):m_dydx{dydx}{}
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};

constexpr auto selfDerivative(const double&)
{
    return d_d_<double,double>(1.0);
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


template<class N,class D>
class Derivative;

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
class Derivative<double,double>: public Primitive<double>, d_d_<double,double>{
public:
    constexpr Derivative(double x, double dx):Primitive<double>{x},
        d_d_<double,double>{dx}{}
    
    
    constexpr Derivative(double x):Primitive<double>{x},
        d_d_<double,double>{selfDerivative(x)}{}
    
};



template <class T>
Derivative(T&&) ->Derivative<std::decay_t<T>,std::decay_t<T>>;



template<class X,class Y>
    requires (std::is_base_of_v<Primitive<X>,Derivative<X,Y>>)
decltype(auto) primitive (const Derivative<X,Y>& d)
{
    return static_cast<Primitive<X> const&>(d)();
}

template<class X,class Y>
    requires (std::is_base_of_v<d_d_<X,Y>,Derivative<X,Y>>)
decltype(auto) derivative (const Derivative<X,Y>& d)
{
    return static_cast<d_d_<X,Y> const&>(d)();
}

template<class>
struct is_derivative: public std::false_type{};

template<class X,class Y>
struct is_derivative<Derivative<X,Y>>: public std::true_type{};

template<class T>
constexpr bool is_derivative_v=is_derivative<T>::value;



struct NoDerivative{};

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



 //auto a=Derivative<double,double>(9.9,9.9);


 //auto b=primitive(a);



template<class var,class T, class X >
    requires (var::is_variable&& std::constructible_from<var,T>)
auto build(Derivative<T,X>&& x){return Derivative<var,X>(std::move(x));}





template<class var, class... T>
    requires (std::constructible_from<var,untransformed_type_t<T>...>&&
             (is_derivative_v<T>||...||false))
auto build(T...x){
    using X=dx_of_dfdx_t<T...>;
    return Derivative<var,X>(std::forward<T>(x)...);}
}



#endif // DERIVATIVE_OPERATOR_H
