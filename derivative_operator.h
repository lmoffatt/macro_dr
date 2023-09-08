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

template<class T>
auto elementary_type(T)
{
    return T{};
}

template<class T>
auto elementary_type(T)-> typename T::value_type
{
    return elementary_type(typename T::value_type{});
}


template<class>
class Primitive;



template<class, class>
class d_d_;


template<class, class>
class Derivative;



template<>
class d_d_<double,double>
{
    double m_dydx;
public:
    constexpr d_d_(double dydx):m_dydx{dydx}{}
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};

template<>
class Primitive<double>
{
    double m_y;
public:
    constexpr Primitive(double y):m_y{y}{}
    constexpr auto& operator()(){return m_y;}
    constexpr auto operator()()const{return m_y;}
};


template<>
class Derivative<double,double>: public Primitive<double>, d_d_<double,double>{
public:
    constexpr Derivative(double x, double dx):Primitive<double>{x}, d_d_<double,double>{dx}{}
};


template<class X,class Y>
constexpr decltype(auto) primitive(const Derivative<X,Y>& d)
{
    return static_cast<Primitive<X> const&>(d)();
}


constexpr auto a=Derivative<double,double>(9.9,9);


constexpr auto b=primitive(a);






 




}



#endif // DERIVATIVE_OPERATOR_H
