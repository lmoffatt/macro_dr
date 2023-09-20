#ifndef MATRIX_DERIVATIVE_H
#define MATRIX_DERIVATIVE_H

#include "derivative_operator.h"
#include "matrix.h"
#include "variables.h"
namespace var {



template<template<class>class T_Matrix>
    requires (T_Matrix<double>::is_Matrix)
class d_d_<T_Matrix<double>,double>
{
    M_der_t<T_Matrix<double>> m_dydx;
public:
    using value_type=T_Matrix<double>;
    template<class aMatrix>
        requires std::constructible_from<M_der_t<T_Matrix<double>>,aMatrix>
    d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    constexpr d_d_(){}
    
    auto& operator()(){return m_dydx;}
    auto& operator()()const{return m_dydx;}
};

template<template<class>class Matrix, class T>
    requires (Matrix<double>::is_Matrix)
class d_d_<T,Matrix<double>>
{
    Matrix<M_der_t<T>> m_dydx;
public:
    using value_type=Matrix<T>;
    
    template<class aMatrix>
        requires std::constructible_from<Matrix<M_der_t<T>>,aMatrix>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    constexpr d_d_(){}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};  



template<template<class>class Matrix>
    requires (Matrix<double>::is_Matrix)
class Derivative<Matrix<double>,double>//: public Primitive<Matrix<double>>, public d_d_<Matrix<double>,double>
{
    using primitive_type=Matrix<double>;
    using derivative_type=d_d_<Matrix<double>,double>;
    primitive_type m_x;
    derivative_type m_d;
public:
    template <class P, class D>
        requires (std::constructible_from<primitive_type,P> &&
                 std::constructible_from<derivative_type,D>)
    Derivative(P&& t_x,D&& t_d): m_x{std::forward<P>(t_x)},m_d{std::forward<D>(t_d)}{}
    Derivative(){}
    
    auto& primitive() {return m_x;}
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}    
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,double>(primitive()(i,j),derivative()(i,j));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,double>(primitive()[i],derivative()[i]);
    }
    
    auto ncols()const {return m_x.ncols();}
    auto nrows()const {return m_x.nrows();}
    auto size()const {return m_x.size();}
    
    
};

template<template<class>class Matrix>
    requires (Matrix<double>::is_Matrix)
class Derivative<double,Matrix<double>>//: public Primitive<double>, public d_d_<double,Matrix<double>>
{
    using primitive_type=double;
    using derivative_type=d_d_<double,Matrix<double>>;
    primitive_type m_x;
    derivative_type m_d;
public:
    template <class D>
        requires (std::constructible_from<derivative_type,D>)
    Derivative(double t_x,D&& t_d): m_x{t_x},m_d{std::forward<D>(t_d)}{}
    Derivative(){}
    
    auto& primitive(){return m_x;}
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}    
    
    
};



template<template<class>class Matrix>
    requires (Matrix<double>::is_Matrix)
class Derivative<Matrix<double>,Matrix<double>>//: public Primitive<Matrix<double>>, public d_d_<Matrix<double>,Matrix<double>>
{
    using primitive_type=Matrix<double>;
    using derivative_type=d_d_<Matrix<double>,Matrix<double>>;
    primitive_type m_x;
    derivative_type m_d;
public:
    auto ncols()const {return m_x.ncols();}
    auto nrows()const {return m_x.nrows();}
    auto size()const {return m_x.size();}
    Derivative(){}
    
    template <class P, class D>
        requires    (std::constructible_from<primitive_type,P> &&
                 std::constructible_from<derivative_type,D>)
    
    Derivative(P&& t_x,D&& t_d): m_x{std::forward<P>(t_x)},m_d{std::forward<D>(t_d)}{}
    
    auto& primitive() {return m_x;}
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}    
    
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,Matrix<double>>(primitive()(i,j),applyMap([i,j](auto const& m){return m(i,j);},derivative()));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,Matrix<double>>(primitive()[i],applyMap([i](auto const& m){return m[i];},derivative()));
    }    
    
    
};




template<class Matrix>
    requires (std::constructible_from<Matrix,std::size_t,std::size_t,std::vector<std::tuple<std::size_t,std::size_t,double>>>)
auto build(std::size_t nrows,std::size_t ncols, std::vector<std::tuple<std::size_t,std::size_t,double>>&& values)
{
    return Matrix(nrows,ncols,std::move(values));  
}






}






#endif // MATRIX_DERIVATIVE_H
