#ifndef PARAMETERS_DERIVATIVE_H
#define PARAMETERS_DERIVATIVE_H

#include "parameters.h"
#include "matrix_derivative.h"

namespace var {

template<class Id>
class d_d_<Parameters<Id>,double>
{
    d_d_<Matrix<double>,double> m_dydx;
public:
    using value_type=d_d_<Matrix<double>,double>;
    template<class aMatrix>
        requires std::constructible_from<d_d_<Matrix<double>,double>,aMatrix>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};

template<class Id>
class d_d_<double,Parameters<Id>>
{
    Matrix<double> m_dydx;
public:
    using value_type=Matrix<double>;
    template<class aMatrix>
        requires std::is_same_v<Matrix<double>,std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
     auto& operator()(){return m_dydx;}
     auto& operator()()const{return m_dydx;}
};

template<class Id>
class d_d_<Matrix<double>,Parameters<Id>>
{
    Matrix<Matrix<double>> m_dydx;
public:
    using value_type=Matrix<Matrix<double>>;
    
    template<class aMatrix>
        requires std::is_same_v<Matrix<Matrix<double>>,std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
    
 };




template<class Id,class Id2>
class d_d_<Parameters<Id>,Parameters<Id2>>
{
    d_d_<Matrix<double>,Parameters<Id2>> m_dydx;
public:
    using value_type=d_d_<Matrix<double>,Parameters<Id2>>;
    
    template<class aMatrix>
        requires std::constructible_from<d_d_<Matrix<double>,Parameters<Id2>>,aMatrix>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};

template<class Id>
d_d_<Parameters<Id>,Parameters<Id>> selfDerivative(const Parameters<Id>& x)
{
    return d_d_<Matrix<double>,Parameters<Id>>(selfDerivative(x())());
}


template<class Id>
class Derivative<double,Parameters<Id>>: public Primitive<double>, public d_d_<double,Parameters<Id>>
{
public:
    template<class aMatrix>
        requires std::constructible_from<d_d_<double,Parameters<Id>>,std::decay_t<aMatrix>>
    Derivative(double x, aMatrix&& dydx): Primitive<double>{x}, d_d_<double,Parameters<Id>>{std::forward<aMatrix>(dydx)}{}    
};


template<class Id>
class Derivative<Matrix<double>,Parameters<Id>>: public Primitive<Matrix<double>>, public d_d_<Matrix<double>,Parameters<Id>>
{
public:
    
    
    template<class aMatrix, class dP>
        requires (std::constructible_from<Matrix<double>,aMatrix>&&
                    std::constructible_from<d_d_<Matrix<double>,Parameters<Id>>,dP>)
    Derivative(aMatrix&& x, dP&& dydx): Primitive<Matrix<double>>{std::forward<aMatrix>(x)}, d_d_<Matrix<double>,Parameters<Id>>{std::forward<dP>(dydx)}{}    
    
    
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,Parameters<Id>>(
            primitive(*this)(i,j),
            applyMap([i,j](auto const& m){return m(i,j);},derivative(*this)()));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,Parameters<Id>>(
            primitive(*this)[i],
            applyMap([i](auto const& m){return m[i];},derivative(*this)));
    }
     
};

template<class Id, class T>
class Derivative<Parameters<Id>,T>{
    Derivative<Matrix<double>,T> m_x;  
    
public:
    
    template<class dParam>
        requires (std::constructible_from<Derivative<Matrix<double>,T>,dParam>)
    Derivative(dParam&& x):m_x{std::forward<dParam>(x)}{}
    
    
    template<class aParam>
        requires (std::constructible_from<Parameters<Id>,aParam>)
    Derivative(aParam&& x): m_x{std::forward<aParam>(x)(), selfDerivative(x)()}{}
    
    auto& operator()()const
    {
        return m_x;
    }
};


    
}



#endif // PARAMETERS_DERIVATIVE_H
