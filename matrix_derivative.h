#ifndef MATRIX_DERIVATIVE_H
#define MATRIX_DERIVATIVE_H

#include "derivative_operator.h"
#include "matrix.h"
#include "variables.h"
namespace var {

template<>
class d_d_<Matrix<double>,double>
{
    Matrix<double> m_dydx;
public:
    using value_type=Matrix<double>;
    template<class aMatrix>
        requires std::is_same_v<Matrix<double>,std::decay_t<aMatrix>>
     d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
     auto& operator()(){return m_dydx;}
     auto& operator()()const{return m_dydx;}
};

template<class T>
class d_d_<T,Matrix<double>>
{
    Matrix<T> m_dydx;
public:
    using value_type=Matrix<T>;
    
    template<class aMatrix>
        requires std::is_same_v<Matrix<T>,std::decay_t<aMatrix>>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};


d_d_<Matrix<double>,Matrix<double>> selfDerivative(const Matrix<double>& x)
{
    auto out=Matrix<Matrix<double>>(x.nrows(),x.ncols(),Matrix<double>(x.nrows(),x.ncols(),0.0));
    for (std::size_t i=0; i<x.size(); i++)
    {
        out[i][i]=1.0;
    }
    return out;
}

template<>
class Derivative<Matrix<double>,double>: public Primitive<Matrix<double>>, public d_d_<Matrix<double>,double>
{
public:
    
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,double>(primitive(*this)(i,j),derivative(*this)(i,j));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,double>(primitive(*this)[i],derivative(*this)[i]);
    }    
    
    
};

template<>
class Derivative<double,Matrix<double>>: public Primitive<double>, public d_d_<double,Matrix<double>>
{
    
    
};



template<>
class Derivative<Matrix<double>,Matrix<double>>: public Primitive<Matrix<double>>, public d_d_<Matrix<double>,Matrix<double>>
{
public:
         
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,Matrix<double>>(primitive(*this)(i,j),applyMap([i,j](auto const& m){return m(i,j);},derivative(*this)));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,Matrix<double>>(primitive(*this)[i],applyMap([i](auto const& m){return m[i];},derivative(*this)));
    }    
    
    
};




template<class Matrix>
    requires (std::constructible_from<Matrix,std::size_t,std::size_t,std::vector<std::tuple<std::size_t,std::size_t,double>>>)
auto build(std::size_t nrows,std::size_t ncols, std::vector<std::tuple<std::size_t,std::size_t,double>>&& values)
{
    return Matrix(nrows,ncols,std::move(values));  
}

template<class Matrix, class T>
    requires std::is_same_v<element_of_Matrix_t<Matrix>,T>
auto build_(std::size_t nrows,std::size_t ncols, std::initializer_list<std::pair<std::size_t,std::size_t>> indexes,
            std::initializer_list<T> values)
{
    Matrix x(nrows,ncols,true);
    auto it_ind=indexes.begin();
    auto it_val=values.begin();
    for (std::size_t k=0; k<values.size(); ++k)
    {
        auto i=it_ind->first;
        auto j=it_ind->second;
        auto value=*it_val;
        x(i,j)=value;
        ++it_ind;
        ++it_val;
    }
    return x;  
}



template<class Matrix, class Der>
    requires(is_derivative_v<Der>)
auto build_(std::size_t nrows,std::size_t ncols, std::initializer_list<std::pair<std::size_t,std::size_t>> indexes,
            std::initializer_list<Der> values)
{
    using X=dx_of_dfdx_t<Der>;
    using T=element_of_Matrix_t<Matrix>;
    
    
    auto df_dX0=derivative(*values.begin());
    auto n=values.size();
    auto dx=applyMap_i([&values,&indexes,n,nrows,ncols](auto e, std::size_t ii){
        auto it_ind=indexes.begin();
        auto it_val=values.begin();
        Matrix ddx(nrows,ncols,true);
        for (std::size_t k=0; k<n; ++k)
        {
            auto i=it_ind->first;
            auto j=it_ind->second;
            auto& Dx=*it_val;
            ddx(i,j)=derivative(Dx)[ii];
            ++it_ind;
            ++it_val;
        }
        return ddx;},df_dX0);
    Matrix x(nrows,ncols,true);
    auto it_ind=indexes.begin();
    auto it_val=values.begin();
    
    for (std::size_t k=0; k<values.size(); ++k)
    {
        auto i=it_ind->first;
        auto j=it_ind->second;
        auto& Dx=*it_val;
        x(i,j)=primitive(Dx);
        ++it_ind;
        ++it_val;
    }
    
    return Derivative<Matrix,X>(x,dx);  
}


}

#endif // MATRIX_DERIVATIVE_H
