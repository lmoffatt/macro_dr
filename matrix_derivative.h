#ifndef MATRIX_DERIVATIVE_H
#define MATRIX_DERIVATIVE_H

#include "derivative_operator.h"
#include "matrix.h"
#include "variables.h"
namespace var {

template<class M>
struct M_der
{
    using type=M;
};

template<>
struct M_der<SymPosDefMatrix<double>>{
    using type=SymmetricMatrix<double>;
};

template<>
struct M_der<DiagPosDetMatrix<double>>{
    using type=DiagonalMatrix<double>;
};

template <class M>
using M_der_t=typename M_der<M>::type;



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
d_d_<Matrix<double>,Matrix<double>> self_derivative(const Matrix<double>& x)
{
    auto out=Matrix<Matrix<double>>(x.nrows(),x.ncols(),Matrix<double>(x.nrows(),x.ncols(),0.0));
    for (std::size_t i=0; i<x.size(); i++)
    {
        out[i][i]=1.0;
    }
    return out;
}

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
            ddx(i,j)=derivative(Dx)()[ii];
            ++it_ind;
            ++it_val;
        }
        return ddx;},df_dX0());
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




template <class T, class X> auto diag(const Derivative<Matrix<T>,X> &a) {
    
    return Derivative<DiagonalMatrix<T>,X>(diag(primitive(a)),
                                            apply_par([](auto const & d){return diag(d);}, derivative(a)));
}


template <class T, class X> auto diagpos(const Derivative<Matrix<T>,X> &a) {
    
    
    
    
    return Derivative<DiagPosDetMatrix<T>,X>(diagpos(primitive(a)),
                                              apply_par([](auto const & d){
                                                  return diag(d);
                                              }, derivative(a)));
}


template<class X>
auto XTX(const Derivative<Matrix<double>,X> &a) {
    
    auto& f=primitive(a);
    return Derivative<SymPosDefMatrix<double>,X>(XTX(f),
                                                   apply_par([&f](auto const & d){return X_plus_XT(tr(d)*f);}, derivative(a)));
}


template<class X>
auto X_plus_XT(const Derivative<Matrix<double>,X> &a) {
    
    auto& f=primitive(a);
    return Derivative<SymmetricMatrix<double>,X>(X_plus_XT(f),
                                                  apply_par([](auto const & d){return X_plus_XT(d);}, derivative(a)));
}


template<class X,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto tr(const Derivative<aMatrix<double>,X> &a) {
    
    auto& f=primitive(a);
    return Derivative<aMatrix<double>,X>(tr(f),
                                                  apply_par([](auto const & d){return tr(d);}, derivative(a)));
}


template<class X,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDiv(const Derivative<aMatrix<double>,X> &a,const Derivative<aMatrix<double>,X> &b) {
    
    return zip([] (auto& x, auto& y){ return x/y;}, a,b);
}

template<class X,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemMult(const Derivative<aMatrix<double>,X> &a,const Derivative<aMatrix<double>,X> &b) {
    
    return zip([] (auto& x, auto& y){ return x*y;}, a,b);
}

template<class X,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>,X> &a, const Derivative<bMatrix<double>,X> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<S,X>(TranspMult(fa,fb),
                                                  zip_par([&fa,&fb](auto const & da, auto const & db)
                                                  {return TranspMult(fa,db)+TranspMult(da,fb);}, derivative(a), derivative(b)));
}


template<class X,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>,X> &a, const bMatrix<double> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<S,X>(TranspMult(fa,fb),
                            apply_par([&fb](auto const & da)
                                    {return TranspMult(da,fb);}, derivative(a)));
}

template<class X,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const aMatrix<double>&a, const Derivative<bMatrix<double>,X> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=a;
    auto& fb=primitive(b);
    
    return Derivative<S,X>(TranspMult(fa,fb),
                            apply_par([&fa](auto const & db)
                                    {return TranspMult(fa,db);}, derivative(b)));
}


template <class X>
Derivative<SymPosDefMatrix<double>,X> AT_B_A(const Derivative<Matrix<double>,X> &a, const Derivative<SymmetricMatrix<double>,X> &b) {
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<SymPosDefMatrix<double>,X>(AT_B_A(fa,fb),
                            zip_par([&fa,&fb](auto const & da, auto const & db)
                                    {return AT_B_A(fa,db)+X_plus_XT(TranspMult(da,fb*fa));}, derivative(a), derivative(b)));
}



template <class X>
Derivative<SymPosDefMatrix<double>,X> AT_B_A(const Derivative<Matrix<double>,X> &a, const SymmetricMatrix<double> &b) {
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<SymPosDefMatrix<double>,X>(AT_B_A(fa,fb),
                                                  apply_par([&fa,&fb](auto const & da)
                                                          {return X_plus_XT(TranspMult(da,fb*fa));}, derivative(a)));
}




template<template<class>class Matrix, class X>
    requires Matrix<double>::is_Matrix
double getvalue(const Derivative<Matrix<double>,X> &x) {
    assert(x.primitive().size() == 1);
    Matrix<double> der(x.derivative()().nrows(),x.derivative()().ncols());
    for (std::size_t i=0; i<der.size(); ++i) der[i]=x.derivative()()[i][0];
    return Derivative<double,X>(x[0],std::move(der)); 
}




template<class X, template<class>class Matrix>
    requires Matrix<double>::is_Matrix
Maybe_error<Derivative<Matrix<double>,X>> inv(const Derivative<Matrix<double>,X>& x)
{
    auto inv_x=inv(x.primitive());
    if (!inv_x)
        return inv_x.error();
    else
    {
        auto dinv=apply([&inv_x](auto const& dx){ return -inv_x.value()*dx*inv_x.value();},x.derivative()());
        return Derivative<Matrix<double>,X>(inv_x.value(),dinv);
    }
}

template<class X>
Maybe_error<std::tuple<Derivative<Matrix<double>,X>,Derivative<DiagonalMatrix<double>,X>,Derivative<Matrix<double>,X>>>
    eigs(const Derivative<Matrix<double>,X> &x, bool does_permutations = true,
          bool does_diagonal_scaling = true,
          bool computes_eigenvalues_condition_numbers = false,
          bool computes_eigenvectors_condition_numbers = false) {
    
    auto res=eigs(x.primitive(),does_permutations,does_diagonal_scaling,
                    computes_eigenvalues_condition_numbers,computes_eigenvectors_condition_numbers);
    
    if (!res)
        return res.error();
    else{
    auto [VR, lambda, VL] = std::move(res).value();
    
    auto derlambda=apply([&VR,&lambda,&VL](auto const & dx){
        auto out=DiagonalMatrix<double>(lambda.nrows(),lambda.ncols());
        for (std::size_t i = 0; i < lambda.size(); ++i) {
            auto vT = VL(i, ":");
            auto u = VR(":", i);
            out[i]= getvalue(vT * dx * u);
        }
        return out;}, x.derivative()());
    
    auto dLambda=Derivative<DiagonalMatrix<double>,X>(lambda,derlambda);
    
    auto derVR=apply([&VR,&lambda,&VL](auto const & dx){
        
        Matrix<double> C(VR.nrows(), VR.ncols());
        for (std::size_t k = 0; k < VR.nrows(); ++k) {
            auto uk = VR(":", k);
            std::size_t m = 0;
            for (std::size_t j = 0; j < VR.ncols(); ++j) {
                if (uk[j]== 1)
                    m = j;
                if (k != j) {
                    auto vTj = VL(j, ":");
                    double dl = lambda[k] - lambda[j];
                    C(k, j) = getvalue(vTj * dx * uk) / dl;
                }
            }
            C(k, k) = 0;
            for (std::size_t j = 0; j < VR.ncols(); ++j) {
                if (k != j)
                    C(k, k) -= VR(m, j) * C(k, j);
            }
        }
        return tr(VR)*C;
    }, x.derivative()());
    auto dVR=Derivative<Matrix<double>,X>(VR,derVR);
    auto dVL=inv(dVR);
    if (! dVL)
        return dVL.error();
    else
        return std::tuple(dVR,dLambda,dVL.value());
    }
}







}






#endif // MATRIX_DERIVATIVE_H
