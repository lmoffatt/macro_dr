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
    constexpr d_d_(){}
    
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
    d_d_(){}
    
    auto& operator()(){return m_dydx;}
    auto& operator()()const{return m_dydx;}
};

template<class Id, template <class> class aMatrix>
    requires (aMatrix<double>::is_Matrix)
class d_d_<aMatrix<double>,Parameters<Id>>
{
    Matrix<aMatrix<double>> m_dydx;
public:
    using value_type=Matrix<aMatrix<double>>;
    d_d_(){}
    
    template<class aaMatrix>
        requires std::constructible_from<Matrix<aMatrix<double>>,std::decay_t<aaMatrix>>
    constexpr d_d_(aaMatrix&& dydx):m_dydx{std::forward<aaMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
    
    template <class F>
    friend auto apply_par(F &&f, d_d_ const &a) {
        using S=std::decay_t<std::invoke_result_t<F,aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(a()[i]);
        return d_d_<S,Parameters<Id>>(x);
    }
    
    template <class F, template<class>class bMatrix>
    friend auto zip_par(F &&f, d_d_ const &a,d_d_<bMatrix<double>,Parameters<Id>>const &b) {
        using S=std::decay_t<std::invoke_result_t<F,aMatrix<double>,bMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(a()[i],b()[i]);
        return d_d_<S,Parameters<Id>>(x);
    }
    
    
    
};

template<class Id,template<class> class notSymmetricMatrix>
    requires ( (notSymmetricMatrix<double>::is_Matrix) && (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,notSymmetricMatrix<double>>)) )

auto inside_out(const d_d_<notSymmetricMatrix<double>,Parameters<Id>>& x)
{
    if (x().size()==0)
        return notSymmetricMatrix<d_d_<double,Parameters<Id>>>{};
    else
    {
        notSymmetricMatrix<d_d_<double,Parameters<Id>>> out(x()[0].nrows(),x()[0].ncols());
        for (std::size_t i=0; i<out.size(); ++i)
        {
            Matrix<double> d_d_par(x().nrows(),x().ncols());
            for (std::size_t j=0; j<d_d_par.size(); ++j)
                d_d_par[j]=x()[j][i];
            out[i]=d_d_<double, Parameters<Id>>(std::move(d_d_par));
        }
        return out;
    }
}

template<class Id,template<class> class aSymmetricMatrix>
    requires (std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,aSymmetricMatrix<double>>)
auto inside_out(const d_d_<aSymmetricMatrix<double>,Parameters<Id>>& x)
{
    aSymmetricMatrix<d_d_<double,Parameters<Id>>> out(x()[0].nrows(),x()[0].ncols());
    for (std::size_t i=0; i<out.nrows(); ++i)
        for (std::size_t j=0; j<=i; ++j)
        {
            Matrix<double> d_d_par(x().nrows(),x().ncols());
            for (std::size_t k=0; k<d_d_par.size(); ++k)
                d_d_par[k]=x()[k](i,j);
            out.set(i,j,d_d_<double, Parameters<Id>>(std::move(d_d_par)));
            
        }
    return out;
}









template<class Id,class Id2>
class d_d_<Parameters<Id>,Parameters<Id2>>
{
    d_d_<Matrix<double>,Parameters<Id2>> m_dydx;
public:
    using value_type=d_d_<Matrix<double>,Parameters<Id2>>;
    d_d_(){}
    
    template<class aMatrix>
        requires std::constructible_from<d_d_<Matrix<double>,Parameters<Id2>>,aMatrix>
    constexpr d_d_(aMatrix&& dydx):m_dydx{std::forward<aMatrix>(dydx)}{}
    
    constexpr auto& operator()(){return m_dydx;}
    constexpr auto operator()()const{return m_dydx;}
};

template<class Id>
Derivative<Parameters<Id>,Parameters<Id>> selfDerivative(const Parameters<Id>& x)
{
    
    return Derivative<Parameters<Id>,Parameters<Id>>(Derivative<Matrix<double>,Parameters<Id>>(x(),self_derivative(x())()));
}


template<class Id>
class Derivative<double,Parameters<Id>>//: public Primitive<double>, public d_d_<double,Parameters<Id>>
{
    using primitive_type=double;
    using d_type=Parameters<Id>;
    using derivative_type=d_d_<primitive_type,d_type>;
    primitive_type m_x;
    derivative_type m_d;
public:
    template <class P, class D>
        requires (std::constructible_from<primitive_type,P> &&
                 std::constructible_from<derivative_type,D>)
    
    Derivative(P t_x,D&& t_d): m_x{t_x},m_d{std::forward<D>(t_d)}{}
    Derivative(){}
    
    Derivative(double t_x): m_x{t_x},m_d{}{}
    
    operator double()const {return m_x;}
    auto& primitive() {return m_x;}
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}    
};


template<class Id, template <class> class T_Matrix>
    requires (T_Matrix<double>::is_Matrix)
class Derivative<T_Matrix<double>,Parameters<Id>>{//: public T_Matrix<double>{
    using primitive_type=T_Matrix<double>;
    using d_type=Parameters<Id>;
    using derivative_type=d_d_<M_der_t<primitive_type>,d_type>;
    primitive_type m_x;
    derivative_type m_d;
public:
    auto ncols()const {return m_x.ncols();}//{return primitive_type::ncols();}
    auto nrows()const{return m_x.nrows();}// {return primitive_type::nrows();}
    auto size()const {return m_x.size();}//{return primitive_type::size();}
    Derivative(){}
    
    //   using gserg=typename primitive_type::sgr;
    //   using gserug=typename derivative_type::sgrdd;
    
    Derivative(std::size_t nrows, std::size_t ncols):m_x(nrows,ncols){}
    
    
    template<class F>
    friend auto apply(F&& f, const Derivative& x)
    {
        return outside_in(apply(std::forward<F>(f),inside_out(x)));
    }
    
    template<class F>
    friend auto zip(F&& f, const Derivative& x,const Derivative& y)
    {
        
        return outside_in(zip(std::forward<F>(f),inside_out(x),inside_out(y)));
    }
    
    template <class P, class D>
        requires (std::constructible_from<primitive_type,P> &&
                 std::constructible_from<derivative_type,D>)
    Derivative(P&& t_x,D&& t_d): m_x{std::forward<P>(t_x)},m_d{std::forward<D>(t_d)}{}
    
    
    
    
    template <class anotherCompatibleMatrix>
    Derivative(Derivative<anotherCompatibleMatrix,d_type>const & t_x): m_x{t_x.primitive()},
        m_d{t_x.derivative()()}{}
    
    
    template <class P>
        requires (std::constructible_from<primitive_type,P>)
    Derivative(P&& t_x):  m_x{std::forward<P>(t_x)},m_d{}{}
    
    auto& primitive(){return m_x;}// {return static_cast<primitive_type&>(*this);}
    auto& primitive() const {return m_x;}//{return static_cast<primitive_type const&>(*this);}
    auto& derivative()const {return m_d;}    
    
    
    auto operator()(std::size_t i, std::size_t j)const
    {
        return Derivative<double,Parameters<Id>>(
            primitive()(i,j),
            applyMap([i,j](auto const& m){return m(i,j);},derivative()()));
    }
    
    auto operator[](std::size_t i)const
    {
        return Derivative<double,Parameters<Id>>(
            primitive()[i],
            applyMap([i](auto const& m){return m[i];},derivative()()));
    }
    
};

template<class Id,template<class> class notSymmetricMatrix>
    requires ( (notSymmetricMatrix<double>::is_Matrix) && (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,notSymmetricMatrix<double>>)) )

auto inside_out(const Derivative<notSymmetricMatrix<double>,Parameters<Id>>& x)
{
    notSymmetricMatrix<Derivative<double,Parameters<Id>>> out(x.nrows(),x.ncols());
    
    auto der=inside_out(x.derivative());
    if (der.size()>0)
        for (std::size_t i=0; i<out.size(); ++i)
            out[i]=Derivative<double, Parameters<Id>>(x.primitive()[i], der[i]);
    return out;
}




template<class Id,template<class> class aSymmetricMatrix>
    requires (std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,aSymmetricMatrix<double>>)
auto inside_out(const Derivative<aSymmetricMatrix<double>,Parameters<Id>>& x)
{
    aSymmetricMatrix<Derivative<double,Parameters<Id>>> out(x.nrows(),x.ncols());
    auto der=inside_out(x.derivative());
    for (std::size_t i=0; i<out.nrows(); ++i)
        for (std::size_t j=0; j<=i; ++j)
            out.set(i,j,Derivative<double, Parameters<Id>>(x.primitive()(i,j), der(i,j)));
    return out;
}






template<class Id,template<class> class notSymmetricMatrix>
    requires ( (notSymmetricMatrix<double>::is_Matrix) && (!(std::is_same_v<SymmetricMatrix<double>, notSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,notSymmetricMatrix<double>>)) )

auto outside_in(const notSymmetricMatrix<Derivative<double,Parameters<Id>>>& x)
{
    auto prim=notSymmetricMatrix<double>(x.nrows(),x.ncols());
    for (std::size_t i=0; i<prim.size(); ++i)
        prim[i]=x[i].primitive();
    
    auto der=Matrix<notSymmetricMatrix<double>>(
        x[0].derivative()().nrows(),x[0].derivative()().ncols(),
        notSymmetricMatrix<double>(x.nrows(),x.ncols()) );
    for (std::size_t i=0; i<der.size(); ++i)
        for (std::size_t j=0; j<der[i].size(); ++j)
            der[i][j]=x[j].derivative()()[i];
    
    return Derivative<notSymmetricMatrix<double>,Parameters<Id>>(std::move(prim),std::move(der));
}




template<class Id,template<class> class aSymmetricMatrix>
    requires (std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,aSymmetricMatrix<double>>)
auto outside_in(const aSymmetricMatrix<Derivative<double,Parameters<Id>>>& x)
{
    auto prim=aSymmetricMatrix<double>(x.nrows(),x.ncols());
    for (std::size_t i=0; i<prim.nrows(); ++i)
        for (std::size_t j=0; j<=i; ++j)
            prim.set(i,j,x(i,j).primitive());
    
    auto der=Matrix<aSymmetricMatrix<double>>(
        x[0].derivative()().nrows(),x[0].derivative()().ncols(),
        aSymmetricMatrix<double>(x.nrows(),x.ncols()) );
    for (std::size_t i=0; i<der.size(); ++i)
        for (std::size_t j1=0; j1<der[i].nrows(); ++j1)
            for (std::size_t j2=0; j2<=j1; ++j2)
                der[i].set(j1,j2,x(j1,j2).derivative()()[i]);
    
    return Derivative<aSymmetricMatrix<double>,Parameters<Id>>(std::move(prim),std::move(der));
}



template<class Id,template<class> class aSymmetricMatrix>
    requires (std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,aSymmetricMatrix<double>>)
void set( Derivative<aSymmetricMatrix<double>,Parameters<Id>>& x, std::size_t i, std::size_t j , double value)
{
    x.primitive().set(i,j,value);
    
}

template<class Id,template<class> class aSymmetricMatrix>
    requires (std::is_same_v<SymmetricMatrix<double>, aSymmetricMatrix<double>>||std::is_same_v<SymPosDefMatrix<double>,aSymmetricMatrix<double>>)
void set( Derivative<aSymmetricMatrix<double>,Parameters<Id>>& x, std::size_t i, std::size_t j , const Derivative<double, Parameters<Id>>& value)
{
    x.primitive().set(i,j,value.primitive());
    if (x.derivative()().empty())
        x.derivative()()=Matrix<aSymmetricMatrix<double>>(value.derivative()().nrows(),value.derivative()().ncols(),
                                                            aSymmetricMatrix<double>(x.nrows(),x.ncols()));
    for (std::size_t k=0; k<x.derivative()().size(); ++k)
        x.derivative()()[k].set(i,j,value.derivative()()[k]);
}










template<class Id, class T>
class Derivative<Parameters<Id>,T>{
    Derivative<Matrix<double>,T> m_x;  
    
public:
    operator Matrix<double>&(){return m_x.primitive();}
    
    operator Id&(){return m_x.primitive();}
    auto ncols()const {return m_x.ncols();}
    auto nrows()const {return m_x.nrows();}
    auto size()const {return m_x.size();}
    
    Derivative(){}
    
    template<class dParam>
        requires (std::constructible_from<Derivative<Matrix<double>,T>,dParam>)
    Derivative(dParam&& x):m_x{std::forward<dParam>(x)}{}
    
    
    template<class aParam>
        requires (std::constructible_from<Parameters<Id>,aParam>)
    Derivative(aParam&& x): m_x{std::forward<aParam>(x)()}{}
    
    auto& primitive()const {return m_x.primitive();}
    auto& derivative()const {return m_x.derivative();}    
    
    
    auto& operator()()const
    {
        return m_x;
    }
};



}



#endif // PARAMETERS_DERIVATIVE_H
