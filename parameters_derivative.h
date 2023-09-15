#ifndef PARAMETERS_DERIVATIVE_H
#define PARAMETERS_DERIVATIVE_H

#include "parameters.h"
//#include "matrix_derivative.h"

namespace var {

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
    
    template <class F>
    friend auto zip_par(F &&f, d_d_ const &a,d_d_<double,Parameters<Id>>const &b) {
        using S=std::decay_t<std::invoke_result_t<F,aMatrix<double>,double>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(a()[i],b()[i]);
        return d_d_<S,Parameters<Id>>(x);
    }
    
    template <class F>
    friend auto zip_par(F &&f, d_d_<double,Parameters<Id>>const &b,d_d_ const &a) {
        using S=std::decay_t<std::invoke_result_t<F,double,aMatrix<double>>>;
        Matrix<S> x(a().nrows(), a().ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(b()[i],a()[i]);
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
    
    auto& primitive() {return m_x;}
    auto& primitive()const {return m_x;}
    auto& derivative()const {return m_d;}
    
    friend auto operator/(const Derivative& x, const Derivative& y){
        auto fx=x.primitive();
        auto fy=x.primitive();
        if ((x.derivative()().size()>0)&&(y.derivative()().size()>0))
            return Derivative(fx/fy,zip([fx,fy](auto dx,auto dy) {return dx/fy-fx/fy/fy*dy;},x.derivative()(),y.derivative()()));
        else if (x.derivative()().size()==0)
            return Derivative(fx/fy,fx/fy/fy*y.derivative()());
        else
            return Derivative(fx/fy,x.derivative()()/fy);
        
        
    }
    friend auto operator*(const Derivative& x, const Derivative& y){
        auto fx=x.primitive();
        auto fy=x.primitive();
        if ((x.derivative()().size()>0)&&(y.derivative()().size()>0))
             return Derivative(fx*fy,zip([fx,fy](auto dx,auto dy) {return dx*fy+fx*dy;},x.derivative()(),y.derivative()()));
        else if (x.derivative()().size()==0)
             return Derivative(fx*fy,fx*y.derivative()());
        else
             return Derivative(fx*fy,x.derivative()()*fy);
             
    }
    
    
    friend auto exp(const Derivative& x){
        auto f=exp(x.primitive());
        return Derivative(f,f*x.derivative()());
    }
    
    friend auto log(const Derivative& x){
        auto f=log(x.primitive());
        return Derivative(f,x.derivative()()*(1.0/x.primitive()));
    }
    
    
    friend auto abs(const Derivative& x){
        auto f=std::abs(x.primitive());
        return Derivative(f,((x.primitive() > 0.0) ? 1.0 : ((x.primitive() < 0) ? -1.0 : 0.0))*x.derivative()());
    }
    
    
    friend bool operator==(Derivative const& one, double val)
    {return one.primitive()==val;}
    
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
    Derivative(std::size_t nrows, std::size_t ncols, double v):m_x(nrows,ncols,v){}
    
    
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
    
    auto friend operator*(const Derivative& x, double y)
    {
        return Derivative(x.primitive()*y,x.derivative()()*y);
    }
    
    
    
    
    
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
    if (x.derivative()().size()==0)
        x.derivative()()=Matrix<aSymmetricMatrix<double>>(value.derivative()().nrows(),value.derivative()().ncols(),
                                                            aSymmetricMatrix<double>(x.nrows(),x.ncols()));
    for (std::size_t k=0; k<x.derivative()().size(); ++k)
        x.derivative()()[k].set(i,j,value.derivative()()[k]);
}










template<class Id, class Id2>
class Derivative<Parameters<Id>,Parameters<Id2>>{
    Derivative<Matrix<double>,Parameters<Id2>> m_x;  
    
public:
    operator Matrix<double>&(){return m_x.primitive();}
    
    operator Id&(){return m_x.primitive();}
    auto ncols()const {return m_x.ncols();}
    auto nrows()const {return m_x.nrows();}
    auto size()const {return m_x.size();}
    
    Derivative(){}
    
    template<class dParam>
        requires (std::constructible_from<Derivative<Matrix<double>,Parameters<Id2>>,dParam>)
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

template<class Id>
Derivative<Parameters<Id>,Parameters<Id>> selfDerivative(const Parameters<Id>& x)
{
    auto out=Matrix<Matrix<double>>(x().nrows(),x().ncols(),Matrix<double>(x().nrows(),x().ncols(),0.0));
    for (std::size_t i=0; i<x().size(); i++)
    {
        out[i][i]=1.0;
    }
    
    return Derivative<Parameters<Id>,Parameters<Id>>(Derivative<Matrix<double>,Parameters<Id>>(x(),d_d_<Matrix<double>,Parameters<Id>>(out)));
}

template <class T, class Id> auto diag(const Derivative<Matrix<T>,Parameters<Id>> &a) {
    
    return Derivative<DiagonalMatrix<T>,Parameters<Id>>(diag(primitive(a)),
                                            apply_par([](auto const & d){return diag(d);}, derivative(a)));
}


template <class T, class Id> auto diagpos(const Derivative<Matrix<T>,Parameters<Id>> &a) {
    
    
    
    
    return Derivative<DiagPosDetMatrix<T>,Parameters<Id>>(diagpos(primitive(a)),
                                              apply_par([](auto const & d){
                                                  return diag(d);
                                              }, derivative(a)));
}


template<class Id>
auto XTX(const Derivative<Matrix<double>,Parameters<Id>> &a) {
    
    auto& f=primitive(a);
    return Derivative<SymPosDefMatrix<double>,Parameters<Id>>(XTX(f),
                                                  apply_par([&f](auto const & d){return X_plus_XT(tr(d)*f);}, derivative(a)));
}


template<class Id>
auto X_plus_XT(const Derivative<Matrix<double>,Parameters<Id>> &a) {
    
    auto& f=primitive(a);
    return Derivative<SymmetricMatrix<double>,Parameters<Id>>(X_plus_XT(f),
                                                  apply_par([](auto const & d){return X_plus_XT(d);}, derivative(a)));
}


template<class Id,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto tr(const Derivative<aMatrix<double>,Parameters<Id>> &a) {
    
    auto& f=primitive(a);
    return Derivative<aMatrix<double>,Parameters<Id>>(tr(f),
                                          apply_par([](auto const & d){return tr(d);}, derivative(a)));
}


template<class Id,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemDiv(const Derivative<aMatrix<double>,Parameters<Id>> &a,const Derivative<aMatrix<double>,Parameters<Id>> &b) {
    
    return zip([] (auto& x, auto& y){ return x/y;}, a,b);
}

template<class Id,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto elemMult(const Derivative<aMatrix<double>,Parameters<Id>> &a,const Derivative<aMatrix<double>,Parameters<Id>> &b) {
    
    return zip([] (auto& x, auto& y){ return x*y;}, a,b);
}

template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>,Parameters<Id>> &a, const Derivative<bMatrix<double>,Parameters<Id>> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<S,Parameters<Id>>(TranspMult(fa,fb),
                            zip_par([&fa,&fb](auto const & da, auto const & db)
                                    {return TranspMult(fa,db)+TranspMult(da,fb);}, derivative(a), derivative(b)));
}


template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const Derivative<aMatrix<double>,Parameters<Id>> &a, const bMatrix<double> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<S,Parameters<Id>>(TranspMult(fa,fb),
                            apply_par([&fb](auto const & da)
                                      {return TranspMult(da,fb);}, derivative(a)));
}

template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto TranspMult(const aMatrix<double>&a, const Derivative<bMatrix<double>,Parameters<Id>> &b) {
    using S=std::decay_t<decltype(TranspMult(aMatrix<double>{},bMatrix<double>{}))>;
    auto& fa=a;
    auto& fb=primitive(b);
    
    return Derivative<S,Parameters<Id>>(TranspMult(fa,fb),
                            apply_par([&fa](auto const & db)
                                      {return TranspMult(fa,db);}, derivative(b)));
}

template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>,Parameters<Id>> &a, const Derivative<bMatrix<double>,Parameters<Id>> &b) {
    using S=std::decay_t<decltype(aMatrix<double>{}*bMatrix<double>{})>;
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    if ((derivative(a)().size()>0)&&(derivative(b)().size()>0))
    
        return Derivative<S,Parameters<Id>>(fa*fb,
                            zip_par([&fa,&fb](auto const & da, auto const & db)
                                    {return fa*db+da*fb;}, derivative(a), derivative(b)));
    else if ((derivative(a)().size()==0)&&(derivative(a)().size()==0))
        return Derivative<S,Parameters<Id>>(fa*fb);
    else if (derivative(a)().size()==0)   
        return Derivative<S,Parameters<Id>>(fa*fb,
                                             apply_par([&fa]( auto const & db)
                                                       {return fa*db;}, derivative(b)));
    else   
        return Derivative<S,Parameters<Id>>(fa*fb,
                                             apply_par([&fb]( auto const & da)
                                                       {return da*fb;}, derivative(a)));
    
}
template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>,Parameters<Id>> &a, bMatrix<double> &b) {
    using S=std::decay_t<decltype(aMatrix<double>{}*bMatrix<double>{})>;
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<S,Parameters<Id>>(fa*fb,
                            apply_par([&fa,&fb](auto const & da)
                                    {return da*fb;}, derivative(a)));
}

template<class Id,template<class> class aMatrix,template<class> class bMatrix >
    requires aMatrix<double>::is_Matrix
auto operator*(const bMatrix<double> &b,const Derivative<aMatrix<double>,Parameters<Id>> &a) {
    using S=std::decay_t<decltype(aMatrix<double>{}*bMatrix<double>{})>;
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<S,Parameters<Id>>(fb*fa,
                            apply_par([&fa,&fb](auto const & da)
                                    {return fb*da;}, derivative(a)));
}


template<class Id,template<class> class aMatrix >
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<aMatrix<double>,Parameters<Id>> &a, const Derivative<double,Parameters<Id>> &b) {
    using S=std::decay_t<decltype(aMatrix<double>{}*double{})>;
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<S,Parameters<Id>>(fa*fb,
                            zip_par([&fa,&fb](auto const & da, auto const & db)
                                    {return fa*db+da*fb;}, derivative(a), derivative(b)));
}


template<class Id,template<class> class aMatrix>
    requires aMatrix<double>::is_Matrix
auto operator*(const Derivative<double,Parameters<Id>> &b,const Derivative<aMatrix<double>,Parameters<Id>> &a) {
    using S=std::decay_t<decltype(double{}*std::declval<aMatrix<double>>())>;
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<S,Parameters<Id>>(fb*fa,
                            zip_par([&fa,&fb](auto const & da, auto const & db)
                                    {return db*fa+fb*da;}, derivative(a), derivative(b)));
}




template <class Id>
Derivative<SymPosDefMatrix<double>,Parameters<Id>> AT_B_A(const Derivative<Matrix<double>,Parameters<Id>> &a, const Derivative<SymmetricMatrix<double>,Parameters<Id>> &b) {
    auto& fa=primitive(a);
    auto& fb=primitive(b);
    
    return Derivative<SymPosDefMatrix<double>,Parameters<Id>>(AT_B_A(fa,fb),
                                                  zip_par([&fa,&fb](auto const & da, auto const & db)
                                                          {return AT_B_A(fa,db)+X_plus_XT(TranspMult(da,fb*fa));}, derivative(a), derivative(b)));
}



template <class Id>
Derivative<SymPosDefMatrix<double>,Parameters<Id>> AT_B_A(const Derivative<Matrix<double>,Parameters<Id>> &a, const SymmetricMatrix<double> &b) {
    auto& fa=primitive(a);
    auto& fb=b;
    
    return Derivative<SymPosDefMatrix<double>,Parameters<Id>>(AT_B_A(fa,fb),
                                                  apply_par([&fa,&fb](auto const & da)
                                                            {return X_plus_XT(TranspMult(da,fb*fa));}, derivative(a)));
}




template<template<class>class Matrix, class Id>
    requires Matrix<double>::is_Matrix
auto getvalue(const Derivative<Matrix<double>,Parameters<Id>> &x) {
    assert(x.primitive().size() == 1);
    Matrix<double> der(x.derivative()().nrows(),x.derivative()().ncols());
    for (std::size_t i=0; i<der.size(); ++i) der[i]=x.derivative()()[i][0];
    return Derivative<double,Parameters<Id>>(x.primitive()[0],std::move(der)); 
}




template<class Id, template<class>class Matrix>
    requires Matrix<double>::is_Matrix
Maybe_error<Derivative<Matrix<double>,Parameters<Id>>> inv(const Derivative<Matrix<double>,Parameters<Id>>& x)
{
    auto inv_x=inv(x.primitive());
    if (!inv_x)
        return inv_x.error();
    else
    {
        auto dinv=apply([&inv_x](auto const& dx){ return -inv_x.value()*dx*inv_x.value();},x.derivative()());
        return Derivative<Matrix<double>,Parameters<Id>>(inv_x.value(),dinv);
    }
}

template<class Id>
Maybe_error<std::tuple<Derivative<Matrix<double>,Parameters<Id>>,Derivative<DiagonalMatrix<double>,Parameters<Id>>,Derivative<Matrix<double>,Parameters<Id>>>>
eigs(const Derivative<Matrix<double>,Parameters<Id>> &x, bool does_permutations = true,
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
        
        auto dLambda=Derivative<DiagonalMatrix<double>,Parameters<Id>>(lambda,derlambda);
        
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
        auto dVR=Derivative<Matrix<double>,Parameters<Id>>(VR,derVR);
        auto dVL=inv(dVR);
        if (! dVL)
            return dVL.error();
        else
            return std::tuple(dVR,dLambda,dVL.value());
    }
}








}



#endif // PARAMETERS_DERIVATIVE_H
