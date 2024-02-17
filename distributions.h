#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "matrix.h"
#include "maybe_error.h"
#include <concepts>
#include <cstddef>
#include <functional>
#include <random>
#include <type_traits>
#include <utility>
#include "random_samplers.h"

template <class Distribution>
concept is_Distribution = requires(Distribution &m, Distribution const& m_const) {
    {
        m(std::declval<mt_64i &>())
    } ;
    
    {
        m_const.logP(m(std::declval<mt_64i &>()))
    }->std::convertible_to<Maybe_error<double>>;
    
};



template <class Distribution, class T>
concept is_Distribution_of = requires(Distribution &m, Distribution const& m_const) {
    {
        m(std::declval<mt_64i &>())
    } ->std::convertible_to<T>;
    
    {
        m_const.logP(m(std::declval<mt_64i &>()))
    }->std::convertible_to<Maybe_error<double>>;
};



template<class Sampler, class ProbabilityFunction>
class Custom_Distribution
{
    std::size_t m_size;
    Sampler m_s;
    ProbabilityFunction m_p;
public:
    using T=std::invoke_result_t<Sampler,mt_64i&>;
    Custom_Distribution(std::size_t t_size,Sampler&& s, ProbabilityFunction&& f):m_size{t_size},m_s{std::move(s)}, m_p{std::move(f)}{}
    
    auto operator()(mt_64i& mt){ return std::invoke(m_s,mt);}
    Maybe_error<double> logP(const T& x)const { return std::invoke(m_p,x);}
    
    auto size()const {return m_size;}
};
    


template <class Distribution>
concept is_vector_sampler = requires(Distribution &m) {
    { m(std::declval<mt_64i &>,0)}->index_accesible;
};

template<class T>
    requires (index_accesible<T>)
auto& get_at(T& v, std::size_t i)
{
    return v[i];
}

template<class T>
auto get_at(T& v, std::size_t i)->std::decay_t<decltype(v[i])>
{
    return v[i];
}
inline double get_at(double const& x,std::size_t i){
    assert(i==0);
    return x;
}


inline double& get_at(double& x,std::size_t i){
    assert(i==0);
    return x;
}

template<class T>
    requires (has_size<T>)
std::size_t size(const T& x) {return x.size();}

template<class T>
    requires (!has_size<T>&&has_size<std::decay_t<decltype(std::declval<T>()())>>)
std::size_t size(const T& x) {return x().size();}


constexpr std::size_t size(double){return 1;}


template<class F, class T>
auto operator|(F&& f, T&& x )->std::invoke_result_t<F,T>
{
    return std::invoke(std::forward<F>(f),std::forward<T>(x));
}


inline auto& append_to(Matrix<double>& m, std::size_t )
{
    return m;
}

template<class X,class ...Xs>
auto& append_to(Matrix<double>& m, std::size_t i, X&& x, Xs&&...xs)
{
    for (std::size_t n=0; n<size(x); ++n)
        m[i+n]=get_at(x,n);
    return append_to(m,i+size(x),std::forward<Xs>(xs)...);
}


template<class ...Xs>
auto concatenate_to_columns(Xs&&...xs)
{
    auto n=(size(xs)+...);
    auto out=Matrix<double>(1,n);
    std::size_t i=0;
    out=append_to(out,i,std::forward<Xs>(xs)...);
    return out;
}


template<class Dist>
concept Multivariate=requires (Dist& d)
{
    {d(std::declval<mt_64i&>())}->std::convertible_to<Matrix<double>>;
};

inline double logP_impl(const Matrix<double>& , std::size_t, double partial_logP)
{
    return partial_logP;
}
template<class Dist, class... Ds>
    requires(Multivariate<Dist>)
Maybe_error<double> logP_impl(const Matrix<double>& x, std::size_t ipos, double partial_logP, const Dist& d, const Ds&...ds)
{
    auto n=d.size();
    auto out=Matrix<double>(1,size(d));
    for (std::size_t i=0; i<n; ++i)
        out[i]=x[ipos+i];
    auto logPi=d.logP(out);
    if(logPi)
        return logP_impl(x,ipos+n,partial_logP+logPi.value(),ds...);
    else
        return logPi.error()+"\n log_impl";
}
template<class Dist, class... Ds>
    requires(!Multivariate<Dist>)
Maybe_error<double> logP_impl(const Matrix<double>& x, std::size_t ipos, double partial_logP, const Dist& d, const Ds&...ds)
{
    auto logPi=d.logP(x[ipos]);
    if (logPi)
        return logP_impl(x,ipos+1,partial_logP+logPi.value(),ds...);
    else
        return logPi.error()+"\n log_impl";
    
}

template<class D>
    requires(is_Distribution<D>)
auto sample(mt_64i& mt, D&d)
{return d(mt);}

template<class D>
    requires(is_Distribution<D>)
auto sampler(const D&d)
{return D(d);}


template<class D>
    requires(is_vector_sampler<D>)
auto sample(mt_64i& mt, D&d, std::size_t n)
{
    return d(mt,n);
}

template<class D>
    requires(is_Distribution<D>)
auto sample(mt_64i& mt, D&&d, std::size_t nrows, std::size_t ncols)
{
    Matrix<double> out(nrows,ncols,false);
    for (std::size_t i=0; i<out.size(); ++i)
        out[i]=std::forward<D>(d)(mt);
    return out;
}






template<class... ds>
    requires(is_Distribution<ds>&&...)
class distributions: public ds...
{
public:
    
    auto operator()(mt_64i& mt)
    {
        return concatenate_to_columns(ds::operator()(mt)...);
    }
    
    Maybe_error<double> logP(const Matrix<double>& x)const
    {
        return logP_impl(x,0ul,0.0,static_cast<ds const&>(*this)...);
    }
    
    explicit distributions(ds&&... d): ds{std::move(d)}...{}
    explicit distributions(ds const &... d): ds{d}...{}
};




#endif // DISTRIBUTIONS_H
