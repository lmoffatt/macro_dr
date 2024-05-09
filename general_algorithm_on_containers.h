#ifndef GENERAL_ALGORITHM_ON_CONTAINERS_H
#define GENERAL_ALGORITHM_ON_CONTAINERS_H

//#include "derivative_operator.h"
#include "maybe_error.h"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
namespace var{

template<class Container>
concept is_Container= requires(Container const&  c)
{{c.size()}->std::convertible_to<std::size_t>;
    {c[std::size_t{}]};};

template<class T>
concept StringLike = std::is_convertible_v<T, std::string_view>;
    
    
auto sum(is_Container auto const& c)
{
    using T=std::decay_t<decltype(c[0])>;
    
    auto out=T{};
    for (std::size_t i=0; i<c.size(); ++i)
        out=out+c[i];
    return out; 
}

auto cumsum(is_Container auto const& c)
{
    using T=std::decay_t<decltype(c[0])>;
        auto out=c;    
    auto sum=T{};
    for (std::size_t i=0; i<c.size(); ++i)
    {
        sum=sum+c[i];
        out[i]=sum;
    }
    return out; 
}

inline double fullsum(double x){return x;}

 double fullsum(std::integral auto x){return x;}
 
 
 
 auto fullsum(is_Container auto const& c)
 {
     double out=0.0;
     for (std::size_t i=0; i<c.size(); ++i)
     {
         out+=fullsum(c[i]);
     }
     return out; 
 }
 



auto i_max(is_Container auto const& c)
{
    std::size_t im=0;
    for (std::size_t i=1; i<c.size(); ++i)
        if (c[im]<c[i])
            im=i;
    return im; 
}

auto i_min(is_Container auto const& c)
{
    std::size_t im=0;
    for (std::size_t i=1; i<c.size(); ++i)
        if (c[im]>c[i])
            im=i;
    return im; 
}

auto min(is_Container auto const& c)
{
    auto out=c[0];
    for (std::size_t i=1; i<c.size(); ++i)
        if (out>c[i])
            out=c[i];
    return out; 
}


auto count_nan(is_Container auto const& c)
{
    auto out=0ul;
    for (std::size_t i=0; i<c.size(); ++i)
        if (std::isnan(c[i]))
            ++out;
    return out; 
}

template<class T>
    requires (std::integral<T>|| std::is_constructible_v<std::string,T>)
inline Maybe_error<bool> compare_contents(T s0, T s1,double =0, double=0,std::size_t =1)
{
    if (s0!=s1)
    {
        std::stringstream ss;
        ss<<"different :\n"<<s0<<"\n"<<s1;
        std::cerr<<"different :\n"<<s0<<"\n"<<s1;
        return error_message(ss.str());
    }
    else
        return true;
}



inline Maybe_error<bool> compare_contents(std::floating_point auto s0, std::floating_point auto s1,double RelError=std::numeric_limits<double>::epsilon()*100, double AbsError=std::numeric_limits<double>::epsilon()*100,std::size_t=1)
{
    if (std::abs(s0-s1)>std::max(AbsError,RelError*std::max(std::abs(s0), std::abs(s1))))
    {
        std::stringstream ss;
        ss<< std::setprecision(std::numeric_limits<double>::digits10 + 1);
        ss<<"relative difference greater than "<<RelError<<":\n"<<s0<<"\n"<<s1;
        return error_message(ss.str());
    }
    else
        return true;
}

inline Maybe_error<bool> compare_contents(is_Container auto const& s0, is_Container auto const& s1,double RelError=std::numeric_limits<double>::epsilon()*100,  double AbsError=std::numeric_limits<double>::epsilon()*100,std::size_t max_errors=10)
{
    if (s0==s1) return true;
    std::size_t n_errors=0;
    std::string message;
    if (s0.size()!=s1.size())
    {
       message+= "differ in size: "+std::to_string(s0.size())+" vs "+std::to_string(s1.size());
        ++n_errors;
    }
    std::size_t i=0;
    auto n=std::min(s0.size(),s1.size());
    while ((n_errors<max_errors)&&(i<n))
    {
        auto Maybe_equal=compare_contents(s0[i], s1[i],RelError,AbsError);
        if (!Maybe_equal)
        {
            message+="\n "+std::to_string(i)+"th element \n"+Maybe_equal.error()()+"\n";
            ++n_errors; 
        }
        ++i;
    }
    if (n_errors>0)
    {
        return error_message(message);
        
    }
    else
        return true;
}



}

#endif // GENERAL_ALGORITHM_ON_CONTAINERS_H
