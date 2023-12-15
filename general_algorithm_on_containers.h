#ifndef GENERAL_ALGORITHM_ON_CONTAINERS_H
#define GENERAL_ALGORITHM_ON_CONTAINERS_H

#include <limits>
#include <type_traits>
namespace var{

template<class Container>
concept is_Container= requires(Container const&  c)
{{c.size()}->std::convertible_to<std::size_t>;
    {c[std::size_t{}]};};
    
    
    
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




auto max(is_Container auto const& c)
{
    using T=std::decay_t<decltype(c[0])>;
    using std::numeric_limits;
    auto out=numeric_limits<T>::min();
    for (std::size_t i=0; i<c.size(); ++i)
        if (out<c[i])
            out=c[i];
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
    using T=std::decay_t<decltype(c[0])>;
    using std::numeric_limits;
    auto out=numeric_limits<T>::max();
    for (std::size_t i=0; i<c.size(); ++i)
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






}

#endif // GENERAL_ALGORITHM_ON_CONTAINERS_H
