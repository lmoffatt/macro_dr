#ifndef TYPE_ALGEBRA_H
#define TYPE_ALGEBRA_H


#include "maybe_error.h"
#include <type_traits>
#include <utility>
#include <variant>

namespace in_progress{
template<class...T> struct P{
    template<class... T2>
    friend auto operator*( P, P<T2...>)
    {
        return P<T...,T2...>{};
    }
    
    
    template<class...aT>
    requires(std::is_convertible_v<aT,T>&&...)
    auto operator()(aT&&...x)const
    {
        return std::tuple<T...>(std::forward<aT>(x)...);
    }
    
 };


template<class... T> struct S{
    
    template<class aT>
        requires(std::is_convertible_v<aT,T>||...)
    auto operator()(aT&& x)const
    {
        return std::variant<T...>(std::forward<aT>(x));
    }
    
    template<class... T2>
    friend auto operator+( S, S<T2...>)
    {
        return S<T...,T2...>{};
    }
    template<class... T2>
    friend auto operator+( S, P<T2...>)
    {
        return S<T...,P<T2...>>{};
    }
    
    template<class... T2>
    friend auto operator+( S, std::variant<T2...>)
    {
        return S<T...,T2...>{};
    }
    
    
    template<class... T2>
    friend auto operator+( P<T2...>,S)
    {
        return S<P<T2...>,T...>{};
    }
    
    
    template<class T2>
    friend auto operator*( S, S<T2>)
    {
        return (S<>{}+...+P<T,T2>{});
    }
    template<class... T2>
    friend auto operator*( S, S<T2...>)
    {
        return (S<>{}+...+(S{}*S<T2>{}));
    }
    template<class... T2>
    friend auto operator*( S, std::variant<T2...>)
    {
        return (S<>{}+...+(S{}*S<T2>{}));
    }
    
};
}

template<auto x>
struct V{};


template<class...T0>

auto prod(std::variant<T0...>const & x)
{
    return x;
}

template<class...T0, class... T1>
             
auto prod_variant(std::tuple<T0...>const & x,const std::variant<T1...>& y)
{
    return std::visit([&x](auto & e)->std::variant<std::tuple<T0...,T1>...>
                      { return std::apply([&e](auto&...t){return std::tuple(e,t...);},x);},y);
}



template<class...Tu, class... T1, class...Vs>
    requires ((is_of_this_template_type_v<Tu,std::tuple>&&...)&&
             (is_of_this_template_type_v<Vs,std::variant>&&...)
             )
auto prod_variant(std::variant<Tu...>const & x,const std::variant<T1...>& y, const Vs&...vs)
{
    return std::visit([&y, &vs...](auto& ex)->std::variant<std::decay_t<decltype(prod_variant(prod_variant(Tu{},y),vs...))>...>
                      {return prod_variant(prod_variant(ex,y),vs...);},x);    
}

template<class F,class...Tus>
    requires (is_of_this_template_type_v<Tus,std::tuple>&&...)
auto Apply_variant(F&& f,std::variant<Tus ...> const&x)
{
    return std::visit([&f](auto& tu)->std::variant<std::decay_t<
                                         decltype(std::apply([f](auto&...ts){return f(ts...);},Tus{}))>...>
                      {
        return std::apply([f](auto&...ts){return f(ts...);},tu);},x);
}


template<class F,class...Vs>
    requires (is_of_this_template_type_v<Vs,std::variant>&&...)
auto Apply_variant(F&& f,std::tuple<Vs ...> const&x)
{
    auto variant_of_tuples=std::apply([](auto&...e){return prod_variant(e...);},x);
    
    return Apply_variant(std::forward<F>(f),variant_of_tuples);      
}





namespace in_progress{
template <class F, class...Ts>
auto Apply(F,P<Ts...>)
{
    return F{}(Ts{}...);
}




template <class F, class...Ps>
auto Map(F,S<Ps...>)
{
    return (S<>{}+...+Apply(F{},Ps{}));
}

template <class...Ts>
auto Variant(S<Ts...>)
{
    return std::variant<Ts...>{};   
}




template <class S>
using Variant_t=decltype(Variant(S{}));



}


#endif // TYPE_ALGEBRA_H
