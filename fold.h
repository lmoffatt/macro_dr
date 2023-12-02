#ifndef FOLD_H
#define FOLD_H



#include "maybe_error.h"
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>

template<template<class...> class Container,
         class T,
         class F>
auto Map (Container<T>const & c, F&& f )
{
    using S=decltype(f(c[0]));
    Container<S> out;
    out.reserve(c.size());
    for (auto & e:c)
    {
        auto v=f(e);
        out.push_back(v);
    }
    return out;    
}




template<class F>
auto Map (F&& f )
{
    return [&f](auto const& c) {return Map(c,std::forward<F>(f));};    
}


template<template<class...> class Container,
         class T,
         class F>
auto MapAdj (Container<T>const & c, F&& f )
{
    using S=std::invoke_result_t<F,T,T>;
    Container<S> out;
    out.reserve(c.size()-1);
    for (std::size_t i=1; i<c.size(); ++i)
    {
        auto v=std::invoke(std::forward<F>(f),c[i-1],c[i]);
        out.push_back(v);
    }
    return out;    
}


template<class C, class F>
auto operator | (C&& c, F&& f)->std::invoke_result_t<F,C>
{
    return std::invoke(std::forward<F>(f),std::forward<C>(c));
}


template<class Container,
         class Init,
         class BinaryReductionOp>
Maybe_error<Init> fold (Container c, Init init,
                       BinaryReductionOp reduce )
{
    auto run=init;
    for (auto& e:c)
    {
        auto v_run=std::invoke(reduce,std::move(run),e);
        if (!v_run) return
                v_run.error();
        else
            run=std::move(v_run.value());
    }
    return run;    
}



template<class Container,
         class Map,
         class BinaryReductionOp>
    requires requires (BinaryReductionOp b,Map f, Container c){{b(f(c[0]),f(c[0]))}->std::convertible_to<decltype(f(c[0]))>;} 
   auto foldMap (Container&& c, Map&& f,BinaryReductionOp&& reduce )
{
    auto  v_run=std::forward<Map>(f)(c[0]);
    for (std::size_t i=1; i<c.size(); ++i)
    {
        v_run=std::invoke(reduce,std::move(v_run),std::forward<Map>(f)(c[i]));
        }
    return v_run;    
}


template<class Container,
         class Map,
         class BinaryReductionOp>
    requires requires (BinaryReductionOp b,Map f, Container c){{b(f(c[0]),f(c[0]))}->std::convertible_to<decltype(f(c[0]))>;} 
auto foldMap (std::size_t i_start, std::size_t i_end,Container&& c, Map&& f,BinaryReductionOp&& reduce )
{
    auto  v_run=std::forward<Map>(f)(c[i_start]);
    for (std::size_t i=i_start+1; i<i_end; ++i)
    {
        v_run=std::invoke(reduce,std::move(v_run),std::forward<Map>(f)(c[i]));
    }
    return v_run;    
}


template<class Init,
         class BinaryReductionOp>
Maybe_error<Init> fold (std::size_t i_start, std::size_t i_end,Init init,
                       BinaryReductionOp reduce )
{
    auto run=init;
    for (auto i=i_start; i<i_end; ++i)
    {
        auto v_run=std::invoke(reduce,std::move(run),i);
        if (!v_run)
            return v_run.error();
        else
            run=std::move(v_run.value());
    }
    return run;    
}




#endif // FOLD_H
