#ifndef FOLD_H
#define FOLD_H



#include "maybe_error.h"
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>




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
