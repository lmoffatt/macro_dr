#ifndef FOLD_H
#define FOLD_H



#include "maybe_error.h"
#include <functional>








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
        if (!v_run) return v_run.error();
        else
            run=std::move(v_run.value());
    }
    return run;    
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
        if (!v_run) return v_run.error();
        else
            run=std::move(v_run.value());
    }
    return run;    
}




#endif // FOLD_H
