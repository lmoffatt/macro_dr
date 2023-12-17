#ifndef CONTINUATION_H
#define CONTINUATION_H
#include "cuevi.h"
#include <cstddef>
namespace cuevi {

template<class ParameterType>
auto continue_evidence(std::string filename, std::size_t max_iter)
{
    std::size_t iter;
    auto Maybe_current =extract_parameters_last<ParameterType>(filename, iter);
    
}
 


} // namespace cuevi


#endif // CONTINUATION_H
