#ifndef INDEXED_H
#define INDEXED_H

#include <cstddef>
namespace var {

template<class,class> class Indexed;


template<class ParameterVariable, class Vector>
requires ParameterVariable::is_Parameter && requires (Vector v, std::size_t i){ v[i];}
class Indexed<ParameterVariable,Vector>
{
     
};

}

#endif // INDEXED_H
