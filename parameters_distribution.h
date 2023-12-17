#ifndef PARAMETERS_DISTRIBUTION_H
#define PARAMETERS_DISTRIBUTION_H

#include "parameters.h"
#include "multivariate_normal_distribution.h"
#include <vector>

namespace var {

template<class Id>
class Parameters_Normal_Distribution: public multivariate_normal_distribution<double,DiagPosDetMatrix<double>>
{
private:
    
    
public:
    using base_type=multivariate_normal_distribution<double,DiagPosDetMatrix<double>>;
    
    
    Parameters<Id> operator()(mt_64i &mt) {
        return base_type::operator ()(mt);
    }
    
    Maybe_error<double> logP(const Matrix<double>& x)const
    {
        return base_type::logP(x); 
    }
    Maybe_error<double> logP(const Parameters<Id>& x)const
    {
        return base_type::logP(x()); 
    }
    
    
    
    
};


template<class Id>
auto prior_around(const Parameters<Id>& x, double error )
{
    assert(error>0);
    auto Maybe_dist=make_multivariate_normal_distribution(x(),DiagPosDetMatrix<double>(x().size(),x().size(),error));
    
    return Parameters_Normal_Distribution<Id>{Maybe_dist.value()};
}

template<class Id>
auto prior_around(const Parameters<Id>& x, std::vector<double> values )
{
    assert(*std::min(values.begin(), values.end())>0);
    auto Maybe_dist=make_multivariate_normal_distribution(x(),DiagPosDetMatrix<double>(values));
    
    return Parameters_Normal_Distribution<Id>{Maybe_dist.value()};
}



}


#endif // PARAMETERS_DISTRIBUTION_H
