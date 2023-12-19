#ifndef MATRIX_RANDOM_H
#define MATRIX_RANDOM_H

#include "matrix.h"
#include "mcmc.h"
#include "multivariate_normal_distribution.h"
#include <iostream>
#include <random>





template<class distribution>
auto random_matrix(mt_64i& mt,distribution&& d,std::size_t nrows,std::size_t ncols)
{
    Matrix<double> out(nrows,ncols,false);
    for (std::size_t i=0; i<out.size(); ++i)
        out[i]=std::forward<distribution>(d)(mt);
    return out;
}


inline auto random_matrix_normal(mt_64i& mt,std::size_t nrows,std::size_t ncols, double mean=0, double stddev=1)
{
    return random_matrix(mt,std::normal_distribution<double>{mean,stddev},nrows,ncols);
}

inline auto random_matrix_exponential(mt_64i& mt,std::size_t nrows,std::size_t ncols, double lambda)
{
    return random_matrix(mt,std::exponential_distribution<double>{lambda},nrows,ncols);
}



//auto random_covariance(mt_64i& mt, std::size_t ndim, std::initializer_list<double> sigmas )
//{
//    auto X=random_matrix_normal(mt,ndim,ndim);
//    auto [Q,R]=qr(X);
//    auto s=diagpos(std::move(sigmas));
//    auto cov=AT_D_A(Q,s);
//    return cov;
//}

inline auto random_covariance(mt_64i& mt, std::size_t ndim, const Matrix<double>& sigmas )
{
    auto X=sample(mt,normal_distribution{0,1},ndim,ndim);
    auto [Q,R]=qr(X);
    auto s=DiagPosDetMatrix(sigmas);
    auto cov=AT_D_A(Q,s);
    return cov;
}


inline auto correlation_matrix(const SymPosDefMatrix<double>& x)
{
    auto out=SymPosDefMatrix<double>(x.nrows(),x.ncols(),false);
    for (std::size_t i=0; i<out.nrows(); ++i)
        for(std::size_t j=0; j<out.ncols(); ++j)
            out.set(i,j,x(i,j)/std::sqrt(x(i,i))/std::sqrt(x(j,j)));
    return out;
}








#endif // MATRIX_RANDOM_H
