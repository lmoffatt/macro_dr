#ifndef MCMC_H
#define MCMC_H

#include "matrix.h"
#include <cassert>
#include <random>
#include <utility>
#include <vector>
#include "distributions.h"
//using Parameters = Matrix<double>;
using Data = Matrix<double>;

using WalkerIndexes = std::vector<std::size_t>;



/*
auto operator+(const Parameters &x, const Parameters &y) {
  assert((x.size() == y.size()) && "sum of vector fields of different sizes");
  auto out = x;
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] + y[i];
  return out;
}
auto operator-(const Parameters &x, const Parameters &y) {
  assert((x.size() == y.size()) && "sum of vector fields of different sizes");
  auto out = x;
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] - y[i];
  return out;
}
*/
auto calc_seed(typename std::mt19937_64::result_type initseed) {
    
    if (initseed == 0) {
        std::random_device rd;
        std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
        
        return useed(rd);
    } else
        return initseed;
}

auto init_mt(typename std::mt19937_64::result_type initseed) {
    initseed = calc_seed(initseed);
    return std::mt19937_64(initseed);
}

template<class D, class T>
    requires(is_Distribution_of<D,T>)
auto logPrior(const D&d, const T& x )
{return d.logP(x);}



template <class Distribution, class Parameters>
concept is_sampler = requires(Distribution & dist) {
    {
        sample(std::declval<std::mt19937_64 &>(),dist)
    } -> std::convertible_to<Parameters>;
};






template <class Prior, class Parameters,class Variables,class DataType>
concept is_prior = requires(Prior const &prior,
                            const Parameters& p,
                            const Variables& var,
                            const DataType& y) {
    {
        sampler(prior)
    } -> is_sampler<Parameters>;
    
    {
        logPrior(prior,p)
    }->std::convertible_to<Maybe_error<double>>;
    
};


template <class FunctionTable,class Likelihood, class Parameters,class Variables,class DataType>
concept is_likelihood_model = requires(FunctionTable&& f,
                                       Likelihood const &lik,
                                       const Parameters& p,
                                       const Variables& var,
                                       const DataType& y) {
    
    {
        simulate(std::declval<std::mt19937_64 &>(),lik,p,var)
    }-> std::convertible_to<DataType>;
    
    {
        logLikelihood(std::forward<FunctionTable>(f),lik,p,y,var)
    }->std::convertible_to<Maybe_error<double>>;
};


struct logLikelihood_f{
    friend constexpr std::string ToString(logLikelihood_f){ return "logLikelihood_f";}
};






template<class Parameters>
struct mcmc {
    Parameters parameter;
    double logP;
    double logL;
};

template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<std::mt19937_64 &>(), std::declval<Prior&>()))>>
 //   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_mcmc(FunctionTable&& f, std::mt19937_64 &mt, Prior const & pr, const Lik& lik,
               const DataType &y, const Variables &x) {
    auto priorsampler=sampler(pr);
    auto par = sample(mt,priorsampler);
    auto logP = logPrior(pr,par);
    auto logL = logLikelihood(f,lik,par, y,x);
    while(!(logP)||!(logL))
    {
        par = sample(mt,priorsampler);
        logP = logPrior(pr,par);
        logL = logLikelihood(f,lik,par, y,x);
        
    }
    return mcmc<Parameters>{std::move(par), logP.value(), logL.value()};
}

#endif // MCMC_H
