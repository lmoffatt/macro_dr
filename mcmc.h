#ifndef MCMC_H
#define MCMC_H

#include "matrix.h"
#include <cassert>
#include <random>
#include <utility>
#include <vector>
#include "distributions.h"
#include "variables.h"
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
inline auto calc_seed(typename mt_64i::result_type initseed) {
    
    if (initseed == 0) {
        std::random_device rd;
        std::uniform_int_distribution<typename mt_64i::result_type> useed;
        
        return useed(rd);
    } else
        return initseed;
}

inline auto init_mt(typename mt_64i::result_type initseed) {
    initseed = calc_seed(initseed);
    return mt_64i(initseed);
}

template<class D, class T>
    requires(is_Distribution_of<D,T>)
auto logPrior(const D&d, const T& x )
{return d.logP(x);}



template <class Distribution, class Parameters>
concept is_sampler = requires(Distribution & dist) {
    {
        sample(std::declval<mt_64i &>(),dist)
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
        simulate(std::declval<mt_64i &>(),lik,p.to_value(),var)
    }-> std::convertible_to<DataType>;
    
    {
        logLikelihood(std::forward<FunctionTable>(f),lik,p.to_value(),y,var)
    }->std::convertible_to<Maybe_error<double>>;
};


struct logLikelihood_f{
    friend constexpr std::string ToString(logLikelihood_f){ return "logLikelihood_f";}
};

struct calc_logLikelihood_f{
    friend constexpr std::string ToString(calc_logLikelihood_f){ return "calc_logLikelihood_f";}
};


class Trial_count : public var::Constant<Trial_count, std::size_t> {};




class Success_count : public var::Constant<Success_count, std::size_t> {};

class Trial_statistics
    : public var::Constant<Trial_statistics,
                           var::Vector_Space<Trial_count, Success_count>> {
public:
    Trial_statistics &operator+=(const Trial_statistics other) {
        get<Trial_count>((*this)())() += get<Trial_count>(other())();
        get<Success_count>((*this)())() += get<Success_count>(other())();
        return *this;
    }
    
    friend void succeeds(Trial_statistics &me) {
        ++get<Trial_count>(me())();
        ++get<Success_count>(me())();
    }
    friend void fails(Trial_statistics &me) { ++get<Trial_count>(me())(); }
    
    void reset() {
        get<Trial_count>((*this)())() = 0;
        get<Success_count>((*this)())() = 0;
    }
    auto count() const { return get<Trial_count>((*this)())(); }
    double rate() const {
        return 1.0 * get<Success_count>((*this)())() /
               get<Trial_count>((*this)())();
    }
};

class emcee_Step_statistics
    : public var::Constant<emcee_Step_statistics, Trial_statistics> {};

class Thermo_Jump_statistics
    : public var::Constant<Thermo_Jump_statistics, Trial_statistics> {};




template<class Parameters>
struct mcmc {
    Parameters parameter;
    double logP;
    double logL;
};

template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<mt_64i &>(), std::declval<Prior&>()))>>
 //   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_mcmc(FunctionTable&& f, mt_64i &mt, Prior const & pr, const Lik& lik,
               const DataType &y, const Variables &x) {
    auto& priorsampler=pr;
    auto par = sample(mt,priorsampler);
    auto logP = logPrior(pr,par);
    auto logL = logLikelihood(std::forward<FunctionTable>(f),lik,par.to_value(), y,x);
    while(!(logP)||!(logL))
    {
        par = sample(mt,priorsampler);
        logP = logPrior(pr,par);
        logL = logLikelihood(f,lik,par.to_value(), y,x);
        
    }
    return mcmc<Parameters>{std::move(par), logP.value(), logL.value()};
}

#endif // MCMC_H
