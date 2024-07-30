#ifndef MCMC_H
#define MCMC_H

#include "matrix.h"
#include <cassert>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include "distributions.h"
#include "maybe_error.h"
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
concept is_likelihood_model = requires(FunctionTable& f,
                                       Likelihood const &lik,
                                       const Parameters& p,
                                       const Variables& var,
                                       const DataType& y) {
    
    {
        simulate(std::declval<mt_64i &>(),lik,p.to_value(),var)
    }-> std::convertible_to<DataType>;
    
    {
        logLikelihood(f,lik,p.to_value(),y,var)
    }->std::convertible_to<Maybe_error<logLs>>;
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
    using base_type=var::Constant<Trial_statistics,
                                    var::Vector_Space<Trial_count, Success_count>>;
    Trial_statistics &operator+=(const Trial_statistics other) {
        get<Trial_count>((*this)())() += get<Trial_count>(other())();
        get<Success_count>((*this)())() += get<Success_count>(other())();
        return *this;
    }
    Trial_statistics(Trial_count n, Success_count s): var::Constant<Trial_statistics,
                        var::Vector_Space<Trial_count, Success_count>>(var::Vector_Space(n,s)){}
    
    Trial_statistics(): base_type{var::Vector_Space(Trial_count(0), Success_count(0))}{}
    
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
    Maybe_error<double> rate() const {
        if (count()>0)
            return 1.0 * get<Success_count>((*this)())() /count();
        else
            return error_message("");
    }
};


struct count: public var::Constant<count, std::size_t>{};

template <class Va>
struct mean: public var::Var<mean<Va>,Va>{};

template <class Va>
struct variance: public var::Var<variance<Va>,Va>{};









template <class Va>
    class Moment_statistics : public var::Var<Moment_statistics<Va>,
                                          var::Vector_Space<mean<Va>,variance<Va>,count>>
{
public:
    using base_type=var::Var<Moment_statistics<Va>,
        var::Vector_Space<mean<Va>,variance<Va>,count>>;
    

public:
    Moment_statistics():base_type{var::Vector_Space{mean<Va>(Va{}),variance<Va>(Va{}), count{0}}}{}
    Moment_statistics(Va x):base_type{var::Vector_Space{mean<Va>(x),variance<Va>(x-x), count{1}}}{}
    Moment_statistics(mean<Va> x,variance<Va> v,count n):base_type{var::Vector_Space{x,v,n}}{}
    
   friend  Moment_statistics operator&(const Moment_statistics& one, const Moment_statistics& other) {
        double n0=get<count>((one)())();
        double n1 =get<count>((other)())();
        
        auto m0=get<mean<Va>>((one)())()();
        auto m1 =get<mean<Va>>((other)())()();
        
        auto v0=get<variance<Va>>((one)())()();
        auto v1 =get<variance<Va>>((other)())()();
        double n=n0+n1;
        
        auto m= m0+ n1/n *(m1-m0);
        
        auto malt= (m0*n0 + m1*n1)/n;
        
      //  auto v = n>1?v0 + (n1-1)/(n-1)* (v1-v0) + n0/(n-1) * sqr_X(m0-m) + n1/(n-1) * sqr_X(m1-m):0;
        auto v = n>1? ((n0-1)* v0 + n0* sqr_X(m0) + (n1-1)* v1 + n1* sqr_X(m1) - n * sqr_X(m))/(n-1):0;          
        Moment_statistics out;
        get<mean<Va>>((out)())()()=m;
        get<variance<Va>>((out)())()()=v;
        get<count>((out)())()=n;
        
        
        return out;
    }
    
     Moment_statistics& operator&=( const Moment_statistics& other)
    {
         *this=  *this & other;
        return *this;
    }
    
    Moment_statistics& operator&=( const Va& x)
    {
        auto n0=get<count>((*this)())();
        
        auto m0=get<mean<Va>>((*this)())()();
        
        auto v0=get<variance<Va>>((*this)())()();
        auto n=n0+1;
        
        auto m= m0+ (x()-m0)/n;
        auto v = n>1?v0 + (sqr_X(x()-m)-v0)/n0 + sqr_X(m0-m):0;
       // auto v = n>1? ((n0-1)* v0 + n0* sqr_X(m0) +   sqr_X(x()) - n * sqr_X(m))/(n-1):0;          
        
        get<count>((*this)())()=n;
        get<variance<Va>>((*this)())()()=std::move(v);
        get<mean<Va>>((*this)())()()=std::move(m);
        return *this;
          
    }
    void reset()
    {
        get<count>((*this)())()=0;
        get<mean<Va>>((*this)())()=Va{};
        get<variance<Va>>((*this)())()=Va{};
    }
    
    
  friend  Moment_statistics operator+(const Moment_statistics& one, const Moment_statistics& other) {
      
      auto n0=get<count>((one)())();
      auto m0=get<mean<Va>>((one)())()();
      auto v0=get<variance<Va>>((one)())()();
      auto n1=get<count>((other)())();
      auto m1=get<mean<Va>>((other)())()();
      auto v1=get<variance<Va>>((other)())()();
      return Moment_statistics(
              mean<Va>(Va(m0+m1)),
            variance<Va>(Va(v0+v1)),
            count(std::min(n0,n1) ));
  }
      
      friend Moment_statistics operator*(double a, Moment_statistics const& one )
      {
          auto n0=get<count>((one)())();
          auto m0=get<mean<Va>>((one)())()();
          auto v0=get<variance<Va>>((one)())()();
          return Moment_statistics(
              mean<Va>(Va(m0*a)),
              variance<Va>(Va(v0*a*a)),
              count(n0 ));
      }
};




class emcee_Step_statistics
    : public var::Constant<emcee_Step_statistics, Trial_statistics> {
    
};

class Thermo_Jump_statistics
    : public var::Constant<Thermo_Jump_statistics, Trial_statistics> {};



class logL_statistics : public var::Constant<logL_statistics, Moment_statistics<logL>>{
 public:
    using var::Constant<logL_statistics, Moment_statistics<logL>>::Constant;
    logL_statistics(logL t_logL): var::Constant<logL_statistics, Moment_statistics<logL>>{Moment_statistics<logL>(t_logL)}{}
    logL_statistics():var::Constant<logL_statistics, Moment_statistics<logL>>{Moment_statistics<logL>{}}{}
    
    
    friend logL_statistics operator+(const logL_statistics& one, const logL_statistics& two)
    {
        return logL_statistics(one()+two());
    }
    friend logL_statistics operator*(double a, const logL_statistics& one)
    {
        return logL_statistics(a * one());
    }
    
};


template<class Parameters>
struct mcmc {
    Parameters parameter;
    double logP;
    logLs logL;
  };
  
  
  template<class Parameters>
  struct error_log
  {
      std::string error;
      Parameters parameter;
  };
  
  

  
  
  
template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<mt_64i &>(), std::declval<Prior&>()))>>
 //   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_mcmc(FunctionTable& f, mt_64i &mt, Prior const & pr, const Lik& lik,
               const DataType &y, const Variables &x) {
    auto& priorsampler=pr;
    auto par = sample(mt,priorsampler);
    auto logP = logPrior(pr,par);
    auto t_logLs = logLikelihood(f,lik,par.to_value(), y,x);
    while(!(logP)||!(t_logLs))
    {
        par = sample(mt,priorsampler);
        logP = logPrior(pr,par);
        t_logLs = logLikelihood(f,lik,par.to_value(), y,x);
        
    }
    return mcmc<Parameters>{std::move(par), logP.value(), t_logLs.value()};
}




template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<mt_64i &>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
Maybe_error<bool> calc_mcmc(FunctionTable& f, Prior const & pr, const Lik& lik,
                 const DataType &y, const Variables &x, mcmc<Parameters>& t_mcmc) {
    auto par = t_mcmc.parameter;
    auto logP = logPrior(pr,par);
    auto t_logLs = logLikelihood(f,lik,par.to_value(), y,x);
    
    if (logP.valid()&& t_logLs.valid())
    {
        t_mcmc.logP=logP.value();
        t_mcmc.logL=t_logLs.value();
        return true;
    }
    else
    {
        t_mcmc.logP=std::numeric_limits<double>::quiet_NaN();
        get<logL>(t_mcmc.logL)()=std::numeric_limits<double>::quiet_NaN();;
    }
    return error_message(logP.error()()+t_logLs.error()());
}

#endif // MCMC_H
