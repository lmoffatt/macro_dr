#ifndef CUEVI_H
#define CUEVI_H
#include "bayesian_linear_regression.h"
#include "experiment.h"
#include "fold.h"
#include "function_measure_verification_and_optimization.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "random_samplers.h"
#include "variables.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <istream>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

class Save_Parameter_every
    : public var::Var<Save_Parameter_every, std::size_t> {};
class Save_Predictions_every
    : public var::Var<Save_Predictions_every, std::size_t> {};
class Save_Likelihood_every
    : public var::Var<Save_Likelihood_every, std::size_t> {};
class Save_Evidence_every
    : public var::Var<Save_Likelihood_every, std::size_t> {};

class Saving_intervals
    : public var::Var<
          Saving_intervals,
          var::Vector_Space<Save_Evidence_every, Save_Likelihood_every,
                            Save_Parameter_every, Save_Predictions_every>> {};

inline void report_model(save_Evidence &, ...) {}

using DataIndexes = std::vector<std::size_t>;

inline auto
generate_random_Indexes(mt_64i &mt, std::size_t num_samples,
                        std::size_t min_num_extra_samples,
                        double num_jumps_per_decade,
                        std::vector<std::size_t> initial_samples = {}) {
    
    std::size_t num_initial_samples = size(initial_samples);
    std::size_t n_jumps = std::max(
        0.0,
        std::floor(num_jumps_per_decade *
                   (std::log10(num_samples) -
                    std::log10(min_num_extra_samples + num_initial_samples))));
    auto indexsizes = DataIndexes(n_jumps + 1);
    
    for (std::size_t i = 0; i < n_jumps + 1; ++i)
        indexsizes[i] = num_samples * std::pow(10.0, -(1.0 * (n_jumps - i)) /
                                                         num_jumps_per_decade);
    auto out = std::vector<DataIndexes>(n_jumps + 1);
    if (n_jumps > 0) {
        auto index = DataIndexes(num_samples);
        std::iota(index.begin(), index.end(), 0u);
        auto it = index.begin();
        if (num_initial_samples > 0) {
            auto new_index = DataIndexes{};
            std::copy(initial_samples.begin(), initial_samples.end(),
                      std::back_inserter(new_index));
            std::set_difference(index.begin(), index.end(), initial_samples.begin(),
                                initial_samples.end(), std::back_inserter(new_index));
            
            std::swap(index, new_index);
            it = index.begin();
            std::advance(it, initial_samples.size());
        }
        it = randomly_extract_n(mt, it, index.end(),
                                indexsizes[0] - num_initial_samples);
        auto res = DataIndexes(index.begin(), it);
        std::sort(res.begin(), res.end());
        out[0] = std::move(res);
        for (auto i = 1u; i < n_jumps + 1; ++i) {
            auto n = (indexsizes[i] - indexsizes[i - 1]);
            it = randomly_extract_n(mt, it, index.end(), n);
            std::sort(index.begin(), it);
            out[i] = DataIndexes(index.begin(), it);
        }
    }
    return out;
}



namespace cuevi {

class Th_Beta : public var::Var<Th_Beta, double> {};
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

class LogPrior : public var::Var<LogPrior, double> {};

class LogLik_value : public var::Var<LogLik_value, double> {
public:
  friend std::ostream &operator<<(std::ostream &os, const LogLik_value &x) {
    return os << x();
  }
};

class Fraction_Index : public var::Constant<Fraction_Index, std::size_t> {
public:
  Fraction_Index friend operator+(Fraction_Index one, std::size_t i) {
    return Fraction_Index(one() + i);
  }
  Fraction_Index friend operator-(Fraction_Index one, std::size_t i) {
    return Fraction_Index(one() - i);
  }
  Fraction_Index &operator++() {
    ++((*this)());
    return *this;
  }
};

class Cuevi_Index : public var::Constant<Cuevi_Index, std::size_t> {
public:
  Cuevi_Index(std::size_t i) : var::Constant<Cuevi_Index, std::size_t>(i) {}
  Cuevi_Index friend operator+(Cuevi_Index one, std::size_t i) {
    return Cuevi_Index(one() + i);
  }
  Cuevi_Index friend operator-(Cuevi_Index one, std::size_t i) {
    return Cuevi_Index(one() - i);
  }
  Cuevi_Index &operator=(std::size_t i) {
    (*this)() = i;
    return *this;
  }
  Cuevi_Index &operator++() {
    ++(*this)();
    return *this;
  }
  friend bool operator<(const Cuevi_Index &one, std::size_t n) {
    return one() < n;
  }
};

class Number_of_Fractions
    : public var::Constant<Number_of_Fractions, std::size_t> {};

class Number_of_samples : public var::Constant<Number_of_samples, std::size_t> {
};

class Walker_id : public var::Constant<Walker_id, std::size_t> {};

class Walker_Index : public var::Constant<Walker_Index, std::size_t> {
public:
  Walker_Index &operator=(std::size_t i) {
    (*this)() = i;
    return *this;
  }
  Walker_Index &operator++() {
    ++(*this)();
    return *this;
  }
  friend bool operator<(const Walker_Index &one, std::size_t n) {
    return one() < n;
  }
};

class LogLik_by_Fraction
    : public var::Var<LogLik_by_Fraction,
                      std::map<Fraction_Index, LogLik_value>> {
public:
  bool has(Fraction_Index i_fra) const {
    return (*this)().find(i_fra) != (*this)().end();
  }

  Maybe_error<double> operator[](Fraction_Index i) const {
    auto it = (*this)().find(i);
    if (it != (*this)().end())
      return it->second();
    else
      return error_message("");
  }
};

class Num_Walkers_Per_Ensemble
    : public var::Constant<Num_Walkers_Per_Ensemble, std::size_t> {};

class Points_per_decade : public var::Constant<Points_per_decade, double> {};
class Points_per_decade_low
    : public var::Constant<Points_per_decade_low, double> {};
class Min_value : public var::Constant<Min_value, double> {};
class Med_value : public var::Constant<Med_value, double> {};
class Includes_zero : public var::Constant<Includes_zero, bool> {};

class Number_trials_until_give_up
    : public var::Constant<Number_trials_until_give_up, std::size_t> {};

class Random_jumps : public var::Constant<Random_jumps, bool> {};

class Thermo_Jumps_every
    : public var::Constant<Thermo_Jumps_every, std::size_t> {

  friend void report_model(save_Evidence &s, Thermo_Jumps_every n) {
    std::ofstream f(s.fname + "_thermo_jumps_every");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << n()
      << "\n";
  }
};

class Th_Beta_Param
    : public var::Constant<
          Th_Beta_Param,
          var::Vector_Space<Includes_zero, Med_value, Points_per_decade,
                            Min_value, Points_per_decade_low>> {
public:
  using var::Constant<
      Th_Beta_Param,
      var::Vector_Space<Includes_zero, Med_value, Points_per_decade, Min_value,
                        Points_per_decade_low>>::Constant;
};
class Fractions_Param
    : public var::Constant<Fractions_Param,
                           var::Vector_Space<Min_value, Points_per_decade>> {};

class Cuevi_temperatures_Param
    : public var::Constant<Cuevi_temperatures_Param,
                           var::Vector_Space<Th_Beta_Param, Fractions_Param>> {
};

class Prior_statistics
    : public var::Constant<Prior_statistics, Trial_statistics> {};

class Likelihood_statistics
    : public var::Constant<Likelihood_statistics,
                           std::map<Fraction_Index, Trial_statistics>> {
public:
  Likelihood_statistics &operator+=(const Likelihood_statistics &other) {
    for (auto &e : other()) {
      (*this)()[e.first] += e.second;
    }
    return *this;
  }
  
  Trial_statistics operator[](Fraction_Index i) const {
      auto it = (*this)().find(i);
      if (it != (*this)().end())
          return it->second;
      else
          return Trial_statistics{};
  }
};

class emcee_Step_statistics
    : public var::Constant<emcee_Step_statistics, Trial_statistics> {};

class Thermo_Jump_statistics
    : public var::Constant<Thermo_Jump_statistics, Trial_statistics> {};

class Cuevi_Jump_statistics
    : public var::Constant<Cuevi_Jump_statistics, Trial_statistics> {};

class Walker_statistics
    : public var::Var<
          Walker_statistics,
          var::Vector_Space<Prior_statistics, Likelihood_statistics,
                            emcee_Step_statistics, Thermo_Jump_statistics,
                            Cuevi_Jump_statistics>> {
public:
    void reset() {
        get<Prior_statistics>((*this)())().reset();
        get<Likelihood_statistics>((*this)())().clear();
        get<emcee_Step_statistics>((*this)())().reset();
        get<Thermo_Jump_statistics>((*this)())().reset();
        get<Cuevi_Jump_statistics>((*this)())().reset();
    }
    
    Walker_statistics &operator+=(const Walker_statistics &other) {
        get<Prior_statistics>((*this)())()+=get<Prior_statistics>(other())();
        get<Likelihood_statistics>((*this)())+=get<Likelihood_statistics>(other());
        get<emcee_Step_statistics>((*this)())()+=get<emcee_Step_statistics>(other())();
        get<Thermo_Jump_statistics>((*this)())()+=get<Thermo_Jump_statistics>(other())();
        get<Cuevi_Jump_statistics>((*this)())()+=get<Cuevi_Jump_statistics>(other())();
        return *this;
       }
    
};

using Walker_statistics_pair =
    std::pair<Walker_statistics &, Walker_statistics &>;

class Init_seed
    : public var::Constant<Init_seed, typename mt_64i::result_type> {};
template <class T> using by_fraction = std::vector<T>;

class Cuevi_temperatures
    : public var::Var<Cuevi_temperatures,
                      std::vector<var::Vector_Space<Th_Beta, Fraction_Index>>> {
};

class Parameter {};
class step_stretch_cuevi_mcmc;
class step_stretch_cuevi_mcmc_per_walker;
class thermo_cuevi_jump_mcmc;



class Cuevi_statistics
    : public var::Var<Cuevi_statistics,
                      ensemble<std::vector<Walker_statistics>>> {
    
public:
    auto operator[](Cuevi_Index icu)
    {
        Walker_statistics out={};
        for (std::size_t i=0; i<(*this)().size(); ++i)
        {
            out+=(*this)()[i][icu()];
        }
        return out;
    }
    
};



template <class ParameterType> class Cuevi_mcmc;

template <class ParameterType> class Cuevi_mcmc {

  using myParameter = var::Var<Parameter, ParameterType>;

  class Walker_value
      : public var::Var<Walker_value,
                        var::Vector_Space<var::Var<Parameter, ParameterType>,
                                          LogPrior, LogLik_by_Fraction>> {
      friend Maybe_error<double> thermo_step(const Walker_value &candidate,
                            const Walker_value &current, Th_Beta beta,
                            Fraction_Index i_fra) {
        
      auto th_ca= get<LogPrior>(candidate())()+
                     beta() * get<LogLik_by_Fraction>(candidate())[i_fra] ;
      auto th_cu =   get<LogPrior>(current())() +
                   beta() *  get<LogLik_by_Fraction>(current())[i_fra];
      
      if (!th_cu.valid()&& th_ca.valid())
          return 100.0;
      else
          return th_ca-th_cu;
    }
    
    friend Maybe_error<double> thermo_jump(Th_Beta ca_beta, Fraction_Index ca_fra,
                            const Walker_value &candidate, Th_Beta cu_beta,
                            Fraction_Index cu_fra,
                            const Walker_value &current) {
        
        auto logL_ca_ca=get<LogLik_by_Fraction>(candidate())[ca_fra];
        auto logL_cu_cu=get<LogLik_by_Fraction>(current())[cu_fra];
        
        auto logL_ca_cu=get<LogLik_by_Fraction>(candidate())[cu_fra];
        auto logL_cu_ca=get<LogLik_by_Fraction>(current())[ca_fra];
                 
        
      auto current_sum =
          ca_beta() * logL_ca_ca +
          cu_beta() * logL_cu_cu;
      auto after_jump_sum =
          cu_beta() * logL_ca_cu +
          ca_beta() * logL_cu_ca;
      
      bool current_is_big_0_is_NaN=(cu_fra()>ca_fra())&&(!logL_cu_ca.valid())&&(logL_cu_ca.valid());
      bool candidate_is_big_0_is_NaN=((cu_fra()<ca_fra())&&(logL_cu_ca.valid())&&(!logL_cu_ca.valid()));
      bool current_sum_is_NaN_after_is_not=(!current_sum.valid()&&after_jump_sum.valid());
      
      if (current_is_big_0_is_NaN || candidate_is_big_0_is_NaN || current_sum_is_NaN_after_is_not)
          return 0.0; // takes exponential---> prob 1
      else      
          return after_jump_sum - current_sum;
    }

    friend Maybe_error<double>
    logEvidence_walker_pos(const Walker_value &w, Cuevi_Index i,
                           const Cuevi_temperatures &t) {
      if (i() + 1 == t().size())
        return 0.0;
      else {
        auto fra_0 = get<Fraction_Index>(t()[i()]);
        auto fra_1 = get<Fraction_Index>(t()[i() + 1]);
        if (fra_0() != fra_1()) {
          auto logL_0 = get<LogLik_by_Fraction>(w())[fra_0];
          auto logL_1 = get<LogLik_by_Fraction>(w())[fra_1];
          return (logL_1 - logL_0) * 0.5;
        } else {
          auto logL_1 = get<LogLik_by_Fraction>(w())[fra_1];
          auto beta_0 = get<Th_Beta>(t()[i()])();
          auto beta_1 = get<Th_Beta>(t()[i() + 1])();
          return (beta_1 - beta_0) * 0.5 * logL_1;
        }
      }
    }

    friend auto logEvidence_walker(const Walker_value &w, Cuevi_Index i,
                                   const Cuevi_temperatures &t) {
      if (i() == 0)
        return logEvidence_walker_pos(w, i, t);
      else
        return logEvidence_walker_pos(w, i, t) +
               logEvidence_walker_pos(w, i - 1, t);
    }

    friend auto cuevi_jump(const Walker_value &w_0, Cuevi_Index i_0,
                           const Walker_value &w_1, Cuevi_Index i_1,
                           const Cuevi_temperatures &t) {
      auto current_sum =
          logEvidence_walker(w_0, i_0, t) + logEvidence_walker(w_1, i_1, t);
      auto after_jump_sum =
          logEvidence_walker(w_0, i_1, t) + logEvidence_walker(w_1, i_0, t);
      return (after_jump_sum - current_sum);
    }
  };

  class Walker
      : public var::Var<Walker, var::Vector_Space<Walker_id, Walker_value,
                                                  Walker_statistics>> {};

  class Walkers_ensemble
      : public var::Var<Walkers_ensemble, ensemble<std::vector<Walker>>> {};
  
  
  template <class FunctionTable, class t_logLikelihood, class Data,
            class Variables>
  static Maybe_error<bool>
  calc_Likelihood(FunctionTable &f, t_logLikelihood &&lik, Walker_value &w,
                  const by_fraction<Data> &y, const by_fraction<Variables> &x,
                  Fraction_Index i_frac, Walker_statistics_pair wa_sta) {
    assert(i_frac() < size(y));

    auto &r_logLikf = get<LogLik_by_Fraction>(w());
    if (r_logLikf.has(i_frac))
      return true;
    else {
      auto &ca_par = get<Parameter>(w());
      auto v_logL = f.f(logLikelihood_f{}, lik, ca_par(),
                                  y[i_frac()], x[i_frac()]);
      if (!v_logL) {
        fails(get<Likelihood_statistics>(wa_sta.first())()[i_frac]);
        fails(get<Likelihood_statistics>(wa_sta.second())()[i_frac]);
        //r_logLikf()[i_frac] = LogLik_value(-std::numeric_limits<double>::infinity());
        return v_logL.error();

      } else {
        succeeds(get<Likelihood_statistics>(wa_sta.first())()[i_frac]);
        succeeds(get<Likelihood_statistics>(wa_sta.second())()[i_frac]);
        r_logLikf()[i_frac] = LogLik_value(v_logL.value());
        return true;
      }
    }
  }
  
  
  
  
  template <class FunctionTable, class t_logLikelihood, class Data,
            class Variables>
  static Maybe_error<bool> calc_Relevant_Likelihoods(
      FunctionTable &f, t_logLikelihood &&lik, const by_fraction<Data> &y,
      const by_fraction<Variables> &x, const Cuevi_temperatures &t,
      Walker_value &w, Cuevi_Index i, Walker_statistics_pair wa_sta) {

    Maybe_error<bool> out = true;
    for (auto i_frac = get<Fraction_Index>(t()[std::max(1ul, i()) - 1]);
         i_frac < get<Fraction_Index>(t()[std::min(t().size() - 1, i() + 1)]);
         ++i_frac) {
      auto Maybe = calc_Likelihood(f, lik, w, y, x,
                                   i_frac, wa_sta);
      if (!Maybe)
        return Maybe.error();
    }
    return true;
  }
  
  
 
  
  template <class FunctionTable, class logLikelihood, class Data,
            class Variables>
  static Maybe_error<double> thermo_jump_logProb(
      FunctionTable &f, logLikelihood &&lik, const by_fraction<Data> &y,
      const by_fraction<Variables> &x, const Cuevi_temperatures &t,
      Cuevi_Index i_0, Walker_value &w_0, Cuevi_Index i_1, Walker_value &w_1,
      Walker_statistics_pair wa_sta_0, Walker_statistics_pair wa_sta_1) {
    auto i_frac_0 = get<Fraction_Index>(t()[i_0()]);
    auto i_frac_1 = get<Fraction_Index>(t()[i_1()]);

    auto Maybe_0 = calc_Likelihood(f, lik, w_0, y,
                                   x, i_frac_1, wa_sta_0);
    auto Maybe_1 = calc_Likelihood(f, lik, w_1, y,
                                   x, i_frac_0, wa_sta_1);
    auto be_0 = get<Th_Beta>(t()[i_0()]);
    auto be_1 = get<Th_Beta>(t()[i_1()]);

    // if (!Maybe_0.valid() || !Maybe_1.valid())
    //   return error_message(Maybe_0.error()() + Maybe_1.error()());
    // else {
      return thermo_jump(be_0, i_frac_0, w_0, be_1, i_frac_1, w_1);
   // }
  }

  template <class FunctionTable, class logLikelihood, class Data,
            class Variables>
  static Maybe_error<double> cuevi_jump_logProb(
      FunctionTable &f, logLikelihood &&lik, const by_fraction<Data> &y,
      const by_fraction<Variables> &x, const Cuevi_temperatures &t,
      Cuevi_Index i_0, Walker_value &w_0, Cuevi_Index i_1, Walker_value &w_1,
      Walker_statistics_pair wa_sta_0, Walker_statistics_pair wa_sta_1) {

    auto Maybe_00 = calc_Relevant_Likelihoods(f,
                                              lik, y, x, t, w_0, i_0, wa_sta_0);
    auto Maybe_01 = calc_Relevant_Likelihoods(f,
                                              lik, y, x, t, w_0, i_1, wa_sta_0);
    auto Maybe_10 = calc_Relevant_Likelihoods(f,
                                              lik, y, x, t, w_1, i_0, wa_sta_1);
    auto Maybe_11 = calc_Relevant_Likelihoods(f,
                                              lik, y, x, t, w_1, i_1, wa_sta_1);

    if (!Maybe_00.valid() || !Maybe_01.valid() || !Maybe_10.valid() ||
        !Maybe_11.valid())
      return error_message(Maybe_00.error()() + Maybe_01.error()() +
                           Maybe_10.error()() + Maybe_11.error()());
    else {
      return cuevi_jump(w_0, i_0, w_1, i_1, t);
    }
  }

  template <class FunctionTable, class Prior, class logLikelihood, class Data,
            class Variables>
  static Maybe_error<Walker_value>
  calc_Walker_value(FunctionTable &f, myParameter &&ca_par, Prior &&p,
                    logLikelihood &&lik, const by_fraction<Data> &y,
                    const by_fraction<Variables> &x,
                    const Cuevi_temperatures &t, Cuevi_Index i_cu,
                    Walker_statistics_pair wa_sta) {
    auto v_logP = logPrior(p, ca_par());
    if (!v_logP) {
      fails(get<Prior_statistics>(wa_sta.first())());
      fails(get<Prior_statistics>(wa_sta.second())());
      return v_logP.error();
    } else {
      succeeds(get<Prior_statistics>(wa_sta.first())());
      succeeds(get<Prior_statistics>(wa_sta.second())());
      Walker_value out(var::Vector_Space(
          std::move(ca_par), LogPrior(v_logP.value()), LogLik_by_Fraction{}));

      auto i_frac = get<Fraction_Index>(t()[i_cu()]);
      auto Maybe_succeed = calc_Likelihood(f, lik,
                                           out, y, x, i_frac, wa_sta);
      if (!Maybe_succeed)
        return Maybe_succeed.error();
      else
        return out;
    }
  }

  template <class FunctionTable, class Sampler, class Prior,
            class logLikelihood, class Data, class Variables>
  static Maybe_error<Walker_value>
  sample_Walker_for(FunctionTable &f, mt_64i &mt, Sampler &&sampler, Prior &&p,
                    logLikelihood &&lik, const by_fraction<Data> &y,
                    const by_fraction<Variables> &x,
                    const Cuevi_temperatures &t, Cuevi_Index i_cu,
                    Walker_statistics_pair wa_sta,
                    Number_trials_until_give_up max_trials) {
    assert(i_cu() < t().size());

    Maybe_error<Walker_value> v_walker(error_message{});
    auto n_trial = 0ul;
    while (!v_walker && n_trial < max_trials) {
      auto ca_par = std::forward<Sampler>(sampler)(mt);
      v_walker = calc_Walker_value(f, ca_par, p,
                                   lik, y, x, t, i_cu, wa_sta);
      ++n_trial;
    }
    if (v_walker.valid())
      return v_walker.value();
    else
      return error_message("more than " + std::to_string(max_trials()) +
                           " and not a single valid sample, last error " +
                           v_walker.error()());
  }

  template <class FunctionTable, class Prior, class logLikelihood, class Data,
            class Variables>
  static auto sample_Walker(FunctionTable &f, mt_64i &mt, Prior &&p,
                            logLikelihood &&lik, const by_fraction<Data> &y,
                            const by_fraction<Variables> &x,
                            const Cuevi_temperatures &t, Cuevi_Index i_cu,
                            Walker_statistics_pair w_sta,
                            Number_trials_until_give_up max_trials) {
    return sample_Walker_for(
        f, mt,
        [&p](mt_64i &mt) { return sample(mt, std::forward<Prior>(p)); }, p, lik,
        y, x, t, i_cu, w_sta, max_trials);
  }

  template <class FunctionTable, class Prior, class logLikelihood, class Data,
            class Variables>
  Maybe_error<Walkers_ensemble> friend init_cuevi(
      FunctionTable &f, ensemble<mt_64i> &mts, Prior &&p, logLikelihood &&lik,
      const by_fraction<Data> &y, const by_fraction<Variables> &x,
      const Cuevi_temperatures &t, Num_Walkers_Per_Ensemble n,
      Cuevi_statistics &sta, std::size_t max_trials_per_sample) {
    auto num_temp = t().size();

    Walkers_ensemble out(n(), std::vector<Walker>(num_temp));

    Maybe_error<bool> succeeds = true;
    for (std::size_t half = 0; half < 2; ++half)
#pragma omp parallel for
      for (std::size_t iiw = 0; iiw < n() / 2; ++iiw) {
        auto iw = half ? iiw + n() / 2 : iiw;
        for (std::size_t i_cu = 0; i_cu < num_temp; ++i_cu) {
          auto &wa_va = out()[iw]()[i_cu];
          Walker_statistics_pair wa_sta(get<Walker_statistics>(wa_va()),
                                        sta()[iw][i_cu]);
          get<Walker_id>(wa_va())() = iw + num_temp * i_cu;
          auto Maybe_Walker_value =
              sample_Walker(f.fork(var::I_thread(iiw)), mts[iiw], p, lik, y, x,
                            t, i_cu, wa_sta, max_trials_per_sample);
          if (Maybe_Walker_value)
            get<Walker_value>(wa_va()) = std::move(Maybe_Walker_value.value());
          else
            succeeds = Maybe_Walker_value.error();
        }
      }
    if (succeeds)
      return out;
    else
      return succeeds.error();
  }

  template <class FunctionTable, class Prior, class logLikelihood, class Data,
            class Variables>
  static Maybe_error<typename Cuevi_mcmc<ParameterType>::Walkers_ensemble>
  init_walkers_ensemble(FunctionTable &f, ensemble<mt_64i> &mts, Prior &&p,
                        logLikelihood &&lik, const by_fraction<Data> &y,
                        const by_fraction<Variables> &x,
                        const Cuevi_temperatures &t, Num_Walkers_Per_Ensemble n,
                        Number_trials_until_give_up max_trials_per_sample) {
    auto num_temp = t().size();
    using Walkers_ensemble =
        typename Cuevi_mcmc<ParameterType>::Walkers_ensemble;
    using Walker = typename Cuevi_mcmc<ParameterType>::Walker;
    using Walker_value = typename Cuevi_mcmc<ParameterType>::Walker_value;
    Cuevi_statistics sta(ensemble<std::vector<Walker_statistics>>(
        n(), std::vector<Walker_statistics>(num_temp)));
    ;
    Walkers_ensemble out(std::vector(n(), std::vector<Walker>(num_temp)));
    
    auto ff=f.fork(n()/2);
    Maybe_error<bool> succeeds = true;
    for (std::size_t half = 0; half < 2; ++half)
      //#pragma omp parallel for
      for (std::size_t iiw = 0; iiw < n() / 2; ++iiw) {
        auto iw = half ? iiw + n() / 2 : iiw;
        for (std::size_t i_cu = 0; i_cu < num_temp; ++i_cu) {
          auto &wa_va = out()[iw][i_cu];
          Walker_statistics_pair wa_sta(get<Walker_statistics>(wa_va()),
                                        sta()[iw][i_cu]);
          get<Walker_id>(wa_va())() = iw + num_temp * i_cu;
          auto Maybe_Walker_value = sample_Walker(
              ff[iiw], mts[iiw], std::forward<Prior>(p), lik,
              y, x, t, Cuevi_Index(i_cu), wa_sta, max_trials_per_sample);
          if (Maybe_Walker_value)
            get<Walker_value>(wa_va()) = std::move(Maybe_Walker_value.value());
          else
            succeeds = Maybe_Walker_value.error();
        }
      }
    f+=ff;
    if (succeeds)
      return out;
    else
      return succeeds.error();
  }

  Cuevi_temperatures m_temperatures;
  Cuevi_statistics m_sta;
  Walkers_ensemble m_data;
  std::size_t m_max_i_frac;

  Cuevi_mcmc(Cuevi_temperatures &&t_temperatures, Cuevi_statistics &&t_sta,
             Walkers_ensemble &&t_data, std::size_t t_max_i_frac)
      : m_temperatures{std::move(t_temperatures)}, m_sta{std::move(t_sta)},
        m_data{std::move(t_data)}, m_max_i_frac{t_max_i_frac} {}

public:
  auto calc_Mean_logLik(Cuevi_Index j) {
    assert(j() < get_Cuevi_Temperatures_Number());
    auto i_frac = get_Fraction(j);
    return foldMap(
               make_Range(Walker_Index(0ul),
                          Walker_Index(get_Walkers_number())),
               [j, i_frac, this](auto i_w) {
                 return get<LogLik_by_Fraction>(
                     get_Walker_Value(i_w, j)())[i_frac];
               },
               [](auto x, auto y) { return x + y; }) /
           get_Walkers_number();
  }

  Maybe_error<double> calc_Mean_logLik_0(Cuevi_Index j) {
    assert(j() < get_Cuevi_Temperatures_Number());

    auto i_frac = get_Fraction(j);
    if (i_frac() == 0)
      return error_message();
    else
      return foldMap(
                 make_Range(Walker_Index(0ul),
                            Walker_Index(get_Walkers_number())),
                 [j, i_frac, this](auto i_w) {
                   return get<LogLik_by_Fraction>(
                       get_Walker_Value(i_w, j)())[i_frac - 1];
                 },
                 [](auto x, auto y) { return x + y; }) /
             get_Walkers_number();
  }

  Maybe_error<double> calc_Mean_logLik_2(Cuevi_Index j) {
    assert(j() < get_Cuevi_Temperatures_Number());
    auto i_frac = get_Fraction(j);
    auto beta = get_Beta(j);
    if (beta < 1.0)
      return error_message();
    else if (i_frac() + 1 >= m_max_i_frac)
      return error_message();
    else
      return foldMap(
                 make_Range(Walker_Index(0ul),
                            Walker_Index(get_Walkers_number())),
                 [j, i_frac, this](auto i_w) {
                   return get<LogLik_by_Fraction>(
                       get_Walker_Value(i_w, j)())[i_frac + 1];
                 },
                 [](auto x, auto y) { return x + y; }) /
             get_Walkers_number();
  }

  auto calc_Mean_logPrior(Cuevi_Index j) {
    assert(j() < get_Cuevi_Temperatures_Number());
    return foldMap(
               make_Range(Walker_Index(0ul),
                          Walker_Index(get_Walkers_number())),
               [j, this](auto i_w) {
                 return get<LogPrior>(get_Walker_Value(i_w, j)())();
               },
               [](auto x, auto y) { return x + y; }) /
           get_Walkers_number();
  }

  auto &get_Walker_Value(Walker_Index i, Cuevi_Index j) {
    assert(i() < get_Walkers_number());
    assert(j() < get_Cuevi_Temperatures_Number());
    return get<Walker_value>(m_data()[i()][j()]());
  }
  auto &get_Walker_Value(Walker_Index i, Cuevi_Index j) const {
    assert(i() < get_Walkers_number());
    assert(j() < get_Cuevi_Temperatures_Number());
    return get<Walker_value>(m_data()[i()][j()]());
  }
  Walker_statistics_pair get_Walker_Statistics(Walker_Index i, Cuevi_Index j) {
    assert(i() < get_Walkers_number());
    assert(j() < get_Cuevi_Temperatures_Number());
    return Walker_statistics_pair(get<Walker_statistics>(m_data()[i()][j()]()),
                                  m_sta()[i()][j()]);
  }

  auto &get_Cuevi_Temperatures() const { return m_temperatures; }

  auto get_Cuevi_Temperatures_Number() const { return m_temperatures().size(); }

  auto get_Walkers_number() const { return m_data().size(); }
  auto get_Parameters_number() const {
    return size(get<Parameter>(get<Walker_value>(m_data()[0][0]())())());
  }
  auto get_Parameter(Walker_Index i, Cuevi_Index j) const {
    assert(i() < get_Walkers_number());
    assert(j() < get_Cuevi_Temperatures_Number());
    return get<Parameter>(get_Walker_Value(i, j)())();
  }
  auto &get_Walker(Walker_Index i, Cuevi_Index j) {
    assert(i() < get_Walkers_number());
    assert(j() < get_Cuevi_Temperatures_Number());

    return m_data()[i()][j()];
  }
  auto &get_Walker(Walker_Index i, Cuevi_Index j) const {
    return m_data()[i()][j()];
  }

  auto get_Beta(Cuevi_Index i_cu) {
    assert(i_cu() < get_Cuevi_Temperatures_Number());
    return get<Th_Beta>(m_temperatures()[i_cu()]);
  }
  
  auto get_Cuevi_Statistics(Cuevi_Index i_cu){
      return m_sta[i_cu];
  }
  
  auto get_Fraction(Cuevi_Index i_cu) {
    assert(i_cu() < get_Cuevi_Temperatures_Number());

    return get<Fraction_Index>(m_temperatures()[i_cu()]);
  }

  template <class DataType>
  static Cuevi_temperatures build_temperatures(const by_fraction<DataType> &y,
                                               const Th_Beta_Param &p) {
    auto n_points_per_decade_beta = get<Points_per_decade>(p());
    auto first_stops_at = get<Med_value>(p());
    auto n_points_per_decade_beta_low = get<Points_per_decade_low>(p());
    auto stops_at = get<Min_value>(p());
    auto includes_zero = get<Includes_zero>(p());

    std::size_t num_beta_high =
        std::ceil(-std::log10(first_stops_at()) * n_points_per_decade_beta());
    std::size_t num_beta_low =
        std::ceil((std::log10(first_stops_at()) - std::log10(stops_at())) *
                  n_points_per_decade_beta_low());

    auto beta_size = 1 + num_beta_high + num_beta_low;
    if (includes_zero()) {
      num_beta_low = num_beta_low + 1;
      beta_size = beta_size + 1;
    }
    auto beta_mid =
        std::pow(10.0, -1.0 * num_beta_high / n_points_per_decade_beta());

    auto fraction_size = y.size();

    auto out = Cuevi_temperatures{};
    auto Cuevi_size = beta_size + fraction_size - 1;

    out().resize(Cuevi_size);
    if (includes_zero()) {
      get<Th_Beta>(out()[0])() = 0;
      get<Fraction_Index>(out()[0])() = 0ul;
    }

    for (std::size_t i = includes_zero() ? 1 : 0; i < num_beta_low; ++i) {
      get<Th_Beta>(out()[i])() =
          beta_mid * std::pow(10.0, -1.0 * (num_beta_low - i) /
                                        n_points_per_decade_beta_low());
      get<Fraction_Index>(out()[i])() = 0ul;
    }
    for (std::size_t i = 0; i < num_beta_high + 1; ++i) {
      get<Th_Beta>(out()[i + num_beta_low])() = std::pow(
          10.0, -1.0 * (num_beta_high - i) / n_points_per_decade_beta());
      get<Fraction_Index>(out()[i + num_beta_low])() = 0ul;
    }
    for (std::size_t i = 1; i < fraction_size; ++i) {
      get<Th_Beta>(out()[i + beta_size - 1])() = 1.0;
      get<Fraction_Index>(out()[i + beta_size - 1])() = i;
    }

    return out;
  }

  friend class step_stretch_cuevi_mcmc;
  friend class step_stretch_cuevi_mcmc_per_walker;
  friend class thermo_cuevi_jump_mcmc;
  friend class thermo_cuevi_jump_mcmc_per_walker;

  template <class FunctionTable, class Prior, class logLikelihood, class Data,
            class Variables>
  Maybe_error<Cuevi_mcmc> static init(
      FunctionTable &f, ensemble<mt_64i> &mts, Prior &&prior,
      logLikelihood &&lik, const by_fraction<Data> &y,
      const by_fraction<Variables> &x, const Th_Beta_Param &beta,
      Num_Walkers_Per_Ensemble num_walkers,
      Number_trials_until_give_up max_trials_per_sample) {
    using Cumc = Cuevi_mcmc<ParameterType>;
    auto t = Cumc::build_temperatures(y, beta);
    auto sta = Cuevi_statistics(
        std::vector(num_walkers(), std::vector<Walker_statistics>(t().size())));
    auto Maybe_walkers =
        Cumc::init_walkers_ensemble(f, mts, std::forward<Prior>(prior), lik, y,
                                    x, t, num_walkers, max_trials_per_sample);
    if (!Maybe_walkers)
      return Maybe_walkers.error();
    else
      return Cumc(std::move(t), std::move(sta),
                  std::move(Maybe_walkers.value()), size(y));
  }

  template <class FunctionTable, class t_logLikelihood, class Data,
            class Variables>
  void calculate_Likelihoods_for_Evidence_calulation(
      FunctionTable &f, t_logLikelihood &&lik, const by_fraction<Data> &y,
      const by_fraction<Variables> &x) {
      auto ff=f.fork(get_Walkers_number() / 2);
    for (std::size_t half = 0; half < 2; ++half)
      // #pragma omp parallel for
      for (std::size_t iiw = 0; iiw < get_Walkers_number() / 2; ++iiw)

      {
        Walker_Index i_w(half ? iiw + get_Walkers_number() / 2 : iiw);
        for (Cuevi_Index i_cu = Cuevi_Index(0ul);
             i_cu < this->get_Cuevi_Temperatures_Number(); ++i_cu) {
          auto i_frac = this->get_Fraction(i_cu)();
          auto beta = this->get_Beta(i_cu)();
          if (beta == 1) {
            Walker_value &wa = this->get_Walker_Value(i_w, i_cu);
            Walker_statistics_pair wa_sta =
                this->get_Walker_Statistics(i_w, i_cu);
            for (std::size_t i = std::max(1ul, i_frac) - 1;
                 i < std::min(i_frac + 2, m_max_i_frac); ++i) {
              Fraction_Index ifrac = Fraction_Index(i);
              calc_Likelihood(ff[iiw], lik, wa, y, x, ifrac,
                              wa_sta);
            }
          }
        }
      }
    f+=ff;
  }

  template <class FunctionTable, class Prior, class t_logLikelihood, class Data,
            class Variables>
  friend void
  report(FunctionTable &f, std::size_t iter, save_likelihood<ParameterType> &s,
         Cuevi_mcmc &data, Prior &&, t_logLikelihood &&lik,
         const by_fraction<Data> &y, const by_fraction<Variables> &x, ...) {

    auto &t = data.get_Cuevi_Temperatures();
    if (iter % s.save_every == 0) {
      data.calculate_Likelihoods_for_Evidence_calulation(f, lik, y, x);
      for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
           ++i_walker) {
        auto iw = Walker_Index(i_walker);
        auto &wav = data.get_Walker_Value(iw, Cuevi_Index(0ul));
        auto logL1 = get<LogLik_by_Fraction>(wav())[Fraction_Index(0)];

        using Maybe_L = decltype(logL1);
        auto beta1 = 0.0;
        Maybe_L log_Evidence = 0.0;
        Maybe_L log_Evidence_no_0 = 0.0;

        for (Cuevi_Index i_cu = 0ul;
             i_cu() < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
          auto i_fra = data.get_Fraction(i_cu);
          auto nsamples = size(y[i_fra()]);
          auto icu = Cuevi_Index(i_cu);
          auto &wa = data.get_Walker(iw, icu);
          auto &wav = data.get_Walker_Value(iw, icu);
          auto &stat = get<Walker_statistics>(wa());
          
          auto &prior_sta = get<Prior_statistics>(stat());
          auto lik_sta_2 = get<Likelihood_statistics>(stat())[i_fra + 1];
          auto lik_sta_1 = get<Likelihood_statistics>(stat())[i_fra];
          auto lik_sta_0 = get<Likelihood_statistics>(stat())[i_fra - 1];
          
          auto &emcee_sta = get<emcee_Step_statistics>(stat());
          auto &th_sta = get<Thermo_Jump_statistics>(stat());
          auto &cuevi_sta = get<Cuevi_Jump_statistics>(stat());
          
          auto logPrior = get<LogPrior>(wav());
          double beta0;
          Maybe_L logL0;
          Maybe_L plog_Evidence;
          Maybe_L logL1_1 = error_message{};
          Maybe_L logL1_0 = error_message{};
          Maybe_L logL1_2 = error_message{};
          Maybe_L logL0_1 = error_message{};
          Maybe_L logL0_0 = error_message{};
          if (i_fra() == 0) {
              beta1 = data.get_Beta(i_cu)();
            logL0 = logL1;
            beta0 = beta1;
            logL1 = get<LogLik_by_Fraction>(wav())[i_fra];
            logL1_1 = logL1;
            if (beta1 == 1.0)
                logL1_2 = get<LogLik_by_Fraction>(wav())[i_fra + 1];
            
            plog_Evidence = (beta1 - beta0) * (logL0 + logL1) / 2;
            log_Evidence = log_Evidence + plog_Evidence;
            if (beta0 == 0)
              log_Evidence_no_0 = log_Evidence_no_0 + beta1 * logL1;
            else
              log_Evidence_no_0 = log_Evidence_no_0 + plog_Evidence;

          } else {
            logL0 = logL1;
            beta0 = beta1;
            logL1_1 = get<LogLik_by_Fraction>(wav())[i_fra];
            logL1_0 = get<LogLik_by_Fraction>(wav())[i_fra - 1];
            logL1_2 = get<LogLik_by_Fraction>(wav())[i_fra + 1];

            auto &wav0 = data.get_Walker_Value(iw, icu - 1);
            logL0_1 = get<LogLik_by_Fraction>(wav0())[i_fra];
            logL0_0 = get<LogLik_by_Fraction>(wav0())[i_fra - 1];
            logL0 = logL0_1 - logL0_0;
            logL1 = logL1_1 - logL1_0;
            beta0 = 0;
            beta1 = data.get_Beta(i_cu)();
            plog_Evidence = (beta1 - beta0) * (logL0 + logL1) / 2;
            log_Evidence = log_Evidence + plog_Evidence;
            log_Evidence_no_0 = log_Evidence_no_0 + plog_Evidence;
          }
          s.f << iter << s.sep << i_cu << s.sep << i_fra() << s.sep << nsamples
              << s.sep << beta1 << s.sep << beta0 << s.sep << i_walker << s.sep
              << get<Walker_id>(wa())() << s.sep << logPrior << s.sep << logL1
              << s.sep << logL0 << s.sep << plog_Evidence << s.sep
              << log_Evidence << s.sep << log_Evidence_no_0 << s.sep << logL1_1
              << s.sep << logL1_0 << s.sep << logL1_2 << s.sep << logL0_1
              << s.sep << logL0_0 << s.sep
              
              << prior_sta().count() << s.sep << prior_sta().rate() << s.sep
              << lik_sta_0.count() << s.sep << lik_sta_0.rate() << s.sep
              << lik_sta_1.count() << s.sep << lik_sta_1.rate() << s.sep
              << lik_sta_2.count() << s.sep << lik_sta_2.rate() << s.sep
              << emcee_sta().count() << s.sep << emcee_sta().rate() << s.sep
              << th_sta().count() << s.sep << th_sta().rate() << s.sep
              << cuevi_sta().count() << s.sep << cuevi_sta().rate() << "\n";
          stat.reset();
        }
      }
    }
  }

  friend void report_title(save_likelihood<ParameterType> &s,
                           Cuevi_mcmc const &...) {

    s.f << "iter" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep << "nsamples"
        << s.sep << "beta1" << s.sep << "beta0" << s.sep << "i_walker" << s.sep
        << "walker_id" << s.sep << "logPrior" << s.sep << "logL1" << s.sep
        << "logL0" << s.sep << "plog_Evidence" << s.sep << "log_Evidence"
        << s.sep << "log_Evidence_no_0" << s.sep << "logL1_1" << s.sep
        << "logL1_0" << s.sep << "logL1_2" << s.sep << "logL0_1" << s.sep
          << "logL0_0" << s.sep
          
          << "prior_sta_count" << s.sep << "prior_sta_rate" << s.sep
          << "lik_sta_0_count" << s.sep << "lik_sta_0_rate" << s.sep
          << "lik_sta_1_count" << s.sep << "lik_sta_1_rate" << s.sep
          << "lik_sta_2_count" << s.sep << "lik_sta_2_rate" << s.sep
          << "emcee_sta_count" << s.sep << "emcee_sta_rate" << s.sep
          << "th_sta_count" << s.sep << "th_sta_rate" << s.sep
          << "cuevi_sta_count" << s.sep << "cuevi_sta_rate"
          << "\n";
  }

  template <class FunctionTable, class Prior, class t_logLikelihood, class Data,
            class Variables, class Finalizer>
  friend void
  report(FunctionTable &, std::size_t iter, save_Parameter<ParameterType> &s,
         Cuevi_mcmc const &data, Prior &&, t_logLikelihood &&,
         const by_fraction<Data> &y, const by_fraction<Variables> &,
         const std::vector<mt_64i> &mts, const Finalizer &mcmc, ...) {

    auto &t = data.get_Cuevi_Temperatures();
    if (iter % s.save_every == 0) {
      for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number();
           ++i_cu) {
        auto icu = Cuevi_Index(i_cu);
        auto i_fra = get<Fraction_Index>(t()[i_cu]);
        for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
             ++i_walker) {
          auto iw = Walker_Index(i_walker);
          auto &wa = data.get_Walker(iw, icu);
          auto &wav = data.get_Walker_Value(iw, icu);

          auto i_mts = iw() % 2;

          for (std::size_t i_par = 0; i_par < data.get_Parameters_number();
               ++i_par) {
            s.f << iter << s.sep << i_cu << s.sep << i_fra() << s.sep
                << size(y[i_fra()]) << s.sep << get<Th_Beta>(t()[i_cu])()
                << s.sep << i_walker << s.sep << get<Walker_id>(wa())() << s.sep
                << i_mts << s.sep << mts[i_mts].pos() << s.sep << i_par << s.sep
                << get<Parameter>(wav())()[i_par];
            report_finalizer_data(s, mcmc);
            s.f << "\n";
          }
        }
      }
    }
  }

  template <bool verifying, class Prior, class t_logLikelihood, class Data,
            class Variables>
  static Maybe_error<std::pair<Cuevi_mcmc, std::vector<unsigned long long>>>
  extract_iter(save_Parameter<ParameterType> &s, std::string &last_line,
               std::size_t &iter0, std::size_t &npar, std::size_t &nwalkers,
               std::size_t &n_cuevi) {
    std::size_t i_cu0;
    std::size_t i_fra0;
    std::size_t nsamples0;
    double beta0;
    std::size_t i_walker0;
    std::size_t walker_id0;

    std::size_t i_mt0;
    std::size_t mt_pos0;
    std::size_t i_par0;
    double par_value0;

    std::size_t iter;
    std::size_t i_cu;
    std::size_t i_fra;
    std::size_t nsamples;
    double beta;
    std::size_t i_walker;
    std::size_t walker_id;
    std::size_t i_mt;
    std::size_t mt_pos;
    std::size_t i_par;
    double par_value;

    if (last_line.empty()) {
      std::getline(s.f, last_line);
      if (!s.f)
        return error_message("no more iterations");
    }

    std::stringstream ss(last_line);

    ss >> iter0 >> s.sep >> i_cu0 >> s.sep >> i_fra0 >> s.sep >> nsamples0 >>
        s.sep >> beta0 >> s.sep >> i_walker0 >> s.sep >> walker_id0 >> s.sep >>
        i_mt0 >> s.sep >> mt_pos0 >> s.sep >> i_par0 >> s.sep >> par_value0;

    iter = iter0;
    i_cu = i_cu0;
    i_walker = i_walker0;
    Cuevi_temperatures t;
    t().reserve(n_cuevi);
    std::vector<std::vector<Walker>> d;
    d.reserve(n_cuevi);
    std::vector<unsigned long long> mts_pos0;
    std::vector<unsigned long long> mts_pos;
    while (iter == iter0) {

      std::vector<Walker> wav;
      mts_pos.reserve(nwalkers / 2);
      wav.reserve(nwalkers);
      while (i_cu == i_cu0) {
        std::vector<double> par;
        par.reserve(npar);
        while (i_walker == i_walker0) {
          std::getline(s.f, last_line);
          std::stringstream ss2(last_line);

          ss2 >> iter >> s.sep >> i_cu >> s.sep >> i_fra >> s.sep >> nsamples >>
              s.sep >> beta >> s.sep >> i_walker >> s.sep >> walker_id >>
              i_mt >> s.sep >> mt_pos >> s.sep >> s.sep >> i_par >> s.sep >>
              par_value;
          if constexpr (verifying) {
            if (par.size() != i_par0)
              return error_message("i_par out of order");
            if (walker_id != walker_id0)
              return error_message("walker_id out of order");
            if (i_fra != i_fra0)
              return error_message("i_fra out of order");
            if (i_cu != i_cu0)
              return error_message("i_cu out of order");
            if (iter != iter0)
              return error_message("iter out of order");
          }
          par.push_back(par_value0);
          par_value0 = par_value;
          i_par0 = i_par;
        }
        npar = par.size();
        auto wa = Walker{};
        get<Walker_id>(wa())() = walker_id0;
        get<Parameter>(get<Walker_value>(wa())())() = ParameterType(par);
        wav.push_back(std::move(wa));
        walker_id0 = walker_id;
        mts_pos.push_back(mt_pos0);
        mt_pos0 = mt_pos;
      }
      nwalkers = wav.size();
      if constexpr (verifying) {
        if (t().size() != i_cu0)
          return error_message("i_cu out of order");
      }
      t().push_back(var::Vector_Space(Th_Beta(beta0), Fraction_Index(i_fra0)));
      d.push_back(wav);
      beta0 = beta;
      i_fra0 = i_fra;
      i_cu0 = i_cu;
      if constexpr (verifying) {
        if (mts_pos0.size() > 0)
          if (mts_pos0 != mts_pos)
            return error_message("mts_pos differ beween i_cusr");
      }
      mts_pos0 = mts_pos;
    }
    n_cuevi = d.size();
    Walkers_ensemble data(
        std::vector(nwalkers, std::vector(n_cuevi, Walker{})));
    for (std::size_t iw = 0; iw < nwalkers; ++iw)
      for (std::size_t icu = 0; icu < n_cuevi; ++icu)
        data()[iw][icu] = d[icu][iw];
    Cuevi_mcmc out;
    out.m_temperatures = std::move(t);
    out.m_data = std::move(data);
    return std::pair(out, mts_pos);
  }

  template <class mcmcm_type>
  friend void report_title(save_Parameter<ParameterType> &s, Cuevi_mcmc const &,
                           const mcmcm_type &mcmc, ...) {

    s.f << "iter" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep << "nsamples"
        << s.sep << "beta" << s.sep << "i_walker" << s.sep << "walker_id"
        << s.sep << "i_mt" << s.sep << "mt_pos" << s.sep << "i_par" << s.sep
        << "parameter_value";
    report_finalizer_title(s, mcmc);
    s.f << "\n";
  }

  template <class FunctionTable, class Prior, class t_logLikelihood, class Data,
            class Variables>
  friend void report(FunctionTable &f, std::size_t iter, save_Evidence &s,
                     Cuevi_mcmc &data, Prior &&, t_logLikelihood &&lik,
                     const by_fraction<Data> &y,
                     const by_fraction<Variables> &x, ...) {

    auto &t = data.get_Cuevi_Temperatures();
    if (iter % s.save_every == 0) {
      data.calculate_Likelihoods_for_Evidence_calulation(f, lik, y, x);
      auto logL1 = data.calc_Mean_logLik(Cuevi_Index(0ul));
      auto logL1_0 = data.calc_Mean_logLik_0(Cuevi_Index(0ul));
      auto logL1_2 = data.calc_Mean_logLik_2(Cuevi_Index(0ul));

      auto beta1 = 0.0;
      Maybe_error<double> log_Evidence = 0.0;
      Maybe_error<double> log_Evidence_no_0 = 0.0;

      for (Cuevi_Index i_cu = 0ul;
           i_cu() < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
        auto i_frac = data.get_Fraction(i_cu);
        auto nsamples = size(y[i_frac()]);
        auto logPrior = data.calc_Mean_logPrior(i_cu);
        double beta0;
        Maybe_error<double> logL0;
        Maybe_error<double> plog_Evidence;
        Maybe_error<double> logL1_1 = error_message{};
        Maybe_error<double> logL1_0 = error_message{};
        Maybe_error<double> logL0_1 = error_message{};
        Maybe_error<double> logL0_0 = error_message{};
        
        auto stat = data.get_Cuevi_Statistics(i_cu);
        
        auto &prior_sta = get<Prior_statistics>(stat());
        auto lik_sta_2 = get<Likelihood_statistics>(stat())[i_frac + 1];
        auto lik_sta_1 = get<Likelihood_statistics>(stat())[i_frac];
        auto lik_sta_0 = get<Likelihood_statistics>(stat())[i_frac - 1];
        
        auto &emcee_sta = get<emcee_Step_statistics>(stat());
        auto &th_sta = get<Thermo_Jump_statistics>(stat());
        auto &cuevi_sta = get<Cuevi_Jump_statistics>(stat());
        
        
        
        if (i_frac() == 0) {
          logL0 = logL1;
          beta0 = beta1;
          logL1 = data.calc_Mean_logLik(i_cu);
          beta1 = data.get_Beta(i_cu)();
          plog_Evidence = (beta1 - beta0) * (logL0 + logL1) / 2;
          log_Evidence = log_Evidence + plog_Evidence;
          if (beta0 == 0)
            log_Evidence_no_0 = log_Evidence_no_0 + beta1 * logL1;
          else
            log_Evidence_no_0 = log_Evidence_no_0 + plog_Evidence;

        } else {
          logL0 = logL1;
          beta0 = beta1;
          logL1_1 = data.calc_Mean_logLik(i_cu);
          logL1_0 = data.calc_Mean_logLik_0(i_cu);
          logL1_2 = data.calc_Mean_logLik_2(i_cu);

          logL0_1 = data.calc_Mean_logLik_2(i_cu - 1);
          logL0_0 = data.calc_Mean_logLik(i_cu - 1);
          logL0 = logL0_1 - logL0_0;
          logL1 = logL1_1 - logL1_0;
          beta0 = 0;
          beta1 = data.get_Beta(i_cu)();
          plog_Evidence = (beta1 - beta0) * (logL0 + logL1) / 2;
          log_Evidence = log_Evidence + plog_Evidence;
          log_Evidence_no_0 = log_Evidence_no_0 + plog_Evidence;
        }
        s.f << iter << s.sep << i_cu << s.sep << i_frac() << s.sep << nsamples
            << s.sep << beta1 << s.sep << beta0 << s.sep << logPrior << s.sep
            << logL1 << s.sep << logL0 << s.sep << plog_Evidence << s.sep
            << log_Evidence << s.sep << log_Evidence_no_0 << s.sep << logL1_1
            << s.sep << logL1_0 << s.sep << logL1_2 << s.sep << logL0_1 << s.sep
            << logL0_0<<s.sep
            << prior_sta().count() << s.sep << prior_sta().rate() << s.sep
            << lik_sta_0.count() << s.sep << lik_sta_0.rate() << s.sep
            << lik_sta_1.count() << s.sep << lik_sta_1.rate() << s.sep
            << lik_sta_2.count() << s.sep << lik_sta_2.rate() << s.sep
            << emcee_sta().count() << s.sep << emcee_sta().rate() << s.sep
            << th_sta().count() << s.sep << th_sta().rate() << s.sep
            << cuevi_sta().count() << s.sep << cuevi_sta().rate()
            << "\n";
      }
    }
  }

  friend void report_title(save_Evidence &s, Cuevi_mcmc const &, ...) {

    s.f << "iter" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep << "nsamples"
        << s.sep << "beta1" << s.sep << "beta0" << s.sep << "logPrior" << s.sep
        << "logL1" << s.sep << "logL0" << s.sep << "plog_Evidence" << s.sep
        << "log_Evidence" << s.sep << "log_Evidence_no_0" << s.sep << "logL1_1"
        << s.sep << "logL1_0" << s.sep << "logL1_2" << s.sep << "logL0_1"
          << s.sep << "logL0_0"<<s.sep
          << "prior_sta_count" << s.sep << "prior_sta_rate" << s.sep
          << "lik_sta_0_count" << s.sep << "lik_sta_0_rate" << s.sep
          << "lik_sta_1_count" << s.sep << "lik_sta_1_rate" << s.sep
          << "lik_sta_2_count" << s.sep << "lik_sta_2_rate" << s.sep
          << "emcee_sta_count" << s.sep << "emcee_sta_rate" << s.sep
          << "th_sta_count" << s.sep << "th_sta_rate" << s.sep
          << "cuevi_sta_count" << s.sep << "cuevi_sta_rate"

        << "\n";
  }
};

template <class ParameterType>
auto extract_parameters_last(std::string const &filename, std::size_t iter) {
  save_Parameter<ParameterType> s(filename, 1);
  std::string last_line;
  std::size_t npar = 0ul;
  std::size_t nwalkers = 0ul;
  std::size_t n_cuevi = 0ul;
  bool cont = true;
  Maybe_error<Cuevi_mcmc<ParameterType>> out = error_message{};
  while (cont) {
    auto outnew = Cuevi_mcmc<ParameterType>::template extract_iter<false>(
        s, last_line, iter, npar, nwalkers, n_cuevi);
    cont = outnew.valid();
    if (outnew)
      out = outnew;
  }
  return out;
}

template <class... saving, class... Ts, class Parameters>
void report_title(save_mcmc<Parameters, saving...> &f,
                  Cuevi_mcmc<Parameters> const &data, const Ts &...ts) {
  (report_title(static_cast<saving &>(f), data, ts...), ..., 1);
}

class step_stretch_cuevi_mcmc_per_walker {
public:
  friend std::string ToString(step_stretch_cuevi_mcmc_per_walker) {
    return "step_stretch_cuevi_mcmc_per_walker";
  }

  template <class FunctionTable, class Prior, class Likelihood, class Variables,
            class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  auto operator()(FunctionTable &f, Cuevi_mcmc<Parameters> &current,
                  ensemble<mt_64i> &mts,
                  std::vector<std::uniform_real_distribution<double>> &rdist,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, std::size_t n_par,
                  std::size_t i, Walker_Index iw, Walker_Index jw,
                  Cuevi_Index i_cu) const {
    auto &wa_i = current.get_Walker_Value(iw, i_cu);
    auto const &wa_j = current.get_Walker_Value(jw, i_cu);

    auto &param_i = get<Parameter>(wa_i())();
    auto &param_j = get<Parameter>(wa_j())();
    auto &t = current.get_Cuevi_Temperatures();

    auto [ca_par, z] = stretch_move(mts[i], rdist[i], param_i, param_j);
    auto wa_sta = current.get_Walker_Statistics(iw, i_cu);
    auto Maybe_Walker = current.calc_Walker_value(
       f, std::move(ca_par), prior, lik, y, x, t,
        i_cu, wa_sta);
    if (Maybe_Walker) {
      auto ca_wa = std::move(Maybe_Walker.value());
      auto &cu_wa = current.get_Walker_Value(iw, i_cu);
      auto dthLogL = thermo_step(ca_wa, cu_wa, current.get_Beta(i_cu),
                                 current.get_Fraction(i_cu));
      if (dthLogL) {
        auto pJump =
            std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL.value()));
        auto r = rdist[i](mts[i]);
        if (pJump < r) {
          fails(get<emcee_Step_statistics>(wa_sta.first())());
          fails(get<emcee_Step_statistics>(wa_sta.second())());
        } else {
          succeeds(get<emcee_Step_statistics>(wa_sta.first())());
          succeeds(get<emcee_Step_statistics>(wa_sta.second())());
          cu_wa = ca_wa;
        }
      }
    }
  }
};

class step_stretch_cuevi_mcmc {
public:
  friend std::string ToString(step_stretch_cuevi_mcmc) {
    return "step_stretch_cuevi_mcmc";
  }

  template <class FunctionTable, class Prior, class Likelihood, class Variables,
            class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &f, Cuevi_mcmc<Parameters> &current,
                  ensemble<mt_64i> &mts, Prior const &prior,
                  Likelihood const &lik, const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  double alpha_stretch = 2) const {
    auto cuevi_temp = current.get_Cuevi_Temperatures();
    auto n_walkers = current.get_Walkers_number();

    auto &mt = mts[0];
    auto n_par = current.get_Parameters_number();

    WalkerIndexes shuffled_walker(n_walkers);
    std::iota(shuffled_walker.begin(), shuffled_walker.end(), 0);
    std::shuffle(shuffled_walker.begin(), shuffled_walker.end(), mt);

    std::uniform_int_distribution<std::size_t> uniform_walker(0, n_walkers / 2 -
                                                                     1);
    std::vector<std::uniform_int_distribution<std::size_t>> udist(
        n_walkers, uniform_walker);

    std::uniform_real_distribution<double> uniform_stretch_zdist(
        1.0 / alpha_stretch, alpha_stretch);
    std::vector<std::uniform_real_distribution<double>> zdist(
        n_walkers, uniform_stretch_zdist);

    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                              uniform_real);
    
    auto ff=f.fork( n_walkers / 2);
    for (bool half : {false, true}) {
      if (half)
        std::shuffle(shuffled_walker.begin(),
                     shuffled_walker.begin() + n_walkers / 2, mt);

#pragma omp parallel for
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto ii = half ? i + n_walkers / 2 : i;
        auto jj = half ? i : i + n_walkers / 2;

        auto iw = Walker_Index(shuffled_walker[ii]);
        auto jw = Walker_Index(shuffled_walker[jj]);

        for (std::size_t i_cu = 0; i_cu < size(cuevi_temp()); ++i_cu) {
          ff[i].f(step_stretch_cuevi_mcmc_per_walker{}, current, mts, rdist,
                 prior, lik, y, x, n_par, i, iw, jw, Cuevi_Index(i_cu));
        }
      }
    }
    f+=ff;
  }
};

class thermo_cuevi_jump_mcmc_per_walker {
public:
  friend std::string ToString(thermo_cuevi_jump_mcmc_per_walker) {
    return "thermo_cuevi_jump_mcmc_per_walker";
  }

  template <class FunctionTable, class Likelihood, class Variables,
            class DataType, class Parameters>

  void operator()(FunctionTable &f, Cuevi_mcmc<Parameters> &current,
                  Likelihood const &lik, const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, double r, Walker_Index i_w,
                  Walker_Index j_w, Cuevi_Index icu, Cuevi_Index jcu,
                  Random_jumps randomize) const {
    auto &t = current.get_Cuevi_Temperatures();
    auto &w_0 = current.get_Walker_Value(i_w, icu);
    auto &w_1 = current.get_Walker_Value(j_w, jcu);
    auto wa_sta_0 = current.get_Walker_Statistics(i_w, icu);
    auto wa_sta_1 = current.get_Walker_Statistics(j_w, jcu);
    auto Maybe_logA = current.thermo_jump_logProb(f, lik, y, x, t, icu, w_0,
                                                  jcu, w_1, wa_sta_0, wa_sta_1);
    if (Maybe_logA.valid()) {
      auto logA = std::move(Maybe_logA.value());
      auto pJump = std::min(1.0, std::exp(logA));
      if (pJump > r) {
        auto &W0 = current.get_Walker(i_w, icu);
        auto &W1 = current.get_Walker(j_w, jcu);
        if (randomize()) {
            succeeds(get<Cuevi_Jump_statistics>(wa_sta_0.first())());
            succeeds(get<Cuevi_Jump_statistics>(wa_sta_0.second())());
            succeeds(get<Cuevi_Jump_statistics>(wa_sta_1.first())());
            succeeds(get<Cuevi_Jump_statistics>(wa_sta_1.second())());
        } else {
            succeeds(get<Thermo_Jump_statistics>(wa_sta_0.first())());
            succeeds(get<Thermo_Jump_statistics>(wa_sta_0.second())());
            succeeds(get<Thermo_Jump_statistics>(wa_sta_1.first())());
            succeeds(get<Thermo_Jump_statistics>(wa_sta_1.second())());
        }

        std::swap(W0, W1);
      } else {
          if (randomize()) {
              fails(get<Cuevi_Jump_statistics>(wa_sta_0.first())());
              fails(get<Cuevi_Jump_statistics>(wa_sta_1.second())());
              fails(get<Cuevi_Jump_statistics>(wa_sta_0.first())());
              fails(get<Cuevi_Jump_statistics>(wa_sta_1.second())());
          } else {
              fails(get<Thermo_Jump_statistics>(wa_sta_0.first())());
              fails(get<Thermo_Jump_statistics>(wa_sta_1.second())());
              fails(get<Thermo_Jump_statistics>(wa_sta_0.first())());
              fails(get<Thermo_Jump_statistics>(wa_sta_1.second())());
          }
      }
    }
  }
};

class thermo_cuevi_jump_mcmc {
public:
  friend std::string ToString(thermo_cuevi_jump_mcmc) {
    return "thermo_cuevi_jump_mcmc";
  }

  template <class FunctionTable, class Prior, class Likelihood, class Variables,
            class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &f, std::size_t iter,
                  Cuevi_mcmc<Parameters> &current, ensemble<mt_64i> &mts,
                  Prior const &, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  Thermo_Jumps_every thermo_jumps_every,
                  Random_jumps randomize) const {
    if (iter % (thermo_jumps_every()) == 0) {
      auto t = current.get_Cuevi_Temperatures();
      auto n_walkers = current.get_Walkers_number();

      std::uniform_real_distribution<double> uniform_real(0, 1);
      auto n_par = current.get_Parameters_number();
      WalkerIndexes shuffled_walker(n_walkers);
      std::iota(shuffled_walker.begin(), shuffled_walker.end(), 0);
      std::shuffle(shuffled_walker.begin(), shuffled_walker.end(), mts[0]);

      std::vector<std::size_t> shuffled_cuevi(t().size());
      std::iota(shuffled_cuevi.begin(), shuffled_cuevi.end(), 0);
      if (randomize())
        std::shuffle(shuffled_cuevi.begin(), shuffled_cuevi.end(), mts[0]);

      std::uniform_int_distribution<std::size_t> uniform_walker(
          0, n_walkers / 2 - 1);
      std::vector<std::uniform_int_distribution<std::size_t>> udist(
          n_walkers, uniform_walker);

      std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                                uniform_real);

      Maybe_error<bool> success = true;
      
      auto ff=f.fork(n_walkers / 2);

#pragma omp parallel for
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {

        auto i_w = Walker_Index(shuffled_walker[i]);
        auto j_w = Walker_Index(shuffled_walker[i + n_walkers / 2]);

        for (std::size_t i_cu = 0; i_cu + 1 < size(t()); ++i_cu) {
          auto icu = Cuevi_Index(shuffled_cuevi[i_cu]);
          auto jcu = Cuevi_Index(shuffled_cuevi[i_cu + 1]);
          double r = rdist[i](mts[i]);
          thermo_cuevi_jump_mcmc_per_walker{}(ff[i], current,
                                              lik, y, x, r, i_w, j_w, icu, jcu,
                                              randomize);
        }
      }
      f+=ff;
    }
  }
};

struct Fractioner {};
struct Reporter {};
struct Finalizer {};

template <class myFractioner, class t_Reporter, class t_Finalizer>
class Cuevi_Algorithm
    : public var::Var<
          Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>,
          var::Vector_Space<
              Num_Walkers_Per_Ensemble, var::Constant<Fractioner, myFractioner>,
              var::Constant<Reporter, t_Reporter>,
              var::Constant<Finalizer, t_Finalizer>, Fractions_Param,
              Th_Beta_Param, Number_trials_until_give_up, Thermo_Jumps_every,
              Random_jumps, Saving_intervals>> {
  using base_type = var::Var<
      Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>,
      var::Vector_Space<Num_Walkers_Per_Ensemble,
                        var::Constant<Fractioner, myFractioner>,
                        var::Constant<Reporter, t_Reporter>,
                        var::Constant<Finalizer, t_Finalizer>, Fractions_Param,
                        Th_Beta_Param, Number_trials_until_give_up,
                        Thermo_Jumps_every, Random_jumps, Saving_intervals>>;

public:
  Cuevi_Algorithm(myFractioner &&frac, t_Reporter &&rep, t_Finalizer &&f,
                  Num_Walkers_Per_Ensemble n, Fractions_Param frac_param,
                  Th_Beta_Param beta,
                  Number_trials_until_give_up max_iter_for_sampling,
                  Thermo_Jumps_every thermo_jumps_every, Random_jumps randomize,
                  Saving_intervals saving_intervals)
      : base_type(var::Vector_Space(
            std::move(n),
            var::Constant<Fractioner, myFractioner>(std::move(frac)),
            var::Constant<Reporter, t_Reporter>(std::move(rep)),
            var::Constant<Finalizer, t_Finalizer>(std::move(f)),
            std::move(frac_param), std::move(beta),
            std::move(max_iter_for_sampling), std::move(thermo_jumps_every),
            std::move(randomize), std::move(saving_intervals))) {}
};


template <class t_Reporter, class t_Finalizer>
class Cuevi_Algorithm_no_Fractioner
    : public var::Var<
          Cuevi_Algorithm_no_Fractioner< t_Reporter, t_Finalizer>,
          var::Vector_Space<
              Num_Walkers_Per_Ensemble, 
              var::Constant<Reporter, t_Reporter>,
              var::Constant<Finalizer, t_Finalizer>, 
              Th_Beta_Param, Number_trials_until_give_up, Thermo_Jumps_every,
              Random_jumps, Saving_intervals>> {
    using base_type = var::Var<
        Cuevi_Algorithm_no_Fractioner< t_Reporter, t_Finalizer>,
        var::Vector_Space<Num_Walkers_Per_Ensemble,
                          var::Constant<Reporter, t_Reporter>,
                          var::Constant<Finalizer, t_Finalizer>, 
                          Th_Beta_Param, Number_trials_until_give_up,
                          Thermo_Jumps_every, Random_jumps, Saving_intervals>>;
    
public:
    Cuevi_Algorithm_no_Fractioner(t_Reporter &&rep, t_Finalizer &&f,
                    Num_Walkers_Per_Ensemble n, 
                    Th_Beta_Param beta,
                    Number_trials_until_give_up max_iter_for_sampling,
                    Thermo_Jumps_every thermo_jumps_every, Random_jumps randomize,
                    Saving_intervals saving_intervals)
        : base_type(var::Vector_Space(
            std::move(n),
            var::Constant<Reporter, t_Reporter>(std::move(rep)),
            var::Constant<Finalizer, t_Finalizer>(std::move(f)),
             std::move(beta),
            std::move(max_iter_for_sampling), std::move(thermo_jumps_every),
            std::move(randomize), std::move(saving_intervals))) {}
};


template <class FunctionTable, class Parameters, class... saving, class... T>
void report_all(FunctionTable &f, std::size_t iter,
                save_mcmc<Parameters, saving...> &s,
                Cuevi_mcmc<Parameters> &data, T const &...ts) {
  (report(f, iter, static_cast<saving &>(s), data, ts...), ..., 1);
}

template <class Parameter, class... saving>
void report_model_all(save_mcmc<Parameter, saving...> &) {}

template <class Parameter, class... saving, class T, class... Ts>
void report_model_all(save_mcmc<Parameter, saving...> &s, T const &t,
                      Ts const &...ts) {
  (report_model(static_cast<saving &>(s), t), ..., report_model_all(s, ts...));
}

template <class myFractioner, class t_Reporter, class t_Finalizer>
Cuevi_Algorithm(myFractioner &&frac, t_Reporter &&rep, t_Finalizer &&f,
                Num_Walkers_Per_Ensemble n, Fractions_Param frac_param,
                Random_jumps random_jumps, Th_Beta_Param beta)
    -> Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>;

template <class FunctionTable, class myFractioner, class t_Reporter,
          class t_Finalizer, class Prior, class Likelihood, class DataType,
          class Variables>
auto evidence_old(FunctionTable &ff,
                  Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer> &&cue,
                  Prior &&prior, Likelihood const &lik, const DataType &y,
                  const Variables &x, const Init_seed init_seed) {
  using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i &>>;

  using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

  using Return_Type =
      Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;

  auto f = ff.fork(var::I_thread(0));
  auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());

  auto mt = init_mt(init_seed());
  auto mts = init_mts(mt, n_walkers() / 2);
  auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
  auto n_points_per_decade_fraction =
      get<Points_per_decade>(get<Fractions_Param>(cue())());

  auto [ys, xs] = get<Fractioner>(cue())()(
      y, x, mt, size(prior) * min_fraction(), n_points_per_decade_fraction());

  auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
      std::forward<FunctionTable>(ff), mts, std::forward<Prior>(prior), lik, ys,
      xs, get<Th_Beta_Param>(cue()), n_walkers,
      get<Number_trials_until_give_up>(cue()));

  if (!Maybe_current)
    return Return_Type(Maybe_current.error());
  else {

    auto current = std::move(Maybe_current.value());
    auto &rep = get<Reporter>(cue())();

    auto a = get<Finalizer>(cue())();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    // report_model(rep, prior, lik, ys, xs);
    report_title(ff, "Iter");
    report_title(rep, current, prior, lik, ys, xs);

    while (!mcmc_run.second) {
    f.f(
      step_stretch_cuevi_mcmc{}, current, mts, prior, lik, ys, xs);
      report_point(ff, iter);

      ++iter;
      f.f(
      thermo_cuevi_jump_mcmc{}, iter, current, mt, mts, prior, lik, ys, xs,
                               get<Thermo_Jumps_every>(cue())(), false);
       f.f(
      thermo_cuevi_jump_mcmc{}, iter + get<Thermo_Jumps_every>(cue())() % 2,
                               current, mt, mts, prior, lik, ys, xs,
                               get<Thermo_Jumps_every>(cue())(), true);

      report_all(f, iter, rep, current, prior, lik, ys, xs);
      mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
    }
    return Return_Type(std::pair(mcmc_run.first, current));
  }
}

template <class FunctionTable, class ParameterType, class t_Reporter,
          class mcmc_type, class Prior, class Likelihood, class DataType,
          class Variables>
auto evidence_loop(FunctionTable &f, std::pair<mcmc_type, bool> &&mcmc_run,
                   t_Reporter &rep, std::size_t iter,
                   Cuevi_mcmc<ParameterType> &current,

                   std::vector<mt_64i> &mts, Random_jumps randomjumps,
                   Thermo_Jumps_every v_thermo_jump_every, Prior &&prior,
                   Likelihood const &lik, const by_fraction<DataType> &ys,
                   const by_fraction<Variables> &xs) {
  using Return_Type =
      Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<ParameterType>>>;
  report_title(rep, current, mcmc_run.first, prior, lik, ys, xs);
 
  while (!mcmc_run.second) {
    f.f(
    step_stretch_cuevi_mcmc{}, current, mts, prior, lik, ys, xs);
    report_point(f, iter);

    ++iter;
    f.f(
    thermo_cuevi_jump_mcmc{}, iter, current, mts, prior, lik, ys, xs,
                             v_thermo_jump_every, Random_jumps(false));
    if (randomjumps())
    f.f(
      thermo_cuevi_jump_mcmc{}, iter + v_thermo_jump_every() % 2, current,
                               mts, prior, lik, ys, xs, v_thermo_jump_every,
                               Random_jumps(true));

    report_all(f, iter, rep, current, prior, lik, ys, xs, mts, mcmc_run.first);
    mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
  }
  return Return_Type(std::pair(mcmc_run.first, current));
}




template <class FunctionTable, class myFractioner, class t_Reporter,
          class t_Finalizer, class Prior, class Likelihood, class DataType,
          class Variables>
auto evidence(FunctionTable &f,
              Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer> &&cue,
              Prior &&prior, Likelihood const &lik, const DataType &y,
              const Variables &x, const Init_seed init_seed) {
  using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i &>>;

  using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

  using Return_Type =
      Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;


  auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
  auto mt = init_mt(init_seed());
  auto mts = init_mts(mt, n_walkers() / 2);
  auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
  auto n_points_per_decade_fraction =
      get<Points_per_decade>(get<Fractions_Param>(cue())());

  auto [ys, xs] = get<Fractioner>(cue())()(
      y, x, mt, size(prior) * min_fraction(), n_points_per_decade_fraction());

  auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
      f, mts, std::forward<Prior>(prior), lik, ys,
      xs, get<Th_Beta_Param>(cue()), n_walkers,
      get<Number_trials_until_give_up>(cue()));

  if (!Maybe_current)
    return Return_Type(Maybe_current.error());
  else {

    auto current = std::move(Maybe_current.value());
    auto &rep = get<Reporter>(cue())();

    auto a = get<Finalizer>(cue())();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    report_title(f, "Iter");
    auto v_random_jumps = get<Random_jumps>(cue());
    auto v_thermo_jump_every = get<Thermo_Jumps_every>(cue());
    report_model_all(rep, v_thermo_jump_every, prior, lik, ys, xs);

    // report_init(rep,v_thermo_jump_every,prior,lik,ys,xs);

    return evidence_loop(f, std::move(mcmc_run),
                         rep, iter, current, mts, v_random_jumps,
                         v_thermo_jump_every, std::forward<Prior>(prior), lik,
                         ys, xs);
  }
}



template <class FunctionTable,  class t_Reporter,
         class t_Finalizer, class Prior, class Likelihood, class DataType,
         class Variables>
auto evidence_fraction(FunctionTable &f,
              Cuevi_Algorithm_no_Fractioner< t_Reporter, t_Finalizer> &&cue,
                       Prior &&prior, Likelihood const &lik, const std::vector<DataType> &ys,
                       const std::vector<Variables> &xs, const Init_seed init_seed) {
    using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i &>>;
    
    using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;
    
    using Return_Type =
        Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;
    
    
    auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
    auto mt = init_mt(init_seed());
    auto mts = init_mts(mt, n_walkers() / 2);
 //   auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
    // auto n_points_per_decade_fraction =
    //     get<Points_per_decade>(get<Fractions_Param>(cue())());
    
    // auto [ys, xs] = get<Fractioner>(cue())()(
    //     y, x, mt, size(prior) * min_fraction(), n_points_per_decade_fraction());
    
    auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
        f, mts, std::forward<Prior>(prior), lik, ys,
        xs, get<Th_Beta_Param>(cue()), n_walkers,
        get<Number_trials_until_give_up>(cue()));
    
    if (!Maybe_current)
        return Return_Type(Maybe_current.error());
    else {
        
        auto current = std::move(Maybe_current.value());
        auto &rep = get<Reporter>(cue())();
        
        auto a = get<Finalizer>(cue())();
        auto mcmc_run = checks_convergence(std::move(a), current);
        
        std::size_t iter = 0;
        report_title(f, "Iter");
        auto v_random_jumps = get<Random_jumps>(cue());
        auto v_thermo_jump_every = get<Thermo_Jumps_every>(cue());
        report_model_all(rep, v_thermo_jump_every, prior, lik, ys, xs);
        
        // report_init(rep,v_thermo_jump_every,prior,lik,ys,xs);
        
        return evidence_loop(f, std::move(mcmc_run),
                             rep, iter, current, mts, v_random_jumps,
                             v_thermo_jump_every, std::forward<Prior>(prior), lik,
                             ys, xs);
    }
}



template <class FunctionTable, class myFractioner, class t_Reporter,
          class t_Finalizer, class Prior, class Likelihood, class DataType,
          class Variables>
auto continue_evidence(
    FunctionTable &f,
    Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer> &&cue, Prior &&prior,
    Likelihood const &lik, const DataType &y, const Variables &x,
    const Init_seed init_seed) {
  using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i &>>;

  using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

  using Return_Type =
      Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;

 
  auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
  auto mt = init_mt(init_seed());
  auto mts = init_mts(mt, n_walkers() / 2);
  auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
  auto n_points_per_decade_fraction =
      get<Points_per_decade>(get<Fractions_Param>(cue())());

  auto [ys, xs] = get<Fractioner>(cue())()(
      y, x, mt, size(prior) * min_fraction(), n_points_per_decade_fraction());

  auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
      f, mts, std::forward<Prior>(prior), lik, ys,
      xs, get<Th_Beta_Param>(cue()), n_walkers,
      get<Number_trials_until_give_up>(cue()));

  if (!Maybe_current)
    return Return_Type(Maybe_current.error());
  else {

    auto current = std::move(Maybe_current.value());
    auto &rep = get<Reporter>(cue())();
    report_title(rep, current, prior, lik, ys, xs);

    auto a = get<Finalizer>(cue())();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    // report_model(rep, prior, lik, ys, xs);
    report_title(f, "Iter");
    std::size_t v_thermo_jump_every = get<Thermo_Jumps_every>(cue());

    return evidence_loop(f, v_thermo_jump_every,
                         mcmc_run, rep, iter, current, mts,
                         std::forward<Prior>(prior), lik, ys, xs);
  }
}

} // namespace cuevi

namespace deprecated{

template <class Parameters> struct mcmc2 : public mcmc<Parameters> {
  double logPa;
};

template <class T> using by_fraction = std::vector<T>;

template <class Parameters> struct cuevi_mcmc {
  by_fraction<std::size_t> nsamples;
  by_fraction<by_beta<double>> beta;
  ensemble<by_fraction<by_beta<mcmc2<Parameters>>>> walkers;
  ensemble<by_fraction<by_beta<std::size_t>>> i_walkers;
  by_fraction<by_beta<int>> is_active;

  auto current_number_of_temperatures() const {
    std::size_t count = 0;
    if (is_active[0][0] == 1)
      ++count;
    for (auto i_frac = 0ul; i_frac < is_active.size(); ++i_frac) {
      for (std::size_t ib = 1; ib < is_active[i_frac].size(); ++ib)
        if (is_active[i_frac][ib] == 1)
          ++count;
    }
    return count;
  }

  friend void report_title(save_likelihood<Parameters> &s,
                           cuevi_mcmc<Parameters> const &, ...) {

    s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
        << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "logPa" << s.sep << "logP" << s.sep
        << "logLik"
        << "\n";
  }
  friend void report_title(save_Evidence &s, cuevi_mcmc const &, ...) {

    s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
        << "nsamples" << s.sep << "beta" << s.sep << "meanPrior" << s.sep
        << "meanLik" << s.sep << "varLik" << s.sep
        << "fraction_Evidence_by_mean" << s.sep << "fraction_Evidence_by_var"
        << s.sep << "Evidence_by_mean" << s.sep << "Evidence_by_var"
        << "\n";
  }

  friend void report_title(save_Parameter<Parameters> &s, cuevi_mcmc const &,
                           ...) {

    s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
        << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "i_par" << s.sep << "par_value"
        << "\n";
  }

  template <class... saving, class... Ts>
  friend void report_title(save_mcmc<Parameters, saving...> &f,
                           cuevi_mcmc const &data, const Ts &...ts) {
    (report_title(static_cast<saving &>(f), data, ts...), ..., 1);
  }
};
template <class... T>
std::ostream &print_i_walkers_title(std::ostream &os, const T &...ts) {
  ((os << ts << ",") << ... << "i_frac");
  os << ","
     << "i_beta"
     << ","
     << "i_walker"
     << ","
     << "id_walker"
     << "\n";
  return os;
}

template <class Parameters, class... T>
std::ostream &print_i_walkers(std::ostream &os,
                              cuevi_mcmc<Parameters> const &data,
                              const T &...ts) {
  os << "i_walker"
     << ","
     << "i_frac"
     << ","
     << "i_beta"
     << ","
     << "id_walker"
     << "\n";
  for (std::size_t i_w = 0; i_w < data.i_walkers.size(); ++i_w) {
    for (std::size_t i_frac = 0; i_frac < data.beta.size(); ++i_frac) {
      for (std::size_t i_beta = 0; i_beta < data.beta[i_frac].size();
           ++i_beta) {
        ((os << ts << ",") << ... << "");
        os << i_w << "," << i_frac << "," << i_beta << ","
           << data.i_walkers[i_w][i_frac][i_beta] << "\n";
      }
    }
  }
  return os;
}

template <class Parameters> auto mean_logL(cuevi_mcmc<Parameters> const &mcmc) {
  auto out = by_fraction<by_beta<double>>{};
  out.reserve(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {
    out.emplace_back(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t iwalker = 0; iwalker < n_walkers; ++iwalker)
      for (std::size_t ibeta = 0; ibeta < size(mcmc.beta[i_frac]); ++ibeta)
        out[i_frac][ibeta] +=
            mcmc.walkers[iwalker][i_frac][ibeta].logL / n_walkers;
  }
  return out;
}
template <class Parameters> auto mean_logP(cuevi_mcmc<Parameters> const &mcmc) {
  auto out = by_fraction<by_beta<double>>{};
  out.reserve(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {
    out.emplace_back(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t iwalker = 0; iwalker < n_walkers; ++iwalker)
      for (std::size_t ibeta = 0; ibeta < size(mcmc.beta[i_frac]); ++ibeta)
        out[i_frac][ibeta] +=
            mcmc.walkers[iwalker][i_frac][ibeta].logP / n_walkers;
  }
  return out;
}

template <class Parameters>
auto var_logL(cuevi_mcmc<Parameters> const &mcmc,
              by_fraction<by_beta<double>> const &mean) {
  auto out = by_fraction<by_beta<double>>{};
  out.reserve(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {
    out.emplace_back(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t iwalker = 0; iwalker < n_walkers; ++iwalker)
      for (std::size_t ibeta = 0; ibeta < size(mcmc.beta[i_frac]); ++ibeta)
        out[i_frac][ibeta] +=
            std::pow(mcmc.walkers[iwalker][i_frac][ibeta].logL -
                         mean[i_frac][ibeta],
                     2) /
            n_walkers;
  }
  return out;
}

template <class Parameters>
auto mean_logL(by_iteration<cuevi_mcmc<Parameters>> const &series) {
  auto const &mcmc = series[0];
  auto out = by_fraction<by_beta<double>>(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  auto n_iters = size(series);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {

    out[i_frac] = by_beta<double>(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t i = 0; i < size(series); ++i)
      for (std::size_t iwalker = 0; iwalker < size(mcmc.walkers); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < size(mcmc.beta[i_frac]); ++ibeta)
          out[i_frac][ibeta] += series[i].walkers[iwalker][i_frac][ibeta].logL /
                                n_iters / n_walkers;
  }
  return out;
}

template <class Parameters>
auto mean_logP(by_iteration<cuevi_mcmc<Parameters>> const &series) {
  auto const &mcmc = series[0];
  auto out = by_fraction<by_beta<double>>(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  auto n_iters = size(series);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {
    out[i_frac] = by_beta<double>(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t i = 0; i < size(series); ++i)
      for (std::size_t iwalker = 0; iwalker < size(mcmc.walkers); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < size(mcmc.beta[i_frac]); ++ibeta)
          out[i_frac][ibeta] += series[i].walkers[iwalker][i_frac][ibeta].logP /
                                n_iters / n_walkers;
  }
  return out;
}

template <class Parameters>
auto var_logL(by_iteration<cuevi_mcmc<Parameters>> const &series,
              by_fraction<by_beta<double>> const &mean) {
  auto const &mcmc = series[0];
  auto out = by_fraction<by_beta<double>>(size(mcmc.beta));
  auto n_walkers = size(mcmc.walkers);
  auto n_iters = size(series);
  for (std::size_t i_frac = 0; i_frac < size(mcmc.beta); ++i_frac) {
    out[i_frac] = by_beta<double>(size(mcmc.beta[i_frac]), 0.0);
    for (std::size_t i = 0; i < size(series); ++i)
      for (std::size_t iwalker = 0; iwalker < size(series[0].walkers);
           ++iwalker)
        for (std::size_t ibeta = 0; ibeta < size(series[0].beta[i_frac]);
             ++ibeta)
          out[i_frac][ibeta] +=
              std::pow(series[i].walkers[iwalker][i_frac][ibeta].logL -
                           mean[i_frac][ibeta],
                       2) /
              n_iters / n_walkers;
  }
  return out;
}

inline by_fraction<double>
calculate_Evidence(by_fraction<by_beta<double>> const &beta,
                   by_fraction<by_beta<double>> const &meanLik) {
  auto nfraction = beta.size();
  auto out = by_fraction<double>(nfraction, 0.0);
  for (std::size_t i_frac = 0; i_frac < nfraction; ++i_frac) {
      out[i_frac] = ::calculate_Evidence(beta[i_frac], meanLik[i_frac]);
  }
  return out;
}

inline by_fraction<double>
calculate_Evidence(by_fraction<by_beta<double>> const &beta,
                   by_fraction<by_beta<double>> const &meanLik,
                   by_fraction<by_beta<double>> const &varLik) {
  auto nfraction = beta.size();
  auto out = by_fraction<double>(nfraction, 0.0);
  for (std::size_t i_frac = 0; i_frac < nfraction; ++i_frac) {
    out[i_frac] =
          ::calculate_Evidence(beta[i_frac], meanLik[i_frac], varLik[i_frac]);
  }
  return out;
}

template <class Cova>
  requires Covariance<double, Cova>
auto cuevi_posterior(
    conjugate,
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const linear_model &, const Matrix<double> &y0, const Matrix<double> &X0,
    const Matrix<double> &y1, const Matrix<double> &X1) {
  auto a_0 = prior.alpha();
  ;
  auto prior_eps_df = 2.0 * a_0;
  auto b_0 = prior.beta();
  auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

  auto L_0 = prior.Gamma();
  auto SSx = XTX(X1) - XTX(X0);
  auto n = y1.nrows() - y0.nrows();
  auto beta_0 = prior.mean();
  SymPosDefMatrix<double> L_n =
      SymPosDefMatrix<double>::I_sware_it_is_possitive(L_0 + SSx);

  auto beta_n =
      tr(inv(L_n) * (tr(X1) * y1 - tr(X0) * y0 + (L_0 * tr(prior.mean()))));

  auto yfit1 = X1 * tr(beta_n);
  auto ydiff1 = y1 - yfit1;
  auto yfit0 = X0 * tr(beta_n);
  auto ydiff0 = y0 - yfit0;
  auto SS = xtx(ydiff1.value()) - xtx(ydiff0.value());

  auto a_n = a_0 + n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

  auto posterior_Normal = make_multivariate_normal_distribution_from_precision(
      std::move(beta_n), std::move(L_n));

  return multivariate_gamma_normal_distribution<double,
                                                SymPosDefMatrix<double>>(
      log_inverse_gamma_distribution(a_n, b_n.value()),
      std::move(posterior_Normal.value()));
}

template <class Cova>
  requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_posterior(
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const Matrix<double> &y0, const Matrix<double> &X0,
    const Matrix<double> &y1, const Matrix<double> &X1) {
  auto a_0 = prior.alpha();
  ;
  auto prior_eps_df = 2.0 * a_0;
  auto b_0 = prior.beta();
  auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

  auto L_0 = prior.Gamma();
  auto SSx = XTX(X1) - XTX(X0);
  auto n = y1.nrows() - y0.nrows();
  auto beta_0 = prior.mean();
  SymPosDefMatrix<double> L_n =
      SymPosDefMatrix<double>::I_sware_it_is_possitive(L_0 + SSx);

  auto beta_n =
      tr(inv(L_n) * (tr(X1) * y1 - tr(X0) * y0 + (L_0 * tr(prior.mean()))));

  auto yfit1 = X1 * tr(beta_n);
  auto ydiff1 = y1 - yfit1;
  auto yfit0 = X0 * tr(beta_n);
  auto ydiff0 = y0 - yfit0;
  auto SS = xtx(ydiff1.value()) - xtx(ydiff0.value());

  auto a_n = a_0 + n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

  auto posterior_Normal = make_multivariate_normal_distribution_from_precision(
      std::move(beta_n), std::move(L_n));

  return multivariate_gamma_normal_distribution<double,
                                                SymPosDefMatrix<double>>(
      log_inverse_gamma_distribution(a_n, b_n.value()),
      std::move(posterior_Normal.value()));
}

template <class Cova>
  requires Covariance<double, Cova>
auto cuevi_evidence(
    conjugate,
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const linear_model &, const Matrix<double> &y0, const Matrix<double> &X0,
    const Matrix<double> &y1, const Matrix<double> &X1) {
  auto a_0 = prior.alpha();
  ;
  auto prior_eps_df = 2.0 * a_0;
  auto b_0 = prior.beta();
  auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

  auto L_0 = prior.Gamma();

  auto SSx = XTX(X1) + XTX(X0) * -1.0;
  auto n = y1.nrows() - y0.nrows();
  auto beta_0 = prior.mean();
  auto L_n = SSx + L_0;
  auto beta_n =
      tr(inv(L_n) * (tr(X1) * y1 - tr(X0) * y0 + (L_0 * tr(prior.mean()))));

  auto yfit1 = X1 * tr(beta_n);
  auto ydiff1 = y1 - yfit1;
  auto yfit0 = X0 * tr(beta_n);
  auto ydiff0 = y0 - yfit0;
  auto SS = xtx(ydiff1.value()) - xtx(ydiff0.value());
  auto a_n = a_0 + n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

  auto logE_n = -0.5 * n * std::log(2 * std::numbers::pi) +
                0.5 * (logdet(L_0) - logdet(L_n)) + a_0 * log(b_0) -
                a_n * log(b_n) + var::lgamma(a_n) - var::lgamma(a_0);
  return logE_n;
}

template <class Cova>
  requires Covariance<double, Cova>
auto cuevi_mean_logLik(
    conjugate,
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const linear_model &, const Matrix<double> &y0, const Matrix<double> &X0,
    const Matrix<double> &y1, const Matrix<double> &X1, double beta0) {
  auto a_0 = prior.alpha();
  auto b_0 = prior.beta();
  auto L_0 = prior.Gamma();
  auto SSx = XTX(X1) + XTX(X0) * -1.0;
  auto n = y1.nrows() - y0.nrows();
  auto beta_0 = prior.mean();
  auto L_n = L_0 + beta0 * SSx;
  auto beta_n = tr(inv(L_n) * (beta0 * (tr(X1) * y1 - tr(X0) * y0) +
                               (L_0 * tr(prior.mean()))));
  auto yfit1 = X1 * tr(beta_n);
  auto ydiff1 = y1 - yfit1;
  auto yfit0 = X0 * tr(beta_n);
  auto ydiff0 = y0 - yfit0;
  auto SS = beta0 * xtx(ydiff1.value()) - beta0 * xtx(ydiff0.value());

  auto a_n = a_0 + beta0 * n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
  double d_a_n = 1.0 * n / 2.0;
  auto d_b_n = 0.5 * xtx(ydiff1.value()) - 0.5 * xtx(ydiff0.value());
  auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                    0.5 * Trace(inv(L_n) * SSx) - a_n / b_n * d_b_n +
                    (digamma(a_n) - log(b_n)) * d_a_n;
  return mean_logLi;
}

template <class Cova>
  requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const Matrix<double> &y0, const Matrix<double> &X0,
    const Matrix<double> &y1, const Matrix<double> &X1, double beta0) {
  auto a_0 = prior.alpha();
  auto b_0 = prior.beta();
  auto L_0 = prior.Gamma();
  auto SSx = XTX(X1) + XTX(X0) * -1.0;
  auto n = y1.nrows() - y0.nrows();
  auto beta_0 = prior.mean();
  auto L_n = L_0 + beta0 * SSx;
  auto beta_n = tr(inv(L_n) * (beta0 * (tr(X1) * y1 - tr(X0) * y0) +
                               (L_0 * tr(prior.mean()))));
  auto yfit1 = X1 * tr(beta_n);
  auto ydiff1 = y1 - yfit1;
  auto yfit0 = X0 * tr(beta_n);
  auto ydiff0 = y0 - yfit0;
  auto SS = beta0 * xtx(ydiff1.value()) - beta0 * xtx(ydiff0.value());

  auto a_n = a_0 + beta0 * n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
  double d_a_n = 1.0 * n / 2.0;
  auto d_b_n = 0.5 * xtx(ydiff1.value()) - 0.5 * xtx(ydiff0.value());
  auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                    0.5 * Trace(inv(L_n) * SSx) - a_n / b_n * d_b_n +
                    (digamma(a_n) - log(b_n)) * d_a_n;
  return mean_logLi;
}

template <class FunctionTable, class Parameters>
void report(FunctionTable &&, std::size_t iter, save_likelihood<Parameters> &s,
            cuevi_mcmc<Parameters> const &data, ...) {
  if (iter % s.save_every == 0)
    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
      for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
        if (data.is_active[i_frac][i_beta] == 1)
          for (std::size_t i_walker = 0; i_walker < size(data.walkers);
               ++i_walker)

            s.f << size(data.beta) << s.sep << size(data.beta[i_frac]) << s.sep
                << iter << s.sep << data.nsamples[i_frac] << s.sep
                << data.beta[i_frac][i_beta] << s.sep << i_walker << s.sep
                << data.i_walkers[i_walker][i_frac][i_beta] << s.sep
                << data.walkers[i_walker][i_frac][i_beta].logPa << s.sep
                << data.walkers[i_walker][i_frac][i_beta].logP << s.sep
                << data.walkers[i_walker][i_frac][i_beta].logL << "\n";
}

template <class FunctionTable, class Parameters>
void report(FunctionTable &&, std::size_t iter, save_Parameter<Parameters> &s,
            cuevi_mcmc<Parameters> const &data, ...) {
  if (iter % s.save_every == 0)
    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
      for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
        if (data.is_active[i_frac][i_beta] == 1)
          for (std::size_t i_walker = 0; i_walker < size(data.walkers);
               ++i_walker)
            for (std::size_t i_par = 0;
                 i_par < size(data.walkers[i_walker][i_frac][i_beta].parameter);
                 ++i_par)

              s.f << size(data.beta) << s.sep << size(data.beta[i_frac])
                  << s.sep << iter << s.sep << data.nsamples[i_frac] << s.sep
                  << data.beta[i_frac][i_beta] << s.sep << i_walker << s.sep
                  << data.i_walkers[i_walker][i_frac][i_beta] << s.sep << i_par
                  << s.sep
                  << data.walkers[i_walker][i_frac][i_beta].parameter[i_par]
                  << "\n";
}

template <class Prior, class Likelihood, class Variables, class DataType>
concept has_conjugate = requires(Prior const &prior, Likelihood const &lik,
                                 const DataType &y, const Variables &x) {
  { posterior(conjugate{}, prior, lik, y[0], x[0]) };
};

template <class Prior, class Likelihood, class Variables, class DataType>
  requires has_conjugate<Prior, Likelihood, Variables, DataType>
void report_model(save_Evidence &s, Prior const &prior, Likelihood const &lik,
                  const DataType &y, const Variables &x,
                  by_fraction<by_beta<double>> const &beta0) {

  by_fraction<by_beta<Maybe_error<double>>> expected_meanLik(size(beta0));

  auto expected_Evidence =
      evidence(conjugate{}, prior, lik, y.back(), x.back());
  // bayesian_linear_regression_calculate_Evidence(prior,lik, y.back(),
  // x.back());

  by_fraction<Maybe_error<double>> partial_expected_evidence(size(beta0));
  by_fraction<Maybe_error<double>> expected_partial_evidence_by_logLik(
      size(beta0));

  partial_expected_evidence[0] = evidence(conjugate{}, prior, lik, y[0], x[0]);
  //   bayesian_linear_regression_calculate_Evidence(prior,lik, y[0], x[0]);

  expected_meanLik[0] = by_beta<Maybe_error<double>>(size(beta0[0]));

  for (std::size_t i_beta = 0; i_beta < size(beta0[0]); ++i_beta)
    expected_meanLik[0][i_beta] =
        mean_logLik(conjugate{}, prior, lik, y[0], x[0], beta0[0][i_beta]);
  auto model_i = posterior(conjugate{}, prior, lik, y[0], x[0]);
  auto model_f = posterior(conjugate{}, prior, lik, y.back(), x.back());

  for (std::size_t i_frac = 1; i_frac < size(beta0); ++i_frac) {
    partial_expected_evidence[i_frac] =
        cuevi_evidence(conjugate{}, model_i, linear_model{}, y[i_frac - 1],
                       x[i_frac - 1], y[i_frac], x[i_frac]);
    expected_meanLik[i_frac] =
        by_beta<Maybe_error<double>>(size(beta0[i_frac]));
    for (std::size_t i_beta = 0; i_beta < size(beta0[i_frac]); ++i_beta)
      expected_meanLik[i_frac][i_beta] = cuevi_mean_logLik(
          conjugate{}, model_i, linear_model{}, y[i_frac - 1], x[i_frac - 1],
          y[i_frac], x[i_frac], beta0[i_frac][i_beta]);
    model_i =
        cuevi_posterior(conjugate{}, model_i, linear_model{}, y[i_frac - 1],
                        x[i_frac - 1], y[i_frac], x[i_frac]);
  }

  Maybe_error<double> sum_partial_Evidence = 0.0;

  for (std::size_t i_frac = 0; i_frac < size(beta0); ++i_frac) {
      expected_partial_evidence_by_logLik[i_frac] = ::calculate_Evidence(
        beta0[i_frac], promote_Maybe_error(expected_meanLik[i_frac]).value());
    sum_partial_Evidence =
        sum_partial_Evidence + expected_partial_evidence_by_logLik[i_frac];
  }

  for (std::size_t i_frac = 0; i_frac < size(beta0); ++i_frac)
    for (std::size_t i_beta = 0; i_beta < size(beta0[i_frac]); ++i_beta)
      if (beta0[i_frac].back() == 1) {
        s.f << size(beta0) << s.sep << size(beta0[i_frac]) << s.sep << 0
            << s.sep << y[i_frac].nrows() << s.sep << beta0[i_frac][i_beta]
            << s.sep << 0 << s.sep << expected_meanLik[i_frac][i_beta] << s.sep
            << 0 << s.sep << partial_expected_evidence[i_frac] << s.sep
            << expected_partial_evidence_by_logLik[i_frac] << s.sep
            << expected_Evidence << s.sep << sum_partial_Evidence << "\n";
      }
}

template <class Prior, class Likelihood, class Variables, class DataType>
  requires(!has_conjugate<Prior, Likelihood, Variables, DataType>)
void report_model(save_Evidence &, Prior const &, Likelihood const &,
                  const DataType &, const Variables &,
                  by_fraction<by_beta<double>> const &) {}
}
namespace deprecated {
template <class FunctionTable, class Parameters, class T>
void report(FunctionTable &&, std::size_t iter, save_Evidence &s,
            cuevi_mcmc<Parameters> const &data, T const &...) {
  if (iter % s.save_every == 0) {

    auto meanLik = mean_logL(data);
    auto meanPrior = mean_logP(data);
    auto varLik = var_logL(data, meanLik);
    auto partial_Evidence2 = calculate_Evidence(data.beta, meanLik, varLik);
    auto partial_Evidence1 = calculate_Evidence(data.beta, meanLik);

    double Evidence1 = 0;
    double Evidence2 = 0;
    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac) {
      Evidence1 += partial_Evidence1[i_frac];
      Evidence2 += partial_Evidence2[i_frac];
    }

    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
      for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
        if (data.beta[i_frac].back() == 1) {
          s.f << size(data.beta) << s.sep << size(data.beta[i_frac]) << s.sep
              << iter << s.sep << data.nsamples[i_frac] << s.sep
              << data.beta[i_frac][i_beta] << s.sep << meanPrior[i_frac][i_beta]
              << s.sep << meanLik[i_frac][i_beta] << s.sep
              << varLik[i_frac][i_beta] << s.sep << partial_Evidence1[i_frac]
              << s.sep << partial_Evidence2[i_frac] << s.sep << Evidence1
              << s.sep << Evidence2 << "\n";
        }
  }
}

template <class FunctionTable, class Parameters, class... saving, class... T>
void report_all(FunctionTable &&f, std::size_t iter,
                save_mcmc<Parameters, saving...> &s,
                cuevi_mcmc<Parameters> const &data, T const &...ts) {
  (report(f, iter, static_cast<saving &>(s), data, ts...), ..., 1);
}

template <class Parameters>
void check_sanity(std::size_t iter, cuevi_mcmc<Parameters> const &data) {
  for (std::size_t iw = 0; iw < data.walkers.size(); ++iw)
    for (std::size_t i_frac = 0; i_frac < data.walkers[iw].size(); ++i_frac)
      for (std::size_t ib = 0; ib < data.walkers[iw][i_frac].size(); ++ib)
        if ((data.walkers[iw][i_frac][ib].logL < -1e8) &&
            (i_frac > 0 || data.beta[0][ib] > 0) && iter % 100 > 10)
          std::cerr << "report " << iter << " " << iw << " " << i_frac << " "
                    << ib << " " << data.walkers[iw][i_frac][ib].logL << "\n";
}

struct step_stretch_cuevi_mcmc_per_walker {
  friend std::string ToString(step_stretch_cuevi_mcmc_per_walker) {
    return "step_stretch_cuevi_mcmc_per_walker";
  }

  template <class FunctionTable, class Observer, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                  Observer &, ensemble<mt_64i> &mt,
                  std::vector<std::uniform_real_distribution<double>> &rdist,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, std::size_t n_par,
                  std::size_t i, std::size_t iw, std::size_t jw, std::size_t ib,
                  std::size_t i_fr) const {
    if (current.is_active[i_fr][ib] == 1) {
      auto r = rdist[i](mt[i]);
      // candidate[ib].walkers[iw].
      auto [ca_par, z] =
          stretch_move(mt[i], rdist[i], current.walkers[iw][i_fr][ib].parameter,
                       current.walkers[jw][i_fr][ib].parameter);

      auto ca_logPa_ = logPrior(prior, ca_par);
      auto ca_logL_0 = i_fr > 0 ? f.f(logLikelihood_f{}, lik, ca_par,
                                      y[i_fr - 1], x[i_fr - 1])
                                : Maybe_error(0.0);
      auto ca_logL_1 = f.f(logLikelihood_f{}, lik, ca_par, y[i_fr], x[i_fr]);
      if (is_valid(ca_logPa_) && is_valid(ca_logL_0) && is_valid(ca_logL_1)) {
        auto ca_logPa = ca_logPa_.value();
        auto ca_logP0 = ca_logPa_.value() + ca_logL_0.value();
        auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();

        auto dthLogL = ca_logP0 - current.walkers[iw][i_fr][ib].logP +
                       current.beta[i_fr][ib] *
                           (ca_logL0 - current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        if (pJump >= r) {
          if (i_fr + 1 < size(current.beta) &&
              (current.beta[i_fr][ib] == 1.0)) {
            auto ca_logL_2 =
                f.f(logLikelihood_f{}, lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
            if ((ca_logL_2)) {
              auto ca_logP1 = ca_logPa + ca_logL_1.value();
              auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
              // if ((current.beta[i_fr][ib]>0)&&(ca_logL0<-1e6||ca_logL1<-1e6))
              // {
              //     std::cerr<<"here guy\n";
              // }
              current.walkers[iw][i_fr][ib].parameter = ca_par;
              current.walkers[iw][i_fr][ib].logPa = ca_logPa;
              current.walkers[iw][i_fr][ib].logP = ca_logP0;
              current.walkers[iw][i_fr][ib].logL = ca_logL0;
              current.walkers[iw][i_fr + 1][0].parameter = ca_par;
              current.walkers[iw][i_fr + 1][0].logPa = ca_logPa;
              current.walkers[iw][i_fr + 1][0].logP = ca_logP1;
              current.walkers[iw][i_fr + 1][0].logL = ca_logL1;
            }
          } else {
            current.walkers[iw][i_fr][ib].parameter = ca_par;
            current.walkers[iw][i_fr][ib].logPa = ca_logPa;
            current.walkers[iw][i_fr][ib].logP = ca_logP0;
            current.walkers[iw][i_fr][ib].logL = ca_logL0;
          }
        }
      }

      else {

        std::cerr << iw << " " << jw << " "
                  << "i_fr=" << i_fr << " "
                  << "ib=" << ib << " " << ca_logPa_.error()()
                  << ca_logL_0.error()() << ca_logL_1.error()() << "\n";
      }
    }
  }
};

struct step_stretch_cuevi_mcmc {
  friend std::string ToString(step_stretch_cuevi_mcmc) {
    return "step_stretch_cuevi_mcmc";
  }

  template <class FunctionTable, class Observer, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                  Observer &obs, ensemble<mt_64i> &mt, Prior const &prior,
                  Likelihood const &lik, const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  double alpha_stretch = 2) const {
    assert(current.beta.size() == num_betas(current));
    auto n_walkers = size(current.walkers);

    auto n_par = size(current.walkers[0][0][0].parameter);

    std::uniform_int_distribution<std::size_t> uniform_walker(0, n_walkers / 2 -
                                                                     1);
    std::vector<std::uniform_int_distribution<std::size_t>> udist(
        n_walkers, uniform_walker);

    std::uniform_real_distribution<double> uniform_stretch_zdist(
        1.0 / alpha_stretch, alpha_stretch);
    std::vector<std::uniform_real_distribution<double>> zdist(
        n_walkers, uniform_stretch_zdist);

    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                              uniform_real);

    for (bool half : {false, true})
#pragma omp parallel for
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto iw = half ? i + n_walkers / 2 : i;
        auto j = udist[i](mt[i]);
        auto jw = half ? j : j + n_walkers / 2;
        if (current.is_active[0][0] == 1)
          f.fork(var::I_thread(i))
              .f(step_stretch_cuevi_mcmc_per_walker{}, current, obs, mt, rdist,
                 prior, lik, y, x, n_par, i, iw, jw, 0, 0);
        for (std::size_t i_fr = 0; i_fr < size(current.beta); ++i_fr) {
          for (std::size_t ib = 1; ib < size(current.beta[i_fr]); ++ib)
            if (current.is_active[i_fr][ib] == 1)
              f.fork(var::I_thread(i))
                  .f(step_stretch_cuevi_mcmc_per_walker{}, current, obs, mt,
                     rdist, prior, lik, y, x, n_par, i, iw, jw, ib, i_fr);
        }
      }
  }

  template <class FunctionTable, class Observer, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)
  
  void operator()(FunctionTable &&, cuevi_mcmc<Parameters> &current, Observer &,
                  ensemble<mt_64i> &mt,
                  std::vector<std::uniform_real_distribution<double>> &rdist,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, std::size_t n_par,
                  std::size_t i, std::size_t iw, std::size_t jw, std::size_t ib,
                  std::size_t i_fr) {
    if (current.is_active[i_fr][ib] == 1) {
      auto r = rdist[i](mt[i]);
      // candidate[ib].walkers[iw].
      auto [ca_par, z] =
          stretch_move(mt[i], rdist[i], current.walkers[iw][i_fr][ib].parameter,
                       current.walkers[jw][i_fr][ib].parameter);

      auto ca_logPa_ = logPrior(prior, ca_par);
      auto ca_logL_0 =
          i_fr > 0 ? logLikelihood(lik, ca_par, y[i_fr - 1], x[i_fr - 1])
                   : Maybe_error(0.0);
      auto ca_logL_1 = logLikelihood(lik, ca_par, y[i_fr], x[i_fr]);
      if (is_valid(ca_logPa_) && is_valid(ca_logL_0) && is_valid(ca_logL_1)) {
        auto ca_logPa = ca_logPa_.value();
        auto ca_logP0 = ca_logPa_.value() + ca_logL_0.value();
        auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();

        auto dthLogL = ca_logP0 - current.walkers[iw][i_fr][ib].logP +
                       current.beta[i_fr][ib] *
                           (ca_logL0 - current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        if (pJump >= r) {
          if (i_fr + 1 < size(current.beta) &&
              (current.beta[i_fr][ib] == 1.0)) {
            auto ca_logL_2 =
                logLikelihood(lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
            if ((ca_logL_2)) {
              auto ca_logP1 = ca_logPa + ca_logL_1.value();
              auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
              // if ((current.beta[i_fr][ib]>0)&&(ca_logL0<-1e6||ca_logL1<-1e6))
              // {
              //     std::cerr<<"here guy\n";
              // }
              current.walkers[iw][i_fr][ib].parameter = ca_par;
              current.walkers[iw][i_fr][ib].logPa = ca_logPa;
              current.walkers[iw][i_fr][ib].logP = ca_logP0;
              current.walkers[iw][i_fr][ib].logL = ca_logL0;
              current.walkers[iw][i_fr + 1][0].parameter = ca_par;
              current.walkers[iw][i_fr + 1][0].logPa = ca_logPa;
              current.walkers[iw][i_fr + 1][0].logP = ca_logP1;
              current.walkers[iw][i_fr + 1][0].logL = ca_logL1;
            }
          } else {
            current.walkers[iw][i_fr][ib].parameter = ca_par;
            current.walkers[iw][i_fr][ib].logPa = ca_logPa;
            current.walkers[iw][i_fr][ib].logP = ca_logP0;
            current.walkers[iw][i_fr][ib].logL = ca_logL0;
          }
        }
      }

      else {

        std::cerr << iw << " " << jw << " "
                  << "i_fr=" << i_fr << " "
                  << "ib=" << ib << " " << ca_logPa_.error()()
                  << ca_logL_0.error()() << ca_logL_1.error()() << "\n";
      }
    }
  }
};

using DataIndexes = std::vector<std::size_t>;

struct fractioner {
  auto operator()(const Matrix<double> &y, const Matrix<double> &x, mt_64i &mt,
                  std::size_t num_parameters, double n_points_per_decade_beta,
                  double n_points_per_decade_fraction, double stops_at,
                  bool includes_zero) const {
    assert(y.nrows() == x.nrows());
    std::size_t num_samples = size(y);
    // generate_random_Indexes differs from repository CumulativeEvidence
    auto indexes = generate_random_Indexes(mt, num_samples, 2 * num_parameters,
                                           n_points_per_decade_fraction);
    auto n_frac = size(indexes);
    by_fraction<Matrix<double>> y_out(n_frac);
    by_fraction<Matrix<double>> x_out(n_frac);
    for (std::size_t i = 0; i < n_frac; ++i) {
      auto n = size(indexes[i]);
      auto ii = indexes[i];
      Matrix<double> yi(n, 1, false);
      Matrix<double> xi(n, x.ncols(), false);

      for (std::size_t j = 0; j < n; ++j) {
        yi[j] = y[ii[j]];
        for (std::size_t k = 0; k < x.ncols(); ++k)
          xi(j, k) = x(ii[j], k);
      }
      y_out[i] = std::move(yi);
      x_out[i] = std::move(xi);
    }

    auto beta0 = get_beta_list(n_points_per_decade_beta,
                               stops_at * num_samples /
                                   (n_frac > 1 ? size(indexes[0]) : 1),
                               includes_zero);
    by_beta<double> betan = {0, 1};
    by_fraction<by_beta<double>> beta(n_frac, betan);
    beta[0] = std::move(beta0);

    return std::tuple(std::move(y_out), std::move(x_out), std::move(beta));
  }
};

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

inline auto init_mcmc2(FunctionTable &&f, mt_64i &mt, const Prior &prior,
                       const Likelihood &lik, const by_fraction<DataType> &y,
                       const by_fraction<Variables> &x) {
  auto prior_sampler = sampler(prior);
  auto par = sample(mt, prior_sampler);
  auto logP = logPrior(prior, par);
  auto logL = logLikelihood(f, lik, par, y[0], x[0]);
  auto logPa = logP;
  while (!(logP) || !(logL)) {
    // std::cerr<<"\npar\ņ"<<par;
    // std::cerr<<"\nlogL\ņ"<<logL;

    par = sample(mt, prior_sampler);
    logP = logPrior(prior, par);
    logL = logLikelihood(f, lik, par, y[0], x[0]);
  }
  return mcmc2<Parameters>{
      mcmc<Parameters>{std::move(par), logP.value(), logL.value()},
      logPa.value()};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

inline auto init_mcmc_resample(FunctionTable &&f, ensemble<mt_64i> &mt,
                               cuevi_mcmc<Parameters> &current,
                               const Prior &prior, const Likelihood &lik,
                               const by_fraction<DataType> &y,
                               const by_fraction<Variables> &x) {
  auto prior_sampler = sampler(prior);
  auto n_walkers = current.walkers.size();
  auto n_frac = current.beta.size();
  for (std::size_t half = 0; half < 2; ++half)
#pragma omp parallel for
    for (std::size_t iiw = 0; iiw < n_walkers / 2; ++iiw) {
      auto iw = iiw + half * n_walkers / 2;
      for (std::size_t i_fr = 0; i_fr < n_frac; ++i_fr)
        for (std::size_t ib = 0; ib < current.beta[i_fr].size(); ++ib) {
          if (i_fr == 0 || ib > 0) {
            auto par = std::decay_t<decltype(sample(mt[iiw], prior_sampler))>{};
            auto logP = Maybe_error<double>(error_message("not_init"));
            auto logL_0 = Maybe_error<double>(error_message("not_init"));
            auto logL_1 = Maybe_error<double>(error_message("not_init"));
            auto logL_2 = Maybe_error<double>(error_message("not_init"));

            bool share_with_next =
                (i_fr + 1 < n_frac) && (ib + 1 == current.beta[i_fr].size());
            while (!(logP.valid() && logL_0.valid() && logL_1.valid() &&
                     logL_2.valid())) {
              par = sample(mt[iiw], prior_sampler);
              logP = logPrior(prior, par);
              if (i_fr > 0)
                logL_0 = f.fork(var::I_thread(iiw))
                             .f(logLikelihood_f{}, lik, par, y[i_fr - 1],
                                x[i_fr - 1]);
              else
                logL_0 = Maybe_error(0.0);

              logL_1 = f.fork(var::I_thread(iiw))
                           .f(logLikelihood_f{}, lik, par, y[i_fr], x[i_fr]);

              if (share_with_next)
                logL_2 = f.fork(var::I_thread(iiw))
                             .f(logLikelihood_f{}, lik, par, y[i_fr + 1],
                                x[i_fr + 1]);
              else
                logL_2 = Maybe_error(0.0);
            }
            auto logPa = logP.value();
            auto logP0 = logP.value() + logL_0.value();
            auto logL0 = logL_1.value() - logL_0.value();
            current.walkers[iw][i_fr][ib].parameter = par;
            current.walkers[iw][i_fr][ib].logPa = logPa;
            current.walkers[iw][i_fr][ib].logP = logP0;
            current.walkers[iw][i_fr][ib].logL = logL0;
            if (share_with_next) {
              auto logP1 = logPa + logL_1.value();
              auto logL1 = logL_2.value() - logL_1.value();
              current.walkers[iw][i_fr + 1][0].parameter = par;
              current.walkers[iw][i_fr + 1][0].logPa = logPa;
              current.walkers[iw][i_fr + 1][0].logP = logP1;
              current.walkers[iw][i_fr + 1][0].logL = logL1;
            }
          }
        }
    }
  return current;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

inline auto init_cuevi_mcmc(FunctionTable &&f, std::size_t n_walkers,
                            by_beta<double> const &beta, ensemble<mt_64i> &mt,
                            Prior const &prior, Likelihood const &lik,
                            const by_fraction<DataType> &y,
                            const by_fraction<Variables> &x) {
  by_fraction<std::size_t> nsamples_out(1, size(y[0]));
  by_fraction<by_beta<double>> beta_out(1, beta);
  auto active_beta = by_beta<int>(beta.size(), 1);
  by_fraction<by_beta<int>> active_out(1, active_beta);

  ensemble<by_fraction<by_beta<std::size_t>>> i_walker(
      n_walkers,
      by_fraction<by_beta<std::size_t>>(1, by_beta<std::size_t>(beta.size())));
  ensemble<by_fraction<by_beta<mcmc2<Parameters>>>> walker(
      n_walkers, by_fraction<by_beta<mcmc2<Parameters>>>(
                     1, by_beta<mcmc2<Parameters>>(beta.size())));

  for (std::size_t half = 0; half < 2; ++half)
#pragma omp parallel for
    for (std::size_t iiw = 0; iiw < n_walkers / 2; ++iiw) {
      auto iw = iiw + half * n_walkers / 2;
      for (std::size_t i = 0; i < beta.size(); ++i) {
        i_walker[iw][0][i] = iw + (beta.size() - i - 1) * n_walkers;
        walker[iw][0][i] =
            init_mcmc2(f.fork(var::I_thread(iiw)), mt[iiw], prior, lik, y, x);
      }
    }

  return cuevi_mcmc<Parameters>{nsamples_out, beta_out, walker, i_walker,
                                active_out};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

inline auto init_cuevi_mcmc_all(FunctionTable &&f, std::size_t n_walkers,
                                by_fraction<by_beta<double>> const &beta_out,
                                ensemble<mt_64i> &mt, Prior const &prior,
                                Likelihood const &lik,
                                const by_fraction<DataType> &y,
                                const by_fraction<Variables> &x) {
  by_fraction<std::size_t> nsamples_out(beta_out.size());
  for (std::size_t i = 0; i < beta_out.size(); ++i)
    nsamples_out[i] = size(y[i]);
  by_fraction<by_beta<int>> active_out(beta_out.size());
  for (std::size_t i = 0; i < beta_out.size(); ++i)
    active_out[i] = by_beta<int>(beta_out[i].size(), 1);

  by_fraction<by_beta<std::size_t>> ii_walker(beta_out.size());
  for (std::size_t i = 0; i < beta_out.size(); ++i)
    ii_walker[i] = by_beta<std::size_t>(beta_out[i].size());
  ensemble<by_fraction<by_beta<std::size_t>>> i_walker(n_walkers, ii_walker);

  by_fraction<by_beta<mcmc2<Parameters>>> walker_i(beta_out.size());
  for (std::size_t i = 0; i < beta_out.size(); ++i)
    walker_i[i] = by_beta<mcmc2<Parameters>>(beta_out[i].size());
  ensemble<by_fraction<by_beta<mcmc2<Parameters>>>> walker(n_walkers, walker_i);

  for (std::size_t iw = 0; iw < n_walkers; ++iw) {
    i_walker[iw][0][0] = iw;
    for (std::size_t ib = 0; ib < beta_out[0].size(); ++ib) {
      i_walker[iw][0][ib] = iw + ib * n_walkers;
    }
    for (std::size_t i_frac = 1; i_frac < beta_out.size(); ++i_frac) {
      i_walker[iw][i_frac][0] = i_walker[iw][i_frac - 1].back();
      for (std::size_t ib = 1; ib < beta_out[i_frac].size(); ++ib) {
        i_walker[iw][i_frac][ib] = i_walker[iw][i_frac][ib - 1] + n_walkers;
      }
    }
  }
  auto current = cuevi_mcmc<Parameters>{nsamples_out, beta_out, walker,
                                        i_walker, active_out};

  return init_mcmc_resample(f, mt, current, prior, lik, y, x);
}

template <class Parameters>
std::size_t next_id_walker(const cuevi_mcmc<Parameters> &c) {
  std::size_t out = 0;
  for (auto &i_w : c.i_walkers)
    for (auto &i_wf : i_w)
      for (auto &i_wfb : i_wf)
        out = std::max(out, i_wfb);
  return out + 1;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

inline Maybe_error<bool>
calculate_Likelihoods_sample(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                             Prior const &prior, Likelihood const &lik,
                             const by_fraction<DataType> &y,
                             const by_fraction<Variables> &x, std::size_t iw,
                             std::size_t i_frac, std::size_t ib) {
  auto const &ca_par = current.walkers[iw][i_frac][ib].parameter;
  auto ca_logPa_ = logPrior(prior, ca_par);
  auto ca_logL_0 =
      i_frac > 0 ? logLikelihood(f, lik, ca_par, y[i_frac - 1], x[i_frac - 1])
                 : Maybe_error(0.0);
  auto ca_logL_1 = logLikelihood(f, lik, ca_par, y[i_frac], x[i_frac]);
  if (!(is_valid(ca_logPa_) && is_valid(ca_logL_0) && is_valid(ca_logL_1))) {
    return error_message(ca_logPa_.error()() + ca_logL_0.error()() +
                         ca_logL_1.error()());
  } else {
    auto ca_logPa = ca_logPa_.value();
    auto ca_logP0 = ca_logPa_.value() + ca_logL_0.value();
    auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();
    if (i_frac + 1 < size(current.walkers[iw]) &&
        (current.beta[i_frac][ib] == 1.0)) {
      auto ca_logL_2 =
          logLikelihood(f, lik, ca_par, y[i_frac + 1], x[i_frac + 1]);
      if (!(ca_logL_2))
        return ca_logL_2.error();
      else {
        // assert(test_equality(ca_par,current.walkers[iw][i_frac +
        // 1][0].parameter, eps));
        auto ca_logP1 = ca_logPa + ca_logL_1.value();
        auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
        current.walkers[iw][i_frac][ib].logPa = ca_logPa;
        current.walkers[iw][i_frac][ib].logP = ca_logP0;
        current.walkers[iw][i_frac][ib].logL = ca_logL0;
        current.walkers[iw][i_frac + 1][0].logPa = ca_logPa;
        current.walkers[iw][i_frac + 1][0].logP = ca_logP1;
        current.walkers[iw][i_frac + 1][0].logL = ca_logL1;
        return true;
      }
    } else {
      current.walkers[iw][i_frac][ib].logPa = ca_logPa;
      current.walkers[iw][i_frac][ib].logP = ca_logP0;
      current.walkers[iw][i_frac][ib].logL = ca_logL0;
      return true;
    }
  }
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
Maybe_error<cuevi_mcmc<Parameters>> calculate_current_Likelihoods(
    FunctionTable &&f, cuevi_mcmc<Parameters> &current, Prior const &prior,
    Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {

  for (std::size_t iw = 0; iw < current.walkers.size(); ++iw) {
    if (current.is_active[0][0] == 1) {
      auto res =
          calculate_Likelihoods_sample(f, current, prior, lik, y, x, iw, 0, 0);
      if (!res)
        return res.error();
    }
    for (std::size_t i_frac = 0; i_frac < current.walkers[iw].size(); ++i_frac)
      for (std::size_t ib = 1; ib < current.walkers[iw][i_frac].size(); ++ib) {
        if (current.is_active[i_frac][ib] == 1) {
          auto res = calculate_Likelihoods_sample(f, current, prior, lik, y, x,
                                                  iw, i_frac, ib);
          if (!res)
            return res.error();
        }
      }
  }
  return current;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
inline auto
create_new_walkers(FunctionTable &&f, const cuevi_mcmc<Parameters> &current,
                   ensemble<mt_64i> &mts, Prior const &prior,
                   Likelihood const &lik, const by_fraction<DataType> &y,
                   const by_fraction<Variables> &x) {

  auto n_walkers = current.walkers.size();
  auto sum_walkers = next_id_walker(current);
  ensemble<mcmc2<Parameters>> new_walkers(n_walkers);
  ensemble<std::size_t> new_i_walkers(n_walkers);
  for (std::size_t half = 0; half < 2; ++half)
    for (std::size_t i = 0; i < n_walkers / 2; ++i) {
      auto iw = i + half * n_walkers / 2;
      new_walkers[iw] = init_mcmc2(f, mts[i], prior, lik, y, x);
      new_i_walkers[iw] = sum_walkers + iw;
    }

  return std::tuple(new_walkers, new_i_walkers);
}

template <class Parameters>
auto get_soon_to_be_inactive_new_walkers_and_its_position(
    const cuevi_mcmc<Parameters> &current) {
  std::size_t i_frac = 0;
  std::size_t ib = 0;
  while (!current.is_active[i_frac][ib] == 1) {
    if (ib + 1 < current.is_active[i_frac].size())
      ++ib;
    else {
      ++i_frac;
      ib = 1;
    }
  }

  auto n_walkers = current.walkers.size();
  ensemble<mcmc2<Parameters>> new_walkers(n_walkers);
  ensemble<std::size_t> new_i_walkers(n_walkers);
  for (std::size_t iw = 0; iw < n_walkers; ++iw) {
    new_walkers[iw] = current.walkers[iw][i_frac][ib];
    new_i_walkers[iw] = current.i_walkers[iw][i_frac][ib];
  }
  return std::tuple(i_frac, ib, new_walkers, new_i_walkers);
}

template <class DataType, class Parameters>
void insert_new_walkers(cuevi_mcmc<Parameters> &current,
                        const by_fraction<by_beta<double>> &final_beta,
                        const by_fraction<DataType> &y,
                        std::size_t insert_i_frac, std::size_t insert_ib,
                        ensemble<mcmc2<Parameters>> &&new_walkers,
                        ensemble<std::size_t> &&new_i_walkers) {

  for (std::size_t iw = 0; iw < current.walkers.size(); ++iw) {
    if (insert_i_frac == 0)
      for (std::size_t ib = insert_ib; ib < current.walkers[iw][0].size();
           ++ib) {
        std::swap(current.walkers[iw][0][ib], new_walkers[iw]);
        std::swap(current.i_walkers[iw][0][ib], new_i_walkers[iw]);
      }
    for (std::size_t i_frac = std::max(1ul, insert_i_frac);
         i_frac < current.walkers[iw].size(); ++i_frac) {
      current.walkers[iw][i_frac][0] = current.walkers[iw][i_frac - 1].back();
      current.i_walkers[iw][i_frac][0] =
          current.i_walkers[iw][i_frac - 1].back();

      for (std::size_t ib = 1; ib < current.walkers[iw][i_frac].size(); ++ib) {
        std::swap(current.walkers[iw][i_frac][ib], new_walkers[iw]);
        std::swap(current.i_walkers[iw][i_frac][ib], new_i_walkers[iw]);
      }
    }
    auto i_frac = current.walkers[iw].size() - 1;
    auto ib = current.walkers[iw][i_frac].size() - 1;

    if (current.beta[i_frac][ib] < 1.0) {
      current.walkers[iw][i_frac].push_back(new_walkers[iw]);
      current.i_walkers[iw][i_frac].push_back(new_i_walkers[iw]);
      if (iw == 0) {
        current.beta[i_frac].push_back(final_beta[i_frac][ib + 1]);
        current.is_active[i_frac].push_back(1);
      }
    } else {
      current.walkers[iw].push_back(by_beta<mcmc2<Parameters>>(2));
      current.i_walkers[iw].push_back(by_beta<std::size_t>(2));
      current.walkers[iw][i_frac + 1][0] = current.walkers[iw][i_frac].back();
      current.i_walkers[iw][i_frac + 1][0] =
          current.i_walkers[iw][i_frac].back();

      current.walkers[iw][i_frac + 1][1] = new_walkers[iw];
      current.i_walkers[iw][i_frac + 1][1] = new_i_walkers[iw];

      if (iw == 0) {
        current.beta.push_back(by_beta<double>{final_beta[i_frac + 1][0],
                                               final_beta[i_frac + 1][1]});
        current.nsamples.push_back(size(y[i_frac + 1]));
        current.is_active.push_back(by_beta<int>{1, 1});
      }
    }
  }
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
Maybe_error<cuevi_mcmc<Parameters>> push_back_new_fraction(
    FunctionTable &&f, const cuevi_mcmc<Parameters> &current_old,
    ensemble<mt_64i> &mts, const by_fraction<by_beta<double>> &final_beta,
    std::size_t max_number_of_simultaneous_temperatures, Prior const &prior,
    Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {

  auto current = current_old;

  if (current.current_number_of_temperatures() <
      max_number_of_simultaneous_temperatures) {
    auto [new_walkers, new_i_walkers] =
        create_new_walkers(f, current, mts, prior, lik, y, x);

    insert_new_walkers(current, final_beta, y, 0ul, 0ul, std::move(new_walkers),
                       std::move(new_i_walkers));
  } else {
    auto [i_frac, ib, new_walkers, new_i_walkers] =
        get_soon_to_be_inactive_new_walkers_and_its_position(current);
    current.is_active[i_frac][ib] = 0;
    if (current.is_active[i_frac].size() == ib + 1)
      current.is_active[i_frac + 1][0] = 0;
    insert_new_walkers(current, final_beta, y, i_frac, ib,
                       std::move(new_walkers), std::move(new_i_walkers));
  }

  return calculate_current_Likelihoods(f, current, prior, lik, y, x);
}

struct thermo_cuevi_jump_mcmc {
  friend std::string ToString(thermo_cuevi_jump_mcmc) {
    return "thermo_cuevi_jump_mcmc";
  }

  template <class FunctionTable, class Observer, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, std::size_t iter,
                  cuevi_mcmc<Parameters> &current, Observer &obs, mt_64i &mt,
                  ensemble<mt_64i> &mts, Prior const &, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  std::size_t thermo_jumps_every) const {
    if (iter % (thermo_jumps_every) == 0) {
      std::uniform_real_distribution<double> uniform_real(0, 1);
      auto n_walkers = mts.size() * 2;
      auto n_par = current.walkers[0][0][0].parameter.size();

      WalkerIndexes shuffled_walker(n_walkers);
      std::iota(shuffled_walker.begin(), shuffled_walker.end(), 0);
      std::shuffle(shuffled_walker.begin(), shuffled_walker.end(), mt);
      std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                                uniform_real);

#pragma omp parallel for // not currently working
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto iw = shuffled_walker[i];
        auto jw = shuffled_walker[i + n_walkers / 2];

        if (size(current.beta) == 1)
          for (std::size_t i_fr = 0; i_fr < 1; ++i_fr) {
            for (std::size_t ib = 0; ib < current.beta[i_fr].size() - 1; ++ib) {
              if (current.is_active[i_fr][ib] == 1) {

                auto r = rdist[i](mts[i]);
                double logA = calc_logA(current.beta[i_fr][ib],
                                        current.beta[i_fr][ib + 1],
                                        current.walkers[iw][i_fr][ib].logL,
                                        current.walkers[jw][i_fr][ib + 1].logL);
                auto pJump = std::min(1.0, std::exp(logA));
                observe_thermo_jump_mcmc(
                    obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                    current.walkers[jw][i_fr][ib + 1].parameter,
                    current.walkers[iw][i_fr][ib].logL,
                    current.walkers[jw][i_fr][ib + 1].logL,
                    -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                    logA, pJump, r, pJump > r);
                if (pJump > r) {
                  std::swap(current.walkers[iw][i_fr][ib],
                            current.walkers[jw][i_fr][ib + 1]);
                  std::swap(current.i_walkers[iw][i_fr][ib],
                            current.i_walkers[jw][i_fr][ib + 1]);
                }
              }
            }
          }
        else {
          for (std::size_t i_fr = 0; i_fr < 1; ++i_fr) {
            for (std::size_t ib = 0; ib + 2 < current.beta[i_fr].size(); ++ib) {
              if (current.is_active[i_fr][ib] == 1) {

                auto r = rdist[i](mts[i]);
                double logA = calc_logA(current.beta[i_fr][ib],
                                        current.beta[i_fr][ib + 1],
                                        current.walkers[iw][i_fr][ib].logL,
                                        current.walkers[jw][i_fr][ib + 1].logL);
                auto pJump = std::min(1.0, std::exp(logA));
                observe_thermo_jump_mcmc(
                    obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                    current.walkers[jw][i_fr][ib + 1].parameter,
                    current.walkers[iw][i_fr][ib].logL,
                    current.walkers[jw][i_fr][ib + 1].logL,
                    -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                    logA, pJump, r, pJump > r);
                if (pJump > r) {
                  std::swap(current.walkers[iw][i_fr][ib],
                            current.walkers[jw][i_fr][ib + 1]);
                  std::swap(current.i_walkers[iw][i_fr][ib],
                            current.i_walkers[jw][i_fr][ib + 1]);
                }
              }
            }
            for (std::size_t ib = current.beta[i_fr].size() - 2;
                 ib < current.beta[i_fr].size() - 1; ++ib) {
              if (current.is_active[i_fr][ib] == 1) {

                auto r = rdist[i](mts[i]);
                double logA = calc_logA(current.beta[i_fr][ib],
                                        current.beta[i_fr][ib + 1],
                                        current.walkers[iw][i_fr][ib].logL,
                                        current.walkers[jw][i_fr][ib + 1].logL);
                auto pJump = std::min(1.0, std::exp(logA));
                observe_thermo_jump_mcmc(
                    obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                    current.walkers[jw][i_fr][ib + 1].parameter,
                    current.walkers[iw][i_fr][ib].logL,
                    current.walkers[jw][i_fr][ib + 1].logL,
                    -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                    logA, pJump, r, pJump > r);
                if (pJump > r) {
                  auto ca_par = current.walkers[iw][i_fr][ib].parameter;
                  auto ca_logL1 =
                      logLikelihood(f.fork(var::I_thread(i)), lik, ca_par,
                                    y[i_fr + 1], x[i_fr + 1]);
                  if (ca_logL1) {
                    auto ca_logPa = current.walkers[iw][i_fr][ib].logPa;
                    auto ca_logP = current.walkers[iw][i_fr][ib].logP;
                    auto ca_logL0 = current.walkers[iw][i_fr][ib].logL;
                    std::swap(current.walkers[iw][i_fr][ib],
                              current.walkers[jw][i_fr][ib + 1]);
                    std::swap(current.i_walkers[iw][i_fr][ib],
                              current.i_walkers[jw][i_fr][ib + 1]);
                    if (size(current.beta) > 1) {
                      current.walkers[jw][i_fr + 1][0].parameter = ca_par;
                      current.walkers[jw][i_fr + 1][0].logPa = ca_logPa;
                      current.walkers[jw][i_fr + 1][0].logP =
                          ca_logP + ca_logL0;
                      current.walkers[jw][i_fr + 1][0].logL =
                          ca_logL1.value() - ca_logL0 - ca_logP + ca_logPa;
                    }
                  }
                }
              }
            }
          }
          for (std::size_t i_fr = 1; i_fr + 1 < current.beta.size(); ++i_fr) {
            if (current.beta[i_fr].size() < 3) {
              for (std::size_t ib = 0; ib + 1 < current.beta[i_fr].size();
                   ++ib) {
                if (current.is_active[i_fr][ib] == 1) {

                  auto r = rdist[i](mts[i]);
                  double logA = calc_logA(
                      current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL);
                  auto pJump = std::min(1.0, std::exp(logA));
                  observe_thermo_jump_mcmc(
                      obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                      current.walkers[jw][i_fr][ib + 1].parameter,
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL,
                      -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                      logA, pJump, r, pJump > r);
                  if (pJump > r) {
                    auto ca_par_1 = current.walkers[iw][i_fr][ib].parameter;
                    auto ca_logL_11 =
                        logLikelihood(f.fork(var::I_thread(i)), lik, ca_par_1,
                                      y[i_fr + 1], x[i_fr + 1]);
                    auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                    auto ca_logL_00 =
                        i_fr == 1
                            ? Maybe_error<double>{0.0}
                            : logLikelihood(f.fork(var::I_thread(i)), lik,
                                            ca_par_0, y[i_fr - 2], x[i_fr - 2]);
                    if (is_valid(ca_logL_11) && is_valid(ca_logL_00)) {
                      auto ca_logPa_1 = current.walkers[iw][i_fr][ib].logPa;
                      auto ca_logP_1 = current.walkers[iw][i_fr][ib].logP;
                      auto ca_logL_1 = current.walkers[iw][i_fr][ib].logL;
                      auto ca_logPa_0 = current.walkers[jw][i_fr][ib + 1].logPa;
                      auto ca_logP_0 = current.walkers[jw][i_fr][ib + 1].logP;
                      auto ca_logL_0 = current.walkers[jw][i_fr][ib + 1].logL;
                      std::swap(current.walkers[iw][i_fr][ib],
                                current.walkers[jw][i_fr][ib + 1]);
                      std::swap(current.i_walkers[iw][i_fr][ib],
                                current.i_walkers[jw][i_fr][ib + 1]);
                      current.walkers[jw][i_fr + 1][0].parameter = ca_par_1;
                      current.walkers[jw][i_fr + 1][0].logPa = ca_logPa_1;
                      current.walkers[jw][i_fr + 1][0].logP =
                          ca_logP_1 + ca_logL_1;
                      current.walkers[jw][i_fr + 1][0].logL =
                          ca_logL_11.value() - ca_logL_1 - ca_logP_1 +
                          ca_logPa_1;
                      auto ib0 = current.beta[i_fr - 1].size() - 1;
                      current.walkers[iw][i_fr - 1][ib0].parameter = ca_par_0;
                      current.walkers[iw][i_fr - 1][ib0].logPa = ca_logPa_0;
                      current.walkers[iw][i_fr - 1][ib0].logP =
                          ca_logPa_0 + ca_logL_00.value();
                      current.walkers[iw][i_fr - 1][ib0].logL =
                          ca_logP_0 - ca_logPa_0 - ca_logL_00.value();
                    }
                  }
                }
              }

            } else {
              for (std::size_t ib = 0; ib < 1; ++ib) {
                if (current.is_active[i_fr][ib] == 1) {

                  auto r = rdist[i](mts[i]);
                  double logA = calc_logA(
                      current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL);
                  auto pJump = std::min(1.0, std::exp(logA));
                  observe_thermo_jump_mcmc(
                      obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                      current.walkers[jw][i_fr][ib + 1].parameter,
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL,
                      -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                      logA, pJump, r, pJump > r);
                  if (pJump > r) {

                    auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                    auto ca_logL_00 =
                        i_fr == 1
                            ? Maybe_error<double>{0.0}
                            : logLikelihood(f.fork(var::I_thread(i)), lik,
                                            ca_par_0, y[i_fr - 2], x[i_fr - 2]);
                    if (ca_logL_00) {
                      auto ca_logPa_0 = current.walkers[jw][i_fr][ib + 1].logPa;
                      auto ca_logP_0 = current.walkers[jw][i_fr][ib + 1].logP;
                      auto ca_logL_0 = current.walkers[jw][i_fr][ib + 1].logL;
                      std::swap(current.walkers[iw][i_fr][ib],
                                current.walkers[jw][i_fr][ib + 1]);
                      std::swap(current.i_walkers[iw][i_fr][ib],
                                current.i_walkers[jw][i_fr][ib + 1]);
                      auto ib0 = current.beta[i_fr - 1].size() - 1;
                      current.walkers[iw][i_fr - 1][ib0].parameter = ca_par_0;
                      current.walkers[iw][i_fr - 1][ib0].logPa = ca_logPa_0;
                      current.walkers[iw][i_fr - 1][ib0].logP =
                          ca_logPa_0 + ca_logL_00.value();
                      current.walkers[iw][i_fr - 1][ib0].logL =
                          ca_logP_0 - ca_logPa_0 - ca_logL_00.value();
                    }
                  }
                }
              }
              for (std::size_t ib = 1; ib + 2 < current.beta[i_fr].size();
                   ++ib) {
                if (current.is_active[i_fr][ib] == 1) {

                  auto r = rdist[i](mts[i]);
                  double logA = calc_logA(
                      current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL);
                  auto pJump = std::min(1.0, std::exp(logA));
                  observe_thermo_jump_mcmc(
                      obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                      current.walkers[jw][i_fr][ib + 1].parameter,
                      current.walkers[iw][i_fr][ib].logL,
                      current.walkers[jw][i_fr][ib + 1].logL,
                      -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                      logA, pJump, r, pJump > r);
                  if (pJump > r) {
                    std::swap(current.walkers[iw][i_fr][ib],
                              current.walkers[jw][i_fr][ib + 1]);
                    std::swap(current.i_walkers[iw][i_fr][ib],
                              current.i_walkers[jw][i_fr][ib + 1]);
                  }
                }
              }
            }
          }
          for (std::size_t i_fr = std::max(1ul, current.beta.size() - 1);
               i_fr < current.beta.size(); ++i_fr) {
            for (std::size_t ib = 0; ib < 1; ++ib) {
              if (current.is_active[i_fr][ib] == 1) {

                auto r = rdist[i](mts[i]);
                double logA = calc_logA(current.beta[i_fr][ib],
                                        current.beta[i_fr][ib + 1],
                                        current.walkers[iw][i_fr][ib].logL,
                                        current.walkers[jw][i_fr][ib + 1].logL);
                auto pJump = std::min(1.0, std::exp(logA));
                observe_thermo_jump_mcmc(
                    obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                    current.walkers[jw][i_fr][ib + 1].parameter,
                    current.walkers[iw][i_fr][ib].logL,
                    current.walkers[jw][i_fr][ib + 1].logL,
                    -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                    logA, pJump, r, pJump > r);
                if (pJump > r) {

                  auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                  auto ca_logL_00 =
                      i_fr == 1
                          ? Maybe_error<double>{0.0}
                          : logLikelihood(f.fork(var::I_thread(i)), lik,
                                          ca_par_0, y[i_fr - 2], x[i_fr - 2]);
                  if (ca_logL_00) {
                    auto ca_logPa_0 = current.walkers[jw][i_fr][ib + 1].logPa;
                    auto ca_logP_0 = current.walkers[jw][i_fr][ib + 1].logP;
                    auto ca_logL_0 = current.walkers[jw][i_fr][ib + 1].logL;
                    std::swap(current.walkers[iw][i_fr][ib],
                              current.walkers[jw][i_fr][ib + 1]);
                    std::swap(current.i_walkers[iw][i_fr][ib],
                              current.i_walkers[jw][i_fr][ib + 1]);
                    auto ib0 = current.beta[i_fr - 1].size() - 1;
                    auto cai_logP = ca_logPa_0 + ca_logL_00.value();
                    auto cai_logL = ca_logP_0 - ca_logPa_0 - ca_logL_00.value();
                    current.walkers[iw][i_fr - 1][ib0] = mcmc2<Parameters>{
                        mcmc<Parameters>{ca_par_0, cai_logP, cai_logL},
                        ca_logPa_0};
                  }
                }
              }
            }
            for (std::size_t ib = 1; ib + 1 < current.beta[i_fr].size(); ++ib) {
              if (current.is_active[i_fr][ib] == 1) {

                auto r = rdist[i](mts[i]);
                double logA = calc_logA(current.beta[i_fr][ib],
                                        current.beta[i_fr][ib + 1],
                                        current.walkers[iw][i_fr][ib].logL,
                                        current.walkers[jw][i_fr][ib + 1].logL);
                auto pJump = std::min(1.0, std::exp(logA));
                observe_thermo_jump_mcmc(
                    obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                    current.walkers[jw][i_fr][ib + 1].parameter,
                    current.walkers[iw][i_fr][ib].logL,
                    current.walkers[jw][i_fr][ib + 1].logL,
                    -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]),
                    logA, pJump, r, pJump > r);
                if (pJump > r) {
                  std::swap(current.walkers[iw][i_fr][ib],
                            current.walkers[jw][i_fr][ib + 1]);
                  std::swap(current.i_walkers[iw][i_fr][ib],
                            current.i_walkers[jw][i_fr][ib + 1]);
                }
              }
            }
          }
        }
      }
    }
  }
};

using FractionIndexes = std::vector<std::size_t>;

struct calculate_cuevi_walker {
  friend std::string ToString(calculate_cuevi_walker) {
    return "calculate_cuevi_walker";
  }
  template <class FunctionTable, class Prior, class Likelihood, class Variables,
            class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  Maybe_error<std::pair<mcmc2<Parameters>, mcmc2<Parameters>>>
  operator()(FunctionTable &&f, mcmc2<Parameters> &current,
             by_fraction<by_beta<double>> const &beta, std::size_t i_frac,

             std::size_t ib, Prior const &, Likelihood const &lik,
             const by_fraction<DataType> &y, const by_fraction<Variables> &x) {
    auto const &ca_par = current.parameter;
    auto ca_logPa_ = current.logPa;
    auto ca_logL_0 = i_frac > 0 ? f.f(logLikelihood_f{}, lik, ca_par,
                                      y[i_frac - 1], x[i_frac - 1])
                                : Maybe_error(0.0);
    auto ca_logL_1 = f.f(logLikelihood_f{}, lik, ca_par, y[i_frac], x[i_frac]);
    if (!(is_valid(ca_logL_0) && is_valid(ca_logL_1))) {
      return error_message(ca_logL_0.error()() + ca_logL_1.error()());
    } else {
      auto ca_logPa = ca_logPa_;
      auto ca_logP0 = ca_logPa_ + ca_logL_0.value();
      auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();
      mcmc2<Parameters> out0;
      mcmc2<Parameters> out1;
      if (beta[i_frac][ib] == 1.0) {
        out0.parameter = ca_par;
        out0.logPa = ca_logPa;
        out0.logP = ca_logP0;
        out0.logL = ca_logL0;

        if (i_frac + 1 < size(y)) {
          auto ca_logL_2 =
              f.f(logLikelihood_f{}, lik, ca_par, y[i_frac + 1], x[i_frac + 1]);
          if (!(ca_logL_2))
            return ca_logL_2.error();
          else {
            auto ca_logP1 = ca_logPa + ca_logL_1.value();
            auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
            out1.parameter = ca_par;
            out1.logPa = ca_logPa;
            out1.logP = ca_logP1;
            out1.logL = ca_logL1;
          }
        }
      } else {
        out1.parameter = ca_par;
        out1.logPa = ca_logPa;
        out1.logP = ca_logP0;
        out1.logL = ca_logL0;
        if ((i_frac > 0) && (beta[i_frac][ib] == 0.0)) {
          auto ca_logL_00 = i_frac > 1 ? f.f(logLikelihood_f{}, lik, ca_par,
                                             y[i_frac - 2], x[i_frac - 2])
                                       : Maybe_error<double>(0.0);
          if (!(ca_logL_00))
            return ca_logL_00.error();
          else {
            auto ca_logP00 = ca_logPa + ca_logL_00.value();
            auto ca_logL00 = ca_logL_0.value() - ca_logL_00.value();
            out0.parameter = ca_par;
            out0.logPa = ca_logPa;
            out0.logP = ca_logP00;
            out0.logL = ca_logL00;
          }
        }
      }
      return std::pair(std::move(out0), std::move(out1));
    }
  }
};
struct thermo_cuevi_randomized_jump_mcmc_per_walker {
  friend std::string ToString(thermo_cuevi_randomized_jump_mcmc_per_walker) {
    return "thermo_cuevi_randomized_jump_mcmc";
  }
  template <class FunctionTable, class Prior, class Likelihood, class Variables,
            class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                  std::size_t iw, std::size_t jw, std::size_t i_frac_origin,
                  std::size_t i_frac_destination, Prior const &prior,
                  Likelihood const &lik, const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, double r) {

    if (i_frac_origin > i_frac_destination)
      (*this)(std::forward<FunctionTable>(f), current, jw, iw,
              i_frac_destination, i_frac_origin, prior, lik, y, x, r);
    else {
      auto ib_origin = i_frac_origin > 0 ? 0ul : current.beta[0].size() - 2;
      auto ib_destination = current.beta[i_frac_destination].size() - 1;
      mcmc2<Parameters> &origin = current.walkers[iw][i_frac_origin][ib_origin];
      mcmc2<Parameters> &destination =
          current.walkers[jw][i_frac_destination][ib_destination];
      auto Maybe_destination_new =
          calculate_cuevi_walker{}(f, origin, current.beta, i_frac_destination,
                                   ib_destination, prior, lik, y, x);
      auto Maybe_origin_new =
          calculate_cuevi_walker{}(f, destination, current.beta, i_frac_origin,
                                   ib_origin, prior, lik, y, x);
      if (Maybe_origin_new.valid() && Maybe_destination_new.valid()) {
        auto [destination_new0, destination_new1] =
            std::move(Maybe_destination_new.value());
        auto [origin_new0, origin_new1] = std::move(Maybe_origin_new.value());
        double beta_origin = current.beta[i_frac_origin][ib_origin];
        double beta_destination =
            current.beta[i_frac_destination][ib_destination];

        auto logA =
            origin_new1.logP - origin.logP + destination_new0.logP -
            destination.logP + beta_origin * (origin_new1.logL - origin.logL) +
            beta_destination * (destination_new0.logL - destination.logL);
        auto pJump = std::min(1.0, std::exp(logA));
        if (pJump > r) {
          std::swap(origin, origin_new1);
          std::swap(destination, destination_new0);
          std::swap(current.i_walkers[iw][i_frac_origin][ib_origin],
                    current.i_walkers[jw][i_frac_destination][ib_destination]);
          if ((i_frac_origin > 0) && (beta_origin == 0)) {
            auto ib_origin0 = current.beta[i_frac_origin - 1].size() - 1;
            current.i_walkers[iw][i_frac_origin - 1][ib_origin0] =
                current.i_walkers[iw][i_frac_origin][ib_origin];
            current.walkers[iw][i_frac_origin - 1][ib_origin0] = origin_new0;
          }
          if ((i_frac_destination + 1 < current.beta.size()) &&
              beta_destination == 1.0) {
            auto ib_destination1 = 0ul;
            current.i_walkers[jw][i_frac_destination + 1][ib_destination1] =
                current.i_walkers[jw][i_frac_destination][ib_destination];
            current.walkers[jw][i_frac_destination + 1][ib_destination1] =
                destination_new1;
          }
        }
      }
    }
  }
};

struct thermo_cuevi_randomized_jump_mcmc {
  friend std::string ToString(thermo_cuevi_randomized_jump_mcmc) {
    return "thermo_cuevi_randomized_jump_mcmc";
  }

  template <class FunctionTable, class Observer, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, std::size_t iter,
                  cuevi_mcmc<Parameters> &current, Observer &, mt_64i &mt,
                  ensemble<mt_64i> &mts, Prior const &prior,
                  Likelihood const &lik, const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  std::size_t thermo_jumps_every) const {
    if (iter % (thermo_jumps_every) == 0) {
      // print_i_walkers(std::cerr, current);
      std::uniform_real_distribution<double> uniform_real(0, 1);
      auto n_walkers = mts.size() * 2;
      auto n_par = current.walkers[0][0][0].parameter.size();
      std::uniform_int_distribution<std::size_t> booldist(0, 1);

      WalkerIndexes shuffled_walker(n_walkers);
      std::iota(shuffled_walker.begin(), shuffled_walker.end(), 0);
      std::shuffle(shuffled_walker.begin(), shuffled_walker.end(), mt);
      std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                                uniform_real);

      auto n_fractions = current.beta.size();
      FractionIndexes landing_fraction(n_fractions);
      std::iota(landing_fraction.begin(), landing_fraction.end(), 0);
      std::shuffle(landing_fraction.begin(), landing_fraction.end(), mt);

#pragma omp parallel for // not currently working
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto iw = shuffled_walker[i];
        auto jw = shuffled_walker[i + n_walkers / 2];

        for (std::size_t i_fr = 0; i_fr < n_fractions; ++i_fr) {
          auto i_frac_origin = i_fr;
          auto i_frac_destination = landing_fraction[i_fr];
          auto r = rdist[i](mts[i]);

          thermo_cuevi_randomized_jump_mcmc_per_walker{}(
              f.fork(var::I_thread(i)), current, iw, jw, i_frac_origin,
              i_frac_destination, prior, lik, y, x, r);
        }
      }
      // print_i_walkers(std::cerr, current);
    }
  }
};

template <class Parameters>
auto derivative_var_ratio(by_fraction<by_beta<double>> const &mean,
                          by_fraction<by_beta<double>> const &var,
                          cuevi_mcmc<Parameters> const &current) {

  auto &beta = current.beta;
  by_fraction<by_beta<double>> out(size(mean));
  for (std::size_t i = 0; i < mean.size(); ++i) {
    out[i] = derivative_var_ratio_beta(mean[i], var[i], beta[i]);
  }
  return out;
}
}

template <>
inline bool compare_to_max_ratio(deprecated::by_fraction<by_beta<double>> const &beta,
                                 deprecated::by_fraction<by_beta<double>> const &mean_logL,
                                 deprecated::by_fraction<by_beta<double>> const &var_ratio,
                                 double max_ratio) {
  for (std::size_t i_frac = 0; i_frac < size(var_ratio); ++i_frac) {
    for (std::size_t ib = 0; ib < size(var_ratio[i_frac]); ++ib) {
      std::cerr << "(" << beta[i_frac][ib] << "[~" << mean_logL[i_frac][ib]
                << "]=> " << var_ratio[i_frac][ib] << ")  ";
      if ((var_ratio[i_frac][ib] > max_ratio) ||
          (var_ratio[i_frac][ib] < 1.0 / max_ratio)) {
        std::cerr << "  FALSE \n";
        return false;
      }
    }
  }
  std::cerr << " TRUE\n";
  return true;
}

namespace deprecated{
template <class Algorithm, class Fractioner, class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, cuevi_mcmc<Parameters>>)
class cuevi_integration {
  Algorithm alg_;
  Fractioner frac_;
  Reporter rep_;
  std::size_t num_scouts_per_ensemble_;
  std::size_t max_num_simultaneous_temperatures_;
  double min_fraction_;
  std::size_t thermo_jumps_every_;
  double n_points_per_decade_beta_;
  double n_points_per_decade_fraction_;
  double stops_at_;
  bool includes_zero_;
  std::size_t initseed_;

public:
  cuevi_integration(Algorithm &&alg, Fractioner &&frac, Reporter &&rep,
                    std::size_t num_scouts_per_ensemble,
                    std::size_t max_num_simultaneous_temperatures,
                    double min_fraction, std::size_t thermo_jumps_every,
                    double n_points_per_decade_beta,
                    double n_points_per_decade_fraction, double stops_at,
                    bool includes_zero, std::size_t initseed)
      : alg_{std::move(alg)}, frac_{std::move(frac)}, rep_{std::move(rep)},
        num_scouts_per_ensemble_{num_scouts_per_ensemble},
        max_num_simultaneous_temperatures_{max_num_simultaneous_temperatures},
        min_fraction_{min_fraction}, thermo_jumps_every_{thermo_jumps_every},
        n_points_per_decade_beta_{n_points_per_decade_beta},
        n_points_per_decade_fraction_{n_points_per_decade_fraction},
        stops_at_{stops_at}, includes_zero_{includes_zero},
        initseed_{initseed} {}

  auto &algorithm() const { return alg_; }
  auto &fractioner() const { return frac_; }
  auto &reporter() { return rep_; }
  auto &min_fraction() const { return min_fraction_; }
  auto &num_scouts_per_ensemble() const { return num_scouts_per_ensemble_; }
  auto &max_num_simultaneous_temperatures() const {
    return max_num_simultaneous_temperatures_;
  }
  auto &thermo_jumps_every() const { return thermo_jumps_every_; }
  auto &n_points_per_decade_beta() const { return n_points_per_decade_beta_; }
  auto &n_points_per_decade_fraction() const {
    return n_points_per_decade_fraction_;
  }
  auto &stops_at() const { return stops_at_; }
  auto &includes_zero() const { return includes_zero_; }
  auto &initseed() const { return initseed_; }
};

template <class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Fractioner, class Reporter,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_Algorithm_conditions<Algorithm, cuevi_mcmc<Parameters>> &&
           is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                                                                                                                                                   DataType>)

auto evidence(FunctionTable &&ff,
              cuevi_integration<Algorithm, Fractioner, Reporter> &&cue,
              Prior const &prior, Likelihood const &lik, const DataType &y,
              const Variables &x, bool all_at_once) {
  auto f = ff.fork(var::I_thread(0));
  auto a = cue.algorithm();
  auto mt = init_mt(cue.initseed());
  auto n_walkers = cue.num_scouts_per_ensemble();
  auto max_num_simultaneous_temperatures =
      cue.max_num_simultaneous_temperatures();
  auto mts = init_mts(mt, cue.num_scouts_per_ensemble() / 2);
  auto [ys, xs, beta_final] = cue.fractioner()(
      y, x, mt, size(prior) * cue.min_fraction(),
      cue.n_points_per_decade_beta(), cue.n_points_per_decade_fraction(),
      cue.stops_at(), cue.includes_zero());
  auto beta_init =
      by_beta<double>(beta_final[0].begin(), beta_final[0].begin() + 2);
  auto current = all_at_once ? init_cuevi_mcmc_all(f, n_walkers, beta_final,
                                                   mts, prior, lik, ys, xs)
                             : init_cuevi_mcmc(f, n_walkers, beta_init, mts,
                                               prior, lik, ys, xs);
  // print_i_walkers(std::cerr, current);

  auto mcmc_run = checks_convergence(std::move(a), current);
  if ((current.beta.size() == beta_final.size()) &&
      (current.beta.back().back() == 1.0))
    mcmc_run.first.we_reach_final_temperature();

  std::size_t iter = 0;
  auto &rep = cue.reporter();
  report_title(rep, current, prior, lik, ys, xs);
  report_model(rep, prior, lik, ys, xs, beta_final);
  report_title(ff, "Iter");

  // auto it_frac = beta_final.begin();
  // auto it_beta = it_frac->begin() + 2;
  while ((current.nsamples.back() < size(ys[size(ys) - 1])) ||
         (size(current.beta.back()) < size(beta_final.back())) ||
         !mcmc_run.second) {
    while (!mcmc_run.second) {
      // print_i_walkers(std::cerr, current);
      f.f(step_stretch_cuevi_mcmc{}, current, rep, mts, prior, lik, ys, xs);
      // print_i_walkers(std::cerr, current);
      // check_sanity(iter,current);
      report_point(ff, iter);

      ++iter;
      // print_i_walkers(std::cerr, current);
      f.f(thermo_cuevi_jump_mcmc{}, iter, current, rep, mt, mts, prior, lik, ys,
          xs, cue.thermo_jumps_every());
      // print_i_walkers(std::cerr, current);
      f.f(thermo_cuevi_randomized_jump_mcmc{},
          iter + cue.thermo_jumps_every() % 2, current, rep, mt, mts, prior,
          lik, ys, xs, cue.thermo_jumps_every());

      // check_sanity(iter,current);

      report_all(f, iter, rep, current, prior, lik, ys, xs);
      mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
    }
    if ((current.nsamples.back() < size(ys[size(ys) - 1])) ||
        (size(current.beta.back()) < size(beta_final.back()))) {
      //   std::cerr<<"\n---walkers!!------------------------\n"<<current.walkers;

      auto is_current = push_back_new_fraction(
          f, current, mts, beta_final, max_num_simultaneous_temperatures, prior,
          lik, ys, xs);
      while (!(is_current)) {
        std::cerr << is_current.error()();
        f.f(step_stretch_cuevi_mcmc{}, current, rep, mts, prior, lik, ys, xs);
        report_point(ff, iter);
        ++iter;
        f.f(thermo_cuevi_jump_mcmc{}, iter, current, rep, mt, mts, prior, lik,
            ys, xs, cue.thermo_jumps_every());
        report_all(f, iter, rep, current, prior, lik, ys, xs);
        is_current = push_back_new_fraction(f, current, mts, beta_final,
                                            max_num_simultaneous_temperatures,
                                            prior, lik, ys, xs);
      }
      current = std::move(is_current.value());
      std::cerr << "\niwalkers!!------------------------\n"
                << current.i_walkers;
      std::cerr << "\n  nsamples=" << current.nsamples.back()
                << "   beta_run=" << current.beta.back().back() << "\n";
      mcmc_run.first.reset();
      if ((current.beta.size() == beta_final.size()) &&
          (current.beta.back().back() == 1.0))
        mcmc_run.first.we_reach_final_temperature();

      mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
    }
  }

  return std::pair(mcmc_run, current);
}

template <class Parameters>
auto cuevi_by_convergence(std::string path, std::string filename,
                          std::size_t num_scouts_per_ensemble,
                          std::size_t max_num_simultaneous_temperatures,
                          double min_fraction, std::size_t thermo_jumps_every,
                          std::size_t max_iter, double max_ratio,
                          double n_points_per_decade_beta,
                          double n_points_per_decade_fraction, double stops_at,
                          bool includes_zero, std::size_t initseed) {
  return cuevi_integration(
      checks_derivative_var_ratio<cuevi_mcmc, Parameters>(max_iter, max_ratio),
      fractioner{},
      save_mcmc<Parameters, save_likelihood<Parameters>,
                save_Parameter<Parameters>, save_Evidence>(path, filename,
                                                           100ul, 100ul, 100ul),
      num_scouts_per_ensemble, max_num_simultaneous_temperatures, min_fraction,
      thermo_jumps_every, n_points_per_decade_beta,
      n_points_per_decade_fraction, stops_at, includes_zero, initseed);
}
} // namespace deprecated
#endif // CUEVI_H
