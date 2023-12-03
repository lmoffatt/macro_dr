#ifndef MICROR_STOCHASTIC_H
#define MICROR_STOCHASTIC_H

#include "experiment.h"
#include "matrix.h"
#include "maybe_error.h"
#include "multivariate_normal_distribution.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "qmodel.h"
#include "variables.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace macrodr {

class mt_64 : public var::Var<mt_64, std::mt19937_64> {};

class N_channel_transition_count
    : public var::Var<N_channel_transition_count, Matrix<std::size_t>> {};

class N_channel_state_XTX_count
    : public var::Var<N_channel_state_XTX_count, SymPosDefMatrix<std::size_t>> {
};

class N_channel_state_count
    : public var::Var<N_channel_state_count, Matrix<std::size_t>> {};

class N_channel_transition_mean
    : public var::Var<N_channel_transition_mean, Matrix<double>> {

  friend bool is_valid(N_channel_transition_mean const &x) {

    return (var::max(x()) > 0) && (var::min(x()) >= 0.0);
  }
};

static auto sample_Multinomial_state(std::mt19937_64 &mt,
                                     P_mean const &t_P_mean, std::size_t N) {
  auto k = t_P_mean().size();
  N_channel_state out(Matrix<double>(1, k));
  std::size_t N_remaining = N;
  double p_remaining = 1;
  for (std::size_t i = 0; i + 1 < k; ++i) {
    auto n = std::binomial_distribution<std::size_t>(
        N_remaining, t_P_mean()[i] / p_remaining)(mt);
    N_remaining -= n;
    p_remaining -= t_P_mean()[i];
    out()[i] = n;
  }
  out()[k - 1] = N_remaining;
  return out;
}

auto N_total_ij__mean_to_count(std::mt19937_64 &mt,
                               N_channel_transition_mean const &nr,
                               double max_eps=std::sqrt(std::numeric_limits<double>::epsilon())) {

  double N_mean = var::sum(nr());
  std::size_t N_count = sample_N(mt, N_mean);

  auto v_Nr = N_channel_transition_mean(nr() * (N_count / N_mean));
  Matrix<std::size_t> out =
      applyMap([max_eps](double N) {
      return std::size_t(std::floor(N+max_eps)); }, v_Nr());

  Matrix<double> sum_r = v_Nr() - out;
  
  assert(N_count >= var::sum(out));
  double Nr = N_count - var::sum(out);
  if (Nr==0)
      return N_channel_transition_count(out);
  else  
  {
      Matrix<double> r = sum_r / Nr;
  

  Matrix<std::size_t> rs = multinomial_distribution(Nr, r)(mt);

  return N_channel_transition_count(out + rs);
  }
}

template <uses_averaging_aproximation averaging>
class Micror_parameters_distribution;

template <uses_averaging_aproximation averaging>
class Micror_state

    : public var::Var<
          Micror_state<averaging>,
          Vector_Space<N_channel_transition_count, N_channel_transition_mean>> {
    
    

public:
    Micror_state(std::mt19937_64 &mt, N_channel_transition_mean const &mean)
        : base_type{Vector_Space(N_total_ij__mean_to_count(mt, mean), mean)} {}
    // Micror_state(const Micror_state&)=default;
    // Micror_state(Micror_state&&)=default;
    // Micror_state& operator=(Micror_state&&)=default;
    // Micror_state& operator=(Micror_state const&)=default;
    // Micror_state()=default;
  using base_type = var::Var<
      Micror_state<averaging>,
      Vector_Space<N_channel_transition_count, N_channel_transition_mean>>;
  using base_type::Var;
  auto size() const {
    return get<N_channel_transition_mean>((*this)())().size();
  }

  friend class Micror_parameters_distribution<averaging>;

  friend auto stretch_move(std::mt19937_64 &mt,
                           std::uniform_real_distribution<double> &rdist,
                           const Micror_state &Xk, const Micror_state &Xj) {
    assert((Xj.size() == Xk.size()) &&
           "sum of vector fields of different sizes");

    auto Nk = get<N_channel_transition_mean>(Xk());
    auto Nj = get<N_channel_transition_mean>(Xj());
    
    N_channel_transition_mean out; 
    double z;
    double is_valid = false;
    while (!is_valid) {
      out= Nj;
      z = std::pow(rdist(mt) + 1, 2) / 2.0;
      is_valid = true;
      for (std::size_t i = 0; i < Nj().size(); ++i) {
        if (Nj()[i] + z * (Nk()[i] - Nj()[i]) >= 0.0) {
          out()[i] += z * (Nk()[i] - Nj()[i]);
        } else {
          is_valid = false;
        }
      }
    }
    return std::tuple(Micror_state(mt, out), z);
  }
};

template <uses_averaging_aproximation averaging>
class Micror_parameters_distribution {
  N_Ch_mean_value m_N;
  dirchlet_distribution P_tra_dist;
  multinomial_distribution m_N_tra_dist;

public:
  Micror_parameters_distribution()=default;
  Micror_state<averaging> operator()(std::mt19937_64 &mt) {
    auto r_Nij_mean = N_channel_transition_mean(P_tra_dist(mt) * m_N());
    return Micror_state<averaging>(mt, r_Nij_mean);
  }
  
  template<class Ptotalij>
      requires (std::is_same_v<Ptotalij,Ptotal_ij>&&averaging.value==2)
  Micror_parameters_distribution(N_Ch_mean_value const &N,
                                 Ptotalij const &Ptotal)
      : m_N{N}, P_tra_dist{Ptotal() * N()}, m_N_tra_dist{N(), Ptotal()} {}
  
  template<class Ptotalj>
      requires (std::is_same_v<Ptotalj,P_mean>&&averaging.value<2)
  Micror_parameters_distribution(N_Ch_mean_value const &N,
                                 Ptotalj const &Ptotal)
      : m_N{N}, P_tra_dist{Ptotal() * N()}, m_N_tra_dist{N(), Ptotal()} {}
  
  Maybe_error<double> logP(const Micror_state<averaging> &s)const  {
    return m_N_tra_dist.logP(get<N_channel_transition_count>(s())());
  }
};

template <uses_averaging_aproximation averaging,
          uses_variance_aproximation variance>
class Micror_parameters_likelihood {
  template <class FuncTable, class Qdt>
  friend double
  logLikelihoodd(FuncTable &&f, Micror_parameters_likelihood,
                Micror_state<averaging> const &p, const Patch_current &y,
                  const Qdt &x) {
    if constexpr (averaging.value==2) {
      auto v_Nij = get<N_channel_transition_count>(p());
      auto y_mean = var::sum(elemMult(get<gmean_ij>(x)(), v_Nij()));
      auto r_var = get<Current_Noise>(x)();
      if constexpr (variance.value)
        r_var = r_var + var::sum(elemMult(get<gvar_ij>(x)(), v_Nij()));
      auto chi2 = sqr(y() - y_mean) / r_var;
      
      return -0.5 * log(2 * std::numbers::pi * r_var) - 0.5 * chi2;
    } else     if constexpr (averaging.value==1) {
      auto v_Nj = get<N_channel_transition_count>(p());
      auto y_mean = var::sum(elemMult(get<gmean_i>(x)(), v_Nj()));

      auto r_var = get<Current_Noise>(x)();
      if constexpr (variance.value)
        r_var = r_var + var::sum(elemMult(get<gvar_i>(x)(), v_Nj()));
      auto chi2 = sqr(y() - y_mean) / r_var;

      return -0.5 * log(2 * std::numbers::pi * r_var) - 0.5 * chi2;
    } else     /*if constexpr (averaging.value==0) */{
        auto v_Nj = get<N_channel_transition_count>(p());
        auto y_mean = var::sum(elemMult(get<g>(x)(), v_Nj()));
        
        auto r_var = get<Current_Noise>(x)();
        auto chi2 = sqr(y() - y_mean) / r_var;
        
        return -0.5 * log(2 * std::numbers::pi * r_var) - 0.5 * chi2;
    }
  }
  template <class FuncTable, class Qdt>
  friend Maybe_error<double>
  logLikelihood(FuncTable &&f, Micror_parameters_likelihood,
                 Micror_state<averaging> const &p, const Patch_current &y,
                 const Qdt &x) {
      auto out=logLikelihoodd(std::forward<FuncTable>(f),Micror_parameters_likelihood{},p,y,x);
      if (std::isfinite(out))
          return out;
      else
          return error_message("not finite");
    
  }
  
  template < class Qdt>
 friend
 Patch_current     simulate(std::mt19937_64 & mt,Micror_parameters_likelihood,Micror_state<averaging> const &p,  const Qdt &x)
  {
      double y_mean;
      double r_var;
     if constexpr (averaging.value==2) {
         auto v_Nij = get<N_channel_transition_count>(p());
          y_mean = var::sum(elemMult(get<gmean_ij>(x)(), v_Nij()));
         r_var = get<Current_Noise>(x)();
         if constexpr (variance.value)
             r_var = r_var + var::sum(elemMult(get<gvar_ij>(x)(), v_Nij()));
         
     } else     if constexpr (averaging.value==1) {
         auto v_Nj = get<N_channel_transition_count>(p());
          y_mean = var::sum(elemMult(get<gmean_i>(x)(), v_Nj()));
         
          r_var = get<Current_Noise>(x)();
         if constexpr (variance.value)
             r_var = r_var + var::sum(elemMult(get<gvar_i>(x)(), v_Nj()));
         } else     /*if constexpr (averaging.value==0) */{
         auto v_Nj = get<N_channel_transition_count>(p());
         y_mean = var::sum(elemMult(get<g>(x)(), v_Nj()));
         
         r_var = get<Current_Noise>(x)();
         
     }
     return Patch_current(std::normal_distribution<>{}(mt)*std::sqrt(r_var)+y_mean);
     
   }
  
};

// the idea is to build a thermodynamic integration procedure to estimate the
// Patch_State

// First we go from Macro -> Micro
//  then we have to define the parameters of the Micro distribution
// then the prior and likelihood of those parameters
// then we build the thermodynamic integration procedure
// we obtain the posterior distibution and the Evidence
// we convert the posterior and evidence into Patch_State

/// so, the parameters of the Micro distribution are defined in such a way
/// that can be used in the emcee algorithm, that is their domain has to be
/// affine, which means among other things unbounded Number of channels or
/// probabilities are not bounded so we need another formulation [I reflected
/// later on the subject and I think is bullshit, there is no problem with
/// hard boundaries in affince sampling, you just put a prior probability of
/// zero, rejecting jumps to prohibited values, so bounded problems are just
/// fine, no need to change the variables]
///
///
/// so, the parameters we have to fit are just the values of the matrix Nij.
///
/// the stetch algorithm being linear has no problem with obtaining new values
/// of Nij, the only thing is that the Nij will become real numbers after
/// multiplied by a real, so we need to convert them back to the discrete
///  Nij, by assigning the fractional channels with a multinomial
///  distribution. This assigment is done just once during the
/// stretch-step, so the calculations of logPrior and logLikelihood are done
/// just once per strech-step
///
/// so, the idea is to generate the Nij using one dirchlet distribution
/// (instead of multinomial) for each state.
///
///
///  Old no longer relevant discussion below
/// the easieast formulation is to just take the logarithm of each value and
/// substract the sum (or the mean) there is technical detail and that is that
/// we cannot take logarithms of zero, so we have to arbitrarily decide a
/// value for the logarithm of zero.
///
/// so, the complete parameterization is given a fixed number of channels,
/// logRi + logRij, the logarithm of the initial ratio of channels in each
/// state and the ratio of transitions.
///
/// as the prior is expressed in terms of multinomial distributions we have to
/// convert the continuous logRi-logRij to the discrete
///  Ni -Nij, by assigning the fractional channels with a multinomial
///  distribution. This assigment is done just once during the
/// stretch-step, so the calculations of logPrior and logLikelihood are done
/// just once per strech-step
///
/// This is an extreme case for thermodynamic integration since we have just
/// one data: per run.
///
/// a version for doing a Micror-stochastic algorithm for the complete series
/// of measurments could be built using a different parameterization  logRij +
/// logRjk  where we just change the state at a given time fixing the starting
/// point and the ending point. In this way we could do metropolis jump for
/// just one measurement, keeping the rest of the random parameters fixed.
/// This would allow to increase the dimensionality of the search enormously,
/// and a more fine grain use of the information. On the other hand if we wish
/// to use a jump that involves the entire set of measurements at the same
/// time, we should have a number of walkers at least as big as the number of
/// random-parameters we have to estimate, which is easy to become impossible.
///
///
///
///

template <uses_averaging_aproximation averaging>
N_channel_state_count calc_Ntotal_j_sum_iter(
    const thermo_mcmc<Micror_state<averaging>> &data,
    std::size_t k, std::size_t i_beta_1) {
    Matrix<std::size_t> u;
    if constexpr (averaging.value>1)
        u=Matrix<std::size_t>(1, k, 1ul);
    else
        u=Matrix<std::size_t>(1, 1, 1ul);
  return foldMap(
      data.walkers,
      [i_beta_1, &u](
          std::vector<mcmc<
              Micror_state<averaging>>> const &walker) {
        return N_channel_state_count(u * (get<N_channel_transition_count>(
                                             walker[i_beta_1].parameter())()));
      },
      [](N_channel_state_count &&one, N_channel_state_count const &two) {
        return N_channel_state_count(std::move(one()) + two());
      });
}

template <uses_averaging_aproximation averaging>
N_channel_state_count
calc_Ntotal_j_sum(std::size_t i_start, std::size_t i_end,
                  const std::vector<thermo_mcmc<Micror_state<averaging>>> &data_chain,
                  std::size_t k, std::size_t i_beta_1) {

  return foldMap(
      i_start, i_end, data_chain,
        [k, i_beta_1](thermo_mcmc<Micror_state<averaging>> const &x) {
        return calc_Ntotal_j_sum_iter(x, k, i_beta_1);
      },
      [](N_channel_state_count &&a, N_channel_state_count &&b) {
        return N_channel_state_count(std::move(a()) + b());
      });
}

template <uses_averaging_aproximation averaging>
N_channel_state_XTX_count
calc_Ntotal_j_XTX_iter(const thermo_mcmc<Micror_state<averaging>> &data, std::size_t k,
                       std::size_t i_beta_1) {
    Matrix<std::size_t> u;
    if constexpr (averaging.value>1)
        u=Matrix<std::size_t>(1, k, 1ul);
    else
        u=Matrix<std::size_t>(1, 1, 1ul);
        
  return foldMap(
      data.walkers,
      [i_beta_1, &u](std::vector<mcmc<Micror_state<averaging>>> const &walker) {
        return N_channel_state_XTX_count(XTX(
            u *
            (get<N_channel_transition_count>(walker[i_beta_1].parameter())())));
      },
      [](N_channel_state_XTX_count &&one,
         N_channel_state_XTX_count const &two) {
        return N_channel_state_XTX_count(std::move(one()) + two());
      });
}

template <uses_averaging_aproximation averaging>
N_channel_state_XTX_count
calc_Ntotal_j_sum_XTX(std::size_t i_start, std::size_t i_end,
                      const std::vector<thermo_mcmc<Micror_state<averaging>>> &data_chain,
                      std::size_t k, std::size_t i_beta_1) {

  return foldMap(
      i_start, i_end, data_chain,
      [k, i_beta_1](thermo_mcmc<Micror_state<averaging>> const &x) {
        return calc_Ntotal_j_XTX_iter(x, k, i_beta_1);
      },
      [](N_channel_state_XTX_count &&a, N_channel_state_XTX_count &&b) {
        return N_channel_state_XTX_count(std::move(a()) + b());
      });
}

template <uses_averaging_aproximation averaging>
Vector_Space<P_mean, P_Cov>
calculate_Pmean_cov(std::size_t i_start, std::size_t i_end,
                    const std::vector<thermo_mcmc<Micror_state<averaging>>> &data_chain,
                    std::size_t k) {
  auto beta = data_chain[0].beta;
  std::size_t i_beta_1 = 0;
  if (beta[i_beta_1] != 1.0)
    i_beta_1 = beta.size() - 1;

  auto N_sum = calc_Ntotal_j_sum(i_start, i_end, data_chain, k, i_beta_1);
  auto N_sum_sqr =
      calc_Ntotal_j_sum_XTX(i_start, i_end, data_chain, k, i_beta_1);
  auto n = data_chain[i_start].walkers.size() * (i_end - i_start);
  auto N_mean = N_sum() * (1.0 / n);
  auto N = var::sum(N_mean);
  auto v_P_mean = N_mean * (1.0 / N);

  auto N_XTX_mean = N_sum_sqr() * 1.0 / n;
  auto v_P_cov = (N_XTX_mean - XTX(N_mean)) * (1.0 / N);

  return Vector_Space(P_mean(v_P_mean), P_Cov(v_P_cov));
}

template <uses_averaging_aproximation averaging>
auto calculate_Evidence_iter(const thermo_mcmc<Micror_state<averaging>> &data) {
  auto meanLik = mean_logL(data);
  return calculate_Evidence(data.beta, meanLik);
}

template <uses_averaging_aproximation averaging>
auto calculate_Evidence_mean(
    std::size_t i_start, std::size_t i_end,
    const std::vector<thermo_mcmc<Micror_state<averaging>>> &data_chain) {
  std::size_t n = i_end - i_start;

  return plogL(foldMap(
                   i_start, i_end, data_chain,
                   [](thermo_mcmc<Micror_state<averaging>> const &data) {
                     return calculate_Evidence_iter(data);
                   },
                   [](auto one, auto two) { return one + two; }) /
               n);
}

template <uses_averaging_aproximation averaging>
auto get_Patch_State(std::size_t i_start, std::size_t i_end,
                     const std::vector<thermo_mcmc<Micror_state<averaging>>> &data_chain,
                     std::size_t k) {
  //    Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean,
  //                 y_var, plogL, eplogL, vplogL>

  auto v_Evi = calculate_Evidence_mean(i_start, i_end, data_chain);
  auto v_P = calculate_Pmean_cov(i_start, i_end, data_chain, k);
  return concatenate(Vector_Space(v_Evi), std::move(v_P));
}
 
template <uses_averaging_aproximation averaging,
          uses_variance_aproximation variance,
         class FunctionTable,
          class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>
  requires(
      /*(U<std::decay_t<C_Patch_State>,
         Patch_State>||U<std::decay_t<C_Patch_State>,
         Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
      U<C_double, double> && U<C_Qdt, Qdt>)

auto
Micror_stochastic(FunctionTable &ftbl, C_Patch_State t_prior,
                  C_Qdt const &t_Qdt, C_Patch_Model const &m,
                  C_double const &Nch, const Patch_current &p_y, double fs,
                  std::size_t myseed,std::size_t number_of_samples,std::vector<double> calculation_intervals_nodes,
                  std::size_t save_every_iter, std::size_t n_points_per_decade, double stops_at) {

  auto &p_P_cov = get<P_Cov>(t_prior);
  auto &p_P_mean = get<P_mean>(t_prior);

  std::size_t num_states = p_P_cov().nrows();
  auto num_par = num_states * (num_states + 1);
  std::string path = "";
  std::size_t num_scouts_per_ensemble =
      std::pow(2, std::ceil(std::log2(num_par)));
  
  auto max_num_simultaneous_temperatures = 10000;
  auto thermo_jumps_every = 4ul;
  auto max_iter_equilibrium = number_of_samples*save_every_iter;
  bool includes_zero = true;
  
  std::vector<std::pair<double, double>> calculation_intervals =
      MapAdj(calculation_intervals_nodes,[](auto one, auto two){return std::pair(one,two);});
  auto tmi = thermo_store_every<Micror_state<averaging>>(
      num_scouts_per_ensemble, max_num_simultaneous_temperatures,
      thermo_jumps_every, save_every_iter, max_iter_equilibrium,
      n_points_per_decade, stops_at, includes_zero, myseed);
  
  
  Micror_parameters_distribution<averaging> r_prior;
  if constexpr (averaging.value==2)
  {
      auto r_Pij = make_Ptotal_ij(diag(get<P_mean>(t_prior)()) * get<P>(t_Qdt)(),get<min_P>(m)());
      
      
      r_prior =
      Micror_parameters_distribution<averaging>(N_Ch_mean_value(Nch), r_Pij.value());
  }
  else
  {
      auto r_Pj = P_mean(get<P_mean>(t_prior)() * get<P>(t_Qdt)());
      r_prior =
          Micror_parameters_distribution<averaging>(N_Ch_mean_value(Nch), r_Pj);
      
      
  }
  auto k = p_P_cov().nrows();

  auto res =
      evidence(ftbl, std::move(tmi), r_prior,
               Micror_parameters_likelihood<averaging, variance>{}, p_y,
               Vector_Space(get<gmean_ij>(t_Qdt), get<gvar_ij>(t_Qdt),
                            get<Current_Noise>(m)));
  std::vector<thermo_mcmc<Micror_state<averaging>>> chains =std::move(res.first.chains());
  
  return Map(calculation_intervals,
             [&chains, k](std::pair<double, double> inter) {
               auto n = chains.size();
               std::size_t i0 = inter.first * n;
               std::size_t i1 = inter.second * n;
               return get_Patch_State(i0, i1, chains, k);
             });
}


template <uses_averaging_aproximation averaging,
         uses_variance_aproximation variance>
struct MicroR{
    friend std::string ToString(MicroR){
        std::string out="MicroR";
        if (averaging.value==2)
            out+="_2";
        else
            out+="__";
        if (variance.value)
            out+="_V";
        else
            out+="_M";
        
        return out; 
    }
};

} // namespace macrodr

#endif // MICROR_STOCHASTIC_H
