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
#include <random>
#include <vector>

namespace macrodr {

class mt_64 : public var::Var<mt_64, std::mt19937_64> {};

class N_channel_transition_count
    : public var::Var<N_channel_transition_count, Matrix<std::size_t>> {};

class N_channel_state_XTX_count
    : public var::Var<N_channel_state_XTX_count, SymPosDefMatrix<std::size_t>> {};

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
                               N_channel_transition_mean const &nr) {

  double N_mean = var::sum(nr());
  std::size_t N_count = sample_N(mt, N_mean);

  auto v_Nr = N_channel_transition_mean(nr() * (N_count / N_mean));
  Matrix<std::size_t> out =
      applyMap([](double N) { return std::size_t(std::floor(N)); }, v_Nr());

  Matrix<double> sum_r = v_Nr() - out;

  assert(N_count > var::sum(out));
  double Nr = N_count - var::sum(out);
  Matrix<double> r = sum_r / Nr;
  ;

  Matrix<std::size_t> rs = multinomial_distribution(Nr, r)(mt);

  return N_channel_transition_count(out + rs);
}

class Micror_parameters_distribution;

class Micror_state
    : public var::Var<Micror_state, Vector_Space<N_channel_transition_count,
                                                 N_channel_transition_mean>> {
  Micror_state(std::mt19937_64 &mt, N_channel_transition_mean const &mean)
      : base_type{Vector_Space(N_total_ij__mean_to_count(mt, mean), mean)} {}

public:
  using base_type =
      var::Var<Micror_state, Vector_Space<N_channel_transition_count,
                                          N_channel_transition_mean>>;
  using base_type::Var;
  auto size() const {
    return get<N_channel_transition_mean>((*this)())().size();
  }

  friend class Micror_parameters_distribution;

  friend auto stretch_move(std::mt19937_64 &mt,
                           std::uniform_real_distribution<double> &rdist,
                           const Micror_state &Xk, const Micror_state &Xj) {
    assert((Xj.size() == Xk.size()) &&
           "sum of vector fields of different sizes");

    auto Nk = get<N_channel_transition_mean>(Xk());
    auto Nj = get<N_channel_transition_mean>(Xj());

    N_channel_transition_mean out = Nj;
    double z;
    double is_valid = false;
    while (!is_valid) {
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

class Micror_parameters_distribution {
  N_Ch_mean_value m_N;
  dirchlet_distribution P_tra_dist;
  multinomial_distribution m_N_tra_dist;

public:
  Micror_state operator()(std::mt19937_64 &mt) {
    auto r_Nij_mean = N_channel_transition_mean(P_tra_dist(mt) * m_N());
    return Micror_state(mt, r_Nij_mean);
  }

  Micror_parameters_distribution(N_Ch_mean_value const &N,
                                 Ptotal_ij const &Ptotal)
      : m_N{N}, P_tra_dist{Ptotal() * N()}, m_N_tra_dist{N(), Ptotal()} {}
  Maybe_error<double> logP(const Micror_state &s) {
    return m_N_tra_dist.logP(get<N_channel_transition_count>(s())());
  }
};

class Micror_parameters_likelihood {
  template <class FuncTable>
  Maybe_error<double>
  logLikelihood(FuncTable &&f, Micror_parameters_likelihood,
                Micror_state const &p, const Patch_current &y,
                const Vector_Space<gmean_ij, gvar_ij, Current_Noise> &x) {
    auto v_Nij = get<N_channel_transition_count>(p());
    auto y_mean = var::sum(elemMult(get<gmean_ij>(x)(), v_Nij()));
    auto y_var = var::sum(elemMult(get<gvar_ij>(x)(), v_Nij()));
    auto r_var = get<Current_Noise>(x)() + y_var;
    auto chi2 = sqr(y() - y_mean) / r_var;

    return -0.5 * log(2 * std::numbers::pi * r_var) - 0.5 * chi2;
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

N_channel_state_count
calc_Ntotal_j_sum_iter(const thermo_mcmc<Micror_state> &data, std::size_t k,
                        std::size_t i_beta_1) {
  Matrix<std::size_t> u(1, k, 1ul);
  return fold(
      data.walkers,
      [i_beta_1, &u](std::vector<mcmc<Micror_state>> const &walker) {
        return N_channel_state_count(u * (get<N_channel_transition_count>(
                                             walker[i_beta_1].parameter())()));
      },
      [](N_channel_state_count &&one, N_channel_state_count const &two) {
        return N_channel_state_count(std::move(one()) + two());
      });
}

N_channel_state_count
calc_Ntotal_j_sum(std::size_t i_start, std::size_t i_end,
                   const std::vector<thermo_mcmc<Micror_state>> &data_chain,
                   std::size_t k, std::size_t i_beta_1) {

  return fold(
      i_start, i_end, data_chain,
      [k, i_beta_1](thermo_mcmc<Micror_state> const &x) {
        return calc_Ntotal_j_sum_iter(x, k, i_beta_1);
      },
      [](N_channel_state_count &&a, N_channel_state_count &&b) {
        return N_channel_state_count(std::move(a()) + b());
      });
}

N_channel_state_XTX_count
calc_Ntotal_j_XTX_iter(const thermo_mcmc<Micror_state> &data,
                            std::size_t k, std::size_t i_beta_1) {
  Matrix<std::size_t> u(1, k, 1ul);
  return fold(
      data.walkers,
      [i_beta_1, &u](std::vector<mcmc<Micror_state>> const &walker) {
          return N_channel_state_XTX_count(XTX(
                  u * (get<N_channel_transition_count>(
                          walker[i_beta_1].parameter())())));
      },
      [](N_channel_state_XTX_count &&one,
         N_channel_state_XTX_count const &two) {
        return N_channel_state_XTX_count(std::move(one()) + two());
      });
}

N_channel_state_XTX_count
calc_Ntotal_j_sum_XTX(std::size_t i_start, std::size_t i_end,
                       const std::vector<thermo_mcmc<Micror_state>> &data_chain,
                       std::size_t k, std::size_t i_beta_1) {

  return fold(
      i_start, i_end, data_chain,
      [k, i_beta_1](thermo_mcmc<Micror_state> const &x) {
        return calc_Ntotal_j_XTX_iter(x, k, i_beta_1);
      },
      [](N_channel_state_XTX_count &&a, N_channel_state_XTX_count &&b) {
        return N_channel_state_XTX_count(std::move(a()) + b());
      });
}

Vector_Space<P_mean,P_Cov> calculate_Pmean_cov(std::size_t i_start, std::size_t i_end,
                       const std::vector<thermo_mcmc<Micror_state>> &data_chain,
                       std::size_t k, std::size_t i_beta_1) {
    auto N_sum =calc_Ntotal_j_sum(i_start,i_end,data_chain,k,i_beta_1);
    auto N_sum_sqr =calc_Ntotal_j_sum_XTX(i_start,i_end,data_chain,k,i_beta_1);
    auto n=data_chain[i_start].walkers.size()*(i_end-i_start);
    auto N_mean=N_sum()*(1.0/n);
    auto N=var::sum(N_mean);
    auto v_P_mean=N_mean*(1.0/N);
    
    auto N_XTX_mean=N_sum_sqr()*1.0/n;
    auto v_P_cov=(N_XTX_mean-XTX(N_mean))*(1.0/N);
    
    return Vector_Space(P_mean(v_P_mean),P_Cov(v_P_cov));
}



auto calculate_Evidence_iter(const thermo_mcmc<Micror_state> &data) {
  auto meanLik = mean_logL(data);
  return calculate_Evidence(data.beta, meanLik);
}

auto calculate_Evidence_mean(std::size_t i_start, std::size_t i_end,
    const std::vector<thermo_mcmc<Micror_state>> &data_chain) {
  std::size_t n = i_end - i_start;

  return fold(
             i_start, i_end, data_chain,
             [](thermo_mcmc<Micror_state> const &data) {
               return calculate_Evidence_iter(data);
             },
             [](auto one, auto two) { return one + two; }) /
         n;
}

auto
get_Patch_State(std::size_t i_start, std::size_t i_end,const std::vector<thermo_mcmc<Micror_state>> &data_chain, std::size_t k, std::size_t i_beta_1) {
  //    Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean,
  //                 y_var, plogL, eplogL, vplogL>
    
    auto v_Evi=calculate_Evidence_mean(i_start,i_end,data_chain);
    auto v_P=calculate_Pmean_cov(i_start,i_end,data_chain,k,i_beta_1);
    return std::tuple(v_Evi,v_P);
}

template <uses_recursive_aproximation recursive,
          uses_averaging_aproximation averaging,
          uses_variance_aproximation variance, class FunctionTable,
          class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>
  requires(
      /*(U<std::decay_t<C_Patch_State>,
         Patch_State>||U<std::decay_t<C_Patch_State>,
         Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
      U<C_double, double> && U<C_Qdt, Qdt>)

Maybe_error<C_Patch_State>
Micror_stochastic(FunctionTable &ftbl, C_Patch_State &&t_prior,
                  C_Qdt const &t_Qdt, C_Patch_Model const &m,
                  C_double const &Nch, const Patch_current &p_y, double fs,
                  std::size_t myseed) {

  auto &p_P_cov = get<P_Cov>(t_prior);
  auto &p_P_mean = get<P_mean>(t_prior);

  std::size_t num_states = p_P_cov().nrows();
  auto num_par = num_states * (num_states + 1);
  std::string path = "";
  std::size_t num_scouts_per_ensemble =
      std::pow(2, std::ceil(std::log2(num_par)));

  auto max_num_simultaneous_temperatures = 100;
  auto thermo_jumps_every = 4ul;
  auto max_iter_warming = 100ul;
  auto max_iter_equilibrium = 10000ul;
  auto n_points_per_decade = 6ul;
  auto stops_at = 1e-2;
  bool includes_zero = true;

  auto tmi = thermo_store_every<Micror_state>(
      num_scouts_per_ensemble, max_num_simultaneous_temperatures,
      thermo_jumps_every, max_iter_warming, max_iter_equilibrium,
      n_points_per_decade, stops_at, includes_zero, myseed);

  auto r_Pij = Ptotal_ij(diag(get<P_mean>(t_prior)()) * get<P>(t_Qdt)());

  auto r_prior = Micror_parameters_distribution(N_Ch_mean_value(Nch), r_Pij);

  std::vector<thermo_mcmc<Micror_state>> chains =
      evidence(ftbl, std::move(tmi), r_prior, Micror_parameters_likelihood{},
               p_y,
               Vector_Space(get<gmean_ij>(t_Qdt), get<gvar_ij>(t_Qdt),
                            get<Current_Noise>(m)))
          .first;
}

} // namespace macrodr

#endif // MICROR_STOCHASTIC_H
