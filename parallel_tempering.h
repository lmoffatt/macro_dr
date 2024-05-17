#ifndef PARALLEL_TEMPERING_H
#define PARALLEL_TEMPERING_H
// #include "bayesian_linear_regression.h"
// #include "bayesian_linear_regression.h"
#include "function_measure_verification_and_optimization.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "multivariate_normal_distribution.h"
#include "parameters.h"
#include "random_samplers.h"
#include "variables.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

struct observer {
  observer() {}
  auto &operator[](std::size_t) const { return *this; }
};

template <class T> using ensemble = std::vector<T>;
template <class T> using by_beta = std::vector<T>;
template <class T> using by_iteration = std::vector<T>;

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

template <class Parameters>
auto stretch_move(mt_64i &, const Parameters &Xk, const Parameters &Xj,
                  double z) {
  assert((Xj.size() == Xk.size()) && "sum of vector fields of different sizes");
  auto out = Xj;
  for (std::size_t i = 0; i < Xj.size(); ++i)
    out[i] += z * (Xk[i] - Xj[i]);
  return out;
}

template <class Parameters>
auto stretch_move(mt_64i &mt, std::uniform_real_distribution<double> &rdist,
                  const Parameters &Xk, const Parameters &Xj, double a=2) {
  assert((Xj.size() == Xk.size()) && "sum of vector fields of different sizes");
    auto z = std::pow((a-1)*rdist(mt) + 1, 2) / a;

  auto out = Xj;
  for (std::size_t i = 0; i < Xj.size(); ++i)
    out[i] += z * (Xk[i] - Xj[i]);
  return std::tuple(out, z);
}

inline auto init_mts(mt_64i &mt, std::size_t n) {
  std::uniform_int_distribution<typename mt_64i::result_type> useed;
  std::vector<mt_64i> out;
  out.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    out.emplace_back(useed(mt));
  return out;
}

inline auto get_beta_list(double n_points_per_decade, double stops_at,
                          bool includes_zero) {
  std::size_t num_beta =
      std::ceil(-std::log10(stops_at) * n_points_per_decade) + 1;

  auto beta_size = num_beta;
  if (includes_zero)
    beta_size = beta_size + 1;

  auto out = std::vector<double>(beta_size, 0.0);
  for (std::size_t i = 0; i < num_beta; ++i)
    out[beta_size - 1 - i] = std::pow(10, -1.0 * i / n_points_per_decade);
  return out;
}

inline auto new_get_beta_list(std::size_t beta_size,
                              std::size_t beta_upper_size,
                              std::size_t beta_medium_size,
                              double beta_upper_value, double beta_medium_value,
                              double stops_at, bool includes_zero) {
  assert(beta_size > beta_upper_size);

  assert(beta_upper_size > beta_medium_size);
  assert(beta_upper_value > beta_medium_value && beta_upper_value < 1);
  assert(beta_medium_value > stops_at && stops_at > 0);

  std::size_t num_beta = includes_zero ? beta_size - 1 : beta_size;
  std::size_t beta_inferior_size =
      num_beta - beta_upper_size - beta_medium_size - 1;
  auto out = std::vector<double>(beta_size, 0.0);

  double n_points_per_decade_high =
      -(1.0 * beta_upper_size / std::log10(beta_upper_value));
  double n_points_per_decade_med =
      -(1.0 * beta_medium_size /
        std::log10(beta_medium_value / beta_upper_value));
  double n_points_per_decade_low =
      -(1.0 * beta_inferior_size / std::log10(stops_at / beta_medium_value));

  for (std::size_t i = 0; i < beta_upper_size; ++i)
    out[beta_size - 1 - i] =
        std::pow(10.0, -(1.0 * i) / n_points_per_decade_high);
  for (std::size_t i = 0; i < beta_medium_size; ++i)
    out[beta_size - 1 - i - beta_upper_size] =
        beta_upper_value * std::pow(10.0, -(1.0 * i) / n_points_per_decade_med);
  for (std::size_t i = 0; i < beta_inferior_size; ++i)
    out[beta_size - 1 - i - beta_upper_size - beta_medium_size] =
        beta_medium_value *
        std::pow(10.0, -(1.0 * i) / n_points_per_decade_low);
  out[beta_size - num_beta] = stops_at;
  assert(includes_zero ? (out[0] == 0) && (out[1] == stops_at)
                       : out[0] == stops_at);
  assert(out.back() == 1.0);
  assert(out[beta_size - beta_upper_size - 1] == beta_upper_value);
  assert(out[beta_size - beta_upper_size - beta_medium_size - 1] ==
         beta_medium_value);

  return out;
}

template <class Parameters>
  requires std::is_assignable_v<Parameters, Parameters const &>
struct thermo_mcmc {
  by_beta<double> beta;
  by_beta<ensemble<mcmc<Parameters>>> walkers;
  by_beta<ensemble<std::size_t>> i_walkers;
  by_beta<emcee_Step_statistics> emcee_stat;
  by_beta<Thermo_Jump_statistics> thermo_stat;
  
  
  void reset_statistics()
  {
      for (auto& e: emcee_stat)
          e().reset();
      for (auto& e: thermo_stat)
          e().reset();
      
  }
  
  auto get_Walkers_number() const { return walkers[0].size(); }
  auto &get_Beta() const { return beta; }
  auto &get_Parameter(std::size_t iw, std::size_t i_b) const {
    return walkers[i_b][iw].parameter;
  }

  auto get_Walker(std::size_t iw, std::size_t i_b) const {
    return i_walkers[i_b][iw];
  }
};

template <template <class> class Thermo_mcmc, class Parameters>
std::size_t num_betas(Thermo_mcmc<Parameters> const &x) {
  return x.beta.size();
}

template <class Parameters>
std::size_t num_Parameters(thermo_mcmc<Parameters> const &x) {
  return x.walkers[0][0].parameter.size();
}

template <class Parameters>
std::size_t num_walkers(thermo_mcmc<Parameters> const &x) {
  return x.walkers[0].size();
}

template <class Parameters>
std::size_t num_samples(by_iteration<thermo_mcmc<Parameters>> const &series) {
  return series.size();
}

template <class Parameters>
auto mean_logL(thermo_mcmc<Parameters> const &mcmc) {
    auto out = by_beta<logLs>(num_betas(mcmc), logLs());
  auto n_walkers = num_walkers(mcmc);
  for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
    for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
      out[ibeta] = out[ibeta]+ mcmc.walkers[ibeta][iwalker].logL / n_walkers;
  return out;
}
template <class Parameters>
auto mean_logP(thermo_mcmc<Parameters> const &mcmc) {
  auto out = by_beta<double>(num_betas(mcmc), 0);
  auto n_walkers = num_walkers(mcmc);
  for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
    for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
      out[ibeta] += mcmc.walkers[ibeta][iwalker].logP / n_walkers;
  return out;
}

template <class Parameters>
auto var_logL(thermo_mcmc<Parameters> const &mcmc,
              by_beta<logLs> const &mean) {
  auto out = by_beta<logLs>(num_betas(mcmc));
  auto n_walkers = num_walkers(mcmc);
  for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
    for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
      out[ibeta] =out[ibeta] +
          pow(mcmc.walkers[ibeta][iwalker].logL - mean[ibeta], 2) /
          n_walkers;
  return out;
}

template <class Parameters>
auto mean_logL(by_iteration<thermo_mcmc<Parameters>> const &series) {
  auto out = by_beta<double>(num_betas(series[0]), 0);
  auto n_walkers = num_walkers(series[0]);
  auto n_iters = num_samples(series);
  for (std::size_t i = 0; i < num_samples(series); ++i)
     for (std::size_t ibeta = 0; ibeta < num_betas(series[0]); ++ibeta)
      for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
        out[ibeta] +=
            series[i].walkers[ibeta][iwalker].logL / n_iters / n_walkers;
  return out;
}

template <class Parameters>
auto mean_logP(by_iteration<thermo_mcmc<Parameters>> const &series) {
  auto out = by_beta<double>(num_betas(series[0]), 0);
  auto n_walkers = num_walkers(series[0]);
  auto n_iters = num_samples(series);
  for (std::size_t i = 0; i < num_samples(series); ++i)
    for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
      for (std::size_t ibeta = 0; ibeta < num_betas(series[0]); ++ibeta)
        out[ibeta] +=
            series[i].walkers[ibeta][iwalker].logP / n_iters / n_walkers;
  return out;
}

template <class Parameters>
auto var_logL(by_iteration<thermo_mcmc<Parameters>> const &series,
              by_beta<logLs> const &mean) {
    auto out = by_beta<logLs>(num_betas(series[0]));
  auto n_walkers = num_walkers(series[0]);
  auto n_iters = num_samples(series);
  for (std::size_t i = 0; i < num_samples(series); ++i)
    for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
      for (std::size_t ibeta = 0; ibeta < num_betas(series[0]); ++ibeta)
        out[ibeta] +=
            pow(series[i].walkers[ibeta][iwalker].logL - mean[ibeta], 2) /
                        (n_iters * n_walkers);
  return out;
}

inline auto derivative_var_ratio_beta(by_beta<double> const &mean,
                                      by_beta<double> const &var,
                                      by_beta<double> const &beta) {
  by_beta<double> out(mean.size() - 1);
  for (std::size_t i = 0; i < mean.size() - 1; ++i) {
    auto dL = (mean[i + 1] - mean[i]) / (beta[i + 1] - beta[i]);
    auto dL0 = var[i + 1];
    auto dL1 = var[i];
    out[i] = (var[i + 1] + var[i]) / (2 * (mean[i + 1] - mean[i])) *
             (beta[i + 1] - beta[i]);
  }
  return out;
}

template <class Parameters>
auto derivative_var_ratio(by_beta<double> const &mean,
                          by_beta<double> const &var,
                          thermo_mcmc<Parameters> const &current) {
  return derivative_var_ratio_beta(mean, var, current.beta);
}

template <class Parameters>
auto mean_logL_walker(by_iteration<thermo_mcmc<Parameters>> const &series) {
  auto out = ensemble<by_beta<double>>(
      nun_walkers(series[0]), by_beta<double>(num_betas(series[0]), 0));
  for (std::size_t i = 0; i < num_samples(series); ++i)
    for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)

      for (std::size_t ibeta = 0; ibeta < num_beta(series[0]); ++ibeta)
        out[iwalker][ibeta] +=
            series[i].walkers[ibeta][iwalker].logL / num_samples(series);
  return out;
}
inline double calcEvidence(double b1, double b2, double L1, double L2) {
  return 0.5 * (L2 + L1) * (b2 - b1);
}
inline double calcEvidence_(double b1, double b2, double L1, double L2) {
  if (b1 == 0)
    return calcEvidence(b1, b2, L1, L2);
  auto db = b2 - b1;
  auto dL = L2 - L1;

  return L1 * db + dL / (std::log(b2) - std::log(b1)) *
                       (b2 * std::log(b2) - b1 * std::log(b1) - b2 + b1 -
                        std::log(b1) * db);
}

inline double calcEvidence(double b1, double b2, double L1, double L2,
                           double dL1, double dL2) {
  auto dL = dL2 - dL1;
  auto db = b2 - b1;
  return L1 * b2 + L2 * b1 + dL * b2 * b1 +
         0.5 * (b2 * b2 - b1 * b1) / db * (L2 - L1 + dL * (b1 + b2)) +
         1.0 / 6.0 * (std::pow(b2, 3) - std::pow(b1, 3)) * dL / db;
}

inline double calcEvidence(double b0, double L0, double dL0) {
  return L0 * b0 - 0.5 * dL0 * b0 * b0;
}
inline double calcEvidence(double b0, double L0) { return L0 * b0; }

inline double calculate_Evidence(by_beta<double> const &beta,
                                 by_beta<double> const &meanLik) {
  auto nb = beta.size();

  if (beta[0] > beta[nb - 1]) {
    double sum = calcEvidence(beta[nb - 1], meanLik[nb - 1]);
    for (std::size_t i = 1; i < beta.size(); ++i)
      sum += calcEvidence(beta[i], beta[i - 1], meanLik[i], meanLik[i - 1]);
    return sum;
  } else {
    double sum = calcEvidence(beta[0], meanLik[0]);
    for (std::size_t i = 1; i < beta.size(); ++i)
      sum += calcEvidence(beta[i - 1], beta[i], meanLik[i - 1], meanLik[i]);
    return sum;
  }
}

inline double calculate_Evidence(by_beta<double> const &beta,
                                 by_beta<double> const &meanLik,
                                 by_beta<double> const &varLik) {
  auto nb = beta.size();
  if (beta[0] > beta[nb - 1]) {
    double sum = calcEvidence(beta[nb - 1], meanLik[nb - 1], varLik[nb - 1]);
    for (std::size_t i = 1; i < beta.size(); ++i)
      sum += calcEvidence(beta[i], beta[i - 1], meanLik[i], meanLik[i - 1],
                          varLik[i], varLik[i - 1]);
    return sum;
  } else {
    double sum = calcEvidence(beta[0], meanLik[0], varLik[0]);
    for (std::size_t i = 1; i < beta.size(); ++i)
      sum += calcEvidence(beta[i - 1], beta[i], meanLik[i - 1], meanLik[i],
                          varLik[i - 1], varLik[i]);
    return sum;
  }
}

template <class Parameters>
auto var_logL_walker(by_iteration<thermo_mcmc<Parameters>> const &series,
                     by_beta<ensemble<double>> const &mean_logL_walker) {
  auto out = by_beta<ensemble<double>>(
      nun_walkers(series[0]), by_beta<double>(num_betas(series[0]), 0));
  for (std::size_t i = 0; i < num_samples(series); ++i)
    for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
      for (std::size_t ibeta = 0; ibeta < num_beta(series[0]); ++ibeta)
        out[iwalker][ibeta] += std::pow(series[i].walkers[ibeta][iwalker].logL -
                                        mean_logL_walker[ibeta][iwalker]) /
                               num_samples(series);
  return out;
}
template <class Parameters>
auto var_mean_logL_walker(by_beta<ensemble<double>> const &mean_logL_walker,
                          by_beta<double> const &mean) {
    auto n_beta = mean_logL_walker[0].size();
    auto n_walkers = mean_logL_walker.size();
  auto out = by_beta<double>(n_beta, 0.0);
  for (std::size_t i = 0; i < n_beta; ++i)
    for (std::size_t j = 0; j < n_walkers; ++j)
      out[i] += std::pow(mean_logL_walker[i][j] - mean[i], 2) / n_walkers;
  return out;
}

template <class Parameters>
auto mean_var_logL_walker(ensemble<by_beta<double>> const &var_logL_walker) {
  auto n_beta = var_logL_walker[0].size();
  auto n_walkers = var_logL_walker.size();
  auto out = by_beta<double>(n_beta, 0.0);
  for (std::size_t i = 0; i < n_beta; ++i)
    for (std::size_t j = 0; j < n_walkers; ++j)
      out[i] += var_logL_walker[j][i] / n_walkers;
  return out;
}

template <class Parameters>
auto mixing_var_ratio(by_beta<double> const &mean_var,
                      by_beta<double> const &var_mean) {
  by_beta<double> out(mean_var.size());
  for (std::size_t i = 0; i < mean_var.size(); ++i)
    out[i] = var_mean[i] / mean_var[i];
  return out;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(
      !is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
auto init_thermo_mcmc(FunctionTable &&f, std::size_t n_walkers,
                      by_beta<double> const &beta, ensemble<mt_64i> &mt,
                      Prior const &prior, Likelihood const &lik,
                      const DataType &y, const Variables &x) {

  ensemble<by_beta<std::size_t>> i_walker(n_walkers,
                                          by_beta<std::size_t>(beta.size()));
  ensemble<by_beta<mcmc<Parameters>>> walker(
      n_walkers, by_beta<mcmc<Parameters>>(beta.size()));
  
  by_beta<emcee_Step_statistics> emcee_stat(beta.size());
  by_beta<Thermo_Jump_statistics> thermo_stat(beta.size()-1);
  
  for (std::size_t half = 0; half < 2; ++half)
#pragma omp parallel for
    for (std::size_t iiw = 0; iiw < n_walkers / 2; ++iiw) {
      auto iw = iiw + half * n_walkers / 2;
      for (std::size_t i = 0; i < beta.size(); ++i) {
        i_walker[iw][i] = iw + i * n_walkers;
        walker[iw][i] =
            init_mcmc(f.fork(var::I_thread(iiw)), mt[iiw], prior, lik, y, x);
      }
    }
  return thermo_mcmc<Parameters>{beta, walker, i_walker,emcee_stat,thermo_stat};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables,
          class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(
      is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
auto init_thermo_mcmc(FunctionTable &&f, std::size_t n_walkers,
                      by_beta<double> const &beta, ensemble<mt_64i> &mt,
                      Prior const &prior, Likelihood const &lik,
                      const DataType &y, const Variables &x) {

  ensemble<by_beta<std::size_t>> i_walker(beta.size(),
                                          by_beta<std::size_t>(n_walkers));
  ensemble<by_beta<mcmc<Parameters>>> walker(beta.size()
      , by_beta<mcmc<Parameters>>(n_walkers));
  by_beta<emcee_Step_statistics> emcee_stat(beta.size());
  by_beta<Thermo_Jump_statistics> thermo_stat(beta.size()-1);
  auto ff = f.fork(omp_get_max_threads());

#pragma omp parallel for //collapse(2)
      for (std::size_t i = 0; i < beta.size(); ++i) {
    for (std::size_t iw = 0; iw < n_walkers ; ++iw) {
        i_walker[i][iw] = iw + i * n_walkers;
        auto i_th = omp_get_thread_num();
        walker[i][iw] = init_mcmc(ff[i_th], mt[i_th], prior, lik, y, x);
      }
    }
  f += ff;
  return thermo_mcmc<Parameters>{beta, walker, i_walker,emcee_stat,thermo_stat};
}

template <class Prior, 
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
auto create_thermo_mcmc(std::size_t n_walkers,
                      by_beta<double> const &beta,mt_64i& mt,Prior const &pr) {
    
    auto& priorsampler=pr;
    auto par = sample(mt,priorsampler);
    
    by_beta<ensemble<std::size_t>> i_walker(beta.size(),
                                            by_beta<std::size_t>(n_walkers));
    by_beta<ensemble<mcmc<Parameters>>> walker(beta.size()
        , by_beta<mcmc<Parameters>>(n_walkers,mcmc<Parameters>{par}));
    by_beta<emcee_Step_statistics> emcee_stat(beta.size());
    by_beta<Thermo_Jump_statistics> thermo_stat(beta.size()-1);
    return thermo_mcmc<Parameters>{beta, walker, i_walker,emcee_stat,thermo_stat};
}


template <class FunctionTable, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
Maybe_error<bool> calc_thermo_mcmc_continuation(FunctionTable &f, std::size_t n_walkers,
                      by_beta<double> const &beta, ensemble<mt_64i> &mt,
                      Prior const &prior, Likelihood const &lik,
                      const DataType &y, const Variables &x, thermo_mcmc<Parameters>& t_mcmc) {
    
    auto ff = f.fork(omp_get_max_threads());
    
    std::string error;
    bool good=true;
  #pragma omp parallel for //collapse(2)
      for (std::size_t i = 0; i < beta.size(); ++i) {
        auto i_th = omp_get_thread_num();
        for (std::size_t iw = 0; iw < n_walkers ; ++iw) {
                auto res=calc_mcmc(ff[i_th],prior, lik, y, x,t_mcmc.walkers[i][iw]);
                if (!res)
                {
                    error+=res.error()();
                    good=false;
                }
            }
        }
    f += ff;
        if (good)
        return good;
    else
            return error_message(error);        
 }
 
 
 
 
 
 
 
template <class Parameters>
std::pair<std::pair<std::size_t, std::size_t>, bool>
check_iterations(std::pair<std::size_t, std::size_t> current_max,
                 const thermo_mcmc<Parameters> &) {
  if (current_max.first >= current_max.second)
    return std::pair(std::pair(0ul, current_max.second), true);
  else
    return std::pair(std::pair(current_max.first + 1, current_max.second),
                     false);
};
struct no_save {

  template <class FunctionTable, class Parameters, class... T>
  friend void report(FunctionTable &&, std::size_t, no_save &,
                     thermo_mcmc<Parameters> const &, T &&...) {}
  template <class Prior, class Likelihood, class Variables, class DataType>
  friend void report_model(no_save &, Prior const &p, Likelihood const &,
                           const DataType &, const Variables &,
                           by_beta<double> const &) {}

  template <class Parameters>
  friend void report_title(no_save &, thermo_mcmc<Parameters> const &, ...) {}
};

template <class Algorithm, class Thermo_mcmc>
concept is_Algorithm_conditions = requires(Algorithm &&a) {
  {
    checks_convergence(std::move(a), std::declval<const Thermo_mcmc &>())
  } -> std::convertible_to<std::pair<Algorithm, bool>>;
};

class less_than_max_iteration {
  std::size_t current_iteration_;
  std::size_t max_iter_warming_;
  std::size_t max_iter_final_;
  std::size_t num_betas_;
  bool is_final_;

public:
  less_than_max_iteration(std::size_t max_iter_warming,
                          std::size_t max_iter_final)
      : current_iteration_{0ul}, max_iter_warming_{max_iter_warming},
        max_iter_final_{max_iter_final}, num_betas_{1ul}, is_final_{false} {}

  less_than_max_iteration &operator++() {
    ++current_iteration_;
    return *this;
  }
  std::size_t current_iteration() const { return current_iteration_; }
  std::size_t max_iteration_warming() const { return max_iter_warming_; }
  std::size_t max_iteration_final() const { return max_iter_final_; }

  std::size_t betas_number() const { return num_betas_; }

  bool is_final() const { return is_final_; }

  bool stop() const {
    if (is_final())
      return current_iteration() >= max_iteration_final();
    else
      return current_iteration() >= max_iteration_warming();
  }
  template <class Anything>
  friend auto checks_convergence(less_than_max_iteration &&c,
                                 const Anything &mcmc) {

    if (c.betas_number() < num_betas(mcmc)) {
      c.current_iteration_ = 0;
      c.num_betas_ = num_betas(mcmc);
    }
    if (c.stop()) {
      return std::pair(std::move(c), true);

    } else {
      ++c;
      return std::pair(std::move(c), false);
    }
  }

  void reset() { current_iteration_ = 0; }

  void we_reach_final_temperature() { is_final_ = true; }
};
// static_assert(is_Algorithm_conditions<less_than_max_iteration, thermo_mcmc>);

template <class mcmc>
  requires std::is_assignable_v<mcmc, mcmc const &>
class store_every_n_iter {
  std::size_t m_save_every;
  std::size_t m_max;
  std::size_t m_current_iter = 0ul;
  std::vector<mcmc> m_values;

public:
  store_every_n_iter(std::size_t save_every, std::size_t max)
      : m_save_every{save_every}, m_max{max} {
    m_values.reserve(max);
  }

  auto &chains() { return m_values; }

  friend auto checks_convergence(store_every_n_iter &&c, const mcmc &t_mcmc) {
    if (c.m_current_iter % c.m_save_every == 0)
      c.m_values.push_back(t_mcmc);
    if (c.m_current_iter < c.m_max) {
      ++c.m_current_iter;
      return std::pair<store_every_n_iter, bool>(std::move(c), false);
    } else {
      return std::pair<store_every_n_iter, bool>(std::move(c), true);
    }
  }
  void reset() {}
};

template <class Beta, class Var_ratio>
bool compare_to_max_ratio(Beta const &beta, Var_ratio const &mean_logL,
                          Var_ratio const &var_ratio, double max_ratio);

inline bool compare_to_max_ratio(by_beta<double> const &beta,
                                 by_beta<double> const &mean_logL,
                                 by_beta<double> const &var_ratio,
                                 double max_ratio) {
  for (std::size_t i = 0; i < var_ratio.size(); ++i) {
    std::cerr << "(" << beta[i] << "[~" << mean_logL[i] << "]=> "
              << var_ratio[i] << ")  ";
    if (var_ratio[i] > max_ratio) {
      std::cerr << "  FALSE \n";
      return false;
    }
  }
  std::cerr << " TRUE\n";
  return true;
}

template <template <class> class Thermo_mcmc, class Parameters>
class checks_derivative_var_ratio {
  std::size_t current_iteration_;
  double max_ratio_;
  by_iteration<Thermo_mcmc<Parameters>> curr_samples_;

public:
  checks_derivative_var_ratio(std::size_t sample_size, double max_ratio = 2)
      : current_iteration_(0ul), max_ratio_{max_ratio},
        curr_samples_{sample_size} {}

  checks_derivative_var_ratio &add(Thermo_mcmc<Parameters> const &x) {
    curr_samples_[current_iteration_ % curr_samples_.size()] = x;
    ++current_iteration_;
    return *this;
  }
  auto &current_samples() const { return curr_samples_; }

  auto get_derivative_var_ratio() const {
    auto m = mean_logL(curr_samples_);
    auto var = var_logL(current_samples(), m);
    return std::tuple(m, derivative_var_ratio(m, var, curr_samples_[0]));
  }

  bool converges() const {
    if (current_iteration_ % current_samples().size() == 0) {
      auto [meanL, var_ratio] = get_derivative_var_ratio();
      return compare_to_max_ratio(curr_samples_[0].beta, meanL, var_ratio,
                                  max_ratio_);
    } else {
      return false;
    }
  }
  friend auto checks_convergence(checks_derivative_var_ratio &&c,
                                 const Thermo_mcmc<Parameters> &mcmc) {
    if ((c.current_iteration_ > 0) &&
        (num_betas(c.current_samples()[0]) < num_betas(mcmc))) {
      c.current_iteration_ = 0;
    }
    c.add(mcmc);
    if (c.converges()) {
      return std::pair(std::move(c), true);
    } else {
      return std::pair(std::move(c), false);
    }
  }

  void reset() {}

  void we_reach_final_temperature() {}
};
// static_assert(is_Algorithm_conditions<checks_derivative_var_ratio<thermo_mcmc>,
//                                       thermo_mcmc>);

template <class Observer, class Parameters>
void observe_step_stretch_thermo_mcmc(
    Observer &obs, std::size_t j, double z, double r, const Parameters &current,
    const Parameters &candidate, double currlogP,
    Maybe_error<double> const calogP, double culogL,
    Maybe_error<double> const &calogL, bool doesChange) {}

template <class Observer, class Parameters>
void observe_thermo_jump_mcmc(Observer &obs, std::size_t jlanding,
                              const Parameters &current,
                              const Parameters &candidate, double culogL,
                              double calogL, double deltabeta, double logA,
                              double pJump, double r, bool doesChange) {}

template <class FunctionTable, class Observer, class Prior, class Likelihood,
          class Variables, class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(!is_of_this_template_type_v<std::decay_t<FunctionTable>,
                                       var::FuncMap_St> &&
           is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
void step_stretch_thermo_mcmc(FunctionTable &&f, std::size_t &iter,
                              thermo_mcmc<Parameters> &current, Observer &obs,
                              const by_beta<double> &beta, ensemble<mt_64i> &mt,
                              Prior const &prior, Likelihood const &lik,
                              const DataType &y, const Variables &x,
                              double alpha_stretch = 2) {
  assert(beta.size() == num_betas(current));
  auto n_walkers = num_walkers(current);
  auto n_beta = beta.size();
  auto n_par = current.walkers[0][0].parameter.size();

  std::uniform_int_distribution<std::size_t> uniform_walker(0,
                                                            n_walkers / 2 - 1);
  std::vector<std::uniform_int_distribution<std::size_t>> udist(n_walkers,
                                                                uniform_walker);

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
      for (std::size_t ib = 0; ib < n_beta; ++ib) {
        // we can try in the outer loop

        auto r = rdist[i](mt[i]);

        // candidate[ib].walkers[iw].
        auto [ca_par, z] =
            stretch_move(mt[i], rdist[i], current.walkers[ib][iw].parameter,
                         current.walkers[ib][jw].parameter);

        auto ca_logP = logPrior(prior, ca_par);
        auto ca_logL =
            logLikelihood(f.fork(var::I_thread(i)), lik, ca_par, y, x);

        if ((ca_logP.valid()) && (ca_logL.valid())) {
            auto dthLogL =std::isfinite(get<logL>(current.walkers[ib][iw].logL)())?
              ca_logP.value() - current.walkers[ib][iw].logP +
                           beta[ib] * (get<logL>(ca_logL.value())() - get<logL>(current.walkers[ib][iw].logL))():1.0;
          auto pJump =
              std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
          if constexpr (!std::is_same_v<Observer, no_save>)
            observe_step_stretch_thermo_mcmc(
                obs[iw][ib], jw, z, r, current.walkers[ib][iw].parameter,
                current.walkers[ib][jw].parameter, current.walkers[ib][iw].logP,
                ca_logP, current.walkers[ib][iw].logL, ca_logL, pJump >= r);
          if (pJump >= r) {
            current.walkers[ib][iw].parameter = std::move(ca_par);
            current.walkers[ib][iw].logP = ca_logP.value();
            current.walkers[ib][iw].logL = ca_logL.value();
          }
        }
      }
    }
  ++iter;
}

template <class FunctionTable, std::size_t N,class Observer, class Prior, class Likelihood,
          class Variables, class DataType,
          class Parameters = std::decay_t<decltype(sample(
              std::declval<mt_64i &>(), std::declval<Prior &>()))>>
  requires(is_of_this_template_type_v<std::decay_t<FunctionTable>,
                                      var::FuncMap_St> &&
           is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
void step_stretch_thermo_mcmc(FunctionTable &&f, std::size_t &iter,var::Event_Timing<N>& dur,
                              thermo_mcmc<Parameters> &current, Observer &obs,
                              const by_beta<double> &beta, ensemble<mt_64i> &mt,
                              Prior const &prior, Likelihood const &lik,
                              const DataType &y, const Variables &x,
                              double alpha_stretch = 2) {
  dur.record("stretch_start");  
  assert(beta.size() == num_betas(current));
  auto n_walkers = num_walkers(current);
  auto n_beta = beta.size();
  auto n_par = current.walkers[0][0].parameter.size();

  std::uniform_int_distribution<std::size_t> uniform_walker(0,
                                                            n_walkers / 2 - 1);
  std::vector<std::uniform_int_distribution<std::size_t>> udist(omp_get_max_threads(),
                                                                uniform_walker);

  std::uniform_real_distribution<double> uniform_stretch_zdist(
      1.0 / alpha_stretch, alpha_stretch);
  std::vector<std::uniform_real_distribution<double>> zdist(
      omp_get_max_threads(), uniform_stretch_zdist);

  std::uniform_real_distribution<double> uniform_real(0, 1);
  std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(),
                                                            uniform_real);

  auto ff = f.fork(omp_get_max_threads());
  std::vector<by_beta<emcee_Step_statistics>> emcee_stat(omp_get_max_threads(), by_beta<emcee_Step_statistics>(n_beta));
  
  dur.record("stretch_before_loop");
  auto num_threads=omp_get_max_threads();
  #pragma omp parallel for //collapse(2)
  for (std::size_t i_thread=0; i_thread<num_threads; ++i_thread)
  {
      std::size_t n_beta_f=std::ceil(n_beta/num_threads);
      for (std::size_t iib = 0; iib < n_beta_f; ++iib) {
      std::size_t ib=i_thread*n_beta_f+iib;    
      dur.record("begin_loop_walker", ib*2);
    //  for (std::size_t iii=0; iii<4; ++iii)
    //      {
      for (bool half : {false, true}){
          for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto i_th = omp_get_thread_num();

        auto iw = half ? i + n_walkers / 2 : i;
        auto j = udist[i_th](mt[i_th]);
        auto jw = half ? j : j + n_walkers / 2;
        // we can try in the outer loop

        auto r = rdist[i_th](mt[i_th]);

        // candidate[ib].walkers[iw].
        auto [ca_par, z] = stretch_move(mt[i_th], rdist[i_th],
                                        current.walkers[ib][iw].parameter,
                                        current.walkers[ib][jw].parameter);

        auto ca_logP = logPrior(prior, ca_par);
        auto ca_logL = logLikelihood(ff[i_th], lik, ca_par.to_value(), y, x);
        
        if (!((ca_logP.valid()) && (ca_logL.valid()))) {
            fails(emcee_stat[i_th][ib]());
        }
        else{
          auto dthLogL =
              ca_logP.value() - current.walkers[ib][iw].logP +
                           beta[ib] * (get<logL>(ca_logL.value())() - get<logL>(current.walkers[ib][iw].logL)());
          auto pJump =
              std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
          if (pJump >= r) {
            current.walkers[ib][iw].parameter = std::move(ca_par);
            current.walkers[ib][iw].logP = ca_logP.value();
            current.walkers[ib][iw].logL = ca_logL.value();
            succeeds(emcee_stat[i_th][ib]());
          }
          else
          {
              fails(emcee_stat[i_th][ib]());
          }
        }
      //}
          
      }
      }
      dur.record("end_loop_walker", ib*2+1);
      }
  }
      dur.advance(n_beta*2);
    dur.record("stretch_after_loop");  
for (std::size_t ib = 0; ib < n_beta; ++ib)
  {
      for(auto& e:emcee_stat)
        current.emcee_stat[ib]()+=e[ib]();
  }
  f += ff;
  ++iter;
  dur.record("stretch_function_end");  
  
}

inline double calc_logA(double betai, double betaj, logLs const& logLi,
                        logLs const &logLj) {
    if (std::isfinite((get<logL>(logLi)() - get<logL>(logLj)())))
         return -(betai - betaj) * (get<logL>(logLi)() - get<logL>(logLj)());
    else{
        if (!std::isfinite(get<logL>(logLi)())&&std::isfinite(get<logL>(logLj)()))
            return 20.0*(betai>betaj);
        if (std::isfinite(get<logL>(logLi)())&&!std::isfinite(get<logL>(logLj)()))
            return -20.0*(betai>betaj);
        else
            return -20.0;
    }
}

template <class Parameters, class Observer>
void thermo_jump_mcmc(std::size_t iter, thermo_mcmc<Parameters> &current,
                      Observer &obs, const by_beta<double> &beta, mt_64i &mt,
                      ensemble<mt_64i> &mts, std::size_t thermo_jumps_every) {
  if (iter % (thermo_jumps_every) == 0) {
    std::uniform_real_distribution<double> uniform_real(0, 1);
    auto n_walkers = current.get_Walkers_number();
    auto n_beta = beta.size();
    auto n_par = current.walkers[0][0].parameter.size();

    WalkerIndexes shuffeld_walkers(n_walkers);
    std::iota(shuffeld_walkers.begin(), shuffeld_walkers.end(), 0);
    std::shuffle(shuffeld_walkers.begin(), shuffeld_walkers.end(), mt);
    std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(),
                                                              uniform_real);
    
    std::vector<by_beta<Thermo_Jump_statistics>> thermo_stat(omp_get_max_threads(), by_beta<Thermo_Jump_statistics>(n_beta-1));
   

#pragma omp parallel for //collapse(2)
    for (std::size_t i = 0; i < n_walkers / 2; ++i) {
      for (std::size_t ib = 0; ib < n_beta - 1; ++ib) {
            auto i_th = omp_get_thread_num();
            auto iw = shuffeld_walkers[i];
            auto jw = shuffeld_walkers[i + n_walkers / 2];

        auto r = rdist[i_th](mts[i_th]);
        double logA =
            calc_logA(beta[ib], beta[ib + 1], current.walkers[ib][iw].logL,
                                current.walkers[ib + 1][jw].logL);
        auto pJump = std::min(1.0, std::exp(logA));
        if (pJump > r) {
          std::swap(current.walkers[ib][iw], current.walkers[ib + 1][jw]);
          std::swap(current.i_walkers[ib][iw], current.i_walkers[ib + 1][jw]);
          succeeds(thermo_stat[i_th][ib]());
        }
        else
        {
            fails(thermo_stat[i_th][ib]());
        }
      }
    }
    
    
    for (std::size_t ib = 0; ib < n_beta-1; ++ib)
    {
        for(auto& e:thermo_stat)
            current.thermo_stat[ib]()+=e[ib]();
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
auto push_back_new_beta(FunctionTable &&f, std::size_t &iter,
                        thermo_mcmc<Parameters> &current, ensemble<mt_64i> &mts,
                        by_beta<double> const &new_beta, Prior const &prior,
                        Likelihood const &lik, const DataType &y,
                        const Variables &x) {
  auto n_walkers = current.walkers[0].size();
  auto n_beta_old = current.walkers.size();
  current.walkers.push_back(ensemble<mcmc<Parameters>>{});
  current.i_walkers.push_back(ensemble<std::size_t>{});
  
  for (std::size_t half = 0; half < 2; ++half)
    for (std::size_t i = 0; i < n_walkers / 2; ++i) {
      auto iw = i + half * n_walkers / 2;
        current.walkers.back().push_back(init_mcmc(f, mts[i], prior, lik, y, x));
      current.i_walkers.back().push_back(n_beta_old * n_walkers + iw);
    }
  iter = 0;
  current.beta = new_beta;
  return current;
}

template <class Parameters> class save_likelihood {

public:
  std::string sep = ",";
  std::ofstream f;
  std::size_t save_every = 1;
  save_likelihood(std::string const &path, std::size_t interval)
      : f{std::ofstream(path + "__i_beta__i_walker.csv")},
        save_every{interval} {
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  }

  friend void report_title(save_likelihood &s, thermo_mcmc<Parameters> const &,
                           ...) {

    s.f << "iter" << s.sep << "iter_time" << s.sep
        << "beta" << s.sep << "i_walker" << s.sep << "id_walker" << s.sep
        << "logP" << s.sep << "logLik"<< s.sep << "elogLik" << s.sep
          << "vlogLik"<< s.sep << "plog_Evidence"<< s.sep
          << "pelog_Evidence"<< s.sep << "pvlog_Evidence"
          << s.sep << "log_Evidence"<< s.sep << "elog_Evidence"<< s.sep << "vlog_Evidence"
        << "\n";
  }
  template <class Prior, class Likelihood, class Variables, class DataType>
  friend void report_model(save_likelihood &, const Prior &, const Likelihood &,
                           const DataType &, const Variables &,
                           by_beta<double> const &) {}

  friend void report_model(save_likelihood &, ...) {}

  template <class FunctionTable, class Duration>
  friend void report(FunctionTable &&, std::size_t iter, const Duration &dur,
                     save_likelihood &s, thermo_mcmc<Parameters> const &data,
                     ...) {
    if (iter % s.save_every == 0) {
      for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker) {
            logLs t_logL = {};
        double beta = 0;
        logLs log_Evidence = var::Vector_Space<logL,elogL,vlogL>(logL(0.0),elogL(0.0),vlogL(0.0));  
        for (std::size_t i_beta = num_betas(data); i_beta > 0; --i_beta) {
          auto logL0 = t_logL;
          double beta0 = beta;
          t_logL = data.walkers[i_beta - 1][i_walker].logL;
          beta = data.beta[i_beta - 1];
          auto plog_Evidence = (beta - beta0) * (logL0 + t_logL) / 2.0;
          log_Evidence =log_Evidence+ plog_Evidence;
          s.f  << iter << s.sep << dur << s.sep
              << beta << s.sep << i_walker << s.sep
              << data.i_walkers[i_beta - 1][i_walker] << s.sep
              << data.walkers[i_beta - 1][i_walker].logP << t_logL.sep(s.sep)
              << plog_Evidence.sep(s.sep) << log_Evidence.sep(s.sep)
              <<"\n";
        }
      }
    }
  }
};
template <class Cova>
  requires Covariance<double, Cova>
Maybe_error<by_beta<double>> bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_gamma_normal_distribution<double, Cova> &prior,
    const Matrix<double> &y, const Matrix<double> &X,
    by_beta<double> const &beta0) {

  by_beta<double> out(size(beta0));
  for (std::size_t i = 0; i < size(out); ++i) {
    auto meanLogLiki =
        bayesian_linear_regression_calculate_mean_logLik(prior, y, X, beta0[i]);
    if (!meanLogLiki)
      return "bayesian_linear_regression_calculate_mean_logLik error for beta "
             "=" +
             std::to_string(beta0[i]) + ":  " + meanLogLiki.error();
    else
      out[i] = meanLogLiki.value();
  }
  return out;
}

template <class Parameters> class save_Parameter {

public:
  class separator : public std::string {
  public:
    using std::string::string;
    //   separator(std::string s):std::string(std::move(s)){}

    std::string operator()() const { return *this; }
    friend std::ostream &operator<<(std::ostream &os, const separator &sep) {
      return os << sep();
    }
    friend std::istream &operator>>(std::istream &is, const separator &sep) {
      std::string ss = sep();
      for (std::size_t i = 0; i < ss.size(); ++i) {
        is.get(ss[i]);
      }
      if (ss != sep())
        is.setstate(std::ios::failbit);
      return is;
    }
  };

  separator sep = ",";
  std::string fname;
  std::ofstream f;
  std::size_t save_every;
  save_Parameter(std::string const &path, std::size_t interval)
      : fname{path}, f{std::ofstream(path + "__i_beta__i_walker__i_par.csv")},
        save_every{interval} {
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  }

  friend void report_title(save_Parameter &s, thermo_mcmc<Parameters> const &,
                           ...) {

    s.f << "iter" << s.sep << "iter_time" << s.sep
        << "beta" << s.sep << "i_walker" << s.sep << "id_walker" << s.sep
        << "i_par" << s.sep << "par_value"
        << "\n";
  }
  template <class Prior, class Likelihood, class Variables, class DataType>
  friend void report_model(save_Parameter &, const Prior &, const Likelihood &,
                           const DataType &, const Variables &,
                           by_beta<double> const &) {}

  friend void report_model(save_Parameter &, ...) {}

  template <class FunctionTable, class Duration>
  friend void report(FunctionTable &&f, std::size_t iter, const Duration &dur,
                     save_Parameter &s, thermo_mcmc<Parameters> const &data,
                     ...) {
    if (iter % s.save_every == 0)
      for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
        for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker)
          for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)

            s.f << iter << s.sep << dur << s.sep
                << data.beta[i_beta] << s.sep << i_walker << s.sep
                << data.i_walkers[i_beta][i_walker] << s.sep << i_par << s.sep
                << data.walkers[i_beta][i_walker].parameter[i_par] << "\n";
  }
  
  
  
  
  
};


template <class Parameters, class Duration>
 bool extract_iter(std::istream& f,std::size_t& iter, Duration &dur,
                          thermo_mcmc<Parameters> &data) {
    for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
        for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker)
            for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)
            {
                double v_beta;
                std::size_t v_i_walker;
                std::size_t v_walker_id;
                std::size_t v_i_par;
                double v_param_value;
                if (load_vars_line(f, iter, dur
                               , v_beta  , v_i_walker 
                               , v_walker_id ,  v_i_par 
                               , v_param_value  ))
                {
                    if (std::abs(v_beta-data.beta[i_beta])>eps*100)
                        return false;
                    if (v_i_walker!=i_walker)
                        return false;
                    if (v_i_par!=i_par)
                        return false;
                    data.i_walkers[i_beta][i_walker]=v_walker_id;
                    data.walkers[i_beta][i_walker].parameter[i_par]=v_param_value;
                }
                else
                {return false;}
            }
    return true;
 }
template <class Parameters, class Duration>
 bool extract_iter(const std::string line,std::size_t& iter, Duration &dur,
                  thermo_mcmc<Parameters> &data) {
    std::stringstream ss(line);
    return extract_iter(ss,iter,dur,data);
} 
 
template <class Parameters, class Duration>
 auto extract_parameters_last(const std::string& fname,std::size_t& iter, Duration &dur,
                  thermo_mcmc<Parameters> &data) {
    auto candidate=data;
    auto candidate_prev=data;
    auto f=std::ifstream(fname);
    std::string line;
    std::getline(f,line);
    
    while (extract_iter(f,iter,dur,candidate))
    {
        std::swap(data,candidate);
        std::swap(candidate_prev,candidate);
    }
    return candidate_prev;
    
}


template <class Parameters> class save_Predictions {

public:
  std::string sep = ",";
  std::ofstream f;
  std::size_t save_every;
  save_Predictions(std::string const &path, std::size_t interval)
      : f{std::ofstream(path + "__i_beta__i_walker__i_x.csv")},
        save_every{interval} {
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  }

  template <class Prior, class Likelihood, class Variables, class DataType>
  friend void report_model(save_Predictions &, const Prior &,
                           const Likelihood &, const DataType &,
                           const Variables &, by_beta<double> const &) {}

  friend void report_model(save_Predictions &, ...) {}
};

namespace macrodr{
void report_title(
    save_Predictions<var::Parameters_transformed> &s,
    thermo_mcmc<var::Parameters_transformed> const &, ...);
}
template <class Parameters, class... saving>
class save_mcmc : public observer, public saving... {

  std::string directory_;
  std::string filename_prefix_;

public:
  template <typename... Size>
    requires((std::integral<Size> && ...) &&
             (sizeof...(saving) == sizeof...(Size)))
  save_mcmc(std::string dir, std::string filename_prefix,
            Size... sampling_intervals)
      : saving{dir + filename_prefix, sampling_intervals}..., directory_{dir},
        filename_prefix_{filename_prefix} {}

  save_mcmc(std::string dir, std::string filename_prefix)
      : saving{dir + filename_prefix, 1ul}..., directory_{dir},
        filename_prefix_{filename_prefix} {}

  template <class FunctionTable, class Duration, class... T>
  friend void report(FunctionTable &&f, std::size_t iter, const Duration &dur,
                     save_mcmc &smcmc, thermo_mcmc<Parameters> const &data,
                     T &&...ts) {
    (report(f, iter, dur, static_cast<saving &>(smcmc), data,
            std::forward<T>(ts)...),
     ..., 1);
  }
  template <class Prior, class Likelihood, class Variables, class DataType>
  friend void report_model(save_mcmc &s, Prior const &prior,
                           Likelihood const &lik, const DataType &y,
                           const Variables &x, by_beta<double> const &beta0) {
    (report_model(static_cast<saving &>(s), prior, lik, y, x, beta0), ..., 1);
  }

  friend void report_title(save_mcmc &f, thermo_mcmc<Parameters> const &data,
                           ...) {
      using namespace macrodr;
    (report_title(static_cast<saving &>(f), data), ..., 1);
  }
};

template <class Parameter, class... saving>
void report_model_all(save_mcmc<Parameter, saving...> &) {}

template <class Parameter, class... saving, class T, class... Ts>
void report_model_all(save_mcmc<Parameter, saving...> &s, T const &t,
                      Ts const &...ts) {
  (report_model(static_cast<saving &>(s), t), ..., report_model_all(s, ts...));
}

class thermo_less_than_max_iteration {
  std::size_t current_iteration_;
  std::size_t max_iter_final_;

public:
  thermo_less_than_max_iteration(std::size_t max_iter_final)
      : current_iteration_{0ul}, max_iter_final_{max_iter_final} {}

  thermo_less_than_max_iteration &operator++() {
    ++current_iteration_;
    return *this;
  }
  std::size_t current_iteration() const { return current_iteration_; }
  std::size_t max_iteration_final() const { return max_iter_final_; }

  bool stop() const { return current_iteration() >= max_iteration_final(); }
  template <class Anything>
  friend auto checks_convergence(thermo_less_than_max_iteration &&c,
                                 const Anything &) {
    if (c.stop()) {
      return std::pair(std::move(c), true);

    } else {
      ++c;
      return std::pair(std::move(c), false);
    }
  }

  void reset() { current_iteration_ = 0; }
  void reset(std::size_t iter) { current_iteration_ = iter; }

  template <class P>
  friend void
  report_finalizer_data(save_Parameter<P> &s,
                        const thermo_less_than_max_iteration &mcmc) {
    s.f << s.sep << mcmc.current_iteration();
  }

  template <class P>
  friend void report_finalizer_title(save_Parameter<P> &s,
                                     const thermo_less_than_max_iteration &) {
    s.f << s.sep << "current_iter";
  }
};

class cuevi_less_than_max_iteration {
  std::size_t current_iteration_;
  std::size_t max_iter_final_;

public:
  cuevi_less_than_max_iteration(std::size_t max_iter_final)
      : current_iteration_{0ul}, max_iter_final_{max_iter_final} {}

  cuevi_less_than_max_iteration &operator++() {
    ++current_iteration_;
    return *this;
  }
  std::size_t current_iteration() const { return current_iteration_; }
  std::size_t max_iteration_final() const { return max_iter_final_; }

  bool stop() const { return current_iteration() >= max_iteration_final(); }
  template <class Anything>
  friend auto checks_convergence(cuevi_less_than_max_iteration &&c,
                                 const Anything &) {
    if (c.stop()) {
      return std::pair(std::move(c), true);

    } else {
      ++c;
      return std::pair(std::move(c), false);
    }
  }

  void reset() { current_iteration_ = 0; }

  template <class P>
  friend void report_finalizer_data(save_Parameter<P> &s,
                                    const cuevi_less_than_max_iteration &mcmc) {
    s.f << s.sep << mcmc.current_iteration();
  }

  template <class P>
  friend void report_finalizer_title(save_Parameter<P> &s,
                                     const cuevi_less_than_max_iteration &) {
    s.f << s.sep << "current_iter";
  }
};

#endif // PARALLEL_TEMPERING_H
