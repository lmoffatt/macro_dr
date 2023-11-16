#ifndef CUEVI_H
#define CUEVI_H
#include "bayesian_linear_regression.h"
#include "function_measure_verification_and_optimization.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "random_samplers.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

template <class T> using by_fraction = std::vector<T>;

template <class Parameters> struct mcmc2 : public mcmc<Parameters> {
  double logPa;
};

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

by_fraction<double>
calculate_Evidence(by_fraction<by_beta<double>> const &beta,
                   by_fraction<by_beta<double>> const &meanLik) {
  auto nfraction = beta.size();
  auto out = by_fraction<double>(nfraction, 0.0);
  for (std::size_t i_frac = 0; i_frac < nfraction; ++i_frac) {
    out[i_frac] = calculate_Evidence(beta[i_frac], meanLik[i_frac]);
  }
  return out;
}

by_fraction<double>
calculate_Evidence(by_fraction<by_beta<double>> const &beta,
                   by_fraction<by_beta<double>> const &meanLik,
                   by_fraction<by_beta<double>> const &varLik) {
  auto nfraction = beta.size();
  auto out = by_fraction<double>(nfraction, 0.0);
  for (std::size_t i_frac = 0; i_frac < nfraction; ++i_frac) {
    out[i_frac] =
        calculate_Evidence(beta[i_frac], meanLik[i_frac], varLik[i_frac]);
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
                a_n * log(b_n) + std::lgamma(a_n) - std::lgamma(a_0);
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
void report(FunctionTable &&f, std::size_t iter, save_likelihood<Parameters> &s,
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
void report(FunctionTable &&f, std::size_t iter, save_Parameter<Parameters> &s,
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
    expected_partial_evidence_by_logLik[i_frac] = calculate_Evidence(
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

template <class FunctionTable, class Parameters, class T>
void report(FunctionTable &&f, std::size_t iter, save_Evidence &s,
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
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                  Observer &obs, ensemble<std::mt19937_64> &mt,
                  std::vector<std::uniform_real_distribution<double>> &rdist,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, std::size_t n_par,
                  std::size_t i, std::size_t iw, std::size_t jw, std::size_t ib,
                  std::size_t i_fr) const {
    if (current.is_active[i_fr][ib] == 1) {
      auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
      auto r = rdist[i](mt[i]);
      // candidate[ib].walkers[iw].
      auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                                 current.walkers[jw][i_fr][ib].parameter, z);

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
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void
  operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current, Observer &obs,
             ensemble<std::mt19937_64> &mt, Prior const &prior,
             Likelihood const &lik, const by_fraction<DataType> &y,
             const by_fraction<Variables> &x, double alpha_stretch = 2) const {
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
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, cuevi_mcmc<Parameters> &current,
                  Observer &obs, ensemble<std::mt19937_64> &mt,
                  std::vector<std::uniform_real_distribution<double>> &rdist,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x, std::size_t n_par,
                  std::size_t i, std::size_t iw, std::size_t jw, std::size_t ib,
                  std::size_t i_fr) {
    if (current.is_active[i_fr][ib] == 1) {
      auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
      auto r = rdist[i](mt[i]);
      // candidate[ib].walkers[iw].
      auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                                 current.walkers[jw][i_fr][ib].parameter, z);

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

auto generate_random_Indexes(std::mt19937_64 &mt, std::size_t num_samples,
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
      auto f = *it;
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

struct fractioner {
  auto operator()(const Matrix<double> &y, const Matrix<double> &x,
                  std::mt19937_64 &mt, std::size_t num_parameters,
                  double n_points_per_decade_beta,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

auto init_mcmc2(FunctionTable &&f, std::mt19937_64 &mt, const Prior &prior,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

auto init_mcmc_resample(FunctionTable &&f, ensemble<std::mt19937_64> &mt,
                        cuevi_mcmc<Parameters> &current, const Prior &prior,
                        const Likelihood &lik, const by_fraction<DataType> &y,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

auto init_cuevi_mcmc(FunctionTable &&f, std::size_t n_walkers,
                     by_beta<double> const &beta, ensemble<std::mt19937_64> &mt,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

auto init_cuevi_mcmc_all(FunctionTable &&f, std::size_t n_walkers,
                         by_fraction<by_beta<double>> const &beta_out,
                         ensemble<std::mt19937_64> &mt, Prior const &prior,
                         Likelihood const &lik, const by_fraction<DataType> &y,
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

  auto r_i_walker = 0ul;
  for (std::size_t i_frac = 0; i_frac < beta_out.size(); ++i_frac)
    for (std::size_t ib = 0; ib < beta_out[i_frac].size(); ++ib)
      for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        i_walker[iw][i_frac][ib] = r_i_walker;
        ++r_i_walker;
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)

Maybe_error<bool>
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
auto create_new_walkers(FunctionTable &&f,
                        const cuevi_mcmc<Parameters> &current,
                        ensemble<std::mt19937_64> &mts, Prior const &prior,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
  requires(is_prior<Prior, Parameters, Variables, DataType> &&
           is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                               DataType>)
Maybe_error<cuevi_mcmc<Parameters>> push_back_new_fraction(
    FunctionTable &&f, const cuevi_mcmc<Parameters> &current_old,
    ensemble<std::mt19937_64> &mts,
    const by_fraction<by_beta<double>> &final_beta,
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
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, std::size_t iter,
                  cuevi_mcmc<Parameters> &current, Observer &obs,
                  std::mt19937_64 &mt, ensemble<std::mt19937_64> &mts,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  std::size_t thermo_jumps_every) const {
    if (iter % (thermo_jumps_every) == 0) {
      std::uniform_real_distribution<double> uniform_real(0, 1);
      auto n_walkers = mts.size() * 2;
      auto n_par = current.walkers[0][0][0].parameter.size();
      std::uniform_int_distribution<std::size_t> booldist(0, 1);
      auto half = booldist(mt) == 1;

      WalkerIndexes landing_walker(n_walkers / 2);
      std::iota(landing_walker.begin(), landing_walker.end(), 0);
      std::shuffle(landing_walker.begin(), landing_walker.end(), mt);
      std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                                uniform_real);

#pragma omp parallel for // not currently working
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto iw = half ? i + n_walkers / 2 : i;
        auto j = landing_walker[i];
        auto jw = half ? j : j + n_walkers / 2;

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
  template <class FunctionTable, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
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
        out0.parameter=ca_par;  
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
            out1.parameter=ca_par;  
            out1.logPa = ca_logPa;
            out1.logP = ca_logP1;
            out1.logL = ca_logL1;
          }
        }
      } else  {
        out1.parameter=ca_par;  
        out1.logPa = ca_logPa;
        out1.logP = ca_logP0;
        out1.logL = ca_logL0;
        if ((i_frac > 0)&&(beta[i_frac][ib] == 0.0)) {
          auto ca_logL_00 =i_frac > 1 ?
              f.f(logLikelihood_f{}, lik, ca_par, y[i_frac - 2], x[i_frac - 2]):
                                  Maybe_error<double>(0.0);
          if (!(ca_logL_00))
            return ca_logL_00.error();
          else {
            auto ca_logP00 = ca_logPa + ca_logL_00.value();
            auto ca_logL00 = ca_logL_0.value() - ca_logL_00.value();
            out0.parameter=ca_par;  
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
  template <class FunctionTable, class Prior, class Likelihood,
            class Variables, class DataType,
            class Parameters = std::decay_t<decltype(sample(
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
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
                std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters,
                                 Variables, DataType>)

  void operator()(FunctionTable &&f, std::size_t iter,
                  cuevi_mcmc<Parameters> &current, Observer &obs,
                  std::mt19937_64 &mt, ensemble<std::mt19937_64> &mts,
                  Prior const &prior, Likelihood const &lik,
                  const by_fraction<DataType> &y,
                  const by_fraction<Variables> &x,
                  std::size_t thermo_jumps_every) const {
    if (iter % (thermo_jumps_every) == 0) {
      std::uniform_real_distribution<double> uniform_real(0, 1);
      auto n_walkers = mts.size() * 2;
      auto n_par = current.walkers[0][0][0].parameter.size();
      std::uniform_int_distribution<std::size_t> booldist(0, 1);
      auto half = booldist(mt) == 1;

      WalkerIndexes landing_walker(n_walkers / 2);
      std::iota(landing_walker.begin(), landing_walker.end(), 0);
      std::shuffle(landing_walker.begin(), landing_walker.end(), mt);
      std::vector<std::uniform_real_distribution<double>> rdist(n_walkers,
                                                                uniform_real);

      auto n_fractions = current.beta.size();
      FractionIndexes landing_fraction(n_fractions);
      std::iota(landing_fraction.begin(), landing_fraction.end(), 0);
      std::shuffle(landing_fraction.begin(), landing_fraction.end(), mt);

#pragma omp parallel for // not currently working
      for (std::size_t i = 0; i < n_walkers / 2; ++i) {
        auto iw = half ? i + n_walkers / 2 : i;
        auto j = landing_walker[i];
        auto jw = half ? j : j + n_walkers / 2;

        for (std::size_t i_fr = 0; i_fr < n_fractions; ++i_fr) {
          auto i_frac_origin = i_fr;
          auto i_frac_destination = landing_fraction[i_fr];
          auto r = rdist[i](mts[i]);

          thermo_cuevi_randomized_jump_mcmc_per_walker{}(
              f.fork(var::I_thread(i)), current, iw, jw, i_frac_origin, i_frac_destination, prior, lik,
              y, x, r);
        }
      }
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

template <>
bool compare_to_max_ratio(by_fraction<by_beta<double>> const &beta,
                          by_fraction<by_beta<double>> const &mean_logL,
                          by_fraction<by_beta<double>> const &var_ratio,
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
              std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
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
      f.f(step_stretch_cuevi_mcmc{}, current, rep, mts, prior, lik, ys, xs);
      // check_sanity(iter,current);
      report_point(ff, iter);

      ++iter;
      f.f(thermo_cuevi_jump_mcmc{}, iter, current, rep, mt, mts, prior, lik, ys,
          xs, cue.thermo_jumps_every());
      f.f(thermo_cuevi_randomized_jump_mcmc{}, iter+cue.thermo_jumps_every()%2, current, rep, mt, mts, prior, lik, ys,
          xs, cue.thermo_jumps_every());
      
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

#endif // CUEVI_H
