#ifndef CUEVI_H
#define CUEVI_H
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "random_samplers.h"
#include <algorithm>
#include <cassert>
#include <random>
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
    friend void report_title(save_likelihood<Parameters> &s, cuevi_mcmc<Parameters> const &,...) {
        
        s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
            << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep
            << "id_walker" << s.sep << "logPa"  << s.sep << "logP" << s.sep << "logLik"
            << "\n";
    }
    friend void report_title(save_Evidence &s, cuevi_mcmc const &,...) {
        
        s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
            << "nsamples" << s.sep << "beta" << s.sep << "meanPrior" << s.sep
            << "meanLik" << s.sep << "varLik" << s.sep
            << "fraction_Evidence_by_mean" << s.sep << "fraction_Evidence_by_var"
            << s.sep << "Evidence_by_mean" << s.sep << "Evidence_by_var"
            << "\n";
    }
    
    friend void report_title(save_Parameter<Parameters> &s, cuevi_mcmc const &,...) {
        
        s.f << "n_fractions" << s.sep
            <<"n_betas" << s.sep
            <<"iter" << s.sep
            << "nsamples" << s.sep
            << "beta" << s.sep
            <<"i_walker" << s.sep
            << "id_walker" << s.sep
            << "i_par" << s.sep
            << "par_value"
            << "\n";
    }
    
    
    
    
    
    template <class... saving, class... Ts>
    friend void report_title(save_mcmc<Parameters,saving...> &f,
                             cuevi_mcmc const &data, const Ts&...ts) {
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

template<class Parameters>
void report(std::size_t iter, save_likelihood<Parameters> &s,
            cuevi_mcmc<Parameters> const &data,...) {
    if (iter % s.save_every == 0)
        for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
            for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
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

template<class Parameters>
void report(std::size_t iter, save_Parameter<Parameters> &s,
            cuevi_mcmc<Parameters> const &data,...) {
    if (iter % s.save_every == 0)
        for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
            for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
                for (std::size_t i_walker = 0; i_walker < size(data.walkers);
                     ++i_walker)
                    for (std::size_t i_par = 0;
                         i_par < size(data.walkers[i_walker][i_frac][i_beta].parameter);
                         ++i_par)
                        
                        
                        
                        s.f << size(data.beta) << s.sep
                            << size(data.beta[i_frac]) << s.sep
                            << iter << s.sep
                            << data.nsamples[i_frac] << s.sep
                            << data.beta[i_frac][i_beta] << s.sep
                            << i_walker << s.sep
                            << data.i_walkers[i_walker][i_frac][i_beta] << s.sep
                            << i_par<< s.sep
                            << data.walkers[i_walker][i_frac][i_beta].parameter[i_par]
                            << "\n";
}






template <class Prior, class Likelihood, class Variables, class DataType>
concept has_conjugate= requires (Prior const &prior, Likelihood const &lik,
                                 const DataType &y, const Variables &x)
{
    {posterior(conjugate{}, prior, lik, y[0], x[0])};
};








template <class Prior, class Likelihood, class Variables, class DataType>
    requires has_conjugate<Prior,Likelihood,Variables,DataType>
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
    requires (!has_conjugate<Prior,Likelihood,Variables,DataType>)
void report_model(save_Evidence &, Prior const &, Likelihood const &,
                  const DataType &, const Variables &,
                  by_fraction<by_beta<double>> const &) {
    
}








template<class Parameters>
void report(std::size_t iter, save_Evidence &s,
            cuevi_mcmc<Parameters> const &data,...) {
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

template <class Parameters, class... saving, class... T>
void report(std::size_t iter, save_mcmc<Parameters,saving...> &f,
            cuevi_mcmc<Parameters> const &data,T const &...ts) {
    (report(iter, static_cast<saving &>(f), data, ts...), ..., 1);
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void step_stretch_cuevi_mcmc(
    cuevi_mcmc<Parameters> &current, Observer &obs,
    ensemble<std::mt19937_64> &mt,
    std::vector<std::uniform_real_distribution<double>> &rdist,
    Prior const &prior, Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x, std::size_t n_par, std::size_t i,
    std::size_t iw, std::size_t jw, std::size_t ib, std::size_t i_fr) {
    auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
    auto r = rdist[i](mt[i]);
    // candidate[ib].walkers[iw].
    auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                               current.walkers[jw][i_fr][ib].parameter, z);
    
    auto ca_logPa_ = logPrior(prior, ca_par);
    auto ca_logL_0 = i_fr > 0
                         ? logLikelihood(lik, ca_par, y[i_fr - 1], x[i_fr - 1])
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
            if (i_fr + 1 < size(current.beta) && (current.beta[i_fr][ib] == 1.0)) {
                auto ca_logL_2 = logLikelihood(lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
                if ((ca_logL_2)) {
                    auto ca_logP1 = ca_logPa + ca_logL_1.value();
                    auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
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
}
template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)
void double_step_stretch_cuevi_mcmc(
    cuevi_mcmc<Parameters> &current, Observer &obs,
    ensemble<std::mt19937_64> &mt,
    std::vector<std::uniform_real_distribution<double>> &rdist,
    Prior const &prior, Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x, std::size_t n_par, std::size_t i,
    std::size_t iw, std::size_t jw, std::size_t ib, std::size_t i_fr)

{
    
    auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
    auto r = rdist[i](mt[i]);
    
    // candidate[ib].walkers[iw].
    auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                               current.walkers[jw][i_fr][ib].parameter, z);
    auto ca_logP = logPrior(prior, ca_par);
    auto ca_logL0 = logLikelihood(lik, ca_par, y[i_fr], x[i_fr]);
    auto ca_logL1 = logLikelihood(lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
    
    if ((ca_logP) && (ca_logL0) && (ca_logL1)) {
        auto dthLogL = ca_logP.value() - current.walkers[iw][i_fr][ib].logP +
                       current.beta[i_fr][ib] *
                           (ca_logL0.value() - current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        observe_step_stretch_thermo_mcmc(
            obs[iw][i_fr][ib], jw, z, r, current.walkers[iw][i_fr][ib].parameter,
            current.walkers[jw][i_fr][ib].parameter,
            current.walkers[iw][i_fr][ib].logP, ca_logP,
            current.walkers[iw][i_fr][ib].logL, ca_logL1, pJump >= r);
        if (pJump >= r) {
            current.walkers[iw][i_fr][ib].parameter = ca_par;
            current.walkers[iw][i_fr][ib].logPa = ca_logP.value();
            current.walkers[iw][i_fr][ib].logP = ca_logP.value();
            current.walkers[iw][i_fr][ib].logL = ca_logL0.value();
            current.walkers[iw][i_fr + 1][0].parameter = std::move(ca_par);
            current.walkers[iw][i_fr + 1][0].logPa = ca_logP.value();
            current.walkers[iw][i_fr + 1][0].logP =
                ca_logP.value() + ca_logL0.value();
            current.walkers[iw][i_fr + 1][0].logL =
                ca_logL1.value() - ca_logL0.value();
        }
    }
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void middle_step_stretch_cuevi_mcmc(
    cuevi_mcmc<Parameters> &current, Observer &obs,
    ensemble<std::mt19937_64> &mt,
    std::vector<std::uniform_real_distribution<double>> &rdist,
    Prior const &prior, Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x, std::size_t n_par, std::size_t i,
    std::size_t iw, std::size_t jw, std::size_t ib, std::size_t i_fr)

{
    
    auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
    auto r = rdist[i](mt[i]);
    
    // candidate[ib].walkers[iw].
    auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                               current.walkers[jw][i_fr][ib].parameter, z);
    auto ca_logP = logPrior(prior, ca_par);
    auto ca_logL0 = logLikelihood(lik, ca_par, y[i_fr - 1], x[i_fr - 1]);
    auto ca_logL1 = logLikelihood(lik, ca_par, y[i_fr], x[i_fr]);
    
    if ((ca_logP) && (ca_logL0) && (ca_logL1)) {
        auto dthLogL =
            ca_logP.value() + ca_logL0.value() -
            current.walkers[iw][i_fr][ib].logP +
            current.beta[i_fr][ib] * (ca_logL1.value() - ca_logL0.value() -
                                      current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        observe_step_stretch_thermo_mcmc(
            obs[iw][ib], jw, z, r, current.walkers[iw][i_fr][ib].parameter,
            current.walkers[jw][i_fr][ib].parameter,
            current.walkers[iw][i_fr][ib].logP, ca_logP,
            current.walkers[iw][i_fr][ib].logL, ca_logL0, pJump >= r);
        if (pJump >= r) {
            current.walkers[iw][i_fr][ib].parameter = std::move(ca_par);
            current.walkers[iw][i_fr][ib].logPa = ca_logP.value();
            current.walkers[iw][i_fr][ib].logP = ca_logP.value() + ca_logL0.value();
            current.walkers[iw][i_fr][ib].logL = ca_logL1.value() - ca_logL0.value();
        }
    }
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void triple_step_stretch_cuevi_mcmc(
    cuevi_mcmc<Parameters> &current, Observer &obs,
    ensemble<std::mt19937_64> &mt,
    std::vector<std::uniform_real_distribution<double>> &rdist,
    Prior const &prior, Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x, std::size_t n_par, std::size_t i,
    std::size_t iw, std::size_t jw, std::size_t ib, std::size_t i_fr) {
    
    auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
    auto r = rdist[i](mt[i]);
    
    // candidate[ib].walkers[iw].
    auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                               current.walkers[jw][i_fr][ib].parameter, z);
    auto ca_logP = logPrior(prior, ca_par);
    auto ca_logL0 = logLikelihood(lik, ca_par, y[i_fr - 1], x[i_fr - 1]);
    auto ca_logL1 = logLikelihood(lik, ca_par, y[i_fr], x[i_fr]);
    auto ca_logL2 = logLikelihood(lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
    
    if ((ca_logP) && (ca_logL0) && (ca_logL1) && (ca_logL2)) {
        auto dthLogL =
            ca_logP.value() + ca_logL0.value() -
            current.walkers[iw][i_fr][ib].logP +
            current.beta[i_fr][ib] * (ca_logL1.value() - ca_logL0.value() -
                                      current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        observe_step_stretch_thermo_mcmc(
            obs[iw][ib], jw, z, r, current.walkers[iw][i_fr][ib].parameter,
            current.walkers[jw][i_fr][ib].parameter,
            current.walkers[iw][i_fr][ib].logP, ca_logP,
            current.walkers[iw][i_fr][ib].logL, ca_logL0, pJump >= r);
        if (pJump >= r) {
            current.walkers[iw][i_fr][ib].parameter = ca_par;
            current.walkers[iw][i_fr][ib].logPa = ca_logP.value();
            current.walkers[iw][i_fr][ib].logP = ca_logP.value() + ca_logL0.value();
            current.walkers[iw][i_fr][ib].logL = ca_logL1.value() - ca_logL0.value();
            current.walkers[iw][i_fr + 1][0].parameter = std::move(ca_par);
            current.walkers[iw][i_fr + 1][0].logPa = ca_logP.value();
            current.walkers[iw][i_fr + 1][0].logP =
                ca_logP.value() + ca_logL1.value();
            current.walkers[iw][i_fr + 1][0].logL =
                ca_logL2.value() - ca_logL1.value();
        }
    }
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void last_step_stretch_cuevi_mcmc(
    cuevi_mcmc<Parameters> &current, Observer &obs,
    ensemble<std::mt19937_64> &mt,
    std::vector<std::uniform_real_distribution<double>> &rdist,
    Prior const &prior, Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x, std::size_t n_par, std::size_t i,
    std::size_t iw, std::size_t jw, std::size_t ib, std::size_t i_fr) {
    
    auto z = std::pow(rdist[i](mt[i]) + 1, 2) / 2.0;
    auto r = rdist[i](mt[i]);
    
    auto ca_par = stretch_move(current.walkers[iw][i_fr][ib].parameter,
                               current.walkers[jw][i_fr][ib].parameter, z);
    auto ca_logP = logPrior(prior, ca_par);
    auto ca_logL0 = logLikelihood(lik, ca_par, y[i_fr - 1], x[i_fr - 1]);
    auto ca_logL1 = logLikelihood(lik, ca_par, y[i_fr], x[i_fr]);
    
    if ((ca_logP) && (ca_logL0) && (ca_logL1)) {
        auto dthLogL =
            ca_logP.value() + ca_logL0.value() -
            current.walkers[iw][i_fr][ib].logP +
            current.beta[i_fr][ib] * (ca_logL1.value() - ca_logL0.value() -
                                      current.walkers[iw][i_fr][ib].logL);
        auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
        observe_step_stretch_thermo_mcmc(
            obs[iw][ib], jw, z, r, current.walkers[iw][i_fr][ib].parameter,
            current.walkers[jw][i_fr][ib].parameter,
            current.walkers[iw][i_fr][ib].logP, ca_logP,
            current.walkers[iw][i_fr][ib].logL, ca_logL0, pJump >= r);
        if (pJump >= r) {
            current.walkers[iw][i_fr][ib].parameter = std::move(ca_par);
            current.walkers[iw][i_fr][ib].logPa = ca_logP.value();
            current.walkers[iw][i_fr][ib].logP = ca_logP.value() + ca_logL0.value();
            current.walkers[iw][i_fr][ib].logL = ca_logL1.value() - ca_logL0.value();
        }
    }
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void step_stretch_cuevi_mcmc(cuevi_mcmc<Parameters> &current, Observer &obs,
                             ensemble<std::mt19937_64> &mt, Prior const &prior,
                             Likelihood const &lik,
                             const by_fraction<DataType> &y,
                             const by_fraction<Variables> &x,
                             double alpha_stretch = 2) {
    assert(current.beta.size() == num_betas(current));
    auto n_walkers = size(current.walkers);
    
    auto n_par = size(current.walkers[0][0][0].parameter);
    
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
            step_stretch_cuevi_mcmc(current, obs, mt, rdist, prior, lik, y, x, n_par,
                                    i, iw, jw, 0, 0);
            for (std::size_t i_fr = 0; i_fr < size(current.beta); ++i_fr) {
                for (std::size_t ib = 1; ib < size(current.beta[i_fr]); ++ib)
                    step_stretch_cuevi_mcmc(current, obs, mt, rdist, prior, lik, y, x,
                                            n_par, i, iw, jw, ib, i_fr);
            }
        }
}

using DataIndexes = std::vector<std::size_t>;

auto generate_random_Indexes(std::mt19937_64 &mt, std::size_t num_samples,
                             std::size_t min_num_extra_samples,
                             double num_jumps_per_decade, std::vector<std::size_t> initial_samples={}) {
    
    std::size_t num_initial_samples= size(initial_samples);
    std::size_t n_jumps = std::max(
        0.0, std::floor(num_jumps_per_decade * (std::log10(num_samples) -
                                           std::log10(min_num_extra_samples+num_initial_samples))));
    auto indexsizes = DataIndexes(n_jumps+1);
    
    for (std::size_t i = 0; i < n_jumps+1; ++i)
        indexsizes[i] = num_samples * std::pow(10.0, -(1.0 * (n_jumps - i)) /
                                                         num_jumps_per_decade);
    auto out = std::vector<DataIndexes>(n_jumps+1);
    if (n_jumps > 0) {
        auto index = DataIndexes(num_samples);
        std::iota(index.begin(), index.end(), 0u);
        auto it=index.begin();
        if (num_initial_samples>0)
        {
            auto new_index=DataIndexes{};
            std::copy(initial_samples.begin(), initial_samples.end(), std::back_inserter(new_index));
            std::set_difference(index.begin(),index.end(),initial_samples.begin(),initial_samples.end(), std::back_inserter(new_index));
            
            std::swap(index,new_index);
            it=index.begin();
            std::advance(it,initial_samples.size());
            auto f=*it;
        }
        it = randomly_extract_n(mt, it, index.end(), indexsizes[0]-num_initial_samples);
        auto res=DataIndexes(index.begin(), it);
        std::sort(res.begin(), res.end());
        out[0] = std::move(res);
        for (auto i = 1u; i < n_jumps+1; ++i) {
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
        auto indexes = generate_random_Indexes(mt, num_samples, 2*num_parameters,
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



template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

auto init_mcmc2(std::mt19937_64 &mt, const Prior &prior, const Likelihood &lik,
                const by_fraction<DataType> &y,
                const by_fraction<Variables> &x) {
    auto prior_sampler = sampler(prior);
    auto par = sample(mt, prior_sampler);
    auto logP = logPrior(prior, par);
    auto logL = logLikelihood(lik, par, y[0], x[0]);
    auto logPa = logP;
    while (!(logP) || !(logL)) {
        // std::cerr<<"\npar\ņ"<<par;
        //std::cerr<<"\nlogL\ņ"<<logL;
        
        par = sample(mt, prior_sampler);
        logP = logPrior(prior, par);
        logL = logLikelihood(lik, par, y[0], x[0]);
    }
    return mcmc2<Parameters>{mcmc<Parameters>{std::move(par), logP.value(), logL.value()}, logPa.value()};
}

template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

auto init_cuevi_mcmc(std::size_t n_walkers, by_beta<double> const &beta,
                     ensemble<std::mt19937_64> &mt, Prior const &prior,
                     Likelihood const &lik, const by_fraction<DataType> &y,
                     const by_fraction<Variables> &x) {
    by_fraction<std::size_t> nsamples_out(1, size(y[0]));
    by_fraction<by_beta<double>> beta_out(1, beta);
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
                i_walker[iw][0][i] = iw + (beta.size()-i-1) * n_walkers;
                walker[iw][0][i] = init_mcmc2(mt[iiw], prior, lik, y, x);
            }
        }
    return cuevi_mcmc<Parameters>{nsamples_out, beta_out, walker, i_walker};
}

template <class Parameters>
std::size_t last_walker(const cuevi_mcmc<Parameters> &c) {
    std::size_t tot = 0;
    for (std::size_t i_frac = 0; i_frac < size(c.walkers[0]); ++i_frac)
        tot += size(c.walkers[0][i_frac]);
    return tot * size(c.walkers);
}


template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

Maybe_error<bool> calculate_Likelihoods_sample(
    cuevi_mcmc<Parameters> &current,
    Prior const &prior,
    Likelihood const &lik,
    const by_fraction<DataType> &y,
    const by_fraction<Variables> &x,
    std::size_t iw,
    std::size_t i_frac,
    std::size_t ib) {
    auto const & ca_par = current.walkers[iw][i_frac][ib].parameter;
    auto ca_logPa_ = logPrior(prior, ca_par);
    auto ca_logL_0 = i_frac > 0
                         ? logLikelihood(lik, ca_par, y[i_frac - 1], x[i_frac - 1])
                         : Maybe_error(0.0);
    auto ca_logL_1 = logLikelihood(lik, ca_par, y[i_frac], x[i_frac]);
    if (!(is_valid(ca_logPa_) && is_valid(ca_logL_0) && is_valid(ca_logL_1)))
    {
        return error_message(ca_logPa_.error()()+ca_logL_0.error()()+ca_logL_1.error()());
    }
    else{
        auto ca_logPa = ca_logPa_.value();
        auto ca_logP0 = ca_logPa_.value() + ca_logL_0.value();
        auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();
        if (i_frac + 1 < size(current.walkers[iw]) && (current.beta[i_frac][ib] == 1.0)) {
            auto ca_logL_2 = logLikelihood(lik, ca_par, y[i_frac + 1], x[i_frac + 1]);
            if (!(ca_logL_2))
                return ca_logL_2.error();
            else    
            {
               // assert(test_equality(ca_par,current.walkers[iw][i_frac + 1][0].parameter, eps));
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


template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)
Maybe_error<cuevi_mcmc<Parameters>> calculate_current_Likelihoods(
    cuevi_mcmc<Parameters> &current,
    Prior const &prior,
    Likelihood const &lik,
    const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {
    
    for (std::size_t iw=0; iw<current.walkers.size(); ++iw)
    {
        auto res=calculate_Likelihoods_sample(current,prior,lik,y,x,iw,0,0);
        if (!res) return res.error();
        for (std::size_t i_frac=0; i_frac<current.walkers[iw].size(); ++i_frac)
            for (std::size_t ib=1; ib<current.walkers[iw][i_frac].size(); ++ib)
            {
                auto res=calculate_Likelihoods_sample(current,prior,lik,y,x,iw,i_frac,ib);
                if (!res) return res.error();
                
            }
    }
    return current;
}


template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)
auto create_new_walkers(
    const cuevi_mcmc<Parameters> &current, ensemble<std::mt19937_64> &mts,
    Prior const &prior,
    Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {
    
    auto n_walkers = current.walkers.size();
    auto sum_walkers = last_walker(current);
    ensemble<mcmc2<Parameters>> new_walkers(n_walkers);
    ensemble<std::size_t> new_i_walkers(n_walkers);
    for (std::size_t half = 0; half < 2; ++half)
        for (std::size_t i = 0; i < n_walkers / 2; ++i) {
            auto iw = i + half * n_walkers / 2;
            new_walkers[iw] = init_mcmc2(mts[i], prior, lik, y, x);
            new_i_walkers[iw] = sum_walkers + iw;
        }
    
    return std::tuple(new_walkers,new_i_walkers);
}

template <class DataType,
         class Parameters>
void insert_new_walkers(
    cuevi_mcmc<Parameters> &current, 
    const by_fraction<by_beta<double>> &final_beta,  const by_fraction<DataType> &y,ensemble<mcmc2<Parameters>>&& new_walkers, ensemble<std::size_t>&& new_i_walkers) {
    for (std::size_t iw=0; iw<current.walkers.size(); ++iw)
    {
        for (std::size_t ib=0; ib<current.walkers[iw][0].size(); ++ib)
        {
            std::swap(current.walkers[iw][0][ib],new_walkers[iw]);
            std::swap(current.i_walkers[iw][0][ib],new_i_walkers[iw]);
        }
        for (std::size_t i_frac=1; i_frac<current.walkers[iw].size(); ++i_frac)
        {
            current.walkers[iw][i_frac][0]=current.walkers[iw][i_frac-1].back();
            current.i_walkers[iw][i_frac][0]=current.i_walkers[iw][i_frac-1].back();
           
            for (std::size_t ib=1; ib<current.walkers[iw][i_frac].size(); ++ib)
            {
                std::swap(current.walkers[iw][i_frac][ib],new_walkers[iw]);
                std::swap(current.i_walkers[iw][i_frac][ib],new_i_walkers[iw]);
            }
        }
        auto i_frac=current.walkers[iw].size()-1;
        auto ib=current.walkers[iw][i_frac].size()-1;
        
        if (current.beta[i_frac][ib]<1.0)
        {
            current.walkers[iw][i_frac].push_back(new_walkers[iw]);
            current.i_walkers[iw][i_frac].push_back(new_i_walkers[iw]);
            if (iw==0)
                current.beta[i_frac].push_back(final_beta[i_frac][ib+1]);
        }
        else
        {
            current.walkers[iw].push_back(by_beta<mcmc2<Parameters>>(2));
            current.i_walkers[iw].push_back(by_beta<std::size_t>(2));
            current.walkers[iw][i_frac+1][0]=current.walkers[iw][i_frac].back();
            current.i_walkers[iw][i_frac+1][0]=current.i_walkers[iw][i_frac].back();
            
            current.walkers[iw][i_frac+1][1]=new_walkers[iw];
            current.i_walkers[iw][i_frac+1][1]=new_i_walkers[iw];
            
            if (iw==0)
            {
                current.beta.push_back(by_beta<double>{final_beta[i_frac+1][0],final_beta[i_frac+1][1]});
                current.nsamples.push_back(size(y[i_frac+1]));
            }
        }
    }
}

template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)
Maybe_error<cuevi_mcmc<Parameters>> push_back_new_fraction(
    const cuevi_mcmc<Parameters> &current_old, ensemble<std::mt19937_64> &mts,
    const by_fraction<by_beta<double>> &final_beta, Prior const &prior,
    Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {
    
    auto current=current_old;
    
    auto [new_walkers,new_i_walkers]=create_new_walkers(current,mts,prior,lik,y,x);
    
    insert_new_walkers(current,final_beta,y,std::move(new_walkers), std::move(new_i_walkers));
    
    return calculate_current_Likelihoods(current,prior,lik,y,x);
    
}





template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)




Maybe_error<cuevi_mcmc<Parameters>> push_back_new_fraction_old(
    const cuevi_mcmc<Parameters> &current_old, ensemble<std::mt19937_64> &mts,
    const by_fraction<by_beta<double>> &final_beta, Prior const &prior,
    Likelihood const &lik, const by_fraction<DataType> &y,
    const by_fraction<Variables> &x) {
    auto current = current_old;
    std::cerr<<"\ncurrent.walkers before\n";
    for (std::size_t i_fr=0; i_fr<size(current.walkers[0]); ++i_fr)
    {
        for (std::size_t i_b=0; i_b<size(current.walkers[0][i_fr]); ++i_b)
            std::cerr<<get_value(current.walkers[0][i_fr][i_b].parameter)[10]<<"\t";
        std::cerr<<"\t\t";
    }
    std::cerr<<"\n";
    auto n_walkers = current.walkers.size();
    auto sum_walkers = last_walker(current);
    ensemble<mcmc2<Parameters>> new_walkers(n_walkers);
    ensemble<std::size_t> new_i_walkers(n_walkers);
    for (std::size_t half = 0; half < 2; ++half)
        for (std::size_t i = 0; i < n_walkers / 2; ++i) {
            auto iw = i + half * n_walkers / 2;
            new_walkers[iw] = init_mcmc2(mts[i], prior, lik, y, x);
            new_i_walkers[iw] = sum_walkers + iw;
        }
    
    
    
    
    
    
    auto swap_walker=[&current,&current_old,&final_beta,&lik,&prior,&y,&x](std::size_t iw,
                                                                                   std::size_t ib,
                                                                                   std::size_t i_frac,
                                                                                   auto& new_walkers,
                                                                                   auto& new_i_walkers)->Maybe_error<bool>
    {
        if ((i_frac<current_old.walkers[iw].size())&&(ib<current_old.walkers[iw][i_frac].size())){  
            std::swap(current.walkers[iw][i_frac][ib], new_walkers[iw]);
            std::swap(current.i_walkers[iw][i_frac][ib], new_i_walkers[iw]);
            
            if ((final_beta[i_frac][ib]==1.0)&&(i_frac+1==current.walkers[iw].size())&&(i_frac+1<final_beta.size()))
            {
                current.i_walkers[iw].push_back(by_beta<std::size_t>(1));
                current.walkers[iw].push_back(by_beta<mcmc2<Parameters>>(1));
                if (iw==0)
                    current.beta.push_back(by_beta<double>{0.0});
            }
        }
        else 
        {
            current.walkers[iw][i_frac].push_back(new_walkers[iw]);
            current.i_walkers[iw][i_frac].push_back(new_i_walkers[iw]);
            if (iw==0)
                current.beta[i_frac].push_back(final_beta[i_frac][ib]);
            
        }
        
        if (final_beta[i_frac][ib]!=1.0){
            return true;
        }
        else
        {
            auto& ca_par = current.walkers[iw][i_frac][ib].parameter;
            auto& i_walker = current.i_walkers[iw][i_frac][ib];
            auto ca_logPa_ = logPrior(prior, ca_par);
            auto ca_logL_0 = i_frac > 0
                                 ? logLikelihood(lik, ca_par, y[i_frac - 1], x[i_frac - 1])
                                 : Maybe_error(0.0);
            auto ca_logL_1 = logLikelihood(lik, ca_par, y[i_frac], x[i_frac]);
            if (!(is_valid(ca_logPa_) && is_valid(ca_logL_0) && is_valid(ca_logL_1)))
            {
                return error_message(ca_logPa_.error()()+ca_logL_0.error()()+ca_logL_1.error()());
            }
            else{
                auto ca_logPa = ca_logPa_.value();
                auto ca_logP0 = ca_logPa_.value() + ca_logL_0.value();
                auto ca_logL0 = ca_logL_1.value() - ca_logL_0.value();
                if (i_frac + 1 < size(current.walkers[iw]) && (current.beta[i_frac][ib] == 1.0)) {
                    auto ca_logL_2 = logLikelihood(lik, ca_par, y[i_frac + 1], x[i_frac + 1]);
                    if (!(ca_logL_2))
                        return ca_logL_2.error();
                    else    
                    {
                        auto ca_logP1 = ca_logPa + ca_logL_1.value();
                        auto ca_logL1 = ca_logL_2.value() - ca_logL_1.value();
                        current.walkers[iw][i_frac][ib].parameter = ca_par;
                        current.walkers[iw][i_frac][ib].logPa = ca_logPa;
                        current.walkers[iw][i_frac][ib].logP = ca_logP0;
                        current.walkers[iw][i_frac][ib].logL = ca_logL0;
                        current.walkers[iw][i_frac + 1][0].parameter = ca_par;
                        current.walkers[iw][i_frac + 1][0].logPa = ca_logPa;
                        current.walkers[iw][i_frac + 1][0].logP = ca_logP1;
                        current.walkers[iw][i_frac + 1][0].logL = ca_logL1;
                        current.i_walkers[iw][i_frac + 1][0] = i_walker;
                        
                        return true;
                    }
                } else {
                    current.walkers[iw][i_frac][ib].parameter = ca_par;
                    current.walkers[iw][i_frac][ib].logPa = ca_logPa;
                    current.walkers[iw][i_frac][ib].logP = ca_logP0;
                    current.walkers[iw][i_frac][ib].logL = ca_logL0;
                    return true;
                }
            }
        }
    };
    
    
    
    for (std::size_t half = 0; half < 2; ++half)
        for (std::size_t i = 0; i < n_walkers / 2; ++i) {
            auto iw = i + half * n_walkers / 2;
            swap_walker(iw,0ul,0ul,new_walkers,new_i_walkers);
            for (std::size_t i_frac=0; i_frac<current_old.walkers[iw].size(); ++i_frac)
            {
                for (std::size_t ib=1; ib<current_old.walkers[iw][i_frac].size(); ++ib)
                {
                    auto res=swap_walker(iw,ib,i_frac,new_walkers,new_i_walkers);
                    if (!res)
                        return res.error();
                }
            }
            auto i_frac=current_old.walkers[iw].size()-1;
            auto i_b=current_old.walkers[iw][i_frac].size()-1;
            if (final_beta[i_frac][i_b]<1.0)
            {
                swap_walker(iw,i_b+1,i_frac,new_walkers,new_i_walkers);     
            }
            else
            {
                swap_walker(iw,1ul,i_frac+1,new_walkers,new_i_walkers);
                if ((iw==0)&& (i_frac+1<size(y)))
                    current.nsamples.push_back(size(y[i_frac+1]));
                
            }    
            
        }
    std::cerr<<"\ncurrent.walkers\n";
    for (std::size_t i_fr=0; i_fr<size(current.walkers[0]); ++i_fr)
    {
        for (std::size_t i_b=0; i_b<size(current.walkers[0][i_fr]); ++i_b)
            std::cerr<<get_value(current.walkers[0][i_fr][i_b].parameter)[10]<<"\t";
        std::cerr<<"\t\t";
    }
    std::cerr<<"\n beta final\n"<<final_beta;
    
    std::cerr<<"\ncurrent.i_walkers[0]-----------------------popo\n"<<current.i_walkers[0];
    std::cerr<<"\ncurrent.nsamples----------------popa\n"<<current.nsamples;
    std::cerr<<"\ncurrent.beta----------------popi\n"<<current.beta;
    return current;    
    
    //    bool we_are_still_filling_the_first_fraction=beta_first < size(final_beta[0]);
    //    bool we_just_completed_the_first_faction=(i_frac_old==0)&&(!we_are_still_filling_the_first_fraction);
    //    if (we_are_still_filling_the_first_fraction) {
    //        for (std::size_t half = 0; half < 2; ++half)
    //            for (std::size_t i = 0; i < n_walkers / 2; ++i) {
    //                auto iw = i + half * n_walkers / 2;
    //                current.walkers[iw][0].push_back(new_walkers[iw]);
    //                current.i_walkers[iw][0].push_back(new_i_walkers[iw]);
    //            }
    //        current.beta[0].push_back(final_beta[0][beta_first]);
    //        return current;
    //    } else {
    
    
    
    //        // we process now the fractions 
    //        for (auto i_frac = 1ul; i_frac < n_frac_old; ++i_frac) {
    //            for (std::size_t i_b = 0; i_b < 1; ++i_b) {
    //                for (std::size_t half = 0; half < 2; ++half)
    //                    for (std::size_t i = 0; i < n_walkers / 2; ++i) {
    //                        auto iw = i + half * n_walkers / 2;
    //                        std::swap(current.walkers[iw][i_frac][i_b], new_walkers[iw]);
    //                        std::swap(current.i_walkers[iw][i_frac][i_b], new_i_walkers[iw]);
    //                        auto &ca_wa = current.walkers[iw][i_frac - 1].back();
    //                        auto ca_logPa = ca_wa.logPa;
    //                        auto ca_par = ca_wa.parameter;
    //                        auto ca_logP = ca_wa.logP + ca_wa.logL;
    //                        auto ca_logL1 = logLikelihood(lik, ca_par, y[i_frac], x[i_frac]);
    //                        if (!(ca_logL1))
    //                            return error_message(ca_logL1.error()() + " push back new fraction at walker " +
    //                                                 std::to_string(iw) + "of fraction " +
    //                                                 std::to_string(i_frac) + " beta 0");
    //                        auto ca_logL = ca_logL1.value() - ca_logP + ca_logPa;
    //                        current.walkers[iw][i_frac][i_b] = {{ca_par, ca_logP, ca_logL},
    //                                                            ca_logPa};
    //                    }
    //            }
    
    //            for (std::size_t i_b = 1; i_b < size(current.beta[i_frac]); ++i_b) {
    //                for (std::size_t half = 0; half < 2; ++half)
    //                    for (std::size_t i = 0; i < n_walkers / 2; ++i) {
    //                        auto iw = i + half * n_walkers / 2;
    //                        std::swap(current.walkers[iw][i_frac][i_b], new_walkers[iw]);
    //                        std::swap(current.i_walkers[iw][i_frac][i_b], new_i_walkers[iw]);
    //                    }
    //            }
    //        }
    //        auto n_beta_current = size(current.beta[i_frac_old]);
    //        auto n_beta_current_final = size(final_beta[i_frac_old]);
    //        auto n_frac_final = size(final_beta);
    //        bool theres_room_for_one_more_beta=(n_beta_current < n_beta_current_final) ;
    //        bool this_is_the_final_fraction_in_a_more_than_2_betas_situation=n_frac_final == n_frac_old;
    
    //        if ( theres_room_for_one_more_beta||this_is_the_final_fraction_in_a_more_than_2_betas_situation
    //            ) {
    //            for (std::size_t half = 0; half < 2; ++half)
    //                for (std::size_t i = 0; i < n_walkers / 2; ++i) {
    //                    auto iw = i + half * n_walkers / 2;
    //                    current.walkers[iw][n_frac_old].push_back(new_walkers[iw]);
    //                    current.i_walkers[iw][n_frac_old].push_back(new_i_walkers[iw]);
    //                }
    //            current.beta[n_frac_old].push_back(
    //                final_beta[n_frac_old][n_beta_current]);
    //            return current;
    //        } else {
    //            for (std::size_t half = 0; half < 2; ++half)
    //                for (std::size_t i = 0; i < n_walkers / 2; ++i) {
    //                    auto iw = i + half * n_walkers / 2;
    //                    current.walkers[iw].push_back(by_beta<mcmc2<Parameters>>(2));
    //                    current.i_walkers[iw].push_back(by_beta<std::size_t>(2));
    
    //                    auto &ca_wa0 = current.walkers[iw][i_frac_old].back();
    //                    auto ca_logPa0 = ca_wa0.logPa;
    //                    auto ca_par0 = ca_wa0.parameter;
    //                    auto ca_logP0 = ca_wa0.logP + ca_wa0.logL;
    //                    auto ca_logL10 =
    //                        logLikelihood(lik, ca_par0, y[i_frac_old + 1], x[i_frac_old + 1]);
    //                    if (!(ca_logL10))
    //                        return error_message(ca_logL10.error()() + " push back new fraction at walker " +
    //                                             std::to_string(iw) + "of fraction " +
    //                                             std::to_string(i_frac_old + 1) + " beta 0");
    //                    auto ca_logL0 = ca_logL10.value() - ca_logP0 + ca_logPa0;
    //                    current.walkers[iw][i_frac_old + 1][0] = {
    //                                                              {ca_par0, ca_logP0, ca_logL0}, ca_logPa0};
    //                    current.i_walkers[iw][i_frac_old + 1][0] =
    //                        current.i_walkers[iw][i_frac_old].back();
    
    //                    auto &ca_wa1 = new_walkers[iw];
    //                    auto ca_logPa1 = ca_wa1.logPa;
    //                    auto ca_par1 = ca_wa1.parameter;
    //                    auto ca_logP1 = ca_wa1.logP + ca_wa1.logL;
    //                    auto ca_logL11 =
    //                        logLikelihood(lik, ca_par1, y[i_frac_old + 1], x[i_frac_old + 1]);
    //                    if (!(ca_logL11)) {
    //                        return error_message(ca_logL11.error()() + " push back new fraction at walker " +
    //                                             std::to_string(iw) + "of fraction " +
    //                                             std::to_string(i_frac_old + 1) + " beta 0");
    //                    } else {
    //                        auto ca_logL1 = ca_logL11.value() - ca_logP1 + ca_logPa1;
    //                        current.walkers[iw][i_frac_old + 1][1] = {
    //                                                                  {ca_par1, ca_logP1, ca_logL1}, ca_logPa1};
    //                        current.i_walkers[iw][i_frac_old + 1][1] = new_i_walkers[iw];
    //                    }
    //                }
    //            auto new_beta = by_beta<double>{0.0, 1.0};
    //            current.beta.push_back(new_beta);
    //            current.nsamples.push_back(size(y[n_frac_old]));
    //            return current;
    //        }
    //    }
}

template <class Observer, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

void thermo_cuevi_jump_mcmc(std::size_t iter, cuevi_mcmc<Parameters> &current,
                            Observer &obs, std::mt19937_64 &mt,
                            ensemble<std::mt19937_64> &mts, Prior const &prior,
                            Likelihood const &lik,
                            const by_fraction<DataType> &y,
                            const by_fraction<Variables> &x,
                            std::size_t thermo_jumps_every) {
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
                        
                        auto r = rdist[i](mts[i]);
                        double logA =
                            calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                      current.walkers[iw][i_fr][ib].logL,
                                      current.walkers[jw][i_fr][ib + 1].logL);
                        auto pJump = std::min(1.0, std::exp(logA));
                        observe_thermo_jump_mcmc(
                            obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                            current.walkers[jw][i_fr][ib + 1].parameter,
                            current.walkers[iw][i_fr][ib].logL,
                            current.walkers[jw][i_fr][ib + 1].logL,
                            -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                            pJump, r, pJump > r);
                        if (pJump > r) {
                            std::swap(current.walkers[iw][i_fr][ib],
                                      current.walkers[jw][i_fr][ib + 1]);
                            std::swap(current.i_walkers[iw][i_fr][ib],
                                      current.i_walkers[jw][i_fr][ib + 1]);
                        }
                    }
                }
            else {
                for (std::size_t i_fr = 0; i_fr < 1; ++i_fr) {
                    for (std::size_t ib = 0; ib < current.beta[i_fr].size() - 2; ++ib) {
                        
                        auto r = rdist[i](mts[i]);
                        double logA =
                            calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                      current.walkers[iw][i_fr][ib].logL,
                                      current.walkers[jw][i_fr][ib + 1].logL);
                        auto pJump = std::min(1.0, std::exp(logA));
                        observe_thermo_jump_mcmc(
                            obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                            current.walkers[jw][i_fr][ib + 1].parameter,
                            current.walkers[iw][i_fr][ib].logL,
                            current.walkers[jw][i_fr][ib + 1].logL,
                            -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                            pJump, r, pJump > r);
                        if (pJump > r) {
                            std::swap(current.walkers[iw][i_fr][ib],
                                      current.walkers[jw][i_fr][ib + 1]);
                            std::swap(current.i_walkers[iw][i_fr][ib],
                                      current.i_walkers[jw][i_fr][ib + 1]);
                        }
                    }
                    for (std::size_t ib = current.beta[i_fr].size() - 2;
                         ib < current.beta[i_fr].size() - 1; ++ib) {
                        
                        auto r = rdist[i](mts[i]);
                        double logA =
                            calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                      current.walkers[iw][i_fr][ib].logL,
                                      current.walkers[jw][i_fr][ib + 1].logL);
                        auto pJump = std::min(1.0, std::exp(logA));
                        observe_thermo_jump_mcmc(
                            obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                            current.walkers[jw][i_fr][ib + 1].parameter,
                            current.walkers[iw][i_fr][ib].logL,
                            current.walkers[jw][i_fr][ib + 1].logL,
                            -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                            pJump, r, pJump > r);
                        if (pJump > r) {
                            auto ca_par = current.walkers[iw][i_fr][ib].parameter;
                            auto ca_logL1 =
                                logLikelihood(lik, ca_par, y[i_fr + 1], x[i_fr + 1]);
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
                                    current.walkers[jw][i_fr + 1][0].logP = ca_logP + ca_logL0;
                                    current.walkers[jw][i_fr + 1][0].logL =
                                        ca_logL1.value() - ca_logL0 - ca_logP + ca_logPa;
                                }
                            }
                        }
                    }
                }
                for (std::size_t i_fr = 1; i_fr + 1 < current.beta.size(); ++i_fr) {
                    if (current.beta[i_fr].size() < 3) {
                        for (std::size_t ib = 0; ib + 1 < current.beta[i_fr].size(); ++ib) {
                            
                            auto r = rdist[i](mts[i]);
                            double logA =
                                calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                          current.walkers[iw][i_fr][ib].logL,
                                          current.walkers[jw][i_fr][ib + 1].logL);
                            auto pJump = std::min(1.0, std::exp(logA));
                            observe_thermo_jump_mcmc(
                                obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                                current.walkers[jw][i_fr][ib + 1].parameter,
                                current.walkers[iw][i_fr][ib].logL,
                                current.walkers[jw][i_fr][ib + 1].logL,
                                -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                                pJump, r, pJump > r);
                            if (pJump > r) {
                                auto ca_par_1 = current.walkers[iw][i_fr][ib].parameter;
                                auto ca_logL_11 =
                                    logLikelihood(lik, ca_par_1, y[i_fr + 1], x[i_fr + 1]);
                                auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                                auto ca_logL_00 = i_fr == 1
                                                      ? Maybe_error<double>{0.0}
                                                      : logLikelihood(lik, ca_par_0,
                                                                      y[i_fr - 2], x[i_fr - 2]);
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
                                    current.walkers[jw][i_fr + 1][0].logP = ca_logP_1 + ca_logL_1;
                                    current.walkers[jw][i_fr + 1][0].logL =
                                        ca_logL_11.value() - ca_logL_1 - ca_logP_1 + ca_logPa_1;
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
                    } else {
                        for (std::size_t ib = 0; ib < 1; ++ib) {
                            
                            auto r = rdist[i](mts[i]);
                            double logA =
                                calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                          current.walkers[iw][i_fr][ib].logL,
                                          current.walkers[jw][i_fr][ib + 1].logL);
                            auto pJump = std::min(1.0, std::exp(logA));
                            observe_thermo_jump_mcmc(
                                obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                                current.walkers[jw][i_fr][ib + 1].parameter,
                                current.walkers[iw][i_fr][ib].logL,
                                current.walkers[jw][i_fr][ib + 1].logL,
                                -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                                pJump, r, pJump > r);
                            if (pJump > r) {
                                
                                auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                                auto ca_logL_00 = i_fr == 1
                                                      ? Maybe_error<double>{0.0}
                                                      : logLikelihood(lik, ca_par_0,
                                                                      y[i_fr - 2], x[i_fr - 2]);
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
                        for (std::size_t ib = 1; ib + 2 < current.beta[i_fr].size(); ++ib) {
                            
                            auto r = rdist[i](mts[i]);
                            double logA =
                                calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                          current.walkers[iw][i_fr][ib].logL,
                                          current.walkers[jw][i_fr][ib + 1].logL);
                            auto pJump = std::min(1.0, std::exp(logA));
                            observe_thermo_jump_mcmc(
                                obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                                current.walkers[jw][i_fr][ib + 1].parameter,
                                current.walkers[iw][i_fr][ib].logL,
                                current.walkers[jw][i_fr][ib + 1].logL,
                                -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                                pJump, r, pJump > r);
                            if (pJump > r) {
                                std::swap(current.walkers[iw][i_fr][ib],
                                          current.walkers[jw][i_fr][ib + 1]);
                                std::swap(current.i_walkers[iw][i_fr][ib],
                                          current.i_walkers[jw][i_fr][ib + 1]);
                            }
                        }
                    }
                }
                for (std::size_t i_fr = std::max(1ul, current.beta.size() - 1);
                     i_fr < current.beta.size(); ++i_fr) {
                    for (std::size_t ib = 0; ib < 1; ++ib) {
                        
                        auto r = rdist[i](mts[i]);
                        double logA =
                            calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                      current.walkers[iw][i_fr][ib].logL,
                                      current.walkers[jw][i_fr][ib + 1].logL);
                        auto pJump = std::min(1.0, std::exp(logA));
                        observe_thermo_jump_mcmc(
                            obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                            current.walkers[jw][i_fr][ib + 1].parameter,
                            current.walkers[iw][i_fr][ib].logL,
                            current.walkers[jw][i_fr][ib + 1].logL,
                            -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                            pJump, r, pJump > r);
                        if (pJump > r) {
                            
                            auto ca_par_0 = current.walkers[jw][i_fr][ib + 1].parameter;
                            auto ca_logL_00 =
                                i_fr == 1
                                    ? Maybe_error<double>{0.0}
                                    : logLikelihood(lik, ca_par_0, y[i_fr - 2], x[i_fr - 2]);
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
                                                                                       mcmc<Parameters>{ca_par_0, cai_logP, cai_logL}, ca_logPa_0};
                            }
                        }
                    }
                    for (std::size_t ib = 1; ib + 1 < current.beta[i_fr].size(); ++ib) {
                        
                        auto r = rdist[i](mts[i]);
                        double logA =
                            calc_logA(current.beta[i_fr][ib], current.beta[i_fr][ib + 1],
                                      current.walkers[iw][i_fr][ib].logL,
                                      current.walkers[jw][i_fr][ib + 1].logL);
                        auto pJump = std::min(1.0, std::exp(logA));
                        observe_thermo_jump_mcmc(
                            obs[iw][ib], jw, current.walkers[iw][i_fr][ib].parameter,
                            current.walkers[jw][i_fr][ib + 1].parameter,
                            current.walkers[iw][i_fr][ib].logL,
                            current.walkers[jw][i_fr][ib + 1].logL,
                            -(current.beta[i_fr][ib] - current.beta[i_fr][ib + 1]), logA,
                            pJump, r, pJump > r);
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
    double min_fraction_;
    std::size_t thermo_jumps_every_;
    double n_points_per_decade_beta_;
    double n_points_per_decade_fraction_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;
    
public:
    cuevi_integration(Algorithm &&alg, Fractioner &&frac, Reporter &&rep,
                      std::size_t num_scouts_per_ensemble, double min_fraction,
                      std::size_t thermo_jumps_every,
                      double n_points_per_decade_beta,
                      double n_points_per_decade_fraction, double stops_at,
                      bool includes_zero, std::size_t initseed)
        : alg_{std::move(alg)}, frac_{std::move(frac)}, rep_{std::move(rep)},
        num_scouts_per_ensemble_{num_scouts_per_ensemble},
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
    auto &thermo_jumps_every() const { return thermo_jumps_every_; }
    auto &n_points_per_decade_beta() const { return n_points_per_decade_beta_; }
    auto &n_points_per_decade_fraction() const {
        return n_points_per_decade_fraction_;
    }
    auto &stops_at() const { return stops_at_; }
    auto &includes_zero() const { return includes_zero_; }
    auto &initseed() const { return initseed_; }
};

template <class Algorithm, class Prior, class Likelihood, class Variables,
         class DataType, class Fractioner, class Reporter,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_Algorithm_conditions<Algorithm, cuevi_mcmc<Parameters>> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<Likelihood, Parameters, Variables, DataType>)

auto evidence(cuevi_integration<Algorithm, Fractioner, Reporter> &&cue,
              Prior const &prior, Likelihood const &lik, const DataType &y,
              const Variables &x) {
    
    auto a = cue.algorithm();
    auto mt = init_mt(cue.initseed());
    auto n_walkers = cue.num_scouts_per_ensemble();
    auto mts = init_mts(mt, cue.num_scouts_per_ensemble() / 2);
    auto [ys, xs, beta_final] = cue.fractioner()(
        y, x, mt, size(prior) * cue.min_fraction(),
        cue.n_points_per_decade_beta(), cue.n_points_per_decade_fraction(),
        cue.stops_at(), cue.includes_zero());
    auto beta_init =
        by_beta<double>(beta_final[0].begin(), beta_final[0].begin() + 2);
    auto current = init_cuevi_mcmc(n_walkers, beta_init, mts, prior, lik, ys, xs);
    auto mcmc_run = checks_convergence(std::move(a), current);
    std::size_t iter = 0;
    auto &rep = cue.reporter();
    report_title(rep, current,prior,lik, ys, xs);
    report_model(rep, prior, lik, ys, xs, beta_final);
    while ((size(current.beta) < size(beta_final)) ||
           (size(current.beta.back()) < size(beta_final.back())) ||
           !mcmc_run.second) {
        while (!mcmc_run.second) {
            step_stretch_cuevi_mcmc(current, rep, mts, prior, lik, ys, xs);
            ++iter;
            thermo_cuevi_jump_mcmc(iter, current, rep, mt, mts, prior, lik, ys, xs,
                                   cue.thermo_jumps_every());
            report(iter, rep, current, prior, lik, ys, xs);
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
        if ((size(current.beta) < size(ys)) ||
            (size(current.beta.back()) < size(beta_final.back()))) {
            //   std::cerr<<"\n---walkers!!------------------------\n"<<current.walkers;
            
            auto is_current =
                push_back_new_fraction(current, mts, beta_final, prior, lik, ys, xs);
            while (!(is_current)) {
                std::cerr<<is_current.error()();
                step_stretch_cuevi_mcmc(current, rep, mts, prior, lik, ys, xs);
                ++iter;
                thermo_cuevi_jump_mcmc(iter, current, rep, mt, mts, prior, lik, ys, xs,
                                       cue.thermo_jumps_every());
                report(iter, rep, current, prior, lik, ys, xs);
                is_current = push_back_new_fraction(current, mts, beta_final, prior,
                                                    lik, ys, xs);
            }
            current = std::move(is_current.value());
            std::cerr<<"\niwalkers!!------------------------\n"<<current.i_walkers;
            std::cerr << "\n  nsamples=" << current.nsamples.back()
                      << "   beta_run=" << current.beta.back().back() << "\n";
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
    }
    
    return std::pair(mcmc_run, current);
}

template<class Parameters>
auto cuevi_by_convergence(std::string path, std::string filename,
                          std::size_t num_scouts_per_ensemble,
                          double min_fraction, std::size_t thermo_jumps_every,
                          std::size_t max_iter, double max_ratio,
                          double n_points_per_decade_beta,
                          double n_points_per_decade_fraction, double stops_at,
                          bool includes_zero, std::size_t initseed) {
    return cuevi_integration(
        checks_derivative_var_ratio<cuevi_mcmc, Parameters>(max_iter, max_ratio),
        fractioner{},
        save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
            path, filename, 100ul, 100ul, 100ul),
        num_scouts_per_ensemble, min_fraction, thermo_jumps_every,
        n_points_per_decade_beta, n_points_per_decade_fraction, stops_at,
        includes_zero, initseed);
}

#endif // CUEVI_H
