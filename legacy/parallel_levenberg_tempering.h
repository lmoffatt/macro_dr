#ifndef PARALLEL_LEVENBERG_TEMPERING_H
#define PARALLEL_LEVENBERG_TEMPERING_H

#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "bayesian_linear_regression.h"
#include "distributions.h"
#include "function_measure_verification_and_optimization.h"
#include "matrix.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "multivariate_normal_distribution.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "random_samplers.h"
#include "variables.h"
#include "variables_derivative.h"

struct levenberg_mcmc;

struct levenberg_Marquart_mcmc;

class levenberg_lambda_adaptive_distribution {
   private:
    std::vector<double> m_lambdas;
    std::vector<Trial_statistics> m_stat;

    std::map<double, std::size_t> m_cumulative;
    std::uniform_real_distribution<double> m_d;

    std::string m_algorithm;

    normal_distribution m_optimal_dist = normal_distribution(logit(0.23), 0.2);

    static auto build_lambdas(std::size_t n) {
        std::vector<double> out(n);
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = std::pow(2.0, i) - 1;
        }
        return out;
    }

    double weight(Trial_statistics const& p, double lambda) const {
        if (m_algorithm == "velocity")
            return p.rate().value() / (lambda + 1.0);
        else {
            return m_optimal_dist.P(logit(p.rate().value())).value();

        }  // beta distribution
        if constexpr (false) {
            double alfa = p.rate().value() * p.count();
            double beta = p.count() - alfa;
            double logBeta = std::lgamma(alfa) + std::lgamma(beta) - std::lgamma(alfa + beta);
            double logP = ((alfa - 1) * std::log(0.23) + (beta - 1) * std::log(1 - 0.23)) - logBeta;
            return std::exp(logP);
        }
    }

   public:
    levenberg_lambda_adaptive_distribution() = default;
    levenberg_lambda_adaptive_distribution(std::size_t n, const std::string& algorithm)
        : m_lambdas{build_lambdas(n)},
          m_stat{n, Trial_statistics(Trial_count(2), Success_count(1))},
          m_cumulative{},
          m_algorithm{algorithm},
          m_optimal_dist{normal_distribution(logit(0.23), 0.2)} {
        update_distribution();
    }

    auto size() const {
        return m_lambdas.size();
    }

    void update_distribution() {
        m_cumulative.clear();
        std::vector<double> sumcum_i(size());
        double sumcum = 0;
        for (std::size_t i = 0; i < size(); ++i) {
            if (m_stat[i].rate().valid())
                sumcum += weight(m_stat[i], m_lambdas[i]);
            sumcum_i[i] = sumcum;
        }
        for (std::size_t i = 0; i < size(); ++i) {
            m_cumulative[sumcum_i[i] / sumcum] = i;
        }
    }

    std::size_t operator()(mt_64i& mt) {
        auto max = m_cumulative.rbegin()->first;
        auto r = m_d(mt) * max;
        auto it = m_cumulative.lower_bound(r);
        return it->second;
    }

    double lambda(std::size_t i) const {
        return m_lambdas[i];
    }

    auto& stat(std::size_t i) {
        return m_stat[i];
    }
    auto& stat(std::size_t i) const {
        return m_stat[i];
    }

    auto& cumulative() const {
        return m_cumulative;
    }
};

struct levenberg_mcmc {
    std::size_t i_walker;
    var::Parameters_transformed m_x;
    dlogPs m_logP;
    dlogLs m_logL;

    logLs get_logLs() const {
        return {get<logL>(m_logL), get<elogL>(m_logL), get<vlogL>(m_logL)};
    }

    Maybe_error<bool> check_dimensions() const {
        auto npar = m_x.size();
        if (!((get<FIM>(m_logP)().nrows() == npar) && (get<FIM>(m_logP)().ncols() == npar) &&
              (get<Grad>(m_logP)().nrows() == npar)))
            return error_message("error in logP gradient or FIM dimensions");
        if (!((get<FIM>(m_logL)().nrows() == npar) && (get<FIM>(m_logL)().ncols() == npar) &&
              (get<Grad>(m_logL)().nrows() == npar)))
            return error_message("error in logL gradient or FIM dimensions");
        return true;
    }

    Maybe_error<bool> is_good() const {
        return ::is_good(m_logP) && ::is_good(m_logL) && m_x.is_good() && check_dimensions();
    }

    double logThP(double beta) const {
        return get<logL>(m_logP)() + beta * get<logL>(m_logL)();
    }

    Maybe_error<multivariate_normal_distribution<double, SymPosDefMatrix<double>>>
        ProposedDistribution(double beta, double lambda) const {
        auto thFim = (get<FIM>(m_logP)() + get<FIM>(m_logL)() * beta);

        auto FIml = SymPosDefMatrix<double>::I_sware_it_is_possitive(thFim + diag(thFim) * lambda);
        auto G = get<Grad>(m_logP)() + get<Grad>(m_logL)() * beta;

        auto& x = m_x;
        auto chol_inv = cholesky(FIml);
        if (!chol_inv)
            return chol_inv.error();
        auto chol = inv(chol_inv.value());
        if (!chol)
            return chol.error();
        auto thFiminv = XXT(tr(chol.value()));

        assert(([&FIml, &thFiminv]() {
            auto Maybe = test_equality(inv(FIml).value(), thFiminv, std::sqrt(eps));
            if (!Maybe) {
                std::cerr << Maybe;
                return false;
            }
            return true;
        }()));

        auto opt = x() + thFiminv * G;
        auto logDetCov = logdet(diag(chol.value()));
        if (!logDetCov)
            return logDetCov.error();
        return multivariate_normal_distribution(std::move(opt), std::move(thFiminv),
                                                std::move(chol.value()), std::move(FIml),
                                                logDetCov.value());
    }
};
class levenberg_Thermo_Jump_statistics
    : public var::Constant<levenberg_Thermo_Jump_statistics, Trial_statistics> {};
class levenberg_Step_statistics
    : public var::Constant<levenberg_Step_statistics, Trial_statistics> {};

class levenberg_log_errors
    : public var::Constant<levenberg_log_errors,
                           std::vector<error_log<var::Parameters_transformed>>> {
   public:
    void save_error(std::string error, var::Parameters_transformed&& val) {
        (*this)().emplace_back(std::move(error), std::move(val));
    }
};

struct levenberg_Marquart_mcmc {
    levenberg_mcmc m_data;
    levenberg_lambda_adaptive_distribution m_lambda;
    double m_beta;
    levenberg_Thermo_Jump_statistics m_thermo_stat;
    levenberg_Step_statistics m_step_stat;
    levenberg_log_errors m_error;
};

class levenberg_lambda_error_statistics
    : public var::Constant<levenberg_lambda_error_statistics, Trial_statistics> {};

class levenberg_likelihood_error_statistics
    : public var::Constant<levenberg_likelihood_error_statistics, Trial_statistics> {};

struct thermo_levenberg_mcmc {
    by_beta<levenberg_Marquart_mcmc> walkers;
    // by_beta<levenberg_lambda_error_statistics> lambda_stat;
    // by_beta<levenberg_likelihood_error_statistics> lik_stat;
    // by_beta<levenberg_Step_statistics> step_stat;
    // by_beta<levenberg_Thermo_Jump_statistics> thermo_stat;

    Maybe_error<bool> is_good(by_beta<double> const& beta) const {
        //      std::cerr<<"iside_current_is_good\n";
        if (walkers.size() != beta.size())
            return error_message("missing walkers " + std::to_string(walkers.size()) + "vs " +
                                 std::to_string(beta.size()));
        //  std::cerr<<"iside_current_is_good mid\n";

        Maybe_error<bool> out(true);
        for (std::size_t i = 0; i < walkers.size(); ++i) {
            //     std::cerr<<"iside_current_is_good walker"<<i<<"\n";
            out = out && walkers[i].m_data.is_good();
        }
        //std::cerr<<"iside_current_is_good end\n";
        return out;
    }

    void reset_statistics() {
        // for (auto &e : step_stat)
        //     e().reset();
        // for (auto &e : lik_stat)
        //     e().reset();
        // for (auto &e : thermo_stat)
        //     e().reset();
    }
};

inline Maybe_error<var::Parameters_transformed> levenberg_step(mt_64i& mt, levenberg_mcmc& current,
                                                               double beta, double lambda) {
    auto Maybe_Prop = current.ProposedDistribution(beta, lambda);
    if (!Maybe_Prop)
        return Maybe_Prop.error();
    auto Prop = std::move(Maybe_Prop.value());
    return current.m_x.create(Prop(mt));
}

inline Maybe_error<dlogPs> dlogPrior(var::Parameters_Normal_Distribution const& prior,
                                     var::Parameters_transformed const& ca) {
    auto r_logP = prior.logP(ca);
    if (!r_logP)
        return r_logP.error();
    return dlogPs{logL(r_logP.value()), Grad(prior.dlogP(ca)),
                  FIM(SymPosDefMatrix<double>(prior.FIM(ca)))};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

Maybe_error<levenberg_mcmc> calculate_next_levenberg_mcmc(FunctionTable& f, Parameters&& ca_par,
                                                          std::size_t i_walker, Prior const& prior,
                                                          Likelihood const& lik, const DataType& y,
                                                          const Variables& x) {
    auto Maybe_ca_logP = dlogPrior(prior, ca_par);
    if (!Maybe_ca_logP)
        return Maybe_ca_logP.error();

    auto Maybe_ca_logL = dlogLikelihood(f, lik, ca_par, y, x);
    if (!Maybe_ca_logL)
        return Maybe_ca_logL.error();
    return levenberg_mcmc{i_walker, std::move(ca_par), std::move(Maybe_ca_logP.value()),
                          std::move(Maybe_ca_logL.value())};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

Maybe_error<levenberg_mcmc> calculate_next_levenberg_mcmc(FunctionTable& f, Parameters&& ca_par,
                                                          std::size_t i_walker, Prior const& prior,
                                                          Likelihood const& lik, const DataType& y,
                                                          const Variables& x, double delta_par) {
    auto Maybe_ca_logP = dlogPrior(prior, ca_par);
    if (!Maybe_ca_logP)
        return Maybe_ca_logP.error();

    auto Maybe_ca_logL = diff_logLikelihood(f, lik, ca_par, y, x, delta_par);

    if (!Maybe_ca_logL)
        return Maybe_ca_logL.error();

    assert(([&ca_par, &f, &lik, &x, &y, &Maybe_ca_logL]() {
        auto ca_logL = logLikelihood(f, lik, ca_par.to_value(), y, x).value();
        if (std::abs(get<logL>(ca_logL)() - get<logL>(Maybe_ca_logL.value())()) >
            std::sqrt(eps) * std::abs(get<logL>(ca_logL)())) {
            std::cerr << get<logL>(ca_logL)() << " vs " << get<logL>(Maybe_ca_logL.value())();
            return false;
        } else
            return true;
    }()));
    return levenberg_mcmc{i_walker, std::move(ca_par), std::move(Maybe_ca_logP.value()),
                          std::move(Maybe_ca_logL.value())};
}

inline Maybe_error<double> levenberg_Acceptance(
    levenberg_mcmc const& current, levenberg_mcmc const& candidate,
    multivariate_normal_distribution<double, SymPosDefMatrix<double>> const& currentProp,
    multivariate_normal_distribution<double, SymPosDefMatrix<double>> const& candidateProp,
    double beta) {
    double logL_current = current.logThP(beta);

    double logL_candidate = candidate.logThP(beta);

    auto Maybe_logP_forward = currentProp.logP(candidate.m_x());
    auto Maybe_logP_backward = candidateProp.logP(current.m_x());
    if (!Maybe_logP_forward)
        return Maybe_logP_forward.error();
    if (!Maybe_logP_backward)
        return Maybe_logP_backward.error();
    double logP_forward = Maybe_logP_forward.value();
    double logP_backward = Maybe_logP_backward.value();

    double logA = (logL_candidate + logP_backward) - (logL_current + logP_forward);
    double A = std::min(1.0, std::exp(logA));

    return A;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
Maybe_error<bool> step_levenberg_thermo_mcmc_i(FunctionTable& f, mt_64i& mt,
                                               std::uniform_real_distribution<double>& rdist,
                                               levenberg_Marquart_mcmc& currentm,
                                               Prior const& prior, Likelihood const& lik,
                                               const DataType& y, const Variables& x, double beta) {
    auto& Lambda = currentm.m_lambda;
    auto i_lambda = Lambda(mt);
    auto lambda = Lambda.lambda(i_lambda);
    auto& current = currentm.m_data;
    auto Maybe_cu_Prop = current.ProposedDistribution(beta, lambda);
    if (!Maybe_cu_Prop) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());

        return Maybe_cu_Prop.error();
    }
    auto cu_Prop = std::move(Maybe_cu_Prop.value());
    auto ca_par = current.m_x.create(cu_Prop(mt));
    auto Maybe_ca_wa =
        calculate_next_levenberg_mcmc(f, std::move(ca_par), current.i_walker, prior, lik, y, x);

    if (!Maybe_ca_wa) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        currentm.m_error.save_error(Maybe_ca_wa.error()(), std::move(ca_par));
        return Maybe_ca_wa.error();
    }
    auto candidate = std::move(Maybe_ca_wa.value());
    auto Maybe_ca_Prop = candidate.ProposedDistribution(beta, lambda);
    if (!Maybe_ca_Prop) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return Maybe_ca_Prop.error();
    }
    auto ca_Prop = std::move(Maybe_ca_Prop.value());
    auto Maybe_A = levenberg_Acceptance(current, candidate, cu_Prop, ca_Prop, beta);
    if (!Maybe_A) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return Maybe_A.error();
    }
    double A = Maybe_A.value();
    double r = rdist(mt);
    if (r < A) {
        succeeds(Lambda.stat(i_lambda));
        succeeds(currentm.m_step_stat());
        current = std::move(candidate);
        return true;
    } else {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return false;
    }
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
Maybe_error<bool> step_levenberg_thermo_mcmc_i(FunctionTable& f, mt_64i& mt,
                                               std::uniform_real_distribution<double>& rdist,
                                               levenberg_Marquart_mcmc& currentm,
                                               Prior const& prior, Likelihood const& lik,
                                               const DataType& y, const Variables& x, double beta,
                                               double delta_par) {
    auto& Lambda = currentm.m_lambda;
    auto i_lambda = Lambda(mt);
    auto lambda = Lambda.lambda(i_lambda);
    auto& current = currentm.m_data;
    auto Maybe_cu_Prop = current.ProposedDistribution(beta, lambda);
    if (!Maybe_cu_Prop) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return Maybe_cu_Prop.error();
    }
    auto cu_Prop = std::move(Maybe_cu_Prop.value());
    auto ca_par = current.m_x.create(cu_Prop(mt));
    auto Maybe_ca_wa = calculate_next_levenberg_mcmc(f, std::move(ca_par), current.i_walker, prior,
                                                     lik, y, x, delta_par);

    if (!Maybe_ca_wa) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        currentm.m_error.save_error(Maybe_ca_wa.error()(), std::move(ca_par));
        return Maybe_ca_wa.error();
    }
    auto candidate = std::move(Maybe_ca_wa.value());
    auto Maybe_ca_Prop = candidate.ProposedDistribution(beta, lambda);
    if (!Maybe_ca_Prop) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return Maybe_ca_Prop.error();
    }
    auto ca_Prop = std::move(Maybe_ca_Prop.value());
    auto Maybe_A = levenberg_Acceptance(current, candidate, cu_Prop, ca_Prop, beta);
    if (!Maybe_A) {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return Maybe_A.error();
    }
    double A = Maybe_A.value();

    double r = rdist(mt);
    if (r < A) {
        succeeds(Lambda.stat(i_lambda));
        succeeds(currentm.m_step_stat());
        current = std::move(candidate);
        return true;
    } else {
        fails(Lambda.stat(i_lambda));
        fails(currentm.m_step_stat());
        return false;
    }
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
void step_levenberg_thermo_mcmc(FunctionTable& f, std::size_t& iter, thermo_levenberg_mcmc& current,
                                const by_beta<double>& beta, ensemble<mt_64i>& mt,
                                Prior const& prior, Likelihood const& lik, const DataType& y,
                                const Variables& x) {
    assert(beta.size() == current.walkers.size());
    auto n_beta = beta.size();
    if (beta.size() != current.walkers.size())
        return;
    auto n_par = current.walkers[0].m_data.m_x.size();
    std::size_t num_threads = omp_get_max_threads();

    auto ff = f.fork(num_threads);

    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(num_threads, uniform_real);
    std::size_t n_beta_f = std::ceil(n_beta / num_threads);

#pragma omp parallel for  // collapse(2)
    for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
        std::size_t ib0 = i_thread * n_beta_f;
        std::size_t ib1 = std::min(n_beta, (i_thread + 1) * n_beta_f);
        for (std::size_t ib = ib0; ib < ib1; ++ib) {
            // dur.record("begin_loop_walker", ib * 2);
            //   for (std::size_t iii=0; iii<4; ++iii)
            //       {
            auto Maybe_trans =
                step_levenberg_thermo_mcmc_i(ff[i_thread], mt[i_thread], rdist[i_thread],
                                             current.walkers[ib], prior, lik, y, x, beta[ib]);
            // if (!Maybe_trans)
            //    std::cerr << "error at iter=" << iter << "beta = " << beta[ib]
            //              << ": " << Maybe_trans.error()() << "\n";
        }
    }
    f += ff;
    ++iter;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
void step_levenberg_thermo_mcmc(FunctionTable& f, std::size_t& iter, thermo_levenberg_mcmc& current,
                                const by_beta<double>& beta, ensemble<mt_64i>& mt,
                                Prior const& prior, Likelihood const& lik, const DataType& y,
                                const Variables& x, double delta_par) {
    assert(beta.size() == current.walkers.size());
    auto n_beta = beta.size();
    auto n_par = current.walkers[0].m_data.m_x.size();
    std::size_t num_threads = omp_get_max_threads();

    auto ff = f.fork(num_threads);

    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(num_threads, uniform_real);
    std::size_t n_beta_f = std::ceil(n_beta / num_threads);

#pragma omp parallel for  // collapse(2)
    for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
        std::size_t ib0 = i_thread * n_beta_f;
        std::size_t ib1 = std::min(n_beta, (i_thread + 1) * n_beta_f);
        for (std::size_t ib = ib0; ib < ib1; ++ib) {
            // dur.record("begin_loop_walker", ib * 2);
            //   for (std::size_t iii=0; iii<4; ++iii)
            //       {
            auto Maybe_trans = step_levenberg_thermo_mcmc_i(ff[i_thread], mt[i_thread],
                                                            rdist[i_thread], current.walkers[ib],
                                                            prior, lik, y, x, beta[ib], delta_par);
            // if (!Maybe_trans)
            //    std::cerr << "error at iter=" << iter << "beta = " << beta[ib]
            //              << ": " << Maybe_trans.error()() << "\n";
        }
    }
    f += ff;
    ++iter;
}

inline double calc_logA(double betai, double betaj, dlogLs const& logLi, dlogLs const& logLj) {
    return -(betai - betaj) * (get<logL>(logLi)() - get<logL>(logLj)());
}

inline void thermo_levenberg_jump_mcmc(std::size_t iter, thermo_levenberg_mcmc& current,
                                       const by_beta<double>& beta, mt_64i& mt,
                                       ensemble<mt_64i>& mts, std::size_t thermo_jumps_every) {
    if (iter % (thermo_jumps_every) == 0) {
        std::uniform_real_distribution<double> uniform_real(0, 1);
        auto n_beta = beta.size();
        auto n_par = current.walkers[0].m_data.m_x.size();

        std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(),
                                                                  uniform_real);

        std::vector<by_beta<Thermo_Jump_statistics>> thermo_stat(
            omp_get_max_threads(), by_beta<Thermo_Jump_statistics>(n_beta - 1));

        std::size_t num_threads = omp_get_max_threads();
        std::size_t n_beta_f = 2ul * std::ceil((1.0 * n_beta) / num_threads / 2.0);

        for (std::size_t odd_or_even = 0; odd_or_even < 2; ++odd_or_even) {
#pragma omp parallel for  // collapse(2)
            for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
                std::size_t ib0 = i_thread * n_beta_f + odd_or_even;
                std::size_t ib1 = std::min(n_beta - 1, (i_thread + 1) * n_beta_f + odd_or_even);
                for (std::size_t ib = ib0; ib < ib1; ib += 2) {
                    auto r = rdist[i_thread](mts[i_thread]);
                    double logA =
                        calc_logA(beta[ib], beta[ib + 1], current.walkers[ib].m_data.m_logL,
                                  current.walkers[ib + 1].m_data.m_logL);
                    auto pJump = std::min(1.0, std::exp(logA));
                    if (pJump > r) {
                        std::swap(current.walkers[ib].m_data, current.walkers[ib + 1].m_data);
                        succeeds(current.walkers[ib].m_thermo_stat());
                    } else {
                        fails(current.walkers[ib].m_thermo_stat());
                    }
                }
            }
        }
    }
}

template <class Parameters>
class save_Levenberg_Lambdas {
   public:
    class separator : public std::string {
       public:
        using std::string::string;
        //   separator(std::string s):std::string(std::move(s)){}

        std::string operator()() const {
            return *this;
        }
        friend std::ostream& operator<<(std::ostream& os, const separator& sep) {
            return os << sep();
        }
        friend std::istream& operator>>(std::istream& is, const separator& sep) {
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
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;

    save_Levenberg_Lambdas(std::string const& path, std::size_t t_sampling_interval,
                           std::size_t t_max_number_of_values_per_iteration)
        : fname{path},
          f{std::ofstream(path + "__i_beta__i_lambda.csv")},
          sampling_interval{t_sampling_interval},
          max_number_of_values_per_iteration{t_max_number_of_values_per_iteration} {
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    }

    friend void report_model(save_Levenberg_Lambdas&, ...) {
    }
};

template <class Parameters>
class save_Levenberg_Errors {
   public:
    class separator : public std::string {
       public:
        using std::string::string;
        //   separator(std::string s):std::string(std::move(s)){}

        std::string operator()() const {
            return *this;
        }
        friend std::ostream& operator<<(std::ostream& os, const separator& sep) {
            return os << sep();
        }
        friend std::istream& operator>>(std::istream& is, const separator& sep) {
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
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;

    save_Levenberg_Errors(std::string const& path, std::size_t t_sampling_interval,
                          std::size_t t_max_number_of_values_per_iteration)
        : fname{path},
          f{std::ofstream(path + "__i_beta__i_errors.csv")},
          sampling_interval{t_sampling_interval},
          max_number_of_values_per_iteration{t_max_number_of_values_per_iteration} {
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    }

    friend void report_model(save_Levenberg_Errors&, ...) {
    }
};

class save_Levenberg_Lambdas_every
    : public var::Var<save_Levenberg_Lambdas_every, std::pair<std::size_t, std::size_t>> {};

class save_Levenberg_Errors_every
    : public var::Var<save_Levenberg_Errors_every, std::pair<std::size_t, std::size_t>> {};

class Saving_Levenberg_intervals
    : public var::Var<Saving_intervals,
                      var::Vector_Space<Save_Evidence_every, Save_Likelihood_every,
                                        Save_Parameter_every, save_Levenberg_Lambdas_every,
                                        save_Levenberg_Errors_every, Save_Predictions_every>> {};

template <class... saving, class... Ts, class Parameters>
void report_title(save_mcmc<Parameters, saving...>& f, thermo_levenberg_mcmc const& data,
                  const Ts&... ts) {
    (report_title(get<saving>(f.m_m), data, ts...), ...);
}

template <class Parameter>
void report_title(save_likelihood<Parameter>& s, thermo_levenberg_mcmc const&...) {
    s.f << "iter" << s.sep << "dur" << s.sep << "beta" << s.sep << "i_walker" << s.sep << "logP"
        << s.sep << " logL" << s.sep << " elogL" << s.sep << " vlogL" << s.sep << "plog_Evidence"
        << s.sep << "eplog_Evidence" << s.sep << "vplog_Evidence" << s.sep << "log_Evidence"
        << s.sep << "elog_Evidence" << s.sep << "vlog_Evidence" << s.sep << "step_count" << s.sep
        << "step_rate" << s.sep << "thermo_count" << s.sep << "thermo_rate"
        << "\n";
}

template <class Parameter, class FunctionTable, class Duration>
void report(FunctionTable&, std::size_t iter, const Duration& dur, save_likelihood<Parameter>& s,
            thermo_levenberg_mcmc& current, by_beta<double> t_beta, ...) {
    std::size_t num_values = 16;
    auto num_beta = t_beta.size();
    std::size_t point_size = num_values * num_beta;
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if (iter % sampling_interval == 0) {
        //std::cerr<<"report save_Likelihood\n";
        logLs t_logL = {};
        double beta = 0;
        logLs log_Evidence =
            var::Vector_Space<logL, elogL, vlogL>(logL(0.0), elogL(0.0), vlogL(0.0));
        for (std::size_t i_beta0 = t_beta.size(); i_beta0 > 0; --i_beta0) {
            auto logL0 = t_logL;
            double beta0 = beta;
            auto i_beta = i_beta0 - 1;
            t_logL = current.walkers[i_beta].m_data.get_logLs();
            beta = t_beta[i_beta];
            auto plog_Evidence = (beta - beta0) * (logL0 + t_logL) / 2.0;
            log_Evidence = log_Evidence + plog_Evidence;
            s.f << iter << s.sep << dur.count() << s.sep << beta << s.sep
                << current.walkers[i_beta].m_data.i_walker << s.sep
                << get<logL>(current.walkers[i_beta].m_data.m_logP) << t_logL.sep(s.sep)
                << plog_Evidence.sep(s.sep) << log_Evidence.sep(s.sep) << s.sep;
            s.f << current.walkers[i_beta].m_step_stat().count() << s.sep
                << current.walkers[i_beta].m_step_stat().rate() << s.sep
                << current.walkers[i_beta].m_thermo_stat().count() << s.sep
                << current.walkers[i_beta].m_thermo_stat().rate() << "\n";
            if (iter % (sampling_interval * 10) == 0) {
                current.walkers[i_beta].m_step_stat().reset();
                current.walkers[i_beta].m_thermo_stat().reset();
            }
        }
        //std::cerr<<"report save_Likelihood end\n";
    }
}

inline std::size_t num_Parameters(thermo_levenberg_mcmc const& x) {
    return x.walkers[0].m_data.m_x.size();
}

inline std::size_t num_lambdas(thermo_levenberg_mcmc const& x) {
    return x.walkers[0].m_lambda.size();
}

template <class Parameter>
void report_title(save_Parameter<Parameter>& s, thermo_levenberg_mcmc const&, ...) {
    s.f << "iter" << s.sep << "dur" << s.sep << "beta" << s.sep << "i_walker" << s.sep << "i_par"
        << s.sep << "par_value" << s.sep << "gradient"
        << "\n";
}

template <class Parameter, class FunctionTable, class Duration>
void report(FunctionTable& f, std::size_t iter, const Duration& dur, save_Parameter<Parameter>& s,
            thermo_levenberg_mcmc const& current, by_beta<double> t_beta, ...) {
    std::size_t num_values = 4;
    auto num_beta = current.walkers.size();
    std::size_t point_size = num_values * num_beta * num_Parameters(current);
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if (iter % sampling_interval == 0) {
        //std::cerr<<"report save_Parameter\n";

        for (std::size_t i_beta = 0; i_beta < current.walkers.size(); ++i_beta)
            for (std::size_t i_par = 0; i_par < current.walkers[i_beta].m_data.m_x.size(); ++i_par)

                s.f << iter << s.sep << dur.count() << s.sep << t_beta[i_beta] << s.sep
                    << current.walkers[i_beta].m_data.i_walker << s.sep << i_par << s.sep
                    << current.walkers[i_beta].m_data.m_x[i_par] << s.sep
                    << get<Grad>(current.walkers[i_beta].m_data.m_logL)()[i_par] << "\n";
        //std::cerr<<"report save_Parameter end\n";
    }
}
template <class Parameter>
void report_title(save_Levenberg_Lambdas<Parameter>& s, thermo_levenberg_mcmc const&, ...) {
    s.f << "iter" << s.sep << "dur" << s.sep << "beta" << s.sep << "lambda" << s.sep << "count"
        << s.sep << "rate" << s.sep << "probability" << s.sep << "cumulative"
        << "\n";
}

template <class Parameter>
void report_title(save_Levenberg_Errors<Parameter>& s, thermo_levenberg_mcmc const& m, ...) {
    s.f << "iter" << s.sep << "dur" << s.sep << "beta" << s.sep << "i_error" << s.sep << "error";
    for (std::size_t ipar = 0; ipar < m.walkers[0].m_data.m_x.size(); ++ipar)
        s.f << s.sep << m.walkers[0].m_data.m_x.parameters().transformed_names()[ipar];
    s.f << "\n";
}

template <class Parameter, class FunctionTable, class Duration>
void report(FunctionTable&, std::size_t iter, const Duration& dur,
            save_Levenberg_Lambdas<Parameter>& s, thermo_levenberg_mcmc& current,
            by_beta<double> t_beta, ...) {
    std::size_t num_values = 8;
    auto num_beta = size(t_beta);
    auto v_num_landa = num_lambdas(current);
    std::size_t point_size = num_values * num_beta * v_num_landa;
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if (iter % sampling_interval == 0) {
        //std::cerr<<"report save_Levenberg_Lambdas\n";
        for (std::size_t i_beta = 0; i_beta < t_beta.size(); ++i_beta) {
            auto it_cum = current.walkers[i_beta].m_lambda.cumulative().begin();
            double p = 0;
            for (std::size_t i_lambda = 0; i_lambda < num_lambdas(current); ++i_lambda) {
                auto p1 = std::max(it_cum->first, p);

                s.f << iter << s.sep << dur.count() << s.sep << t_beta[i_beta] << s.sep
                    << current.walkers[i_beta].m_lambda.lambda(i_lambda) << s.sep
                    << current.walkers[i_beta].m_lambda.stat(i_lambda).count() << s.sep
                    << current.walkers[i_beta].m_lambda.stat(i_lambda).rate() << s.sep << p1 - p
                    << s.sep << it_cum->first << "\n";
                p = p1;
                ++it_cum;
            }
            current.walkers[i_beta].m_lambda.update_distribution();
        }
        //std::cerr<<"report save_Levenberg_Lambdas end\n";
    }
}

template <class Parameter, class FunctionTable, class Duration>
void report(FunctionTable&, std::size_t iter, const Duration& dur,
            save_Levenberg_Errors<Parameter>& s, thermo_levenberg_mcmc& current,
            by_beta<double> t_beta, ...) {
    std::size_t num_values = 8;
    auto num_beta = size(t_beta);
    auto v_num_error = current.walkers[0].m_error().size();
    std::size_t point_size = num_values * num_beta * v_num_error;
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if (iter % sampling_interval == 0) {
        for (std::size_t i_beta = 0; i_beta < t_beta.size(); ++i_beta) {
            for (std::size_t i_error = 0; i_error < current.walkers[i_beta].m_error().size();
                 ++i_error) {
                s.f << iter << s.sep << dur.count() << s.sep << t_beta[i_beta] << s.sep << i_error
                    << s.sep << current.walkers[i_beta].m_error()[i_error].error;
                for (std::size_t ipar = 0;
                     ipar < current.walkers[i_beta].m_error()[i_error].parameter.size(); ++ipar) {
                    s.f << s.sep << current.walkers[i_beta].m_error()[i_error].parameter[ipar];
                }
                s.f << "\n";
            }
            current.walkers[i_beta].m_error().clear();
        }
    }
}

template <class FunctionTable, class Duration, class Parameters, class... saving, class... T>
void report_all(FunctionTable& f, std::size_t iter, const Duration& dur,
                save_mcmc<Parameters, saving...>& s, thermo_levenberg_mcmc& data, T const&... ts) {
    //std::cerr<<"in report_all\n";
    (report(f, iter, dur, get<saving>(s.m_m), data, ts...), ..., 1);
    //std::cerr<<"after report_all\n";
}

template <class Algorithm, class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> )
class thermodynamic_levenberg_integration {
    Algorithm alg_;
    Reporter rep_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t thermo_jumps_every_;
    std::size_t beta_size_;
    std::size_t beta_upper_size_;
    std::size_t beta_medium_size_;
    std::size_t m_n_lambdas;
    std::string m_lambda_adaptive_algorithm;
    double beta_upper_value_;
    double beta_medium_value_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;
    double diff_par_;

   public:
    thermodynamic_levenberg_integration(Algorithm&& alg, Reporter&& rep,
                                        std::size_t num_scouts_per_ensemble,
                                        std::size_t thermo_jumps_every, std::size_t beta_size,
                                        std::size_t beta_upper_size, std::size_t beta_medium_size,
                                        double beta_upper_value, double beta_medium_value,
                                        std::size_t n_lambdas,
                                        std::string lambda_adaptive_algorithm, double stops_at,
                                        bool includes_zero, std::size_t initseed, double dp)
        : alg_{std::move(alg)},
          rep_{std::move(rep)},
          num_scouts_per_ensemble_{num_scouts_per_ensemble},
          thermo_jumps_every_{thermo_jumps_every},
          beta_size_{beta_size},
          beta_upper_size_{beta_upper_size},
          beta_medium_size_{beta_medium_size},
          beta_upper_value_{beta_upper_value},
          beta_medium_value_{beta_medium_value},
          m_n_lambdas{n_lambdas},
          m_lambda_adaptive_algorithm{lambda_adaptive_algorithm},
          stops_at_{stops_at},
          includes_zero_{includes_zero},
          initseed_{initseed},
          diff_par_{dp} {
    }

    auto& algorithm() const {
        return alg_;
    }
    auto& reporter() {
        return rep_;
    }
    auto& num_scouts_per_ensemble() const {
        return num_scouts_per_ensemble_;
    }

    auto& thermo_jumps_every() const {
        return thermo_jumps_every_;
    }
    auto& beta_size() const {
        return beta_size_;
    }
    auto& beta_upper_size() const {
        return beta_upper_size_;
    }
    auto& beta_medium_size() const {
        return beta_medium_size_;
    }
    auto& beta_upper_value() const {
        return beta_upper_value_;
    }
    auto& beta_medium_value() const {
        return beta_medium_value_;
    }
    auto& n_lambdas() const {
        return m_n_lambdas;
    }
    auto& lambda_adaptive_algorithm() const {
        return m_lambda_adaptive_algorithm;
    }
    auto& stops_at() const {
        return stops_at_;
    }
    auto& includes_zero() const {
        return includes_zero_;
    }
    auto& initseed() const {
        return initseed_;
    }
    auto& delta_par() const {
        return diff_par_;
    }
};

template <class FunctionTable, class Algorithm, class Prior, class Likelihood, class Variables,
          class DataType, class Reporter, class mcmc, class timepoint>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_levenberg_evidence_loop(
    FunctionTable& f, thermodynamic_levenberg_integration<Algorithm, Reporter>&& therm,
    Prior const& prior, Likelihood const& lik, const DataType& y, const Variables& x, mcmc mcmc_run,
    std::size_t iter, thermo_levenberg_mcmc& current, Reporter& rep, const by_beta<double> beta_run,
    mt_64i& mt, std::vector<mt_64i>& mts, const timepoint& start) {
    var::Event_Timing<200> even_dur(start);
    std::ofstream event_file(f.file() + "_event_timing.csv");

    auto delta_par = therm.delta_par();

    while (!mcmc_run.second) {
        //std::cerr << "main_loop_start"<< "\n";
        even_dur.record("main_loop_start");

        if (delta_par > 0)
            step_levenberg_thermo_mcmc(f, iter, current, beta_run, mts, prior, lik, y, x,
                                       delta_par);
        else
            step_levenberg_thermo_mcmc(f, iter, current, beta_run, mts, prior, lik, y, x);
        // auto check_current=current.is_good(beta_run);
        // if (!check_current)
        //std::cerr<<check_current.error()();
        even_dur.record("befor_thermo_jump");
        // std::cerr << "befor_thermo_jump"
        //          << "\n";

        thermo_levenberg_jump_mcmc(iter, current, beta_run, mt, mts, therm.thermo_jumps_every());
        even_dur.record("after_thermo_jump");
        //std::cerr << "after_thermo_jump" << "\n";
        //  check_current=current.is_good(beta_run);
        //std::cerr << "current.is_good"<< "\n";
        // if (!check_current)
        //     std::cerr<<check_current.error()();

        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start);
        report_all(f, iter, dur, rep, current, beta_run, prior, lik, y, x, mts, mcmc_run.first);
        even_dur.record("after_report_all");
        // std::cerr << "after_report_all"<< "\n";

        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        even_dur.record("after_checks_convergence");
        // std::cerr << "after_checks_convergence"<< "\n";
        if (iter == 1)
            even_dur.report_title(event_file);
        even_dur.report_iter(event_file, iter);
        if (iter % 10 == 0)
            event_file.flush();
        //      std::cerr << "end of loop"                  << "\n";
    }

    return std::pair(std::move(mcmc_run.first), current);
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&&
//   is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_levenberg_mcmc(FunctionTable& f, mt_64i& mt, std::size_t i_walker, Prior const& pr,
                         const Lik& lik, const DataType& y, const Variables& x) {
    auto& priorsampler = pr;
    auto par = sample(mt, priorsampler);
    auto logP = dlogPrior(pr, par);
    auto t_logLs = dlogLikelihood(f, lik, par, y, x);
    while (!(logP) || !(t_logLs)) {
        par = sample(mt, priorsampler);
        logP = dlogPrior(pr, par);
        t_logLs = dlogLikelihood(f, lik, par, y, x);
    }
    return levenberg_mcmc{i_walker, std::move(par), logP.value(), t_logLs.value()};
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&&
//   is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_levenberg_mcmc(FunctionTable& f, mt_64i& mt, std::size_t i_walker, Prior const& pr,
                         const Lik& lik, const DataType& y, const Variables& x, double delta_par) {
    auto& priorsampler = pr;
    auto par = sample(mt, priorsampler);
    auto logP = dlogPrior(pr, par);
    auto t_logLs = diff_logLikelihood(f, lik, par, y, x, delta_par);
    while (!(logP) || !(t_logLs)) {
        par = sample(mt, priorsampler);
        logP = dlogPrior(pr, par);
        t_logLs = diff_logLikelihood(f, lik, par, y, x, delta_par);
    }
    return levenberg_mcmc{i_walker, std::move(par), logP.value(), t_logLs.value()};
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&&
//   is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_levenberg_marquardt_mcmc(FunctionTable& f, mt_64i& mt, std::size_t i_walker,
                                   Prior const& pr, const Lik& lik, const DataType& y,
                                   const Variables& x, double beta, std::size_t n_lambdas,
                                   std::string lambda_adaptive_algorithm) {
    return levenberg_Marquart_mcmc{
        init_levenberg_mcmc(f, mt, i_walker, pr, lik, y, x),
        levenberg_lambda_adaptive_distribution(n_lambdas, lambda_adaptive_algorithm), beta};
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&&
//   is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_levenberg_marquardt_mcmc(FunctionTable& f, mt_64i& mt, std::size_t i_walker,
                                   Prior const& pr, const Lik& lik, const DataType& y,
                                   const Variables& x, double beta, std::size_t n_lambdas,
                                   std::string lambda_adaptive_algorithm, double delta_par) {
    return levenberg_Marquart_mcmc{
        init_levenberg_mcmc(f, mt, i_walker, pr, lik, y, x, delta_par),
        levenberg_lambda_adaptive_distribution(n_lambdas, lambda_adaptive_algorithm), beta};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
auto init_levenberg_thermo_mcmc(FunctionTable& f, by_beta<double> const& beta,
                                std::size_t n_lambdas, const std::string& lambda_adaptive_algorithm,
                                ensemble<mt_64i>& mt, Prior const& prior, Likelihood const& lik,
                                const DataType& y, const Variables& x, double delta_par) {
    by_beta<levenberg_Marquart_mcmc> walker(beta.size());
    by_beta<emcee_Step_statistics> emcee_stat(beta.size());
    by_beta<Thermo_Jump_statistics> thermo_stat(beta.size() - 1);
    auto ff = f.fork(omp_get_max_threads());

#pragma omp parallel for  // collapse(2)
    for (std::size_t ib = 0; ib < beta.size(); ++ib) {
        auto i_th = omp_get_thread_num();
        if (delta_par > 0)
            walker[ib] =
                init_levenberg_marquardt_mcmc(ff[i_th], mt[i_th], ib, prior, lik, y, x, beta[ib],
                                              n_lambdas, lambda_adaptive_algorithm, delta_par);
        else
            walker[ib] =
                init_levenberg_marquardt_mcmc(ff[i_th], mt[i_th], ib, prior, lik, y, x, beta[ib],
                                              n_lambdas, lambda_adaptive_algorithm);
    }

    f += ff;
    return thermo_levenberg_mcmc{walker};
}

template <class FunctionTable, class Algorithm, class Prior, class Likelihood, class Variables,
          class DataType, class Reporter>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_levenberg_evidence(FunctionTable& f,
                               thermodynamic_levenberg_integration<Algorithm, Reporter>&& therm,
                               Prior const& prior, Likelihood const& lik, const DataType& y,
                               const Variables& x) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto mts = init_mts(mt, omp_get_max_threads());
    auto beta =
        new_get_beta_list(therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
                          therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
                          therm.includes_zero());
    auto n_lambdas = therm.n_lambdas();
    auto lambda_adaptive_algorithm = therm.lambda_adaptive_algorithm();
    auto it_beta_run_begin = beta.rend() - beta.size();
    auto it_beta_run_end = beta.rend();
    auto beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);

    auto current = init_levenberg_thermo_mcmc(f, beta_run, n_lambdas, lambda_adaptive_algorithm,
                                              mts, prior, lik, y, x, therm.delta_par());
    // auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    const auto start = std::chrono::high_resolution_clock::now();
    auto& rep = therm.reporter();
    report_title(rep, current, lik, y, x, beta);
    // report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta);

    return thermo_levenberg_evidence_loop(
        f, std::forward<thermodynamic_levenberg_integration<Algorithm, Reporter>>(therm), prior,
        lik, y, x, mcmc_run, iter, current, rep, beta_run, mt, mts, start);
}

#endif  // PARALLEL_LEVENBERG_TEMPERING_H
