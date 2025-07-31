#ifndef PARALLEL_TEMPERING_LINEAR_REGRESSION_H
#define PARALLEL_TEMPERING_LINEAR_REGRESSION_H
#include <chrono>
#include <cstddef>
#include <fstream>
#include <type_traits>
#include <utility>
#include <vector>

#include "bayesian_linear_regression.h"
#include "function_measure_verification_and_optimization.h"
#include "mcmc.h"
#include "multivariate_normal_distribution.h"
#include "parallel_tempering.h"

class save_Iter {
   public:
    std::string sep = ",";
    std::string fname;
    std::ofstream f;
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;
    std::chrono::duration<double> m_prev_dur;
    save_Iter(std::string const& path, std::size_t t_sampling_interval,
              std::size_t t_max_number_of_values_per_iteration)
        : fname{path},
          f{std::ofstream(path + "_iter_time.csv")},
          sampling_interval{t_sampling_interval},
          max_number_of_values_per_iteration{t_max_number_of_values_per_iteration} {
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    }

    template <class mcmc>
    friend void report_title(save_Iter& s, mcmc const&, ...) {
        s.f << "iter" << s.sep << "iter_time" << s.sep << "iter_dur"
            << "\n";
    }

    template <class FunctionTable>
    friend void report(FunctionTable&&, std::size_t iter, const std::chrono::duration<double> dur,
                       save_Iter& s, ...) {
        std::size_t point_size = 2;
        std::size_t sampling_interval =
            std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

        if (iter % sampling_interval == 0) {
            s.f << iter << s.sep << dur.count() << s.sep << (dur - s.m_prev_dur).count() << "\n";
            s.m_prev_dur = dur;
        }
        if (iter % 10 == 0)
            s.f.flush();
    }

    friend void report_model(save_Iter&, ...) {
    }
};

class save_Evidence {
   public:
    std::string sep = ",";
    std::string fname;
    std::ofstream f;
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;
    save_Evidence(std::string const& path, std::size_t t_sampling_interval,
                  std::size_t t_max_number_of_values_per_iteration)
        : fname{path},
          f{std::ofstream(path + "__i_iter.csv")},
          sampling_interval{t_sampling_interval},
          max_number_of_values_per_iteration{t_max_number_of_values_per_iteration} {
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    }

    template <class Parameters>
    friend void report_title(save_Evidence& s, thermo_mcmc<Parameters> const&, ...) {
        s.f << "iter" << s.sep << "iter_time" << s.sep << "i_beta" << s.sep << "num_beta" << s.sep
            << "beta" << s.sep << "sample_size"

            << s.sep << "effective_sample_size"

            << s.sep << "meanPrior" << s.sep << "logL" << s.sep << "elogL" << s.sep << "vlogL"

            << s.sep << "var_logL" << s.sep << "var_elogL" << s.sep << "var_vlogL"

            << s.sep << "mean_logL_across" << s.sep << "var_logL_across" << s.sep
            << "count_logL_across"

            << s.sep << "var_logL_within"

            << s.sep << "plog_Evidence" << s.sep << "eplog_Evidence" << s.sep << "vplog_Evidence"

            << s.sep << "log_Evidence" << s.sep << "elog_Evidence" << s.sep << "vlog_Evidence"

            << s.sep << "mean_logL" << s.sep << "var_logL" << s.sep << "count_logL"

            << s.sep << "mean_plog_Evidence" << s.sep << "var_plog_Evidence" << s.sep
            << "count_plog_Evidence" << s.sep << "mean_log_Evidence" << s.sep << "var_log_Evidence"
            << s.sep << "count_log_Evidence"

            << s.sep << "deltaEvidence_variance" << s.sep << "Acceptance_variance" << s.sep
            << "emcee_stat_count" << s.sep << "emcee_stat_rate" << s.sep << "thermo_jump_stat_count"
            << s.sep << "thermo_jump_rate" << s.sep << "deltaBeta_deltalogL"
            << "\n";
    }

    template <class FunctionTable, class Duration, class Parameters>
    friend void report(FunctionTable& f, std::size_t iter, const Duration& dur, save_Evidence& s,
                       thermo_mcmc<Parameters> const& data, ...) {
        std::size_t num_values = 32;
        std::size_t point_size = num_values * num_betas(data);
        std::size_t sampling_interval =
            std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
        if ((iter > 0) && (data.num_samples() > 0) && (iter % sampling_interval == 0)) {
            auto across = calculate_across_sta(data.walkers_sta);
            auto within = calculate_within_sta(data.walkers_sta);
            auto eff_size = calculate_effective_sample_size(across, within);

            auto deltaEvidence_variance = calculate_delta_Evidence_variance(data, data.beta);
            auto Acceptance_variance = calculate_Acceptance_variance(deltaEvidence_variance, data);

            auto dBdL = calculate_deltaBeta_deltaL(data);

            auto meanLik = mean_logL(data);
            auto meanPrior = mean_logP(data);

            auto varLik = var_logL(data, meanLik);
            logLs r_logL = {};
            double beta = 0;
            logL_statistics m_logL = {};
            logLs log_Evidence = {};
            logL_statistics m_log_Evidence = {};
            for (std::size_t i_beta = 0; i_beta < data.get_Beta().size(); ++i_beta) {
                auto logL0 = r_logL;
                auto m_logL0 = m_logL;
                double beta0 = beta;
                beta = data.beta[i_beta];
                r_logL = meanLik[i_beta];
                m_logL = across[i_beta];
                auto plog_Evidence = (beta - beta0) / 2.0 * (logL0 + r_logL);
                auto m_plog_Evidence = (beta - beta0) / 2.0 * (m_logL0 + m_logL);
                if (beta0 > 0) {
                    log_Evidence = log_Evidence + plog_Evidence;

                    m_log_Evidence = m_log_Evidence + m_plog_Evidence;
                }
                auto emcee_count = data.emcee_stat[i_beta]().count();
                auto emcee_rate = data.emcee_stat[i_beta]().rate();
                auto thermo_count = data.thermo_stat[std::max(1ul, i_beta) - 1ul]().count();
                auto thermo_rate = data.thermo_stat[std::max(1ul, i_beta) - 1ul]().rate();

                s.f << iter << s.sep << dur.count() << s.sep << i_beta << s.sep << data.beta.size()
                    << s.sep << beta << s.sep << calculate_sample_size(data.walkers_sta, i_beta)
                    << s.sep << eff_size[i_beta]

                    << s.sep << meanPrior[i_beta] << r_logL.sep(s.sep) << varLik[i_beta].sep(s.sep)
                    << across[i_beta]()().sep(s.sep) << s.sep << within[i_beta]()()
                    << plog_Evidence.sep(s.sep) << log_Evidence.sep(s.sep) << m_logL()().sep(s.sep)
                    << m_plog_Evidence()().sep(s.sep) << m_log_Evidence()().sep(s.sep) << s.sep
                    << deltaEvidence_variance[std::max(1ul, i_beta) - 1] << s.sep
                    << Acceptance_variance[std::max(1ul, i_beta) - 1] << s.sep << emcee_count
                    << s.sep << emcee_rate << s.sep << thermo_count << s.sep << thermo_rate << s.sep
                    << dBdL[std::max(1ul, i_beta) - 1] << "\n";
            }
        }
    }

    template <class FunctionTable, class Parameters>
    friend void report_old(FunctionTable&& f, std::size_t iter, save_Evidence& s,
                           thermo_mcmc<Parameters> const& data, ...) {
        std::size_t num_values = 10;
        std::size_t point_size = num_values * num_betas(data);
        if ((iter > 0) && (data.num_samples() > 0) &&
            (iter %
                 std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration) ==
             0)) {
            auto meanLik = mean_logL(data);
            auto meanPrior = mean_logP(data);

            auto varLik = var_logL(data, meanLik);
            if (data.beta[0] == 1) {
                auto Evidence2 = calculate_Evidence(data.beta, meanLik, varLik);
                auto Evidence1 = calculate_Evidence(data.beta, meanLik);
                for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
                    s.f << num_betas(data) << s.sep << iter << s.sep << data.beta[i_beta] << s.sep
                        << meanPrior[i_beta] << s.sep << meanLik[i_beta] << s.sep << varLik[i_beta]
                        << s.sep << Evidence1 << s.sep << Evidence2 << "\n";
            }
        }
    }
    template <class Prior, class Likelihood, class Variables, class DataType>
        requires requires(Prior const& prior, Likelihood const& lik, const DataType& y,
                          const Variables& x) {
            { evidence(conjugate{}, prior, lik, y, x) };
        }
    friend void report_model(save_Evidence& s, Prior const& prior, Likelihood const& lik,
                             const DataType& y, const Variables& x, by_beta<double> const& beta0) {
        auto expected_Evidence = evidence(conjugate{}, prior, lik, y, x);

        auto meanLik = mean_logLik(conjugate{}, prior, lik, y, x, beta0);
        if (is_valid(meanLik) && is_valid(expected_Evidence))
            for (std::size_t i_beta = 0; i_beta < size(beta0); ++i_beta)
                s.f << 0 << s.sep << 0 << s.sep << beta0[i_beta] << s.sep << 0 << s.sep
                    << meanLik.value()[i_beta] << s.sep << 0 << s.sep << expected_Evidence << s.sep
                    << expected_Evidence << "\n";
    }

    template <class Prior, class Likelihood, class Variables, class DataType>
    friend void report_model(save_Evidence& s, Prior const& prior, Likelihood const& lik,
                             const DataType& y, const Variables& x, by_beta<double> const& beta0) {
        for (std::size_t i_beta = 0; i_beta < size(beta0); ++i_beta)
            s.f << 0 << s.sep << 0 << s.sep << beta0[i_beta] << "\n";
    }
};

inline void report_model(save_Evidence&, ...) {
}
template <class Parameters, class Algorithm, class Reporter>
    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>>)
class thermo {
    Algorithm alg;
    Reporter rep;
    std::size_t num_scouts_per_ensemble;
    std::size_t max_num_simultaneous_temperatures;
    std::size_t thermo_jumps_every;
    double n_points_per_decade;
    double stops_at;
    bool includes_zero;
    std::size_t initseed;

    template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType>
        requires(is_prior<Prior, Parameters, Variables, DataType> &&
                 is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

    auto logEvidence(FunctionTable&& f, const Prior prior, const Likelihood& lik, const DataType& y,
                     const Variables& x) {
        auto a = alg;
        auto mt = init_mt(initseed);
        auto n_walkers = num_scouts_per_ensemble;
        auto mts = init_mts(mt, num_scouts_per_ensemble / 2);
        auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);

        auto beta_run = by_beta<double>(beta.rend() - 2, beta.rend());

        auto current = init_thermo_mcmc(f, n_walkers, beta_run, mts, prior, lik, y, x);
        auto n_par = current.walkers[0][0].parameter.size();
        auto mcmc_run = checks_convergence(std::move(a), current);
        std::size_t iter = 0;
        const auto start = std::chrono::high_resolution_clock::now();
        report_title(rep, current, lik, y, x);
        report_model(rep, prior, lik, y, x, beta);

        while (beta_run.size() < beta.size() || !mcmc_run.second) {
            while (!mcmc_run.second) {
                step_stretch_thermo_mcmc(f, iter, current, rep, beta_run, mts, prior, lik, y, x);
                thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts, thermo_jumps_every);
                const auto end = std::chrono::high_resolution_clock::now();
                auto dur = std::chrono::duration<double>(end - start);
                report(f, iter, dur, rep, current, prior, lik, y, x);
                mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
            }
            if (beta_run.size() < beta.size()) {
                if (beta_run.size() < max_num_simultaneous_temperatures) {
                    beta_run.insert(beta_run.begin(), beta[beta_run.size()]);
                    current = push_back_new_beta(f, iter, current, mts, beta_run, prior, lik, y, x);
                }

                reset(mcmc_run.first);
                if (current.beta.back() == 1.0)
                    mcmc_run.first.we_reach_final_temperature();

                mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
            }
        }

        return std::pair(mcmc_run, current);
    }
};

template <class FunctionTable, class Algorithm, class Prior, class Likelihood, class Variables,
          class DataType, class Reporter,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

auto thermo_impl(FunctionTable&& f, const Algorithm& alg, Prior const& prior, Likelihood const& lik,
                 const DataType& y, const Variables& x, Reporter rep,
                 std::size_t num_scouts_per_ensemble, std::size_t max_num_simultaneous_temperatures,
                 std::size_t thermo_jumps_every, double n_points_per_decade, double stops_at,
                 bool includes_zero, std::size_t initseed) {
    auto a = alg;
    auto mt = init_mt(initseed);
    auto n_walkers = num_scouts_per_ensemble;
    auto mts = init_mts(mt, num_scouts_per_ensemble / 2);
    auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);

    auto beta_run = by_beta<double>(beta.rend() - 2, beta.rend());
    auto current = init_thermo_mcmc(f, n_walkers, beta_run, mts, prior, lik, y, x);
    auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);
    std::size_t iter = 0;
    const auto start = std::chrono::high_resolution_clock::now();

    report_title(rep, current, lik, y, x);
    report_model(rep, prior, lik, y, x, beta);

    while (beta_run.size() < beta.size() || !mcmc_run.second) {
        while (!mcmc_run.second) {
            step_stretch_thermo_mcmc(f, iter, current, rep, beta_run, mts, prior, lik, y, x);
            thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts, thermo_jumps_every);
            const auto end = std::chrono::high_resolution_clock::now();
            auto dur = std::chrono::duration<double>(end - start);
            report(f, iter, dur, rep, current);
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
        if (beta_run.size() < beta.size()) {
            beta_run.insert(beta_run.begin(), beta[beta_run.size()]);
            current = push_back_new_beta(f, iter, current, mts, beta_run, prior, lik, y, x);
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
    }

    return std::pair(mcmc_run, current);
}

template <class Algorithm, class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> )
class thermodynamic_integration {
    Algorithm alg_;
    Reporter rep_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t max_num_simultaneous_temperatures_;
    std::size_t thermo_jumps_every_;
    double n_points_per_decade_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;

   public:
    thermodynamic_integration(Algorithm&& alg, Reporter&& rep, std::size_t num_scouts_per_ensemble,
                              std::size_t max_num_simultaneous_temperatures,
                              std::size_t thermo_jumps_every, double n_points_per_decade,
                              double stops_at, bool includes_zero, std::size_t initseed)
        : alg_{std::move(alg)},
          rep_{std::move(rep)},
          num_scouts_per_ensemble_{num_scouts_per_ensemble},
          max_num_simultaneous_temperatures_{max_num_simultaneous_temperatures},
          thermo_jumps_every_{thermo_jumps_every},
          n_points_per_decade_{n_points_per_decade},
          stops_at_{stops_at},
          includes_zero_{includes_zero},
          initseed_{initseed} {
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
    auto& max_num_simultaneous_temperatures() const {
        return max_num_simultaneous_temperatures_;
    }
    auto& thermo_jumps_every() const {
        return thermo_jumps_every_;
    }
    auto& n_points_per_decade() const {
        return n_points_per_decade_;
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
};

template <class Algorithm, class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> )
class new_thermodynamic_integration {
    Algorithm alg_;
    Reporter rep_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t thermo_jumps_every_;
    std::size_t beta_size_;
    std::size_t beta_upper_size_;
    std::size_t beta_medium_size_;
    double beta_upper_value_;
    double beta_medium_value_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;
    std::size_t adapt_beta_every_;
    std::string adapt_beta_equalizer_;
    std::string adapt_beta_controler_;
    std::string adapt_beta_variance_;
    double adapt_beta_nu_;
    double adapt_beta_t0_;
    double adapt_beta_threshold_;
    bool adjust_beta_;
    double acceptance_upper_limit_;
    double acceptance_lower_limit_;
    double desired_acceptance_;

   public:
    new_thermodynamic_integration(Algorithm&& alg, Reporter&& rep,
                                  std::size_t num_scouts_per_ensemble,
                                  std::size_t thermo_jumps_every, std::size_t beta_size,
                                  std::size_t beta_upper_size, std::size_t beta_medium_size,
                                  double beta_upper_value, double beta_medium_value,
                                  double stops_at, bool includes_zero, std::size_t initseed)
        : alg_{std::move(alg)},
          rep_{std::move(rep)},
          num_scouts_per_ensemble_{num_scouts_per_ensemble},
          thermo_jumps_every_{thermo_jumps_every},
          beta_size_{beta_size},
          beta_upper_size_{beta_upper_size},
          beta_medium_size_{beta_medium_size},
          beta_upper_value_{beta_upper_value},
          beta_medium_value_{beta_medium_value},

          stops_at_{stops_at},
          includes_zero_{includes_zero},
          initseed_{initseed} {
    }

    new_thermodynamic_integration(Algorithm&& alg, Reporter&& rep,
                                  std::size_t num_scouts_per_ensemble,
                                  std::size_t thermo_jumps_every, std::size_t beta_size,
                                  std::size_t initseed, std::size_t t_adapt_beta_every,
                                  std::string t_adapt_beta_equalizer,
                                  std::string t_adapt_beta_constroler,
                                  std::string t_adapt_beta_variance, double t_adapt_beta_nu,
                                  double t_adapt_beta_t0, double t_adapt_beta_threshold,
                                  bool t_adjust_beta, double t_acceptance_upper_limit,
                                  double t_acceptance_lower_limit, double t_desired_acceptance)
        : alg_{std::move(alg)},
          rep_{std::move(rep)},
          num_scouts_per_ensemble_{num_scouts_per_ensemble},
          thermo_jumps_every_{thermo_jumps_every},
          beta_size_{beta_size},
          beta_upper_size_{},
          beta_medium_size_{},
          beta_upper_value_{},
          beta_medium_value_{},

          stops_at_{},
          includes_zero_{},
          initseed_{initseed},
          adapt_beta_every_{t_adapt_beta_every},
          adapt_beta_equalizer_{t_adapt_beta_equalizer},
          adapt_beta_controler_{t_adapt_beta_constroler},
          adapt_beta_variance_{t_adapt_beta_variance},
          adapt_beta_nu_{t_adapt_beta_nu},
          adapt_beta_t0_{t_adapt_beta_t0},
          adapt_beta_threshold_{t_adapt_beta_threshold},
          adjust_beta_{t_adjust_beta},
          acceptance_upper_limit_{t_acceptance_upper_limit},
          acceptance_lower_limit_{t_acceptance_lower_limit},
          desired_acceptance_{t_desired_acceptance} {
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
    auto& stops_at() const {
        return stops_at_;
    }
    auto& includes_zero() const {
        return includes_zero_;
    }
    auto& initseed() const {
        return initseed_;
    }
    auto& adapt_beta_every() const {
        return adapt_beta_every_;
    }
    auto& adapt_beta_equalizer() const {
        return adapt_beta_equalizer_;
    }
    auto& adapt_beta_controler() const {
        return adapt_beta_controler_;
    }
    auto& adapt_beta_variance() const {
        return adapt_beta_variance_;
    }
    auto& adapt_beta_nu() const {
        return adapt_beta_nu_;
    }
    auto& adapt_beta_t0() const {
        return adapt_beta_t0_;
    }
    auto& adapt_beta_threshold() const {
        return adapt_beta_threshold_;
    }
    bool adjust_beta() const {
        return adjust_beta_;
    }
    auto& desired_acceptance() const {
        return desired_acceptance_;
    }
    auto& acceptance_lower_limit() const {
        return acceptance_lower_limit_;
    }
    auto& acceptance_upper_limit() const {
        return acceptance_upper_limit_;
    }
};

template <class FunctionTable, class Duration, class Parameters, class... saving, class... T>
void report_all(FunctionTable& f, std::size_t iter, const Duration& dur,
                save_mcmc<Parameters, saving...>& s, thermo_mcmc<Parameters>& data,
                T const&... ts) {
    (report(f, iter, dur, get<saving>(s.m_m), data, ts...), ..., 1);
}

template <class FunctionTable, class Algorithm, class Prior, class Likelihood, class Variables,
          class DataType, class Reporter>
    requires(!is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto evidence(FunctionTable&& ff, thermodynamic_integration<Algorithm, Reporter>&& therm,
              Prior const& prior, Likelihood const& lik, const DataType& y, const Variables& x) {
    auto f = ff.fork(var::I_thread(0));
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, therm.num_scouts_per_ensemble() / 2);
    auto beta = get_beta_list(therm.n_points_per_decade(), therm.stops_at(), therm.includes_zero());

    auto it_beta_run_begin = beta.rend() - beta.size();
    auto it_beta_run_end = beta.rend();
    auto beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);

    auto current = init_thermo_mcmc(f, n_walkers, beta_run, mts, prior, lik, y, x);
    // auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    const auto start = std::chrono::high_resolution_clock::now();
    auto& rep = therm.reporter();
    report_title(rep, current, lik, y, x);
    report_model(rep, prior, lik, y, x, beta);

    while (it_beta_run_begin != beta.rbegin() || !mcmc_run.second) {
        while (!mcmc_run.second) {
            step_stretch_thermo_mcmc(f, iter, current, rep, beta_run, mts, prior, lik, y, x);
            thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts, therm.thermo_jumps_every());
            const auto end = std::chrono::high_resolution_clock::now();
            auto dur = std::chrono::duration<double>(end - start);
            report(f, iter, dur, rep, current);
            // using geg=typename
            // decltype(checks_convergence(std::move(mcmc_run.first), current))::eger;
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
        if (it_beta_run_begin != beta.rbegin()) {
            --it_beta_run_begin;
            if (beta_run.size() < therm.max_num_simultaneous_temperatures()) {
                beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);
                current = push_back_new_beta(f, iter, current, mts, beta_run, prior, lik, y, x);
            } else {
                --it_beta_run_end;
                beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);
                current.beta = beta_run;
            }
            std::cerr << "\n  beta_run=" << beta_run[0] << "\n";
            mcmc_run.first.reset();
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
    }

    return std::pair(std::move(mcmc_run.first), current);
}

template <class FunctionTable, class Algorithm, class Prior, class Likelihood, class Variables,
          class DataType, class Reporter>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_evidence_(FunctionTable&& f, new_thermodynamic_integration<Algorithm, Reporter>&& therm,
                      Prior const& prior, Likelihood const& lik, const DataType& y,
                      const Variables& x) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());
    auto beta =
        new_get_beta_list(therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
                          therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
                          therm.includes_zero());

    auto it_beta_run_begin = beta.rend() - beta.size();
    auto it_beta_run_end = beta.rend();
    auto beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);

    auto current = init_thermo_mcmc(f, n_walkers, beta_run, mts, prior, lik, y, x);
    // auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 0;
    const auto start = std::chrono::high_resolution_clock::now();

    auto& rep = therm.reporter();
    report_title(rep, current, lik, y, x);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta);

    var::Event_Timing<200> even_dur(start);
    std::ofstream event_file(f.file() + "event.csv");

    while (!mcmc_run.second) {
        even_dur.record("main_loop_start");
        step_stretch_thermo_mcmc(f, iter, even_dur, current, rep, beta_run, mts, prior, lik, y, x);

        even_dur.record("befor_thermo_jump");
        thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts, therm.thermo_jumps_every());
        even_dur.record("after_thermo_jump");

        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start);
        report_all(f, iter, dur, rep, current, prior, lik, y, x, mts, mcmc_run.first);

        even_dur.record("after_report_all");
        if (iter == 1)
            even_dur.report_title(event_file);
        even_dur.report_iter(event_file, iter);

        // report_point(f, iter);

        // using geg=typename
        // decltype(checks_convergence(std::move(mcmc_run.first), current))::eger;
        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
    }

    return std::pair(std::move(mcmc_run.first), current);
}

template <bool Adapt_beta, class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Reporter, class mcmc, class Parameters,
          class timepoint>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_evidence_loop(FunctionTable&& f,
                          new_thermodynamic_integration<Algorithm, Reporter>&& therm,
                          Prior const& prior, Likelihood const& lik, const DataType& y,
                          const Variables& x, mcmc mcmc_run, std::size_t iter,
                          thermo_mcmc<Parameters>& current, Reporter& rep,
                          by_beta<double>& beta_run, mt_64i& mt, std::vector<mt_64i>& mts,
                          const timepoint& start,
                          const std::chrono::duration<double>& previous_duration) {
    var::Event_Timing<200> even_dur(start);
    std::ofstream event_file(f.file() + "_event_timing.csv");

    while (!mcmc_run.second) {
        even_dur.record("main_loop_start");
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start) + previous_duration;

        report_all(f, iter, dur, rep, current, prior, lik, y, x, mts, mcmc_run.first);
        if constexpr (Adapt_beta) {
            adapt_beta(iter, current, beta_run, therm.adapt_beta_every(),
                       therm.adapt_beta_equalizer(), therm.adapt_beta_controler(),
                       therm.adapt_beta_variance(), therm.desired_acceptance(),
                       therm.adapt_beta_nu(), therm.adapt_beta_t0());
            if (therm.adjust_beta())
                adjust_beta(f, iter, therm.adapt_beta_every(), therm.acceptance_upper_limit(),
                            therm.acceptance_lower_limit(), current, beta_run, mts, prior, lik, y,
                            x);
            if (iter % therm.adapt_beta_every() == 0)
                current.reset_statistics();
        }
        step_stretch_thermo_mcmc(f, iter, even_dur, current, rep, beta_run, mts, prior, lik, y, x);

        even_dur.record("befor_thermo_jump");

        thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts, therm.thermo_jumps_every());

        even_dur.record("after_thermo_jump");

        even_dur.record("after_report_all");
        // report_point(f, iter);

        // using geg=typename
        // decltype(checks_convergence(std::move(mcmc_run.first), current))::eger;
        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        even_dur.record("after_checks_convergence");
        if (iter == 2)
            even_dur.report_title(event_file);
        even_dur.report_iter(event_file, iter);
        if (iter % 10 == 0)
            event_file.flush();
    }
    return std::pair(std::move(mcmc_run.first), current);
}

template <bool Adapt_beta, class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Reporter>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_evidence(FunctionTable&& f, new_thermodynamic_integration<Algorithm, Reporter>&& therm,
                     Prior const& prior, Likelihood const& lik, const DataType& y,
                     const Variables& x) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());
    by_beta<double> beta_run;

    if constexpr (Adapt_beta) {
        beta_run = by_beta<double>(therm.beta_size(), 0);
    } else {
        auto beta =
            new_get_beta_list(therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
                              therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
                              therm.includes_zero());

        auto it_beta_run_begin = beta.rend() - beta.size();
        auto it_beta_run_end = beta.rend();
        beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);
    }
    auto current = init_thermo_mcmc(f, n_walkers, beta_run, mts, prior, lik, y, x);
    // auto n_par = current.walkers[0][0].parameter.size();

    if constexpr (Adapt_beta) {
        beta_run = initial_beta_dts(current);
    }
    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 1;
    const auto start = std::chrono::high_resolution_clock::now();
    auto& rep = therm.reporter();
    report_title(rep, current, lik, y, x);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta_run);
    std::chrono::duration<double> previous_duration(0.0);
    return thermo_evidence_loop<Adapt_beta>(
        f, std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm), prior, lik, y,
        x, mcmc_run, iter, current, rep, beta_run, mt, mts, start, previous_duration);
}

template <bool Adapt_beta, class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Reporter>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_evidence_continuation(const std::string& idName, FunctionTable& f,
                                  new_thermodynamic_integration<Algorithm, Reporter>&& therm,
                                  Prior const& prior, Likelihood const& lik, const DataType& y,
                                  const Variables& x) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());

    by_beta<double> beta;

    if constexpr (Adapt_beta) {
        beta = by_beta<double>(therm.beta_size(), 0);
    } else {
        auto beta_ =
            new_get_beta_list(therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
                              therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
                              therm.includes_zero());

        auto it_beta_run_begin = beta_.rend() - beta_.size();
        auto it_beta_run_end = beta_.rend();
        beta = by_beta<double>(it_beta_run_begin, it_beta_run_end);
    }

    auto current = create_thermo_mcmc(n_walkers, beta, mt, prior);
    auto& rep = therm.reporter();

    auto fname = idName + "__i_beta__i_walker__i_par.csv";
    std::size_t iter = 0;
    auto start = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::ratio<1>> duration;

    current = extract_parameters_last(fname, iter, duration, current);

    a.reset(iter);
    beta = current.beta;
    auto res = calc_thermo_mcmc_continuation(f, n_walkers, beta, mts, prior, lik, y, x, current);

    auto mcmc_run = checks_convergence(std::move(a), current);
    //using return_type=Maybe_error<decltype(std::pair(std::move(mcmc_run.first), current))>;
    report_title(rep, current, lik, y, x);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta);

    return thermo_evidence_loop<Adapt_beta>(
        f, std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm), prior, lik, y,
        x, mcmc_run, iter, current, rep, beta, mt, mts, start, duration);
}

class thermo_max {
    std::string path_;
    std::string filename_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t max_num_simultaneous_temperatures_;
    std::size_t thermo_jumps_every_;
    std::size_t max_iter_warming_;
    std::size_t max_iter_final_;

    double n_points_per_decade_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;

   public:
    thermo_max(std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
               std::size_t max_num_simultaneous_temperatures,

               std::size_t thermo_jumps_every, std::size_t max_iter_warming,
               std::size_t max_iter_equilibrium, double n_points_per_decade, double stops_at,
               bool includes_zero, std::size_t initseed)
        : path_{path},
          filename_{filename},
          num_scouts_per_ensemble_{num_scouts_per_ensemble},
          max_num_simultaneous_temperatures_{max_num_simultaneous_temperatures},
          thermo_jumps_every_{thermo_jumps_every},
          max_iter_warming_{max_iter_warming},
          max_iter_final_{max_iter_equilibrium},
          n_points_per_decade_{n_points_per_decade},
          stops_at_{stops_at},
          includes_zero_{includes_zero},
          initseed_{initseed} {
    }

    template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
              class Parameters =
                  std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
        requires(is_prior<Prior, Parameters, Variables, DataType> &&
                 is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
    auto operator()(FunctionTable&& f, const Prior& prior, const Likelihood& lik, const DataType& y,
                    const Variables& x) {
        return thermo_impl(
            less_than_max_iteration(max_iter_warming_, max_iter_final_), prior, lik, y, x,
            save_mcmc<Parameters, save_likelihood<Parameters>, save_Parameter<Parameters>,
                      save_Evidence>(path_, filename_, 10ul, 10ul, 10ul),
            num_scouts_per_ensemble_, max_num_simultaneous_temperatures_, thermo_jumps_every_,
            n_points_per_decade_, stops_at_, includes_zero_, initseed_);
    }
};

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
auto thermo_max_iter(FunctionTable&& f, const Prior& prior, const Likelihood& lik,
                     const DataType& y, const Variables& x, std::string path, std::string filename,
                     std::size_t num_scouts_per_ensemble,
                     std::size_t max_num_simultaneous_temperatures, std::size_t thermo_jumps_every,
                     std::size_t max_iter_warming, std::size_t max_iter_equilibrium,
                     double n_points_per_decade, double stops_at, bool includes_zero,
                     std::size_t initseed) {
    return thermo_impl(
        f, less_than_max_iteration(max_iter_warming, max_iter_equilibrium), prior, lik, y, x,
        save_mcmc<Parameters, save_Iter, save_likelihood<Parameters>, save_Parameter<Parameters>,
                  save_Evidence, save_Predictions<Parameters>>(path, filename, 1ul, 10ul, 10ul,
                                                               10ul, 100ul),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class Parameters>
auto thermo_by_max_iter(std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
                        std::size_t max_num_simultaneous_temperatures,
                        std::size_t thermo_jumps_every, std::size_t max_iter_warming,
                        std::size_t max_iter_equilibrium, double n_points_per_decade,
                        double stops_at, bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(
        less_than_max_iteration(max_iter_warming, max_iter_equilibrium),
        save_mcmc<Parameters, save_likelihood<Parameters>, save_Parameter<Parameters>,
                  save_Evidence>(path, filename, 10ul, 10ul, 10ul),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class Parameters>
auto thermo_store_every(std::size_t num_scouts_per_ensemble,
                        std::size_t max_num_simultaneous_temperatures,
                        std::size_t thermo_jumps_every, std::size_t save_every_iter,
                        std::size_t max_iter_equilibrium, double n_points_per_decade,
                        double stops_at, bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(
        store_every_n_iter<thermo_mcmc<Parameters>>(save_every_iter, max_iter_equilibrium),
        no_save{}, num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class Parameters>
auto thermo_store_every_and_report(const std::string path, const std::string filename,
                                   std::size_t num_scouts_per_ensemble,
                                   std::size_t max_num_simultaneous_temperatures,
                                   std::size_t thermo_jumps_every, std::size_t save_every_iter,
                                   std::size_t max_iter_equilibrium, double n_points_per_decade,
                                   double stops_at, bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(
        store_every_n_iter<thermo_mcmc<Parameters>>(save_every_iter, max_iter_equilibrium),
        save_mcmc<Parameters, save_likelihood<Parameters>, save_Parameter<Parameters>,
                  save_Evidence>(path, filename, 10ul, 1000ul, 10ul),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
auto thermo_convergence(FunctionTable&& f, const Prior& prior, const Likelihood& lik,
                        const DataType& y, const Variables& x, std::string path,
                        std::string filename, std::size_t num_scouts_per_ensemble,
                        std::size_t max_num_simultaneous_temperatures,
                        std::size_t thermo_jumps_every, std::size_t max_iter,
                        double n_points_per_decade, double stops_at, bool includes_zero,
                        std::size_t initseed) {
    return thermo_impl(
        f, checks_derivative_var_ratio<thermo_mcmc, Parameters>(max_iter * prior.size()), prior,
        lik, y, x,
        save_mcmc<Parameters, save_likelihood<Parameters>, save_Parameter<Parameters>,
                  save_Evidence>(path, filename, 10ul, 100ul, 10ul),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class Parameters>
auto thermo_by_convergence(std::string path, std::string filename,
                           std::size_t num_scouts_per_ensemble,
                           std::size_t max_num_simultaneous_temperatures,
                           std::size_t thermo_jumps_every, std::size_t max_iter,
                           double n_points_per_decade, double stops_at, bool includes_zero,
                           std::size_t initseed) {
    return thermodynamic_integration(
        checks_derivative_var_ratio<thermo_mcmc, Parameters>(max_iter),
        save_mcmc<Parameters, save_likelihood<Parameters>, save_Parameter<Parameters>,
                  save_Evidence>(path, filename, 10ul, 100ul, 10ul),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

template <class Cova>
    requires Covariance<double, Cova>
Maybe_error<by_beta<double>> mean_logLik(
    conjugate, const multivariate_gamma_normal_distribution<double, Cova>& prior,
    const linear_model&, const Matrix<double>& y, const Matrix<double>& X,
    by_beta<double> const& beta) {
    auto a_0 = prior.alpha();
    ;
    auto b_0 = prior.beta();
    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    by_beta<double> mean_logLik_(size(beta));
    auto beta_0 = prior.mean();
    for (std::size_t i = 0; i < size(beta); ++i) {
        auto beta0 = beta[i];
        auto L_n = L_0 + beta0 * SSx;
        auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));
        auto yfit = X * tr(beta_n);
        auto ydiff = y - yfit;
        auto SS = beta0 * xtx(ydiff.value());

        auto a_n = a_0 + beta0 * n / 2.0;
        auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
        double d_a_n = 1.0 * n / 2.0;
        auto d_b_n = 0.5 * xtx(ydiff.value());
        auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                          a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;
        if (mean_logLi)
            mean_logLik_[i] = mean_logLi.value();
        else
            return error_message("Error at beta=" + std::to_string(beta0) + " " +
                                 mean_logLi.error()());
    }
    return mean_logLik_;
}

template <class Parameters, class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_normal_distribution<double, Cova>& prior, double prior_eps_df,
    double prior_eps_variance, const Matrix<double>& y, const Matrix<double>& X,
    by_beta<double> const& beta0) {
    auto beta = beta0;
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    by_beta<std::tuple<Maybe_error<double>, double, Parameters>> mean_logLik;
    mean_logLik.reserve(beta.size());

    //  by_beta<Maybe_error<double>> mean_Ev;
    //  mean_Ev.reserve(beta.size());
    //  by_beta<Maybe_error<double>> mean_logLik_diff;
    //  mean_logLik_diff.reserve(beta.size());

    //  by_beta<Maybe_error<double>> diff_logdetLn;
    //  diff_logdetLn.reserve(beta.size());

    //  by_beta<Maybe_error<double>> diff_alogb;
    //  diff_alogb.reserve(beta.size());

    //  by_beta<Maybe_error<double>> diff_lgamma;
    //  diff_lgamma.reserve(beta.size());

    for (std::size_t i = 0; i < beta.size(); ++i) {
        auto L_n = L_0 + beta[i] * SSx;

        auto beta_n = tr(inv(L_n) * (beta[i] * (tr(X) * y) + (L_0 * tr(prior.mean()))));

        auto yfit = X * tr(beta_n);
        auto ydiff = y - yfit;
        auto SS = beta[i] * xtx(ydiff.value());
        std::cerr << "SS\n" << SS << "\n";

        auto a_n = a_0 + beta[i] * n / 2.0;
        auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

        auto logE_n = -0.5 * beta[i] * n * std::log(2 * std::numbers::pi) +
                      0.5 * (logdet(L_0) - logdet(L_n)) + a_0 * log(b_0) - a_n * log(b_n) +
                      var::lgamma(a_n) - var::lgamma(a_0);
        //    std::cerr<<beta[i]<<"\t"<<"Ev"<<logE_n<<"\n";
        //    std::cerr<<beta[i]<<"\tL_0\t"<<L_0 <<"\n";
        //    std::cerr<<beta[i]<<"\tL_n\t"<<L_n<<"\n";
        //    std::cerr<<beta[i]<<"\tlogdet(L_0) - logdet(L_n)\t"<<logdet(L_0) -
        //    logdet(L_n)<<"\n"; std::cerr<<beta[i]<<"\t+ a_0 * log(b_0) -a_n *
        //    log(b_n)\t"<<+ a_0 * log(b_0) -a_n * log(b_n)<<"\n";
        //    std::cerr<<beta[i]<<"\t+ var::lgamma(a_n) - var::lgamma(a_0)\t"<<+
        //    var::lgamma(a_n) - var::lgamma(a_0)<<"\n";

        double d_a_n = 1.0 * n / 2.0;
        auto d_b_n = 0.5 * xtx(ydiff.value());

        auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                          a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;

        mean_logLik.push_back(std::tuple(mean_logLi, std::log(b_n.value() / a_n), beta_n.value()));
        //    mean_Ev.push_back(logE_n);
        //    diff_logdetLn.push_back(logdet(L_n));
        //    diff_alogb.push_back(a_n * log(b_n));
        //    diff_lgamma.push_back(var::lgamma(a_n));
    }
    /*
for (std::size_t i = 0; i < beta.size(); ++i) {
  beta[i]=std::max(beta[i]*(1+1e-6),beta[i]+1e-9);
  auto L_n = beta[i] * SSx + L_0;

  auto beta_n =
      tr(inv(L_n) * (beta[i] * (tr(X) * y) + (L_0 * tr(prior.mean()))));

  auto yfit = X * tr(beta_n);
  auto ydiff = y - yfit;
  auto SS = beta[i] * xtx(ydiff.value());
  std::cerr<<"SS\n"<<SS<<"\n";

  auto a_n = a_0 + beta[i] * n / 2.0;
  auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);


  auto logE_n = -0.5 * beta[i] * n * std::log(2 * std::numbers::pi) +
                0.5 * (logdet(L_0) - logdet(L_n)) +
                a_0 * log(b_0) - a_n * log(b_n) +
                var::lgamma(a_n) - var::lgamma(a_0);
  double d_a_n = 1.0 * n / 2.0;
  auto d_b_n = 0.5 * xtx(ydiff.value());
  std::cerr<<beta[i]<<"\t"<<"Ev"<<logE_n<<"\n\n";
  std::cerr<<beta[i]<<"\t\t diff
logdet(L_n))\t"<<-(diff_logdetLn[i]-logdet(L_n))/(beta[i]-beta0[i])<<"\n";
  std::cerr<<std::setprecision(12)<<"L_n\n"<<L_n<<"\n";
  std::cerr<<std::setprecision(12)<<"SSx\n"<<SSx<<"\n";
  std::cerr<<beta[i]<<"\t\t inv(L_n) * SSx\t"<<inv(L_n) * SSx<<"\n\n";
  std::cerr<<beta[i]<<"\t\t inv(L_n) \t"<<inv(L_n) <<"\n\n";

  std::cerr<<beta[i]<<"\t\tTrace(inv(L_n) * SSx)\t"<<Trace(inv(L_n) *
SSx)<<"\n\n";
  std::cerr<<beta[i]<<"\t\tTrace(inv(static_cast<Matrix<double>const&>(L_n)) *
SSx)\t"<<Trace(inv(static_cast<Matrix<double>const&>(L_n)) * SSx)<<"\n\n";

  std::cerr<<beta[i]<<"\t\tdiff_alogb[i]-a_n * log(b_n))\t"<<-(diff_alogb[i]-a_n
* log(b_n))/(beta[i]-beta0[i])<<"\n"; std::cerr<<beta[i]<<"\t\t(a_n / b_n *
d_b_n +log(b_n) * d_a_n)\t"<<(a_n / b_n * d_b_n +log(b_n) * d_a_n)<<"\n\n";


  std::cerr<<beta[i]<<"\t\t-(diff_lgamma[i]-var::lgamma(a_n))\t"<<-(diff_lgamma[i]-var::lgamma(a_n))/(beta[i]-beta0[i])<<"\n";
  std::cerr<<beta[i]<<"\t\tdigamma(a_n)  * d_a_n\t"<<+digamma(a_n)  *
d_a_n<<"\n\n";



  auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                    0.5 * Trace(inv(L_n) * SSx) -
                    (a_n / b_n * d_b_n +log(b_n) * d_a_n)+
                    digamma(a_n)  * d_a_n;
  mean_logLik_diff.push_back((logE_n-mean_Ev[i])/(beta[i]-beta0[i]));

}

*/

    //  for (std::size_t i=0; i<beta0.size(); ++i)
    //    std::cerr<<"beta= "<<beta0[i]<<"\tmean_logLi"<<mean_logLik[i]<<"\tdiff:
    //    "<<mean_logLik_diff[i]<<"\n";

    return mean_logLik;
}

#endif  // PARALLEL_TEMPERING_LINEAR_REGRESSION_H
