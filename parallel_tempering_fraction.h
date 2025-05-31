#ifndef PARALLEL_TEMPERING_FRACTION_H
#define PARALLEL_TEMPERING_FRACTION_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "derivative_operator.h"
#include "distributions.h"
#include "general_algorithm_on_containers.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"

template <class Parameters>
struct fraction_mcmc {
    Parameters parameter;
    double logP;
    std::map<std::size_t, logLs> logL_map;

    auto& logL(std::size_t i_frac) const {
        auto it = logL_map.find(i_frac);
        assert(it != logL_map.end());
        return it->second;
    }
};

auto i_fract_beta_to_global_beta(std::vector<std::size_t> const& i_fracs,
                                 std::vector<double> const& beta,
                                 std::vector<std::size_t> num_samples, std::size_t max_samples) {
    std::vector<double> global_beta(i_fracs.size());
    global_beta[0] = 0;

    for (std::size_t i = 1; i < global_beta.size(); ++i) {
        auto n0 = i_fracs[i] > 0 ? num_samples[i_fracs[i] - 1] : 0;
        global_beta[i] = (n0 + (num_samples[i_fracs[i]] - n0) * beta[i]) / max_samples;
    }
    assert(global_beta.back() == 1);
    return global_beta;
}

auto global_beta_to_S(const std::vector<double>& beta) {
    auto iB = beta.back() == 1 ? 1 : 0;
    auto S = std::vector<double>(beta.size() - iB);
    auto n = S.size() - 1;
    S[n] = std::log(1 / beta[n] - 1);
    for (std::size_t i = 1; i < S.size(); ++i)
        S[n - i] = std::log(1 / beta[n - i] - 1 / beta[n - i + 1]);
    return S;
}

auto S_t0_global_beta(const std::vector<double>& S) {
    auto beta = std::vector<double>(S.size());
    auto n = beta.size() - 1;
    beta[n] = 1 / (1 + std::exp(S[n]));
    for (std::size_t i = 1; i < S.size(); ++i)
        beta[n - i] = 1.0 / (1.0 / beta[n - i + 1] + std::exp(S[n - i]));
    return beta;
}

template <class Parameters>
    requires std::is_assignable_v<Parameters, Parameters const&>
struct thermo_fraction_mcmc  //:public thermo_mcmc<Parameters>
{
    std::vector<std::size_t> samples_size;
    std::size_t max_samples;
    by_beta<std::size_t> i_fraction;
    by_beta<std::vector<std::size_t>> i_fractions;
    by_beta<double> beta;
    by_beta<double> global_beta;
    by_beta<ensemble<fraction_mcmc<Parameters>>> walkers;
    by_beta<ensemble<std::map<std::size_t, logL_statistics>>> walkers_sta;

    by_beta<ensemble<std::size_t>> i_walkers;
    by_beta<emcee_Step_statistics> emcee_stat;
    by_beta<Thermo_Jump_statistics> thermo_stat;

    auto num_fractions() const {
        return i_fraction.size();
    }
    void reset_statistics() {
        for (std::size_t i_beta = 0; i_beta < walkers_sta.size(); ++i_beta) {
            for (auto& ee : walkers_sta[i_beta]) {
                for (auto& eee : ee) {
                    eee.second().reset();
                }
            }
        }
        for (auto& e : emcee_stat) e().reset();
        for (auto& e : thermo_stat) e().reset();
    }

    auto get_Walkers_number() const {
        return walkers[0].size();
    }
    auto& get_Beta() const {
        return beta;
    }
    auto& get_Parameter(std::size_t iw, std::size_t i_b) const {
        return walkers[i_b][iw].parameter;
    }

    auto get_Walker(std::size_t iw, std::size_t i_b) const {
        return i_walkers[i_b][iw];
    }

    auto num_samples() const {
        std::size_t min_samples = get<count>(walkers_sta[0][0].begin()->second()())();
        for (std::size_t i = 0; i < walkers_sta.size(); ++i)
            for (std::size_t j = 0; j < walkers_sta[i].size(); ++j)
                if (get<count>(walkers_sta[0][0].begin()->second()())() < min_samples)
                    min_samples = get<count>(walkers_sta[0][0].begin()->second()())();
        return min_samples;
    }

    auto last_walker() const {
        std::size_t last = 0;
        for (auto& e : i_walkers)
            for (auto iw : e)
                if (iw > last)
                    last = iw;
        return last;
    }
};

template <class Parameters>
std::size_t num_walkers(thermo_fraction_mcmc<Parameters> const& x) {
    return x.walkers[0].size();
}

template <class Parameters>
std::size_t num_Parameters(thermo_fraction_mcmc<Parameters> const& x) {
    return x.walkers[0][0].parameter.size();
}

template <class Parameters>
auto calculate_initial_logL0(thermo_fraction_mcmc<Parameters> const& initial_data) {
    auto logL_0 = std::map<std::size_t, logL_statistics>{};
    for (std::size_t i = 0; i < initial_data.walkers_sta.size(); ++i) {
        for (std::size_t j = 0; j < initial_data.walkers_sta[i].size(); ++j) {
            for (auto& e : initial_data.walkers_sta[i][j]) {
                logL_0[e.first]() &= e.second();
            }
        }
    }

    return logL_0;
}
inline double calculate_logL_mean(ensemble<std::map<std::size_t, logL_statistics>> const& sta,
                                  std::size_t i_fr) {
    double sum = 0.0;
    for (std::size_t i = 0; i < sta.size(); ++i) {
        auto it = sta[i].find(i_fr);
        if (it != sta[i].end())
            sum += get<mean<logL>>(it->second()())()();
        else {
            std::cerr << "error in " << i << " i_fr" << i_fr << "\n";
            abort();
        }
    }
    return sum / sta.size();
    ;
}

inline std::map<std::size_t, double> calculate_logL_mean(
    ensemble<std::map<std::size_t, logL_statistics>> const& sta,
    const std::vector<std::size_t>& i_fracs) {
    std::map<std::size_t, double> sum;
    for (std::size_t i = 0; i < sta.size(); ++i)
        for (auto i_fr : i_fracs) {
            auto it = sta[i].find(i_fr);
            if (it != sta[i].end())
                sum[i_fr] += get<mean<logL>>(it->second()())()();
            else {
                std::cerr << "error in " << i << " i_fr" << i_fr << "\n";
                abort();
            }
        }
    for (auto i_fr : i_fracs) sum[i_fr] = sum[i_fr] / sta.size();

    return sum;
}
inline by_beta<std::map<std::size_t, double>> calculate_logL_mean(
    by_beta<ensemble<std::map<std::size_t, logL_statistics>>> const& sta,
    by_beta<std::vector<std::size_t>> const& i_fracs) {
    by_beta<std::map<std::size_t, double>> out(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i) out[i] = calculate_logL_mean(sta[i], i_fracs[i]);
    return out;
}

template <class Parameters>
inline std::map<std::size_t, double> calculate_logL_mean(
    ensemble<fraction_mcmc<Parameters>> const& wa, const std::vector<std::size_t>& i_fracs) {
    std::map<std::size_t, double> sum;
    for (std::size_t i = 0; i < wa.size(); ++i)
        for (auto i_fr : i_fracs) {
            sum[i_fr] += get<logL>(wa[i].logL(i_fr))();
        }
    for (auto i_fr : i_fracs) sum[i_fr] = sum[i_fr] / wa.size();

    return sum;
}

template <class Parameters>
inline by_beta<std::map<std::size_t, double>> calculate_logL_mean(
    by_beta<ensemble<fraction_mcmc<Parameters>>> const& wa,
    by_beta<std::vector<std::size_t>> const& i_fracs) {
    by_beta<std::map<std::size_t, double>> out(wa.size());
    for (std::size_t i = 0; i < wa.size(); ++i) out[i] = calculate_logL_mean(wa[i], i_fracs[i]);
    return out;
}

inline std::map<std::size_t, logL_statistics> calculate_across_sta(
    ensemble<std::map<std::size_t, logL_statistics>>& sta, std::vector<std::size_t> i_frac) {
    std::map<std::size_t, logL_statistics> across;
    for (std::size_t iw = 0; iw < sta.size(); ++iw)
        for (auto i_fr : i_frac) {
            auto it = sta[iw].find(i_fr);
            assert(it != sta[iw].end());
            across[i_fr]() &= get<mean<logL>>(it->second()())();
        }
    return across;
}
inline auto calculate_across_sta(by_beta<ensemble<std::map<std::size_t, logL_statistics>>>& sta,
                                 by_beta<std::vector<std::size_t>> i_frac) {
    by_beta<std::map<std::size_t, logL_statistics>> across(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i)
        across[i] = calculate_across_sta(sta[i], i_frac[i]);
    return across;
}

inline auto calculate_sample_size(by_beta<ensemble<std::map<std::size_t, logL_statistics>>>& sta,
                                  std::size_t i_beta, std::size_t i_frac) {
    auto it = sta[i_beta][0].find(i_frac);
    assert(it != sta[i_beta][0].end());
    return get<count>(it->second()());
}

inline variance<logL> calculate_within_sta(ensemble<std::map<std::size_t, logL_statistics>>& sta,
                                           std::size_t i_frac) {
    variance<logL> within(logL(0.0));
    count df(count(0));
    for (std::size_t i = 0; i < sta.size(); ++i) {
        auto it = sta[i].find(i_frac);
        assert(it != sta[i].end());
        auto r_df = get<count>(it->second()())() - 1;
        within()() += r_df * get<variance<logL>>(it->second()())()();
        df() += r_df;
    }
    within()() /= df();
    return within;
}

inline auto calculate_within_sta(by_beta<ensemble<std::map<std::size_t, logL_statistics>>>& sta,
                                 by_beta<std::size_t> i_frac) {
    by_beta<variance<logL>> within(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i)
        within[i] = calculate_within_sta(sta[i], i_frac[i]);
    return within;
}

template <class Parameters>
auto calculate_deltaBeta_deltaL(const thermo_fraction_mcmc<Parameters>& current) {
    std::vector<double> dBdL(current.walkers.size() - 1);

    auto L = calculate_logL_mean(current.walkers_sta, current.i_fractions);

    for (std::size_t i = 0; i < dBdL.size(); ++i) {
        auto i_frac1 = current.i_fraction[i + 1];
        auto i_frac = current.i_fraction[i];
        auto deltabeta =
            i_frac1 == i_frac ? current.beta[i + 1] - current.beta[i] : current.beta[i + 1];
        auto deltaL = i_frac1 > 0 ? L[i + 1][i_frac1] - L[i][i_frac1] - L[i + 1][i_frac1 - 1] +
                                        L[i][i_frac1 - 1]
                                  : L[i + 1][i_frac1] - L[i][i_frac1];

        dBdL[i] = deltaL * deltabeta;
    }
    return dBdL;
}

template <class Parameters>
auto calculate_partial_Evidence_sta(const thermo_fraction_mcmc<Parameters>& current) {
    std::vector<logEv_statistics> out(current.beta.size() - 1);

    for (std::size_t i = 0; i < out.size(); ++i) {
        auto i_frac1 = current.i_fraction[i + 1];
        auto i_frac = current.i_fraction[i];
        auto deltabeta =
            i_frac1 == i_frac ? current.beta[i + 1] - current.beta[i] : current.beta[i + 1];
        if (i_frac1 > 0) {
            auto L1 = current.walkers[i + 1][0].logL(i_frac1) -
                      current.walkers[i + 1][0].logL(i_frac1 - 1);
            auto L0 = current.walkers[i][0].logL(i_frac1) - current.walkers[i][0].logL(i_frac1 - 1);

            out[i] = logEv_statistics(logEv(deltabeta * get<logL>(L1 - L0)()));
        } else {
            auto L1 = current.walkers[i + 1][0].logL(i_frac1);
            auto L0 = current.walkers[i][0].logL(i_frac1);

            out[i] = logEv_statistics(logEv(deltabeta * get<logL>(L1 - L0)()));
        }
        for (auto i_w = 1ul; i_w < num_walkers(current); ++i_w) {
            if (i_frac1 > 0) {
                auto L1 = current.walkers[i + 1][i_w].logL(i_frac1) -
                          current.walkers[i + 1][i_w].logL(i_frac1 - 1);
                auto L0 = current.walkers[i][i_w].logL(i_frac1) -
                          current.walkers[i][i_w].logL(i_frac1 - 1);

                out[i]() &= logEv_statistics(logEv(deltabeta * get<logL>(L1 - L0)()))();
            } else {
                auto L1 = current.walkers[i + 1][i_w].logL(i_frac1);
                auto L0 = current.walkers[i][i_w].logL(i_frac1);

                out[i]() &= logEv_statistics(logEv(deltabeta * get<logL>(L1 - L0)()))();
            }
        }
    }
    return out;
}

template <class Parameters>
auto mean_logL(by_iteration<thermo_fraction_mcmc<Parameters>>& series) {
    auto out = by_beta<std::map<std::size_t, double>>(num_betas(series[0]));
    auto n_walkers = num_walkers(series[0]);
    auto n_iters = num_samples(series);
    for (std::size_t i = 0; i < num_samples(series); ++i)
        for (std::size_t ibeta = 0; ibeta < num_betas(series[0]); ++ibeta)
            for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
                for (auto& e : series[i].walkers[ibeta][iwalker].logL)
                    out[ibeta][e.first] += e.second / n_iters / n_walkers;
    return out;
}

template <class Parameters>
auto mean_logL_per_fraction(thermo_fraction_mcmc<Parameters> const& mcmc) {
    auto out = by_beta<std::map<std::size_t, logLs>>(num_betas(mcmc));
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            for (auto i_fr : mcmc.i_fractions[ibeta])
                out[ibeta][i_fr] =
                    out[ibeta][i_fr] + mcmc.walkers[ibeta][iwalker].logL(i_fr) / n_walkers;
    return out;
}

template <class Parameters>
auto mean_logP(thermo_fraction_mcmc<Parameters> const& mcmc) {
    auto out = by_beta<double>(num_betas(mcmc), 0);
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            out[ibeta] += mcmc.walkers[ibeta][iwalker].logP / n_walkers;
    return out;
}

template <class Parameters>
auto var_logL(thermo_fraction_mcmc<Parameters> const& mcmc,
              by_beta<std::map<std::size_t, logLs>>& mean) {
    auto out = by_beta<std::map<std::size_t, logLs>>(num_betas(mcmc));
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            for (auto i_fr : mcmc.i_fractions[ibeta])
                out[ibeta][i_fr] =
                    out[ibeta][i_fr] +
                    pow(mcmc.walkers[ibeta][iwalker].logL(i_fr) - mean[ibeta][i_fr], 2) / n_walkers;
    return out;
}

auto S_to_i_fract_beta(std::vector<double> const& S, std::vector<std::size_t> num_samples,
                       std::size_t max_samples) {
    auto frac_samples =
        var::apply_to([max_samples](auto x) { return 1.0 * x / max_samples; }, num_samples);
    std::vector<std::size_t> frac_i(S.size());
    std::vector<double> beta(S.size());
    auto n = S.size() - 1;
    auto nf = num_samples.size() - 1;
    std::size_t i_frac = nf;
    auto i_frac0 = i_frac > 0 ? i_frac - 1 : 0;
    // qué onda? da para que sea la ultima fraccion?
    // el threshold debería ser a mitad de camino entre 0 y la proxima fraccion.

    auto threshold_s = -std::log(frac_samples[i_frac] / 2.0 + frac_samples[i_frac0] / 2.0);

    while ((i_frac > 0) && (S[n] > threshold_s)) {
        --i_frac;
        threshold_s = -std::log(frac_samples[i_frac] / 2.0 + frac_samples[i_frac - 1] / 2.0);
    }
    // tenemos el primer i_frac! y el beta (debe ser 1)
    beta[n - 0] = 1.0;
    frac_i[n - 0] = i_frac;
    double current_T = 1.0 / frac_samples[i_frac];
    // el proximo threshold esta dado por
    for (std::size_t i = 1; i < S.size(); ++i) {
        if ((i_frac > 0) &&
            (current_T + std::exp(S[n - i] - S[n - i + 1]) > 1.0 / frac_samples[i_frac - 1])) {
            --i_frac;
            beta[n - i] = 1.0;
            frac_i[n - i] = i_frac;
            current_T = 1.0 / frac_samples[i_frac];

        } else {
            if (std::isfinite(S[n - i])) {
                current_T = current_T + std::exp(S[n - i] - S[n - i + 1]);
                auto gbeta = 1 / current_T;
                auto fraction_below = i_frac > 0 ? frac_samples[i_frac - 1] : 0.0;

                auto r_beta =
                    (gbeta - fraction_below) / (1.0 / frac_samples[i_frac] - fraction_below);
                beta[n - i] = r_beta;
            } else {
                beta[n - i] = 0;
            }
            frac_i[n - i] = i_frac;
        }
    }
    return std::tuple(std::move(frac_i), std::move(beta));
}

auto global_beta_to_i_fract_beta(std::vector<double> const& global_beta,
                                 std::vector<std::size_t> p_i_frac,
                                 std::vector<std::size_t> num_samples, std::size_t max_samples,
                                 double threshold) {
    auto samples_global_beta =
        var::apply_to([max_samples](auto x) { return 1.0 * x / max_samples; }, num_samples);
    std::vector<std::size_t> frac_i(global_beta.size());
    std::vector<double> beta(global_beta.size());
    frac_i[0] = 0;
    beta[0] = 0;
    frac_i.back() = num_samples.size() - 1;
    beta.back() = 1;
    // first re-assign i_fracs
    for (std::size_t i = 1; i + 1 < global_beta.size(); ++i) {
        auto dgb = threshold * (global_beta[i] - global_beta[i - 1]);
        if (frac_i[i - 1] == p_i_frac[i] && p_i_frac[i] + 1 < num_samples.size() &&
            samples_global_beta[p_i_frac[i]] + dgb < global_beta[i]) {
            frac_i[i] = p_i_frac[i] + 1;
        } else if (p_i_frac[i + 1] == p_i_frac[i] && p_i_frac[i] > 0 &&
                   (samples_global_beta[p_i_frac[i] - 1] > dgb + global_beta[i])) {
            frac_i[i] = p_i_frac[i] - 1;
        } else {
            frac_i[i] = p_i_frac[i];
        }
    }
    // find the global_betas junctures
    std::vector<double> global_beta_junctures(num_samples.size() + 1);
    std::vector<double> global_beta_segments(num_samples.size());
    global_beta_junctures[0] = 0;
    global_beta_junctures.back() = 1;
    for (std::size_t i = 1; i + 1 < global_beta.size(); ++i) {
        if (frac_i[i + 1] > frac_i[i])
            global_beta_junctures[frac_i[i + 1]] = global_beta[i];
    }
    for (std::size_t i = 0; i < global_beta_segments.size(); ++i) {
        global_beta_segments[i] = global_beta_junctures[i + 1] - global_beta_junctures[i];
    }
    // now calculate the betas
    for (std::size_t i = 1; i + 1 < global_beta.size(); ++i) {
        beta[i] =
            (global_beta[i] - global_beta_junctures[frac_i[i]]) / global_beta_segments[frac_i[i]];
    }
    return std::tuple(std::move(frac_i), std::move(beta));
}

auto i_frac_to_i_fracs(const std::vector<std::size_t>& i_frac) {
    std::vector<std::vector<std::size_t>> i_fracs(i_frac.size());
    for (std::size_t i = 0; i < i_frac.size(); ++i) {
        if (i_frac[i] == 0) {
            if ((i + 1) < i_frac.size() && i_frac[i + 1] > i_frac[i])
                i_fracs[i] = {i_frac[i], i_frac[i] + 1};
            else
                i_fracs[i] = {i_frac[i]};
        } else if (i + 1 < i_frac.size() && i_frac[i + 1] > i_frac[i])
            i_fracs[i] = {i_frac[i] - 1, i_frac[i], i_frac[i] + 1};
        else
            i_fracs[i] = {i_frac[i] - 1, i_frac[i]};
    }
    return i_fracs;
}

template <class Parameters>
ensemble<std::map<std::size_t, logL_statistics>> make_logL_statistics(
    ensemble<fraction_mcmc<Parameters>> const& walkers) {
    ensemble<std::map<std::size_t, logL_statistics>> out(walkers.size());
    for (auto i = 0ul; i < walkers.size(); ++i) {
        for (auto const& e : walkers[i].logL_map)
            out[i][e.first] = logL_statistics(get<logL>(e.second));
    }
    return out;
}

template <class Parameters>
by_beta<ensemble<std::map<std::size_t, logL_statistics>>> make_logL_statistics(
    by_beta<ensemble<fraction_mcmc<Parameters>>> const& walkers) {
    by_beta<ensemble<std::map<std::size_t, logL_statistics>>> out;
    out.reserve(walkers.size());
    for (auto i = 0ul; i < walkers.size(); ++i) {
        out.push_back(ensemble<std::map<std::size_t, logL_statistics>>(walkers[i].size()));

        for (std::size_t j = 0; j < walkers[i].size(); ++j) {
            for (auto const& e : walkers[i][j].logL_map)
                out[i][j][e.first] = logL_statistics(get<logL>(e.second));
        }
    }
    return out;
}

template <class LikelihoodModel, class FuncTable, class Parameters, class Variables, class DataType>
Maybe_error<std::map<std::size_t, logLs>> logLikelihoods(FuncTable& f, const LikelihoodModel& lik,
                                                         Parameters const& p,
                                                         const std::vector<DataType>& y,
                                                         const std::vector<Variables>& var,
                                                         std::vector<std::size_t> i_fracs) {
    std::map<std::size_t, logLs> out;
    for (auto i_frac : i_fracs) {
        auto Maybe_lik = logLikelihood(f, lik, p, y[i_frac], var[i_frac]);
        if (!Maybe_lik.valid())
            return Maybe_lik.error();
        else
            out[i_frac] = std::move(Maybe_lik.value());
    }
    return out;
}

template <class LikelihoodModel, class FuncTable, class Parameters, class Variables, class DataType>
Maybe_error<bool> update_logLikelihoods(FuncTable& f, const LikelihoodModel& lik,
                                        Parameters const& p, const std::vector<DataType>& y,
                                        const std::vector<Variables>& var,
                                        std::vector<std::size_t> i_fracs,
                                        std::map<std::size_t, logLs>& rlogL,
                                        std::map<std::size_t, logL_statistics>& rlog_sta) {
    Maybe_error<bool> out(true);
    for (auto i_frac : i_fracs) {
        if (rlogL.find(i_frac) == rlogL.end()) {
            auto Maybe_lik = logLikelihood(f, lik, p, y[i_frac], var[i_frac]);
            if (!Maybe_lik.valid()) {
                rlogL[i_frac] = nan_logL();
                out = error_message(out.error()() + Maybe_lik.error()());
            } else {
                rlog_sta[i_frac] = logL_statistics(get<logL>(Maybe_lik.value()));
                rlogL[i_frac] = std::move(Maybe_lik.value());
            }
        }
    }
    return out;
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_fraction_mcmc(FunctionTable& f, mt_64i& mt, Prior const& pr, const Lik& lik,
                        const std::vector<DataType>& y, const std::vector<Variables>& x,
                        std::vector<std::size_t> i_frac) {
    auto& priorsampler = pr;
    auto par = sample(mt, priorsampler);
    auto logP = logPrior(pr, par);
    std::map<std::size_t, logLs> t_logsLs;
    auto logL = logLikelihoods(f, lik, par.to_value(), y, x, i_frac);
    while (!(logP) || !(logL)) {
        par = sample(mt, priorsampler);
        logP = logPrior(pr, par);
        logL = logLikelihoods(f, lik, par.to_value(), y, x, i_frac);
    }
    return fraction_mcmc<Parameters>{std::move(par), logP.value(), std::move(logL.value())};
}

template <class Parameters>
auto calculate_controler_step(const thermo_fraction_mcmc<Parameters>& current,
                              std::string equalizing_paramter, double desired_acceptance,
                              std::string variance_approximation) {
    if (equalizing_paramter == "Acceptance_vfm") {
        auto A = calculate_Acceptance(current);
        auto d = std::vector<double>(A.size() - 1);
        for (std::size_t i = 0; i < d.size(); ++i) {
            d[i] = std::max(std::min(1.0, A[i + 1] - A[i]), -1.0);
        }

        return d;
    } else if (equalizing_paramter == "deltaBeta_deltaL_vfm") {
        auto dBdL = calculate_deltaBeta_deltaL(current);
        auto d = std::vector<double>(dBdL.size());

        for (std::size_t i = 0; i + 1 < dBdL.size(); ++i) {
            d[i] = std::max(std::min(1.0, std::exp(std::min(0.0, -dBdL[i + 1])) -
                                              std::exp(std::min(0.0, -dBdL[i]))),
                            -1.0);
        }
        d.back() = -std::exp(std::min(0.0, -dBdL.back()));
        return d;
    } else
        return std::vector<double>{};
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
void insert_high_temperture_beta(FunctionTable& f, std::size_t tested_index,
                                 thermo_fraction_mcmc<Parameters>& current, double new_beta,
                                 ensemble<mt_64i>& mt, Prior const& prior, Likelihood const& lik,
                                 const DataType& y, const Variables& x) {
    auto n_walkers = current.walkers[0].size();
    ensemble<std::size_t> i_walker(n_walkers);
    ensemble<fraction_mcmc<Parameters>> walker(n_walkers);
    emcee_Step_statistics emcee_stat(emcee_Step_statistics{});
    Thermo_Jump_statistics thermo_stat(Thermo_Jump_statistics{});
    auto ff = f.fork(omp_get_max_threads());

    auto last_walker = current.last_walker();

    auto i_fracts_0 = std::vector<std::size_t>{0};

#pragma omp parallel for  // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        i_walker[iw] = iw + last_walker + 1;
        auto i_th = omp_get_thread_num();
        walker[iw] = init_fraction_mcmc(ff[i_th], mt[i_th], prior, lik, y, x, i_fracts_0);
    }
    f += ff;
    // update S, i_frac etc
    auto walker_sta = make_logL_statistics(walker);
    double new_global_beta =
        std::sqrt(current.global_beta[tested_index + 1] * current.global_beta[tested_index]);
    ;
    current.global_beta.insert(current.global_beta.begin() + tested_index + 1, new_global_beta);

    current.beta.insert(current.beta.begin() + tested_index + 1, new_beta);
    current.i_fractions.insert(current.i_fractions.begin() + tested_index + 1, {0});
    current.i_fraction.insert(current.i_fraction.begin() + tested_index + 1, 0);
    current.i_walkers.insert(current.i_walkers.begin(), i_walker);
    current.walkers.insert(current.walkers.begin(), walker);
    current.walkers_sta.insert(current.walkers_sta.begin(), walker_sta);
    current.emcee_stat.insert(current.emcee_stat.begin(), emcee_stat);
    current.thermo_stat.insert(current.thermo_stat.begin(), thermo_stat);
}

template <class Parameters>
void remove_high_temperture_beta(thermo_fraction_mcmc<Parameters>& current,
                                 std::size_t tested_index) {
    current.beta.erase(current.beta.begin() + tested_index + 1);
    current.i_walkers.erase(current.i_walkers.begin());
    current.walkers.erase(current.walkers.begin());
    current.walkers_sta.erase(current.walkers_sta.begin());
    current.emcee_stat.erase(current.emcee_stat.begin());
    current.thermo_stat.erase(current.thermo_stat.begin());
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
void adjust_fraction_beta(FunctionTable& f, std::size_t iter, std::size_t adapt_beta_every,
                          double acceptance_upper_limit, double acceptance_lower_limit,
                          thermo_fraction_mcmc<Parameters>& current, ensemble<mt_64i>& mt,
                          Prior const& prior, Likelihood const& lik, const DataType& ys,
                          const Variables& xs) {
    if ((iter > 0) && (current.num_samples() > 0) && (iter % adapt_beta_every == 0)) {
        assert(current.beta[current.beta.size() - 1] = 1);
        std::size_t tested_index = 1;
        // auto A=calculate_Acceptance(current);
        update_likelihoods_all(f, current, lik, ys, xs);
        auto A = calculate_deltaBeta_deltaL(current);
        if (false && std::exp(std::min(0.0, -A[tested_index])) > acceptance_upper_limit)
            remove_high_temperture_beta(current, tested_index);
        else
            //             if (std::reduce(A.begin()+tested_index, A.end(),0.0,
            // [acceptance_lower_limit](double a,double b){
            //          if ((a==1)||(b<acceptance_lower_limit))
            //               return 1.0;
            //          else
            // return 0.0;})==1.0)
            if (std::exp(std::min(0.0, -A[tested_index + 1])) < acceptance_lower_limit) {
                double new_beta;
                if (current.beta[tested_index] == 0) {
                    new_beta =
                        current.beta[tested_index + 1] *
                        std::sqrt(current.beta[tested_index + 1] / current.beta[tested_index + 2]);
                } else {
                    new_beta =
                        std::sqrt(current.beta[tested_index + 1] * current.beta[tested_index]);
                }
                insert_high_temperture_beta(f, tested_index, current, new_beta, mt, prior, lik, ys,
                                            xs);
                update_likelihoods_all(f, current, lik, ys, xs);
            }
    }
}

template <class FunctionTable, class Likelihood, class Variables, class DataType, class Parameters>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

void adapt_fraction_beta(FunctionTable& f, std::size_t iter,
                         thermo_fraction_mcmc<Parameters>& current, std::size_t adapt_beta_every,
                         std::string equalizing_paramter, std::string variance_approximation,
                         double desired_acceptance, double nu, double t0,
                         double beta_adapt_threshold, Likelihood const& lik, const DataType& ys,
                         const Variables& xs) {
    if ((iter > 0) && (current.num_samples() > 0) && (iter % adapt_beta_every == 0)) {
        // assert(beta[beta.size() - 1] = 1);
        double kappa = 1.0 / nu * t0 / (t0 + iter);

        update_likelihoods_all(f, current, lik, ys, xs);

        auto d = calculate_controler_step(current, equalizing_paramter, desired_acceptance,
                                          variance_approximation);

        auto beta = current.beta;
        std::vector<double> T(beta.size() - (beta[0] == 0 ? 1 : 0));
        for (std::size_t i = 0; i < T.size(); ++i) {
            T[i] = 1.0 / beta[i + (beta[0] == 0 ? 1 : 0)];
        }
        std::vector<double> S(T.size() - 1);

        for (std::size_t i = 0; i < S.size(); ++i) {
            S[S.size() - 1 - i] = std::log(T[T.size() - 2 - i] - T[T.size() - 1 - i]);
        }
        if (equalizing_paramter.ends_with("vfm")) {
            for (std::size_t i = 0; i + 1 < d.size(); ++i) {
                S[i] += kappa * d[i];
            }
        } else {
            auto dbds = calculate_d_beta_d_s(beta);
            for (std::size_t i = 0; i < S.size(); ++i) {
                S[i] += kappa * d[i] / dbds[i];
            }
        }
        std::size_t tested_index = 1;
        for (std::size_t i = 0; i < S.size() - tested_index; ++i) {
            if (std::isfinite(S[S.size() - 1 - i])) {
                T[T.size() - 2 - i] = T[T.size() - 1 - i] + std::exp(S[S.size() - 1 - i]);
            }
        }
        for (std::size_t i = 0; i < T.size(); ++i) beta[i + (beta[0] == 0 ? 1 : 0)] = 1.0 / T[i];

        current.beta = beta;
    }
}

template <class FunctionTable, class Likelihood, class Variables, class DataType, class Parameters>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

void adapt_fraction_i_fraction(FunctionTable& f, std::size_t iter,
                               thermo_fraction_mcmc<Parameters>& current,
                               std::size_t adapt_beta_every, std::string equalizing_paramter,
                               std::string variance_approximation, double desired_acceptance,
                               double nu, double t0, double beta_adapt_threshold,
                               double beta_upper_limit, Likelihood const& lik, const DataType& ys,
                               const Variables& xs) {
    if ((iter > 0) && (current.num_samples() > 0) && (iter % adapt_beta_every == 0)) {
        // assert(beta[beta.size() - 1] = 1);
        double kappa = 1.0 / nu * t0 / (t0 + iter);

        update_likelihoods_all(f, current, lik, ys, xs);

        auto d = calculate_controler_step(current, equalizing_paramter, desired_acceptance,
                                          variance_approximation);

        auto i_fraction = current.i_fraction;
        auto new_i_fraction = i_fraction;

        auto beta = current.beta;
        auto new_beta = beta;
        std::size_t tested_index = 1;

        bool something_changed = false;
        for (std::size_t i = tested_index; i + 1 < beta.size(); ++i) {
            if (i_fraction[i] < i_fraction[i + 1]) {
                if (false && (beta[i + 1] < 1.0) && (d[i - tested_index] > beta_adapt_threshold)) {
                    something_changed = true;
                    new_i_fraction[i + 1] = i_fraction[i];

                    new_beta[i + 1] = 1;
                    if (beta[i - 1] < 1)
                        new_beta[i] = std::sqrt(beta[i - 1]);
                    else
                        beta[i] =
                            -1.0 / (calculate_logL_mean(current.walkers_sta[i], i_fraction[i]) -
                                    calculate_logL_mean(current.walkers_sta[i], i_fraction[i - 1]));

                } else if ((beta[i - 1] < 1.0) && (-d[i - tested_index] > beta_adapt_threshold)) {
                    something_changed = true;
                    new_i_fraction[i] = i_fraction[i + 1];
                    new_beta[i - 1] = 1;
                    if (beta[i + 1] < 1)
                        new_beta[i] = beta[i + 1] * beta[i + 1];
                    else
                        new_beta[i] =
                            -1.0 / (calculate_logL_mean(current.walkers_sta[i], i_fraction[i + 1]) -
                                    calculate_logL_mean(current.walkers_sta[i], i_fraction[i]));
                    ;
                }
            }
        }
        auto n = beta.size() - 1;
        if (i_fraction.back() + 1 < current.samples_size.size() && (beta[beta.size() - 2] < 1.0) &&
            (-d[n - tested_index] > beta_upper_limit)) {
            something_changed = true;
            new_i_fraction[n] = i_fraction[n] + 1;
            new_beta[n - 1] = 1;
            new_beta[n] = 1;
        }

        if (something_changed) {
            current.i_fraction = std::move(new_i_fraction);
            current.beta = std::move(new_beta);
            current.i_fractions = i_frac_to_i_fracs(current.i_fraction);
            update_likelihoods_all(f, current, lik, ys, xs);
        }
    }
}

template <class Parameters>
void update(by_beta<ensemble<std::map<std::size_t, logL_statistics>>>& walkers_sta,
            by_beta<ensemble<fraction_mcmc<Parameters>>> const& walkers,
            std::vector<std::vector<std::size_t>> const& i_fractions) {
    for (auto i = 0ul; i < walkers_sta.size(); ++i)
        for (std::size_t j = 0; j < walkers_sta[i].size(); ++j) {
            for (auto i_fr : i_fractions[i]) {
                walkers_sta[i][j][i_fr]() &= get<logL>(walkers[i][j].logL(i_fr));
            }
        }
}

template <class Walker, class logP, class logLik>
double calc_emcee_jump_delta(Walker const& current, logP const& ca_logP, logLik const& ca_logL,
                             double beta, std::size_t i_frac) {
    auto delta_P = ca_logP - current.logP;
    auto it = ca_logL.find(i_frac);
    assert(it != ca_logL.end());
    auto& calogL = it->second;
    auto delta_lik = beta * (get<logL>(calogL)() - get<logL>(current.logL(i_frac))());
    if (i_frac > 0) {
        auto it = ca_logL.find(i_frac - 1);
        assert(it != ca_logL.end());
        auto& calogL0 = it->second;

        delta_P =
            delta_P + (1 - beta) * (get<logL>(calogL0)() - get<logL>(current.logL(i_frac - 1))());
    }
    auto out = delta_P + delta_lik;
    // it is possible that current is not finite
    if (std::isfinite(out))
        return out;
    else
        return +20.0;
}

auto get_samples_size(auto const& y) {
    std::vector<std::size_t> out(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        out[i] = y[i]().size();
    }
    return out;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
auto init_thermo_fraction_mcmc(FunctionTable& f, std::size_t n_walkers, std::size_t beta_size,
                               ensemble<mt_64i>& mt, Prior const& prior, Likelihood const& lik,
                               const DataType& y, const Variables& x) {
    by_beta<ensemble<std::size_t>> i_walker(beta_size, by_beta<std::size_t>(n_walkers));
    by_beta<ensemble<fraction_mcmc<Parameters>>> walker(
        beta_size, by_beta<fraction_mcmc<Parameters>>(n_walkers));
    by_beta<emcee_Step_statistics> emcee_stat(beta_size, emcee_Step_statistics{});
    by_beta<Thermo_Jump_statistics> thermo_stat(beta_size - 1, Thermo_Jump_statistics{});
    auto ff = f.fork(omp_get_max_threads());

    by_beta<std::size_t> i_fraction(beta_size, 0ul);

    auto i_fractions = i_frac_to_i_fracs(i_fraction);

    double sum_prior_logLikelihood = 0;
    std::size_t count_l = 0;
#pragma omp parallel for  // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        for (std::size_t i = 0; i < beta_size; ++i) {
            i_walker[i][iw] = iw + i * n_walkers;
            auto i_th = omp_get_thread_num();
            walker[i][iw] =
                init_fraction_mcmc(ff[i_th], mt[i_th], prior, lik, y, x, i_fractions[i]);
        }
    }
    f += ff;
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        for (std::size_t i = 0; i < beta_size; ++i) {
            if (i_fractions[i][0] == 0) {
                sum_prior_logLikelihood += get<logL>(walker[i][iw].logL(0))();
                ++count_l;
            }
        }
    }
    std::vector<std::size_t> samples_size(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        samples_size[i] = y[i]().size();
    }
    std::size_t max_samples = var::max(samples_size);
    double min_beta = -1.0 / sum_prior_logLikelihood * count_l;
    double log_dgb = -std::log(min_beta) / (beta_size - 2);
    std::vector<double> beta(beta_size);
    beta[0] = 0;
    for (std::size_t i = 1; i + 1 < beta_size; ++i)
        beta[i] = min_beta * std::exp(log_dgb * (i - 1));
    beta[beta_size - 1] = 1;

    std::vector<double> globalbeta(beta_size);
    for (std::size_t i = 0; i < beta_size; ++i) {
        auto i_fr = i_fraction[i];
        auto sample_size0 = i_fr > 0 ? samples_size[i_fr - 1] : 0;
        globalbeta[i] =
            (beta[i] * (samples_size[i_fr] - sample_size0) + sample_size0) / max_samples;
    }
    auto walker_sta = make_logL_statistics(walker);
    return thermo_fraction_mcmc<Parameters>{samples_size, max_samples, i_fraction, i_fractions,
                                            beta,         globalbeta,  walker,     walker_sta,
                                            i_walker,     emcee_stat,  thermo_stat};
}

template <class FunctionTable, class Likelihood, class Variables, class DataType, class Parameters>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
void update_likelihoods_all(FunctionTable& f, thermo_fraction_mcmc<Parameters>& current,
                            Likelihood const& lik, const DataType& y, const Variables& x) {
    auto n_walkers = current.walkers[0].size();
    auto n_beta = current.beta.size();
    auto n_par = current.walkers[0][0].parameter.size();

    auto ff = f.fork(omp_get_max_threads());
    std::size_t num_threads = omp_get_max_threads();
    auto total_walkers = n_beta * n_walkers;
    std::size_t walkers_per_thread = std::ceil(1.0 * total_walkers / num_threads);

#pragma omp parallel for  // collapse(2)
    for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
        for (std::size_t iwb = i_thread * walkers_per_thread;
             iwb < std::min(walkers_per_thread * (i_thread + 1), total_walkers); ++iwb) {
            //  dur.record("begin_loop_walker", iwb * 2);
            //   std::size_t ib = iwb / (n_walkers / 2);
            //   std::size_t i = iwb - ib * n_walkers / 2;
            std::size_t iw = iwb / (n_beta);
            std::size_t ib = iwb - iw * n_beta;

            auto i_th = omp_get_thread_num();
            auto& cu_par = current.walkers[ib][iw].parameter;
            auto Maybe_i = update_logLikelihoods(
                ff[i_th], lik, cu_par.to_value(), y, x, current.i_fractions[ib],
                current.walkers[ib][iw].logL_map, current.walkers_sta[ib][iw]);
            assert(Maybe_i.valid());
        }
    }
    f += ff;
}

template <class FunctionTable, std::size_t N, class Observer, class Prior, class Likelihood,
          class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)
void step_stretch_thermo_fraction_mcmc(FunctionTable& f, std::size_t& iter,
                                       var::Event_Timing<N>& dur,
                                       thermo_fraction_mcmc<Parameters>& current, Observer& obs,
                                       ensemble<mt_64i>& mt, Prior const& prior,
                                       Likelihood const& lik, const DataType& y, const Variables& x,
                                       double alpha_stretch = 2) {
    dur.record("stretch_start");
    auto n_walkers = current.walkers[0].size();
    auto n_beta = current.beta.size();
    auto n_par = current.walkers[0][0].parameter.size();

    std::uniform_int_distribution<std::size_t> uniform_walker(0, n_walkers / 2 - 1);
    std::vector<std::uniform_int_distribution<std::size_t>> udist(omp_get_max_threads(),
                                                                  uniform_walker);

    std::uniform_real_distribution<double> uniform_stretch_zdist(1.0 / alpha_stretch,
                                                                 alpha_stretch);
    std::vector<std::uniform_real_distribution<double>> zdist(omp_get_max_threads(),
                                                              uniform_stretch_zdist);

    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(), uniform_real);

    auto ff = f.fork(omp_get_max_threads());
    std::vector<by_beta<emcee_Step_statistics>> emcee_stat(omp_get_max_threads(),
                                                           by_beta<emcee_Step_statistics>(n_beta));

    dur.record("stretch_before_loop");
    std::size_t num_threads = omp_get_max_threads();

    auto total_walkers_per_half = n_beta * n_walkers / 2;
    std::size_t walkers_per_thread = std::ceil(1.0 * total_walkers_per_half / num_threads);

    // std::size_t n_beta_f = std::ceil(1.0*n_beta / num_threads);

    for (bool half : {false, true}) {
#pragma omp parallel for  // collapse(2)
        for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
            for (std::size_t iwb = i_thread * walkers_per_thread;
                 iwb < std::min(walkers_per_thread * (i_thread + 1), total_walkers_per_half);
                 ++iwb) {
                //  dur.record("begin_loop_walker", iwb * 2);
                //    std::size_t ib = iwb / (n_walkers / 2);
                //    std::size_t i = iwb - ib * n_walkers / 2;

                std::size_t i = iwb / (n_beta);
                std::size_t ib = iwb - i * n_beta;

                auto i_th = omp_get_thread_num();

                auto iw = half ? i + n_walkers / 2 : i;
                auto j = udist[i_th](mt[i_th]);
                auto jw = half ? j : j + n_walkers / 2;
                // we can try in the outer loop

                auto r = rdist[i_th](mt[i_th]);
                auto& cu_par = current.walkers[ib][iw].parameter;
                auto Maybe_i = update_logLikelihoods(
                    ff[i_th], lik, cu_par.to_value(), y, x, current.i_fractions[ib],
                    current.walkers[ib][iw].logL_map, current.walkers_sta[ib][iw]);

                // candidate[ib].walkers[iw].
                auto [ca_par, z] =
                    stretch_move(mt[i_th], rdist[i_th], current.walkers[ib][iw].parameter,
                                 current.walkers[ib][jw].parameter);

                auto ca_logP = logPrior(prior, ca_par);
                auto i_fracts = current.i_fractions[ib];

                auto ca_logL = logLikelihoods(ff[i_th], lik, ca_par.to_value(), y, x, i_fracts);

                if (!((ca_logP.valid()) && (ca_logL.valid()))) {
                    fails(emcee_stat[i_th][ib]());
                } else {
                    auto dthLogL = calc_emcee_jump_delta(current.walkers[ib][iw], ca_logP.value(),
                                                         ca_logL.value(), current.beta[ib],
                                                         current.i_fraction[ib]);
                    auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
                    if (pJump >= r) {
                        current.walkers[ib][iw].parameter = std::move(ca_par);
                        current.walkers[ib][iw].logP = ca_logP.value();
                        current.walkers[ib][iw].logL_map = ca_logL.value();
                        succeeds(emcee_stat[i_th][ib]());
                    } else {
                        fails(emcee_stat[i_th][ib]());
                    }
                }
                //}

                //      dur.record("end_loop_walker", ib * 2 + 1);
            }
        }
    }

    dur.advance(n_beta * 2);
    dur.record("stretch_after_loop");
#pragma omp parallel for
    for (std::size_t ib = 0; ib < n_beta; ++ib) {
        for (auto& e : emcee_stat) current.emcee_stat[ib]() += e[ib]();
    }
    f += ff;
    ++iter;
    update(current.walkers_sta, current.walkers, current.i_fractions);

    dur.record("stretch_function_end");
}

template <class LikelihoodModel, class FuncTable, class Parameters, class Variables, class DataType,
          class Observer>
void thermo_fraction_jump_mcmc(FuncTable& f, const LikelihoodModel& lik,
                               const std::vector<DataType>& y, const std::vector<Variables>& x,
                               std::size_t iter, thermo_fraction_mcmc<Parameters>& current,
                               Observer& obs, mt_64i& mt, ensemble<mt_64i>& mts,
                               std::size_t thermo_jumps_every) {
    if ((iter > 0) && (current.num_samples() > 0) && (iter % (thermo_jumps_every) == 0)) {
        auto ff = f.fork(omp_get_max_threads());

        std::uniform_real_distribution<double> uniform_real(0, 1);
        auto n_walkers = current.get_Walkers_number();
        auto n_beta = current.beta.size();
        auto n_par = current.walkers[0][0].parameter.size();

        WalkerIndexes shuffeld_walkers(n_walkers);
        std::iota(shuffeld_walkers.begin(), shuffeld_walkers.end(), 0);
        std::shuffle(shuffeld_walkers.begin(), shuffeld_walkers.end(), mt);
        std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(),
                                                                  uniform_real);

        std::vector<by_beta<Thermo_Jump_statistics>> thermo_stat(
            omp_get_max_threads(), by_beta<Thermo_Jump_statistics>(n_beta - 1));

        std::size_t num_threads = omp_get_max_threads();

        auto total_walkers_per_half = (n_beta - 1) * n_walkers / 2;
        std::size_t walkers_per_thread = std::ceil(1.0 * total_walkers_per_half / num_threads);

        // std::size_t n_beta_f = std::ceil(1.0*n_beta / num_threads);

#pragma omp parallel for  // collapse(2)
        for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
            for (std::size_t iwb = i_thread * walkers_per_thread;
                 iwb < std::min(walkers_per_thread * (i_thread + 1), total_walkers_per_half);
                 ++iwb) {
                //  dur.record("begin_loop_walker", iwb * 2);
                //                std::size_t ib = iwb / (n_walkers / 2);
                //                std::size_t i = iwb - ib * n_walkers / 2;

                std::size_t i = iwb / (n_beta - 1);
                std::size_t ib = iwb - i * (n_beta - 1);

                auto i_th = omp_get_thread_num();
                auto iw = shuffeld_walkers[i];
                auto jw = shuffeld_walkers[i + n_walkers / 2];

                auto r = rdist[i_th](mts[i_th]);
                auto beta0 = current.beta[ib] == 1 ? 0 : current.beta[ib];
                auto i_frac = current.i_fraction[ib + 1];
                auto Lik_i_0 =
                    i_frac > 0 ? get<logL>(current.walkers[ib][iw].logL(i_frac - 1))() : 0;
                auto Liki = get<logL>(current.walkers[ib][iw].logL(i_frac))() - Lik_i_0;
                auto Lik_j_0 =
                    i_frac > 0 ? get<logL>(current.walkers[ib + 1][jw].logL(i_frac - 1))() : 0;
                auto Likj = get<logL>(current.walkers[ib + 1][jw].logL(i_frac))() - Lik_j_0;

                double logA = calc_logA(beta0, current.beta[ib + 1], Liki, Likj);
                auto pJump = std::min(1.0, std::exp(logA));
                if (pJump > r) {
                    auto Maybe_i = update_logLikelihoods(
                        ff[i_th], lik, current.walkers[ib][iw].parameter.to_value(), x, y,
                        current.i_fractions[ib + 1], current.walkers[ib][iw].logL_map,
                        current.walkers_sta[ib][iw]);

                    auto Maybe_j = update_logLikelihoods(
                        ff[i_th], lik, current.walkers[ib + 1][jw].parameter.to_value(), x, y,
                        current.i_fractions[ib], current.walkers[ib + 1][jw].logL_map,
                        current.walkers_sta[ib + 1][jw]);
                    if (Maybe_i.valid() && Maybe_j.valid()) {
                        std::swap(current.walkers[ib][iw], current.walkers[ib + 1][jw]);
                        std::swap(current.i_walkers[ib][iw], current.i_walkers[ib + 1][jw]);
                        succeeds(thermo_stat[i_th][ib]());
                    } else {
                        fails(thermo_stat[i_th][ib]());
                    }
                } else {
                    fails(thermo_stat[i_th][ib]());
                }
            }
        }

        for (std::size_t ib = 0; ib < n_beta - 1; ++ib) {
            for (auto& e : thermo_stat) current.thermo_stat[ib]() += e[ib]();
        }
    }
}

template <class FunctionTable, class Prior, class Lik, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
Maybe_error<bool> calc_mcmc(FunctionTable& f, Prior const& pr, const Lik& lik, const DataType& y,
                            const Variables& x, fraction_mcmc<Parameters>& t_mcmc,
                            std::vector<std::size_t> const& i_fracs) {
    auto par = t_mcmc.parameter;
    auto logP = logPrior(pr, par);
    auto t_logLs = logLikelihoods(f, lik, par.to_value(), y, x, i_fracs);

    if (logP.valid() && t_logLs.valid()) {
        t_mcmc.logP = logP.value();
        t_mcmc.logL_map = t_logLs.value();
        return true;
    } else {
        t_mcmc.logP = std::numeric_limits<double>::quiet_NaN();
        for (auto& e : t_mcmc.logL_map)
            get<logL>(e.second)() = std::numeric_limits<double>::quiet_NaN();
        ;
    }
    return error_message(logP.error()() + t_logLs.error()());
}

template <bool Adapt_beta, class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Reporter, class mcmc, class Parameters,
          class timepoint>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_fraction_evidence_loop(
    FunctionTable&& f, new_thermodynamic_integration<Algorithm, Reporter>&& therm,
    Prior const& prior, Likelihood const& lik, const std::vector<DataType>& ys,
    const std::vector<Variables>& xs, mcmc mcmc_run, std::size_t iter,
    thermo_fraction_mcmc<Parameters>& current, Reporter& rep, mt_64i& mt, std::vector<mt_64i>& mts,
    const timepoint& start, const std::chrono::duration<double>& previous_duration) {
    var::Event_Timing<200> even_dur(start);
    std::ofstream event_file(f.file() + "_event_timing.csv");

    while (!mcmc_run.second) {
        even_dur.record("main_loop_start");
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start) + previous_duration;

        report_all(f, iter, dur, rep, current, prior, lik, ys, xs, mts, mcmc_run.first);
        if constexpr (Adapt_beta) {
            adapt_fraction_beta(f, iter, current, therm.adapt_beta_every(),
                                therm.adapt_beta_equalizer(), therm.adapt_beta_variance(),
                                therm.desired_acceptance(), therm.adapt_beta_nu(),
                                therm.adapt_beta_t0(), therm.adapt_beta_threshold(), lik, ys, xs);

            if (therm.adjust_beta()) {
                adjust_fraction_beta(f, iter, therm.adapt_beta_every(),
                                     therm.acceptance_upper_limit(), therm.acceptance_lower_limit(),
                                     current, mts, prior, lik, ys, xs);
                adapt_fraction_i_fraction(f, iter, current, therm.adapt_beta_every(),
                                          therm.adapt_beta_equalizer(), therm.adapt_beta_variance(),
                                          therm.desired_acceptance(), therm.adapt_beta_nu(),
                                          therm.adapt_beta_t0(), therm.adapt_beta_threshold(),
                                          therm.acceptance_upper_limit(), lik, ys, xs);
            }
            if (iter % therm.adapt_beta_every() == 0)
                current.reset_statistics();
        }

        step_stretch_thermo_fraction_mcmc(f, iter, even_dur, current, rep, mts, prior, lik, ys, xs);

        even_dur.record("befor_thermo_jump");

        thermo_fraction_jump_mcmc(f, lik, xs, ys, iter, current, rep, mt, mts,
                                  therm.thermo_jumps_every());

        even_dur.record("after_thermo_jump");

        even_dur.record("after_report_all");
        report_point(f, iter);

        // using geg=typename
        // decltype(checks_convergence(std::move(mcmc_run.first), current))::eger;
        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        even_dur.record("after_checks_convergence");
        if (iter == 1)
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

auto thermo_fraction_evidence(FunctionTable&& f,
                              new_thermodynamic_integration<Algorithm, Reporter>&& therm,
                              Prior const& prior, Likelihood const& lik,
                              const std::vector<DataType>& ys, const std::vector<Variables>& xs) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());

    auto current =
        init_thermo_fraction_mcmc(f, n_walkers, therm.beta_size(), mts, prior, lik, ys, xs);
    // auto n_par = current.walkers[0][0].parameter.size();

    auto mcmc_run = checks_convergence(std::move(a), current);

    std::size_t iter = 1;
    const auto start = std::chrono::high_resolution_clock::now();
    auto& rep = therm.reporter();
    report_title(rep, current, lik, ys, xs);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, ys, xs, current.beta);
    std::chrono::duration<double> previous_duration(0.0);
    return thermo_fraction_evidence_loop<Adapt_beta>(
        f, std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm), prior, lik, ys,
        xs, mcmc_run, iter, current, rep, mt, mts, start, previous_duration);
}

template <class Prior, class Parameters = std::decay_t<decltype(sample(std::declval<mt_64i&>(),
                                                                       std::declval<Prior&>()))>>
auto create_thermo_fraction_mcmc(std::size_t n_walkers, by_beta<double> const& beta, mt_64i& mt,
                                 Prior const& pr, std::vector<std::size_t> samples_size,
                                 std::size_t max_samples) {
    auto& priorsampler = pr;
    auto par = sample(mt, priorsampler);

    by_beta<ensemble<std::size_t>> i_walker(beta.size(), by_beta<std::size_t>(n_walkers));
    by_beta<ensemble<fraction_mcmc<Parameters>>> walker(
        beta.size(), by_beta<fraction_mcmc<Parameters>>(n_walkers, fraction_mcmc<Parameters>{par}));
    by_beta<ensemble<std::map<std::size_t, logL_statistics>>> walker_sta(
        beta.size(), by_beta<std::map<std::size_t, logL_statistics>>(n_walkers));

    by_beta<emcee_Step_statistics> emcee_stat(beta.size());
    by_beta<Thermo_Jump_statistics> thermo_stat(beta.size() - 1);

    by_beta<std::size_t> i_fraction(beta.size());
    by_beta<std::vector<std::size_t>> i_fractions(beta.size());
    by_beta<double> S(beta.size());
    by_beta<ensemble<std::map<std::size_t, logL_statistics>>> walkers_sta;

    return thermo_fraction_mcmc<Parameters>{samples_size, max_samples, i_fraction, i_fractions,
                                            beta,         S,           walker,     walker_sta,
                                            i_walker,     emcee_stat,  thermo_stat};
}

template <class Parameters, class Duration>
bool extract_iter(std::istream& f, std::size_t& iter, Duration& dur,
                  thermo_fraction_mcmc<Parameters>& data) {
    std::size_t i_beta = 0;
    std::size_t i_walker = 0;
    std::size_t i_par = 0;

    std::size_t v_i_beta;
    std::size_t v_num_beta = 0;
    double v_beta;
    std::size_t v_i_walker;
    std::size_t v_walker_id;

    std::size_t v_i_par;
    double v_param_value;

    bool not_all = true;
    bool from_begining = false;

    double v_dur;

    while (not_all && load_vars_line(f, iter, v_dur, v_i_beta, v_num_beta, v_beta, v_i_walker,
                                     v_walker_id, v_i_par, v_param_value)) {
        dur = std::chrono::duration<double>(v_dur);
        if (!from_begining) {
            from_begining = (v_i_beta == 0) && (v_i_walker == 0) && (v_i_par == 0);
        }
        if (from_begining) {
            if ((v_i_beta < data.beta.size())) {
                data.beta[v_i_beta] = v_beta;
                data.i_walkers[v_i_beta][v_i_walker] = v_walker_id;
                data.walkers[v_i_beta][v_i_walker].parameter[v_i_par] = v_param_value;
            } else if (v_i_beta == data.beta.size()) {
                data.beta.push_back(v_beta);
                data.i_walkers.push_back(data.i_walkers[v_i_beta - 1]);
                data.i_walkers[v_i_beta][v_i_walker] = v_walker_id;
                data.walkers.push_back(data.walkers[v_i_beta - 1]);
                data.walkers[v_i_beta][v_i_walker].parameter[v_i_par] = v_param_value;

                data.walkers_sta.push_back(data.walkers_sta[v_i_beta - 1]);
                data.emcee_stat.push_back(data.emcee_stat[v_i_beta - 1]);
                data.thermo_stat.push_back(data.thermo_stat[v_i_beta - 2]);

                i_beta = v_i_beta;
            }
            if (v_i_par != i_par)
                return false;
            if (v_i_walker != i_walker)
                return false;
            if (v_i_beta != i_beta)
                return false;

            if (i_par + 1 < num_Parameters(data)) {
                ++i_par;
            } else {
                i_par = 0;
                if (i_walker + 1 < data.i_walkers[i_beta].size()) {
                    ++i_walker;
                } else {
                    i_walker = 0;
                    if (i_beta + 1 < data.walkers.size()) {
                        ++i_beta;
                    } else {
                        i_beta = 0;
                        not_all = !(v_num_beta == data.walkers.size());
                    }
                }
            }
        }
    }

    return v_num_beta == data.walkers.size();
}

template <class Parameters, class Duration>
bool extract_iter(const std::string line, std::size_t& iter, Duration& dur,
                  thermo_fraction_mcmc<Parameters>& data) {
    std::stringstream ss(line);
    return extract_iter(ss, iter, dur, data);
}

template <class Parameters, class Duration>
auto extract_parameters_last(const std::string& fname, std::size_t& iter, Duration& dur,
                             thermo_fraction_mcmc<Parameters>& data) {
    auto candidate = data;
    auto f = std::ifstream(fname);
    std::string line;
    std::getline(f, line);
    auto iter_prev = iter;
    auto dur_prev = dur;
    while (extract_iter(f, iter_prev, dur_prev, candidate)) {
        iter = iter_prev;
        dur = dur_prev;
        std::swap(data, candidate);
    }
    return data;
}

template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
          class Parameters =
              std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
Maybe_error<bool> calc_thermo_fraction_mcmc_continuation(FunctionTable& f, std::size_t n_walkers,
                                                         ensemble<mt_64i>& mt, Prior const& prior,
                                                         Likelihood const& lik, const DataType& y,
                                                         const Variables& x,
                                                         thermo_fraction_mcmc<Parameters>& t_mcmc) {
    auto ff = f.fork(omp_get_max_threads());

    std::string error;
    bool good = true;
#pragma omp parallel for  // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        for (std::size_t i = 0; i < t_mcmc.beta.size(); ++i) {
            auto i_th = omp_get_thread_num();
            auto res =
                calc_mcmc(ff[i_th], prior, lik, y, x, t_mcmc.walkers[i][iw], t_mcmc.i_fractions[i]);
            if (!res) {
                error += res.error()();
                good = false;
            }
        }
    }
    f += ff;
    if (good)
        return good;
    else
        return error_message(error);
}

template <bool Adapt_beta, class FunctionTable, class Algorithm, class Prior, class Likelihood,
          class Variables, class DataType, class Reporter>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_fraction_evidence_continuation(
    const std::string& idName, FunctionTable& f,
    new_thermodynamic_integration<Algorithm, Reporter>&& therm, Prior const& prior,
    Likelihood const& lik, const std::vector<DataType>& y, const std::vector<Variables>& x) {
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

    std::vector<std::size_t> samples_size = get_samples_size(y);
    std::size_t max_samples = var::max(samples_size);

    auto current =
        create_thermo_fraction_mcmc(n_walkers, beta, mt, prior, samples_size, max_samples);
    auto& rep = therm.reporter();

    auto fname = idName + "__i_beta__i_walker__i_par.csv";
    std::size_t iter = 0;
    auto start = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::ratio<1>> duration;

    current = extract_parameters_last(fname, iter, duration, current);

    a.reset(iter);
    auto res = calc_thermo_fraction_mcmc_continuation(f, n_walkers, mts, prior, lik, y, x, current);

    auto mcmc_run = checks_convergence(std::move(a), current);
    //using return_type=Maybe_error<decltype(std::pair(std::move(mcmc_run.first), current))>;
    report_title(rep, current, lik, y, x);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta);

    return thermo_fraction_evidence_loop<Adapt_beta>(
        f, std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm), prior, lik, y,
        x, mcmc_run, iter, current, rep, mt, mts, start, duration);
}

#endif  // PARALLEL_TEMPERING_FRACTION_H
