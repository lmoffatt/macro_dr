#ifndef CUEVI_H
#define CUEVI_H
#include <omp.h>

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

using DataIndexes = std::vector<std::size_t>;

inline auto generate_random_Indexes(mt_64i& mt, std::size_t num_samples,
                                    std::size_t min_num_extra_samples, double num_jumps_per_decade,
                                    std::vector<std::size_t> initial_samples = {}) {
    std::size_t num_initial_samples = size(initial_samples);
    std::size_t n_jumps =
        std::max(0.0, std::floor(num_jumps_per_decade *
                                 (std::log10(num_samples) -
                                  std::log10(min_num_extra_samples + num_initial_samples))));
    auto indexsizes = DataIndexes(n_jumps + 1);

    for (std::size_t i = 0; i < n_jumps + 1; ++i)
        indexsizes[i] = num_samples * std::pow(10.0, -(1.0 * (n_jumps - i)) / num_jumps_per_decade);
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
        it = randomly_extract_n(mt, it, index.end(), indexsizes[0] - num_initial_samples);
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

class LogPrior : public var::Var<LogPrior, double> {};

class LogLik_value : public var::Var<LogLik_value, logLs> {
   public:
    friend std::ostream& operator<<(std::ostream& os, const LogLik_value& x) {
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
    Fraction_Index& operator++() {
        ++((*this)());
        return *this;
    }
};

class Cuevi_Index : public var::Constant<Cuevi_Index, std::size_t> {
   public:
    Cuevi_Index(std::size_t i) : var::Constant<Cuevi_Index, std::size_t>(i) {
    }
    Cuevi_Index friend operator+(Cuevi_Index one, std::size_t i) {
        return Cuevi_Index(one() + i);
    }
    Cuevi_Index friend operator-(Cuevi_Index one, std::size_t i) {
        return Cuevi_Index(one() - i);
    }
    Cuevi_Index& operator=(std::size_t i) {
        (*this)() = i;
        return *this;
    }
    Cuevi_Index& operator++() {
        ++(*this)();
        return *this;
    }
    friend bool operator<(const Cuevi_Index& one, std::size_t n) {
        return one() < n;
    }
};

class Number_of_Fractions : public var::Constant<Number_of_Fractions, std::size_t> {};

class Number_of_samples : public var::Constant<Number_of_samples, std::size_t> {};

class Walker_id : public var::Constant<Walker_id, std::size_t> {};

class Walker_Index : public var::Constant<Walker_Index, std::size_t> {
   public:
    Walker_Index& operator=(std::size_t i) {
        (*this)() = i;
        return *this;
    }
    Walker_Index& operator++() {
        ++(*this)();
        return *this;
    }
    friend bool operator<(const Walker_Index& one, std::size_t n) {
        return one() < n;
    }
};

class LogLik_by_Fraction
    : public var::Var<LogLik_by_Fraction, std::map<Fraction_Index, LogLik_value>> {
   public:
    bool has(Fraction_Index i_fra) const {
        return (*this)().find(i_fra) != (*this)().end();
    }

    Maybe_error<logLs> operator[](Fraction_Index i) const {
        auto it = (*this)().find(i);
        if (it != (*this)().end())
            return it->second();
        else
            return error_message("");
    }
};

class Num_Walkers_Per_Ensemble : public var::Constant<Num_Walkers_Per_Ensemble, std::size_t> {};

class Points_per_decade : public var::Constant<Points_per_decade, double> {};
class Points_per_decade_low : public var::Constant<Points_per_decade_low, double> {};
class Min_value : public var::Constant<Min_value, double> {};
class Med_value : public var::Constant<Med_value, double> {};
class Includes_zero : public var::Constant<Includes_zero, bool> {};

class Number_trials_until_give_up : public var::Constant<Number_trials_until_give_up, std::size_t> {
};

class Random_jumps : public var::Constant<Random_jumps, bool> {};

class Thermo_Jumps_every : public var::Constant<Thermo_Jumps_every, std::size_t> {
    friend void report_model(save_Evidence& s, Thermo_Jumps_every n) {
        std::ofstream f(s.fname + "_thermo_jumps_every");
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << n() << "\n";
    }
};

class Th_Beta_Param
    : public var::Constant<Th_Beta_Param,
                           var::Vector_Space<Includes_zero, Med_value, Points_per_decade, Min_value,
                                             Points_per_decade_low>> {
   public:
    using var::Constant<Th_Beta_Param,
                        var::Vector_Space<Includes_zero, Med_value, Points_per_decade, Min_value,
                                          Points_per_decade_low>>::Constant;
};
class Fractions_Param
    : public var::Constant<Fractions_Param, var::Vector_Space<Min_value, Points_per_decade>> {};

class Cuevi_temperatures_Param
    : public var::Constant<Cuevi_temperatures_Param,
                           var::Vector_Space<Th_Beta_Param, Fractions_Param>> {};

class Prior_statistics : public var::Constant<Prior_statistics, Trial_statistics> {};

class Likelihood_statistics
    : public var::Constant<Likelihood_statistics, std::map<Fraction_Index, Trial_statistics>> {
   public:
    Likelihood_statistics& operator+=(const Likelihood_statistics& other) {
        for (auto& e : other()) {
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

class Cuevi_Jump_statistics : public var::Constant<Cuevi_Jump_statistics, Trial_statistics> {};

class Walker_statistics
    : public var::Var<
          Walker_statistics,
          var::Vector_Space<Prior_statistics, Likelihood_statistics, emcee_Step_statistics,
                            Thermo_Jump_statistics, Cuevi_Jump_statistics>> {
   public:
    void reset() {
        get<Prior_statistics>((*this)())().reset();
        get<Likelihood_statistics>((*this)())().clear();
        get<emcee_Step_statistics>((*this)())().reset();
        get<Thermo_Jump_statistics>((*this)())().reset();
        get<Cuevi_Jump_statistics>((*this)())().reset();
    }

    Walker_statistics& operator+=(const Walker_statistics& other) {
        get<Prior_statistics>((*this)())() += get<Prior_statistics>(other())();
        get<Likelihood_statistics>((*this)()) += get<Likelihood_statistics>(other());
        get<emcee_Step_statistics>((*this)())() += get<emcee_Step_statistics>(other())();
        get<Thermo_Jump_statistics>((*this)())() += get<Thermo_Jump_statistics>(other())();
        get<Cuevi_Jump_statistics>((*this)())() += get<Cuevi_Jump_statistics>(other())();
        return *this;
    }
};

using Walker_statistics_pair = std::pair<Walker_statistics&, Walker_statistics&>;

class Init_seed : public var::Constant<Init_seed, typename mt_64i::result_type> {};
template <class T>
using by_fraction = std::vector<T>;

class Cuevi_temperatures
    : public var::Var<Cuevi_temperatures, std::vector<var::Vector_Space<Th_Beta, Fraction_Index>>> {
};

class Parameter {};
class step_stretch_cuevi_mcmc;
class step_stretch_cuevi_mcmc_per_walker;
class thermo_cuevi_jump_mcmc;

class Cuevi_statistics
    : public var::Var<Cuevi_statistics, ensemble<std::vector<Walker_statistics>>> {
   public:
    auto operator[](Cuevi_Index icu) {
        Walker_statistics out = {};
        for (std::size_t i = 0; i < (*this)().size(); ++i) {
            out += (*this)()[i][icu()];
        }
        return out;
    }
};

template <class ParameterType>
class Cuevi_mcmc;

template <class ParameterType>
class Cuevi_mcmc {
    using myParameter = var::Var<Parameter, ParameterType>;

    class Walker_value
        : public var::Var<Walker_value, var::Vector_Space<var::Var<Parameter, ParameterType>,
                                                          LogPrior, LogLik_by_Fraction>> {
        friend Maybe_error<double> thermo_step(const Walker_value& candidate,
                                               const Walker_value& current, Th_Beta beta,
                                               Fraction_Index i_fra) {
            auto th_ca = get<LogLik_by_Fraction>(candidate())[i_fra];
            auto th_cu = get<LogLik_by_Fraction>(current())[i_fra];

            if (!th_cu.valid() && th_ca.valid())
                return 100.0;
            else
                return get<LogPrior>(candidate())() + beta() * get<logL>(th_ca.value())() -
                       get<LogPrior>(current())() - beta() * get<logL>(th_cu.value())();
        }
        friend Maybe_error<double> thermo_jump(Th_Beta ca_beta, Fraction_Index ca_fra,
                                               const Walker_value& candidate, Th_Beta cu_beta,
                                               Fraction_Index cu_fra, const Walker_value& current) {
            auto logL_ca_ca = get<LogLik_by_Fraction>(candidate())[ca_fra];
            auto logL_cu_cu = get<LogLik_by_Fraction>(current())[cu_fra];

            auto logL_ca_cu = get<LogLik_by_Fraction>(candidate())[cu_fra];
            auto logL_cu_ca = get<LogLik_by_Fraction>(current())[ca_fra];

            auto current_sum = ca_beta() * logL_ca_ca + cu_beta() * logL_cu_cu;
            auto after_jump_sum = cu_beta() * logL_ca_cu + ca_beta() * logL_cu_ca;

            bool current_is_big_0_is_NaN =
                (cu_fra() > ca_fra()) && (!logL_cu_ca.valid()) && (logL_cu_ca.valid());
            bool candidate_is_big_0_is_NaN =
                ((cu_fra() < ca_fra()) && (logL_cu_ca.valid()) && (!logL_cu_ca.valid()));
            bool current_sum_is_NaN_after_is_not = (!current_sum.valid() && after_jump_sum.valid());

            if (current_is_big_0_is_NaN || candidate_is_big_0_is_NaN ||
                current_sum_is_NaN_after_is_not)
                return 0.0;  // takes exponential---> prob 1
            else
                return get<logL>((after_jump_sum - current_sum).value())();
        }

        friend Maybe_error<double> logEvidence_walker_pos(const Walker_value& w, Cuevi_Index i,
                                                          const Cuevi_temperatures& t) {
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

        friend auto logEvidence_walker(const Walker_value& w, Cuevi_Index i,
                                       const Cuevi_temperatures& t) {
            if (i() == 0)
                return logEvidence_walker_pos(w, i, t);
            else
                return logEvidence_walker_pos(w, i, t) + logEvidence_walker_pos(w, i - 1, t);
        }

        friend auto cuevi_jump(const Walker_value& w_0, Cuevi_Index i_0, const Walker_value& w_1,
                               Cuevi_Index i_1, const Cuevi_temperatures& t) {
            auto current_sum = logEvidence_walker(w_0, i_0, t) + logEvidence_walker(w_1, i_1, t);
            auto after_jump_sum = logEvidence_walker(w_0, i_1, t) + logEvidence_walker(w_1, i_0, t);
            return (after_jump_sum - current_sum);
        }
    };

    class Walker
        : public var::Var<Walker, var::Vector_Space<Walker_id, Walker_value, Walker_statistics>> {};

    class Walkers_ensemble : public var::Var<Walkers_ensemble, ensemble<std::vector<Walker>>> {};

    template <class FunctionTable, class t_logLikelihood, class Data, class Variables>
    static Maybe_error<bool> calc_Likelihood(FunctionTable& f, t_logLikelihood&& lik,
                                             Walker_value& w, const by_fraction<Data>& y,
                                             const by_fraction<Variables>& x, Fraction_Index i_frac,
                                             Walker_statistics_pair wa_sta) {
        assert(i_frac() < size(y));

        auto& r_logLikf = get<LogLik_by_Fraction>(w());
        if (r_logLikf.has(i_frac))
            return true;
        else {
            auto& ca_par = get<Parameter>(w());
            auto v_logL =
                f.f(logLikelihood_f{}, lik, ca_par().to_value(), y[i_frac()], x[i_frac()]);
            if (!v_logL) {
                fails(get<Likelihood_statistics>(wa_sta.first())()[i_frac]);
                fails(get<Likelihood_statistics>(wa_sta.second())()[i_frac]);
                // r_logLikf()[i_frac] =
                // LogLik_value(-std::numeric_limits<double>::infinity());
                return v_logL.error();

            } else {
                succeeds(get<Likelihood_statistics>(wa_sta.first())()[i_frac]);
                succeeds(get<Likelihood_statistics>(wa_sta.second())()[i_frac]);
                r_logLikf()[i_frac] = LogLik_value(v_logL.value());
                return true;
            }
        }
    }

    template <class FunctionTable, class t_logLikelihood, class Data, class Variables>
    static Maybe_error<bool> calc_Relevant_Likelihoods(FunctionTable& f, t_logLikelihood&& lik,
                                                       const by_fraction<Data>& y,
                                                       const by_fraction<Variables>& x,
                                                       const Cuevi_temperatures& t, Walker_value& w,
                                                       Cuevi_Index i,
                                                       Walker_statistics_pair wa_sta) {
        Maybe_error<bool> out = true;
        for (auto i_frac = get<Fraction_Index>(t()[std::max(1ul, i()) - 1]);
             i_frac < get<Fraction_Index>(t()[std::min(t().size() - 1, i() + 1)]); ++i_frac) {
            auto Maybe = calc_Likelihood(f, lik, w, y, x, i_frac, wa_sta);
            if (!Maybe)
                return Maybe.error();
        }
        return true;
    }

    template <class FunctionTable, class logLikelihood, class Data, class Variables>
    static Maybe_error<double> thermo_jump_logProb(
        FunctionTable& f, logLikelihood&& lik, const by_fraction<Data>& y,
        const by_fraction<Variables>& x, const Cuevi_temperatures& t, Cuevi_Index i_0,
        Walker_value& w_0, Cuevi_Index i_1, Walker_value& w_1, Walker_statistics_pair wa_sta_0,
        Walker_statistics_pair wa_sta_1) {
        auto i_frac_0 = get<Fraction_Index>(t()[i_0()]);
        auto i_frac_1 = get<Fraction_Index>(t()[i_1()]);

        auto Maybe_0 = calc_Likelihood(f, lik, w_0, y, x, i_frac_1, wa_sta_0);
        auto Maybe_1 = calc_Likelihood(f, lik, w_1, y, x, i_frac_0, wa_sta_1);
        auto be_0 = get<Th_Beta>(t()[i_0()]);
        auto be_1 = get<Th_Beta>(t()[i_1()]);

        // if (!Maybe_0.valid() || !Maybe_1.valid())
        //   return error_message(Maybe_0.error()() + Maybe_1.error()());
        // else {
        return thermo_jump(be_0, i_frac_0, w_0, be_1, i_frac_1, w_1);
        // }
    }

    template <class FunctionTable, class logLikelihood, class Data, class Variables>
    static Maybe_error<double> cuevi_jump_logProb(
        FunctionTable& f, logLikelihood&& lik, const by_fraction<Data>& y,
        const by_fraction<Variables>& x, const Cuevi_temperatures& t, Cuevi_Index i_0,
        Walker_value& w_0, Cuevi_Index i_1, Walker_value& w_1, Walker_statistics_pair wa_sta_0,
        Walker_statistics_pair wa_sta_1) {
        auto Maybe_00 = calc_Relevant_Likelihoods(f, lik, y, x, t, w_0, i_0, wa_sta_0);
        auto Maybe_01 = calc_Relevant_Likelihoods(f, lik, y, x, t, w_0, i_1, wa_sta_0);
        auto Maybe_10 = calc_Relevant_Likelihoods(f, lik, y, x, t, w_1, i_0, wa_sta_1);
        auto Maybe_11 = calc_Relevant_Likelihoods(f, lik, y, x, t, w_1, i_1, wa_sta_1);

        if (!Maybe_00.valid() || !Maybe_01.valid() || !Maybe_10.valid() || !Maybe_11.valid())
            return error_message(Maybe_00.error()() + Maybe_01.error()() + Maybe_10.error()() +
                                 Maybe_11.error()());
        else {
            return cuevi_jump(w_0, i_0, w_1, i_1, t);
        }
    }

    template <class FunctionTable, class Prior, class logLikelihood, class Data, class Variables>
    static Maybe_error<Walker_value> calc_Walker_value(
        FunctionTable& f, myParameter&& ca_par, Prior&& p, logLikelihood&& lik,
        const by_fraction<Data>& y, const by_fraction<Variables>& x, const Cuevi_temperatures& t,
        Cuevi_Index i_cu, Walker_statistics_pair wa_sta) {
        auto v_logP = logPrior(p, ca_par());
        if (!v_logP) {
            fails(get<Prior_statistics>(wa_sta.first())());
            fails(get<Prior_statistics>(wa_sta.second())());
            return v_logP.error();
        } else {
            succeeds(get<Prior_statistics>(wa_sta.first())());
            succeeds(get<Prior_statistics>(wa_sta.second())());
            Walker_value out(var::Vector_Space(std::move(ca_par), LogPrior(v_logP.value()),
                                               LogLik_by_Fraction{}));

            auto i_frac = get<Fraction_Index>(t()[i_cu()]);
            auto Maybe_succeed = calc_Likelihood(f, lik, out, y, x, i_frac, wa_sta);
            if (!Maybe_succeed)
                return Maybe_succeed.error();
            else
                return out;
        }
    }

    template <class FunctionTable, class Sampler, class Prior, class logLikelihood, class Data,
              class Variables>
    static Maybe_error<Walker_value> sample_Walker_for(
        FunctionTable& f, mt_64i& mt, Sampler&& sampler, Prior&& p, logLikelihood&& lik,
        const by_fraction<Data>& y, const by_fraction<Variables>& x, const Cuevi_temperatures& t,
        Cuevi_Index i_cu, Walker_statistics_pair wa_sta, Number_trials_until_give_up max_trials) {
        assert(i_cu() < t().size());

        Maybe_error<Walker_value> v_walker(error_message{});
        auto n_trial = 0ul;
        while (!v_walker && n_trial < max_trials) {
            auto ca_par = std::forward<Sampler>(sampler)(mt);
            v_walker = calc_Walker_value(f, ca_par, p, lik, y, x, t, i_cu, wa_sta);
            ++n_trial;
        }
        if (v_walker.valid())
            return v_walker.value();
        else
            return error_message("more than " + std::to_string(max_trials()) +
                                 " and not a single valid sample, last error " +
                                 v_walker.error()());
    }

    template <class FunctionTable, class Prior, class logLikelihood, class Data, class Variables>
    static auto sample_Walker(FunctionTable& f, mt_64i& mt, Prior&& p, logLikelihood&& lik,
                              const by_fraction<Data>& y, const by_fraction<Variables>& x,
                              const Cuevi_temperatures& t, Cuevi_Index i_cu,
                              Walker_statistics_pair w_sta,
                              Number_trials_until_give_up max_trials) {
        return sample_Walker_for(
            f, mt, [&p](mt_64i& mt) { return sample(mt, std::forward<Prior>(p)); }, p, lik, y, x, t,
            i_cu, w_sta, max_trials);
    }

    template <class FunctionTable, class Prior, class logLikelihood, class Data, class Variables>
    static Maybe_error<typename Cuevi_mcmc<ParameterType>::Walkers_ensemble> init_walkers_ensemble(
        FunctionTable& f, ensemble<mt_64i>& mts, Prior&& p, logLikelihood&& lik,
        const by_fraction<Data>& y, const by_fraction<Variables>& x, const Cuevi_temperatures& t,
        Num_Walkers_Per_Ensemble n, Number_trials_until_give_up max_trials_per_sample) {
        auto num_temp = t().size();
        using Walkers_ensemble = typename Cuevi_mcmc<ParameterType>::Walkers_ensemble;
        using Walker = typename Cuevi_mcmc<ParameterType>::Walker;
        using Walker_value = typename Cuevi_mcmc<ParameterType>::Walker_value;
        Cuevi_statistics sta(ensemble<std::vector<Walker_statistics>>(
            n(), std::vector<Walker_statistics>(num_temp)));
        ;
        Walkers_ensemble out(std::vector(n(), std::vector<Walker>(num_temp)));

        auto ff = f.fork(omp_get_max_threads());
        std::vector<Maybe_error<bool>> succeeds(omp_get_max_threads());
#pragma omp parallel for  //collapse(2)
        for (std::size_t iw = 0; iw < n(); ++iw) {
            for (std::size_t i_cu = 0; i_cu < num_temp; ++i_cu) {
                auto i_th = omp_get_thread_num();
                auto& wa_va = out()[iw][i_cu];
                Walker_statistics_pair wa_sta(get<Walker_statistics>(wa_va()), sta()[iw][i_cu]);
                get<Walker_id>(wa_va())() = iw + num_temp * i_cu;
                auto Maybe_Walker_value =
                    sample_Walker(ff[i_th], mts[i_th], std::forward<Prior>(p), lik, y, x, t,
                                  Cuevi_Index(i_cu), wa_sta, max_trials_per_sample);
                if (Maybe_Walker_value) {
                    succeeds[i_th] = true;
                    get<Walker_value>(wa_va()) = std::move(Maybe_Walker_value.value());
                } else
                    succeeds[i_th] =
                        error_message(succeeds[i_th].error()() + Maybe_Walker_value.error()());
            }
        }
        f += ff;
        auto success = consolidate(succeeds);

        if (success)
            return out;
        else
            return success.error();
    }

    Cuevi_temperatures m_temperatures;
    Cuevi_statistics m_sta;
    Walkers_ensemble m_data;
    std::size_t m_max_i_frac;

    Cuevi_mcmc(Cuevi_temperatures&& t_temperatures, Cuevi_statistics&& t_sta,
               Walkers_ensemble&& t_data, std::size_t t_max_i_frac)
        : m_temperatures{std::move(t_temperatures)},
          m_sta{std::move(t_sta)},
          m_data{std::move(t_data)},
          m_max_i_frac{t_max_i_frac} {
    }

   public:
    auto calc_Mean_logLik(Cuevi_Index j) {
        assert(j() < get_Cuevi_Temperatures_Number());
        auto i_frac = get_Fraction(j);
        return foldMap(
                   make_Range(Walker_Index(0ul), Walker_Index(get_Walkers_number())),
                   [j, i_frac, this](auto i_w) {
                       return get<LogLik_by_Fraction>(get_Walker_Value(i_w, j)())[i_frac];
                   },
                   [](auto x, auto y) { return x + y; }) /
               get_Walkers_number();
    }

    Maybe_error<logLs> calc_Mean_logLik_0(Cuevi_Index j) {
        assert(j() < get_Cuevi_Temperatures_Number());

        auto i_frac = get_Fraction(j);
        if (i_frac() == 0)
            return error_message();
        else
            return foldMap(
                       make_Range(Walker_Index(0ul), Walker_Index(get_Walkers_number())),
                       [j, i_frac, this](auto i_w) {
                           return get<LogLik_by_Fraction>(get_Walker_Value(i_w, j)())[i_frac - 1];
                       },
                       [](auto x, auto y) { return x + y; }) /
                   get_Walkers_number();
    }

    Maybe_error<logLs> calc_Mean_logLik_2(Cuevi_Index j) {
        assert(j() < get_Cuevi_Temperatures_Number());
        auto i_frac = get_Fraction(j);
        auto beta = get_Beta(j);
        if (beta < 1.0)
            return error_message();
        else if (i_frac() + 1 >= m_max_i_frac)
            return error_message();
        else
            return foldMap(
                       make_Range(Walker_Index(0ul), Walker_Index(get_Walkers_number())),
                       [j, i_frac, this](auto i_w) {
                           return get<LogLik_by_Fraction>(get_Walker_Value(i_w, j)())[i_frac + 1];
                       },
                       [](auto x, auto y) { return x + y; }) /
                   get_Walkers_number();
    }

    auto calc_Mean_logPrior(Cuevi_Index j) {
        assert(j() < get_Cuevi_Temperatures_Number());
        return foldMap(
                   make_Range(Walker_Index(0ul), Walker_Index(get_Walkers_number())),
                   [j, this](auto i_w) { return get<LogPrior>(get_Walker_Value(i_w, j)())(); },
                   [](auto x, auto y) { return x + y; }) /
               get_Walkers_number();
    }

    auto& get_Walker_Value(Walker_Index i, Cuevi_Index j) {
        assert(i() < get_Walkers_number());
        assert(j() < get_Cuevi_Temperatures_Number());
        return get<Walker_value>(m_data()[i()][j()]());
    }
    auto& get_Walker_Value(Walker_Index i, Cuevi_Index j) const {
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

    auto& get_Cuevi_Temperatures() const {
        return m_temperatures;
    }

    auto get_Cuevi_Temperatures_Number() const {
        return m_temperatures().size();
    }

    std::size_t get_Walkers_number() const {
        return m_data().size();
    }
    auto get_Parameters_number() const {
        return size(get<Parameter>(get<Walker_value>(m_data()[0][0]())())());
    }
    auto get_Parameter(Walker_Index i, Cuevi_Index j) const {
        assert(i() < get_Walkers_number());
        assert(j() < get_Cuevi_Temperatures_Number());
        return get<Parameter>(get_Walker_Value(i, j)())();
    }
    auto& get_Walker(Walker_Index i, Cuevi_Index j) {
        assert(i() < get_Walkers_number());
        assert(j() < get_Cuevi_Temperatures_Number());

        return m_data()[i()][j()];
    }
    auto& get_Walker(Walker_Index i, Cuevi_Index j) const {
        return m_data()[i()][j()];
    }

    auto get_Beta(Cuevi_Index i_cu) {
        assert(i_cu() < get_Cuevi_Temperatures_Number());
        return get<Th_Beta>(m_temperatures()[i_cu()]);
    }

    auto get_Cuevi_Statistics(Cuevi_Index i_cu) {
        return m_sta[i_cu];
    }

    auto get_Fraction(Cuevi_Index i_cu) {
        assert(i_cu() < get_Cuevi_Temperatures_Number());

        return get<Fraction_Index>(m_temperatures()[i_cu()]);
    }

    template <class DataType>
    static Cuevi_temperatures build_temperatures(const by_fraction<DataType>& y,
                                                 const Th_Beta_Param& p) {
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
        auto beta_mid = std::pow(10.0, -1.0 * num_beta_high / n_points_per_decade_beta());

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
                beta_mid *
                std::pow(10.0, -1.0 * (num_beta_low - i) / n_points_per_decade_beta_low());
            get<Fraction_Index>(out()[i])() = 0ul;
        }
        for (std::size_t i = 0; i < num_beta_high + 1; ++i) {
            get<Th_Beta>(out()[i + num_beta_low])() =
                std::pow(10.0, -1.0 * (num_beta_high - i) / n_points_per_decade_beta());
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

    template <class FunctionTable, class Prior, class logLikelihood, class Data, class Variables>
    Maybe_error<Cuevi_mcmc> static init(FunctionTable& f, ensemble<mt_64i>& mts, Prior&& prior,
                                        logLikelihood&& lik, const by_fraction<Data>& y,
                                        const by_fraction<Variables>& x, const Th_Beta_Param& beta,
                                        Num_Walkers_Per_Ensemble num_walkers,
                                        Number_trials_until_give_up max_trials_per_sample) {
        using Cumc = Cuevi_mcmc<ParameterType>;
        auto t = Cumc::build_temperatures(y, beta);
        auto sta = Cuevi_statistics(
            std::vector(num_walkers(), std::vector<Walker_statistics>(t().size())));
        auto Maybe_walkers = Cumc::init_walkers_ensemble(f, mts, std::forward<Prior>(prior), lik, y,
                                                         x, t, num_walkers, max_trials_per_sample);
        if (!Maybe_walkers)
            return Maybe_walkers.error();
        else
            return Cumc(std::move(t), std::move(sta), std::move(Maybe_walkers.value()), size(y));
    }

    template <class FunctionTable, class t_logLikelihood, class Data, class Variables>
    void calculate_Likelihoods_for_Evidence_calulation(FunctionTable& f, t_logLikelihood&& lik,
                                                       const by_fraction<Data>& y,
                                                       const by_fraction<Variables>& x) {
        auto ff = f.fork(omp_get_max_threads());
        std::size_t num_walk = this->get_Walkers_number();

#pragma omp parallel for  // collapse(2)
        for (std::size_t iwa = 0; iwa < num_walk; ++iwa) {
            for (std::size_t icu = 0; icu < this->get_Cuevi_Temperatures_Number(); ++icu) {
                Cuevi_Index i_cu = Cuevi_Index(icu);
                auto i_th = omp_get_thread_num();
                Walker_Index i_w(iwa);
                auto i_frac = this->get_Fraction(i_cu)();
                auto beta = this->get_Beta(i_cu)();
                if (beta == 1) {
                    Walker_value& wa = this->get_Walker_Value(i_w, i_cu);
                    Walker_statistics_pair wa_sta = this->get_Walker_Statistics(i_w, i_cu);
                    for (std::size_t i = std::max(1ul, i_frac) - 1;
                         i < std::min(i_frac + 2, m_max_i_frac); ++i) {
                        Fraction_Index ifrac = Fraction_Index(i);
                        calc_Likelihood(ff[i_th], lik, wa, y, x, ifrac, wa_sta);
                    }
                }
            }
        }
        f += ff;
    }

    template <class FunctionTable, class Duration, class Prior, class t_logLikelihood, class Data,
              class Variables>
    friend void report(FunctionTable& f, std::size_t iter, const Duration& dur,
                       save_likelihood<ParameterType>& s, Cuevi_mcmc& data, Prior&&,
                       t_logLikelihood&& lik, const by_fraction<Data>& y,
                       const by_fraction<Variables>& x, ...) {
        auto& t = data.get_Cuevi_Temperatures();
        auto num_states = lik.m.number_of_states();
        auto num_fr = data.get_Cuevi_Temperatures_Number();
        std::size_t num_values = 32;
        std::size_t point_size = num_values * data.get_Walkers_number() * num_fr;
        std::size_t sampling_interval =
            std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

        if ((iter > 0) && (iter % sampling_interval == 0)) {
            data.calculate_Likelihoods_for_Evidence_calulation(f, lik, y, x);
            for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number(); ++i_walker) {
                auto iw = Walker_Index(i_walker);
                auto& wav = data.get_Walker_Value(iw, Cuevi_Index(0ul));
                auto logL1 = get<LogLik_by_Fraction>(wav())[Fraction_Index(0)];

                using Maybe_L = decltype(logL1);
                auto beta1 = 0.0;
                double beta0 = 0.0;
                Maybe_L log_Evidence = {};
                Maybe_L log_Evidence_no_0 = {};

                for (Cuevi_Index i_cu = 0ul; i_cu() < data.get_Cuevi_Temperatures_Number();
                     ++i_cu) {
                    auto i_fra = data.get_Fraction(i_cu);
                    auto nsamples = size(y[i_fra()]);
                    auto icu = Cuevi_Index(i_cu);
                    auto& wa = data.get_Walker(iw, icu);
                    auto& wav = data.get_Walker_Value(iw, icu);
                    auto& stat = get<Walker_statistics>(wa());

                    auto& prior_sta = get<Prior_statistics>(stat());
                    auto lik_sta_2 = get<Likelihood_statistics>(stat())[i_fra + 1];
                    auto lik_sta_1 = get<Likelihood_statistics>(stat())[i_fra];
                    auto lik_sta_0 = get<Likelihood_statistics>(stat())[i_fra - 1];

                    auto& emcee_sta = get<emcee_Step_statistics>(stat());
                    auto& th_sta = get<Thermo_Jump_statistics>(stat());
                    auto& cuevi_sta = get<Cuevi_Jump_statistics>(stat());

                    auto logPrior = get<LogPrior>(wav());
                    Maybe_L logL0;
                    Maybe_L plog_Evidence;
                    Maybe_L logL1_1 = error_message{};
                    Maybe_L logL1_0 = error_message{};
                    Maybe_L logL1_2 = error_message{};
                    Maybe_L logL0_1 = error_message{};
                    Maybe_L logL0_0 = error_message{};
                    if (i_fra() == 0) {
                        logL0 = logL1;
                        beta0 = beta1;
                        beta1 = data.get_Beta(i_cu)();
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

                        auto& wav0 = data.get_Walker_Value(iw, icu - 1);
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
                    s.f << iter << s.sep << dur.count() << s.sep << i_cu << s.sep << i_fra()
                        << s.sep << nsamples << s.sep << beta1 << s.sep << beta0 << s.sep
                        << i_walker << s.sep << get<Walker_id>(wa())() << s.sep << logPrior
                        << sep(logL1, s.sep) << s.sep << getv<logL>(logL0)
                        << sep(plog_Evidence, s.sep) << sep(log_Evidence, s.sep) << s.sep
                        << getv<logL>(log_Evidence_no_0) << sep(logL1_1, s.sep) << s.sep
                        << getv<logL>(logL1_0) << s.sep << getv<logL>(logL1_2) << s.sep
                        << getv<logL>(logL0_1) << s.sep << getv<logL>(logL0_0) << s.sep

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

    friend void report_title(save_likelihood<ParameterType>& s, Cuevi_mcmc const&...) {
        s.f << "iter" << s.sep << "iter_time" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep
            << "nsamples" << s.sep << "beta1" << s.sep << "beta0" << s.sep << "i_walker" << s.sep
            << "walker_id" << s.sep << "logPrior" << s.sep << "logL1" << s.sep << "elogL1" << s.sep
            << "vlogL1" << s.sep << "logL0" << s.sep << "plog_Evidence" << s.sep << "eplog_Evidence"
            << s.sep << "vplog_Evidence" << s.sep << "log_Evidence" << s.sep << "elog_Evidence"
            << s.sep << "vlog_Evidence" << s.sep << "log_Evidence_no_0" << s.sep << "logL1_1"
            << s.sep << "elogL1_1" << s.sep << "vlogL1_1" << s.sep << "logL1_0" << s.sep
            << "logL1_2" << s.sep << "logL0_1" << s.sep << "logL0_0"

            << s.sep << "prior_sta_count" << s.sep << "prior_sta_rate" << s.sep << "lik_sta_0_count"
            << s.sep << "lik_sta_0_rate" << s.sep << "lik_sta_1_count" << s.sep << "lik_sta_1_rate"
            << s.sep << "lik_sta_2_count" << s.sep << "lik_sta_2_rate" << s.sep << "emcee_sta_count"
            << s.sep << "emcee_sta_rate" << s.sep << "th_sta_count" << s.sep << "th_sta_rate"
            << s.sep << "cuevi_sta_count" << s.sep << "cuevi_sta_rate"
            << "\n";
    }

    template <class FunctionTable, class Duration, class Prior, class t_logLikelihood, class Data,
              class Variables, class Finalizer>
    friend void report(FunctionTable&, std::size_t iter, const Duration& dur,
                       save_Parameter<ParameterType>& s, Cuevi_mcmc const& data, Prior&&,
                       t_logLikelihood&&, const by_fraction<Data>& y, const by_fraction<Variables>&,
                       const std::vector<mt_64i>& mts, const Finalizer& mcmc, ...) {
        auto& t = data.get_Cuevi_Temperatures();

        auto num_fr = data.get_Cuevi_Temperatures_Number();
        std::size_t num_values = 4;
        std::size_t point_size =
            num_values * data.get_Walkers_number() * num_fr * data.get_Parameters_number();
        std::size_t sampling_interval =
            std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

        if ((iter > 0) && (iter % sampling_interval == 0)) {
            for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
                auto icu = Cuevi_Index(i_cu);
                auto i_fra = get<Fraction_Index>(t()[i_cu]);
                for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number(); ++i_walker) {
                    auto iw = Walker_Index(i_walker);
                    auto& wa = data.get_Walker(iw, icu);
                    auto& wav = data.get_Walker_Value(iw, icu);

                    auto i_mts = iw() % 2;

                    for (std::size_t i_par = 0; i_par < data.get_Parameters_number(); ++i_par) {
                        s.f << iter << s.sep << dur.count() << s.sep << i_cu << s.sep << i_fra()
                            << s.sep << size(y[i_fra()]) << s.sep << get<Th_Beta>(t()[i_cu])()
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

    template <bool verifying, class Prior, class t_logLikelihood, class Data, class Variables>
    static Maybe_error<std::pair<Cuevi_mcmc, std::vector<unsigned long long>>> extract_iter(
        save_Parameter<ParameterType>& s, std::string& last_line, std::size_t& iter0,
        std::size_t& npar, std::size_t& nwalkers, std::size_t& n_cuevi) {
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

        ss >> iter0 >> s.sep >> i_cu0 >> s.sep >> i_fra0 >> s.sep >> nsamples0 >> s.sep >> beta0 >>
            s.sep >> i_walker0 >> s.sep >> walker_id0 >> s.sep >> i_mt0 >> s.sep >> mt_pos0 >>
            s.sep >> i_par0 >> s.sep >> par_value0;

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

                    ss2 >> iter >> s.sep >> i_cu >> s.sep >> i_fra >> s.sep >> nsamples >> s.sep >>
                        beta >> s.sep >> i_walker >> s.sep >> walker_id >> i_mt >> s.sep >>
                        mt_pos >> s.sep >> s.sep >> i_par >> s.sep >> par_value;
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
        Walkers_ensemble data(std::vector(nwalkers, std::vector(n_cuevi, Walker{})));
        for (std::size_t iw = 0; iw < nwalkers; ++iw)
            for (std::size_t icu = 0; icu < n_cuevi; ++icu) data()[iw][icu] = d[icu][iw];
        Cuevi_mcmc out;
        out.m_temperatures = std::move(t);
        out.m_data = std::move(data);
        return std::pair(out, mts_pos);
    }

    template <class mcmcm_type>
    friend void report_title(save_Parameter<ParameterType>& s, Cuevi_mcmc const&,
                             const mcmcm_type& mcmc, ...) {
        s.f << "iter" << s.sep << "iter_time" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep
            << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep << "walker_id" << s.sep
            << "i_mt" << s.sep << "mt_pos" << s.sep << "i_par" << s.sep << "parameter_value";
        report_finalizer_title(s, mcmc);
        s.f << "\n";
    }

    template <class mcmcm_type>
    friend void report_title(save_RateParameter<ParameterType>&, Cuevi_mcmc const&,
                             const mcmcm_type&, ...) {
    }

    template <class FunctionTable, class Duration, class Prior, class t_logLikelihood, class Data,
              class Variables>
    friend void report(FunctionTable& f, std::size_t iter, const Duration& dur, save_Evidence& s,
                       Cuevi_mcmc& data, Prior&&, t_logLikelihood&& lik, const by_fraction<Data>& y,
                       const by_fraction<Variables>& x, ...) {
        auto& t = data.get_Cuevi_Temperatures();
        std::size_t num_values = 32;
        std::size_t point_size = num_values * data.get_Cuevi_Temperatures_Number();
        if ((iter > 0) && (iter % std::max(s.sampling_interval,
                                           point_size / s.max_number_of_values_per_iteration) ==
                           0)) {
            data.calculate_Likelihoods_for_Evidence_calulation(f, lik, y, x);
            auto logL1 = data.calc_Mean_logLik(Cuevi_Index(0ul));
            auto logL1_0 = data.calc_Mean_logLik_0(Cuevi_Index(0ul));
            auto logL1_2 = data.calc_Mean_logLik_2(Cuevi_Index(0ul));

            auto beta1 = 0.0;
            Maybe_error<logLs> log_Evidence = {};
            Maybe_error<logLs> log_Evidence_no_0 = {};

            for (Cuevi_Index i_cu = 0ul; i_cu() < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
                auto i_frac = data.get_Fraction(i_cu);
                auto nsamples = size(y[i_frac()]);
                auto logPrior = data.calc_Mean_logPrior(i_cu);
                double beta0;
                Maybe_error<logLs> logL0;
                Maybe_error<logLs> plog_Evidence;
                Maybe_error<logLs> logL1_1 = error_message{};
                Maybe_error<logLs> logL1_0 = error_message{};
                Maybe_error<logLs> logL0_1 = error_message{};
                Maybe_error<logLs> logL0_0 = error_message{};

                auto stat = data.get_Cuevi_Statistics(i_cu);

                auto& prior_sta = get<Prior_statistics>(stat());
                auto lik_sta_2 = get<Likelihood_statistics>(stat())[i_frac + 1];
                auto lik_sta_1 = get<Likelihood_statistics>(stat())[i_frac];
                auto lik_sta_0 = get<Likelihood_statistics>(stat())[i_frac - 1];

                auto& emcee_sta = get<emcee_Step_statistics>(stat());
                auto& th_sta = get<Thermo_Jump_statistics>(stat());
                auto& cuevi_sta = get<Cuevi_Jump_statistics>(stat());

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
                s.f << iter << s.sep << dur.count() << s.sep << i_cu << s.sep << i_frac() << s.sep
                    << nsamples << s.sep << beta1 << s.sep << beta0 << s.sep << logPrior
                    << sep(logL1, s.sep) << sep(logL0, s.sep) << sep(plog_Evidence, s.sep)
                    << sep(log_Evidence, s.sep) << s.sep << getv<logL>(log_Evidence_no_0)
                    << sep(logL1_1, s.sep) << s.sep << getv<logL>(logL1_0) << s.sep
                    << getv<logL>(logL1_2) << s.sep << getv<logL>(logL0_1) << s.sep
                    << getv<logL>(logL0_0) << s.sep << prior_sta().count() << s.sep
                    << prior_sta().rate() << s.sep << lik_sta_0.count() << s.sep << lik_sta_0.rate()
                    << s.sep << lik_sta_1.count() << s.sep << lik_sta_1.rate() << s.sep
                    << lik_sta_2.count() << s.sep << lik_sta_2.rate() << s.sep
                    << emcee_sta().count() << s.sep << emcee_sta().rate() << s.sep
                    << th_sta().count() << s.sep << th_sta().rate() << s.sep << cuevi_sta().count()
                    << s.sep << cuevi_sta().rate() << "\n";
            }
        }
    }

    friend void report_title(save_Evidence& s, Cuevi_mcmc const&, ...) {
        s.f << "iter" << s.sep << "iter_time" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep
            << "nsamples" << s.sep << "beta1" << s.sep << "beta0" << s.sep << "logPrior" << s.sep

            << "logL1" << s.sep << "elogL1" << s.sep << "vlogL1" << s.sep << "logL0" << s.sep
            << "elogL0" << s.sep << "vlogL0" << s.sep << "plog_Evidence" << s.sep
            << "eplog_Evidence" << s.sep << "vplog_Evidence" << s.sep << "log_Evidence" << s.sep
            << "elog_Evidence" << s.sep << "vlog_Evidence" << s.sep << "log_Evidence_no_0" << s.sep
            << "logL1_1" << s.sep << "elogL1_1" << s.sep << "vlogL1_1" << s.sep

            << "logL1_0" << s.sep << "logL1_2" << s.sep << "logL0_1" << s.sep << "logL0_0" << s.sep
            << "prior_sta_count" << s.sep << "prior_sta_rate" << s.sep << "lik_sta_0_count" << s.sep
            << "lik_sta_0_rate" << s.sep << "lik_sta_1_count" << s.sep << "lik_sta_1_rate" << s.sep
            << "lik_sta_2_count" << s.sep << "lik_sta_2_rate" << s.sep << "emcee_sta_count" << s.sep
            << "emcee_sta_rate" << s.sep << "th_sta_count" << s.sep << "th_sta_rate" << s.sep
            << "cuevi_sta_count" << s.sep << "cuevi_sta_rate"

            << "\n";
    }
};

template <class ParameterType>
auto extract_parameters_last(std::string const& filename, std::size_t iter) {
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
void report_title(save_mcmc<Parameters, saving...>& f, Cuevi_mcmc<Parameters> const& data,
                  const Ts&... ts) {
    (report_title(get<saving>(f.m_m), data, ts...), ..., 1);
}

class step_stretch_cuevi_mcmc_per_walker {
   public:
    friend std::string ToString(step_stretch_cuevi_mcmc_per_walker) {
        return "step_stretch_cuevi_mcmc_per_walker";
    }

    template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
              class Parameters =
                  std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
        requires(is_prior<Prior, Parameters, Variables, DataType> &&
                 is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

    auto operator()(FunctionTable& f, Cuevi_mcmc<Parameters>& current, ensemble<mt_64i>& mts,
                    std::vector<std::uniform_real_distribution<double>>& rdist, Prior const& prior,
                    Likelihood const& lik, const by_fraction<DataType>& y,
                    const by_fraction<Variables>& x, std::size_t n_par, std::size_t i_th,
                    Walker_Index iw, Walker_Index jw, Cuevi_Index i_cu) const {
        auto& wa_i = current.get_Walker_Value(iw, i_cu);
        auto const& wa_j = current.get_Walker_Value(jw, i_cu);

        auto& param_i = get<Parameter>(wa_i())();
        auto& param_j = get<Parameter>(wa_j())();
        auto& t = current.get_Cuevi_Temperatures();

        auto [ca_par, z] = stretch_move(mts[i_th], rdist[i_th], param_i, param_j);
        auto wa_sta = current.get_Walker_Statistics(iw, i_cu);
        auto Maybe_Walker =
            current.calc_Walker_value(f, std::move(ca_par), prior, lik, y, x, t, i_cu, wa_sta);
        if (Maybe_Walker) {
            auto ca_wa = std::move(Maybe_Walker.value());
            auto& cu_wa = current.get_Walker_Value(iw, i_cu);
            auto dthLogL =
                thermo_step(ca_wa, cu_wa, current.get_Beta(i_cu), current.get_Fraction(i_cu));
            if (dthLogL) {
                auto pJump = std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL.value()));
                auto r = rdist[i_th](mts[i_th]);
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

    template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
              class Parameters =
                  std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
        requires(is_prior<Prior, Parameters, Variables, DataType> &&
                 is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

    void operator()(FunctionTable& f, Cuevi_mcmc<Parameters>& current, ensemble<mt_64i>& mts,
                    Prior const& prior, Likelihood const& lik, const by_fraction<DataType>& y,
                    const by_fraction<Variables>& x, double alpha_stretch = 2) const {
        auto cuevi_temp = current.get_Cuevi_Temperatures();
        auto n_walkers = current.get_Walkers_number();

        auto& mt = mts[0];
        auto n_par = current.get_Parameters_number();

        WalkerIndexes i_shuffled_walker(n_walkers);
        std::iota(i_shuffled_walker.begin(), i_shuffled_walker.end(), 0);
        std::shuffle(i_shuffled_walker.begin(), i_shuffled_walker.end(), mt);

        std::uniform_int_distribution<std::size_t> uniform_walker(0, n_walkers / 2 - 1);
        std::vector<std::uniform_int_distribution<std::size_t>> udist(n_walkers, uniform_walker);

        std::uniform_real_distribution<double> uniform_stretch_zdist(1.0 / alpha_stretch,
                                                                     alpha_stretch);
        std::vector<std::uniform_real_distribution<double>> zdist(n_walkers, uniform_stretch_zdist);

        std::uniform_real_distribution<double> uniform_real(0, 1);
        std::vector<std::uniform_real_distribution<double>> rdist(n_walkers, uniform_real);

        auto j_shuffled_walker = i_shuffled_walker;
        std::shuffle(j_shuffled_walker.begin(), j_shuffled_walker.begin() + n_walkers / 2, mt);

        auto ff = f.fork(omp_get_max_threads());
        for (bool half : {false, true}) {
#pragma omp parallel for  //collapse(2)
            for (std::size_t i = 0; i < n_walkers / 2; ++i) {
                for (std::size_t i_cu = 0; i_cu < size(cuevi_temp()); ++i_cu) {
                    auto i_th = omp_get_thread_num();

                    auto ii = half ? i + n_walkers / 2 : i;
                    auto jj = half ? i : i + n_walkers / 2;

                    auto iw = Walker_Index(i_shuffled_walker[ii]);
                    auto jw = Walker_Index(j_shuffled_walker[jj]);
                    ff[i_th].f(step_stretch_cuevi_mcmc_per_walker{}, current, mts, rdist, prior,
                               lik, y, x, n_par, i_th, iw, jw, Cuevi_Index(i_cu));
                }
            }
        }
        f += ff;
    }
};

class thermo_cuevi_jump_mcmc_per_walker {
   public:
    friend std::string ToString(thermo_cuevi_jump_mcmc_per_walker) {
        return "thermo_cuevi_jump_mcmc_per_walker";
    }

    template <class FunctionTable, class Likelihood, class Variables, class DataType,
              class Parameters>

    void operator()(FunctionTable& f, Cuevi_mcmc<Parameters>& current, Likelihood const& lik,
                    const by_fraction<DataType>& y, const by_fraction<Variables>& x, double r,
                    Walker_Index i_w, Walker_Index j_w, Cuevi_Index icu, Cuevi_Index jcu,
                    Random_jumps randomize) const {
        auto& t = current.get_Cuevi_Temperatures();
        auto& w_0 = current.get_Walker_Value(i_w, icu);
        auto& w_1 = current.get_Walker_Value(j_w, jcu);
        auto wa_sta_0 = current.get_Walker_Statistics(i_w, icu);
        auto wa_sta_1 = current.get_Walker_Statistics(j_w, jcu);
        auto Maybe_logA =
            current.thermo_jump_logProb(f, lik, y, x, t, icu, w_0, jcu, w_1, wa_sta_0, wa_sta_1);
        if (Maybe_logA.valid()) {
            auto logA = std::move(Maybe_logA.value());
            auto pJump = std::min(1.0, std::exp(logA));
            if (pJump > r) {
                auto& W0 = current.get_Walker(i_w, icu);
                auto& W1 = current.get_Walker(j_w, jcu);
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

    template <class FunctionTable, class Prior, class Likelihood, class Variables, class DataType,
              class Parameters =
                  std::decay_t<decltype(sample(std::declval<mt_64i&>(), std::declval<Prior&>()))>>
        requires(is_prior<Prior, Parameters, Variables, DataType> &&
                 is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables, DataType>)

    void operator()(FunctionTable& f, std::size_t iter, Cuevi_mcmc<Parameters>& current,
                    ensemble<mt_64i>& mts, Prior const&, Likelihood const& lik,
                    const by_fraction<DataType>& y, const by_fraction<Variables>& x,
                    Thermo_Jumps_every thermo_jumps_every, Random_jumps randomize) const {
        if ((iter > 0) && (iter % (thermo_jumps_every()) == 0)) {
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

            std::uniform_int_distribution<std::size_t> uniform_walker(0, n_walkers / 2 - 1);

            std::vector<std::uniform_real_distribution<double>> rdist(omp_get_max_threads(),
                                                                      uniform_real);

            Maybe_error<bool> success = true;

            auto ff = f.fork(omp_get_max_threads());

#pragma omp parallel for
            for (std::size_t i = 0; i < n_walkers / 2; ++i) {
                for (std::size_t i_cu = 0; i_cu + 1 < size(t()); ++i_cu) {
                    auto i_th = omp_get_thread_num();
                    auto i_w = Walker_Index(shuffled_walker[i]);
                    auto j_w = Walker_Index(shuffled_walker[i + n_walkers / 2]);
                    auto icu = Cuevi_Index(shuffled_cuevi[i_cu]);
                    auto jcu = Cuevi_Index(shuffled_cuevi[i_cu + 1]);
                    double r = rdist[i_th](mts[i_th]);
                    thermo_cuevi_jump_mcmc_per_walker{}(ff[i_th], current, lik, y, x, r, i_w, j_w,
                                                        icu, jcu, randomize);
                }
            }
            f += ff;
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
          var::Vector_Space<Num_Walkers_Per_Ensemble, var::Constant<Fractioner, myFractioner>,
                            var::Constant<Reporter, t_Reporter>,
                            var::Constant<Finalizer, t_Finalizer>, Fractions_Param, Th_Beta_Param,
                            Number_trials_until_give_up, Thermo_Jumps_every, Random_jumps,
                            Saving_intervals>> {
    using base_type =
        var::Var<Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>,
                 var::Vector_Space<
                     Num_Walkers_Per_Ensemble, var::Constant<Fractioner, myFractioner>,
                     var::Constant<Reporter, t_Reporter>, var::Constant<Finalizer, t_Finalizer>,
                     Fractions_Param, Th_Beta_Param, Number_trials_until_give_up,
                     Thermo_Jumps_every, Random_jumps, Saving_intervals>>;

   public:
    Cuevi_Algorithm(myFractioner&& frac, t_Reporter&& rep, t_Finalizer&& f,
                    Num_Walkers_Per_Ensemble n, Fractions_Param frac_param, Th_Beta_Param beta,
                    Number_trials_until_give_up max_iter_for_sampling,
                    Thermo_Jumps_every thermo_jumps_every, Random_jumps randomize,
                    Saving_intervals saving_intervals)
        : base_type(var::Vector_Space(
              std::move(n), var::Constant<Fractioner, myFractioner>(std::move(frac)),
              var::Constant<Reporter, t_Reporter>(std::move(rep)),
              var::Constant<Finalizer, t_Finalizer>(std::move(f)), std::move(frac_param),
              std::move(beta), std::move(max_iter_for_sampling), std::move(thermo_jumps_every),
              std::move(randomize), std::move(saving_intervals))) {
    }
};

template <class t_Reporter, class t_Finalizer>
class Cuevi_Algorithm_no_Fractioner
    : public var::Var<
          Cuevi_Algorithm_no_Fractioner<t_Reporter, t_Finalizer>,
          var::Vector_Space<Num_Walkers_Per_Ensemble, var::Constant<Reporter, t_Reporter>,
                            var::Constant<Finalizer, t_Finalizer>, Th_Beta_Param,
                            Number_trials_until_give_up, Thermo_Jumps_every, Random_jumps,
                            Saving_intervals>> {
    using base_type =
        var::Var<Cuevi_Algorithm_no_Fractioner<t_Reporter, t_Finalizer>,
                 var::Vector_Space<Num_Walkers_Per_Ensemble, var::Constant<Reporter, t_Reporter>,
                                   var::Constant<Finalizer, t_Finalizer>, Th_Beta_Param,
                                   Number_trials_until_give_up, Thermo_Jumps_every, Random_jumps,
                                   Saving_intervals>>;

   public:
    Cuevi_Algorithm_no_Fractioner(t_Reporter&& rep, t_Finalizer&& f, Num_Walkers_Per_Ensemble n,
                                  Th_Beta_Param beta,
                                  Number_trials_until_give_up max_iter_for_sampling,
                                  Thermo_Jumps_every thermo_jumps_every, Random_jumps randomize,
                                  Saving_intervals saving_intervals)
        : base_type(var::Vector_Space(
              std::move(n), var::Constant<Reporter, t_Reporter>(std::move(rep)),
              var::Constant<Finalizer, t_Finalizer>(std::move(f)), std::move(beta),
              std::move(max_iter_for_sampling), std::move(thermo_jumps_every), std::move(randomize),
              std::move(saving_intervals))) {
    }
};

template <class FunctionTable, class Duration, class Parameters, class... saving, class... T>
void report_all(FunctionTable& f, std::size_t iter, const Duration& dur,
                save_mcmc<Parameters, saving...>& s, Cuevi_mcmc<Parameters>& data, T const&... ts) {
    (report(f, iter, dur, get<saving>(s.m_m), data, ts...), ..., 1);
}

template <class FunctionTable, class ParameterType, class t_Reporter, class mcmc_type, class Prior,
          class Likelihood, class DataType, class Variables>
auto evidence_loop(FunctionTable& f, std::pair<mcmc_type, bool>&& mcmc_run, t_Reporter& rep,
                   std::size_t iter, Cuevi_mcmc<ParameterType>& current,

                   std::vector<mt_64i>& mts, Random_jumps randomjumps,
                   Thermo_Jumps_every v_thermo_jump_every, Prior&& prior, Likelihood const& lik,
                   const by_fraction<DataType>& ys, const by_fraction<Variables>& xs) {
    using Return_Type = Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<ParameterType>>>;
    report_title(rep, current, mcmc_run.first, prior, lik, ys, xs);
    const auto start = std::chrono::high_resolution_clock::now();

    while (!mcmc_run.second) {
        f.f(step_stretch_cuevi_mcmc{}, current, mts, prior, lik, ys, xs);
        report_point(f, iter);

        ++iter;
        f.f(thermo_cuevi_jump_mcmc{}, iter, current, mts, prior, lik, ys, xs, v_thermo_jump_every,
            Random_jumps(false));
        if (randomjumps())
            f.f(thermo_cuevi_jump_mcmc{}, iter + v_thermo_jump_every() % 2, current, mts, prior,
                lik, ys, xs, v_thermo_jump_every, Random_jumps(true));

        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start);
        report_all(f, iter, dur, rep, current, prior, lik, ys, xs, mts, mcmc_run.first);
        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
    }
    return Return_Type(std::pair(mcmc_run.first, current));
}

template <class FunctionTable, class myFractioner, class t_Reporter, class t_Finalizer, class Prior,
          class Likelihood, class DataType, class Variables>
auto evidence(FunctionTable& f, Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>&& cue,
              Prior&& prior, Likelihood const& lik, const DataType& y, const Variables& x,
              const Init_seed init_seed) {
    using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i&>>;

    using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

    using Return_Type = Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;

    auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
    auto mt = init_mt(init_seed());
    auto mts = init_mts(mt, omp_get_max_threads());
    auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
    auto n_points_per_decade_fraction = get<Points_per_decade>(get<Fractions_Param>(cue())());

    auto [ys, xs] = get<Fractioner>(cue())()(y, x, mt, size(prior) * min_fraction(),
                                             n_points_per_decade_fraction());

    auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
        f, mts, std::forward<Prior>(prior), lik, ys, xs, get<Th_Beta_Param>(cue()), n_walkers,
        get<Number_trials_until_give_up>(cue()));

    if (!Maybe_current)
        return Return_Type(Maybe_current.error());
    else {
        auto current = std::move(Maybe_current.value());
        auto& rep = get<Reporter>(cue())();

        auto a = get<Finalizer>(cue())();
        auto mcmc_run = checks_convergence(std::move(a), current);

        std::size_t iter = 0;
        report_title(f, "Iter");
        auto v_random_jumps = get<Random_jumps>(cue());
        auto v_thermo_jump_every = get<Thermo_Jumps_every>(cue());
        report_model_all(rep, v_thermo_jump_every, prior, lik, ys, xs);

        // report_init(rep,v_thermo_jump_every,prior,lik,ys,xs);

        return evidence_loop(f, std::move(mcmc_run), rep, iter, current, mts, v_random_jumps,
                             v_thermo_jump_every, std::forward<Prior>(prior), lik, ys, xs);
    }
}

template <class FunctionTable, class t_Reporter, class t_Finalizer, class Prior, class Likelihood,
          class DataType, class Variables>
auto evidence_fraction(FunctionTable& f,
                       Cuevi_Algorithm_no_Fractioner<t_Reporter, t_Finalizer>&& cue, Prior&& prior,
                       Likelihood const& lik, const std::vector<DataType>& ys,
                       const std::vector<Variables>& xs, const Init_seed init_seed) {
    using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i&>>;

    using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

    using Return_Type = Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;

    auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
    auto mt = init_mt(init_seed());
    auto mts = init_mts(mt, omp_get_max_threads());
    //   auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
    // auto n_points_per_decade_fraction =
    //     get<Points_per_decade>(get<Fractions_Param>(cue())());

    // auto [ys, xs] = get<Fractioner>(cue())()(
    //     y, x, mt, size(prior) * min_fraction(),
    //     n_points_per_decade_fraction());

    auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
        f, mts, std::forward<Prior>(prior), lik, ys, xs, get<Th_Beta_Param>(cue()), n_walkers,
        get<Number_trials_until_give_up>(cue()));

    if (!Maybe_current)
        return Return_Type(Maybe_current.error());
    else {
        auto current = std::move(Maybe_current.value());
        auto& rep = get<Reporter>(cue())();

        auto a = get<Finalizer>(cue())();
        auto mcmc_run = checks_convergence(std::move(a), current);

        std::size_t iter = 0;
        report_title(f, "Iter");
        auto v_random_jumps = get<Random_jumps>(cue());
        auto v_thermo_jump_every = get<Thermo_Jumps_every>(cue());
        report_model_all(rep, v_thermo_jump_every, prior, lik, ys, xs);

        // report_init(rep,v_thermo_jump_every,prior,lik,ys,xs);

        return evidence_loop(f, std::move(mcmc_run), rep, iter, current, mts, v_random_jumps,
                             v_thermo_jump_every, std::forward<Prior>(prior), lik, ys, xs);
    }
}

template <class FunctionTable, class myFractioner, class t_Reporter, class t_Finalizer, class Prior,
          class Likelihood, class DataType, class Variables>
auto continue_evidence(FunctionTable& f,
                       Cuevi_Algorithm<myFractioner, t_Reporter, t_Finalizer>&& cue, Prior&& prior,
                       Likelihood const& lik, const DataType& y, const Variables& x,
                       const Init_seed init_seed) {
    using Parameter_Type = std::decay_t<std::invoke_result_t<Prior, mt_64i&>>;

    using mcmc_type = std::decay_t<decltype(get<Finalizer>(cue())())>;

    using Return_Type = Maybe_error<std::pair<mcmc_type, Cuevi_mcmc<Parameter_Type>>>;

    auto n_walkers = get<Num_Walkers_Per_Ensemble>(cue());
    auto mt = init_mt(init_seed());
    auto mts = init_mts(mt, omp_get_max_threads());
    auto min_fraction = get<Min_value>(get<Fractions_Param>(cue())());
    auto n_points_per_decade_fraction = get<Points_per_decade>(get<Fractions_Param>(cue())());

    auto [ys, xs] = get<Fractioner>(cue())()(y, x, mt, size(prior) * min_fraction(),
                                             n_points_per_decade_fraction());

    auto Maybe_current = Cuevi_mcmc<Parameter_Type>::init(
        f, mts, std::forward<Prior>(prior), lik, ys, xs, get<Th_Beta_Param>(cue()), n_walkers,
        get<Number_trials_until_give_up>(cue()));

    if (!Maybe_current)
        return Return_Type(Maybe_current.error());
    else {
        auto current = std::move(Maybe_current.value());
        auto& rep = get<Reporter>(cue())();
        report_title(rep, current, prior, lik, ys, xs);

        auto a = get<Finalizer>(cue())();
        auto mcmc_run = checks_convergence(std::move(a), current);

        std::size_t iter = 0;
        // report_model(rep, prior, lik, ys, xs);
        report_title(f, "Iter");
        std::size_t v_thermo_jump_every = get<Thermo_Jumps_every>(cue());

        return evidence_loop(f, v_thermo_jump_every, mcmc_run, rep, iter, current, mts,
                             std::forward<Prior>(prior), lik, ys, xs);
    }
}

}  // namespace cuevi

#endif  // CUEVI_H
