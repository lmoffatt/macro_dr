#ifndef QMODEL_H
#define QMODEL_H
//#include "cuevi.h"
// #include "experiment.h"
#include <derivative_fwd.h>
#include <distributions.h>

#include "experiment.h"
#include "fold.h"
#include "function_memoization.h"
#include "matrix.h"
#include "parallel_levenberg_tempering.h"
#include "parallel_tempering_fraction.h"
#include "parameters_derivative.h"
// Ensure LAPACK-backed helpers are available to Macro_DMR and friends
#include "lapack_headers.h"
// Schur+Parlett path for Qdt — stable at clustered eigenvalues where the eig
// path's eigenvector inverse becomes ill-conditioned (lifted micro Q at
// Nch ≥ ~50, k=2). Provides expm_schur_parlett, frechet_integral_schur,
// frechet_double_integral_schur as the building blocks of calc_Qdt_schur.
#include "schur_parlett.h"
// #include "models_MoffattHume_linear.h"
#include <macrodr/dsl/type_name.h>
#include <moment_statistics.h>

#include <atomic>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "derivative_operator.h"
#include "derivative_test.h"
#include "exponential_matrix.h"
#include "function_measure_verification_and_optimization.h"
#include "general_algorithm_on_containers.h"
#include "general_output_operator.h"
#include "gsl_integrate.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "type_algebra.h"
#include "variables.h"
#include "qmodel_types.h"
namespace macrodr {
/*
class State_Model;


Maybe_error<State_Model>
to_State_Model(std::size_t t_number_of_states,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_transition_rates,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_agonist_transition_rates,
std::map<std::size_t, std::string> t_conductances);

class State_Model {
std::size_t m_number_of_states;

std::map<std::pair<std::size_t, std::size_t>, std::string> m_transition_rates;
std::map<std::pair<std::size_t, std::size_t>, std::string>
m_agonist_transition_rates;
std::map<std::size_t, std::string> m_conductances;

State_Model(std::size_t t_number_of_states,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_transition_rates,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_agonist_transition_rates,
std::map<std::size_t, std::string> t_conductances)
: m_number_of_states{t_number_of_states},
m_transition_rates{t_transition_rates},
m_agonist_transition_rates{t_agonist_transition_rates},
m_conductances{t_conductances} {}

public:
friend Maybe_error<State_Model>
to_State_Model(std::size_t t_number_of_states,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_transition_rates,
std::map<std::pair<std::size_t, std::size_t>, std::string>
&t_agonist_transition_rates,
std::map<std::size_t, std::string> t_conductances) {
for (auto &elem : t_transition_rates) {
if (elem.first.first >= t_number_of_states)
return error_message(
"transition start state greater than number of "
"states; number_of_states = " +
std::to_string(t_number_of_states) +
" start state= " + std::to_string(elem.first.first));
else if (elem.first.second >= t_number_of_states)
return error_message(
"transition end state greater than number of "
"states; number_of_states = " +
std::to_string(t_number_of_states) +
" end state= " + std::to_string(elem.first.second));
else if (elem.first.second == t_number_of_states)
return error_message(
"transition start state same as end state;  start state: " +
std::to_string(elem.first.first) +
" end state= " + std::to_string(elem.first.second));
}
//
for (auto &elem : t_agonist_transition_rates) {
if (elem.first.first >= t_number_of_states)
return error_message(
"agonist transition start state greater than number of "
"states; number_of_states = " +
std::to_string(t_number_of_states) +
" start state= " + std::to_string(elem.first.first));
else if (elem.first.second >= t_number_of_states)
return error_message(
"agonist transition end state greater than number of "
"states; number_of_states = " +
std::to_string(t_number_of_states) +
" end state= " + std::to_string(elem.first.second));
else if (elem.first.second == t_number_of_states)
return error_message(
"agonist transition start state same as end state;  start state: " +
std::to_string(elem.first.first) +
" end state= " + std::to_string(elem.first.second));
}

std::size_t i_cond = 0;
for (auto &elem : t_conductances) {
if (elem.first >= t_number_of_states)
return error_message(
"state conductance number greater than number of states:"
" proposed= " +
std::to_string(elem.first) +
" number_of_states = " + std::to_string(t_number_of_states));

if (elem.first != i_cond)
return error_message("state conductance skipped: current is" +
     std::to_string(i_cond) +
     " proposed= " + std::to_string(elem.first));
++i_cond;
}
if (i_cond != t_number_of_states)
return error_message(
"state conductance missing: number of proposed states=" +
std::to_string(i_cond) +
" number_of_states = " + std::to_string(t_number_of_states));

return State_Model(t_number_of_states, t_transition_rates,
t_agonist_transition_rates, t_conductances);
}
};
*/


/*
class Transition_rate_resting
: public var::Var<Transition_rate_resting, Matrix<double>,
            Power<var::s, -1>> {};
            
class Transition_rate_agonist
: public var::Var<Transition_rate_agonist, Matrix<double>,
            Power<var::Product<var::s, var::microMolar>, -1>> {};
            
class Transition_rate
: public var::Var<Transition_rate, Matrix<double>, Power<var::s, -1>> {};

class State_unitary_current
: public var::Var<State_unitary_current, Matrix<double>,
            Product<var::pA, Power<var::number, -1>>> {};
            
class Number_of_States
: public var::Var<Number_of_States, std::size_t, var::number> {};

class Number_of_Channels_mean
: public var::Var<Number_of_Channels_mean, double, var::number> {};

class Number_of_Channels_stddev
: public var::Var<Number_of_Channels_stddev, double, var::number> {};

class Current_noise : public var::Var<Current_noise, double, var::pA> {};

class State_Probability
: public var::Var<State_Probability, Matrix<double>, var::prob> {};

class State_Probability_Cov
: public var::Var<State_Probability_Cov, Matrix<double>,
            Power<var::prob, 2>> {};
            
class Transition_rate_landa
: public var::Var<Transition_rate_landa, DiagonalMatrix<double>,
            Power<var::s, -1>> {};
            
class Transition_rate_V
: public var::Var<Transition_rate_V, Matrix<double>, var::dimensionless> {};
class Transition_rate_W
: public var::Var<Transition_rate_W, Matrix<double>, var::dimensionless> {};
*/


// class plogL : public var::Var<plogL, double> {
//     friend std::string className(plogL) { return "plogL"; }
// };
// class eplogL : public var::Var<eplogL, double> {
//     friend std::string className(eplogL) { return "eplogL"; }
// };
// class vplogL : public var::Constant<vplogL, double> {
//     friend std::string className(vplogL) { return "vplogL"; }
// };


/*

template <class T>
inline auto get_mean_Probits(std::vector<Evolution_of<T>>const & bootstrap_estimates,
                             const std::set<double>& cis) {
    assert(!bootstrap_estimates.empty());

    std::vector<std::vector<T>> values;
    values.reserve(bootstrap_estimates.size());
    for (auto& estimate : bootstrap_estimates)
        values.push_back(estimate());

    auto [mean_value, probits] = get_mean_Probits(values, cis);

    Evolution_of<T> mean_out(std::move(mean_value));
    std::map<double, Evolution_of<T>> probits_out;
    for (auto& [level, value] : probits)
        probits_out.emplace(level, Evolution_of<T>(std::move(value)));

    return std::make_pair(std::move(mean_out), std::move(probits_out));
}
*/
// has_var_c moved to qmodel_types.h.


template <class C_Patch_Model, class C_double>
C_Patch_Model add_Patch_inactivation(C_Patch_Model&& m, C_double const& deactivation_rate) {
    using Transf = transformation_type_t<C_Patch_Model>;
    auto Nst = get<N_St>(m)() + 1;
    Op_t<Transf, Q0> v_Q0 = Q0(Matrix<double>(Nst, Nst, 0.0));
    for (std::size_t i = 0; i + 1 < Nst; ++i) {
        for (std::size_t j = 0; j + 1 < Nst; ++j) set(v_Q0(), i, j, get<Q0>(m)()(i, j));
        set(v_Q0(), i, Nst - 1, deactivation_rate);
    }
    Op_t<Transf, Qa> v_Qa = Qa(Matrix<double>(Nst, Nst, 0.0));
    for (std::size_t i = 0; i + 1 < Nst; ++i) {
        for (std::size_t j = 0; j + 1 < Nst; ++j) set(v_Qa(), i, j, get<Qa>(m)()(i, j));
    }
    Op_t<Transf, g> v_g = g(Matrix<double>(Nst, 1, 0.0));
    for (std::size_t i = 0; i + 1 < Nst; ++i) {
        set(v_g(), i, 0, get<g>(m)()[i]);
    }
    Op_t<Transf, P_initial> v_Pini = P_initial(Matrix<double>(1, Nst, 0.0));
    for (std::size_t i = 0; i + 1 < Nst; ++i) {
        set(v_Pini(), 0, i, get<P_initial>(m)()[i]);
    }

    get<N_St>(m)() = Nst;
    get<Qa>(m) = v_Qa;
    get<Q0>(m) = v_Q0;
    get<g>(m) = v_g;
    get<P_initial>(m) = v_Pini;
    return std::forward<C_Patch_Model>(m);
}

template <class C_Patch_Model>
std::pair<C_Patch_Model, double> remove_Patch_inactivation(C_Patch_Model const& mi) {
    auto m = mi;

    using Transf = transformation_type_t<C_Patch_Model>;
    auto Nst = get<N_St>(m)() - 1;
    Op_t<Transf, Q0> v_Q0 = Q0(Matrix<double>(Nst, Nst, 0.0));
    auto deactivation_rate = get<Q0>(m)()(0ul, Nst);
    for (std::size_t i = 0; i < Nst; ++i) {
        for (std::size_t j = 0; j < Nst; ++j) set(v_Q0(), i, j, get<Q0>(m)()(i, j));
        assert(deactivation_rate == get<Q0>(m)()(i, Nst));
    }

    Op_t<Transf, Qa> v_Qa = Qa(Matrix<double>(Nst, Nst, 0.0));
    for (std::size_t i = 0; i < Nst; ++i) {
        for (std::size_t j = 0; j < Nst; ++j) set(v_Qa(), i, j, get<Qa>(m)()(i, j));
    }
    Op_t<Transf, g> v_g = g(Matrix<double>(Nst, 1, 0.0));
    for (std::size_t i = 0; i < Nst; ++i) {
        set(v_g(), i, 0, get<g>(m)()[i]);
    }
    Op_t<Transf, P_initial> v_Pini = P_initial(Matrix<double>(1, Nst, 0.0));
    for (std::size_t i = 0; i < Nst; ++i) {
        set(v_Pini(), 0, i, get<P_initial>(m)()[i]);
    }
    v_Pini() = v_Pini() / var::sum(v_Pini());

    get<N_St>(m)() = Nst;
    get<Qa>(m) = v_Qa;
    get<Q0>(m) = v_Q0;
    get<g>(m) = v_g;
    get<P_initial>(m) = v_Pini;
    return std::pair(std::move(m), deactivation_rate);
}


template <typename Simulate_tag>
void save_simulation(std::string const& filename, const std::string& separator,
                     Simulated_Recording<Simulate_tag> const& sim) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    f << "i_step" << separator << "patch_current";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << separator << "i_state" << separator << "N_state";
    f << "\n";

    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        auto N = get<N_Ch_State_Evolution>(sim());
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
            for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
                f << i_step << separator << y()[i_step]() << separator << i_state << separator
                  << N()[i_step]()[i_state] << "\n";
            }
        }
    } else {
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
            f << i_step << separator << y()[i_step]() << "\n";
        }
    }
}

template <typename Simulate_tag>
Maybe_error<bool> load_simulation(std::string const& fname, std::string separator,
                                  Simulated_Recording<Simulate_tag>& r) {
    std::ifstream f(fname);
    if (!f)
        return error_message("cannot open file " + fname);
    std::string line;
    std::getline(f, line);
    std::stringstream ss(line);

    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        if (!(ss >> septr("i_step") >> septr(separator) >> septr("patch_current") >>
              septr(separator) >> septr("i_state") >> septr(separator) >> septr("N_state")))
            return error_message("titles are wrong : expected  \ni_step" + separator +
                                 "patch_current" + separator + "i_state" + separator +
                                 "N_state\n found:\n" + line);
    } else {
        if (!(ss))
            return error_message("titles are wrong : expected  i_state" + separator +
                                 "N_state\n found:\n" + line);
        if (!(ss >> septr("i_step") >> septr(separator) >> septr("patch_current")))
            return error_message("titles are wrong : expected  i_step:" + separator +
                                 "patch_current; found:" + line);
    }
    std::getline(f, line);
    ss = std::stringstream(line);
    std::size_t i_step;
    std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();

    double val;
    auto& e = get<Recording>(r());

    if constexpr (!var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        while (extract_double(ss >> i_step >> septr(separator), val, separator[0])) {
            if (i_step_prev != i_step) {
                if (i_step != e().size())
                    return error_message("i_step missmatch expected" + std::to_string(e().size()) +
                                         " found:" + std::to_string(i_step));
                e().push_back(Patch_current(val));
                i_step_prev = i_step;
            }
            std::getline(f, line);
            ss = std::stringstream(line);
        }
        return true;
    } else {
        std::size_t i_state;
        double N;
        std::vector<double> N_state;
        auto n_channel_states = N_state.size();
        while ((extract_double(ss >> i_step >> septr(separator), val, separator[0]) >>
                septr(separator) >> i_state >> septr(separator) >> N)) {
            if (i_step_prev >= i_step) {
                if (i_state != N_state.size())
                    return error_message(
                        "i_state missmatch expected =" + std::to_string(N_state.size()) +
                        " found " + std::to_string(i_state));
                N_state.push_back(N);
                if (i_step_prev > i_step)
                    e().push_back(Patch_current(val));

            } else {
                if (i_step != e().size())
                    return error_message("i_step missmatch expected" + std::to_string(e().size()) +
                                         " found:" + std::to_string(i_step));
                e().push_back(Patch_current(val));
                if (n_channel_states == 0)
                    n_channel_states = N_state.size();
                else if (n_channel_states != N_state.size())
                    return error_message("n_channel_states missmatch expected" +
                                         std::to_string(n_channel_states) +
                                         " found:" + std::to_string(N_state.size()));
                auto& Ns = get<N_Ch_State_Evolution>(r());
                Ns().emplace_back(N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
                N_state.clear();
                N_state.push_back(N);
            }
            i_step_prev = i_step;

            std::getline(f, line);
            ss = std::stringstream(line);
        }
        if (n_channel_states != N_state.size())
            return error_message("n_channel_states missmatch expected" +
                                 std::to_string(n_channel_states) +
                                 " found:" + std::to_string(N_state.size()));
        auto& Ns = get<N_Ch_State_Evolution>(r());
        Ns().emplace_back(N_channel_state(Matrix<double>(1, n_channel_states, N_state)));

        return true;
    }
}

template <typename Simulate_tag>
void save_fractioned_simulation(std::string filename, std::string separator,
                                std::vector<Simulated_Recording<Simulate_tag>> const& vsim) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    f << "i_frac" << separator << "i_step" << separator << "patch_current";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << separator << "i_state" << separator << "N_state";
    f << "\n";

    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        for (std::size_t i_frac = 0; i_frac < vsim.size(); ++i_frac) {
            auto& sim = vsim[i_frac];
            auto N = get<N_Ch_State_Evolution>(sim());
            auto y = get<Recording>(sim());
            for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
                for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
                    f << i_frac << separator << i_step << separator << y()[i_step]() << separator
                      << i_state << separator << N()[i_step]()[i_state] << "\n";
                }
            }
        }
    } else {
        for (std::size_t i_frac = 0; i_frac < vsim.size(); ++i_frac) {
            auto& sim = vsim[i_frac];
            auto y = get<Recording>(sim());
            for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
                f << i_frac << separator << i_step << separator << y()[i_step]() << "\n";
            }
        }
    }
}

template <typename Simulate_tag>
Maybe_error<bool> load_fractioned_simulation(std::string const& fname, std::string separator,
                                             std::vector<Simulated_Recording<Simulate_tag>>& rs) {
    std::ifstream f(fname);
    if (!f)
        return error_message("cannot open file " + fname);
    std::string line;
    std::getline(f, line);
    std::stringstream ss(line);

    if (!(ss >> septr("i_frac") >> septr(separator) >> septr("i_step") >> septr(separator) >>
          septr("patch_current")))
        return error_message("titles are wrong : expected  i_step:" + separator +
                             "patch_current; found:" + line);
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        if (!(ss >> septr(separator) >> septr("i_state") >> septr(separator) >> septr("N_state")))
            return error_message("titles are wrong : expected  i_state:" + separator +
                                 "N_state; found:" + line);

    std::getline(f, line);
    ss = std::stringstream(line);
    std::size_t i_frac;
    std::size_t i_frac_prev = std::numeric_limits<std::size_t>::max();
    std::size_t i_step;
    std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();

    double val;
    Simulated_Recording<Simulate_tag> r;

    if constexpr (!var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        while (extract_double(ss >> i_frac >> septr(separator) >> i_step >> septr(separator), val,
                              separator[0])) {
            if (i_frac_prev != i_frac) {
                if (i_frac != rs.size())
                    return error_message("i_step missmatch expected " + std::to_string(rs.size()) +
                                         " found:" + std::to_string(i_step));
                if (i_frac > 0) {
                    rs.push_back(r);
                    r = Simulated_Recording<Simulate_tag>{};
                }
                i_step_prev = std::numeric_limits<std::size_t>::max();
            }
            auto& e = get<Recording>(r());
            if (i_step_prev != i_step) {
                if (i_step != e().size())
                    return error_message("i_step missmatch expected " + std::to_string(e().size()) +
                                         " found:" + std::to_string(i_step));
                e().push_back(Patch_current(val));
                i_step_prev = i_step;
            }
            std::getline(f, line);
            ss = std::stringstream(line);
        }
        return true;
    } else {
        double val_prev = std::numeric_limits<double>::max();
        std::size_t i_state;
        double N;
        std::vector<double> N_state;
        auto n_channel_states = N_state.size();
        while ((extract_double(ss >> i_frac >> septr(separator) >> i_step >> septr(separator), val,
                               separator[0]) >>
                septr(separator) >> i_state >> septr(separator) >> N)) {
            if (i_frac_prev != i_frac) {
                if ((i_frac > 0) && (i_frac != rs.size() + 1))
                    return error_message("i_frac missmatch expected: " + std::to_string(rs.size()) +
                                         " found: " + std::to_string(i_frac));
                if (i_frac > 0) {
                    auto& Ns = get<N_Ch_State_Evolution>(r());
                    Ns().emplace_back(
                        N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
                    rs.push_back(r);
                    N_state.clear();
                    r = Simulated_Recording<Simulate_tag>{};
                }
                i_frac_prev = i_frac;
            }

            auto& e = get<Recording>(r());
            if (i_step_prev == i_step) {
                if ((val != val_prev) && !(std::isnan(val) && std::isnan(val_prev)))
                    return error_message(
                        "change patch current in same i_step val=" + std::to_string(val) +
                        " prev_val=" + std::to_string(val_prev));
                N_state.push_back(N);

            } else {
                if (i_step != e().size())
                    return error_message(
                        "i_step missmatch expected: " + std::to_string(e().size()) +
                        " found: " + std::to_string(i_step));
                e().push_back(Patch_current(val));
                val_prev = val;
                if (i_step > i_step_prev) {
                    if (n_channel_states == 0)
                        n_channel_states = N_state.size();
                    else if (n_channel_states != N_state.size())
                        return error_message("n_channel_states missmatch expected: " +
                                             std::to_string(n_channel_states) +
                                             " found: " + std::to_string(N_state.size()));
                    auto& Ns = get<N_Ch_State_Evolution>(r());
                    Ns().emplace_back(
                        N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
                    N_state.clear();
                }
                N_state.push_back(N);
                i_step_prev = i_step;
            }
            std::getline(f, line);
            ss = std::stringstream(line);
        }
        if (n_channel_states != N_state.size())
            return error_message("n_channel_states missmatch expected" +
                                 std::to_string(n_channel_states) +
                                 " found:" + std::to_string(N_state.size()));
        auto& Ns = get<N_Ch_State_Evolution>(r());
        Ns().emplace_back(N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
        rs.push_back(r);

        return true;
    }
}

inline Maybe_error<bool> load_simulation(std::string const& fname, std::string separator,
                                         v_Simulated_Recording& r) {
    Simulated_Recording<var::please_include<N_Ch_State_Evolution>> try_N_state;
    auto Maybe_N_state = load_simulation(fname, separator, try_N_state);
    if (Maybe_N_state.valid()) {
        r = try_N_state;
        return true;
    }
    r = Simulated_Recording<var::please_include<>>{};
    return load_simulation(fname, separator, r);
}

template <typename Simulate_tag>
void save_simulation(std::string name, std::vector<Simulated_Recording<Simulate_tag>> const& r) {
    std::ofstream f(name);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    f << "nrep"
      << ","
      << "i_step"
      << ","
      << "patch_current";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        f << ","
          << "i_state"
          << ","
          << "i_state"
          << "N";
    }
    f << "\n";

    for (std::size_t i = 0; i < r.size(); ++i)
        for (std::size_t n = 0; n < r[i]()().size(); ++n) {
            f << i << "," << n << "," << r[i]()()[n]() << "\n";
        }
}

template <typename uses_recursive, typename uses_averaging, typename uses_variance,
          typename uses_variance_form = uses_variance_form_aproximation<variance_total>>
    requires(var::is_this_constexpr_Var_v<uses_recursive, bool, uses_recursive_aproximation> &&
             var::is_this_constexpr_Var_v<uses_averaging, int, uses_averaging_aproximation> &&
             var::is_this_constexpr_Var_v<uses_variance, bool, uses_variance_aproximation> &&
             uses_variance_form_aproximation_c<uses_variance_form>)
struct MacroR {
    friend std::string ToString(MacroR) {
        std::string out = "MacroR";

        out += uses_recursive::value ? "_R" : "_NR";
        out += (uses_averaging::value == 2) ? "_2" : "__";
        out += uses_variance::value ? "_V" : "_M";
        if constexpr (uses_variance_form::value == variance_residual)
            out += "_res";

        return out;
    }
};

struct Calc_Qdt {
    friend std::string ToString(Calc_Qdt) { return "Calc_Qdt"; }
};

struct Calc_Qdt_step {
    friend std::string ToString(Calc_Qdt_step) { return "Calc_Qdt_step"; }
};

struct Calc_Qdtm_step {
    friend std::string ToString(Calc_Qdtm_step) { return "Calc_Qdtm_step"; }
};

struct Calc_Qdtg_step {
    friend std::string ToString(Calc_Qdtg_step) { return "Calc_Qdtg_step"; }
};

struct Calc_Qx {
    friend std::string ToString(Calc_Qx) { return "Calc_Qx"; }
};

struct Calc_eigen {
    friend std::string ToString(Calc_eigen) { return "Calc_eigen"; }
};

template <class Policy = StabilizerPolicyEnabled, typename CQx>
    requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>> full_expm(const CQx& x) {
    assert(x.ncols() == x.nrows());
    assert(x.size() > 0);

    // Scale A by power of 2 so that its norm is < 1/2 .
    std::size_t s = log2_norm(primitive(x())) + 1;

    auto A = x * (1.0 / std::pow(2.0, int(s)));

    // Pade approximation for exp(A)
    auto eE = expm_pade(A);

    if (!eE)
        return eE.error();
    else {
        auto Maybe_E = to_Transition_Probability<Policy>(eE.value());
        if (!Maybe_E)
            return Maybe_E.error();
        else {
            auto E = std::move(Maybe_E.value());
            // Undo scaling by repeated squaring
            for (std::size_t k = 0; k < s; k++) {
                Maybe_E = to_Transition_Probability<Policy>(E() * E());
                if (!Maybe_E)
                    return Maybe_E.error();
                else
                    E = std::move(Maybe_E.value());
            }
            return E;
        }
    }
}

template <class Policy = StabilizerPolicyEnabled, typename CQx>
    requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>> expm_taylor_scaling_squaring(const CQx& x,
                                                                 std::size_t order = 6) {
    {
        double max = maxAbs(primitive(x()));
        double desired = 0.125 / 8.0;
        int k = std::ceil(std::log2(max / desired));
        int n = std::max(0, k);
        double scale = std::pow(2, -n);
        auto dx = x() * scale;
        auto expm_dx = expm_taylor(dx, order);
        auto Maybe_expm_run = to_Transition_Probability<Policy>(expm_dx);
        if (!Maybe_expm_run)
            return Maybe_expm_run.error();
        else {
            auto expm_run = std::move(Maybe_expm_run.value());
            for (std::size_t i = 0; i < n; ++i) {
                Maybe_expm_run = to_Transition_Probability<Policy>(expm_run() * expm_run());
                if (!Maybe_expm_run)
                    return Maybe_expm_run.error();
                else
                    expm_run = std::move(Maybe_expm_run.value());
            }
            return expm_run;
        }
    }
}

template <class Policy = StabilizerPolicyEnabled, typename CQx>
    requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>> expm_sure(const CQx& x) {
    auto Maybe_expm = full_expm<Policy>(x);
    if (Maybe_expm)
        return Maybe_expm.value();
    else
        return expm_taylor_scaling_squaring<Policy>(x);
}

template <class recursive, class averaging, class variance, class variance_correction>
    requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
             uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
struct MacroR2;
class Macro_DMR {
    template <class C_double>
        requires U<C_double, double>
    static C_double E1(C_double const& x, C_double const& exp_x,
                       double eps = std::numeric_limits<double>::epsilon() * 100) {
        if (abs(primitive(x) * primitive(x)) < eps)
            return 1.0;
        else
            return (exp_x - 1.0) / x;
    }

    static double E1(double x) {
        if (std::abs(x) < std::numeric_limits<double>::epsilon() * 100)
            return 1.0;
        else if (std::abs(x) < 1e-2)
            return std::expm1(x) / x;
        else
            return (std::exp(x) - 1.0) / x;
    }

    static double E2(double x, double y) {
        const double eps = std::numeric_limits<double>::epsilon();
        if (x * x < eps) {
            if (y * y < eps)
                return 0.5;
            else
                return (E1(y) - 1.0) / y;
        } else if (y * y < eps)
            return (E1(x) - 1.0) / x;
        else if ((y - x) * (y - x) < eps)
            return (std::exp(x) - E1(x)) / x;
        else
            return (E1(y) - E1(x)) / (y - x);
    }

    template <class C_double>
        requires U<C_double, double>
    static C_double E2(C_double const& x, C_double const& y, C_double const& exp_x,
                       C_double const& exp_y, double eps = std::numeric_limits<double>::epsilon()) {
        if (primitive(x) * primitive(x) < eps) {
            if (primitive(y) * primitive(y) < eps)
                return 0.5;
            else
                return (E1(y, exp_y, eps) - 1.0) / y;
        } else if (primitive(y) * primitive(y) < eps)
            return (E1(x, exp_x) - 1.0) / x;
        else if ((primitive(y) - primitive(x)) * (primitive(y) - primitive(x)) < eps)
            return (exp_x - E1(x, exp_x)) / x;
        else
            return (E1(y, exp_y) - E1(x, exp_x)) / (y - x);
    }

    template <class C_double>
        requires U<C_double, double>
    static C_double Ee(C_double const& x, C_double const& y, C_double const& exp_x,
                       C_double const& exp_y, double eps = std::numeric_limits<double>::epsilon()) {
        if (sqr(primitive(x) - primitive(y)) < eps)
            return exp_x;
        else
            return (exp_x - exp_y) / (x - y);
    };

    template <class C_double>
        requires U<C_double, double>
    static C_double EX_111(C_double const& x, C_double const& y, C_double const& z,
                           C_double const& exp_x) {
        return exp_x / ((x - y) * (x - z));
    }

    template <class C_double>
        requires U<C_double, double>
    static C_double E111(C_double const& x, C_double const& y, C_double const& z,
                         C_double const& exp_x, C_double const& exp_y, C_double const& exp_z) {
        return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) + EX_111(z, y, x, exp_z);
    }
    template <class C_double>
        requires U<C_double, double>
    static C_double E12(C_double const& x, C_double const& y, C_double const& exp_x,
                        C_double const& exp_y) {
        return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
    }

    template <class C_double>
        requires U<C_double, double>
    static C_double E3(C_double const& x, C_double const& y, C_double const& z,
                       C_double const& exp_x, C_double const& exp_y, C_double const& exp_z,
                       double eps = std::numeric_limits<double>::epsilon()) {
        auto x_ = primitive(x);
        auto y_ = primitive(y);
        auto z_ = primitive(z);

        if (sqr(x_ - y_) < eps)  // x==y
        {
            if (sqr(y_ - z_) < eps)  // y==z
                return exp_x / 2.0;  // x==y==z
            else
                return E12(z, x, exp_z, exp_x);  // x==y!=z
        } else if (sqr(y_ - z_) < eps)           // x!=y==z
        {
            return E12(x, y, exp_x, exp_y);
        } else if (sqr(x_ - z_) < eps)  // y!=z==x!=y
        {
            return E12(y, x, exp_y, exp_x);
        } else
            return E111(x, y, z, exp_x, exp_y, exp_z);  // x!=y!=z!=x
    }

    // P_Cov is stored as the *bare centered covariance* of the one-hot
    // channel-state indicator: Cov(X) = diag(P_mean) − P_mean·P_meanᵀ.
    // Sanity check: deterministic state P = e_i ⇒ Cov(X) = 0 (zero matrix —
    // no variance, the channel is known). Row sums = 0 by simplex constraint
    // on the underlying X. The legacy `test(C_P_Cov, tolerance)` below
    // validates exactly this — restored after a misadventure where I
    // mistakenly reinterpreted the convention as `bare + diag(P_mean)` and
    // briefly removed it.

    template <bool output>
    static Maybe_error_t<bool> test_Probability_value(double e,
                                                      Probability_error_tolerance tolerance) {
        if (!std::isfinite(e)) {
            if constexpr (output)
                return error_message(" not finite value=" + std::to_string(e) + "\n");
            else
                return error_message("");
        } else if (e + tolerance() < 0) {
            if constexpr (output)
                return error_message(" negative prob=" + std::to_string(e) + "\n");
            else
                return error_message("");
        } else if (e - tolerance() > 1) {
            if constexpr (output)
                return error_message("  prob greater than one" + std::to_string(e) +
                                     " 1- prob=" + std::to_string(1 - e) + "\n");
            else
                return error_message("");
        } else
            return true;
    }

    template <bool output, class C_P_mean>
        requires(U<C_P_mean, P_mean>)
    static Maybe_error<bool> test(const C_P_mean& pp, Probability_error_tolerance tolerance) {
        auto p = primitive(pp());
        double sum = 0;
        for (std::size_t i = 0; i < p.size(); ++i) {
            auto Maybe_prob_value = test_Probability_value<output>(p[i], tolerance);
            if (!Maybe_prob_value)
                return Maybe_prob_value.error();
            sum += p[i];
        }
        if (std::abs(primitive(sum) - 1.0) < tolerance())
            return true;
        else if constexpr (output)
            return error_message("sum test sum=" + std::to_string(sum));
        else
            return error_message("");
    }

    template <bool output, class C_P_Cov>
        requires(U<C_P_Cov, P_Cov>)
    static Maybe_error<bool> test(const C_P_Cov& t_p, Probability_error_tolerance tolerance) {
        auto& p = primitive(t_p());

        for (std::size_t i = 0; i < p.nrows(); ++i) {
            for (std::size_t j = 0; j < p.ncols(); ++j) {
                if (auto pijt = test_Probability_value<output>(p(i, i), tolerance); !pijt) {
                    if constexpr (output)
                        return error_message(" at Pcov(" + std::to_string(i) + "," +
                                             std::to_string(j) + "):  " + pijt.error()());
                    else
                        return pijt.error();
                }
            }
            double sum = 0;
            for (std::size_t j = 0; j < p.ncols(); ++j) {
                if (i != j) {
                    if ((p(i, i) * p(j, j) - sqr(p(i, j)) + tolerance() < 0) &&
                        (p(i, i) > tolerance() * tolerance()) &&
                        (p(j, j) > tolerance() * tolerance())) {
                        if constexpr (output) {
                            double corr = sqr(p(i, j)) / p(i, i) / p(j, j);
                            std::stringstream ss;
                            ss << "tolerance=" << tolerance << "\n";
                            ss << " pcov=\n"
                               << p << "\n i=" << i << "j=" << j << " pcov(i,j)=" << p(i, j)
                               << " corr=" << corr << " pcov(i,i)" << p(i, i) << " pcov(j,j)"
                               << p(j, j) << "\n";
                            return error_message(ss.str());
                        }
                        return error_message("");
                    } else
                        sum += p(i, j);
                }
            }
            if (std::abs(p(i, i) + sum) > tolerance()) {
                if constexpr (output) {
                    std::stringstream ss;
                    ss << "tolerance=" << tolerance << "\n";
                    ss << " p=\n"
                       << p << "\n i=" << i << " p(i,j)=" << p(i, i) << " sum=" << sum << "\n";
                    return error_message(ss.str());
                }
                return error_message("");
            }
        }
        return true;
    }

    template <bool output, class C_P_mean, class C_P_Cov>
        requires(U<C_P_Cov, P_Cov> && U<C_P_mean, P_mean>)
    static Maybe_error<bool> test(const C_P_mean& t_P_mean, const C_P_Cov& t_P_cov,
                                  Probability_error_tolerance tolerance) {
        auto ck_mean = test<true>(t_P_mean, tolerance);
        auto ck_cov = test<true>(t_P_cov, tolerance);
        if (ck_mean && ck_cov)
            return true;
        else if constexpr (output) {
            std::stringstream ss;
            ss << " Pmean test: " << ck_mean.error()() << " Pcov test: " << ck_cov.error()();
            return error_message(ss.str());
        } else
            return error_message("");
    }

   public:
    static bool crude_Qx_violations(Qx const& q) {
        auto Qu = q() * Matrix<double>(q().ncols(), 1ul, 1.0);
        if (maxAbs(Qu) > 1e-7 * norm_inf(q())) {
            return true;
        } else
            return false;
    }

    template <class C_Patch_Model>
        requires U<C_Patch_Model, Patch_Model>
    auto calc_Qx(const C_Patch_Model& m, Agonist_concentration x)
        -> Transfer_Op_to<C_Patch_Model, Qx> {
        auto v_Qx = build<Qx>(get<Q0>(m)() + get<Qa>(m)() * x.value());
        Matrix<double> u(v_Qx().ncols(), 1, 1.0);
        v_Qx() = v_Qx() - diag(v_Qx() * u);
        assert(!crude_Qx_violations(primitive(v_Qx)));
        //         std::cerr<<"Qx violation\n";
        return v_Qx;
    }
    template <class C_Q0, class C_Qa>
        requires U<C_Q0, Q0>
    auto calc_Qx(const C_Q0& t_Q0, const C_Qa& t_Qa, Agonist_concentration x)
        -> Transfer_Op_to<C_Q0, Qx> {
        auto v_Qx = build<Qx>(t_Q0() + t_Qa() * x.value());
        Matrix<double> u(v_Qx().ncols(), 1, 1.0);
        v_Qx() = v_Qx() - diag(v_Qx() * u);
        assert(!crude_Qx_violations(primitive(v_Qx)));
        //     std::cerr<<"Qx violation\n";
        return v_Qx;
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Q0, class C_Qa>
        requires U<C_Q0, Q0>
    Maybe_error<Transfer_Op_to<C_Q0, P_initial>> calc_Pinitial(const C_Q0& t_Q0, const C_Qa& t_Qa,
                                                               Agonist_concentration x,
                                                               N_St nstates) {
        auto p0 = Matrix<double>(1ul, nstates(), 1.0 / nstates());
        auto t_Qx = calc_Qx(t_Q0, t_Qa, x);
        auto v_eig_Qx = calc_eigen(t_Qx);
        if (v_eig_Qx) {
            auto& landa = get<lambda>(v_eig_Qx.value())();
            auto& Vv = get<V>(v_eig_Qx.value())();
            auto& Wv = get<W>(v_eig_Qx.value())();
            auto i_landa = var::i_max(primitive(landa));
            auto ladt = primitive(landa) - primitive(landa);
            ladt[i_landa] = 1.0;
            auto p0Vv = p0 * Vv;
            auto p0Vladt = p0Vv * ladt;
            auto p0Vvladt = p0Vladt * Wv;

            auto Maybe_P = to_Probability(p0 * Vv * ladt * Wv);
            if (!Maybe_P)
                return Maybe_P.error();
            return build<P_initial>(std::move(Maybe_P.value()));

        } else {
            auto eP = expm_sure(t_Qx());
            if (!eP)
                return eP.error();
            auto Maybe_P = to_Transition_Probability<Policy>(eP.value());
            if (!Maybe_P)
                return Maybe_P.error();
            else {
                auto P = std::move(Maybe_P.value());
                auto Maybe_P2 = to_Transition_Probability<Policy>(P() * P());
                while (Maybe_P2.valid() && maxAbs(primitive(P() - Maybe_P2.value()())) > 1e-6) {
                    P = std::move(Maybe_P2.value());
                    Maybe_P2 = to_Transition_Probability<Policy>(P() * P());
                }
                if (!Maybe_P2)
                    return Maybe_P2.error();
                auto Maybe_P3 = to_Probability(p0 * Maybe_P2.value()());
                if (!Maybe_P3)
                    return Maybe_P3.error();
                return build<P_initial>(std::move(Maybe_P3.value()));
            }
        }
    }

    template <class C_Qx>
        requires U<C_Qx, Qx>
    auto calc_eigen(const C_Qx& v_Qx) -> Maybe_error<Transfer_Op_to<C_Qx, Eigs>> {
        // Compute eigen-decomposition of Qx
        auto maybe_eig = eigs(v_Qx());
        if (!maybe_eig) {
            return maybe_eig.error();
        }

        auto [r_lambda, r_V, r_VL] = std::move(maybe_eig.value());
        auto r_W = tr(std::move(r_VL));

        // Frobenius condition number κ_F(V) = ‖V‖_F · ‖W‖_F. Primitive only
        // (no derivative payload). Sets the FP-noise floor of P_ij after the
        // V·diag(exp(λdt))·W reconstruction; consumed downstream as the
        // pseudo-count strength for the Bayesian-shrinkage regularizer. See
        // theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md
        const auto& V_p = primitive(r_V);
        const auto& W_p = primitive(r_W);
        double V_norm_sq = 0.0;
        double W_norm_sq = 0.0;
        for (std::size_t i = 0; i < V_p.size(); ++i)
            V_norm_sq += V_p[i] * V_p[i];
        for (std::size_t i = 0; i < W_p.size(); ++i)
            W_norm_sq += W_p[i] * W_p[i];
        double r_kappa = std::sqrt(V_norm_sq) * std::sqrt(W_norm_sq);

        return build<Eigs>(build<lambda>(std::move(r_lambda)), build<V>(std::move(r_V)),
                           build<W>(std::move(r_W)), kappa_V(r_kappa));
    }

    template <class C_Patch_Model>
        requires U<C_Patch_Model, Patch_Model>
    auto calc_eigen(const C_Patch_Model& m, Agonist_concentration x)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Eigs>> {
        return calc_eigen(calc_Qx(m, x));
    }

    // ---------------------------------------------------------------------
    // Bayesian shrinkage of conditional conductance moments with v_g
    // (per-state, "Poisson endpoint") priors.
    //
    // Full derivation:
    //   theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md
    //
    // The raw conditional estimators
    //     gmean_ij = gtotal_ij     / P_ij,
    //     gsqr_ij  = gtotal_sqr_ij / P_ij
    // are well-posed when P_ij is well above the FP noise floor but become
    // dominated by rounding noise for rare (small-Q_ij·dt, multi-hop)
    // transitions where both numerator and denominator approach zero.
    //
    // v_g prior: in the no-data limit, the rare i→j transit averages the
    // two endpoint conductances. For a moment matrix M:
    //     prior_gmean_ij[i,j] = (g[i]  + g[j])  / 2
    //     prior_gsqr_ij[i,j]  = (g[i]² + g[j]²) / 2.
    // (The gvar prior emerges algebraically as ((g[i] − g[j])/2)².)
    //
    // Chosen over EB-from-marginals because the marginals are themselves
    // sums of the eigendecomposition-reconstructed gtotal_ij. When V is
    // ill-conditioned (large κ(V)), gtotal_ij is FP-contaminated exactly
    // when the prior is supposed to dominate — EB marginals would then
    // feed contaminated data back into the shrinkage target, defeating the
    // regularization. v_g priors depend only on per-state model parameters
    // and stay clean regardless of κ(V).
    //
    // The shrunken conditional moment is the conjugate-prior posterior:
    //     gmean_ij_shrunk = (gtotal_ij     + prior_gmean · ε) / (P + ε),
    //     gsqr_ij_shrunk  = (gtotal_sqr_ij + prior_gsqr  · ε) / (P + ε),
    // where ε = min_P_prior is the pseudo-count:
    //   - eig paths:          ε = ε_mach · κ_F(V)   (FP-noise floor of P)
    //   - taylor/schur paths: ε = 10·√N·ε_mach
    //
    // P_ij ≫ ε: identical to the raw quotient.
    // P_ij ≪ ε: shrinks smoothly to the marginal average.
    // Continuous, monotonic, derivative-clean throughout.
    //
    // Conditional variance follows from the second-moment identity
    //     gvar_ij = gsqr_ij_shrunk − gmean_ij_shrunk²
    // applied to the shrunken moments. No explicit gvar prior needed.
    // ---------------------------------------------------------------------

    // v_g-based ("endpoint") prior on gmean_ij and gsqr_ij. In the Poisson
    // limit a rare transit i→j averages the bracketing endpoint
    // conductances:
    //     prior_gmean_ij[i,j] = (g[i] + g[j]) / 2
    //     prior_gsqr_ij[i,j]  = (g[i]² + g[j]²) / 2
    // The gvar prior emerges algebraically:
    //     prior_gvar = prior_gsqr − prior_gmean² = ((g[i] − g[j])/2)²
    //
    // Chosen over EB-from-marginals because it does not depend on the
    // eigendecomposition's reconstruction of gtotal_ij: when V is
    // ill-conditioned (large κ(V)) the marginals of gtotal_ij are
    // FP-contaminated exactly when the prior is supposed to dominate,
    // which defeats the regularization. v_g priors stay clean regardless of
    // κ(V) — they depend only on per-state model parameters.
    template <class C_g>
        requires U<C_g, g>
    static auto gmean_ij_prior(const C_g& t_g) {
        const std::size_t N = t_g().size();
        Matrix<double> uT(1, N, 1.0);
        auto g_outer = t_g() * uT;                       // (i,j) = g[i]
        return build<gmean_ij>((g_outer + tr(g_outer)) * 0.5);
    }

    // Typed gvar_ij to fit calc_g_ij_bayes's (gtotal_sqr_ij, gvar_ij)
    // overload; semantically a prior on the second moment, not on the variance.
    template <class C_g>
        requires U<C_g, g>
    static auto gsqr_ij_prior(const C_g& t_g) {
        const std::size_t N = t_g().size();
        Matrix<double> uT(1, N, 1.0);
        auto g_sq = elemMult(t_g(), t_g());              // N×1, entries g[i]²
        auto gsq_outer = g_sq * uT;                       // (i,j) = g[i]²
        return build<gvar_ij>((gsq_outer + tr(gsq_outer)) * 0.5);
    }

    // Conjugate-prior posterior under pseudo-count form
    //     posterior = (data + prior · ε) / (P + ε).
    // Overloads:  (gtotal_ij,     gmean_ij)  →  gmean_ij_shrunk,
    //             (gtotal_sqr_ij, gvar_ij)   →  gsqr_ij_shrunk.
    template <class C_gtotal_ij, class C_P, class C_prior>
        requires (((U<C_gtotal_ij, gtotal_ij>     && U<C_prior, gmean_ij>) ||
                   (U<C_gtotal_ij, gtotal_sqr_ij> && U<C_prior, gvar_ij>))
                  && U<C_P, P>)
    static auto calc_g_ij_bayes(const C_gtotal_ij& t_gtotal_ij,
                                   const C_P& t_P,
                                   const C_prior& t_prior, double min_P_prior) {
        const std::size_t N = t_P().nrows();
        Matrix<double> UU(N, N, 1.0);
        auto num = t_gtotal_ij() + t_prior() * min_P_prior;
        auto den = t_P() + UU * min_P_prior;
        auto ratio = zip([](auto const& x, auto const& y) { return x / y; }, num, den);
        // build<T> needs a concrete tag type; dispatch on the prior's role.
        if constexpr (U<C_prior, gmean_ij>)
            return build<gmean_ij>(std::move(ratio));
        else
            return build<gvar_ij>(std::move(ratio));
    }

    template <class C_P_mean>
        requires U<C_P_mean, P_mean>

    static C_P_mean normalize(C_P_mean&& pp, double t_min_p) {
        using Trans = transformation_type_t<C_P_mean>;
        auto p = var::inside_out(pp());
        for (std::size_t i = 0; i < p.nrows(); ++i) {
            Op_t<Trans, double> sum = 0;
            for (std::size_t j = 0; j < p.ncols(); ++j) {
                if (primitive(p(i, j)) > 1.0 - t_min_p) {
                    for (std::size_t k = 0; k < p.ncols(); ++k) {
                        p(i, k) = (j == k) ? 1.0 + p(i, k) - p(i, k) : p(i, k) - p(i, k);
                    }
                    return C_P_mean(var::outside_in(p, pp.dx()));
                } else if (primitive(p(i, j)) < t_min_p)
                    p(i, j) = p(i, j) - p(i, j) + 0.0;
                else
                    sum = sum + p(i, j);
            }
            if (primitive(sum) != 1.0)
                for (std::size_t j = 0; j < p.ncols(); ++j) p(i, j) = p(i, j) / sum;
        }
        return C_P_mean(var::outside_in(p, pp.dx()));
    }

    static auto sample_Multinomial(mt_64i& mt, P_mean const t_P_mean, std::size_t N) {
        auto k = t_P_mean().size();
        N_channel_state out(Matrix<double>(1, k));
        std::size_t N_remaining = N;
        auto cumP = var::cumsum_reverse(t_P_mean());
        std::size_t i = 0;
        while (N_remaining > 0 && i + 1 < k)  // in case of numerical errors
        {
            const auto denominator = primitive(cumP[i]);
            const auto probability =
                denominator > 0.0 ? std::clamp(primitive(t_P_mean()[i]) / denominator, 0.0, 1.0)
                                  : 0.0;
            auto n = std::binomial_distribution<std::size_t>(N_remaining, probability)(mt);
            N_remaining -= n;
            out()[i] = static_cast<double>(n);
            i = i + 1;
        }
        out()[k - 1] = static_cast<double>(N_remaining);
        return out;
    }

    static Matrix<double> cumSum_reverse(Matrix<double> const& p) {
        Matrix<double> out = p;
        for (std::size_t i = 0; i < p.nrows(); ++i) {
            for (std::size_t j = p.ncols() > 0 ? p.ncols() - 1 : 0; j > 0; --j) {
                out(i, j - 1) += out(i, j);
            }
        }
        return out;
    }

    static auto sample_Multinomial(mt_64i& mt, P const t_P, N_channel_state N) {
        assert(t_P().nrows() == t_P().ncols());
        auto k = N().size();
        auto cumP = cumSum_reverse(t_P());
        N_channel_state out(Matrix<double>(1, k, 0.0));
        for (std::size_t i = 0; i < k; ++i) {
            auto N_remaining = static_cast<std::size_t>(std::round(N()[i]));
            for (std::size_t j = 0; j + 1 < k; ++j) {
                if (N_remaining > 0) {
                    const auto denominator = cumP(i, j);
                    const auto probability =
                        denominator > 0.0 ? std::clamp(t_P()(i, j) / denominator, 0.0, 1.0) : 0.0;
                    auto n = std::binomial_distribution<std::size_t>(N_remaining, probability)(mt);
                    N_remaining -= n;
                    out()[j] += static_cast<double>(n);
                }
            }
            out()[k - 1] += static_cast<double>(N_remaining);
        }
        return out;
    }

    template <class C_P_Cov>
        requires U<C_P_Cov, P_Cov>
    static C_P_Cov normalize(C_P_Cov&& p, double t_min_p) {
        for (std::size_t i = 0; i < primitive(p)().nrows(); ++i) {
            if (primitive(p)()(i, i) < t_min_p) {
                for (std::size_t j = 0; j < primitive(p)().ncols(); ++j) {
                    set(p(), i, j, 0.0);
                }
            }
        }
        return std::move(p);
    }

    template <class C_P>
        requires U<C_P, P>
    static C_P normalize(C_P&& p, double t_min_p) {
        using Trans = transformation_type_t<C_P>;

        for (std::size_t i = 0; i < p().nrows(); ++i) {
            Op_t<Trans, double> sumP = 0;
            for (std::size_t j = 0; j < p().ncols(); ++j)
                if (primitive(p()(i, j)) < t_min_p)
                    p()(i, j) = 0;
                else
                    sumP = sumP + p()(i, j);
            for (std::size_t j = 0; j < p().ncols(); ++j) p()(i, j) = p()(i, j) / sumP;
        }
        // std::cerr<<p;
        return std::move(p);
    }

    template <class Vs, class Patch_Model>
        requires Vs::is_vector_map_space
    Maybe_error<Eigs const*> get_eigen(Vs& buffer_calc, const Patch_Model& m,
                                       Agonist_concentration x) {
        auto Maybe_eigen =
            buffer_calc[var::Vector_Map<Eigs>{}][var::Vector_Space<Agonist_concentration>(x)];
        if (Maybe_eigen)
            return Maybe_eigen;
        else {
            auto Maybe_new_eigen = calc_eigen(m, x);
            if (Maybe_new_eigen) {
                buffer_calc[var::Vector_Map<Eigs>{}].emplace(x, std::move(Maybe_new_eigen.value()));
                return get_eigen(buffer_calc, m, x);
            } else
                return Maybe_new_eigen.error();
        }
    }

    template <class C_Patch_Model, class C_Qx_eig>
        requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx_eig, Eigs>)
    auto calc_Peq_(C_Qx_eig const& t_Eigs, const C_Patch_Model& m)
        -> Transfer_Op_to<C_Patch_Model, P_mean> {
        auto nstates = get<N_St>(m).value();
        auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);

        auto& landa = get<lambda>(t_Eigs)();
        auto& Vv = get<V>(t_Eigs)();
        auto& Wv = get<W>(t_Eigs)();
        auto ladt = get<lambda>(t_Eigs)() * 1e8;
        auto laexp = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            ladt);

        if constexpr (false) {
            std::cerr << "\np0\n" << p0;
            std::cerr << "\nlanda\n" << landa;
            std::cerr << "\nVv\n" << Vv;
            std::cerr << "\nWv\n" << Wv;
            std::cerr << "\nlaexp\n" << laexp;
            std::cerr << "\nWv*Vv\n" << Wv * Vv;
            std::cerr << "\nVv*landa*Wv\n" << Vv * landa * Wv;
            //  std::cerr<<"\nWv*landa*Vv\n"<<Wv*landa*Vv;
            std::cerr << "\nQx\n" << get<Qx>(t_Eigs);
        }

        return build<P_mean>(p0 * Vv * laexp * Wv);
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Qx>
        requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx, Qx>)
    auto calc_Peq(C_Qx const& t_Eigs, N_St nstates) -> Transfer_Op_to<C_Qx, P_mean> {
        auto p0 = Matrix<double>(1ul, nstates(), 1.0 / nstates());
        auto v_eig_Qx = calc_eigen<Policy>(t_Eigs);
        if (v_eig_Qx) {
            auto& landa = get<lambda>(v_eig_Qx.value())();
            auto& Vv = get<V>(v_eig_Qx.value())();
            auto& Wv = get<W>(v_eig_Qx.value())();
            auto ladt = get<lambda>(v_eig_Qx.value())() * 1e8;

            auto laexp = apply(
                [](auto const& x) {
                    using std::exp;
                    return exp(x);
                },
                ladt);
            if constexpr (false) {
                std::cerr << "\np0\n" << p0;
                std::cerr << "\nlanda\n" << landa;
                std::cerr << "\nVv\n" << Vv;
                std::cerr << "\nWv\n" << Wv;
                std::cerr << "\nlaexp\n" << laexp;
                std::cerr << "\nWv*Vv\n" << Wv * Vv;
                std::cerr << "\nVv*landa*Wv\n" << Vv * landa * Wv;
                //  std::cerr<<"\nWv*landa*Vv\n"<<Wv*landa*Vv;
                std::cerr << "\nQx\n" << get<Qx>(t_Eigs);
            }

            return build<P_mean>(p0 * Vv * laexp * Wv);

        } else {
            // std::cerr << "uses expm_sure\n";
            auto P = expm_sure<Policy>(t_Eigs());
            auto P2 = P * P;
            while (maxAbs(primitive(P - P2)) > 1e-9) {
                P = P2;
                P2 = P * P;
            }
            return build<P_mean>(p0 * P2);
        }
    }

    // template <class C_Patch_Model, class C_Qx_eig>
    //     requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx_eig, Eigs>)
    // auto calc_Peq(C_Qx_eig const &t_Eigs
    // , const C_Patch_Model &m)
    //     -> Transfer_Op_to<C_Patch_Model, P_mean> {
    //     auto nstates = get<N_St>(m).value();
    //     auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);

    //     auto &landa = get<lambda>(t_Eigs
    // )();
    //     auto &Vv = get<V>(t_Eigs
    // )();
    //     auto &Wv = get<W>(t_Eigs
    // )();
    //     auto laexp = DiagonalMatrix<double>(nstates, nstates, 0.0);
    //     for (std::size_t i = 0; i < nstates; ++i) {
    //         if (landa(i, i) == 0.0)
    //             laexp[i] = 1.0;
    //     }
    //     if constexpr (false) {
    //         std::cerr << "\np0\n" << p0;
    //         std::cerr << "\nlanda\n" << landa;
    //         std::cerr << "\nVv\n" << Vv;
    //         std::cerr << "\nWv\n" << Wv;
    //         std::cerr << "\nlaexp\n" << laexp;
    //         std::cerr << "\nWv*Vv\n" << Wv * Vv;
    //         std::cerr << "\nVv*landa*Wv\n" << Vv * landa * Wv;
    //         //  std::cerr<<"\nWv*landa*Vv\n"<<Wv*landa*Vv;
    //         std::cerr << "\nQx\n" << get<Qx>(t_Eigs
    // );
    //     }

    //     return build<P_mean>(p0 * Vv * laexp * Wv);
    // }

    template <class Policy = StabilizerPolicyEnabled, class Patch_Model>
    auto calc_P(const Patch_Model& m, const Eigs& t_Eigs, double dt, double t_min_P) {
        auto ladt = get<lambda>(t_Eigs)() * dt;

        auto exp_ladt = apply([](double x) { return std::exp(x); }, ladt);
        //    return normalize(P(get<V>(t_Eigs
        // )() * exp_ladt * get<W>(t_Eigs
        // )()),
        //    t_min_P);
        return to_Transition_Probability<Policy>(get<V>(t_Eigs)() * exp_ladt * get<W>(t_Eigs)());
    }

    template <class Policy = StabilizerPolicyEnabled, class Patch_Model>
    Maybe_error<Transfer_Op_to<Patch_Model, P>> calc_P(const Patch_Model& m, const Qx& t_Eigs,
                                                       double dt, double t_min_P) {
        auto t_eigenQx = calc_eigen(t_Eigs);
        if (t_eigenQx) {
            auto Maybe_P = calc_P<Policy>(m, t_eigenQx.value(), dt, t_min_P);
            if (Maybe_P)
                return Maybe_P;
        }
        //      return normalize(P(expm_sure(t_Eigs
        // () * dt)), t_min_P);
        auto Maybe_eP = expm_sure(t_Eigs() * dt);
        if (!Maybe_eP)
            return Maybe_eP.error();
        return to_Transition_Probability<Policy>(Maybe_eP.value());
    }

    template <class Patch_Model>
    auto calc_Qdt_old(const Patch_Model& m, const Eigs& t_Eigs, number_of_samples ns, double dt) {
        auto t_min_P = get<min_P>(m)();
        auto& v_g = get<g>(m);

        std::size_t N = t_Eigs[Var<V>{}]().ncols();

        auto ladt = t_Eigs[Var<lambda>{}]() * dt;

        auto exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            ladt);
        auto v_P = P(get<V>(t_Eigs)() * exp_ladt * get<W>(t_Eigs)());

        SymmetricMatrix<double> E2m(N, N);
        SymmetricMatrix<double> E2mb(N, N);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < i + 1; ++j)
                E2m.set(i, j, Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j], t_min_P));

        // build E2
        Matrix<double> WgV_E2(N, N);
        Matrix<double> WgV = get<W>(t_Eigs)() * diag(get<g>(m)()) * get<V>(t_Eigs)();

        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) WgV_E2(i, j) = WgV(i, j) * E2m(i, j);

        auto v_gtotal_ij = gtotal_ij(t_Eigs[Var<V>{}]() * WgV_E2 * t_Eigs[Var<W>{}]());

        Matrix<double> WgV_E3(N, N, 0.0);
        for (std::size_t n1 = 0; n1 < N; n1++)
            for (std::size_t n3 = 0; n3 < N; n3++)
                for (std::size_t n2 = 0; n2 < N; n2++) {
                    WgV_E3(n1, n3) += WgV(n1, n2) * WgV(n2, n3) *
                                      E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1], exp_ladt[n2],
                                         exp_ladt[n3], t_min_P);  // optimizable
                }

        auto v_gtotal_sqr_ij =
            gtotal_sqr_ij(t_Eigs[Var<V>{}]() * WgV_E3 * t_Eigs[Var<W>{}]() * 2.0);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                if (v_P()(i, j) == 0) {
                    v_gtotal_ij()(i, j) = 0;
                    v_gtotal_sqr_ij()(i, j) = 0;
                }

        auto U = Matrix<double>(1, N, 1.0);
        auto UU = Matrix<double>(N, N, 1.0);
        auto gmean_ij_p = X_plus_XT(v_g() * U) * (0.5);
        auto gvar_ij_p =
            (v_g() * U - apply([](double x) { return std::abs(x); }, tr(v_g() * U))) * (0.5);

        std::cerr << "\ngmean_ij_p=\n" << gmean_ij_p << "\ngvar_ij_p=\n" << gvar_ij_p << "\n";
        // std::cerr<<"\n UU=ņ"<<UU<<"\n";
        auto gmean_ij_tot = v_gtotal_ij() + gmean_ij_p * t_min_P;
        auto P_p = v_P() + UU * t_min_P;
        auto v_gmean_ij = gmean_ij(zip([](auto x, auto y) { return x / y; }, gmean_ij_tot, P_p));
        auto v_gtotal_var_ij =
            gtotal_var_ij(v_gtotal_sqr_ij() -
                          zip([](auto x, auto y) { return x * y; }, v_gtotal_ij(), v_gmean_ij()));
        auto gvar_ij_tot = v_gtotal_var_ij() + gvar_ij_p * t_min_P;
        auto v_gvar_ij = gvar_ij(zip([](auto x, auto y) { return x / y; }, gvar_ij_tot, P_p));
        Matrix<double> u(N, 1, 1.0);
        auto v_gmean_i = gmean_i(v_gtotal_ij() * u);
        auto v_gsqr_i = gsqr_i(v_gtotal_sqr_ij() * u);
        auto v_gvar_i = gvar_i(v_gtotal_var_ij() * u);

        return Qdt(ns, min_P(t_min_P), std::move(v_P), std::move(v_gmean_i), std::move(v_gtotal_ij),
                   std::move(v_gmean_ij), std::move(v_gtotal_sqr_ij), std::move(v_gsqr_i),
                   std::move(v_gvar_i), std::move(v_gtotal_var_ij), std::move(v_gvar_ij));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model, class C_Qx_eig>
   requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Eigs>)
     Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> calc_Qdtg_eig(FunctionTable&&,
                                                                   const C_Patch_Model& m,
                                                                   const C_Qx_eig& t_Eigs,
                                                                   number_of_samples ns,
                                                                   double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        // const double eps=std::numeric_limits<double>::epsilon();
        const auto& t_V = get<V>(t_Eigs);
        const auto& t_W = get<W>(t_Eigs);
        const auto& t_landa = get<lambda>(t_Eigs);
        const auto& t_g = get<g>(m);
        const auto t_min_P = get<min_P>(m);
        auto v_ladt = t_landa() * dt *0.5;
        auto v_exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            v_ladt);

        auto Psum = t_V() * v_exp_ladt * t_W();  // zeroed
        auto Maybe_r_P = to_Transition_Probability<Policy>(Psum);
        if (!Maybe_r_P) {
            return Maybe_r_P.error();
        }
        auto r_P = build<P_half>(std::move(Maybe_r_P.value())());
        if (std::isnan(primitive(var::max(r_P())))) {
            return error_message("nan P");
        }

        return build<Qdtg>(ns, min_P(t_min_P), std::move(r_P), t_g);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model,
              class C_Qx_eig>
        requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Eigs>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> calc_Qdtm_eig(FunctionTable&&,
                                                                   const C_Patch_Model& m,
                                                                   const C_Qx_eig& t_Qx,
                                                                   number_of_samples ns,
                                                                   double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        auto& t_V = get<V>(t_Qx);
        auto& t_landa = get<lambda>(t_Qx);
        auto& t_W = get<W>(t_Qx);
        auto& t_g = get<g>(m);
        auto t_min_P = get<min_P>(m);
        // Pseudo-count for source-level Bayesian shrinkage. ε_mach · κ_F(V) is
        // the FP noise floor of P_ij after V·diag(exp(λdt))·W reconstruction.
        const double min_P_prior = std::numeric_limits<double>::epsilon() * get<kappa_V>(t_Qx)();

        auto v_ladt = t_landa() * dt;
        auto v_exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            v_ladt);

        auto Maybe_r_P = to_Transition_Probability<Policy>(t_V() * v_exp_ladt * t_W());
        if (!Maybe_r_P)
            return Maybe_r_P.error();
        auto r_P = std::move(Maybe_r_P.value());
        std::size_t N = r_P().ncols();

        // Padé divided-difference Ee/E3 take a near-equality threshold; default
        // (ε_mach) is fine. The future Opitz/expm1 rewrite (Task 6) drops it.
        SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i + 1; ++j) {
                set(E2m, i, j, Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j]));
            }
        }

        Matrix<Op_t<Trans, double>> WgV_E2(N, N);
        auto v_WgV = t_W() * diag(t_g()) * t_V();
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

        auto r_gtotal_ij_raw = build<gtotal_ij>(t_V() * WgV_E2 * t_W());

        Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
        for (std::size_t n1 = 0; n1 < N; n1++)
            for (std::size_t n3 = 0; n3 < N; n3++)
                for (std::size_t n2 = 0; n2 < N; n2++) {
                    WgV_E3(n1, n3) = WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                                          E3(v_ladt[n1], v_ladt[n2], v_ladt[n3],
                                                             v_exp_ladt[n1], v_exp_ladt[n2],
                                                             v_exp_ladt[n3]);
                }

        auto r_gtotal_sqr_ij_raw = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);

        // v_g-based priors (Poisson endpoint, robust against ill-conditioned V).
        auto r_gmean_ij_prior = gmean_ij_prior(t_g);
        auto r_gsqr_ij_prior  = gsqr_ij_prior (t_g);

        auto r_gmean_ij = calc_g_ij_bayes(r_gtotal_ij_raw,     r_P, r_gmean_ij_prior, min_P_prior);
        auto r_gsqr_ij  = calc_g_ij_bayes(r_gtotal_sqr_ij_raw, r_P, r_gsqr_ij_prior,  min_P_prior);

        // Back-convert gtotal_ij so identity gtotal_ij = P · gmean_ij holds
        // exactly. Qdtm does not store gtotal_sqr_ij / gtotal_var_ij / gvar_ij,
        // but we materialize gtotal_sqr_ij transiently to (a) feed the canary
        // and (b) sum into gsqr_i.
        auto r_gtotal_ij     = build<gtotal_ij>    (elemMult(r_gmean_ij(), r_P()));
        auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(elemMult(r_gsqr_ij(),  r_P()));

        // Source-level canaries.
        if (auto chk = require_gtotal_ij_in_range(r_gtotal_ij, r_P, t_g); !chk)
            return chk.error();
        if (auto chk = require_gtotal_sqr_ij_in_range(r_gtotal_sqr_ij, r_P, t_g); !chk)
            return chk.error();

        // Row-sum marginals. r_gsqr_i is the P-weighted sum of the shrunken
        // gsqr_ij; together with the back-converted gmean_i it gives the full
        // LTV gvar_i = gsqr_i − gmean_i² expected by micro_* consumers.
        Matrix<double> u(N, 1, 1.0);
        auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
        auto r_gsqr_i  = build<gsqr_i> (r_gtotal_sqr_ij() * u);
        auto r_gvar_i  = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        if (std::isnan(primitive(var::max(r_P()))))
            return error_message("nan P");

        return build<Qdtm>(ns, t_min_P, std::move(r_P), std::move(r_gmean_i),
                           std::move(r_gtotal_ij), std::move(r_gmean_ij),
                           std::move(r_gsqr_i), std::move(r_gvar_i));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model,
              class C_Qx_eig>
        requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Eigs>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> calc_Qdt_eig(FunctionTable&&,
                                                                 const C_Patch_Model& m,
                                                                 const C_Qx_eig& t_Qx,
                                                                 number_of_samples ns, double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        auto& t_V = get<V>(t_Qx);
        auto& t_landa = get<lambda>(t_Qx);
        auto& t_W = get<W>(t_Qx);
        auto& t_g = get<g>(m);
        auto t_min_P = get<min_P>(m);
        // Pseudo-count for source-level Bayesian shrinkage. ε_mach · κ_F(V) is
        // the FP noise floor of P_ij after V·diag(exp(λdt))·W reconstruction.
        // See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        const double min_P_prior = std::numeric_limits<double>::epsilon() * get<kappa_V>(t_Qx)();

        auto v_ladt = t_landa() * dt;
        auto v_exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            v_ladt);

        auto Maybe_r_P = to_Transition_Probability<Policy>(t_V() * v_exp_ladt * t_W());
        if (!Maybe_r_P)
            return Maybe_r_P.error();
        auto r_P = std::move(Maybe_r_P.value());

        std::size_t N = r_P().ncols();

        // Ee / E3 fall back to their default eps = ε_mach for the
        // near-equality branch. min_P is not a dimensionally meaningful
        // threshold here; future Opitz/expm1 rewrite (Task 6) drops the
        // argument entirely.
        SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i + 1; ++j) {
                set(E2m, i, j,
                    Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j]));
            }
        }

        Matrix<Op_t<Trans, double>> WgV_E2(N, N);
        auto v_WgV = t_W() * diag(t_g()) * t_V();
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

        auto r_gtotal_ij_raw = build<gtotal_ij>(t_V() * WgV_E2 * t_W());

        Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
        for (std::size_t n1 = 0; n1 < N; n1++)
            for (std::size_t n3 = 0; n3 < N; n3++)
                for (std::size_t n2 = 0; n2 < N; n2++) {
                    WgV_E3(n1, n3) = WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                                          E3(v_ladt[n1], v_ladt[n2], v_ladt[n3],
                                                             v_exp_ladt[n1], v_exp_ladt[n2],
                                                             v_exp_ladt[n3]);
                }

        auto r_gtotal_sqr_ij_raw = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);

        // v_g-based priors (Poisson endpoint, robust against ill-conditioned V).
        auto r_gmean_ij_prior = gmean_ij_prior(t_g);
        auto r_gsqr_ij_prior  = gsqr_ij_prior (t_g);

        auto r_gmean_ij = calc_g_ij_bayes(r_gtotal_ij_raw,     r_P, r_gmean_ij_prior, min_P_prior);
        auto r_gsqr_ij  = calc_g_ij_bayes(r_gtotal_sqr_ij_raw, r_P, r_gsqr_ij_prior,  min_P_prior);

        // Conditional variance from the second-moment identity.
        auto r_gvar_ij = build<gvar_ij>(r_gsqr_ij() - elemMult(r_gmean_ij(), r_gmean_ij()));

        // Back-convert gtotal_* fields so identity gtotal = P · gmoment holds
        // exactly. Downstream IRT/MRT depends on this for
        //     sigma2_i = Σ_j P_ij · gvar_ij    (= gsqr_i − Σ_j P_ij · gmean_ij²).
        auto r_gtotal_ij     = build<gtotal_ij>    (elemMult(r_gmean_ij(), r_P()));
        auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(elemMult(r_gsqr_ij(),  r_P()));
        auto r_gtotal_var_ij = build<gtotal_var_ij>(elemMult(r_gvar_ij(),  r_P()));

        // Source-level canaries: gtotal_ij ∈ [gmin·P, gmax·P],
        // gtotal_sqr_ij ∈ [gmin²·P, gmax²·P]. Warn-band logs to stderr;
        // error-band aborts the Qdt construction.
        if (auto chk = require_gtotal_ij_in_range(r_gtotal_ij, r_P, t_g); !chk)
            return chk.error();
        if (auto chk = require_gtotal_sqr_ij_in_range(r_gtotal_sqr_ij, r_P, t_g); !chk)
            return chk.error();

        // Row-sum marginals over j (P-weighted, since gtotal_* carries P).
        Matrix<double> u(N, 1, 1.0);
        auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
        auto r_gsqr_i  = build<gsqr_i> (r_gtotal_sqr_ij() * u);
        // Full law-of-total-variance: Var(Ā|i) = E[Ā²|i] − E[Ā|i]² = gsqr_i − gmean_i².
        // This is what micro_full and micro_monoid expect when they read
        // get<gvar_i>(t_Qdt) / get<gvar_i>(t_Qdtm).
        auto r_gvar_i  = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        return build<Qdt>(ns, t_min_P, std::move(r_P), std::move(r_gmean_i),
                          std::move(r_gtotal_ij), std::move(r_gmean_ij),
                          std::move(r_gtotal_sqr_ij), std::move(r_gsqr_i), std::move(r_gvar_i),
                          std::move(r_gtotal_var_ij), std::move(r_gvar_ij));
    }

#if 0  // calc_Qdtm_eig_codex / calc_Qdt_eig_codex parked: alternative codex-
       // generated paths superseded by calc_Qdtm_eig / calc_Qdt_eig with v_g
       // priors and back-converted gtotal_*.
    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model,
              class C_Eigs>
        requires(U<C_Patch_Model, Patch_Model> && U<C_Eigs, Eigs>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> calc_Qdtm_eig_codex(FunctionTable&&,
                                                                         const C_Patch_Model& m,
                                                                         const C_Eigs& t_Eigs,
                                                                         number_of_samples ns,
                                                                         double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        // const double eps=std::numeric_limits<double>::epsilon();

        auto& Vmat = get<V>(t_Eigs);
        auto& Wmat = get<W>(t_Eigs);
        auto& Ldiag = get<lambda>(t_Eigs);
        auto& t_g = get<g>(m);
        // Bayesian-shrinkage pseudo-count for gmean_ij regularization.
        // See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        const double min_P_prior = std::numeric_limits<double>::epsilon() * get<kappa_V>(t_Eigs)();
        using GType = std::decay_t<decltype(t_g)>;
        using DX = var::dx_of_dfdx_t<GType>;
        auto const& dx = [&]() -> const DX& {
            if constexpr (var::is_derivative_v<GType>) {
                return var::get_dx_of_dfdx(t_g);
            } else {
                static const DX no{};
                return no;
            }
        }();
        auto t_min_P = get<min_P>(m);
        auto v_ladt = Ldiag() * dt;
        auto v_exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            v_ladt);

        auto Pacc = Vmat() * v_exp_ladt * Wmat();  // zeroed
        auto Maybe_r_P = to_Transition_Probability<Policy>(Pacc);
        if (!Maybe_r_P) {
            return Maybe_r_P.error();
        }
        auto r_P = std::move(Maybe_r_P.value());

        std::size_t N = r_P().ncols();

        SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i + 1; ++j) {
                set(E2m, i, j, Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j], t_min_P()));
            }
        }

        SymmetricMatrix<Op_t<Trans, double>> E2_sqr(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i + 1; ++j) {
                set(E2_sqr, i, j, E2(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j]));
            }
        }

        Matrix<Op_t<Trans, double>> WgV_E2(N, N);

        auto Win = var::inside_out(Wmat());
        auto Vin = var::inside_out(Vmat());
        Matrix<Op_t<Trans, double>> v_WgV(N, N, Op_t<Trans, double>(0.0));
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                Op_t<Trans, double> acc(0.0);
                for (std::size_t k = 0; k < N; ++k) {
                    acc = acc + Win(i, k) * t_g()(k, std::size_t(0)) * Vin(k, j);
                }
                v_WgV(i, j) = acc;
                WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);
            }
        }

        auto r_gtotal_ij = [&]() {
            using VMatType = std::decay_t<decltype(Vmat)>;
            auto mul_mm = [&](const auto& A, const auto& B) {
                Matrix<Op_t<Trans, double>> out(A.nrows(), B.ncols(), Op_t<Trans, double>(0.0));
                for (std::size_t i = 0; i < A.nrows(); ++i)
                    for (std::size_t j = 0; j < B.ncols(); ++j) {
                        Op_t<Trans, double> acc(0.0);
                        for (std::size_t k = 0; k < A.ncols(); ++k) acc = acc + A(i, k) * B(k, j);
                        out(i, j) = acc;
                    }
                return out;
            };
            if constexpr (var::is_derivative_v<VMatType>) {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E2);
                auto prod_in = mul_mm(tmp, Win_in);
                auto base = build<gtotal_ij>(var::outside_in(prod_in, dx));
                if constexpr (Policy::enforce_gmean_bounds)
                    return force_gmean_in_range(std::move(base), t_g);
                else
                    return base;
            } else {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E2);
                auto prod = mul_mm(tmp, Win_in);
                auto base = build<gtotal_ij>(prod);
                if constexpr (Policy::enforce_gmean_bounds)
                    return force_gmean_in_range(std::move(base), t_g);
                else
                    return base;
            }
        }();

        if constexpr (Policy::mask_probability) {
            constexpr double probability_eps = 1e-12;
            auto probability_mask = apply(
                [probability_eps](auto prob_entry) {
                    return prob_entry / (prob_entry + decltype(prob_entry)(probability_eps));
                },
                r_P());
            r_gtotal_ij() = elemMult(r_gtotal_ij(), probability_mask);
        }

        // Source-level Bayesian shrinkage at gtotal_ij.
        // theory/macroir/notes/Gmean_ij_gvarij/handoff_state.md
        auto gtotal_ij_prior_mat = gtotal_ij_endpoint_prior(t_g);
        r_gtotal_ij() = r_gtotal_ij() + gtotal_ij_prior_mat * min_P_prior;

        auto r_gmean_ij = [&]() {
            auto base = build<gmean_ij>(
                divide_by_P_plus_eps(r_gtotal_ij(), r_P(), min_P_prior));
            if constexpr (Policy::enforce_gmean_bounds)
                return force_gmean_in_range(std::move(base), t_g);
            else
                return base;
        }();
        // Smoothly clamp negative variance to preserve invariants while keeping derivatives
        Matrix<double> u(N, 1, 1.0);
        auto r_gmean_i = [&]() {
            auto base = build<gmean_i>(r_gtotal_ij() * u);
            if constexpr (Policy::enforce_gmean_bounds)
                return force_gmean_in_range(std::move(base), t_g);
            else
                return base;
        }();
        if (crude_gmean_violation(primitive(r_gmean_i), primitive(get<g>(m))))
            return error_message("gmean_violation");

        /**
 M_Matrix<double> WgV_Wg_E2(k_u,k_u);
    for (std::size_t k0=0; k0<k_u; k0++)
    {
	double rladt=Qx_v.landa[k0]*xdt.dt();
	if (tol.isEqual(rladt*rladt,0.0))
	{
	    for (std::size_t k2=0; k2<k_u; k2++)
	    {
		double rla2dt=Qx_v.landa[k2]*xdt.dt();
		if (tol.isEqual(rla2dt*rla2dt,0.0))
		    WgV_Wg_E2(k0,k2)=Qx_v.WgV(k0,k2)*
				     Qx_v.Wg[k2]*0.5;
		else
		    WgV_Wg_E2(k0,k2)=Qx_v.WgV(k0,k2)*Qx_v.Wg[k2]*
				     (exp(rla2dt)-rla2dt-1.0)/rla2dt/rla2dt;
	    }
	}   */

        Matrix<Op_t<Trans, double>> WgV_Wg_E2(N, 1, 0.0);

        auto v_Wg = Wmat() * t_g();

        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                WgV_Wg_E2(i, 0) = WgV_Wg_E2(i, 0) + v_WgV(i, j) * E2_sqr(i, j) * v_Wg[j];

        Matrix<Op_t<Trans, double>> rgsqr_i(N, 1, 0.0);

        for (std::size_t i = 0; i < N; i++) {
            for (std::size_t k0 = 0; k0 < N; k0++) {
                auto term = Vmat()(i, k0) * WgV_Wg_E2[k0];
                rgsqr_i[i] = rgsqr_i[i] + term + term;
            }
        }

        auto r_gsqr_i = build<gsqr_i>(var::outside_in(rgsqr_i, dx));

        auto r_gvar_i = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        /* truncate is not derivative safe yet*/

        if constexpr (Policy::clamp_variance) {
            r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
        }

        if constexpr (Policy::mask_probability) {
            constexpr double probability_eps = 1e-12;
            auto probability_mask = apply(
                [probability_eps](auto prob_entry) {
                    return prob_entry / (prob_entry + decltype(prob_entry)(probability_eps));
                },
                r_P());
            using MaskEntry = std::decay_t<decltype(probability_mask(0, 0))>;
            Matrix<MaskEntry> row_indicator(r_P().nrows(), 1, MaskEntry(0.0));
            for (std::size_t i_row = 0; i_row < r_P().nrows(); ++i_row) {
                MaskEntry accum(0.0);
                for (std::size_t j_col = 0; j_col < r_P().ncols(); ++j_col)
                    accum = accum + probability_mask(i_row, j_col);
                auto denom = accum + MaskEntry(probability_eps);
                row_indicator(i_row, 0) = accum / denom;
            }
            auto gsqr_inside = var::inside_out(r_gsqr_i());
            auto gvar_inside = var::inside_out(r_gvar_i());
            for (std::size_t i_row = 0; i_row < gsqr_inside.nrows(); ++i_row) {
                auto indicator = row_indicator(i_row, 0);
                gsqr_inside.set(i_row, 0, gsqr_inside(i_row, 0) * indicator);
                gvar_inside.set(i_row, 0, gvar_inside(i_row, 0) * indicator);
            }
            r_gsqr_i() = var::outside_in(std::move(gsqr_inside), dx);
            r_gvar_i() = var::outside_in(std::move(gvar_inside), dx);
        }
        //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
        //   primitive(t_g())); auto test_g_mean =
        //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
        // if (!test_g_mean || !test_g_var)
        //     return error_message(test_g_var.error()()+test_g_mean.error()());
        if (std::isnan(primitive(var::max(r_P()))))
            return error_message("nan P");

        return build<Qdtm>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                           std::move(r_gtotal_ij), std::move(r_gmean_ij), std::move(r_gsqr_i),
                           std::move(r_gvar_i));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model,
              class C_Qx_eig>
        requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Eigs>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> calc_Qdt_eig_codex(FunctionTable&&,
                                                                       const C_Patch_Model& m,
                                                                       const C_Qx_eig& t_Eigs,
                                                                       number_of_samples ns,
                                                                       double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        // const double eps=std::numeric_limits<double>::epsilon();
        auto& Vmat = get<V>(t_Eigs);
        auto& Ldiag = get<lambda>(t_Eigs);
        auto& Wmat = get<W>(t_Eigs);
        auto& t_g = get<g>(m);
        // Bayesian-shrinkage pseudo-count for gmean_ij / gvar_ij.
        // See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        const double min_P_prior = std::numeric_limits<double>::epsilon() * get<kappa_V>(t_Eigs)();
        using GType = std::decay_t<decltype(t_g)>;
        using DX = var::dx_of_dfdx_t<GType>;
        auto const& dx = [&]() -> const DX& {
            if constexpr (var::is_derivative_v<GType>) {
                return var::get_dx_of_dfdx(t_g);
            } else {
                static const DX no{};
                return no;
            }
        }();
        auto t_min_P = get<min_P>(m);
        auto v_ladt = Ldiag() * dt;
        auto v_exp_ladt = apply(
            [](auto const& x) {
                using std::exp;
                return exp(x);
            },
            v_ladt);

        // Build P via block masks
        auto Ncols = Vmat().ncols();
        auto Pacc = Vmat() * v_exp_ladt * Wmat();  // zeroed

        auto Maybe_r_P = to_Transition_Probability<Policy>(Pacc);
        if (!Maybe_r_P) {
            return Maybe_r_P.error();
        }

        auto r_P = std::move(Maybe_r_P.value());

        std::size_t N = r_P().ncols();

        SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i + 1; ++j) {
                set(E2m, i, j, Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j], t_min_P()));
            }
        }

        Matrix<Op_t<Trans, double>> WgV_E2(N, N);

        auto Win = var::inside_out(Wmat());
        auto Vin = var::inside_out(Vmat());
        Matrix<Op_t<Trans, double>> v_WgV(N, N, Op_t<Trans, double>(0.0));
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                Op_t<Trans, double> acc(0.0);
                for (std::size_t k = 0; k < N; ++k) {
                    acc = acc + Win(i, k) * t_g()(k, std::size_t(0)) * Vin(k, j);
                }
                v_WgV(i, j) = acc;
                WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);
            }
        }

        auto r_gtotal_ij = [&]() {
            using VMatType = std::decay_t<decltype(Vmat)>;
            auto mul_mm = [&](const auto& A, const auto& B) {
                Matrix<Op_t<Trans, double>> out(A.nrows(), B.ncols(), Op_t<Trans, double>(0.0));
                for (std::size_t i = 0; i < A.nrows(); ++i)
                    for (std::size_t j = 0; j < B.ncols(); ++j) {
                        Op_t<Trans, double> acc(0.0);
                        for (std::size_t k = 0; k < A.ncols(); ++k) acc = acc + A(i, k) * B(k, j);
                        out(i, j) = acc;
                    }
                return out;
            };
            if constexpr (var::is_derivative_v<VMatType>) {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E2);
                auto prod_in = mul_mm(tmp, Win_in);
                auto base = build<gtotal_ij>(var::outside_in(prod_in, dx));
                if constexpr (Policy::enforce_gmean_bounds)
                    return force_gmean_in_range(std::move(base), t_g);
                else
                    return base;
            } else {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E2);
                auto prod = mul_mm(tmp, Win_in);
                auto base = build<gtotal_ij>(prod);
                if constexpr (Policy::enforce_gmean_bounds)
                    return force_gmean_in_range(std::move(base), t_g);
                else
                    return base;
            }
        }();

        Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
        for (std::size_t n1 = 0; n1 < N; n1++)
            for (std::size_t n3 = 0; n3 < N; n3++)
                for (std::size_t n2 = 0; n2 < N; n2++) {
                    //      std::cerr<<"\t"<<WgV_E3(n1, n3);

                    WgV_E3(n1, n3) =
                        WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                             E3(v_ladt[n1], v_ladt[n2], v_ladt[n3], v_exp_ladt[n1],
                                                v_exp_ladt[n2], v_exp_ladt[n3], t_min_P());
                }

        auto r_gtotal_sqr_ij = [&]() {
            using VMatType = std::decay_t<decltype(Vmat)>;
            auto mul_mm = [&](const auto& A, const auto& B) {
                Matrix<Op_t<Trans, double>> out(A.nrows(), B.ncols(), Op_t<Trans, double>(0.0));
                for (std::size_t i = 0; i < A.nrows(); ++i)
                    for (std::size_t j = 0; j < B.ncols(); ++j) {
                        Op_t<Trans, double> acc(0.0);
                        for (std::size_t k = 0; k < A.ncols(); ++k) acc = acc + A(i, k) * B(k, j);
                        out(i, j) = acc;
                    }
                return out;
            };
            if constexpr (var::is_derivative_v<VMatType>) {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E3);
                auto prod3_in = mul_mm(tmp, Win_in);
                // scale by 2.0
                for (std::size_t i = 0; i < prod3_in.nrows(); ++i)
                    for (std::size_t j = 0; j < prod3_in.ncols(); ++j)
                        prod3_in.set(i, j, prod3_in(i, j) + prod3_in(i, j));
                return build<gtotal_sqr_ij>(var::outside_in(prod3_in, dx));
            } else {
                auto Vin_in = var::inside_out(Vmat());
                auto Win_in = var::inside_out(Wmat());
                auto tmp = mul_mm(Vin_in, WgV_E3);
                auto prod3 = mul_mm(tmp, Win_in);
                for (std::size_t i = 0; i < prod3.nrows(); ++i)
                    for (std::size_t j = 0; j < prod3.ncols(); ++j)
                        prod3.set(i, j, prod3(i, j) + prod3(i, j));
                return build<gtotal_sqr_ij>(prod3);
            }
        }();

        if constexpr (false) {
            std::cerr << "\nr_gtotal_sqr_ij\n" << r_gtotal_sqr_ij;
            std::cerr << "\nvar::outside_in(var::inside_out(r_gtotal_sqr_ij))\n"
                      << var::outside_in(var::inside_out(r_gtotal_sqr_ij()));

            std::cerr << "\nvar::inside_out(r_gtotal_sqr_ij)\n"
                      << var::inside_out(r_gtotal_sqr_ij());
        }

        if constexpr (Policy::clamp_variance) {
            r_gtotal_sqr_ij() = truncate_negative_variance(std::move(r_gtotal_sqr_ij()));
        }

        if constexpr (Policy::mask_probability) {
            constexpr double probability_eps = 1e-12;
            auto apply_probability_mask = [&](auto& matrix) {
                using MatrixType = std::decay_t<decltype(matrix)>;
                if constexpr (var::is_derivative_v<MatrixType>) {
                    auto matrix_inside = var::inside_out(matrix);
                    auto probability_inside = var::inside_out(r_P());
                    for (std::size_t i_mask = 0; i_mask < matrix_inside.nrows(); ++i_mask)
                        for (std::size_t j_mask = 0; j_mask < matrix_inside.ncols(); ++j_mask) {
                            auto entry = matrix_inside(i_mask, j_mask);
                            auto prob_entry = probability_inside(i_mask, j_mask);
                            auto denom = prob_entry + decltype(prob_entry)(probability_eps);
                            auto indicator = prob_entry / denom;
                            matrix_inside.set(i_mask, j_mask, entry * indicator);
                        }
                    matrix = var::outside_in(std::move(matrix_inside), dx);
                } else {
                    auto indicator = apply(
                        [probability_eps](auto prob_value) {
                            return prob_value / (prob_value + probability_eps);
                        },
                        r_P());
                    matrix = elemMult(matrix, indicator);
                }
            };
            apply_probability_mask(r_gtotal_ij());
            apply_probability_mask(r_gtotal_sqr_ij());
        }

        // Source-level Bayesian shrinkage at gtotal_ij and gtotal_sqr_ij.
        // gvar_ij prior emerges algebraically. See
        // theory/macroir/notes/Gmean_ij_gvarij/handoff_state.md
        auto gtotal_ij_prior_mat = gtotal_ij_endpoint_prior(t_g);
        r_gtotal_ij() = r_gtotal_ij() + gtotal_ij_prior_mat * min_P_prior;
        auto gtotal_sqr_ij_prior_mat = gtotal_sqr_ij_endpoint_prior(t_g);
        r_gtotal_sqr_ij() = r_gtotal_sqr_ij() + gtotal_sqr_ij_prior_mat * min_P_prior;

        auto r_gmean_ij = [&]() {
            auto base = build<gmean_ij>(
                divide_by_P_plus_eps(r_gtotal_ij(), r_P(), min_P_prior));
            if constexpr (Policy::enforce_gmean_bounds)
                return force_gmean_in_range(std::move(base), t_g);
            else
                return base;
        }();
        auto r_gtotal_var_ij =
            build<gtotal_var_ij>(r_gtotal_sqr_ij() - elemMult(r_gtotal_ij(), r_gmean_ij()));

        // Smoothly clamp negative variance to preserve invariants while keeping derivatives
        if constexpr (Policy::clamp_variance) {
            r_gtotal_var_ij() = truncate_negative_variance(std::move(r_gtotal_var_ij()));
        }

        auto r_gvar_ij = build<gvar_ij>(
            divide_by_P_plus_eps(r_gtotal_var_ij(), r_P(), min_P_prior));

        Matrix<double> u(N, 1, 1.0);
        auto r_gmean_i = [&]() {
            auto base = build<gmean_i>(r_gtotal_ij() * u);
            if constexpr (Policy::enforce_gmean_bounds)
                return force_gmean_in_range(std::move(base), t_g);
            else
                return base;
        }();
        if (crude_gmean_violation(primitive(r_gmean_i), primitive(get<g>(m))))
            return error_message("gmean_violation");

        auto r_gsqr_i = build<gsqr_i>(r_gtotal_sqr_ij() * u);
        auto r_gvar_i = build<gvar_i>(r_gtotal_var_ij() * u);
        if constexpr (Policy::clamp_variance) {
            r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
        }
        //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
        //   primitive(t_g())); auto test_g_mean =
        //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
        // if (!test_g_mean || !test_g_var)
        //     return error_message(test_g_var.error()()+test_g_mean.error()());
        // else {

        return build<Qdt>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                          std::move(r_gtotal_ij), std::move(r_gmean_ij), std::move(r_gtotal_sqr_ij),
                          std::move(r_gsqr_i), std::move(r_gvar_i), std::move(r_gtotal_var_ij),
                          std::move(r_gvar_ij));
        // }
    }
#endif // calc_Qdtm_eig_codex / calc_Qdt_eig_codex parked

    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(/*U<C_Patch_Model, Patch_Model> && */ U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> calc_Qdt_taylor(const C_Patch_Model& m,
                                                                    const C_Qx& t_Eigs,
                                                                    number_of_samples ns, double dt,
                                                                    std::size_t order = 6ul) {
        auto v_Qrun = t_Eigs() * dt;
        double max = maxAbs(primitive(v_Qrun));
        double desired = 0.125 / 4.0;
        int k = std::ceil(std::log2(max / desired));
        int n = std::max(0, k);
        double scale = std::pow(2, -n);
        auto t_Qrun_sub = v_Qrun * scale;
        auto Maybe_P_sub = to_Transition_Probability<Policy>(expm_taylor(t_Qrun_sub, order));
        if (!Maybe_P_sub) {
            return Maybe_P_sub.error();
        } else {
            auto sub_ns = number_of_samples(ns() * scale);
            auto P_sub = std::move(Maybe_P_sub.value());
            auto r_Qn = get_Qn_via_taylor_integrals<Policy>(t_Qrun_sub, P_sub, get<g>(m),
                                                             sub_ns, get<min_P>(m), order);
            for (std::size_t i = 0; i < n; ++i) {
                auto Maybe_r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
                if (!Maybe_r_Qn)
                    return Maybe_r_Qn.error();
                r_Qn = std::move(Maybe_r_Qn.value());
            }
            assert(get<number_of_samples>(r_Qn) == ns);
            return Qn_to_Qdt(r_Qn, m);
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(/*U<C_Patch_Model, Patch_Model> && */
                 U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> calc_Qdtg_taylor(const C_Patch_Model& m,
                                                                      const C_Qx& t_Eigs,
                                                                      number_of_samples ns,
                                                                      double dt,
                                                                      std::size_t order = 5ul) {
        auto v_Qrun = t_Eigs() * dt*0.5;
        double max = maxAbs(primitive(v_Qrun));
        double desired = 0.125 / 4.0;
        int n;
        if (max > 0) {
            int k = std::ceil(std::log2(max / desired));
            n = std::max(0, k);
        } else
            n = 0;
        double scale = std::pow(2, -n);
        auto t_Qrun_sub = v_Qrun * scale;
        auto Maybe_P_sub = to_Transition_Probability<Policy>(expm_taylor(t_Qrun_sub, order));
        if (!Maybe_P_sub) {
            return Maybe_P_sub.error();
        }
        // Square n times to recover P_half = expm(Q·dt/2) from the scaled
        // sub-step P_sub = expm(Q·dt/2 · 2⁻ⁿ). Same scaling-and-squaring
        // recipe used by expm_sure (line ~1665). Without these squarings the
        // returned P_half is (P_half)^(1/2ⁿ), wrong whenever n>0 — which
        // turns the lifted micro_R likelihood completely wrong.
        auto P_half_mat = std::move(Maybe_P_sub.value());
        for (int i = 0; i < n; ++i) {
            auto Maybe_P2 = to_Transition_Probability<Policy>(P_half_mat() * P_half_mat());
            if (!Maybe_P2) return Maybe_P2.error();
            P_half_mat = std::move(Maybe_P2.value());
        }
        return build<Qdtg>(ns, get<min_P>(m), build<P_half>(std::move(P_half_mat())), get<g>(m));
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(/*U<C_Patch_Model, Patch_Model> && */ U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> calc_Qdtm_taylor(const C_Patch_Model& m,
                                                                      const C_Qx& t_Eigs,
                                                                      number_of_samples ns,
                                                                      double dt,
                                                                      std::size_t order = 5ul) {
        auto v_Qrun = t_Eigs() * dt;
        double max = maxAbs(primitive(v_Qrun));
        double desired = 0.125 / 4.0;
        int n;
        if (max > 0) {
            int k = std::ceil(std::log2(max / desired));
            n = std::max(0, k);
        } else
            n = 0;
        double scale = std::pow(2, -n);
        auto t_Qrun_sub = v_Qrun * scale;
        auto Maybe_P_sub = to_Transition_Probability<Policy>(expm_taylor(t_Qrun_sub, order));
        if (!Maybe_P_sub) {
            return Maybe_P_sub.error();
        } else {
            auto sub_ns = number_of_samples(ns() * scale);
            auto P_sub = std::move(Maybe_P_sub.value());
            auto r_Qn = get_Qn_via_taylor_integrals<Policy>(t_Qrun_sub, P_sub, get<g>(m),
                                                             sub_ns, get<min_P>(m), order);
            for (std::size_t i = 0; i < n; ++i) {
                auto Maybe_r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
                if (!Maybe_r_Qn)
                    return Maybe_r_Qn.error();
                r_Qn = std::move(Maybe_r_Qn.value());
            }
            assert(get<number_of_samples>(r_Qn) == ns);
            return Qn_to_Qdtm(r_Qn, m);
        }
    }

    // -------------------------------------------------------------------
    // Schur+Parlett path. Plain-double only at step 7; Derivative-aware
    // overload lives further below (step 8). Mirrors the (m, t_Qx, ns, dt)
    // signature of calc_Qdt_eig / calc_Qdt_taylor so dispatch can swap
    // implementations transparently.
    // -------------------------------------------------------------------
    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> calc_Qdtg_schur(const C_Patch_Model& m,
                                                                      const C_Qx& t_Qx,
                                                                      number_of_samples ns,
                                                                      double dt) {
        using GType = std::decay_t<decltype(get<g>(m))>;
        if constexpr (var::is_derivative_v<GType>) {
            // Derivative path: delegate to calc_Qdtg_taylor. With the
            // squaring fix in calc_Qdtg_taylor (qmodel.h:3146-3151), Taylor
            // produces correct Derivative<Qdtg> using N×N Derivative<Matrix>
            // arithmetic — much cheaper than Pade scaling-and-squaring.
            return calc_Qdtg_taylor<Policy>(m, t_Qx, ns, dt);
        } else {
            auto Maybe_P_half = lapack::expm_schur_parlett(t_Qx() * (dt * 0.5));
            if (!Maybe_P_half) return Maybe_P_half.error();
            auto Maybe_r_P = to_Transition_Probability<Policy>(Maybe_P_half.value());
            if (!Maybe_r_P) return Maybe_r_P.error();
            auto r_P = build<P_half>(std::move(Maybe_r_P.value())());
            return build<Qdtg>(ns, get<min_P>(m), std::move(r_P), get<g>(m));
        }
    }

    // Helper: assemble the full Qdt struct from P, gtotal_ij, gtotal_sqr_ij
    // already computed via Schur+VanLoan. Mirrors the post-spectral algebra
    // in calc_Qdt_eig (lines 2576-2621): gmean_ij = gtotal_ij / P,
    // gtotal_var_ij = gtotal_sqr_ij - gtotal_ij·gmean_ij, gvar_ij = .../P,
    // and the row-sum {gmean_i, gsqr_i, gvar_i}. Templated on the matrix
    // type so plain double and Derivative<Matrix<double>> share the same
    // post-expm algebra — the codebase's overloaded operators
    // (elemDivSafe, elemMult, *) propagate derivatives transparently.
    template <class Policy, class C_Patch_Model, class C_Mat>
    auto assemble_Qdt_from_moments(const C_Patch_Model& m, number_of_samples ns,
                                    C_Mat r_P_mat, C_Mat r_gtotal_ij_mat,
                                    C_Mat r_gtotal_sqr_ij_mat) {
        auto t_g = get<g>(m);
        auto t_min_P = get<min_P>(m);

        auto Maybe_r_P = to_Transition_Probability<Policy>(std::move(r_P_mat));
        if (!Maybe_r_P)
            return Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>(Maybe_r_P.error());
        auto r_P = std::move(Maybe_r_P.value());

        auto r_gtotal_ij     = build<gtotal_ij>    (std::move(r_gtotal_ij_mat));
        auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(std::move(r_gtotal_sqr_ij_mat));

        // Bayesian-shrinkage pseudo-count for gmean_ij / gvar_ij. Schur+Padé
        // path is backward-stable on its own orthogonal basis (no κ(V)
        // amplification); the noise floor is the random-walk FP accumulation
        // 10·√N·ε_mach across the N×N matmul chain. Parameter-free.
        // See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        const std::size_t N = r_P().nrows();
        const double min_P_prior = 10.0 * std::sqrt(static_cast<double>(N)) *
                                   std::numeric_limits<double>::epsilon();

        // v_g-based Bayesian shrinkage on both conditional moments via
        // calc_g_ij_bayes, then back-convert gtotal_* = gmoment · P so
        // downstream identities hold. Matches calc_Qdt_eig.
        // See theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        auto r_gmean_ij_prior = gmean_ij_prior(t_g);
        auto r_gsqr_ij_prior  = gsqr_ij_prior (t_g);

        auto r_gmean_ij = calc_g_ij_bayes(r_gtotal_ij,     r_P, r_gmean_ij_prior, min_P_prior);
        auto r_gsqr_ij  = calc_g_ij_bayes(r_gtotal_sqr_ij, r_P, r_gsqr_ij_prior,  min_P_prior);

        auto r_gvar_ij = build<gvar_ij>(r_gsqr_ij() - elemMult(r_gmean_ij(), r_gmean_ij()));

        r_gtotal_ij()     = elemMult(r_gmean_ij(), r_P());
        r_gtotal_sqr_ij() = elemMult(r_gsqr_ij(),  r_P());
        auto r_gtotal_var_ij = build<gtotal_var_ij>(elemMult(r_gvar_ij(), r_P()));

        // Source-level canaries.
        if (auto chk = require_gtotal_ij_in_range(r_gtotal_ij, r_P, t_g); !chk)
            return Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>(chk.error());
        if (auto chk = require_gtotal_sqr_ij_in_range(r_gtotal_sqr_ij, r_P, t_g); !chk)
            return Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>(chk.error());

        Matrix<double> u(r_P().nrows(), 1, 1.0);
        auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
        auto r_gsqr_i  = build<gsqr_i> (r_gtotal_sqr_ij() * u);
        // Full law-of-total-variance: Var(Ā|i) = gsqr_i − gmean_i².
        auto r_gvar_i  = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        return Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>(
            build<Qdt>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                        std::move(r_gtotal_ij), std::move(r_gmean_ij),
                        std::move(r_gtotal_sqr_ij), std::move(r_gsqr_i), std::move(r_gvar_i),
                        std::move(r_gtotal_var_ij), std::move(r_gvar_ij)));
    }

    // Build diag(g_vec) as a full N×N matrix, templated on the inner matrix
    // type. Pass the *unwrapped* matrix (g_vec is Matrix<double> or
    // Derivative<Matrix<double>>, not the g/Derivative<g> wrapper). For
    // Derivative inputs it uses inside_out / outside_in to repack — keeps
    // the Frechet integral builder agnostic to plain vs Derivative.
    template <class C_Mat>
        requires(var::U<C_Mat, Matrix<double>>)
    auto build_g_as_diag_matrix(C_Mat const& g_vec) {
        auto g_inner = var::inside_out(g_vec);
        using scalar_t = std::decay_t<decltype(g_inner(std::size_t{0}, std::size_t{0}))>;
        const std::size_t N = g_inner.nrows();
        scalar_t zero = g_inner(std::size_t{0}, std::size_t{0}) -
                         g_inner(std::size_t{0}, std::size_t{0});
        Matrix<scalar_t> D(N, N, zero);
        for (std::size_t i = 0; i < N; ++i) D(i, i) = g_inner(i, std::size_t{0});

        if constexpr (var::is_derivative_v<C_Mat>) {
            return var::outside_in(D, g_vec.dx());
        } else {
            return D;
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> calc_Qdt_schur(const C_Patch_Model& m,
                                                                    const C_Qx& t_Qx,
                                                                    number_of_samples ns,
                                                                    double dt) {
        using GType = std::decay_t<decltype(get<g>(m))>;
        if constexpr (var::is_derivative_v<GType>) {
            // Derivative path: delegate to calc_Qdt_taylor. With the seed
            // fix in get_Qn_via_taylor_integrals (qmodel.h:3215+), Taylor
            // produces correct Derivative<Qdt> — validated by
            // [macroir][parity][derivative]. Critically, Taylor does its
            // moment-integral arithmetic on N×N Derivative<Matrix> rather
            // than on the 3N×3N Van Loan augmented matrix that Pade-VanLoan
            // requires. For N=51 (Nch=50, k=2), p=6 parameters, that's a
            // ~10× speedup with zero accuracy loss — every Derivative
            // matmul on the 3N×3N augmented matrix is 27× more expensive
            // than the same operation on N×N (cost scales as N³ × (1+2p)).
            return calc_Qdt_taylor<Policy>(m, t_Qx, ns, dt);
        } else {
            // Plain double path. ONE 3N×3N Pade scaling-and-squaring on the
            // Van Loan augmented matrix [[Q,G,0],[0,Q,G],[0,0,Q]]·dt yields
            // P, A and B simultaneously (the (0,0), (0,1), (0,2) blocks of
            // the resulting expm). Replaces the three separate expm calls
            // (N×N for P, 2N×2N for A, 3N×3N for B) with one — saves the
            // N×N and 2N×2N expm work entirely (~25%).
            const auto& t_g = get<g>(m);
            auto G_diag = build_g_as_diag_matrix(t_g());

            auto Maybe_PAB = lapack::frechet_p_a_b_combined(t_Qx(), G_diag, dt);
            if (!Maybe_PAB)
                return Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>(Maybe_PAB.error());
            auto& [P_mat, A_mat, B_mat] = Maybe_PAB.value();

            // gtotal_ij = A / dt; gtotal_sqr_ij = 2·B / dt². The "·2.0"
            // matches the symmetrization factor calc_Qdt_eig applies at
            // qmodel.h:2591, since our B is a single-ordering triangle
            // integral.
            auto r_gtotal_ij_mat = A_mat * (1.0 / dt);
            auto r_gtotal_sqr_ij_mat = B_mat * (2.0 / (dt * dt));

            return assemble_Qdt_from_moments<Policy>(m, ns, std::move(P_mat),
                                                      std::move(r_gtotal_ij_mat),
                                                      std::move(r_gtotal_sqr_ij_mat));
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Patch_Model, class C_Qx>
        requires(U<C_Qx, Qx>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> calc_Qdtm_schur(const C_Patch_Model& m,
                                                                      const C_Qx& t_Qx,
                                                                      number_of_samples ns,
                                                                      double dt) {
        // Compute the full Qdt then project to Qdtm (= Qdt minus the *_ij
        // variance fields). Cheaper to share the full computation than to
        // reimplement; the dropped fields are post-extraction so the cost
        // delta is negligible.
        auto Maybe_full = calc_Qdt_schur<Policy>(m, t_Qx, ns, dt);
        if (!Maybe_full) return Maybe_full.error();
        auto& full = Maybe_full.value();
        return build<Qdtm>(get<number_of_samples>(full), get<min_P>(full), get<P>(full),
                            get<gmean_i>(full), get<gtotal_ij>(full), get<gmean_ij>(full),
                            get<gsqr_i>(full), get<gvar_i>(full));
    }

    template <class C_Qdt>
        requires(U<C_Qdt, Qdt>)
    auto get_Qn(const C_Qdt& x) {
        auto n = get<number_of_samples>(x)();
        return build<Qn>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                         build<PG_n>(get<gtotal_ij>(x)() * n),
                         build<PGG_n>(get<gtotal_sqr_ij>(x)() * (n * n * 0.5)));
    }

    template <class C_P, class C_g>
        requires(U<C_P, P> && U<C_g, g>)
    auto get_Qn(const C_P& t_P, C_g const& t_g, number_of_samples n, min_P t_minP) {
        auto N = t_P().nrows();
        auto u = Matrix<double>(1, N, 1.0);
        auto G = t_g() * u;
        auto GT = tr(G);
        auto Gmean = 0.5 * G + 0.5 * GT;

        auto Gvar = elemMult(G, GT) - elemMult(Gmean, Gmean);

        return build<Qn>(n, t_minP, t_P, build<PG_n>(elemMult(t_P(), Gmean) * n()),
                         build<PGG_n>(elemMult(t_P(), Gvar) * (n() * n() * 0.5)));
    }

    // Seed Qn for the Taylor scaling-and-squaring path with the *correct*
    // small-step time-integrals — replacement for the placeholder
    // endpoint-average get_Qn(P_sub, g, …) which produced wrong gmean*/gvar*
    // (see project_qdt_taylor_seed_bug). Computes
    //
    //   A(dt) = ∫₀^dt P(t) · diag(g) · P(dt − t) dt
    //   B(dt) = ∫₀^dt ∫₀^t1 P(s) · diag(g) · P(t1 − s) · diag(g) · P(dt − t1) ds dt1
    //
    // via Taylor series in Q_sub = Q · dt_sub. The recurrences
    //   M₀ = G,   M_m = Q_sub · M_{m−1} + G · Q_sub^m
    //   N₀ = G², N_p = Q_sub · N_{p−1} + G · M_p
    // give A/dt_sub = Σ M_m / (m+1)!  and  B/dt_sub² = Σ N_p / (p+2)!.
    // Then PG_n = sub_ns · A/dt_sub and PGG_n = sub_ns² · B/dt_sub², which
    // composes correctly under sum_Qn (the same composition rule as for the
    // eig-derived seed, with PG_n and PGG_n carrying the per-time scaling).
    template <class Policy = StabilizerPolicyEnabled, class C_Q_sub, class C_P, class C_g>
        requires(U<C_P, P> && U<C_g, g>)
    auto get_Qn_via_taylor_integrals(const C_Q_sub& Q_sub, const C_P& t_P,
                                      C_g const& t_g, number_of_samples sub_ns,
                                      min_P t_minP, std::size_t order = 6ul) {
        auto N_sz = t_P().nrows();
        auto u = Matrix<double>(1, N_sz, 1.0);
        auto G_row = t_g() * u;  // [i,j] = g[i]

        Matrix<double> I_NxN(N_sz, N_sz, 0.0);
        for (std::size_t i = 0; i < N_sz; ++i) I_NxN(i, i) = 1.0;

        // M₀ = diag(g) as a full N×N matrix: row-scale identity by g[i].
        auto M_prev = elemMult(G_row, I_NxN);

        // PG accumulator = Σ M_m / (m+1)!  (m=0 term is M₀ / 1!).
        auto PG_acc = M_prev;

        using M_t = std::decay_t<decltype(M_prev)>;
        std::vector<M_t> M_store;
        M_store.reserve(order);
        M_store.push_back(M_prev);

        // Q_sub^1 = Q_sub.
        auto Q_pow = Q_sub;

        double inv_fact = 1.0;
        for (std::size_t m = 1; m < order; ++m) {
            auto M_cur = Q_sub * M_prev + elemMult(G_row, Q_pow);
            inv_fact /= static_cast<double>(m + 1);
            PG_acc = PG_acc + M_cur * inv_fact;
            M_store.push_back(M_cur);
            M_prev = std::move(M_cur);
            if (m + 1 < order)
                Q_pow = Q_sub * Q_pow;  // advance to Q_sub^(m+1)
        }

        auto PG_n_val = PG_acc * static_cast<double>(sub_ns());

        // PGG accumulator = Σ N_p / (p+2)!.  N₀ = G·G = diag(g²) materialized
        // via row-scaling M₀ by g.
        auto N_cur = elemMult(G_row, M_store[0]);

        double inv_fact2 = 0.5;  // 1 / 2! for p=0
        auto PGG_acc = N_cur * inv_fact2;

        for (std::size_t p = 1; p + 1 < order; ++p) {
            N_cur = Q_sub * N_cur + elemMult(G_row, M_store[p]);
            inv_fact2 /= static_cast<double>(p + 2);
            PGG_acc = PGG_acc + N_cur * inv_fact2;
        }

        auto PGG_n_val =
            PGG_acc * (static_cast<double>(sub_ns()) * static_cast<double>(sub_ns()));

        return build<Qn>(sub_ns, t_minP, t_P,
                         build<PG_n>(std::move(PG_n_val)),
                         build<PGG_n>(std::move(PGG_n_val)));
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Qn>
        requires(U<C_Qn, Qn>)
    static Maybe_error<C_Qn> sum_Qn(C_Qn&& one, const C_Qn& two) {
        auto n1 = get<number_of_samples>(two)();
        get<PGG_n>(one)() = (get<PGG_n>(one)() * get<P>(two)()) +
                            (get<PG_n>(one)() * get<PG_n>(two)()) +
                            (get<P>(one)() * get<PGG_n>(two)());
        get<PG_n>(one)() = (get<PG_n>(one)() * get<P>(two)()) + (get<P>(one)() * get<PG_n>(two)());
        auto Maybe_P = to_Transition_Probability<Policy>(get<P>(one)() * get<P>(two)());
        if (!Maybe_P)
            return Maybe_P.error();
        get<P>(one) = std::move(Maybe_P.value());
        get<number_of_samples>(one)() = get<number_of_samples>(one)() + n1;
        return one;
    }

    template <class Policy = StabilizerPolicyEnabled, class C_Qn, class C_Qdt>
        requires(U<C_Qn, Qn> && U<C_Qdt, Qdt>)
    static Maybe_error<C_Qn> sum_Qdt(C_Qn&& one, const C_Qdt& two) {
        auto n1 = get<number_of_samples>(two)();
        get<PGG_n>(one)() = (get<PGG_n>(one)() * get<P>(two)()) +
                            (get<PG_n>(one)() * get<gtotal_ij>(two)()) * n1 +
                            (get<P>(one)() * get<gtotal_sqr_ij>(two)()) * (0.5 * n1 * n1);
        get<PG_n>(one)() =
            (get<PG_n>(one)() * get<P>(two)()) + (get<P>(one)() * get<gtotal_ij>(two)()) * n1;

        auto Maybe_P = to_Transition_Probability<Policy>(get<P>(one)() * get<P>(two)());
        if (!Maybe_P)
            return Maybe_P.error();
        get<P>(one) = std::move(Maybe_P.value());

        get<number_of_samples>(one)() = get<number_of_samples>(one)() + n1;
        return one;
    }

    static bool is_Binomial_Approximation_valid(double N, double p, double q, double Np_min) {
        if (N * p < Np_min)
            return false;
        else if (N * q < Np_min)
            return false;
        else
            return true;
    }

    static y_mean max_possible_value_of_ymean(N_Ch_mean_value t_N, const g& t_g,
                                              Current_Baseline b) {
        return y_mean(t_N() * var::max(t_g()) + b());
    }

    static y_mean min_possible_value_of_ymean(N_Ch_mean_value t_N, g t_g, Current_Baseline b) {
        return y_mean(t_N() * var::min(t_g()) + b());
    }

    static bool crude_gmean_violation(gmean_i const& v_gm, const g& v_g) {
        auto max_g_m = var::max(v_gm());
        auto min_g_m = var::min(v_gm());
        auto max_g = var::max(v_g());
        auto min_g = var::min(v_g());
        // Absolute tolerance: numerical noise in the eigendecomposition-driven
        // gmean_ij computation grows with matrix size. For the lifted micro
        // path (M×M matrices with M up to ~100) overshoots of ~1e-13 are
        // routine. Allow a tolerance proportional to the dynamic range.
        double tol = 1e-9 * std::max(std::abs(max_g - min_g), 1.0);
        if (!((max_g_m <= max_g + tol) && (min_g_m >= min_g - tol))) {
            std::cerr << "max_g_m=" << max_g_m << " max_g=" << max_g << " min_g_m=" << min_g_m
                      << " min_g=" << min_g << "\n";
            return true;
        } else
            return false;
    }

    // ---------------------------------------------------------------------
    // In-range canaries for the conductance-derived quantities.
    //
    // With Bayesian shrinkage on gmean_ij and the row-stochastic invariant
    // on P, the values are in their physical ranges by construction. These
    // canaries verify the construction held — they don't clamp.
    //
    //   check_*_in_range    returns Maybe_error<bool>: true=clean,
    //                       false=warn-band excursion (stderr-logged),
    //                       error_message=hard out-of-range.
    //
    //   require_*_in_range  thin wrapper, returns Maybe_error<void>: hard
    //                       failure on error, silent on warn-band.
    //
    // Physical ranges:
    //   gmean_ij  ∈ [gmin, gmax]              (conditional mean conductance)
    //   gmean_i   ∈ [gmin, gmax]              (marginal mean conductance)
    //   gtotal_ij ∈ [gmin·P_ij, gmax·P_ij]    (per-entry, P-aware)
    //
    // Replaces the old force_*_in_range clamping helpers.
    // ---------------------------------------------------------------------
    template <class C_gmean_ij, class C_g>
        requires(U<C_gmean_ij, gmean_ij> && U<C_g, g>)
    static Maybe_error<bool> check_gmean_ij_in_range(const C_gmean_ij& q, const C_g& v_g) {
        const double gmax = primitive(var::max(v_g()));
        const double gmin = primitive(var::min(v_g()));
        const double range = std::max(std::abs(gmax - gmin), 1.0);
        constexpr double warn_band_rel  = 1e-9;
        constexpr double error_band_rel = 1e-3;
        const double warn_band  = range * warn_band_rel;
        const double error_band = range * error_band_rel;

        auto const& m = primitive(q());
        double max_exc = 0.0;
        std::size_t worst_i = 0, worst_j = 0;
        double worst_val = 0.0;
        bool worst_over_max = false;
        for (std::size_t i = 0; i < m.nrows(); ++i) {
            for (std::size_t j = 0; j < m.ncols(); ++j) {
                const double v = m(i, j);
                const double over  = v - gmax;
                const double under = gmin - v;
                const double exc = std::max(0.0, std::max(over, under));
                if (exc > max_exc) {
                    max_exc = exc;
                    worst_i = i;
                    worst_j = j;
                    worst_val = v;
                    worst_over_max = (over >= under);
                }
            }
        }
        if (max_exc > error_band) {
            std::ostringstream ss;
            ss << "check_gmean_ij_in_range: [" << worst_i << "," << worst_j << "] = "
               << std::scientific << std::setprecision(3) << worst_val
               << (worst_over_max ? " > gmax = " : " < gmin = ")
               << (worst_over_max ? gmax : gmin)
               << " (excursion " << max_exc << ", band " << error_band << ")";
            return error_message(ss.str());
        }
        if (max_exc > warn_band) {
            std::cerr << "[warn] check_gmean_ij_in_range: [" << worst_i << "," << worst_j << "] = "
                      << std::scientific << std::setprecision(3) << worst_val
                      << (worst_over_max ? " > gmax = " : " < gmin = ")
                      << (worst_over_max ? gmax : gmin) << "\n";
            return Maybe_error<bool>(false);
        }
        return Maybe_error<bool>(true);
    }

    template <class C_gmean_ij, class C_g>
        requires(U<C_gmean_ij, gmean_ij> && U<C_g, g>)
    static Maybe_error<void> require_gmean_ij_in_range(const C_gmean_ij& q, const C_g& v_g) {
        auto r = check_gmean_ij_in_range(q, v_g);
        if (!r) return r.error();
        return Maybe_error<void>{};
    }

    template <class C_gmean_i, class C_g>
        requires(U<C_gmean_i, gmean_i> && U<C_g, g>)
    static Maybe_error<bool> check_gmean_i_in_range(const C_gmean_i& q, const C_g& v_g) {
        const double gmax = primitive(var::max(v_g()));
        const double gmin = primitive(var::min(v_g()));
        const double range = std::max(std::abs(gmax - gmin), 1.0);
        constexpr double warn_band_rel  = 1e-9;
        constexpr double error_band_rel = 1e-3;
        const double warn_band  = range * warn_band_rel;
        const double error_band = range * error_band_rel;

        auto const& v = primitive(q());
        double max_exc = 0.0;
        std::size_t worst_i = 0;
        double worst_val = 0.0;
        bool worst_over_max = false;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const double x = v[i];
            const double over  = x - gmax;
            const double under = gmin - x;
            const double exc = std::max(0.0, std::max(over, under));
            if (exc > max_exc) {
                max_exc = exc;
                worst_i = i;
                worst_val = x;
                worst_over_max = (over >= under);
            }
        }
        if (max_exc > error_band) {
            std::ostringstream ss;
            ss << "check_gmean_i_in_range: [" << worst_i << "] = "
               << std::scientific << std::setprecision(3) << worst_val
               << (worst_over_max ? " > gmax = " : " < gmin = ")
               << (worst_over_max ? gmax : gmin)
               << " (excursion " << max_exc << ", band " << error_band << ")";
            return error_message(ss.str());
        }
        if (max_exc > warn_band) {
            std::cerr << "[warn] check_gmean_i_in_range: [" << worst_i << "] = "
                      << std::scientific << std::setprecision(3) << worst_val
                      << (worst_over_max ? " > gmax = " : " < gmin = ")
                      << (worst_over_max ? gmax : gmin) << "\n";
            return Maybe_error<bool>(false);
        }
        return Maybe_error<bool>(true);
    }

    template <class C_gmean_i, class C_g>
        requires(U<C_gmean_i, gmean_i> && U<C_g, g>)
    static Maybe_error<void> require_gmean_i_in_range(const C_gmean_i& q, const C_g& v_g) {
        auto r = check_gmean_i_in_range(q, v_g);
        if (!r) return r.error();
        return Maybe_error<void>{};
    }

    template <class C_gtotal_ij, class C_P, class C_g>
        requires(U<C_gtotal_ij, gtotal_ij> && U<C_P, P> && U<C_g, g>)
    static Maybe_error<bool> check_gtotal_ij_in_range(const C_gtotal_ij& q,
                                                     const C_P& v_P,
                                                     const C_g& v_g) {
        const double gmax = primitive(var::max(v_g()));
        const double gmin = primitive(var::min(v_g()));
        const double range = std::max(std::abs(gmax - gmin), 1.0);
        constexpr double warn_band_rel  = 1e-9;
        constexpr double error_band_rel = 1e-3;

        auto const& m = primitive(q());
        auto const& p = primitive(v_P());
        double max_exc = 0.0;
        std::size_t worst_i = 0, worst_j = 0;
        double worst_val = 0.0;
        bool worst_over_max = false;
        for (std::size_t i = 0; i < m.nrows(); ++i) {
            for (std::size_t j = 0; j < m.ncols(); ++j) {
                const double v = m(i, j);
                const double Pij = p(i, j);
                const double upper = gmax * Pij;
                const double lower = gmin * Pij;
                const double over  = v - upper;
                const double under = lower - v;
                const double exc = std::max(0.0, std::max(over, under));
                if (exc > max_exc) {
                    max_exc = exc;
                    worst_i = i;
                    worst_j = j;
                    worst_val = v;
                    worst_over_max = (over >= under);
                }
            }
        }
        if (max_exc > range * error_band_rel) {
            std::ostringstream ss;
            ss << "check_gtotal_ij_in_range: [" << worst_i << "," << worst_j << "] = "
               << std::scientific << std::setprecision(3) << worst_val
               << (worst_over_max ? " > gmax·P_ij " : " < gmin·P_ij ")
               << " (excursion " << max_exc << ")";
            return error_message(ss.str());
        }
        if (max_exc > range * warn_band_rel) {
            std::cerr << "[warn] check_gtotal_ij_in_range: [" << worst_i << "," << worst_j << "] = "
                      << std::scientific << std::setprecision(3) << worst_val
                      << (worst_over_max ? " > gmax·P_ij" : " < gmin·P_ij")
                      << " (excursion " << max_exc << ")\n";
            return Maybe_error<bool>(false);
        }
        return Maybe_error<bool>(true);
    }

    template <class C_gtotal_ij, class C_P, class C_g>
        requires(U<C_gtotal_ij, gtotal_ij> && U<C_P, P> && U<C_g, g>)
    static Maybe_error<void> require_gtotal_ij_in_range(const C_gtotal_ij& q,
                                                       const C_P& v_P,
                                                       const C_g& v_g) {
        auto r = check_gtotal_ij_in_range(q, v_P, v_g);
        if (!r) return r.error();
        return Maybe_error<void>{};
    }

    // Per-entry bound for the back-converted second moment:
    //     gsqr_min · P_ij  ≤  gtotal_sqr_ij[i,j]  ≤  gsqr_max · P_ij,
    // where gsqr_min / gsqr_max are the **min / max of g[k]² over states k**
    // (not gmin² / gmax² — the conductance vector may span zero, e.g. after
    // baseline subtraction, in which case gmin² is the *upper* end of g²
    // and the natural lower bound on E[g²|i→j] is 0).
    // Same two-tier scheme as check_gtotal_ij_in_range.
    template <class C_gtotal_sqr_ij, class C_P, class C_g>
        requires(U<C_gtotal_sqr_ij, gtotal_sqr_ij> && U<C_P, P> && U<C_g, g>)
    static Maybe_error<bool> check_gtotal_sqr_ij_in_range(const C_gtotal_sqr_ij& q,
                                                          const C_P& v_P,
                                                          const C_g& v_g) {
        auto const& g_prim = primitive(v_g());
        double gsqr_max = 0.0;
        double gsqr_min = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < g_prim.size(); ++k) {
            const double g2 = g_prim[k] * g_prim[k];
            if (g2 > gsqr_max) gsqr_max = g2;
            if (g2 < gsqr_min) gsqr_min = g2;
        }
        if (!std::isfinite(gsqr_min)) gsqr_min = 0.0;
        const double range = std::max(std::abs(gsqr_max - gsqr_min), 1.0);
        constexpr double warn_band_rel  = 1e-9;
        constexpr double error_band_rel = 1e-3;

        auto const& m = primitive(q());
        auto const& p = primitive(v_P());
        double max_exc = 0.0;
        std::size_t worst_i = 0, worst_j = 0;
        double worst_val = 0.0;
        bool worst_over_max = false;
        for (std::size_t i = 0; i < m.nrows(); ++i) {
            for (std::size_t j = 0; j < m.ncols(); ++j) {
                const double v = m(i, j);
                const double Pij = p(i, j);
                const double upper = gsqr_max * Pij;
                const double lower = gsqr_min * Pij;
                const double over  = v - upper;
                const double under = lower - v;
                const double exc = std::max(0.0, std::max(over, under));
                if (exc > max_exc) {
                    max_exc = exc;
                    worst_i = i;
                    worst_j = j;
                    worst_val = v;
                    worst_over_max = (over >= under);
                }
            }
        }
        if (max_exc > range * error_band_rel) {
            std::ostringstream ss;
            ss << "check_gtotal_sqr_ij_in_range: [" << worst_i << "," << worst_j << "] = "
               << std::scientific << std::setprecision(3) << worst_val
               << (worst_over_max ? " > gsqr_max·P_ij " : " < gsqr_min·P_ij ")
               << " (excursion " << max_exc << ")";
            return error_message(ss.str());
        }
        if (max_exc > range * warn_band_rel) {
            std::cerr << "[warn] check_gtotal_sqr_ij_in_range: [" << worst_i << "," << worst_j
                      << "] = " << std::scientific << std::setprecision(3) << worst_val
                      << (worst_over_max ? " > gsqr_max·P_ij" : " < gsqr_min·P_ij")
                      << " (excursion " << max_exc << ")\n";
            return Maybe_error<bool>(false);
        }
        return Maybe_error<bool>(true);
    }

    template <class C_gtotal_sqr_ij, class C_P, class C_g>
        requires(U<C_gtotal_sqr_ij, gtotal_sqr_ij> && U<C_P, P> && U<C_g, g>)
    static Maybe_error<void> require_gtotal_sqr_ij_in_range(const C_gtotal_sqr_ij& q,
                                                            const C_P& v_P,
                                                            const C_g& v_g) {
        auto r = check_gtotal_sqr_ij_in_range(q, v_P, v_g);
        if (!r) return r.error();
        return Maybe_error<void>{};
    }

    template <class C_Qn, class C_Patch_Model>
        requires(U<C_Qn, Qn>)
    static auto Qn_to_Qdt(const C_Qn& x, const C_Patch_Model& m)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        auto& t_P = get<P>(x);
        auto& t_g = get<g>(m);
        auto n = get<number_of_samples>(x)();
        const std::size_t N = t_P().nrows();
        Matrix<double> u(N, 1, 1.0);

        // Taylor / uniformization path: no eigendecomposition, no κ(V).
        // Pseudo-count = 10·√N·ε_mach (random-walk FP accumulation across the
        // N×N matmul chain). See
        // theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md.
        const double min_P_prior = 10.0 * std::sqrt(static_cast<double>(N)) *
                                   std::numeric_limits<double>::epsilon();

        auto r_gtotal_ij_raw     = build<gtotal_ij>    (get<PG_n>(x)() * (1.0 / n));
        auto r_gtotal_sqr_ij_raw = build<gtotal_sqr_ij>(get<PGG_n>(x)() * (2.0 / (n * n)));

        // v_g-based priors (Poisson endpoint, robust against ill-conditioned V).
        auto r_gmean_ij_prior = gmean_ij_prior(t_g);
        auto r_gsqr_ij_prior  = gsqr_ij_prior (t_g);

        auto r_gmean_ij = calc_g_ij_bayes(r_gtotal_ij_raw,     t_P, r_gmean_ij_prior, min_P_prior);
        auto r_gsqr_ij  = calc_g_ij_bayes(r_gtotal_sqr_ij_raw, t_P, r_gsqr_ij_prior,  min_P_prior);

        auto r_gvar_ij = build<gvar_ij>(r_gsqr_ij() - elemMult(r_gmean_ij(), r_gmean_ij()));

        auto r_gtotal_ij     = build<gtotal_ij>    (elemMult(r_gmean_ij(), t_P()));
        auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(elemMult(r_gsqr_ij(),  t_P()));
        auto r_gtotal_var_ij = build<gtotal_var_ij>(elemMult(r_gvar_ij(),  t_P()));

        // Source-level canaries.
        if (auto chk = require_gtotal_ij_in_range(r_gtotal_ij, t_P, t_g); !chk)
            return chk.error();
        if (auto chk = require_gtotal_sqr_ij_in_range(r_gtotal_sqr_ij, t_P, t_g); !chk)
            return chk.error();

        auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
        auto r_gsqr_i  = build<gsqr_i> (r_gtotal_sqr_ij() * u);
        // Full LTV: Var(Ā|i) = gsqr_i − gmean_i².
        auto r_gvar_i  = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        return build<Qdt>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                          std::move(r_gmean_i), std::move(r_gtotal_ij), std::move(r_gmean_ij),
                          std::move(r_gtotal_sqr_ij), std::move(r_gsqr_i), std::move(r_gvar_i),
                          std::move(r_gtotal_var_ij), std::move(r_gvar_ij));
    }

    template <class C_Qn, class C_Patch_Model>
        requires(U<C_Qn, Qn>)
    static auto Qn_to_Qdtm(const C_Qn& x, const C_Patch_Model& m)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
        auto& t_P = get<P>(x);
        auto& t_g = get<g>(m);
        auto n = get<number_of_samples>(x)();
        if (n <= 0)
            std::cerr << " nana here";
        const std::size_t N = t_P().nrows();
        Matrix<double> u(N, 1, 1.0);

        // Taylor / uniformization path: no eigendecomposition, no κ(V).
        const double min_P_prior = 10.0 * std::sqrt(static_cast<double>(N)) *
                                   std::numeric_limits<double>::epsilon();

        auto r_gtotal_ij_raw     = build<gtotal_ij>    (get<PG_n>(x)() * (1.0 / n));
        auto r_gtotal_sqr_ij_raw = build<gtotal_sqr_ij>(get<PGG_n>(x)() * (2.0 / (n * n)));

        auto r_gmean_ij_prior = gmean_ij_prior(t_g);
        auto r_gsqr_ij_prior  = gsqr_ij_prior (t_g);

        auto r_gmean_ij = calc_g_ij_bayes(r_gtotal_ij_raw,     t_P, r_gmean_ij_prior, min_P_prior);
        auto r_gsqr_ij  = calc_g_ij_bayes(r_gtotal_sqr_ij_raw, t_P, r_gsqr_ij_prior,  min_P_prior);

        // gvar_ij not stored in Qdtm; gtotal_sqr_ij materialized transiently
        // for the canary check and the gsqr_i row-sum.
        auto r_gtotal_ij     = build<gtotal_ij>    (elemMult(r_gmean_ij(), t_P()));
        auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(elemMult(r_gsqr_ij(),  t_P()));

        // Source-level canaries.
        if (auto chk = require_gtotal_ij_in_range(r_gtotal_ij, t_P, t_g); !chk)
            return chk.error();
        if (auto chk = require_gtotal_sqr_ij_in_range(r_gtotal_sqr_ij, t_P, t_g); !chk)
            return chk.error();

        auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
        auto r_gsqr_i  = build<gsqr_i> (r_gtotal_sqr_ij() * u);
        auto r_gvar_i  = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

        return build<Qdtm>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                           std::move(r_gmean_i), std::move(r_gtotal_ij), std::move(r_gmean_ij),
                           std::move(r_gsqr_i), std::move(r_gvar_i));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtg_agonist_step(FunctionTable& f, const C_Patch_Model& m,
                                const Agonist_step& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> {
        auto dt = get<number_of_samples>(t_step)() / fs;
        auto t_Qeig = f.fstop(Calc_eigen{}, m, get<Agonist_concentration>(t_step));

        if (t_Qeig) {
            auto Maybe_Qdt =
                calc_Qdtg_eig<Policy>(f, m, t_Qeig.value(), get<number_of_samples>(t_step), dt);
            if (Maybe_Qdt)

                return Maybe_Qdt.value();
        }
        auto t_Eigs = build<Qx>(calc_Qx(m, get<Agonist_concentration>(t_step)));
        return calc_Qdtg_taylor<Policy>(m, t_Eigs, get<number_of_samples>(t_step), dt);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtm_agonist_step(FunctionTable& f, const C_Patch_Model& m,
                                const Agonist_step& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
        auto dt = get<number_of_samples>(t_step)() / fs;
        auto t_Qeig = [&f, &t_step, &m]() {
            if constexpr (var::has_it_defined<Calc_eigen, decltype(f)>())
                return f.fstop(Calc_eigen{}, m, get<Agonist_concentration>(t_step));
            else
                return Macro_DMR{}.calc_eigen(m, get<Agonist_concentration>(t_step));
        }();

        if (t_Qeig) {
            auto r_Qeig = std::move(t_Qeig.value());
            auto Maybe_Qdt =
                calc_Qdtm_eig<Policy>(f, m, r_Qeig, get<number_of_samples>(t_step), dt);
            if (Maybe_Qdt)

                return Maybe_Qdt.value();
        }
        auto t_Eigs = build<Qx>(calc_Qx(m, get<Agonist_concentration>(t_step)));
        return calc_Qdtm_taylor(m, t_Eigs, get<number_of_samples>(t_step), dt);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
        requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdt_agonist_step(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step,
                               double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        auto dt = get<number_of_samples>(t_step)() / fs;

        auto t_Qeig = [this, &f, &m, &t_step]() {
            if constexpr (var::has_it_defined<Calc_eigen, decltype(f)>())
                return f.fstop(Calc_eigen{}, m, get<Agonist_concentration>(t_step));
            else
                return calc_eigen(m, get<Agonist_concentration>(t_step));
        }();

        if (t_Qeig) {
            auto Maybe_Qdt =
                calc_Qdt_eig<Policy>(f, m, t_Qeig.value(), get<number_of_samples>(t_step), dt);
            if (Maybe_Qdt)
                return Maybe_Qdt;
        }
        auto t_Qx = build<Qx>(calc_Qx(m, get<Agonist_concentration>(t_step)));
        return calc_Qdt_taylor<Policy>(m, t_Qx, get<number_of_samples>(t_step), dt);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdt(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        if constexpr (std::is_same_v<Policy, StabilizerPolicyEnabled>) {
            if constexpr (std::is_same_v<Nothing, decltype(f[Calc_Qdt_step{}])>)
                return calc_Qdt_agonist_step<Policy>(f, m, t_step, fs);
            else
                return f.f(Calc_Qdt_step{}, m, t_step, fs);
        } else {
            return calc_Qdt_agonist_step<Policy>(f, m, t_step, fs);
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtm(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
        if constexpr (std::is_same_v<Policy, StabilizerPolicyEnabled>) {
            if constexpr (std::is_same_v<Nothing, decltype(f[Calc_Qdtg_step{}])>)
                return calc_Qdtm_agonist_step<Policy>(f, m, t_step, fs);
            else
                return f.f(Calc_Qdtm_step{}, m, t_step, fs);
        } else {
            return calc_Qdtm_agonist_step<Policy>(f, m, t_step, fs);
        }
    }

    template <class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtg(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> {
        if constexpr (std::is_same_v<Nothing, decltype(f[Calc_Qdtg_step{}])>)
            return calc_Qdtg_agonist_step(f, m, t_step, fs);
        else
            return f.f(Calc_Qdtg_step{}, m, t_step, fs);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qn_bisection(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step,
                           double fs, int order) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qn>> {
        auto dt = get<number_of_samples>(t_step)() / fs;
        auto ns = get<number_of_samples>(t_step);
        auto t_Eigs = f.fstop(Calc_eigen{}, m, get<Agonist_concentration>(t_step));

        if (!t_Eigs)
            return t_Eigs.error();
        else {
            double scale = std::pow(2.0, -1.0 * order);

            number_of_samples n_ss(ns() * scale);
            double sdt = dt * scale;
            auto t_Psub = calc_P<Policy>(m, t_Eigs.value(), sdt, get<min_P>(m)() * scale);
            auto r_Qn = get_Qn(t_Psub, get<g>(m), n_ss, min_P(get<min_P>(m)() * scale));
            for (std::size_t i = 0; i < order; ++i) {
                r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
            }
            assert(get<number_of_samples>(r_Qn)() == ns());
            return r_Qn;
        }
    }

    template <class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdt_bisection(FunctionTable& f, const C_Patch_Model& m, const Agonist_step& t_step,
                            double fs, std::size_t order)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        auto maybe_Qn = calc_Qn_bisection(f, m, t_step, fs, order);
        if (!maybe_Qn)
            return maybe_Qn.error();
        else {
            return Qn_to_Qdt(maybe_Qn.value());
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_Qdt(FunctionTable& f, const C_Patch_Model& m, const std::vector<Agonist_step>& t_step,
                  double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        if (t_step.empty())
            return error_message("Emtpy agonist step");
        else {
            auto v_Qdt0 = calc_Qdt<Policy>(f, m, t_step[0], fs);
            if (!v_Qdt0)
                return v_Qdt0.error();
            else {
                auto v_Qrun = get_Qn(v_Qdt0.value());
                for (std::size_t i = 1; i < t_step.size(); ++i) {
                    auto v_Qdti = calc_Qdt<Policy>(f, m, t_step[i], fs);
                    if (!v_Qdti)
                        return v_Qdti.error();
                    else {
                        auto Maybe_Qrun = sum_Qdt(std::move(v_Qrun), v_Qdti.value());
                        if (!Maybe_Qrun)
                            return Maybe_Qrun.error();
                        v_Qrun = std::move(Maybe_Qrun.value());
                    }
                }
                return Qn_to_Qdt(v_Qrun, m);
            }
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_Qdtg(FunctionTable& f, const C_Patch_Model& m,
                   const std::vector<Agonist_step>& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> {
        if (t_step.empty())
            return error_message("Emtpy agonist step");
        if (t_step.size() == 1)
            return calc_Qdtg(f, m, t_step[0], fs);

        auto v_Qdt0 = calc_Qdtg(f, m, t_step[0], fs);
        if (!v_Qdt0)
            return v_Qdt0.error();
        auto v_Prun = get<P_half>(v_Qdt0.value());
        auto v_ns = get<number_of_samples>(v_Qdt0.value())();
        for (std::size_t i = 1; i < t_step.size(); ++i) {
            auto v_Qdti = calc_Qdtg(f, m, t_step[i], fs);
            if (!v_Qdti)
                return v_Qdti.error();
            else {
                auto Maybe_v_Prun =
                    to_Transition_Probability<Policy>(v_Prun() * get<P_half>(v_Qdti.value())());
                if (!Maybe_v_Prun)
                    return Maybe_v_Prun.error();

                v_Prun() = std::move(Maybe_v_Prun.value()());
                v_ns = v_ns + get<number_of_samples>(v_Qdti.value())();
            }
        }
        return build<Qdtg>(number_of_samples(v_ns), get<min_P>(m), std::move(v_Prun), get<g>(m));
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_Qdtm(FunctionTable& f, const C_Patch_Model& m,
                   const std::vector<Agonist_step>& t_step, double fs)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
        if (t_step.empty())
            return error_message("Emtpy agonist step");
        if (t_step.size() == 1)
            return calc_Qdtm<Policy>(f, m, t_step[0], fs);

        std::size_t i0 = 0;
        while ((get<number_of_samples>(t_step[i0]) == 0) && (i0 < t_step.size())) ++i0;
        auto v_Qdt0 = calc_Qdt<Policy>(f, m, t_step[i0], fs);
        if (!v_Qdt0)
            return v_Qdt0.error();
        auto v_Qrun = get_Qn(v_Qdt0.value());
        for (std::size_t i = i0 + 1; i < t_step.size(); ++i) {
            auto v_Qdti = calc_Qdt<Policy>(f, m, t_step[i], fs);
            if (!v_Qdti)
                return v_Qdti.error();
            else {
                auto Maybe_v_Qrun = sum_Qdt(std::move(v_Qrun), v_Qdti.value());
                if (!Maybe_v_Qrun)
                    return Maybe_v_Qrun.error();

                v_Qrun = std::move(Maybe_v_Qrun.value());
            }
        }
        return Qn_to_Qdtm(v_Qrun, m);
    }

    template <class FunctionTable, class C_Patch_Model>
    // requires(U<C_Patch_Model, Patch_Model> )
    auto calc_Qdt_bisection(FunctionTable& f, const C_Patch_Model& m,
                            const std::vector<Agonist_step>& t_step, double fs, std::size_t order)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        if (t_step.empty())
            return error_message("Emtpy agonist step");
        else {
            auto v_Qn0 = calc_Qn_bisection(f, m, t_step[0], fs, order);
            if (!v_Qn0)
                return v_Qn0.error();
            else {
                auto v_Qrun = v_Qn0.value();
                for (std::size_t i = 1; i < t_step.size(); ++i) {
                    auto v_Qni = calc_Qn_bisection(f, m, t_step[i], fs, order);
                    if (!v_Qni)
                        return v_Qni.error();
                    else
                        v_Qrun = sum_Qn(std::move(v_Qrun), v_Qni.value());
                }
                return Qn_to_Qdt(v_Qrun);
            }
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdt(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                  double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        return calc_Qdt<Policy>(f, m, t_step(), fs);
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtm(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                   double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
        return calc_Qdtm<Policy>(f, m, t_step(), fs);
    }

 template <class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdtg(FunctionTable& f, const C_Patch_Model& m, const Agonist_evolution& t_step,
                   double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtg>> {
        return calc_Qdtg(f, m, t_step(), fs);
    }


    template <class FunctionTable, class C_Patch_Model>
    //   requires(U<C_Patch_Model, Patch_Model>)
    auto calc_Qdt_bisection(FunctionTable& f, const C_Patch_Model& m,
                            const Agonist_evolution& t_step, double fs, std::size_t order)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
        return calc_Qdt_bisection(f, m, t_step(), fs, order);
    }

    Maybe_error<bool> test_conductance_mean(const Matrix<double> gmean, const Matrix<double> g) {
        auto max_g = var::max(g);
        auto min_g = var::min(g);
        Maybe_error<bool> out = true;
        return reduce(
            [max_g, min_g](Maybe_error<bool> succeeds, auto e) -> Maybe_error<bool> {
                if (e < min_g) {
                    return error_message("conductance too negative");
                } else if (e > max_g) {
                    return error_message("conductance too positve");
                } else {
                    return succeeds;
                }
            },
            out, gmean);
    }

    Maybe_error<bool> test_conductance_variance(const Matrix<double> gvar, const Matrix<double> g) {
        auto max_abs_g = maxAbs(g) * 0.5;
        Maybe_error<bool> out = true;
        return reduce(
            [max_abs_g](Maybe_error<bool> succeeds, auto e) -> Maybe_error<bool> {
                if (e < 0) {
                    return error_message("variance  negative");
                } else if (e > max_abs_g) {
                    return error_message("variance too big");
                } else {
                    return succeeds;
                }
            },
            out, gvar);
    }

    template <class C_Matrix>
        requires U<C_Matrix, Matrix<double>>
    Maybe_error<bool> test_conductance_variance(const C_Matrix& var,
                                                Conductance_variance_error_tolerance tol) {
        if (var.ncols() == var.nrows()) {
            for (std::size_t i = 0; i < var.nrows(); ++i)
                for (std::size_t j = 0; j < var.ncols(); ++j)
                    if (primitive(var(i, j)) + tol() < 0) {
                        std::stringstream ss;
                        ss << " negative diagonal variance at i=" << i << ", j= " << j << "\n"
                           << var;
                        return error_message(ss.str());
                    }
            return true;
        } else
            for (std::size_t i = 0; i < var.size(); ++i)
                if (primitive(var[i]) + tol() < 0) {
                    std::stringstream ss;
                    ss << " negative  variance at i=" << i << "\n" << var;
                    return error_message(ss.str());
                }
        return true;
    }

    template <class... Variances, class C_Vector_Space>
    Maybe_error<bool> test_conductance_variances(const C_Vector_Space& q,
                                                 Conductance_variance_error_tolerance tol) {
        return ((type_name<Variances>() >> test_conductance_variance(get<Variances>(q)(), tol)) &&
                ...);
    }

    template <class C_Qdt>
        requires U<C_Qdt, Qdt>
    Maybe_error<bool> test(const C_Qdt& q, Conductance_variance_error_tolerance tol) {
        return "fails Qdt test\n" >>
               test_conductance_variances<gmean_i, gtotal_ij, gmean_ij, gtotal_sqr_ij, gsqr_i,
                                          gvar_i, gtotal_var_ij, gvar_ij>(q, tol);
    }

    template <class C_Matrix>
        requires U<std::decay_t<C_Matrix>, Matrix<double>>
    C_Matrix truncate_negative_variance(C_Matrix var) {
        return apply(
            [](auto const& value) {
                using var::max;
                using std::max;
                return max(value * 0.0, value);
            },
            var);
    }

    /*
            template<uses_recursive_aproximation
         recursive,uses_averaging_aproximation averaging,
       uses_variance_aproximation variance> auto run_old(const Patch_State
       &t_prior, Qdt const &t_Qdt, Patch_Model const &m, const Experiment_step
     &p, double fs) const { auto &p_y = get<Patch_current>(p); auto &p_P_mean =
       get<P_mean>(t_prior); auto &p_P_Cov = get<P_Cov>(t_prior);
       
            double e =
                get<Current_Noise>(m).value() * fs/
         get<number_of_samples>(p).value(); double N = get<N_Ch_mean>(m)(); auto
         N_states = p_P_mean().nrows(); Matrix<double> u(N_states, 1, 1.0);
         
            auto SmD = p_P_Cov() - diag(p_P_mean());
            double gSg = xtAx(get<gmean_i>(t_Qdt)(), SmD) +
                         getvalue(p_P_mean() *
                                  zip([](auto x, auto y) { return x * y; },
                                      get<gtotal_ij>(t_Qdt)(),
         get<gmean_ij>(t_Qdt)())
           * u);
           
            double ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
            
            auto e_mu = e + N * ms;
            auto v_y_mean = y_mean(N * getvalue(p_P_mean() *
         get<gmean_i>(t_Qdt)())); auto v_y_var = y_var(e_mu + N * gSg); if
         (std::isnan(p_y.value())) { auto v_vplogL = vplogL(0.0); auto v_plogL =
         plogL(std::numeric_limits<double>::quiet_NaN()); auto v_eplogL =
         eplogL(std::numeric_limits<double>::quiet_NaN()); auto v_P_cov =
         P_Cov(AT_B_A(get<P>(t_Qdt)(), SmD)); auto v_P_mean = P_mean(p_P_mean()
     * get<P>(t_Qdt)()); v_P_cov() = v_P_cov() + diag(v_P_mean());
     
              return Patch_State(logL(get<logL>(t_prior)()),v_P_mean, v_P_cov,
           v_y_mean, v_y_var, v_plogL, v_eplogL, v_vplogL);
              // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
              // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
              //      auto test = mp_state_information::test(P_mean, P__cov,
              //      tolerance_); if (test.has_value())
              //        return
     Op(mp_state_information::adjust(std::move(P_mean),
              // std::move(P__cov),
              //                                               y_mean, y_var,
       plogL,
              //                                               eplogL,
              // vplogL,Q_dt.min_P(), e));
              //      else
              //        return Op(false, "fails at intertrace prediction!!: " +
              //        test.error());
            }
            auto dy = p_y.value() - v_y_mean();
            auto chi = dy / v_y_var();
            auto v_P_cov = P_Cov(AT_B_A(get<P>(t_Qdt)(), SmD));
            auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)());
            v_P_cov() = v_P_cov() + diag(v_P_mean());
            
            auto chi2 = dy * chi;
            
            auto v_plogL = plogL(0);
            if (v_y_var() > 0)
              v_plogL() = -0.5 * log(2 * std::numbers::pi * v_y_var()) - 0.5 *
       chi2; else v_plogL() = std::numeric_limits<double>::infinity();
       
            auto v_eplogL = eplogL(-0.5 * log(2 * std::numbers::pi * v_y_var())
     - 0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)" vplogL v_vplogL(0.5);
            // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
            
            //    auto test = mp_state_information::test(P_mean, P__cov, y_mean,
           y_var,
            //    plogL,
            //                                           eplogL, e,
     tolerance());
            //    if (!test) {
            //      std::stringstream ss;
            
            //      ss << "\nP_mean \n" << P_mean;
            //      ss << "\nPcov \n" << P__cov;
            //      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;
            
            //      return Op(false, "\nfails in trace!!!; error=" +
       test.error()()
         +
            //      ss.str());
            //    } else
            return Patch_State(logL(get<logL>(t_prior)()+v_plogL()),v_P_mean,
         v_P_cov, v_y_mean, v_y_var, v_plogL, v_eplogL, v_vplogL);
            }
            
            */

    //    Maybe_error<Patch_State> DVR(const Patch_State &t_prior, Qdt const
    //    &t_Qdt,
    //                                 Patch_Model const &m, const Agonist_step &p,
    //                                 const Patch_current &p_y, double fs)
    //                                 const
    //                                 {
    //        //  auto &p_y = get<Patch_current>(p);
    //        auto &p_P_mean = get<P_mean>(t_prior);
    //        auto &p_P_Cov = get<P_Cov>(t_prior);

    //        double e =
    //            get<Current_Noise>(m).value() * fs/
    //            get<number_of_samples>(p).value();
    //        double N = get<N_Ch_mean>(m)();

    //        auto N_states = p_P_mean().ncols();
    //        Matrix<double> u(N_states, 1, 1.0);

    //        auto SmD = p_P_Cov() - diag(p_P_mean());

    //        if (std::isnan(p_y.value())) {
    //            auto v_P_cov = P_Cov(AT_B_A(get<P>(t_Qdt)(), SmD));
    //            auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)());
    //            v_P_cov() = v_P_cov() + diag(v_P_mean());

    //            return Patch_State(
    //                logL(get<logL>(t_prior)()), elogL(get<elogL>(t_prior)()),
    //                vlogL(get<vlogL>(t_prior)()), v_P_mean, v_P_cov,
    //                y_mean(NaN), y_var(NaN), plogL(NaN), eplogL(NaN),
    //                vplogL(NaN));
    //            // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
    //            // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
    //            //      auto test = mp_state_information::test(P_mean, P__cov,
    //            //      tolerance_); if (test.has_value())
    //            //        return
    //            Op(mp_state_information::adjust(std::move(P_mean),
    //            // std::move(P__cov),
    //            //                                               y_mean,
    //            y_var, plogL,
    //            //                                               eplogL,
    //            // vplogL,Q_dt.min_P(), e));
    //            //      else
    //            //        return Op(false, "fails at intertrace prediction!!:
    //            "
    //            +
    //            //        test.error());
    //        }
    //        double gSg = xtAx(get<gmean_i>(t_Qdt)(), SmD) +
    //                     getvalue(p_P_mean() *
    //                              zip([](auto x, auto y) { return x * y; },
    //                                  get<gtotal_ij>(t_Qdt)(),
    //                                  get<gmean_ij>(t_Qdt)()) *
    //                              u);

    //        double sSg =
    //            xtAy(get<gvar_i>(t_Qdt)(), SmD, get<gmean_i>(t_Qdt)()) +
    //            getvalue(p_P_mean() *
    //                     zip([](auto x, auto y) { return x * y; },
    //                         get<gtotal_var_ij>(t_Qdt)(),
    //                         get<gmean_ij>(t_Qdt)()) *
    //                     u);

    //        double sSs =
    //            xtAx(get<gvar_i>(t_Qdt)(), SmD) +
    //            getvalue(p_P_mean() *
    //                     zip([](auto x, auto y) { return x * y; },
    //                         get<gtotal_var_ij>(t_Qdt)(),
    //                         get<gvar_ij>(t_Qdt)())
    //                         *
    //                     u);

    //        auto sS = tr(get<gvar_i>(t_Qdt)()) * SmD * get<P>(t_Qdt)() +
    //                  p_P_mean() * get<gtotal_var_ij>(t_Qdt)();

    //        auto gS = tr(get<gmean_i>(t_Qdt)()) * SmD * get<P>(t_Qdt)() +
    //                  p_P_mean() * get<gtotal_ij>(t_Qdt)();

    //        double ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

    //        double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    //        double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    //        auto e_mu = e + N * ms0;

    //        auto v_y_mean = y_mean(N * getvalue(p_P_mean() *
    //        get<gmean_i>(t_Qdt)()) -
    //                               N * 0.5 / e_mu * sSg);

    //        auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    //        auto v_y_var = y_var(std::max(e_mu + N * gSg - N * zeta *
    //        sqr(sSg), e)); auto dy = p_y.value() - v_y_mean();

    //        auto chi = dy / v_y_var();

    //        auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)() + chi * gS -
    //                               (chi * zeta * sSg + 0.5 / e_mu) * sS);

    //        auto v_P_cov = build<P_Cov>(
    //            AT_B_A(get<P>(t_Qdt)(), SmD) + diagpos(v_P_mean()) -
    //            (zeta + N / v_y_var() * sqr(zeta * sSg)) * XTX(sS) +
    //            (2.0 * N / v_y_var() * zeta * sSg) * X_plus_XT(tr(sS) * gS) -
    //            (N / v_y_var()) * XTX(gS));

    //        auto chi2 = dy * chi;

    //        auto v_plogL = plogL(0.0);
    //        if (v_y_var() > 0)
    //            v_plogL =
    //                plogL(-0.5 * log(2 * std::numbers::pi * v_y_var()) - 0.5 *
    //                chi2);
    //        else
    //            v_plogL = plogL(std::numeric_limits<double>::infinity());

    //        auto v_eplogL = eplogL(-0.5 * log(2 * std::numbers::pi *
    //        v_y_var())
    //        -
    //                               0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    //        vplogL v_vplogL(0.5);
    //        // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
    //        std::cerr << get<number_of_samples>(p).value() << "\t" << v_P_mean
    //        << "\n"; return Patch_State(logL(get<logL>(t_prior)() +
    //        v_plogL()),
    //                           elogL(get<elogL>(t_prior)() + v_eplogL()),
    //                           vlogL(get<vlogL>(t_prior)() + v_vplogL()),
    //                           v_P_mean, v_P_cov, v_y_mean, v_y_var, v_plogL,
    //                           v_eplogL, v_vplogL);
    //    }

    auto safely_calculate_error(auto const& e, auto const& N, auto const& gSg,
                                auto const& ms) const {
        if (std::isfinite(primitive(gSg)) && primitive(gSg) > 0) {
            if (isfinite(primitive(ms)) && primitive(ms) > 0) {
                return build<y_var>(e + N * gSg + N * ms);
            } else {
                return build<y_var>(e + N * gSg);
            }
        } else {
            if (isfinite(primitive(ms)) && primitive(ms) > 0)
                return build<y_var>(e + N * ms);
            else
                return build<y_var>(e);
        }
    }

// =============================================================================
// DEAD CODE — kept for reference only.
// safely_calculate_y_mean_yvar_Pmean_PCov (and its sister overload below) is
// the previous deepseek-style closed-form Taylor block (sSg / sSs / e_mu /
// zeta). Both overloads are unreachable: only their own internal recursive
// fallbacks call them, no external caller exists. Their helper
// safely_calculate_Algo_Pmean_Pcov has 12 call sites in this dead family but
// NO definition anywhere in the codebase — confirming the family is never
// instantiated. The current IRT/MRT path lives in
// safely_calculate_Algo_State_recursive (line ~4143). Kept inside #if 0 to
// preserve the earlier algebra for reference, while removing the noise from
// grep / IDE navigation.
// =============================================================================
#if 0
    template <class recursive, class averaging, class variance, class variance_correction,
              class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 (U<C_Patch_State, Algo_State> || U<C_Patch_State, Algo_State_Dynamic>))
    Maybe_error<C_Patch_State> safely_calculate_y_mean_yvar_Pmean_PCov(
        C_Patch_State const& t_prior, C_Qdt const& t_Qdt, C_Patch_Model const& m, C_double const& N,
        const Patch_current& p_y, double fs) const {
        constexpr bool PoissonDif = true;
        using Transf = transformation_type_t<C_Qdt>;

        auto& p_P_mean = get<P_mean>(t_prior());
        auto SmD = get<P_Cov>(t_prior())() - diag(p_P_mean());
        auto& y = p_y.value();
        bool is_y_nan = std::isnan(y);
        auto y_baseline = get<Current_Baseline>(m);
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value();
        Matrix<double> u(p_P_mean().size(), 1, 1.0);
        if constexpr (averaging::value == 0) {
            auto& t_g = get<g>(m);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_g()) + y_baseline());
            auto gSg = getvalue(TranspMult(t_g(), SmD) * t_g());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);

            auto r_y_var = build<y_var>(e);
            if (primitive(gSg) > 0)
                r_y_var() = r_y_var() + N * gSg;

            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        } else if constexpr (averaging::value == 2) {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
            auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
            if constexpr (variance_correction::value) {
                auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
                auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                auto& t_gvar_i = get<gvar_i>(t_Qdt);
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

                auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
                auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

                auto delta_emu = sqr(ms + e / N) - 2.0 / N * sSs;
                if (!std::isfinite(primitive(delta_emu)) || primitive(delta_emu) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

                auto e_mu = e + N * ms0;
                auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) -
                                              N * 0.5 / e_mu * sSg + y_baseline());
                auto zeta = N / (2 * sqr(e_mu) + N * sSs);

                if (!std::isfinite(primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg))) ||
                    primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg)) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto r_y_var = build<y_var>(e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
                auto& t_P = get<P>(t_Qdt);
                auto dy = y - r_y_mean();
                auto chi = dy / r_y_var();
                auto chi2 = dy * chi;

                auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
                auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();

                auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS -
                                                     (chi * zeta * sSg + 0.5 / e_mu) * sS);

                if (!Maybe_r_P_mean)
                    return Maybe_r_P_mean.error();
                auto r_P_mean = std::move(Maybe_r_P_mean.value());

                auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD) + diag(r_P_mean * t_P()) -
                                            (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
                                            (2.0 * N / r_y_var() * zeta * sSg) *
                                                X_plus_XT(TranspMult(sS, gS)) -
                                            (N / r_y_var()) * XTX(gS));

                if (!all_Probability_elements(primitive(r_P_mean)) ||
                    !all_Covariance_elements(primitive(r_P_cov())))
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                else {
                    auto r_macro_algo = macror_algorithm(
                        ToString(MacroR2<recursive, averaging, variance, variance_correction>{}));
                    return std::tuple(std::move(r_y_mean), std::move(r_y_var), std::move(r_P_mean),
                                      std::move(r_P_cov), std::move(chi2), std::move(r_macro_algo));
                }

            } else {
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
                if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, uses_averaging_aproximation<1>, variance, variance_correction>(
                        t_prior, t_Qdt, m, N, p_y, fs);
                auto r_y_mean =
                    build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
                if constexpr (PoissonDif)
                    e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
                else
                    e = e + get<Proportional_Noise>(m).value() * abs(y);
                auto r_y_var = build<y_var>(e + N * gSg);
                if constexpr (variance::value) {
                    auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                    if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                        r_y_var() = r_y_var() + N * ms;
                    } else {
                        return safely_calculate_Algo_Pmean_Pcov<recursive, averaging,
                                                                uses_variance_aproximation<false>,
                                                                variance_correction>(
                            is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                            p_y, fs, SmD);
                    }
                }
                return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                        variance_correction>(
                    is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y,
                    fs, SmD);
            }
        } else /* if constexpr (averaging::value == 1) */ {
            auto& t_gmean_i = [&t_Qdt](){ if constexpr(averaging::value>0) {return get<gmean_i>(t_Qdt);} else {return get<g>(t_Qdt);} }();    
            auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i());
            if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                return safely_calculate_y_mean_yvar_Pmean_PCov<
                    recursive, uses_averaging_aproximation<0>, uses_variance_aproximation<false>,
                    uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N, p_y,
                                                                         fs);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);
            auto r_y_var = build<y_var>(e + N * gSg);
            if constexpr (variance::value) {
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                if (std::isfinite(primitive(ms)) && primitive(ms) > 0) {
                    r_y_var() = r_y_var() + N * ms;
                } else {
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, uses_variance_aproximation<false>,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                }
            }
            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        }
    }

    template <class recursive, class averaging, class variance, class variance_correction,
              class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction>)
    Maybe_error<std::tuple<Transfer_Op_to<C_Qdt, y_mean>, Transfer_Op_to<C_Qdt, y_var>,
                           Transfer_Op_to<C_Qdt, P_mean>, Transfer_Op_to<C_Qdt, P_Cov>,
                           Transfer_Op_to<C_Qdt, double>, macror_algorithm>>
        safely_calculate_y_mean_yvar_Pmean_PCov(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                                C_Patch_Model const& m, C_double const& N,
                                                const Patch_current& p_y, double fs) const {
        constexpr bool PoissonDif = true;
        using Transf = transformation_type_t<C_Qdt>;

        auto& p_P_mean = get<P_mean>(t_prior);
        auto SmD = get<P_Cov>(t_prior)() - diag(p_P_mean());
        auto& y = p_y.value();
        bool is_y_nan = std::isnan(y);
        auto y_baseline = get<Current_Baseline>(m);
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value();
        Matrix<double> u(p_P_mean().size(), 1, 1.0);
        if constexpr (averaging::value == 0) {
            auto& t_g = get<g>(m);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_g()) + y_baseline());
            auto gSg = getvalue(TranspMult(t_g(), SmD) * t_g());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);

            auto r_y_var = build<y_var>(e);
            if (primitive(gSg) > 0)
                r_y_var() = r_y_var() + N * gSg;

            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        } else if constexpr (averaging::value == 2) {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
            auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
            if constexpr (variance_correction::value) {
                auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
                auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                auto& t_gvar_i = get<gvar_i>(t_Qdt);
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

                auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
                auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

                auto delta_emu = sqr(ms + e / N) - 2.0 / N * sSs;
                if (!std::isfinite(primitive(delta_emu)) || primitive(delta_emu) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

                auto e_mu = e + N * ms0;
                auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) -
                                              N * 0.5 / e_mu * sSg + y_baseline());
                auto zeta = N / (2 * sqr(e_mu) + N * sSs);

                if (!std::isfinite(primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg))) ||
                    primitive(N * ms0 + N * gSg - N * zeta * sqr(sSg)) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N,
                                                                             p_y, fs);

                auto r_y_var = build<y_var>(e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
                auto& t_P = get<P>(t_Qdt);
                auto dy = y - r_y_mean();
                auto chi = dy / r_y_var();
                auto chi2 = dy * chi;

                auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
                auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();

                auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS -
                                                     (chi * zeta * sSg + 0.5 / e_mu) * sS);

                if (!Maybe_r_P_mean)
                    return Maybe_r_P_mean.error();
                auto r_P_mean = std::move(Maybe_r_P_mean.value());

                auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD) + diag(r_P_mean * t_P()) -
                                            (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
                                            (2.0 * N / r_y_var() * zeta * sSg) *
                                                X_plus_XT(TranspMult(sS, gS)) -
                                            (N / r_y_var()) * XTX(gS));

                if (!all_Probability_elements(primitive(r_P_mean)) ||
                    !all_Covariance_elements(primitive(r_P_cov())))
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, variance,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                else {
                    auto r_macro_algo = macror_algorithm(
                        ToString(MacroR2<recursive, averaging, variance, variance_correction>{}));
                    return std::tuple(std::move(r_y_mean), std::move(r_y_var), std::move(r_P_mean),
                                      std::move(r_P_cov), std::move(chi2), std::move(r_macro_algo));
                }

            } else {
                auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                           getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
                if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                    return safely_calculate_y_mean_yvar_Pmean_PCov<
                        recursive, uses_averaging_aproximation<1>, variance, variance_correction>(
                        t_prior, t_Qdt, m, N, p_y, fs);
                auto r_y_mean =
                    build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
                if constexpr (PoissonDif)
                    e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
                else
                    e = e + get<Proportional_Noise>(m).value() * abs(y);
                auto r_y_var = build<y_var>(e + N * gSg);
                if constexpr (variance::value) {
                    auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                    if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                        r_y_var() = r_y_var() + N * ms;
                    } else {
                        return safely_calculate_Algo_Pmean_Pcov<recursive, averaging,
                                                                uses_variance_aproximation<false>,
                                                                variance_correction>(
                            is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                            p_y, fs, SmD);
                    }
                }
                return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                        variance_correction>(
                    is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y,
                    fs, SmD);
            }
        } else /* if constexpr (averaging::value == 1) */ {
            auto& t_gmean_i = get<gmean_i>(t_Qdt);
            auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i());
            if (!std::isfinite(primitive(gSg)) || primitive(gSg) <= 0)
                return safely_calculate_y_mean_yvar_Pmean_PCov<
                    recursive, uses_averaging_aproximation<0>, uses_variance_aproximation<false>,
                    uses_taylor_variance_correction_aproximation<false>>(t_prior, t_Qdt, m, N, p_y,
                                                                         fs);
            auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
            if constexpr (PoissonDif)
                e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
            else
                e = e + get<Proportional_Noise>(m).value() * abs(y);
            auto r_y_var = build<y_var>(e + N * gSg);
            if constexpr (variance::value) {
                auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
                if (std::isfinite(primitive(ms)) && primitive(ms) > 0) {
                    r_y_var() = r_y_var() + N * ms;
                } else {
                    return safely_calculate_Algo_Pmean_Pcov<
                        recursive, averaging, uses_variance_aproximation<false>,
                        uses_taylor_variance_correction_aproximation<false>>(
                        is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N,
                        p_y, fs, SmD);
                }
            }
            return safely_calculate_Algo_Pmean_Pcov<recursive, averaging, variance,
                                                    variance_correction>(
                is_y_nan, std::move(r_y_mean), std::move(r_y_var), t_prior, t_Qdt, m, N, p_y, fs,
                SmD);
        }
    }
#endif // dead safely_calculate_y_mean_yvar_Pmean_PCov family
    // =========================================================================

    template <class C_y_var>
    auto calculate_logL(bool y_is_nan, C_y_var const& r_y_var, auto const& chi2,
                         auto& m) const -> Transfer_Op_to<C_y_var, logL> {
        using DX = var::dx_of_dfdx_t<C_y_var>;
        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<logL>(std::move(base));
        }
        if (get<Proportional_Noise>(m).value() == 0) {
            return build<logL>(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
                               0.5 * chi2());
        } else {
            return build<logL>(0.5 * chi2() - log(var::Poisson_noise_normalization(
                                                  primitive(r_y_var()),
                                                  primitive(get<Proportional_Noise>(m).value()))));
        }
    }

	    template <class C_y_var>
	    auto calculate_elogL(bool y_is_nan, C_y_var const& r_y_var, auto& m) const {
	        using DX = var::dx_of_dfdx_t<C_y_var>;
	        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<elogL>(std::move(base));
        }
	        if (get<Proportional_Noise>(m).value() == 0) {
	            return build<elogL>(-0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5);
	        } else {
	            const auto pn_value = get<Proportional_Noise>(m).value();
	            if constexpr (var::is_derivative_v<std::decay_t<decltype(pn_value)>>) {
	                return build<elogL>(var::Poisson_noise_expected_logL(r_y_var() ,
	                                                                    pn_value));
	            } else {
	                auto pn = var::init_with_dx<DX>(pn_value, dx);
	                return build<elogL>(
	                    var::Poisson_noise_expected_logL(r_y_var() , pn));
	            }
	        }
	    }

    // ────────────────────────────────────────────────────────────────────
    // Smooth (C∞) lower-bound of two values.
    //
    //     softmin(a, b; ε) = (a + b − √((a − b)² + ε²)) / 2
    //
    // Properties:
    //   - softmin(a, b; 0) = min(a, b)                   (recovers hard min)
    //   - softmin(a, b; ε) ≤ min(a, b)                   (always ≤ true min)
    //   - smooth in (a, b): no kink, no Heaviside in ∂/∂(a, b)
    //   - softmin(1, x; ε) = x          when  x ≪ 1 − ε
    //                     ≈ 1 − ε / 2  when  x ≈ 1
    //                     = 1          when  x ≫ 1 + ε
    //
    // The trust coefficient uses softmin instead of `min(1, factor·alfa_p)` to
    // eliminate the kink at `factor·alfa_p = 1`. Each kink in α(θ) introduces a
    // step in ∂α/∂θ that randomizes across realizations and pumps variance into
    // the score; softmin removes the step. ε controls the smoothing band — too
    // small recovers the kink, too large biases α below `min(1, factor·alfa_p)`
    // by ε / 2 even far from the boundary. ε = 1e-4 is a reasonable default;
    // the residual bias 5e-5 is tiny vs. the `(1 − factor) = 0.1` margin.
    // ────────────────────────────────────────────────────────────────────
    template <class A, class B>
    auto softmin(A const& a, B const& b, double eps) const {
        using std::sqrt;
        auto diff = a - b;
        auto sum  = a + b;
        auto root = sqrt(diff * diff + eps * eps);
        return (sum - root) * 0.5;
    }


    // ────────────────────────────────────────────────────────────────────
    // LogSumExp soft-min of the simplex trust bounds, capped at 1 — a C∞
    // alternative to calculate_trust_coefficient's hard min + d_i sign-branch.
    //
    // NOT CURRENTLY WIRED IN. The active trust path is calculate_trust_coefficient
    // (softmin-fold over candidates); this is kept as a reference implementation
    // of the LogSumExp formulation we may switch to.
    //
    //   α = −(1/k)·log( exp(−k·1) + Σ_i exp(−k·factor·bound_i) )
    //
    // with the per-component binding bound (the largest α keeping μ_i = p_i + α·d_i
    // inside [0,1]):
    //   d_i > 0 :  bound_i = (1 − p_i)/d_i        (upper face)
    //   d_i < 0 :  bound_i = (0 − p_i)/d_i        (lower face; > 0 since d_i < 0)
    //   d_i = 0 :  no constraint                  (skipped)
    //
    // Why this is smooth where the hard min is not: as d_i → 0 the forward bound
    // → +∞ from BOTH sides, so its term exp(−k·factor·bound_i) → 0 — the component
    // drops out of the sum smoothly, the d_i=0 face-switch is invisible (both faces
    // give +∞ → zero weight), and no 1/d_i blow-up survives into α. The max-shift m
    // is taken over PRIMITIVES only, so it cancels in the ratio and derivatives pass
    // through exactly while the exponentials stay bounded. k = LSE sharpness (→ hard
    // min as k→∞; larger k is sharper but has more curvature near the knee).
    // ────────────────────────────────────────────────────────────────────
    template <class C_Matrix1, class C_Matrix2>
    auto log_Sum_Exp_min_1(C_Matrix1 const& t_pmean, C_Matrix2 const& d,
                           double factor) const {
        using std::exp;
        using std::log;
        constexpr double k = 1e4;  // LSE sharpness

        auto d_p = var::inside_out(d);          // per-element scalars (Derivative-aware)
        auto r_pmean = var::inside_out(t_pmean);

        // Pass 1 (primitives only): the smallest exponent argument, for the shift.
        // x_cap = 1 ; x_i = factor·bound_i ; m = min over the cap and valid bounds.
        double m = 1.0;
        for (std::size_t i = 0; i < d_p.size(); ++i) {
            const double di = primitive(d_p[i]);
            if (di == 0.0)
                continue;
            const double pi = primitive(r_pmean[i]);
            const double bound = (di > 0.0) ? (1.0 - pi) / di : (0.0 - pi) / di;
            m = std::min(m, factor * bound);
        }

        // Pass 2 (Derivative-aware): Σ exp(−k·(x − m)), starting with the cap term.
        auto sum = exp(-k * (1.0 - m)) + 0.0 * d_p[0];  // cap term + dx carrier (zero deriv)
        for (std::size_t i = 0; i < d_p.size(); ++i) {
            const double di = primitive(d_p[i]);
            if (di == 0.0)
                continue;
            auto bound_i = (di > 0.0) ? (1.0 - r_pmean[i]) / d_p[i]   // upper face, derivatives intact
                                      : (0.0 - r_pmean[i]) / d_p[i];  // lower face
            sum = sum + exp(-k * (factor * bound_i - m));
        }
        auto alpha = m - (1.0 / k) * log(sum);
        return build<trust_coefficient>(alpha);
    }


    // ────────────────────────────────────────────────────────────────────
    // Trust coefficient for the recursive POSTERIOR-MEAN update.
    //
    // The mean update has the form  μ_new = μ + α · d  where d = chi · gS.
    // For μ_new to remain a probability vector (μ_new_i ∈ [0, 1]) we need
    //     d_i > 0 :  α  ≤  (1 − μ_i) / d_i
    //     d_i < 0 :  α  ≤  −μ_i      / d_i
    // and α ∈ (0, 1]. The function returns the largest α satisfying every
    // per-index bound, scaled by `factor` (∈ (0, 1)) for safety when α < 1.
    //
    // This bound only constrains the mean. The recursive COVARIANCE update
    //     Σ_new = Σ_pre  −  (α · N / y_var) · gSᵀ gS
    // has its own (often tighter) constraint — see calculate_psd_trust_coefficient
    // and docs/math/trust_coefficient.md.
    // ────────────────────────────────────────────────────────────────────
    // Derivative-aware, and SMOOTH IN THE INDEX SELECTION. The constraint is
    // α ≤ min_i alfa_i with alfa_i = (1−p_i)/d_i (d_i>0) or −p_i/d_i (d_i<0).
    //
    // The previous version took the min via a HARD argmin over primitives and
    // evaluated the binding expression at that single index. That keeps the
    // VALUE continuous but makes ∂α/∂θ DISCONTINUOUS: when the argmin switches
    // between components (a tie, reachable under any θ-perturbation) the
    // derivative jumps from ∂alfa_a/∂θ to ∂alfa_b/∂θ. The per-sample detailed
    // dump localized exactly this jump in ∂(trust_coefficient)/∂θ as the source
    // of the numerical-Fisher (FD-of-AD) instability.
    //
    // Fix: fold a softmin over ALL candidates instead — α = softmin(1,
    // factor·alfa_0, factor·alfa_1, …), derivative-aware on every component. Then
    // ∂α/∂θ = Σ_i w_i·∂alfa_i/∂θ with smooth softmax weights → continuous across
    // ties. The softmin also smooths the 1-vs-alfa cap (no kink there either).
    //
    // For a 2-state simplex this is degenerate (p_C+p_O=1 and the tangent
    // displacement d_C=−d_O make the two candidates identical in value AND
    // derivative, so the old argmin tie was already harmless) — the fix matters
    // for the PSD constraint and for any model with >2 states. Kept here for
    // uniformity. ε = trust_softmin_eps (1e-4); the fold adds ~ε bias per
    // candidate, tiny vs the (1 − factor) = 0.1 margin. No candidates → α = 1
    // with zero derivative (same as the old no-binding branch).
    template <class C_Matrix1, class C_Matrix2>
    auto calculate_trust_coefficient(C_Matrix1 const& t_pmean, C_Matrix2 const& d,
                                     double factor) const {
        constexpr double trust_softmin_eps = 1e-4;
        auto d_p = var::inside_out(d);
        auto r_pmean= var::inside_out(t_pmean);

        auto alpha = 1.0 + 0.0 * d_p[0];  // cap at 1, with (zero) derivatives carried
        // α_μ is now DECOUPLED (applied only to the mean update, never to Σ). At a
        // likelihood-residual zero-crossing d=chi·gS→0 the bound (1−p)/d_i is at its
        // 1/d_i singularity, but the mean step α_μ·chi·gS≡0 there for any α, so the
        // residual irregularity is harmless — it no longer reaches the covariance.
        // The d=0 sign-branch is kept as-is (the singular region is exactly the
        // non-binding region, disjoint from where α_μ actually constrains).
        for (std::size_t i = 0; i < d_p.size(); ++i) {
            const auto&  d_i = d_p[i];
            if (primitive(d_i) > 0) {
                auto alfa_i = (1.0 - r_pmean[i]) / d_p[i];  // bound (1−p_i)/d_i, derivatives intact
                alpha = softmin(alpha, alfa_i * factor, trust_softmin_eps);
            } else if (primitive(d_i) < 0) {
                auto alfa_i = (0.0 - r_pmean[i]) / d_p[i];  // bound −p_i/d_i, derivatives intact
                alpha = softmin(alpha, alfa_i * factor, trust_softmin_eps);
            }
        }
        return build<trust_coefficient>(alpha);
    }

    // ────────────────────────────────────────────────────────────────────
    // Trust coefficient for the recursive POSTERIOR-COVARIANCE update.
    //
    // The Σ down-date is
    //     Σ_new = Σ_pre  −  (α · β) · gSᵀ gS,    β = N / y_var
    // with the rank-1 contribution diagonal at index i equal to (α·β)·gS_i².
    // For each diagonal of Σ_new to stay non-negative (a NECESSARY condition
    // for PSD) we need
    //     Σ_pre_(i,i)  −  (α · β) · gS_i²  ≥  0
    // ⇒   α  ≤  Σ_pre_(i,i)  /  (β · gS_i²)        whenever gS_i ≠ 0.
    //
    // This is a per-index bound; the function returns the smallest such α
    // across i, scaled by `factor` for safety when < 1. Returns 0 if any
    // diagonal of Σ_pre is non-positive at an index where gS_i ≠ 0 — in that
    // case the information update would push an already-degenerate diagonal
    // negative, and the caller should refuse the recursive step.
    //
    // Combined use:  α = min(α_μ, α_Σ) — the more restrictive bound applies
    // to BOTH updates so they are scaled consistently. See
    // docs/math/trust_coefficient.md for the full derivation, the off-diagonal
    // PSD argument, and notes on when a Joseph-form update would be a
    // stronger alternative.
    // ────────────────────────────────────────────────────────────────────
    // Derivative-aware, and SMOOTH in the diagonal selection (same fix as
    // calculate_trust_coefficient): instead of a hard argmin over diagonals it
    // folds a softmin over all valid per-diagonal bounds, so ∂α/∂θ stays
    // continuous when the binding diagonal switches. beta = N / y_var carries
    // derivatives in AD context. See the in-body comment for why this overload
    // is the one that drove the numerical-Fisher instability.
    //
    // The trust coefficient's role is to REDUCE gS so the rank-1 down-date
    // preserves diagonal positivity — never to abort the step with α = 0.
    // Indices where Σ_pre(i,i) ≤ 0 are SKIPPED: those diagonals are already
    // degenerate and no α > 0 fixes them, so the constraint at i contributes no
    // information to α. Other indices set the bound.
    template <class C_Sigma, class C_gS, class C_beta>
    auto calculate_psd_trust_coefficient(C_Sigma const& Sigma_pre,
                                         C_gS const& gS,
                                         C_beta const& beta,
                                         double factor) const {
        constexpr double trust_softmin_eps = 1e-4;
        auto const& Sigma_p = primitive(Sigma_pre);
        auto const& gS_p = primitive(gS);

        // SMOOTH over diagonals (see calculate_trust_coefficient): α = softmin(1,
        // factor·alfa_0, factor·alfa_1, …) with alfa_i = Σ_ii/(β·gS_i²),
        // derivative-aware on every valid diagonal — NO hard argmin. The PSD
        // bounds are INDEPENDENT across i (no simplex degeneracy linking them),
        // so the old argmin made ∂α/∂θ discontinuous whenever the binding
        // diagonal switched (Σ_00·gS_1² = Σ_11·gS_0²). That jump — localized by
        // the per-sample detailed dump — is the numerical-Fisher instability;
        // the softmin fold removes it (∂α/∂θ = Σ_i w_i·∂alfa_i/∂θ, continuous).
        auto alpha = 1.0 + 0.0 * gS[0];  // cap at 1, with (zero) derivatives carried
        for (std::size_t i = 0; i < gS_p.size(); ++i) {
            if (gS_p[i] == 0.0 || Sigma_p(i, i) <= 0.0)
                continue;  // unconstrained / already-degenerate diagonal — skip (as before)
            auto alfa_i = Sigma_pre(i, i) / (beta * gS[i] * gS[i]);  // Σ_ii/(β·gS_i²), derivatives intact
            alpha = softmin(alpha, alfa_i * factor, trust_softmin_eps);
        }
        return build<trust_coefficient>(alpha);
    }

    template <class C_Sigma, class C_Downdate>
    auto calculate_psd_trust_coefficient(C_Sigma const& Sigma_pre,
                                         C_Downdate const& downdate,
                                         double factor) const {
        constexpr double trust_softmin_eps = 1e-4;
        auto const& Sigma_p = primitive(Sigma_pre);
        auto const& down_p = primitive(downdate);

        // SMOOTH over diagonals (same fix as the other PSD overload): softmin
        // over {1} ∪ {factor·Σ_ii/downdate_ii}, derivative-aware, no hard argmin
        // — so ∂α/∂θ stays continuous when the binding diagonal switches.
        auto alpha = 1.0 + 0.0 * downdate(0, 0);  // cap at 1, derivatives (zero) carried
        for (std::size_t i = 0; i < Sigma_p.nrows(); ++i) {
            if (down_p(i, i) <= 0.0 || Sigma_p(i, i) <= 0.0)
                continue;
            auto alfa_i = Sigma_pre(i, i) / downdate(i, i);  // Σ_ii/downdate_ii, derivatives intact
            alpha = softmin(alpha, alfa_i * factor, trust_softmin_eps);
        }
        return build<trust_coefficient>(alpha);
    }

    template <bool dynamic, class averaging, class variance, class C_Patch_State, class C_Qdt,
              class C_Patch_Model, class C_double>
        requires(uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> && (U<C_Patch_State, Patch_State>))
    auto safely_calculate_Algo_State_non_recursive(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                                   C_Patch_Model const& m, C_double const& N,
                                                   const Patch_current& p_y, double fs) const
        -> Maybe_error<Transfer_Op_to<
            C_Patch_State, std::conditional_t<dynamic, Algo_State_Dynamic, Algo_State>>> {
        constexpr bool PoissonDif = true;
        auto& y = p_y.value();

        auto const& t_P = [&t_Qdt]() { 
            if constexpr(averaging::value > 0) 
               {return get<P>(t_Qdt);} 
            else {return get<P_half>(t_Qdt); }}
            ();
        auto  p_P_mean = get<P_mean>(t_prior());
        auto  p_P_Cov = get<P_Cov>(t_prior());
        if constexpr (averaging::value == 0) {
            auto maybe_p_P_mean = to_Probability(p_P_mean() * t_P());
            if (!maybe_p_P_mean.valid())
                return maybe_p_P_mean.error();
            p_P_mean() = std::move(maybe_p_P_mean.value());
            auto maybe_p_P_Cov = to_Covariance_Probability(AT_B_A(t_P(), p_P_Cov() - diag(p_P_mean())) + diag(p_P_mean() * t_P()));
            if (!maybe_p_P_Cov.valid())
                return maybe_p_P_Cov.error();
            p_P_Cov() = std::move(maybe_p_P_Cov.value());
        }


        auto SmD = get<P_Cov>(t_prior())() - diag(p_P_mean());
        auto y_baseline = get<Current_Baseline>(m);
        auto const& t_gmean_i = [&t_Qdt, &m]() -> auto const& {
            if constexpr (averaging::value > 0) {
                return get<gmean_i>(t_Qdt)();
            } else {
                return get<g>(m)();
            }
        }();
        auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i) + y_baseline());
        auto gSg = getvalue(TranspMult(t_gmean_i, get<P_Cov>(t_prior())()) * t_gmean_i);
        auto dy = y - r_y_mean();
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value() +
                 get<Proportional_Noise>(m).value() * abs(y - r_y_mean());

        auto r_y_var = build<y_var>(e);
        r_y_var() = r_y_var() + N * gSg;
        using std::sqrt;
        auto r_r_std= build<r_std>(dy/sqrt(r_y_var()));
        auto chi = dy / r_y_var();
        auto chi2 = build<Chi2>(dy * chi);
        auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
        if (!Maybe_r_P_mean.valid())
            return Maybe_r_P_mean.error();

        auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

        auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
        auto Maybe_r_P_cov = to_Covariance_Probability(r_P_cov() + diag(r_P_mean()));

        if (!Maybe_r_P_cov.valid())
            return Maybe_r_P_cov.error();
        r_P_cov() = std::move(Maybe_r_P_cov.value());

        auto alfa = trust_coefficient(1.0);
        //auto r_logL = calculate_logL(std::isnan(y), r_y_var, chi2, alfa, m);

        if constexpr (!dynamic) {
            Transfer_Op_to<C_Patch_State, Algo_State> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<r_std>(out())= std::move(r_r_std);
            get<trust_coefficient>(out()) = alfa;
            get<taylor_trust_coefficient>(out()) = taylor_trust_coefficient(1.0);
            get<taylor_vSv>(out()) = taylor_vSv(0.0);
            get<taylor_strength>(out()) = taylor_strength(0.0);

            get<Chi2>(out()) = std::move(chi2);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())() = std::move(r_P_cov());
            return out;
        } else {
            Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<r_std>(out())= std::move(r_r_std);
            get<trust_coefficient>(out()) = alfa;
            get<taylor_trust_coefficient>(out()) = taylor_trust_coefficient(1.0);
            get<taylor_vSv>(out()) = taylor_vSv(0.0);
            get<taylor_strength>(out()) = taylor_strength(0.0);

            get<Chi2>(out()) = std::move(chi2);
            get<P_mean_t2_y0>(out())() = std::move(r_P_mean());
            get<P_Cov_t2_y0>(out())() = std::move(r_P_cov());
            return out;
        }
    }

    template <bool dynamic, class averaging, class variance, class variance_correction,
              class variance_form, class C_Patch_State, class C_Qdt,
              class C_Patch_Model, class C_double>

        requires(uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 uses_variance_form_aproximation_c<variance_form> &&
                 U<C_Patch_State, Patch_State>)
    auto safely_calculate_Algo_State_recursive(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                               C_Patch_Model const& m, C_double const& N,
                                               const Patch_current& p_y, double fs) const
        -> Maybe_error<Transfer_Op_to<
            C_Patch_State, std::conditional_t<dynamic, Algo_State_Dynamic, Algo_State>>> {
        constexpr const double trust_multiplying_factor = 0.9;
        auto& y = p_y.value();
        if (std::isnan(y)) {
            return safely_calculate_Algo_State_non_recursive<
                dynamic, averaging, uses_variance_aproximation<false>>(
                t_prior, t_Qdt, m, N, p_y, fs);
        }
        auto const& t_P = [&t_Qdt](){
            if constexpr(averaging::value==0){
                return get<P_half>(t_Qdt);
            }
            else{
                return get<P>(t_Qdt);
            }
        }();   

        auto p_P_mean = get<P_mean>(t_prior());
        if constexpr (averaging::value==0) {
            auto Maybe_p_P_mean = to_Probability(p_P_mean() * t_P());
            if (!Maybe_p_P_mean.valid())
                {return Maybe_p_P_mean.error();}
            p_P_mean() = std::move(Maybe_p_P_mean.value());
        }
        
        auto p_P_Cov = get<P_Cov>(t_prior());
        if constexpr (averaging::value==0) {
            auto Maybe_p_P_Cov = to_Covariance_Probability(AT_B_A(t_P(), get<P_Cov>(t_prior())() - diag(get<P_mean>(t_prior())())) + diag(p_P_mean()));
            if (!Maybe_p_P_Cov.valid())
                {return Maybe_p_P_Cov.error();
            } 
            p_P_Cov() = std::move(Maybe_p_P_Cov.value());
        }   
        
        auto SmD = p_P_Cov()  - diag(p_P_mean());

        auto y_baseline = get<Current_Baseline>(m);
        constexpr bool PoissonDif = true;
        auto& t_gmean_i = [&t_Qdt, &m]() -> decltype(auto) {
            if constexpr (averaging::value > 0)
                return get<gmean_i>(t_Qdt);
            else
                return get<g>(t_Qdt);
        }();

        auto gSg = [&t_gmean_i,&p_P_Cov, &SmD, &p_P_mean, &t_Qdt]() {
            if constexpr (averaging::value == 2) {
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
                Matrix<double> u(p_P_mean().size(), 1, 1.0);
        
                return getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                       getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
                       
            } else {
                return getvalue(TranspMult(t_gmean_i(), p_P_Cov()) * t_gmean_i());
            }
        }();
        auto r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());
        auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
                 get<Pink_Noise>(m).value() +
                 get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
        auto r_y_var = build<y_var>(e + N * gSg);
        if constexpr (variance::value&& averaging::value>0) {
            // Compute the gvar_i flavor that matches the gSg above, on the fly,
            // so the result is independent of which Qdt(m)-flavor t_Qdt provides
            // (their `gvar_i` fields carry different mathematical objects):
            //   av=2 (I family): residual = gsqr_i − (gtotal_ij ∘ gmean_ij)·𝟏
            //                    = E_j[Var(Ā|i,j)|i]. The variance-of-conditional
            //                    -mean piece is already inside gSg's second term,
            //                    so the residual avoids double-counting.
            //   av=1 (M family): total = gsqr_i − gmean_i² = Var(Ā|X₀=i). gSg has
            //                    no boundary cross-cov term, so the full variance
            //                    per starting state is what's needed.
            // See theory/macroir/notes/gvar_i_overcount_audit.md.
            auto ms = [&]() {
                auto& t_gsqr_i = get<gsqr_i>(t_Qdt);
                if constexpr (averaging::value == 2 ||
                              variance_form::value == variance_residual) {
                    auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                    auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
                    Matrix<double> u(p_P_mean().size(), 1, 1.0);
                    auto gvar_i_residual =
                        t_gsqr_i() - elemMult(t_gtotal_ij(), t_gmean_ij()) * u;
                    return getvalue(p_P_mean() * gvar_i_residual);
                } else {
                    auto gvar_i_total =
                        t_gsqr_i() - elemMult(t_gmean_i(), t_gmean_i());
                    return getvalue(p_P_mean() * gvar_i_total);
                }
            }();
            if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                r_y_var() = r_y_var() + N * ms;
            } else {
                return error_message("invalid channel noise", ms);
            }
        }

        auto dy = y - r_y_mean();
        auto chi = dy / r_y_var();
        using std::sqrt; 
        auto r_r_std=build<r_std>( dy/sqrt(r_y_var()));
        auto chi2 = build<Chi2>(dy * chi);
        // gS is the endpoint-frame Bayes gain row vector: gS = Cov(X_end, y).
        // For all averaging values the rank-1 down-date and α_Σ check live at
        // the endpoint frame (sigma_pre), so gS must be there too.
        //   avg=0:  obs at midpoint via instantaneous g.   gS = gᵀ · Σ_mid · P_half
        //   avg=1:  obs depends on start state via gmean_i. gS = gmean_iᵀ · Σ_start · P
        //   avg=2:  obs integrated over interval.          gS = gmean_iᵀ·SmD·P + p·gtotal_ij
        // For avg=0/1 the trailing  · t_P()  propagates the start/mid-frame
        // gain through the remaining Markov dynamics so that XTX(gS) lives in
        // the same frame as sigma_pre. (Pre-fix the trailing · t_P was missing,
        // making the down-date frame-mismatched; manifested as a Distortion-
        // Induced-Bias spike in macro_R at long intervals / large Num_ch.)
        auto gS = [&t_gmean_i, &p_P_Cov,&SmD, &p_P_mean, &t_Qdt, &t_P]() {
            if constexpr (averaging::value == 2) {
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);

                return TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
            } else {
                return TranspMult(t_gmean_i(), p_P_Cov()) * t_P();
            }
        }();
        // gS is a probability displacement: μ·t_P + α·chi·gS must remain a
        // probability, which forces gS·u = 0. Surface any drift in p_P_Cov /
        // gtotal_ij upstream rather than letting it propagate silently into
        // the posterior-mean update.
        if (auto Maybe_gS_check = to_Probability_displacement(gS); !Maybe_gS_check)
            return Maybe_gS_check.error();

        // Σ_pre = post-Markov, pre-Bayesian-update covariance. Computed once and
        // reused for the PSD trust check and for the rank-1 down-date below.
        auto sigma_pre = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P());

        // ===========================================================
        // IRT (av=2) / MRT (av=1) rank-1 quasi-Laplace branch.
        // -----------------------------------------------------------
        // Per theory/macroir/docs/Macro_IRT/macroirt_supplement.tex.
        // Drops the standard Kalman update in favour of a Newton step in
        // direction (γ̃ᵀΣ + ṽᵀΣ), where v = γ̄₀ + (δ/V)·σ̄²₀ absorbs the
        // heteroscedastic open-channel noise contribution.
        //
        // Self-contained: builds and returns its own Algo_State /
        // Algo_State_Dynamic. The standard Kalman block below runs only
        // when this branch is not selected.
        //
        // TODO (dynamic outputs): the derivative-path extras (d_GS,
        // P_mean_t11_y0, P_mean_t10_y1, etc.) below use the standard
        // Kalman intermediates (chi, gS) with the IRT-derived `alfa`.
        // For accurate derivative tracking under vc=true, derive
        // IRT-specific versions of these extras.
        // ===========================================================
        if constexpr (variance_correction::value && averaging::value > 0) {
            constexpr double trust_softmin_eps_irt = 1e-4;

            if constexpr (averaging::value == 1) {
                // Exact rank-2 MacroMRT update.  The observation model is
                // start-state based, so Newton iterations live in the start
                // frame and the final posterior is propagated through P.
                constexpr std::size_t max_newton_steps_mrt = 5;
                // Exit when the applied step is below tolerance: at that point
                // the iterate is stationary, so the final in-loop K coincides
                // with the post-loop K_final (derivative included) and the
                // mean/covariance derivative coupling is restored. See
                // handoff_state.md Task 9.
                constexpr double newton_tol_mrt = 1e-10;

                auto sigma2_i =
                    get<gsqr_i>(t_Qdt)() - elemMult(t_gmean_i(), t_gmean_i());
                auto gammaS0 = TranspMult(t_gmean_i(), p_P_Cov());
                auto sigmaS0 = TranspMult(sigma2_i, p_P_Cov());
                if (auto Maybe_gammaS0_check = to_Probability_displacement(gammaS0);
                    !Maybe_gammaS0_check)
                    return Maybe_gammaS0_check.error();
                if (auto Maybe_sigmaS0_check = to_Probability_displacement(sigmaS0);
                    !Maybe_sigmaS0_check)
                    return Maybe_sigmaS0_check.error();

                auto p_iter = p_P_mean();
                auto alfa_mu_exact = build<trust_coefficient>(1.0 + 0.0 * p_iter[0]);
                auto final_delta = dy;
                auto final_V = r_y_var() - N * gSg;
                auto final_vSv = gSg;

                for (std::size_t iter = 0; iter < max_newton_steps_mrt; ++iter) {
                    auto y_mean_iter =
                        N * getvalue(p_iter * t_gmean_i()) + y_baseline();
                    auto delta_iter = y - y_mean_iter;
                    auto V_iter = e + N * getvalue(p_iter * sigma2_i);
                    if (!(std::isfinite(primitive(V_iter))) || primitive(V_iter) <= 0.0)
                        return error_message("invalid MacroMRT rank-2 V", V_iter);

                    auto v_iter = t_gmean_i() + (delta_iter / V_iter) * sigma2_i;
                    auto vS0 = TranspMult(v_iter, p_P_Cov());
                    auto vSv = getvalue(vS0 * v_iter);
                    auto b = getvalue(vS0 * sigma2_i);
                    auto c = getvalue(sigmaS0 * sigma2_i);

                    auto m11 = V_iter / N + vSv;
                    auto m12 = b;
                    auto m22 = c - 2.0 * V_iter * V_iter / N;
                    auto det = m11 * m22 - m12 * m12;
                    if (!(std::isfinite(primitive(det))) ||
                        std::abs(primitive(det)) <= 1e-30)
                        return error_message("singular MacroMRT rank-2 Woodbury matrix", det);

                    auto k11 = m22 / det;
                    auto k12 = (0.0 - m12) / det;
                    auto k22 = m11 / det;

                    auto delta_p = p_iter - p_P_mean();
                    auto delta_v = getvalue(delta_p * v_iter);
                    auto delta_s = getvalue(delta_p * sigma2_i);

                    auto q_gamma = 2.0 * delta_iter / V_iter;
                    auto q_sigma = delta_iter * delta_iter / (V_iter * V_iter) -
                                   1.0 / V_iter;
                    auto qS0 = q_gamma * gammaS0 + q_sigma * sigmaS0;
                    auto qv = getvalue(qS0 * v_iter);
                    auto qs = getvalue(qS0 * sigma2_i);

                    auto delta_term =
                        (delta_v * k11 + delta_s * k12) * vS0 +
                        (delta_v * k12 + delta_s * k22) * sigmaS0;
                    auto qS_post =
                        qS0 - (qv * k11 + qs * k12) * vS0 -
                        (qv * k12 + qs * k22) * sigmaS0;
                    auto p_candidate = p_P_mean() + delta_term + 0.5 * qS_post;
                    auto mean_step = p_candidate - p_iter;

                    alfa_mu_exact = calculate_trust_coefficient(
                        p_iter, mean_step, trust_multiplying_factor);

                    auto Maybe_next_p =
                        to_Probability(p_iter + alfa_mu_exact() * mean_step);
                    if (!Maybe_next_p)
                        return Maybe_next_p.error();
                    p_iter = std::move(Maybe_next_p.value());
                    final_delta = delta_iter;
                    final_V = V_iter;
                    final_vSv = vSv;
                    const double step_inf = var::max(apply(
                        [](auto const& v) { return std::abs(v); },
                        primitive(alfa_mu_exact() * mean_step)));
                    if (step_inf < newton_tol_mrt)
                        break;
                }

                auto y_mean_final =
                    N * getvalue(p_iter * t_gmean_i()) + y_baseline();
                final_delta = y - y_mean_final;
                final_V = e + N * getvalue(p_iter * sigma2_i);
                if (!(std::isfinite(primitive(final_V))) || primitive(final_V) <= 0.0)
                    return error_message("invalid final MacroMRT rank-2 V", final_V);

                auto v_final = t_gmean_i() + (final_delta / final_V) * sigma2_i;
                auto vS0_final = TranspMult(v_final, p_P_Cov());
                auto vSv_final = getvalue(vS0_final * v_final);
                auto b_final = getvalue(vS0_final * sigma2_i);
                auto c_final = getvalue(sigmaS0 * sigma2_i);

                auto m11_final = final_V / N + vSv_final;
                auto m12_final = b_final;
                auto m22_final = c_final - 2.0 * final_V * final_V / N;
                auto det_final =
                    m11_final * m22_final - m12_final * m12_final;
                if (!(std::isfinite(primitive(det_final))) ||
                    std::abs(primitive(det_final)) <= 1e-30)
                    return error_message("singular final MacroMRT rank-2 Woodbury matrix",
                                         det_final);

                auto k11_final = m22_final / det_final;
                auto k12_final = (0.0 - m12_final) / det_final;
                auto k22_final = m11_final / det_final;
                final_vSv = vSv_final;

                auto vS_end = vS0_final * t_P();
                auto sigmaS_end = sigmaS0 * t_P();
                auto cov_downdate_start =
                    k11_final * XTX(vS0_final) +
                    k12_final * X_plus_XT(TranspMult(vS0_final, sigmaS0)) +
                    k22_final * XTX(sigmaS0);
                auto cov_downdate_end =
                    k11_final * XTX(vS_end) +
                    k12_final * X_plus_XT(TranspMult(vS_end, sigmaS_end)) +
                    k22_final * XTX(sigmaS_end);

                auto alfa_sigma_exact = calculate_psd_trust_coefficient(
                    sigma_pre, cov_downdate_end, trust_multiplying_factor);
                // Start- and end-frame downdates have different diagonals, so
                // each needs its own trust coefficient — same pattern as the
                // IRT branch below. Re-using α_end on the start-frame downdate
                // can over-shoot when sigma_pre[i,i] > p_P_Cov[i,i], pushing
                // P_cov_t1_y1's diagonal substantially negative.
                auto alfa_sigma_start = calculate_psd_trust_coefficient(
                    p_P_Cov(), cov_downdate_start, trust_multiplying_factor);

                auto Maybe_r_P_mean = to_Probability(p_iter * t_P());
                if (!Maybe_r_P_mean)
                    return Maybe_r_P_mean.error();
                auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

                auto Maybe_r_P_cov = to_Covariance_Probability(
                    sigma_pre - alfa_sigma_exact() * cov_downdate_end);
                if (!Maybe_r_P_cov)
                    return Maybe_r_P_cov.error();
                auto r_P_cov = build<P_Cov>(std::move(Maybe_r_P_cov.value()));

                auto Maybe_P_cov_t1_y1 = to_Covariance_Probability(
                    p_P_Cov() - alfa_sigma_start() * cov_downdate_start);
                if (!Maybe_P_cov_t1_y1)
                    return Maybe_P_cov_t1_y1.error();

                if (!all_Probability_elements(primitive(r_P_mean())) ||
                    !all_Covariance_elements(primitive(r_P_cov()))) {
                    return error_message("error in P_mean or P_cov (rank-2 MRT)");
                }

                double taylor_strength_p = [&]() {
                    double sigma_sq = std::abs(primitive(
                        getvalue(TranspMult(sigma2_i, sigma2_i))));
                    double gamma_sq = std::abs(primitive(
                        getvalue(TranspMult(t_gmean_i(), t_gmean_i()))));
                    return std::abs(primitive(final_delta / final_V)) *
                           std::sqrt(sigma_sq) /
                           std::max(std::sqrt(gamma_sq), 1e-30);
                }();

                if constexpr (!dynamic) {
                    Transfer_Op_to<C_Patch_State, Algo_State> out;
                    get<y_mean>(out()) = std::move(r_y_mean);
                    get<y_var>(out()) = std::move(r_y_var);
                    get<r_std>(out()) = std::move(r_r_std);
                    get<Chi2>(out()) = std::move(chi2);
                    get<P_mean>(out())() = std::move(r_P_mean());
                    get<P_Cov>(out())() = std::move(r_P_cov());
                    get<trust_coefficient>(out()) = alfa_mu_exact;
                    get<taylor_trust_coefficient>(out()) =
                        build<taylor_trust_coefficient>(alfa_sigma_exact());
                    get<taylor_vSv>(out()) = build<taylor_vSv>(final_vSv);
                    get<taylor_strength>(out()) =
                        taylor_strength(taylor_strength_p);
                    return out;
                } else {
                    Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;
                    get<y_mean>(out()) = std::move(r_y_mean);
                    get<y_var>(out()) = std::move(r_y_var);
                    get<Chi2>(out()) = std::move(chi2);
                    get<r_std>(out()) = std::move(r_r_std);
                    get<trust_coefficient>(out()) = alfa_mu_exact;
                    get<taylor_trust_coefficient>(out()) =
                        build<taylor_trust_coefficient>(alfa_sigma_exact());
                    get<taylor_vSv>(out()) = build<taylor_vSv>(final_vSv);
                    get<taylor_strength>(out()) =
                        taylor_strength(taylor_strength_p);
                    get<P>(out()) = get<P>(t_Qdt);
                    get<gmean_i>(out()) = get<gmean_i>(t_Qdt);
                    get<gvar_i>(out()) = get<gvar_i>(t_Qdt);

                    auto Maybe_r_P_mean_t2_y0 = to_Probability(p_P_mean() * t_P());
                    if (!Maybe_r_P_mean_t2_y0)
                        return Maybe_r_P_mean_t2_y0.error();
                    auto Maybe_r_P_cov_t2_y0 = to_Covariance_Probability(sigma_pre);
                    if (!Maybe_r_P_cov_t2_y0)
                        return Maybe_r_P_cov_t2_y0.error();

                    auto r_P_mean_0t_y0 = diag(p_P_mean()) * t_P();
                    auto r_P_mean_0t_y1 = diag(p_iter) * t_P();
                    auto r_P_cross_cov_0t_y0 =
                        SmD * t_P() + diag(p_P_mean()) * t_P();
                    auto SmD1 = Maybe_P_cov_t1_y1.value() - diag(p_iter);
                    auto r_P_cross_cov_0t_y1 =
                        SmD1 * t_P() + diag(p_iter) * t_P();

                    get<P_mean_0t_y0>(out())() = std::move(r_P_mean_0t_y0);
                    get<P_mean_0t_y1>(out())() = std::move(r_P_mean_0t_y1);
                    get<P_cross_cov_0t_y0>(out())() =
                        std::move(r_P_cross_cov_0t_y0);
                    get<P_cross_cov_0t_y1>(out())() =
                        std::move(r_P_cross_cov_0t_y1);
                    get<P_Cov_t2_y0>(out())() =
                        std::move(Maybe_r_P_cov_t2_y0.value());
                    get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                    get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                    get<P_mean_t2_y0>(out())() =
                        std::move(Maybe_r_P_mean_t2_y0.value());
                    get<P_mean_t1_y1>(out())() = std::move(p_iter);
                    get<d_gS>(out())() = std::move(gS);
                    get<P_Cov_t1_y1>(out())() =
                        std::move(Maybe_P_cov_t1_y1.value());
                    return out;
                }
            } else {
                {
                    // Exact rank-2 MacroIRT update.  Newton iterations keep the
                    // observation state in the start frame, while the final
                    // posterior mean/covariance are lifted to the endpoint frame
                    // through the IR tilde contractions.
                    constexpr std::size_t max_newton_steps_irt = 5;
                    constexpr double newton_tol_irt = 1e-10;

                    auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
                    auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
                    auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
                    auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
                    Matrix<double> u_irt(p_P_mean().size(), 1, 1.0);

                    auto sigma2_i =
                        get<gsqr_i>(t_Qdt)() -
                        elemMult(t_gtotal_ij(), t_gmean_ij()) * u_irt;

                    auto gammaS0 = TranspMult(t_gmean_i(), p_P_Cov());
                    auto sigmaS0 = TranspMult(sigma2_i, p_P_Cov());
                    if (auto Maybe_gammaS0_check = to_Probability_displacement(gammaS0);
                        !Maybe_gammaS0_check)
                        return Maybe_gammaS0_check.error();
                    if (auto Maybe_sigmaS0_check = to_Probability_displacement(sigmaS0);
                        !Maybe_sigmaS0_check)
                        return Maybe_sigmaS0_check.error();

                    auto sigmaS =
                        TranspMult(sigma2_i, SmD) * t_P() +
                        p_P_mean() * t_gtotal_var_ij();
                    if (auto Maybe_sigmaS_check = to_Probability_displacement(sigmaS);
                        !Maybe_sigmaS_check)
                        return Maybe_sigmaS_check.error();

                    auto gamma_sigma0 = getvalue(gammaS0 * sigma2_i);
                    auto sigma_sigma0 = getvalue(sigmaS0 * sigma2_i);
                    auto gamma_gamma0 = getvalue(gammaS0 * t_gmean_i());

                    auto gamma_sigma =
                        getvalue(TranspMult(t_gmean_i(), SmD) * sigma2_i) +
                        getvalue(p_P_mean() *
                                 (elemMult(t_gtotal_ij(), t_gvar_ij()) * u_irt));
                    auto sigma_sigma =
                        getvalue(TranspMult(sigma2_i, SmD) * sigma2_i) +
                        getvalue(p_P_mean() *
                                 (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u_irt));

                    auto p_iter = p_P_mean();
                    auto alfa_mu_start = build<trust_coefficient>(1.0 + 0.0 * p_iter[0]);
                    auto final_delta = dy;
                    auto final_V = r_y_var() - N * gSg;
                    auto final_beta = final_delta / final_V;
                    auto final_vSv = gSg;
                    auto final_start_k11 = final_V / N;
                    auto final_start_k12 = 0.0 * final_start_k11;
                    auto final_start_k22 = final_start_k11;

                    for (std::size_t iter = 0; iter < max_newton_steps_irt; ++iter) {
                        auto y_mean_iter =
                            N * getvalue(p_iter * t_gmean_i()) + y_baseline();
                        auto delta_iter = y - y_mean_iter;
                        auto V_iter = e + N * getvalue(p_iter * sigma2_i);
                        if (!(std::isfinite(primitive(V_iter))) ||
                            primitive(V_iter) <= 0.0)
                            return error_message("invalid MacroIRT rank-2 V", V_iter);

                        auto beta_iter = delta_iter / V_iter;
                        auto v_iter = t_gmean_i() + beta_iter * sigma2_i;
                        auto vS0 = gammaS0 + beta_iter * sigmaS0;
                        auto vSv0 = gamma_gamma0 + 2.0 * beta_iter * gamma_sigma0 +
                                     beta_iter * beta_iter * sigma_sigma0;
                        auto b0 = gamma_sigma0 + beta_iter * sigma_sigma0;
                        auto c0 = sigma_sigma0;

                        auto m11 = V_iter / N + vSv0;
                        auto m12 = b0;
                        auto m22 = c0 - 2.0 * V_iter * V_iter / N;
                        auto det = m11 * m22 - m12 * m12;
                        if (!(std::isfinite(primitive(det))) ||
                            std::abs(primitive(det)) <= 1e-30)
                            return error_message("singular MacroIRT rank-2 start Woodbury matrix",
                                                 det);

                        auto k11 = m22 / det;
                        auto k12 = (0.0 - m12) / det;
                        auto k22 = m11 / det;

                        auto delta_p = p_iter - p_P_mean();
                        auto delta_v = getvalue(delta_p * v_iter);
                        auto delta_s = getvalue(delta_p * sigma2_i);

                        auto q_gamma = 2.0 * delta_iter / V_iter;
                        auto q_sigma = delta_iter * delta_iter / (V_iter * V_iter) -
                                       1.0 / V_iter;
                        auto qS0 = q_gamma * gammaS0 + q_sigma * sigmaS0;
                        auto qv = getvalue(qS0 * v_iter);
                        auto qs = getvalue(qS0 * sigma2_i);

                        auto delta_term =
                            (delta_v * k11 + delta_s * k12) * vS0 +
                            (delta_v * k12 + delta_s * k22) * sigmaS0;
                        auto qS_post =
                            qS0 - (qv * k11 + qs * k12) * vS0 -
                            (qv * k12 + qs * k22) * sigmaS0;
                        auto p_candidate = p_P_mean() + delta_term + 0.5 * qS_post;
                        auto mean_step = p_candidate - p_iter;

                        alfa_mu_start = calculate_trust_coefficient(
                            p_iter, mean_step, trust_multiplying_factor);

                        auto Maybe_next_p =
                            to_Probability(p_iter + alfa_mu_start() * mean_step);
                        if (!Maybe_next_p)
                            return Maybe_next_p.error();
                        p_iter = std::move(Maybe_next_p.value());
                        final_delta = delta_iter;
                        final_V = V_iter;
                        final_beta = beta_iter;
                        final_vSv = vSv0;
                        final_start_k11 = k11;
                        final_start_k12 = k12;
                        final_start_k22 = k22;
                        const double step_inf = var::max(apply(
                            [](auto const& v) { return std::abs(v); },
                            primitive(alfa_mu_start() * mean_step)));
                        if (step_inf < newton_tol_irt)
                            break;
                    }

                    auto y_mean_final =
                        N * getvalue(p_iter * t_gmean_i()) + y_baseline();
                    final_delta = y - y_mean_final;
                    final_V = e + N * getvalue(p_iter * sigma2_i);
                    if (!(std::isfinite(primitive(final_V))) ||
                        primitive(final_V) <= 0.0)
                        return error_message("invalid final MacroIRT rank-2 V", final_V);

                    final_beta = final_delta / final_V;
                    auto v_final = t_gmean_i() + final_beta * sigma2_i;
                    auto vS0_final = gammaS0 + final_beta * sigmaS0;
                    auto vS_final = gS + final_beta * sigmaS;
                    auto vSv0_final =
                        gamma_gamma0 + 2.0 * final_beta * gamma_sigma0 +
                        final_beta * final_beta * sigma_sigma0;
                    auto b0_final = gamma_sigma0 + final_beta * sigma_sigma0;
                    auto c0_final = sigma_sigma0;

                    auto m11_start_final = final_V / N + vSv0_final;
                    auto m12_start_final = b0_final;
                    auto m22_start_final = c0_final - 2.0 * final_V * final_V / N;
                    auto det_start_final =
                        m11_start_final * m22_start_final -
                        m12_start_final * m12_start_final;
                    if (!(std::isfinite(primitive(det_start_final))) ||
                        std::abs(primitive(det_start_final)) <= 1e-30)
                        return error_message("singular final MacroIRT start Woodbury matrix",
                                             det_start_final);
                    final_start_k11 = m22_start_final / det_start_final;
                    final_start_k12 = (0.0 - m12_start_final) / det_start_final;
                    final_start_k22 = m11_start_final / det_start_final;

                    auto vSv_final = gSg + 2.0 * final_beta * gamma_sigma +
                                     final_beta * final_beta * sigma_sigma;
                    auto b_final = gamma_sigma + final_beta * sigma_sigma;
                    auto c_final = sigma_sigma;

                    auto m11_final = final_V / N + vSv_final;
                    auto m12_final = b_final;
                    auto m22_final = c_final - 2.0 * final_V * final_V / N;
                    auto det_final =
                        m11_final * m22_final - m12_final * m12_final;
                    if (!(std::isfinite(primitive(det_final))) ||
                        std::abs(primitive(det_final)) <= 1e-30)
                        return error_message("singular final MacroIRT rank-2 Woodbury matrix",
                                             det_final);

                    auto k11_final = m22_final / det_final;
                    auto k12_final = (0.0 - m12_final) / det_final;
                    auto k22_final = m11_final / det_final;
                    final_vSv = vSv_final;

                    auto delta_p = p_iter - p_P_mean();
                    auto delta_v = getvalue(delta_p * v_final);
                    auto delta_s = getvalue(delta_p * sigma2_i);

                    auto q_gamma = 2.0 * final_delta / final_V;
                    auto q_sigma = final_delta * final_delta / (final_V * final_V) -
                                   1.0 / final_V;
                    auto qS = q_gamma * gS + q_sigma * sigmaS;
                    auto gamma_v = gSg + final_beta * gamma_sigma;
                    auto sigma_v = gamma_sigma + final_beta * sigma_sigma;
                    auto qv = q_gamma * gamma_v + q_sigma * sigma_v;
                    auto qs = q_gamma * gamma_sigma + q_sigma * sigma_sigma;

                    auto delta_term =
                        (delta_v * k11_final + delta_s * k12_final) * vS_final +
                        (delta_v * k12_final + delta_s * k22_final) * sigmaS;
                    auto qS_post =
                        qS - (qv * k11_final + qs * k12_final) * vS_final -
                        (qv * k12_final + qs * k22_final) * sigmaS;
                    auto mu_prior_end = p_P_mean() * t_P();
                    auto p_candidate_end = mu_prior_end + delta_term + 0.5 * qS_post;
                    auto mean_step_end = p_candidate_end - mu_prior_end;
                    auto alfa_mu_exact = calculate_trust_coefficient(
                        mu_prior_end, mean_step_end, trust_multiplying_factor);

                    auto cov_downdate_end =
                        k11_final * XTX(vS_final) +
                        k12_final * X_plus_XT(TranspMult(vS_final, sigmaS)) +
                        k22_final * XTX(sigmaS);
                    auto cov_downdate_start =
                        final_start_k11 * XTX(vS0_final) +
                        final_start_k12 * X_plus_XT(TranspMult(vS0_final, sigmaS0)) +
                        final_start_k22 * XTX(sigmaS0);

                    auto alfa_sigma_exact = calculate_psd_trust_coefficient(
                        sigma_pre, cov_downdate_end, trust_multiplying_factor);
                    auto alfa_sigma_start = calculate_psd_trust_coefficient(
                        p_P_Cov(), cov_downdate_start, trust_multiplying_factor);

                    auto Maybe_r_P_mean =
                        to_Probability(mu_prior_end + alfa_mu_exact() * mean_step_end);
                    if (!Maybe_r_P_mean)
                        return Maybe_r_P_mean.error();
                    auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

                    auto Maybe_r_P_cov = to_Covariance_Probability(
                        sigma_pre - alfa_sigma_exact() * cov_downdate_end);
                    if (!Maybe_r_P_cov)
                        return Maybe_r_P_cov.error();
                    auto r_P_cov = build<P_Cov>(std::move(Maybe_r_P_cov.value()));

                    auto Maybe_P_cov_t10_y1 = to_Covariance_Probability(
                        p_P_Cov() - alfa_sigma_start() * cov_downdate_start);
                    if (!Maybe_P_cov_t10_y1)
                        return Maybe_P_cov_t10_y1.error();

                    if (!all_Probability_elements(primitive(r_P_mean())) ||
                        !all_Covariance_elements(primitive(r_P_cov()))) {
                        return error_message("error in P_mean or P_cov (rank-2 IRT)");
                    }

                    double taylor_strength_p = [&]() {
                        double sigma_sq = std::abs(primitive(
                            getvalue(TranspMult(sigma2_i, sigma2_i))));
                        double gamma_sq = std::abs(primitive(
                            getvalue(TranspMult(t_gmean_i(), t_gmean_i()))));
                        return std::abs(primitive(final_beta)) *
                               std::sqrt(sigma_sq) /
                               std::max(std::sqrt(gamma_sq), 1e-30);
                    }();

                    if constexpr (!dynamic) {
                        Transfer_Op_to<C_Patch_State, Algo_State> out;
                        get<y_mean>(out()) = std::move(r_y_mean);
                        get<y_var>(out()) = std::move(r_y_var);
                        get<r_std>(out()) = std::move(r_r_std);
                        get<Chi2>(out()) = std::move(chi2);
                        get<P_mean>(out())() = std::move(r_P_mean());
                        get<P_Cov>(out())() = std::move(r_P_cov());
                        get<trust_coefficient>(out()) = alfa_mu_exact;
                        get<taylor_trust_coefficient>(out()) =
                            build<taylor_trust_coefficient>(alfa_sigma_exact());
                        get<taylor_vSv>(out()) = build<taylor_vSv>(final_vSv);
                        get<taylor_strength>(out()) =
                            taylor_strength(taylor_strength_p);
                        return out;
                    } else {
                        Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;
                        get<y_mean>(out()) = std::move(r_y_mean);
                        get<y_var>(out()) = std::move(r_y_var);
                        get<Chi2>(out()) = std::move(chi2);
                        get<r_std>(out()) = std::move(r_r_std);
                        get<trust_coefficient>(out()) = alfa_mu_exact;
                        get<taylor_trust_coefficient>(out()) =
                            build<taylor_trust_coefficient>(alfa_sigma_exact());
                        get<taylor_vSv>(out()) = build<taylor_vSv>(final_vSv);
                        get<taylor_strength>(out()) =
                            taylor_strength(taylor_strength_p);
                        get<P>(out()) = get<P>(t_Qdt);
                        get<gmean_i>(out()) = get<gmean_i>(t_Qdt);
                        get<gvar_i>(out()) = get<gvar_i>(t_Qdt);
                        get<gtotal_ij>(out()) = get<gtotal_ij>(t_Qdt);
                        get<gmean_ij>(out()) = get<gmean_ij>(t_Qdt);

                        auto Maybe_r_P_mean_t11_y0 = to_Probability(mu_prior_end);
                        if (!Maybe_r_P_mean_t11_y0)
                            return Maybe_r_P_mean_t11_y0.error();
                        auto Maybe_r_P_mean_t10_y1 = to_Probability(p_iter);
                        if (!Maybe_r_P_mean_t10_y1)
                            return Maybe_r_P_mean_t10_y1.error();
                        auto Maybe_r_P_cov_t11_y0 =
                            to_Covariance_Probability(sigma_pre);
                        if (!Maybe_r_P_cov_t11_y0)
                            return Maybe_r_P_cov_t11_y0.error();

                        auto r_P_mean_0t_y0 = diag(p_P_mean()) * t_P();
                        auto r_P_mean_0t_y1 = diag(p_iter) * t_P();
                        auto r_P_cross_cov_0t_y0 =
                            SmD * t_P() + diag(p_P_mean()) * t_P();
                        auto SmD1 = Maybe_P_cov_t10_y1.value() - diag(p_iter);
                        auto r_P_cross_cov_0t_y1 =
                            SmD1 * t_P() + diag(p_iter) * t_P();

                        auto GS =
                            diag(TranspMult(t_gmean_i(), SmD)) * t_P() +
                            diag(p_P_mean()) * t_gtotal_ij();
                        if (auto Maybe_GS_check = to_Probability_displacement(GS);
                            !Maybe_GS_check)
                            return Maybe_GS_check.error();
                        get<d_GS>(out())() = std::move(GS);

                        get<P_mean_t11_y0>(out())() =
                            std::move(Maybe_r_P_mean_t11_y0.value());
                        get<P_mean_t10_y1>(out())() =
                            std::move(Maybe_r_P_mean_t10_y1.value());
                        get<P_mean_t20_y1>(out())() = std::move(r_P_mean());
                        get<P_Cov_t20_y1>(out())() = std::move(r_P_cov());
                        get<P_mean_0t_y0>(out())() = std::move(r_P_mean_0t_y0);
                        get<P_mean_0t_y1>(out())() = std::move(r_P_mean_0t_y1);
                        get<P_cross_cov_0t_y0>(out())() =
                            std::move(r_P_cross_cov_0t_y0);
                        get<P_cross_cov_0t_y1>(out())() =
                            std::move(r_P_cross_cov_0t_y1);
                        get<P_Cov_t11_y0>(out())() =
                            std::move(Maybe_r_P_cov_t11_y0.value());
                        get<P_Cov_t10_y1>(out())() =
                            std::move(Maybe_P_cov_t10_y1.value());
                        return out;
                    }
                }

            // V_obs = ε² + N·μ·σ̄² (measurement-noise piece at r=μ_prior).
            // This is the supplement's V in the energy / Hessian / Newton-step
            // derivation. NOT the predictive variance r_y_var (= V_pred =
            // V_obs + N·γ̃ᵀΣγ̃). All β = δ/V_obs scalings below derive from
            // ∂E/∂r at r=μ; using δ/V_pred (as the supplement's pseudo-code
            // does) under-counts the σ² Taylor correction by a factor
            // V_obs/V_pred, which is small precisely in the high-N regime
            // where channel noise dominates.
            auto V_obs = r_y_var() - N * gSg;

            // gvar_i flavor for v: residual for av=2, total for av=1.
            // (Recomputed here as a vector; the +N·μ·gvar_i scalar above
            //  uses its own local lambda — see the in-place fix block.)
            auto gvar_i_for_v = [&]() {
                auto& t_gsqr_i_loc = get<gsqr_i>(t_Qdt);
                if constexpr (averaging::value == 2) {
                    auto& t_gtotal_ij_loc = get<gtotal_ij>(t_Qdt);
                    auto& t_gmean_ij_loc = get<gmean_ij>(t_Qdt);
                    Matrix<double> u_loc(p_P_mean().size(), 1, 1.0);
                    return t_gsqr_i_loc() -
                           elemMult(t_gtotal_ij_loc(), t_gmean_ij_loc()) * u_loc;
                } else {
                    return t_gsqr_i_loc() - elemMult(t_gmean_i(), t_gmean_i());
                }
            }();

            // -----------------------------------------------------------
            // α_vSv (Taylor trust coefficient): largest α ∈ [0,1] keeping
            //   vSv_eff(α) = a + 2·α·β·b + α²·β²·c  ≥ ε·a > 0
            // where β = δ/V and the three tilde scalars are
            //   a = γ̃ᵀΣγ      (= existing gSg)
            //   b = γ̃ᵀΣσ      (cross tilde scalar, NEW)
            //   c = σ̃ᵀΣσ      (σ²-direction tilde scalar, NEW)
            // Cauchy-Schwarz gives b² ≤ a·c, so vSv_eff(α) ≥ 0 in exact
            // arithmetic and α_vSv = 1. Numerical cancellation can produce
            // b² > a·c slightly; in that case the closed-form quadratic
            // root gives the largest α below the dip.
            //
            // α_vSv plays the role of variance inflation by 1/α_vSv on the
            // σ² contribution to v: v_eff = γ̄₀ + α_vSv·(δ/V)·σ̄²₀.  Smoothly
            // interpolates between full IRT (α_vSv=1) and pure MacroIR
            // (α_vSv=0; σ² Taylor correction disabled).
            // -----------------------------------------------------------
            // X̃ᵀΣỸ closed form for tilde scalars on conditional quantities X_ij, Y_ij:
            //   X̃ᵀΣỸ = −(p·X̄)·(p·Ȳ) + Σᵢⱼ pᵢ·P_ij·X_ij·Y_ij
            //         = TranspMult(X̄,SmD)·Ȳ + p · (gtotal_X_ij ∘ Y_ij) · u
            // Note the SECOND factor in elemMult is the CONDITIONAL Y_ij
            // (gvar_ij here), not the joint gtotal_var_ij — using the joint
            // would inject an extra P_ij and make the second term too small,
            // letting the negative SmD piece dominate and violating
            // Cauchy-Schwarz (b² ≤ a·c). Earlier draft had this bug; figure_2
            // macro_IRT was producing negative vSv → invalid posterior.
            auto b_tilde = [&]() {
                if constexpr (averaging::value == 2) {
                    auto& t_gtotal_ij_loc = get<gtotal_ij>(t_Qdt);
                    auto& t_gvar_ij_loc = get<gvar_ij>(t_Qdt);
                    Matrix<double> u_loc(p_P_mean().size(), 1, 1.0);
                    return getvalue(TranspMult(t_gmean_i(), SmD) * gvar_i_for_v) +
                           getvalue(p_P_mean() *
                                    (elemMult(t_gtotal_ij_loc(),
                                              t_gvar_ij_loc()) * u_loc));
                } else {
                    return getvalue(TranspMult(t_gmean_i(), p_P_Cov()) * gvar_i_for_v);
                }
            }();
            auto c_tilde = [&]() {
                if constexpr (averaging::value == 2) {
                    auto& t_gtotal_var_ij_loc = get<gtotal_var_ij>(t_Qdt);
                    auto& t_gvar_ij_loc = get<gvar_ij>(t_Qdt);
                    Matrix<double> u_loc(p_P_mean().size(), 1, 1.0);
                    return getvalue(TranspMult(gvar_i_for_v, SmD) * gvar_i_for_v) +
                           getvalue(p_P_mean() *
                                    (elemMult(t_gtotal_var_ij_loc(),
                                              t_gvar_ij_loc()) * u_loc));
                } else {
                    return getvalue(TranspMult(gvar_i_for_v, p_P_Cov()) * gvar_i_for_v);
                }
            }();

            // Closed-form α_vSv via the quadratic
            //   c·β²·α² + 2·b·β·α + (a − ε·a) ≥ 0
            // Cases on primitives (binding decision); evaluate the binding
            // expression with derivatives intact. We use the
            //   `value + 0.0 * (derivative-aware quantity)`
            // pattern (cf. calculate_trust_coefficient) to construct
            // derivative-aware constants of the right type.
            constexpr double eps_vSv_relative = 1e-6;
            auto beta_vSv = dy / V_obs;
            using std::sqrt;
            auto alfa_vSv = [&]() {
                auto a_p = primitive(gSg);
                auto b_p = primitive(b_tilde);
                auto c_p = primitive(c_tilde);
                auto beta_p = primitive(beta_vSv);
                double eps_floor = eps_vSv_relative * std::max(a_p, 1e-30);
                // Derivative-aware "constants" of the same type as beta_vSv:
                //   make_alfa(1) → derivative-aware 1.0 with zero derivative
                //   make_alfa(0) → derivative-aware 0.0 with zero derivative
                auto zero_alfa = 0.0 * beta_vSv;
                auto one_alfa  = 1.0 + zero_alfa;
                if (!(c_p > 0) || !(beta_p * beta_p * c_p > 0)) {
                    // No σ²-direction contribution → quadratic degenerates;
                    // vSv_eff(α) ≡ a stays positive → full IRT.
                    return one_alfa;
                }
                double disc_p = b_p * b_p * beta_p * beta_p
                              - c_p * beta_p * beta_p * (a_p - eps_floor);
                if (disc_p <= 0.0) {
                    // Cauchy-Schwarz holds → quadratic stays ≥ ε·a → full IRT.
                    return one_alfa;
                }
                // Quadratic dips below ε·a between two real roots.
                // Roots: α± = (−β·b ± √disc) / (β²·c). We want the smaller
                // positive root (largest α below the dip).
                double sqrt_disc_p = std::sqrt(disc_p);
                double inv_two_a_p = 1.0 / (beta_p * beta_p * c_p);
                double r1_p = (-beta_p * b_p - sqrt_disc_p) * inv_two_a_p;
                double r2_p = (-beta_p * b_p + sqrt_disc_p) * inv_two_a_p;
                double root_low_p = std::min(r1_p, r2_p);
                double root_high_p = std::max(r1_p, r2_p);
                if (root_high_p <= 0.0) {
                    // Dip is entirely at α < 0 → safe for all α ≥ 0.
                    return one_alfa;
                }
                if (root_low_p <= 0.0) {
                    // Dip extends from negative to root_high_p > 0 → no safe
                    // α below the dip; fall back to no Taylor correction.
                    return zero_alfa;
                }
                if (root_low_p >= 1.0) {
                    return one_alfa;
                }
                // Evaluate the binding root derivative-aware: re-compute via
                // the closed form on the derivative-tracked b, c, β.
                auto disc = b_tilde * b_tilde * beta_vSv * beta_vSv
                          - c_tilde * beta_vSv * beta_vSv * (gSg - eps_floor);
                auto sqrt_disc = sqrt(disc);
                auto inv_two_a = 1.0 / (beta_vSv * beta_vSv * c_tilde);
                // Pick the same branch the primitive selected (smaller root).
                auto neg_beta_b = (0.0 - beta_vSv) * b_tilde;
                auto root_low_d  = (neg_beta_b - sqrt_disc) * inv_two_a;
                auto root_high_d = (neg_beta_b + sqrt_disc) * inv_two_a;
                if (primitive(root_low_d) > primitive(root_high_d))
                    return root_high_d;
                return root_low_d;
            }();

            // Effective direction v = γ̄₀ + α_vSv·(δ/V_obs)·σ̄²₀  (per-i₀, K-dim)
            auto v = t_gmean_i() + alfa_vSv * (dy / V_obs) * gvar_i_for_v;

            // Tilde scalar vSv = ṽᵀΣv  (mirrors gSg computation; with α_vSv-shrunk v)
            // vSv = ṽᵀΣv tilde — same X̃ᵀΣỸ closed form as b_tilde/c_tilde:
            //   −(p·v̄)² + Σᵢⱼ pᵢ·P_ij·v_ij²
            //   = TranspMult(v̄,SmD)·v̄ + p · (gtotal_v_ij ∘ v_ij) · u
            // Need v_ij = gmean_ij + α_vSv·β·gvar_ij at the (i,j) level so
            // elemMult is joint × conditional, not joint × joint.
            auto vSv = [&]() {
                if constexpr (averaging::value == 2) {
                    auto& t_gtotal_ij_loc = get<gtotal_ij>(t_Qdt);
                    auto& t_gtotal_var_ij_loc = get<gtotal_var_ij>(t_Qdt);
                    auto& t_gmean_ij_loc = get<gmean_ij>(t_Qdt);
                    auto& t_gvar_ij_loc = get<gvar_ij>(t_Qdt);
                    auto gtotal_v_ij = t_gtotal_ij_loc() +
                                       alfa_vSv * (dy / V_obs) *
                                           t_gtotal_var_ij_loc();
                    auto v_ij = t_gmean_ij_loc() +
                                alfa_vSv * (dy / V_obs) *
                                    t_gvar_ij_loc();
                    Matrix<double> u_loc(p_P_mean().size(), 1, 1.0);
                    return getvalue(TranspMult(v, SmD) * v) +
                           getvalue(p_P_mean() *
                                    (elemMult(gtotal_v_ij, v_ij) * u_loc));
                } else {
                    return getvalue(TranspMult(v, p_P_Cov()) * v);
                }
            }();

            // Tilde vector vS = ṽᵀΣ in endpoint frame  (mirrors gS computation; with α_vSv-shrunk v)
            auto vS = [&]() {
                if constexpr (averaging::value == 2) {
                    auto& t_gtotal_ij_loc = get<gtotal_ij>(t_Qdt);
                    auto& t_gtotal_var_ij_loc = get<gtotal_var_ij>(t_Qdt);
                    auto gtotal_v_ij = t_gtotal_ij_loc() +
                                       alfa_vSv * (dy / V_obs) *
                                           t_gtotal_var_ij_loc();
                    return TranspMult(v, SmD) * t_P() + p_P_mean() * gtotal_v_ij;
                } else {
                    return TranspMult(v, p_P_Cov()) * t_P();
                }
            }();

            // Newton step from r₀ = μ_prior. Start-frame Newton step
            // (eq mu_post_MRT, supplement section 4):
            //     Δμ_start_col = (δ/(2V_obs))·Σ_p_post·(γ̄+v)
            // Propagating to endpoint frame as a row vector and expanding
            // the SM down-date in Σ_p_post:
            //     Δμ_end = (δ/(2V_obs))·(γ̄+v)ᵀ·Σ_p_post·P
            //            = (δ/(2V_obs))·[(γ̄+v)ᵀΣP − sm·(γ̄+v)ᵀΣv · vᵀΣP]
            //            = (δ/(2V_obs))·[(gS+vS) − sm·(b' + vSv)·vS]
            // where b' := γ̃ᵀΣv = gSg + α·β·b_tilde and β = δ/V_obs.
            //
            // V_obs is computed above (= ε² + N·μ·σ̄², the supplement's V
            // in the energy/Hessian derivation, not the predictive r_y_var).
            auto sm_factor = N / (V_obs + N * vSv);
            auto cov_downdate = sm_factor * XTX(vS);

            // b' = γ̃ᵀΣv tilde scalar. From γ̃ᵀΣv = γ̃ᵀΣγ̄ + α·β·γ̃ᵀΣσ̄²
            //   b' = gSg + (α_vSv · δ / V_obs) · b_tilde
            auto bv_tilde = gSg + alfa_vSv * (dy / V_obs) * b_tilde;
            auto mean_dir = (dy / (2 * V_obs)) *
                            ((gS + vS) - sm_factor * (bv_tilde + vSv) * vS);

            // Trust region with the IRT/MRT directions
            auto alfa_mu_irt = calculate_trust_coefficient(
                p_P_mean() * t_P(), mean_dir, trust_multiplying_factor);
            auto alfa_sigma_irt = calculate_psd_trust_coefficient(
                sigma_pre, vS, sm_factor, trust_multiplying_factor);
            auto alfa = build<trust_coefficient>(softmin(
                alfa_mu_irt(), alfa_sigma_irt(), trust_softmin_eps_irt));

            // DECOUPLED trust: mean uses α_μ alone, covariance uses α_σ alone
            // (independent — no shared softmin). α_σ is RETAINED here (unlike the
            // standard rank-1 branch) because this Taylor/IRT down-date is not
            // self-limiting. `alfa` (the combined value) is kept only for the
            // diagnostic dump / error logs below.
            auto Maybe_r_P_mean =
                to_Probability(p_P_mean() * t_P() + alfa_mu_irt() * mean_dir);
            if (!Maybe_r_P_mean) {
                std::cerr << "[IRT-fail Pmean] av=" << averaging::value
                          << " N=" << primitive(N)
                          << " alpha_vSv=" << primitive(alfa_vSv)
                          << " gSg=" << primitive(gSg)
                          << " b_tilde=" << primitive(b_tilde)
                          << " c_tilde=" << primitive(c_tilde)
                          << " vSv=" << primitive(vSv)
                          << " sm_factor=" << primitive(sm_factor)
                          << " dy=" << primitive(dy)
                          << " r_y_var=" << primitive(r_y_var())
                          << " alfa_mu=" << primitive(alfa_mu_irt())
                          << " alfa_sigma=" << primitive(alfa_sigma_irt())
                          << " alfa_final=" << primitive(alfa())
                          << " err=" << Maybe_r_P_mean.error()()
                          << "\n";
                return Maybe_r_P_mean.error();
            }
            auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

            auto Maybe_r_P_cov =
                to_Covariance_Probability(sigma_pre - alfa_sigma_irt() * cov_downdate);  // decoupled: Σ uses α_σ
            if (!Maybe_r_P_cov) {
                std::cerr << "[IRT-fail Pcov] av=" << averaging::value
                          << " N=" << primitive(N)
                          << " alpha_vSv=" << primitive(alfa_vSv)
                          << " gSg=" << primitive(gSg)
                          << " b_tilde=" << primitive(b_tilde)
                          << " c_tilde=" << primitive(c_tilde)
                          << " vSv=" << primitive(vSv)
                          << " sm_factor=" << primitive(sm_factor)
                          << " dy=" << primitive(dy)
                          << " r_y_var=" << primitive(r_y_var())
                          << " alfa_mu=" << primitive(alfa_mu_irt())
                          << " alfa_sigma=" << primitive(alfa_sigma_irt())
                          << " alfa_final=" << primitive(alfa())
                          << " err=" << Maybe_r_P_cov.error()()
                          << "\n";
                return Maybe_r_P_cov.error();
            }
            auto r_P_cov = build<P_Cov>(std::move(Maybe_r_P_cov.value()));

            if (!all_Probability_elements(primitive(r_P_mean())) ||
                !all_Covariance_elements(primitive(r_P_cov()))) {
                return error_message("error in P_mean or P_cov (IRT/MRT)");
            }
            auto r_logL = calculate_logL(false, r_y_var, chi2, m);
            auto r_elogL = calculate_elogL(false, r_y_var, m);

            // Diagnostic: relative magnitude of the σ² perturbation in v vs.
            // the baseline γ̄ direction.
            //   taylor_strength = α_vSv · |β| · ‖σ̄²‖₂ / ‖γ̄‖₂      (β = δ/V)
            // Computed primitive-only — it's a diagnostic, doesn't need to
            // carry derivatives.
            double taylor_strength_p = [&]() {
                double sigma_sq = std::abs(primitive(
                    getvalue(TranspMult(gvar_i_for_v, gvar_i_for_v))));
                double gamma_sq = std::abs(primitive(
                    getvalue(TranspMult(t_gmean_i(), t_gmean_i()))));
                double beta_p = primitive(dy / r_y_var());
                double alfa_p = primitive(alfa_vSv);
                return std::abs(alfa_p * beta_p) *
                       std::sqrt(sigma_sq) /
                       std::max(std::sqrt(gamma_sq), 1e-30);
            }();

            if constexpr (!dynamic) {
                Transfer_Op_to<C_Patch_State, Algo_State> out;
                get<y_mean>(out()) = std::move(r_y_mean);
                get<y_var>(out()) = std::move(r_y_var);
                get<r_std>(out()) = std::move(r_r_std);
                get<Chi2>(out()) = std::move(chi2);
                get<P_mean>(out())() = std::move(r_P_mean());
                get<P_Cov>(out())() = std::move(r_P_cov());
                get<trust_coefficient>(out()) = alfa;
                get<taylor_trust_coefficient>(out()) =
                    build<taylor_trust_coefficient>(alfa_vSv);
                get<taylor_vSv>(out()) = build<taylor_vSv>(vSv);
                get<taylor_strength>(out()) = taylor_strength(taylor_strength_p);
                return out;
            } else {
                // Dynamic case: main outputs are IRT/MRT-specific; derivative
                // -path extras (d_GS, P_mean_t11_y0, etc.) reuse the standard
                // Kalman intermediates (chi, gS) with the IRT alfa as
                // placeholders. See TODO in the block header.
                Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;
                get<y_mean>(out()) = std::move(r_y_mean);
                get<y_var>(out()) = std::move(r_y_var);
                get<Chi2>(out()) = std::move(chi2);
                get<r_std>(out()) = std::move(r_r_std);
                get<trust_coefficient>(out()) = alfa;
                get<taylor_trust_coefficient>(out()) =
                    build<taylor_trust_coefficient>(alfa_vSv);
                get<taylor_vSv>(out()) = build<taylor_vSv>(vSv);
                get<taylor_strength>(out()) = taylor_strength(taylor_strength_p);
                get<P>(out()) = get<P>(t_Qdt);
                get<gmean_i>(out()) = get<gmean_i>(t_Qdt);
                get<gvar_i>(out()) = get<gvar_i>(t_Qdt);
                if constexpr (averaging::value == 2) {
                    auto& t_gmean_i_loc = get<gmean_i>(t_Qdt);
                    auto gS0 = TranspMult(t_gmean_i_loc(), SmD) +
                               elemMult(p_P_mean(), t_gmean_i_loc());
                    if (auto Maybe_gS0_check = to_Probability_displacement(gS0);
                        !Maybe_gS0_check)
                        return Maybe_gS0_check.error();
                    auto Maybe_r_P_mean_t11_y0 = to_Probability(p_P_mean() * t_P());
                    auto Maybe_r_P_mean_t10_y1 =
                        to_Probability(p_P_mean() + alfa() * chi * gS0);
                    if (!Maybe_r_P_mean_t11_y0)
                        return Maybe_r_P_mean_t11_y0.error();
                    if (!Maybe_r_P_mean_t10_y1)
                        return Maybe_r_P_mean_t10_y1.error();
                    auto r_P_mean_0t_y0 = diag(p_P_mean()) * t_P();
                    auto& t_gtotal_ij_loc = get<gtotal_ij>(t_Qdt);
                    auto GS = diag(TranspMult(t_gmean_i_loc(), SmD)) * t_P() +
                              diag(p_P_mean()) * t_gtotal_ij_loc();
                    if (auto Maybe_GS_check = to_Probability_displacement(GS);
                        !Maybe_GS_check)
                        return Maybe_GS_check.error();
                    auto r_P_mean_0t_y1 = r_P_mean_0t_y0 + alfa() * chi * GS;
                    auto r_P_cross_cov_0t_y0 =
                        SmD * t_P() + diag(p_P_mean()) * t_P();
                    auto r_P_cross_cov_0t_y1 = r_P_cross_cov_0t_y0 -
                        (alfa() * N / r_y_var()) * TranspMult(gS0, gS);
                    get<d_GS>(out())() = std::move(GS);
                    get<P_mean_t11_y0>(out())() =
                        std::move(Maybe_r_P_mean_t11_y0.value());
                    get<P_mean_t10_y1>(out())() =
                        std::move(Maybe_r_P_mean_t10_y1.value());
                    get<P_mean_t20_y1>(out())() = std::move(r_P_mean());
                    get<P_Cov_t20_y1>(out())() = std::move(r_P_cov());
                    get<gtotal_ij>(out()) = get<gtotal_ij>(t_Qdt);
                    get<gmean_ij>(out()) = get<gmean_ij>(t_Qdt);
                    get<P_mean_0t_y0>(out())() = std::move(r_P_mean_0t_y0);
                    get<P_mean_0t_y1>(out())() = std::move(r_P_mean_0t_y1);
                    get<P_cross_cov_0t_y0>(out())() = std::move(r_P_cross_cov_0t_y0);
                    get<P_cross_cov_0t_y1>(out())() = std::move(r_P_cross_cov_0t_y1);
                    auto Maybe_r_P_cov_t11_y0 = to_Covariance_Probability(
                        AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));
                    auto Maybe_r_P_cov_t10_y1 = to_Covariance_Probability(
                        get<P_Cov>(t_prior())() -
                        (alfa() * N / r_y_var()) * XTX(gS0));
                    if (!Maybe_r_P_cov_t11_y0)
                        return Maybe_r_P_cov_t11_y0.error();
                    if (!Maybe_r_P_cov_t10_y1)
                        return Maybe_r_P_cov_t10_y1.error();
                    get<P_Cov_t11_y0>(out())() =
                        std::move(Maybe_r_P_cov_t11_y0.value());
                    get<P_Cov_t10_y1>(out())() =
                        std::move(Maybe_r_P_cov_t10_y1.value());
                } else {  // averaging::value == 1
                    auto Maybe_r_P_mean_t2_y0 = to_Probability(p_P_mean() * t_P());
                    if (!Maybe_r_P_mean_t2_y0)
                        return Maybe_r_P_mean_t2_y0.error();
                    auto Maybe_r_P_mean_t1_y1 =
                        to_Probability(p_P_mean() + alfa() * chi * gS);
                    if (!Maybe_r_P_mean_t1_y1)
                        return Maybe_r_P_mean_t1_y1.error();
                    auto Maybe_r_P_cov_t2_y0 = to_Covariance_Probability(
                        AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));
                    if (!Maybe_r_P_cov_t2_y0)
                        return Maybe_r_P_cov_t2_y0.error();
                    auto Maybe_r_P_cov_t1_y1 = to_Covariance_Probability(
                        get<P_Cov>(t_prior())() -
                        (alfa() * N / r_y_var()) * XTX(gS));
                    if (!Maybe_r_P_cov_t1_y1)
                        return Maybe_r_P_cov_t1_y1.error();
                    auto r_P_mean_0t_y0 = diag(p_P_mean()) * t_P();
                    auto r_P_mean_0t_y1 =
                        diag(Maybe_r_P_mean_t1_y1.value()) * t_P();
                    auto r_P_cross_cov_0t_y0 =
                        SmD * t_P() + diag(p_P_mean()) * t_P();
                    auto SmD1 = Maybe_r_P_cov_t1_y1.value() -
                                diag(Maybe_r_P_mean_t1_y1.value());
                    auto r_P_cross_cov_0t_y1 =
                        SmD1 * t_P() +
                        diag(Maybe_r_P_mean_t1_y1.value()) * t_P();
                    get<P_mean_0t_y0>(out())() = std::move(r_P_mean_0t_y0);
                    get<P_mean_0t_y1>(out())() = std::move(r_P_mean_0t_y1);
                    get<P_cross_cov_0t_y0>(out())() = std::move(r_P_cross_cov_0t_y0);
                    get<P_cross_cov_0t_y1>(out())() = std::move(r_P_cross_cov_0t_y1);
                    get<P_Cov_t2_y0>(out())() =
                        std::move(Maybe_r_P_cov_t2_y0.value());
                    get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                    get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                    get<P_mean_t2_y0>(out())() =
                        std::move(Maybe_r_P_mean_t2_y0.value());
                    get<P_mean_t1_y1>(out())() =
                        std::move(Maybe_r_P_mean_t1_y1.value());
                    get<d_gS>(out())() = std::move(gS);
                    get<P_Cov_t1_y1>(out())() =
                        std::move(Maybe_r_P_cov_t1_y1.value());
                }
                return out;
            }
            }
        }

        // ===========================================================
        // Standard rank-1 Kalman branch (vc=false, or unhandled av=0).
        // ===========================================================

        // α_μ : largest α keeping (μ·t_P + α·chi·gS) on the simplex. Endpoint
        // frame for all averaging values, since gS is now endpoint-frame too.
        // Pass derivative-aware inputs — the AD-aware overload of
        // calculate_trust_coefficient finds the binding index via primitives
        // (cheap) and evaluates the binding expression with derivatives intact.
        auto alfa_mu = calculate_trust_coefficient(
            p_P_mean() * t_P(), chi * gS, trust_multiplying_factor);

        // α_Σ : per-diagonal PSD trust. REDUNDANT in this rank-1 branch — the
        // undamped (α=1) down-date is PSD by construction (Markov-step Löwner margin;
        // see docs/math/trust_coefficient.md and the α_σ-necessity analysis). It is
        // NO LONGER APPLIED to Σ; kept only for the diagnostic dump. The real PSD
        // guard is the to_Covariance_Probability canary below. (It is genuinely needed
        // only in the variance_correction Taylor branches, which are handled elsewhere.)
        [[maybe_unused]] auto alfa_sigma = calculate_psd_trust_coefficient(
            sigma_pre, gS, N / r_y_var(), trust_multiplying_factor);

        // DECOUPLED trust: the mean uses α_μ alone; the covariance down-date is left
        // undamped (α=1, see the Σ updates below). α_μ no longer reaches Σ, so its
        // residual irregularity at residual zero-crossings is harmless.
        auto alfa = build<trust_coefficient>(alfa_mu());


        // Mean update at endpoint frame, unified across all averaging values.
        // Mathematically equivalent to the prior  (μ + α·chi·gS_obs) · t_P
        // formulation for avg=0/1, since  gS = gS_obs · t_P  by construction.
        auto Maybe_r_P_mean =
            to_Probability(p_P_mean() * t_P() + alfa() * chi * gS);
        if (!Maybe_r_P_mean)
            return Maybe_r_P_mean.error();
        auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

        auto Maybe_r_P_cov = to_Covariance_Probability(
            sigma_pre - (N / r_y_var()) * XTX(gS));  // decoupled: Σ undamped (α=1, α_σ dropped — PSD by construction)
        if (!Maybe_r_P_cov)
            return Maybe_r_P_cov.error();

        auto r_P_cov = build<P_Cov>(std::move(Maybe_r_P_cov.value()));

        if (!all_Probability_elements(primitive(r_P_mean())) ||
            !all_Covariance_elements(primitive(r_P_cov()))) {
            return error_message("error in P_mean or P_cov");
        }
        auto r_logL = calculate_logL(false, r_y_var, chi2,  m);

        auto r_elogL = calculate_elogL(false, r_y_var,  m);

        if constexpr (!dynamic) {
            Transfer_Op_to<C_Patch_State, Algo_State> out;
            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<r_std>(out())= std::move(r_r_std);

            get<Chi2>(out())= std::move(chi2);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())() = std::move(r_P_cov());
            get<r_std>(out())= std::move(r_r_std);

            get<trust_coefficient>(out()) = alfa;
            // DIAGNOSTIC: keep α_μ and α_Σ separately so the per-sample detailed
            // dump can tell which constraint's derivative jumps (mean-simplex vs
            // PSD) at the FD-instability spike.
            get<trust_mean_coefficient>(out())() = alfa_mu();
            get<trust_psd_coefficient>(out())() = alfa_sigma();
            // Standard Kalman branch — Taylor not active.
            get<taylor_trust_coefficient>(out()) = taylor_trust_coefficient(1.0);
            get<taylor_vSv>(out()) = taylor_vSv(0.0);
            get<taylor_strength>(out()) = taylor_strength(0.0);
            return out;
        } else {
            Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;
            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<Chi2>(out()) = std::move(chi2);
            get<r_std>(out())= std::move(r_r_std);
            
            get<trust_coefficient>(out()) = alfa;
            // Standard Kalman branch — Taylor not active.
            get<taylor_trust_coefficient>(out()) = taylor_trust_coefficient(1.0);
            get<taylor_vSv>(out()) = taylor_vSv(0.0);
            get<taylor_strength>(out()) = taylor_strength(0.0);
            if constexpr (averaging::value == 0) {
                get<P_half>(out())= get<P_half>(t_Qdt);
            } else{
                     get<P>(out())= get<P>(t_Qdt);
            }
            if constexpr (averaging::value>0) {
                get<gmean_i>(out())=  get<gmean_i>(t_Qdt);
                get<gvar_i>(out()) =  get<gvar_i>(t_Qdt);
            }
            if constexpr (averaging::value == 2) {


                auto& t_gmean_i = get<gmean_i>(t_Qdt);

                auto gS0 = TranspMult(t_gmean_i(), SmD )+
                          elemMult( p_P_mean(),t_gmean_i()) ;
                // gS0 is the start-side probability displacement: μ + α·chi·gS0 must
                // remain a probability, so gS0·u = 0.
                if (auto Maybe_gS0_check = to_Probability_displacement(gS0); !Maybe_gS0_check)
                    return Maybe_gS0_check.error();

                auto Maybe_r_P_mean_t11_y0 = to_Probability(p_P_mean() * t_P());
                auto Maybe_r_P_mean_t10_y1 = to_Probability(p_P_mean() + alfa() * chi * gS0);

                if (!Maybe_r_P_mean_t11_y0) {
                    return Maybe_r_P_mean_t11_y0.error();
                }
                if (!Maybe_r_P_mean_t10_y1) {
                    return Maybe_r_P_mean_t10_y1.error();
                }

                
                
                // Matrix form of the boundary-state prior mean:
                // M^-_{0,t}(i0,it) = P(X0=i0, Xt=it) = mu_0(i0) P_{i0->it}(t).
                auto r_P_mean_0t_y0= diag(p_P_mean())*t_P();
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);

                auto GS=  diag(TranspMult(t_gmean_i(), SmD)) * t_P() + diag(p_P_mean())* t_gtotal_ij();
                // GS is the joint-displacement matrix: M⁻ + α·chi·GS must remain a
                // joint distribution (total sum 1), so total sum of GS must be 0.
                // Stronger row/column identities (uᵀ·GS = gS, GS·u = gS0) are checked
                // by the asserts a few lines below.
                if (auto Maybe_GS_check = to_Probability_displacement(GS); !Maybe_GS_check)
                    return Maybe_GS_check.error();

                auto r_P_mean_0t_y1= r_P_mean_0t_y0 + alfa()*chi*GS;

                // Reduced start/end cross-covariance:
                // C^-_{0,t} = Cov(x_0, x_t) = Sigma_0 P(t).
                auto r_P_cross_cov_0t_y0= SmD*t_P()+diag(p_P_mean())*t_P(); 

                
                // Posterior reduced cross-covariance after conditioning on y_{0->t}.
                auto r_P_cross_cov_0t_y1= r_P_cross_cov_0t_y0 - (N/r_y_var())* TranspMult(gS0,gS);  // decoupled: α=1
                Matrix<double> uT(1ul,p_P_mean().size(), 1.0);
         
               
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y0).value(), Maybe_r_P_mean_t11_y0.value()));
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y1).value(), r_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y0)).value(), p_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y1)).value(), Maybe_r_P_mean_t10_y1.value()));

                assert(var::test_equality(uT*GS,gS));
                assert(var::test_equality(MultTransp(uT,GS),gS0));
                

                get<d_GS>(out())() = std::move(GS);
                get<P_mean_t11_y0>(out())() = std::move(Maybe_r_P_mean_t11_y0.value());
                get<P_mean_t10_y1>(out())() = std::move(Maybe_r_P_mean_t10_y1.value());
                
                get<P_mean_t20_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t20_y1>(out())() = std::move(r_P_cov());
                get<gtotal_ij>(out()) = get<gtotal_ij>(t_Qdt);
                get<gmean_ij>(out()) = get<gmean_ij>(t_Qdt);

                get<P_mean_0t_y0>(out())()= std::move(r_P_mean_0t_y0);
                get<P_mean_0t_y1>(out())()= std::move(r_P_mean_0t_y1);
                get<P_cross_cov_0t_y0>(out())()= std::move(r_P_cross_cov_0t_y0);
                get<P_cross_cov_0t_y1>(out())()= std::move(r_P_cross_cov_0t_y1);  


                
                
                auto Maybe_r_P_cov_t11_y0 =
                    to_Covariance_Probability(AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));
                auto Maybe_r_P_cov_t10_y1 = to_Covariance_Probability(
                    get<P_Cov>(t_prior())() - (N / r_y_var()) * XTX(gS0));  // decoupled: α=1

                if (!Maybe_r_P_cov_t11_y0) {
                    return Maybe_r_P_cov_t11_y0.error();
                }
                if (!Maybe_r_P_cov_t10_y1) {
                    return Maybe_r_P_cov_t10_y1.error();
                }

                get<P_Cov_t11_y0>(out())() = std::move(Maybe_r_P_cov_t11_y0.value());
                get<P_Cov_t10_y1>(out())() = std::move(Maybe_r_P_cov_t10_y1.value());

            } else  if constexpr (averaging::value ==1)
            {
        
                auto Maybe_r_P_mean_t2_y0 = to_Probability(p_P_mean() * t_P());
                if (!Maybe_r_P_mean_t2_y0) {
                    return Maybe_r_P_mean_t2_y0.error();
                }
                auto Maybe_r_P_mean_t1_y1 = to_Probability(p_P_mean() + alfa() *chi * gS);

                if (!Maybe_r_P_mean_t1_y1) {
                    return Maybe_r_P_mean_t1_y1.error();
                }

                auto Maybe_r_P_cov_t2_y0 =
                    to_Covariance_Probability(AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));
                if (!Maybe_r_P_cov_t2_y0) {
                    return Maybe_r_P_cov_t2_y0.error();  
                }
                auto Maybe_r_P_cov_t1_y1 = to_Covariance_Probability(
                    get<P_Cov>(t_prior())() - (N / r_y_var()) * XTX(gS));  // decoupled: α=1

                if (!Maybe_r_P_cov_t1_y1) {
                    return Maybe_r_P_cov_t1_y1.error();
                }
                
                // Matrix form of the boundary-state prior mean over the interval.
                auto r_P_mean_0t_y0= diag(p_P_mean())*t_P();
                
                auto r_P_mean_0t_y1= diag(Maybe_r_P_mean_t1_y1.value())*t_P();

                // Reduced start/end cross-covariance before conditioning on y_{0->t}.
                auto r_P_cross_cov_0t_y0= SmD*t_P()+diag(p_P_mean())*t_P(); 

                auto SmD1= Maybe_r_P_cov_t1_y1.value() - diag(Maybe_r_P_mean_t1_y1.value());
                // Reduced start/end cross-covariance after conditioning on y_{0->t}.
                auto r_P_cross_cov_0t_y1= SmD1*t_P()+diag(Maybe_r_P_mean_t1_y1.value())*t_P(); 

                

                Matrix<double> uT(1UL,p_P_mean().size(), 1.0);
         
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y0).value(), Maybe_r_P_mean_t2_y0.value()));
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y1).value(), r_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y0)).value(), p_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y1)).value(), Maybe_r_P_mean_t1_y1.value()));
                get<P_mean_0t_y0>(out())()= std::move(r_P_mean_0t_y0);
                get<P_mean_0t_y1>(out())()= std::move(r_P_mean_0t_y1);
                get<P_cross_cov_0t_y0>(out())()= std::move(r_P_cross_cov_0t_y0);
                get<P_cross_cov_0t_y1>(out())()= std::move(r_P_cross_cov_0t_y1);  


                get<P_Cov_t2_y0>(out())() = std::move(Maybe_r_P_cov_t2_y0.value());    
                get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                get<P_mean_t2_y0>(out())() = std::move(Maybe_r_P_mean_t2_y0.value());
                get<P_mean_t1_y1>(out())() = std::move(Maybe_r_P_mean_t1_y1.value());

                get<d_gS>(out())() = std::move(gS);


                get<P_Cov_t1_y1>(out())() = std::move(Maybe_r_P_cov_t1_y1.value());
            }
            else{
                static_assert(averaging::value == 0);
                get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                auto gS0 = TranspMult(t_gmean_i(), p_P_Cov());
                // gS0 is the start-side probability displacement: μ + chi·gS0 must
                // remain a probability, so gS0·u = 0.
                if (auto Maybe_gS0_check = to_Probability_displacement(gS0); !Maybe_gS0_check)
                    return Maybe_gS0_check.error();

                auto Maybe_r_P_mean_t15_y1 = to_Probability(p_P_mean() + chi * gS0);

                if (!Maybe_r_P_mean_t15_y1) {
                    return Maybe_r_P_mean_t15_y1.error();
                }

                get<P_mean_t15_y1>(out())() = std::move(Maybe_r_P_mean_t15_y1.value());
                auto Maybe_r_P_cov_t15_y1 = to_Covariance_Probability(
                    p_P_Cov() - (N / r_y_var()) * XTX(gS0));  // decoupled: α=1

                if (!Maybe_r_P_cov_t15_y1) {
                    return Maybe_r_P_cov_t15_y1.error();
                }

                get<P_Cov_t15_y1>(out())() = std::move(Maybe_r_P_cov_t15_y1.value());
                get<P_mean_t15_y0>(out())() = std::move(p_P_mean());
                get<P_Cov_t15_y0>(out())() = std::move(p_P_Cov());
                get<d_gS>(out())() = std::move(gS);

                

            }
            return out;
        }
    }

    template <bool dynamic, class recursive, class averaging, class variance,
              class variance_correction, class variance_form, class C_Patch_State,
              class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 uses_variance_form_aproximation_c<variance_form> &&
                 (U<C_Patch_State, Patch_State>))
    auto safely_calculate_Algo_State(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                     C_Patch_Model const& m, C_double const& N,
                                     const Patch_current& p_y, double fs) const {
        if constexpr (!recursive::value) {
            // Non-recursive (NaN/gap) path doesn't apply variance_correction.
            return safely_calculate_Algo_State_non_recursive<dynamic, averaging, variance>(
                t_prior, t_Qdt, m, N, p_y, fs);
        }
        return safely_calculate_Algo_State_recursive<dynamic, averaging, variance,
                                                     variance_correction, variance_form>(
            t_prior, t_Qdt, m, N, p_y, fs);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<Macro_State<vVars...>> update_macro_state(Macro_State<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, logL const& t_logL, elogL const& t_elogL,... ) const 
                                              {
        // Update patch state for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            } else if constexpr (requires { algo.get_P_mean(); algo.get_P_Cov(); }) {
                get<P_mean>(ps())() = algo.get_P_mean();
                get<P_Cov>(ps())() = algo.get_P_Cov();
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        // Optional accumulators when both state and Algo_State provide them.
        if constexpr (has_var_c<Macro_State<vVars...>&, elogL> ) {
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        }
        if constexpr (has_var_c<Macro_State<vVars...>&, vlogL> ) {
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;
        }

        if constexpr (has_var_c<Macro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            // Avoid uninitialized Constant defaults.
            if constexpr (has_var_c<Element&, elogL>) {
                get<elogL>(el)() = t_elogL();   
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                get<vlogL>(el)() = 0.5;
            }

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            // Common predicted fields.
            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<trust_mean_coefficient>{});
            copy_component(std::type_identity<trust_psd_coefficient>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            copy_component(std::type_identity<Chi2>{});

            if constexpr (has_var_c<Element&, Algo_State_Dynamic>) {
                if constexpr (requires { get<Algo_State_Dynamic>(el) = algo; }) {
                    get<Algo_State_Dynamic>(el) = algo;
                }
            }

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    
    template <class... vVars, class C_Algo_State, class C_logL, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<Vector_Space<vVars...>> update_macro_state(Vector_Space<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, C_logL const& t_logL,C_elogL const& t_elogL,...) const {
        // If this Vector_Space carries a Patch_State, keep recursion semantics consistent with
        // Macro_State/dMacro_State/ddMacro_State updates.
        if constexpr (has_var_c<Vector_Space<vVars...>&, Patch_State>) {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            } else if constexpr (requires { algo.get_P_mean(); algo.get_P_Cov(); }) {
                get<P_mean>(ps())() = algo.get_P_mean();
                get<P_Cov>(ps())() = algo.get_P_Cov();
            }
        }

        get<logL>(t_prior_all) ()= get<logL>(t_prior_all)() + t_logL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, elogL>)
            get<elogL>(t_prior_all) ()= get<elogL>(t_prior_all)() + t_elogL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, vlogL>)
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (var::has_it_v<Vector_Space<vVars...>, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all);

            evo.emplace_back(algo);
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_elogL >
    Maybe_error<var::Derivative<Vector_Space<vVars...>, var::Parameters_transformed>>
    update(var::Derivative<Vector_Space<vVars...>, var::Parameters_transformed>&& t_prior_all,
           C_Algo_State&& algo,
           var::Derivative<logL, var::Parameters_transformed> const& t_logL, C_elogL const& t_elogL,...) const {
        get<logL>(t_prior_all) ()= get<logL>(t_prior_all)() + t_logL();
        if constexpr (var::has_it_v<Macro_State<vVars...>, elogL> )
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        if constexpr (var::has_it_v<Macro_State<vVars...>, vlogL> )
            get<vlogL>(t_prior_all) ()= get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (var::has_it_v<Vector_Space<vVars...>, Evolution> &&
                      var::has_it_v<std::decay_t<C_Algo_State>, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all);

            evo.emplace_back(algo);
        }
        return std::move(t_prior_all);
    }


    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<dMacro_State<vVars...>> update_macro_state(
        dMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL ,
        C_elogL const& t_elogL, 
        var::Derivative<y_mean, var::Parameters_transformed> const& t_ymean, 
        var::Derivative<y_var, var::Parameters_transformed> const& t_yvar) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }
        auto d_y_mean= t_ymean.derivative()();
        auto d_y_var= t_yvar.derivative()();

        auto r_y_var=t_yvar.primitive()();
        
        
        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();
       if constexpr (var::has_it_v<dMacro_State<vVars...>, covariance<Grad>> ){
        auto t_CovGradient = covariance<Grad>(
            parameter_spd_payload(XXT(t_logL.derivative()()), var::get_dx_of_dfdx(t_logL)));
        get<covariance<Grad>>(t_prior_all)() = get<covariance<Grad>>(t_prior_all)() + t_CovGradient();
       }

       // NOTE: use `has_var_c` (concept-based, follows inheritance via Var<>
       // subscript) instead of `var::has_it_v` (which is template-pattern based
       // and does NOT see slots through the dMacro_State → Vector_Space derivation).
       // The Hessian/covariance<Grad>/elogL/vlogL siblings using has_it_v in
       // this function are dead branches by the same root cause — they should
       // be migrated to has_var_c in a follow-up cleanup.
       if constexpr (has_var_c<dMacro_State<vVars...>&, Gaussian_Fisher_Information> ){
        // Per-step Gaussian Fisher Information block accumulated across
        // timesteps. G_lik = E[-∂²ℓ/∂θ²] under the moment-matched Gaussian
        // observation model. Used as PSD curvature by MLE optimizers.
        auto t_GFI = parameter_spd_payload(
            XXT(d_y_mean) / r_y_var + XXT(d_y_var) / (2 * r_y_var * r_y_var),
            var::get_dx_of_dfdx(t_logL));
         auto& current_GFI = get<Gaussian_Fisher_Information>(t_prior_all)();
         // First step: slot is default-constructed (empty SymPosDef); initialise
         // to t_GFI directly. Subsequent steps accumulate.
         if (current_GFI.value().nrows() == 0) {
             current_GFI = t_GFI;
         } else {
             current_GFI = current_GFI + t_GFI;
         }
       }
       if constexpr (var::has_it_v<dMacro_State<vVars...>, elogL> )
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL ();
        if constexpr (var::has_it_v<dMacro_State<vVars...>, vlogL> )
            get<vlogL>(t_prior_all) ()= get<vlogL>(t_prior_all)() + 0.5;

        if constexpr (has_var_c<dMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            if constexpr (has_var_c<Element&, elogL>) {
                get<elogL>(el) = t_elogL;
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                if constexpr (requires { get<vlogL>(el)() = 0.5; }) {
                    get<vlogL>(el)() = 0.5;
                } else {
                    get<vlogL>(el) = 0.5;
                }
            }



            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<trust_mean_coefficient>{});
            copy_component(std::type_identity<trust_psd_coefficient>{});
            copy_component(std::type_identity<r_std>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update_macro_state(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,...) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        if constexpr (has_var_c<ddMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            auto seed_zero = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id>) {
                    using Comp = std::decay_t<decltype(get<Id>(el))>;
                    if constexpr (var::is_derivative_v<Comp>) {
                        if constexpr (std::constructible_from<Comp, decltype(t_logL.dx()) const&>) {
                            get<Id>(el) = Comp(t_logL.dx());
                        }
                    } else if constexpr (requires { get<Id>(el)() = 0.0; }) {
                        get<Id>(el)() = 0.0;
                    }
                }
            };

            seed_zero(std::type_identity<elogL>{});
            seed_zero(std::type_identity<vlogL>{});

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<trust_mean_coefficient>{});
            copy_component(std::type_identity<trust_psd_coefficient>{});

            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update_macro_state(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,
        C_elogL const& t_elogL,...) const {
        // Update patch state (including derivatives) for recursion.
        {
            auto& ps = get<Patch_State>(t_prior_all);
            if constexpr (has_var_c<decltype(algo()) const&, P_mean> &&
                          has_var_c<decltype(algo()) const&, P_Cov>) {
                get<P_mean>(ps()) = get<P_mean>(algo());
                get<P_Cov>(ps()) = get<P_Cov>(algo());
            }
        }

        get<logL>(t_prior_all)() = get<logL>(t_prior_all)() + t_logL();

        if constexpr (var::has_it_v<ddMacro_State<vVars...>, elogL>) {
            get<elogL>(t_prior_all)() = get<elogL>(t_prior_all)() + t_elogL();
        }
        if constexpr (var::has_it_v<ddMacro_State<vVars...>, vlogL>) {
            get<vlogL>(t_prior_all)() = get<vlogL>(t_prior_all)() + 0.5;
        }

        if constexpr (has_var_c<ddMacro_State<vVars...>&, Evolution>) {
            auto& evo = get<Evolution>(t_prior_all)();
            using Evo = std::decay_t<decltype(get<Evolution>(t_prior_all))>;
            using Element = typename Evo::element_type;

            Element el{};

            if constexpr (has_var_c<Element&, logL>) {
                get<logL>(el) = t_logL;
            }
            if constexpr (has_var_c<Element&, elogL>) {
                if constexpr (requires { get<elogL>(el) = t_elogL; }) {
                    get<elogL>(el) = t_elogL;
                } else if constexpr (requires { get<elogL>(el) = var::primitive(t_elogL); }) {
                    get<elogL>(el) = var::primitive(t_elogL);
                } else if constexpr (requires { get<elogL>(el)() = var::primitive(t_elogL)(); }) {
                    get<elogL>(el)() = var::primitive(t_elogL)();
                }
            }
            if constexpr (has_var_c<Element&, vlogL>) {
                if constexpr (requires { get<vlogL>(el)() = 0.5; }) {
                    get<vlogL>(el)() = 0.5;
                }
            }

            auto copy_component = [&](auto type_tag) {
                using Id = typename decltype(type_tag)::type;
                if constexpr (has_var_c<Element&, Id> && has_var_c<decltype(algo()) const&, Id>) {
                    if constexpr (requires { get<Id>(el) = get<Id>(algo()); }) {
                        get<Id>(el) = get<Id>(algo());
                    } else if constexpr (requires { get<Id>(el) = var::primitive(get<Id>(algo())); }) {
                        get<Id>(el) = var::primitive(get<Id>(algo()));
                    }
                }
            };

            copy_component(std::type_identity<y_mean>{});
            copy_component(std::type_identity<y_var>{});
            copy_component(std::type_identity<trust_coefficient>{});
            copy_component(std::type_identity<trust_mean_coefficient>{});
            copy_component(std::type_identity<trust_psd_coefficient>{});
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<r_std>{});
            
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class recursive, class averaging, class variance, class variance_correction,
              class variance_form, class FunctionTable, class C_Macro_State, class C_Qdt, class C_Patch_Model,
              class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 uses_variance_form_aproximation_c<variance_form> &&
                 /*(U<std::decay_t<C_Patch_State>,
                                                         Patch_State>||U<std::decay_t<C_Patch_State>,
                                                         Patch_State_and_Evolution>
               )&& U<C_Patch_Model, Patch_Model>
                                                   &&*/
                 U<C_double, double> && (U<C_Qdt, Qdt> || U<C_Qdt, Qdtm> || U<C_Qdt, Qdtg>))

        Maybe_error<C_Macro_State> Macror(FunctionTable&, C_Macro_State&& t_prior_all,
                                        C_Qdt const& t_Qdt, C_Patch_Model const& m,
                                        C_double const& Nch, const Patch_current& p_y,
                                        double fs) const {
        using Transf = transformation_type_t<C_Qdt>;

        auto& t_prior = get<Patch_State>(t_prior_all);

        auto Maybe_Algo =
            safely_calculate_Algo_State<is_Algo_dynamic<C_Macro_State>(), recursive, averaging,
                                        variance, variance_correction, variance_form>(
                t_prior, t_Qdt, m, Nch, p_y, fs);
        if (!Maybe_Algo)
            return Maybe_Algo.error();

        auto r_Algo_state = std::move(Maybe_Algo.value());
        auto r_y_mean = get<y_mean>(r_Algo_state());
        auto r_y_var = get<y_var>(r_Algo_state());
        auto r_r_std = get<r_std>(r_Algo_state());
        
        auto r_chi2 = get<Chi2>(r_Algo_state());
        auto y = p_y.value();
        bool y_is_nan = std::isnan(y);
        auto r_logL = calculate_logL(y_is_nan, r_y_var, r_chi2, m);
        auto r_elogL = calculate_elogL(y_is_nan, r_y_var, m);
        
       
        auto r_prior_all =
            update_macro_state(std::move(t_prior_all), std::move(r_Algo_state), std::move(r_logL), std::move(r_elogL), r_y_mean, r_y_var);
        return {std::move(r_prior_all)};
    }

    // template <class recursive, class averaging, class variance, class variance_correction,
    //           class FunctionTable, class C_Patch_State, class C_Qdt, class C_Patch_Model,
    //           class C_double>

    //     requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
    //              uses_variance_aproximation_c<variance> &&
    //              uses_taylor_variance_correction_aproximation_c<variance_correction> &&
    //              /*(U<std::decay_t<C_Patch_State>,
    //                                                          Patch_State>||U<std::decay_t<C_Patch_State>,
    //                                                          Patch_State_and_Evolution>
    //                )&& U<C_Patch_Model, Patch_Model>
    //                                                    &&*/
    //              U<C_double, double> && (U<C_Qdt, Qdt> || U<C_Qdt, Qdtm> || U<C_Qdt, Qdtg>))

    // Maybe_error<C_Patch_State> Macror_old(FunctionTable&, C_Patch_State&& t_prior, C_Qdt const& t_Qdt,
    //                                       C_Patch_Model const& m, C_double const& Nch,
    //                                       const Patch_current& p_y, double fs) const {
    //     get<macror_algorithm>(t_prior)() =
    //         ToString(MacroR2<recursive, averaging, variance, variance_correction>{});
    //     using Transf = transformation_type_t<C_Qdt>;

    //     auto& p_P_mean = get<P_mean>(t_prior);
    //     auto SmD = get<P_Cov>(t_prior)() - diag(p_P_mean());
    //     auto& y = p_y.value();
    //     auto& t_tolerance = get<Probability_error_tolerance>(m);
    //     auto& t_min_P = get<min_P>(m);
    //     auto y_baseline = get<Current_Baseline>(m);
    //     auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value() +
    //              get<Pink_Noise>(m).value();

    //     auto N = Nch;
    //     Matrix<double> u(p_P_mean().size(), 1, 1.0);

    //     auto N_states = p_P_mean().ncols();

    //     Op_t<Transf, double> ms = 0;
    //     if constexpr (variance::value)
    //         ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

    //     auto& t_gmean_i = get<gmean_i>(t_Qdt);
    //     auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    //     auto& t_gmean_ij = get<gmean_ij>(t_Qdt);
    //     auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
    //                getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    //     Op_t<Transf, y_mean> r_y_mean;
    //     Op_t<Transf, y_var> r_y_var;

    //     Op_t<Transf, double> sSg;
    //     auto t_P = get<P>(t_Qdt);

    //     r_y_mean = build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());

    //     if (std::isnan(y)) {
    //         get<macror_algorithm>(t_prior)() =
    //             ToString(MacroR2<uses_recursive_aproximation<false>, averaging, variance,
    //                              variance_correction>{});

    //         auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();

    //         auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

    //         auto Maybe_P_cov = to_Covariance_Probability(r_P_cov() + diag(r_P_mean()));
    //         if (!Maybe_P_cov)
    //             return Maybe_P_cov.error();
    //         r_P_cov() = std::move(Maybe_P_cov.value());
    //         if constexpr (U<C_Patch_State, Patch_State_and_Evolution>) {
    //             auto& ev = get<Macro_State_Evolution>(t_prior);
    //             ev().push_back(build<Patch_State>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), r_P_mean, r_P_cov, r_y_mean, r_y_var,
    //                 plogL(NaN), eplogL(NaN), vplogL(NaN), get<macror_algorithm>(t_prior)));
    //             return build<Patch_State_and_Evolution>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), std::move(ev));
    //         } else if constexpr (U<C_Patch_State, Patch_State_and_y_Evolution>) {
    //             auto& yev = get<ymean_Evolution>(t_prior);
    //             yev().push_back(r_y_mean);
    //             auto& yvev = get<yvar_Evolution>(t_prior);
    //             yvev().push_back(r_y_var);

    //             return build<Patch_State_and_y_Evolution>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), std::move(yev), std::move(yvev));
    //         } else if constexpr (var::is_derivative_v<C_Patch_State>) {
    //             return build<Patch_State_and_Hessian>(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior), get<FIM>(t_prior));

    //         } else
    //             return Patch_State(
    //                 build<logL>(get<logL>(t_prior)()), build<elogL>(get<elogL>(t_prior)()),
    //                 build<vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean), std::move(r_P_cov),
    //                 std::move(r_y_mean), std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN),
    //                 get<macror_algorithm>(t_prior));
    //     }

    //     constexpr bool PoissonDif = true;
    //     using std::abs;
    //     if constexpr (PoissonDif)
    //         e = e + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
    //     else
    //         e = e + get<Proportional_Noise>(m).value() * abs(y);

    //     auto r_y_mean_max = max_possible_value_of_ymean(N_Ch_mean_value(primitive(Nch)),
    //                                                     primitive(get<g>(m)), primitive(y_baseline));

    //     auto r_y_mean_min = min_possible_value_of_ymean(N_Ch_mean_value(primitive(Nch)),
    //                                                     primitive(get<g>(m)), primitive(y_baseline));

    //     // if ((primitive(r_y_mean()) - r_y_mean_max()) >
    //     //     std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_max())) * 1e-3)
    //     //     std::cerr << "\n max violation" << r_y_mean() << "  vs  max: " << r_y_mean_max();
    //     // if ((r_y_mean_min() - primitive(r_y_mean())) >
    //     //     std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_min())) *
    //     //         1e-1)
    //     // std::cerr << "\n min violation\n"
    //     //          << r_y_mean() << "  vs  min: " << r_y_mean_min();

    //     if (std::isfinite(primitive(gSg)) && primitive(gSg) > 0) {
    //         if (isfinite(primitive(ms)) && primitive(ms) > 0) {
    //             r_y_var = build<y_var>(e + N * gSg + N * ms);
    //         } else {
    //             r_y_var = build<y_var>(e + N * gSg);
    //         }
    //     } else {
    //         if (isfinite(primitive(ms)) && primitive(ms) > 0)
    //             r_y_var = build<y_var>(e + N * ms);
    //         else
    //             r_y_var = build<y_var>(e);
    //     }

    //     auto dy = y - r_y_mean();
    //     auto chi = dy / r_y_var();
    //     Op_t<Transf, P_mean> r_P_mean;
    //     Op_t<Transf, P_Cov> r_P_cov;

    //     if constexpr (!recursive::value) {
    //         r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));

    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean.valid())
    //             return Maybe_r_P_mean.error();

    //         r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

    //         auto Maybe_r_P_cov = to_Covariance_Probability(r_P_cov() + diag(r_P_mean()));

    //         if (!Maybe_r_P_cov.valid())
    //             return Maybe_r_P_cov.error();
    //         r_P_cov() = std::move(Maybe_r_P_cov.value());

    //     } else if constexpr (!variance_correction::value) {
    //         auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();

    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P() + chi * gS);
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();
    //         r_P_mean() = std::move(Maybe_r_P_mean.value());

    //         auto Maybe_r_P_cov = to_Covariance_Probability(
    //             AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) - (N / r_y_var()) * XTX(gS));
    //         if (!Maybe_r_P_cov)
    //             return Maybe_r_P_cov.error();

    //         r_P_cov() = std::move(Maybe_r_P_cov.value());
    //     } else {
    //         auto& t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
    //         auto& t_gvar_ij = get<gvar_ij>(t_Qdt);
    //         auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    //         auto& t_gvar_i = get<gvar_i>(t_Qdt);
    //         auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    //         auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
    //         auto sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
    //                    getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));

    //         auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    //         auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    //         auto e_mu = e + N * ms0;
    //         r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) - N * 0.5 / e_mu * sSg + y_baseline();
    //         auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    //         r_y_var() = var::max(e, e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
    //         auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
    //         auto sS = TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();
    //         r_P_mean() =
    //             to_Probability(p_P_mean() * t_P() + chi * gS - (chi * zeta * sSg + 0.5 / e_mu) * sS);

    //         r_P_cov() = AT_B_A(t_P(), SmD) + diag(r_P_mean() * t_P()) -
    //                     (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
    //                     (2.0 * N / r_y_var() * zeta * sSg) * X_plus_XT(TranspMult(sS, gS)) -
    //                     (N / r_y_var()) * XTX(gS);
    //     }

    //     if (!all_Probability_elements(primitive(r_P_mean())) ||
    //         !all_Covariance_elements(primitive(r_P_cov()))) {
    //         auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
    //         if (!Maybe_r_P_mean)
    //             return Maybe_r_P_mean.error();

    //         r_P_mean() = Maybe_r_P_mean.value();
    //         auto Maybe_r_P_cov =
    //             to_Covariance_Probability(AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));

    //         if (!Maybe_r_P_cov)
    //             return Maybe_r_P_cov.error();

    //         r_P_cov() = std::move(Maybe_r_P_cov.value());
    //         get<macror_algorithm>(t_prior)() =
    //             ToString(MacroR2<uses_recursive_aproximation<false>, averaging, variance,
    //                              variance_correction>{});
    //     }

    //     auto chi2 = dy * chi;

    //     Op_t<Transf, plogL> r_plogL;
    //     Op_t<Transf, eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5);
    //     if (primitive(r_y_var()) > 0.0) {
    //         if (get<Proportional_Noise>(m).value() == 0) {
    //             r_plogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2;
    //             if constexpr (var::is_derivative_v<std::decay_t<decltype(r_plogL)>>)
    //                 if (std::isnan(var::derivative(r_plogL())()[0]))
    //                     std::cerr << "nan der\n";
    //             r_eplogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5;
    //         } else {
    //             r_plogL() = -log(var::Poisson_noise_normalization(
    //                             primitive(r_y_var()), primitive(get<Proportional_Noise>(m).value()))) -
    //                         0.5 * chi2;
    //             r_eplogL() = var::Poisson_noise_expected_logL(
    //                 primitive(r_y_var()), primitive(get<Proportional_Noise>(m).value()));
    //         }
    //     } else {
    //         std::stringstream ss;
    //         ss << "Negative variance!!\n";
    //         ss << "\nr_y_var=\t" << r_y_var;
    //         ss << "\ngSg=\t" << gSg;
    //         return error_message(ss.str());
    //     }

    //     vplogL r_vlogL(0.5);
    //     if (std::isnan(primitive(r_plogL()))) {
    //         std::stringstream ss;
    //         ss << "likelihood is nan \n patch current=";
    //         print(ss, p_y) << "Qdt";
    //         print(ss, primitive(t_Qdt)) << "tprior";
    //         print(ss, primitive(t_prior));
    //         return error_message(ss.str());
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_Evolution>) {
    //         auto& ev = get<Macro_State_Evolution>(t_prior);
    //         ev().push_back(build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
    //                                           build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //                                           build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), r_P_mean,
    //                                           r_P_cov, r_y_mean, r_y_var, r_plogL, r_eplogL, r_vlogL,
    //                                           get<macror_algorithm>(t_prior)));
    //         return build<Patch_State_and_Evolution>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), std::move(ev));
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_y_Evolution>) {
    //         auto& yev = get<ymean_Evolution>(t_prior);
    //         yev().push_back(r_y_mean);
    //         auto& yvev = get<yvar_Evolution>(t_prior);
    //         yvev().push_back(r_y_var);
    //         return build<Patch_State_and_y_Evolution>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), std::move(yev), std::move(yvev));
    //     } else if constexpr (U<C_Patch_State, Patch_State_and_Hessian>) {
    //         auto r_J = derivative(r_y_mean)();
    //         auto r_JS = derivative(r_y_var)();

    //         auto r_FIM =
    //             XXT(r_J) / primitive(r_y_var()) + XXT(r_JS) / (2.0 * sqr(primitive(r_y_var())));

    //         return build<Patch_State_and_Hessian>(
    //             build<logL>(get<logL>(t_prior)() + r_plogL()),
    //             build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //             build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
    //             std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //             get<macror_algorithm>(t_prior), FIM(get<FIM>(t_prior)() + r_FIM));
    //     } else {
    //         return build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
    //                                   build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
    //                                   build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
    //                                   std::move(r_P_mean), std::move(r_P_cov), std::move(r_y_mean),
    //                                   std::move(r_y_var), r_plogL, r_eplogL, r_vlogL,
    //                                   get<macror_algorithm>(t_prior));
    //     }
    // }

    template <class C_Patch_Model>

    auto init(const C_Patch_Model& m) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Patch_State>> {
        Transfer_Op_to<C_Patch_Model, Patch_State> out;
        auto r_P_mean = build<P_mean>(get<P_initial>(m)());
        // P_Cov is the bare centered covariance of the one-hot channel-state
        // indicator: Cov(X) = E[XXᵀ] − E[X]E[X]ᵀ = diag(P) − Pᵀ·P.
        // Sanity check: deterministic state P = e_i ⇒ Cov(X) = 0.
        auto r_P_cov = build<P_Cov>(diagpos(r_P_mean()) - XTX(r_P_mean.value()));
        auto r_test = test<true>(r_P_mean, r_P_cov, get<Probability_error_tolerance>(m));
        if (!r_test) {
            return error_message("fails at init: " + r_test.error()());
        }
        get<P_mean>(out()) = std::move(r_P_mean);
        get<P_Cov>(out()) = std::move(r_P_cov);
        return out;
    }

    template <class adaptive, class recursive, class averaging, class variance,
              class variance_correction, class variance_form, class MacroState, class FuncTable, class C_Parameters,
              class Model>

        requires(uses_adaptive_aproximation_c<adaptive> &&
                 uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 uses_variance_form_aproximation_c<variance_form> &&
                 // what_to_include_c<predictions> &&
                 is_of_this_template_type_v<FuncTable, FuncMap_St>)
    auto log_Likelihood(FuncTable& f, const Model& model, const C_Parameters& par,
                        const Recording& y, const Experiment& e) -> Maybe_error<MacroState>
    // -> Maybe_error<std::conditional_t<
    //     var::is_derivative_v<C_Parameters>,
    //     std::conditional_t<predictions::value == 2,
    //                        Transfer_Op_to<C_Parameters, Macro_State_Evolution>,
    //                        Vector_Space<logL, elogL, vlogL, Grad, FIM>>,
    //     std::conditional_t<predictions::value == 2, Macro_State_Evolution,
    //                        std::conditional_t<predictions::value == 1, logL_y_yvar,
    //                                           Vector_Space<logL, elogL, vlogL>>>>>
    {
        //    using v_adaptive = ::V<adaptive>;
        //    using v_recursive = ::V<recursive>;
        //    using v_averaging = ::V<averaging>;
        //    using v_variance = ::V<variance>;
        //    using v_variance_correction = ::V<variance_correction>;

        using Transf = transformation_type_t<C_Parameters>;

        auto Maybe_m = model(par);
        if (!is_valid(Maybe_m)) {
            return get_error(Maybe_m);
        }

        auto m = std::move(get_value(Maybe_m));
        auto fs = get<Frequency_of_Sampling>(e).value();
        auto ini = init(m);
        auto f_local = f.create("_lik");
        constexpr bool test_derivative = false;
        // Flip to true to run a Clarke FD vs. analytic derivative check on
        // MacroR2 at each step. The test lifts t_prior's Derivative<logL,Pt>
        // and Derivative<Patch_State,Pt> components into a
        // Derivative<Vector_Space<logL,Patch_State>,Pt>, so test_derivative_clarke
        // can perturb all four MacroR2 inputs (c_prior, Qdtm, m, Nch) uniformly.
        // Off by default since FD calls are expensive on figure_2 sweeps;
        // each block also has a static counter capping firings at max_tests.
        constexpr bool test_macroir_derivative = true;
        if (!ini) {
            return ini.error();
        }

        using DX = var::dx_of_dfdx_t<C_Parameters>;
        if constexpr (!std::is_same_v<DX, var::NoDerivative>) {
            MACRODR_DX_ASSERT(var::has_dx(par) &&
                              "log_Likelihood: derivative parameters missing dx()");
            MACRODR_DX_ASSERT(var::has_dx(ini.value()) &&
                              "log_Likelihood: init(m) result missing dx in derivative mode");
        }

        auto t_macro = MacroState(std::move(ini.value()));

        if constexpr (!std::is_same_v<DX, var::NoDerivative>) {
            MACRODR_DX_ASSERT(var::has_dx(t_macro) &&
                              "log_Likelihood: MacroState constructed without dx in derivative mode");
        }

        auto Maybe_run = fold(
            0UL, y().size(), std::move(t_macro),
            [this, &f_local, &m, fs, &e, &y](MacroState&& t_prior, std::size_t i_step) {
                Agonist_evolution const& t_step =
                    get<Agonist_evolution>(get<Recording_conditions>(e)()[i_step]);

                auto time = get<Time>(get<Recording_conditions>(e)()[i_step])();
                auto time_segment = get<N_Ch_mean_time_segment_duration>(m)();
                auto Nchs = get<N_Ch_mean>(m)();
                std::size_t i_segment =
                    std::min(Nchs.size() - 1.0, std::floor(time / time_segment));
                auto j_segment = std::min(Nchs.size() - 1, i_segment + 1);
                auto r = std::max(1.0, time / time_segment - i_segment);
                auto Nch = Nchs[i_segment] * (1 - r) + r * Nchs[j_segment];
                constexpr bool test_eigen = true;
                constexpr bool test_Qx = true;
                if constexpr (test_Qx) {
                    const auto h = 1e-7;
                    auto test_der_Qx = var::test_derivative_clarke(
                        [this, &t_step](auto l_m) {
                            return calc_Qx(l_m, get<Agonist_concentration>(t_step()[0]));
                        },
                        h, m);
                    if (!test_der_Qx) {
                        std::cerr << "\n -----start Qx error-----\n";

                        std::cerr << test_der_Qx.error()();
                        std::cerr << "\n -----end Qx test error-----\n";
                        return Maybe_error<MacroState>(test_der_Qx.error());
                    }
                }

                if constexpr (test_eigen && false) {
                    const auto h = 1e-7;
                    auto test_der_eigen = var::test_derivative_clarke(
                        [this, &t_step](auto l_m) {
                            return calc_eigen(l_m, get<Agonist_concentration>(t_step()[0]));
                        },
                        h, m);
                    if (!test_der_eigen) {
                        std::cerr << "\n -----start eigen test error-----\n";

                        std::cerr << test_der_eigen.error()();
                        std::cerr << "\n -----end eigen test error-----\n";
                        return Maybe_error<MacroState>(test_der_eigen.error());
                    }
                }
                if constexpr (!adaptive::value) {
                    if constexpr (!variance_correction::value) {
                        // Qdt-flavor by `vc` (vc=false → Qdtm/Qdtg, vc=true → Qdt below).
                        // For av=2 + vc=false (macro_IR), the residual gvar_i is
                        // recomputed in-place at line ~4084 from Qdtm's
                        // gsqr_i/gtotal_ij/gmean_ij — these are algebraically the
                        // same as Qdt's (E3(x,y,0) ≡ E2(x,y)), so we get the correct
                        // residual without paying calc_Qdt's E3 cost.
                        auto Maybe_t_Qdtm = [this, &f_local,&m,&t_step,&fs]()
                        {
                            if constexpr( averaging::value>0){
                                return calc_Qdtm(f_local, m, t_step, fs);}
                             else{
                                return calc_Qdtg(f_local, m, t_step, fs);
                             }
                        }();
                        if (!Maybe_t_Qdtm)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | calc_Qdt(m) | " +
                                Maybe_t_Qdtm.error()()));
                        auto t_Qdtm = std::move(Maybe_t_Qdtm.value());
                        if constexpr (test_derivative) {
                            const auto h = 1e-7;
                            auto f_no_memoi = f_local.to_bare_functions();

                            auto test_der_t_Qdtm = test_derivative_clarke(
                                [this, &t_step, &fs, &f_no_memoi](auto const& l_m)
                                    -> Maybe_error<Transfer_Op_to<
                                        std::decay_t<decltype(l_m)>,
                                        var::Vector_Space<P, gmean_i, gtotal_ij, gsqr_i, gvar_i>>> {
                                    auto maybe_res = calc_Qdtm<StabilizerPolicyDisabled>(
                                        f_no_memoi, l_m, t_step, fs);
                                    if (!maybe_res.valid())
                                        return maybe_res.error();
                                    return select<P, gmean_i, gtotal_ij, gsqr_i, gvar_i>(
                                        std::move(maybe_res.value()));
                                },
                                h, m);
                            if (true && !test_der_t_Qdtm) {
                                std::cerr << test_der_t_Qdtm.error()();
                                auto Maybe_t_Qdtm_no_stab =
                                    calc_Qdtm<StabilizerPolicyDisabled>(f_no_memoi, m, t_step, fs);
                                if (Maybe_t_Qdtm_no_stab) {
                                    auto t_Qdtm_no_stab = std::move(Maybe_t_Qdtm_no_stab.value());
                                    auto f_no_memoi = f_local.to_bare_functions();

                                    auto test_no_stab = var::test_derivative_clarke<true>(
                                        [this, &t_step, &fs, &f_no_memoi](auto const& l_m) {
                                            return calc_Qdtm<StabilizerPolicyDisabled>(
                                                f_no_memoi, l_m, t_step, fs);
                                        },
                                        h, m);
                                    auto test_no_stab_text = var::test_derivative_clarke<true>(
                                        [this, &t_step, &fs, &f_no_memoi](auto const& l_m) {
                                            return calc_Qdtm<StabilizerPolicyDisabled>(
                                                f_no_memoi, l_m, t_step, fs);
                                        },
                                        h, m);
                                    if (!test_no_stab) {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled also fails:\n"
                                                  << test_no_stab.error()();
                                    } else {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled passes Taylor "
                                                  << test_no_stab_text.error()();
                                    }
                                }
                                std::cerr << "\nt_step\n" << t_step;
                                return Maybe_error<MacroState>(test_der_t_Qdtm.error());
                            }
                        }

                        // Clarke derivative test on MacroR2's logL/Patch_State gradients.
                        //
                        // Lifting (Option A): t_prior is a dMacro_State<...> (a
                        // Vector_Space whose components are Derivative<logL,Pt> /
                        // Derivative<Patch_State,Pt>), which is *not* itself a
                        // Derivative<>. test_derivative_clarke needs primitive() and
                        // Taylor_first() over each input; both are missing for
                        // dMacro_State. We rebuild a minimal
                        // Derivative<Vector_Space<logL,Patch_State>,Pt> from t_prior's
                        // two derivative-bearing components — that lifted type already
                        // has the framework's primitive()/Taylor_first() overloads
                        // (variables_derivative.h:179, derivative_fwd.h:200), so all
                        // four args (c_prior, Qdtm, m, Nch) get perturbed uniformly.
                        //
                        // The lambda then has to reshape on the way in and out:
                        //   derivative pass → input is Derivative<...>: reconstruct
                        //     dMacro_State<>(get<DLogL_t>, get<DPatch_t>) for
                        //     MacroR2; back out the two Derivative components.
                        //   primitive/Taylor pass → input is Vector_Space<logL,
                        //     Patch_State>: reconstruct Macro_State<>(get<logL>,
                        //     get<Patch_State>) for MacroR2; back out the two
                        //     primitive components.
                        // The two branches return different Maybe_error<> types — the
                        // lambda is templated on auto inputs, so each instantiation
                        // has its own return type. test_derivative_clarke stays
                        // happy because Y/Yp/Yn are all Vector_Space<logL,Patch_State>
                        // and dY is Derivative<Vector_Space<logL,Patch_State>,Pt>,
                        // which test_clarke_brackets dispatches on
                        // (derivative_test.h:799).
                        //
                        // Static atomic counter caps firings at max_tests since
                        // calc_dlikelihood_predictions runs this on every (sim ×
                        // algo × step).
                        if constexpr (test_macroir_derivative &&
                                      !std::is_same_v<DX, var::NoDerivative>) {
                            static std::atomic<int> test_macroir_count{0};
                            constexpr int max_tests = 200;
                            int prev_count = test_macroir_count.fetch_add(
                                1, std::memory_order_relaxed);
                            if (prev_count == 0) {
                                std::cerr << "[test_macroir_derivative] FIRST FIRE algo="
                                          << ToString(MacroR2<recursive, averaging,
                                                              variance, variance_correction,
                                                              variance_form>{})
                                          << " (R=" << recursive::value
                                          << " av=" << averaging::value
                                          << " V=" << variance::value
                                          << " vc=" << variance_correction::value
                                          << ", budget=" << max_tests << ")\n";
                            }
                            if (prev_count < max_tests) {
                                const auto h = 1e-7;
                                auto f_no_memoi = f_local.to_bare_functions();

                                using DLogL_t =
                                    var::Derivative<logL, var::Parameters_transformed>;
                                using DPatch_t =
                                    var::Derivative<Patch_State,
                                                    var::Parameters_transformed>;
                                using Combined_t =
                                    var::Derivative<var::Vector_Space<logL, Patch_State>,
                                                    var::Parameters_transformed>;

                                auto c_prior_lifted = Combined_t(
                                    DLogL_t(get<DLogL_t>(t_prior)),
                                    DPatch_t(get<DPatch_t>(t_prior)));

                                auto test_der_macroir = test_derivative_clarke<false>(
                                    [this, &fs, &f_no_memoi, &y, i_step](
                                        auto&& l_t_prior, auto const& l_Qdtm,
                                        auto const& l_m, auto const& l_Nch) {
                                        using TPrior =
                                            std::decay_t<decltype(l_t_prior)>;
                                        if constexpr (var::is_derivative_v<TPrior>) {
                                            using DL =
                                                var::Derivative<logL,
                                                                var::Parameters_transformed>;
                                            using DP =
                                                var::Derivative<Patch_State,
                                                                var::Parameters_transformed>;
                                            using Combined =
                                                var::Derivative<var::Vector_Space<logL,
                                                                                  Patch_State>,
                                                                var::Parameters_transformed>;
                                            dMacro_State<> dprior(
                                                DL(get<logL>(l_t_prior)),
                                                DP(get<Patch_State>(l_t_prior)));
                                            auto Maybe_res =
                                                MacroR2<recursive, averaging, variance,
                                                        variance_correction, variance_form>{}(
                                                    f_no_memoi, std::move(dprior),
                                                    l_Qdtm, l_m, l_Nch, y()[i_step], fs);
                                            if (!Maybe_res)
                                                return Maybe_error<Combined>(
                                                    Maybe_res.error());
                                            auto& res = Maybe_res.value();
                                            return Maybe_error<Combined>(
                                                Combined(DL(get<DL>(res)),
                                                         DP(get<DP>(res))));
                                        } else {
                                            using Result =
                                                var::Vector_Space<logL, Patch_State>;
                                            Macro_State<> mprior(
                                                logL(get<logL>(l_t_prior)),
                                                Patch_State(
                                                    get<Patch_State>(l_t_prior)));
                                            auto Maybe_res =
                                                MacroR2<recursive, averaging, variance,
                                                        variance_correction, variance_form>{}(
                                                    f_no_memoi, std::move(mprior),
                                                    l_Qdtm, l_m, l_Nch, y()[i_step], fs);
                                            if (!Maybe_res)
                                                return Maybe_error<Result>(
                                                    Maybe_res.error());
                                            auto& res = Maybe_res.value();
                                            return Maybe_error<Result>(
                                                Result(logL(get<logL>(res)),
                                                       Patch_State(get<Patch_State>(res))));
                                        }
                                    },
                                    h, c_prior_lifted, t_Qdtm, m, Nch);
                                if (!test_der_macroir) {
                                    std::cerr << "[test_macroir_derivative] FAIL"
	                                              << " algo=" << ToString(MacroR2<recursive,
	                                                  averaging, variance, variance_correction,
	                                                  variance_form>{})
                                              << " i_step=" << i_step
                                              << " y=" << y()[i_step]
                                              << " t_step=" << t_step << "\n"
                                              << test_der_macroir.error()() << "\n";
                                }
                            }
                        }

                        auto Maybe_macror = MacroR2<recursive, averaging, variance, variance_correction,
                                                    variance_form>{}(
                            f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);
                        if (!Maybe_macror)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | MacroR2 | " +
                                Maybe_macror.error()()));
                        return Maybe_error<MacroState>(std::move(Maybe_macror.value()));

                    } else {
                        auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                        if (!Maybe_t_Qdt)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | calc_Qdt | " +
                                Maybe_t_Qdt.error()()));
                        auto t_Qdt = std::move(Maybe_t_Qdt.value());

                        if constexpr (test_derivative) {
                            const auto h = 1e-7;
                            auto f_no_memoi = f_local.to_bare_functions();

                            auto test_der_t_Qdt = var::test_derivative_clarke(
                                [this, &t_step, &fs, &f_no_memoi](auto const& l_m,
                                                                  auto const& /*l_Qx*/) {
                                    return calc_Qdt<StabilizerPolicyDisabled>(f_no_memoi, l_m,
                                                                              t_step, fs);
                                },
                                h, m, t_Qdt);
                            if (true && !test_der_t_Qdt) {
                                std::cerr << test_der_t_Qdt.error()();
                                auto Maybe_t_Qdt_no_stab = calc_Qdt(f_no_memoi, m, t_step, fs);
                                if (Maybe_t_Qdt_no_stab) {
                                    auto t_Qdt_no_stab = std::move(Maybe_t_Qdt_no_stab.value());

                                    auto test_no_stab = var::test_derivative_clarke<true>(
                                        [this, &t_step, &fs, &f_no_memoi](auto const& l_m) {
                                            return calc_Qdt<StabilizerPolicyDisabled>(
                                                f_no_memoi, l_m, t_step, fs);
                                        },
                                        h, m);
                                    if (!test_no_stab) {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled also fails:\n"
                                                  << test_no_stab.error()();
                                    } else {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled passes Taylor "
                                                     "test.\n";
                                    }
                                }
                                std::cerr << "\nt_step\n" << t_step;
                                return Maybe_error<MacroState>(test_der_t_Qdt.error());
                            }
                        }

                        // Clarke derivative test on the variance_correction=true path
                        // (e.g. macro_IRT in figure_2). Same lift-and-reshape pattern
                        // as the !variance_correction branch above — see that comment
                        // for the rationale.
                        if constexpr (test_macroir_derivative &&
                                      !std::is_same_v<DX, var::NoDerivative>) {
                            static std::atomic<int> test_macroir_count_vc{0};
                            constexpr int max_tests = 200;
                            int prev_count = test_macroir_count_vc.fetch_add(
                                1, std::memory_order_relaxed);
                            if (prev_count == 0) {
                                std::cerr << "[test_macroir_derivative] FIRST FIRE algo="
                                          << ToString(MacroR2<recursive, averaging,
                                                              variance, variance_correction,
                                                              variance_form>{})
                                          << " (R=" << recursive::value
                                          << " av=" << averaging::value
                                          << " V=" << variance::value
                                          << " vc=" << variance_correction::value
                                          << ", budget=" << max_tests << ")\n";
                            }
                            if (prev_count < max_tests) {
                                const auto h = 1e-7;
                                auto f_no_memoi = f_local.to_bare_functions();

                                using DLogL_t =
                                    var::Derivative<logL, var::Parameters_transformed>;
                                using DPatch_t =
                                    var::Derivative<Patch_State,
                                                    var::Parameters_transformed>;
                                using Combined_t =
                                    var::Derivative<var::Vector_Space<logL, Patch_State>,
                                                    var::Parameters_transformed>;

                                auto c_prior_lifted = Combined_t(
                                    DLogL_t(get<DLogL_t>(t_prior)),
                                    DPatch_t(get<DPatch_t>(t_prior)));

                                auto test_der_macroir = test_derivative_clarke<false>(
                                    [this, &fs, &f_no_memoi, &y, i_step](
                                        auto&& l_t_prior, auto const& l_Qdt,
                                        auto const& l_m, auto const& l_Nch) {
                                        using TPrior =
                                            std::decay_t<decltype(l_t_prior)>;
                                        if constexpr (var::is_derivative_v<TPrior>) {
                                            using DL =
                                                var::Derivative<logL,
                                                                var::Parameters_transformed>;
                                            using DP =
                                                var::Derivative<Patch_State,
                                                                var::Parameters_transformed>;
                                            using Combined =
                                                var::Derivative<var::Vector_Space<logL,
                                                                                  Patch_State>,
                                                                var::Parameters_transformed>;
                                            dMacro_State<> dprior(
                                                DL(get<logL>(l_t_prior)),
                                                DP(get<Patch_State>(l_t_prior)));
                                            auto Maybe_res =
                                                MacroR2<recursive, averaging, variance,
                                                        variance_correction, variance_form>{}(
                                                    f_no_memoi, std::move(dprior),
                                                    l_Qdt, l_m, l_Nch, y()[i_step], fs);
                                            if (!Maybe_res)
                                                return Maybe_error<Combined>(
                                                    Maybe_res.error());
                                            auto& res = Maybe_res.value();
                                            return Maybe_error<Combined>(
                                                Combined(DL(get<DL>(res)),
                                                         DP(get<DP>(res))));
                                        } else {
                                            using Result =
                                                var::Vector_Space<logL, Patch_State>;
                                            Macro_State<> mprior(
                                                logL(get<logL>(l_t_prior)),
                                                Patch_State(
                                                    get<Patch_State>(l_t_prior)));
                                            auto Maybe_res =
                                                MacroR2<recursive, averaging, variance,
                                                        variance_correction, variance_form>{}(
                                                    f_no_memoi, std::move(mprior),
                                                    l_Qdt, l_m, l_Nch, y()[i_step], fs);
                                            if (!Maybe_res)
                                                return Maybe_error<Result>(
                                                    Maybe_res.error());
                                            auto& res = Maybe_res.value();
                                            return Maybe_error<Result>(
                                                Result(logL(get<logL>(res)),
                                                       Patch_State(get<Patch_State>(res))));
                                        }
                                    },
                                    h, c_prior_lifted, t_Qdt, m, Nch);
                                if (!test_der_macroir) {
                                    std::cerr << "[test_macroir_derivative] FAIL"
	                                              << " algo=" << ToString(MacroR2<recursive,
	                                                  averaging, variance, variance_correction,
	                                                  variance_form>{})
                                              << " i_step=" << i_step
                                              << " y=" << y()[i_step]
                                              << " t_step=" << t_step << "\n"
                                              << test_der_macroir.error()() << "\n";
                                }
                            }
                        }

                        auto Maybe_macror = MacroR2<recursive, averaging, variance, variance_correction,
                                                    variance_form>{}(
                            f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                        if (!Maybe_macror)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | MacroR2 | " +
                                Maybe_macror.error()()));
                        return Maybe_error<MacroState>(std::move(Maybe_macror.value()));
                    }
                } else {
                    if constexpr (!variance_correction::value) {
                        // Same vc-branched Qdt-flavor rule as the !adaptive branch.
                        // av=2 residual gvar_i is recomputed in-place from Qdtm.
                        auto Maybe_t_Qdtm = [this,& f_local, &m, &t_step, &fs]() {
                           if constexpr(averaging::value>0){
                            return calc_Qdtm(f_local, m, t_step, fs);}
                           else{
                            return calc_Qdtg(f_local, m, t_step, fs);}
                        }();
                        if (!Maybe_t_Qdtm)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | calc_Qdt(m) | " +
                                Maybe_t_Qdtm.error()()));
                        auto t_Qdtm = std::move(Maybe_t_Qdtm.value());

                        if constexpr (test_derivative && averaging::value>0 ) {
                            const auto h = 1e-7;
                            auto f_no_memoi = f_local.to_bare_functions();

                            auto test_der_t_Qdtm = var::test_derivative_clarke(
                                [this, &t_step, &fs, &f_no_memoi](auto const& l_m) {
                                    return calc_T_Qdtm(f_no_memoi, l_m, t_step, fs);
                                },
                                h, m, t_Qdtm);
                            if (true && !test_der_t_Qdtm) {
                                std::cerr << test_der_t_Qdtm.error()();
                                auto Maybe_t_Qdtm_no_stab =
                                    calc_Qdtm<StabilizerPolicyDisabled>(f_no_memoi, m, t_step, fs);
                                if (Maybe_t_Qdtm_no_stab) {
                                    auto t_Qdtm_no_stab = std::move(Maybe_t_Qdtm_no_stab.value());
                                    auto test_no_stab = var::test_derivative_clarke(
                                        [this, &t_step, &fs, &f_no_memoi](auto const& l_m) {
                                            return calc_Qdtm<StabilizerPolicyDisabled>(
                                                f_no_memoi, l_m, t_step, fs);
                                        },
                                        h, m);
                                    if (!test_no_stab) {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled also fails:\n"
                                                  << test_no_stab.error()();
                                    } else {
                                        std::cerr << "\n[diagnostic] "
                                                     "StabilizersDisabled passes Taylor "
                                                     "test.\n";
                                    }
                                }
                                std::cerr << "\nt_step\n" << t_step;
                                return Maybe_error<MacroState>(test_der_t_Qdtm.error());
                            }
                        }

                        auto& t_gmean_i = [&t_Qdtm]() ->auto &{
                            if constexpr (averaging::value > 0) {
                                return get<gmean_i>(t_Qdtm);
                            } else {
                                return get<g>(t_Qdtm);
                            }
                        }();
                        double mg = getvalue(primitive(get<P_mean>(get<Patch_State>(t_prior)())()) *
                                             primitive(t_gmean_i)());
                        double g_max = var::max(primitive(t_gmean_i)());
                        double g_min = var::min(primitive(t_gmean_i)());
                        double g_range = g_max - g_min;
                        auto N = primitive(Nch);
                        auto p_bi = (g_max - mg) / g_range;
                        auto q_bi = (mg - g_min) / g_range;
                        bool test_Binomial = is_Binomial_Approximation_valid(
                            N, p_bi, q_bi, get<Binomial_magical_number>(m)());
                        if (test_Binomial) {
                            // using egsr=typename decltype(f.f(MacroR<recursive,
                            // averaging, variance>{}))::ege;
                            //  auto r=egsr();
                            // return f_local.f(
                            //     MacroR2<v_recursive, v_averaging, v_variance,
                            //             v_variance_correction>{},
                            //     std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                            auto Maybe_macror = MacroR2<recursive, averaging, variance, variance_correction,
                                                        variance_form>{}(
                                f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);
                            if (!Maybe_macror)
                                return Maybe_error<MacroState>(error_message(
                                    "k=" + std::to_string(i_step) + " | MacroR2 (Binomial) | " +
                                    Maybe_macror.error()()));
                            return Maybe_error<MacroState>(std::move(Maybe_macror.value()));

                        } else {
                            auto Maybe_macror =
                                MacroR2<uses_recursive_aproximation<false>, averaging, variance,
                                        uses_taylor_variance_correction_aproximation<false>,
                                        variance_form>{}(
                                    f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);
                            if (!Maybe_macror)
                                return Maybe_error<MacroState>(error_message(
                                    "k=" + std::to_string(i_step) + " | MacroR2 (non-Binomial fallback) | " +
                                    Maybe_macror.error()()));
                            return Maybe_error<MacroState>(std::move(Maybe_macror.value()));
                        }
                    } else {
                        auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                        if (!Maybe_t_Qdt)
                            return Maybe_error<MacroState>(error_message(
                                "k=" + std::to_string(i_step) + " | calc_Qdt | " +
                                Maybe_t_Qdt.error()()));
                        auto t_Qdt = std::move(Maybe_t_Qdt.value());
                        if constexpr (test_derivative) {
                            const auto dx = 1e-6;
                            const auto eps = 1e-2;

                            auto test_der_t_Qdt = var::test_Derivative(
                                [this, &t_step, &fs, &f_local](auto const& l_m, auto const& l_Qx) {
                                    return calc_Qdt(f_local, l_m, t_step, fs);
                                },
                                dx, eps, m, t_Qdt);
                            if (true && !test_der_t_Qdt) {
                                std::cerr << test_der_t_Qdt.error()();
                                std::cerr << "\nt_step\n" << t_step;
                                return Maybe_error<MacroState>(test_der_t_Qdt.error());
                            }
                        }

                        double mg = getvalue(primitive(get<P_mean>(get<Patch_State>(t_prior)())()) *
                                             primitive(get<gmean_i>(t_Qdt)()));
                        double g_max = var::max(get<gmean_i>(primitive(t_Qdt))());
                        double g_min = var::min(get<gmean_i>(primitive(t_Qdt))());
                        double g_range = g_max - g_min;
                        auto N = primitive(Nch);
                        auto p_bi = (g_max - mg) / g_range;
                        auto q_bi = (mg - g_min) / g_range;
                        bool test_Binomial = is_Binomial_Approximation_valid(
                            N, p_bi, q_bi, get<Binomial_magical_number>(m)());
                        if (test_Binomial) {
                            // using egsr=typename decltype(f.f(MacroR<recursive,
                            // averaging, variance>{}))::ege;
                            //  auto r=egsr();
                            // return f_local.f(
                            //     MacroR2<v_recursive, v_averaging, v_variance,
                            //             v_variance_correction>{},
                            //     std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                            auto Maybe_macror = MacroR2<recursive, averaging, variance, variance_correction,
                                                        variance_form>{}(
                                f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                            if (!Maybe_macror)
                                return Maybe_error<MacroState>(error_message(
                                    "k=" + std::to_string(i_step) + " | MacroR2 (Binomial, taylor-correction) | " +
                                    Maybe_macror.error()()));
                            return Maybe_error<MacroState>(std::move(Maybe_macror.value()));
                        } else {
                            auto Maybe_macror = MacroR2<uses_recursive_aproximation<false>, averaging, variance,
                                           uses_taylor_variance_correction_aproximation<false>,
                                           variance_form>{}(
                                f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                            if (!Maybe_macror)
                                return Maybe_error<MacroState>(error_message(
                                    "k=" + std::to_string(i_step) + " | MacroR2 (non-Binomial fallback, taylor-correction) | " +
                                    Maybe_macror.error()()));
                            return Maybe_error<MacroState>(std::move(Maybe_macror.value()));
                        }
                    }
                }
            });
        f += f_local;
        if (!Maybe_run)
            return Maybe_run.error();
        return std::move(Maybe_run.value());
    }

    // ────────────────────────────────────────────────────────────────────
    // nonlinearsqr_logLikelihood — lean classical nonlinear least-squares
    // (Moffatt & Hume 2007 JGP) driver. Marginalizes the noise scale σ²
    // (Jeffreys); computes ONLY the deterministic mean
    //     μ_i = N·(P_mean_i · gmean_i) + baseline
    // per interval — NO P_Cov, NO y_var, NO gSg, NO Kalman down-date. The
    // within-interval kinetics (calc_Qdt / calc_Qdtg) are shared verbatim with
    // the macro driver; only the reduction differs.
    //
    //   SSE = Σ_{finite i} (y_i − μ_i)²,   n = #finite y_i
    //   logL = lgamma(n/2) − (n/2)·log(π) − (n/2)·log(SSE)   [Jeffreys marginal]
    //   score = ∂logL/∂θ  (carried for free by the AD chain riding SSE)
    //   Gaussian_Fisher_Information = (n / SSE)·Σ_i (∂μ_i/∂θ)(∂μ_i/∂θ)ᵀ
    //
    // AD-trap (user flag): the Fisher is assembled from the FIRST-order AD
    // Jacobian rows with the primitive(SSE) prefactor — NEVER as
    // derivative(derivative(logL)) (the 2nd-order AD Hessian carries residual-
    // curvature + score-outer terms that are wrong off the optimum).
    //
    // Returns the SAME reduced states as the macro driver:
    //   dMacro_State_Hessian_minimal (Derivative<logL> + Gaussian_Fisher_Information)
    //     for the MLE / derivative path, and
    //   Macro_State_reg (logL/elogL/vlogL) for the value path.
    // ALSO produces the per-interval Evolution when MacroState carries one
    // (dMacro_State_Ev_gradient_all → fig 4 cumulative J_T/F_T, fig 5 distortion):
    // see the pass-2 block in the finalize and nonlinearsqr_cpp_spec.md §H. The
    // remaining Evolution-carrying states (predictions / diagnostic / detailed)
    // stay guarded at the visit sites. `adaptive`, `recursive`, `variance` and
    // `variance_correction` are accepted and ignored (the marginalized LSE mean
    // has no recursion, no emission variance — the [VAR-BLOCK] behaviour).
    template <class adaptive, class recursive, class averaging, class variance,
              class variance_correction, class MacroState, class FuncTable, class C_Parameters,
              class Model>

        requires(uses_adaptive_aproximation_c<adaptive> &&
                 uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                 is_of_this_template_type_v<FuncTable, FuncMap_St>)
    auto nonlinearsqr_logLikelihood(FuncTable& f, const Model& model, const C_Parameters& par,
                                    const Recording& y, const Experiment& e)
        -> Maybe_error<MacroState> {
        using DX = var::dx_of_dfdx_t<C_Parameters>;
        constexpr bool is_deriv = !std::is_same_v<DX, var::NoDerivative>;

        auto Maybe_m = model(par);
        if (!is_valid(Maybe_m)) {
            return get_error(Maybe_m);
        }
        auto m = std::move(get_value(Maybe_m));
        auto fs = get<Frequency_of_Sampling>(e).value();

        auto ini = init(m);
        if (!ini) {
            return ini.error();
        }
        auto f_local = f.create("_lik");

        // Threaded patch state. Only P_mean is propagated (the deterministic
        // mean trajectory); P_Cov rides untouched and is discarded downstream —
        // the reduced result states never read it.
        auto t_patch = std::move(ini.value());
        auto const& dx = var::get_dx_of_dfdx(t_patch);

        // Accumulators. SSE rides as a Derivative<> scalar so the AD chain
        // carries the score; Fisher_acc is the plain Gauss-Newton JᵀJ block.
        auto SSE = var::init_with_dx<DX>(0.0, dx);
        SymPosDefMatrix<double> Fisher_acc;
        std::size_t n = 0;

        // Per-interval retention for the Evolution (fig 4 / fig 5). The (n/SSE)
        // prefactor is global — unknown until the fold ends — so pass 1 only
        // RETAINS (μ_t with dμ_t, y_t) and pass 2 fills the Evolution in the
        // finalize. That post-pass is O(T) and does NOT re-run calc_Qdt.
        // NOTE: not gated on is_deriv — the DERIVATIVE path feeds figures 3/4
        // (dMacro_State_Ev_gradient_all) and the VALUE path feeds figure 1
        // (Macro_State_Ev_diagnostic). Both fills live in the finalize below.
        constexpr bool wants_evo = has_var_c<MacroState&, Evolution>;
        using YMeanT = std::decay_t<decltype(build<y_mean>(var::init_with_dx<DX>(0.0, dx)))>;
        std::vector<YMeanT> ev_ymean;
        std::vector<double> ev_y;
        std::vector<Matrix<double>> ev_pmean;  // open-loop prior mean (figure 1 only)
        if constexpr (wants_evo) {
            ev_ymean.reserve(y().size());
            ev_y.reserve(y().size());
            ev_pmean.reserve(y().size());
        }

        for (std::size_t i_step = 0; i_step < y().size(); ++i_step) {
            Agonist_evolution const& t_step =
                get<Agonist_evolution>(get<Recording_conditions>(e)()[i_step]);

            auto time = get<Time>(get<Recording_conditions>(e)()[i_step])();
            auto time_segment = get<N_Ch_mean_time_segment_duration>(m)();
            auto Nchs = get<N_Ch_mean>(m)();
            std::size_t i_segment =
                std::min(Nchs.size() - 1.0, std::floor(time / time_segment));
            auto j_segment = std::min(Nchs.size() - 1, i_segment + 1);
            auto r = std::max(1.0, time / time_segment - i_segment);
            auto Nch = Nchs[i_segment] * (1 - r) + r * Nchs[j_segment];

            auto y_baseline = get<Current_Baseline>(m);
            double yi = y()[i_step].value();

            // Per-interval mean + P_mean propagation, sharing the accumulators.
            // t_P is the endpoint transition (P for av>0, P_half for av=0);
            // t_gmean_i is the interval mean single-channel current (gmean_i for
            // av>0, the bare conductance g for av=0) — mirroring
            // safely_calculate_Algo_State_non_recursive's r_y_mean.
            auto process = [&](auto const& t_P, auto const& t_gmean_i) -> Maybe_error<bool> {
                auto& p_P_mean = get<P_mean>(t_patch());
                auto r_y_mean =
                    build<y_mean>(Nch * getvalue(p_P_mean() * t_gmean_i) + y_baseline());
                // Pass 1: retain every interval (NaN ones included, so the
                // Evolution stays index-aligned with the recording).
                if constexpr (wants_evo) {
                    ev_ymean.push_back(r_y_mean);
                    ev_y.push_back(yi);
                    // the PRIOR (pre-propagation) mean, i.e. the open-loop state that
                    // produced mu_i — figure 1's "prior open probability" row.
                    ev_pmean.push_back(var::primitive(p_P_mean()));
                }
                if (!std::isnan(yi)) {
                    auto dy = yi - r_y_mean();
                    SSE = SSE + dy * dy;
                    if constexpr (is_deriv) {
                        auto blk = XXT(var::derivative(r_y_mean)());
                        if (Fisher_acc.nrows() == 0)
                            Fisher_acc = std::move(blk);
                        else
                            Fisher_acc = Fisher_acc + blk;
                    }
                    ++n;
                }
                // Propagate P_mean even across NaN (pre-agonist) intervals: the
                // state must advance so every post-prefix μ_i is correct.
                auto Maybe_r_P_mean = to_Probability(p_P_mean() * t_P());
                if (!Maybe_r_P_mean.valid())
                    return Maybe_r_P_mean.error();
                p_P_mean() = std::move(Maybe_r_P_mean.value());
                return true;
            };

            if constexpr (averaging::value > 0) {
                auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                if (!Maybe_t_Qdt)
                    return Maybe_error<MacroState>(error_message(
                        "k=" + std::to_string(i_step) + " | calc_Qdt | " +
                        Maybe_t_Qdt.error()()));
                auto t_Qdt = std::move(Maybe_t_Qdt.value());
                auto ok = process(get<P>(t_Qdt), get<gmean_i>(t_Qdt)());
                if (!ok)
                    return Maybe_error<MacroState>(error_message(
                        "k=" + std::to_string(i_step) + " | nonlinearsqr step | " + ok.error()()));
            } else {
                auto Maybe_t_Qdtg = calc_Qdtg(f_local, m, t_step, fs);
                if (!Maybe_t_Qdtg)
                    return Maybe_error<MacroState>(error_message(
                        "k=" + std::to_string(i_step) + " | calc_Qdtg | " +
                        Maybe_t_Qdtg.error()()));
                auto t_Qdtg = std::move(Maybe_t_Qdtg.value());
                auto ok = process(get<P_half>(t_Qdtg), get<g>(m)());
                if (!ok)
                    return Maybe_error<MacroState>(error_message(
                        "k=" + std::to_string(i_step) + " | nonlinearsqr step | " + ok.error()()));
            }
        }
        f += f_local;

        if (n == 0)
            return Maybe_error<MacroState>(
                error_message("nonlinearsqr_logLikelihood: no finite observations"));
        double SSE_v = var::primitive(SSE);
        if (!std::isfinite(SSE_v) || SSE_v <= 0)
            return Maybe_error<MacroState>(error_message(
                "nonlinearsqr_logLikelihood: SSE not finite/positive: " + std::to_string(SSE_v)));

        double n_over_2 = 0.5 * static_cast<double>(n);
        double logL_const = std::lgamma(n_over_2) - n_over_2 * std::log(std::numbers::pi);

        if constexpr (is_deriv) {
            // logL = C − (n/2)·log(SSE): the AD chain riding SSE carries the
            // score = ∂logL/∂θ into the Derivative<logL> the optimizer reads.
            auto C_deriv = var::init_with_dx<DX>(logL_const, dx);
            auto t_dlogL = build<logL>(C_deriv + log(SSE) * (-n_over_2));
            // Gaussian_Fisher_Information = (n / SSE)·JᵀJ, σ̂² = SSE/n (the GN
            // plug-in NLS curvature). The 1/σ̂² prefactor uses primitive(SSE) —
            // NOT differentiated (AD-trap above).
            double sigma2_inv = static_cast<double>(n) / SSE_v;
            Gaussian_Fisher_Information t_GFI{};
            t_GFI() = parameter_spd_payload(sigma2_inv * Fisher_acc,
                                            var::get_dx_of_dfdx(t_dlogL));

            if constexpr (wants_evo) {
                // ---- Pass 2: per-interval Evolution (fig 4 / fig 5) -------------
                // sigma_hat^2 is known only here. Setting y_var == sigma_hat^2 with
                // an identically ZERO derivative makes every EXISTING downstream
                // consumer's  XXT(dmu)/v + XXT(dv)/(2 v^2)  collapse to
                // (n/SSE)*dmu_t dmu_t^T = the LSE per-interval Fisher, and makes the
                // per-interval scores sum to the scalar total. Nothing downstream
                // changes (no new type, no writer/R/battery edit). See
                // nonlinearsqr_cpp_spec.md section H.
                MacroState out{};
                get<logL>(out) = t_dlogL;
                if constexpr (has_var_c<MacroState&, Patch_State>)
                    get<Patch_State>(out) = t_patch;
                if constexpr (has_var_c<MacroState&, Gaussian_Fisher_Information>)
                    get<Gaussian_Fisher_Information>(out) = t_GFI;

                const double sigma2 = SSE_v / static_cast<double>(n);
                const double sigma = std::sqrt(sigma2);
                const double half_log = 0.5 * std::log(2.0 * std::numbers::pi * sigma2);

                auto& evo = get<Evolution>(out)();
                using Evo = std::decay_t<decltype(get<Evolution>(out))>;
                using Element = typename Evo::element_type;
                evo.reserve(ev_ymean.size());
                for (std::size_t i = 0; i < ev_ymean.size(); ++i) {
                    Element el{};
                    const bool finite = !std::isnan(ev_y[i]);

                    // residual as a Derivative; zero (value and derivative) on a
                    // NaN interval, mirroring calculate_logL's hard zero.
                    auto d = var::init_with_dx<DX>(0.0, dx);
                    if (finite)
                        d = ev_y[i] - ev_ymean[i]();

                    // y_mean = mu_t with dmu_t. On a NaN interval the derivative is
                    // ZEROED so sum_t F_t still equals the returned total GFI (the
                    // fold skips NaN intervals in Fisher_acc).
                    if constexpr (has_var_c<Element&, y_mean>) {
                        if (finite)
                            get<y_mean>(el) = ev_ymean[i];
                        else
                            get<y_mean>(el) = build<y_mean>(
                                var::init_with_dx<DX>(var::primitive(ev_ymean[i]()), dx));
                    }
                    // y_var == sigma_hat^2 with EXPLICIT shaped zeros. A default
                    // (0x0) derivative would emit no y_var/derivative rows and
                    // silently empty figure 4's panel instead of erroring.
                    if constexpr (has_var_c<Element&, y_var>)
                        get<y_var>(el) = build<y_var>(var::init_with_dx<DX>(sigma2, dx));
                    // logL_t = -0.5*log(2*pi*s2) - 0.5*d^2/s2. s2 is a plug-in
                    // constant, so the AD chain yields d(logL_t)/dtheta =
                    // (n/SSE)*d_t*dmu_t exactly, and sum_t equals the scalar score.
                    if constexpr (has_var_c<Element&, logL>)
                        get<logL>(el) = build<logL>(
                            var::init_with_dx<DX>(finite ? -half_log : 0.0, dx) +
                            (d * d) * (finite ? -0.5 / sigma2 : 0.0));
                    if constexpr (has_var_c<Element&, elogL>)
                        get<elogL>(el) = build<elogL>(
                            var::init_with_dx<DX>(finite ? (-half_log - 0.5) : 0.0, dx));
                    // NOTE: sum_t r_std_t^2 == n IDENTICALLY for LSE, so any panel
                    // reading r2_std reads "calibrated" by construction.
                    if constexpr (has_var_c<Element&, r_std>)
                        get<r_std>(el) = build<r_std>(d * (1.0 / sigma));
                    if constexpr (has_var_c<Element&, trust_coefficient>) {
                        if constexpr (requires { get<trust_coefficient>(el)() = 1.0; })
                            get<trust_coefficient>(el)() = 1.0;
                    }
                    evo.emplace_back(std::move(el));
                }
                return Maybe_error<MacroState>(std::move(out));
            } else {
                return Maybe_error<MacroState>(
                    MacroState(std::move(t_dlogL), std::move(t_patch), std::move(t_GFI)));
            }
        } else {
            double logL_v = logL_const - n_over_2 * std::log(SSE_v);

            if constexpr (wants_evo) {
                // ---- Value-path Evolution: the figure-1 per-interval diagnostic ----
                // Same substitution as the derivative path, without AD: the per-interval
                // predictive variance IS the global plug-in sigma_hat^2 = SSE/n, so the
                // LSE's band is FLAT where the macro ones breathe with the gating
                // (e + N*gSg). That contrast is the point of putting the LSE in figure 1.
                //
                // The diagnostic element is NESTED: please_include<logL, elogL, vlogL,
                // Algo_State_Dynamic>, so the per-interval quantities go INSIDE the
                // Algo_State_Dynamic space, not as flat slots.
                //
                // The intra-interval Kalman snapshots (P_mean_t2_y*, P_Cov_t*, d_gS,
                // d_GS, ...) are left EMPTY on purpose: the LSE is non-recursive, so
                // figure 1 should show it as "no update (open loop)", exactly as it
                // already marks NR and MNR. The open-loop prior mean goes into
                // P_mean_t20_y1 because that is the slot Algo_State_Dynamic::get_P_mean()
                // falls back to when the conditioned ones are empty.
                MacroState out{};
                get<logL>(out) = logL(logL_v);
                if constexpr (has_var_c<MacroState&, Patch_State>)
                    get<Patch_State>(out) = t_patch;
                if constexpr (has_var_c<MacroState&, elogL>)
                    get<elogL>(out) = elogL(0.0);
                if constexpr (has_var_c<MacroState&, vlogL>)
                    get<vlogL>(out) = vlogL(0.0);

                const double sigma2 = SSE_v / static_cast<double>(n);
                const double sigma = std::sqrt(sigma2);
                const double half_log = 0.5 * std::log(2.0 * std::numbers::pi * sigma2);

                auto& evo = get<Evolution>(out)();
                using Evo = std::decay_t<decltype(get<Evolution>(out))>;
                using Element = typename Evo::element_type;
                evo.reserve(ev_ymean.size());
                for (std::size_t i = 0; i < ev_ymean.size(); ++i) {
                    Element el{};
                    const bool finite = !std::isnan(ev_y[i]);
                    const double mu = var::primitive(ev_ymean[i]());
                    const double d = finite ? (ev_y[i] - mu) : 0.0;

                    // per-interval PROFILED Gaussian term. NOTE it sums to the profiled
                    // logL, NOT to the reported Jeffreys marginal (they differ by a
                    // function of n only). Label the figure-1 cumulative row accordingly.
                    if constexpr (has_var_c<Element&, logL>)
                        get<logL>(el)() = finite ? (-half_log - 0.5 * d * d / sigma2) : 0.0;
                    if constexpr (has_var_c<Element&, elogL>)
                        get<elogL>(el)() = finite ? (-half_log - 0.5) : 0.0;
                    if constexpr (has_var_c<Element&, vlogL>)
                        get<vlogL>(el)() = 0.5;

                    if constexpr (has_var_c<Element&, Algo_State_Dynamic>) {
                        auto& algo = get<Algo_State_Dynamic>(el)();
                        get<y_mean>(algo)() = mu;
                        get<y_var>(algo)() = sigma2;  // FLAT: the LSE homoscedastic band
                        get<r_std>(algo)() = finite ? d / sigma : 0.0;
                        get<Chi2>(algo)() = finite ? d * d / sigma2 : 0.0;
                        if (i < ev_pmean.size())
                            get<P_mean_t20_y1>(algo)() = ev_pmean[i];
                    }
                    evo.emplace_back(std::move(el));
                }
                return Maybe_error<MacroState>(std::move(out));
            } else {
                // NOTE the `else` is load-bearing: `if constexpr` only discards the
                // branch not taken, so a return placed AFTER the block would still be
                // instantiated for Evolution-carrying states, whose Macro_State ctor
                // takes one more argument (the Evolution) than this 4-arg form.
                return Maybe_error<MacroState>(
                    MacroState(logL(logL_v), std::move(t_patch), elogL(0.0), vlogL(0.0)));
            }
        }
    }

    template <typename Simulate_tag, class Policy = StabilizerPolicyEnabled, class Patch_Model>
    Maybe_error<Simulated_Sub_Step_t<Simulate_tag>> sub_sub_sample(
        mt_64i& mt, Simulated_Sub_Step_t<Simulate_tag> const& t_sim_step, const Patch_Model& m,
        const Agonist_step& t_s, std::size_t n_sub_dt, double fs) {
        auto r_sim_step = t_sim_step;
        auto& t_g = get<g>(m);
        auto N = get<N_channel_state>(t_sim_step);
        double ysum = get<y_sum>(t_sim_step)();
        auto sum_samples = get<number_of_samples>(t_sim_step)();

        auto n_samples = get<number_of_samples>(t_s)();
        auto tQx = calc_Qx(m, get<Agonist_concentration>(t_s));

        auto dt = n_samples / fs;
        auto sub_dt = dt / n_sub_dt;

        double sub_sample = 1.0 * n_samples / n_sub_dt;
        auto Maybe_t_P = calc_P<Policy>(m, tQx, sub_dt, get<min_P>(m)());
        if (!Maybe_t_P) {
            return Maybe_t_P.error();
        }
        auto t_P = std::move(Maybe_t_P.value());
        for (std::size_t i = 0; i < n_sub_dt; ++i) {
            N = sample_Multinomial(mt, t_P, N);
            auto y = getvalue(N() * t_g());
            ysum += y * sub_sample;
            if constexpr (var::has_it_v<Simulate_tag, Only_Ch_Curent_Sub_Evolution>) {
                get<Only_Ch_Curent_Sub_Evolution>(r_sim_step)().emplace_back(y);
            }
            if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Sub_Evolution>) {
                get<N_Ch_State_Sub_Evolution>(r_sim_step)().emplace_back(N);
            }
        }
        sum_samples += n_samples;
        //  std::cerr << N << sum_samples << "  " << ysum << "  "
        //            << ysum / sum_samples << "\n";
        get<N_channel_state>(r_sim_step) = N;
        get<number_of_samples>(r_sim_step)() = sum_samples;
        get<y_sum>(r_sim_step)() = ysum;
        return r_sim_step;
    }

    struct Uniformization_Segment {
        Matrix<double> embedded_transition;
        double lambda = 0.0;
        std::size_t n_samples = 0;
    };

    template <class Patch_Model>
    static std::vector<std::size_t> expand_channel_states(const Patch_Model& m,
                                                          const N_channel_state& counts) {
        auto n_states = get<g>(m)().size();
        std::vector<std::size_t> channel_states;
        std::size_t n_channels = 0;
        for (std::size_t i = 0; i < n_states; ++i) {
            n_channels += static_cast<std::size_t>(std::llround(counts()[i]));
        }
        channel_states.reserve(n_channels);
        for (std::size_t i = 0; i < n_states; ++i) {
            auto n_i = static_cast<std::size_t>(std::llround(counts()[i]));
            channel_states.insert(channel_states.end(), n_i, i);
        }
        return channel_states;
    }

    template <class Patch_Model>
    static N_channel_state collapse_channel_states(const Patch_Model& m,
                                                   const std::vector<std::size_t>& channel_states) {
        auto n_states = get<g>(m)().size();
        N_channel_state out(Matrix<double>(1, n_states, 0.0));
        for (auto state : channel_states) {
            assert(state < n_states);
            out()[state] += 1.0;
        }
        return out;
    }

    template <class Patch_Model>
    static std::vector<double> state_conductances(const Patch_Model& m) {
        auto n_states = get<g>(m)().size();
        std::vector<double> conductances(n_states, 0.0);
        for (std::size_t i = 0; i < n_states; ++i) {
            conductances[i] = primitive(get<g>(m)()[i]);
        }
        return conductances;
    }

    template <class Qx>
    static Uniformization_Segment make_uniformization_segment(const Qx& t_Qx,
                                                              std::size_t n_samples) {
        auto n_states = t_Qx().nrows();
        Uniformization_Segment segment{Matrix<double>(n_states, n_states, 0.0), 0.0, n_samples};

        for (std::size_t i = 0; i < n_states; ++i) {
            segment.lambda = std::max(segment.lambda, -primitive(t_Qx()(i, i)));
        }

        if (segment.lambda <= 0.0) {
            return segment;
        }

        for (std::size_t i = 0; i < n_states; ++i) {
            double row_sum = 0.0;
            for (std::size_t j = 0; j < n_states; ++j) {
                if (i == j) {
                    continue;
                }
                auto p = std::max(0.0, primitive(t_Qx()(i, j)) / segment.lambda);
                segment.embedded_transition(i, j) = p;
                row_sum += p;
            }

            auto diagonal = std::max(0.0, 1.0 + primitive(t_Qx()(i, i)) / segment.lambda);
            segment.embedded_transition(i, i) = diagonal;
            row_sum += diagonal;

            if (row_sum > 0.0 && std::abs(row_sum - 1.0) > 1e-12) {
                for (std::size_t j = 0; j < n_states; ++j) {
                    segment.embedded_transition(i, j) /= row_sum;
                }
            }
        }

        return segment;
    }

    static std::size_t sample_embedded_transition(mt_64i& mt, const Matrix<double>& transition,
                                                  std::size_t state) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const auto u = uniform(mt);
        double cumulative = 0.0;
        auto fallback = state;
        for (std::size_t j = 0; j < transition.ncols(); ++j) {
            cumulative += transition(state, j);
            if (u <= cumulative || j + 1 == transition.ncols()) {
                return j;
            }
        }
        return fallback;
    }

    static double sample_uniformization_segment(mt_64i& mt, std::size_t& channel_state,
                                                const std::vector<double>& conductances,
                                                const Uniformization_Segment& segment, double fs) {
        if (segment.n_samples == 0) {
            return 0.0;
        }

        const auto duration = segment.n_samples / fs;
        if (segment.lambda <= 0.0) {
            return conductances[channel_state] * segment.n_samples;
        }

        std::exponential_distribution<double> waiting_time(segment.lambda);
        double elapsed = 0.0;
        double integrated_current = 0.0;

        while (elapsed < duration) {
            auto wait = waiting_time(mt);
            auto next_elapsed = elapsed + wait;
            if (next_elapsed >= duration) {
                integrated_current += conductances[channel_state] * (duration - elapsed) * fs;
                break;
            }

            integrated_current += conductances[channel_state] * wait * fs;
            elapsed = next_elapsed;
            channel_state =
                sample_embedded_transition(mt, segment.embedded_transition, channel_state);
        }

        return integrated_current;
    }

    template <typename Simulate_tag, class Patch_Model>
    Maybe_error<Simulated_Sub_Step_t<Simulate_tag>> sub_sub_sample(
        mt_64i& mt, Simulated_Sub_Step_t<Simulate_tag>&& t_sim_step, const Patch_Model& m,
        const std::vector<Agonist_step>& t_s, std::size_t n_sub_dt, double fs) {
        for (std::size_t i = 0; i < t_s.size(); ++i) {
            auto Maybe_sub_step =
                sub_sub_sample<Simulate_tag>(mt, std::move(t_sim_step), m, t_s[i], n_sub_dt, fs);
            if (!Maybe_sub_step)
                return Maybe_sub_step.error();
            else
                t_sim_step = std::move(Maybe_sub_step.value());
        }
        return t_sim_step;
    }

    template <typename Simulate_tag, class Patch_Model>
    Maybe_error<Simulated_Step<Simulate_tag>> uniformization_sample(
        mt_64i& mt, Simulated_Step<Simulate_tag>&& t_sim_step, const Patch_Model& m,
        const Agonist_evolution& t_s, double fs) {
        auto channel_states = expand_channel_states(m, get<N_channel_state>(t_sim_step()));
        auto conductances = state_conductances(m);
        double ysum = 0.0;
        std::size_t total_samples = 0;

        for (auto const& segment_step : t_s()) {
            const auto segment_samples = get<number_of_samples>(segment_step)();
            const auto tQx = calc_Qx(m, get<Agonist_concentration>(segment_step));
            const auto segment = make_uniformization_segment(tQx, segment_samples);
            total_samples += segment_samples;
            for (auto& channel_state : channel_states) {
                ysum +=
                    sample_uniformization_segment(mt, channel_state, conductances, segment, fs);
            }
        }

        if (total_samples == 0) {
            return error_message("uniformization simulation requires positive samples per step");
        }

        auto y_mean = ysum / total_samples;
        get<N_channel_state>(t_sim_step()) = collapse_channel_states(m, channel_states);

        auto& t_e_step = get<Recording>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())());
        double e =
            get<Current_Noise>(m)() * fs / total_samples + get<Pink_Noise>(m).value();

        auto y_baseline = get<Current_Baseline>(m);
        auto y = y_mean + y_baseline() + std::normal_distribution<double>()(mt) * std::sqrt(e);
        auto ey = get<Proportional_Noise>(m).value() * std::abs(y);
        if (ey > 0) {
            y = y + std::normal_distribution<double>()(mt) * std::sqrt(ey);
        }

        t_e_step().emplace_back(Patch_current(y));
        if constexpr (var::has_it_v<Simulate_tag, Only_Ch_Curent_Evolution>) {
            get<Only_Ch_Curent_Evolution>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())())()
                .emplace_back(y_mean);
        }
        if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
            get<N_Ch_State_Evolution>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())())()
                .push_back(get<N_channel_state>(t_sim_step()));
        }

        return t_sim_step;
    }

    template <typename Simulate_tag, class Patch_Model>
    Maybe_error<Simulated_Step<Simulate_tag>> sub_sample(mt_64i& mt,
                                                         Simulated_Step<Simulate_tag>&& t_sim_step,
                                                         const Patch_Model& m,
                                                         const Agonist_evolution& t_s,
                                                         std::size_t n_sub_dt, double fs) {
        // auto &N = get<N_channel_state>(t_sim_step);

        // std::cerr<<N();

        auto t_sub_step =
            Simulated_Sub_Step_build<Simulate_tag>(get<N_channel_state>(t_sim_step()));

        auto Maybe_t_sub_step =
            sub_sub_sample<Simulate_tag>(mt, std::move(t_sub_step), m, t_s(), n_sub_dt, fs);

        if (!Maybe_t_sub_step) {
            return Maybe_t_sub_step.error();
        }
        t_sub_step = std::move(Maybe_t_sub_step.value());
        double y_mean = get<y_sum>(t_sub_step)() / get<number_of_samples>(t_sub_step)();
        get<N_channel_state>(t_sim_step()) = get<N_channel_state>(t_sub_step);

        auto& t_e_step = get<Recording>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())());
        double e = get<Current_Noise>(m)() * fs / get<number_of_samples>(t_sub_step)() +
                   get<Pink_Noise>(m).value();

        auto y_baseline = get<Current_Baseline>(m);

        auto y = y_mean + y_baseline() + std::normal_distribution<double>()(mt) * std::sqrt(e);
        auto ey = get<Proportional_Noise>(m).value() * std::abs(y);
        if (ey > 0)
            y = y + std::normal_distribution<double>()(mt) * std::sqrt(ey);
        t_e_step().emplace_back(Patch_current(y));
        if constexpr (var::has_it_v<Simulate_tag, Only_Ch_Curent_Evolution>) {
            get<Only_Ch_Curent_Evolution>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())())()
                .emplace_back(y_mean);
        }
        if constexpr (var::has_it_v<Simulate_tag, Only_Ch_Curent_Sub_Evolution>) {
            auto& c = get<Only_Ch_Curent_Sub_Evolution>(
                get<Simulated_Recording<Simulate_tag>>(t_sim_step())())();
            auto& n = get<Only_Ch_Curent_Sub_Evolution>(t_sub_step);
            c.insert(c.end(), n().begin(), n().end());
        }
        if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Sub_Evolution>) {
            auto& c = get<N_Ch_State_Sub_Evolution>(
                get<Simulated_Recording<Simulate_tag>>(t_sim_step())())();
            auto& n = get<N_Ch_State_Sub_Evolution>(t_sub_step);
            c.insert(c.end(), n().begin(), n().end());
        }

        if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
            get<N_Ch_State_Evolution>(get<Simulated_Recording<Simulate_tag>>(t_sim_step())())()
                .push_back(get<N_channel_state>(t_sim_step()));
        }

        return t_sim_step;
    }

    template <typename Simulate_tag, class Patch_Model>
    Simulated_Step<Simulate_tag> init_sim(mt_64i& mt, const Patch_Model& m, const Experiment& e) {
        auto initial_x = get<initial_agonist_concentration>(e);
        auto v_Qx = calc_Qx(m, initial_x());
        auto r_P_mean = P_mean(get<P_initial>(m)());
        auto N = static_cast<std::size_t>(std::llround(get<N_Ch_mean>(m)()[0]));
        auto sim = Simulated_Recording<Simulate_tag>{};
        get<SeedNumber>(sim())() = mt.initial_seed();
        auto N_state = sample_Multinomial(mt, r_P_mean, N);
        return Simulated_Step<Simulate_tag>(Vector_Space(std::move(N_state), std::move(sim)));
    }

    template <typename Simulate_tag>
    static Simulated_Recording<Simulate_tag> copy_NaNs(Simulated_Recording<Simulate_tag>&& sim,
                                                       const Recording& r) {
        for (std::size_t i = 0; i < size(r()); ++i)
            if (std::isnan(r()[i]()))
                get<Recording>(sim())()[i]() = r()[i]();

        return std::move(sim);
    }

    template <typename Simulate_tag, class Model>

    Maybe_error<Simulated_Recording<Simulate_tag>> sample_(mt_64i& mt, const Model& model,
                                                           const var::Parameters_values& par,
                                                           const Experiment& e,
                                                           const Simulation_Parameters& sim,
                                                           const Recording& r = Recording{}) {
        auto Maybe_m = model(par);
        if (!Maybe_m)
            return Maybe_m.error();
        else {
            auto m = std::move(Maybe_m.value());

            auto n_sub_dt = get<Simulation_n_sub_dt>(sim);
            auto simulation_mode = get<Simulation_Mode>(sim).value();
            auto fs = get<Frequency_of_Sampling>(e).value();

            auto ini = init_sim<Simulate_tag>(mt, m, e);
            Maybe_error<Simulated_Step<Simulate_tag>> run =
                error_message("unknown simulation mode ", simulation_mode);

            if (simulation_mode == simulation_algorithm_substeps_name) {
                if (n_sub_dt() == 0) {
                    return error_message(
                        "number_of_substeps must be greater than zero for substep simulation");
                }

                run = fold(get<Recording_conditions>(e)(), ini,
                           [this, &m, fs, n_sub_dt,
                            &mt](Simulated_Step<Simulate_tag>&& t_sim_step,
                                 Experiment_step const& t_step) {
                               return Maybe_error<Simulated_Step<Simulate_tag>>(sub_sample(
                                   mt, std::move(t_sim_step), m, t_step, n_sub_dt(), fs));
                           });
            } else if (simulation_mode == simulation_algorithm_uniformization_name) {
                if constexpr (var::has_it_v<Simulate_tag, Only_Ch_Curent_Sub_Evolution> ||
                              var::has_it_v<Simulate_tag, N_Ch_State_Sub_Evolution>) {
                    return error_message(
                        "simulate_with_sub_intervals does not support "
                        "simulation_algorithm=\"uniformization\"; use "
                        "number_of_substeps instead");
                } else {
                    run = fold(get<Recording_conditions>(e)(), ini,
                               [this, &m, fs,
                                &mt](Simulated_Step<Simulate_tag>&& t_sim_step,
                                     Experiment_step const& t_step) {
                                   return Maybe_error<Simulated_Step<Simulate_tag>>(
                                       uniformization_sample(mt, std::move(t_sim_step), m, t_step,
                                                             fs));
                               });
                }
            }
            if (!run)
                return run.error();
            else {
                return copy_NaNs(std::move(get<Simulated_Recording<Simulate_tag>>(run.value()())),
                                 r);
            }
        }
    }
    template <class Model>
    Maybe_error<Simulated_Recording<var::please_include<>>> sample(
        mt_64i& mt, const Model& model, const var::Parameters_values& par, const Experiment& e,
        const Simulation_Parameters& sim, const Recording& r = Recording{}) {
        return sample_<var::please_include<>>(mt, model, par, e, sim, r);
    }

    template <class Model>
    Maybe_error<Simulated_Recording<var::please_include<N_Ch_State_Evolution>>> sample_N(
        mt_64i& mt, const Model& model, const var::Parameters_values& par, const Experiment& e,
        const Simulation_Parameters& sim, const Recording& r = Recording{}) {
        return sample_<var::please_include<N_Ch_State_Evolution>>(mt, model, par, e, sim, r);
    }

    template <class Model>
    Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
        sample_sub_y(mt_64i& mt, const Model& model, const var::Parameters_values& par,
                     const Experiment& e, const Simulation_Parameters& sim,
                     const Recording& r = Recording{}) {
        return sample_<var::please_include<Only_Ch_Curent_Sub_Evolution>>(mt, model, par, e, sim,
                                                                          r);
    }
};

// Pull in the unified dMacro_State_Ev_gradient_* aliases now that Macro_DMR is
// fully defined. micro_types.h transitively includes the (deprecated)
// micro_full.h which references Macro_DMR, so the include must come *after*
// Macro_DMR's closing brace above.
}  // namespace macrodr (briefly, so the include re-enters at file scope)
#include "micro_types.h"
namespace macrodr {

template <class recursive, class averaging, class variance, class variance_correction,
          class variance_form = uses_variance_form_aproximation<variance_total>>

    requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
             uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_variance_form_aproximation_c<variance_form>)
struct MacroR2 {
    friend std::string ToString(MacroR2) {
        std::string out = "MacroR";
        if (recursive::value)
            out += "_R";
        else
            out += "_NR";
        if (averaging::value == 2)
            out += "_2";
        else
            out += "__";
        if (variance::value)
            out += "_V";
        else
            out += "_M";
        if (variance_correction::value)
            out += "_V";
        else
            out += "_M";
        if constexpr (variance_form::value == variance_residual)
            out += "_res";

        return out;
    }

    template <class T, class... Ts>
    auto operator()(T&& x, Ts&&... xs) {
        auto m = Macro_DMR{};
        auto l1 = m.Macror<recursive, averaging, variance, variance_correction, variance_form>(
            std::forward<T>(x), std::forward<Ts>(xs)...);

        return std::move(l1);
    }
};

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class Model,
          class qdt_method = uses_qdt_method<0>,
          class variance_form = uses_variance_form_aproximation<variance_total>>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
struct Likelihood_Model_constexpr {
    // Expose the template parameters as nested aliases so std::visit lambdas
    // can reach them via typename decltype(modelLikelihood)::micro_type, etc.
    // family_type is the stored 3-valued selector; micro_type / nonlinearsqr_type
    // are DERIVED bool_constants so the many visit sites reading micro_type::value
    // keep working while family carries the third (nonlinearsqr) branch.
    using adaptive_type = adaptive;
    using recursive_type = recursive;
    using averaging_type = averaging;
    using variance_type = variance;
    using variance_correction_type = variance_correction;
    using family_type = family;
    using micro_type = std::bool_constant<family_type::value == family_micro>;
    using nonlinearsqr_type = std::bool_constant<family_type::value == family_nonlinearsqr>;
    using qdt_method_type = qdt_method;
    using variance_form_type = variance_form;

    Model m;
    Simulation_n_sub_dt n_sub_dt;
    Likelihood_Model_constexpr(const Model& model, Simulation_n_sub_dt n_sub_dt)
        : m{model}, n_sub_dt{n_sub_dt} {}
    Likelihood_Model_constexpr(adaptive, recursive, averaging, variance, variance_correction, family,
                               const Model& model, Simulation_n_sub_dt n_sub_dt,
                               qdt_method = qdt_method{}, variance_form = variance_form{})
        : m{model}, n_sub_dt{n_sub_dt} {}

    template <class Parameter>
    friend void report_model(save_Parameter<Parameter>& s, Likelihood_Model_constexpr const& d) {
        std::ofstream f(s.fname + "_likelihood_model.csv");
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        f << "adaptive: " << adaptive::value << "\n";
        f << "recursive: " << recursive::value << "\n";
        f << "averaging: " << averaging::value << "\n";
        f << "variance: " << variance::value << "\n";
        f << "variance_correction: " << variance_correction::value << "\n";
        f << "family: " << family::value << "\n";
        f << "taylor_qdt: " << qdt_method::value << "\n";
        f << "variance_form: " << variance_form::value << "\n";
        f << "Simulation_n_sub_dt: " << d.n_sub_dt << "\n";
        report_model(s, d.m);
    }
};

template <class adaptive_range, class recursive_range, class averaging_range, class variance_range,
          class taylor_variance_correction_range, class family_range, class Model,
          class qdt_method_range =
              var::constexpr_Var_domain<int, uses_qdt_method, 0>,
          class variance_form_range =
              var::constexpr_Var_domain<int, uses_variance_form_aproximation, variance_total>>
    requires(
        var::is_this_constexpr_Var_domain_c<adaptive_range, bool, uses_adaptive_aproximation> &&
        var::is_this_constexpr_Var_domain_c<recursive_range, bool, uses_recursive_aproximation> &&
        var::is_this_constexpr_Var_domain_c<averaging_range, int, uses_averaging_aproximation> &&
        var::is_this_constexpr_Var_domain_c<variance_range, bool, uses_variance_aproximation> &&
        var::is_this_constexpr_Var_domain_c<taylor_variance_correction_range, bool,
                                            uses_taylor_variance_correction_aproximation> &&
        var::is_this_constexpr_Var_domain_c<family_range, int, uses_family_aproximation> &&
        var::is_this_constexpr_Var_domain_c<qdt_method_range, int, uses_qdt_method> &&
        var::is_this_constexpr_Var_domain_c<variance_form_range, int,
                                            uses_variance_form_aproximation>)
struct Likelihood_Model_regular {
    Model m;
    Simulation_n_sub_dt n_sub_dt;
    uses_adaptive_aproximation_value adaptive;
    uses_recursive_aproximation_value recursive;
    uses_averaging_aproximation_value averaging;
    uses_variance_aproximation_value variance;
    uses_taylor_variance_correction_aproximation_value taylor_variance_correction;
    uses_family_aproximation_value family;
    uses_qdt_method_value qdt_method;
    uses_variance_form_aproximation_value variance_form;

    static constexpr adaptive_range range_adaptive = {};
    static constexpr recursive_range range_recursive = {};
    static constexpr averaging_range range_averaging = {};
    static constexpr variance_range range_variance = {};
    static constexpr taylor_variance_correction_range range_variance_correction = {};
    static constexpr family_range range_family = {};
    static constexpr qdt_method_range range_qdt_method = {};
    static constexpr variance_form_range range_variance_form = {};

    Likelihood_Model_regular(
        const Model& model, Simulation_n_sub_dt n_sub_dt, uses_adaptive_aproximation_value adaptive,
        uses_recursive_aproximation_value recursive, uses_averaging_aproximation_value averaging,
        uses_variance_aproximation_value variance,
        uses_taylor_variance_correction_aproximation_value taylor_variance_correction,
        uses_family_aproximation_value family,
        uses_qdt_method_value qdt_method = uses_qdt_method_value{},
        uses_variance_form_aproximation_value variance_form =
            uses_variance_form_aproximation_value(variance_total))
        : m{model},
          n_sub_dt{n_sub_dt},
          adaptive{adaptive},
          recursive{recursive},
          averaging{averaging},
          variance{variance},
          taylor_variance_correction{taylor_variance_correction},
          family{family},
          qdt_method{qdt_method},
          variance_form{variance_form} {}

    using cartesian = algebra_2<std::tuple, std::variant>::P_constexpr<
        typename adaptive_range::variant_type, typename recursive_range::variant_type,
        typename averaging_range::variant_type, typename variance_range::variant_type,
        typename taylor_variance_correction_range::variant_type,
        typename family_range::variant_type,
        typename qdt_method_range::variant_type,
        typename variance_form_range::variant_type>;

    template <class, class>
    struct Likelihood_Model_variant_impl;
    template <class... adaptive, class... recursive, class... averaging, class... variance,
              class... taylor_variance_correction, class... family, class... qdt_method,
              class... variance_form, class M>

        requires(
            (uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<taylor_variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>) &&
            ...)
    struct Likelihood_Model_variant_impl<
        std::variant<
            std::tuple<adaptive, recursive, averaging, variance, taylor_variance_correction,
                       family, qdt_method, variance_form>...>,
        M> {
        // Likelihood_Model_constexpr template arg order is
        // (adaptive, recursive, averaging, variance, variance_correction, family,
        //  Model, taylor_qdt, variance_form) — Model precedes the optional flags to keep
        // existing 7-arg instantiations backward-compatible.
        using type =
            std::variant<Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                                    taylor_variance_correction, family, M,
                                                    qdt_method, variance_form>...>;
    };

    using Likelihood_Model_variant = typename Likelihood_Model_variant_impl<cartesian, Model>::type;

    Maybe_error<Likelihood_Model_variant> get_variant() const {
        auto car = promote_Maybe_error(std::make_tuple(
            adaptive_range::to_variant(adaptive.value),
            recursive_range::to_variant(recursive.value),
            averaging_range::to_variant(averaging.value),
            variance_range::to_variant(variance.value),
            taylor_variance_correction_range::to_variant(taylor_variance_correction.value),
            family_range::to_variant(family.value),
            qdt_method_range::to_variant(qdt_method.value),
            variance_form_range::to_variant(variance_form.value)));

        if (!car) {
            return car.error();
        }
        return std::apply(
            [this](auto... t) {
                return std::visit(
                    [this](auto... v) -> Likelihood_Model_variant {
                        // Tuple element order: (adaptive, recursive, averaging,
                        // variance, variance_correction, family, taylor_qdt,
                        // variance_form). Likelihood_Model_constexpr expects
                        // the first 6 then Model then the optional flags.
                        auto reorder = [this](auto a, auto r, auto av, auto va, auto vc, auto fa,
                                              auto tq, auto vf) -> Likelihood_Model_variant {
                            using LMC = Likelihood_Model_constexpr<
                                std::decay_t<decltype(a)>, std::decay_t<decltype(r)>,
                                std::decay_t<decltype(av)>, std::decay_t<decltype(va)>,
                                std::decay_t<decltype(vc)>, std::decay_t<decltype(fa)>, Model,
                                std::decay_t<decltype(tq)>, std::decay_t<decltype(vf)>>;
                            return Likelihood_Model_variant{LMC(m, n_sub_dt)};
                        };
                        return reorder(v...);
                    },
                    t...);
            },
            car.value());
    }
};


template<class ...Vs, class...Ws>
auto merge_Maybe_variant( Maybe_error<std::variant<Vs...>>&& v1, Maybe_error<std::variant<Ws...>>&& v2)-> Maybe_error<std::variant<Vs...,Ws...>>   
{
    using merged_variant = std::variant<Vs..., Ws...>;
    if (v1.valid())
        return std::visit(
            [](auto&& x) -> Maybe_error<merged_variant> {
                return merged_variant{std::forward<decltype(x)>(x)};
            },
            std::move(v1).value());
    if (v2.valid())
        return std::visit(
            [](auto&& x) -> Maybe_error<merged_variant> {
                return merged_variant{std::forward<decltype(x)>(x)};
            },
            std::move(v2).value());
    return {error_message(v1.error()(), " ", v2.error()())} ;
}

// struct Likelihood_Model_v_all {
//     using v_uses_adaptive_aproximation =
//         std::variant<::V<uses_adaptive_aproximation<false>>, ::V<uses_adaptive_aproximation<true>>>;

//     using v_uses_recursive_aproximation = std::variant<::V<uses_recursive_aproximation<false>>,
//                                                        ::V<uses_recursive_aproximation>true>>>;

//     using v_uses_averaging_aproximation =
//         std::variant<::V<uses_averaging_aproximation<0>>, ::V<uses_averaging_aproximation<1>>,
//                      ::V<uses_averaging_aproximation<2>>>;

//     using v_uses_variance_aproximation =
//         std::variant<::V<uses_variance_aproximation<false>>, ::V<uses_variance_aproximation<true>>>;

//     using v_uses_variance_correction_aproximation =
//         std::variant<::V<uses_taylor_variance_correction_aproximation<false>>,
//                      ::V<uses_taylor_variance_correction_aproximation<true>>>;

//     template <uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
//               uses_averaging_aproximation averaging, uses_variance_aproximation variance,
//               uses_taylor_variance_correction_aproximation variance_correction, class Model>
//     auto template_op(::V<adaptive>, ::V<recursive>, ::V<averaging>, ::V<variance>,
//                      ::V<variance_correction>, const Model& model,
//                      Simulation_n_sub_dt n_sub_dt) const {
//         return Likelihood_Model<adaptive, recursive, averaging, variance, variance_correction,
//                                 Model>(model, n_sub_dt);
//     }

//     template <class Model>
//     auto variant_op(v_uses_adaptive_aproximation t_adaptive,
//                     v_uses_recursive_aproximation t_recursive,
//                     v_uses_averaging_aproximation t_averaging,
//                     v_uses_variance_aproximation t_variance,
//                     v_uses_variance_correction_aproximation t_var_corr, Model const& model,
//                     Simulation_n_sub_dt n_sub_dt) const {
//         auto tu = std::tuple(t_adaptive, t_recursive, t_averaging, t_variance, t_var_corr);
//         return Apply_variant(
//             [this, &model, n_sub_dt](auto const&... x) {
//                 // using m11=decltype(template_op(x...,model,
//                 // n_sub_dt))::llego_aquiÑ;
//                 //   return m11{};
//                 return this->template_op(x..., model, n_sub_dt);
//             },
//             tu);
//     }

//     template <class Model>
//     auto bool_op(uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
//                  uses_averaging_aproximation averaging, uses_variance_aproximation variance,
//                  uses_taylor_variance_correction_aproximation variance_correction, const Model& model,
//                  Simulation_n_sub_dt n_sub_dt) const {
//         v_uses_adaptive_aproximation t_adaptive;
//         if (adaptive::value)
//             t_adaptive = ::V<uses_adaptive_aproximation<true>>{};
//         else
//             t_adaptive = ::V<uses_adaptive_aproximation<false>>{};

//         v_uses_recursive_aproximation t_recursive;
//         if (recursive::value)
//             t_recursive = ::V<uses_recursive_aproximation<true>>{};
//         else
//             t_recursive = ::V<uses_recursive_aproximation<false>>{};

//         v_uses_averaging_aproximation t_averaging;
//         if (averaging::value == 0)
//             t_averaging = ::V<uses_averaging_aproximation<0>>{};
//         else if (averaging::value == 1)

//             t_averaging = ::V<uses_averaging_aproximation<1>>{};
//         else
//             t_averaging = ::V<uses_averaging_aproximation<2>>{};

//         v_uses_variance_aproximation t_variance;
//         if (variance::value)
//             t_variance = ::V<uses_variance_aproximation<true>>{};
//         else
//             t_variance = ::V<uses_variance_aproximation<false>>{};

//         v_uses_variance_correction_aproximation t_var_corr;
//         if (variance_correction::value)
//             t_var_corr = ::V<uses_taylor_variance_correction_aproximation<true>>{};
//         else
//             t_var_corr = ::V<uses_taylor_variance_correction_aproximation<false>>{};

//         return this->variant_op(t_adaptive, t_recursive, t_averaging, t_variance, t_var_corr, model,
//                                 n_sub_dt);
//     }
// };

// struct Likelihood_Model_v {
//     using v_uses_adaptive_aproximation =
//         std::variant<::V<uses_adaptive_aproximation<false>>, ::V<uses_adaptive_aproximation<true>>>;

//     using v_uses_recursive_aproximation = std::variant<::V<uses_recursive_aproximation<false>>,
//                                                        ::V<uses_recursive_aproximation<true>>>;

//     using v_uses_averaging_aproximation =
//         std::variant<::V<uses_averaging_aproximation<0>>, ::V<uses_averaging_aproximation<1>>,
//                      ::V<uses_averaging_aproximation<2>>>;

//     using v_uses_variance_aproximation = std::variant<::V<uses_variance_aproximation<true>>>;

//     using v_uses_variance_correction_aproximation =
//         std::variant<::V<uses_taylor_variance_correction_aproximation<false>>>;

//     template <uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
//               uses_averaging_aproximation averaging, uses_variance_aproximation variance,
//               uses_taylor_variance_correction_aproximation variance_correction, class Model>
//     auto template_op(::V<adaptive>, ::V<recursive>, ::V<averaging>, ::V<variance>,
//                      ::V<variance_correction>, const Model& model,
//                      Simulation_n_sub_dt n_sub_dt) const {
//         return Likelihood_Model<adaptive, recursive, averaging, variance, variance_correction,
//                                 Model>(model, n_sub_dt);
//     }

//     template <class Model>
//     auto variant_op(v_uses_adaptive_aproximation t_adaptive,
//                     v_uses_recursive_aproximation t_recursive,
//                     v_uses_averaging_aproximation t_averaging,
//                     v_uses_variance_aproximation t_variance,
//                     v_uses_variance_correction_aproximation t_var_corr, Model const& model,
//                     Simulation_n_sub_dt n_sub_dt) const {
//         auto tu = std::tuple(t_adaptive, t_recursive, t_averaging, t_variance, t_var_corr);
//         return Apply_variant(
//             [this, &model, n_sub_dt](auto const&... x) {
//                 // using m11=decltype(template_op(x...,model,
//                 // n_sub_dt))::llego_aquiÑ;
//                 //   return m11{};
//                 return this->template_op(x..., model, n_sub_dt);
//             },
//             tu);
//     }

//     template <class Model>
//     auto bool_op(uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
//                  uses_averaging_aproximation averaging, uses_variance_aproximation variance,
//                  uses_taylor_variance_correction_aproximation variance_correction, const Model& model,
//                  Simulation_n_sub_dt n_sub_dt) const {
//         v_uses_adaptive_aproximation t_adaptive;
//         if (adaptive::value)
//             t_adaptive = ::V<uses_adaptive_aproximation<true>>{};
//         else
//             t_adaptive = ::V<uses_adaptive_aproximation<false>>{};

//         v_uses_recursive_aproximation t_recursive;
//         if (recursive::value)
//             t_recursive = ::V<uses_recursive_aproximation<true>>{};
//         else
//             t_recursive = ::V<uses_recursive_aproximation<false>>{};

//         v_uses_averaging_aproximation t_averaging;
//         if (averaging::value == 0)
//             t_averaging = ::V<uses_averaging_aproximation<0>>{};
//         else if (averaging::value == 0)
//             t_averaging = ::V<uses_averaging_aproximation<1>>{};
//         else
//             t_averaging = ::V<uses_averaging_aproximation<2>>{};

//         v_uses_variance_aproximation t_variance;
//         t_variance = ::V<uses_variance_aproximation<true>>{};

//         v_uses_variance_correction_aproximation t_var_corr;
//         t_var_corr = ::V<uses_taylor_variance_correction_aproximation<false>>{};

//         return this->variant_op(t_adaptive, t_recursive, t_averaging, t_variance, t_var_corr, model,
//                                 n_sub_dt);
//     }
// };

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class Model,
          class family = uses_family_aproximation<family_macro>,
          class variance_form = uses_variance_form_aproximation<variance_total>>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             uses_variance_form_aproximation_c<variance_form>)
auto make_Likelihood_Model(const Model& m, Simulation_n_sub_dt n_sub_dt) {
    return Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                      family, Model, uses_qdt_method<0>, variance_form>(m, n_sub_dt);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Parameters,
          class Recoding, class Experiment>

    requires(uses_adaptive_aproximation_c<adaptive>,
             uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
                     uses_variance_aproximation_c<variance> &&
                     uses_taylor_variance_correction_aproximation_c<variance_correction> &&
                     uses_family_aproximation_c<family> &&
                     var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
                     uses_variance_form_aproximation_c<variance_form>)
Maybe_error<logLs> logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const Recoding& y, const Experiment& e) {
    auto v_logL = Macro_DMR{}
                      .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                                      variance_form, Macro_State_reg>(f, lik.m, p, y, e);
    if (!v_logL)
        return v_logL.error();
    else
        return logLs(get<logL>(v_logL.value()), get<elogL>(v_logL.value()),
                     get<vlogL>(v_logL.value()));
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<dMacro_State_Hessian_minimal> dlogLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        dMacro_State_Hessian_minimal>(f, lik.m, dpp, y, var);
}

// ── nonlinearsqr (classical LSE) entrypoints ──────────────────────────────
// Two overloads under one name (mirroring the macro logLikelihood/dlogLikelihood
// pair) routed from the likelihood.cpp visit arms:
//   value path (@761)   → nonlinearsqr_logLikelihood(f, lik, par_values, ...)  → logLs
//   MLE / derivative (@787, via calculate_mdlikelihood_nonlinearsqr_impl)
//                       → nonlinearsqr_logLikelihood(f, lik, Parameters_transformed, ...)
//                         → dMacro_State_Hessian_minimal {Derivative<logL>, GFI}
// The concrete Parameters_transformed 3rd param + the value overload's
// !same_as<Parameters, Parameters_transformed> constraint keep the two
// unambiguous.
template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Parameters,
          class Recoding, class Experiment>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form> &&
             !std::is_same_v<std::decay_t<Parameters>, var::Parameters_transformed>)
Maybe_error<logLs> nonlinearsqr_logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const Recoding& y, const Experiment& e) {
    auto v_logL = Macro_DMR{}
                      .nonlinearsqr_logLikelihood<adaptive, recursive, averaging, variance,
                                                  variance_correction, Macro_State_reg>(
                          f, lik.m, p, y, e);
    if (!v_logL)
        return v_logL.error();
    else
        return logLs(get<logL>(v_logL.value()), get<elogL>(v_logL.value()),
                     get<vlogL>(v_logL.value()));
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<dMacro_State_Hessian_minimal> nonlinearsqr_logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .nonlinearsqr_logLikelihood<adaptive, recursive, averaging, variance, variance_correction,
                                    dMacro_State_Hessian_minimal>(f, lik.m, dpp, y, var);
}

// Evolution-producing entrypoint (figure 4 cumulative J_T/F_T, figure 5 distortion).
// Same lean fold as above, but MacroState carries an Evolution, which switches on
// the driver's pass-2 post-fill: per-interval y_var == sigma_hat^2 = SSE/n with an
// identically ZERO derivative, so every existing downstream consumer's
// XXT(dmu)/v + XXT(dv)/(2 v^2) collapses to (n/SSE)*dmu dmu^T (the LSE per-interval
// Fisher) and the per-interval scores sum to the scalar total. Nothing downstream
// changes. Distinct NAME (not an overload) because the signature matches the
// dMacro_State_Hessian_minimal one and only the return type differs.
// Routed from likelihood.cpp calculate_mdlikelihood_predictions_visit.
template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<dMacro_State_Ev_gradient_all> nonlinearsqr_dlogLikelihoodPredictions(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .nonlinearsqr_logLikelihood<adaptive, recursive, averaging, variance, variance_correction,
                                    dMacro_State_Ev_gradient_all>(f, lik.m, dpp, y, var);
}

// VALUE-path Evolution entrypoint (figure 1 per-interval diagnostic). Same lean fold,
// MacroState = Macro_State_Ev_diagnostic, which switches on the driver's value-path
// Evolution fill: per interval mu_t, y_var == sigma_hat^2 (FLAT band), r_std, Chi2 and
// the open-loop prior mean. The intra-interval Kalman snapshots stay empty — the LSE is
// non-recursive, so figure 1 shows it as "no update (open loop)" like NR and MNR.
// Takes Parameters_values (the diagnostics visit passes par.to_value()).
template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Parameters,
          class Recoding, class Experiment>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<Macro_State_Ev_diagnostic> nonlinearsqr_logLikelihoodDiagnostic(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const Recoding& y, const Experiment& e) {
    return Macro_DMR{}
        .nonlinearsqr_logLikelihood<adaptive, recursive, averaging, variance, variance_correction,
                                    Macro_State_Ev_diagnostic>(f, lik.m, p, y, e);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<diff_Macro_State_Gradient_Hessian> diff_logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var,
    double delta_par) {
    auto v_p = p.to_value();
    auto Maybe_MacroEv =
        Macro_DMR{}
            .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                            variance_form,
                            Macro_State_Ev_predictions>(f, lik.m, v_p, y, var);
    if (!Maybe_MacroEv)
        return Maybe_MacroEv.error();

    Macro_State_Ev_predictions v_MacroEv = std::move(Maybe_MacroEv.value());
    auto npar = p.size();
    auto ny = y().size();

    Matrix<double> J_y(ny, npar);
    Matrix<double> J_v(ny, npar);
    Matrix<double> G(npar, 1ul);
    DiagPosDetMatrix<double> t_yvar(ny, ny);
    DiagPosDetMatrix<double> t_yvar_2(ny, ny);

    auto& r_Ev = get<Evolution>(v_MacroEv);
    for (std::size_t i = 0; i < ny; ++i) {
        auto r_yvar = get<y_var>(r_Ev()[i])();
        auto r_ymean = get<y_mean>(r_Ev()[i])();
        if (std::isfinite(r_yvar) && std::isfinite(r_ymean) && (r_yvar > 0)) {
            t_yvar[i] = 1.0 / r_yvar;
            t_yvar_2[i] = 1.0 / (2.0 * r_yvar * r_yvar);
        } else {
            t_yvar[i] = 0;
            t_yvar_2[i] = 0;
        }
    }

    for (std::size_t ipar = 0; ipar < npar; ++ipar) {
        auto pi = p;
        pi[ipar] = pi[ipar] + delta_par;
        auto v_pi = pi.to_value();
        auto Maybe_Ev_dif =
            Macro_DMR{}
                .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                                variance_form,
                                Macro_State_Ev_predictions>(f, lik.m, v_pi, y, var);
        if (!Maybe_Ev_dif)
            return Maybe_Ev_dif.error();
        auto v_MacroEvi = std::move(Maybe_Ev_dif.value());
        auto& r_Evi = get<Evolution>(v_MacroEvi);

        G[ipar] = (get<logL>(v_MacroEvi)() - get<logL>(v_MacroEv)()) / delta_par;
        for (std::size_t i = 0; i < ny; ++i) {
            J_y(i, ipar) = (get<y_mean>(r_Evi()[i])() - get<y_mean>(r_Ev()[i])()) / delta_par;
            J_v(i, ipar) = (get<y_var>(r_Evi()[i])() - get<y_var>(r_Ev()[i])()) / delta_par;
        }
    }

    auto r_FIM = AT_D_A_Safe(J_y, t_yvar) + AT_D_A_Safe(J_v, t_yvar_2);

    return diff_Macro_State_Gradient_Hessian(
        std::move(get<logL>(v_MacroEv)), std::move(get<Patch_State>(v_MacroEv)),
        std::move(get<elogL>(v_MacroEv)), std::move(get<vlogL>(v_MacroEv)),
        Grad(std::move(G), p), Gaussian_Fisher_Information(std::move(r_FIM), p));
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<Macro_State_Ev_predictions> logLikelihoodPredictions(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        Macro_State_Ev_predictions>(f, lik.m, p, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<Macro_State_Ev_diagnostic> logLikelihoodDiagnostic(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        Macro_State_Ev_diagnostic>(f, lik.m, p, y, var);
}







template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<dMacro_State_Ev_gradient_all> dlogLikelihoodPredictions(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        dMacro_State_Ev_gradient_all>(f, lik.m, dpp, y, var);
}

// Detailed-state sibling: same evaluation, but the per-step Evolution carries
// the rich detailed_element (P_mean/P_Cov/y_mean/y_var/trust_coefficient/logL
// as Derivatives). Used for the FD-instability localization diagnostic.
template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FuncTable, class Model, class Variables, class DataType>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<dMacro_State_Ev_detailed> dlogLikelihoodPredictionsDetailed(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        dMacro_State_Ev_detailed>(f, lik.m, dpp, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<Macro_State_Ev_diagnostic> logLikelihood_Diagnostic(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction, variance_form,
                        Macro_State_Ev_diagnostic>(f, lik.m, p, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
Maybe_error<std::vector<Macro_State_Ev_predictions>> fractioned_logLikelihoodPredictions(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     family, Model, qdt_method, variance_form>& lik,
    Parameters const& p, const std::vector<Variables>& y, const std::vector<DataType>& var) {
    // std::vector<Macro_State_Evolution> out;
    // for (std::size_t i = 0; i < y.size(); ++i) {
    //     auto Maybe_e = Macro_DMR{}
    //                        .log_Likelihood<adaptive, recursive, averaging, variance,
    //                                        variance_correction, Macro_State_Ev_predictions>(
    //                            f, lik.m, p, get<Recording>(y[i]()), var[i]);
    //     if (Maybe_e)
    //         out.push_back(Maybe_e.value());
    //     else
    //         return Maybe_e.error();
    // }
    // return out;
}

inline auto get_num_samples(const Agonist_step& e) {
    return get<number_of_samples>(e)();
}

inline auto get_num_samples(const std::vector<Agonist_step>& e) {
    double out = 0;
    for (auto& elem : e) out += get_num_samples(elem);
    return out;
}

inline auto get_num_samples(const Agonist_evolution& e) {
    return get_num_samples(e());
}

inline Agonist_evolution average_agonist_step(Agonist_step const& x) {
    return std::vector<Agonist_step>(1, x);
}

inline Agonist_evolution average_agonist_step(std::vector<Agonist_step> const& x) {
    if (x.empty())
        return x;
    else if (x.size() == 1)
        return x;
    else {
        auto out = std::vector<Agonist_step>{x[0]};
        for (std::size_t i = 1; i < x.size(); ++i) {
            if ((get<Agonist_concentration>(out.back())() != get<Agonist_concentration>(x[i])()) &&
                ((get<Agonist_concentration>(out.back())() == 0) ||
                 (get<Agonist_concentration>(x[i])() == 0))) {
                out.push_back(x[i]);
            } else {
                auto n1 = get<number_of_samples>(out.back())();
                auto n2 = get<number_of_samples>(x[i])();
                auto new_agonist = (get<Agonist_concentration>(out.back())() * n1 +
                                    get<Agonist_concentration>(x[i])() * n2) /
                                   (n1 + n2);
                get<number_of_samples>(out.back())() = n1 + n2;
                get<Agonist_concentration>(out.back())() = new_agonist;
            }
        }
        return out;
    }
}

// static Agonist_evolution average_agonist_step(Agonist_evolution const & x){
//     return std::visit([] (auto const& a){return
//     average_agonist_step(a);},x());}

inline Agonist_evolution average_agonist_step(Agonist_evolution const& x,
                                              bool average_the_evolution) {
    if (average_the_evolution)
        return average_agonist_step(x());
    else
        return x;
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class Model, class Parameters, class Variables>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form> &&
             !is_of_this_template_type_v<Variables, std::vector>)
auto simulate(mt_64i& mt,
              const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                               variance_correction, family, Model, qdt_method,
                                               variance_form>& lik,
              Parameters const& p, const Variables& var) {
    return Macro_DMR{}
        .sample(mt, lik.m, p, var, make_substep_simulation_parameters(lik.n_sub_dt()))
        .value()();
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class family, class qdt_method, class variance_form,
          class Model, class Parameters, class Variables>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             uses_family_aproximation_c<family> &&
             var::is_this_constexpr_Var_v<qdt_method, int, uses_qdt_method> &&
             uses_variance_form_aproximation_c<variance_form>)
auto simulate(mt_64i& mt,
              const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                               variance_correction, family, Model, qdt_method,
                                               variance_form>& lik,
              Parameters const& p, const std::vector<Variables>& var) {
    std::vector<Recording> out(var.size());
    for (std::size_t i = 0; i < var.size(); ++i)
        out[i] = get<Recording>(simulate(mt, lik, p, var[i]));
    return out;
}

static Experiment_step average_Experimental_step(Experiment_step const& x,
                                                 bool average_the_evolution) {
    return Experiment_step(get<Time>(x),
                           average_agonist_step(get<Agonist_evolution>(x), average_the_evolution));
}

static Agonist_evolution add_agonist_step_i(std::vector<Agonist_step>&& x,
                                            std::vector<Agonist_step>&& y) {
    if (x.empty())
        return std::move(y);
    else {
        for (std::size_t i = 0; i < y.size(); ++i) {
            if (get<Agonist_concentration>(x.back())() == get<Agonist_concentration>(y[i])()) {
                get<number_of_samples>(x.back())() += get<number_of_samples>(y[i])();
            } else {
                x.push_back(y[i]);
            }
        }
        return std::move(x);
    }
}

static Agonist_evolution add_agonist_step_i(Agonist_step&& x, std::vector<Agonist_step>&& y) {
    return add_agonist_step_i(std::vector<Agonist_step>{std::move(x)}, std::move(y));
}

static Agonist_evolution add_agonist_step_i(std::vector<Agonist_step>&& x, Agonist_step&& e) {
    if (x.empty()) {
        x.push_back(e);
        return x;
    } else if (get<Agonist_concentration>(x.back())() == get<Agonist_concentration>(e)()) {
        get<number_of_samples>(x.back())() += get<number_of_samples>(e)();
    } else {
        x.push_back(std::move(e));
    }
    return x;
}
static Agonist_evolution add_agonist_step_i(Agonist_step&& x, Agonist_step&& e) {
    if (get<Agonist_concentration>(x)() == get<Agonist_concentration>(e)()) {
        get<number_of_samples>(x)() += get<number_of_samples>(e)();
        return std::vector<Agonist_step>{std::move(x)};
    } else {
        return std::vector<Agonist_step>{std::move(x), std::move(e)};
    }
}

static Agonist_evolution add_agonist_step(Agonist_evolution&& x, Agonist_evolution&& e) {
    return add_agonist_step_i(std::move(x()), std::move(e()));
}

inline void report_title(save_Predictions<var::Parameters_transformed>& s,
                         thermo_mcmc<var::Parameters_transformed> const&, ...) {
    s.f << "iter" << s.sep << "iter_time" << s.sep << "i_beta" << s.sep << "num_beta" << s.sep
        << "beta" << s.sep << "i_walker" << s.sep << "id_walker" << s.sep << "i_step" << s.sep
        << "time" << s.sep << "num_samples" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep
        << "Y_var" << s.sep << "plogL" << s.sep << "pelogL"
        << "\n";

    s.g << "iter" << s.sep << "iter_time" << s.sep << "i_beta" << s.sep << "num_beta" << s.sep
        << "beta" << s.sep << "i_walker" << s.sep << "id_walker" << s.sep << "i_step" << s.sep
        << "i_state" << s.sep << "j_state" << s.sep << "moment" << s.sep << "value"
        << "\n";
}

inline void report_title(save_Predictions<Matrix<double>>& s, thermo_mcmc<Matrix<double>> const&,
                         ...) {}

template <class FunctionTable, class Duration, class Prior, class t_logLikelihood, class Data,
          class Variables>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, FuncMap_St>)
void report(FunctionTable& f, std::size_t iter, const Duration& dur,
            save_Predictions<var::Parameters_transformed>& s,
            thermo_mcmc<var::Parameters_transformed> const& data, Prior const&,
            t_logLikelihood const& lik, const Data& y, const Variables& x, ...) {
    auto num_samples = size(y);
    auto num_states = lik.m.number_of_states();
    auto num_states_a = std::pow(2, std::round(std::log2(num_states)));
    std::size_t num_values = 4;
    std::size_t num_beta_portions = 2;

    std::size_t num_samples_a = std::pow(2, std::round(std::log2(size(y))));
    std::size_t point_size =
        num_values * num_beta_portions * data.get_Walkers_number() * num_samples_a;
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    std::size_t state_num_values = 1;
    std::size_t num_moments = 2;
    std::size_t state_point_size = state_num_values * num_beta_portions *
                                   data.get_Walkers_number() * num_samples_a * num_states_a *
                                   num_moments;
    std::size_t state_sampling_interval =
        std::max(s.sampling_interval, state_point_size / s.max_number_of_values_per_iteration);

    if ((iter == 0) || (iter % sampling_interval != 0))
        return;
    auto ff = f.fork(omp_get_max_threads());

    auto all_Predictions = std::vector<std::vector<std::decay_t<decltype(logLikelihoodPredictions(
        ff[0], lik, data.get_Parameter(0, 0), y, x))>>>(data.get_Walkers_number());

    auto beta = data.get_Beta();
    auto num_beta = beta.size();

#pragma omp parallel for
    for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number(); ++i_walker) {
        for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
            if ((beta[i_b] == 1) || (iter % (beta.size() * sampling_interval) == 0)) {
                auto i_th = omp_get_thread_num();

                auto par = data.get_Parameter(i_walker, i_b);
                auto walker_id = data.get_Walker(i_walker, i_b);
                all_Predictions[i_walker].push_back(
                    logLikelihoodPredictions(ff[i_th], lik, par.to_value(), y, x));
            }
        }
    }
    f += ff;
    for (std::size_t half = 0; half < 2; ++half) {
        for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
            auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
            for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
                if ((beta[i_b] == 1) || (iter % (beta.size() * sampling_interval) == 0)) {
                    auto par = data.get_Parameter(i_walker, i_b);
                    auto walker_id = data.get_Walker(i_walker, i_b);
                    auto prediction = (iter % (beta.size() * sampling_interval) == 0)
                                          ? all_Predictions[i_walker][i_b]
                                          : all_Predictions[i_walker][0];
                    if (is_valid(prediction)) {
                        auto& Macro_predictions = prediction.value();
                        auto& predictions = get<Evolution>(Macro_predictions);
                        for (std::size_t i_step = 0; i_step < size(y); ++i_step) {
                            auto v_ev =
                                get<Agonist_evolution>(get<Recording_conditions>(x)()[i_step]);

                            auto time = get<Time>(get<Recording_conditions>(x)()[i_step]);
                            auto num_smples = get_num_samples(v_ev);

                            s.f << iter << s.sep << dur.count() << s.sep << i_b << s.sep
                                << beta.size() << s.sep << beta[i_b] << s.sep << i_walker << s.sep
                                << walker_id << s.sep << i_step << s.sep << time << s.sep
                                << num_samples << s.sep << y()[i_step]() << s.sep
                                << get<y_mean>(predictions()[i_step]) << s.sep
                                << get<y_var>(predictions()[i_step]) << s.sep;
                            //  << get<plogL>(predictions()[i_step]) << s.sep
                            //  << get<eplogL>(predictions()[i_step]) << "\n";

                            if (((iter % state_sampling_interval == 0) && (beta[i_b] == 1)) ||
                                (iter % (state_sampling_interval * num_beta) == 0)) {
                                auto& v_P = get<P_mean>(predictions()[i_step]);
                                auto& v_Pc = get<P_Cov>(predictions()[i_step]);
                                auto n = v_P().size();
                                for (std::size_t i_state = 0; i_state < n; ++i_state) {
                                    s.g << iter << s.sep << dur.count() << s.sep << i_b << s.sep
                                        << beta[i_b] << s.sep << i_walker << s.sep << walker_id
                                        << s.sep << i_step << s.sep << i_state << s.sep << 0
                                        << s.sep << "mean" << s.sep << v_P()[i_state] << "\n";
                                    if ((iter % (state_sampling_interval * num_states) == 0) &&
                                            ((beta[i_b] == 1)) ||
                                        (iter % (state_sampling_interval * num_states * num_beta) ==
                                         0))
                                        for (std::size_t j_state = 0; j_state <= i_state;
                                             ++j_state) {
                                            s.g << iter << s.sep << dur.count() << s.sep << i_b
                                                << s.sep << beta[i_b] << s.sep << i_walker << s.sep
                                                << walker_id << s.sep << i_step << s.sep << i_state
                                                << s.sep << j_state << s.sep << "Cov" << s.sep
                                                << v_Pc()(i_state, j_state) << "\n";
                                        }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <class Parameters>
void report_title(save_Predictions<Parameters>& s, thermo_levenberg_mcmc const&, ...) {
    s.f << "iter" << s.sep << "iter_time" << s.sep << "beta" << s.sep << "walker_id" << s.sep
        << "i_step" << s.sep << "time" << s.sep << "num_samples" << s.sep << "agonist" << s.sep
        << "Agonist_evolution" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_var" << s.sep
        << "plogL" << s.sep << "pelogL"
        << "\n";
}

template <class FunctionTable, class Duration, class Prior, class t_logLikelihood, class Data,
          class Variables>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, FuncMap_St>)
void report(FunctionTable& f, std::size_t iter, const Duration& dur,
            save_Predictions<var::Parameters_transformed>& s, thermo_levenberg_mcmc const& data,
            by_beta<double> t_beta, Prior&&, t_logLikelihood&& lik, const Data& y,
            const Variables& x, ...) {
    std::size_t num_values = 8;
    std::size_t point_size = num_values * size(y) * t_beta.size();

    if ((iter == 0) ||
        (iter % std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration) !=
         0))
        return;
    // std::cerr<<"report save_Predictions\n";

    auto ff = f.fork(omp_get_max_threads());

    auto all_Predictions = std::vector<std::decay_t<decltype(logLikelihoodPredictions(
        ff[0], lik, data.walkers[0].m_data.m_x, y, x))>>(t_beta.size());

    auto num_samples = size(y);
#pragma omp parallel for
    for (std::size_t i_b = 0; i_b < t_beta.size(); ++i_b) {
        auto i_th = omp_get_thread_num();

        auto par = data.walkers[i_b].m_data.m_x;
        auto walker_id = data.walkers[i_b].m_data.i_walker;

        all_Predictions[i_b] = logLikelihoodPredictions(ff[i_th], lik, par.to_value(), y, x);
    }

    f += ff;
    for (std::size_t i_b = 0; i_b < t_beta.size(); ++i_b) {
        auto par = data.walkers[i_b].m_data.m_x;
        auto walker_id = data.walkers[i_b].m_data.i_walker;
        auto prediction = all_Predictions[i_b];
        if (is_valid(prediction)) {
            auto& predictions = prediction.value();
            for (std::size_t i_step = 0; i_step < size(y); ++i_step) {
                auto v_ev = get<Agonist_evolution>(get<Recording_conditions>(x)()[i_step]);

                auto time = get<Time>(get<Recording_conditions>(x)()[i_step]);
                auto num_smples = get_num_samples(v_ev);

                s.f << iter << s.sep << dur.count() << s.sep << t_beta[i_b] << s.sep << walker_id
                    << s.sep << i_step << s.sep << time << s.sep << num_samples << s.sep
                    << ToString(average_agonist_step(v_ev, true)) << s.sep << ToString(v_ev)
                    << s.sep << y()[i_step]() << s.sep

                    << get<y_mean>(predictions()[i_step]) << s.sep
                    << get<y_var>(predictions()[i_step]) << s.sep;
                //  << get<plogL>(predictions()[i_step]) << s.sep
                //  << get<eplogL>(predictions()[i_step]) << "\n";
            }
        }
    }
    //   std::cerr<<"report save_Predictions end\n";
}

inline std::string ToString(const Agonist_step& ev) {
    return std::to_string(get<number_of_samples>(ev)()) + "  " +
           std::to_string(get<Agonist_concentration>(ev)());
}

inline std::string ToString(const std::vector<Agonist_step>& ev) {
    std::string out;
    for (auto& e : ev) {
        out += ToString(e) + "  ";
    }
    return out;
}

inline std::string ToString(const Agonist_evolution& ev) {
    return ToString(ev());
}

template <typename Simulate_tag>

void report(std::string filename, const Macro_State_Ev_predictions& predictions,
            const Simulated_Recording<Simulate_tag>& y, const Experiment& xs) {
    auto& ys = get<Recording>(y());
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "i_step"
      << ","
      << "time"
      << ","
      << "num_samples"
      << ","
      << "Agonist_step"
      << ","
      << "v_ev"
      << ","
      << "y"
      << ","
      << "y_mean"
      << ","
      << "y_var"
      << ","
      << "plogL"
      << ","
      << "eplogL"
      << ","
      << "logL"
      << ","
      << "i"
      << ","
      << "P_mean"
      << ","
      << "j"
      << ","
      << "P_Cov";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        f << ","
          << "N";
    }
    f << "\n";
    for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
        auto v_ev = get<Agonist_evolution>(get<Recording_conditions>(xs)()[i_step]);
        // for (std::size_t i = 0; i < get<P_Cov>(predictions()[i_step])().nrows(); ++i) {
        //     for (std::size_t j = 0; j < get<P_Cov>(predictions()[i_step])().ncols(); ++j) {
        //         f << i_step << "," << get<Time>(get<Recording_conditions>(xs)()[i_step]) << ","
        //           << get_num_samples(v_ev) << "," << ToString(average_agonist_step(v_ev, true))
        //           << "," << ToString(v_ev) << "," << ys()[i_step]() << ","
        //           << get<y_mean>(predictions()[i_step]) << "," << get<y_var>(predictions()[i_step])
        //           << "," << get<plogL>(predictions()[i_step]) << ","
        //           << get<eplogL>(predictions()[i_step]) << "," << get<logL>(predictions()[i_step])
        //           << "," << i << "," << get<P_mean>(predictions()[i_step])()[i] << "," << j << ","
        //           << get<P_Cov>(predictions()[i_step])()(i, j);
        //         if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        //             f << "," << get<N_Ch_State_Evolution>(y())()[i_step]()[i] << "\n";
        //         else
        //             f << "\n";
        //     }
        // }
    }
}

template <class FunctionTable, class Duration, class Prior, class t_logLikelihood>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, FuncMap_St>)
void report(FunctionTable&, std::size_t iter, const Duration& dur,
            save_RateParameter<var::Parameters_transformed>& s,
            thermo_mcmc<var::Parameters_transformed> const& data, Prior const&,
            t_logLikelihood const& lik, ...) {
    auto num_states = lik.m.number_of_states();
    std::size_t num_values = 1;
    std::size_t num_beta_portions = 2;
    std::size_t point_size =
        num_values * num_beta_portions * data.get_Walkers_number() * num_states * num_states;
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if ((iter == 0) || (iter % sampling_interval != 0))
        return;
    auto& model = lik.m;
    auto beta = data.get_Beta();

    for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number(); ++i_walker) {
        for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
            if (beta[i_b] == 1) {
                auto par = data.get_Parameter(i_walker, i_b);
                auto walker_id = data.get_Walker(i_walker, i_b);
                auto Maybe_mo = model(par.to_value());
                if (is_valid(Maybe_mo)) {
                    auto& mo = Maybe_mo.value();
                    auto v_Q0 = get<Q0>(mo);
                    auto v_Qa = get<Qa>(mo);

                    for (std::size_t i_from = 0; i_from < v_Q0().nrows(); ++i_from)
                        for (std::size_t i_to = 0; i_to < v_Q0().ncols(); ++i_to) {
                            if (v_Qa()(i_from, i_to) > 0)
                                s.f << iter << s.sep << dur.count() << s.sep << i_b << s.sep
                                    << beta.size() << s.sep << beta[i_b] << s.sep << i_walker
                                    << s.sep << walker_id << s.sep << "agonist" << s.sep << i_from
                                    << s.sep << i_to << s.sep << v_Qa()(i_from, i_to) << "\n";
                            if (v_Q0()(i_from, i_to) > 0)
                                s.f << iter << s.sep << dur.count() << s.sep << i_b << s.sep
                                    << beta.size() << s.sep << beta[i_b] << s.sep << i_walker
                                    << s.sep << walker_id << s.sep << "no_agonist" << s.sep
                                    << i_from << s.sep << i_to << s.sep << v_Q0()(i_from, i_to)
                                    << "\n";
                        }
                }
            }
        }
    }
}

template <typename Simulate_tag>

void save_Likelihood_Predictions(std::string filename, const Macro_State_Ev_diagnostic& predictions,
                                 const Simulated_Recording<Simulate_tag>& y, const Experiment& xs) {
    auto& ys = get<Recording>(y());
    std::ofstream f(filename + "_patch_state_evolution.csv");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "i_step"
      << ","
      << "time"
      << ","
      << "num_samples"
      << ","
      << "Agonist_step"
      << ","
      << "v_ev"
      << ","
      << "y"
      << ","
      << "y_mean"
      << ","
      << "y_var"
      << ","
      << "plogL"
      << ","
      << "eplogL"
      << ","
      << "logL"
      << ","
      << "i"
      << ","
      << "P_mean"
      << ","
      << "j"
      << ","
      << "P_Cov";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << ","
          << "N"
          << "\n";
    else
        f << "\n";
    // for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
    //     auto v_ev = get<Agonist_evolution>(get<Recording_conditions>(xs)()[i_step]);
    //     for (std::size_t i = 0; i < get<P_Cov>(predictions()[i_step])().nrows(); ++i) {
    //         for (std::size_t j = 0; j < get<P_Cov>(predictions()[i_step])().ncols(); ++j) {
    //             f << i_step << "," << get<Time>(get<Recording_conditions>(xs)()[i_step]) << ","
    //               << get_num_samples(v_ev) << "," << ToString(average_agonist_step(v_ev, true))
    //               << "," << ToString(v_ev) << "," << ys()[i_step]() << ","
    //               << get<y_mean>(predictions()[i_step]) << "," << get<y_var>(predictions()[i_step])
    //               << "," << get<plogL>(predictions()[i_step]) << ","
    //               << get<eplogL>(predictions()[i_step]) << "," << get<logL>(predictions()[i_step])
    //               << "," << i << "," << get<P_mean>(predictions()[i_step])()[i] << "," << j << ","
    //               << get<P_Cov>(predictions()[i_step])()(i, j);
    //             if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
    //                 f << "," << get<N_Ch_State_Evolution>(y())()[i_step]()[i] << "\n";
    //             else
    //                 f << "\n";
    //         }
    //     }
    // }
}

inline void save_Likelihood_Predictions(std::string filename,
                                        const Macro_State_Ev_predictions& predictions,
                                        const v_Simulated_Recording& sim_y, const Experiment& xs) {
    std::visit([&](auto& y) { save_Likelihood_Predictions(filename, predictions, y, xs); }, sim_y);
}

template <typename Simulate_tag>

void save_Likelihood_Predictions(std::string filename,
                                 const Macro_State_Ev_predictions& predictions,
                                 const Simulated_Recording<Simulate_tag>& y, const Experiment& xs) {
    auto& ys = get<Recording>(y());
    std::ofstream flogL(filename + "_logL.csv");
    flogL << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "logL" << ","
          << "elogL" << "," << "vlogL" << "\n";
    flogL << get<logL>(predictions)() << "," << get<elogL>(predictions)() << ","
          << get<vlogL>(predictions)() << "\n";
    std::ofstream f(filename + "_logL_evolution.csv");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "i_step"
      << ","
      << "time"
      << ","
      << "num_samples"
      << ","
      << "Agonist_step"
      << ","
      << "v_ev"
      << ","
      << "y"
      << ","
      << "y_mean"
      << ","
      << "y_var";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << ","
          << "N"
          << "\n";
    else
        f << "\n";
    for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
        auto v_ev = get<Agonist_evolution>(get<Recording_conditions>(xs)()[i_step]);
        f << i_step << "," << get<Time>(get<Recording_conditions>(xs)()[i_step]) << ","
          << get_num_samples(v_ev) << "," << ToString(average_agonist_step(v_ev, true)) << ","
          << ToString(v_ev) << "," << ys()[i_step]() << "," /*
          << get<ymean_Evolution>(predictions)()[i_step] << ","
          << get<yvar_Evolution>(predictions)()[i_step] */
          << "\n";
    }
}
inline auto new_thermo_Model_by_max_iter_dts(
    std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
    std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium, std::size_t beta_size,
    Saving_intervals sint, std::size_t initseed, std::size_t t_adapt_beta_every,
    std::string t_adapt_beta_equalizer, std::string t_adapt_beta_constroler,
    std::string t_adapt_beta_variance, double t_adapt_beta_nu, double t_adapt_beta_t0,
    double t_adapt_beta_threshold, bool t_adjust_beta, double t_acceptance_upper_limit,
    double t_acceptance_lower_limit, double t_desired_acceptance) {
    return new_thermodynamic_integration(
        thermo_less_than_max_iteration(max_iter_equilibrium),
        save_mcmc<var::Parameters_transformed, save_Iter,
                  save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>,
                  save_RateParameter<var::Parameters_transformed>, save_Evidence,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, std::pair(1ul, 1ul), get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<Save_RateParameter_every>(sint())(),
            get<Save_Evidence_every>(sint())(), get<Save_Predictions_every>(sint())()),
        num_scouts_per_ensemble, thermo_jumps_every, beta_size, initseed, t_adapt_beta_every,
        t_adapt_beta_equalizer, t_adapt_beta_constroler, t_adapt_beta_variance, t_adapt_beta_nu,
        t_adapt_beta_t0, t_adapt_beta_threshold, t_adjust_beta, t_acceptance_upper_limit,
        t_acceptance_lower_limit, t_desired_acceptance);
}

#ifdef ZOMBIE
namespace zombie {

template <typename Simulate_tag>

void save_fractioned_Likelihood_Predictions(
    std::string filename, const std::string& sep,
    const std::vector<Macro_State_Evolution>& v_predictions,
    const std::vector<Simulated_Recording<Simulate_tag>>& v_y,
    const std::vector<Experiment>& v_xs) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    f << "i_frac" << sep << "i_step" << sep << "t_ini" << sep << "time" << sep << "i_sub_step"
      << sep << "number_of_samples" << sep << "Agonist_concentration" << sep << "macror_algorithm"
      << sep << "y" << sep << "y_mean" << sep << "y_var" << sep << "plogL" << sep << "eplogL" << sep
      << "logL" << sep << "i_state" << sep << "P_mean" << sep << "j_state" << sep << "P_Cov";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << ","
          << "N"
          << "\n";
    else
        f << "\n";
    for (std::size_t i_frac = 0; i_frac < v_y.size(); ++i_frac) {
        double t_ini = 0;
        auto& y = v_y[i_frac];
        auto& ys = get<Recording>(y());
        auto& xs = v_xs[i_frac];
        auto& predictions = v_predictions[i_frac];
        for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
            auto v_ev = get<Agonist_evolution>(get<Recording_conditions>(xs)()[i_step]);
            for (std::size_t i_sub_step = 0; i_sub_step < v_ev.size(); ++i_sub_step) {
                for (std::size_t i_state = 0; i_state < get<P_Cov>(predictions()[i_step])().nrows();
                     ++i_state) {
                    for (std::size_t j_state = 0;
                         j_state < get<P_Cov>(predictions()[i_step])().ncols(); ++j_state) {
                        f << i_frac << sep << i_step << sep << t_ini << sep
                          << get<Time>(get<Recording_conditions>(xs)()[i_step]) << sep << i_sub_step
                          << sep << get<number_of_samples>(v_ev[i_sub_step]) << sep
                          << get<Agonist_concentration>(v_ev[i_sub_step]) << sep
                          << get<macror_algorithm>(predictions()[i_step]) << sep << ys()[i_step]()
                          << sep << get<y_mean>(predictions()[i_step]) << sep
                          << get<y_var>(predictions()[i_step]) << sep
                          << get<plogL>(predictions()[i_step]) << sep
                          << get<eplogL>(predictions()[i_step]) << sep
                          << get<logL>(predictions()[i_step]) << sep << i_state << sep
                          << get<P_mean>(predictions()[i_step])()[i_state] << sep << j_state << sep
                          << get<P_Cov>(predictions()[i_step])()(i_state, j_state);
                        if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
                            f << sep << get<N_Ch_State_Evolution>(y())()[i_step]()[i_state] << "\n";
                        else
                            f << "\n";
                    }
                }
                t_ini +=
                    get<number_of_samples>(v_ev[i_sub_step])() / get<Frequency_of_Sampling>(xs)();
            }
        }
    }
}

template <class ParameterType, class FunctionTable, class Duration, class Prior,
          class t_logLikelihood, class Data, class Variables>
    requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
void report(FunctionTable& f, std::size_t iter, const Duration& dur,
            save_Predictions<ParameterType>& s, cuevi::Cuevi_mcmc<ParameterType>& data,
            Prior const&, t_logLikelihood const& lik, const cuevi::by_fraction<Data>& ys,
            const cuevi::by_fraction<Variables>& xs, ...) {
    auto& t = data.get_Cuevi_Temperatures();
    auto num_fr = data.get_Cuevi_Temperatures_Number();
    std::size_t num_values = 8;
    std::size_t point_size =
        num_values * data.get_Walkers_number() * num_fr * data.get_Parameters_number();
    std::size_t sampling_interval =
        std::max(s.sampling_interval, point_size / s.max_number_of_values_per_iteration);

    if ((iter == 0) || (iter % sampling_interval != 0)) {
        return;
    }
    data.calculate_Likelihoods_for_Evidence_calulation(f, lik, ys, xs);

    auto ff = f.fork(omp_get_max_threads());
    using Predition_type = std::decay_t<decltype(logLikelihoodPredictions(
        ff[0], lik, data.get_Parameter(cuevi::Walker_Index(0), cuevi::Cuevi_Index(0)), ys[0],
        xs[0]))>;
    auto allPredictions = std::vector<std::vector<Predition_type>>(
        data.get_Cuevi_Temperatures_Number(),
        std::vector<Predition_type>(data.get_Walkers_number()));
#pragma omp parallel for  // collapse(2)
    for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number(); ++i_walker) {
        for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
            auto iw = cuevi::Walker_Index(i_walker);
            auto i_th = omp_get_thread_num();

            auto icu = cuevi::Cuevi_Index(i_cu);
            auto i_frac = data.get_Fraction(i_cu);
            auto beta = data.get_Beta(icu);
            auto nsamples = size(ys[i_frac()]);
            auto& wa = data.get_Walker(iw, icu);
            auto& wav = data.get_Walker_Value(iw, icu);
            auto par = data.get_Parameter(iw, i_cu);
            auto prediction =
                logLikelihoodPredictions(ff[i_th], lik, par.to_value(), ys[i_frac()], xs[i_frac()]);
            allPredictions[i_cu][i_walker] = std::move(prediction);
        }
    }
    f += ff;
    for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number(); ++i_cu) {
        auto icu = cuevi::Cuevi_Index(i_cu);
        auto i_frac = data.get_Fraction(i_cu);
        auto beta = data.get_Beta(icu);
        auto nsamples = size(ys[i_frac()]);
        for (std::size_t half = 0; half < 2; ++half)
            for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
                auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
                auto iw = cuevi::Walker_Index(i_walker);
                auto& wa = data.get_Walker(iw, icu);
                auto& wav = data.get_Walker_Value(iw, icu);
                auto par = data.get_Parameter(iw, i_cu);
                auto& prediction = allPredictions[i_cu][i_walker];
                if (is_valid(prediction)) {
                    auto& predictions = prediction.value();
                    for (std::size_t i_step = 0; i_step < size(ys[i_frac()]); ++i_step) {
                        auto v_ev = get<Agonist_evolution>(
                            get<Recording_conditions>(xs[i_frac()])()[i_step]);

                        s.f << iter << s.sep << dur.count() << s.sep << i_cu << s.sep << i_frac()
                            << s.sep << nsamples << s.sep << beta() << s.sep << i_walker << s.sep
                            << get<cuevi::Walker_id>(wa())() << s.sep << i_step << s.sep
                            << get<Time>(get<Recording_conditions>(xs[i_frac()])()[i_step]) << s.sep
                            << get_num_samples(v_ev) << s.sep
                            << ToString(average_agonist_step(v_ev, true)) << s.sep << ToString(v_ev)
                            << s.sep << ys[i_frac()]()[i_step]() << s.sep
                            << get<y_mean>(predictions()[i_step]) << s.sep
                            << get<y_var>(predictions()[i_step]) << s.sep
                            << get<plogL>(predictions()[i_step]) << s.sep
                            << get<eplogL>(predictions()[i_step]) << "\n";
                    }
                }
            }
    }
}
template <class ParameterType>
void report_title(save_Predictions<ParameterType>& s, cuevi::Cuevi_mcmc<ParameterType> const&,
                  ...) {
    s.f << "iter" << s.sep << "iter_time" << s.sep << "i_cu" << s.sep << "i_frac" << s.sep
        << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep << "Walker_id" << s.sep
        << "i_step" << s.sep << "Time" << s.sep << "num_samples" << s.sep << "average_agonist_step"
        << s.sep << "v_ev" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_var" << s.sep
        << "plogL" << s.sep << "eplogL"
        << "\n";
}

/*
auto cuevi_Model_by_convergence(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments, bool average_the_agonist_evolution,
    
    std::size_t num_scouts_per_ensemble, double min_fraction,
    std::size_t thermo_jumps_every, std::size_t max_iter, double max_ratio,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double stops_at, bool includes_zero, std::size_t initseed) {
  return deprecated::cuevi_integration(
      checks_derivative_var_ratio<deprecated::cuevi_mcmc,
                                  var::Parameters_transformed>(max_iter,
                                                                   max_ratio),
      experiment_fractioner(t_segments, average_the_agonist_evolution),
      save_mcmc<var::Parameters_transformed,
                save_likelihood<var::Parameters_transformed>,
                save_Parameter<var::Parameters_transformed>, save_Evidence,
                save_Predictions<var::Parameters_transformed>>(
          path, filename, 100ul, 100ul, 100ul, 100ul),
      num_scouts_per_ensemble, min_fraction, thermo_jumps_every,
      n_points_per_decade_beta, n_points_per_decade_fraction, stops_at,
      includes_zero, initseed);
}
*/

cuevi::Cuevi_Algorithm<
    experiment_fractioner,
    save_mcmc<var::Parameters_transformed, save_likelihood<var::Parameters_transformed>,
              save_Parameter<var::Parameters_transformed>, save_Evidence,
              save_Predictions<var::Parameters_transformed>>,
    cuevi_less_than_max_iteration>
    new_cuevi_Model_by_iteration(
        std::string path, std::string filename, const std::vector<std::size_t>& t_segments,
        bool average_the_agonist_evolution, std::size_t num_scouts_per_ensemble,
        std::size_t number_trials_until_give_up, double min_fraction,
        std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium,
        double n_points_per_decade_beta, double n_points_per_decade_fraction, double medium_beta,
        double stops_at, bool includes_the_zero, Saving_intervals sint, bool random_jumps) {
    return cuevi::Cuevi_Algorithm(
        experiment_fractioner(t_segments, average_the_agonist_evolution),
        save_mcmc<var::Parameters_transformed, save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>, save_Evidence,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<Save_Evidence_every>(sint())(),
            get<Save_Predictions_every>(sint())()),
        cuevi_less_than_max_iteration(max_iter_equilibrium),
        cuevi::Num_Walkers_Per_Ensemble(num_scouts_per_ensemble),
        cuevi::Fractions_Param(
            Vector_Space(cuevi::Min_value(min_fraction),
                         cuevi::Points_per_decade(n_points_per_decade_fraction))),
        cuevi::Th_Beta_Param(Vector_Space(  // Includes_zero,
            // Med_value,Points_per_decade,Min_value,
            // Points_per_decade_low
            cuevi::Includes_zero(includes_the_zero), cuevi::Med_value(medium_beta),
            cuevi::Points_per_decade(n_points_per_decade_fraction), cuevi::Min_value(stops_at),
            cuevi::Points_per_decade_low(n_points_per_decade_beta))),
        cuevi::Number_trials_until_give_up(number_trials_until_give_up),
        cuevi::Thermo_Jumps_every(thermo_jumps_every), cuevi::Random_jumps(random_jumps),
        std::move(sint));
}

cuevi::Cuevi_Algorithm_no_Fractioner<
    save_mcmc<var::Parameters_transformed, save_likelihood<var::Parameters_transformed>,
              save_Parameter<var::Parameters_transformed>, save_Evidence,
              save_Predictions<var::Parameters_transformed>>,
    cuevi_less_than_max_iteration>
    new_cuevi_Model_already_fraction_by_iteration(
        std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
        std::size_t number_trials_until_give_up, std::size_t thermo_jumps_every,
        std::size_t max_iter_equilibrium, double n_points_per_decade_beta_high,
        double n_points_per_decade_beta_low, double medium_beta, double stops_at,
        bool includes_the_zero, Saving_intervals sint, bool random_jumps) {
    return cuevi::Cuevi_Algorithm_no_Fractioner(
        save_mcmc<var::Parameters_transformed, save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>, save_Evidence,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<Save_Evidence_every>(sint())(),
            get<Save_Predictions_every>(sint())()),
        cuevi_less_than_max_iteration(max_iter_equilibrium),
        cuevi::Num_Walkers_Per_Ensemble(num_scouts_per_ensemble),
        cuevi::Th_Beta_Param(Vector_Space(  // Includes_zero,
            // Med_value,Points_per_decade,Min_value,
            // Points_per_decade_low
            cuevi::Includes_zero(includes_the_zero), cuevi::Med_value(medium_beta),
            cuevi::Points_per_decade(n_points_per_decade_beta_high), cuevi::Min_value(stops_at),
            cuevi::Points_per_decade_low(n_points_per_decade_beta_low))),
        cuevi::Number_trials_until_give_up(number_trials_until_give_up),
        cuevi::Thermo_Jumps_every(thermo_jumps_every), cuevi::Random_jumps(random_jumps),
        std::move(sint));
}

auto new_thermo_Model_by_max_iter(std::string path, std::string filename,
                                  std::size_t num_scouts_per_ensemble,
                                  std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium,
                                  std::size_t beta_size, std::size_t beta_upper_size,
                                  std::size_t beta_medium_size, double beta_upper_value,
                                  double beta_medium_value,

                                  double stops_at, bool includes_zero, Saving_intervals sint,
                                  std::size_t initseed) {
    return new_thermodynamic_integration(
        thermo_less_than_max_iteration(max_iter_equilibrium),
        save_mcmc<var::Parameters_transformed, save_Iter,
                  save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>,
                  save_RateParameter<var::Parameters_transformed>, save_Evidence,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, std::pair(1ul, 1ul), get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<Save_RateParameter_every>(sint())(),
            get<Save_Evidence_every>(sint())(), get<Save_Predictions_every>(sint())()),
        num_scouts_per_ensemble, thermo_jumps_every, beta_size, beta_upper_size, beta_medium_size,
        beta_upper_value, beta_medium_value, stops_at, includes_zero, initseed);
}
}  // namespace zombie
namespace zombie {

auto thermo_levenberg_Model_by_max_iter(
    std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
    std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium, std::size_t beta_size,
    std::size_t beta_upper_size, std::size_t beta_medium_size, double beta_upper_value,
    double beta_medium_value, std::size_t n_lambdas, std::string lambda_adaptive_algorithm,
    double stops_at, bool includes_zero, Saving_Levenberg_intervals sint, std::size_t initseed,
    double dp) {
    return thermodynamic_levenberg_integration(
        thermo_less_than_max_iteration(max_iter_equilibrium),
        save_mcmc<var::Parameters_transformed, save_Iter,
                  save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>,
                  save_Levenberg_Lambdas<var::Parameters_transformed>,
                  save_Levenberg_Errors<var::Parameters_transformed>,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, std::pair(1ul, 1ul), get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<save_Levenberg_Lambdas_every>(sint())(),
            get<save_Levenberg_Errors_every>(sint())(), get<Save_Predictions_every>(sint())()),
        num_scouts_per_ensemble, thermo_jumps_every, beta_size, beta_upper_size, beta_medium_size,
        beta_upper_value, beta_medium_value, n_lambdas, lambda_adaptive_algorithm, stops_at,
        includes_zero, initseed, dp);
}

auto thermo_Model_by_max_iter(std::string path, std::string filename,
                              std::size_t num_scouts_per_ensemble,
                              std::size_t max_num_simultaneous_temperatures,
                              std::size_t thermo_jumps_every, std::size_t max_iter_warming,
                              std::size_t max_iter_equilibrium, double n_points_per_decade,
                              double stops_at, bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(
        less_than_max_iteration(max_iter_warming, max_iter_equilibrium),
        save_mcmc<var::Parameters_transformed, save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>, save_Evidence,
                  save_Predictions<var::Parameters_transformed>>(
            path, filename, std::pair(10ul, 200ul), std::pair(10ul, 200ul), std::pair(10ul, 200ul),
            std::pair(10ul, 200ul)),
        num_scouts_per_ensemble, max_num_simultaneous_temperatures, thermo_jumps_every,
        n_points_per_decade, stops_at, includes_zero, initseed);
}

}  // namespace zombie
#endif

template <class Parameter, class Simulate_tag>
void report_model(save_Parameter<Parameter>& s, Simulated_Recording<Simulate_tag> const& sim) {
    std::ofstream f(s.fname + "_simulation.csv");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    f << "i_step" << s.sep << "patch_current";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << s.sep << "i_state" << s.sep << "N_state";
    f << "\n";

    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        auto N = get<N_Ch_State_Evolution>(sim());
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
            for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_step) {
                f << i_step << s.sep << y()[i_step]() << i_state << s.sep << N()[i_step]()[i_state]
                  << "\n";
            }
        }
    } else {
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
            f << i_step << s.sep << y()[i_step]() << "\n";
        }
    }
}

template <class Simulate_tag>
void save_Simulated_Recording(std::string const& filename, const std::string& separator,
                              Simulated_Recording<Simulate_tag> const& sim) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    f << "i_step" << separator << "patch_current";
    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
        f << separator << "i_state" << separator << "N_state";
    f << "\n";

    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>) {
        auto N = get<N_Ch_State_Evolution>(sim());
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
            for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
                f << i_step << separator << y()[i_step]() << separator << i_state << separator
                  << N()[i_step]()[i_state] << "\n";
            }
        }
    } else {
        auto y = get<Recording>(sim());
        for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
            f << i_step << separator << y()[i_step]() << "\n";
        }
    }
}

inline void report_model(save_Parameter<var::Parameters_transformed>& s,
                         Parameters_Transformations const& m) {
    std::ofstream f(s.fname + "_parameter.csv");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto n = m.size();
    f << "i_par" << s.sep << "moment" << s.sep << "value"
      << "\n";
    for (auto i_par = 0ul; i_par < n; ++i_par)
        f << i_par << s.sep << "mean" << s.sep << m.standard_values()[i_par] << "\n";
}

// Maybe_error<Parameters_Transformations>
// load_Parameters(save_Parameter<var::Parameters_transformed> &s) {
//   return load_Parameters(s.fname + "_parameter.csv", s.sep);
//}
namespace cmd {
#ifdef ZOMBIE
namespace zombie {

inline auto set_Fraction_algorithm(double min_fraction, double n_points_per_decade_fraction,
                                   std::string segments) {
    return std::tuple(min_fraction, n_points_per_decade_fraction, segments);
}

using fraction_algo_type =
    typename return_type<std::decay_t<decltype(&set_Fraction_algorithm)>>::type;

inline Maybe_error<std::vector<std::size_t>> load_segments_length_for_fractioning(
    const std::string& filename, std::string sep) {
    std::ifstream f(filename);
    if (!f)
        return error_message(filename + " cannot be opened");
    std::string line;
    std::getline(f, line);
    if (!f)
        return error_message(filename + " has no data");

    std::stringstream ss(line);
    std::vector<std::size_t> out;
    std::size_t number_of_samples;

    while (ss >> number_of_samples) {
        out.push_back(number_of_samples);
        ss >> septr(sep);
    }

    return out;
}

inline Maybe_error<std::tuple<std::string, std::string, double, double>> calc_experiment_fractions(
    std::string save_name, std::string recording, macrodr::cmd::experiment_type experiment,
    fraction_algo_type fraction_algo, Maybe_error<std::size_t> Maybe_num_param,
    std::size_t i_seed) {
    auto myseed = calc_seed(i_seed);

    if (!Maybe_num_param)
        return Maybe_num_param.error();
    auto num_param = Maybe_num_param.value();

    auto init_seed = calc_seed(i_seed);
    mt_64i mt(init_seed);

    auto [min_fraction, n_points_per_decade_fraction, segments] = std::move(fraction_algo);

    auto Maybe_segments = load_segments_length_for_fractioning(segments, ",");
    if (!Maybe_segments)
        return Maybe_segments.error();
    macrodr::Recording y;
    auto Maybe_y = load_Recording_Data(recording, ",", y);

    macrodr::experiment_fractioner frac(Maybe_segments.value(), 0);

    auto [ys, xs] = frac(y, experiment, mt, num_param * min_fraction, n_points_per_decade_fraction);
    auto filename = save_name + "_" + std::to_string(myseed);
    save_fractioned_experiment(filename + "_experiment.csv", ",", xs);
    save_fractioned_Recording(filename + "_recording.csv", ",", ys);
    return std::tuple(filename + "_experiment.csv", filename + "_recording.csv",
                      get<Frequency_of_Sampling>(experiment)(),
                      get<initial_agonist_concentration>(experiment)()());
}

inline Maybe_error<std::tuple<std::string, std::string, double, double>> calc_simulation_fractions(
    std::string save_name, std::string simulation, ::macrodr::cmd::experiment_type experiment,
    fraction_algo_type fraction_algo, Maybe_error<std::size_t> Maybe_num_param,
    std::size_t i_seed) {
    if (!Maybe_num_param)
        return Maybe_num_param.error();
    auto num_param = Maybe_num_param.value();
    auto myseed = calc_seed(i_seed);

    auto init_seed = calc_seed(i_seed);
    mt_64i mt(init_seed);

    auto [min_fraction, n_points_per_decade_fraction, segments] = std::move(fraction_algo);

    auto Maybe_segments = load_segments_length_for_fractioning(segments, ",");
    if (!Maybe_segments)
        return Maybe_segments.error();
    Simulated_Recording<var::please_include<N_Ch_State_Evolution>> y;
    auto Maybe_y = load_simulation(simulation, ",", y);
    if (!Maybe_y)
        return Maybe_y.error();

    macrodr::experiment_fractioner frac(Maybe_segments.value(), 0);

    auto [ys, xs] = frac(y, experiment, mt, num_param * min_fraction, n_points_per_decade_fraction);
    auto filename = save_name + "_" + std::to_string(myseed);
    save_fractioned_experiment(filename + "_frac_experiment.csv", ",", xs);
    save_fractioned_simulation(filename + "_frac_recording.csv", ",", ys);
    return std::tuple(filename + "_frac_experiment.csv", filename + "_frac_recording.csv",
                      get<Frequency_of_Sampling>(experiment)(),
                      get<initial_agonist_concentration>(experiment)()());
}

using fractioned_experiment_type =
    typename return_type<std::decay_t<decltype(&calc_experiment_fractions)>>::type;

using fractioned_simulation_type =
    typename return_type<std::decay_t<decltype(&calc_simulation_fractions)>>::type;
}  // namespace zombie
#endif  //ZOMBIE
}  // namespace cmd

}  // namespace macrodr

// Include CLI function table after all qmodel types are defined to avoid cycles
#ifndef QMODEL_NO_CLI
#include <CLI_function_table.h>
#endif

#endif  // QMODEL_H
