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
// #include "models_MoffattHume_linear.h"
#include <macrodr/dsl/type_name.h>

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
namespace macrodr {
using dsl::type_name;
using var::Parameters_Transformations;
using var::Power;
using var::Product;
using var::Var;
using var::Vector_Space;

using var::build;
using var::Constant;
using var::Fun;
using var::Op_t;
using var::primitive;
using var::Transfer_Op_to;
using var::transformation_type_t;

using var::U;

using std::abs;
using std::exp;
using std::log10;
using std::max;

using var::F;
using var::FuncMap_St;
using var::Time_it_st;

using dsl::type_name;
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

class Q0 : public var::Var<Q0, Matrix<double>> {};

inline std::size_t get_max_state(
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::string>> const& new_formulas) {
    std::size_t out = 0;
    for (auto& e : new_formulas) {
        if (e.first.first > out)
            out = e.first.first;
        if (e.first.second > out)
            out = e.first.second;
    }
    return out;
}
class Q0_formula : public var::Var<Q0_formula, std::vector<std::vector<std::string>>> {
   public:
    using var::Var<Q0_formula, std::vector<std::vector<std::string>>>::Var;

    Q0_formula(std::size_t N)
        : var::Var<Q0_formula, std::vector<std::vector<std::string>>>{
              std::vector<std::vector<std::string>>{N, std::vector<std::string>{N, ""}}} {}

    friend std::ostream& operator<<(std::ostream& os, Q0_formula const& x) {
        os << "Q0 formula"
           << "\n";
        for (std::size_t i = 0; i < x().size(); ++i) {
            for (std::size_t j = 0; j < x()[i].size(); ++j) {
                if (x()[i][j].size() > 0)
                    os << "Q(" << i << "," << j << ")->" << x()[i][j] << "\t";
            }
            os << "\n";
        }
        return os;
    }
};

class Qa : public var::Var<Qa, Matrix<double>> {};
class Qa_formula : public var::Var<Qa_formula, std::vector<std::vector<std::string>>> {
   public:
    using base_type = var::Var<Qa_formula, std::vector<std::vector<std::string>>>;
    using base_type::Var;
    Qa_formula(std::size_t N)
        : var::Var<Qa_formula, std::vector<std::vector<std::string>>>{
              std::vector<std::vector<std::string>>{N, std::vector<std::string>{N, ""}}} {}
    friend std::ostream& operator<<(std::ostream& os, Qa_formula const& x) {
        os << "Q0 formula"
           << "\n";
        for (std::size_t i = 0; i < x().size(); ++i) {
            for (std::size_t j = 0; j < x()[i].size(); ++j) {
                if (x()[i][j].size() > 0)
                    os << "Q(" << i << "," << j << ")->" << x()[i][j] << "\t";
            }
            os << "\n";
        }
        return os;
    }
};
template <class Q_formula>
    requires(std::is_same_v<Q_formula, Q0_formula> || std::is_same_v<Q_formula, Qa_formula>)
auto change_states_number(const Q_formula& f, std::size_t N) {
    Q_formula out(N);
    for (std::size_t i = 0; i < f().size(); ++i)
        for (std::size_t j = 0; j < f().size(); ++j) out()[i][j] = f()[i][j];

    return out;
}

template <class Q_formula>
    requires(std::is_same_v<Q_formula, Q0_formula> || std::is_same_v<Q_formula, Qa_formula>)
auto insert_new_formula(const Q_formula& f, std::size_t i_ini, std::size_t i_end,
                        std::string&& formula) {
    auto N = std::max(i_ini, i_end) + 1;
    auto out = change_states_number(f, N);
    out()[i_ini][i_end] = std::move(formula);

    return out;
}

class Qx : public var::Var<Qx, Matrix<double>> {};
class P_initial : public var::Var<P_initial, Matrix<double>> {};

class g : public var::Var<g, Matrix<double>> {};
class g_formula : public var::Var<g_formula, std::vector<std::string>> {};

inline auto change_states_number(const g_formula& f, std::size_t N) {
    g_formula out(std::vector<std::string>{N, ""});
    for (std::size_t i = 0; i < f().size(); ++i) out()[i] = f()[i];
    return out;
}

class N_St : public var::Constant<N_St, std::size_t> {};

class N_Ch_mean : public var::Var<N_Ch_mean, Matrix<double>> {};
class N_Ch_mean_value : public var::Var<N_Ch_mean, double> {};

class N_Ch_mean_time_segment_duration
    : public var::Constant<N_Ch_mean_time_segment_duration, double> {};

class N_Ch_init : public var::Var<N_Ch_init, double> {};
class N_Ch_eq : public var::Var<N_Ch_eq, double> {};
class N_Ch_tau : public var::Var<N_Ch_tau, double> {};

class Binomial_magical_number : public var::Constant<Binomial_magical_number, double> {};

class min_P : public var::Constant<min_P, double> {};

class N_Ch_std : public var::Var<N_Ch_std, double> {};

class Current_Noise : public var::Var<Current_Noise, double> {};

class Pink_Noise : public var::Var<Pink_Noise, double> {};

class Proportional_Noise : public var::Var<Proportional_Noise, double> {};

class Current_Baseline : public var::Var<Current_Baseline, double> {};

template <class C_Matrix>
auto to_Probability(C_Matrix const& x) -> Maybe_error<C_Matrix> {
    using std::max;
    using std::min;
    if (!isfinite(var::fullsum(x)))
        return error_message("error in Probability");

    for (std::size_t i = 0; i < x.size(); ++i)
        if (std::isnan(primitive(x[i])))

            if constexpr (var::is_derivative_v<C_Matrix>)
                for (std::size_t i = 0; i < x.size(); ++i)
                    if (std::isnan(derivative(x)()[i][0]))
                        return error_message("error in derivative");

    constexpr double smooth_eps = 0;
    auto out = apply(
        [smooth_eps](auto const& value) {
            using ValueType = std::decay_t<decltype(value)>;
            using std::sqrt;
            auto half = 0.5;
            auto eps_like = smooth_eps*smooth_eps;
            auto term = value * value + eps_like;
            auto sqrt_term = sqrt(term);
            return half * (value + sqrt_term);
        },
        x);
    if (var::min(primitive(out)) < 0)
        return error_message("how did negative Probability arise?");

    auto s = var::sum(out);
    if (s == 0)
        return error_message(" cero probability");
    return out * (1.0 / s);
}

template <class C_Matrix>
auto to_Covariance_Probability(C_Matrix const& x) -> Maybe_error<C_Matrix> {
    using std::max;
    using std::min;
    if (!isfinite(var::fullsum(x)))
        return error_message("nan Cov");
    auto out = x;
    for (auto i = 0ul; i < x.nrows(); ++i) {
        if (primitive(out(i, i)) <= 0.0) {
            for (auto j = 0ul; j < x.ncols(); ++j) {
                set(out, i, j, 0.0);
            }
        }
    }

    return out;
}

inline bool all_Probability_elements(Matrix<double> const& x) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(primitive(x[i])))
            return false;
        else if (x[i] > 1.0)
            return false;
        else if (x[i] < 0.0)
            return false;
    }
    return true;
}

template <class C_Matrix>
inline bool all_Covariance_elements(C_Matrix const& x) {
    if (x.ncols() != x.nrows())
        return false;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (!std::isfinite(primitive(x(i, i))) || x(i, i) < 0)
            return false;
        // for (std::size_t j = 0; j < i; ++j)
        //  if (std::isfinite(primitive(x(i, j))) && (x(i, i) > 0) && (x(j, j) > 0))
        //  {
        //      auto r=x(i, j) * x(i, j) / x(i, i) / x(j, j);
        //      if ( r> 1.0)
        //      return false;
        //  } else
        //{
        //   if (x(i, j) != 0)
        //     return false;
        // }
    }
    return true;
}

inline bool crude_lambda_violations(DiagonalMatrix<double> const& l) {
    if (var::max(l) > 1e-2)
        return true;
    else
        return false;
}

template <class C_Matrix>
auto to_Transition_Probability_Eigenvalues(C_Matrix&& lambda) {
    // if (crude_lambda_violations(primitive(lambda)))
    //     std::cerr<<"crude lambda violations\n";

    auto i_max = var::i_max(primitive(lambda));
    lambda.set(i_max, 0.0);
    auto j_max = var::i_max(primitive(lambda));
    while (primitive(lambda[j_max]) > 0.0) {
        lambda.set(j_max, 0.0);
        j_max = var::i_max(primitive(lambda));
    }

    return lambda;
}

struct StabilizerPolicyEnabled {
    static constexpr bool clamp_variance = true;
    static constexpr bool mask_probability = true;
    static constexpr bool enforce_gmean_bounds = true;
    static constexpr bool project_transition_probability = true;
    static constexpr bool sanitize_eigenvalues = true;
};

struct StabilizerPolicyDisabled {
    static constexpr bool clamp_variance = false;
    static constexpr bool mask_probability = false;
    static constexpr bool enforce_gmean_bounds = false;
    static constexpr bool project_transition_probability = false;
    static constexpr bool sanitize_eigenvalues = false;
};

class N_channel_state : public var::Var<N_channel_state, Matrix<double>> {
    using var::Var<N_channel_state, Matrix<double>>::Var;
};

class y_sum : public var::Var<y_sum, double> {};

class t_sum : public var::Var<t_sum, double> {};
class P_mean : public var::Var<P_mean, Matrix<double>> {
   public:
    friend std::string className(P_mean) { return "P_mean"; }
};

class P_mean_t2_y0 : public var::Var<P_mean_t2_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t2_y0) { return "P_mean_t2_y0"; }
};

class P_mean_t2_y1 : public var::Var<P_mean_t2_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t2_y1) { return "P_mean_t2_y1"; }
};

class P_mean_t15_y0 : public var::Var<P_mean_t15_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t15_y0) { return "P_mean_t15_y0"; }
};

class P_mean_t15_y1 : public var::Var<P_mean_t15_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t15_y1) { return "P_mean_t15_y1"; }
};


class P_mean_t1_y1 : public var::Var<P_mean_t1_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t1_y1) { return "P_mean_t1_y1"; }
};

class P_mean_t20_y1 : public var::Var<P_mean_t20_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t20_y1) { return "P_mean_t20_y1"; }
};

class P_mean_t11_y0 : public var::Var<P_mean_t11_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_t11_y0) { return "P_mean_t11_y0"; }
};

class P_mean_t10_y1 : public var::Var<P_mean_t10_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_t10_y1) { return "P_mean_t10_y1"; }
};

class P_mean_0t_y0 : public var::Var<P_mean_0t_y0, Matrix<double>> {
   public:
    friend std::string className(P_mean_0t_y0) { return "P_mean_0t_y0"; }
};

class P_mean_0t_y1 : public var::Var<P_mean_0t_y1, Matrix<double>> {
   public:
    friend std::string className(P_mean_0t_y1) { return "P_mean_0t_y1"; }
};


class P_Cov : public var::Var<P_Cov, SymmetricMatrix<double>> {
    friend std::string className(P_Cov) { return "P_Cov"; }
};

class P_var_ii_0t_y0 : public var::Var<P_var_ii_0t_y0, Matrix<double>> {
    friend std::string className(P_var_ii_0t_y0) { return "P_var_ii_0t_y0"; }
};


class P_var_ii_0t_y1 : public var::Var<P_var_ii_0t_y1, Matrix<double>> {
    friend std::string className(P_var_ii_0t_y1) { return "P_var_ii_0t_y1"; }
};


class P_Cov_t2_y0 : public var::Var<P_Cov_t2_y0, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t2_y0) { return "P_Cov_t2_y0"; }
};

class P_Cov_t2_y1 : public var::Var<P_Cov_t2_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t2_y1) { return "P_Cov_t2_y1"; }
};

class P_Cov_t15_y0 : public var::Var<P_Cov_t15_y0, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t15_y0) { return "P_Cov_t15_y0"; }
};

class P_Cov_t15_y1 : public var::Var<P_Cov_t15_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t15_y1) { return "P_Cov_t15_y1"; }
};



class P_Cov_t1_y1 : public var::Var<P_Cov_t1_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t1_y1) { return "P_Cov_t1_y1"; }
};

class P_Cov_t20_y1 : public var::Var<P_Cov_t20_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t20_y1) { return "P_Cov_t20_y1"; }
};

class P_Cov_t11_y0 : public var::Var<P_Cov_t11_y0, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t11_y0) { return "P_Cov_t11_y0"; }
};

class P_Cov_t10_y1 : public var::Var<P_Cov_t10_y1, SymmetricMatrix<double>> {
   public:
    friend std::string className(P_Cov_t10_y1) { return "P_Cov_t10_y1"; }
};

class lambda : public var::Var<lambda, DiagonalMatrix<double>> {};

class V : public var::Var<V, Matrix<double>> {};
class W : public var::Var<W, Matrix<double>> {};
// Block partition of the spectrum: rows [begin, end) per block
class Blocks : public var::Constant<Blocks, Matrix<std::size_t>> {};

template <bool b>
class uses_variance_aproximation : public var::constexpr_Var<bool, uses_variance_aproximation, b> {
};
class uses_variance_aproximation_value
    : public var::constexpr_Var_value<bool, uses_variance_aproximation> {};

template <bool b>
class uses_taylor_variance_correction_aproximation
    : public var::constexpr_Var<bool, uses_taylor_variance_correction_aproximation, b> {};
class uses_taylor_variance_correction_aproximation_value
    : public var::constexpr_Var_value<bool, uses_taylor_variance_correction_aproximation> {};

template <bool b>
class uses_adaptive_aproximation : public var::constexpr_Var<bool, uses_adaptive_aproximation, b> {
};

class uses_adaptive_aproximation_value
    : public var::constexpr_Var_value<bool, uses_adaptive_aproximation> {};

template <bool b>
class uses_recursive_aproximation
    : public var::constexpr_Var<bool, uses_recursive_aproximation, b> {};
class uses_recursive_aproximation_value
    : public var::constexpr_Var_value<bool, uses_recursive_aproximation> {};

template <int b>
class uses_averaging_aproximation : public var::constexpr_Var<int, uses_averaging_aproximation, b> {
};
class uses_averaging_aproximation_value
    : public var::constexpr_Var_value<int, uses_averaging_aproximation> {};

template <int b>
class return_predictions : public var::constexpr_Var<int, return_predictions, b> {};
class return_predictions_value : public var::constexpr_Var_value<int, return_predictions> {};

template <class T>
concept uses_averaging_aproximation_c =
    var::is_this_constexpr_Var_c<T, int, uses_averaging_aproximation>;

template <class T>
concept uses_adaptive_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_adaptive_aproximation>;

template <class T>
concept uses_recursive_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_recursive_aproximation>;

template <class T>
concept what_to_include_c = is_of_this_template_type_v<T, var::please_include>;

template <class T>
concept uses_variance_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_variance_aproximation>;

template <class T>
concept uses_taylor_variance_correction_aproximation_c =
    var::is_this_constexpr_Var_c<T, bool, uses_taylor_variance_correction_aproximation>;

class Probability_error_tolerance : public var::Constant<Probability_error_tolerance, double> {};

class Conductance_variance_error_tolerance
    : public var::Constant<Conductance_variance_error_tolerance, double> {};

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

class P : public Var<P, Matrix<double>> {
    friend std::string className(const P&) { return "P_ij"; }
};

class P_half : public Var<P_half, Matrix<double>> {
    friend std::string className(const P_half&) { return "Ph_ij"; }
};

template <class Policy = StabilizerPolicyEnabled, class C_Matrix>
Maybe_error<Transfer_Op_to<C_Matrix, P>> to_Transition_Probability(C_Matrix const& x) {
    if constexpr (!Policy::project_transition_probability) {
        return build<P>(x);
    }
    constexpr double smooth_eps = 1e-12;
    auto out = apply(
        [smooth_eps](auto const& value) {
            using ValueType = std::decay_t<decltype(value)>;
            using std::sqrt;
            auto half = ValueType(0.5);
            auto eps_like = ValueType(smooth_eps);
            auto term = value * value + eps_like;
            auto sqrt_term = sqrt(term);
            return half * (value + sqrt_term);
        },
        x);
    auto sumP = out * Matrix<double>(out.ncols(), 1ul, 1.0);
    auto s = inv(diag(sumP));

    for (std::size_t i = 0; i < sumP.size(); ++i)
        if (std::isnan(primitive(sumP[i]))) {
            //  std::cerr << "rro";
            return error_message("not transition prob");
        }
    if (s)
        // auto test=s*out*Matrix<double>(out.ncols(),1ul, 1.0);
        return build<P>(s.value() * out);
    else
        return s.error();
}

class Ptotal_ij : public Var<Ptotal_ij, Matrix<double>> {
    friend std::string className(Ptotal_ij) { return "Ptotal_ij"; }

   public:
};
inline Maybe_error<Ptotal_ij> make_Ptotal_ij(Matrix<double>&& x, double max_dt) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        if ((x[i] < -max_dt) || (x[i] > max_dt + 1))
            return error_message(std::to_string(i) + "= " + std::to_string(x[i]) +
                                 " istoo big or small ");
        else
            x[i] = std::min(1.0, std::max(0.0, x[i]));
    }
    return Ptotal_ij(x / var::sum(x));
}

class gmean_i : public Var<gmean_i, Matrix<double>> {
    friend std::string className(gmean_i) { return "gmean_i"; }
};
class gtotal_ij : public Var<gtotal_ij, Matrix<double>> {
    friend std::string className(gtotal_ij) { return "gtotal_ij"; }
};
class gmean_ij : public Var<gmean_ij, Matrix<double>> {
    friend std::string className(gmean_ij) { return "gmean_ij"; }
};
class gtotal_sqr_ij : public Var<gtotal_sqr_ij, Matrix<double>> {
    friend std::string className(gtotal_sqr_ij) { return "gtotal_sqr_ij"; }
};
class gsqr_i : public Var<gsqr_i, Matrix<double>> {
    friend std::string className(gsqr_i) { return "gsqr_i"; }
};
class gvar_i : public Var<gvar_i, Matrix<double>> {
    friend std::string className(gvar_i) { return "gvar_i"; }
};
class gtotal_var_ij : public Var<gtotal_var_ij, Matrix<double>> {
    friend std::string className(gtotal_var_ij) { return "gtotal_var_ij"; }
};
class gvar_ij : public Var<gvar_ij, Matrix<double>> {
    friend std::string className(gvar_ij) { return "gvar_ij"; }
};

class y_mean : public var::Var<y_mean, double> {
    friend std::string className(y_mean) { return "y_mean"; }
};
class y_var : public var::Var<y_var, double> {
    friend std::string className(y_var) { return "y_var"; }
};
class trust_coefficient : public var::Var<trust_coefficient, double> {
    friend std::string className(trust_coefficient) { return "trust_coefficient"; }
};

class Chi2 : public var::Var<Chi2, double> {
    friend std::string className(Chi2) { return "Chi2"; }
};

// class plogL : public var::Var<plogL, double> {
//     friend std::string className(plogL) { return "plogL"; }
// };
// class eplogL : public var::Var<eplogL, double> {
//     friend std::string className(eplogL) { return "eplogL"; }
// };
// class vplogL : public var::Constant<vplogL, double> {
//     friend std::string className(vplogL) { return "vplogL"; }
// };

class macror_algorithm : public var::Constant<macror_algorithm, std::string> {
    using var::Constant<macror_algorithm, std::string>::Constant;
    friend std::string className(macror_algorithm) { return "macror_algorithm"; }
};

class PGn : public var::Var<PGn, Matrix<double>> {};
class PGG_n : public var::Var<PGG_n, Matrix<double>> {};
class PG_n : public var::Var<PG_n, Matrix<double>> {};
// class PPn : public var::Var<PPn, Matrix<double>> {};

using Qn = Vector_Space<number_of_samples, min_P, P, PG_n, PGG_n>;

using Eigs = Vector_Space<lambda, V, W>;

using Qdtg = Vector_Space<number_of_samples, min_P, P_half, g>;

using Qdtm =
    Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij, gsqr_i, gvar_i>;

using Qdt = Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij, gtotal_sqr_ij,
                         gsqr_i, gvar_i, gtotal_var_ij, gvar_ij>;

template <class recursive, class averaging, class variance, class variance_correction>

    requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
             uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
struct Qdt_u {
    using type = double;
};

using Patch_Model = Vector_Space<N_St, Q0, Qa, P_initial, g, N_Ch_mean, Current_Noise, Pink_Noise,
                                 Proportional_Noise, Current_Baseline,
                                 N_Ch_mean_time_segment_duration, Binomial_magical_number, min_P,
                                 Probability_error_tolerance, Conductance_variance_error_tolerance>;

inline void save(const std::string name, const Patch_Model& m) {
    std::ofstream f_Q0(name + "_Q0.txt");
    f_Q0 << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<Q0>(m) << "\n";
    std::ofstream f_Qa(name + "_Qa.txt");
    f_Qa << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<Qa>(m) << "\n";
    std::ofstream f_g(name + "_g.txt");
    f_g << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get<g>(m) << "\n";
}

 

class Algo_State_Dynamic
    : public var::Var<
          Algo_State_Dynamic,
          Vector_Space<y_mean, y_var, trust_coefficient, Chi2, P,P_half,gmean_i,gvar_i,gmean_ij,gtotal_ij,P_mean_t2_y0, P_mean_t2_y1,P_mean_t15_y0, P_mean_t15_y1,
                       P_mean_t1_y1, P_mean_t20_y1, P_mean_t11_y0, P_mean_t10_y1, 
                       P_mean_0t_y0,P_mean_0t_y1,P_var_ii_0t_y0,P_var_ii_0t_y1,
                       P_Cov_t2_y0,
                       P_Cov_t2_y1,P_Cov_t15_y0,
                       P_Cov_t15_y1, P_Cov_t1_y1, P_Cov_t20_y1, P_Cov_t11_y0, P_Cov_t10_y1>> {
   public:
    Matrix<double> const& get_P_mean() const {
        if (get<P_mean_t2_y1>((*this)())().size() > 0) {
            return get<P_mean_t2_y1>((*this)())();
        }
        if (get<P_mean_t2_y0>((*this)())().size() > 0) {
            return get<P_mean_t2_y0>((*this)())();
        }
        
        return get<P_mean_t20_y1>((*this)())();
    }
    SymmetricMatrix<double> const& get_P_Cov() const {
        if (get<P_Cov_t2_y1>((*this)())().size() > 0) {
            return get<P_Cov_t2_y1>((*this)())();
        }if (get<P_Cov_t2_y0>((*this)())().size() > 0) {
            return get<P_Cov_t2_y0>((*this)())();
        }
        
        return get<P_Cov_t20_y1>((*this)())();
    }
};

class Algo_State
    : public var::Var<Algo_State,
                      Vector_Space<y_mean, y_var, trust_coefficient, Chi2, P_mean, P_Cov>> {
   public:
    using base_type =
        var::Var<Algo_State, Vector_Space<y_mean, y_var, trust_coefficient, Chi2, P_mean, P_Cov>>;
    Algo_State(const Algo_State_Dynamic& p)
        : base_type{Vector_Space(get<y_mean>(p()), get<y_var>(p()), get<trust_coefficient>(p()),
                                 get<Chi2>(p()),

                                 P_mean(p.get_P_mean()), P_Cov(p.get_P_Cov()))} {}

    using base_type::Var;
};

struct Patch_State : public var::Var<Patch_State, Vector_Space<P_mean, P_Cov>> {};

struct Evolution {};

template <class T>
class Evolution_of : public Var<Evolution_of<T>, std::vector<T>> {
   public:
    using element_type = T;
    auto& operator[](var::Var<Evolution>) { return *this; }

    auto const& operator[](var::Var<Evolution>) const { return *this; }

    friend std::string className(Evolution_of) { return "Macro_State_Evolution"; }
    using value_type = std::vector<T>;
};

template <class V, class Id>
concept has_var_c = requires(V&& v) {
    std::forward<V>(v)[var::Var<Id>{}];
};

template <typename... Vars>
struct Macro_State : public Vector_Space<logL, Patch_State, Vars...> {
    Macro_State() = default;
    Macro_State(Patch_State&& ps) { get<Patch_State>((*this)) = std::move(ps); }
    Macro_State(logL&& l, Patch_State&& ps, Vars&&... vars)
        : Vector_Space<logL, Patch_State, Vars...>(std::move(l), std::move(ps),
                                                   std::forward<Vars>(vars)...) {}
};

template <typename... Vars>
struct dMacro_State
    : public Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                          var::Derivative<Patch_State, var::Parameters_transformed>, Vars...> {
    dMacro_State(var::Derivative<Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);

        // Seed the prior patch state (carries dx).
        get<var::Derivative<Patch_State, var::Parameters_transformed>>(*this) = std::move(dps);
        // Accumulators start at zero/empty but share the same dx.
        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<var::Derivative<logL, var::Parameters_transformed>>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    dMacro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                 var::Derivative<Patch_State, var::Parameters_transformed>&& dps, Vars&&... vars)
        : Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                       var::Derivative<Patch_State, var::Parameters_transformed>, Vars...>(
              std::move(dl), std::move(dps), std::forward<Vars>(vars)...) {}
    dMacro_State() = default;   
    
                        
    
};

template <typename... Vars>
struct ddMacro_State
    : public var::Derivative<Vector_Space<logL, Patch_State, Vars...>, var::Parameters_transformed> {
    ddMacro_State(var::Derivative<Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);
        get<Patch_State>(*this) = std::move(dps);

        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<logL>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    ddMacro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                  var::Derivative<Patch_State, var::Parameters_transformed>&& dps,
                  var::Derivative_t<Vars, var::Parameters_transformed>&&... vars)
        : var::Derivative<Vector_Space<logL, Patch_State, Vars...>, var::Parameters_transformed>(
              std::move(dl), std::move(dps), std::forward<Vars>(vars)...) {}
};

}  // namespace macrodr

// Teach the generic derivative machinery that ddMacro_State lives in the
// Parameters_transformed derivative universe, so that Transfer_Op_to and
// related helpers propagate derivatives correctly.
namespace var {
template <class... Vars>
struct transformation_type<macrodr::ddMacro_State<Vars...>> {
    using type = Derivative_Op<Parameters_transformed>;
};

template <class... Vars>
struct is_derivative<macrodr::ddMacro_State<Vars...>> : std::true_type {};

template <class... Vars, class... Ds>
struct dx_of_dfdx<macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class G, class... Vars, class... Ds>
    requires(!is_derivative_v<G>)
struct dx_of_dfdx<G, macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class F, class... Vars, class... Ds>
struct dx_of_dfdx<Derivative<F, Parameters_transformed>, macrodr::ddMacro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};
}  // namespace var

namespace macrodr {

using predictions_element =
    var::please_include<logL, elogL, vlogL, y_mean, y_var, P_mean, P_Cov, trust_coefficient>;

using diagnostic_element = var::please_include<logL, elogL, vlogL, Algo_State_Dynamic>;

using gradient_minimal_element =
    var::please_include<var::Derivative<logL, var::Parameters_transformed>, elogL, y_mean, y_var,
                        trust_coefficient>;

using gradient_all_element =
    var::please_include<var::Derivative<logL, var::Parameters_transformed>,
                        var::Derivative<elogL, var::Parameters_transformed>,
                        var::Derivative<y_mean, var::Parameters_transformed>,
                        var::Derivative<y_var, var::Parameters_transformed>, trust_coefficient>;

using Macro_State_minimal = Macro_State<>;

using Macro_State_reg = add_t<Macro_State_minimal, var::please_include<elogL, vlogL>>;

using dMacro_State_Hessian_minimal =
    add_t<dMacro_State<>, var::please_include<FIM>>;

using diff_Macro_State_Gradient_Hessian =
    add_t<Macro_State<>, var::please_include<elogL, vlogL, Grad, FIM>>;

using Macro_State_Ev_predictions =
    add_t<Macro_State_reg,
          var::please_include<Evolution_of<add_t<Vector_Space<>, predictions_element>>>>;

using Macro_State_Ev_diagnostic =
    add_t<Macro_State_reg,
          var::please_include<Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>>;

using dMacro_State_Ev_gradient_minimal =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, gradient_minimal_element>>>>;

using dMacro_State_Ev_gradient_all =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>>;

template <class VS>
constexpr bool is_Algo_dynamic() {
    if constexpr (has_var_c<VS const&, Evolution>) {
        using Evo = std::decay_t<decltype(std::declval<VS const&>()[var::Var<Evolution>{}])>;
        using El = typename Evo::element_type;
        return has_var_c<El const&, Algo_State_Dynamic>;
    } else {
        return false;
    }
}

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

class Simulation_n_sub_dt : public Var<Simulation_n_sub_dt, std::size_t> {};

class N_Ch_State_Evolution : public Var<N_Ch_State_Evolution, std::vector<N_channel_state>> {};

class Only_Ch_Curent_Evolution : public Var<Only_Ch_Curent_Evolution, std::vector<Patch_current>> {
};

class N_Ch_State_Sub_Evolution
    : public Var<N_Ch_State_Sub_Evolution, std::vector<N_channel_state>> {};

class Only_Ch_Curent_Sub_Evolution
    : public Var<Only_Ch_Curent_Sub_Evolution, std::vector<Patch_current>> {};

template <typename Simulate_tag>
class Simulated_Recording
    : public Var<Simulated_Recording<Simulate_tag>, add_t<Vector_Space<Recording>, Simulate_tag>> {
   public:
    constexpr static const bool includes_N = var::has_it_v<Simulate_tag, N_Ch_State_Evolution>;
    using Var<Simulated_Recording<Simulate_tag>, add_t<Vector_Space<Recording>, Simulate_tag>>::Var;
};

using Simulated_recording = Simulated_Recording<var::please_include<>>;

using v_Simulated_Recording =
    std::variant<Simulated_Recording<var::please_include<>>,
                 Simulated_Recording<var::please_include<Only_Ch_Curent_Evolution>>,
                 Simulated_Recording<var::please_include<N_Ch_State_Evolution>>>;

template <typename Simulate_tag>
class Simulated_Step
    : public Var<Simulated_Step<Simulate_tag>,
                 Vector_Space<N_channel_state, Simulated_Recording<Simulate_tag>>> {
    // using Vector_Space<N_channel_state,
    // Simulated_Recording<Simulate_tag>>::Vector_Space;
};

template <typename Simulate_tag>
using Simulated_Sub_Step_t =
    add_if_present_t<Vector_Space<N_channel_state, number_of_samples, y_sum>, Simulate_tag,
                     N_Ch_State_Sub_Evolution, Only_Ch_Curent_Sub_Evolution>;

template <typename Simulate_tag>
Simulated_Sub_Step_t<Simulate_tag> Simulated_Sub_Step_build(N_channel_state N) {
    Simulated_Sub_Step_t<Simulate_tag> out;
    get<N_channel_state>(out) = N;
    return out;
}

using Simulation_Parameters = Vector_Space<Simulation_n_sub_dt>;

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

template <typename uses_recursive, typename uses_averaging, typename uses_variance>
    requires(var::is_this_constexpr_Var_v<uses_recursive, bool, uses_recursive_aproximation> &&
             var::is_this_constexpr_Var_v<uses_averaging, int, uses_averaging_aproximation> &&
             var::is_this_constexpr_Var_v<uses_variance, bool, uses_variance_aproximation>)
struct MacroR {
    friend std::string ToString(MacroR) {
        std::string out = "MacroR";

        out += uses_recursive::value ? "_R" : "_NR";
        out += (uses_averaging::value == 2) ? "_2" : "__";
        out += uses_variance::value ? "_V" : "_M";

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
        return build<Eigs>(build<lambda>(std::move(r_lambda)), build<V>(std::move(r_V)),
                           build<W>(std::move(r_W)));
    }

    template <class C_Patch_Model>
        requires U<C_Patch_Model, Patch_Model>
    auto calc_eigen(const C_Patch_Model& m, Agonist_concentration x)
        -> Maybe_error<Transfer_Op_to<C_Patch_Model, Eigs>> {
        return calc_eigen(calc_Qx(m, x));
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
            auto n = std::binomial_distribution<std::size_t>(
                N_remaining, cumP[i] > 0 ? t_P_mean()[i] / cumP[i] : 0)(mt);
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
                    auto n = std::binomial_distribution<std::size_t>(N_remaining,
                                                                     t_P()(i, j) / cumP(i, j))(mt);
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
        // const double eps=std::numeric_limits<double>::epsilon();
        auto& t_V = get<V>(t_Qx);
        auto& t_landa = get<lambda>(t_Qx);
        auto& t_W = get<W>(t_Qx);
        auto& t_g = get<g>(m);
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
        if constexpr (!Policy::enforce_gmean_bounds) {
            t_min_P() = eps;
        }
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
        else {
            auto r_P = std::move(Maybe_r_P.value());

            std::size_t N = r_P().ncols();

            SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < i + 1; ++j) {
                    set(E2m, i, j,
                        Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j], t_min_P()));
                }
            }

            SymmetricMatrix<Op_t<Trans, double>> E2_sqr(N, N);
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < i + 1; ++j) {
                    set(E2_sqr, i, j, E2(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j]));
                }
            }

            Matrix<Op_t<Trans, double>> WgV_E2(N, N);

            auto v_WgV = var::inside_out(t_W() * diag(t_g()) * t_V());

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < N; ++j) {
                    WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);
                }
            }

            auto t_WgV_E2 = var::outside_in(std::move(WgV_E2), dx);
            auto r_gtotal_ij =
                force_gtotal_in_range<Policy>(build<gtotal_ij>(t_V() * t_WgV_E2 * t_W()), t_g, r_P);

            auto r_gmean_ij = build<gmean_ij>(elemDivSafe(r_gtotal_ij(), r_P(), t_min_P()));
            /* truncate is not derivative safe yet*/

            Matrix<double> u(N, 1, 1.0);
            auto r_gmean_i = force_gmean_in_range(build<gmean_i>(r_gtotal_ij() * u), t_g);
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

            auto v_Wg = t_W() * t_g();

            auto v_eps = eps;
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j)
                    WgV_Wg_E2(i, 0) = WgV_Wg_E2(i, 0) + v_WgV(i, j) * E2_sqr(i, j) * v_Wg[j];

            Matrix<Op_t<Trans, double>> rgsqr_i(N, 1, 0.0);

            for (std::size_t i = 0; i < N; i++) {
                for (std::size_t k0 = 0; k0 < N; k0++)
                    rgsqr_i[i] = rgsqr_i[i] + 2 * t_V()(i, k0) * WgV_Wg_E2[k0];
            }

            auto r_gsqr_i = build<gsqr_i>(var::outside_in(rgsqr_i, dx));

            auto r_gvar_i = build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

            /* truncate is not derivative safe yet*/

            if constexpr (StabilizerPolicyEnabled::clamp_variance) {
                r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
            }
            //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
            //   primitive(t_g())); auto test_g_mean =
            //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
            // if (!test_g_mean || !test_g_var)
            //     return error_message(test_g_var.error()()+test_g_mean.error()());
            // else {

            if (std::isnan(primitive(var::max(r_P()))))
                return error_message("nan P");

            return build<Qdtm>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                               std::move(r_gtotal_ij), std::move(r_gmean_ij), std::move(r_gsqr_i),
                               std::move(r_gvar_i));
            // }
        }
    }

    template <class Policy = StabilizerPolicyEnabled, class FunctionTable, class C_Patch_Model,
              class C_Qx_eig>
        requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Eigs>)
    Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> calc_Qdt_eig(FunctionTable&&,
                                                                 const C_Patch_Model& m,
                                                                 const C_Qx_eig& t_Qx,
                                                                 number_of_samples ns, double dt) {
        using Trans = transformation_type_t<C_Patch_Model>;
        // const double eps=std::numeric_limits<double>::epsilon();
        auto& t_V = get<V>(t_Qx);
        auto& t_landa = get<lambda>(t_Qx);
        auto& t_W = get<W>(t_Qx);
        auto& t_g = get<g>(m);
        auto t_min_P = get<min_P>(m);
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
        else {
            auto r_P = std::move(Maybe_r_P.value());

            std::size_t N = r_P().ncols();

            SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < i + 1; ++j) {
                    set(E2m, i, j,
                        Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j], t_min_P()));
                }
            }

            Matrix<Op_t<Trans, double>> WgV_E2(N, N);

            auto v_WgV = t_W() * diag(t_g()) * t_V();

            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j) WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

            auto r_gtotal_ij =
                force_gtotal_in_range<Policy>(build<gtotal_ij>(t_V() * WgV_E2 * t_W()), t_g, r_P);

            Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
            for (std::size_t n1 = 0; n1 < N; n1++)
                for (std::size_t n3 = 0; n3 < N; n3++)
                    for (std::size_t n2 = 0; n2 < N; n2++) {
                        //      std::cerr<<"\t"<<WgV_E3(n1, n3);

                        WgV_E3(n1, n3) = WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                                              E3(v_ladt[n1], v_ladt[n2], v_ladt[n3],
                                                                 v_exp_ladt[n1], v_exp_ladt[n2],
                                                                 v_exp_ladt[n3], t_min_P());
                    }

            auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);

            auto r_gmean_ij = build<gmean_ij>(elemDivSafe(r_gtotal_ij(), r_P(), eps));

            auto r_gtotal_var_ij = force_gtotal_var_in_range<Policy>(
                build<gtotal_var_ij>(r_gtotal_sqr_ij() - elemMult(r_gtotal_ij(), r_gmean_ij())),
                t_g, r_P);

            /* truncate is not derivative safe yet*/

            auto r_gvar_ij = build<gvar_ij>(elemDivSafe(r_gtotal_var_ij(), r_P(), eps));

            Matrix<double> u(N, 1, 1.0);
            auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
            if (crude_gmean_violation(primitive(r_gmean_i), primitive(get<g>(m))))
                return error_message("gmean_violation");

            auto r_gsqr_i = build<gsqr_i>(r_gtotal_sqr_ij() * u);
            auto r_gvar_i = build<gvar_i>(r_gtotal_var_ij() * u);

            //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
            //   primitive(t_g())); auto test_g_mean =
            //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
            // if (!test_g_mean || !test_g_var)
            //     return error_message(test_g_var.error()()+test_g_mean.error()());
            // else {

            return build<Qdt>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                              std::move(r_gtotal_ij), std::move(r_gmean_ij),
                              std::move(r_gtotal_sqr_ij), std::move(r_gsqr_i), std::move(r_gvar_i),
                              std::move(r_gtotal_var_ij), std::move(r_gvar_ij));
            // }
        }
    }

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

        auto r_gmean_ij = [&]() {
            auto base = build<gmean_ij>(elemDivSafe(r_gtotal_ij(), r_P(), t_min_P()));
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

        auto r_gmean_ij = [&]() {
            auto base = build<gmean_ij>(elemDivSoftAbs(r_gtotal_ij(), r_P(), t_min_P()));
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

        auto r_gvar_ij = build<gvar_ij>(elemDivSoftAbs(r_gtotal_var_ij(), r_P(), t_min_P()));

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
            auto r_Qn = get_Qn(P_sub, get<g>(m), sub_ns, get<min_P>(m));
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
        } else {
            return build<Qdtg>(ns, get<min_P>(m), build<P_half>(std::move(Maybe_P_sub.value()())), get<g>(m));
        }
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
            auto r_Qn = get_Qn(P_sub, get<g>(m), sub_ns, get<min_P>(m));
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
        if (!((max_g_m <= max_g) && (min_g_m >= min_g))) {
            std::cerr << "max_g_m=" << max_g_m << " max_g=" << max_g << " min_g_m=" << min_g_m
                      << " min_g=" << min_g << "\n";
            return true;
        } else
            return false;
    }

    template <class Policy = StabilizerPolicyEnabled, class C_gmean, class C_g>
        requires((U<C_gmean, gmean_i> || U<C_gmean, gtotal_ij> || U<C_gmean, gmean_ij>) &&
                 U<C_g, g>)
    static auto force_gmean_in_range(C_gmean&& g_mean, const C_g& v_g) {
        if constexpr (StabilizerPolicyEnabled::enforce_gmean_bounds) {
            auto gmax = var::max(v_g());

            auto gmin = var::min(v_g());
            g_mean() = apply(
                [&gmax, &gmin](auto const& value) {
                    using std::max;
                    using var::max;
                    using std::min;
                    using var::min;
                    return max(min(value, gmax), gmin);
                },
                std::forward<C_gmean>(g_mean)());
        }
        return std::forward<C_gmean>(g_mean);
    }

    template <class Policy = StabilizerPolicyEnabled, class C_gtotal, class C_g, class C_P>
        requires((U<C_gtotal, gtotal_ij>) && (U<C_P, P>) && U<C_g, g>)
    static auto force_gtotal_in_range(C_gtotal&& g_total, const C_g& v_g, const C_P& v_P) {
        if constexpr (StabilizerPolicyEnabled::enforce_gmean_bounds) {
            using std::max;
            using std::min;
            using var::max;
            using var::min;
            auto gmax = var::max(v_g());

            auto gmin = var::min(v_g());
            g_total() = zip(
                [&gmax, &gmin](auto const& value, auto const& Pij) {
                    return max(min(value, gmax * Pij), gmin * Pij);
                },
                g_total(), v_P());
        }
        return std::forward<C_gtotal>(g_total);
    }

    template <class Policy = StabilizerPolicyEnabled, class C_gtotal_var_ij, class C_g, class C_P>
        requires((U<C_gtotal_var_ij, gtotal_var_ij>) && (U<C_P, P>) && U<C_g, g>)
    static auto force_gtotal_var_in_range(C_gtotal_var_ij&& g_totalvar, const C_g& v_g,
                                          const C_P& v_P) {
        if constexpr (StabilizerPolicyEnabled::enforce_gmean_bounds) {
            using std::max;
            using std::min;
            using var::max;
            using var::min;
            auto gmax = var::max(v_g());
            auto gmin = var::min(v_g());
            auto gmaxvar = (gmax - gmin) / 2;
            auto gminvar = gmin - gmin;
            g_totalvar() = zip(
                [&gmaxvar, &gminvar](auto const& value, auto const& Pij) {
                    return max(min(value, gmaxvar * Pij), gminvar * Pij);
                },
                g_totalvar(), v_P());
        }
        return std::forward<C_gtotal_var_ij>(g_totalvar);
    }

    template <class C_Qn, class C_Patch_Model>
        requires(U<C_Qn, Qn>)
    static auto Qn_to_Qdt(const C_Qn& x, const C_Patch_Model& m) {
        auto u = Matrix<double>(get<P>(x)().ncols(), 1ul, 1.0);
        auto r_P = get<P>(x)();
        auto n = get<number_of_samples>(x)();
        auto r_gtotal_sqr_ij = get<PGG_n>(x)() * (2.0 / (n * n));
        auto r_gtotal_ij = get<PG_n>(x)() * (1.0 / n);
        auto b_gtotal_ij = force_gmean_in_range(build<gtotal_ij>(r_gtotal_ij), get<g>(m));

        auto r_gmean_ij = elemDivSoftAbs(b_gtotal_ij(), get<P>(x)(), get<min_P>(x)());
        auto b_gmean_ij = force_gmean_in_range(build<gmean_ij>(r_gmean_ij), get<g>(m));

        auto r_gtotal_var_ij = r_gtotal_sqr_ij - elemMult(b_gtotal_ij(), b_gmean_ij());
        auto r_gmean_i = r_gtotal_ij * u;
        auto b_gmean_i = force_gmean_in_range(build<gmean_i>(r_gmean_i), get<g>(m));
        auto r_gsqr_i = r_gtotal_sqr_ij * u;
        auto r_gvar_ij = elemDivSoftAbs(r_gtotal_var_ij, r_P, get<min_P>(x)());
        auto r_gvar_i = r_gtotal_var_ij * u;

        return build<Qdt>(get<number_of_samples>(x), get<min_P>(x), get<P>(x), std::move(b_gmean_i),
                          std::move(b_gtotal_ij), std::move(b_gmean_ij),
                          build<gtotal_sqr_ij>(r_gtotal_sqr_ij), build<gsqr_i>(r_gsqr_i),
                          build<gvar_i>(r_gvar_i), build<gtotal_var_ij>(r_gtotal_var_ij),
                          build<gvar_ij>(r_gvar_ij));
    }

    template <class C_Qn, class C_Patch_Model>
        requires(U<C_Qn, Qn>)
    static auto Qn_to_Qdtm(const C_Qn& x, const C_Patch_Model& m) {
        auto u = Matrix<double>(get<P>(x)().ncols(), 1ul, 1.0);
        auto r_P = get<P>(x)();
        auto n = get<number_of_samples>(x)();
        if (n <= 0)
            std::cerr << " nana here";

        auto r_gtotal_sqr_ij = get<PGG_n>(x)() * (2.0 / (n * n));
        auto r_gtotal_ij = get<PG_n>(x)() * (1.0 / n);
        auto b_gtotal_ij = force_gmean_in_range(build<gtotal_ij>(r_gtotal_ij), get<g>(m));

        auto r_gmean_ij = elemDivSoftAbs(b_gtotal_ij(), get<P>(x)(), get<min_P>(x)());
        auto b_gmean_ij = force_gmean_in_range(build<gmean_ij>(r_gmean_ij), get<g>(m));
        auto r_gtotal_var_ij = r_gtotal_sqr_ij - elemMult(b_gtotal_ij(), b_gmean_ij());
        auto r_gmean_i = r_gtotal_ij * u;
        auto b_gmean_i = force_gmean_in_range(build<gmean_i>(r_gmean_i), get<g>(m));

        auto r_gsqr_i = r_gtotal_sqr_ij * u;
        auto r_gvar_i = r_gtotal_var_ij * u;

        auto out = build<Qdtm>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                               std::move(b_gmean_i), std::move(b_gtotal_ij), std::move(b_gmean_ij),
                               build<gsqr_i>(r_gsqr_i), build<gvar_i>(r_gvar_i));
        // if (!is_finite(out))
        //   std::cerr<<" nana here";
        return out;
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

    template <class C_y_var>
    auto calculate_logL(bool y_is_nan, C_y_var const& r_y_var, auto const& chi2,
                        auto const& alfa, auto& m) const -> Transfer_Op_to<C_y_var, logL> {
        using DX = var::dx_of_dfdx_t<C_y_var>;
        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<logL>(std::move(base));
        }
        if (get<Proportional_Noise>(m).value() == 0) {
            return build<logL>(-0.5 * log(2 * std::numbers::pi * r_y_var() * alfa()) -
                               0.5 * chi2());
        } else {
            return build<logL>(0.5 * chi2() - log(var::Poisson_noise_normalization(
                                                  primitive(r_y_var() * alfa()),
                                                  primitive(get<Proportional_Noise>(m).value()))));
        }
    }

	    template <class C_y_var>
	    auto calculate_elogL(bool y_is_nan, C_y_var const& r_y_var, auto const& alfa, auto& m) const {
	        using DX = var::dx_of_dfdx_t<C_y_var>;
	        auto const& dx = var::get_dx_of_dfdx(r_y_var);

        if (y_is_nan) {
            auto base = var::init_with_dx<DX>(0.0, dx);
            return build<elogL>(std::move(base));
        }
	        if (get<Proportional_Noise>(m).value() == 0) {
	            return build<elogL>(-0.5 * log(2 * std::numbers::pi * r_y_var()*alfa()) - 0.5);
	        } else {
	            const auto pn_value = get<Proportional_Noise>(m).value();
	            if constexpr (var::is_derivative_v<std::decay_t<decltype(pn_value)>>) {
	                return build<elogL>(var::Poisson_noise_expected_logL(r_y_var() * alfa(),
	                                                                    pn_value));
	            } else {
	                auto pn = var::init_with_dx<DX>(pn_value, dx);
	                return build<elogL>(
	                    var::Poisson_noise_expected_logL(r_y_var() * alfa(), pn));
	            }
	        }
	    }

    auto calculate_trust_coefficient(Matrix<double> const& t_pmean, Matrix<double> const& d,
                                     double factor) const {
        auto alfa = 1.0;
        for (std::size_t i = 0; i < d.size(); ++i) {
            double d_i = d[i];
            double p_i = t_pmean[i];
            if (d_i > 0) {
                double alfa_i = (1.0 - p_i) / d_i;
                alfa = std::min(alfa_i, alfa);
            } else if (d_i < 0) {
                double alfa_i = -(p_i) / d_i;
                alfa = std::min(alfa_i, alfa);
            }
        }
        if (alfa < 1)
            alfa = alfa * factor;
        return trust_coefficient(alfa);
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
            p_P_mean() = p_P_mean() * t_P();
            p_P_Cov() = AT_B_A(t_P(), p_P_Cov() - diag(p_P_mean())) + diag(p_P_mean() * t_P());
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
            get<trust_coefficient>(out()) = alfa;

            get<Chi2>(out()) = std::move(chi2);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())() = std::move(r_P_cov());
            return out;
        } else {
            Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;

            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<trust_coefficient>(out()) = alfa;

            get<Chi2>(out()) = std::move(chi2);
            get<P_mean_t2_y0>(out())() = std::move(r_P_mean());
            get<P_Cov_t2_y0>(out())() = std::move(r_P_cov());
            return out;
        }
    }

    template <bool dynamic, class averaging, class variance, class C_Patch_State, class C_Qdt,
              class C_Patch_Model, class C_double>

        requires(uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> && U<C_Patch_State, Patch_State>)
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
        auto p_P_mean =[&t_prior, &t_P]()->decltype(auto){
            if constexpr(averaging::value>0){
                return get<P_mean>(t_prior());
            }
            else{
            return build<P_mean>(to_Probability(get<P_mean>(t_prior())()* t_P()).value());
            }   
        }(); 
        auto p_P_Cov =[&t_prior, &t_P, &p_P_mean]()->decltype(auto){
            if constexpr(averaging::value>0){
                return get<P_Cov>(t_prior());
            }
            else{
            return build<P_Cov>(AT_B_A(t_P(), get<P_Cov>(t_prior())() -diag(get<P_mean>(t_prior())())) + diag(p_P_mean()));
            }   
        }(); 
        
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
            auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
            if (std::isfinite(primitive(ms)) && primitive(ms) >= 0) {
                r_y_var() = r_y_var() + N * ms;
            } else {
                return error_message("invalid channel noise", ms);
            }
        }

        auto dy = y - r_y_mean();
        auto chi = dy / r_y_var();
        auto chi2 = build<Chi2>(dy * chi);
        auto gS = [&t_gmean_i, &p_P_Cov,&SmD, &p_P_mean, &t_Qdt, &t_P]() {
            if constexpr (averaging::value == 2) {
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);

                return TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
                ;
            } else {
                return TranspMult(t_gmean_i(), p_P_Cov());
            }
        }();
        auto alfa = [&](){
            if constexpr (averaging::value==2)  {
            return  calculate_trust_coefficient(primitive(p_P_mean()*t_P()), primitive(chi) * primitive(gS),
                                                trust_multiplying_factor);
                                            
            }
            else {
            return  calculate_trust_coefficient(primitive(p_P_mean()), primitive(chi) * primitive(gS),
                                                trust_multiplying_factor);
            }}();


        auto Maybe_r_P_mean = [&](){ 
            if constexpr (averaging::value==2){
            return to_Probability( p_P_mean() * t_P() + alfa() * chi * gS);}
            else{
            return to_Probability((p_P_mean() + alfa() * chi * gS) * t_P());
        }
        }();
        if (!Maybe_r_P_mean)
            return Maybe_r_P_mean.error();
        auto r_P_mean = build<P_mean>(std::move(Maybe_r_P_mean.value()));

        auto Maybe_r_P_cov = to_Covariance_Probability(
            AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) - (alfa() * N / r_y_var()) * XTX(gS));
        if (!Maybe_r_P_cov)
            return Maybe_r_P_cov.error();

        auto r_P_cov = build<P_Cov>(std::move(Maybe_r_P_cov.value()));

        if (!all_Probability_elements(primitive(r_P_mean())) ||
            !all_Covariance_elements(primitive(r_P_cov()))) {
            return error_message("error in P_mean or P_cov");
        }
        auto r_logL = calculate_logL(false, r_y_var, chi2, alfa, m);

        auto r_elogL = calculate_elogL(false, r_y_var, alfa,  m);

        if constexpr (!dynamic) {
            Transfer_Op_to<C_Patch_State, Algo_State> out;
            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<Chi2>(out())= std::move(chi2);
            get<P_mean>(out())() = std::move(r_P_mean());
            get<P_Cov>(out())() = std::move(r_P_cov());
            get<trust_coefficient>(out()) = alfa;
            return out;
        } else {
            Transfer_Op_to<C_Patch_State, Algo_State_Dynamic> out;
            get<y_mean>(out()) = std::move(r_y_mean);
            get<y_var>(out()) = std::move(r_y_var);
            get<Chi2>(out()) = std::move(chi2);
            get<trust_coefficient>(out()) = alfa;
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

                auto Maybe_r_P_mean_t11_y0 = to_Probability(p_P_mean() * t_P());
                auto Maybe_r_P_mean_t10_y1 = to_Probability(p_P_mean() + alfa() * chi * gS0);

                if (!Maybe_r_P_mean_t11_y0) {
                    return Maybe_r_P_mean_t11_y0.error();
                }
                if (!Maybe_r_P_mean_t10_y1) {
                    return Maybe_r_P_mean_t10_y1.error();
                }

                
                
                auto r_P_mean_0t_y0= diag(p_P_mean())*t_P();
                auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);

                auto GS=  diag(TranspMult(t_gmean_i(), SmD)) * t_P() + diag(p_P_mean())* t_gtotal_ij();
                
                auto r_P_mean_0t_y1= r_P_mean_0t_y0 + alfa()*chi*GS;

                auto r_P_var_ii_0t_y0= diag(p_P_Cov())*elemMult(t_P(),t_P())+diag(p_P_mean())*(t_P()-elemMult(t_P(),t_P()));
                auto r_P_var_ii_0t_y1= r_P_var_ii_0t_y0 - (alfa()*N/r_y_var())* elemMult(GS,GS);
                Matrix<double> uT(1ul,p_P_mean().size(), 1.0);
         
               
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y0).value(), Maybe_r_P_mean_t11_y0.value()));
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y1).value(), r_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y0)).value(), p_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y1)).value(), Maybe_r_P_mean_t10_y1.value()));

                get<P_mean_t11_y0>(out())() = std::move(Maybe_r_P_mean_t11_y0.value());
                get<P_mean_t10_y1>(out())() = std::move(Maybe_r_P_mean_t10_y1.value());
                
                get<P_mean_t20_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t20_y1>(out())() = std::move(r_P_cov());
                get<gtotal_ij>(out()) = get<gtotal_ij>(t_Qdt);
                get<gmean_ij>(out()) = get<gmean_ij>(t_Qdt);

                get<P_mean_0t_y0>(out())()= std::move(r_P_mean_0t_y0);
                get<P_mean_0t_y1>(out())()= std::move(r_P_mean_0t_y1);
                get<P_var_ii_0t_y0>(out())()= std::move(r_P_var_ii_0t_y0);
                get<P_var_ii_0t_y1>(out())()= std::move(r_P_var_ii_0t_y1);  


                
                
                auto Maybe_r_P_cov_t11_y0 =
                    to_Covariance_Probability(AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()));
                auto Maybe_r_P_cov_t10_y1 = to_Covariance_Probability(
                    get<P_Cov>(t_prior())() - (alfa() * N / r_y_var()) * XTX(gS0));

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
                    get<P_Cov>(t_prior())() - (alfa() * N / r_y_var()) * XTX(gS));

                if (!Maybe_r_P_cov_t1_y1) {
                    return Maybe_r_P_cov_t1_y1.error();
                }
                
                auto r_P_mean_0t_y0= diag(p_P_mean())*t_P();
                
                auto r_P_mean_0t_y1= diag(Maybe_r_P_mean_t1_y1.value())*t_P();

                auto r_P_var_ii_0t_y0= diag(p_P_Cov())*elemMult(t_P(),t_P())+
                diag(p_P_mean())*(t_P()-elemMult(t_P(),t_P()));

                auto r_P_var_ii_0t_y1= diag(Maybe_r_P_cov_t1_y1.value())*elemMult(t_P(),t_P())+
                diag(Maybe_r_P_mean_t1_y1.value())*(t_P()-elemMult(t_P(),t_P()));


                Matrix<double> uT(1UL,p_P_mean().size(), 1.0);
         
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y0).value(), Maybe_r_P_mean_t2_y0.value()));
                assert(var::test_equality(to_Probability(uT*r_P_mean_0t_y1).value(), r_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y0)).value(), p_P_mean()));
                assert(var::test_equality(to_Probability(MultTransp(uT,r_P_mean_0t_y1)).value(), Maybe_r_P_mean_t1_y1.value()));
                get<P_mean_0t_y0>(out())()= std::move(r_P_mean_0t_y0);
                get<P_mean_0t_y1>(out())()= std::move(r_P_mean_0t_y1);
                get<P_var_ii_0t_y0>(out())()= std::move(r_P_var_ii_0t_y0);
                get<P_var_ii_0t_y1>(out())()= std::move(r_P_var_ii_0t_y1);  


                get<P_Cov_t2_y0>(out())() = std::move(Maybe_r_P_cov_t2_y0.value());    
                get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                get<P_mean_t2_y0>(out())() = std::move(Maybe_r_P_mean_t2_y0.value());
                get<P_mean_t1_y1>(out())() = std::move(Maybe_r_P_mean_t1_y1.value());



                get<P_Cov_t1_y1>(out())() = std::move(Maybe_r_P_cov_t1_y1.value());
            }
            else{
                static_assert(averaging::value == 0);
                get<P_mean_t2_y1>(out())() = std::move(r_P_mean());
                get<P_Cov_t2_y1>(out())() = std::move(r_P_cov());
                auto gS0 = TranspMult(t_gmean_i(), p_P_Cov());

                auto Maybe_r_P_mean_t15_y1 = to_Probability(p_P_mean() + chi * gS0);

                if (!Maybe_r_P_mean_t15_y1) {
                    return Maybe_r_P_mean_t15_y1.error();
                }

                get<P_mean_t15_y1>(out())() = std::move(Maybe_r_P_mean_t15_y1.value());
                auto Maybe_r_P_cov_t15_y1 = to_Covariance_Probability(
                    p_P_Cov() - (alfa() * N / r_y_var()) * XTX(gS0));

                if (!Maybe_r_P_cov_t15_y1) {
                    return Maybe_r_P_cov_t15_y1.error();
                }

                get<P_Cov_t15_y1>(out())() = std::move(Maybe_r_P_cov_t15_y1.value());
                get<P_mean_t15_y0>(out())() = std::move(p_P_mean());
                get<P_Cov_t15_y0>(out())() = std::move(p_P_Cov());

                

            }
            return out;
        }
    }

    template <bool dynamic, class recursive, class averaging, class variance, class C_Patch_State,
              class C_Qdt, class C_Patch_Model, class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> && (U<C_Patch_State, Patch_State>))
    auto safely_calculate_Algo_State(C_Patch_State const& t_prior, C_Qdt const& t_Qdt,
                                     C_Patch_Model const& m, C_double const& N,
                                     const Patch_current& p_y, double fs) const {
        if constexpr (!recursive::value) {
            return safely_calculate_Algo_State_non_recursive<dynamic, averaging, variance>(
                t_prior, t_Qdt, m, N, p_y, fs);
        }
        return safely_calculate_Algo_State_recursive<dynamic, averaging, variance>(t_prior, t_Qdt,
                                                                                   m, N, p_y, fs);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<Macro_State<vVars...>> update(Macro_State<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, logL const& t_logL, elogL const& t_elogL ) const 
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
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
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
    Maybe_error<Vector_Space<vVars...>> update(Vector_Space<vVars...>&& t_prior_all,
                                              C_Algo_State&& algo, C_logL const& t_logL,C_elogL const& t_elogL) const {
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
           var::Derivative<logL, var::Parameters_transformed> const& t_logL, C_elogL const& t_elogL) const {
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
    Maybe_error<dMacro_State<vVars...>> update(
        dMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL, C_elogL const& t_elogL) const {
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
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL) const {
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
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class... vVars, class C_Algo_State, class C_elogL>
        requires(U<C_Algo_State, Algo_State> || U<C_Algo_State, Algo_State_Dynamic>)
    Maybe_error<ddMacro_State<vVars...>> update(
        ddMacro_State<vVars...>&& t_prior_all, C_Algo_State&& algo,
        var::Derivative<logL, var::Parameters_transformed> const& t_logL,
        C_elogL const& t_elogL) const {
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
            copy_component(std::type_identity<P_mean>{});
            copy_component(std::type_identity<P_Cov>{});
            copy_component(std::type_identity<Chi2>{});

            evo.emplace_back(std::move(el));
        }
        return std::move(t_prior_all);
    }

    template <class recursive, class averaging, class variance, class variance_correction,
              class FunctionTable, class C_Macro_State, class C_Qdt, class C_Patch_Model,
              class C_double>

        requires(uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
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
                                        variance>(t_prior, t_Qdt, m, Nch, p_y, fs);
        if (!Maybe_Algo)
            return Maybe_Algo.error();

        auto r_Algo_state = std::move(Maybe_Algo.value());
        auto r_y_mean = get<y_mean>(r_Algo_state());
        auto r_y_var = get<y_var>(r_Algo_state());
        auto r_chi2 = get<Chi2>(r_Algo_state());
        auto r_trust_coefficient = get<trust_coefficient>(r_Algo_state());
        auto y = p_y.value();
        bool y_is_nan = std::isnan(y);
        auto r_logL = calculate_logL(y_is_nan, r_y_var, r_chi2, r_trust_coefficient, m);
        auto r_elogL = calculate_elogL(y_is_nan, r_y_var, r_trust_coefficient, m);
        
       
        auto r_prior_all =
            update(std::move(t_prior_all), std::move(r_Algo_state), std::move(r_logL), std::move(r_elogL));
        return std::move(r_prior_all);
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
              class variance_correction, class MacroState, class FuncTable, class C_Parameters,
              class Model>

        requires(uses_adaptive_aproximation_c<adaptive> &&
                 uses_recursive_aproximation_c<recursive> &&
                 uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction> &&
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
        constexpr bool test_derivative = true;
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
            0ul, y().size(), std::move(t_macro),
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
                        auto Maybe_t_Qdtm = [this, &f_local,&m,&t_step,&fs]()
                        {
                            if constexpr( averaging::value>0){
                                return calc_Qdtm(f_local, m, t_step, fs);}
                             else{
                                return calc_Qdtg(f_local, m, t_step, fs);
                             }
                        }();
                        if (!Maybe_t_Qdtm)
                            return Maybe_error<MacroState>(Maybe_t_Qdtm.error());
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

                        if constexpr (test_macroir_derivative &&
                                      var::has_it_v<MacroState, y_mean>) {
                            const auto h = 1e-7;
                            auto f_no_memoi = f_local.to_bare_functions();
                            auto c_prior = t_prior;

                            auto tt_prior =
                                select<logL, Patch_State, y_mean, y_var>(std::move(c_prior));

                            auto test_der_macroir = test_derivative_clarke<false>(
                                [this, &fs, &f_no_memoi, &y, i_step](
                                    auto l_t_prior, auto const& l_Qdtm, auto const& l_m,
                                    auto const& l_Nch)
                                    -> Maybe_error<Transfer_Op_to<
                                        std::decay_t<decltype(l_t_prior)>,
                                        var::Vector_Space<logL, Patch_State, y_mean, y_var>>> {
                                    auto Maybe_res = MacroR2<recursive, averaging, variance,
                                                             variance_correction>{}(
                                        f_no_memoi, std::move(l_t_prior), l_Qdtm, l_m, l_Nch,
                                        y()[i_step], fs);
                                    if (!Maybe_res)
                                        return Maybe_res.error();
                                    return select<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var>(
                                        std::move(Maybe_res.value()));
                                },
                                h, tt_prior, t_Qdtm, m, Nch);
                            if (true && !test_der_macroir) {
                                std::cerr << "\nError on i_step: " << i_step
                                          << "\ty: " << y()[i_step] << " t_step: " << t_step
                                          << "\n";
                                std::cerr << test_der_macroir.error()();
                            }
                        }

                        return MacroR2<recursive, averaging, variance, variance_correction>{}(
                            f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);

                    } else {
                        auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                        if (!Maybe_t_Qdt)
                            return Maybe_error<MacroState>(Maybe_t_Qdt.error());
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
                                        [this, &t_step, &fs, &f_no_memoi](auto const& l_m,
                                                                          auto const& /*l_Qx*/) {
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

                        return MacroR2<recursive, averaging, variance, variance_correction>{}(
                            f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                    }
                } else {
                    if constexpr (!variance_correction::value) {
                        auto Maybe_t_Qdtm = [this,& f_local, &m, &t_step, &fs]() {
                           if constexpr(averaging::value>0){
                            return calc_Qdtm(f_local, m, t_step, fs);}
                           else{
                            return calc_Qdtg(f_local, m, t_step, fs);}
                        }();
                        if (!Maybe_t_Qdtm)
                            return Maybe_error<MacroState>(Maybe_t_Qdtm.error());
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
                            return MacroR2<recursive, averaging, variance, variance_correction>{}(
                                f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);

                        } else {
                            return
                                //   f_local.f(
                                // MacroR2<::V<uses_recursive_aproximation<false>>,
                                //         v_averaging,
                                //         ::V<uses_variance_aproximation<false>>,
                                //         ::V<uses_taylor_variance_correction_aproximation(
                                //             false)>>{},
                                // std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);

                                MacroR2<uses_recursive_aproximation<false>, averaging, variance,
                                        uses_taylor_variance_correction_aproximation<false>>{}(
                                    f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step], fs);
                        }
                    } else {
                        auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                        if (!Maybe_t_Qdt)
                            return Maybe_error<MacroState>(Maybe_t_Qdt.error());
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
                            return MacroR2<recursive, averaging, variance, variance_correction>{}(
                                f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                        } else {
                            return MacroR2<uses_recursive_aproximation<false>, averaging, variance,
                                           uses_taylor_variance_correction_aproximation<false>>{}(
                                f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                        }
                    }
                }
            });
        f += f_local;
        if (!Maybe_run)
            return Maybe_run.error();
        return std::move(Maybe_run.value());
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
        auto N = get<N_Ch_mean>(m)()[0];
        auto sim = Simulated_Recording<Simulate_tag>{};
        auto N_state = sample_Multinomial(mt, r_P_mean, N);
        return Simulated_Step<Simulate_tag>(
            Vector_Space(std::move(N_state), Simulated_Recording<Simulate_tag>{}));
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
            auto fs = get<Frequency_of_Sampling>(e).value();
            auto sim_recording = Recording{};

            auto ini = init_sim<Simulate_tag>(mt, m, e);
            auto run = fold(get<Recording_conditions>(e)(), ini,
                            [this, &m, fs, n_sub_dt, &mt](Simulated_Step<Simulate_tag>&& t_sim_step,
                                                          Experiment_step const& t_step) {
                                return Maybe_error<Simulated_Step<Simulate_tag>>(sub_sample(
                                    mt, std::move(t_sim_step), m, t_step, n_sub_dt(), fs));
                            });
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

template <class recursive, class averaging, class variance, class variance_correction>

    requires(uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
             uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
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

        return out;
    }

    template <class T, class... Ts>
    auto operator()(T&& x, Ts&&... xs) {
        auto m = Macro_DMR{};
        auto l1 = m.Macror<recursive, averaging, variance, variance_correction>(
            std::forward<T>(x), std::forward<Ts>(xs)...);

        return std::move(l1);
    }
};

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class Model>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
struct Likelihood_Model_constexpr {
    Model m;
    Simulation_n_sub_dt n_sub_dt;
    Likelihood_Model_constexpr(const Model& model, Simulation_n_sub_dt n_sub_dt)
        : m{model}, n_sub_dt{n_sub_dt} {}
    Likelihood_Model_constexpr(adaptive, recursive, averaging, variance, variance_correction,
                               const Model& model, Simulation_n_sub_dt n_sub_dt)
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
        f << "Simulation_n_sub_dt: " << d.n_sub_dt << "\n";
        report_model(s, d.m);
    }
};

template <class adaptive_range, class recursive_range, class averaging_range, class variance_range,
          class taylor_variance_correction_range, class Model>
    requires(
        var::is_this_constexpr_Var_domain_c<adaptive_range, bool, uses_adaptive_aproximation> &&
        var::is_this_constexpr_Var_domain_c<recursive_range, bool, uses_recursive_aproximation> &&
        var::is_this_constexpr_Var_domain_c<averaging_range, int, uses_averaging_aproximation> &&
        var::is_this_constexpr_Var_domain_c<variance_range, bool, uses_variance_aproximation> &&
        var::is_this_constexpr_Var_domain_c<taylor_variance_correction_range, bool,
                                            uses_taylor_variance_correction_aproximation>)
struct Likelihood_Model_regular {
    Model m;
    Simulation_n_sub_dt n_sub_dt;
    uses_adaptive_aproximation_value adaptive;
    uses_recursive_aproximation_value recursive;
    uses_averaging_aproximation_value averaging;
    uses_variance_aproximation_value variance;
    uses_taylor_variance_correction_aproximation_value taylor_variance_correction;

    static constexpr adaptive_range range_adaptive = {};
    static constexpr recursive_range range_recursive = {};
    static constexpr averaging_range range_averaging = {};
    static constexpr variance_range range_variance = {};
    static constexpr taylor_variance_correction_range range_variance_correction = {};

    Likelihood_Model_regular(
        const Model& model, Simulation_n_sub_dt n_sub_dt, uses_adaptive_aproximation_value adaptive,
        uses_recursive_aproximation_value recursive, uses_averaging_aproximation_value averaging,
        uses_variance_aproximation_value variance,
        uses_taylor_variance_correction_aproximation_value taylor_variance_correction)
        : m{model},
          n_sub_dt{n_sub_dt},
          adaptive{adaptive},
          recursive{recursive},
          averaging{averaging},
          variance{variance},
          taylor_variance_correction{taylor_variance_correction} {}

    using cartesian = algebra_2<std::tuple, std::variant>::P_constexpr<
        typename adaptive_range::variant_type, typename recursive_range::variant_type,
        typename averaging_range::variant_type, typename variance_range::variant_type,
        typename taylor_variance_correction_range::variant_type>;

    template <class, class>
    struct Likelihood_Model_variant_impl;
    template <class... adaptive, class... recursive, class... averaging, class... variance,
              class... taylor_variance_correction, class M>

        requires(
            (uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<taylor_variance_correction>) &&
            ...)
    struct Likelihood_Model_variant_impl<
        std::variant<
            std::tuple<adaptive, recursive, averaging, variance, taylor_variance_correction>...>,
        M> {
        using type =
            std::variant<Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                                    taylor_variance_correction, M>...>;
    };

    using Likelihood_Model_variant = typename Likelihood_Model_variant_impl<cartesian, Model>::type;

    Maybe_error<Likelihood_Model_variant> get_variant() const {
        auto car = promote_Maybe_error(std::make_tuple(
            adaptive_range::to_variant(adaptive.value),
            recursive_range::to_variant(recursive.value),
            averaging_range::to_variant(averaging.value),
            variance_range::to_variant(variance.value),
            taylor_variance_correction_range::to_variant(taylor_variance_correction.value)));

        if (!car) {
            return car.error();
        }
        return std::apply(
            [this](auto... t) {
                return std::visit(
                    [this](auto... v) -> Likelihood_Model_variant {
                        auto l = Likelihood_Model_constexpr<std::decay_t<decltype(v)>..., Model>(
                            m, n_sub_dt);
                        return Likelihood_Model_variant{std::move(l)};
                    },
                    t...);
            },
            car.value());
    }
};

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
          class variance_correction, class Model>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
auto make_Likelihood_Model(const Model& m, Simulation_n_sub_dt n_sub_dt) {
    return Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                      Model>(m, n_sub_dt);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FuncTable, class Model, class Parameters, class Recoding,
          class Experiment>

    requires(uses_adaptive_aproximation_c<adaptive>,
             uses_recursive_aproximation_c<recursive> && uses_averaging_aproximation_c<averaging> &&
                 uses_variance_aproximation_c<variance> &&
                 uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<logLs> logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    Parameters const& p, const Recoding& y, const Experiment& e) {
    auto v_logL = Macro_DMR{}
                      .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                                      Macro_State_reg>(f, lik.m, p, y, e);
    if (!v_logL)
        return v_logL.error();
    else
        return logLs(get<logL>(v_logL.value()), get<elogL>(v_logL.value()),
                     get<vlogL>(v_logL.value()));
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FuncTable, class Model, class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<dMacro_State_Hessian_minimal> dlogLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                        dMacro_State_Hessian_minimal>(f, lik.m, dpp, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FuncTable, class Model, class Variables, class DataType>
    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<diff_Macro_State_Gradient_Hessian> diff_logLikelihood(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var,
    double delta_par) {
    auto v_p = p.to_value();
    auto Maybe_MacroEv =
        Macro_DMR{}
            .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
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
            t_yvar_2[i] = 2.0 * t_yvar[i] * t_yvar[i];
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
        std::move(get<elogL>(v_MacroEv)), std::move(get<vlogL>(v_MacroEv)), Grad(std::move(G)),
        FIM(std::move(r_FIM)));
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<Macro_State_Ev_predictions> logLikelihoodPredictions(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                        Macro_State_Ev_predictions>(f, lik.m, p, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<Macro_State_Ev_diagnostic> logLikelihoodDiagnostic(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                        Macro_State_Ev_diagnostic>(f, lik.m, p, y, var);
}







template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FuncTable, class Model, class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<dMacro_State_Ev_gradient_all> dlogLikelihoodPredictions(
    FuncTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    var::Parameters_transformed const& p, const DataType& y, const Variables& var) {
    auto dp = var::selfDerivative(p);
    auto dpp = dp.to_value();
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                        dMacro_State_Ev_gradient_all>(f, lik.m, dpp, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<Macro_State_Ev_diagnostic> logLikelihood_Diagnostic(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
    Parameters const& p, const DataType& y, const Variables& var) {
    return Macro_DMR{}
        .log_Likelihood<adaptive, recursive, averaging, variance, variance_correction,
                        Macro_State_Ev_diagnostic>(f, lik.m, p, y, var);
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class FunctionTable, class Model, class Parameters,
          class Variables, class DataType>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
Maybe_error<std::vector<Macro_State_Ev_predictions>> fractioned_logLikelihoodPredictions(
    FunctionTable& f,
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, variance_correction,
                                     Model>& lik,
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
          class variance_correction, class Model, class Parameters, class Variables>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction> &&
             !is_of_this_template_type_v<Variables, std::vector>)
auto simulate(mt_64i& mt,
              const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                               variance_correction, Model>& lik,
              Parameters const& p, const Variables& var) {
    return Macro_DMR{}.sample(mt, lik.m, p, var, lik.n_sub_dt).value()();
}

template <class adaptive, class recursive, class averaging, class variance,
          class variance_correction, class Model, class Parameters, class Variables>

    requires(uses_adaptive_aproximation_c<adaptive> && uses_recursive_aproximation_c<recursive> &&
             uses_averaging_aproximation_c<averaging> && uses_variance_aproximation_c<variance> &&
             uses_taylor_variance_correction_aproximation_c<variance_correction>)
auto simulate(mt_64i& mt,
              const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance,
                                               variance_correction, Model>& lik,
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
#ifdef ZOMBIE
namespace zombie {

class experiment_fractioner {
    std::vector<std::size_t> segments = {73, 33, 22, 22, 4};
    bool average_agonist_evolution = 10;

   public:
    experiment_fractioner(const std::vector<std::size_t>& t_segments,
                          bool average_the_agonist_evolution)
        : segments{t_segments}, average_agonist_evolution{average_the_agonist_evolution} {}

    static auto average_Recording(const Recording_conditions& e, const Recording& y,
                                  const std::vector<std::size_t>& indexes0,
                                  const std::vector<std::size_t>& indexes1,
                                  bool average_the_evolution) {
        assert(size(y()) == size(indexes0));
        assert(size(e()) == size(y()));
        auto out_size = indexes1.size();
        std::vector<Patch_current> out_y(out_size);
        std::vector<Experiment_step> out_x(out_size);
        double sum_y = 0;
        std::size_t sum_samples = 0;
        std::size_t ii = 0;

        Agonist_evolution v_agonist = std::vector<Agonist_step>{};

        for (std::size_t i = 0; i < size(y()); ++i) {
            if (indexes0[i] == indexes1[ii]) {
                if (sum_samples == 0) {
                    out_y[ii] = y()[i];
                    out_x[ii] = average_Experimental_step(e()[i], average_the_evolution);
                } else {
                    auto n_samples = get_num_samples(e()[i]);
                    sum_y += y()[i]() * n_samples;
                    sum_samples += n_samples;
                    v_agonist = add_agonist_step(
                        std::move(v_agonist), average_agonist_step(e()[i], average_the_evolution));

                    out_y[ii] = Patch_current(sum_y / sum_samples);
                    out_x[ii] = Experiment_step(get<Time>(e()[i]), v_agonist);
                    sum_y = 0;
                    sum_samples = 0;
                    v_agonist = std::vector<Agonist_step>{};
                }
                ++ii;
            } else {
                assert(indexes0[i] < indexes1[ii] && "indexes fails");
                auto n_samples = get_num_samples(e()[i]);
                sum_y += y()[i]() * n_samples;
                sum_samples += n_samples;
                v_agonist = add_agonist_step(std::move(v_agonist),
                                             average_agonist_step(e()[i], average_the_evolution));
            }
        }

        return std::tuple(Recording_conditions(out_x), Recording(out_y));
    }

    template <typename Simulate_tag>

    static auto average_Recording(const Recording_conditions& e,
                                  const Simulated_Recording<Simulate_tag>& sim,
                                  const std::vector<std::size_t>& indexes0,
                                  const std::vector<std::size_t>& indexes1,
                                  bool average_the_evolution) {
        auto& yr = get<Recording>(sim());
        assert(size(yr()) == size(indexes0));
        assert(size(e()) == size(yr()));
        auto out_size = indexes1.size();
        std::vector<Patch_current> out_y(out_size);
        std::vector<N_channel_state> out_N(out_size);

        std::vector<Experiment_step> out_x(out_size);
        double sum_y = 0;
        std::size_t sum_samples = 0;
        std::size_t ii = 0;

        Agonist_evolution v_agonist = std::vector<Agonist_step>{};

        for (std::size_t i = 0; i < size(yr()); ++i) {
            if (indexes0[i] == indexes1[ii]) {
                if (sum_samples == 0) {
                    out_y[ii] = yr()[i];
                    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
                        out_N[ii] = get<N_Ch_State_Evolution>(sim())()[i];
                    out_x[ii] = average_Experimental_step(e()[i], average_the_evolution);
                } else {
                    auto n_samples = get_num_samples(e()[i]);
                    sum_y += yr()[i]() * n_samples;
                    sum_samples += n_samples;
                    v_agonist = add_agonist_step(
                        std::move(v_agonist), average_agonist_step(e()[i], average_the_evolution));

                    out_y[ii] = Patch_current(sum_y / sum_samples);

                    out_x[ii] = Experiment_step(get<Time>(e()[i]), v_agonist);
                    if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
                        out_N[ii] = get<N_Ch_State_Evolution>(sim())()[i];

                    sum_y = 0;
                    sum_samples = 0;
                    v_agonist = std::vector<Agonist_step>{};
                }
                ++ii;
            } else {
                assert(indexes0[i] < indexes1[ii] && "indexes fails");
                auto n_samples = get_num_samples(e()[i]);
                sum_y += yr()[i]() * n_samples;
                sum_samples += n_samples;
                v_agonist = add_agonist_step(std::move(v_agonist),
                                             average_agonist_step(e()[i], average_the_evolution));
            }
        }
        Simulated_Recording<Simulate_tag> out_sim;
        get<Recording>(out_sim())() = std::move(out_y);
        if constexpr (var::has_it_v<Simulate_tag, N_Ch_State_Evolution>)
            get<N_Ch_State_Evolution>(out_sim())() = std::move(out_N);

        return std::tuple(Recording_conditions(out_x), std::move(out_sim));
    }

    auto operator()(const Recording& y, const Experiment& x, mt_64i& mt, std::size_t num_parameters,
                    double n_points_per_decade_beta, double n_points_per_decade_fraction,
                    double stops_at, bool includes_zero) const {
        assert(size(y()) == size(get<Recording_conditions>(x)()));
        assert(size(y()) == var::sum(segments));

        std::size_t num_samples = size(y());
        std::size_t max_num_samples_per_segment = var::max(segments);

        auto cum_segments = var::cumsum(segments);

        auto indexes =
            generate_random_Indexes(mt, num_samples, 1, n_points_per_decade_fraction, cum_segments);
        // std::cerr <<
        // "\nindexes\n**************************************************"
        //              "*************************\n";
        // std::cerr << indexes;
        //  std::abort();
        auto n_frac = size(indexes);
        cuevi::by_fraction<Recording> y_out(n_frac);
        cuevi::by_fraction<Experiment> x_out(
            n_frac, Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                               get<initial_agonist_concentration>(x)));
        y_out[n_frac - 1] = y;
        x_out[n_frac - 1] = x;

        for (std::size_t i = n_frac - 1; i > 0; --i) {
            std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
                average_Recording(get<Recording_conditions>(x_out[i]), y_out[i], indexes[i],
                                  indexes[i - 1], average_agonist_evolution);
        }

        // std::abort();
        auto beta0 = get_beta_list(n_points_per_decade_beta,
                                   stops_at * num_samples / (n_frac > 1 ? size(indexes[0]) : 1),
                                   includes_zero);
        by_beta<double> betan = {0, 1};
        cuevi::by_fraction<by_beta<double>> beta(n_frac, betan);
        beta[0] = std::move(beta0);

        return std::tuple(std::move(y_out), std::move(x_out), std::move(beta));
    }

    auto operator()(const Recording& y, const Experiment& x, mt_64i& mt, std::size_t num_parameters,
                    double n_points_per_decade_fraction) const {
        assert(size(y()) == size(get<Recording_conditions>(x)()));
        assert(size(y()) == var::sum(segments));

        std::size_t num_samples = size(y());
        //  std::size_t max_num_samples_per_segment = var::max(segments);

        auto cum_segments = var::cumsum(segments);

        auto indexes = generate_random_Indexes(mt, num_samples, num_parameters,
                                               n_points_per_decade_fraction, cum_segments);
        // std::cerr <<
        // "\nindexes\n**************************************************"
        //              "*************************\n";
        // std::cerr << indexes;
        // std::abort();
        auto n_frac = size(indexes);
        cuevi::by_fraction<Recording> y_out(n_frac);
        cuevi::by_fraction<Experiment> x_out(
            n_frac, Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                               get<initial_agonist_concentration>(x)));
        y_out[n_frac - 1] = y;
        x_out[n_frac - 1] = x;

        for (std::size_t i = n_frac - 1; i > 0; --i) {
            std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
                average_Recording(get<Recording_conditions>(x_out[i]), y_out[i], indexes[i],
                                  indexes[i - 1], average_agonist_evolution);
        }
        return std::tuple(std::move(y_out), std::move(x_out));
    }

    template <typename Simulate_tag>

    auto operator()(const Simulated_Recording<Simulate_tag>& sim, const Experiment& x, mt_64i& mt,
                    std::size_t num_parameters, double n_points_per_decade_fraction) const {
        auto& yr = get<Recording>(sim());
        assert(size(yr()) == size(get<Recording_conditions>(x)()));
        assert(size(segments) == 0 || size(yr()) == var::sum(segments));

        std::size_t num_samples = size(yr());
        //  std::size_t max_num_samples_per_segment = var::max(segments);

        auto cum_segments = var::cumsum(segments);

        auto indexes = generate_random_Indexes(mt, num_samples, num_parameters,
                                               n_points_per_decade_fraction, cum_segments);
        // std::cerr <<
        // "\nindexes\n**************************************************"
        //              "*************************\n";
        // std::cerr << indexes;
        // std::abort();
        auto n_frac = size(indexes);
        cuevi::by_fraction<Simulated_Recording<Simulate_tag>> y_out(n_frac);
        cuevi::by_fraction<Experiment> x_out(
            n_frac, Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                               get<initial_agonist_concentration>(x)));
        y_out[n_frac - 1] = sim;
        x_out[n_frac - 1] = x;

        for (std::size_t i = n_frac - 1; i > 0; --i) {
            std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
                average_Recording(get<Recording_conditions>(x_out[i]), y_out[i], indexes[i],
                                  indexes[i - 1], average_agonist_evolution);
        }
        return std::tuple(std::move(y_out), std::move(x_out));
    }
};
}  // namespace zombie
#endif

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
