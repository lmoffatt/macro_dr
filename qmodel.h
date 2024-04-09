#pragma once
#include "cuevi.h"
// #include "experiment.h"
#include "experiment.h"
#include "fold.h"
#include "function_memoization.h"
#include "matrix.h"
// #include "models_MoffattHume_linear.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>
#ifndef QMODEL_H
#define QMODEL_H
#include <map>
#include <string>

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

template <class T> T sqr(T x) { return x * x; }

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

inline std::size_t
get_max_state(std::vector<std::pair<std::pair<std::size_t, std::size_t>,
                                    std::string>> const &new_formulas) {
  std::size_t out = 0;
  for (auto &e : new_formulas) {
    if (e.first.first > out)
      out = e.first.first;
    if (e.first.second > out)
      out = e.first.second;
  }
  return out;
}
class Q0_formula
    : public var::Var<Q0_formula, std::vector<std::vector<std::string>>> {
public:
  using var::Var<Q0_formula, std::vector<std::vector<std::string>>>::Var;

  Q0_formula(std::size_t N)
      : var::Var<Q0_formula, std::vector<std::vector<std::string>>>{
            std::vector<std::vector<std::string>>{
                N, std::vector<std::string>{N, ""}}} {}

  friend std::ostream &operator<<(std::ostream &os, Q0_formula const &x) {
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
class Qa_formula
    : public var::Var<Qa_formula, std::vector<std::vector<std::string>>> {
public:
  using base_type = var::Var<Qa_formula, std::vector<std::vector<std::string>>>;
  using base_type::Var;
  Qa_formula(std::size_t N)
      : var::Var<Qa_formula, std::vector<std::vector<std::string>>>{
            std::vector<std::vector<std::string>>{
                N, std::vector<std::string>{N, ""}}} {}
  friend std::ostream &operator<<(std::ostream &os, Qa_formula const &x) {
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
  requires(std::is_same_v<Q_formula, Q0_formula> ||
           std::is_same_v<Q_formula, Qa_formula>)
auto change_states_number(const Q_formula &f, std::size_t N) {
  Q_formula out(N);
  for (std::size_t i = 0; i < f().size(); ++i)
    for (std::size_t j = 0; j < f().size(); ++j)
      out()[i][j] = f()[i][j];

  return out;
}

template <class Q_formula>
  requires(std::is_same_v<Q_formula, Q0_formula> ||
           std::is_same_v<Q_formula, Qa_formula>)
auto insert_new_formula(const Q_formula &f, std::size_t i_ini,
                        std::size_t i_end, std::string &&formula) {
  auto N = std::max(i_ini, i_end) + 1;
  auto out = change_states_number(f, N);
  out()[i_ini][i_end] = std::move(formula);

  return out;
}

class Qx : public var::Var<Qx, Matrix<double>> {};
class P_initial : public var::Var<P_initial, Matrix<double>> {};

class g : public var::Var<g, Matrix<double>> {};
class g_formula : public var::Var<g_formula, std::vector<std::string>> {};

inline auto change_states_number(const g_formula &f, std::size_t N) {
  g_formula out(std::vector<std::string>{N, ""});
  for (std::size_t i = 0; i < f().size(); ++i)
    out()[i] = f()[i];
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

class Binomial_magical_number
    : public var::Constant<Binomial_magical_number, double> {};

class min_P : public var::Constant<min_P, double> {};

class N_Ch_std : public var::Var<N_Ch_std, double> {};

class Current_Noise : public var::Var<Current_Noise, double> {};

class Pink_Noise : public var::Var<Pink_Noise, double> {};

class Proportional_Noise : public var::Var<Proportional_Noise, double> {};

class Current_Baseline : public var::Var<Current_Baseline, double> {};

class P_mean : public var::Var<P_mean, Matrix<double>> {};

template <class C_Matrix> auto to_Probability(C_Matrix const &x) {

  using std::max;
  using std::min;
  auto out = apply([](auto e) { return min(max(e, 0.0), 1.0); }, x);
  auto s = var::sum(out);
  return out * (1.0 / s);
}

inline bool all_Probability_elements(Matrix<double> const &x) {

  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] > 1.0)
      return false;
    else if (x[i] < 0.0)
      return false;
  }
  return true;
}

template <class C_Matrix>
inline bool all_Covariance_elements(C_Matrix const &x) {
  if (x.ncols() != x.nrows())
    return false;
  for (std::size_t i = 0; i < x.nrows(); ++i) {
    if (x(i, i) < 0)
      return false;
    for (std::size_t j = 0; j < i; ++j)
      if ((x(i, i) > 0) && (x(j, j) > 0)) {
        if (x(i, j) * x(i, j) / x(i, i) / x(j, j) > 1.0)
          return false;
      } else {
        if (x(i, j) != 0)
          return false;
      }
  }
  return true;
}

inline bool crude_lambda_violations(DiagonalMatrix<double> const &l) {

  if (var::max(l) > 1e-2)
    return true;
  else
    return false;
}

template <class C_Matrix>
auto to_Transition_Probability_Eigenvalues(C_Matrix &&lambda) {

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

class N_channel_state : public var::Var<N_channel_state, Matrix<double>> {
  using var::Var<N_channel_state, Matrix<double>>::Var;
};

class y_sum : public var::Var<y_sum, double> {};

class t_sum : public var::Var<t_sum, double> {};

class P_Cov : public var::Var<P_Cov, SymmetricMatrix<double>> {};

class lambda : public var::Var<lambda, DiagonalMatrix<double>> {};

class V : public var::Var<V, Matrix<double>> {};
class W : public var::Var<W, Matrix<double>> {};

class uses_variance_aproximation
    : public var::struct_Var<uses_variance_aproximation, bool> {};

class uses_variance_correction_aproximation
    : public var::struct_Var<uses_variance_correction_aproximation, bool> {};

class uses_adaptive_aproximation
    : public var::struct_Var<uses_adaptive_aproximation, bool> {};
class uses_recursive_aproximation
    : public var::struct_Var<uses_recursive_aproximation, bool> {};
class uses_averaging_aproximation
    : public var::struct_Var<uses_averaging_aproximation, int> {};

class return_predictions : public var::struct_Var<return_predictions, bool> {};

class Probability_error_tolerance
    : public var::Constant<Probability_error_tolerance, double> {};

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
  friend std::string className(P) { return "P_ij"; }
};

template <class C_Matrix>
Maybe_error<Transfer_Op_to<C_Matrix, P>>
to_Transition_Probability(C_Matrix const &x) {
  using std::max;
  using std::min;
  auto out = apply([](auto e) { return min(max(e, 0.0), 1.0); }, x);
  auto sumP = out * Matrix<double>(out.ncols(), 1ul, 1.0);
  auto s = inv(diag(sumP));
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
inline Maybe_error<Ptotal_ij> make_Ptotal_ij(Matrix<double> &&x,
                                             double max_dt) {
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
class plogL : public var::Var<plogL, double> {
  friend std::string className(plogL) { return "plogL"; }
};
class eplogL : public var::Var<eplogL, double> {
  friend std::string className(eplogL) { return "eplogL"; }
};
class vplogL : public var::Var<vplogL, double> {
  friend std::string className(vplogL) { return "vplogL"; }
};

class macror_algorithm : public var::Constant<macror_algorithm, std::string> {
  using var::Constant<macror_algorithm, std::string>::Constant;
  friend std::string className(macror_algorithm) { return "macror_algorithm"; }
};

class PGn : public var::Var<PGn, Matrix<double>> {};
class PGG_n : public var::Var<PGG_n, Matrix<double>> {};
class PG_n : public var::Var<PG_n, Matrix<double>> {};
// class PPn : public var::Var<PPn, Matrix<double>> {};

using Qn = Vector_Space<number_of_samples, min_P, P, PG_n, PGG_n>;

using Qx_eig = Vector_Space<Qx, V, lambda, W>;

using Qdtm = Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij,
                          gmean_ij, gsqr_i, gvar_i>;

using Qdt =
    Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij,
                 gtotal_sqr_ij, gsqr_i, gvar_i, gtotal_var_ij, gvar_ij>;

using Patch_Model =
    Vector_Space<N_St, Q0, Qa, P_initial, g, N_Ch_mean, Current_Noise,
                 Pink_Noise, Proportional_Noise, Current_Baseline,
                 N_Ch_mean_time_segment_duration, Binomial_magical_number,
                 min_P, Probability_error_tolerance,
                 Conductance_variance_error_tolerance>;

inline void save(const std::string name, const Patch_Model &m) {
  std::ofstream f_Q0(name + "_Q0.txt");
  f_Q0 << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << get<Q0>(m) << "\n";
  std::ofstream f_Qa(name + "_Qa.txt");
  f_Qa << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << get<Qa>(m) << "\n";
  std::ofstream f_g(name + "_g.txt");
  f_g << std::setprecision(std::numeric_limits<double>::digits10 + 1)
      << get<g>(m) << "\n";
}

using Patch_State =
    Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL,
                 eplogL, vplogL, macror_algorithm>;

template <class C_Patch_Model, class C_double>
C_Patch_Model add_Patch_inactivation(C_Patch_Model &&m,
                                     C_double const &deactivation_rate) {
  using Transf = transformation_type_t<C_Patch_Model>;
  auto Nst = get<N_St>(m)() + 1;
  Op_t<Transf, Q0> v_Q0 = Q0(Matrix<double>(Nst, Nst, 0.0));
  for (std::size_t i = 0; i + 1 < Nst; ++i) {
    for (std::size_t j = 0; j + 1 < Nst; ++j)
      set(v_Q0(), i, j, get<Q0>(m)()(i, j));
    set(v_Q0(), i, Nst - 1, deactivation_rate);
  }
  Op_t<Transf, Qa> v_Qa = Qa(Matrix<double>(Nst, Nst, 0.0));
  for (std::size_t i = 0; i + 1 < Nst; ++i) {
    for (std::size_t j = 0; j + 1 < Nst; ++j)
      set(v_Qa(), i, j, get<Qa>(m)()(i, j));
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
  return std::move(m);
}

template <class C_Patch_Model>
std::pair<C_Patch_Model, double> remove_Patch_inactivation(C_Patch_Model const &mi) {
  auto m=mi;
  
  
  using Transf = transformation_type_t<C_Patch_Model>;
  auto Nst = get<N_St>(m)() - 1;
  Op_t<Transf, Q0> v_Q0 = Q0(Matrix<double>(Nst, Nst, 0.0));
  auto deactivation_rate=get<Q0>(m)()(0ul,Nst);
  for (std::size_t i = 0; i  < Nst; ++i) {
    for (std::size_t j = 0; j  < Nst; ++j)
      set(v_Q0(), i, j, get<Q0>(m)()(i, j));
    assert(deactivation_rate==get<Q0>(m)()(i,Nst));
  }
  
  Op_t<Transf, Qa> v_Qa = Qa(Matrix<double>(Nst, Nst, 0.0));
  for (std::size_t i = 0; i  < Nst; ++i) {
    for (std::size_t j = 0; j < Nst; ++j)
      set(v_Qa(), i, j, get<Qa>(m)()(i, j));
  }
  Op_t<Transf, g> v_g = g(Matrix<double>(Nst, 1, 0.0));
  for (std::size_t i = 0; i  < Nst; ++i) {
    set(v_g(), i, 0, get<g>(m)()[i]);
  }
  Op_t<Transf, P_initial> v_Pini = P_initial(Matrix<double>(1, Nst, 0.0));
  for (std::size_t i = 0; i  < Nst; ++i) {
    set(v_Pini(), 0, i, get<P_initial>(m)()[i]);
  }
  v_Pini()=v_Pini()/var::sum(v_Pini());
  
  get<N_St>(m)() = Nst;
  get<Qa>(m) = v_Qa;
  get<Q0>(m) = v_Q0;
  get<g>(m) = v_g;
  get<P_initial>(m) = v_Pini;
  return std::pair(std::move(m),deactivation_rate);
}



class Patch_State_Evolution
    : public Var<Patch_State_Evolution, std::vector<Patch_State>> {};

using Patch_State_and_Evolution =
    Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL,
                 eplogL, vplogL, macror_algorithm, Patch_State_Evolution>;

class Simulation_n_sub_dt : public Var<Simulation_n_sub_dt, std::size_t> {};

class includes_N_state_evolution
    : public var::struct_Var<includes_N_state_evolution, bool> {};
class N_Ch_State_Evolution
    : public Var<N_Ch_State_Evolution, std::vector<N_channel_state>> {};

template <includes_N_state_evolution keep_N_state>
class Simulated_Recording
    : public Var<
          Simulated_Recording<keep_N_state>,
          std::conditional_t<keep_N_state.value,
                             Vector_Space<Recording, N_Ch_State_Evolution>,
                             Vector_Space<Recording>>> {
public:
  constexpr static const bool includes_N = keep_N_state.value;
  using Var<Simulated_Recording<keep_N_state>,
            std::conditional_t<keep_N_state.value,
                               Vector_Space<Recording, N_Ch_State_Evolution>,
                               Vector_Space<Recording>>>::Var;
};

using v_Simulated_Recording =
    std::variant<Simulated_Recording<includes_N_state_evolution(false)>,
                 Simulated_Recording<includes_N_state_evolution(true)>>;

template <includes_N_state_evolution keep_N_state>
class Simulated_Step
    : public Var<
          Simulated_Step<keep_N_state>,
          Vector_Space<N_channel_state, Simulated_Recording<keep_N_state>>> {
  // using Vector_Space<N_channel_state,
  // Simulated_Recording<keep_N_state>>::Vector_Space;
};

using Simulated_Sub_Step =
    Vector_Space<N_channel_state, number_of_samples, y_sum>;

using Simulation_Parameters = Vector_Space<Simulation_n_sub_dt>;

template <includes_N_state_evolution includesN>
void save_simulation(std::string const &filename, const std::string &separator,
                     Simulated_Recording<includesN> const &sim) {
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_step" << separator << "patch_current";
  if constexpr (includesN.value)
    f << separator << "i_state" << separator << "N_state";
  f << "\n";

  if constexpr (includesN.value) {
    auto N = get<N_Ch_State_Evolution>(sim());
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
      for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
        f << i_step << separator << y()[i_step]() << separator << i_state
          << separator << N()[i_step]()[i_state] << "\n";
      }
    }
  } else {
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
      f << i_step << separator << y()[i_step]() << "\n";
    }
  }
}

template <includes_N_state_evolution keep_N_state>
Maybe_error<bool> load_simulation(std::string const &fname,
                                  std::string separator,
                                  Simulated_Recording<keep_N_state> &r) {
  std::ifstream f(fname);
  if (!f)
    return error_message("cannot open file " + fname);
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);

  if constexpr (keep_N_state.value) {
    if (!(ss >> septr("i_step") >> septr(separator) >> septr("patch_current") >>
          septr(separator) >> septr("i_state") >> septr(separator) >>
          septr("N_state")))
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
  auto &e = get<Recording>(r());

  if constexpr (!keep_N_state.value) {
    while (
        extract_double(ss >> i_step >> septr(separator), val, separator[0])) {
      if (i_step_prev != i_step) {
        if (i_step != e().size())
          return error_message("i_step missmatch expected" +
                               std::to_string(e().size()) +
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
    while (
        (extract_double(ss >> i_step >> septr(separator), val, separator[0]) >>
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
          return error_message("i_step missmatch expected" +
                               std::to_string(e().size()) +
                               " found:" + std::to_string(i_step));
        e().push_back(Patch_current(val));
        if (n_channel_states == 0)
          n_channel_states = N_state.size();
        else if (n_channel_states != N_state.size())
          return error_message("n_channel_states missmatch expected" +
                               std::to_string(n_channel_states) +
                               " found:" + std::to_string(N_state.size()));
        auto &Ns = get<N_Ch_State_Evolution>(r());
        Ns().emplace_back(
            N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
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
    auto &Ns = get<N_Ch_State_Evolution>(r());
    Ns().emplace_back(
        N_channel_state(Matrix<double>(1, n_channel_states, N_state)));

    return true;
  }
}

template <includes_N_state_evolution includesN>
void save_fractioned_simulation(
    std::string filename, std::string separator,
    std::vector<Simulated_Recording<includesN>> const &vsim) {
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_frac" << separator << "i_step" << separator << "patch_current";
  if constexpr (includesN.value)
    f << separator << "i_state" << separator << "N_state";
  f << "\n";

  if constexpr (includesN.value) {
    for (std::size_t i_frac = 0; i_frac < vsim.size(); ++i_frac) {
      auto &sim = vsim[i_frac];
      auto N = get<N_Ch_State_Evolution>(sim());
      auto y = get<Recording>(sim());
      for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
        for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
          f << i_frac << separator << i_step << separator << y()[i_step]()
            << separator << i_state << separator << N()[i_step]()[i_state]
            << "\n";
        }
      }
    }
  } else {
    for (std::size_t i_frac = 0; i_frac < vsim.size(); ++i_frac) {
      auto &sim = vsim[i_frac];
      auto y = get<Recording>(sim());
      for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
        f << i_frac << separator << i_step << separator << y()[i_step]()
          << "\n";
      }
    }
  }
}

template <includes_N_state_evolution keep_N_state>
Maybe_error<bool>
load_fractioned_simulation(std::string const &fname, std::string separator,
                           std::vector<Simulated_Recording<keep_N_state>> &rs) {
  std::ifstream f(fname);
  if (!f)
    return error_message("cannot open file " + fname);
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);

  if (!(ss >> septr("i_frac") >> septr(separator) >> septr("i_step") >>
        septr(separator) >> septr("patch_current")))
    return error_message("titles are wrong : expected  i_step:" + separator +
                         "patch_current; found:" + line);
  if constexpr (keep_N_state.value)
    if (!(ss >> septr(separator) >> septr("i_state") >> septr(separator) >>
          septr("N_state")))
      return error_message("titles are wrong : expected  i_state:" + separator +
                           "N_state; found:" + line);

  std::getline(f, line);
  ss = std::stringstream(line);
  std::size_t i_frac;
  std::size_t i_frac_prev = std::numeric_limits<std::size_t>::max();
  std::size_t i_step;
  std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();

  double val;
  Simulated_Recording<keep_N_state> r;

  if constexpr (!keep_N_state.value) {
    while (extract_double(ss >> i_frac >> septr(separator) >> i_step >>
                              septr(separator),
                          val, separator[0])) {
      if (i_frac_prev != i_frac) {
        if (i_frac != rs.size())
          return error_message("i_step missmatch expected " +
                               std::to_string(rs.size()) +
                               " found:" + std::to_string(i_step));
        if (i_frac > 0) {
          rs.push_back(r);
          r = Simulated_Recording<keep_N_state>{};
        }
        i_step_prev = std::numeric_limits<std::size_t>::max();
      }
      auto &e = get<Recording>(r());
      if (i_step_prev != i_step) {
        if (i_step != e().size())
          return error_message("i_step missmatch expected " +
                               std::to_string(e().size()) +
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
    while ((extract_double(ss >> i_frac >> septr(separator) >> i_step >>
                               septr(separator),
                           val, separator[0]) >>
            septr(separator) >> i_state >> septr(separator) >> N)) {
      if (i_frac_prev != i_frac) {
        if ((i_frac > 0) && (i_frac != rs.size() + 1))
          return error_message(
              "i_frac missmatch expected: " + std::to_string(rs.size()) +
              " found: " + std::to_string(i_frac));
        if (i_frac > 0) {
          auto &Ns = get<N_Ch_State_Evolution>(r());
          Ns().emplace_back(
              N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
          rs.push_back(r);
          N_state.clear();
          r = Simulated_Recording<keep_N_state>{};
        }
        i_frac_prev = i_frac;
      }

      auto &e = get<Recording>(r());
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
          auto &Ns = get<N_Ch_State_Evolution>(r());
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
    auto &Ns = get<N_Ch_State_Evolution>(r());
    Ns().emplace_back(
        N_channel_state(Matrix<double>(1, n_channel_states, N_state)));
    rs.push_back(r);

    return true;
  }
}

inline Maybe_error<bool> load_simulation(std::string const &fname,
                                         std::string separator,
                                         v_Simulated_Recording &r) {

  Simulated_Recording<includes_N_state_evolution(true)> try_N_state;
  auto Maybe_N_state = load_simulation(fname, separator, try_N_state);
  if (Maybe_N_state.valid()) {
    r = try_N_state;
    return true;
  }
  r = Simulated_Recording<includes_N_state_evolution(false)>{};
  return load_simulation(fname, separator, r);
}

template <includes_N_state_evolution keep_N_state>
void save_simulation(std::string name,
                     std::vector<Simulated_Recording<keep_N_state>> const &r) {
  std::ofstream f(name);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  f << "nrep"
    << ","
    << "i_step"
    << ","
    << "patch_current";
  if constexpr (keep_N_state.value) {
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

template <uses_recursive_aproximation recursive,
          uses_averaging_aproximation averaging,
          uses_variance_aproximation variance>
struct MacroR {
  friend std::string ToString(MacroR) {
    std::string out = "MacroR";
    if (recursive.value)
      out += "_R";
    else
      out += "_NR";
    if (averaging.value == 2)
      out += "_2";
    else
      out += "__";
    if (variance.value)
      out += "_V";
    else
      out += "_M";

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

struct Calc_Qx {
  friend std::string ToString(Calc_Qx) { return "Calc_Qx"; }
};

struct Calc_eigen {
  friend std::string ToString(Calc_eigen) { return "Calc_eigen"; }
};

template <typename CQx>
  requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>> full_expm(const CQx &x) {
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

    auto Maybe_E = to_Transition_Probability(eE.value());
    if (!Maybe_E)
      return Maybe_E.error();
    else {

      auto E = std::move(Maybe_E.value());
      // Undo scaling by repeated squaring
      for (std::size_t k = 0; k < s; k++) {
        Maybe_E = to_Transition_Probability(E() * E());
        if (!Maybe_E)
          return Maybe_E.error();
        else
          E = std::move(Maybe_E.value());
      }
      return E;
    }
  }
}

template <typename CQx>
  requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>>
expm_taylor_scaling_squaring(const CQx &x, std::size_t order = 6) {
  {
    double max = maxAbs(primitive(x()));
    double desired = 0.125;
    int k = std::ceil(std::log2(max / desired));
    std::size_t n = std::max(0, k);
    double scale = std::pow(2, -n);
    auto dx = x() * scale;
    auto expm_dx = expm_taylor(dx, order);
    auto Maybe_expm_run = to_Transition_Probability(expm_dx);
    if (!Maybe_expm_run)
      return Maybe_expm_run.error();
    else {
      auto expm_run = std::move(Maybe_expm_run.value());
      for (std::size_t i = 0; i < n; ++i) {
        Maybe_expm_run = to_Transition_Probability(expm_run() * expm_run());
        if (!Maybe_expm_run)
          return Maybe_expm_run.error();
        else
          expm_run = std::move(Maybe_expm_run.value());
      }
      return expm_run;
    }
  }
}

template <typename CQx>
  requires(var::U<CQx, Qx>)
Maybe_error<Transfer_Op_to<CQx, P>> expm_sure(const CQx &x) {
  auto Maybe_expm = full_expm(x);
  if (Maybe_expm)
    return Maybe_expm.value();
  else
    return expm_taylor_scaling_squaring(x);
}

template <class recursive, class averaging, class variance,
          class variance_correction>
struct MacroR2;
class Macro_DMR {
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
  static C_double Ee(C_double const &x, C_double const &y,
                     C_double const &exp_x, C_double const &exp_y,
                     double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(primitive(x) - primitive(y)) < eps)
      return exp_x;
    else
      return (exp_x - exp_y) / (x - y);
  };

  template <class C_double>
    requires U<C_double, double>
  static C_double EX_111(C_double const &x, C_double const &y,
                         C_double const &z, C_double const &exp_x) {
    return exp_x / ((x - y) * (x - z));
  }

  template <class C_double>
    requires U<C_double, double>
  static C_double E111(C_double const &x, C_double const &y, C_double const &z,
                       C_double const &exp_x, C_double const &exp_y,
                       C_double const &exp_z) {
    return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
           EX_111(z, y, x, exp_z);
  }
  template <class C_double>
    requires U<C_double, double>
  static C_double E12(C_double const &x, C_double const &y,
                      C_double const &exp_x, C_double const &exp_y) {
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
  }

  template <class C_double>
    requires U<C_double, double>
  static C_double E3(C_double const &x, C_double const &y, C_double const &z,
                     C_double const &exp_x, C_double const &exp_y,
                     C_double const &exp_z,
                     double eps = std::numeric_limits<double>::epsilon()) {
    auto x_ = primitive(x);
    auto y_ = primitive(y);
    auto z_ = primitive(z);

    if (sqr(x_ - y_) < eps) // x==y
    {
      if (sqr(y_ - z_) < eps) // y==z
        return exp_x / 2.0;   // x==y==z
      else
        return E12(z, x, exp_z, exp_x); // x==y!=z
    } else if (sqr(y_ - z_) < eps)      // x!=y==z
    {
      return E12(x, y, exp_x, exp_y);
    } else if (sqr(x_ - z_) < eps) // y!=z==x!=y
    {
      return E12(y, x, exp_y, exp_x);
    } else
      return E111(x, y, z, exp_x, exp_y, exp_z); // x!=y!=z!=x
  }

  template <bool output>
  static Maybe_error_t<bool>
  test_Probability_value(double e, Probability_error_tolerance tolerance) {
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
  static Maybe_error<bool> test(const C_P_mean &pp,
                                Probability_error_tolerance tolerance) {
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
  static Maybe_error<bool> test(const C_P_Cov &t_p,
                                Probability_error_tolerance tolerance) {
    auto &p = primitive(t_p());

    for (std::size_t i = 0; i < p.nrows(); ++i) {
      for (std::size_t j = 0; j < p.ncols(); ++j) {
        if (auto pijt = test_Probability_value<output>(p(i, i), tolerance);
            !pijt) {
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
             << p << "\n i=" << i << " p(i,j)=" << p(i, i) << " sum=" << sum
             << "\n";
          return error_message(ss.str());
        }
        return error_message("");
      }
    }
    return true;
  }

  template <bool output, class C_P_mean, class C_P_Cov>
    requires(U<C_P_Cov, P_Cov> && U<C_P_mean, P_mean>)
  static Maybe_error<bool> test(const C_P_mean &t_P_mean,
                                const C_P_Cov &t_P_cov,
                                Probability_error_tolerance tolerance) {
    auto ck_mean = test<true>(t_P_mean, tolerance);
    auto ck_cov = test<true>(t_P_cov, tolerance);
    if (ck_mean && ck_cov)
      return true;
    else if constexpr (output) {
      std::stringstream ss;
      ss << " Pmean test: " << ck_mean.error()()
         << " Pcov test: " << ck_cov.error()();
      return error_message(ss.str());
    } else
      return error_message("");
  }

public:
  static bool crude_Qx_violations(Qx const &q) {
    auto Qu = q() * Matrix<double>(q().ncols(), 1ul, 1.0);
    if (maxAbs(Qu) > 1e-7 * norm_inf(q())) {
      return true;
    } else
      return false;
  }

  template <class C_Patch_Model>
    requires U<C_Patch_Model, Patch_Model>
  auto calc_Qx(const C_Patch_Model &m, ATP_concentration x)
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
  auto calc_Qx(const C_Q0 &t_Q0, const C_Qa &t_Qa, ATP_concentration x)
      -> Transfer_Op_to<C_Q0, Qx> {
    auto v_Qx = build<Qx>(t_Q0() + t_Qa() * x.value());
    Matrix<double> u(v_Qx().ncols(), 1, 1.0);
    v_Qx() = v_Qx() - diag(v_Qx() * u);
    assert(!crude_Qx_violations(primitive(v_Qx)));
    //     std::cerr<<"Qx violation\n";
    return v_Qx;
  }

  template <class C_Q0, class C_Qa>
    requires U<C_Q0, Q0>
  Maybe_error<Transfer_Op_to<C_Q0, P_initial>>
  calc_Pinitial(const C_Q0 &t_Q0, const C_Qa &t_Qa, ATP_concentration x,
                N_St nstates) {
    auto p0 = Matrix<double>(1ul, nstates(), 1.0 / nstates());
    auto t_Qx = calc_Qx(t_Q0, t_Qa, x);
    auto v_eig_Qx = calc_eigen(t_Qx);
    if (v_eig_Qx) {

      auto &landa = get<lambda>(v_eig_Qx.value())();
      auto &Vv = get<V>(v_eig_Qx.value())();
      auto &Wv = get<W>(v_eig_Qx.value())();
      auto i_landa = var::i_max(primitive(landa));
      auto ladt = primitive(landa) - primitive(landa);
      ladt[i_landa] = 1.0;

      return build<P_initial>(to_Probability(p0 * Vv * ladt * Wv));

    } else {
      //  std::cerr << "uses expm_sure\n";
      auto Maybe_P = to_Transition_Probability(expm_sure(t_Qx()));
      if (!Maybe_P)
        return Maybe_P.error();
      else {
        auto P = std::move(Maybe_P.value());
        auto Maybe_P2 = to_Transition_Probability(P() * P());
        while (Maybe_P2.valid() &&
               maxAbs(primitive(P() - Maybe_P2.value()())) > 1e-6) {
          P = std::move(Maybe_P2.value());
          Maybe_P2 = to_Transition_Probability(P() * P());
        }
        if (!Maybe_P2)
          return Maybe_P2.error();
        else
          return build<P_initial>(to_Probability(p0 * Maybe_P2.value()()));
      }
    }
  }

  template <class C_Qx>
  auto calc_eigen(const C_Qx &v_Qx)
      -> Maybe_error<Transfer_Op_to<C_Qx, Qx_eig>> {
    auto maybe_eig = eigs(v_Qx());
    if (maybe_eig) {
      auto [v_V, v_l, ignore] = std::move(maybe_eig.value());
      auto Maybe_W = inv(v_V);
      if (!Maybe_W)
        return Maybe_W.error();
      else {
        auto v_W = std::move(Maybe_W.value());
        // auto &la = primitive(v_l);
        // auto i_lambda = var::i_max(la);
        // auto la_dt_inf = DiagonalMatrix<double>(la.nrows(), la.ncols(),
        // 0.0); la_dt_inf[i_lambda] = 1.0; auto Peq =
        // to_Transition_Probability(v_V * la_dt_inf * v_W); if (!Peq)
        //   return Peq.error();
        // else
        return build<Qx_eig>(
            std::move(v_Qx), build<V>(std::move(v_V)),
            build<lambda>(
                to_Transition_Probability_Eigenvalues(std::move(v_l))),
            build<W>(std::move(v_W)));
      }
    } else
      return maybe_eig.error();
  }

  template <class C_Patch_Model>
    requires U<C_Patch_Model, Patch_Model>
  auto calc_eigen(const C_Patch_Model &m, ATP_concentration x)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qx_eig>> {
    return calc_eigen(calc_Qx(m, x));
  }

  template <class C_P_mean>
    requires U<C_P_mean, P_mean>

  static C_P_mean normalize(C_P_mean &&pp, double t_min_p) {
    using Trans = transformation_type_t<C_P_mean>;
    auto p = var::inside_out(pp());
    for (std::size_t i = 0; i < p.nrows(); ++i) {
      Op_t<Trans, double> sum = 0;
      for (std::size_t j = 0; j < p.ncols(); ++j) {
        if (primitive(p(i, j)) > 1.0 - t_min_p) {
          for (std::size_t k = 0; k < p.ncols(); ++k) {
            p(i, k) = (j == k) ? 1.0 + p(i, k) - p(i, k) : p(i, k) - p(i, k);
          }
          return C_P_mean(var::outside_in(p));
        } else if (primitive(p(i, j)) < t_min_p)
          p(i, j) = p(i, j) - p(i, j) + 0.0;
        else
          sum = sum + p(i, j);
      }
      if (primitive(sum) != 1.0)
        for (std::size_t j = 0; j < p.ncols(); ++j)
          p(i, j) = p(i, j) / sum;
    }
    return C_P_mean(var::outside_in(p));
  }

  static auto sample_Multinomial(mt_64i &mt, P_mean const t_P_mean,
                                 std::size_t N) {
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

  static auto sample_Multinomial(mt_64i &mt, P const t_P, N_channel_state N) {
    assert(t_P().nrows() == t_P().ncols());
    auto k = N().size();
    N_channel_state out(Matrix<double>(1, k, 0.0));
    for (std::size_t i = 0; i < k; ++i) {
      std::size_t N_remaining = N()[i];
      double p_remaining = 1;
      for (std::size_t j = 0; j + 1 < k; ++j) {
        if (N_remaining > 0) {
          auto n = std::binomial_distribution<std::size_t>(
              N_remaining, t_P()(i, j) / p_remaining)(mt);
          N_remaining -= n;
          p_remaining -= t_P()(i, j);
          out()[j] += n;
        }
      }
      out()[k - 1] += N_remaining;
    }
    return out;
  }

  template <class C_P_Cov>
    requires U<C_P_Cov, P_Cov>
  static C_P_Cov normalize(C_P_Cov &&p, double t_min_p) {
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
  static C_P normalize(C_P &&p, double t_min_p) {
    using Trans = transformation_type_t<C_P>;

    for (std::size_t i = 0; i < p().nrows(); ++i) {
      Op_t<Trans, double> sumP = 0;
      for (std::size_t j = 0; j < p().ncols(); ++j)
        if (primitive(p()(i, j)) < t_min_p)
          p()(i, j) = 0;
        else
          sumP = sumP + p()(i, j);
      for (std::size_t j = 0; j < p().ncols(); ++j)
        p()(i, j) = p()(i, j) / sumP;
    }
    // std::cerr<<p;
    return std::move(p);
  }

  template <class Vs, class Patch_Model>
    requires Vs::is_vector_map_space
  Maybe_error<Qx_eig const *> get_eigen(Vs &buffer_calc, const Patch_Model &m,
                                        ATP_concentration x) {
    auto Maybe_eigen = buffer_calc[var::Vector_Map<Qx_eig>{}]
                                  [var::Vector_Space<ATP_concentration>(x)];
    if (Maybe_eigen)
      return Maybe_eigen;
    else {
      auto Maybe_new_eigen = calc_eigen(m, x);
      if (Maybe_new_eigen) {
        buffer_calc[var::Vector_Map<Qx_eig>{}].emplace(
            x, std::move(Maybe_new_eigen.value()));
        return get_eigen(buffer_calc, m, x);
      } else
        return Maybe_new_eigen.error();
    }
  }

  template <class C_Patch_Model, class C_Qx_eig>
    requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx_eig, Qx_eig>)
  auto calc_Peq_(C_Qx_eig const &t_Qx, const C_Patch_Model &m)
      -> Transfer_Op_to<C_Patch_Model, P_mean> {
    auto nstates = get<N_St>(m).value();
    auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);

    auto &landa = get<lambda>(t_Qx)();
    auto &Vv = get<V>(t_Qx)();
    auto &Wv = get<W>(t_Qx)();
    auto ladt = get<lambda>(t_Qx)() * 1e8;
    auto laexp = apply(
        [](auto const &x) {
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
      std::cerr << "\nQx\n" << get<Qx>(t_Qx);
    }

    return build<P_mean>(p0 * Vv * laexp * Wv);
  }

  template <class C_Qx>
    requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx, Qx>)
  auto calc_Peq(C_Qx const &t_Qx, N_St nstates)
      -> Transfer_Op_to<C_Qx, P_mean> {
    auto p0 = Matrix<double>(1ul, nstates(), 1.0 / nstates());
    auto v_eig_Qx = calc_eigen(t_Qx);
    if (v_eig_Qx) {

      auto &landa = get<lambda>(v_eig_Qx.value())();
      auto &Vv = get<V>(v_eig_Qx.value())();
      auto &Wv = get<W>(v_eig_Qx.value())();
      auto ladt = get<lambda>(v_eig_Qx.value())() * 1e8;

      auto laexp = apply(
          [](auto const &x) {
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
        std::cerr << "\nQx\n" << get<Qx>(t_Qx);
      }

      return build<P_mean>(p0 * Vv * laexp * Wv);

    } else {
      // std::cerr << "uses expm_sure\n";
      auto P = expm_sure(t_Qx());
      auto P2 = P * P;
      while (maxAbs(primitive(P - P2)) > 1e-9) {
        P = P2;
        P2 = P * P;
      }
      return build<P_mean>(p0 * P2);
    }
  }

  // template <class C_Patch_Model, class C_Qx_eig>
  //     requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx_eig, Qx_eig>)
  // auto calc_Peq(C_Qx_eig const &t_Qx, const C_Patch_Model &m)
  //     -> Transfer_Op_to<C_Patch_Model, P_mean> {
  //     auto nstates = get<N_St>(m).value();
  //     auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);

  //     auto &landa = get<lambda>(t_Qx)();
  //     auto &Vv = get<V>(t_Qx)();
  //     auto &Wv = get<W>(t_Qx)();
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
  //         std::cerr << "\nQx\n" << get<Qx>(t_Qx);
  //     }

  //     return build<P_mean>(p0 * Vv * laexp * Wv);
  // }

  template <class Patch_Model>
  auto calc_P(const Patch_Model &m, const Qx_eig &t_Qx, double dt,
              double t_min_P) {
    auto ladt = get<lambda>(t_Qx)() * dt;

    auto exp_ladt = apply([](double x) { return std::exp(x); }, ladt);
    //    return normalize(P(get<V>(t_Qx)() * exp_ladt * get<W>(t_Qx)()),
    //    t_min_P);
    return to_Transition_Probability(get<V>(t_Qx)() * exp_ladt *
                                     get<W>(t_Qx)());
  }

  template <class Patch_Model>
  auto calc_P(const Patch_Model &m, const Qx &t_Qx, double dt, double t_min_P) {
    auto t_eigenQx = calc_eigen(t_Qx);
    if (t_eigenQx) {
      auto Maybe_P = calc_P(m, t_eigenQx.value(), dt, t_min_P);
      if (Maybe_P)
        return Maybe_P;
    }
    //      return normalize(P(expm_sure(t_Qx() * dt)), t_min_P);
    return to_Transition_Probability(expm_sure(t_Qx() * dt));
  }

  template <class Patch_Model>
  auto calc_Qdt_old(const Patch_Model &m, const Qx_eig &t_Qx,
                    number_of_samples ns, double dt) {

    auto t_min_P = get<min_P>(m)();
    auto &v_g = get<g>(m);

    std::size_t N = t_Qx[Var<Qx>{}]().ncols();

    auto ladt = t_Qx[Var<lambda>{}]() * dt;

    auto exp_ladt = apply(
        [](auto const &x) {
          using std::exp;
          return exp(x);
        },
        ladt);
    auto v_P = P(get<V>(t_Qx)() * exp_ladt * get<W>(t_Qx)());

    SymmetricMatrix<double> E2m(N, N);
    SymmetricMatrix<double> E2mb(N, N);
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < i + 1; ++j)
        E2m.set(i, j, Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j], t_min_P));

    // build E2
    Matrix<double> WgV_E2(N, N);
    Matrix<double> WgV = get<W>(t_Qx)() * diag(get<g>(m)()) * get<V>(t_Qx)();

    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = WgV(i, j) * E2m(i, j);

    auto v_gtotal_ij = gtotal_ij(t_Qx[Var<V>{}]() * WgV_E2 * t_Qx[Var<W>{}]());

    Matrix<double> WgV_E3(N, N, 0.0);
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
          WgV_E3(n1, n3) +=
              WgV(n1, n2) * WgV(n2, n3) *
              E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1], exp_ladt[n2],
                 exp_ladt[n3], t_min_P); // optimizable
        }

    auto v_gtotal_sqr_ij =
        gtotal_sqr_ij(t_Qx[Var<V>{}]() * WgV_E3 * t_Qx[Var<W>{}]() * 2.0);
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        if (v_P()(i, j) == 0) {
          v_gtotal_ij()(i, j) = 0;
          v_gtotal_sqr_ij()(i, j) = 0;
        }

    auto U = Matrix<double>(1, N, 1.0);
    auto UU = Matrix<double>(N, N, 1.0);
    auto gmean_ij_p = X_plus_XT(v_g() * U) * (0.5);
    auto gvar_ij_p = (v_g() * U - apply([](double x) { return std::abs(x); },
                                        tr(v_g() * U))) *
                     (0.5);

    std::cerr << "\ngmean_ij_p=\n"
              << gmean_ij_p << "\ngvar_ij_p=\n"
              << gvar_ij_p << "\n";
    // std::cerr<<"\n UU=ņ"<<UU<<"\n";
    auto gmean_ij_tot = v_gtotal_ij() + gmean_ij_p * t_min_P;
    auto P_p = v_P() + UU * t_min_P;
    auto v_gmean_ij =
        gmean_ij(zip([](auto x, auto y) { return x / y; }, gmean_ij_tot, P_p));
    auto v_gtotal_var_ij = gtotal_var_ij(
        v_gtotal_sqr_ij() -
        zip([](auto x, auto y) { return x * y; }, v_gtotal_ij(), v_gmean_ij()));
    auto gvar_ij_tot = v_gtotal_var_ij() + gvar_ij_p * t_min_P;
    auto v_gvar_ij =
        gvar_ij(zip([](auto x, auto y) { return x / y; }, gvar_ij_tot, P_p));
    Matrix<double> u(N, 1, 1.0);
    auto v_gmean_i = gmean_i(v_gtotal_ij() * u);
    auto v_gsqr_i = gsqr_i(v_gtotal_sqr_ij() * u);
    auto v_gvar_i = gvar_i(v_gtotal_var_ij() * u);

    return Qdt(ns, min_P(t_min_P), std::move(v_P), std::move(v_gmean_i),
               std::move(v_gtotal_ij), std::move(v_gmean_ij),
               std::move(v_gtotal_sqr_ij), std::move(v_gsqr_i),
               std::move(v_gvar_i), std::move(v_gtotal_var_ij),
               std::move(v_gvar_ij));
  }

  template <class FunctionTable, class C_Patch_Model, class C_Qx_eig>
    requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Qx_eig>)
  Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>>
  calc_Qdtm_eig(FunctionTable &&, const C_Patch_Model &m, const C_Qx_eig &t_Qx,
                number_of_samples ns, double dt) {
    using Trans = transformation_type_t<C_Patch_Model>;
    // const double eps=std::numeric_limits<double>::epsilon();
    auto &t_V = get<V>(t_Qx);
    auto &t_landa = get<lambda>(t_Qx);
    auto &t_W = get<W>(t_Qx);
    auto &t_g = get<g>(m);
    auto t_min_P = get<min_P>(m);
    auto v_ladt = t_landa() * dt;
    auto v_exp_ladt = apply(
        [](auto const &x) {
          using std::exp;
          return exp(x);
        },
        v_ladt);

    auto Maybe_r_P = to_Transition_Probability(t_V() * v_exp_ladt * t_W());
    if (!Maybe_r_P)
      return Maybe_r_P.error();
    else {
      auto r_P = std::move(Maybe_r_P.value());

      std::size_t N = r_P().ncols();

      SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
      for (std::size_t i = 0; i < N; ++i) {

        for (std::size_t j = 0; j < i + 1; ++j) {
          set(E2m, i, j,
              Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j],
                 t_min_P()));
        }
      }

      Matrix<Op_t<Trans, double>> WgV_E2(N, N);

      auto v_WgV = t_W() * diag(t_g()) * t_V();

      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
          WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

      auto r_gtotal_ij =
          force_gmean_in_range(build<gtotal_ij>(t_V() * WgV_E2 * t_W()), t_g);

      auto r_gmean_ij = force_gmean_in_range(
          build<gmean_ij>(elemDivSafe(r_gtotal_ij(), r_P(), t_min_P())), t_g);
      /* truncate is not derivative safe yet*/

      Matrix<double> u(N, 1, 1.0);
      auto r_gmean_i =
          force_gmean_in_range(build<gmean_i>(r_gtotal_ij() * u), t_g);
      // if (crude_gmean_violation(primitive(r_gmean_i),
      // primitive(get<g>(m))))
      //     std::cerr<<"gmean_violation\n";

      Matrix<Op_t<Trans, double>> WgV_Wg_E2(N, N, 0.0);

      auto v_Wg = t_W() * t_g();

      auto v_eps = eps;
      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
          WgV_Wg_E2(i, j) = v_WgV(i, j) * E2m(i, j) * v_Wg[j];
      for (std::size_t k0 = 0; k0 < N; k0++) {
        auto rladt = v_ladt[k0];
        if (rladt * rladt > v_eps) {
          for (std::size_t k2 = 0; k2 < N; k2++) {
            auto rla2dt = v_ladt[k2];
            if (rla2dt * rla2dt > v_eps)
              WgV_Wg_E2(k0, k2) = v_WgV(k0, k2) * v_Wg[k2] * 0.5;
            else
              WgV_Wg_E2(k0, k2) = v_WgV(k0, k2) * v_Wg[k2] *
                                  (v_exp_ladt[k2] - rla2dt - 1.0) / rla2dt /
                                  rla2dt;
          }
        } else {
          for (std::size_t k2 = 0; k2 < N; k2++) {
            double rla2dt = v_ladt[k2];
            if (rla2dt * rla2dt > v_eps) {
              WgV_Wg_E2(k0, k2) = v_WgV(k0, k2) * v_Wg[k2] *
                                  (v_exp_ladt[k0] - rla2dt - 1.0) / rla2dt /
                                  rla2dt;
            } else if ((rla2dt - rladt) * (rla2dt - rladt) >
                       v_eps) // comparing squared difference
            {
              WgV_Wg_E2(k0, k2) = v_WgV(k0, k2) * v_Wg[k2] *
                                  (1.0 - v_exp_ladt[k0] * (1.0 - rladt)) /
                                  rladt / rladt;
            } else {
              WgV_Wg_E2(k0, k2) = v_WgV(k0, k2) * v_Wg[k2] *
                                  (1.0 / rladt / rla2dt +
                                   v_exp_ladt[k2] / rla2dt / (rla2dt - rladt) +
                                   v_exp_ladt[k0] / rladt / (rladt - rla2dt));
            }
          }
        }
      }

      Matrix<Op_t<Trans, double>> rgsqr_i(N, 1, 0.0);

      for (std::size_t i = 0; i < N; i++) {
        for (std::size_t k0 = 0; k0 < N; k0++)
          for (std::size_t k2 = 0; k2 < N; k2++)
            rgsqr_i[i] += 2 * t_V()(i, k0) * WgV_Wg_E2(k0, k2);
      }

      auto r_gsqr_i = build<gsqr_i>(rgsqr_i);

      auto r_gvar_i =
          build<gvar_i>(r_gsqr_i() - elemMult(r_gmean_i(), r_gmean_i()));

      /* truncate is not derivative safe yet*/

      if constexpr (true) {
        r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
      }
      //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
      //   primitive(t_g())); auto test_g_mean =
      //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
      // if (!test_g_mean || !test_g_var)
      //     return error_message(test_g_var.error()()+test_g_mean.error()());
      // else {

      return build<Qdtm>(ns, min_P(t_min_P), std::move(r_P),
                         std::move(r_gmean_i), std::move(r_gtotal_ij),
                         std::move(r_gmean_ij), std::move(r_gsqr_i),
                         std::move(r_gvar_i));
      // }
    }
  }

  template <class FunctionTable, class C_Patch_Model, class C_Qx_eig>
    requires(U<C_Patch_Model, Patch_Model> && U<C_Qx_eig, Qx_eig>)
  Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>
  calc_Qdt_eig(FunctionTable &&, const C_Patch_Model &m, const C_Qx_eig &t_Qx,
               number_of_samples ns, double dt) {
    using Trans = transformation_type_t<C_Patch_Model>;
    // const double eps=std::numeric_limits<double>::epsilon();
    auto &t_V = get<V>(t_Qx);
    auto &t_landa = get<lambda>(t_Qx);
    auto &t_W = get<W>(t_Qx);
    auto &t_g = get<g>(m);
    auto t_min_P = get<min_P>(m);
    auto v_ladt = t_landa() * dt;
    auto v_exp_ladt = apply(
        [](auto const &x) {
          using std::exp;
          return exp(x);
        },
        v_ladt);

    auto Maybe_r_P = to_Transition_Probability(t_V() * v_exp_ladt * t_W());
    if (!Maybe_r_P)
      return Maybe_r_P.error();
    else {
      auto r_P = std::move(Maybe_r_P.value());

      std::size_t N = r_P().ncols();

      SymmetricMatrix<Op_t<Trans, double>> E2m(N, N);
      for (std::size_t i = 0; i < N; ++i) {

        for (std::size_t j = 0; j < i + 1; ++j) {
          set(E2m, i, j,
              Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j],
                 t_min_P()));
        }
      }

      Matrix<Op_t<Trans, double>> WgV_E2(N, N);

      auto v_WgV = t_W() * diag(t_g()) * t_V();

      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
          WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

      auto r_gtotal_ij =
          force_gmean_in_range(build<gtotal_ij>(t_V() * WgV_E2 * t_W()), t_g);

      Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
      for (std::size_t n1 = 0; n1 < N; n1++)
        for (std::size_t n3 = 0; n3 < N; n3++)
          for (std::size_t n2 = 0; n2 < N; n2++) {
            //      std::cerr<<"\t"<<WgV_E3(n1, n3);

            WgV_E3(n1, n3) = WgV_E3(n1, n3) +
                             v_WgV(n1, n2) * v_WgV(n2, n3) *
                                 E3(v_ladt[n1], v_ladt[n2], v_ladt[n3],
                                    v_exp_ladt[n1], v_exp_ladt[n2],
                                    v_exp_ladt[n3], t_min_P()); 
      }

      auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);

      if constexpr (false) {
        std::cerr << "\nr_gtotal_sqr_ij\n" << r_gtotal_sqr_ij;
        std::cerr << "\nvar::outside_in(var::inside_out(r_gtotal_sqr_ij))\n"
                  << var::outside_in(var::inside_out(r_gtotal_sqr_ij()));

        std::cerr << "\nvar::inside_out(r_gtotal_sqr_ij)\n"
                  << var::inside_out(r_gtotal_sqr_ij());
      }

      if constexpr (true) {
        r_gtotal_sqr_ij() =
            truncate_negative_variance(std::move(r_gtotal_sqr_ij()));
        for (std::size_t i = 0; i < N; ++i)
          for (std::size_t j = 0; j < N; ++j)
            if (r_P()(i, j) == 0) {
              r_gtotal_ij().set(i, j, 0.0);
              r_gtotal_sqr_ij().set(i, j, 0.0);
            }
      }

      auto r_gmean_ij = force_gmean_in_range(
          build<gmean_ij>(elemDivSafe(r_gtotal_ij(), r_P(), t_min_P())), t_g);
      auto r_gtotal_var_ij = build<gtotal_var_ij>(
          r_gtotal_sqr_ij() - elemMult(r_gtotal_ij(), r_gmean_ij()));

      /* truncate is not derivative safe yet*/
      if constexpr (true) {
        r_gtotal_var_ij() =
            truncate_negative_variance(std::move(r_gtotal_var_ij()));
      }

      auto r_gvar_ij =
          build<gvar_ij>(elemDivSafe(r_gtotal_var_ij(), r_P(), t_min_P()));

      Matrix<double> u(N, 1, 1.0);
      auto r_gmean_i =
          force_gmean_in_range(build<gmean_i>(r_gtotal_ij() * u), t_g);
      // if (crude_gmean_violation(primitive(r_gmean_i),
      // primitive(get<g>(m))))
      //     std::cerr<<"gmean_violation\n";

      auto r_gsqr_i = build<gsqr_i>(r_gtotal_sqr_ij() * u);
      auto r_gvar_i = build<gvar_i>(r_gtotal_var_ij() * u);
      if constexpr (true) {
        r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
      }
      //   auto test_g_var = test_conductance_variance(primitive(r_gvar_i()),
      //   primitive(t_g())); auto test_g_mean =
      //   test_conductance_mean(primitive(r_gmean_i()), primitive(t_g()));
      // if (!test_g_mean || !test_g_var)
      //     return error_message(test_g_var.error()()+test_g_mean.error()());
      // else {

      return build<Qdt>(ns, min_P(t_min_P), std::move(r_P),
                        std::move(r_gmean_i), std::move(r_gtotal_ij),
                        std::move(r_gmean_ij), std::move(r_gtotal_sqr_ij),
                        std::move(r_gsqr_i), std::move(r_gvar_i),
                        std::move(r_gtotal_var_ij), std::move(r_gvar_ij));
      // }
    }
  }

  template <class C_Patch_Model, class C_Qx>
    requires(/*U<C_Patch_Model, Patch_Model> && */ U<C_Qx, Qx>)
  Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>>
  calc_Qdt_taylor(const C_Patch_Model &m, const C_Qx &t_Qx,
                  number_of_samples ns, double dt, std::size_t order = 5ul) {
    auto v_Qrun = t_Qx() * dt;
    double max = maxAbs(primitive(v_Qrun));
    double desired = 0.125;
    int k = std::ceil(std::log2(max / desired));
    std::size_t n = std::max(0, k);
    double scale = std::pow(2, -n);
    auto t_Qrun_sub = v_Qrun * scale;
    auto Maybe_P_sub =
        to_Transition_Probability(expm_taylor(t_Qrun_sub, order));
    if (!Maybe_P_sub) {
      return Maybe_P_sub.error();
    } else {
      auto P_sub = std::move(Maybe_P_sub.value());
      auto r_Qn = get_Qn(P_sub, get<g>(m), ns, get<min_P>(m));
      for (std::size_t i = 0; i < n; ++i) {
        r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
      }
      get<number_of_samples>(r_Qn) = ns;
      return Qn_to_Qdt(r_Qn);
    }
  }

  template <class C_Patch_Model, class C_Qx>
    requires(/*U<C_Patch_Model, Patch_Model> && */ U<C_Qx, Qx>)
  Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>>
  calc_Qdtm_taylor(const C_Patch_Model &m, const C_Qx &t_Qx,
                   number_of_samples ns, double dt, std::size_t order = 5ul) {
    auto v_Qrun = t_Qx() * dt;
    double max = maxAbs(primitive(v_Qrun));
    double desired = 0.125;
    int k = std::ceil(std::log2(max / desired));
    std::size_t n = std::max(0, k);
    double scale = std::pow(2, -n);
    auto t_Qrun_sub = v_Qrun * scale;
    auto Maybe_P_sub =
        to_Transition_Probability(expm_taylor(t_Qrun_sub, order));
    if (!Maybe_P_sub) {
      return Maybe_P_sub.error();
    } else {
      auto P_sub = std::move(Maybe_P_sub.value());
      auto r_Qn = get_Qn(P_sub, get<g>(m), ns, get<min_P>(m));
      for (std::size_t i = 0; i < n; ++i) {
        r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
      }
      get<number_of_samples>(r_Qn) = ns;
      return Qn_to_Qdtm(r_Qn);
    }
  }

  template <class C_Qdt>
    requires(U<C_Qdt, Qdt>)
  auto get_Qn(const C_Qdt &x) {
    auto n = get<number_of_samples>(x)();
    return build<Qn>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                     build<PG_n>(get<gtotal_ij>(x)() * n),
                     build<PGG_n>(get<gtotal_sqr_ij>(x)() * (n * n * 0.5)));
  }

  template <class C_P, class C_g>
    requires(U<C_P, P> && U<C_g, g>)
  auto get_Qn(const C_P &t_P, C_g const &t_g, number_of_samples n,
              min_P t_minP) {

    auto N = t_P().nrows();
    auto u = Matrix<double>(1, N, 1.0);
    auto G = t_g() * u;
    auto GT = tr(G);
    auto Gmean = 0.5 * G + 0.5 * GT;

    auto Gvar = elemMult(G, GT) - elemMult(Gmean, Gmean);

    return build<Qn>(n, t_minP, t_P, build<PG_n>(elemMult(t_P(), Gmean) * n()),
                     build<PGG_n>(elemMult(t_P(), Gvar) * (n() * n() * 0.5)));
  }

  template <class C_Qn>
    requires(U<C_Qn, Qn>)
  static C_Qn sum_Qn(C_Qn &&one, const C_Qn &two) {
    auto n1 = get<number_of_samples>(two)();
    get<PGG_n>(one)() = (get<PGG_n>(one)() * get<P>(two)()) +
                        (get<PG_n>(one)() * get<PG_n>(two)()) +
                        (get<P>(one)() * get<PGG_n>(two)());
    get<PG_n>(one)() =
        (get<PG_n>(one)() * get<P>(two)()) + (get<P>(one)() * get<PG_n>(two)());
    get<P>(one) =
        to_Transition_Probability(get<P>(one)() * get<P>(two)()).value();
    get<number_of_samples>(one)() = get<number_of_samples>(one)() + n1;
    return one;
  }

  template <class C_Qn, class C_Qdt>
    requires(U<C_Qn, Qn> && U<C_Qdt, Qdt>)
  static C_Qn sum_Qdt(C_Qn &&one, const C_Qdt &two) {
    auto n1 = get<number_of_samples>(two)();
    get<PGG_n>(one)() =
        (get<PGG_n>(one)() * get<P>(two)()) +
        (get<PG_n>(one)() * get<gtotal_ij>(two)()) * n1 +
        (get<P>(one)() * get<gtotal_sqr_ij>(two)()) * (0.5 * n1 * n1);
    get<PG_n>(one)() = (get<PG_n>(one)() * get<P>(two)()) +
                       (get<P>(one)() * get<gtotal_ij>(two)()) * n1;
    get<P>(one) =
        to_Transition_Probability(get<P>(one)() * get<P>(two)()).value();
    get<number_of_samples>(one)() = get<number_of_samples>(one)() + n1;
    return one;
  }

  static bool is_Binomial_Approximation_valid(double N, double p, double q,
                                              double Np_min) {
    if (N * p < Np_min)
      return false;
    else if (N * q < Np_min)
      return false;
    else
      return true;
  }

  static y_mean max_possible_value_of_ymean(N_Ch_mean_value t_N, const g &t_g,
                                            Current_Baseline b) {
    return y_mean(t_N() * var::max(t_g()) + b());
  }

  static y_mean min_possible_value_of_ymean(N_Ch_mean_value t_N, g t_g,
                                            Current_Baseline b) {
    return y_mean(t_N() * var::min(t_g()) + b());
  }

  static bool crude_gmean_violation(gmean_i const &v_gm, const g &v_g) {
    auto max_g_m = var::max(v_gm());
    auto min_g_m = var::min(v_gm());
    auto max_g = var::max(v_g());
    auto min_g = var::min(v_g());
    if ((max_g_m <= max_g) && (min_g_m >= min_g))
      return false;
    else
      return true;
  }

  template <class C_gmean, class C_g>
    requires((U<C_gmean, gmean_i> || U<C_gmean, gtotal_ij> ||
              U<C_gmean, gmean_ij>) &&
             U<C_g, g>)
  auto force_gmean_in_range(C_gmean &&g_mean, const C_g &v_g) {
    auto gmax = var::max(primitive(v_g()));
    auto gmin = var::min(primitive(v_g()));
    for (std::size_t i = 0; i < g_mean().size(); ++i) {
      if (primitive(g_mean()[i]) < gmin)
        g_mean().set(i, gmin);
      else if (primitive(g_mean()[i]) > gmax)
        g_mean().set(i, gmax);
    }
    return std::move(g_mean);
  }

  template <class C_Qn>
    requires(U<C_Qn, Qn>)
  static auto Qn_to_Qdt(const C_Qn &x) {
    auto u = Matrix<double>(get<P>(x)().ncols(), 1ul, 1.0);
    auto r_P = get<P>(x)();
    auto n = get<number_of_samples>(x)();
    auto r_gtotal_sqr_ij = get<PGG_n>(x)() * (2.0 / (n * n));
    auto r_gtotal_ij = get<PG_n>(x)() * (1.0 / n);
    auto r_gmean_ij = elemDivSafe(r_gtotal_ij, get<P>(x)(), get<min_P>(x)());
    auto r_gtotal_var_ij = r_gtotal_sqr_ij - elemMult(r_gtotal_ij, r_gmean_ij);
    auto r_gmean_i = r_gtotal_ij * u;
    auto r_gsqr_i = r_gtotal_sqr_ij * u;
    auto r_gvar_ij = elemDivSafe(r_gtotal_var_ij, r_P, get<min_P>(x)());
    auto r_gvar_i = r_gtotal_var_ij * u;

    return build<Qdt>(
        get<number_of_samples>(x), get<min_P>(x), get<P>(x),
        build<gmean_i>(r_gmean_i), build<gtotal_ij>(r_gtotal_ij),
        build<gmean_ij>(r_gmean_ij), build<gtotal_sqr_ij>(r_gtotal_sqr_ij),
        build<gsqr_i>(r_gsqr_i), build<gvar_i>(r_gvar_i),
        build<gtotal_var_ij>(r_gtotal_var_ij), build<gvar_ij>(r_gvar_ij));
  }

  template <class C_Qn>
    requires(U<C_Qn, Qn>)
  static auto Qn_to_Qdtm(const C_Qn &x) {
    auto u = Matrix<double>(get<P>(x)().ncols(), 1ul, 1.0);
    auto r_P = get<P>(x)();
    auto n = get<number_of_samples>(x)();
    auto r_gtotal_sqr_ij = get<PGG_n>(x)() * (2.0 / (n * n));
    auto r_gtotal_ij = get<PG_n>(x)() * (1.0 / n);
    auto r_gmean_ij = elemDivSafe(r_gtotal_ij, get<P>(x)(), get<min_P>(x)());
    auto r_gtotal_var_ij = r_gtotal_sqr_ij - elemMult(r_gtotal_ij, r_gmean_ij);
    auto r_gmean_i = r_gtotal_ij * u;
    auto r_gsqr_i = r_gtotal_sqr_ij * u;
    auto r_gvar_i = r_gtotal_var_ij * u;

    return build<Qdtm>(get<number_of_samples>(x), get<min_P>(x), get<P>(x),
                       build<gmean_i>(r_gmean_i), build<gtotal_ij>(r_gtotal_ij),
                       build<gmean_ij>(r_gmean_ij), build<gsqr_i>(r_gsqr_i),
                       build<gvar_i>(r_gvar_i));
  }

  template <class FunctionTable, class C_Patch_Model>
    requires(!is_of_this_template_type_v<FunctionTable, FuncMap_St>)
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_ATP_step(FunctionTable &&f, const C_Patch_Model &m,
                         const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    auto dt = get<number_of_samples>(t_step)() / fs;
    auto t_Qeig = f.fstop(Calc_eigen{}, m, get<ATP_concentration>(t_step));
    if constexpr (false) {
      auto test_der_eigen = var::test_Derivative(
          [this, &t_step](auto l_m) {
            return calc_eigen(l_m, get<ATP_concentration>(t_step));
          },
          1, 1e-9, m);
      if (!test_der_eigen) {
        std::cerr << test_der_eigen.error()();
        return test_der_eigen.error();
      }
    }

    if (t_Qeig) {
      auto Maybe_Qdt = calc_Qdt_eig(f, m, t_Qeig.value(),
                                    get<number_of_samples>(t_step), dt);
      if (Maybe_Qdt)
        return Maybe_Qdt;
    }
    auto t_Qx = build<Qx>(calc_Qx(m, get<ATP_concentration>(t_step)));
    return calc_Qdt_taylor(m, t_Qx, get<number_of_samples>(t_step), dt);
  }

  template <class FunctionTable, class C_Patch_Model>
    requires(!is_of_this_template_type_v<FunctionTable, FuncMap_St>)
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdtm_ATP_step(FunctionTable &&f, const C_Patch_Model &m,
                          const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
    auto dt = get<number_of_samples>(t_step)() / fs;
    auto t_Qeig = f.fstop(Calc_eigen{}, m, get<ATP_concentration>(t_step));

    if (t_Qeig) {
      auto Maybe_Qdt = calc_Qdtm_eig(f, m, t_Qeig.value(),
                                     get<number_of_samples>(t_step), dt);
      if (Maybe_Qdt)
        return Maybe_Qdt.value();
    }
    auto t_Qx = build<Qx>(calc_Qx(m, get<ATP_concentration>(t_step)));
    return calc_Qdtm_taylor(m, t_Qx, get<number_of_samples>(t_step), dt);
  }

  template <class FunctionTable, class C_Patch_Model>
    requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_ATP_step(FunctionTable &f, const C_Patch_Model &m,
                         const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    auto dt = get<number_of_samples>(t_step)() / fs;
    auto t_Qeig = f.fstop(Calc_eigen{}, m, get<ATP_concentration>(t_step));
    if constexpr (false) {
      auto test_der_eigen = var::test_Derivative(
          [this, &t_step](auto l_m) {
            return calc_eigen(l_m, get<ATP_concentration>(t_step));
          },
          1, 1e-9, m);
      if (!test_der_eigen) {
        std::cerr << test_der_eigen.error()();
        return test_der_eigen.error();
      }
    }

    if (t_Qeig) {
      auto Maybe_Qdt = calc_Qdt_eig(f, m, t_Qeig.value(),
                                    get<number_of_samples>(t_step), dt);
      if (Maybe_Qdt)
        return Maybe_Qdt;
    }
    auto t_Qx = build<Qx>(calc_Qx(m, get<ATP_concentration>(t_step)));
    return calc_Qdt_taylor(m, t_Qx, get<number_of_samples>(t_step), dt);
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt(FunctionTable &&f, const C_Patch_Model &m,
                const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    if constexpr (std::is_same_v<Nothing, decltype(f[Calc_Qdt_step{}])>)
      return calc_Qdt_ATP_step(f, m, t_step, fs);
    else
      return f.f(Calc_Qdt_step{}, m, t_step, fs);
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdtm(FunctionTable &&f, const C_Patch_Model &m,
                 const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
    if constexpr (std::is_same_v<Nothing, decltype(f[Calc_Qdt_step{}])>)
      return calc_Qdtm_ATP_step(f, m, t_step, fs);
    else
      return f.f(Calc_Qdtm_step{}, m, t_step, fs);
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qn_bisection(FunctionTable &&f, const C_Patch_Model &m,
                         const ATP_step &t_step, double fs, std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qn>> {
    auto dt = get<number_of_samples>(t_step)() / fs;
    auto ns = get<number_of_samples>(t_step);
    auto t_Qx = f.fstop(Calc_eigen{}, m, get<ATP_concentration>(t_step));

    if (!t_Qx)
      return t_Qx.error();
    else {
      double scale = std::pow(2.0, -1.0 * order);

      number_of_samples n_ss(ns() * scale);
      double sdt = dt * scale;
      auto t_Psub = calc_P(m, t_Qx.value(), sdt, get<min_P>(m)() * scale);
      auto r_Qn =
          get_Qn(t_Psub, get<g>(m), n_ss, min_P(get<min_P>(m)() * scale));
      for (std::size_t i = 0; i < order; ++i) {
        r_Qn = sum_Qn(std::move(r_Qn), r_Qn);
      }
      assert(get<number_of_samples>(r_Qn)() == ns());
      return r_Qn;
    }
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_bisection(FunctionTable &&f, const C_Patch_Model &m,
                          const ATP_step &t_step, double fs, std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    auto maybe_Qn =
        calc_Qn_bisection(std::forward<FunctionTable>(f), m, t_step, fs, order);
    if (!maybe_Qn)
      return maybe_Qn.error();
    else {
      return Qn_to_Qdt(maybe_Qn.value());
    }
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model> )
  auto calc_Qdt(FunctionTable &&f, const C_Patch_Model &m,
                const std::vector<ATP_step> &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    if (t_step.empty())
      return error_message("Emtpy ATP step");
    else {
      auto v_Qdt0 = calc_Qdt(std::forward<FunctionTable>(f), m, t_step[0], fs);
      if (!v_Qdt0)
        return v_Qdt0.error();
      else {
        auto v_Qrun = get_Qn(v_Qdt0.value());
        for (std::size_t i = 1; i < t_step.size(); ++i) {
          auto v_Qdti =
              calc_Qdt(std::forward<FunctionTable>(f), m, t_step[i], fs);
          if (!v_Qdti)
            return v_Qdti.error();
          else
            v_Qrun = sum_Qdt(std::move(v_Qrun), v_Qdti.value());
        }
        return Qn_to_Qdt(v_Qrun);
      }
    }
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model> )
  auto calc_Qdtm(FunctionTable &&f, const C_Patch_Model &m,
                 const std::vector<ATP_step> &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
    if (t_step.empty())
      return error_message("Emtpy ATP step");
    if (t_step.size() == 1)
      return calc_Qdtm(std::forward<FunctionTable>(f), m, t_step[0], fs);

    auto v_Qdt0 = calc_Qdt(std::forward<FunctionTable>(f), m, t_step[0], fs);
    if (!v_Qdt0)
      return v_Qdt0.error();
    auto v_Qrun = get_Qn(v_Qdt0.value());
    for (std::size_t i = 1; i < t_step.size(); ++i) {
      auto v_Qdti = calc_Qdt(std::forward<FunctionTable>(f), m, t_step[i], fs);
      if (!v_Qdti)
        return v_Qdti.error();
      else
        v_Qrun = sum_Qdt(std::move(v_Qrun), v_Qdti.value());
    }
    return Qn_to_Qdtm(v_Qrun);
  }

  template <class FunctionTable, class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model> )
  auto calc_Qdt_bisection(FunctionTable &&f, const C_Patch_Model &m,
                          const std::vector<ATP_step> &t_step, double fs,
                          std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    if (t_step.empty())
      return error_message("Emtpy ATP step");
    else {
      auto v_Qn0 = calc_Qn_bisection(std::forward<FunctionTable>(f), m,
                                     t_step[0], fs, order);
      if (!v_Qn0)
        return v_Qn0.error();
      else {
        auto v_Qrun = v_Qn0.value();
        for (std::size_t i = 1; i < t_step.size(); ++i) {
          auto v_Qni = calc_Qn_bisection(std::forward<FunctionTable>(f), m,
                                         t_step[i], fs, order);
          if (!v_Qni)
            return v_Qni.error();
          else
            v_Qrun = sum_Qn(std::move(v_Qrun), v_Qni.value());
        }
        return Qn_to_Qdt(v_Qrun);
      }
    }
  }

  template <class FunctionTable, class C_Patch_Model>
  //   requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt(FunctionTable &&f, const C_Patch_Model &m,
                const ATP_evolution &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    return calc_Qdt(std::forward<FunctionTable>(f), m, t_step(), fs);
  }

  template <class FunctionTable, class C_Patch_Model>
  //   requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdtm(FunctionTable &&f, const C_Patch_Model &m,
                const ATP_evolution &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdtm>> {
    return calc_Qdtm(std::forward<FunctionTable>(f), m, t_step(), fs);
  }
  
  template <class FunctionTable, class C_Patch_Model>
  //   requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_bisection(FunctionTable &&f, const C_Patch_Model &m,
                          const ATP_evolution &t_step, double fs,
                          std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    return calc_Qdt_bisection(std::forward<FunctionTable>(f), m, t_step(), fs,
                              order);
  }

  Maybe_error<bool> test_conductance_mean(const Matrix<double> gmean,
                                          const Matrix<double> g) {
    auto max_g = var::max(g);
    auto min_g = var::min(g);
    Maybe_error<bool> out = true;
    return reduce(
        [max_g, min_g](Maybe_error<bool> succeeds,
                       auto e) -> Maybe_error<bool> {
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

  Maybe_error<bool> test_conductance_variance(const Matrix<double> gvar,
                                              const Matrix<double> g) {
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
  Maybe_error<bool>
  test_conductance_variance(const C_Matrix &var,
                            Conductance_variance_error_tolerance tol) {
    if (var.ncols() == var.nrows()) {
      for (std::size_t i = 0; i < var.nrows(); ++i)
        for (std::size_t j = 0; j < var.ncols(); ++j)
          if (primitive(var(i, j)) + tol() < 0) {
            std::stringstream ss;
            ss << " negative diagonal variance at i=" << i << ", j= " << j
               << "\n"
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
  Maybe_error<bool>
  test_conductance_variances(const C_Vector_Space &q,
                             Conductance_variance_error_tolerance tol) {
    return ((typeid(Variances).name() >>
             test_conductance_variance(get<Variances>(q)(), tol)) &&
            ...);
  }

  template <class C_Qdt>
    requires U<C_Qdt, Qdt>
  Maybe_error<bool> test(const C_Qdt &q,
                         Conductance_variance_error_tolerance tol) {
    return "fails Qdt test\n" >>
           test_conductance_variances<gmean_i, gtotal_ij, gmean_ij,
                                      gtotal_sqr_ij, gsqr_i, gvar_i,
                                      gtotal_var_ij, gvar_ij>(q, tol);
  }

  template <class C_Matrix>
    requires U<C_Matrix, Matrix<double>>

  C_Matrix truncate_negative_variance(C_Matrix &&var) {
    for (std::size_t i = 0; i < var.nrows(); ++i)
      for (std::size_t j = 0; j < var.ncols(); ++j)
        var.set(i, j, max(0.0, var(i, j)));
    return var;
  }

  /*
  template<uses_recursive_aproximation recursive,uses_averaging_aproximation
  averaging, uses_variance_aproximation variance> auto run_old(const
  Patch_State &t_prior, Qdt const &t_Qdt, Patch_Model const &m, const
  Experiment_step &p, double fs) const { auto &p_y = get<Patch_current>(p);
  auto &p_P_mean = get<P_mean>(t_prior); auto &p_P_Cov = get<P_Cov>(t_prior);

  double e =
      get<Current_Noise>(m).value() * fs/ get<number_of_samples>(p).value();
  double N = get<N_Ch_mean>(m)();
  auto N_states = p_P_mean().nrows();
  Matrix<double> u(N_states, 1, 1.0);

  auto SmD = p_P_Cov() - diag(p_P_mean());
  double gSg = xtAx(get<gmean_i>(t_Qdt)(), SmD) +
               getvalue(p_P_mean() *
                        zip([](auto x, auto y) { return x * y; },
                            get<gtotal_ij>(t_Qdt)(), get<gmean_ij>(t_Qdt)()) *
                        u);

  double ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

  auto e_mu = e + N * ms;
  auto v_y_mean = y_mean(N * getvalue(p_P_mean() * get<gmean_i>(t_Qdt)()));
  auto v_y_var = y_var(e_mu + N * gSg);
  if (std::isnan(p_y.value())) {
    auto v_vplogL = vplogL(0.0);
    auto v_plogL = plogL(std::numeric_limits<double>::quiet_NaN());
    auto v_eplogL = eplogL(std::numeric_limits<double>::quiet_NaN());
    auto v_P_cov = P_Cov(AT_B_A(get<P>(t_Qdt)(), SmD));
    auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)());
    v_P_cov() = v_P_cov() + diag(v_P_mean());

    return Patch_State(logL(get<logL>(t_prior)()),v_P_mean, v_P_cov, v_y_mean,
  v_y_var, v_plogL, v_eplogL, v_vplogL);
    // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
    // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
    //      auto test = mp_state_information::test(P_mean, P__cov,
    //      tolerance_); if (test.has_value())
    //        return Op(mp_state_information::adjust(std::move(P_mean),
    //                                               std::move(P__cov),
    //                                               y_mean, y_var, plogL,
    //                                               eplogL,
    //                                               vplogL,Q_dt.min_P(), e));
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
    v_plogL() = -0.5 * log(2 * std::numbers::pi * v_y_var()) - 0.5 * chi2;
  else
    v_plogL() = std::numeric_limits<double>::infinity();

  auto v_eplogL = eplogL(-0.5 * log(2 * std::numbers::pi * v_y_var()) -
                         0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
  vplogL v_vplogL(0.5);
  // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

  //    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var,
  //    plogL,
  //                                           eplogL, e, tolerance());
  //    if (!test) {
  //      std::stringstream ss;

  //      ss << "\nP_mean \n" << P_mean;
  //      ss << "\nPcov \n" << P__cov;
  //      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

  //      return Op(false, "\nfails in trace!!!; error=" + test.error()() +
  //      ss.str());
  //    } else
  return Patch_State(logL(get<logL>(t_prior)()+v_plogL()),v_P_mean, v_P_cov,
  v_y_mean, v_y_var, v_plogL, v_eplogL, v_vplogL);
  }

  */

  //    Maybe_error<Patch_State> DVR(const Patch_State &t_prior, Qdt const
  //    &t_Qdt,
  //                                 Patch_Model const &m, const ATP_step &p,
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

  template <uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance, class FunctionTable,
            class C_Patch_State, class C_Qdt, class C_Patch_Model,
            class C_double>
    requires(
        /*(U<std::decay_t<C_Patch_State>,
           Patch_State>||U<std::decay_t<C_Patch_State>,
           Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
        U<C_double, double> && U<C_Qdt, Qdt>)

  Maybe_error<C_Patch_State>
  Macror_old(FunctionTable &, C_Patch_State &&t_prior, C_Qdt const &t_Qdt,
             C_Patch_Model const &m, C_double const &Nch,
             const Patch_current &p_y, double fs) const {

    using Transf = transformation_type_t<C_Qdt>;
    auto &p_P_cov = get<P_Cov>(t_prior);
    auto &p_P_mean = get<P_mean>(t_prior);
    //    auto &y = get<Patch_current>(p).value();
    auto &y = p_y.value();

    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto e = get<Current_Noise>(m).value() * fs /
                 get<number_of_samples>(t_Qdt).value() +
             get<Pink_Noise>(m).value() +
             get<Proportional_Noise>(m).value() * std::abs(y);
    ;
    auto y_baseline = get<Current_Baseline>(m);
    auto N = Nch;
    Matrix<double> u(p_P_mean().size(), 1, 1.0);

    auto SmD = p_P_cov() - diag(p_P_mean());

    auto N_states = p_P_mean().ncols();

    auto &t_gmean_i = get<gmean_i>(t_Qdt);
    auto &t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    auto &t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
    auto &t_gmean_ij = get<gmean_ij>(t_Qdt);
    auto &t_gvar_i = get<gvar_i>(t_Qdt);
    auto gSg =
        getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
        getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    if constexpr (false) {
      auto test_gSg = var::test_Derivative(
          [this, &u](auto const &t_gmean_i, auto const &t_gmean_ij,
                     auto const &SmD, auto const &p_P_mean,
                     auto const &t_gtotal_ij) {
            return getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                   getvalue(p_P_mean() *
                            (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
          },
          1e-4, 1e-6, t_gmean_i, t_gmean_ij, SmD, p_P_mean, t_gtotal_ij);
      if (!test_gSg) {
        std::cerr << "\n error in test_gSg!!\n" << test_gSg.error()();
        return Maybe_error<C_Patch_State>(test_gSg.error());
      }
    }
    auto ms = getvalue(p_P_mean() * t_gvar_i());

    Op_t<Transf, double> e_mu;
    Op_t<Transf, y_mean> r_y_mean;
    Op_t<Transf, y_var> r_y_var;

    Op_t<Transf, double> sSg;
    Op_t<Transf, double> sSs;
    Op_t<Transf, double> zeta;
    auto t_P = get<P>(t_Qdt);

    if constexpr ((!variance.value) && (!recursive.value)) {
      e_mu = e + max(0.0, N * ms);
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline();
      r_y_var() = e + max(0.0, N * ms + N * gSg);
      if (!(primitive(r_y_var()) > 0.0)) {
        std::stringstream ss;
        ss << "Negative variance!!\n";
        ss << "\nr_y_var=\t" << r_y_var;
        ss << "\ne_mu=\t" << e_mu;
        ss << "\ne=\t" << e;
        ss << "\nN=\t" << N;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\ngSg=\t" << gSg;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\n"
           << "p_P_mean()"
           << "=\t" << p_P_mean();
        ss << "\n"
           << "t_gvar_i()"
           << "=\t" << t_gvar_i();

        return error_message(ss.str());
      }

    } else if constexpr (!variance.value && recursive.value) {
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();

      auto ms = getvalue(p_P_mean() * t_gvar_i());

      e_mu = e + max(N * ms, 0.0);
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline();
      r_y_var() = e + max(0.0, N * ms + N * gSg);
      if (!(primitive(r_y_var()) > 0)) {
        std::stringstream ss;
        ss << "Negative variance!!\n";
        ss << "\nr_y_var=\t" << r_y_var;
        ss << "\ne_mu=\t" << e_mu;
        ss << "\ne=\t" << e;
        ss << "\nN=\t" << N;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\ngSg=\t" << gSg;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\n"
           << "p_P_mean()"
           << "=\t" << p_P_mean();
        ss << "\n"
           << "t_gvar_i()"
           << "=\t" << t_gvar_i();

        return error_message(ss.str());
      }

    } else // (variance && (recursive || !recursive))
    {
      auto &t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
      auto &t_gvar_ij = get<gvar_ij>(t_Qdt);

      sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
            getvalue(p_P_mean() *
                     (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
      sSs =
          getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
          getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));

      auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
      auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

      e_mu = e + N * ms0;
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) -
                   N * 0.5 / e_mu * sSg + y_baseline();
      zeta = N / (2 * sqr(e_mu) + N * sSs);
      r_y_var() = var::max(e, e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
      if (!(primitive(r_y_var()) > 0)) {
        std::stringstream ss;
        ss << "Negative variance!!\n";
        ss << "\nr_y_var=\t" << r_y_var;
        ss << "\ne_mu=\t" << e_mu;
        ss << "\ne=\t" << e;
        ss << "\nN=\t" << N;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\ngSg=\t" << gSg;
        ss << "\n"
           << "ms"
           << "=\t" << ms;
        ss << "\n"
           << "p_P_mean()"
           << "=\t" << p_P_mean();
        ss << "\n"
           << "t_gvar_i()"
           << "=\t" << t_gvar_i();

        return error_message(ss.str());
      }
    }
    if (std::isnan(y)) {

      auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
      auto r_P_mean = build<P_mean>(to_Probability(p_P_mean() * t_P()));
      r_P_cov() = r_P_cov() + diag(r_P_mean());
      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto r_test = test<true>(r_P_mean, r_P_cov, t_tolerance);
      if constexpr (true || r_test)
        if constexpr (U<C_Patch_State, Patch_State>)
          return Op_t<Transf, Patch_State>(
              Op_t<Transf, logL>(get<logL>(t_prior)()),
              Op_t<Transf, elogL>(get<elogL>(t_prior)()),
              Op_t<Transf, vlogL>(get<vlogL>(t_prior)()),
              build<P_mean>(to_Probability(std::move(r_P_mean))),
              normalize(std::move(r_P_cov), t_min_P()), std::move(r_y_mean),
              std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN));
        else {
          auto &ev = get<Patch_State_Evolution>(t_prior);
          r_P_mean = normalize(std::move(r_P_mean), t_min_P());
          r_P_cov = normalize(std::move(r_P_cov), t_min_P());
          ev().push_back(Op_t<Transf, Patch_State>(
              Op_t<Transf, logL>(get<logL>(t_prior)()),
              Op_t<Transf, elogL>(get<elogL>(t_prior)()),
              Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), r_P_mean, r_P_cov,
              r_y_mean, r_y_var, plogL(NaN), eplogL(NaN), vplogL(NaN)));
          return Op_t<Transf, Patch_State_and_Evolution>(
              Op_t<Transf, logL>(get<logL>(t_prior)()),
              Op_t<Transf, elogL>(get<elogL>(t_prior)()),
              Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean),
              std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var),
              plogL(NaN), eplogL(NaN), vplogL(NaN), std::move(ev));
        }
      else
        return error_message("fails at intertrace prediction!!: " +
                             r_test.error()());
    }

    auto dy = y - r_y_mean();
    auto chi = dy / r_y_var();
    Op_t<Transf, P_mean> r_P_mean;
    Op_t<Transf, P_Cov> r_P_cov;
    if constexpr (!recursive.value) {
      r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));

      r_P_mean = build<P_mean>(p_P_mean() * t_P());
      r_P_cov() = r_P_cov() + diag(r_P_mean());
    } else if constexpr (!variance.value) {
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
      auto gseg = chi * gS;

      r_P_mean() = p_P_mean() * t_P() + chi * gS;

      r_P_cov() = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) -
                  (N / r_y_var()) * XTX(gS);

    } else {
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
      auto sS =
          TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();
      r_P_mean() =
          p_P_mean() * t_P() + chi * gS - (chi * zeta * sSg + 0.5 / e_mu) * sS;

      r_P_cov() =
          AT_B_A(t_P(), SmD) + diag(r_P_mean() * t_P()) -
          (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
          (2.0 * N / r_y_var() * zeta * sSg) * X_plus_XT(TranspMult(sS, gS)) -
          (N / r_y_var()) * XTX(gS);
    }

    auto chi2 = dy * chi;

    Op_t<Transf, plogL> r_plogL;
    if (primitive(r_y_var()) > 0.0)
      r_plogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2;
    else {
      std::stringstream ss;
      ss << "Negative variance!!\n";
      ss << "\nr_y_var=\t" << r_y_var;
      ss << "\ne_mu=\t" << e_mu;
      ss << "\ngSg=\t" << gSg;
      return error_message(ss.str());
    }

    if constexpr (false) {
      auto test_plogL = var::test_Derivative(
          [](auto const &r_y_var, auto const &chi2, auto const &t_prior) {
            return -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2 +
                   get<logL>(t_prior)();
          },
          1e-6, 1e-8, r_y_var, chi2, t_prior);
      if (!test_plogL) {
        std::cerr << "\n error in test_plogL!!\n" << test_plogL.error()();
        return Maybe_error<C_Patch_State>(test_plogL.error());
      }
    }

    //    std::cerr<<p<<"\n";
    //    std::cerr<<"r_plogL\n"<<r_plogL<<"\n";
    //    std::cerr<<"r_P_mean\n"<<r_P_mean<<"\n";
    //    std::cerr<<"r_y_mean\n"<<r_y_mean<<"\n";
    //    std::cerr<<"r_y_var\n"<<r_y_var<<"\n";

    //    if (get<Time>(p)()>1)
    //      std::abort();

    Op_t<Transf, eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
                                  0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    vplogL r_vlogL(0.5);
    auto r_test = test<true>(r_P_mean, r_P_cov, t_tolerance);
    if constexpr (false) {
      if (!r_test) {
        std::stringstream ss;

        ss << "\nP_mean \n" << r_P_mean;
        ss << "\nPcov \n" << r_P_cov;
        // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

        return error_message("\nfails in trace!!!; error=" + r_test.error()() +
                             ss.str());
      }
    } else if (std::isnan(primitive(r_plogL())))
      return error_message("likelihood is nan");
    else if constexpr (U<C_Patch_State, Patch_State>)
      return build<Patch_State>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
          normalize(std::move(r_P_mean), t_min_P()),
          normalize(std::move(r_P_cov), t_min_P()), std::move(r_y_mean),
          std::move(r_y_var), r_plogL, r_eplogL, r_vlogL);
    else {
      auto &ev = get<Patch_State_Evolution>(t_prior);
      r_P_mean = normalize(std::move(r_P_mean), t_min_P());
      r_P_cov = normalize(std::move(r_P_cov), t_min_P());
      ev().push_back(build<Patch_State>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), r_P_mean, r_P_cov,
          r_y_mean, r_y_var, r_plogL, r_eplogL, r_vlogL));
      return build<Patch_State_and_Evolution>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
          std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL,
          r_eplogL, r_vlogL, std::move(ev));
    }
  }

  template <uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
            uses_variance_correction_aproximation variance_correction,
            class FunctionTable, class C_Patch_State, class C_Qdt,
            class C_Patch_Model, class C_double>
    requires(
        /*(U<std::decay_t<C_Patch_State>,
         Patch_State>||U<std::decay_t<C_Patch_State>,
         Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
        U<C_double, double> && (U<C_Qdt, Qdt>||U<C_Qdt, Qdtm>))

  Maybe_error<C_Patch_State> Macror(FunctionTable &, C_Patch_State &&t_prior,
                                    C_Qdt const &t_Qdt, C_Patch_Model const &m,
                                    C_double const &Nch,
                                    const Patch_current &p_y, double fs) const {

    get<macror_algorithm>(t_prior)() =
        ToString(MacroR2<::V<recursive>, ::V<averaging>, ::V<variance>,
                         ::V<variance_correction>>{});
    using Transf = transformation_type_t<C_Qdt>;

    auto &p_P_mean = get<P_mean>(t_prior);
    auto SmD = get<P_Cov>(t_prior)() - diag(p_P_mean());
    auto &y = p_y.value();

    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto y_baseline = get<Current_Baseline>(m);
    auto e = get<Current_Noise>(m).value() * fs /
                 get<number_of_samples>(t_Qdt).value() +
             get<Pink_Noise>(m).value();

    Op_t<Transf, double> ms = 0;
    if constexpr (variance.value)
      ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

    auto N = Nch;
    Matrix<double> u(p_P_mean().size(), 1, 1.0);

    auto N_states = p_P_mean().ncols();

    auto &t_gmean_i = get<gmean_i>(t_Qdt);
    auto &t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    auto &t_gmean_ij = get<gmean_ij>(t_Qdt);
    auto gSg =
        getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
        getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

    Op_t<Transf, y_mean> r_y_mean;
    Op_t<Transf, y_var> r_y_var;

    Op_t<Transf, double> sSg;
    auto t_P = get<P>(t_Qdt);

    r_y_mean =
        build<y_mean>(N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline());

    if (std::isnan(y)) {
      get<macror_algorithm>(t_prior)() = ToString(
          MacroR2<::V<uses_recursive_aproximation(false)>, ::V<averaging>,
                  ::V<variance>, ::V<variance_correction>>{});

      auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
      auto r_P_mean = build<P_mean>(to_Probability(p_P_mean() * t_P()));
      r_P_cov() = r_P_cov() + diag(r_P_mean());
      if constexpr (U<C_Patch_State, Patch_State>)
        return Op_t<Transf, Patch_State>(
            Op_t<Transf, logL>(get<logL>(t_prior)()),
            Op_t<Transf, elogL>(get<elogL>(t_prior)()),
            Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean),
            std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var),
            plogL(NaN), eplogL(NaN), vplogL(NaN),
            get<macror_algorithm>(t_prior));
      else {
        auto &ev = get<Patch_State_Evolution>(t_prior);
        ev().push_back(Op_t<Transf, Patch_State>(
            Op_t<Transf, logL>(get<logL>(t_prior)()),
            Op_t<Transf, elogL>(get<elogL>(t_prior)()),
            Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), r_P_mean, r_P_cov,
            r_y_mean, r_y_var, plogL(NaN), eplogL(NaN), vplogL(NaN),
            get<macror_algorithm>(t_prior)));
        return Op_t<Transf, Patch_State_and_Evolution>(
            Op_t<Transf, logL>(get<logL>(t_prior)()),
            Op_t<Transf, elogL>(get<elogL>(t_prior)()),
            Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean),
            std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var),
            plogL(NaN), eplogL(NaN), vplogL(NaN),
            get<macror_algorithm>(t_prior), std::move(ev));
      }
    }

    constexpr bool PoissonDif = true;

    if constexpr (PoissonDif)
      e = e + get<Proportional_Noise>(m).value() * std::abs(y - r_y_mean());
    else
      e = e + get<Proportional_Noise>(m).value() * std::abs(y);

    auto r_y_mean_max = max_possible_value_of_ymean(
        N_Ch_mean_value(primitive(Nch)), primitive(get<g>(m)),
        primitive(y_baseline));

    auto r_y_mean_min = min_possible_value_of_ymean(
        N_Ch_mean_value(primitive(Nch)), primitive(get<g>(m)),
        primitive(y_baseline));

    if ((primitive(r_y_mean()) - r_y_mean_max()) >
        std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_max())) *
            1e-3)
      std::cerr << "\n max violation" << r_y_mean()
                << "  vs  max: " << r_y_mean_max();
    // if ((r_y_mean_min() - primitive(r_y_mean())) >
    //     std::max(std::abs(primitive(r_y_mean())), std::abs(r_y_mean_min())) *
    //         1e-1)
    // std::cerr << "\n min violation\n"
    //          << r_y_mean() << "  vs  min: " << r_y_mean_min();

    if (primitive(gSg) > 0) {
      if (primitive(ms) > 0) {
        r_y_var = build<y_var>(e + N * gSg + N * ms);
      } else {
        r_y_var = build<y_var>(e + N * gSg);
      }
    } else {
      if (primitive(ms) > 0)
        r_y_var = build<y_var>(e + N * ms);
      else
        r_y_var = build<y_var>(e);
    }

    auto dy = y - r_y_mean();
    auto chi = dy / r_y_var();
    Op_t<Transf, P_mean> r_P_mean;
    Op_t<Transf, P_Cov> r_P_cov;

    if constexpr (!recursive.value) {
      r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
      r_P_mean = build<P_mean>(to_Probability(p_P_mean() * t_P()));
      r_P_cov() = r_P_cov() + diag(r_P_mean());

    } else if constexpr (!variance_correction.value) {
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();

      r_P_mean() = p_P_mean() * t_P() + chi * gS;

      r_P_cov() = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) -
                  (N / r_y_var()) * XTX(gS);
    } else {
      auto &t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
      auto &t_gvar_ij = get<gvar_ij>(t_Qdt);
      auto &t_gtotal_ij = get<gtotal_ij>(t_Qdt);
      auto &t_gvar_i = get<gvar_i>(t_Qdt);
      auto gSg =
          getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
          getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));

      auto sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
                 getvalue(p_P_mean() *
                          (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
      auto sSs =
          getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
          getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));

      auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
      auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

      auto e_mu = e + N * ms0;
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) -
                   N * 0.5 / e_mu * sSg + y_baseline();
      auto zeta = N / (2 * sqr(e_mu) + N * sSs);
      r_y_var() = var::max(e, e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
      auto sS =
          TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();
      r_P_mean() =
          p_P_mean() * t_P() + chi * gS - (chi * zeta * sSg + 0.5 / e_mu) * sS;

      r_P_cov() =
          AT_B_A(t_P(), SmD) + diag(r_P_mean() * t_P()) -
          (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
          (2.0 * N / r_y_var() * zeta * sSg) * X_plus_XT(TranspMult(sS, gS)) -
          (N / r_y_var()) * XTX(gS);
    }

    if (!all_Probability_elements(primitive(r_P_mean())) ||
        !all_Covariance_elements(primitive(r_P_cov()))) {
      r_P_mean() = p_P_mean() * t_P();

      r_P_cov() = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P());

      get<macror_algorithm>(t_prior)() = ToString(
          MacroR2<::V<uses_recursive_aproximation(false)>, ::V<averaging>,
                  ::V<variance>, ::V<variance_correction>>{});
    }

    auto chi2 = dy * chi;

    Op_t<Transf, plogL> r_plogL;
    Op_t<Transf, eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
                                  0.5);
    if (primitive(r_y_var()) > 0.0) {
      if (get<Proportional_Noise>(m).value() == 0) {
        r_plogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2;
        r_eplogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5;
      } else {
        r_plogL() = -log(Poisson_noise_normalization(
                        r_y_var(), get<Proportional_Noise>(m).value())) -
                    0.5 * chi2;
        r_eplogL() = Poisson_noise_expected_logL(
            r_y_var(), get<Proportional_Noise>(m).value());
      }
    } else {
      std::stringstream ss;
      ss << "Negative variance!!\n";
      ss << "\nr_y_var=\t" << r_y_var;
      ss << "\ngSg=\t" << gSg;
      return error_message(ss.str());
    }

    vplogL r_vlogL(0.5);
    if (std::isnan(primitive(r_plogL()))) {
      std::stringstream ss;
      ss << "likelihood is nan \n patch current=";
      print(ss, p_y) << "Qdt";
      print(ss, primitive(t_Qdt)) << "tprior";
      print(ss, primitive(t_prior));
      return error_message(ss.str());
    } else if constexpr (U<C_Patch_State, Patch_State>)
      return build<Patch_State>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
          std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL,
          r_eplogL, r_vlogL, get<macror_algorithm>(t_prior));
    else {
      auto &ev = get<Patch_State_Evolution>(t_prior);
      ev().push_back(build<Patch_State>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), r_P_mean, r_P_cov,
          r_y_mean, r_y_var, r_plogL, r_eplogL, r_vlogL,
          get<macror_algorithm>(t_prior)));
      return build<Patch_State_and_Evolution>(
          build<logL>(get<logL>(t_prior)() + r_plogL()),
          build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
          build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
          std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL,
          r_eplogL, r_vlogL, get<macror_algorithm>(t_prior), std::move(ev));
    }
  }

  template <return_predictions predictions, class C_Patch_Model>
    requires U<C_Patch_Model, Patch_Model>
  auto init(const C_Patch_Model &m, initial_ATP_concentration initial_x)
      -> Maybe_error<Transfer_Op_to<
          C_Patch_Model,
          std::conditional_t<predictions.value, Patch_State_and_Evolution,
                             Patch_State>>> {
    auto v_Qx = calc_Qx(m, initial_x());
    auto r_P_mean = build<P_mean>(get<P_initial>(m)());
    auto r_P_cov = build<P_Cov>(diagpos(r_P_mean()) - XTX(r_P_mean.value()));
    auto r_test =
        test<true>(r_P_mean, r_P_cov, get<Probability_error_tolerance>(m));
    if (r_test) {
      auto t_min_P = get<min_P>(m);
      if (false) {
        std::cerr << "initial\n";
        std::cerr << "r_P_mean" << r_P_mean;
        std::cerr << "r_P_cov" << r_P_cov;
        //  std::cerr<<"normalized r_P_cov"<<normalize(std::move(r_P_cov),
        //  t_min_P());
      }
      if constexpr (!predictions.value)

        return Transfer_Op_to<C_Patch_Model, Patch_State>(
            logL(0.0), elogL(0.0), vlogL(0.0), std::move(r_P_mean),
            std::move(r_P_cov), y_mean(NaN), y_var(NaN), plogL(NaN),
            eplogL(NaN), vplogL(NaN), macror_algorithm(""));
      else
        return Transfer_Op_to<C_Patch_Model, Patch_State_and_Evolution>(
            logL(0.0), elogL(0.0), vlogL(0.0), std::move(r_P_mean),
            std::move(r_P_cov), y_mean(NaN), y_var(NaN), plogL(NaN),
            eplogL(NaN), vplogL(NaN), macror_algorithm(""),
            Patch_State_Evolution());

    } else {
      return error_message("fails at init: " + r_test.error()());
    }
  }

  template <uses_adaptive_aproximation adaptive,
            uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
            uses_variance_correction_aproximation variance_correction,
            return_predictions predictions, class FuncTable, class C_Parameters,
            class Model>
    requires(!is_of_this_template_type_v<FuncTable, FuncMap_St>)
  auto log_Likelihood(FuncTable &&f, const Model &model,
                      const C_Parameters &par, const Experiment &e,
                      const Recording &y)
      -> Maybe_error<Transfer_Op_to<
          C_Parameters,
          std::conditional_t<predictions.value, Patch_State_Evolution,
                             Vector_Space<logL, elogL, vlogL>>>> {

    f.clear();

    using Transf = transformation_type_t<C_Parameters>;
    using C_Patch_State =
        Op_t<Transf,
             std::conditional_t<predictions.value, Patch_State_and_Evolution,
                                Patch_State>>;
    auto Maybe_m = model(par);
    if (!is_valid(Maybe_m))
      return get_error(Maybe_m);
    else {
      auto m = std::move(get_value(Maybe_m));
      auto fs = get<Frequency_of_Sampling>(e).value();
      auto ini = init<predictions>(m, get<initial_ATP_concentration>(e));

      auto gege = 0;
      if (!ini)
        return ini.error();
      else {
        auto run = fold(
            0ul, y().size(), std::move(ini).value(),
            [this, &f, &m, fs, &e, &y, &gege](C_Patch_State &&t_prior,
                                              std::size_t i_step) {
              ATP_evolution const &t_step =
                  get<ATP_evolution>(get<Recording_conditions>(e)()[i_step]);

              auto time = get<Time>(get<Recording_conditions>(e)()[i_step])();
              auto time_segment = get<N_Ch_mean_time_segment_duration>(m)();
              auto Nchs = get<N_Ch_mean>(m)();
              std::size_t i_segment =
                  std::min(Nchs.size() - 1.0, std::floor(time / time_segment));
              auto j_segment = std::min(Nchs.size() - 1, i_segment + 1);
              auto r = std::max(1.0, time / time_segment - i_segment);
              auto Nch = Nchs[i_segment] * (1 - r) + r * Nchs[j_segment];

              auto Maybe_t_Qdt = calc_Qdt(f, m, t_step, fs);
              if (!Maybe_t_Qdt)
                return Maybe_error<C_Patch_State>(Maybe_t_Qdt.error());
              else {
                auto t_Qdt = std::move(Maybe_t_Qdt.value());

                //
                if constexpr (false) {
                  auto test_der_t_Qdt = var::test_Derivative(
                      [this, &t_step, &fs, &gege, &f](auto const &l_m,
                                                      auto const &l_Qx) {
                        return calc_Qdt(f, l_m, t_step, fs);
                      },
                      1e-6, 1e-2, m, t_Qdt);
                  if (true && !test_der_t_Qdt) {
                    std::cerr << test_der_t_Qdt.error()();
                    std::cerr << "\nt_step\n" << t_step;
                    return Maybe_error<C_Patch_State>(test_der_t_Qdt.error());
                  }
                }
                if constexpr (false) {
                  if (gege < 10)
                    ++gege;
                  else
                    abort();
                }
                if (false) {
                  auto test_Qdt =
                      test(t_Qdt, get<Conductance_variance_error_tolerance>(m));

                  if (!test_Qdt)
                    return Maybe_error<C_Patch_State>(test_Qdt.error());
                }

                if constexpr (false) {
                  auto test_der_Macror = var::test_Derivative(
                      [this, &t_step, &fs, &f](auto const &l_m,
                                               auto const &l_prior,
                                               auto const &l_Qdt) {
                        return f.ff(
                            MacroR<uses_recursive_aproximation(false),
                                   uses_averaging_aproximation(1),
                                   uses_variance_aproximation(false)>{})(
                            l_prior, l_Qdt, l_m, fs);
                      },
                      1e-7, 1e-14, m, t_prior, t_Qdt);
                  if (!test_der_Macror) {
                    std::cerr << "\nt_step\n" << t_step;
                    std::cerr << test_der_Macror.error()();
                    std::cerr << "\nt_step\n" << t_step;
                    //   return
                    //   Maybe_error<C_Patch_State>(test_der_Macror.error());
                  }
                }

                if constexpr (false) {

                  std::cerr << "\nplogL\n" << get<plogL>(t_prior);
                  std::cerr << "\nlogL\n" << get<logL>(t_prior);
                  std::cerr << "\nt_step\n" << t_step << "\n";
                }
                if constexpr (!adaptive.value) {
                  return f.f(MacroR<recursive, averaging, variance>{},
                             std::move(t_prior), t_Qdt, m, Nch, y()[i_step],
                             fs);
                } else {
                  double mg = getvalue(primitive(get<P_mean>(t_prior)()) *
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
                    return f.f(MacroR<recursive, averaging, variance>{},
                               std::move(t_prior), t_Qdt, m, Nch, y()[i_step],
                               fs);
                  } else {
                    return f.f(
                        MacroR<uses_recursive_aproximation(false), averaging,
                               uses_variance_aproximation(false)>{},
                        std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                  }
                }
              }
            });
        if (!run)
          return run.error();
        else if constexpr (predictions.value)
          return get<Patch_State_Evolution>(run.value());
        else
          return build<Vector_Space<logL, elogL, vlogL>>(
              get<logL>(run.value()), get<elogL>(run.value()),
              get<vlogL>(run.value()));
      }
    }
  }

  template <uses_adaptive_aproximation adaptive,
            uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
            uses_variance_correction_aproximation variance_correction,
            return_predictions predictions, class FuncTable, class C_Parameters,
            class Model>
    requires(is_of_this_template_type_v<FuncTable, FuncMap_St>)
  auto log_Likelihood(FuncTable &f, const Model &model, const C_Parameters &par,
                      const Experiment &e, const Recording &y)
      -> Maybe_error<Transfer_Op_to<
          C_Parameters,
          std::conditional_t<predictions.value, Patch_State_Evolution,
                             Vector_Space<logL, elogL, vlogL>>>> {

    using v_adaptive = ::V<adaptive>;
    using v_recursive = ::V<recursive>;
    using v_averaging = ::V<averaging>;
    using v_variance = ::V<variance>;
    using v_variance_correction = ::V<variance_correction>;

    using Transf = transformation_type_t<C_Parameters>;
    using C_Patch_State =
        Op_t<Transf,
             std::conditional_t<predictions.value, Patch_State_and_Evolution,
                                Patch_State>>;

    auto Maybe_m = model(par);
    if (!is_valid(Maybe_m))
      return get_error(Maybe_m);
    else {
      auto m = std::move(get_value(Maybe_m));
      auto fs = get<Frequency_of_Sampling>(e).value();
      auto ini = init<predictions>(m, get<initial_ATP_concentration>(e));

      auto gege = 0;
      auto f_local = f.create("_lik");
      if (!ini)
        return ini.error();
      else {
        auto run = fold(
            0ul, y().size(), std::move(ini).value(),
            [this, &f_local, &m, fs, &e, &y, &gege](C_Patch_State &&t_prior,
                                                    std::size_t i_step) {
              ATP_evolution const &t_step =
                  get<ATP_evolution>(get<Recording_conditions>(e)()[i_step]);

              auto time = get<Time>(get<Recording_conditions>(e)()[i_step])();
              auto time_segment = get<N_Ch_mean_time_segment_duration>(m)();
              auto Nchs = get<N_Ch_mean>(m)();
              std::size_t i_segment =
                  std::min(Nchs.size() - 1.0, std::floor(time / time_segment));
              auto j_segment = std::min(Nchs.size() - 1, i_segment + 1);
              auto r = std::max(1.0, time / time_segment - i_segment);
              auto Nch = Nchs[i_segment] * (1 - r) + r * Nchs[j_segment];

              
               
                if constexpr (!adaptive.value) {
                  if constexpr(!variance_correction.value){
                    auto Maybe_t_Qdtm = calc_Qdtm(f_local, m, t_step, fs);
                    if (!Maybe_t_Qdtm)
                      return Maybe_error<C_Patch_State>(Maybe_t_Qdtm.error());
                    auto t_Qdtm = std::move(Maybe_t_Qdtm.value());
                  return MacroR2<v_recursive, v_averaging, v_variance,
                                 v_variance_correction>{}(
                      f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step],
                      fs);
                    
                  }
               else{
                    auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                    if (!Maybe_t_Qdt)
                      return Maybe_error<C_Patch_State>(Maybe_t_Qdt.error());
                    auto t_Qdt = std::move(Maybe_t_Qdt.value());
                  return MacroR2<v_recursive, v_averaging, v_variance,
                                 v_variance_correction>{}(
                      f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step],
                      fs);
                  }
                } else {if (!variance_correction.value){
                      auto Maybe_t_Qdtm = calc_Qdtm(f_local, m, t_step, fs);
                      if (!Maybe_t_Qdtm)
                        return Maybe_error<C_Patch_State>(Maybe_t_Qdtm.error());
                      auto t_Qdtm = std::move(Maybe_t_Qdtm.value());

                  double mg = getvalue(primitive(get<P_mean>(t_prior)()) *
                                       primitive(get<gmean_i>(t_Qdtm)()));
                  double g_max = var::max(get<gmean_i>(primitive(t_Qdtm))());
                  double g_min = var::min(get<gmean_i>(primitive(t_Qdtm))());
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
                    return MacroR2<v_recursive, v_averaging, v_variance,
                                   v_variance_correction>{}(
                        f_local, std::move(t_prior), t_Qdtm, m, Nch, y()[i_step],
                        fs);
                  } else {
                    return
                        //   f_local.f(
                        // MacroR2<::V<uses_recursive_aproximation(false)>,
                        //         v_averaging,
                        //         ::V<uses_variance_aproximation(false)>,
                        //         ::V<uses_variance_correction_aproximation(
                        //             false)>>{},
                        // std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);

                        MacroR2<::V<uses_recursive_aproximation(false)>,
                                v_averaging,
                                ::V<uses_variance_aproximation(false)>,
                                ::V<uses_variance_correction_aproximation(
                                    false)>>{}(f_local, std::move(t_prior),
                                               t_Qdtm, m, Nch, y()[i_step], fs);
                  }
                }
                    else  {
                      auto Maybe_t_Qdt = calc_Qdt(f_local, m, t_step, fs);
                      if (!Maybe_t_Qdt)
                        return Maybe_error<C_Patch_State>(Maybe_t_Qdt.error());
                      auto t_Qdt = std::move(Maybe_t_Qdt.value());
                      
                      double mg = getvalue(primitive(get<P_mean>(t_prior)()) *
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
                        return MacroR2<v_recursive, v_averaging, v_variance,
                                       v_variance_correction>{}(
                            f_local, std::move(t_prior), t_Qdt, m, Nch, y()[i_step],
                            fs);
                      } else {
                        return MacroR2<::V<uses_recursive_aproximation(false)>,
                                    v_averaging,
                                    ::V<uses_variance_aproximation(false)>,
                                    ::V<uses_variance_correction_aproximation(
                                        false)>>{}(f_local, std::move(t_prior),
                                                   t_Qdt, m, Nch, y()[i_step], fs);
                      }
                    }
                  
            }});
        f += f_local;
        if (!run)
          return run.error();
        else if constexpr (predictions.value)
          return get<Patch_State_Evolution>(run.value());
        else
          return build<Vector_Space<logL, elogL, vlogL>>(
              get<logL>(run.value()), get<elogL>(run.value()),
              get<vlogL>(run.value()));
      }
    }
  }

  template <class Patch_Model>
  Maybe_error<Simulated_Sub_Step>
  sub_sub_sample(mt_64i &mt, Simulated_Sub_Step const &t_sim_step,
                 const Patch_Model &m, const ATP_step &t_s,
                 std::size_t n_sub_dt, double fs) {
    auto &t_g = get<g>(m);
    auto N = get<N_channel_state>(t_sim_step);
    double ysum = get<y_sum>(t_sim_step)();
    auto sum_samples = get<number_of_samples>(t_sim_step)();

    auto n_samples = get<number_of_samples>(t_s)();
    auto tQx = calc_Qx(m, get<ATP_concentration>(t_s));

    auto dt = n_samples / fs;
    auto sub_dt = dt / n_sub_dt;

    double sub_sample = 1.0 * n_samples / n_sub_dt;
    auto Maybe_t_P = calc_P(m, tQx, sub_dt, get<min_P>(m)());
    if (!Maybe_t_P)
      return Maybe_t_P.error();
    else {

      auto t_P = std::move(Maybe_t_P.value());
      for (std::size_t i = 0; i < n_sub_dt; ++i) {
        N = sample_Multinomial(mt, t_P, N);
        ysum += getvalue(N() * t_g()) * sub_sample;
      }
      sum_samples += n_samples;
      //  std::cerr << N << sum_samples << "  " << ysum << "  "
      //            << ysum / sum_samples << "\n";
      return Simulated_Sub_Step(N_channel_state(N),
                                number_of_samples(sum_samples), y_sum(ysum));
    }
  }

  template <class Patch_Model>
  Maybe_error<Simulated_Sub_Step>
  sub_sub_sample(mt_64i &mt, Simulated_Sub_Step &&t_sim_step,
                 const Patch_Model &m, const std::vector<ATP_step> &t_s,
                 std::size_t n_sub_dt, double fs) {

    for (std::size_t i = 0; i < t_s.size(); ++i) {
      auto Maybe_sub_step =
          sub_sub_sample(mt, std::move(t_sim_step), m, t_s[i], n_sub_dt, fs);
      if (!Maybe_sub_step)
        return Maybe_sub_step.error();
      else
        t_sim_step = std::move(Maybe_sub_step.value());
    }
    return t_sim_step;
  }

  template <includes_N_state_evolution keep_N_state, class Patch_Model>
  Maybe_error<Simulated_Step<keep_N_state>>
  sub_sample(mt_64i &mt, Simulated_Step<keep_N_state> &&t_sim_step,
             const Patch_Model &m, const ATP_evolution &t_s,
             std::size_t n_sub_dt, double fs) {
    // auto &N = get<N_channel_state>(t_sim_step);

    // std::cerr<<N();

    auto t_sub_step = Simulated_Sub_Step(get<N_channel_state>(t_sim_step()),
                                         number_of_samples(0ul), y_sum(0.0));

    auto Maybe_t_sub_step =
        sub_sub_sample(mt, std::move(t_sub_step), m, t_s(), n_sub_dt, fs);

    if (!Maybe_t_sub_step)
      return Maybe_t_sub_step.error();
    else {
      t_sub_step = std::move(Maybe_t_sub_step.value());
      double y_mean =
          get<y_sum>(t_sub_step)() / get<number_of_samples>(t_sub_step)();
      get<N_channel_state>(t_sim_step()) = get<N_channel_state>(t_sub_step);

      auto &t_e_step = get<Recording>(
          get<Simulated_Recording<keep_N_state>>(t_sim_step())());
      double e =
          get<Current_Noise>(m)() * fs / get<number_of_samples>(t_sub_step)() +
          get<Pink_Noise>(m).value();
      ;
      auto y_baseline = get<Current_Baseline>(m);

      auto y = y_mean + y_baseline() +
               std::normal_distribution<double>()(mt) * std::sqrt(e);
      auto ey = get<Proportional_Noise>(m).value() * std::abs(y);
      if (ey > 0)
        y = y + std::normal_distribution<double>()(mt) * std::sqrt(ey);
      t_e_step().emplace_back(Patch_current(y));
      if constexpr (keep_N_state.value) {
        get<N_Ch_State_Evolution>(
            get<Simulated_Recording<keep_N_state>>(t_sim_step())())()
            .push_back(get<N_channel_state>(t_sim_step()));
      }

      return t_sim_step;
    }
  }

  template <includes_N_state_evolution keep_N_state, class Patch_Model>
  Simulated_Step<keep_N_state> init_sim(mt_64i &mt, const Patch_Model &m,
                                        const Experiment &e) {
    auto initial_x = get<initial_ATP_concentration>(e);
    auto v_Qx = calc_Qx(m, initial_x());
    auto r_P_mean = P_mean(get<P_initial>(m)());
    auto N = get<N_Ch_mean>(m)()[0];
    auto sim = Simulated_Recording<keep_N_state>{};
    auto N_state = sample_Multinomial(mt, r_P_mean, N);
    return Simulated_Step<keep_N_state>(
        Vector_Space(std::move(N_state), Simulated_Recording<keep_N_state>{}));
  }

  template <includes_N_state_evolution keep_N_state>
  static Simulated_Recording<keep_N_state>
  copy_NaNs(Simulated_Recording<keep_N_state> &&sim, const Recording &r) {
    for (std::size_t i = 0; i < size(r()); ++i)
      if (std::isnan(r()[i]()))
        get<Recording>(sim())()[i]() = r()[i]();

    return std::move(sim);
  }

  template <includes_N_state_evolution keep_N_state, class Model, class Id>
  Maybe_error<Simulated_Recording<keep_N_state>>
  sample_(mt_64i &mt, const Model &model, const var::Parameters_values<Id> &par,
          const Experiment &e, const Simulation_Parameters &sim,
          const Recording &r = Recording{}) {

    auto Maybe_m = model(par);
    if (!Maybe_m)
      return Maybe_m.error();
    else {

      auto m = std::move(Maybe_m.value());

      auto n_sub_dt = get<Simulation_n_sub_dt>(sim);
      auto fs = get<Frequency_of_Sampling>(e).value();
      auto sim_recording = Recording{};

      auto ini = init_sim<keep_N_state>(mt, m, e);
      auto run =
          fold(get<Recording_conditions>(e)(), ini,
               [this, &m, fs, n_sub_dt,
                &mt](Simulated_Step<keep_N_state> &&t_sim_step,
                     Experiment_step const &t_step) {
                 return Maybe_error<Simulated_Step<keep_N_state>>(sub_sample(
                     mt, std::move(t_sim_step), m, t_step, n_sub_dt(), fs));
               });
      if (!run)
        return run.error();
      else {
        return copy_NaNs(
            std::move(get<Simulated_Recording<keep_N_state>>(run.value()())),
            r);
      }
    }
  }
  template <class Model, class Id>
  Maybe_error<Simulated_Recording<includes_N_state_evolution(false)>>
  sample(mt_64i &mt, const Model &model, const var::Parameters_values<Id> &par,
         const Experiment &e, const Simulation_Parameters &sim,
         const Recording &r = Recording{}) {
    return sample_<includes_N_state_evolution(false)>(mt, model, par, e, sim,
                                                      r);
  }

  template <class Model, class Id>
  Maybe_error<Simulated_Recording<includes_N_state_evolution(true)>>
  sample_N(mt_64i &mt, const Model &model,
           const var::Parameters_values<Id> &par, const Experiment &e,
           const Simulation_Parameters &sim, const Recording &r = Recording{}) {
    return sample_<includes_N_state_evolution(true)>(mt, model, par, e, sim, r);
  }
};

template <class recursive, class averaging, class variance,
          class variance_correction>
struct MacroR2 {
  friend std::string ToString(MacroR2) {
    std::string out = "MacroR";
    if (recursive{}.value.value)
      out += "_R";
    else
      out += "_NR";
    if (averaging{}.value.value == 2)
      out += "_2";
    else
      out += "__";
    if (variance{}.value.value)
      out += "_V";
    else
      out += "_M";
    if (variance_correction{}.value.value)
      out += "_V";
    else
      out += "_M";

    return out;
  }

  template <class T, class... Ts> auto operator()(T &&x, Ts &&...xs) {
    auto m = Macro_DMR{};

    return m.Macror<recursive{}.value, averaging{}.value, variance{}.value,
                    variance_correction{}.value>(std::forward<T>(x),
                                                 std::forward<Ts>(xs)...);
  }
};

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction, class Model>
struct Likelihood_Model {
  Model m;
  Simulation_n_sub_dt n_sub_dt;
  Likelihood_Model(const Model &model, Simulation_n_sub_dt n_sub_dt)
      : m{model}, n_sub_dt{n_sub_dt} {}

  template <class Parameter>
  friend void report_model(save_Parameter<Parameter> &s,
                           Likelihood_Model const &d) {
    std::ofstream f(s.fname + "_likelihood_model.csv");
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    f << "adaptive: " << adaptive << "\n";
    f << "recursive: " << recursive << "\n";
    f << "averaging: " << averaging << "\n";
    f << "variance: " << variance << "\n";
    f << "variance_correction: " << variance_correction << "\n";
    f << "Simulation_n_sub_dt: " << d.n_sub_dt << "\n";
    report_model(s, d.m);
  }
};

struct Likelihood_Model_v_all {
  using v_uses_adaptive_aproximation =
      std::variant<::V<uses_adaptive_aproximation(false)>,
                   ::V<uses_adaptive_aproximation(true)>>;

  using v_uses_recursive_aproximation =
      std::variant<::V<uses_recursive_aproximation(false)>,
                   ::V<uses_recursive_aproximation(true)>>;

  using v_uses_averaging_aproximation =
      std::variant<::V<uses_averaging_aproximation(0)>,
                   ::V<uses_averaging_aproximation(1)>,
                   ::V<uses_averaging_aproximation(2)>>;

  using v_uses_variance_aproximation =
      std::variant<::V<uses_variance_aproximation(false)>,
                   ::V<uses_variance_aproximation(true)>>;

  using v_uses_variance_correction_aproximation =
      std::variant<::V<uses_variance_correction_aproximation(false)>,
                   ::V<uses_variance_correction_aproximation(true)>>;

  template <uses_adaptive_aproximation adaptive,
            uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
            uses_variance_correction_aproximation variance_correction,
            class Model>
  auto template_op(::V<adaptive>, ::V<recursive>, ::V<averaging>, ::V<variance>,
                   ::V<variance_correction>, const Model &model,
                   Simulation_n_sub_dt n_sub_dt) const {
    return Likelihood_Model<adaptive, recursive, averaging, variance,
                            variance_correction, Model>(model, n_sub_dt);
  }

  template <class Model>
  auto variant_op(v_uses_adaptive_aproximation t_adaptive,
                  v_uses_recursive_aproximation t_recursive,
                  v_uses_averaging_aproximation t_averaging,
                  v_uses_variance_aproximation t_variance,
                  v_uses_variance_correction_aproximation t_var_corr,
                  Model const &model, Simulation_n_sub_dt n_sub_dt) const {

    auto tu = std::tuple(t_adaptive, t_recursive, t_averaging, t_variance,
                         t_var_corr);
    return Apply_variant(
        [this, &model, n_sub_dt](auto const &...x) {
          // using m11=decltype(template_op(x...,model,
          // n_sub_dt))::llego_aquiÑ;
          //   return m11{};
          return this->template_op(x..., model, n_sub_dt);
        },
        tu);
  }

  template <class Model>
  auto bool_op(uses_adaptive_aproximation adaptive,
               uses_recursive_aproximation recursive,
               uses_averaging_aproximation averaging,
               uses_variance_aproximation variance,
               uses_variance_correction_aproximation variance_correction,
               const Model &model, Simulation_n_sub_dt n_sub_dt) const {
    v_uses_adaptive_aproximation t_adaptive;
    if (adaptive.value)
      t_adaptive = ::V<uses_adaptive_aproximation(true)>{};
    else
      t_adaptive = ::V<uses_adaptive_aproximation(false)>{};

    v_uses_recursive_aproximation t_recursive;
    if (recursive.value)
      t_recursive = ::V<uses_recursive_aproximation(true)>{};
    else
      t_recursive = ::V<uses_recursive_aproximation(false)>{};

    v_uses_averaging_aproximation t_averaging;
    if (averaging.value == 0)
      t_averaging = ::V<uses_averaging_aproximation(0)>{};
    else if (averaging.value == 1)

      t_averaging = ::V<uses_averaging_aproximation(1)>{};
    else
      t_averaging = ::V<uses_averaging_aproximation(2)>{};

    v_uses_variance_aproximation t_variance;
    if (variance.value)
      t_variance = ::V<uses_variance_aproximation(true)>{};
    else
      t_variance = ::V<uses_variance_aproximation(false)>{};

    v_uses_variance_correction_aproximation t_var_corr;
    if (variance_correction.value)
      t_var_corr = ::V<uses_variance_correction_aproximation(true)>{};
    else
      t_var_corr = ::V<uses_variance_correction_aproximation(false)>{};

    return this->variant_op(t_adaptive, t_recursive, t_averaging, t_variance,
                            t_var_corr, model, n_sub_dt);
  }
};

struct Likelihood_Model_v {
  using v_uses_adaptive_aproximation =
      std::variant<::V<uses_adaptive_aproximation(false)>>;

  using v_uses_recursive_aproximation =
      std::variant<::V<uses_recursive_aproximation(false)>,
                   ::V<uses_recursive_aproximation(true)>>;

  using v_uses_averaging_aproximation =
      std::variant<::V<uses_averaging_aproximation(2)>>;

  using v_uses_variance_aproximation =
      std::variant<::V<uses_variance_aproximation(true)>>;

  using v_uses_variance_correction_aproximation =
      std::variant<::V<uses_variance_correction_aproximation(false)>>;

  template <uses_adaptive_aproximation adaptive,
            uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
            uses_variance_correction_aproximation variance_correction,
            class Model>
  auto template_op(::V<adaptive>, ::V<recursive>, ::V<averaging>, ::V<variance>,
                   ::V<variance_correction>, const Model &model,
                   Simulation_n_sub_dt n_sub_dt) const {
    return Likelihood_Model<adaptive, recursive, averaging, variance,
                            variance_correction, Model>(model, n_sub_dt);
  }

  template <class Model>
  auto variant_op(v_uses_adaptive_aproximation t_adaptive,
                  v_uses_recursive_aproximation t_recursive,
                  v_uses_averaging_aproximation t_averaging,
                  v_uses_variance_aproximation t_variance,
                  v_uses_variance_correction_aproximation t_var_corr,
                  Model const &model, Simulation_n_sub_dt n_sub_dt) const {

    auto tu = std::tuple(t_adaptive, t_recursive, t_averaging, t_variance,
                         t_var_corr);
    return Apply_variant(
        [this, &model, n_sub_dt](auto const &...x) {
          // using m11=decltype(template_op(x...,model,
          // n_sub_dt))::llego_aquiÑ;
          //   return m11{};
          return this->template_op(x..., model, n_sub_dt);
        },
        tu);
  }

  template <class Model>
  auto bool_op(uses_adaptive_aproximation adaptive,
               uses_recursive_aproximation recursive,
               uses_averaging_aproximation averaging,
               uses_variance_aproximation variance,
               uses_variance_correction_aproximation variance_correction,
               const Model &model, Simulation_n_sub_dt n_sub_dt) const {
    v_uses_adaptive_aproximation t_adaptive;
    t_adaptive = ::V<uses_adaptive_aproximation(false)>{};

    v_uses_recursive_aproximation t_recursive;
    if (recursive.value)
      t_recursive = ::V<uses_recursive_aproximation(true)>{};
    else
      t_recursive = ::V<uses_recursive_aproximation(false)>{};

    v_uses_averaging_aproximation t_averaging;
    t_averaging = ::V<uses_averaging_aproximation(2)>{};

    v_uses_variance_aproximation t_variance;
    t_variance = ::V<uses_variance_aproximation(true)>{};

    v_uses_variance_correction_aproximation t_var_corr;
    t_var_corr = ::V<uses_variance_correction_aproximation(false)>{};

    return this->variant_op(t_adaptive, t_recursive, t_averaging, t_variance,
                            t_var_corr, model, n_sub_dt);
  }
};

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction, class Model>
auto make_Likelihood_Model(const Model &m, Simulation_n_sub_dt n_sub_dt) {
  return Likelihood_Model<adaptive, recursive, averaging, variance,
                          variance_correction, Model>(m, n_sub_dt);
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction, class FuncTable,
    class Model, class Parameters, class Variables, class DataType>
Maybe_error<logLs>
logLikelihood(FuncTable &&f,
              const Likelihood_Model<adaptive, recursive, averaging, variance,
                                     variance_correction, Model> &lik,
              Parameters const &p, const Variables &var, const DataType &y) {
  auto v_logL =
      Macro_DMR{}
          .log_Likelihood<adaptive, recursive, averaging, variance,
                          variance_correction, return_predictions(false)>(
              std::forward<FuncTable>(f), lik.m, p, y, var);
  if (!v_logL)
    return v_logL.error();
  else
    return logLs(get<logL>(v_logL.value()), get<elogL>(v_logL.value()),
                 get<vlogL>(v_logL.value()));
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction,
    class FunctionTable, class Model, class Parameters, class Variables,
    class DataType>
Maybe_error<Patch_State_Evolution> logLikelihoodPredictions(
    FunctionTable &&f,
    const Likelihood_Model<adaptive, recursive, averaging, variance,
                           variance_correction, Model> &lik,
    Parameters const &p, const Variables &var, const DataType &y) {
  return Macro_DMR{}
      .log_Likelihood<adaptive, recursive, averaging, variance,
                      variance_correction, return_predictions(true)>(
          std::forward<FunctionTable>(f), lik.m, p, y, var);
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction,
    class FunctionTable, class Model, class Parameters, class Variables,
    class DataType>
Maybe_error<std::vector<Patch_State_Evolution>>
fractioned_logLikelihoodPredictions(
    FunctionTable &&f,
    const Likelihood_Model<adaptive, recursive, averaging, variance,
                           variance_correction, Model> &lik,
    Parameters const &p, const std::vector<Variables> &y,
    const std::vector<DataType> &var) {
  std::vector<Patch_State_Evolution> out;
  for (std::size_t i = 0; i < y.size(); ++i) {
    auto Maybe_e =
        Macro_DMR{}
            .log_Likelihood<adaptive, recursive, averaging, variance,
                            variance_correction, return_predictions(true)>(
                std::forward<FunctionTable>(f), lik.m, p, var[i],
                get<Recording>(y[i]()));
    if (Maybe_e)
      out.push_back(Maybe_e.value());
    else
      return Maybe_e.error();
  }
  return out;
}

inline auto get_num_samples(const ATP_step &e) {
  return get<number_of_samples>(e)();
}

inline auto get_num_samples(const std::vector<ATP_step> &e) {
  double out = 0;
  for (auto &elem : e)
    out += get_num_samples(elem);
  return out;
}

inline auto get_num_samples(const ATP_evolution &e) {
  return get_num_samples(e());
}

inline ATP_evolution average_ATP_step(ATP_step const &x) {
  return std::vector<ATP_step>(1, x);
}

inline ATP_evolution average_ATP_step(std::vector<ATP_step> const &x) {
  if (x.empty())
    return x;
  else if (x.size() == 1)
    return x;
  else {
    auto out = std::vector<ATP_step>{x[0]};
    for (std::size_t i = 1; i < x.size(); ++i) {

      if ((get<ATP_concentration>(out.back())() !=
           get<ATP_concentration>(x[i])()) &&
          ((get<ATP_concentration>(out.back())() == 0) ||
           (get<ATP_concentration>(x[i])() == 0))) {
        out.push_back(x[i]);
      } else {
        auto n1 = get<number_of_samples>(out.back())();
        auto n2 = get<number_of_samples>(x[i])();
        auto new_ATP = (get<ATP_concentration>(out.back())() * n1 +
                        get<ATP_concentration>(x[i])() * n2) /
                       (n1 + n2);
        get<number_of_samples>(out.back())() = n1 + n2;
        get<ATP_concentration>(out.back())() = new_ATP;
      }
    }
    return out;
  }
}

// static ATP_evolution average_ATP_step(ATP_evolution const & x){
//     return std::visit([] (auto const& a){return
//     average_ATP_step(a);},x());}

inline ATP_evolution average_ATP_step(ATP_evolution const &x,
                                      bool average_the_evolution) {
  if (average_the_evolution)
    return average_ATP_step(x());
  else
    return x;
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    uses_variance_correction_aproximation variance_correction, class Model,
    class Parameters, class Variables>
auto simulate(mt_64i &mt,
              const Likelihood_Model<adaptive, recursive, averaging, variance,
                                     variance_correction, Model> &lik,
              Parameters const &p, const Variables &var) {
  return Macro_DMR{}.sample(mt, lik.m, p, var, lik.n_sub_dt).value()();
}

static Experiment_step average_Experimental_step(Experiment_step const &x,
                                                 bool average_the_evolution) {
  return Experiment_step(get<Time>(x), average_ATP_step(get<ATP_evolution>(x),
                                                        average_the_evolution));
}

static ATP_evolution add_ATP_step_i(std::vector<ATP_step> &&x,
                                    std::vector<ATP_step> &&y) {
  if (x.empty())
    return std::move(y);
  else {
    for (std::size_t i = 0; i < y.size(); ++i) {
      if (get<ATP_concentration>(x.back())() ==
          get<ATP_concentration>(y[i])()) {
        get<number_of_samples>(x.back())() += get<number_of_samples>(y[i])();
      } else {
        x.push_back(y[i]);
      }
    }
    return std::move(x);
  }
}

static ATP_evolution add_ATP_step_i(ATP_step &&x, std::vector<ATP_step> &&y) {
  return add_ATP_step_i(std::vector<ATP_step>{std::move(x)}, std::move(y));
}

static ATP_evolution add_ATP_step_i(std::vector<ATP_step> &&x, ATP_step &&e) {
  if (x.empty()) {
    x.push_back(e);
    return x;
  } else if (get<ATP_concentration>(x.back())() ==
             get<ATP_concentration>(e)()) {
    get<number_of_samples>(x.back())() += get<number_of_samples>(e)();
  } else {
    x.push_back(std::move(e));
  }
  return x;
}
static ATP_evolution add_ATP_step_i(ATP_step &&x, ATP_step &&e) {
  if (get<ATP_concentration>(x)() == get<ATP_concentration>(e)()) {
    get<number_of_samples>(x)() += get<number_of_samples>(e)();
    return std::vector<ATP_step>{std::move(x)};
  } else {
    return std::vector<ATP_step>{std::move(x), std::move(e)};
  }
}

static ATP_evolution add_ATP_step(ATP_evolution &&x, ATP_evolution &&e) {
  return add_ATP_step_i(std::move(x()), std::move(e()));
}

class experiment_fractioner {
  std::vector<std::size_t> segments = {73, 33, 22, 22, 4};
  bool average_ATP_evolution = 10;

public:
  experiment_fractioner(const std::vector<std::size_t> &t_segments,
                        bool average_the_ATP_evolution)
      : segments{t_segments}, average_ATP_evolution{average_the_ATP_evolution} {
  }

  static auto average_Recording(const Recording_conditions &e,
                                const Recording &y,
                                const std::vector<std::size_t> &indexes0,
                                const std::vector<std::size_t> &indexes1,
                                bool average_the_evolution) {

    assert(size(y()) == size(indexes0));
    assert(size(e()) == size(y()));
    auto out_size = indexes1.size();
    std::vector<Patch_current> out_y(out_size);
    std::vector<Experiment_step> out_x(out_size);
    double sum_y = 0;
    std::size_t sum_samples = 0;
    std::size_t ii = 0;

    ATP_evolution v_ATP = std::vector<ATP_step>{};

    for (std::size_t i = 0; i < size(y()); ++i) {
      if (indexes0[i] == indexes1[ii]) {
        if (sum_samples == 0) {
          out_y[ii] = y()[i];
          out_x[ii] = average_Experimental_step(e()[i], average_the_evolution);
        } else {

          auto n_samples = get_num_samples(e()[i]);
          sum_y += y()[i]() * n_samples;
          sum_samples += n_samples;
          v_ATP = add_ATP_step(std::move(v_ATP),
                               average_ATP_step(e()[i], average_the_evolution));

          out_y[ii] = Patch_current(sum_y / sum_samples);
          out_x[ii] = Experiment_step(get<Time>(e()[i]), v_ATP);
          sum_y = 0;
          sum_samples = 0;
          v_ATP = std::vector<ATP_step>{};
        }
        ++ii;
      } else {
        assert(indexes0[i] < indexes1[ii] && "indexes fails");
        auto n_samples = get_num_samples(e()[i]);
        sum_y += y()[i]() * n_samples;
        sum_samples += n_samples;
        v_ATP = add_ATP_step(std::move(v_ATP),
                             average_ATP_step(e()[i], average_the_evolution));
      }
    }

    return std::tuple(Recording_conditions(out_x), Recording(out_y));
  }

  template <includes_N_state_evolution keep_N_state>
  static auto average_Recording(const Recording_conditions &e,
                                const Simulated_Recording<keep_N_state> &sim,
                                const std::vector<std::size_t> &indexes0,
                                const std::vector<std::size_t> &indexes1,
                                bool average_the_evolution) {
    auto &yr = get<Recording>(sim());
    assert(size(yr()) == size(indexes0));
    assert(size(e()) == size(yr()));
    auto out_size = indexes1.size();
    std::vector<Patch_current> out_y(out_size);
    std::vector<N_channel_state> out_N(out_size);

    std::vector<Experiment_step> out_x(out_size);
    double sum_y = 0;
    std::size_t sum_samples = 0;
    std::size_t ii = 0;

    ATP_evolution v_ATP = std::vector<ATP_step>{};

    for (std::size_t i = 0; i < size(yr()); ++i) {
      if (indexes0[i] == indexes1[ii]) {
        if (sum_samples == 0) {
          out_y[ii] = yr()[i];
          if constexpr (keep_N_state.value)
            out_N[ii] = get<N_Ch_State_Evolution>(sim())()[i];
          out_x[ii] = average_Experimental_step(e()[i], average_the_evolution);
        } else {

          auto n_samples = get_num_samples(e()[i]);
          sum_y += yr()[i]() * n_samples;
          sum_samples += n_samples;
          v_ATP = add_ATP_step(std::move(v_ATP),
                               average_ATP_step(e()[i], average_the_evolution));

          out_y[ii] = Patch_current(sum_y / sum_samples);

          out_x[ii] = Experiment_step(get<Time>(e()[i]), v_ATP);
          if constexpr (keep_N_state.value)
            out_N[ii] = get<N_Ch_State_Evolution>(sim())()[i];

          sum_y = 0;
          sum_samples = 0;
          v_ATP = std::vector<ATP_step>{};
        }
        ++ii;
      } else {
        assert(indexes0[i] < indexes1[ii] && "indexes fails");
        auto n_samples = get_num_samples(e()[i]);
        sum_y += yr()[i]() * n_samples;
        sum_samples += n_samples;
        v_ATP = add_ATP_step(std::move(v_ATP),
                             average_ATP_step(e()[i], average_the_evolution));
      }
    }
    Simulated_Recording<keep_N_state> out_sim;
    get<Recording>(out_sim())() = std::move(out_y);
    if constexpr (keep_N_state.value)
      get<N_Ch_State_Evolution>(out_sim())() = std::move(out_N);

    return std::tuple(Recording_conditions(out_x), std::move(out_sim));
  }

  auto operator()(const Recording &y, const Experiment &x, mt_64i &mt,
                  std::size_t num_parameters, double n_points_per_decade_beta,
                  double n_points_per_decade_fraction, double stops_at,
                  bool includes_zero) const {
    assert(size(y()) == size(get<Recording_conditions>(x)()));
    assert(size(y()) == var::sum(segments));

    std::size_t num_samples = size(y());
    std::size_t max_num_samples_per_segment = var::max(segments);

    auto cum_segments = var::cumsum(segments);

    auto indexes = generate_random_Indexes(
        mt, num_samples, 1, n_points_per_decade_fraction, cum_segments);
    // std::cerr <<
    // "\nindexes\n**************************************************"
    //              "*************************\n";
    // std::cerr << indexes;
    //  std::abort();
    auto n_frac = size(indexes);
    deprecated::by_fraction<Recording> y_out(n_frac);
    deprecated::by_fraction<Experiment> x_out(
        n_frac,
        Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                   get<initial_ATP_concentration>(x)));
    y_out[n_frac - 1] = y;
    x_out[n_frac - 1] = x;

    for (std::size_t i = n_frac - 1; i > 0; --i) {
      std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
          average_Recording(get<Recording_conditions>(x_out[i]), y_out[i],
                            indexes[i], indexes[i - 1], average_ATP_evolution);
    }

    // std::abort();
    auto beta0 = get_beta_list(n_points_per_decade_beta,
                               stops_at * num_samples /
                                   (n_frac > 1 ? size(indexes[0]) : 1),
                               includes_zero);
    by_beta<double> betan = {0, 1};
    deprecated::by_fraction<by_beta<double>> beta(n_frac, betan);
    beta[0] = std::move(beta0);

    return std::tuple(std::move(y_out), std::move(x_out), std::move(beta));
  }

  auto operator()(const Recording &y, const Experiment &x, mt_64i &mt,
                  std::size_t num_parameters,
                  double n_points_per_decade_fraction) const {
    assert(size(y()) == size(get<Recording_conditions>(x)()));
    assert(size(y()) == var::sum(segments));

    std::size_t num_samples = size(y());
    //  std::size_t max_num_samples_per_segment = var::max(segments);

    auto cum_segments = var::cumsum(segments);

    auto indexes =
        generate_random_Indexes(mt, num_samples, num_parameters,
                                n_points_per_decade_fraction, cum_segments);
    // std::cerr <<
    // "\nindexes\n**************************************************"
    //              "*************************\n";
    // std::cerr << indexes;
    // std::abort();
    auto n_frac = size(indexes);
    deprecated::by_fraction<Recording> y_out(n_frac);
    deprecated::by_fraction<Experiment> x_out(
        n_frac,
        Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                   get<initial_ATP_concentration>(x)));
    y_out[n_frac - 1] = y;
    x_out[n_frac - 1] = x;

    for (std::size_t i = n_frac - 1; i > 0; --i) {
      std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
          average_Recording(get<Recording_conditions>(x_out[i]), y_out[i],
                            indexes[i], indexes[i - 1], average_ATP_evolution);
    }
    return std::tuple(std::move(y_out), std::move(x_out));
  }

  template <includes_N_state_evolution keep_N_state>
  auto operator()(const Simulated_Recording<keep_N_state> &sim,
                  const Experiment &x, mt_64i &mt, std::size_t num_parameters,
                  double n_points_per_decade_fraction) const {
    auto &yr = get<Recording>(sim());
    assert(size(yr()) == size(get<Recording_conditions>(x)()));
    assert(size(segments) == 0 || size(yr()) == var::sum(segments));

    std::size_t num_samples = size(yr());
    //  std::size_t max_num_samples_per_segment = var::max(segments);

    auto cum_segments = var::cumsum(segments);

    auto indexes =
        generate_random_Indexes(mt, num_samples, num_parameters,
                                n_points_per_decade_fraction, cum_segments);
    // std::cerr <<
    // "\nindexes\n**************************************************"
    //              "*************************\n";
    // std::cerr << indexes;
    // std::abort();
    auto n_frac = size(indexes);
    deprecated::by_fraction<Simulated_Recording<keep_N_state>> y_out(n_frac);
    deprecated::by_fraction<Experiment> x_out(
        n_frac,
        Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                   get<initial_ATP_concentration>(x)));
    y_out[n_frac - 1] = sim;
    x_out[n_frac - 1] = x;

    for (std::size_t i = n_frac - 1; i > 0; --i) {
      std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
          average_Recording(get<Recording_conditions>(x_out[i]), y_out[i],
                            indexes[i], indexes[i - 1], average_ATP_evolution);
    }
    return std::tuple(std::move(y_out), std::move(x_out));
  }
};

template <class Id, class... Ts>
void report_title(
    save_Predictions<var::Parameters_transformed<Id>> &s,
    deprecated::cuevi_mcmc<var::Parameters_transformed<Id>> const &,
    const Ts &...t) {

  s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
      << "iter_time" << s.sep << "nsamples" << s.sep << "beta" << s.sep
      << "i_walker" << s.sep << "id_walker" << s.sep << "i_step" << s.sep
      << "time" << s.sep << "num_samples" << s.sep << "ATP" << s.sep
      << "ATPevol" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_var"
      << s.sep << "plogL" << s.sep << "pelogL"
      << "\n";
}

template <class Id>
void report_title(save_Predictions<var::Parameters_transformed<Id>> &s,
                  thermo_mcmc<var::Parameters_transformed<Id>> const &, ...) {

  s.f << "n_betas" << s.sep << "iter" << s.sep << "iter_time" << s.sep << "beta"
      << s.sep << "i_walker" << s.sep << "id_walker" << s.sep << "i_step"
      << s.sep << "time" << s.sep << "num_samples" << s.sep << "ATP" << s.sep
      << "ATP_evolution" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep
      << "Y_var" << s.sep << "plogL" << s.sep << "pelogL"
      << "\n";
}

inline void report_title(save_Predictions<Matrix<double>> &s,
                         thermo_mcmc<Matrix<double>> const &, ...) {}

template <class Id, class FunctionTable, class Duration>
  requires(!is_of_this_template_type_v<std::decay_t<FunctionTable>, FuncMap_St>)
void report(FunctionTable &&, std::size_t iter, const Duration &dur,
            save_Predictions<var::Parameters_transformed<Id>> &s,
            thermo_mcmc<var::Parameters_transformed<Id>> const &data, ...) {
  if (iter % s.save_every == 0)
    for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
      for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker)
        for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)

          s.f << num_betas(data) << s.sep << iter << s.sep << dur << s.sep
              << data.beta[i_beta] << s.sep << i_walker << s.sep
              << data.i_walkers[i_walker][i_beta] << s.sep << i_par << s.sep
              << data.walkers[i_walker][i_beta].parameter[i_par] << "\n";
}

template <class Id, class FunctionTable, class Duration, class Prior,
          class t_logLikelihood, class Data, class Variables>
  requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, FuncMap_St>)
void report(FunctionTable &f, std::size_t iter, const Duration &dur,
            save_Predictions<var::Parameters_transformed<Id>> &s,
            thermo_mcmc<var::Parameters_transformed<Id>> const &data, Prior &&,
            t_logLikelihood &&lik, const Data &y, const Variables &x, ...) {

  if (iter % s.save_every != 0)
    return;
  auto ff = f.fork(omp_get_max_threads());

  auto all_Predictions =
      std::vector<std::vector<std::decay_t<decltype(logLikelihoodPredictions(
          ff[0], lik, data.get_Parameter(0, 0), y, x))>>>(
          data.get_Walkers_number());

  auto beta = data.get_Beta();
  auto num_samples = size(y);
#pragma omp parallel for
  for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
       ++i_walker) {
    for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
      auto i_th = omp_get_thread_num();

      auto par = data.get_Parameter(i_walker, i_b);
      auto walker_id = data.get_Walker(i_walker, i_b);
      all_Predictions[i_walker].push_back(
          logLikelihoodPredictions(ff[i_th], lik, par.to_value(), y, x));
    }
  }
  f += ff;
  for (std::size_t half = 0; half < 2; ++half) {
    for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
      auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
      for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
        auto par = data.get_Parameter(i_walker, i_b);
        auto walker_id = data.get_Walker(i_walker, i_b);
        auto prediction = all_Predictions[i_walker][i_b];
        if (is_valid(prediction)) {
          auto &predictions = prediction.value();
          for (std::size_t i_step = 0; i_step < size(y); ++i_step) {
            auto v_ev =
                get<ATP_evolution>(get<Recording_conditions>(x)()[i_step]);

            auto time = get<Time>(get<Recording_conditions>(x)()[i_step]);
            auto num_smples = get_num_samples(v_ev);

            s.f << beta.size() << s.sep << iter << s.sep << dur << s.sep
                << beta[i_b] << s.sep << i_walker << s.sep << walker_id << s.sep
                << i_step << s.sep << time << s.sep << num_samples << s.sep
                << ToString(average_ATP_step(v_ev, true)) << s.sep
                << ToString(v_ev) << s.sep << y()[i_step]() << s.sep
                << get<y_mean>(predictions()[i_step]) << s.sep
                << get<y_var>(predictions()[i_step]) << s.sep
                << get<plogL>(predictions()[i_step]) << s.sep
                << get<eplogL>(predictions()[i_step]) << "\n";
          }
        }
      }
    }
  }
}

inline std::string ToString(const ATP_step &ev) {
  return std::to_string(get<number_of_samples>(ev)()) + "  " +
         std::to_string(get<ATP_concentration>(ev)());
}

inline std::string ToString(const std::vector<ATP_step> &ev) {
  std::string out;
  for (auto &e : ev) {
    out += ToString(e) + "  ";
  }
  return out;
}

inline std::string ToString(const ATP_evolution &ev) { return ToString(ev()); }

template <includes_N_state_evolution keep_N_state>
void report(std::string filename, const Patch_State_Evolution &predictions,
            const Simulated_Recording<keep_N_state> &y, const Experiment &xs) {
  auto &ys = get<Recording>(y());
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "i_step"
    << ","
    << "time"
    << ","
    << "num_samples"
    << ","
    << "ATP_step"
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
  if constexpr (keep_N_state.value) {
    f << ","
      << "N";
  }
  f << "\n";
  for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
    auto v_ev = get<ATP_evolution>(get<Recording_conditions>(xs)()[i_step]);
    for (std::size_t i = 0; i < get<P_Cov>(predictions()[i_step])().nrows();
         ++i) {
      for (std::size_t j = 0; j < get<P_Cov>(predictions()[i_step])().ncols();
           ++j) {
        f << i_step << "," << get<Time>(get<Recording_conditions>(xs)()[i_step])
          << "," << get_num_samples(v_ev) << ","
          << ToString(average_ATP_step(v_ev, true)) << "," << ToString(v_ev)
          << "," << ys()[i_step]() << "," << get<y_mean>(predictions()[i_step])
          << "," << get<y_var>(predictions()[i_step]) << ","
          << get<plogL>(predictions()[i_step]) << ","
          << get<eplogL>(predictions()[i_step]) << ","
          << get<logL>(predictions()[i_step]) << "," << i << ","
          << get<P_mean>(predictions()[i_step])()[i] << "," << j << ","
          << get<P_Cov>(predictions()[i_step])()(i, j);
        if constexpr (keep_N_state.value)
          f << "," << get<N_Ch_State_Evolution>(y())()[i_step]()[i] << "\n";
        else
          f << "\n";
      }
    }
  }
}

template <includes_N_state_evolution keep_N_state>
void save_Likelihood_Predictions(std::string filename,
                                 const Patch_State_Evolution &predictions,
                                 const Simulated_Recording<keep_N_state> &y,
                                 const Experiment &xs) {
  auto &ys = get<Recording>(y());
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "i_step"
    << ","
    << "time"
    << ","
    << "num_samples"
    << ","
    << "ATP_step"
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
  if constexpr (keep_N_state.value)
    f << ","
      << "N"
      << "\n";
  else
    f << "\n";
  for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
    auto v_ev = get<ATP_evolution>(get<Recording_conditions>(xs)()[i_step]);
    for (std::size_t i = 0; i < get<P_Cov>(predictions()[i_step])().nrows();
         ++i) {
      for (std::size_t j = 0; j < get<P_Cov>(predictions()[i_step])().ncols();
           ++j) {
        f << i_step << "," << get<Time>(get<Recording_conditions>(xs)()[i_step])
          << "," << get_num_samples(v_ev) << ","
          << ToString(average_ATP_step(v_ev, true)) << "," << ToString(v_ev)
          << "," << ys()[i_step]() << "," << get<y_mean>(predictions()[i_step])
          << "," << get<y_var>(predictions()[i_step]) << ","
          << get<plogL>(predictions()[i_step]) << ","
          << get<eplogL>(predictions()[i_step]) << ","
          << get<logL>(predictions()[i_step]) << "," << i << ","
          << get<P_mean>(predictions()[i_step])()[i] << "," << j << ","
          << get<P_Cov>(predictions()[i_step])()(i, j);
        if constexpr (keep_N_state.value)
          f << "," << get<N_Ch_State_Evolution>(y())()[i_step]()[i] << "\n";
        else
          f << "\n";
      }
    }
  }
}

inline void save_Likelihood_Predictions(
    std::string filename, const Patch_State_Evolution &predictions,
    const v_Simulated_Recording &sim_y, const Experiment &xs) {
  std::visit(
      [&](auto &y) {
        save_Likelihood_Predictions(filename, predictions, y, xs);
      },
      sim_y);
}

template <includes_N_state_evolution keep_N_state>
void save_fractioned_Likelihood_Predictions(
    std::string filename, const std::string &sep,
    const std::vector<Patch_State_Evolution> &v_predictions,
    const std::vector<Simulated_Recording<keep_N_state>> &v_y,
    const std::vector<Experiment> &v_xs) {
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  f << "i_frac" << sep << "i_step" << sep << "t_ini" << sep << "time" << sep
    << "i_sub_step" << sep << "number_of_samples" << sep << "ATP_concentration"
    << sep << "macror_algorithm" << sep << "y" << sep << "y_mean" << sep
    << "y_var" << sep << "plogL" << sep << "eplogL" << sep << "logL" << sep
    << "i_state" << sep << "P_mean" << sep << "j_state" << sep << "P_Cov";
  if constexpr (keep_N_state.value)
    f << ","
      << "N"
      << "\n";
  else
    f << "\n";
  for (std::size_t i_frac = 0; i_frac < v_y.size(); ++i_frac) {
    double t_ini = 0;
    auto &y = v_y[i_frac];
    auto &ys = get<Recording>(y());
    auto &xs = v_xs[i_frac];
    auto &predictions = v_predictions[i_frac];
    for (std::size_t i_step = 0; i_step < size(ys); ++i_step) {
      auto v_ev = get<ATP_evolution>(get<Recording_conditions>(xs)()[i_step]);
      for (std::size_t i_sub_step = 0; i_sub_step < v_ev.size(); ++i_sub_step) {
        for (std::size_t i_state = 0;
             i_state < get<P_Cov>(predictions()[i_step])().nrows(); ++i_state) {
          for (std::size_t j_state = 0;
               j_state < get<P_Cov>(predictions()[i_step])().ncols();
               ++j_state) {
            f << i_frac << sep << i_step << sep << t_ini << sep
              << get<Time>(get<Recording_conditions>(xs)()[i_step]) << sep
              << i_sub_step << sep << get<number_of_samples>(v_ev[i_sub_step])
              << sep << get<ATP_concentration>(v_ev[i_sub_step]) << sep
              << get<macror_algorithm>(predictions()[i_step]) << sep
              << ys()[i_step]() << sep << get<y_mean>(predictions()[i_step])
              << sep << get<y_var>(predictions()[i_step]) << sep
              << get<plogL>(predictions()[i_step]) << sep
              << get<eplogL>(predictions()[i_step]) << sep
              << get<logL>(predictions()[i_step]) << sep << i_state << sep
              << get<P_mean>(predictions()[i_step])()[i_state] << sep << j_state
              << sep << get<P_Cov>(predictions()[i_step])()(i_state, j_state);
            if constexpr (keep_N_state.value)
              f << sep << get<N_Ch_State_Evolution>(y())()[i_step]()[i_state]
                << "\n";
            else
              f << "\n";
          }
        }
        t_ini += get<number_of_samples>(v_ev[i_sub_step])() /
                 get<Frequency_of_Sampling>(xs)();
      }
    }
  }
}
template <class FunctionTable, class Duration, class Prior, class Likelihood,
          class Variables, class DataType, class Parameters>
void report(FunctionTable &&f, std::size_t iter, const Duration &dur,
            save_Predictions<Parameters> &s,
            deprecated::cuevi_mcmc<Parameters> const &data, Prior const &prior,
            Likelihood const &lik, const DataType &ys, const Variables &xs,
            ...) {
  if (iter % s.save_every == 0)
    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
      for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
        for (std::size_t i_walker = 0; i_walker < size(data.walkers);
             ++i_walker) {
          Maybe_error<Patch_State_Evolution> prediction =
              logLikelihoodPredictions(
                  f.fork(var::I_thread(i_walker)), lik,
                  data.walkers[i_walker][i_frac][i_beta].parameter, ys[i_frac],
                  xs[i_frac]);
          if (is_valid(prediction)) {
            auto &predictions = prediction.value();
            for (std::size_t i_step = 0; i_step < size(ys[i_frac]); ++i_step) {
              auto v_ev = get<ATP_evolution>(
                  get<Recording_conditions>(xs[i_frac])()[i_step]);

              s.f << size(data.beta) << s.sep << size(data.beta[i_frac])
                  << s.sep << iter << s.sep << dur << s.sep
                  << data.nsamples[i_frac] << s.sep << data.beta[i_frac][i_beta]
                  << s.sep << i_walker << s.sep
                  << data.i_walkers[i_walker][i_frac][i_beta] << s.sep << i_step
                  << s.sep
                  << get<Time>(get<Recording_conditions>(xs[i_frac])()[i_step])
                  << s.sep << get_num_samples(v_ev) << s.sep
                  << ToString(average_ATP_step(v_ev)) << s.sep << ToString(v_ev)
                  << s.sep << ys[i_frac]()[i_step]() << s.sep
                  << get<y_mean>(predictions()[i_step]) << s.sep
                  << get<y_var>(predictions()[i_step]) << s.sep
                  << get<plogL>(predictions()[i_step]) << s.sep
                  << get<eplogL>(predictions()[i_step]) << "\n";
            }
          }
        }
}

template <class ParameterType, class FunctionTable, class Duration, class Prior,
          class t_logLikelihood, class Data, class Variables>
  requires(!is_of_this_template_type_v<FunctionTable, FuncMap_St>)
void report(FunctionTable &&f, std::size_t iter, const Duration &dur,
            save_Predictions<ParameterType> &s,
            cuevi::Cuevi_mcmc<ParameterType> &data, Prior &&,
            t_logLikelihood &&lik, const deprecated::by_fraction<Data> &ys,
            const deprecated::by_fraction<Variables> &xs, ...) {

  auto &t = data.get_Cuevi_Temperatures();
  if (iter % s.save_every == 0) {
    data.calculate_Likelihoods_for_Evidence_calulation(f, lik, ys, xs);
    for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number();
         ++i_cu) {
      auto icu = cuevi::Cuevi_Index(i_cu);
      auto i_frac = data.get_Fraction(i_cu);
      auto beta = data.get_Beta(icu);
      auto nsamples = size(ys[i_frac()]);
      for (std::size_t half = 0; half < 2; ++half)
        // #pragma omp parallel for
        for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
          auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
          auto iw = cuevi::Walker_Index(i_walker);
          auto &wa = data.get_Walker(iw, icu);
          auto &wav = data.get_Walker_Value(iw, icu);
          auto par = data.get_Parameter(iw, i_cu);
          auto prediction = logLikelihoodPredictions(
              f.fork(var::I_thread(iiw)), lik, par, ys[i_frac()], xs[i_frac()]);
          if (is_valid(prediction)) {
            auto &predictions = prediction.value();
            for (std::size_t i_step = 0; i_step < size(ys[i_frac()]);
                 ++i_step) {
              auto v_ev = get<ATP_evolution>(
                  get<Recording_conditions>(xs[i_frac()])()[i_step]);

              s.f << iter << s.sep << dur << s.sep << i_cu << s.sep << i_frac()
                  << s.sep << nsamples << s.sep << beta() << s.sep << i_walker
                  << s.sep << get<cuevi::Walker_id>(wa())() << s.sep << i_step
                  << s.sep
                  << get<Time>(
                         get<Recording_conditions>(xs[i_frac()])()[i_step])
                  << s.sep << get_num_samples(v_ev) << s.sep
                  << ToString(average_ATP_step(v_ev)) << s.sep << ToString(v_ev)
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
}

template <class ParameterType, class FunctionTable, class Duration, class Prior,
          class t_logLikelihood, class Data, class Variables>
  requires(is_of_this_template_type_v<FunctionTable, FuncMap_St>)
void report(FunctionTable &f, std::size_t iter, const Duration &dur,
            save_Predictions<ParameterType> &s,
            cuevi::Cuevi_mcmc<ParameterType> &data, Prior &&,
            t_logLikelihood &&lik, const cuevi::by_fraction<Data> &ys,
            const cuevi::by_fraction<Variables> &xs, ...) {

  auto &t = data.get_Cuevi_Temperatures();
  if (iter % s.save_every != 0) {
    return;
  }
  data.calculate_Likelihoods_for_Evidence_calulation(f, lik, ys, xs);

  auto ff = f.fork(omp_get_max_threads());
  using Predition_type = std::decay_t<decltype(logLikelihoodPredictions(
      ff[0], lik,
      data.get_Parameter(cuevi::Walker_Index(0), cuevi::Cuevi_Index(0)), ys[0],
      xs[0]))>;
  auto allPredictions = std::vector<std::vector<Predition_type>>(
      data.get_Cuevi_Temperatures_Number(),
      std::vector<Predition_type>(data.get_Walkers_number()));
#pragma omp parallel for //collapse(2)
  for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
       ++i_walker) {
    for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number();
         ++i_cu) {
      auto iw = cuevi::Walker_Index(i_walker);
      auto i_th = omp_get_thread_num();

      auto icu = cuevi::Cuevi_Index(i_cu);
      auto i_frac = data.get_Fraction(i_cu);
      auto beta = data.get_Beta(icu);
      auto nsamples = size(ys[i_frac()]);
      auto &wa = data.get_Walker(iw, icu);
      auto &wav = data.get_Walker_Value(iw, icu);
      auto par = data.get_Parameter(iw, i_cu);
      auto prediction = logLikelihoodPredictions(ff[i_th], lik, par.to_value(),
                                                 ys[i_frac()], xs[i_frac()]);
      allPredictions[i_cu][i_walker] = std::move(prediction);
    }
  }
  f += ff;
  for (std::size_t i_cu = 0; i_cu < data.get_Cuevi_Temperatures_Number();
       ++i_cu) {
    auto icu = cuevi::Cuevi_Index(i_cu);
    auto i_frac = data.get_Fraction(i_cu);
    auto beta = data.get_Beta(icu);
    auto nsamples = size(ys[i_frac()]);
    for (std::size_t half = 0; half < 2; ++half)
      for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
        auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
        auto iw = cuevi::Walker_Index(i_walker);
        auto &wa = data.get_Walker(iw, icu);
        auto &wav = data.get_Walker_Value(iw, icu);
        auto par = data.get_Parameter(iw, i_cu);
        auto &prediction = allPredictions[i_cu][i_walker];
        if (is_valid(prediction)) {
          auto &predictions = prediction.value();
          for (std::size_t i_step = 0; i_step < size(ys[i_frac()]); ++i_step) {
            auto v_ev = get<ATP_evolution>(
                get<Recording_conditions>(xs[i_frac()])()[i_step]);

            s.f << iter << s.sep << dur << s.sep << i_cu << s.sep << i_frac()
                << s.sep << nsamples << s.sep << beta() << s.sep << i_walker
                << s.sep << get<cuevi::Walker_id>(wa())() << s.sep << i_step
                << s.sep
                << get<Time>(get<Recording_conditions>(xs[i_frac()])()[i_step])
                << s.sep << get_num_samples(v_ev) << s.sep
                << ToString(average_ATP_step(v_ev, true)) << s.sep
                << ToString(v_ev) << s.sep << ys[i_frac()]()[i_step]() << s.sep
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
void report_title(save_Predictions<ParameterType> &s,
                  cuevi::Cuevi_mcmc<ParameterType> const &, ...) {

  s.f << "iter" << s.sep << "iter_time" << s.sep << "i_cu" << s.sep << "i_frac"
      << s.sep << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep
      << "Walker_id" << s.sep << "i_step" << s.sep << "Time" << s.sep
      << "num_samples" << s.sep << "average_ATP_step" << s.sep << "v_ev"
      << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_var" << s.sep
      << "plogL" << s.sep << "eplogL"
      << "\n";
}

template <class Id>
auto cuevi_Model_by_convergence(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments, bool average_the_ATP_evolution,

    std::size_t num_scouts_per_ensemble, double min_fraction,
    std::size_t thermo_jumps_every, std::size_t max_iter, double max_ratio,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double stops_at, bool includes_zero, std::size_t initseed) {
  return cuevi_integration(
      checks_derivative_var_ratio<deprecated::cuevi_mcmc,
                                  var::Parameters_transformed<Id>>(max_iter,
                                                                   max_ratio),
      experiment_fractioner(t_segments, average_the_ATP_evolution),
      save_mcmc<var::Parameters_transformed<Id>,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, 100ul, 100ul, 100ul, 100ul),
      num_scouts_per_ensemble, min_fraction, thermo_jumps_every,
      n_points_per_decade_beta, n_points_per_decade_fraction, stops_at,
      includes_zero, initseed);
}

template <class Id>
auto cuevi_Model_by_iteration(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments, bool average_the_ATP_evolution,

    std::size_t num_scouts_per_ensemble,
    std::size_t max_number_of_simultaneous_temperatures, double min_fraction,
    std::size_t thermo_jumps_every, std::size_t max_iter_warming,
    std::size_t max_iter_equilibrium, double max_ratio,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double stops_at, bool includes_zero, std::size_t initseed) {
  return deprecated::cuevi_integration(
      less_than_max_iteration(max_iter_warming, max_iter_equilibrium),
      experiment_fractioner(t_segments, average_the_ATP_evolution),
      save_mcmc<var::Parameters_transformed<Id>,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, 10ul, 100ul, 10ul, 100ul),
      num_scouts_per_ensemble, max_number_of_simultaneous_temperatures,
      min_fraction, thermo_jumps_every, n_points_per_decade_beta,
      n_points_per_decade_fraction, stops_at, includes_zero, initseed);
}

template <class Id>
cuevi::Cuevi_Algorithm<
    experiment_fractioner,
    save_mcmc<var::Parameters_transformed<Id>,
              save_likelihood<var::Parameters_transformed<Id>>,
              save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
              save_Predictions<var::Parameters_transformed<Id>>>,
    cuevi_less_than_max_iteration>
new_cuevi_Model_by_iteration(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments, bool average_the_ATP_evolution,
    std::size_t num_scouts_per_ensemble,
    std::size_t number_trials_until_give_up, double min_fraction,
    std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double medium_beta, double stops_at, bool includes_the_zero,
    Saving_intervals sint, bool random_jumps) {
  return cuevi::Cuevi_Algorithm(
      experiment_fractioner(t_segments, average_the_ATP_evolution),
      save_mcmc<var::Parameters_transformed<Id>,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, get<Save_Likelihood_every>(sint())(),
          get<Save_Parameter_every>(sint())(),
          get<Save_Evidence_every>(sint())(),
          get<Save_Predictions_every>(sint())()),
      cuevi_less_than_max_iteration(max_iter_equilibrium),
      cuevi::Num_Walkers_Per_Ensemble(num_scouts_per_ensemble),
      cuevi::Fractions_Param(
          Vector_Space(cuevi::Min_value(min_fraction),
                       cuevi::Points_per_decade(n_points_per_decade_fraction))),
      cuevi::Th_Beta_Param(
          Vector_Space( // Includes_zero,
                        // Med_value,Points_per_decade,Min_value,
                        // Points_per_decade_low
              cuevi::Includes_zero(includes_the_zero),
              cuevi::Med_value(medium_beta),
              cuevi::Points_per_decade(n_points_per_decade_fraction),
              cuevi::Min_value(stops_at),
              cuevi::Points_per_decade_low(n_points_per_decade_beta))),
      cuevi::Number_trials_until_give_up(number_trials_until_give_up),
      cuevi::Thermo_Jumps_every(thermo_jumps_every),
      cuevi::Random_jumps(random_jumps), std::move(sint));
}

template <class Id>
cuevi::Cuevi_Algorithm_no_Fractioner<
    save_mcmc<var::Parameters_transformed<Id>,
              save_likelihood<var::Parameters_transformed<Id>>,
              save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
              save_Predictions<var::Parameters_transformed<Id>>>,
    cuevi_less_than_max_iteration>
new_cuevi_Model_already_fraction_by_iteration(
    std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
    std::size_t number_trials_until_give_up, std::size_t thermo_jumps_every,
    std::size_t max_iter_equilibrium, double n_points_per_decade_beta_high,
    double n_points_per_decade_beta_low, double medium_beta, double stops_at,
    bool includes_the_zero, Saving_intervals sint, bool random_jumps) {
  return cuevi::Cuevi_Algorithm_no_Fractioner(
      save_mcmc<var::Parameters_transformed<Id>,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, get<Save_Likelihood_every>(sint())(),
          get<Save_Parameter_every>(sint())(),
          get<Save_Evidence_every>(sint())(),
          get<Save_Predictions_every>(sint())()),
      cuevi_less_than_max_iteration(max_iter_equilibrium),
      cuevi::Num_Walkers_Per_Ensemble(num_scouts_per_ensemble),
      cuevi::Th_Beta_Param(
          Vector_Space( // Includes_zero,
                        // Med_value,Points_per_decade,Min_value,
                        // Points_per_decade_low
              cuevi::Includes_zero(includes_the_zero),
              cuevi::Med_value(medium_beta),
              cuevi::Points_per_decade(n_points_per_decade_beta_high),
              cuevi::Min_value(stops_at),
              cuevi::Points_per_decade_low(n_points_per_decade_beta_low))),
      cuevi::Number_trials_until_give_up(number_trials_until_give_up),
      cuevi::Thermo_Jumps_every(thermo_jumps_every),
      cuevi::Random_jumps(random_jumps), std::move(sint));
}
template <class Id>
auto new_thermo_Model_by_max_iter(
    std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
    std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium,
    std::size_t beta_size, std::size_t beta_upper_size,
    std::size_t beta_medium_size, double beta_upper_value,
    double beta_medium_value,

    double stops_at, bool includes_zero, Saving_intervals sint,
    std::size_t initseed) {
  return new_thermodynamic_integration(
      thermo_less_than_max_iteration(max_iter_equilibrium),
      save_mcmc<var::Parameters_transformed<Id>, save_Iter,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, 1ul, get<Save_Likelihood_every>(sint())(),
          get<Save_Parameter_every>(sint())(),
          get<Save_Evidence_every>(sint())(),
          get<Save_Predictions_every>(sint())()),
      num_scouts_per_ensemble, thermo_jumps_every, beta_size, beta_upper_size,
      beta_medium_size, beta_upper_value, beta_medium_value, stops_at,
      includes_zero, initseed);
}

template <class Id>
auto thermo_Model_by_max_iter(std::string path, std::string filename,
                              std::size_t num_scouts_per_ensemble,
                              std::size_t max_num_simultaneous_temperatures,
                              std::size_t thermo_jumps_every,
                              std::size_t max_iter_warming,
                              std::size_t max_iter_equilibrium,
                              double n_points_per_decade, double stops_at,
                              bool includes_zero, std::size_t initseed) {
  return thermodynamic_integration(
      less_than_max_iteration(max_iter_warming, max_iter_equilibrium),
      save_mcmc<var::Parameters_transformed<Id>,
                save_likelihood<var::Parameters_transformed<Id>>,
                save_Parameter<var::Parameters_transformed<Id>>, save_Evidence,
                save_Predictions<var::Parameters_transformed<Id>>>(
          path, filename, 10ul, 10ul, 10ul, 100ul),
      num_scouts_per_ensemble, max_num_simultaneous_temperatures,
      thermo_jumps_every, n_points_per_decade, stops_at, includes_zero,
      initseed);
}

template <class Parameter, includes_N_state_evolution includesN>
void report_model(save_Parameter<Parameter> &s,
                  Simulated_Recording<includesN> const &sim) {
  std::ofstream f(s.fname + "_simulation.csv");
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_step" << s.sep << "patch_current";
  if constexpr (includesN.value)
    f << s.sep << "i_state" << s.sep << "N_state";
  f << "\n";

  if constexpr (includesN.value) {
    auto N = get<N_Ch_State_Evolution>(sim());
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
      for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_step) {
        f << i_step << s.sep << y()[i_step]() << i_state << s.sep
          << N()[i_step]()[i_state] << "\n";
      }
    }
  } else {
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
      f << i_step << s.sep << y()[i_step]() << "\n";
    }
  }
}

template <includes_N_state_evolution includesN>
void save_Simulated_Recording(std::string const &filename,
                              const std::string &separator,
                              Simulated_Recording<includesN> const &sim) {
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_step" << separator << "patch_current";
  if constexpr (includesN.value)
    f << separator << "i_state" << separator << "N_state";
  f << "\n";

  if constexpr (includesN.value) {
    auto N = get<N_Ch_State_Evolution>(sim());
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < N().size(); ++i_step) {
      for (auto i_state = 0ul; i_state < N()[i_step]().size(); ++i_state) {
        f << i_step << separator << y()[i_step]() << separator << i_state
          << separator << N()[i_step]()[i_state] << "\n";
      }
    }
  } else {
    auto y = get<Recording>(sim());
    for (auto i_step = 0ul; i_step < y().size(); ++i_step) {
      f << i_step << separator << y()[i_step]() << "\n";
    }
  }
}

template <class Id>
void report_model(save_Parameter<var::Parameters_transformed<Id>> &s,
                  Parameters_Transformations<Id> const &m) {
  std::ofstream f(s.fname + "_parameter.csv");
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  auto n = m.size();
  f << "i_par" << s.sep << "moment" << s.sep << "value"
    << "\n";
  for (auto i_par = 0ul; i_par < n; ++i_par)
    f << i_par << s.sep << "mean" << s.sep << m[i_par] << "\n";
}
template <class Id>
Maybe_error<Parameters_Transformations<Id>>
load_Parameters(save_Parameter<var::Parameters_transformed<Id>> &s) {
  return load_Parameters<Id>(s.fname + "_parameter.csv", s.sep);
}

namespace cmd {
inline auto set_Fraction_algorithm(double min_fraction,
                                   double n_points_per_decade_fraction,
                                   std::string segments) {
  return std::tuple(min_fraction, n_points_per_decade_fraction, segments);
}

using fraction_algo_type =
    typename return_type<std::decay_t<decltype(&set_Fraction_algorithm)>>::type;

inline Maybe_error<std::vector<std::size_t>>
load_segments_length_for_fractioning(const std::string &filename,
                                     std::string sep) {

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

inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_experiment_fractions(std::string save_name, std::string recording,
                          experiment_type experiment,
                          fraction_algo_type fraction_algo,
                          Maybe_error<std::size_t> Maybe_num_param,
                          std::size_t i_seed) {
  auto myseed = calc_seed(i_seed);

  if (!Maybe_num_param)
    return Maybe_num_param.error();
  auto num_param = Maybe_num_param.value();

  auto init_seed = calc_seed(i_seed);
  mt_64i mt(init_seed);

  auto [min_fraction, n_points_per_decade_fraction, segments] =
      std::move(fraction_algo);

  auto Maybe_segments = load_segments_length_for_fractioning(segments, ",");
  if (!Maybe_segments)
    return Maybe_segments.error();
  Recording y;
  auto Maybe_y = load_Recording_Data(recording, ",", y);

  macrodr::experiment_fractioner frac(Maybe_segments.value(), 0);

  auto [ys, xs] = frac(y, experiment, mt, num_param * min_fraction,
                       n_points_per_decade_fraction);
  auto filename = save_name + "_" + std::to_string(myseed);
  save_fractioned_experiment(filename + "_experiment.csv", ",", xs);
  save_fractioned_Recording(filename + "_recording.csv", ",", ys);
  return std::tuple(filename + "_experiment.csv", filename + "_recording.csv",
                    get<Frequency_of_Sampling>(experiment)(),
                    get<initial_ATP_concentration>(experiment)()());
}

inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_simulation_fractions(std::string save_name, std::string simulation,
                          experiment_type experiment,
                          fraction_algo_type fraction_algo,
                          Maybe_error<std::size_t> Maybe_num_param,
                          std::size_t i_seed) {
  if (!Maybe_num_param)
    return Maybe_num_param.error();
  auto num_param = Maybe_num_param.value();
  auto myseed = calc_seed(i_seed);

  auto init_seed = calc_seed(i_seed);
  mt_64i mt(init_seed);

  auto [min_fraction, n_points_per_decade_fraction, segments] =
      std::move(fraction_algo);

  auto Maybe_segments = load_segments_length_for_fractioning(segments, ",");
  if (!Maybe_segments)
    return Maybe_segments.error();
  Simulated_Recording<includes_N_state_evolution(true)> y;
  auto Maybe_y = load_simulation(simulation, ",", y);
  if (!Maybe_y)
    return Maybe_y.error();

  macrodr::experiment_fractioner frac(Maybe_segments.value(), 0);

  auto [ys, xs] = frac(y, experiment, mt, num_param * min_fraction,
                       n_points_per_decade_fraction);
  auto filename = save_name + "_" + std::to_string(myseed);
  save_fractioned_experiment(filename + "_frac_experiment.csv", ",", xs);
  save_fractioned_simulation(filename + "_frac_recording.csv", ",", ys);
  return std::tuple(filename + "_frac_experiment.csv",
                    filename + "_frac_recording.csv",
                    get<Frequency_of_Sampling>(experiment)(),
                    get<initial_ATP_concentration>(experiment)()());
}

using fractioned_experiment_type = typename return_type<
    std::decay_t<decltype(&calc_experiment_fractions)>>::type;

using fractioned_simulation_type = typename return_type<
    std::decay_t<decltype(&calc_simulation_fractions)>>::type;

} // namespace cmd

} // namespace macrodr

#endif // QMODEL_H
