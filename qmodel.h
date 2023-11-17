#pragma once
#include "cuevi.h"
#include "experiment.h"
#include "fold.h"
#include "function_memoization.h"
#include "matrix.h"
#include <cstddef>
#include <functional>
#include <numeric>
#include <random>
#include <set>
#include <type_traits>
#ifndef QMODEL_H
#define QMODEL_H
#include <map>
#include <string>

#include "derivative_operator.h"
#include "derivative_test.h"
#include "general_algorithm_on_containers.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parameters.h"
#include "variables.h"
#include "exponential_matrix.h"
#include "function_measure_verification_and_optimization.h"
namespace macrodr {

using var::Parameters;
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

using var::FuncMap;
using var::Time_it;
using var::F;

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
class Q0_formula
    : public var::Var<Q0_formula, std::vector<std::vector<std::string>>> {
public:
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

class Qx : public var::Var<Qx, Matrix<double>> {};
class P_initial : public var::Var<P_initial, Matrix<double>> {};

class g : public var::Var<g, Matrix<double>> {};
class g_formula : public var::Var<g_formula, std::vector<std::string>> {};

class N_St : public var::Constant<N_St, std::size_t> {};

class N_Ch_mean : public var::Var<N_Ch_mean, Matrix<double>> {};

class N_Ch_mean_time_segment_duration : public var::Constant<N_Ch_mean_time_segment_duration, double> {};

class N_Ch_init : public var::Var<N_Ch_init, double> {};
class N_Ch_eq : public var::Var<N_Ch_eq, double> {};
class N_Ch_tau : public var::Var<N_Ch_tau, double> {};

class Binomial_magical_number
    : public var::Constant<Binomial_magical_number, double> {};

class min_P : public var::Constant<min_P, double> {};

class N_Ch_std : public var::Var<N_Ch_std, double> {};

class Current_Noise : public var::Var<Current_Noise, double> {};

class Current_Baseline : public var::Var<Current_Baseline, double> {};

class P_mean : public var::Var<P_mean, Matrix<double>> {};

class N_channel_state : public var::Var<N_channel_state, Matrix<double>> {};

class y_sum : public var::Var<y_sum, double> {};

class t_sum : public var::Var<t_sum, double> {};

class P_Cov : public var::Var<P_Cov, SymmetricMatrix<double>> {};

class lambda : public var::Var<lambda, DiagonalMatrix<double>> {};

class V : public var::Var<V, Matrix<double>> {};
class W : public var::Var<W, Matrix<double>> {};

class uses_variance_aproximation
    : public var::struct_Var<uses_variance_aproximation, bool> {};

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

class P : public Var<P, Matrix<double>> {};

class gmean_i : public Var<gmean_i, Matrix<double>> {};
class gtotal_ij : public Var<gtotal_ij, Matrix<double>> {};
class gmean_ij : public Var<gmean_ij, Matrix<double>> {};
class gtotal_sqr_ij : public Var<gtotal_sqr_ij, Matrix<double>> {};
class gsqr_i : public Var<gsqr_i, Matrix<double>> {};
class gvar_i : public Var<gvar_i, Matrix<double>> {};
class gtotal_var_ij : public Var<gtotal_var_ij, Matrix<double>> {};
class gvar_ij : public Var<gvar_ij, Matrix<double>> {};

class y_mean : public var::Var<y_mean, double> {};
class y_var : public var::Var<y_var, double> {};
class plogL : public var::Var<plogL, double> {};
class eplogL : public var::Var<eplogL, double> {};
class vplogL : public var::Var<vplogL, double> {};

class logL : public var::Var<logL, double> {};
class elogL : public var::Var<elogL, double> {};
class vlogL : public var::Var<vlogL, double> {};

class PGn : public var::Var<PGn, Matrix<double>> {};
class PGG_n : public var::Var<PGG_n, Matrix<double>> {};
class PG_n : public var::Var<PG_n, Matrix<double>> {};
// class PPn : public var::Var<PPn, Matrix<double>> {};

using Qn = Vector_Space<number_of_samples, min_P, P, PG_n, PGG_n>;

using Qx_eig = Vector_Space<Qx, V, lambda, W>;

using Qdt =
    Vector_Space<number_of_samples, min_P, P, gmean_i, gtotal_ij, gmean_ij,
                 gtotal_sqr_ij, gsqr_i, gvar_i, gtotal_var_ij, gvar_ij>;

using Patch_Model =
    Vector_Space<N_St, Q0, Qa, P_initial,g, N_Ch_mean,
                                 Current_Noise,Current_Baseline, N_Ch_mean_time_segment_duration,
                 Binomial_magical_number, min_P, Probability_error_tolerance,
                 Conductance_variance_error_tolerance>;

template <class Id> struct Model_Patch {
  template <class F> class Model {
    std::tuple<F, Parameters<Id>, typename Parameters<Id>::Names, Q0_formula,
               Qa_formula, g_formula>
        m_f;

  public:
    static constexpr bool is_Model_Patch = true;
    template <class G> Model(G &&t_g) : m_f{std::forward<G>(t_g)()} {}

    auto &names() const {
      return std::get<typename Parameters<Id>::Names>(m_f);
    }
    auto &parameters() const { return std::get<Parameters<Id>>(m_f); }

    auto &get_Q0_formula() const { return std::get<Q0_formula>(m_f); }
    auto &get_Qa_formula() const { return std::get<Qa_formula>(m_f); }
    auto &get_g_formula() const { return std::get<g_formula>(m_f); }

    template <class P>
      requires std::is_same_v<var::untransformed_type_t<P>, Parameters<Id>>
    auto operator()(const P &t_p) const {
      return std::invoke(std::get<F>(m_f), t_p);
    }
  };
  template <class F>
  Model(F &&f)->Model_Patch<Id>::
      Model<std::tuple_element_t<0, decltype(std::declval<F &&>()())>>;
};

using Patch_State = Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean,
                                 y_var, plogL, eplogL, vplogL>;

template<class C_Patch_Model, class C_double>
C_Patch_Model add_Patch_inactivation(C_Patch_Model&& m, C_double const& deactivation_rate)
{
    using Transf = transformation_type_t<C_Patch_Model>;
    auto Nst=get<N_St>(m)()+1;
    Op_t<Transf, Q0> v_Q0=Q0(Matrix<double>(Nst,Nst,0.0));
    for (std::size_t i=0; i+1<Nst; ++i)
    {
        for (std::size_t j=0; j+1<Nst; ++j)
            set(v_Q0(),i,j,get<Q0>(m)()(i,j));
        set(v_Q0(),i,Nst-1, deactivation_rate);
    }
    Op_t<Transf, Qa> v_Qa=Qa(Matrix<double>(Nst,Nst,0.0));
    for (std::size_t i=0; i+1<Nst; ++i)
    {
        for (std::size_t j=0; j+1<Nst; ++j)
            set(v_Qa(),i,j,get<Qa>(m)()(i,j));
    }
    Op_t<Transf, g> v_g=g(Matrix<double>(Nst,1,0.0));
    for (std::size_t i=0; i+1<Nst; ++i)
    {
        set(v_g(),i,0,get<g>(m)()[i]);
    }
    Op_t<Transf, P_initial> v_Pini=P_initial(Matrix<double>(1,Nst,0.0));
    for (std::size_t i=0; i+1<Nst; ++i)
    {
        set(v_Pini(),0,i,get<P_initial>(m)()[i]);
    }
    
    get<N_St>(m)()=Nst;
    get<Qa>(m)=v_Qa;
    get<Q0>(m)=v_Q0;
    get<g>(m)=v_g;
    get<P_initial>(m)=v_Pini;
    return std::move(m);
 }


class Patch_State_Evolution
    : public Var<Patch_State_Evolution, std::vector<Patch_State>> {};

using Patch_State_and_Evolution =
    Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL,
                 eplogL, vplogL, Patch_State_Evolution>;

class Number_of_simulation_sub_steps
    : public Var<Number_of_simulation_sub_steps, std::size_t> {};

class Simulated_Recording : public Var<Simulated_Recording, Recording> {};

using Simulated_Step = Vector_Space<N_channel_state, Simulated_Recording>;

using Simulated_Sub_Step =
    Vector_Space<N_channel_state, number_of_samples, y_sum>;

using Simulation_Parameters = Vector_Space<Number_of_simulation_sub_steps>;


template <uses_recursive_aproximation recursive,
         uses_averaging_aproximation averaging,
         uses_variance_aproximation variance>
struct MacroR{
    friend std::string ToString(MacroR){
        std::string out="MacroR";
        if (recursive.value)
            out+="_R";
        else
            out+="_NR";
        if (averaging.value==2)
            out+="_2";
        else
            out+="__";
        if (variance.value)
            out+="_V";
        else
            out+="_M";
        
        return out; 
        }
};


struct Calc_Qdt{
    friend std::string ToString(Calc_Qdt){return "Calc_Qdt";}
};

struct Calc_Qdt_step{
    friend std::string ToString(Calc_Qdt_step){return "Calc_Qdt_step";}
};


struct Calc_Qx{
    friend std::string ToString(Calc_Qx){return "Calc_Qx";}
};

struct Calc_eigen{
    friend std::string ToString(Calc_eigen){return "Calc_eigen";}
};



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
  template <class C_Patch_Model>
       requires U<C_Patch_Model, Patch_Model>
  auto calc_Qx(const C_Patch_Model &m, ATP_concentration x)
      -> Transfer_Op_to<C_Patch_Model, Qx> {
    using Trans = transformation_type_t<C_Patch_Model>;
    auto v_Qx = build<Qx>(get<Q0>(m)() + get<Qa>(m)() * x.value());
    Matrix<double> u(v_Qx().ncols(), 1, 1.0);
    v_Qx() = v_Qx() - diag(v_Qx() * u);
    return v_Qx;
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
      else
        return build<Qx_eig>(std::move(v_Qx), build<V>(std::move(v_V)),
                             build<lambda>(std::move(v_l)),
                             build<W>(std::move(Maybe_W.value())));
    } else
      return maybe_eig.error();
  }
  
  template <class C_Patch_Model>
       requires U<C_Patch_Model, Patch_Model>
  auto calc_eigen(const C_Patch_Model &m, ATP_concentration x)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model,  Qx_eig>> {
      return calc_eigen(calc_Qx(m,x));
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

  static auto sample_Multinomial(std::mt19937_64 &mt, P_mean const t_P_mean,
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

  static auto sample_Multinomial(std::mt19937_64 &mt, P const t_P,
                                 N_channel_state N) {
    assert(t_P().nrows() == t_P().ncols());
    auto k = N().size();
    N_channel_state out(Matrix<double>(1, k, 0.0));
    for (std::size_t i = 0; i < k; ++i) {
      std::size_t N_remaining = N()[i];
      double p_remaining = 1;
      for (std::size_t j = 0; j + 1 < k; ++j) {
        auto n = std::binomial_distribution<std::size_t>(
            N_remaining, t_P()(i, j) / p_remaining)(mt);
        N_remaining -= n;
        p_remaining -= t_P()(i, j);
        out()[j] += n;
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
    auto laexp = apply([](auto const & x) { using std::exp; return exp(x); }, ladt);
    
    
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

  template <class C_Patch_Model, class C_Qx>
    requires(/*U<C_Patch_Model, Patch_Model> &&*/ U<C_Qx, Qx>)
  auto calc_Peq(C_Qx const &t_Qx, const C_Patch_Model &m)
      -> Transfer_Op_to<C_Patch_Model, P_mean> {
    auto nstates = get<N_St>(m).value();
    auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);
    auto v_eig_Qx = calc_eigen(t_Qx);
    if (v_eig_Qx) {
        
      auto &landa = get<lambda>(v_eig_Qx.value())();
      auto &Vv = get<V>(v_eig_Qx.value())();
      auto &Wv = get<W>(v_eig_Qx.value())();
      auto ladt = get<lambda>(v_eig_Qx.value())() * 1e8;
     
      auto laexp = apply([](auto const & x) { using std::exp; return exp(x); }, ladt);
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
        std::cerr<<"uses expm_sure\n";
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
    return normalize(P(get<V>(t_Qx)() * exp_ladt * get<W>(t_Qx)()), t_min_P);
  }
  
  template <class Patch_Model>
  auto calc_P(const Patch_Model &m, const Qx &t_Qx, double dt,
              double t_min_P) {
      auto t_eigenQx=calc_eigen(t_Qx);
      if (t_eigenQx)
          return calc_P(m,t_eigenQx.value(),dt,t_min_P);
      else
      
      return normalize(P(expm_sure(t_Qx()*dt)), t_min_P);
  }
  
  
  template <class Patch_Model>
  auto calc_Qdt_old(const Patch_Model &m, const Qx_eig &t_Qx,
                    number_of_samples ns, double dt) {

    auto t_min_P = get<min_P>(m)();
    auto &v_g = get<g>(m);

    std::size_t N = t_Qx[Var<Qx>{}]().ncols();

    auto ladt = t_Qx[Var<lambda>{}]() * dt;
    
    auto exp_ladt = apply([](auto const& x) { using std::exp; return exp(x); }, ladt);
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

  template <class FunctionTable,class C_Patch_Model, class C_Qx_eig>
    requires(U<C_Patch_Model, Patch_Model> &&  U<C_Qx_eig, Qx_eig>)
  auto calc_Qdt( FunctionTable &&f,const C_Patch_Model &m, const C_Qx_eig &t_Qx,
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

    auto r_P = build<P>(t_V() * v_exp_ladt * t_W());

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
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);

    auto r_gtotal_ij = build<gtotal_ij>(t_V() * WgV_E2 * t_W());

    Matrix<Op_t<Trans, double>> WgV_E3(N, N, Op_t<Trans, double>(0.0));
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
          //      std::cerr<<"\t"<<WgV_E3(n1, n3);

          WgV_E3(n1, n3) =
              WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                   E3(v_ladt[n1], v_ladt[n2], v_ladt[n3],
                                      v_exp_ladt[n1], v_exp_ladt[n2],
                                      v_exp_ladt[n3], t_min_P()); // optimizable
        }

    auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);

    if constexpr (false) {
      std::cerr << "\nr_gtotal_sqr_ij\n" << r_gtotal_sqr_ij;
      std::cerr << "\nvar::outside_in(var::inside_out(r_gtotal_sqr_ij))\n"
                << var::outside_in(var::inside_out(r_gtotal_sqr_ij()));

      std::cerr << "\nvar::inside_out(r_gtotal_sqr_ij)\n"
                << var::inside_out(r_gtotal_sqr_ij());
    }

    if constexpr (false) {
      auto test_r_gtotal_sqr_ij = var::test_Derivative(
          [this, &N, &t_min_P](const auto &t_V, const auto &t_W,
                               const auto &v_WgV, const auto &v_ladt,
                               const auto &v_exp_ladt) {
            using Trans2 =
                transformation_type_t<std::decay_t<decltype(v_ladt)>>;

            Matrix<Op_t<Trans2, double>> WgV_E3(N, N,
                                                Op_t<Trans2, double>(0.0));
            for (std::size_t n1 = 0; n1 < N; n1++)
              for (std::size_t n3 = 0; n3 < N; n3++)
                for (std::size_t n2 = 0; n2 < N; n2++) {
                  //      std::cerr<<"\t"<<WgV_E3(n1, n3);

                  WgV_E3(n1, n3) =
                      WgV_E3(n1, n3) + v_WgV(n1, n2) * v_WgV(n2, n3) *
                                           E3(v_ladt[n1], v_ladt[n2],
                                              v_ladt[n3], v_exp_ladt[n1],
                                              v_exp_ladt[n2], v_exp_ladt[n3],
                                              t_min_P()); // optimizable
                }

            //    return var::outside_in(WgV_E3);
            return build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);
          },
          1e-4, 1e-6, t_V, t_W, v_WgV, v_ladt, v_exp_ladt);
      if (!test_r_gtotal_sqr_ij) {
        std::cerr << "\n error in test_r_gtotal_sqr_ij!!\n"
                  << test_r_gtotal_sqr_ij.error()();
        std::abort();
      }
    }

    if constexpr (false) {
      assert(test_conductance_variance(
          r_gtotal_sqr_ij(), get<Conductance_variance_error_tolerance>(m)));
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
    Matrix<double> U(1, t_g().size(), 1.0);
    Matrix<double> UU(t_g().size(), t_g().size(), 1.0);
    auto gmean_ij_p = X_plus_XT(t_g() * U) * (0.5);

    auto gvar_ij_p =
        apply([](auto x) { return abs(x); }, t_g() * U - tr(t_g() * U)) * (0.5);

    auto gmean_ij_tot = r_gtotal_ij(); // + gmean_ij_p * t_min_P();
    // auto P_p = r_P() + UU * t_min_P();
    auto r_gmean_ij =
        build<gmean_ij>(elemDivSafe(gmean_ij_tot, r_P(), t_min_P()));
    auto r_gtotal_var_ij = build<gtotal_var_ij>(
        r_gtotal_sqr_ij() - elemMult(r_gtotal_ij(), r_gmean_ij()));

    if constexpr (false) {
      auto test_elemMu = var::test_Derivative(
          [&U](const auto &r_gtotal_sqr_ij, const auto &r_gtotal_ij,
               const auto &r_gmean_ij) {
            return r_gtotal_sqr_ij() - elemMult(r_gtotal_ij(), r_gmean_ij());
          },
          1e-4, 1e-6, r_gtotal_sqr_ij, r_gtotal_ij, r_gmean_ij);
      if (!test_elemMu) {
        std::cerr << "\n error in test_elemMu!!\n" << test_elemMu.error()();
        std::abort();
      }
    }

    //    std::cerr<<"\n---------------------fin---------------------------------------------------\n";
    //    std::cerr<<"\n elemDiv(gmean_ij_tot,
    //    P_p)"<<primitive(elemDiv(gmean_ij_tot, P_p)); std::cerr<<"\n
    //    gmean_ij_tot"<<primitive(gmean_ij_tot); std::cerr<<"\n
    //    P_p"<<primitive(P_p); std::cerr<<"\n elemDiv(primitive(gmean_ij_tot),
    //    primitive(P_p))"<<elemDiv(primitive(gmean_ij_tot), primitive(P_p));

    //    std::cerr<<"\n------------------------------------------------------------------------\n";

    if constexpr (false) {
      assert(test_conductance_variance(
          r_gtotal_var_ij(), get<Conductance_variance_error_tolerance>(m)));
    }

    /* truncate is not derivative safe yet*/
    if constexpr (true) {
      r_gtotal_var_ij() =
          truncate_negative_variance(std::move(r_gtotal_var_ij()));
    }

    auto gvar_ij_tot = r_gtotal_var_ij(); // + gvar_ij_p * t_min_P();
    auto r_gvar_ij = build<gvar_ij>(elemDivSafe(gvar_ij_tot, r_P(), t_min_P()));
    Matrix<double> u(N, 1, 1.0);
    auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
    auto r_gsqr_i = build<gsqr_i>(r_gtotal_sqr_ij() * u);
    auto r_gvar_i = build<gvar_i>(r_gtotal_var_ij() * u);
    if constexpr (true) {
      r_gvar_i() = truncate_negative_variance(std::move(r_gvar_i()));
    }

    return build<Qdt>(ns, min_P(t_min_P), std::move(r_P), std::move(r_gmean_i),
                      std::move(r_gtotal_ij), std::move(r_gmean_ij),
                      std::move(r_gtotal_sqr_ij), std::move(r_gsqr_i),
                      std::move(r_gvar_i), std::move(r_gtotal_var_ij),
                      std::move(r_gvar_ij));
  }
  
  
  
  template <class C_Patch_Model, class C_Qx>
      requires(/*U<C_Patch_Model, Patch_Model> && */ U<C_Qx, Qx>)
  auto calc_Qdt_taylor(const C_Patch_Model &m, const C_Qx &t_Qx,
                     number_of_samples ns, double dt, std::size_t order) {
      auto v_Qrun=t_Qx()*dt;
      double max=maxAbs(primitive(v_Qrun));
      double desired=0.125;
      int k=std::ceil(std::log2(max/desired));
      std::size_t n=std::max(0,k);
      double scale=std::pow(2,-n);
      auto t_Qrun_sub=v_Qrun*scale;
      auto P_sub=build<P>(expm_taylor(t_Qrun_sub,order));
      auto r_Qn=get_Qn( P_sub,get<g>(m),ns,get<min_P>(m));
      for (std::size_t i=0; i<n; ++n)
      {
          r_Qn=sum_Qn(std::move(r_Qn),r_Qn);
      }
      get<number_of_samples>(r_Qn)()=ns;
      return Qn_to_Qdt(r_Qn);    
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
      requires(U<C_P, P>&&U<C_g,g>)
  auto get_Qn(const C_P &t_P, C_g const & t_g,number_of_samples n, min_P t_minP) {
      
      auto N=t_P().nrows();
      auto t_g2=apply([](auto x){return x*x;}, t_g());
      auto u=Matrix<double>(1,N,1.0);
      auto G=t_g()*u;
      auto GT=tr(G);
      auto G2=t_g2*u;
      auto G2T=tr(G2);
      auto Gmean=0.5*G+0.5*GT;
      
      auto Gvar=0.5*G2+0.5*G2T-Gmean;
      
      return build<Qn>(n, t_minP, t_P,
                       build<PG_n>(t_P() *Gmean* n()),
                       build<PGG_n>(t_P()* Gvar * (n() * n() * 0.5)));
  }
  
  template <class C_Qn>
      requires(U<C_Qn, Qn>)
  static C_Qn sum_Qn(C_Qn &&one, const C_Qn &two) {
      auto n1 = get<number_of_samples>(two)();
      get<PGG_n>(one)() =
          (get<PGG_n>(one)() * get<P>(two)()) +
          (get<PG_n>(one)() * get<PG_n>(two)()) +
          (get<P>(one)() * get<PGG_n>(two)());
      get<PG_n>(one)() = (get<PG_n>(one)() * get<P>(two)()) +
                         (get<P>(one)() * get<PG_n>(two)());
      get<P>(one) =
          normalize(build<P>(get<P>(one)() * get<P>(two)()), get<min_P>(two)());
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
        normalize(build<P>(get<P>(one)() * get<P>(two)()), get<min_P>(two)());
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

  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_ATP_step(FunctionTable && f,const C_Patch_Model &m, const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    auto dt = get<number_of_samples>(t_step)() / fs;
    auto t_Qx = f.fstop(Calc_eigen{},m, get<ATP_concentration>(t_step));

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

    if (!t_Qx)
      return t_Qx.error();
    else {
      return calc_Qdt(f,m, t_Qx.value(), get<number_of_samples>(t_step), dt);
    }
  }
  
  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt(FunctionTable && f,const C_Patch_Model &m, const ATP_step &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
      if constexpr(std::is_same_v<Nothing,decltype(f[Calc_Qdt_step{}])>)
          return calc_Qdt_ATP_step(f,m,t_step,fs);
      else
          return f.f(Calc_Qdt_step{},m,t_step,fs);
  }
  
  
  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qn_bisection(FunctionTable && f,const C_Patch_Model &m, const ATP_step &t_step, double fs, std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qn>> {
      auto dt = get<number_of_samples>(t_step)() / fs;
      auto ns =get<number_of_samples>(t_step) ;
      auto tQx=f.fstop(Calc_Qx{},m, get<ATP_concentration>(t_step));
      auto t_Qx = f.fstop(Calc_eigen{},tQx);
      
      if (!t_Qx)
          return t_Qx.error();
      else {
          double scale=std::pow(2.0,-1.0*order);
          
          number_of_samples n_ss(ns()*scale);
          double sdt=dt*scale;
          auto t_Psub=calc_P(m,t_Qx.value(),sdt,get<min_P>(m)());
          auto r_Qn=get_Qn( t_Psub,get<g>(m),n_ss,get<min_P>(m));
          for (std::size_t i=0; i<order; ++i)
          {
              r_Qn=sum_Qn(std::move(r_Qn),r_Qn);
          }
          assert(get<number_of_samples>(r_Qn)()==ns());
          return r_Qn;    
      }
  }
  
  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_bisection(FunctionTable && f,const C_Patch_Model &m, const ATP_step &t_step, double fs, std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
      auto maybe_Qn=calc_Qn_bisection(f,m,t_step,fs,order);
      if (!maybe_Qn) return maybe_Qn.error();
      else
      {
          return Qn_to_Qdt(maybe_Qn.value());
      }
  }
  
  
  
  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model> )
  auto calc_Qdt(FunctionTable&&f,const C_Patch_Model &m, const std::vector<ATP_step> &t_step,
                double fs) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    if (t_step.empty())
      return error_message("Emtpy ATP step");
    else {
      auto v_Qdt0 = calc_Qdt(f,m, t_step[0], fs);
      if (!v_Qdt0)
        return v_Qdt0.error();
      else {
        auto v_Qrun = get_Qn(v_Qdt0.value());
        for (std::size_t i = 1; i < t_step.size(); ++i) {
          auto v_Qdti = calc_Qdt(f,m, t_step[i], fs);
          if (!v_Qdti)
            return v_Qdti.error();
          else
            v_Qrun = sum_Qdt(std::move(v_Qrun), v_Qdti.value());
        }
        return Qn_to_Qdt(v_Qrun);
      }
    }
  }
  
  
  template <class FunctionTable,class C_Patch_Model>
  // requires(U<C_Patch_Model, Patch_Model> )
  auto calc_Qdt_bisection(FunctionTable&&f,const C_Patch_Model &m, const std::vector<ATP_step> &t_step,
                          double fs, std::size_t order) -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
      if (t_step.empty())
          return error_message("Emtpy ATP step");
      else {
          auto v_Qn0 = calc_Qn_bisection(f,m, t_step[0], fs,order);
          if (!v_Qn0)
              return v_Qn0.error();
          else {
              auto v_Qrun = v_Qn0.value();
              for (std::size_t i = 1; i < t_step.size(); ++i) {
                  auto v_Qni = calc_Qn_bisection(f,m, t_step[i], fs,order);
                  if (!v_Qni)
                      return v_Qni.error();
                  else
                      v_Qrun = sum_Qn(std::move(v_Qrun), v_Qni.value());
              }
              return Qn_to_Qdt(v_Qrun);
          }
      }
  }
  
  
  
  
  template <class FunctionTable,class C_Patch_Model>
  //   requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt(FunctionTable&&f,const C_Patch_Model &m, const ATP_evolution &t_step, double fs)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
    return std::visit([this, &m, &f,fs](auto &&a) { return calc_Qdt(f,m, a, fs); },
                      t_step());
  }
  
  template <class FunctionTable,class C_Patch_Model>
  //   requires(U<C_Patch_Model, Patch_Model>)
  auto calc_Qdt_bisection(FunctionTable&&f,const C_Patch_Model &m, const ATP_evolution &t_step, double fs, std::size_t order)
      -> Maybe_error<Transfer_Op_to<C_Patch_Model, Qdt>> {
      return std::visit([this, &m, &f,fs,order](auto &&a) { return calc_Qdt_bisection(f,m, a, fs,order  ); },
                        t_step());
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
averaging, uses_variance_aproximation variance> auto run_old(const Patch_State
&t_prior, Qdt const &t_Qdt, Patch_Model const &m, const Experiment_step &p,
double fs) const { auto &p_y = get<Patch_current>(p); auto &p_P_mean =
get<P_mean>(t_prior); auto &p_P_Cov = get<P_Cov>(t_prior);

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
  //                                 const Patch_current &p_y, double fs) const
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
  //            //                                               y_mean, y_var,
  //            plogL,
  //            //                                               eplogL,
  //            // vplogL,Q_dt.min_P(), e));
  //            //      else
  //            //        return Op(false, "fails at intertrace prediction!!: "
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
  //                         get<gtotal_var_ij>(t_Qdt)(), get<gvar_ij>(t_Qdt)())
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
  //        auto v_y_var = y_var(std::max(e_mu + N * gSg - N * zeta * sqr(sSg),
  //        e)); auto dy = p_y.value() - v_y_mean();

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

  //        auto v_eplogL = eplogL(-0.5 * log(2 * std::numbers::pi * v_y_var())
  //        -
  //                               0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
  //        vplogL v_vplogL(0.5);
  //        // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
  //        std::cerr << get<number_of_samples>(p).value() << "\t" << v_P_mean
  //        << "\n"; return Patch_State(logL(get<logL>(t_prior)() + v_plogL()),
  //                           elogL(get<elogL>(t_prior)() + v_eplogL()),
  //                           vlogL(get<vlogL>(t_prior)() + v_vplogL()),
  //                           v_P_mean, v_P_cov, v_y_mean, v_y_var, v_plogL,
  //                           v_eplogL, v_vplogL);
  //    }

  template <uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
           class FunctionTable,
           class C_Patch_State,
            class C_Qdt, class C_Patch_Model, class C_double>
    requires(
        /*(U<std::decay_t<C_Patch_State>,
           Patch_State>||U<std::decay_t<C_Patch_State>,
           Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
        U<C_double, double> &&
        U<C_Qdt, Qdt>)

  Maybe_error<C_Patch_State>
  Macror(FunctionTable&,C_Patch_State &&t_prior, C_Qdt const &t_Qdt, C_Patch_Model const &m,
         C_double const &Nch, const Patch_current &p_y, double fs)const  {

    using Transf = transformation_type_t<C_Qdt>;
    auto &p_P_cov = get<P_Cov>(t_prior);
    auto &p_P_mean = get<P_mean>(t_prior);
    //    auto &y = get<Patch_current>(p).value();
    auto &y = p_y.value();

    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto e = get<Current_Noise>(m).value() * fs /
             get<number_of_samples>(t_Qdt).value();
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
      auto r_P_mean = build<P_mean>(p_P_mean() * t_P());
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
              normalize(std::move(r_P_mean), t_min_P()),
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
            logL(0.0), elogL(0.0), vlogL(0.0),
            normalize(std::move(r_P_mean), t_min_P()),
            normalize(std::move(r_P_cov), t_min_P()), y_mean(NaN), y_var(NaN),
            plogL(NaN), eplogL(NaN), vplogL(NaN));
      else
        return Transfer_Op_to<C_Patch_Model, Patch_State_and_Evolution>(
            logL(0.0), elogL(0.0), vlogL(0.0),
            normalize(std::move(r_P_mean), t_min_P()),
            normalize(std::move(r_P_cov), t_min_P()), y_mean(NaN), y_var(NaN),
            plogL(NaN), eplogL(NaN), vplogL(NaN), Patch_State_Evolution());

    } else {
      return error_message("fails at init: " + r_test.error()());
    }
  }

  template <uses_adaptive_aproximation adaptive,
            uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance, return_predictions predictions,
            class FuncTable, class C_Parameters, class Model>
  auto log_Likelihood(FuncTable&& f,const Model &model, const C_Parameters &par,
                      const Experiment &e, const Recording &y)
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
    auto m = model(par);
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto ini = init<predictions>(m, get<initial_ATP_concentration>(e));

    auto gege = 0;
    if (!ini)
      return ini.error();
    else {
      auto run = fold(
          0ul, y().size(), std::move(ini).value(),
          [this, &f,&m, fs, &e, &y, &gege](C_Patch_State &&t_prior,
                                        std::size_t i_step) {
            ATP_evolution const &t_step =
                get<ATP_evolution>(get<Recording_conditions>(e)()[i_step]);
                
            auto time=get<Time>(get<Recording_conditions>(e)()[i_step])();
            auto time_segment=get<N_Ch_mean_time_segment_duration>(m)();
            auto Nchs = get<N_Ch_mean>(m)();
            std::size_t i_segment=std::floor(time/time_segment);
            auto j_segment= std::min(Nchs.size()-1,i_segment+1);
            auto r= std::max(1.0,time/time_segment-i_segment);
            auto Nch = Nchs[i_segment]*(1-r)+r*Nchs[j_segment];
            
            auto Maybe_t_Qdt = calc_Qdt(f,m, t_step, fs);
            if (!Maybe_t_Qdt)
              return Maybe_error<C_Patch_State>(Maybe_t_Qdt.error());
            else {
              auto t_Qdt = std::move(Maybe_t_Qdt.value());

              //
              if constexpr (false) {
                auto test_der_t_Qdt = var::test_Derivative(
                    [this, &t_step, &fs, &gege,&f](auto const &l_m,
                                                auto const &l_Qx) {
                      return f.f(Calc_Qdt{},l_m, t_step, fs);
                    },
                    1e-6, 1e-2, m);
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
                    [this, &t_step, &fs,&f](auto const &l_m, auto const &l_prior,
                                         auto const &l_Qdt) {
                      return f.ff(MacroR<uses_recursive_aproximation(false),
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
              if constexpr (!adaptive.value)
              {
                  return f.f(MacroR<recursive, averaging, variance>{},
                    std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
              }
              else {
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
                   // using egsr=typename decltype(f.f(MacroR<recursive, averaging, variance>{}))::ege;
                  //  auto r=egsr();
                   return f.f(MacroR<recursive, averaging, variance>{},
                      std::move(t_prior), t_Qdt, m, Nch, y()[i_step], fs);
                } else {
                  return f.f(MacroR<uses_recursive_aproximation(false), averaging,
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
        return build<Vector_Space<logL, elogL, vlogL>>(get<logL>(run.value()),
                                                       get<elogL>(run.value()),
                                                       get<vlogL>(run.value()));
    }
  }

  template <class Patch_Model>
  Maybe_error<Simulated_Sub_Step>
  sub_sub_sample(std::mt19937_64 &mt, Simulated_Sub_Step &&t_sim_step,
                 const Patch_Model &m, const ATP_step &t_s, std::size_t n_sub,
                 double fs) {
    auto &t_g = get<g>(m);
    auto &N = get<N_channel_state>(t_sim_step);
    double &ysum = get<y_sum>(t_sim_step)();
    auto &sum_samples = get<number_of_samples>(t_sim_step)();

    auto n_samples = get<number_of_samples>(t_s)();
    auto tQx=calc_Qx(m, get<ATP_concentration>(t_s));
    
    auto dt = n_samples / fs;
    auto sub_dt = dt / n_sub;

    double sub_sample = 1.0 * n_samples / n_sub;
    auto t_P = calc_P(m, tQx, sub_dt, get<min_P>(m)());
    for (std::size_t i = 0; i < n_sub; ++i) {
      N = sample_Multinomial(mt, t_P, N);
      ysum += getvalue(N() * t_g()) * sub_sample;
    }
    sum_samples += n_samples;
    return t_sim_step;
  }

  template <class Patch_Model>
  Maybe_error<Simulated_Sub_Step>
  sub_sub_sample(std::mt19937_64 &mt, Simulated_Sub_Step &&t_sim_step,
                 const Patch_Model &m, const std::vector<ATP_step> &t_s,
                 std::size_t n_sub, double fs) {

    for (std::size_t i = 0; i < t_s.size(); ++i) {
      auto Maybe_sub_step =
          sub_sub_sample(mt, std::move(t_sim_step), m, t_s[i], n_sub, fs);
      if (!Maybe_sub_step)
        return Maybe_sub_step.error();
      else
        t_sim_step = std::move(Maybe_sub_step.value());
    }
    return t_sim_step;
  }

  template <class Patch_Model>
  Maybe_error<Simulated_Step>
  sub_sample(std::mt19937_64 &mt, Simulated_Step &&t_sim_step,
             const Patch_Model &m, const ATP_evolution &t_s, std::size_t n_sub,
             double fs) {
    auto &N = get<N_channel_state>(t_sim_step);

    // std::cerr<<N();

    auto t_sub_step = Simulated_Sub_Step(get<N_channel_state>(t_sim_step),
                                         number_of_samples(0ul), y_sum(0.0));

    auto Maybe_t_sub_step = std::visit(
        [this, &mt, &m, n_sub, &t_sub_step, fs](auto const &a) {
          return sub_sub_sample(mt, std::move(t_sub_step), m, a, n_sub, fs);
        },
        t_s());

    if (!Maybe_t_sub_step)
      return Maybe_t_sub_step.error();
    else {
      t_sub_step = std::move(Maybe_t_sub_step.value());
      double y_mean =
          get<y_sum>(t_sub_step)() / get<number_of_samples>(t_sub_step)();
      get<N_channel_state>(t_sim_step) = get<N_channel_state>(t_sub_step);
      auto &t_e_step = get<Simulated_Recording>(t_sim_step);
      double e =
          get<Current_Noise>(m)() * fs / get<number_of_samples>(t_sub_step)();
      auto y_baseline = get<Current_Baseline>(m);

      t_e_step()().emplace_back(
          Patch_current(y_mean + y_baseline() +
                        std::normal_distribution<double>()(mt) * std::sqrt(e)));
      return t_sim_step;
    }
  }

  template <class Patch_Model>
  Simulated_Step init_sim(std::mt19937_64 &mt, const Patch_Model &m,
                          const Experiment &e) {
    auto initial_x = get<initial_ATP_concentration>(e);
    auto v_Qx = calc_Qx(m, initial_x());
    auto r_P_mean = P_mean(get<P_initial>(m)());
    auto N = get<N_Ch_mean>(m)()[0];
    auto sim = Simulated_Recording(Recording{});
    auto N_state = sample_Multinomial(mt, r_P_mean, N);
    return Simulated_Step(std::move(N_state), std::move(sim));
  }

  static Simulated_Recording copy_NaNs(Simulated_Recording &&sim,
                                       const Recording &r) {
    for (std::size_t i = 0; i < size(r()); ++i)
      if (std::isnan(r()[i]()))
        sim()()[i]() = r()[i]();

    return std::move(sim);
  }

  template <class Model, class Id>
  Maybe_error<Simulated_Recording>
  sample(std::mt19937_64 &mt, const Model &model, const Parameters<Id> &par,
         const Experiment &e, const Simulation_Parameters &sim,
         const Recording &r = Recording{}) {

    auto m = model(par);
    auto n_sub = get<Number_of_simulation_sub_steps>(sim);
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto sim_recording = Recording{};

    auto ini = init_sim(mt, m, e);
      auto run =
          fold(get<Recording_conditions>(e)(), ini,
               [this, &m, fs, n_sub, &mt](Simulated_Step &&t_sim_step,
                                          Experiment_step const &t_step) {
                 return Maybe_error<Simulated_Step>(sub_sample(
                     mt, std::move(t_sim_step), m, t_step, n_sub(), fs));
               });
      if (!run)
        return run.error();
      else {
        return copy_NaNs(std::move(get<Simulated_Recording>(run.value())), r);
      }
    
  }
};

struct Model0 : public Model_Patch<Model0> {};
struct Model1 : public Model_Patch<Model1> {};

struct Allost1 : public Model_Patch<Allost1> {};

template <uses_adaptive_aproximation adaptive,
          uses_recursive_aproximation recursive,
          uses_averaging_aproximation averaging,
          uses_variance_aproximation variance, class Model>
struct Likelihood_Model {
  Model m;
  Number_of_simulation_sub_steps number_of_simulation_sub_steps;
  Likelihood_Model(const Model &model, Number_of_simulation_sub_steps n)
      : m{model}, number_of_simulation_sub_steps{n} {}
};

template <uses_adaptive_aproximation adaptive,
          uses_recursive_aproximation recursive,
          uses_averaging_aproximation averaging,
          uses_variance_aproximation variance, class Model>
auto make_Likelihood_Model(const Model &m, Number_of_simulation_sub_steps n) {
  return Likelihood_Model<adaptive, recursive, averaging, variance, Model>(m,
                                                                           n);
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    class FuncTable,class Model, class Parameters, class Variables, class DataType>
Maybe_error<double>
logLikelihood(FuncTable&& f,const Likelihood_Model<adaptive, recursive, averaging, variance,
                                     Model> &lik,
              Parameters const &p, const Variables &var, const DataType &y) {
  auto v_logL =
      Macro_DMR{}
          .log_Likelihood<adaptive, recursive, averaging, variance,
                          return_predictions(false)>(f,lik.m, p, y, var);
  if (!v_logL)
    return v_logL.error();
  else
    return get<logL>(v_logL.value())();
}

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    class FunctionTable,class Model, class Parameters, class Variables, class DataType>
Maybe_error<Patch_State_Evolution> logLikelihoodPredictions(FunctionTable&& f,
    const Likelihood_Model<adaptive, recursive, averaging, variance, Model>
        &lik,
    Parameters const &p, const Variables &var, const DataType &y) {
  return Macro_DMR{}
      .log_Likelihood<adaptive, recursive, averaging, variance,
                      return_predictions(true)>(f,lik.m, p, y, var);
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
  return std::visit([](auto const &a) { return get_num_samples(a); }, e());
}

inline ATP_evolution average_ATP_step(ATP_step const &x) { return x; }

inline ATP_evolution average_ATP_step(std::vector<ATP_step> const &x) {
  if (x.empty())
    return x;
  else if (x.size() == 1)
    return x[0];
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
    if (out.size() == 1)
      return std::move(out[0]);
    else
      return out;
  }
}

// static ATP_evolution average_ATP_step(ATP_evolution const & x){
//     return std::visit([] (auto const& a){return average_ATP_step(a);},x());}

inline ATP_evolution average_ATP_step(ATP_evolution const &x) { return x; }

template <
    uses_adaptive_aproximation adaptive, uses_recursive_aproximation recursive,
    uses_averaging_aproximation averaging, uses_variance_aproximation variance,
    class Model, class Parameters, class Variables>
auto simulate(std::mt19937_64 &mt,
              const Likelihood_Model<adaptive, recursive, averaging, variance,
                                     Model> &lik,
              Parameters const &p, const Variables &var) {
  return Macro_DMR{}
      .sample(mt, lik.m, p, var, lik.number_of_simulation_sub_steps)
      .value()();
}

static Experiment_step average_Experimental_step(Experiment_step const &x) {
  return Experiment_step(get<Time>(x), average_ATP_step(get<ATP_evolution>(x)));
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
    return e;
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
    return x;
  } else {
    return std::vector<ATP_step>{std::move(x), std::move(e)};
  }
}

static ATP_evolution add_ATP_step(ATP_evolution &&x, ATP_evolution &&e) {
  return std::visit(
      [&x](auto &&a) {
        return std::visit(
            [&a](auto &&ax) {
              return add_ATP_step_i(std::move(ax), std::move(a));
            },
            x());
      },
      e());
}

class experiment_fractioner {
  std::vector<std::size_t> segments = {73, 33, 22, 22, 4};
  std::size_t min_number_of_samples = 10;

public:
  experiment_fractioner(const std::vector<std::size_t> &t_segments,
                        std::size_t t_min_number_of_samples)
      : segments{t_segments}, min_number_of_samples{t_min_number_of_samples} {}

  static auto average_Recording(const Recording_conditions &e,
                                const Recording &y,
                                const std::vector<std::size_t> &indexes0,
                                const std::vector<std::size_t> &indexes1) {

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
          out_x[ii] = average_Experimental_step(e()[i]);
        } else {

          auto n_samples = get_num_samples(e()[i]);
          sum_y += y()[i]() * n_samples;
          sum_samples += n_samples;
          v_ATP = add_ATP_step(std::move(v_ATP), average_ATP_step(e()[i]));

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
        v_ATP = add_ATP_step(std::move(v_ATP), average_ATP_step(e()[i]));
      }
    }
    std::cerr << "\nout_x\n****************************************************"
                 "***********************\n";
    std::cerr << out_x;
    std::cerr << "\nout_y\n****************************************************"
                 "***********************\n";
    std::cerr << out_y;

    return std::tuple(Recording_conditions(out_x), Recording(out_y));
  }

  auto operator()(const Recording &y, const Experiment &x, std::mt19937_64 &mt,
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
    std::cerr << "\nindexes\n**************************************************"
                 "*************************\n";
    std::cerr << indexes;
    // std::abort();
    auto n_frac = size(indexes);
    by_fraction<Recording> y_out(n_frac);
    by_fraction<Experiment> x_out(
        n_frac,
        Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                   get<initial_ATP_concentration>(x)));
    y_out[n_frac - 1] = y;
    x_out[n_frac - 1] = x;

    for (std::size_t i = n_frac - 1; i > 0; --i) {
      std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
          average_Recording(get<Recording_conditions>(x_out[i]), y_out[i],
                            indexes[i], indexes[i - 1]);
    }

    // std::abort();
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

template <class Id, class... Ts>
void report_title(save_Predictions<Parameters<Id>> &s,
                  cuevi_mcmc<Parameters<Id>> const &, const Ts &...t) {

  s.f << "n_fractions" << s.sep << "n_betas" << s.sep << "iter" << s.sep
      << "nsamples" << s.sep << "beta" << s.sep << "i_walker" << s.sep
      << "id_walker" << s.sep << "i_x" << s.sep << "time" << s.sep
      << "num_samples" << s.sep << "ATP" << s.sep << "ATPevol" << s.sep
      << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_std" << s.sep << "plogL"
      << s.sep << "pelogL"
      << "\n";
}

template <class Id>
void report_title(save_Predictions<Parameters<Id>> &s,
                  thermo_mcmc<Parameters<Id>> const &, ...) {

  s.f << "n_betas" << s.sep << "iter" << s.sep << "beta" << s.sep << "i_walker"
      << s.sep << "id_walker" << s.sep << "i_x" << s.sep << "time" << s.sep
      << "num_samples" << s.sep << "ATP" << s.sep << "ATP_evolution" << s.sep
      << "Y_obs" << s.sep << "Y_pred" << s.sep << "Y_std" << s.sep << "plogL"
      << s.sep << "pelogL"
      << "\n";
}

void report_title(save_Predictions<Matrix<double>> &s,
                  thermo_mcmc<Matrix<double>> const &, ...) {}

template <class Id,class FunctionTable>
void report(FunctionTable&& f,std::size_t iter, save_Predictions<Parameters<Id>> &s,
            thermo_mcmc<Parameters<Id>> const &data, ...) {
  if (iter % s.save_every == 0)
    for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
      for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker)
        for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)

          s.f << num_betas(data) << s.sep << iter << s.sep << data.beta[i_beta]
              << s.sep << i_walker << s.sep << data.i_walkers[i_walker][i_beta]
              << s.sep << i_par << s.sep
              << data.walkers[i_walker][i_beta].parameter[i_par] << "\n";
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

inline std::string ToString(const ATP_evolution &ev) {
  return std::visit([](auto const &a) { return ToString(a); }, ev());
}

template <class FunctionTable,class Prior, class Likelihood, class Variables, class DataType,
          class Parameters>
void report(FunctionTable &&f, std::size_t iter, save_Predictions<Parameters> &s,
            cuevi_mcmc<Parameters> const &data, Prior const &prior,
            Likelihood const &lik, const DataType &ys, const Variables &xs) {
  if (iter % s.save_every == 0)
    for (std::size_t i_frac = 0; i_frac < size(data.beta); ++i_frac)
      for (std::size_t i_beta = 0; i_beta < size(data.beta[i_frac]); ++i_beta)
        for (std::size_t i_walker = 0; i_walker < size(data.walkers);
             ++i_walker) {
          Maybe_error<Patch_State_Evolution> prediction =
                logLikelihoodPredictions(f.fork( var::I_thread(i_walker)),
                  lik, data.walkers[i_walker][i_frac][i_beta].parameter,
                  ys[i_frac], xs[i_frac]);
          if (is_valid(prediction)) {
            auto &predictions = prediction.value();
            for (std::size_t i_x = 0; i_x < size(ys[i_frac]); ++i_x) {
              auto v_ev = get<ATP_evolution>(
                  get<Recording_conditions>(xs[i_frac])()[i_x]);

              s.f << size(data.beta) << s.sep << size(data.beta[i_frac])
                  << s.sep << iter << s.sep << data.nsamples[i_frac] << s.sep
                  << data.beta[i_frac][i_beta] << s.sep << i_walker << s.sep
                  << data.i_walkers[i_walker][i_frac][i_beta] << s.sep << i_x
                  << s.sep
                  << get<Time>(get<Recording_conditions>(xs[i_frac])()[i_x])
                  << s.sep << get_num_samples(v_ev) << s.sep
                  << ToString(average_ATP_step(v_ev)) << s.sep << ToString(v_ev)
                  << s.sep << ys[i_frac]()[i_x]() << s.sep
                  << get<y_mean>(predictions()[i_x]) << s.sep
                  << get<y_var>(predictions()[i_x]) << s.sep
                  << get<plogL>(predictions()[i_x]) << s.sep
                  << get<eplogL>(predictions()[i_x]) << "\n";
            }
          }
        }
}

template <class Id>
auto cuevi_Model_by_convergence(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments,
    std::size_t t_min_number_of_samples,

    std::size_t num_scouts_per_ensemble, double min_fraction,
    std::size_t thermo_jumps_every, std::size_t max_iter, double max_ratio,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double stops_at, bool includes_zero, std::size_t initseed) {
  return cuevi_integration(
      checks_derivative_var_ratio<cuevi_mcmc, Parameters<Id>>(max_iter,
                                                              max_ratio),
      experiment_fractioner(t_segments, t_min_number_of_samples),
      save_mcmc<Parameters<Id>, save_likelihood<Parameters<Id>>,
                save_Parameter<Parameters<Id>>, save_Evidence,
                save_Predictions<Parameters<Id>>>(path, filename, 100ul, 100ul,
                                                  100ul, 100ul),
      num_scouts_per_ensemble, min_fraction, thermo_jumps_every,
      n_points_per_decade_beta, n_points_per_decade_fraction, stops_at,
      includes_zero, initseed);
}

template <class Id>
auto cuevi_Model_by_iteration(
    std::string path, std::string filename,
    const std::vector<std::size_t> &t_segments,
    std::size_t t_min_number_of_samples,

    std::size_t num_scouts_per_ensemble,
    std::size_t max_number_of_simultaneous_temperatures,
    double min_fraction,
    std::size_t thermo_jumps_every,
    std::size_t max_iter_warming,
    std::size_t max_iter_equilibrium,
    double max_ratio,
    double n_points_per_decade_beta, double n_points_per_decade_fraction,
    double stops_at, bool includes_zero, std::size_t initseed) {
  return cuevi_integration(
      less_than_max_iteration(max_iter_warming, max_iter_equilibrium),
      experiment_fractioner(t_segments, t_min_number_of_samples),
      save_mcmc<Parameters<Id>, save_likelihood<Parameters<Id>>,
                save_Parameter<Parameters<Id>>, save_Evidence,
                save_Predictions<Parameters<Id>>>(path, filename, 10ul, 100ul,
                                                  10ul, 100ul),
      num_scouts_per_ensemble, max_number_of_simultaneous_temperatures,min_fraction, thermo_jumps_every,
      n_points_per_decade_beta, n_points_per_decade_fraction, stops_at,
      includes_zero, initseed);
}

template <class Id>
auto thermo_Model_by_max_iter(std::string path, std::string filename,
                              std::size_t num_scouts_per_ensemble,
                              std::size_t max_num_simultaneous_temperatures,
                              std::size_t thermo_jumps_every,
                              std::size_t max_iter_warming,
                              std::size_t max_iter_equilibrium,
                              double n_points_per_decade,
                              double stops_at, bool includes_zero,
                              std::size_t initseed) {
  return thermodynamic_integration(
      less_than_max_iteration(max_iter_warming,max_iter_equilibrium),
      save_mcmc<Parameters<Id>, save_likelihood<Parameters<Id>>,
                save_Parameter<Parameters<Id>>, save_Evidence,
                save_Predictions<Parameters<Id>>>(path, filename, 10ul, 10ul,
                                                  10ul, 100ul),
      num_scouts_per_ensemble, max_num_simultaneous_temperatures,
      thermo_jumps_every, n_points_per_decade, stops_at, includes_zero,
      initseed);
}

} // namespace macrodr

#endif // QMODEL_H
