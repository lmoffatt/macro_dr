#pragma once
#include "experiment.h"
#include "fold.h"
#include "matrix.h"
#include <functional>
#include <numeric>
#include <random>
#include <set>
#ifndef QMODEL_H
#define QMODEL_H
#include <map>
#include <string>

#include "maybe_error.h"
#include "variables.h"
#include "parameters.h"
#include "derivative_operator.h"
namespace macrodr {

using var::Power;
using var::Product;
using var::Var;
using var::Vector_Space;
using var::Parameters;

using var::Constant;
using var::build;
using var::Transfer_Op_to;
using var::transformation_type_t;
using var::Op_t;
using var::primitive;

using var::U;

using std::exp;
using std::max;
using std::abs;

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

class Qa : public var::Var<Qa, Matrix<double>> {};

class Qx : public var::Var<Qx, Matrix<double>> {};

class g : public var::Var<g, Matrix<double>> {};

class N_St : public var::Constant<N_St, std::size_t> {};

class N_Ch_mean : public var::Var<N_Ch_mean, double> {};

class min_P : public var::Constant<min_P, double> {};

class N_Ch_std : public var::Var<N_Ch_std, double> {};

class curr_noise : public var::Var<curr_noise, double> {};

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
class uses_recursive_aproximation
    : public var::struct_Var<uses_recursive_aproximation, bool> {};
class uses_averaging_aproximation
    : public var::struct_Var<uses_averaging_aproximation, int> {};

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

using Qx_eig = Vector_Space<Qx, V, lambda, W>;

using Qdt = Vector_Space<P, gmean_i, gtotal_ij, gmean_ij, gtotal_sqr_ij, gsqr_i,
                         gvar_i, gtotal_var_ij, gvar_ij>;

using Patch_Model = Vector_Space<N_St, Q0, Qa, g, N_Ch_mean,
    curr_noise, min_P, Probability_error_tolerance,
                                 Conductance_variance_error_tolerance>;


template<class Id>
struct Model_Patch
{
    template<class F>
    class Model{
    F m_f;
public:
    static constexpr bool is_Model_Patch=true;
    Model(F t_f):m_f{t_f}{}
    
    
    template<class P>
     requires std::is_same_v<var::untransformed_type_t<P>,Parameters<Id>>
    auto operator()(const P& t_p)const{
        return std::invoke(m_f,t_p);
    }
};
    template<class F>
    Model(F f)->Model_Patch<Id>::Model<F>;
    
};




using Patch_State = Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean,
                                 y_var, plogL, eplogL, vplogL>;

class Number_of_simulation_sub_steps
    : public Var<Number_of_simulation_sub_steps, std::size_t> {};

class Simulated_Experiment : public Var<Simulated_Experiment, Experiment> {};

using Simulated_Step = Vector_Space<N_channel_state, Simulated_Experiment>;

using Simulated_Sub_Step = Vector_Space<N_channel_state, y_sum, t_sum>;

using Simulation_Parameters = Vector_Space<Number_of_simulation_sub_steps>;

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
  
  template<class C_double>
      requires U<C_double,double>
  static C_double Ee(C_double const& x, C_double const& y, C_double const& exp_x, C_double const& exp_y,
                   double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(primitive(x) - primitive(y)) < eps)
      return exp_x;
    else
      return (exp_x - exp_y) / (x - y);
  };
  
  template<class C_double>
      requires U<C_double,double>
  static C_double EX_111(C_double const& x, C_double const& y, C_double const& z, C_double const& exp_x) {
    return exp_x / ((x - y) * (x - z));
  }
  
  template<class C_double>
      requires U<C_double,double>
  static C_double E111(C_double const& x, C_double const& y, C_double const& z, C_double const& exp_x, C_double const& exp_y,
                     C_double const& exp_z) {
    return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
           EX_111(z, y, x, exp_z);
  }
  template<class C_double>
      requires U<C_double,double>
  static C_double E12(C_double const& x, C_double const& y, C_double const& exp_x, C_double const& exp_y) {
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
  }
  
  template<class C_double>
      requires U<C_double,double>
  static C_double E3(C_double const& x, C_double const& y, C_double const& z, C_double const& exp_x, C_double const& exp_y,
                   C_double const& exp_z,
                   double eps = std::numeric_limits<double>::epsilon()) {
    auto x_=primitive(x);
    auto y_=primitive(y);
    auto z_=primitive(z);
    
    
    if (sqr(x_ - y_) < eps) // x==y
    {
      if (sqr(y_ - z_) < eps) // y==z
        return exp_x / 2.0; // x==y==z
      else
        return E12(z, x, exp_z, exp_x); // x==y!=z
    } else if (sqr(y_ - z_) < eps)        // x!=y==z
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

  template <bool output,class C_P_mean>
               requires (U<C_P_mean,P_mean>)
  static Maybe_error<bool> test(const C_P_mean &pp,
                                Probability_error_tolerance tolerance) {
    auto p=primitive(pp());
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
  
  template <bool output,class C_P_Cov>
      requires (U<C_P_Cov,P_Cov>)
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
  
  template <bool output,class C_P_mean,class C_P_Cov>
      requires (U<C_P_Cov,P_Cov>&&U<C_P_mean,P_mean>)
  static Maybe_error<bool> test(const C_P_mean &t_P_mean, const C_P_Cov &t_P_cov,
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
  template<class C_Patch_Model>
        requires U<C_Patch_Model,Patch_Model>
    auto calc_eigen(const C_Patch_Model &m, ATP_concentration x) ->Maybe_error<Transfer_Op_to<C_Patch_Model,Qx_eig>>
    {
    using Trans=transformation_type_t<C_Patch_Model>;
    auto v_Qx=build<Qx> (get<Q0>(m)() + get<Qa>(m)() * x.value());
    Matrix<double> u(v_Qx().ncols(), 1, 1.0);
    v_Qx() = v_Qx() - diag(v_Qx() * u);
    auto maybe_eig = eigs(v_Qx());
    if (maybe_eig) {
      auto [v_V, v_l, v_W] = maybe_eig.value();
      return build<Qx_eig>(std::move(v_Qx), build<V>(std::move(v_V)), build<lambda>(std::move(v_l)),
                           build<W>(std::move(v_W)));
    } else
      return maybe_eig.error();
  }
    
    template<class C_P_mean>
        requires U<C_P_mean,P_mean>
  static C_P_mean normalize(C_P_mean &&pp, double t_min_p) {
    auto& p=primitive(pp());
    for (std::size_t i = 0; i < p.nrows(); ++i) {
      double sum = 0;
      for (std::size_t j = 0; j < p.ncols(); ++j) {
        if (p(i, j) > 1 - t_min_p) {
          for (std::size_t k = 0; k < p.ncols(); ++k) {
            p(i, k) = (j == k) ? 1.0 : 0.0;
          }
          return std::move(pp);
        } else if (p(i, j) < t_min_p)
          p(i, j) = 0;
        else
          sum += p(i, j);
      }
      if (sum != 1)
        for (std::size_t j = 0; j < p.ncols(); ++j)
          p(i, j) = p(i, j) / sum;
    }
    return std::move(pp);
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
  
  template<class C_P_Cov>
      requires U<C_P_Cov,P_Cov>
    static C_P_Cov normalize(C_P_Cov &&pp, double t_min_p) {
    auto& p=primitive(pp());
    for (std::size_t i = 0; i < p.nrows(); ++i) {
      if (p(i, i) < t_min_p) {
        for (std::size_t j = 0; j < p.ncols(); ++j) {
          set(p,i, j, 0.0);
        }
      }
    }
    return std::move(pp);
  }
  
  
  
  static P normalize(P &&p, double t_min_p) {
    // std::cerr<<p;
    for (std::size_t i = 0; i < p().nrows(); ++i) {
      double sumP = 0;
      for (std::size_t j = 0; j < p().ncols(); ++j)
        if (p()(i, j) < t_min_p)
          p()(i, j) = 0;
        else
          sumP += p()(i, j);
      for (std::size_t j = 0; j < p().ncols(); ++j)
        p()(i, j) = p()(i, j) / sumP;
    }
    // std::cerr<<p;
    return std::move(p);
  }

  template <class Vs>
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
  
  template<class C_Patch_Model, class C_Qx_eig>
      requires (U<C_Patch_Model,Patch_Model>&&U<C_Qx_eig,Qx_eig>)
  auto calc_Peq(C_Qx_eig const &t_Qx, const C_Patch_Model &m)->Transfer_Op_to<C_Patch_Model,P_mean> {
    auto nstates = get<N_St>(m).value();
    auto p0 = Matrix<double>(1ul, nstates, 1.0 / nstates);
    
    auto &landa = get<lambda>(t_Qx)();
    auto &Vv = get<V>(t_Qx)();
    auto &Wv = get<W>(t_Qx)();
    auto laexp = DiagonalMatrix<double>(nstates, nstates, 0.0);
    for (std::size_t i = 0; i < nstates; ++i) {
      if (landa(i, i) == 0.0)
        laexp[i] = 1.0;
    }
    if (false) {
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

  auto calc_P(const Patch_Model &m, const Qx_eig &t_Qx, double dt,
              double t_min_P) {
    auto ladt = get<lambda>(t_Qx)() * dt;

    auto exp_ladt = apply([](double x) { return std::exp(x); }, ladt);
    return normalize(P(get<V>(t_Qx)() * exp_ladt * get<W>(t_Qx)()), t_min_P);
  }

  auto calc_Qdt_old(const Patch_Model &m, const Qx_eig &t_Qx, double dt) {

    auto t_min_P = get<min_P>(m)();
    auto &v_g = get<g>(m);

    std::size_t N = t_Qx[Var<Qx>{}]().ncols();

    auto ladt = t_Qx[Var<lambda>{}]() * dt;

    auto exp_ladt = apply([](double x) { return std::exp(x); }, ladt);
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

    return Qdt(std::move(v_P), std::move(v_gmean_i), std::move(v_gtotal_ij),
               std::move(v_gmean_ij), std::move(v_gtotal_sqr_ij),
               std::move(v_gsqr_i), std::move(v_gvar_i),
               std::move(v_gtotal_var_ij), std::move(v_gvar_ij));
  }
  
  template<class C_Patch_Model, class C_Qx_eig>
      requires (U<C_Patch_Model,Patch_Model>&& U<C_Qx_eig,Qx_eig>)
  auto calc_Qdt(const C_Patch_Model &m, const C_Qx_eig &t_Qx, double dt) {
    using Trans=transformation_type_t<C_Patch_Model>;
    // const double eps=std::numeric_limits<double>::epsilon();
    auto &t_V = get<V>(t_Qx);
    auto &t_landa = get<lambda>(t_Qx);
    auto &t_W = get<W>(t_Qx);
    auto &t_g = get<g>(m);
    auto t_min_P = get<min_P>(m);
    auto v_ladt = t_landa() * dt;
    auto v_exp_ladt = apply([](auto const & x) { using std::exp;
return exp(x); }, v_ladt);
    
    
    auto r_P = build<P>(t_V() * v_exp_ladt * t_W());

    std::size_t N = r_P().ncols();
    
    SymmetricMatrix<Op_t<Trans,double>> E2m(N, N);
    for (std::size_t i = 0; i < N; ++i)
    {
      
      for (std::size_t j = 0; j < i + 1; ++j)
      {
        set(E2m,
            i, j,
            Ee(v_ladt[i], v_ladt[j], v_exp_ladt[i], v_exp_ladt[j], t_min_P()));
        
      }
    }
    
    
    
    Matrix<Op_t<Trans,double>> WgV_E2(N, N);

    auto v_WgV = t_W() * diag(t_g()) * t_V();

    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = v_WgV(i, j) * E2m(i, j);
    
    
    
    auto r_gtotal_ij = build<gtotal_ij>(t_V() * WgV_E2 * t_W());
    
    
    Matrix<Op_t<Trans,double>> WgV_E3(N, N,Op_t<Trans,double>(0.0));
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
    //      std::cerr<<"\t"<<WgV_E3(n1, n3);
          
          WgV_E3(n1, n3) = WgV_E3(n1, n3)+
              v_WgV(n1, n2) * v_WgV(n2, n3) *
              E3(v_ladt[n1], v_ladt[n2], v_ladt[n3], v_exp_ladt[n1],
                 v_exp_ladt[n2], v_exp_ladt[n3], t_min_P()); // optimizable
        }
    
    auto r_gtotal_sqr_ij = build<gtotal_sqr_ij>(t_V() * WgV_E3 * t_W() * 2.0);
    
    
   // assert(test_conductance_variance(
   //     r_gtotal_sqr_ij(), get<Conductance_variance_error_tolerance>(m)));
    r_gtotal_sqr_ij() =
        truncate_negative_variance(std::move(r_gtotal_sqr_ij()));
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        if (r_P()(i, j) == 0) {
          r_gtotal_ij()(i, j) = 0;
          r_gtotal_sqr_ij()(i, j) = 0;
        }

    Matrix<double> U(1, t_g().size(), 1.0);
    Matrix<double> UU(t_g().size(), t_g().size(), 1.0);
    auto gmean_ij_p = X_plus_XT(t_g() * U) * (0.5);
    
    
    auto gvar_ij_p =
        apply([](auto x) { return abs(x); }, t_g() * U - tr(t_g() * U)) *
        (0.5);

    auto gmean_ij_tot = r_gtotal_ij() + gmean_ij_p * t_min_P();
    auto P_p = r_P() + UU * t_min_P();
    auto r_gmean_ij = build<gmean_ij>(elemDiv(gmean_ij_tot, P_p));
    auto r_gtotal_var_ij = build<gtotal_var_ij>(r_gtotal_sqr_ij() -
                                         elemMult(r_gtotal_ij(), r_gmean_ij()));
    
//    std::cerr<<"\n---------------------fin---------------------------------------------------\n";
//    std::cerr<<"\n elemDiv(gmean_ij_tot, P_p)"<<primitive(elemDiv(gmean_ij_tot, P_p));
//    std::cerr<<"\n gmean_ij_tot"<<primitive(gmean_ij_tot);
//    std::cerr<<"\n P_p"<<primitive(P_p);
//    std::cerr<<"\n elemDiv(primitive(gmean_ij_tot), primitive(P_p))"<<elemDiv(primitive(gmean_ij_tot), primitive(P_p));
    
//    std::cerr<<"\n------------------------------------------------------------------------\n";
    
    
    
    //assert(test_conductance_variance(
    //    r_gtotal_var_ij(), get<Conductance_variance_error_tolerance>(m)));
    r_gtotal_var_ij() =
        truncate_negative_variance(std::move(r_gtotal_var_ij()));

    auto gvar_ij_tot = r_gtotal_var_ij() + gvar_ij_p * t_min_P();
    auto r_gvar_ij = build<gvar_ij>(elemDiv(gvar_ij_tot, P_p));
    Matrix<double> u(N, 1, 1.0);
    auto r_gmean_i = build<gmean_i>(r_gtotal_ij() * u);
    auto r_gsqr_i = build<gsqr_i>(r_gtotal_sqr_ij() * u);
    auto r_gvar_i = build<gvar_i>(r_gtotal_var_ij() * u);
    
    return build<Qdt>(std::move(r_P), std::move(r_gmean_i), std::move(r_gtotal_ij),
               std::move(r_gmean_ij), std::move(r_gtotal_sqr_ij),
               std::move(r_gsqr_i), std::move(r_gvar_i),
               std::move(r_gtotal_var_ij), std::move(r_gvar_ij));
  }
  
  
  template<class C_Matrix>
      requires U<C_Matrix,Matrix<double>>
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
  
  template<class C_Qdt>
      requires U<C_Qdt,Qdt>
  Maybe_error<bool> test(const C_Qdt &q,
                         Conductance_variance_error_tolerance tol) {
    return "fails Qdt test" >>
           test_conductance_variances<gmean_i, gtotal_ij, gmean_ij,
                                      gtotal_sqr_ij, gsqr_i, gvar_i,
                                      gtotal_var_ij, gvar_ij>(q, tol);
  }
  
  template<class C_Matrix>
      requires U<C_Matrix,Matrix<double>>
  
  C_Matrix truncate_negative_variance(C_Matrix &&var) {
    for (std::size_t i = 0; i < var.size(); ++i)
      var[i] = max(0.0, var[i]);
    return var;
  }

  /*
  template<uses_recursive_aproximation recursive,uses_averaging_aproximation
  averaging, uses_variance_aproximation variance> auto run_old(const Patch_State
  &t_prior, Qdt const &t_Qdt, Patch_Model const &m, const Experiment_step &p,
  double fs) const { auto &p_y = get<Patch_current>(p); auto &p_P_mean =
  get<P_mean>(t_prior); auto &p_P_Cov = get<P_Cov>(t_prior);

    double e =
        get<curr_noise>(m).value() * get<number_of_samples>(p).value() / fs;
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

    //      return Op(false, "\nfails in trace!!!; error=" + test.error() +
    //      ss.str());
    //    } else
    return Patch_State(logL(get<logL>(t_prior)()+v_plogL()),v_P_mean, v_P_cov,
  v_y_mean, v_y_var, v_plogL, v_eplogL, v_vplogL);
  }

  */

  Maybe_error<Patch_State> DVR(const Patch_State &t_prior, Qdt const &t_Qdt,
                               Patch_Model const &m, const Experiment_step &p,
                               double fs) const {
    auto &p_y = get<Patch_current>(p);
    auto &p_P_mean = get<P_mean>(t_prior);
    auto &p_P_Cov = get<P_Cov>(t_prior);

    double e =
        get<curr_noise>(m).value() * get<number_of_samples>(p).value() / fs;
    double N = get<N_Ch_mean>(m)();

    auto N_states = p_P_mean().ncols();
    Matrix<double> u(N_states, 1, 1.0);

    auto SmD = p_P_Cov() - diag(p_P_mean());

    if (std::isnan(p_y.value())) {
      auto v_P_cov = P_Cov(AT_B_A(get<P>(t_Qdt)(), SmD));
      auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)());
      v_P_cov() = v_P_cov() + diag(v_P_mean());

      return Patch_State(
          logL(get<logL>(t_prior)()), elogL(get<elogL>(t_prior)()),
          vlogL(get<vlogL>(t_prior)()), v_P_mean, v_P_cov, y_mean(NaN),
          y_var(NaN), plogL(NaN), eplogL(NaN), vplogL(NaN));
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
    double gSg = xtAx(get<gmean_i>(t_Qdt)(), SmD) +
                 getvalue(p_P_mean() *
                          zip([](auto x, auto y) { return x * y; },
                              get<gtotal_ij>(t_Qdt)(), get<gmean_ij>(t_Qdt)()) *
                          u);

    double sSg =
        xtAy(get<gvar_i>(t_Qdt)(), SmD, get<gmean_i>(t_Qdt)()) +
        getvalue(p_P_mean() *
                 zip([](auto x, auto y) { return x * y; },
                     get<gtotal_var_ij>(t_Qdt)(), get<gmean_ij>(t_Qdt)()) *
                 u);

    double sSs =
        xtAx(get<gvar_i>(t_Qdt)(), SmD) +
        getvalue(p_P_mean() *
                 zip([](auto x, auto y) { return x * y; },
                     get<gtotal_var_ij>(t_Qdt)(), get<gvar_ij>(t_Qdt)()) *
                 u);

    auto sS = tr(get<gvar_i>(t_Qdt)()) * SmD * get<P>(t_Qdt)() +
              p_P_mean() * get<gtotal_var_ij>(t_Qdt)();

    auto gS = tr(get<gmean_i>(t_Qdt)()) * SmD * get<P>(t_Qdt)() +
              p_P_mean() * get<gtotal_ij>(t_Qdt)();

    double ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());

    double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    auto e_mu = e + N * ms0;

    auto v_y_mean = y_mean(N * getvalue(p_P_mean() * get<gmean_i>(t_Qdt)()) -
                           N * 0.5 / e_mu * sSg);

    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto v_y_var = y_var(std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e));
    auto dy = p_y.value() - v_y_mean();

    auto chi = dy / v_y_var();

    auto v_P_mean = P_mean(p_P_mean() * get<P>(t_Qdt)() + chi * gS -
                           (chi * zeta * sSg + 0.5 / e_mu) * sS);

    auto v_P_cov =
        build<P_Cov>(AT_B_A(get<P>(t_Qdt)(), SmD) + diagpos(v_P_mean()) -
              (zeta + N / v_y_var() * sqr(zeta * sSg)) * XTX(sS) +
              (2.0 * N / v_y_var() * zeta * sSg) * X_plus_XT(tr(sS) * gS) -
              (N / v_y_var()) * XTX(gS));

    auto chi2 = dy * chi;

    auto v_plogL = plogL(0);
    if (v_y_var() > 0)
      v_plogL =
          plogL(-0.5 * log(2 * std::numbers::pi * v_y_var()) - 0.5 * chi2);
    else
      v_plogL = plogL(std::numeric_limits<double>::infinity());

    auto v_eplogL = eplogL(-0.5 * log(2 * std::numbers::pi * v_y_var()) -
                           0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    vplogL v_vplogL(0.5);
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
    std::cerr << get<Time>(p).value() << "\t" << v_P_mean << "\n";
    return Patch_State(logL(get<logL>(t_prior)() + v_plogL()),
                       elogL(get<elogL>(t_prior)() + v_eplogL()),
                       vlogL(get<vlogL>(t_prior)() + v_vplogL()), v_P_mean,
                       v_P_cov, v_y_mean, v_y_var, v_plogL, v_eplogL, v_vplogL);
  }

  template <uses_recursive_aproximation recursive,
            uses_averaging_aproximation averaging,
            uses_variance_aproximation variance,
           class C_Patch_State, class C_Qdt,class C_Patch_Model>
      requires (U<C_Patch_State,Patch_State>&&U<C_Patch_Model,Patch_Model>&&U<C_Qdt,Qdt>)
  
  
  Maybe_error<C_Patch_State> Macror(const C_Patch_State &t_prior, C_Qdt const &t_Qdt,
                                  C_Patch_Model const &m,
                                    const Experiment_step &p, double fs) {
    
    using Transf=transformation_type_t<C_Qdt>;
    auto &p_P_cov = get<P_Cov>(t_prior);
    auto &p_P_mean = get<P_mean>(t_prior);
    auto &y = get<Patch_current>(p).value();

    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto e =
        get<curr_noise>(m).value() * get<number_of_samples>(p).value() / fs;
    auto N = get<N_Ch_mean>(m)();
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

    auto ms = getvalue(p_P_mean() * t_gvar_i());
    
    Op_t<Transf,double> e_mu;
    Op_t<Transf,y_mean> r_y_mean;
    Op_t<Transf,y_var> r_y_var;

    double sSg;
    double sSs;
    double zeta;
    auto t_P = get<P>(t_Qdt);

    if constexpr ((!variance.value) && (!recursive.value)) {
      e_mu = e + N * ms;
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i());
      r_y_var() = e_mu + N * gSg;
      if (!(r_y_var() > 0)) {
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

      auto  ms = getvalue(p_P_mean() * t_gvar_i());

      e_mu = e + N * ms;
      r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i());
      r_y_var() = e_mu + N * gSg;
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

      auto delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
      auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

      e_mu = e + N * ms0;
      r_y_mean() =
          N * getvalue(p_P_mean() * t_gmean_i()) - N * 0.5 / e_mu * sSg;
      zeta = N / (2 * sqr(e_mu) + N * sSs);
      r_y_var() = std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e);
      if (!(r_y_var() > 0)) {
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
      if (r_test)
        return Op_t<Transf,Patch_State>(
            Op_t<Transf,logL>(get<logL>(t_prior)()), Op_t<Transf,elogL>(get<elogL>(t_prior)()),
            Op_t<Transf,vlogL>(get<vlogL>(t_prior)()),
            normalize(std::move(r_P_mean), t_min_P()),
            normalize(std::move(r_P_cov), t_min_P()), std::move(r_y_mean),
            std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN));
      else
        return error_message("fails at intertrace prediction!!: " +
                             r_test.error()());
    }

    auto dy = y - r_y_mean();
    auto chi = dy / r_y_var();
    Op_t<Transf,P_mean> r_P_mean;
    Op_t<Transf,P_Cov> r_P_cov;
    if constexpr (!recursive.value) {
      r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
      
      r_P_mean = build<P_mean>(p_P_mean() * t_P());
      r_P_cov() = r_P_cov() + diag(r_P_mean());
    } else if constexpr (!variance.value) {
      auto gS =
          TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
      auto gseg=chi*gS;
      
      
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
    
    Op_t<Transf,plogL> r_plogL;
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
    
    Op_t<Transf,eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
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
    } else
      return build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
                                build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
                                build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
                         normalize(std::move(r_P_mean), t_min_P()),
                         normalize(std::move(r_P_cov), t_min_P()),
                         std::move(r_y_mean), std::move(r_y_var),
                         r_plogL, r_eplogL,
                         r_vlogL);
  }
  
  template<class C_Patch_Model>
      requires U<C_Patch_Model,Patch_Model>
  auto init(const C_Patch_Model &m,
            initial_ATP_concentration initial_x)->Maybe_error<Transfer_Op_to<C_Patch_Model,Patch_State>> {
    auto v_Qx = calc_eigen(m, initial_x());
    if (v_Qx) {
      auto r_P_mean = calc_Peq(v_Qx.value(), m);
      auto r_P_cov = build<P_Cov>(diagpos(r_P_mean()) - XTX(r_P_mean.value()));
      auto r_test =
          test<true>(r_P_mean, r_P_cov, get<Probability_error_tolerance>(m));
      if (r_test)
      {
        auto t_min_P = get<min_P>(m);
        if (true) {
          std::cerr << "initial\n";
          std::cerr << "r_P_mean" << r_P_mean;
          std::cerr << "r_P_cov" << r_P_cov;
          //  std::cerr<<"normalized r_P_cov"<<normalize(std::move(r_P_cov),
          //  t_min_P());
        }
        return Transfer_Op_to<C_Patch_Model,Patch_State>(logL(0.0), elogL(0.0), vlogL(0.0),
                           normalize(std::move(r_P_mean), t_min_P()),
                           normalize(std::move(r_P_cov), t_min_P()),
                           y_mean(NaN), y_var(NaN), plogL(NaN), eplogL(NaN),
                           vplogL(NaN));
      }
      else
      {
        return error_message("fails at init: " + r_test.error()());
      }
    } else
      return v_Qx.error();
  }
  
  
  template<class C_Parameters, class Model>
  auto log_Likelihood(const Model &model,const C_Parameters& par, const Experiment &e)->Maybe_error<Transfer_Op_to<C_Parameters,Vector_Space<logL,elogL,vlogL>>>  {
    
    
    using Transf=transformation_type_t<C_Parameters>;
    using C_Patch_State=Op_t<Transf,Patch_State>;
    auto m=model(par);
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto ini = init(m, get<initial_ATP_concentration>(e));
    
    
    auto gege=0;
    if (!ini)
      return ini.error();
    else {
      auto run = fold(
          get<Recording>(e)(), ini.value(),
          [this, &m, fs,&gege]( C_Patch_State const &t_prior,
                         Experiment_step const &t_step) {
            auto t_Qx = calc_eigen(m, get<ATP_concentration>(t_step));

            if (!t_Qx)
                return Maybe_error<C_Patch_State>(t_Qx.error());
            // print(std::cerr,t_Qx.value());
            auto t_Qdt = calc_Qdt(m, t_Qx.value(),
                                  get<number_of_samples>(t_step).value() / fs);
            
//            print(std::cerr,t_prior);
//            if (gege<10) ++gege;
//            else abort();
            
            if (false){
            auto test_Qdt =
                test(t_Qdt, get<Conductance_variance_error_tolerance>(m));

            if (!test_Qdt)
              return Maybe_error<C_Patch_State>(test_Qdt.error());
            } 
            return Macror<uses_recursive_aproximation(true),
                          uses_averaging_aproximation(2),
                          uses_variance_aproximation(false)>(t_prior, t_Qdt, m,
                                                             t_step, fs);
          });
      if (!run)
        return run.error();
      else
        return build<Vector_Space<logL,elogL,vlogL>>(get<logL>(run.value()),get<elogL>(run.value()),get<vlogL>(run.value()));
    }
  }

  Simulated_Step sub_sample(std::mt19937_64 &mt, Simulated_Step &&t_sim_step,
                            const Patch_Model &m, const Experiment_step &t_s,
                            P t_P, std::size_t n_sub, double e) {
    auto &t_g = get<g>(m);
    auto &N = get<N_channel_state>(t_sim_step);
    double ysum = 0;
    for (std::size_t i = 0; i < n_sub; ++i) {
      N = sample_Multinomial(mt, t_P, N);
      ysum += getvalue(N() * t_g());
    }
    auto t_e_step = t_s;
    
    get<Patch_current>(t_e_step) = Patch_current(ysum / n_sub+std::normal_distribution<double>()(mt)*std::sqrt(e));
    get<Recording>(get<Simulated_Experiment>(t_sim_step)())().push_back(
        t_e_step);
    // std::cerr<<t_e_step;
    std::cerr << N;
    return t_sim_step;
  }

  Maybe_error<Simulated_Step>
  init_sim(std::mt19937_64 &mt, const Patch_Model &m, const Experiment &e) {
    auto initial_x = get<initial_ATP_concentration>(e);
    auto v_Qx = calc_eigen(m, initial_x());
    if (!v_Qx)
      return v_Qx.error();
    auto r_P_mean = calc_Peq(v_Qx.value(), m);
    auto N = get<N_Ch_mean>(m);
    auto sim = Simulated_Experiment(
        Experiment(Recording{}, get<Frequency_of_Sampling>(e),
                   get<initial_ATP_concentration>(e)));
    auto N_state = sample_Multinomial(mt, r_P_mean, N());
    return Simulated_Step(std::move(N_state), std::move(sim));
  }
  
  template<class Model, class Id>
  Maybe_error<Simulated_Experiment> sample(std::mt19937_64 &mt,
                                           const Model &model,const Parameters<Id>& par, 
                                           const Experiment &e,
                                           const Simulation_Parameters &sim) {
    
    auto m=model(par);
    auto n_sub = get<Number_of_simulation_sub_steps>(sim);
    auto fs = get<Frequency_of_Sampling>(e).value();
    auto sim_recording = Recording{};

    auto ini = init_sim(mt, m, e);
    if (!ini)
      return ini.error();
    else {
      auto run = fold(
          get<Recording>(e)(), ini.value(),
          [this, &m, fs, n_sub, &mt](Simulated_Step &&t_sim_step,
                                     Experiment_step const &t_step) {
            auto t_Qx = calc_eigen(m, get<ATP_concentration>(t_step));

            if (!t_Qx)
              return Maybe_error<Simulated_Step>(t_Qx.error());
            // print(std::cerr,t_Qx.value());
            
            auto dt = get<number_of_samples>(t_step).value() / fs;
            auto sub_dt = dt / n_sub();
            auto t_min_P = get<min_P>(m);
            auto t_P = calc_P(m, t_Qx.value(), sub_dt, t_min_P());
            double e =
                get<curr_noise>(m).value() * get<number_of_samples>(t_step).value() / fs;
            return Maybe_error<Simulated_Step>(
                sub_sample(mt, std::move(t_sim_step), m, t_step, t_P, n_sub(),e));
          });
      if (!run)
        return run.error();
      else
        return get<Simulated_Experiment>(run.value());
    }
  }
};

struct Model1: public Model_Patch<Model1>{};


} // namespace macrodr

#endif // QMODEL_H
