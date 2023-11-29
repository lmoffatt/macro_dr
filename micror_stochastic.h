#ifndef MICROR_STOCHASTIC_H
#define MICROR_STOCHASTIC_H

#include "matrix.h"
#include "maybe_error.h"
#include "multivariate_normal_distribution.h"
#include "qmodel.h"
#include "variables.h"
#include <cmath>
#include <cstddef>
#include <random>

namespace macrodr {

class mt_64 : public var::Var<mt_64, std::mt19937_64> {};

class N_channel_transition_count
    : public var::Var<N_channel_transition_count, Matrix<std::size_t>> {};

class N_channel_transition_mean
    : public var::Var<N_channel_transition_mean, Matrix<double>> {

  friend bool is_valid(N_channel_transition_mean const &x) {

    return (var::max(x()) > 0) && (var::min(x()) >= 0.0);
  }
};

static auto sample_Multinomial_state(std::mt19937_64 &mt,
                                     P_mean const &t_P_mean, std::size_t N) {
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

auto N_total_ij__mean_to_count(std::mt19937_64 &mt,
                       N_channel_transition_mean const &nr) {
    
  double N_mean=var::sum(nr());
  std::size_t N_count=sample_N(mt,N_mean);
  
  auto v_Nr=N_channel_transition_mean(nr()*(N_count/N_mean));
  Matrix<std::size_t> out=applyMap([](double N){return std::size_t(std::floor(N));},v_Nr());
  
  Matrix<double> sum_r= v_Nr()-out;
  
  
  assert(N_count>var::sum(out));
  double Nr=N_count-var::sum(out);
  Matrix<double> r= sum_r/Nr;;
  
  Matrix<std::size_t> rs=multinomial_distribution(Nr,r)(mt);
  
  
  return N_channel_transition_count(out+rs);
}


class Micror_parameters_distribution;

class Micror_state
    : public var::Var<Micror_state, Vector_Space<N_channel_transition_count,N_channel_transition_mean>> {
  Micror_state(std::mt19937_64 &mt, N_channel_transition_mean const &mean)
        : base_type{Vector_Space(N_total_ij__mean_to_count(mt, mean),mean)} {}

public:
  using base_type =
      var::Var<Micror_state, Vector_Space<N_channel_transition_count,N_channel_transition_mean>>;
  using base_type::Var;
  auto size() const {
    return get<N_channel_transition_mean>((*this)())().size();
  }
  
  friend class Micror_parameters_distribution;
  
  friend auto stretch_move(std::mt19937_64 &mt,
                           std::uniform_real_distribution<double> &rdist,
                           const Micror_state &Xk, const Micror_state &Xj) {
    assert((Xj.size() == Xk.size()) &&
           "sum of vector fields of different sizes");

    auto Nk = get<N_channel_transition_mean>(Xk());
    auto Nj = get<N_channel_transition_mean>(Xj());

    N_channel_transition_mean out = Nj;
    double z;
    double is_valid = false;
    while (!is_valid) {
      z = std::pow(rdist(mt) + 1, 2) / 2.0;
      is_valid = true;
      for (std::size_t i = 0; i < Nj().size(); ++i){
          if (Nj()[i] + z * (Nk()[i] - Nj()[i]) >= 0.0){
              out()[i] += z * (Nk()[i] - Nj()[i]);}
          else{
              is_valid = false;}
      }
    }
    return std::tuple(Micror_state(mt, out), z);
  }
};

class Micror_parameters_distribution {
  N_Ch_mean_value m_N;
  dirchlet_distribution P_tra_dist;
  multinomial_distribution m_N_tra_dist;
  
public:
  Micror_state operator()(std::mt19937_64 &mt) {
    auto r_Nij_mean = N_channel_transition_mean(P_tra_dist(mt) * m_N());
      return Micror_state(mt,r_Nij_mean);
  }
  
  Micror_parameters_distribution(N_Ch_mean_value const& N, Ptotal_ij const& Ptotal ):
      m_N{N},P_tra_dist{Ptotal()*N()},m_N_tra_dist{N(),Ptotal()}{}
  Maybe_error<double> logP(const Micror_state &s) {
    return m_N_tra_dist.logP(get<N_channel_transition_count>(s())());
  }
};


class Micror_parameters_likelihood{
    template <
        class FuncTable,class Variables, class DataType>
    Maybe_error<double>
    logLikelihood(FuncTable&& f,const Micror_parameters_likelihood &lik,
                  Micror_state const &p, const Variables &var, const DataType &y)
    {
        
    }
};



template <uses_recursive_aproximation recursive,
          uses_averaging_aproximation averaging,
          uses_variance_aproximation variance, class FunctionTable,
          class C_Patch_State, class C_Qdt, class C_Patch_Model, class C_double>
  requires(
      /*(U<std::decay_t<C_Patch_State>,
         Patch_State>||U<std::decay_t<C_Patch_State>,
         Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
      U<C_double, double> && U<C_Qdt, Qdt>)

Maybe_error<C_Patch_State>
Micror_stochastic(FunctionTable &ftbl, C_Patch_State &&t_prior,
                  C_Qdt const &t_Qdt, C_Patch_Model const &m,
                  C_double const &Nch, const Patch_current &p_y, double fs,
                  std::size_t myseed) {

  // the idea is to build a thermodynamic integration procedure to estimate the
  // Patch_State

  // First we go from Macro -> Micro
  //  then we have to define the parameters of the Micro distribution
  // then the prior and likelihood of those parameters
  // then we build the thermodynamic integration procedure
  // we obtain the posterior distibution and the Evidence
  // we convert the posterior and evidence into Patch_State

  /// so, the parameters of the Micro distribution are defined in such a way
  /// that can be used in the emcee algorithm, that is their domain has to be
  /// affine, which means among other things unbounded Number of channels or
  /// probabilities are not bounded so we need another formulation [I reflected
  /// later on the subject and I think is bullshit, there is no problem with
  /// hard boundaries in affince sampling, you just put a prior probability of
  /// zero, rejecting jumps to prohibited values, so bounded problems are just
  /// fine, no need to change the variables]
  ///
  ///
  /// so, the parameters we have to fit are just the values of the matrix Nij.
  ///
  /// the stetch algorithm being linear has no problem with obtaining new values
  /// of Nij, the only thing is that the Nij will become real numbers after
  /// multiplied by a real, so we need to convert them back to the discrete
  ///  Nij, by assigning the fractional channels with a multinomial
  ///  distribution. This assigment is done just once during the
  /// stretch-step, so the calculations of logPrior and logLikelihood are done
  /// just once per strech-step
  ///
  /// so, the idea is to generate the Nij using one dirchlet distribution
  /// (instead of multinomial) for each state.
  ///
  ///
  ///  Old no longer relevant discussion below
  /// the easieast formulation is to just take the logarithm of each value and
  /// substract the sum (or the mean) there is technical detail and that is that
  /// we cannot take logarithms of zero, so we have to arbitrarily decide a
  /// value for the logarithm of zero.
  ///
  /// so, the complete parameterization is given a fixed number of channels,
  /// logRi + logRij, the logarithm of the initial ratio of channels in each
  /// state and the ratio of transitions.
  ///
  /// as the prior is expressed in terms of multinomial distributions we have to
  /// convert the continuous logRi-logRij to the discrete
  ///  Ni -Nij, by assigning the fractional channels with a multinomial
  ///  distribution. This assigment is done just once during the
  /// stretch-step, so the calculations of logPrior and logLikelihood are done
  /// just once per strech-step
  ///
  /// This is an extreme case for thermodynamic integration since we have just
  /// one data: per run.
  ///
  /// a version for doing a Micror-stochastic algorithm for the complete series
  /// of measurments could be built using a different parameterization  logRij +
  /// logRjk  where we just change the state at a given time fixing the starting
  /// point and the ending point. In this way we could do metropolis jump for
  /// just one measurement, keeping the rest of the random parameters fixed.
  /// This would allow to increase the dimensionality of the search enormously,
  /// and a more fine grain use of the information. On the other hand if we wish
  /// to use a jump that involves the entire set of measurements at the same
  /// time, we should have a number of walkers at least as big as the number of
  /// random-parameters we have to estimate, which is easy to become impossible.
  ///
  ///
  auto &p_P_cov = get<P_Cov>(t_prior);
  auto &p_P_mean = get<P_mean>(t_prior);

  std::size_t num_states = p_P_cov().nrows();
  auto num_par = num_states * (num_states + 1);
  std::string path = "";
  std::size_t num_scouts_per_ensemble =
      std::pow(2, std::ceil(std::log2(num_par)));

  auto max_num_simultaneous_temperatures = 100;
  auto thermo_jumps_every = 4ul;
  auto max_iter_warming = 100ul;
  auto max_iter_equilibrium = 10000ul;
  auto n_points_per_decade = 6ul;
  auto stops_at = 1e-2;
  bool includes_zero = true;

  auto tmi = thermo_by_max_iter<Matrix<double>>(
      path, "Iteri", num_scouts_per_ensemble, max_num_simultaneous_temperatures,
      thermo_jumps_every, max_iter_warming, max_iter_equilibrium,
      n_points_per_decade, stops_at, includes_zero, myseed);

  auto opt = evidence(ftbl, std::move(tmi), my_linear_model.prior(),
                      my_linear_model.likelihood(), y, X);

  using Transf = transformation_type_t<C_Qdt>;
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
  auto gSg = getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
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
    auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();

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

    sSg =
        getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
        getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
    sSs = getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
          getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));

    auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    e_mu = e + N * ms0;
    r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) - N * 0.5 / e_mu * sSg +
                 y_baseline();
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
    auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
    auto gseg = chi * gS;

    r_P_mean() = p_P_mean() * t_P() + chi * gS;

    r_P_cov() = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) -
                (N / r_y_var()) * XTX(gS);

  } else {
    auto gS = TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
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
    return build<Patch_State>(build<logL>(get<logL>(t_prior)() + r_plogL()),
                              build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
                              build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
                              normalize(std::move(r_P_mean), t_min_P()),
                              normalize(std::move(r_P_cov), t_min_P()),
                              std::move(r_y_mean), std::move(r_y_var), r_plogL,
                              r_eplogL, r_vlogL);
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

} // namespace macrodr

#endif // MICROR_STOCHASTIC_H
