#include "CLI_macro_dr.h"
#include "distributions.h"
#include "experiment.h"
#include "function_measure_verification_and_optimization.h"
#include "function_memoization.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "maybe_error.h"
#include "parameters_derivative.h"
#include "parameters_distribution.h"
#include "variables.h"
#include "variables_derivative.h"
// #include "multivariate_normal_distribution.h"
#include "allosteric_models.h"
#include "cuevi.h"
#include "micror_stochastic.h"
#include "parallel_tempering.h"
#include "qmodel.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
using namespace macrodr;



int main(int argc, char **argv) {

  constexpr bool test_Bayesian_Linear_Regression = false;
  if constexpr (test_Bayesian_Linear_Regression) {
    /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
    auto myseed = 9762841416869310605ul;
    //  auto myseed = 0;

    /**
     * @brief npar number of parameters of the linear model
     */
    auto npar = 40ul;

    /**
     * @brief nsamples number of samples of the linear model
     */
    auto nsamples = 5000ul;

    /**
     * @brief std_log10_std_par  standard deviation of the log10 of the standard
     * deviation of the parameters values
     */
    auto std_log10_std_par = 1.0;

    /**
     * @brief mean_mean_par mean of the mean of the parameters values
     */
    auto mean_mean_par = 0.0;

    /**
     * @brief std_mean_par standard deviation of the mean of the parameters
     * values
     */
    auto std_mean_par = 10.0;

    /**
     * @brief mean_b mean of the linear coefficient of all parameters
     */
    auto mean_b = 0.0;

    /**
     * @brief std_b standard deviation of the paramters linear coefficient
     */
    auto std_b = 10.0;

    /**
     * @brief prior_error_factor error factor between the prior and real
     * parameter coefficient covariance value
     */
    auto prior_error_factor = 1;

    /**
     * @brief prior_eps_df prior value for the degrees of freedom used for
     * estimating the value of the variance of the data point error
     */
    double prior_eps_df = 1.0;

    /**
     * @brief prior_eps_variance prior value of the variance of the data point
     * error
     */
    double prior_eps_variance = 1.0;

    myseed = calc_seed(myseed);
    std::cerr << "myseed =" << myseed;

    /**
     * @brief mt random engine
     */
    auto mt = init_mt(myseed);

    /**
     * @brief a_0 parameter alpha of the beta distribution for the prior of the
     * data point error variance
     */
    auto a_0 = prior_eps_df / 2.0;

    /**
     * @brief b_0 parameter beta of the beta distribution for the prior of the
     * data point error variance
     */
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;

    /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
    std::size_t num_scouts_per_ensemble = 64;

    /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
    std::size_t max_num_simultaneous_temperatures = 4;

    /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */
    double n_points_per_decade = 3;

    /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
    double n_points_per_decade_fraction = 3;

    /**
     * @brief stops_at minimum value of beta greater than zero
     */
    double stops_at = 1e-7;

    /**
     * @brief includes_zero considers also beta equal zero
     */
    bool includes_zero = true;

    /**
     * @brief max_iter_warming maximum number of iterations on each warming step
     */
    std::size_t max_iter_warming = 200;

    /**
     * @brief max_iter_equilibrium maximum number of iterations on the
     * equilibrium step
     */
    std::size_t max_iter_equilibrium = 10000;

    /**
     * @brief path directory for the output
     */
    std::string path = "";

    /**
     * @brief filename prefix filename for the output
     */
    std::string filename = "A";

    /**
     * @brief checks_derivative_every_model_size number of steps before every
     * check of the derivative against the beta thermo parameter for stopping
     */
    std::size_t checks_derivative_every_model_size = 5000;

    /**
     * @brief max_ratio maximimum tolerated ratio for the beta derivative method
     */
    double max_ratio = 8;

    /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
    double min_fraction = 2;

    /**
     * @brief my_linear_model bayesian linear model a bayesian linear model with
     * all the prior information
     */
    auto my_linear_model =
        make_bayesian_linear_model(prior_eps_df, prior_eps_variance, npar,
                                   mean_b, std_b, prior_error_factor)
            .value();

    /**
     * @brief X random generated independent variables
     * @brief mean_par random generated mean of the parameters linear
     * coefficients
     * @brief cov_par random generated covariance of the parameters linear
     * coefficients
     */
    auto [X, mean_par, cov_par] = independent_variable_model(
        npar, std_log10_std_par, mean_mean_par, std_mean_par)(mt, nsamples);
    //  auto y = X * tr(b) + eps;

    /**
     * @brief b random generated linear coefficient values
     */
    auto b = sample(mt, my_linear_model);

    /**
     * @brief y simulated data using the constructed linear model, the
     * independent variables and the linear coefficient values
     */
    auto y = simulate(mt, my_linear_model, b, X);

    /**
     * @brief beta values for the thermodynamic parameter beta (that ranges from
     * 0 -only prior- to 1 -posterior likelihood
     */
    auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);

    /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
    std::size_t thermo_jumps_every = my_linear_model.size() * 1e0;
    auto ModelName = "Model";

    auto ftbl = FuncMap(
        path + ModelName + std::to_string(myseed) + "_" + time_now(),
        Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                  thermo_cuevi_randomized_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                       step_stretch_cuevi_mcmc_per_walker{}),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(logLikelihood_f{},
                       [](auto &&...x) {
                         return logLikelihood(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false), uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false), uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_Qdt{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_Qx{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qx(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_eigen{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_eigen(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2));

    if (false) {

      /**
       * @brief tmi classical thermodynamic algorithm ends by maximum iteration
       */
      auto tmi = thermo_by_max_iter<Matrix<double>>(
          path, "Iteri", num_scouts_per_ensemble,
          max_num_simultaneous_temperatures, thermo_jumps_every,
          max_iter_warming, max_iter_equilibrium, n_points_per_decade, stops_at,
          includes_zero, myseed);
      auto opt = evidence(ftbl, std::move(tmi), my_linear_model.prior(),
                          my_linear_model.likelihood(), y, X);
    }
    if (false) {
      /**
       * @brief tbc classical thermodynamic algorithm, ends using convergence
       * criteria
       */
      auto tbc = thermo_by_convergence<Matrix<double>>(
          path, "exp_thermo", num_scouts_per_ensemble,
          max_num_simultaneous_temperatures, thermo_jumps_every,
          checks_derivative_every_model_size, n_points_per_decade, stops_at,
          includes_zero, myseed);

      auto opt2 = evidence(ftbl, std::move(tbc), my_linear_model.prior(),
                           my_linear_model.likelihood(), y, X);

      // std::cout<<y;
    }
    if (true) {

      /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
      auto cbc = cuevi_by_convergence<Matrix<double>>(
          path, "exp_cuevi_40", num_scouts_per_ensemble,
          max_num_simultaneous_temperatures, min_fraction, thermo_jumps_every,
          checks_derivative_every_model_size, max_ratio, n_points_per_decade,
          n_points_per_decade_fraction, stops_at, includes_zero, myseed);

      bool all_at_once = false;
      auto opt3 = evidence(ftbl, std::move(cbc), my_linear_model.prior(),
                           my_linear_model.likelihood(), y, X, all_at_once);
    }
  }

  constexpr bool test_dynamice_command_line_interprester = false;

  if constexpr (test_dynamice_command_line_interprester) {
    auto cm = dcli::Compiler{};

    cm.push_function("load_experiment",
                     dcli::to_typed_function<std::string, double, double>(
                         &macrodr::load_experiment, "filename",
                         "frequency_of_sampling", "initial_ATP"));

    auto filename = "../macro_dr/run_script.txt";
    std::ifstream f(filename);

    std::string s;
    while (f) {
      std::string line;
      std::getline(f, line);
      s += line + "\n";
    }
    std::cout << "\ntest file \n" << s << "\n";
    auto p = dcli::extract_program(s);

    std::cerr << p;

    if (p) {
      auto c = dcli::compile_program(cm, p.value());
      if (c) {
        auto exec = c.value().run();
      }
    }

    if (p) {
      auto ss = p.value().str();
      std::cerr << ss;
    } else
      std::cerr << p.error()();
  }

  auto Efilename_all = "../macro_dr/Moffatt_Hume_2007_ATP_time.txt";

  auto [recording_conditions_all, recording_all] =
      macrodr::load_recording(Efilename_all);

  auto experiment_all = Experiment(
      std::move(recording_conditions_all), Frequency_of_Sampling(50e3),
      initial_ATP_concentration(ATP_concentration(0.0)));

  auto Efilename_7 = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt";
  auto Efilename_7_const_dt = "../macro_dr/Moffatt_Hume_2007_ATP_time_7_constant_dt.txt";

  auto [recording_conditions_7, recording_7] =
      macrodr::load_recording(Efilename_7);
  
  auto [recording_conditions_7_dt, recording_7_dt] =
      macrodr::load_recording(Efilename_7_const_dt);
  
  auto experiment_7 =
      Experiment(std::move(recording_conditions_7), Frequency_of_Sampling(50e3),
                 initial_ATP_concentration(ATP_concentration(0.0)));
  
  auto experiment_7_dt =
      Experiment(std::move(recording_conditions_7_dt), Frequency_of_Sampling(50e3),
                 initial_ATP_concentration(ATP_concentration(0.0)));
  
  auto &experiment = experiment_7;
  auto &recording = recording_7;

  struct Model0 : public Model_Patch<Model0> {};
  struct Model1 : public Model_Patch<Model1> {};

  struct Allost1 : public Model_Patch<Allost1> {};

  auto model00_7 = Model0::Model([]() {
    auto names_model = std::vector<std::string>{"kon",
                                                "koff",
                                                "gatin_on",
                                                "gating_off",
                                                "inactivation_rate",
                                                "unitary_current"};
    auto names_other = std::vector<std::string>{
        "Current_Noise", "Current_Baseline", "Num_ch_mean"};

    std::size_t N = 6ul;

    auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Q0_formula()[1][0] = "koff";
    v_Q0_formula()[2][1] = "2*koff";
    v_Q0_formula()[3][2] = "3*koff";
    v_Q0_formula()[3][4] = "gatin_on";
    v_Q0_formula()[4][3] = "gating_off";
    v_Q0_formula()[0][5] = "inactivation_rate";

    auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Qa_formula()[0][1] = "3*kon";
    v_Qa_formula()[1][2] = "2*kon";
    v_Qa_formula()[2][3] = "kon";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[4] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(),
                       names_other.end());
    auto p = Parameters<Model0>(
        std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1e-4, 1, 1000});

    auto logp =
        Parameters<Model0>(apply([](auto x) { return std::log10(x); }, p()));

    return std::tuple(
        [](const auto &logp) {
          using std::pow;
          auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
          auto kon = p[0];
          auto koff = p[1];
          auto gating_on = p[2];
          auto gating_off = p[3];
          auto inactivation_rate = p[4];
          auto v_unitary_current = p[5] * -1.0;
          auto Npar = 6ul;
          auto v_curr_noise = p[Npar];
          auto v_baseline = logp()[Npar + 1];
          //  auto v_Num_ch_mean=p[Npar+2];
          //  auto v_std_log_Num_ch=p[Npar+3];

          auto Npar2 = Npar + 1;
          auto v_N0 = p[std::pair(Npar2 + 1, Npar2 + 1)];
          return build<Patch_Model>(
              N_St(6),
              build<Q0>(var::build_<Matrix<double>>(
                  6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                  {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                   gating_off})),
              build<Qa>(var::build_<Matrix<double>>(
                  6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
              build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                   {v_unitary_current})),
              build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(12100),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        logp, typename Parameters<Model0>::Names(names_model),
        std::move(v_Q0_formula), std::move(v_Qa_formula),
        std::move(v_g_formula));
  });

  auto prior_model00_7 = Custom_Distribution(
      9ul,
      [](std::mt19937_64 &mt) {
        double kon = 10.0;
        double koff = 100;
        double gating_on = 100;
        double gating_off = 100;
        double inactivation_rate = 1e-4;
        double unitary_current = 0.5;
        double current_noise = 1e-4;
        double current_baseline = 1; // zero after log10
        double Num_ch_mean = 1000;
        double global_std_log10_par = 2.0;

        auto parv = std::vector<double>{kon,
                                        koff,
                                        gating_on,
                                        gating_off,
                                        inactivation_rate,
                                        unitary_current,
                                        current_noise,
                                        current_baseline,
                                        Num_ch_mean};

        auto n_parv = parv.size();

        auto par = Matrix<double>(parv.size(), 1ul, parv);
        auto logpar = apply([](auto x) { return std::log10(x); }, par);

        auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);

        auto covPar = DiagPosDetMatrix<double>(s_logpar);

        auto distPar = make_multivariate_normal_distribution(logpar, covPar);

        auto sample_par = distPar.value()(mt);

        // now recalculate the N_ch_numbers according to segements

        return Parameters<Model0>(sample_par);
      },
      [](Parameters<Model0> const &logp) {
        double kon = 10.0;
        double koff = 100;
        double gating_on = 100;
        double gating_off = 100;
        double inactivation_rate = 1e-4;
        double unitary_current = 0.5;
        double current_noise = 1e-4;
        double current_baseline = 1; // zero after log10
        double Num_ch_mean = 1000;
        double global_std_log10_par = 2.0;

        auto parv = std::vector<double>{kon,
                                        koff,
                                        gating_on,
                                        gating_off,
                                        inactivation_rate,
                                        unitary_current,
                                        current_noise,
                                        current_baseline,
                                        Num_ch_mean};

        auto n_parv = parv.size();
        auto par = Matrix<double>(parv.size(), 1ul, parv);
        auto logpar = apply([](auto x) { return std::log10(x); }, par);

        auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);

        auto covPar = DiagPosDetMatrix<double>(s_logpar);

        auto distPar = make_multivariate_normal_distribution(logpar, covPar);

        return distPar.value().logP(logp());
      });

  auto model00 = Model0::Model([]() {
    auto names_model = std::vector<std::string>{"kon",
                                                "koff",
                                                "gatin_on",
                                                "gating_off",
                                                "inactivation_rate",
                                                "unitary_current"};
    auto names_other = std::vector<std::string>{
        "Current_Noise", "Current_Baseline", "Num_ch_mean", "Num_ch_stddev",
        "Num_ch_0",      "Num_ch_1",         "Num_ch_2",    "Num_ch_3",
        "Num_ch_4",      "Num_ch_5",         "Num_ch_6"};

    std::size_t N = 6ul;

    auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Q0_formula()[1][0] = "koff";
    v_Q0_formula()[2][1] = "2*koff";
    v_Q0_formula()[3][2] = "3*koff";
    v_Q0_formula()[3][4] = "gatin_on";
    v_Q0_formula()[4][3] = "gating_off";
    v_Q0_formula()[0][5] = "inactivation_rate";

    auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Qa_formula()[0][1] = "3*kon";
    v_Qa_formula()[1][2] = "2*kon";
    v_Qa_formula()[2][3] = "kon";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[4] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(),
                       names_other.end());
    auto p = Parameters<Model0>(
        std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1e-4, 1, 1000, 100,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000});

    auto logp =
        Parameters<Model0>(apply([](auto x) { return std::log10(x); }, p()));

    return std::tuple(
        [](const auto &logp) {
          using std::pow;
          auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
          auto kon = p[0];
          auto koff = p[1];
          auto gating_on = p[2];
          auto gating_off = p[3];
          auto inactivation_rate = p[4];
          auto v_unitary_current = p[5] * -1.0;
          auto Npar = 6ul;
          auto v_curr_noise = p[Npar];
          auto v_baseline = logp()[Npar + 1];
          //  auto v_Num_ch_mean=p[Npar+2];
          //  auto v_std_log_Num_ch=p[Npar+3];

          auto Npar2 = Npar + 3;
          auto v_N0 = p[std::pair(Npar2 + 1, Npar2 + 7)];
          return build<Patch_Model>(
              N_St(6),
              build<Q0>(var::build_<Matrix<double>>(
                  6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                  {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                   gating_off})),
              build<Qa>(var::build_<Matrix<double>>(
                  6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
              build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                   {v_unitary_current})),
              build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(121),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        logp, typename Parameters<Model0>::Names(names_model),
        std::move(v_Q0_formula), std::move(v_Qa_formula),
        std::move(v_g_formula));
  });

  auto prior_model00 = Custom_Distribution(
      17ul,
      [](std::mt19937_64 &mt) {
        double kon = 10.0;
        double koff = 100;
        double gating_on = 100;
        double gating_off = 100;
        double inactivation_rate = 1e-4;
        double unitary_current = 0.5;
        double current_noise = 1e-4;
        double current_baseline = 1; // zero after log10
        double Num_ch_mean = 1000;
        double std_log_Num_ch = 5e-2;
        double Num_ch_factor_n = 1.0;
        std::size_t N_ch_segments = 7;
        double global_std_log10_par = 2.0;

        auto parv = std::vector<double>{kon,
                                        koff,
                                        gating_on,
                                        gating_off,
                                        inactivation_rate,
                                        unitary_current,
                                        current_noise,
                                        current_baseline,
                                        Num_ch_mean,
                                        std_log_Num_ch};

        auto n_parv = parv.size();

        parv.insert(parv.end(), N_ch_segments, Num_ch_factor_n);

        auto par = Matrix<double>(parv.size(), 1ul, parv);
        auto logpar = apply([](auto x) { return std::log10(x); }, par);

        auto s_logpar = std::vector<double>(n_parv - 1, global_std_log10_par);
        s_logpar.insert(s_logpar.end(), N_ch_segments + 1, 1);

        auto covPar = DiagPosDetMatrix<double>(s_logpar);

        auto distPar = make_multivariate_normal_distribution(logpar, covPar);

        auto sample_par = distPar.value()(mt);

        // now recalculate the N_ch_numbers according to segements

        auto s_log_Num_ch_mean = sample_par[n_parv - 2];
        auto s_std_log_Num_ch_cv = std::pow(10.0, sample_par[n_parv - 1]);

        for (std::size_t i = n_parv; i < par.size(); ++i)
          sample_par[i] =
              sample_par[i] / global_std_log10_par * s_std_log_Num_ch_cv +
              s_log_Num_ch_mean;

        return Parameters<Model0>(sample_par);
      },
      [](Parameters<Model0> const &logp) {
        double kon = 10.0;
        double koff = 100;
        double gating_on = 100;
        double gating_off = 100;
        double inactivation_rate = 1e-4;
        double unitary_current = 0.5;
        double current_noise = 1e-4;
        double current_baseline = 1; // zero after log10
        double Num_ch_mean = 1000;
        double std_logNum_ch = 5e-2;
        double Num_ch_factor_n = 1.0; //
        std::size_t N_ch_segments = 7;
        double global_std_log10_par = 2.0;

        auto parv = std::vector<double>{kon,
                                        koff,
                                        gating_on,
                                        gating_off,
                                        inactivation_rate,
                                        unitary_current,
                                        current_noise,
                                        current_baseline,
                                        Num_ch_mean,
                                        std_logNum_ch};

        auto n_parv = parv.size();
        auto p_Num_ch_mean = std::pow(10.0, logp()[n_parv - 2]);
        auto p_std_logNum_ch = std::pow(10.0, logp()[n_parv - 1]);

        parv.insert(parv.end(), N_ch_segments, p_Num_ch_mean);

        auto par = Matrix<double>(parv.size(), 1ul, parv);
        auto logpar = apply([](auto x) { return std::log10(x); }, par);

        auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);
        s_logpar.insert(s_logpar.end(), N_ch_segments, p_std_logNum_ch);

        auto covPar = DiagPosDetMatrix<double>(s_logpar);

        auto distPar = make_multivariate_normal_distribution(logpar, covPar);

        return distPar.value().logP(logp());
      });

  auto model01 = Model0::Model([]() {
    auto names_model = std::vector<std::string>{"kon",
                                                "koff",
                                                "gatin_on",
                                                "gating_off",
                                                "inactivation_rate",
                                                "unitary_current"};
    auto names_other =
        std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

    std::size_t N = 6ul;

    auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Q0_formula()[1][0] = "koff";
    v_Q0_formula()[2][1] = "2*koff";
    v_Q0_formula()[3][2] = "3*koff";
    v_Q0_formula()[3][4] = "gatin_on";
    v_Q0_formula()[4][3] = "gating_off";
    v_Q0_formula()[0][5] = "inactivation_rate";

    auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Qa_formula()[0][1] = "3*kon";
    v_Qa_formula()[1][2] = "2*kon";
    v_Qa_formula()[2][3] = "kon";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[4] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(),
                       names_other.end());
    auto p = Parameters<Model0>(
        std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1000, 1e-4, 1});

    auto logp =
        Parameters<Model0>(apply([](auto x) { return std::log10(x); }, p()));

    return std::tuple(
        [](const auto &logp) {
          using std::pow;
          auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
          auto kon = p[0];
          auto koff = p[1];
          auto gating_on = p[2];
          auto gating_off = p[3];
          auto inactivation_rate = p[4];
          auto v_unitary_current = p[5] * -1.0;
          auto Npar = 6ul;
          auto v_N0 = p[Npar];
          auto v_curr_noise = p[Npar + 1];
          auto v_baseline = logp()[Npar + 2];
          return build<Patch_Model>(
              N_St(6),
              build<Q0>(var::build_<Matrix<double>>(
                  6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                  {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                   gating_off})),
              build<Qa>(var::build_<Matrix<double>>(
                  6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
              build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                   {v_unitary_current})),
              build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline), Binomial_magical_number(5.0),
              min_P(1e-7), Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        logp, typename Parameters<Model0>::Names(names_model),
        std::move(v_Q0_formula), std::move(v_Qa_formula),
        std::move(v_g_formula));
  });

  auto model4 = Model0::Model([]() {
    auto names_model = std::vector<std::string>{"k01",
                                                "k10",
                                                "k12",
                                                "k21",
                                                "k23",
                                                "k32",
                                                "k34",
                                                "k43",
                                                "k45",
                                                "k54",
                                                "k46",
                                                "k64",
                                                "k57",
                                                "k75",
                                                "k67",
                                                "k08",
                                                "unitary_current"};
    auto names_other =
        std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

    std::size_t N = 9ul;

    auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Q0_formula()[0][8] = "k08";
    v_Q0_formula()[1][0] = "k10";
    v_Q0_formula()[2][1] = "k21";
    v_Q0_formula()[3][2] = "k32";
    v_Q0_formula()[3][4] = "k34";
    v_Q0_formula()[4][3] = "k43";
    v_Q0_formula()[4][5] = "k45";
    v_Q0_formula()[5][4] = "k54";
    v_Q0_formula()[4][6] = "k46";
    v_Q0_formula()[6][4] = "k64";
    v_Q0_formula()[5][7] = "k57";
    v_Q0_formula()[7][5] = "k75";
    v_Q0_formula()[6][7] = "k67";
    v_Q0_formula()[7][6] = "(k75 * k54 * k46 * k67)/ (k64 * k45 * k57)";
    auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Qa_formula()[0][1] = "k01";
    v_Qa_formula()[1][2] = "k12";
    v_Qa_formula()[2][3] = "k23";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[5] = "unitary_current";
    v_g_formula()[6] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(),
                       names_other.end());

    auto p_kinetics = std::vector<double>(16, 100.0);
    p_kinetics[0] = 10;
    p_kinetics[2] = 10;
    p_kinetics[4] = 10;
    p_kinetics[12] = 100;
    p_kinetics[15] = 1e-3;

    auto p_k_MH2007 = std::vector<double>{
        15.98, 0.019, 16.3,  380,   11.6,  6822, 3718, 43.54,
        540,   1088,  0.033, 0.246, 31.16, 79.0, 4.53, 1e-5};
    auto p_other = std::vector<double>{1, 4800, 1e-3, 1};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    p_k_MH2007.insert(p_k_MH2007.end(), p_other.begin(), p_other.end());

    //    auto p = Parameters<Model0>(p_kinetics);
    auto p = Parameters<Model0>(p_k_MH2007);

    auto logp =
        Parameters<Model0>(apply([](auto x) { return std::log10(x); }, p()));

    return std::tuple(
        [](const auto &logp) {
          using std::pow;
          auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
          auto Nst = 9ul;
          auto k01 = p[0];
          auto k10 = p[1];
          auto k12 = p[2];
          auto k21 = p[3];
          auto k23 = p[4];
          auto k32 = p[5];
          auto k34 = p[6];
          auto k43 = p[7];
          auto k45 = p[8];
          auto k54 = p[9];
          auto k46 = p[10];
          auto k64 = p[11];
          auto k57 = p[12];
          auto k75 = p[13];
          auto k67 = p[14];
          auto k76 = (k75 * k54 * k46 * k67) / (k64 * k45 * k57);
          auto k08 = p[15];
          auto v_g = p[16] * -1.0;
          auto Npar = 17ul;
          auto v_N0 = p[std::pair(Npar, Npar)];
          auto v_curr_noise = p[Npar + 1];
          auto v_baseline = logp()[Npar + 2];

          return build<Patch_Model>(
              N_St(Nst),
              build<Q0>(var::build_<Matrix<double>>(Nst, Nst,
                                                    {{1, 0},
                                                     {2, 1},
                                                     {3, 2},
                                                     {3, 4},
                                                     {4, 3},
                                                     {4, 5},
                                                     {5, 4},
                                                     {4, 6},
                                                     {6, 4},
                                                     {5, 7},
                                                     {7, 5},
                                                     {6, 7},
                                                     {7, 6},
                                                     {0, 8}},
                                                    {k10, k21, k32, k34, k43,
                                                     k45, k54, k46, k64, k57,
                                                     k75, k67, k76, k08})),
              build<Qa>(var::build_<Matrix<double>>(
                  Nst, Nst, {{0, 1}, {1, 2}, {2, 3}}, {k01, k12, k23})),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
              build<g>(var::build_<Matrix<double>>(Nst, 1, {{5, 0}, {6, 0}},
                                                   {v_g, v_g})),
              build<N_Ch_mean>(v_N0),

              build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(121000),
              Binomial_magical_number(5.0), min_P(1e-14),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        std::move(logp), std::move(names_model), v_Q0_formula, v_Qa_formula,
        v_g_formula);
  });

  auto model4_g_lin = Model0::Model([]() {
    auto names_model = std::vector<std::string>{"k01",
                                                "k10",
                                                "k12",
                                                "k21",
                                                "k23",
                                                "k32",
                                                "k34",
                                                "k43",
                                                "k45",
                                                "k54",
                                                "k46",
                                                "k64",
                                                "k57",
                                                "k75",
                                                "k67",
                                                "k08",
                                                "unitary_current"};
    auto names_other =
        std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

    std::size_t N = 9ul;

    auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Q0_formula()[0][8] = "k08";
    v_Q0_formula()[1][0] = "k10";
    v_Q0_formula()[2][1] = "k21";
    v_Q0_formula()[3][2] = "k32";
    v_Q0_formula()[3][4] = "k34";
    v_Q0_formula()[4][3] = "k43";
    v_Q0_formula()[4][5] = "k45";
    v_Q0_formula()[5][4] = "k54";
    v_Q0_formula()[4][6] = "k46";
    v_Q0_formula()[6][4] = "k64";
    v_Q0_formula()[5][7] = "k57";
    v_Q0_formula()[7][5] = "k75";
    v_Q0_formula()[6][7] = "k67";
    v_Q0_formula()[7][6] = "(k75 * k54 * k46 * k67)/ (k64 * k45 * k57)";
    auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
        N, std::vector<std::string>(N, "")));
    v_Qa_formula()[0][1] = "k01";
    v_Qa_formula()[1][2] = "k12";
    v_Qa_formula()[2][3] = "k23";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[5] = "unitary_current";
    v_g_formula()[6] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(),
                       names_other.end());

    auto p_kinetics = std::vector<double>(16, 100.0);
    p_kinetics[0] = 10;
    p_kinetics[2] = 10;
    p_kinetics[4] = 10;
    p_kinetics[12] = 100;
    p_kinetics[15] = 1e-3;

    auto p_other = std::vector<double>{0.1, 1000, 1e-3, 1};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Model0>(p_kinetics);

    auto logp =
        Parameters<Model0>(apply([](auto x) { return std::log10(x); }, p()));

    return std::tuple(
        [](const auto &logp) {
          using std::pow;
          auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
          auto Nst = 9ul;
          auto k01 = p[0];
          auto k10 = p[1];
          auto k12 = p[2];
          auto k21 = p[3];
          auto k23 = p[4];
          auto k32 = p[5];
          auto k34 = p[6];
          auto k43 = p[7];
          auto k45 = p[8];
          auto k54 = p[9];
          auto k46 = p[10];
          auto k64 = p[11];
          auto k57 = p[12];
          auto k75 = p[13];
          auto k67 = p[14];
          auto k76 = (k75 * k54 * k46 * k67) / (k64 * k45 * k57);
          auto k08 = p[15];
          auto v_g = logp()[16];
          auto Npar = 17ul;
          auto v_N0 = p[std::pair(Npar, Npar)];
          auto v_curr_noise = p[Npar + 1];
          auto v_baseline = logp()[Npar + 2];

          return build<Patch_Model>(
              N_St(Nst),
              build<Q0>(var::build_<Matrix<double>>(Nst, Nst,
                                                    {{1, 0},
                                                     {2, 1},
                                                     {3, 2},
                                                     {3, 4},
                                                     {4, 3},
                                                     {4, 5},
                                                     {5, 4},
                                                     {4, 6},
                                                     {6, 4},
                                                     {5, 7},
                                                     {7, 5},
                                                     {6, 7},
                                                     {7, 6},
                                                     {0, 8}},
                                                    {k10, k21, k32, k34, k43,
                                                     k45, k54, k46, k64, k57,
                                                     k75, k67, k76, k08})),
              build<Qa>(var::build_<Matrix<double>>(
                  Nst, Nst, {{0, 1}, {1, 2}, {2, 3}}, {k01, k12, k23})),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
              build<g>(var::build_<Matrix<double>>(Nst, 1, {{5, 0}, {6, 0}},
                                                   {v_g, v_g})),
              build<N_Ch_mean>(v_N0),

              build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(121000),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        std::move(logp), std::move(names_model), v_Q0_formula, v_Qa_formula,
        v_g_formula);
  });

  auto model6 = Allost1::Model([]() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto v_gating = Conformational_change_label{"Gating"};

    auto mo = make_Conformational_model(
        Agonist_dependency_map{

            std::map<Conformational_change_label, Agonist_dependency>{
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}},
                {v_gating, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_binding, v_binding, v_rocking, v_gating},
        std::vector<Conformational_interaction>{
            {Vector_Space{
                 Conformational_interaction_label{"BR"},
                 Conformational_interaction_players{{v_binding, v_rocking}},
                 Conformational_interaction_positions{
                     {{0, 3}, {1, 3}, {2, 3}}}},
             Vector_Space{
                 Conformational_interaction_label{"BG"},
                 Conformational_interaction_players{{v_binding, v_gating}},
                 Conformational_interaction_positions{
                     {{0, 4}, {1, 4}, {2, 4}}}},
             Vector_Space{
                 Conformational_interaction_label{"RG"},
                 Conformational_interaction_players{{v_rocking, v_gating}},
                 Conformational_interaction_positions{{{3, 4}}}}

            }},
        std::vector<Conductance_interaction>{
            Vector_Space{Conductance_interaction_label{"Gating_Current"},
                         Conductance_interaction_players{{v_gating}},
                         Conductance_interaction_positions{{{{4}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec = std::vector<std::string>{
        "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
        "Gating_on",  "Gating_off",  "BR",         "BR_0",
        "BR_1",       "BG",          "BG_0",       "BG_1",
        "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
    auto names_other = std::vector<std::string>{
        "Inactivation_rate", "Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics = std::vector<double>{
        10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
    auto p_other = std::vector<double>{1e-3, 1000, 1e-3, 1};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));

          p()[names["BR_0"].value()] =
              p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
          p()[names["BR_1"].value()] =
              p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

          p()[names["BG_0"].value()] =
              p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
          p()[names["BG_1"].value()] =
              p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

          p()[names["RG_0"].value()] =
              p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
          p()[names["RG_1"].value()] =
              p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

          p()[names["Gating_Current"].value()] =
              p()[names["Gating_Current"].value()] * -1.0;
          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          auto v_Inac_rate = p()[Npar];
          auto v_N0 = p()[std::pair{Npar + 1, Npar + 1}];
          auto v_curr_noise = p()[Npar + 2];
          auto v_baseline = logp()[Npar + 3];
          auto Nst = get<N_St>(m())();
          return add_Patch_inactivation(
              build<Patch_Model>(N_St(get<N_St>(m())), std::move(a_Q0),
                                 std::move(a_Qa),
                                 build<P_initial>(var::build_<Matrix<double>>(
                                     1, Nst, {{0, 0}}, {1.0})),
                                 std::move(a_g), build<N_Ch_mean>(v_N0),
                                 build<Current_Noise>(v_curr_noise),
                                 build<Current_Baseline>(v_baseline),
                                 N_Ch_mean_time_segment_duration(120000),
                                 Binomial_magical_number(5.0), min_P(1e-7),
                                 Probability_error_tolerance(1e-2),
                                 Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  auto model6_no_inactivation = Allost1::Model([]() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto v_gating = Conformational_change_label{"Gating"};

    auto mo = make_Conformational_model(
        Agonist_dependency_map{

            std::map<Conformational_change_label, Agonist_dependency>{
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}},
                {v_gating, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_binding, v_binding, v_rocking, v_gating},
        std::vector<Conformational_interaction>{
            {Vector_Space{
                 Conformational_interaction_label{"BR"},
                 Conformational_interaction_players{{v_binding, v_rocking}},
                 Conformational_interaction_positions{
                     {{0, 3}, {1, 3}, {2, 3}}}},
             Vector_Space{
                 Conformational_interaction_label{"BG"},
                 Conformational_interaction_players{{v_binding, v_gating}},
                 Conformational_interaction_positions{
                     {{0, 4}, {1, 4}, {2, 4}}}},
             Vector_Space{
                 Conformational_interaction_label{"RG"},
                 Conformational_interaction_players{{v_rocking, v_gating}},
                 Conformational_interaction_positions{{{3, 4}}}}

            }},
        std::vector<Conductance_interaction>{
            Vector_Space{Conductance_interaction_label{"Gating_Current"},
                         Conductance_interaction_players{{v_gating}},
                         Conductance_interaction_positions{{{{4}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec = std::vector<std::string>{
        "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
        "Gating_on",  "Gating_off",  "BR",         "BR_0",
        "BR_1",       "BG",          "BG_0",       "BG_1",
        "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
    auto names_other =
        std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics = std::vector<double>{
        10, 1000, 1000, 100000, 1, 100, 100, 1, 1, 1, 1, 1, 100, 1, 1, 1};
    auto p_Moffatt_Hume_transformed = std::vector<double>{
        9.28,     1871, 2547.88,   295207, 0.220378, 150.312, 74.865, 0.0323846,
        0.187903, 1.77, -0.457748, 1,      123,      1,       1.3411, 1};
    auto p_other = std::vector<double>{1000, 1e-3, 1};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));

          p()[names["BR_0"].value()] =
              p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
          p()[names["BR_1"].value()] =
              p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

          p()[names["BG_0"].value()] =
              p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
          p()[names["BG_1"].value()] =
              p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

          p()[names["RG_0"].value()] =
              p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
          p()[names["RG_1"].value()] =
              p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

          p()[names["Gating_Current"].value()] =
              p()[names["Gating_Current"].value()] * -1.0;
          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          // auto v_Inac_rate = p()[Npar];
          auto v_N0 = p()[std::pair{Npar, Npar}];
          auto v_curr_noise = p()[Npar + 1];
          auto v_baseline = logp()[Npar + 2];
          auto Nst = get<N_St>(m())();
          return build<Patch_Model>(
              N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
              std::move(a_g), build<N_Ch_mean>(v_N0),
              build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(120000),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  auto model6_Eff_no_inactivation = Allost1::Model([]() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto v_gating = Conformational_change_label{"Gating"};

    auto mo = make_Conformational_model(
        Agonist_dependency_map{

            std::map<Conformational_change_label, Agonist_dependency>{
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}},
                {v_gating, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_binding, v_binding, v_rocking, v_gating},
        std::vector<Conformational_interaction>{
            {Vector_Space{
                 Conformational_interaction_label{"BR"},
                 Conformational_interaction_players{{v_binding, v_rocking}},
                 Conformational_interaction_positions{
                     {{0, 3}, {1, 3}, {2, 3}}}},
             Vector_Space{
                 Conformational_interaction_label{"BG"},
                 Conformational_interaction_players{{v_binding, v_gating}},
                 Conformational_interaction_positions{
                     {{0, 4}, {1, 4}, {2, 4}}}},
             Vector_Space{
                 Conformational_interaction_label{"RG"},
                 Conformational_interaction_players{{v_rocking, v_gating}},
                 Conformational_interaction_positions{{{3, 4}}}}

            }},
        std::vector<Conductance_interaction>{
            Vector_Space{Conductance_interaction_label{"Gating_Current"},
                         Conductance_interaction_players{{v_gating}},
                         Conductance_interaction_positions{{{{4}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec_untransformed = std::vector<std::string>{
        "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
        "Gating_on",  "Gating_off",  "BR",         "BR_0",
        "BR_1",       "BG",          "BG_0",       "BG_1",
        "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
    auto names_vec = std::vector<std::string>{"Binding_Act_on",
                                              "Binding_Act_off",
                                              "Rocking_on_B",
                                              "Rocking_off_B",
                                              "Gating_on_BR",
                                              "Gating_off_BR",
                                              "BR",
                                              "BR_Bon",
                                              "BR_Ron",
                                              "BG",
                                              "BG_Bon",
                                              "BG_Gon",
                                              "RG",
                                              "RG_Ron",
                                              "RG_Gon",
                                              "Gating_Current"};

    auto names_other =
        std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics =
        std::vector<double>{9.28, 1871, 3875, 1.07, 914, 776, 65.1 * 1.15, 1.15,
                            33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
    auto p_other = std::vector<double>{4800, 1e-3, 1};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec_untransformed);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto tr_p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));
          auto Binding_on = tr_p()[0];
          auto Binding_off = tr_p()[1];
          auto Rocking_on_B = tr_p()[2];
          auto Rocking_off_B = tr_p()[3];

          auto Gating_on_BR = tr_p()[4];
          auto Gating_off_BR = tr_p()[5];
          auto BR = tr_p()[6];
          auto BR_Bon = tr_p()[7];
          auto BR_Ron = tr_p()[8];
          auto BG = tr_p()[9];
          auto BG_Bon = tr_p()[10];
          auto BG_Gon = tr_p()[11];
          auto RG = tr_p()[12];
          auto RG_Ron = tr_p()[13];
          auto RG_Gon = tr_p()[14];
          auto Gating_Current = tr_p()[15] * (-1.0);

          auto Rocking_on = Rocking_on_B / pow(BR_Bon, 3);
          auto Rocking_off = Rocking_off_B * pow(BR / BR_Bon, 3);

          auto Gating_on = Gating_off_BR / pow(BG_Gon, 3) / RG_Gon;
          auto Gating_off = Gating_off_BR * pow(BG / BG_Gon, 3) * RG / RG_Gon;

          auto BR_0 = log(BR_Bon) / log(BR);
          auto BR_1 = log(BR_Ron) / log(BR);

          auto BG_0 = log(BG_Bon) / log(BG);
          auto BG_1 = log(BG_Gon) / log(BG);

          auto RG_0 = log(RG_Ron) / log(RG);
          auto RG_1 = log(RG_Gon) / log(RG);

          auto p = tr_p;
          p()[2] = Rocking_on;
          p()[3] = Rocking_off;

          p()[4] = Gating_on;
          p()[5] = Gating_off;

          p()[7] = BR_0;
          p()[8] = BR_1;

          p()[10] = BG_0;
          p()[11] = BG_1;

          p()[13] = RG_0;
          p()[14] = RG_1;
          p()[15] = Gating_Current;
          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          // std::cerr<<"parameters\n"<<p();

          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          // auto v_Inac_rate = p()[Npar];
          auto v_N0 = tr_p()[std::pair{Npar, Npar}];
          auto v_curr_noise = tr_p()[Npar + 1];
          auto v_baseline = logp()[Npar + 2];
          auto Nst = get<N_St>(m())();
          return build<Patch_Model>(
              N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
              build<P_initial>(
                  var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
              std::move(a_g), build<N_Ch_mean>(v_N0),
              build<Current_Noise>(v_curr_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(120000),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  auto sd_log_model6_eff = std::vector<double>(15, 1e-5);
  sd_log_model6_eff.insert(sd_log_model6_eff.end(), 4, 2);
  sd_log_model6_eff[18] = 50;

  auto prior_model6_Eff_no_inactivation = var::prior_around(
      model6_Eff_no_inactivation.parameters(), sd_log_model6_eff);

  auto model7 = Allost1::Model([]() {
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_gating = Conformational_change_label{"Gating"};

    auto mo = make_Conformational_model(
        Agonist_dependency_map{

            std::map<Conformational_change_label, Agonist_dependency>{
                {v_rocking, Agonist_dependency{}},
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_gating, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_rocking, v_binding, v_rocking, v_binding, v_rocking,
            v_gating},
        std::vector<Conformational_interaction>{
            {Vector_Space{
                 Conformational_interaction_label{"BR"},
                 Conformational_interaction_players{{v_binding, v_rocking}},
                 Conformational_interaction_positions{
                     {{0, 1}, {2, 3}, {4, 5}}}},
             Vector_Space{
                 Conformational_interaction_label{"BG"},
                 Conformational_interaction_players{{v_binding, v_gating}},
                 Conformational_interaction_positions{
                     {{0, 6}, {2, 6}, {4, 6}}}},
             Vector_Space{
                 Conformational_interaction_label{"RG"},
                 Conformational_interaction_players{{v_rocking, v_gating}},
                 Conformational_interaction_positions{{{1, 6}, {3, 6}, {5, 6}}}}

            }},
        std::vector<Conductance_interaction>{
            Vector_Space{Conductance_interaction_label{"Gating_Current"},
                         Conductance_interaction_players{{v_gating}},
                         Conductance_interaction_positions{{{{6}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec = std::vector<std::string>{
        "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
        "Gating_on",  "Gating_off",  "BR",         "BR_0",
        "BR_1",       "BG",          "BG_0",       "BG_1",
        "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
    auto names_other = std::vector<std::string>{
        "Inactivation_rate", "Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics = std::vector<double>{
        10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
    auto p_other = std::vector<double>{1, 1, 100, 20, 1000};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));

          p()[names["BR_0"].value()] =
              p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
          p()[names["BR_1"].value()] =
              p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

          p()[names["BG_0"].value()] =
              p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
          p()[names["BG_1"].value()] =
              p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

          p()[names["RG_0"].value()] =
              p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
          p()[names["RG_1"].value()] =
              p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

          p()[names["Gating_Current"].value()] =
              p()[names["Gating_Current"].value()] * -1.0;
          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          auto v_Inac_rate = p()[Npar];
          auto v_N0 = p()[Npar + 1];
          auto v_curr_noise = p()[Npar + 2];
          auto v_baseline = logp()[Npar + 3];
          auto Nst = get<N_St>(m())();
          return add_Patch_inactivation(
              build<Patch_Model>(N_St(get<N_St>(m())), std::move(a_Q0),
                                 std::move(a_Qa),
                                 build<P_initial>(var::build_<Matrix<double>>(
                                     1, Nst, {{0, 0}}, {1.0})),
                                 std::move(a_g), build<N_Ch_mean>(v_N0),
                                 build<Current_Noise>(v_curr_noise),
                                 build<Current_Baseline>(v_baseline),
                                 Binomial_magical_number(5.0), min_P(1e-7),
                                 Probability_error_tolerance(1e-2),
                                 Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  auto model8 = Allost1::Model([]() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto mo = make_Conformational_model(
        Agonist_dependency_map{
            std::map<Conformational_change_label, Agonist_dependency>{
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_rocking, v_binding, v_rocking, v_binding, v_rocking},
        std::vector<Conformational_interaction>{
            {Vector_Space{Conformational_interaction_label{"RBR"},
                          Conformational_interaction_players{
                              {v_rocking, v_binding, v_rocking}},
                          Conformational_interaction_positions{{{5, 0, 1},
                                                                {1, 0, 5},
                                                                {1, 2, 3},
                                                                {3, 2, 1},
                                                                {3, 4, 5},
                                                                {5, 4, 3}}}}}},
        std::vector<Conductance_interaction>{Vector_Space{
            Conductance_interaction_label{"Rocking_Current_factor"},
            Conductance_interaction_players{{v_rocking}},
            Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec = std::vector<std::string>{
        "Binding_on",  "Binding_off", "Rocking_on",
        "Rocking_off", "RBR",         "RBR_0",
        "RBR_1",       "RBR_2",       "Rocking_Current_factor"}; //--> 8
    auto names_other =
        std::vector<std::string>{"Inactivation_rate", "Leaking_current",
                                 "Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics =
        std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 1.0, 100};
    auto p_other = std::vector<double>{1, 1, 100, 20, 1000, 1e-3};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));
          p()[names["RBR_0"].value()] =
              p()[names["RBR_0"].value()] / (1.0 + p()[names["RBR_0"].value()]);
          p()[names["RBR_1"].value()] =
              p()[names["RBR_1"].value()] / (1.0 + p()[names["RBR_1"].value()]);
          p()[names["RBR_2"].value()] =
              p()[names["RBR_2"].value()] / (1.0 + p()[names["RBR_2"].value()]);

          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          auto v_Inac_rate = p()[Npar];
          auto v_leaking_current = p()[Npar + 1];
          auto v_N0 = p()[Npar + 2];
          auto v_curr_noise = p()[Npar + 3];
          auto v_baseline = logp()[Npar + 4];

          auto v_g = build<g>(apply(
              [&v_leaking_current](const auto &x) {
                return v_leaking_current * pow(10.0, x) * (-1.0);
              },
              a_g()));
          auto Nst = get<N_St>(m())();
          return add_Patch_inactivation(
              build<Patch_Model>(N_St(get<N_St>(m())), std::move(a_Q0),
                                 std::move(a_Qa),
                                 build<P_initial>(var::build_<Matrix<double>>(
                                     1, Nst, {{0, 0}}, {1.0})),
                                 std::move(v_g), build<N_Ch_mean>(v_N0),
                                 build<Current_Noise>(v_curr_noise),
                                 build<Current_Baseline>(v_baseline),
                                 Binomial_magical_number(5.0), min_P(1e-7),
                                 Probability_error_tolerance(1e-2),
                                 Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  auto model9 = Allost1::Model([]() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto mo = make_Conformational_model(
        Agonist_dependency_map{
            std::map<Conformational_change_label, Agonist_dependency>{
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{
            v_binding, v_rocking, v_binding, v_rocking, v_binding, v_rocking},
        std::vector<Conformational_interaction>{{Vector_Space{
            Conformational_interaction_label{"RB"},
            Conformational_interaction_players{{v_rocking, v_binding}},
            Conformational_interaction_positions{
                {{5, 0}, {1, 0}, {1, 2}, {3, 2}, {3, 4}, {5, 4}}}}}},
        std::vector<Conductance_interaction>{Vector_Space{
            Conductance_interaction_label{"Rocking_Current_factor"},
            Conductance_interaction_players{{v_rocking}},
            Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

    assert(mo);
    auto m = std::move(mo.value());

    auto names = make_ModelNames<Allost1>(m);

    auto names_vec = std::vector<std::string>{
        "Binding_on", "Binding_off",
        "Rocking_on", "Rocking_off",
        "RB",         "RB_0",
        "RB_1",       "Rocking_Current_factor"}; //--> 8
    auto names_other = std::vector<std::string>{
        "Inactivation_rate", "Num_ch", "Current_Noise", "Current_Baseline"};

    auto p_kinetics =
        std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 100};
    auto p_other = std::vector<double>{1, 1, 100, 20, 1000, 1e-3};

    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Parameters<Allost1>(p_kinetics);

    auto logp =
        Parameters<Allost1>(apply([](auto x) { return std::log10(x); }, p()));

    assert(names() == names_vec);

    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &logp) {
          using std::pow;
          auto p = build<Parameters<Allost1>>(
              apply([](const auto &x) { return pow(10.0, x); }, logp()));
          p()[names["RB_0"].value()] =
              p()[names["RB_0"].value()] / (1.0 + p()[names["RB_0"].value()]);
          p()[names["RB_1"].value()] =
              p()[names["RB_1"].value()] / (1.0 + p()[names["RB_1"].value()]);

          auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
          assert(Maybe_Q0Qag);
          auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
          auto Npar = names().size();

          auto v_Inac_rate = p()[Npar];
          auto v_leaking_current = p()[Npar + 1];
          auto v_N0 = p()[Npar + 2];
          auto v_curr_noise = p()[Npar + 3];
          auto v_baseline = logp()[Npar + 4];

          auto v_g = build<g>(apply(
              [&v_leaking_current](const auto &x) {
                return v_leaking_current * pow(10.0, x) * (-1.0);
              },
              a_g()));
          auto Nst = get<N_St>(m())();
          return add_Patch_inactivation(
              build<Patch_Model>(N_St(get<N_St>(m())), std::move(a_Q0),
                                 std::move(a_Qa),
                                 build<P_initial>(var::build_<Matrix<double>>(
                                     1, Nst, {{0, 0}}, {1.0})),
                                 std::move(v_g), build<N_Ch_mean>(v_N0),
                                 build<Current_Noise>(v_curr_noise),
                                 build<Current_Baseline>(v_baseline),
                                 Binomial_magical_number(5.0), min_P(1e-7),
                                 Probability_error_tolerance(1e-2),
                                 Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        },
        logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
  });

  //    auto egrw=var::Derivative_t<decltype(Fun(Var<N_Ch_mean>{},[](auto
  //    N,auto...){return N;},9.0 )),Parameters<Model0>>;
  // using tetw=typename dsge::ik;

  auto param11 = Parameters<Model0>(apply(
      [](auto x) { return std::log10(x); },
      Matrix<double>(1, 14,
                     std::vector<double>{18, 12, 6, 210, 420, 630, 1680, 54,
                                         0.5, 100, 50, 1000, 1e-4, 1.0})));

  auto param4 = model4.parameters();
  auto param4Names = model4.names();
  auto param8 = model8.parameters();
  auto param8Names = model8.names();
  auto param9 = model9.parameters();
  auto param9Names = model9.names();
  auto param6 = model6.parameters();
  auto param6Names = model6.names();
  auto param7 = model7.parameters();
  auto param7Names = model7.names();
  auto param00 = model00.parameters();
  auto param00Names = model00.names();

  auto param00_7 = model00_7.parameters();
  auto param00_7Names = model00_7.names();

  auto param11Names = std::vector<std::string>{
      "k01",         "k12",          "k23",   "k10",         "k21",
      "k32",         "k34",          "k43",   "conductance", "Num_Chan_0",
      "Num_Chan_eq", "Num_Chan_tau", "noise", "baseline"};

  auto &model0 = model4;
  auto &param1Names = model0.names();
  auto &param1 = model0.parameters();
  std::string ModelName = "model4";
  using MyModel = Model0;
  
  auto &model0_alt = model6_Eff_no_inactivation;
  auto &param1Names_alt = model0_alt.names();
  auto &param1_alt = model0_alt.parameters();
  std::string ModelName_alt = "model6_eff_no_inactivation";
  using MyModel_alt = Allost1;
  
  
  assert(param1Names().size() == param1.size());
  auto dparam1 = var::selfDerivative(param1);

  auto NNN = build<Fun>(
      Var<N_Ch_mean>{}, [](auto N, auto...) { return N; }, param1[0]);
  auto dNNN = build<Fun>(
      Var<N_Ch_mean>{}, [](auto N, auto...) { return N; }, dparam1()[0]);

  auto Npar = 5ul;
  auto v_N0 = dparam1()[Npar + 2];
  auto v_Neq = dparam1()[Npar + 3];
  auto v_Ntau = dparam1()[Npar] + 4;

  auto ggg = build<Fun>(
      Var<N_Ch_mean>{},
      [](auto &N0, auto &Neq, auto &Ntau, auto &time) {
        return Neq + (N0 - Neq) * (1.0 - exp(-time() / Ntau));
      },
      v_N0, v_Neq, v_Ntau);

  // auto ggb = ggg(Time(10.0)); does not compile with N_Ch_mean being a Matrix

  // using jger=typename decltype(NNN)::kgerge;
  // using jger2=typename decltype(dNNN)::kgerge;

  // std::cerr<<"dparam1\n"<<dparam1;
  // auto n = build<N_Ch_mean>(dparam1()[0]);

  auto dp0 = dparam1()[0];

  auto qq =
      var::build_<Matrix<double>>(5, 5, {{0, 1}, {1, 2}, {2, 3}},
                                  {dparam1()[2], dparam1()[3], dparam1()[4]});

  auto m = model0(param1);
  auto dm = model0(dparam1);

  // std::cerr<<"\n!--------------------------------------------------------------------------!\n";
  //  print(std::cerr,dm);
  //  print(std::cerr,m);
  auto fs = get<Frequency_of_Sampling>(experiment).value();
  auto dini = macrodr::Macro_DMR{}.init<return_predictions(false)>(
      dm, get<initial_ATP_concentration>(experiment));
  auto ini = macrodr::Macro_DMR{}.init<return_predictions(false)>(
      m, get<initial_ATP_concentration>(experiment));

  auto t_step = get<Recording_conditions>(experiment)()[0];
  auto tQx = macrodr::Macro_DMR{}.calc_Qx(m, ATP_concentration(100.0));
  auto t_Qx = macrodr::Macro_DMR{}.calc_eigen(tQx);

  //  auto dt_Qx = macrodr::Macro_DMR{}.calc_eigen(dm,
  //  get<ATP_concentration>(t_step));

  // auto dt_Qdt = macrodr::Macro_DMR{}.calc_Qdt(m, t_Qx.value(),
  //                                             get<number_of_samples>(t_step).value()
  //                                             / fs);

  //  auto t_Qdt = macrodr::Macro_DMR{}.calc_Qdt(m, t_Qx.value(),1e-3);
  auto num_scouts_per_ensemble = 16ul;
  auto ftbl = FuncMap(
      "_" + time_now(),
      Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
              num_scouts_per_ensemble / 2),
      Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
              num_scouts_per_ensemble / 2),
      Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                thermo_cuevi_randomized_jump_mcmc{}),
              num_scouts_per_ensemble / 2),
      var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                     step_stretch_cuevi_mcmc_per_walker{}),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(logLikelihood_f{},
                     [](auto &&...x) {
                       return logLikelihood(std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false)>{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.Macror<uses_recursive_aproximation(true),
                                       uses_averaging_aproximation(2),
                                       uses_variance_aproximation(false),
                                       uses_variance_correction_aproximation(false)>(
                           std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false)>{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.Macror<uses_recursive_aproximation(false),
                                       uses_averaging_aproximation(2),
                                       uses_variance_aproximation(false), uses_variance_correction_aproximation(false)>(
                           std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(Calc_Qdt{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(Calc_Qx{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qx(std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2),
      var::Time_it(F(Calc_eigen{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_eigen(std::forward<decltype(x)>(x)...);
                     }),
                   num_scouts_per_ensemble / 2));

  if constexpr (false) {
    auto t_Qdt2 = macrodr::Macro_DMR{}.calc_Qdt(
        ftbl, m, ATP_step(number_of_samples(50), ATP_concentration(100.0)),
        50e3);
    auto t_Qdt3 = macrodr::Macro_DMR{}.calc_Qdt(
        ftbl, m,
        std::vector<ATP_step>{
            ATP_step(number_of_samples(10.0), ATP_concentration(100.0)),
            ATP_step(number_of_samples(40.0), ATP_concentration(100.0))},
        50e3);

    auto test = test_equality(t_Qdt2.value(), t_Qdt3.value(), eps);
    std::cerr << test;
    // std::abort();
  }

  std::random_device rd;
  typename std::mt19937_64::result_type seed = rd();

  std::mt19937_64 mt(seed);

  /**
   * @brief analyze_and_test_Qdt
   *
   * the objective of this section is to do two things
   *
   * 1. a sanity check of the coding of Qdt by comparing it to Qdt_by_bisection
   * both should give comparable results.
   *
   * 2. provide data to build some nice graphs ilustrating the concept of Qdt
   *
   *
   *  with all this in mind we need to show that both functions are comparable
   * in a wide range of concentration and dt and for current model. (model0)
   *
   */
  constexpr bool analyze_and_test_Qdt = false;
  if constexpr (analyze_and_test_Qdt) {

    auto orders = std::vector<std::size_t>{1, 2, 3, 5, 8, 12, 16};

    // auto xs = std::vector<ATP_concentration>{
    //     ATP_concentration(0),    ATP_concentration(100),
    //     ATP_concentration(200),  ATP_concentration(500),
    //     ATP_concentration(1000), ATP_concentration(2000),
    //     ATP_concentration(5000), ATP_concentration(10000)};

    auto xs = std::vector<ATP_concentration>{ATP_concentration(0.0),
                                             ATP_concentration(10000.0)};

    auto dts = std::vector<double>();
    double dt_min = 1e-6;
    double dt_max = 1;
    double n_per_decade_dt = 12;

    double log10_dt_run = std::log10(dt_min);
    while (log10_dt_run <= std::log10(dt_max)) {
      dts.push_back(std::pow(10.0, log10_dt_run));
      log10_dt_run += 1.0 / n_per_decade_dt;
    }

    std::string filename = "Qdt_comparison_";

    std::ofstream fo_i(filename + ModelName + "_istate_.txt");
    std::ofstream fo_ij(filename + ModelName + "_istate_jstate.txt");

    fo_i << "ATP"
         << ","
         << "dt"
         << ","
         << "algoritm"
         << ","
         << "order"
         << ","
         << "variable"
         << ","
         << "i_state"
         << ","
         << "value"
         << "\n";
    fo_ij << "ATP"
          << ","
          << "dt"
          << ","
          << "algoritm"
          << ","
          << "order"
          << ","
          << "variable"
          << ","
          << "i_state"
          << ","
          << "j_state"
          << ","
          << "value"
          << "\n";

    auto m = model0(param1);
    auto fs = 50e3;
    for (auto x : xs) {
      for (auto dt : dts) {
        //, , ,
        //  , , , ,
        auto n = number_of_samples(dt * fs);
        auto step = ATP_step(n, x);
        Qdt Qdtd = macrodr::Macro_DMR{}
                       .calc_Qdt(ftbl.fork(var::I_thread(0)), m, step, fs)
                       .value();
        for (std::size_t i = 0; i < get<gmean_i>(Qdtd)().size(); ++i) {
          std::apply(
              [x, dt, i, &fo_i](auto const &...g_i) {
                ((fo_i << x.value() << "," << dt << ","
                       << "Qdt_eig"
                       << "," << 0 << "," << className(g_i) << "," << i << ","
                       << g_i()[i] << "\n"),
                 ...);
              },
              std::forward_as_tuple(get<gmean_i>(Qdtd), get<gsqr_i>(Qdtd),
                                    get<gvar_i>(Qdtd)));
        }
        for (std::size_t i = 0; i < get<gtotal_ij>(Qdtd)().nrows(); ++i) {
          for (std::size_t j = 0; j < get<gtotal_ij>(Qdtd)().ncols(); ++j) {

            std::apply(
                [x, dt, &fo_ij, i, j](auto const &...g_ij) {
                  ((fo_ij << x.value() << "," << dt << ","
                          << "Qdt_eig"
                          << "," << 0 << "," << className(g_ij) << "," << i
                          << "," << j << "," << g_ij()(i, j) << "\n"),
                   ...);
                },
                std::forward_as_tuple(
                    get<P>(Qdtd), get<gtotal_ij>(Qdtd), get<gmean_ij>(Qdtd),
                    get<gtotal_sqr_ij>(Qdtd), get<gtotal_var_ij>(Qdtd),
                    get<gvar_ij>(Qdtd)));
          }
        }

        for (auto order : orders) {
          Qdt Q_dtb = macrodr::Macro_DMR{}
                          .calc_Qdt_bisection(ftbl.fork(var::I_thread(0)), m,
                                              step, fs, order)
                          .value();

          for (std::size_t i = 0; i < get<gmean_i>(Qdtd)().size(); ++i) {
            std::apply(
                [x, dt, i, &fo_i, order](auto const &...g_i) {
                  ((fo_i << x.value() << "," << dt << ","
                         << "Qdt_bis"
                         << "," << order << "," << className(g_i) << "," << i
                         << "," << g_i()[i] << "\n"),
                   ...);
                },
                std::forward_as_tuple(get<gmean_i>(Q_dtb), get<gsqr_i>(Q_dtb),
                                      get<gvar_i>(Q_dtb)));
          }
          for (std::size_t i = 0; i < get<gtotal_ij>(Q_dtb)().nrows(); ++i) {
            for (std::size_t j = 0; j < get<gtotal_ij>(Q_dtb)().ncols(); ++j) {

              std::apply(
                  [x, dt, &fo_ij, i, j, order](auto const &...g_ij) {
                    ((fo_ij << x.value() << "," << dt << ","
                            << "Qdt_bis"
                            << "," << order << "," << className(g_ij) << ","
                            << i << "," << j << "," << g_ij()(i, j) << "\n"),
                     ...);
                  },
                  std::forward_as_tuple(
                      get<P>(Q_dtb), get<gtotal_ij>(Q_dtb),
                      get<gmean_ij>(Q_dtb), get<gtotal_sqr_ij>(Q_dtb),
                      get<gtotal_var_ij>(Q_dtb), get<gvar_ij>(Q_dtb)));
            }
          }
        }
      }
    }
  }

    constexpr bool analyze_and_test_Macror = false;
  if constexpr (analyze_and_test_Macror) {

    /// the idea is to run several Macror algorithms with their Micro_stochastic
    /// counterparts. the Micro_stochastic algorithm has several parameters that
    /// determine the sampling effort we want to see that the sampling effort is
    /// big enough. the dimensions of sampling effort are the following
    /// 1. warming up period
    /// 2. number of samples
    /// 3. betas:
    /// 3. a number per decade
    /// 3. b minimum beta
    /// 4. number of walkers

    std::string filename = "micro_stochastic_48_";
    std::ofstream fo_i(filename + ModelName +  time_now()+".txt");
    //  std::ofstream fo_ij(filename + ModelName + "_istate_jstate.txt");

    fo_i << "i_t"
         << ","
         << "algoritm"
         << ","
         << "number_of_samples"
         << ","
         << "interval"
         << ","
         << "n_points_per_decade"
         << ","
         << "stops_at"
         << ","
         << "variable"
         << ","
         << "i_state"
         << ","
         << "j_state"
         << ","
         << "value"
         << "\n";
    constexpr const auto averaging = uses_averaging_aproximation(2);
    constexpr const auto variance = uses_variance_aproximation(false);
    constexpr const auto variance_correction = uses_variance_correction_aproximation(false);

    std::size_t current_i = 0;
    std::size_t current_i_s = 0;

    auto save = [](auto &fo_i, std::size_t i_t, const std::string &algorithm,
                   std::size_t number_of_samples, double interval,
                   std::size_t n_points_per_decade, double stops_at,
                   auto const &patch) {
      auto r_PCov = get<P_Cov>(patch);
      auto r_Pmean = get<P_mean>(patch);
      auto r_logL = get<plogL>(patch);

      fo_i << i_t << "," << algorithm << "," << number_of_samples << ","
           << interval << "," << n_points_per_decade << "," << stops_at << ","
           << "plogL"
           << "," << 0 << "," << 0 << "," << r_logL << "\n";
      for (std::size_t i = 0; i < r_PCov().nrows(); ++i)
        fo_i << i_t << "," << algorithm << "," << number_of_samples << ","
             << interval << "," << n_points_per_decade << "," << stops_at << ","
             << "Pmean"
             << "," << i << "," << 0 << "," << r_Pmean()[i] << "\n";
      for (std::size_t i = 0; i < r_PCov().nrows(); ++i)
        for (std::size_t j = i; j < r_PCov().ncols(); ++j)
          fo_i << i_t << "," << algorithm << "," << number_of_samples << ","
               << interval << "," << n_points_per_decade << "," << stops_at
               << ","
               << "PCov"
               << "," << i << "," << j << "," << r_PCov()(i, j) << "\n";
      fo_i.flush();
    };

    auto ftbl2 = insert(
        "micror_stoch_48", ftbl,
        var::Time_it(
            F(MacroR<uses_recursive_aproximation(true), averaging, variance>{},
              [&fo_i, &current_i, &current_i_s, &mt, &save, averaging,
               variance](auto &ft, Patch_State &&t_prior, Qdt const &t_Qdt,
                         Patch_Model const &mo, double const &Nch,
                         const Patch_current &p_y, double fs) {
                std::vector<std::size_t> sampled_i = { 120,130,135};
                if (current_i == sampled_i[current_i_s]) {
                  ++current_i_s;
                    std::vector<std::size_t> v_number_of_samples = { 1000};
                  std::vector<double> calculation_intervals = {0.01,0.02,0.05,0.1,0.2,0.5,  0.75,
                                                                1};
                  std::size_t save_every_iter = 50;
                  std::vector<std::size_t> v_n_points_per_decade = {12};
                  std::vector<double> v_stops_at = {1e-4};

                  std::uniform_int_distribution<
                      typename std::mt19937_64::result_type>
                      useed;

                  std::string algorithm=ToString(MicroR<averaging, variance>{});

                  for (auto number_of_samples : v_number_of_samples) {
                    for (auto n_points_per_decade : v_n_points_per_decade)
                      for (auto stops_at : v_stops_at) {

                          /*
                           * FunctionTable &ftbl, C_Patch_State &&t_prior,
                  C_Qdt const &t_Qdt, C_Patch_Model const &m,
                  C_double const &Nch, const Patch_current &p_y, double fs,
                  std::size_t myseed,std::size_t number_of_samples,std::vector<double> calculation_intervals,
                  std::size_t save_every_iter, std::size_t n_points_per_decade, double stops_at
                           * */

                           auto myseed = useed(mt);

                          auto micro =
                            Micror_stochastic<averaging, variance>(
                                ft, t_prior, t_Qdt, mo, Nch, p_y, fs, myseed,
                                number_of_samples, calculation_intervals,
                                save_every_iter, n_points_per_decade, stops_at)
                                ;
                           std::cerr<<"done\n"<<number_of_samples<<" npoint "<<n_points_per_decade<<" sto at"<<stops_at<<"\n";
                        for (std::size_t i_interval = 0;
                             i_interval < micro.size(); ++i_interval) {
                          save(fo_i, current_i, algorithm, number_of_samples,
                               calculation_intervals[i_interval],
                               n_points_per_decade, stops_at,
                               micro[i_interval]);
                        }
                      }
                  }
                }
                auto m = Macro_DMR{};
                auto out_Macro = m.Macror<uses_recursive_aproximation(true),
                                          averaging, variance,variance_correction>(
                    ft, std::move(t_prior), t_Qdt, mo, Nch, p_y, fs);
                save(fo_i, current_i,
                     ToString(MacroR<uses_recursive_aproximation(true),
                                     averaging, variance>{}),
                     0, 0, 0, 0, out_Macro.value());
                ++current_i;
                return m.Macror<uses_recursive_aproximation(true),
                                averaging, variance,variance_correction>(
                    ft, std::move(t_prior), t_Qdt, mo, Nch, p_y, fs);
              }),
            num_scouts_per_ensemble / 2));
    
    auto sim_7 = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
        recording_7);
    auto lik_7 = Macro_DMR{}
                   .log_Likelihood<uses_adaptive_aproximation(false),
                                   uses_recursive_aproximation(true), averaging,
                                   variance, variance_correction,return_predictions(false)>(
                        ftbl2.fork(var::I_thread(0)), model0, param1, experiment_7,
                       sim_7.value()());
  }

  constexpr bool test_derivative = false;

  if constexpr (test_derivative) {
    auto number_replicates = 100;
    auto outputfilename = "../macro_dr/output";

    std::string algorithm = "_new_MacroRC100_log_N100";
    std::ofstream fo(outputfilename + algorithm + ".txt");

    Matrix<double> mean_dlogL;
    SymPosDefMatrix<double> Cov_dlogL;

    Matrix<double> mean_edlogL;
    SymPosDefMatrix<double> Cov_edlogL;

    std::string path = "derivative";
    auto num_scouts_per_ensemble = 2ul;

    for (std::size_t i = 0; i < number_replicates; ++i) {
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
          recording);
      auto lik = Macro_DMR{}
                     .log_Likelihood<uses_adaptive_aproximation(false),
                                     uses_recursive_aproximation(true),
                                     uses_averaging_aproximation(2),
                                     uses_variance_aproximation(false),
                                     uses_variance_correction_aproximation(false),
                                     return_predictions(false)>(
                         ftbl.fork(var::I_thread(0)), model0, dparam1,
                         experiment, sim.value()());
      // std::cerr << "\n" << i << "th likelihood!!\n" << lik;
      if (lik) {
        auto v_logL = get<logL>(lik.value());
        auto v_elogL = get<elogL>(lik.value());
        auto v_vlogL = get<vlogL>(lik.value());

        auto v_dlogL = derivative(v_logL);
        auto v_delogL = derivative(v_elogL);

        if (i == 0) {
          fo << "logL,elogL,vlogL";
          for (std::size_t j = 0; j < v_dlogL().size(); ++j)
            fo << ","
               << "dlogL_d" << param1Names[j];
          // for (std::size_t j= 0; j<v_delogL().size(); ++j)
          //     fo<<","<<"delogL_d"<<param1Names[j];
          fo << "\n";
          mean_dlogL =
              Matrix<double>(v_dlogL().nrows(), v_dlogL().ncols(), 0.0);
          Cov_dlogL =
              SymPosDefMatrix<double>(v_dlogL().size(), v_dlogL().size(), 0.0);

          mean_edlogL =
              Matrix<double>(v_delogL().nrows(), v_delogL().ncols(), 0.0);
          Cov_edlogL = SymPosDefMatrix<double>(v_delogL().size(),
                                               v_delogL().size(), 0.0);
        }
        fo << primitive(v_logL) << ",";
        fo << primitive(v_elogL) << ",";
        fo << primitive(v_vlogL);
        for (std::size_t j = 0; j < v_dlogL().size(); ++j)
          fo << "," << v_dlogL()[j];
        // for (std::size_t j= 0; j<v_delogL().size(); ++j)
        //     fo<<"\t"<<v_delogL()[j];
        fo << "\n";
        mean_dlogL = mean_dlogL + v_dlogL();
        Cov_dlogL = Cov_dlogL + XXT(v_dlogL());
        mean_edlogL = mean_edlogL + v_delogL();
        Cov_edlogL = Cov_edlogL + XXT(v_delogL());
      }
    }
    mean_dlogL = mean_dlogL * (1.0 / number_replicates);
    std::ofstream foa(outputfilename + algorithm + "ana.txt");

    std::cerr << "\nmean_dlogL\n" << mean_dlogL << "\n";
    foa << "\nmean_dlogL\n" << mean_dlogL << "\n";
    Cov_dlogL = Cov_dlogL * (1.0 / number_replicates) - XXT(mean_dlogL);
    std::cerr << "\nCov_dlogL\n" << Cov_dlogL << "\n";
    foa << "\nCov_dlogL\n" << Cov_dlogL << "\n";

    auto Cov_inv = inv(Cov_dlogL);
    std::cerr << "\nCov_inv\n" << Cov_inv;
    foa << "\nCov_inv\n" << Cov_inv;
    if (Cov_inv) {
      std::cerr << "\nparameters\n" << param1() << "\n";

      auto accuracy = mean_dlogL * inv(Cov_dlogL).value();
      auto sensitivity =
          apply([](auto x) { return std::sqrt(x); }, diag(Cov_inv.value()));
      std::cerr << "\naccuracy\n" << accuracy << "\n";
      std::cerr << "\nsensitivityy\n" << sensitivity << "\n";

      std::cerr << "\naccuracy rate\n" << elemDiv(accuracy, param1()) << "\n";
      std::cerr << "\nsensitivityy\n"
                << sensitivity * inv(diag(param1())).value() << "\n";
      std::cerr << "\naccuracy  over se \n"
                << accuracy *
                       inv(sensitivity * (1.0 / std::sqrt(number_replicates)))
                << "\n";

      foa << "\naccuracy\n" << accuracy << "\n";
      foa << "\nsensitivityy\n" << sensitivity << "\n";

      foa << "\naccuracy rate\n" << elemDiv(accuracy, param1()) << "\n";
      foa << "\nsensitivity rate\n"
          << sensitivity * inv(diag(param1())).value() << "\n";

      foa << "\naccuracy  over se \n"
          << accuracy * inv(sensitivity * (1.0 / std::sqrt(number_replicates)))
          << "\n";
    }
  }

  constexpr bool test_simulation_model = false;
  if constexpr (test_simulation_model) {
    auto &model0 = model4;
    auto &param1Names = model0.names();
    auto &param1 = model0.parameters();
    std::string ModelName = "model4";
    using MyModel = Allost1;

    auto m = model0(param1);

    save(ModelName, m);

    std::string fname = "simulate_model";
    std::size_t nrep = 20;
    std::vector<Simulated_Recording<includes_N_state_evolution(false)>> out(20);
    for (std::size_t i = 0; i < nrep; ++i) {
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
          recording);
      out[i] = sim.value();
    }
    save(ModelName, out);
  }

  constexpr bool test_derivative_per_algorithm = false;

  if constexpr (test_derivative_per_algorithm) {
    auto number_replicates = 100;
    auto outputfilename = "test_derivative_" + ModelName;

    std::string fname = "_new_MacroRC100_log_N100";

    std::vector<std::string> algo = {"NR", "R", "VR", "aNR", "aR", "aVR"};

    std::vector<std::ofstream> fo;
    for (auto e : algo) {
      fo.push_back(std::ofstream(outputfilename + e + fname + ".txt"));
    }

    std::vector<Matrix<double>> mean_dlogL(algo.size());
    std::vector<SymPosDefMatrix<double>> Cov_dlogL(algo.size());

    // std::vector<Matrix<double>> mean_edlogL(algo.size());
    // std::vector<SymPosDefMatrix<double>> Cov_edlogL(algo.size());

    std::string path = "derivative";
    auto ftb = FuncMap(
        ModelName,
        Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                  thermo_cuevi_randomized_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                       step_stretch_cuevi_mcmc_per_walker{}),
                     num_scouts_per_ensemble / 2),

        var::Time_it(F(logLikelihood_f{},
                       [](auto &&...x) {
                         return logLikelihood(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(true)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(true),
                                         uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(
            F(MacroR<uses_recursive_aproximation(false),
                     uses_averaging_aproximation(2),
                     uses_variance_aproximation(false)>{},
              [](auto &&...x) {
                auto m = Macro_DMR{};
                return m.template Macror<uses_recursive_aproximation(false),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false)>(
                    std::forward<decltype(x)>(x)...);
              }),
            num_scouts_per_ensemble / 2),
        // var::Thread_Memoizer(
        //     var::F(Calc_Qdt_step{},
        //            [](auto &&...x) {
        //              auto m = Macro_DMR{};
        //                auto bisection_order=16ul;
        //              return m.calc_Qdt_bisection(
        //                  std::forward<decltype(x)>(x)...,bisection_order);
        //            }),
        //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
        //     num_scouts_per_ensemble / 2),
        var::Time_it(
            var::F(Calc_Qdt_step{},
                   [](auto &&...x) {
                     auto m = Macro_DMR{};
                     return m.calc_Qdt_ATP_step(
                         std::forward<decltype(x)>(x)...);
                   }),
            // var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
            num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_Qx{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qx(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_eigen{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_eigen(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2)

    );
    auto logLik_by_algo = [&ftb, &model0, &dparam1,
                           &experiment](auto const &sim, std::string algo) {
      auto ftbl = ftb.fork(var::I_thread(0ul));
      if (algo == "NR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(false),
                            uses_recursive_aproximation(false),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "R")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(false),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "VR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(false),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(true),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "aNR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(false),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "aR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else // if(algo=="aVR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(true),
                              uses_variance_correction_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
    };

    for (std::size_t i = 0; i < number_replicates; ++i) {
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
          recording);
      for (std::size_t j = 0; j < algo.size(); ++j) {
        auto lik = logLik_by_algo(sim, algo[j]);
        std::cerr << "\n"
                  << i << "th likelihood for " << algo[j] << " !!\n"
                  << lik;
        if (lik) {
          auto v_logL = get<logL>(lik.value());
          auto v_elogL = get<elogL>(lik.value());
          auto v_vlogL = get<vlogL>(lik.value());

          auto v_dlogL = derivative(v_logL);
          //   auto v_delogL = derivative(v_elogL);

          if (i == 0) {
            fo[j] << "logL,elogL,vlogL";
            for (std::size_t jj = 0; jj < v_dlogL().size(); ++jj)
              fo[j] << ","
                    << "dlogL_d" << param1Names[jj];
            // for (std::size_t j= 0; j<v_delogL().size(); ++j)
            //     fo<<","<<"delogL_d"<<param1Names[j];
            fo[j] << "\n";
            mean_dlogL[j] =
                Matrix<double>(v_dlogL().nrows(), v_dlogL().ncols(), 0.0);
            Cov_dlogL[j] = SymPosDefMatrix<double>(v_dlogL().size(),
                                                   v_dlogL().size(), 0.0);

            //   mean_edlogL[j] = Matrix<double>(v_delogL().nrows(),
            //   v_delogL().ncols(), 0.0); Cov_edlogL[j] =
            //   SymPosDefMatrix<double>(v_delogL().size(),
            //                                          v_delogL().size(), 0.0);
          }
          fo[j] << primitive(v_logL) << ",";
          fo[j] << primitive(v_elogL) << ",";
          fo[j] << primitive(v_vlogL);
          for (std::size_t jj = 0; jj < v_dlogL().size(); ++jj)
            fo[j] << "," << v_dlogL()[jj];
          // for (std::size_t j= 0; j<v_delogL().size(); ++j)
          //     fo<<"\t"<<v_delogL()[j];
          fo[j] << "\n";
          mean_dlogL[j] = mean_dlogL[j] + v_dlogL();
          Cov_dlogL[j] = Cov_dlogL[j] + XXT(v_dlogL());
          // mean_edlogL[j] = mean_edlogL[j] + v_delogL();
          // Cov_edlogL[j] = Cov_edlogL[j] + XXT(v_delogL());
        }
      }
    }
    for (std::size_t j = 0; j < algo.size(); ++j) {
      std::cerr << "\n......... algo " << algo[j] << "...........\n"
                << mean_dlogL[j] << "\n";

      mean_dlogL[j] = mean_dlogL[j] * (1.0 / number_replicates);
      std::ofstream foa(outputfilename + fname + algo[j] + "ana.txt");

      std::cerr << "\nmean_dlogL\n" << mean_dlogL[j] << "\n";
      foa << "\nmean_dlogL\n" << mean_dlogL[j] << "\n";
      Cov_dlogL[j] =
          Cov_dlogL[j] * (1.0 / number_replicates) - XXT(mean_dlogL[j]);
      std::cerr << "\nCov_dlogL\n" << Cov_dlogL[j] << "\n";
      foa << "\nCov_dlogL\n" << Cov_dlogL[j] << "\n";

      auto Cov_inv = inv(Cov_dlogL[j]);
      std::cerr << "\nCov_inv\n" << Cov_inv;
      foa << "\nCov_inv\n" << Cov_inv;
      if (Cov_inv) {
        std::cerr << "\nparameters\n" << param1() << "\n";

        auto accuracy = mean_dlogL[j] * inv(Cov_dlogL[j]).value();
        auto sensitivity =
            apply([](auto x) { return std::sqrt(x); }, diag(Cov_inv.value()));
        std::cerr << "\naccuracy\n" << accuracy << "\n";
        std::cerr << "\nsensitivityy\n" << sensitivity << "\n";

        std::cerr << "\naccuracy rate\n" << elemDiv(accuracy, param1()) << "\n";
        std::cerr << "\nsensitivityy\n"
                  << sensitivity * inv(diag(param1())).value() << "\n";
        std::cerr << "\naccuracy  over se \n"
                  << accuracy *
                         inv(sensitivity * (1.0 / std::sqrt(number_replicates)))
                  << "\n";

        foa << "\naccuracy\n" << accuracy << "\n";
        foa << "\nsensitivityy\n" << sensitivity << "\n";

        foa << "\naccuracy rate\n" << elemDiv(accuracy, param1()) << "\n";
        foa << "\nsensitivity rate\n"
            << sensitivity * inv(diag(param1())).value() << "\n";

        foa << "\naccuracy  over se \n"
            << accuracy *
                   inv(sensitivity * (1.0 / std::sqrt(number_replicates)))
            << "\n";
      }
    }
  }

  constexpr bool thermo_int_by_max_iter = false;

  if constexpr (thermo_int_by_max_iter) {
    /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
    // auto myseed = 9762841416869310605ul;
    auto myseed = 0;

    myseed = calc_seed(myseed);
    std::cerr << "myseed =" << myseed;

    /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
    std::size_t num_scouts_per_ensemble = 16;
    /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
    std::size_t max_num_simultaneous_temperatures = 4;

    /**

* @brief n_points_per_decade number of points per 10 times increment in beta
thermodynamic parameter
*/
    double n_points_per_decade = 3;

    /**
     * @brief stops_at minimum value of beta greater than zero
     */
    double stops_at = 1e-3;

    /**
     * @brief includes_zero considers also beta equal zero
     */
    bool includes_zero = true;

    /**
     * @brief max_iter maximum number of iterations on each warming step
     */
    std::size_t max_iter_warming = 200;

    /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
    std::size_t max_iter_equilibrium = 50000;

    /**
     * @brief path directory for the output
     */
    std::string path = "";

    /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
    std::size_t thermo_jumps_every = param1().size() * 10e0;

    double prior_error = 2;

    auto param1_prior = var::prior_around(param1, prior_error);

    /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
    auto tmi = thermo_Model_by_max_iter<MyModel>(
        path, "newNaN_thermo_R__2000", num_scouts_per_ensemble,
        max_num_simultaneous_temperatures, thermo_jumps_every, max_iter_warming,
        max_iter_equilibrium, n_points_per_decade, stops_at, includes_zero,
        myseed);

    auto modelLikelihood = make_Likelihood_Model<
        uses_adaptive_aproximation(false), uses_recursive_aproximation(true),
        uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)
>(
        model0, Simulation_n_sub_dt(1000ul));

    auto sim = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
        recording);
    auto ftbl = FuncMap(
        path,
        Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                  thermo_cuevi_randomized_jump_mcmc{}),
                num_scouts_per_ensemble / 2),

        var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                       step_stretch_cuevi_mcmc_per_walker{}),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(logLikelihood_f{},
                       [](auto &&...x) {
                         return logLikelihood(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(false),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_Qdt{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_Qx{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qx(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2),
        var::Time_it(F(Calc_eigen{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_eigen(std::forward<decltype(x)>(x)...);
                       }),
                     num_scouts_per_ensemble / 2)

    );

    if (sim)
      auto opt = evidence(ftbl, std::move(tmi), param1_prior, modelLikelihood,
                          sim.value()(), experiment);
  }
  
  
  
  constexpr bool cuevi_by_max_iter = false;
  if (cuevi_by_max_iter) {
      /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
      //   auto myseed = 9762841416869310605ul;
      //    auto myseed = 2555984001541913735ul;
      auto myseed = 0ul;
      
      myseed = calc_seed(myseed);
      std::cerr << "myseed =" << myseed << "\n";
      
      /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
      std::size_t num_scouts_per_ensemble = 16;
      
      /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
      std::size_t max_num_simultaneous_temperatures = 1e5;
      
      /**
     * @brief stops_at minimum value of beta greater than zero
     */
      double stops_at = 0.0001;
      
      /**
     * @brief includes_zero considers also beta equal zero
     */
      bool includes_zero = true;
      
      /**
     * @brief max_iter maximum number of iterations on each warming step
     */
      std::size_t max_iter_warming = 50;
      
      /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
      std::size_t max_iter_equilibrium = 8000;
      
      /**
     * @brief path directory for the output
     */
      std::string path = "";
      
      /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
      double min_fraction = 2;
      
      /**
     * @brief checks_derivative_every_model_size number of steps before every
     * check of the derivative against the beta thermo parameter for stopping
     */
      std::size_t checks_derivative_every_model_size = 10;
      
      /**
     * @brief max_ratio maximimum tolerated ratio for the beta derivative method
     */
      double max_ratio = 8000e16;
      
      /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */
      
      double n_points_per_decade = 6;
      /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
      double n_points_per_decade_fraction = 6;
      
      /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
      std::size_t thermo_jumps_every = param1().size() / 4;
      
      double prior_error = 2;
      
      auto param1_prior = var::prior_around(param1, prior_error);
      
      // auto& param1_prior = prior_model00_7;
      
      /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
      
      auto modelLikelihood = make_Likelihood_Model<
          uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
          uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
          model0, Simulation_n_sub_dt(1000ul));
      
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
          recording);
      
      if (sim) {
          std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
          auto number_of_traces = 7;
          auto number_of_segments = t_segments.size();
          t_segments.reserve(number_of_traces * t_segments.size());
          
          for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
              std::copy_n(t_segments.begin(), number_of_segments,
                          std::back_inserter(t_segments));
          std::cerr << "t_segments\n" << t_segments;
          std::cerr << "cum t_segments\n" << var::cumsum(t_segments);
          
          std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};
          
          std::size_t t_min_number_of_samples = 10;
          
          /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
          std::size_t bisection_count = 2ul;
          std::string filename_bisection =
              ModelName + "_bisection_" + std::to_string(bisection_count) + "_" +
              std::to_string(myseed) + "_" + time_now();
          
          bool all_at_once = true;
          
          std::string all_at_once_str =
              all_at_once ? "_all_at_once_" : "_progressive_";
          
          std::string n_points_per_decade_str =
              "_" + std::to_string(n_points_per_decade) + "_";
          
          std::string filename = ModelName + "_sim_eig_4800ch_MRAdap_only_7_" +
                                 all_at_once_str + "_randomized_jump_" +
                                 n_points_per_decade_str + time_now() + "_" +
                                 // std::to_string(bisection_count) + "_" +
                                 std::to_string(myseed);
          
          auto &t_segments_used = t_segments_7;
          
          auto cbc = cuevi_Model_by_iteration<MyModel>(
              path, filename, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed);
          
          // auto opt3 = evidence(std::move(cbc), param1_prior, modelLikelihood,
          //                      sim.value()(), experiment);
          auto ftbl3 = FuncMap(
              path + filename,
              Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                             step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          auto lik = Macro_DMR{}
                         .log_Likelihood<uses_adaptive_aproximation(false),
                                         uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false),
                                         return_predictions(true)>(
                             ftbl3.fork(var::I_thread(0)), model0, param1,
                             experiment, sim.value()());
          report(filename+"_lik.csv",lik.value(),sim.value(), experiment);
          if (false)
              auto opt3 = evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                                   sim.value()(), experiment, all_at_once);
      }
  }
  
  constexpr bool new_cuevi_by_max_iter = true;
  if (new_cuevi_by_max_iter) {
     
      /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
      //   auto myseed = 9762841416869310605ul;
      //    auto myseed = 2555984001541913735ul;
      auto myseed = 0ul;
      
      myseed = calc_seed(myseed);
      std::cerr << "myseed =" << myseed << "\n";
      
      /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
      std::size_t num_scouts_per_ensemble = 16;
      
      /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
      std::size_t max_num_simultaneous_temperatures = 1e5;
      
      /**
     * @brief stops_at minimum value of beta greater than zero
     */
      double stops_at = 1e-9;
      
      /**
     * @brief includes_zero considers also beta equal zero
     */
      bool includes_zero = true;
      
      /**
     * @brief max_iter maximum number of iterations on each warming step
     */
      std::size_t max_iter_warming = 50;
      
      /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
      std::size_t max_iter_equilibrium = 50000;
      
      /**
     * @brief path directory for the output
     */
      std::string path = "";
      
      /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
      double min_fraction = 2;
      
      /**
     * @brief checks_derivative_every_model_size number of steps before every
     * check of the derivative against the beta thermo parameter for stopping
     */
      std::size_t checks_derivative_every_model_size = 10;
      
      /**
     * @brief max_ratio maximimum tolerated ratio for the beta derivative method
     */
      double max_ratio = 8000e16;
      
      /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */
      
      double n_points_per_decade = 3;
      /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
      double n_points_per_decade_fraction = 6;
      
      /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
      std::size_t thermo_jumps_every = param1().size() / 4;
      
      double prior_error = 2;
      
      auto param1_prior = var::prior_around(param1, prior_error);
      
      // auto& param1_prior = prior_model00_7;
      
      /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
      
      auto modelLikelihood = make_Likelihood_Model<
          uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
          uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
          model0, Simulation_n_sub_dt(1000ul));
      
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
          recording);
      
      if (sim) {
          std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
          auto number_of_traces = 7;
          auto number_of_segments = t_segments.size();
          t_segments.reserve(number_of_traces * t_segments.size());
          
          for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
              std::copy_n(t_segments.begin(), number_of_segments,
                          std::back_inserter(t_segments));
          std::cerr << "t_segments\n" << t_segments;
          std::cerr << "cum t_segments\n" << var::cumsum(t_segments);
          
          std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};
          
          std::size_t t_min_number_of_samples = 10;
          
          /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
          std::size_t bisection_count = 2ul;
          std::string filename_bisection =
              ModelName + "_bisection_" + std::to_string(bisection_count) + "_" +
              std::to_string(myseed) + "_" + time_now();
          
          bool all_at_once = true;
          
          std::string all_at_once_str =
              all_at_once ? "_all_at_once_" : "_progressive_";
          
          std::string n_points_per_decade_str =
              "_" + std::to_string(n_points_per_decade) + "_";
          
          std::string filename = ModelName + "_new_cuevi_sim_eig_4800ch_only_7_" +
                                 n_points_per_decade_str + time_now() + "_" +
                                 // std::to_string(bisection_count) + "_" +
                                 std::to_string(myseed);
          
          auto &t_segments_used = t_segments_7;
          
          auto saving_itervals=Saving_intervals(Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),Save_Likelihood_every(num_scouts_per_ensemble),Save_Parameter_every(num_scouts_per_ensemble*10),Save_Predictions_every(num_scouts_per_ensemble*20)));
          
          auto cbc = new_cuevi_Model_by_iteration<MyModel>(
              path, filename, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed, saving_itervals);
          
          // auto opt3 = evidence(std::move(cbc), param1_prior, modelLikelihood,
          //                      sim.value()(), experiment);
          auto ftbl3 = FuncMap(
              path + filename,
              Time_it(F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                             cuevi::step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          auto lik = Macro_DMR{}
                         .log_Likelihood<uses_adaptive_aproximation(false),
                                         uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false),
                                         return_predictions(true)>(
                             ftbl3.fork(var::I_thread(0)), model0, param1,
                             experiment, sim.value()());
          report(filename+"_lik.csv",lik.value(),sim.value(), experiment);
          if (true)
              auto opt3 = cuevi::evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                                          sim.value()(), experiment, cuevi::Init_seed(seed));
      }
  }
  
  
  constexpr bool test_partial_logL = false;
  if (test_partial_logL) {
    /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
    //   auto myseed = 9762841416869310605ul;
    //    auto myseed = 2555984001541913735ul;
    auto myseed = 0ul;

    myseed = calc_seed(myseed);
    std::cerr << "myseed =" << myseed << "\n";


    /**
     * @brief path directory for the output
     */
    std::string path = "";



    /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */

    auto modelLikelihood = make_Likelihood_Model<
        uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
        uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
        model0, Simulation_n_sub_dt(1000));

    auto sim = Macro_DMR{}.sample_N(
        mt, model0, param1, experiment,
        Simulation_Parameters(Simulation_n_sub_dt(1000)),
        recording);

    if (sim) {

      std::string filename = ModelName + "_partial_logL_1e-9" +
                              time_now() + "_" +
                             std::to_string(myseed);

     
      auto ftbl3 = FuncMap(
          path + filename,
          Time_it(F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{}),
                  num_scouts_per_ensemble / 2),
          Time_it(F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{}),
                  num_scouts_per_ensemble / 2),
          Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                    thermo_cuevi_randomized_jump_mcmc{}),
                  num_scouts_per_ensemble / 2),
          var::Time_it(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                         cuevi::step_stretch_cuevi_mcmc_per_walker{}),
                       num_scouts_per_ensemble / 2),
          var::Time_it(F(logLikelihood_f{},
                         [](auto &&...x) {
                           return logLikelihood(
                               std::forward<decltype(x)>(x)...);
                         }),
                       num_scouts_per_ensemble / 2),
          var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                uses_averaging_aproximation(2),
                                uses_variance_aproximation(true)>{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.Macror<uses_recursive_aproximation(true),
                                           uses_averaging_aproximation(2),
                                           uses_variance_aproximation(true),
                                           uses_variance_correction_aproximation(false)>(
                               std::forward<decltype(x)>(x)...);
                         }),
                       num_scouts_per_ensemble / 2),
          var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                uses_averaging_aproximation(2),
                                uses_variance_aproximation(false)>{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.Macror<uses_recursive_aproximation(true),
                                           uses_averaging_aproximation(2),
                                           uses_variance_aproximation(false),
                                           uses_variance_correction_aproximation(false)>(
                               std::forward<decltype(x)>(x)...);
                         }),
                       num_scouts_per_ensemble / 2),
          var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                uses_averaging_aproximation(2),
                                uses_variance_aproximation(false)>{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.Macror<uses_recursive_aproximation(false),
                                           uses_averaging_aproximation(2),
                                           uses_variance_aproximation(false),
                                           uses_variance_correction_aproximation(false)>(
                               std::forward<decltype(x)>(x)...);
                         }),
                       num_scouts_per_ensemble / 2),
          // var::Thread_Memoizer(
          //     var::F(Calc_Qdt_step{},
          //            [](auto &&...x) {
          //              auto m = Macro_DMR{};
          //                auto bisection_order=16ul;
          //              return m.calc_Qdt_bisection(
          //                  std::forward<decltype(x)>(x)...,bisection_order);
          //            }),
          //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
          //     num_scouts_per_ensemble / 2),
          var::Thread_Memoizer(
              var::F(Calc_Qdt_step{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt_ATP_step(
                           std::forward<decltype(x)>(x)...);
                     }),
              var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              num_scouts_per_ensemble / 2),
          // var::Time_it(
          //     var::F(Calc_Qdt_step{},
          //            [](auto &&...x) {
          //              auto m = Macro_DMR{};
          //              return
          //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
          //            })),

          var::F(Calc_Qdt{},
                 [](auto &&...x) {
                   auto m = Macro_DMR{};
                   return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                 }),
          F(Calc_Qx{},
            [](auto &&...x) {
              auto m = Macro_DMR{};
              return m.calc_Qx(std::forward<decltype(x)>(x)...);
            }),
          var::Thread_Memoizer(
              F(Calc_eigen{},
                [](auto &&...x) {
                  auto m = Macro_DMR{};
                  return m.calc_eigen(std::forward<decltype(x)>(x)...);
                }),
              var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
              num_scouts_per_ensemble / 2)
          // var::Time_it(
          //     F(Calc_eigen{},
          //       [](auto &&...x) {
          //           auto m = Macro_DMR{};
          //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
          //       }))

      );
      auto lik = Macro_DMR{}
                         .log_Likelihood<uses_adaptive_aproximation(false),
                                         uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false),
                                         return_predictions(true)>(
                             ftbl3.fork(var::I_thread(0)), model0, param1,
                             experiment, sim.value()());
      report(filename+"_lik.csv",lik.value(),sim.value(), experiment);
    }
  }
  
  constexpr bool cuevi_by_max_iter_cross_model = false;
  if (cuevi_by_max_iter_cross_model) {
      /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
      //   auto myseed = 9762841416869310605ul;
      //    auto myseed = 2555984001541913735ul;
      auto myseed = 0ul;
      
      myseed = calc_seed(myseed);
      std::cerr << "myseed =" << myseed << "\n";
      
      /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
      std::size_t num_scouts_per_ensemble = 16;
      
      /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
      std::size_t max_num_simultaneous_temperatures = 1e5;
      
      /**
     * @brief stops_at minimum value of beta greater than zero
     */
      double stops_at = 0.0001;
      
      /**
     * @brief includes_zero considers also beta equal zero
     */
      bool includes_zero = true;
      
      /**
     * @brief max_iter maximum number of iterations on each warming step
     */
      std::size_t max_iter_warming = 50;
      
      /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
      std::size_t max_iter_equilibrium = 50000;
      
      /**
     * @brief path directory for the output
     */
      std::string path = "";
      
      /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
      double min_fraction = 2;
      
      /**
     * @brief checks_derivative_every_model_size number of steps before every
     * check of the derivative against the beta thermo parameter for stopping
     */
      std::size_t checks_derivative_every_model_size = 10;
      
      /**
     * @brief max_ratio maximimum tolerated ratio for the beta derivative method
     */
      double max_ratio = 8000e16;
      
      /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */
      
      double n_points_per_decade = 6;
      /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
      double n_points_per_decade_fraction = 6;
      
      /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
      std::size_t thermo_jumps_every = param1().size() / 4;
      
      double prior_error = 2;
      
      auto param1_prior = var::prior_around(param1, prior_error);
      auto param_alt_prior = var::prior_around(param1_alt, prior_error);
      
      // auto& param1_prior = prior_model00_7;
      
      /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
      
      auto modelLikelihood = make_Likelihood_Model<
          uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
          uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
          model0, Simulation_n_sub_dt(1000));
      
      auto modelLikelihood_alt = make_Likelihood_Model<
          uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
          uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
          model0_alt, Simulation_n_sub_dt(1000));
      
      auto n_sub_dt=Simulation_n_sub_dt(10000);
      
      std::string min_dt_st="_n_sub_dt_"+std::to_string(n_sub_dt())+"_";
      
      auto sim_N = Macro_DMR{}.sample_N(
          mt, model0, param1, experiment,
          Simulation_Parameters(n_sub_dt),
          recording);
      
      auto sim=get<Recording>(sim_N.value()());
      
      auto sim_alt_N = Macro_DMR{}.sample_N(
          mt, model0_alt, param1_alt, experiment,
          Simulation_Parameters(n_sub_dt),
          recording);
      auto sim_alt=get<Recording>(sim_alt_N.value()());
      
      if (sim_N.valid()&&sim_alt_N.valid()) {
          std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
          auto number_of_traces = 7;
          auto number_of_segments = t_segments.size();
          t_segments.reserve(number_of_traces * t_segments.size());
          
          for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
              std::copy_n(t_segments.begin(), number_of_segments,
                          std::back_inserter(t_segments));
          std::cerr << "t_segments\n" << t_segments;
          std::cerr << "cum t_segments\n" << var::cumsum(t_segments);
          
          std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};
          
          std::size_t t_min_number_of_samples = 10;
          
          /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
          std::size_t bisection_count = 2ul;
          std::string filename_bisection =
              ModelName + "_bisection_" + std::to_string(bisection_count) + "_" +
              std::to_string(myseed) + "_" + time_now();
          
          bool all_at_once = true;
          
          std::string all_at_once_str =
              all_at_once ? "_all_at_once_" : "_progressive_";
          
          std::string n_points_per_decade_str =
              "_" + std::to_string(n_points_per_decade) + "_";
          
          std::string prior_error_str =
              "_prior_error_" + std::to_string(prior_error) + "_";
          
          
          
          std::string filename_0_0 ="sim_"+min_dt_st+ ModelName +"_model_"+ModelName+ 
                                 prior_error_str+n_points_per_decade_str + time_now() + "_" +
                                 // std::to_string(bisection_count) + "_" +
                                 std::to_string(myseed);
          std::string filename_0_alt ="sim_"+min_dt_st+ ModelName +"_model_"+ModelName_alt+ 
                                    prior_error_str+ n_points_per_decade_str + time_now() + "_" +
                                     // std::to_string(bisection_count) + "_" +
                                     std::to_string(myseed);
          std::string filename_alt_0 ="sim_"+min_dt_st+ ModelName_alt +"_model_"+ModelName+ 
                                     prior_error_str+n_points_per_decade_str + time_now() + "_" +
                                     // std::to_string(bisection_count) + "_" +
                                     std::to_string(myseed);
          std::string filename_alt_alt ="sim_"+min_dt_st+ ModelName_alt +"_model_"+ModelName_alt+ 
                                     prior_error_str+n_points_per_decade_str + time_now() + "_" +
                                     // std::to_string(bisection_count) + "_" +
                                     std::to_string(myseed);
          
          
          auto &t_segments_used = t_segments_7;
          
          auto cbc_0_0 = cuevi_Model_by_iteration<MyModel>(
              path, filename_0_0, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed);
          
          auto cbc_0_alt = cuevi_Model_by_iteration<MyModel_alt>(
              path, filename_0_alt, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed);
          auto cbc_alt_0 = cuevi_Model_by_iteration<MyModel>(
              path, filename_alt_0, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed);
          auto cbc_alt_alt = cuevi_Model_by_iteration<MyModel_alt>(
              path, filename_alt_alt, t_segments_used, t_min_number_of_samples,
              num_scouts_per_ensemble, max_num_simultaneous_temperatures,
              min_fraction, thermo_jumps_every, max_iter_warming,
              max_iter_equilibrium, max_ratio, n_points_per_decade,
              n_points_per_decade_fraction, stops_at, includes_zero, myseed);
          
          // auto opt3 = evidence(std::move(cbc), param1_prior, modelLikelihood,
          //                      sim.value()(), experiment);
          auto ftbl3_0_0 = FuncMap(
              path + filename_0_0,
              Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                             step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          auto ftbl3_0_alt = FuncMap(
              path + filename_0_alt,
              Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                             step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          auto ftbl3_alt_0 = FuncMap(
              path + filename_alt_0,
              Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                             step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          auto ftbl3_alt_alt = FuncMap(
              path + filename_alt_alt,
              Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                        thermo_cuevi_randomized_jump_mcmc{}),
                      num_scouts_per_ensemble / 2),
              var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                             step_stretch_cuevi_mcmc_per_walker{}),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(logLikelihood_f{},
                             [](auto &&...x) {
                                 return logLikelihood(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(true)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(true),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(true),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false)>{},
                             [](auto &&...x) {
                                 auto m = Macro_DMR{};
                                 return m.Macror<uses_recursive_aproximation(false),
                                                 uses_averaging_aproximation(2),
                                                 uses_variance_aproximation(false),
                                                 uses_variance_correction_aproximation(false)>(
                                     std::forward<decltype(x)>(x)...);
                             }),
                           num_scouts_per_ensemble / 2),
              // var::Thread_Memoizer(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //                auto bisection_order=16ul;
              //              return m.calc_Qdt_bisection(
              //                  std::forward<decltype(x)>(x)...,bisection_order);
              //            }),
              //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
              //     num_scouts_per_ensemble / 2),
              var::Thread_Memoizer(
                  var::F(Calc_Qdt_step{},
                         [](auto &&...x) {
                             auto m = Macro_DMR{};
                             return m.calc_Qdt_ATP_step(
                                 std::forward<decltype(x)>(x)...);
                         }),
                  var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                  num_scouts_per_ensemble / 2),
              // var::Time_it(
              //     var::F(Calc_Qdt_step{},
              //            [](auto &&...x) {
              //              auto m = Macro_DMR{};
              //              return
              //              m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
              //            })),
              
              var::F(Calc_Qdt{},
                     [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              F(Calc_Qx{},
                [](auto &&...x) {
                    auto m = Macro_DMR{};
                    return m.calc_Qx(std::forward<decltype(x)>(x)...);
                }),
              var::Thread_Memoizer(
                  F(Calc_eigen{},
                    [](auto &&...x) {
                        auto m = Macro_DMR{};
                        return m.calc_eigen(std::forward<decltype(x)>(x)...);
                    }),
                  var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{},
                  num_scouts_per_ensemble / 2)
              // var::Time_it(
              //     F(Calc_eigen{},
              //       [](auto &&...x) {
              //           auto m = Macro_DMR{};
              //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
              //       }))
              
              );
          
          auto lik_0_0 = Macro_DMR{}
                         .log_Likelihood<uses_adaptive_aproximation(false),
                                         uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false),
                                         uses_variance_correction_aproximation(false),
                                         return_predictions(true)>(
                             ftbl3_0_0.fork(var::I_thread(0)), model0, param1,
                             experiment, sim);
          report(filename_0_0+"_lik.csv",lik_0_0.value(),sim_N.value(), experiment);
          auto lik_1_1 = Macro_DMR{}
                             .log_Likelihood<uses_adaptive_aproximation(false),
                                             uses_recursive_aproximation(true),
                                             uses_averaging_aproximation(2),
                                             uses_variance_aproximation(false),
                                             uses_variance_correction_aproximation(false),
                                             return_predictions(true)>(
                                 ftbl3_alt_alt.fork(var::I_thread(0)), model0_alt, param1_alt,
                                 experiment, sim_alt);
          report(filename_alt_alt+"_lik.csv",lik_1_1.value(),sim_alt_N.value(), experiment);
          
          if (true)
          {
#pragma omp parallel for
          for (auto i=0; i<4; ++i)
          {
              if (i==0){
                                    
                  auto opt3 = evidence(ftbl3_0_0, std::move(cbc_0_0), param1_prior, modelLikelihood,
                               sim, experiment, all_at_once);
              }
              else if (i==1)
              {   auto opt3 = evidence(ftbl3_0_alt, std::move(cbc_0_alt), param_alt_prior, modelLikelihood_alt,
                                       sim, experiment, all_at_once);
              }else if (i==2)
              {   auto opt3 = evidence(ftbl3_alt_0, std::move(cbc_alt_0), param1_prior, modelLikelihood,
                                       sim_alt, experiment, all_at_once);
              }else if (i==3)
              {   auto opt3 = evidence(ftbl3_alt_alt, std::move(cbc_alt_alt), param_alt_prior, modelLikelihood_alt,
                                       sim_alt, experiment, all_at_once);
              }}
          
      }
      }
  }
  
  
  
  if (false) {
    std::string ModelName = "test_der_Likelihood";
    std::string path = "";
    auto num_scouts_per_ensemble = 2ul;
    auto myseed = 0ul;

    auto sim = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Simulation_n_sub_dt(1000)));

    auto test_der_Likelihood = var::test_Derivative(
        [&model0, &sim, &experiment, &ftbl](auto const &dparam1) {
          return Macro_DMR{}
              .log_Likelihood<uses_adaptive_aproximation(false),
                              uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false),
                                uses_variance_correction_aproximation(false),
                              return_predictions(false)>(
                  ftbl.fork(var::I_thread(0)), model0, dparam1, experiment,
                  sim.value()());
        },
        1, 1e-10, dparam1);
    if (!test_der_Likelihood) {
      std::cerr << test_der_Likelihood.error()();
    }
  }
}
