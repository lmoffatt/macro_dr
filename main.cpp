#include "CLI_macro_dr.h"
#include "experiment.h"
#include "function_measure_verification_and_optimization.h"
#include "function_memoization.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "maybe_error.h"
#include "parameters_derivative.h"
#include "parameters_distribution.h"
#include "variables_derivative.h"
// #include "multivariate_normal_distribution.h"
#include "allosteric_models.h"
#include "cuevi.h"
#include "parallel_tempering.h"
#include "qmodel.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
using namespace macrodr;

inline std::string leadingZero(int i) {
  if (i == 0)
    return "00";
  else if (i < 10)
    return "0" + std::to_string(i);
  else
    return std::to_string(i);
}

inline std::string time_now() {
  auto tc = std::chrono::system_clock::now();
  std::time_t rawtime = std::chrono::system_clock::to_time_t(tc);

  auto tcount = (tc.time_since_epoch().count() / 1000) % 1000000;

  struct std::tm *t;
  time(&rawtime);
  t = localtime(&rawtime);
  return leadingZero(t->tm_hour) + leadingZero(t->tm_min) +
         leadingZero(t->tm_sec) + "s" + std::to_string(tcount);
}

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
        Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{})),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{})),
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
                                         uses_variance_aproximation(false)>(
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
                                         uses_variance_aproximation(false)>(
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
                       })),
        var::Time_it(F(Calc_eigen{}, [](auto &&...x) {
          auto m = Macro_DMR{};
          return m.calc_eigen(std::forward<decltype(x)>(x)...);
        })));

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

      auto opt3 = evidence(ftbl, std::move(cbc), my_linear_model.prior(),
                           my_linear_model.likelihood(), y, X);
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

  auto Efilename = "../macro_dr/Moffatt_Hume_2007_ATP_time.txt";

  auto [recording_conditions, recording] = macrodr::load_recording(Efilename);

  auto experiment =
      Experiment(std::move(recording_conditions), Frequency_of_Sampling(50e3),
                 initial_ATP_concentration(ATP_concentration(0.0)));

  struct Model0 : public Model_Patch<Model0> {};
  struct Model1 : public Model_Patch<Model1> {};

  struct Allost1 : public Model_Patch<Allost1> {};

  auto model00 = Model0::Model([]() {
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

    auto p_other = std::vector<double>{1, 1000, 1e-3, 1};

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
          auto v_g = p[16] * -1.0;
          auto Npar = 17ul;
          auto v_N0 = p[Npar];
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
              build<Current_Baseline>(v_baseline), Binomial_magical_number(5.0),
              min_P(1e-7), Probability_error_tolerance(1e-2),
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

  auto param11Names = std::vector<std::string>{
      "k01",         "k12",          "k23",   "k10",         "k21",
      "k32",         "k34",          "k43",   "conductance", "Num_Chan_0",
      "Num_Chan_eq", "Num_Chan_tau", "noise", "baseline"};

  auto &model0 = model6;
  auto &param1Names = param6Names;
  auto &param1 = param6;
  std::string ModelName = "Model6";
  using MyModel = Allost1;

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

  auto ggb = ggg(Time(10.0));

  // using jger=typename decltype(NNN)::kgerge;
  // using jger2=typename decltype(dNNN)::kgerge;

  // std::cerr<<"dparam1\n"<<dparam1;
  auto n = build<N_Ch_mean>(dparam1()[0]);

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

  auto ftbl = FuncMap(
      "_" + time_now(),
      Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{})),
      Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{})),
      var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                     step_stretch_cuevi_mcmc_per_walker{})),
      var::Time_it(F(logLikelihood_f{},
                     [](auto &&...x) {
                       return logLikelihood(std::forward<decltype(x)>(x)...);
                     })),
      var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false)>{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.Macror<uses_recursive_aproximation(true),
                                       uses_averaging_aproximation(2),
                                       uses_variance_aproximation(false)>(
                           std::forward<decltype(x)>(x)...);
                     })),
      var::Time_it(F(MacroR<uses_recursive_aproximation(false),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false)>{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.Macror<uses_recursive_aproximation(false),
                                       uses_averaging_aproximation(2),
                                       uses_variance_aproximation(false)>(
                           std::forward<decltype(x)>(x)...);
                     })),
      var::Time_it(F(Calc_Qdt{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     })),
      var::Time_it(F(Calc_Qx{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qx(std::forward<decltype(x)>(x)...);
                     })),
      var::Time_it(F(Calc_eigen{}, [](auto &&...x) {
        auto m = Macro_DMR{};
        return m.calc_eigen(std::forward<decltype(x)>(x)...);
      })));

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
          Simulation_Parameters(Number_of_simulation_sub_steps(100ul)),
          recording);
      auto lik = Macro_DMR{}
                     .log_Likelihood<uses_adaptive_aproximation(false),
                                     uses_recursive_aproximation(true),
                                     uses_averaging_aproximation(2),
                                     uses_variance_aproximation(false),
                                     return_predictions(false)>(
                         ftbl, model0, dparam1, experiment, sim.value()());
      std::cerr << "\n" << i << "th likelihood!!\n" << lik;
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
          mean_dlogL = Matrix<double>(1, v_delogL().size(), 0.0);
          Cov_dlogL = SymPosDefMatrix<double>(v_delogL().size(),
                                              v_delogL().size(), 0.0);

          mean_edlogL = Matrix<double>(1, v_delogL().size(), 0.0);
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
        Cov_dlogL = Cov_dlogL + XTX(v_dlogL());
        mean_edlogL = mean_edlogL + v_delogL();
        Cov_edlogL = Cov_edlogL + XTX(v_delogL());
      }
    }
    mean_dlogL = mean_dlogL * (1.0 / number_replicates);
    std::ofstream foa(outputfilename + algorithm + "ana.txt");

    std::cerr << "\nmean_dlogL\n" << mean_dlogL << "\n";
    foa << "\nmean_dlogL\n" << mean_dlogL << "\n";
    Cov_dlogL = Cov_dlogL * (1.0 / number_replicates) - XTX(mean_dlogL);
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

  constexpr bool test_derivative_per_algorithm = false;

  if constexpr (test_derivative_per_algorithm) {
    auto number_replicates = 100;
    auto outputfilename = "../macro_dr/output";

    std::string fname = "_new_MacroRC100_log_N100";

    std::vector<std::string> algo = {"NR", "R", "VR", "aNR", "aR", "aVR"};

    std::vector<std::ofstream> fo;
    for (auto e : algo) {
      fo.push_back(std::ofstream(outputfilename + e + fname + ".txt"));
    }

    std::vector<Matrix<double>> mean_dlogL(algo.size());
    std::vector<SymPosDefMatrix<double>> Cov_dlogL(algo.size());

    std::vector<Matrix<double>> mean_edlogL(algo.size());
    std::vector<SymPosDefMatrix<double>> Cov_edlogL(algo.size());

    std::string path = "derivative";
    auto ftb = FuncMap(
        path, Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{})),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{})),
        var::Time_it(F(step_stretch_cuevi_mcmc_per_walker{},
                       step_stretch_cuevi_mcmc_per_walker{})),

        var::Time_it(F(logLikelihood_f{},
                       [](auto &&...x) {
                         return logLikelihood(std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it(F(MacroR<uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(true)>{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.Macror<uses_recursive_aproximation(true),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(true)>(
                             std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it(
            F(MacroR<uses_recursive_aproximation(false),
                     uses_averaging_aproximation(2),
                     uses_variance_aproximation(false)>{},
              [](auto &&...x) {
                auto m = Macro_DMR{};
                return m.template Macror<uses_recursive_aproximation(false),
                                         uses_averaging_aproximation(2),
                                         uses_variance_aproximation(false)>(
                    std::forward<decltype(x)>(x)...);
              })),
        var::Time_it(F(Calc_Qdt{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it(F(Calc_Qx{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qx(std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it(F(Calc_eigen{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_eigen(std::forward<decltype(x)>(x)...);
                       }))

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
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "R")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(false),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "VR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(false),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(true),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "aNR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(false),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else if (algo == "aR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(false),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
      else // if(algo=="aVR")
        return Macro_DMR{}
            .log_Likelihood<uses_adaptive_aproximation(true),
                            uses_recursive_aproximation(true),
                            uses_averaging_aproximation(2),
                            uses_variance_aproximation(true),
                            return_predictions(false)>(
                ftbl, model0, dparam1, experiment, sim.value()());
    };

    for (std::size_t i = 0; i < number_replicates; ++i) {
      auto sim = Macro_DMR{}.sample(
          mt, model0, param1, experiment,
          Simulation_Parameters(Number_of_simulation_sub_steps(100ul)),
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
          auto v_delogL = derivative(v_elogL);

          if (i == 0) {
            fo[j] << "logL,elogL,vlogL";
            for (std::size_t jj = 0; jj < v_dlogL().size(); ++jj)
              fo[j] << ","
                    << "dlogL_d" << param1Names[jj];
            // for (std::size_t j= 0; j<v_delogL().size(); ++j)
            //     fo<<","<<"delogL_d"<<param1Names[j];
            fo[j] << "\n";
            mean_dlogL[j] = Matrix<double>(1, v_delogL().size(), 0.0);
            Cov_dlogL[j] = SymPosDefMatrix<double>(v_delogL().size(),
                                                   v_delogL().size(), 0.0);

            mean_edlogL[j] = Matrix<double>(1, v_delogL().size(), 0.0);
            Cov_edlogL[j] = SymPosDefMatrix<double>(v_delogL().size(),
                                                    v_delogL().size(), 0.0);
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
          Cov_dlogL[j] = Cov_dlogL[j] + XTX(v_dlogL());
          mean_edlogL[j] = mean_edlogL[j] + v_delogL();
          Cov_edlogL[j] = Cov_edlogL[j] + XTX(v_delogL());
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
          Cov_dlogL[j] * (1.0 / number_replicates) - XTX(mean_dlogL[j]);
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

  if (thermo_int_by_max_iter) {
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
    std::size_t max_iter_equilibrium = 10000;

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
        uses_averaging_aproximation(2), uses_variance_aproximation(false)>(
        model0, Number_of_simulation_sub_steps(10ul));

    auto sim = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Number_of_simulation_sub_steps(100ul)),
        recording);
    auto ftbl = FuncMap(
        path, Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{})),
        Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{})),
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
                                         uses_variance_aproximation(false)>(
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
                                         uses_variance_aproximation(false)>(
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
                       })),
        var::Time_it(F(Calc_eigen{},
                       [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_eigen(std::forward<decltype(x)>(x)...);
                       }))

    );

    if (sim)
      auto opt = evidence(ftbl, std::move(tmi), param1_prior, modelLikelihood,
                          sim.value()(), experiment);
  }

  constexpr bool cuevi_by_max_iter = true;
  if (cuevi_by_max_iter) {
    /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
    //   auto myseed = 9762841416869310605ul;
    auto myseed = 0ul;

    myseed = calc_seed(myseed);
    std::cerr << "myseed =" << myseed << "\n";

    /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
    std::size_t num_scouts_per_ensemble = 8;

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
    std::size_t max_iter_warming = 200;

    /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
    std::size_t max_iter_equilibrium = 40000;

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

    double n_points_per_decade = 24;
    /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
    double n_points_per_decade_fraction = 24;

    /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
    std::size_t thermo_jumps_every = param1().size() / 4;

    double prior_error = 2;

    auto param1_prior = var::prior_around(param1, prior_error);

    /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */

    auto modelLikelihood = make_Likelihood_Model<
        uses_adaptive_aproximation(false), uses_recursive_aproximation(true),
        uses_averaging_aproximation(2), uses_variance_aproximation(false)>(
        model0, Number_of_simulation_sub_steps(10ul));

    auto sim = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Number_of_simulation_sub_steps(100ul)),
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

      std::size_t t_min_number_of_samples = 10;

      /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
      std::size_t bisection_count = 2ul;
      std::string filename_bisection = ModelName + "_bisection_" +
                             std::to_string(bisection_count) + "_" +
                             std::to_string(myseed) + "_" + time_now();
      
      std::string filename = ModelName + "_memoization_" +
                                       std::to_string(bisection_count) + "_" +
                                       std::to_string(myseed) + "_" + time_now();

      auto cbc = cuevi_Model_by_iteration<MyModel>(
          path, filename, t_segments, t_min_number_of_samples,
          num_scouts_per_ensemble, max_num_simultaneous_temperatures,
          min_fraction, thermo_jumps_every, max_iter_warming,
          max_iter_equilibrium, max_ratio, n_points_per_decade,
          n_points_per_decade_fraction, stops_at, includes_zero, myseed);

      // auto opt3 = evidence(std::move(cbc), param1_prior, modelLikelihood,
      //                      sim.value()(), experiment);
      auto ftbl3 = FuncMap(
          path + filename,
          Time_it(F(step_stretch_cuevi_mcmc{}, step_stretch_cuevi_mcmc{})),
          Time_it(F(thermo_cuevi_jump_mcmc{}, thermo_cuevi_jump_mcmc{})),
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
                                uses_variance_aproximation(false)>{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.Macror<uses_recursive_aproximation(true),
                                           uses_averaging_aproximation(2),
                                           uses_variance_aproximation(false)>(
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
                                           uses_variance_aproximation(false)>(
                               std::forward<decltype(x)>(x)...);
                         }),
                       num_scouts_per_ensemble / 2),
          var::Time_it(var::Thread_Memoizer(
              var::F(Calc_Qdt_step{},
                     [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
                     }),
              var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{})),
           var::Time_it(
               var::F(Calc_Qdt{},
                      [](auto &&...x) {
                         auto m = Macro_DMR{};
                         return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                     }),
              num_scouts_per_ensemble / 2),
          var::Time_it(F(Calc_Qx{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.calc_Qx(std::forward<decltype(x)>(x)...);
                         })),
          var::Time_it(var::Thread_Memoizer(
                       F(Calc_eigen{},
                         [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.calc_eigen(std::forward<decltype(x)>(x)...);
                         }),              var::Memoiza_all_values<Maybe_error<Qx_eig>, Qx>{}))

      );

      auto opt3 = evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                           recording, experiment);
    }
  }

  if (false) {
    std::string ModelName = "test_der_Likelihood";
    std::string path = "";
    auto num_scouts_per_ensemble = 2ul;
    auto myseed = 0ul;

    auto sim = Macro_DMR{}.sample(
        mt, model0, param1, experiment,
        Simulation_Parameters(Number_of_simulation_sub_steps(100ul)));

    auto test_der_Likelihood = var::test_Derivative(
        [&model0, &sim, &experiment, &ftbl](auto const &dparam1) {
          return Macro_DMR{}
              .log_Likelihood<uses_adaptive_aproximation(false),
                              uses_recursive_aproximation(true),
                              uses_averaging_aproximation(2),
                              uses_variance_aproximation(false),
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
