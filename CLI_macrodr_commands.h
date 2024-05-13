#ifndef CLI_MACRODR_COMMANDS_H
#define CLI_MACRODR_COMMANDS_H

#include "experiment.h"
#include "function_measure_verification_and_optimization.h"
#include "function_memoization.h"
#include "matrix.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "random_samplers.h"
#include "variables.h"
// #include "multivariate_normal_distribution.h"
#include "cuevi.h"
// #include "models.h"
#include "models_MoffattHume_linear.h"
#include "qmodel.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iterator>
#include <string>
#include <vector>

namespace macrodr {

namespace cmd {

namespace deprecated {

namespace deprecated {
using namespace var::deprecated;
inline auto get_function_Table(std::string filename,
                               std::size_t num_scouts_per_ensemble) {
  using namespace macrodr;
  using namespace var::deprecated;

  return var::deprecated::FuncMap(
      filename,
      Time_it(
          F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{}),
          num_scouts_per_ensemble / 2),
      Time_it(
          F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{}),
          num_scouts_per_ensemble / 2),
      Time_it(F(::deprecated::thermo_cuevi_randomized_jump_mcmc{},
                ::deprecated::thermo_cuevi_randomized_jump_mcmc{}),
              num_scouts_per_ensemble / 2),
      Time_it(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                cuevi::step_stretch_cuevi_mcmc_per_walker{}),
              num_scouts_per_ensemble / 2),
      Time_it(F(logLikelihood_f{},
                [](auto &&...x) {
                  return logLikelihood(std::forward<decltype(x)>(x)...);
                }),
              num_scouts_per_ensemble / 2),
      Time_it(F(MacroR<uses_recursive_aproximation(true),
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
      Time_it(F(MacroR<uses_recursive_aproximation(true),
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
      Time_it(F(MacroR<uses_recursive_aproximation(false),
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
      //            [](auto &&.   ..x) {
      //              auto m = Macro_DMR{};
      //                auto bisection_order=16ul;
      //              return m.calc_Qdt_bisection(
      //                  std::forward<decltype(x)>(x)...,bisection_order);
      //            }),
      //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
      //     num_scouts_per_ensemble / 2),
      Thread_Memoizer(
          var::F(Calc_Qdt_step{},
                 [](auto &&...x) {
                   auto m = Macro_DMR{};
                   return m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
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
      Thread_Memoizer(
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
}
} // namespace deprecated

} // namespace deprecated


namespace deprecated {
using namespace var::deprecated;

inline auto get_function_Table_maker(std::string filename,
                                     std::size_t num_scouts_per_ensemble) {
  using namespace macrodr;
  return [filename, num_scouts_per_ensemble]() {
    return var::deprecated::FuncMap(
        filename,
        Time_it(F(cuevi::step_stretch_cuevi_mcmc{},
                  cuevi::step_stretch_cuevi_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(
            F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{}),
            num_scouts_per_ensemble / 2),
        Time_it(F(::deprecated::thermo_cuevi_randomized_jump_mcmc{},
                  ::deprecated::thermo_cuevi_randomized_jump_mcmc{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                  cuevi::step_stretch_cuevi_mcmc_per_walker{}),
                num_scouts_per_ensemble / 2),
        Time_it(F(logLikelihood_f{},
                  [](auto &&...x) {
                    return logLikelihood(std::forward<decltype(x)>(x)...);
                  }),
                num_scouts_per_ensemble / 2),
        Time_it(
            F(MacroR<uses_recursive_aproximation(true),
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
        Time_it(
            F(MacroR<uses_recursive_aproximation(true),
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
        Time_it(
            F(MacroR<uses_recursive_aproximation(false),
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
        //            [](auto &&.   ..x) {
        //              auto m = Macro_DMR{};
        //                auto bisection_order=16ul;
        //              return m.calc_Qdt_bisection(
        //                  std::forward<decltype(x)>(x)...,bisection_order);
        //            }),
        //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
        //     num_scouts_per_ensemble / 2),
        Thread_Memoizer(
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
        Thread_Memoizer(
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
  };
}

} // namespace deprecated


/*
make_Likelihood_Model<uses_adaptive_aproximation(true),
                      uses_recursive_aproximation(true),
                      uses_averaging_aproximation(2),
                      uses_variance_aproximation(false),
                      uses_variance_correction_aproximation(
                          false)>

*/

  



// evidence(prior= prior_model, likelihoodModel= likelihood, data = simulation,
// experiment = experiment, algorithm= algorithm, function_table
// =function_table, init_seed =0)




/**
   * @brief fuck
   * @param num_scouts_per_ensemble number of scouts per ensemble in the
   * affine ensemble mcmc model
     number_trials_until_give_up when the number of parallel
   * temepratures reaches this number, it stops growing, the same scouts
   * drifts on temperature
   * @param number_trials_until_give_up for trials each sampled intial value
   * @param stops_at minimum value of beta greater than zero
   * @param medium_beta beta where the slope changes
   * @param
   * @param
   *
   */

/**
 * @brief stops_at
 */

template <class TableFunction>
inline auto calc_evidence_old(
    TableFunction ftbl3, typename mt_64i::result_type myseed,
    std::string Efilename_7 = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
    std::size_t num_scouts_per_ensemble = 16,
    std::size_t number_trials_until_give_up = 1e5, double stops_at = 1e-15,
    double medium_beta = 1e-2,
    /**
     * @brief includes_zero considers also beta equal zero
     */
    bool includes_zero = true,

    /**
     * @brief randomly tries thermodynamic jumps
     */
    bool random_jumps = true,

    /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
    std::size_t max_iter_equilibrium = 50000,

    /**
     * @brief path directory for the output
     */
    std::string path = "",

    /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
    double min_fraction = 4,
    /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */

    double n_points_per_decade = 1,
    /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
    double n_points_per_decade_fraction = 6,
    /**
     * @brief thermo_jump_factor factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
    double thermo_jump_factor = 0.25,

    /**
     * @brief prior_error error around the prior
     *
     */

    double prior_error = 2) {
  using namespace macrodr;

  myseed = calc_seed(myseed);
  mt_64i mt(myseed);

  auto [recording_conditions_7, recording_7] =
      macrodr::load_recording(Efilename_7);

  auto experiment_7 =
      Experiment(std::move(recording_conditions_7), Frequency_of_Sampling(50e3),
                 initial_ATP_concentration(ATP_concentration(0.0)));
  auto &experiment = experiment_7;
  auto &recording = recording_7;

  auto &model0 = model6_Eff_no_inactivation;
  auto &param1 = model0.parameters();
  std::string ModelName = "model6_Eff_no_inactivation";
  using MyModel = Allost1;

  std::size_t thermo_jumps_every = param1().size() * thermo_jump_factor;

  auto param1_prior = var::prior_around(param1, prior_error);

  auto modelLikelihood = make_Likelihood_Model<
      uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
      uses_averaging_aproximation(2), uses_variance_aproximation(false),
      uses_variance_correction_aproximation(false)>(
      model0, Simulation_n_sub_dt(1000ul));

  auto sim = Macro_DMR{}.sample(
      mt, model0, param1, experiment,
      Simulation_Parameters(Simulation_n_sub_dt(1000ul)), recording);

  if (sim) {
    std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
    auto number_of_traces = 7;
    auto number_of_segments = t_segments.size();
    t_segments.reserve(number_of_traces * t_segments.size());

    for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
      std::copy_n(t_segments.begin(), number_of_segments,
                  std::back_inserter(t_segments));

    std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};

    bool average_the_ATP_evolution = true;

    /**
     * @brief cbc cumulative evidence algorithm, ends using convergence
     * criteria
     */
    std::size_t bisection_count = 2ul;
    std::string filename_bisection = ModelName + "_bisection_" +
                                     std::to_string(bisection_count) + "_" +
                                     std::to_string(myseed) + "_" + time_now();

    std::string n_points_per_decade_str =
        "_" + std::to_string(n_points_per_decade) + "_";

    std::string filename = ModelName + "_new_cuevi_sim_eig_4800ch_only_7_" +
                           n_points_per_decade_str + time_now() + "_" +
                           // std::to_string(bisection_count) + "_" +
                           std::to_string(myseed);

    auto &t_segments_used = t_segments_7;

    auto saving_itervals = Saving_intervals(
        Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),
                     Save_Likelihood_every(num_scouts_per_ensemble),
                     Save_Parameter_every(num_scouts_per_ensemble),
                     Save_Predictions_every(num_scouts_per_ensemble * 50)));

    auto cbc = new_cuevi_Model_by_iteration(
        path, filename, t_segments_used, average_the_ATP_evolution,
        num_scouts_per_ensemble, number_trials_until_give_up, min_fraction,
        thermo_jumps_every, max_iter_equilibrium, n_points_per_decade,
        n_points_per_decade_fraction, medium_beta, stops_at, includes_zero,
        saving_itervals, random_jumps);

    auto lik = Macro_DMR{}
                   .log_Likelihood<uses_adaptive_aproximation(false),
                                   uses_recursive_aproximation(true),
                                   uses_averaging_aproximation(2),
                                   uses_variance_aproximation(false),
                                   uses_variance_correction_aproximation(false),
                                   return_predictions(true)>(
                       ftbl3.fork(var::I_thread(0)), model0, param1, experiment,
                       sim.value()());
    report(filename + "_lik.csv", lik.value(), sim.value(), experiment);
    auto opt3 =
        cuevi::evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                        sim.value()(), experiment, cuevi::Init_seed(myseed));
    // auto opt4= cuevi::continue_evidence(filename, 2* maxiter);
  }
}
} // namespace cmd

} // namespace macrodr

#endif // CLI_MACRODR_COMMANDS_H
