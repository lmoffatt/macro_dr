#ifndef CLI_MACRODR_COMMANDS_H
#define CLI_MACRODR_COMMANDS_H

#include "CLI_macro_dr.h"
#include "distributions.h"
#include "experiment.h"
#include "function_measure_verification_and_optimization.h"
#include "function_memoization.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parameters.h"
#include "parameters_derivative.h"
#include "parameters_distribution.h"
#include "random_samplers.h"
#include "variables.h"
#include "variables_derivative.h"
// #include "multivariate_normal_distribution.h"
#include "allosteric_models.h"
#include "cuevi.h"
#include "micror_stochastic.h"
// #include "models.h"
#include "models_MoffattHume_linear.h"
#include "parallel_tempering.h"
#include "qmodel.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace macrodr {

inline auto get_random_id(std::string prefix) {

  return prefix + "_" + std::to_string(calc_seed(0ul));
}

inline void write_text(std::string filename, std::string s) {
  std::ofstream f(filename);
  f << s;
}

static constexpr auto temp_script_file = "temp_script_file.txt";
inline void write_script(std::string program_name) {
  rename(temp_script_file, program_name.c_str());
}

inline auto load_Parameter_value(const std::string filename, std::string separator) {
  return std::pair(filename, separator);
}

auto load_Prior_value(std::string filename, std::string separator = ",") {
  return std::pair(filename, separator);
}

inline auto load_Recording_value(std::string filename, std::string separator = ",") {
  return std::pair(filename, separator);
}

namespace deprecated{


namespace deprecated{
using namespace var::deprecated;
auto get_function_Table(std::string filename,
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
}

}
auto get_function_Table_St(std::string filename) {
    using namespace macrodr;
    return var::FuncMap_St(
        filename,
        Time_it_st(
            F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{})),
        Time_it_st(
            F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{})),
    //    Time_it_st(F(::deprecated::thermo_cuevi_randomized_jump_mcmc{},
    //              ::deprecated::thermo_cuevi_randomized_jump_mcmc{})),
        var::Time_it_st(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                       cuevi::step_stretch_cuevi_mcmc_per_walker{})),
        var::Time_it_st(F(logLikelihood_f{},
                       [](auto &&...x) {
                           return logLikelihood(std::forward<decltype(x)>(x)...);
                       })),
        var::Time_it_st(
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
              })),
        var::Time_it_st(
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
              })),
        var::Time_it_st(
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
              })),
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
        var::Single_Thread_Memoizer(
            var::F(Calc_Qdt_step{},
                   [](auto &&...x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt_ATP_step(std::forward<decltype(x)>(x)...);
                   }),
            var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{}),
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
        var::Single_Thread_Memoizer(
            F(Calc_eigen{},
              [](auto &&...x) {
                  auto m = Macro_DMR{};
                  return m.calc_eigen(std::forward<decltype(x)>(x)...);
              }),
            var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{})
        // var::Time_it(
        //     F(Calc_eigen{},
        //       [](auto &&...x) {
        //           auto m = Macro_DMR{};
        //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
        //       }))
        
        );
}



auto get_function_Table_maker_value(std::string filename,
                                    std::size_t num_scouts_per_ensemble) {
  return std::pair(filename, num_scouts_per_ensemble);
}

namespace deprecated{
using namespace var::deprecated;


auto get_function_Table_maker(std::string filename,
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

}
auto get_function_Table_maker_value_St(std::string filename) {
    return filename;
}
auto get_function_Table_maker_St(std::string filename) {
    using namespace macrodr;
    return [filename]() {
        return var::FuncMap_St(
            filename,
            Time_it_st(F(cuevi::step_stretch_cuevi_mcmc{},
                      cuevi::step_stretch_cuevi_mcmc{})),
            Time_it_st(
                F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{})),
            // Time_it_st(F(::deprecated::thermo_cuevi_randomized_jump_mcmc{},
            //           ::deprecated::thermo_cuevi_randomized_jump_mcmc{})),
            var::Time_it_st(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                           cuevi::step_stretch_cuevi_mcmc_per_walker{})),
            var::Time_it_st(F(logLikelihood_f{},
                           [](auto &&...x) {
                               return logLikelihood(std::forward<decltype(x)>(x)...);
                           })),
            var::Time_it_st(
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
                  })),
            var::Time_it_st(
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
                  })),
            var::Time_it_st(
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
                  })),
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
            var::Single_Thread_Memoizer(
                var::F(Calc_Qdt_step{},
                       [](auto &&...x) {
                           auto m = Macro_DMR{};
                           return m.calc_Qdt_ATP_step(
                               std::forward<decltype(x)>(x)...);
                       }),
                var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{}),
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
            var::Single_Thread_Memoizer(
                F(Calc_eigen{},
                  [](auto &&...x) {
                      auto m = Macro_DMR{};
                      return m.calc_eigen(std::forward<decltype(x)>(x)...);
                  }),
                var::Memoiza_all_values<Maybe_error<Qx_eig>, ATP_concentration>{})
            // var::Time_it(
            //     F(Calc_eigen{},
            //       [](auto &&...x) {
            //           auto m = Macro_DMR{};
            //           return m.calc_eigen(std::forward<decltype(x)>(x)...);
            //       }))
            
            );
    };
}




auto get_Experiment(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
    double frequency_of_sampling = 50e3, double initial_ATP = 0) {
  using namespace macrodr;
  auto [recording_conditions, recording] = macrodr::load_recording(filename);

  return Experiment(std::move(recording_conditions),
                    Frequency_of_Sampling(frequency_of_sampling),
                    initial_ATP_concentration(ATP_concentration(initial_ATP)));
}

auto get_Observations(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt") {
  auto [recording_conditions, recording] = macrodr::load_recording(filename);
  std::string out = filename.substr(0, filename.size() - 4) + "_recording.txt";
  save_Recording(out, ",", recording);
  return out;
}

auto get_Prior(double prior_error, const std::string modelname) {

  return std::pair(prior_error, modelname);
}

auto get_Likelihood(const std::string &model, bool adaptive_aproximation,
                    bool recursive_approximation, int averaging_approximation,
                    bool variance_correction_approximation,
                    bool variance_correction, std::size_t n_sub_dt) {

  return std::tuple(model, adaptive_aproximation, recursive_approximation,
                    averaging_approximation, variance_correction_approximation,
                    variance_correction, n_sub_dt);
}

auto get_CueviAlgorithm(
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

    std::size_t t_min_number_of_samples = 20, std::string filename = "haha",
    std::size_t thermo_jumps_every = 10) {
  using namespace macrodr;

  std::vector<std::size_t> t_segments_used = {73, 33, 22, 22};

  auto saving_itervals = Saving_intervals(
      Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),
                   Save_Likelihood_every(num_scouts_per_ensemble),
                   Save_Parameter_every(num_scouts_per_ensemble),
                   Save_Predictions_every(num_scouts_per_ensemble * 20)));

  return std::tuple(path, filename, t_segments_used, t_min_number_of_samples,
                    num_scouts_per_ensemble, number_trials_until_give_up,
                    min_fraction, thermo_jumps_every, max_iter_equilibrium,
                    n_points_per_decade, n_points_per_decade_fraction,
                    medium_beta, stops_at, includes_zero, saving_itervals,
                    random_jumps);
}

template <class T> struct return_type;

template <class R, class... A> struct return_type<R (*)(A...)> {
  using type = R;
};



using tablefun_value_type = typename return_type<
    std::decay_t<decltype(&get_function_Table_maker_value)>>::type;

using algo_type =
    typename return_type<std::decay_t<decltype(&get_CueviAlgorithm)>>::type;
using experiment_type =
    typename return_type<std::decay_t<decltype(&get_Experiment)>>::type;
using recording_type =
    typename return_type<std::decay_t<decltype(&get_Observations)>>::type;
using recording_value_type =
    typename return_type<std::decay_t<decltype(&load_Recording_value)>>::type;
using prior_value_type =
    typename return_type<std::decay_t<decltype(&load_Prior_value)>>::type;
using parameters_value_type =
    typename return_type<std::decay_t<decltype(&load_Parameter_value)>>::type;
using likelihood_type =
    typename return_type<std::decay_t<decltype(&get_Likelihood)>>::type;

// evidence(prior= prior_model, likelihoodModel= likelihood, data = simulation,
// experiment = experiment, algorithm= algorithm, function_table
// =function_table, init_seed =0)

void calc_likelihood(std::string outfilename, parameters_value_type par,
                     likelihood_type likelihood, recording_value_type recording,
                     experiment_type experiment, algo_type algorithm,
                     tablefun_value_type ft, std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto ftbl3 = get_function_Table_maker_St(filename)();

  auto Maybe_model_v = get_model(std::get<0>(likelihood));
  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [outfilename, &par, &ftbl3, &experiment, &recording, &likelihood,
         &algorithm, &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);

          auto [path, filename, t_segments_used, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                min_fraction, thermo_jump_factor, max_iter_equilibrium,
                n_points_per_decade, n_points_per_decade_fraction, medium_beta,
                stops_at, includes_zero, saving_itervals, random_jumps] =
              std::move(algorithm);

          auto [model, adaptive_aproximation, recursive_approximation,
                averaging_approximation, variance_correction_approximation,
                variance_correction, n_sub_dt] = likelihood;

          using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

          auto Maybe_param1 = var::load_Parameters<MyModel>(
              par.first, par.second, model0.model_name(), model0.names());
          Simulated_Recording<includes_N_state_evolution(true)> y;
          auto Maybe_y =
              load_Simulated_Recording(recording.first, recording.second, y);
          if (!Maybe_param1.valid() || !Maybe_y.valid()) {
            std::cerr << Maybe_param1.error()() << Maybe_y.error()();
          } else {

            auto param1 = std::move(Maybe_param1.value());

            std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
            auto number_of_traces = 7;
            auto number_of_segments = t_segments.size();
            t_segments.reserve(number_of_traces * t_segments.size());

            for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
              std::copy_n(t_segments.begin(), number_of_segments,
                          std::back_inserter(t_segments));

            std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};
            auto modelLikelihood =
                make_Likelihood_Model<uses_adaptive_aproximation(true),
                                      uses_recursive_aproximation(true),
                                      uses_averaging_aproximation(2),
                                      uses_variance_aproximation(false),
                                      uses_variance_correction_aproximation(
                                          false)>(
                    model0, Simulation_n_sub_dt(n_sub_dt));
            auto lik =
                Macro_DMR{}
                    .log_Likelihood<uses_adaptive_aproximation(false),
                                    uses_recursive_aproximation(true),
                                    uses_averaging_aproximation(2),
                                    uses_variance_aproximation(false),
                                    uses_variance_correction_aproximation(
                                        false),
                                    return_predictions(true)>(
                        ftbl3, model0, param1,
                        experiment, get<Recording>(y()));
            if (lik)
              save_Likelihood_Predictions(outfilename, lik.value(), y,
                                          experiment);
            else
              std::cerr << lik.error()();
          }
        },
        model_v);
  }
}

inline Maybe_error<std::vector<std::size_t>> load_segments_length_for_fractioning(const std::string& filename, std::string sep)
{
    std::ifstream f(filename);
    if (!f)
        return error_message(filename+ " cannot be opened");
    std::string line;
    std::getline(f,line);
    if (!f)
        return error_message(filename+ " has no data");
    
    std::stringstream ss(line);
    std::vector<std::size_t> out;
    std::size_t number_of_samples;
    
    while (ss>>number_of_samples)
    {
        out.push_back(number_of_samples);
        ss>>septr(sep);
    }
    
    return out;
    
}


void calc_evidence(prior_value_type prior, likelihood_type likelihood,
                   std::string recording, experiment_type experiment,
                   std::string segments,
                   algo_type algorithm, tablefun_value_type ft,
                   std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto ftbl3 = get_function_Table_maker_St(filename)();

  auto Maybe_model_v = get_model(std::get<0>(likelihood));
  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&ftbl3, &experiment, &recording,&segments, &prior, &likelihood, &algorithm,
         &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);

          auto [path, file_name, t_segments_used, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                min_fraction, thermo_jump_factor, max_iter_equilibrium,
                n_points_per_decade, n_points_per_decade_fraction, medium_beta,
                stops_at, includes_zero, saving_itervals, random_jumps] =
              std::move(algorithm);

          auto [model, adaptive_aproximation, recursive_approximation,
                averaging_approximation, variance_correction_approximation,
                variance_correction, n_sub_dt] = likelihood;

          using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

          std::string ModelName = model0.model_name();

          auto Maybe_param1_prior = var::load_Prior<MyModel>(
              prior.first, prior.second, model0.model_name(), model0.names());
          if (!Maybe_param1_prior) {
            std::cerr << Maybe_param1_prior.error()();
          } else {
            auto param1_prior = std::move(Maybe_param1_prior.value());
            std::size_t thermo_jumps_every =
                param1_prior.size() * thermo_jump_factor;

            Recording y;
            auto Maybe_y = load_Recording_Data(recording, ",", y);
            if (Maybe_y) {
              std::vector<std::size_t> t_segments = {73, 33, 22, 22,
                                                     1,  1,  1,  1};
              auto number_of_traces = 7;
              auto number_of_segments = t_segments.size();
              t_segments.reserve(number_of_traces * t_segments.size());

              for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
                std::copy_n(t_segments.begin(), number_of_segments,
                            std::back_inserter(t_segments));

              std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};

              std::size_t t_min_number_of_samples = 20;

              /**
               * @brief cbc cumulative evidence algorithm, ends using
               * convergence criteria
               */

              std::string filename = file_name + "_" + ModelName + "_" +
                                     time_now() + "_" + std::to_string(myseed);
              
              auto Maybe_t_segments_used = load_segments_length_for_fractioning(segments,",");
              
              if (Maybe_t_segments_used)
              {
                  auto t_segments_used=std::move(Maybe_t_segments_used.value());

              auto saving_itervals = Saving_intervals(Vector_Space(
                  Save_Evidence_every(num_scouts_per_ensemble),
                  Save_Likelihood_every(num_scouts_per_ensemble),
                  Save_Parameter_every(num_scouts_per_ensemble),
                  Save_Predictions_every(num_scouts_per_ensemble * 20)));

              auto cbc = new_cuevi_Model_by_iteration<MyModel>(
                  path, filename, t_segments_used, t_min_number_of_samples,
                  num_scouts_per_ensemble, number_trials_until_give_up,
                  min_fraction, thermo_jumps_every, max_iter_equilibrium,
                  n_points_per_decade, n_points_per_decade_fraction,
                  medium_beta, stops_at, includes_zero, saving_itervals,
                  random_jumps);

              auto modelLikelihood =
                  make_Likelihood_Model<uses_adaptive_aproximation(true),
                                        uses_recursive_aproximation(true),
                                        uses_averaging_aproximation(2),
                                        uses_variance_aproximation(false),
                                        uses_variance_correction_aproximation(
                                            false)>(
                      model0, Simulation_n_sub_dt(n_sub_dt));
              auto opt3 = cuevi::evidence(ftbl3, std::move(cbc), param1_prior,
                                          modelLikelihood, y, experiment,
                                          cuevi::Init_seed(myseed));
            }
            }
          }
        },
        model_v);
  }
}

void calc_evidence_continuation(prior_value_type prior,
                                likelihood_type likelihood,
                                std::string recording,
                                experiment_type experiment, algo_type algorithm,
                                tablefun_value_type ft, std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto ftbl3 = get_function_Table_maker_St(filename)();

  auto Maybe_model_v = get_model(std::get<0>(likelihood));
  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&ftbl3, &experiment, &recording, &prior, &likelihood, &algorithm,
         &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);

          auto [path, file_name, t_segments_used, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                min_fraction, thermo_jump_factor, max_iter_equilibrium,
                n_points_per_decade, n_points_per_decade_fraction, medium_beta,
                stops_at, includes_zero, saving_itervals, random_jumps] =
              std::move(algorithm);

          auto [model, adaptive_aproximation, recursive_approximation,
                averaging_approximation, variance_correction_approximation,
                variance_correction, n_sub_dt] = likelihood;

          using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

          std::string ModelName = model0.model_name();

          auto Maybe_param1_prior = var::load_Prior<MyModel>(
              prior.first, prior.second, model0.model_name(), model0.names());
          if (!Maybe_param1_prior) {
            std::cerr << Maybe_param1_prior.error()();
          } else {
            auto param1_prior = std::move(Maybe_param1_prior.value());
            std::size_t thermo_jumps_every =
                param1_prior.size() * thermo_jump_factor;

            Recording y;
            auto Maybe_y = load_Recording_Data(recording, ",", y);
            if (Maybe_y) {
              std::vector<std::size_t> t_segments = {73, 33, 22, 22,
                                                     1,  1,  1,  1};
              auto number_of_traces = 7;
              auto number_of_segments = t_segments.size();
              t_segments.reserve(number_of_traces * t_segments.size());

              for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
                std::copy_n(t_segments.begin(), number_of_segments,
                            std::back_inserter(t_segments));

              std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};

              std::size_t t_min_number_of_samples = 20;

              /**
               * @brief cbc cumulative evidence algorithm, ends using
               * convergence criteria
               */

              std::string filename = file_name + "_" + ModelName + "_" +
                                     time_now() + "_" + std::to_string(myseed);

              auto &t_segments_used = t_segments_7;

              auto saving_itervals = Saving_intervals(Vector_Space(
                  Save_Evidence_every(num_scouts_per_ensemble),
                  Save_Likelihood_every(num_scouts_per_ensemble),
                  Save_Parameter_every(num_scouts_per_ensemble),
                  Save_Predictions_every(num_scouts_per_ensemble * 20)));

              auto cbc = new_cuevi_Model_by_iteration<MyModel>(
                  path, filename, t_segments_used, t_min_number_of_samples,
                  num_scouts_per_ensemble, number_trials_until_give_up,
                  min_fraction, thermo_jumps_every, max_iter_equilibrium,
                  n_points_per_decade, n_points_per_decade_fraction,
                  medium_beta, stops_at, includes_zero, saving_itervals,
                  random_jumps);

              auto modelLikelihood =
                  make_Likelihood_Model<uses_adaptive_aproximation(true),
                                        uses_recursive_aproximation(true),
                                        uses_averaging_aproximation(2),
                                        uses_variance_aproximation(false),
                                        uses_variance_correction_aproximation(
                                            false)>(
                      model0, Simulation_n_sub_dt(n_sub_dt));
              auto opt3 = cuevi::evidence(ftbl3, std::move(cbc), param1_prior,
                                          modelLikelihood, y, experiment,
                                          cuevi::Init_seed(myseed));
            }
          }
        },
        model_v);
  }
}

std::string run_simulation(std::string filename, recording_type recording_file,
                           experiment_type experiment, std::size_t myseed,
                           std::string modelName,
                           parameters_value_type parameter_files, bool includeN,
                           std::size_t n_sub_dt) {
  using namespace macrodr;

  auto Maybe_recording = load_Recording(recording_file, std::string(","));
  if (!Maybe_recording)
    return Maybe_recording.error()();
  else {
    auto recording = std::move(Maybe_recording.value());
    auto Maybe_model_v = get_model(modelName);
    if (!Maybe_model_v)
      return Maybe_model_v.error()();
    else {
      auto model_v = std::move(Maybe_model_v.value());

      return std::visit(
          [&filename, &experiment, &recording, &myseed, &parameter_files,
           n_sub_dt, includeN](auto model0ptr) -> std::string {
            auto &model0 = *model0ptr;
            myseed = calc_seed(myseed);
            mt_64i mt(myseed);
            using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
            auto Maybe_parameter_values = var::load_Parameters<MyModel>(
                parameter_files.first, parameter_files.second,
                model0.model_name(), model0.names());
            if (!Maybe_parameter_values)
              return Maybe_parameter_values.error()();
            else {
              auto param1 = std::move(Maybe_parameter_values.value());
              save_Parameter<Parameters<MyModel>> s(filename, 1);
              if (includeN) {
                auto sim = Macro_DMR{}.sample(
                    mt, model0, param1, experiment,
                    Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)),
                    recording);
                if (!sim)
                  return "";
                else {
                  report_model(s, sim.value());
                  return filename + "_simulation.csv";
                }
              } else {
                auto sim = Macro_DMR{}.sample_N(
                    mt, model0, param1, experiment,
                    Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)),
                    recording);
                if (!sim)
                  return sim.error()();
                else {
                  report_model(s, sim.value());
                  return filename + "_simulation.csv";
                }
              }
            }
          },
          model_v);
    }
  }
}

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
auto calc_evidence_old(
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

    std::size_t t_min_number_of_samples = 20;

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
                     Save_Predictions_every(num_scouts_per_ensemble * 20)));

    auto cbc = new_cuevi_Model_by_iteration<MyModel>(
        path, filename, t_segments_used, t_min_number_of_samples,
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
} // namespace macrodr

#endif // CLI_MACRODR_COMMANDS_H
