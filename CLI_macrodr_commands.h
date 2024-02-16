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
#include "type_algebra.h"
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

template <class T> struct return_type;

template <class R, class... A> struct return_type<R (*)(A...)> {
  using type = R;
};

inline auto get_random_id(std::string prefix) {

  return prefix + "_" + std::to_string(calc_seed(0ul));
}

inline void write_text(std::string filename, std::string s) {
  std::ofstream f(filename);
  f << s;
}

static constexpr inline auto temp_script_file = "temp_script_file.txt";
inline void write_script(std::string program_name) {
  rename(temp_script_file, program_name.c_str());
}

inline auto load_Parameter_value(const std::string filename,
                                 std::string separator) {
  return std::pair(filename, separator);
}

inline auto load_Prior_value(std::string filename,
                             std::string separator = ",") {
  return std::pair(filename, separator);
}

inline auto load_Recording_value(std::string filename,
                                 std::string separator = ",") {
  return std::pair(filename, separator);
}


inline auto load_fractioned_Recording_value(std::string filename,
                                 std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_fractioned_Experiment(std::string filename,
                                            std::string separator = ",") {
    return std::pair(filename, separator);
}



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
inline auto get_function_Table_St(std::string filename,
                                  std::size_t save_every) {
  using namespace macrodr;
  return var::FuncMap_St(
      filename, save_every,
      Time_it_st(F(cuevi::step_stretch_cuevi_mcmc{},
                   cuevi::step_stretch_cuevi_mcmc{})),
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

inline auto
get_function_Table_maker_value(std::string filename,
                               std::size_t num_scouts_per_ensemble) {
  return std::pair(filename, num_scouts_per_ensemble);
}

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
inline auto get_function_Table_maker_value_St(std::string filename) {
  return filename;
}
inline auto get_function_Table_maker_St(std::string filename,
                                        std::size_t save_every) {
  using namespace macrodr;
  return [filename, save_every]() {
    return var::FuncMap_St(
               filename, save_every,
               Time_it_st(F(cuevi::step_stretch_cuevi_mcmc{},
                            cuevi::step_stretch_cuevi_mcmc{})),
               Time_it_st(F(cuevi::thermo_cuevi_jump_mcmc{},
                            cuevi::thermo_cuevi_jump_mcmc{})),
               // Time_it_st(F(::deprecated::thermo_cuevi_randomized_jump_mcmc{},
               //           ::deprecated::thermo_cuevi_randomized_jump_mcmc{})),
               var::Time_it_st(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                                 cuevi::step_stretch_cuevi_mcmc_per_walker{})),
               var::Time_it_st(F(logLikelihood_f{},
                                 [](auto &&...x) {
                                   return logLikelihood(
                                       std::forward<decltype(x)>(x)...);
                                 })),
               var::Time_it_st(F(
                   MacroR<uses_recursive_aproximation(true),
                          uses_averaging_aproximation(2),
                          uses_variance_aproximation(true)>{},
                   [](auto &&...x) {
                     auto m = Macro_DMR{};
                     return m
                         .Macror<uses_recursive_aproximation(true),
                                 uses_averaging_aproximation(2),
                                 uses_variance_aproximation(true),
                                 uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                   })),
               var::Time_it_st(F(
                   MacroR<uses_recursive_aproximation(true),
                          uses_averaging_aproximation(2),
                          uses_variance_aproximation(false)>{},
                   [](auto &&...x) {
                     auto m = Macro_DMR{};
                     return m
                         .Macror<uses_recursive_aproximation(true),
                                 uses_averaging_aproximation(2),
                                 uses_variance_aproximation(false),
                                 uses_variance_correction_aproximation(false)>(
                             std::forward<decltype(x)>(x)...);
                   })),
               var::Time_it_st(F(
                   MacroR<uses_recursive_aproximation(false),
                          uses_averaging_aproximation(2),
                          uses_variance_aproximation(false)>{},
                   [](auto &&...x) {
                     auto m = Macro_DMR{};
                     return m
                         .Macror<uses_recursive_aproximation(false),
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
               //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step,
               //     double>{}, num_scouts_per_ensemble / 2),
               var::Single_Thread_Memoizer(
                   var::F(Calc_Qdt_step{},
                          [](auto &&...x) {
                            auto m = Macro_DMR{};
                            return m.calc_Qdt_ATP_step(
                                std::forward<decltype(x)>(x)...);
                          }),
                   var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step,
                                           double>{}),
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
                   var::Memoiza_all_values<Maybe_error<Qx_eig>,
                                           ATP_concentration>{})
               // var::Time_it(
               //     F(Calc_eigen{},
               //       [](auto &&...x) {
               //           auto m = Macro_DMR{};
               //           return
               //           m.calc_eigen(std::forward<decltype(x)>(x)...);
               //       }))

               )
        .append_Fs<MacroR2>(
            in_progress::P<
                in_progress::S<::V<uses_recursive_aproximation(false)>,
                               ::V<uses_recursive_aproximation(true)>>,
                in_progress::S<::V<uses_averaging_aproximation(0)>,
                               ::V<uses_averaging_aproximation(1)>,
                               ::V<uses_averaging_aproximation(2)>>,
                in_progress::S<::V<uses_variance_aproximation(false)>,
                               ::V<uses_variance_aproximation(true)>>,
                in_progress::S<
                    ::V<uses_variance_correction_aproximation(false)>,
                    ::V<uses_variance_correction_aproximation(true)>>>{});
  };
}

inline auto get_Experiment(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
    double frequency_of_sampling = 50e3, double initial_ATP = 0) {
  using namespace macrodr;
  auto [recording_conditions, recording] = macrodr::load_recording(filename);

  return Experiment(std::move(recording_conditions),
                    Frequency_of_Sampling(frequency_of_sampling),
                    initial_ATP_concentration(ATP_concentration(initial_ATP)));
}

inline auto get_Observations(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt") {
  auto [recording_conditions, recording] = macrodr::load_recording(filename);
  std::string out = filename.substr(0, filename.size() - 4) + "_recording.txt";
  save_Recording(out, ",", recording);
  return out;
}

/*
make_Likelihood_Model<uses_adaptive_aproximation(true),
                      uses_recursive_aproximation(true),
                      uses_averaging_aproximation(2),
                      uses_variance_aproximation(false),
                      uses_variance_correction_aproximation(
                          false)>

*/

inline auto get_Likelihood_Algorithm(
    bool use_adaptive_approx, bool use_recursive_approx,
    std::size_t n_point_averaging, bool variance_approximation,
    bool uses_variance_approximation_correction, std::string model, double) {}

inline auto get_Prior(double prior_error, const std::string modelname) {

  return std::pair(prior_error, modelname);
}

inline auto set_Likelihood_algorithm(bool adaptive_aproximation,
                                     bool recursive_approximation,
                                     int averaging_approximation,
                                     bool variance_correction_approximation,
                                     bool variance_correction,
                                     std::size_t n_sub_dt) {

  return std::tuple(adaptive_aproximation, recursive_approximation,
                    averaging_approximation, variance_correction_approximation,
                    variance_correction, n_sub_dt);
}

inline auto set_Fraction_algorithm(double min_fraction,
                                   double n_points_per_decade_fraction,
                                   std::string segments) {
  return std::tuple(min_fraction, n_points_per_decade_fraction, segments);
}
using fraction_algo_type =
    typename return_type<std::decay_t<decltype(&set_Fraction_algorithm)>>::type;

inline auto set_CueviAlgorithm(
    std::size_t num_scouts_per_ensemble = 16,
    std::size_t number_trials_until_give_up = 1e5,
    double stops_at = 1e-15,
    double medium_beta = 1e-2,
    bool includes_zero = true,
    bool random_jumps = true,
    std::size_t max_iter_equilibrium = 50000,
    std::string path = "",
    double n_points_per_decade = 1,
    std::size_t t_min_number_of_samples = 20,
    std::string filename = "haha",
    std::size_t thermo_jumps_every = 10) {
  using namespace macrodr;

  auto saving_itervals = Saving_intervals(
      Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),
                   Save_Likelihood_every(num_scouts_per_ensemble),
                   Save_Parameter_every(num_scouts_per_ensemble),
                   Save_Predictions_every(num_scouts_per_ensemble * 20)));

  return std::tuple(path, filename, t_min_number_of_samples,
                    num_scouts_per_ensemble, number_trials_until_give_up,
                    thermo_jumps_every, max_iter_equilibrium,
                    n_points_per_decade, medium_beta, stops_at, includes_zero,
                    saving_itervals, random_jumps);
}

using tablefun_value_type = typename return_type<
    std::decay_t<decltype(&get_function_Table_maker_value)>>::type;

using cuevi_algo_type =
    typename return_type<std::decay_t<decltype(&set_CueviAlgorithm)>>::type;
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
calc_experiment_fractions(std::string save_name,std::string recording, experiment_type experiment,                                  fraction_algo_type fraction_algo, std::string model,std::size_t i_seed )
{
    auto myseed = calc_seed(i_seed);
    
    auto init_seed=calc_seed(i_seed);
    mt_64i mt(init_seed);
    
    auto [min_fraction, n_points_per_decade_fraction, segments] =
        std::move(fraction_algo);
    
    auto Maybe_segments=load_segments_length_for_fractioning(segments,",");
    if (!Maybe_segments)
        return Maybe_segments.error();
    Recording y;
    auto Maybe_y = load_Recording_Data(recording, ",", y);
    
    macrodr::experiment_fractioner frac(Maybe_segments.value(),0);
    
    auto maybe_model=get_model(model);
    if (!maybe_model)
        return maybe_model.error();
    
    auto param_size=  std::visit([&](auto model0ptr) {
        return model0ptr->parameters().size();},maybe_model.value());
    
    auto [ys, xs] = frac(
        y, experiment, mt, param_size * min_fraction, n_points_per_decade_fraction);
    auto filename=save_name+"_"+std::to_string(myseed);
    save_fractioned_experiment(filename+"_experiment.csv",",",xs);
    save_fractioned_Recording(filename+"_recording.csv",",",ys);
    return std::tuple(filename+"_experiment.csv",filename+"_recording.csv", get<Frequency_of_Sampling>(experiment)(),get<initial_ATP_concentration>(experiment)()());
    
    
}


inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_simulation_fractions(std::string save_name,std::string simulation, experiment_type experiment,                                  fraction_algo_type fraction_algo, std::string model,std::size_t i_seed )
{
    auto myseed = calc_seed(i_seed);
    
    auto init_seed=calc_seed(i_seed);
    mt_64i mt(init_seed);
    
    auto [min_fraction, n_points_per_decade_fraction, segments] =
        std::move(fraction_algo);
    
    auto Maybe_segments=load_segments_length_for_fractioning(segments,",");
    if (!Maybe_segments)
        return Maybe_segments.error();
    Simulated_Recording<includes_N_state_evolution(true)> y;
    auto Maybe_y = load_Simulated_Recording(simulation, ",", y);
    if (!Maybe_y)
        return Maybe_y.error();
    
    macrodr::experiment_fractioner frac(Maybe_segments.value(),0);
    
    auto maybe_model=get_model(model);
    if (!maybe_model)
        return maybe_model.error();
    
    auto param_size=  std::visit([&](auto model0ptr) {
        return model0ptr->parameters().size();},maybe_model.value());
    
    auto [ys, xs] = frac(
        y, experiment, mt, param_size * min_fraction, n_points_per_decade_fraction);
    auto filename=save_name+"_"+std::to_string(myseed);
    save_fractioned_experiment(filename+"_frac_experiment.csv",",",xs);
    save_fractioned_Recording(filename+"_frac_recording.csv",",",ys);
    return std::tuple(filename+"_frac_experiment.csv",filename+"_frac_recording.csv", get<Frequency_of_Sampling>(experiment)(),get<initial_ATP_concentration>(experiment)()());
    
    
}


using fractioned_experiment_type =
    typename return_type<std::decay_t<decltype(&calc_experiment_fractions)>>::type;


using fractioned_simulation_type =
    typename return_type<std::decay_t<decltype(&calc_simulation_fractions)>>::type;

using likelihood_algo_type = typename return_type<
    std::decay_t<decltype(&set_Likelihood_algorithm)>>::type;

// evidence(prior= prior_model, likelihoodModel= likelihood, data = simulation,
// experiment = experiment, algorithm= algorithm, function_table
// =function_table, init_seed =0)






inline void calc_likelihood(std::string outfilename,
                            std::string model,
                            parameters_value_type par,
                            likelihood_algo_type likelihood,
                            recording_value_type recording,
                            experiment_type experiment, cuevi_algo_type algorithm,
                            tablefun_value_type ft) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto save_every = num_scouts_per_ensemble;
  auto ftbl3 = get_function_Table_maker_St(filename, save_every)();

  auto Maybe_model_v = get_model(model);
  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [outfilename, &par, &ftbl3, &experiment, &recording, &likelihood,
         &algorithm](auto model0ptr) {
          auto &model0 = *model0ptr;

          auto [path, filename, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                thermo_jump_factor, max_iter_equilibrium, n_points_per_decade,
                medium_beta, stops_at, includes_zero, saving_itervals,
                random_jumps] = std::move(algorithm);

          auto [adaptive_aproximation, recursive_approximation,
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
                        ftbl3, model0, param1, experiment, get<Recording>(y()));
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


inline std::string calc_fraction_likelihood( std::string file_name,
                                            std::string model,
                                     parameters_value_type par,
                                     fractioned_simulation_type  Maybe_frac_simulation,
                                     likelihood_algo_type likelihood_algo,
                                     tablefun_value_type ft) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto save_every=num_scouts_per_ensemble;
  
  auto ftbl3 = get_function_Table_maker_St(filename, save_every)();
  
  auto Maybe_model_v = get_model(model);
  if (!Maybe_model_v)
    return Maybe_model_v.error()();
  auto model_v = std::move(Maybe_model_v.value());
  if (!Maybe_frac_simulation)
      return Maybe_frac_simulation.error()();
  auto frac_simulation=std::move(Maybe_frac_simulation.value());
  return std::visit(
      [&ftbl3, &frac_simulation, &likelihood_algo, &par,
       &file_name](auto model0ptr) {
        auto &model0 = *model0ptr;
          auto [experiment,simulation, fs,iniATP]=frac_simulation;

        auto [adaptive_aproximation, recursive_approximation,
              averaging_approximation, variance_correction,
              variance_correction_approximation, n_sub_dt] = likelihood_algo;

        using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

        std::string ModelName = model0.model_name();
        auto Maybe_param1 = var::load_Parameters<MyModel>(
            par.first, par.second, model0.model_name(), model0.names());
        
        std::vector<Experiment> xs;
        std::vector<Simulated_Recording<includes_N_state_evolution(true)>> ys;
        auto Maybe_e = load_fractioned_experiment(experiment, ",",
                                                  fs,iniATP,xs);
        
        auto Maybe_y = load_fractioned_simulation(simulation, ",",
                                                  ys);
        
        if (xs.size()!=ys.size())
            return std::string()+"number of fractions mismatch between recordings and experiments";
        
        if (!(Maybe_e.valid()&&Maybe_y.valid() && Maybe_param1.valid()))
            return Maybe_e.error()()+Maybe_y.error()() + Maybe_param1.error()();
        
        
        auto param1 = std::move(Maybe_param1.value());

        std::string filename = file_name + "_" + ModelName + "_" + time_now();

        auto modelLikelihood_v = Likelihood_Model_v{}.bool_op(
            uses_adaptive_aproximation(adaptive_aproximation),
            uses_recursive_aproximation(recursive_approximation),
            uses_averaging_aproximation(averaging_approximation),
            uses_variance_aproximation(variance_correction),
            uses_variance_correction_aproximation(
                variance_correction_approximation),
            model0, Simulation_n_sub_dt(n_sub_dt));

        return std::visit(
            [&ftbl3, &param1, &ys, &xs,
             &filename](auto &modelLikelihood) {
              auto Maybe_lik = fractioned_logLikelihoodPredictions(ftbl3, modelLikelihood,
                                                        param1, ys, xs);
              if (!Maybe_lik)
                return Maybe_lik.error()();
              else
                  {
                 save_fractioned_Likelihood_Predictions(filename+"frac_likelihood.csv", Maybe_lik.value(), ys,
                                            xs);
                  return filename;
              }
            },
            modelLikelihood_v);
      },
      model_v);
}

inline void calc_fraction_evidence(std::string model,
                          prior_value_type prior,
                          likelihood_algo_type likelihood,
                          fractioned_simulation_type Maybe_frac_experiment,
                          fraction_algo_type fraction_algo,
                          cuevi_algo_type cuevi_algorithm,
                          tablefun_value_type ft, std::size_t myseed) {
    using namespace macrodr;
    if (! Maybe_frac_experiment)
        return;
    auto frac_experiment=Maybe_frac_experiment.value();
    
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    
    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every)();
    
    auto Maybe_model_v = get_model(model);
    
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &frac_experiment,   &prior, &likelihood,
             &fraction_algo, &cuevi_algorithm, &myseed](auto model0ptr) {
                auto &model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                auto [fname_experiment,fname_simulation, fs,iniATP]=frac_experiment;
                std::vector<Experiment> xs;
                std::vector<Recording> ys;
                auto Maybe_e = load_fractioned_experiment(fname_experiment, ",",
                                                          fs,iniATP,xs);
                
                auto Maybe_ys = load_fractioned_Recording(fname_simulation, ",",
                                                          ys);
                
                auto [min_fraction, n_points_per_decade_fraction, segments] =
                    std::move(fraction_algo);
                auto [path, file_name, t_min_number_of_samples,
                      num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium, n_points_per_decade,
                      medium_beta, stops_at, includes_zero, saving_itervals,
                      random_jumps] = std::move(cuevi_algorithm);
                
                auto [adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction,
                      variance_correction_approximation, n_sub_dt] = likelihood;
                
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
                    
                    if (Maybe_ys.valid()&& Maybe_e.valid()) {
                        
                        std::string filename = file_name + "_" + ModelName + "_" +
                                               time_now() + "_" + std::to_string(myseed);
                        
                        auto Maybe_t_segments_used =
                            load_segments_length_for_fractioning(segments, ",");
                        
                        if (Maybe_t_segments_used) {
                            auto t_segments_used = std::move(Maybe_t_segments_used.value());
                            
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
                            
                            auto modelLikelihood_v = Likelihood_Model_v{}.bool_op(
                                uses_adaptive_aproximation(adaptive_aproximation),
                                uses_recursive_aproximation(recursive_approximation),
                                uses_averaging_aproximation(averaging_approximation),
                                uses_variance_aproximation(variance_correction),
                                uses_variance_correction_aproximation(
                                    variance_correction_approximation),
                                model0, Simulation_n_sub_dt(n_sub_dt));
                            // using m2=typename
                            // decltype(modelLikelihood_v)::paseModelLikelihoodv;
                            
                            std::visit(
                                [&ftbl3, &cbc, &param1_prior, &ys, &xs,
                                 &myseed](auto &modelLikelihood) {
                                    auto opt3 = cuevi::evidence_fraction(
                                        ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                                        ys, xs, cuevi::Init_seed(myseed));
                                },
                                modelLikelihood_v);
                        }
                    } else
                        std::cerr << Maybe_ys.error()();
                }
            },
            model_v);
    }
}



inline void calc_evidence(std::string model,
                          prior_value_type prior,
                          likelihood_algo_type likelihood,
                          std::string recording,
                          experiment_type experiment,
                          fraction_algo_type fraction_algo,
                          cuevi_algo_type cuevi_algorithm,
                          tablefun_value_type ft, std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  
  auto save_every = num_scouts_per_ensemble;
  auto ftbl3 = get_function_Table_maker_St(filename, save_every)();

  auto Maybe_model_v = get_model(model);

  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&ftbl3, &experiment, &recording,  &prior, &likelihood,
        &fraction_algo, &cuevi_algorithm, &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);
          
          auto [min_fraction, n_points_per_decade_fraction, segments] =
              std::move(fraction_algo);
          auto [path, file_name, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                thermo_jump_factor, max_iter_equilibrium, n_points_per_decade,
                medium_beta, stops_at, includes_zero, saving_itervals,
                random_jumps] = std::move(cuevi_algorithm);

          auto [adaptive_aproximation, recursive_approximation,
                averaging_approximation, variance_correction,
                variance_correction_approximation, n_sub_dt] = likelihood;

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

              std::string filename = file_name + "_" + ModelName + "_" +
                                     time_now() + "_" + std::to_string(myseed);

              auto Maybe_t_segments_used =
                  load_segments_length_for_fractioning(segments, ",");

              if (Maybe_t_segments_used) {
                auto t_segments_used = std::move(Maybe_t_segments_used.value());

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

                auto modelLikelihood_v = Likelihood_Model_v{}.bool_op(
                    uses_adaptive_aproximation(adaptive_aproximation),
                    uses_recursive_aproximation(recursive_approximation),
                    uses_averaging_aproximation(averaging_approximation),
                    uses_variance_aproximation(variance_correction),
                    uses_variance_correction_aproximation(
                        variance_correction_approximation),
                    model0, Simulation_n_sub_dt(n_sub_dt));
                // using m2=typename
                // decltype(modelLikelihood_v)::paseModelLikelihoodv;

                std::visit(
                    [&ftbl3, &cbc, &param1_prior, &y, &experiment,
                     &myseed](auto &modelLikelihood) {
                      auto opt3 = cuevi::evidence(
                          ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                          y, experiment, cuevi::Init_seed(myseed));
                    },
                    modelLikelihood_v);
              }
            } else
              std::cerr << Maybe_y.error()();
          }
        },
        model_v);
  }
}

inline void calc_evidence_continuation(std::string model,
    prior_value_type prior, likelihood_algo_type likelihood,
    std::string recording, experiment_type experiment,fraction_algo_type fraction_algo, cuevi_algo_type algorithm,
    tablefun_value_type ft, std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);
  auto save_every = num_scouts_per_ensemble;
  auto ftbl3 = get_function_Table_maker_St(filename, save_every)();

  auto Maybe_model_v = get_model(model);
  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&ftbl3, &experiment, &recording, &prior, &likelihood, &fraction_algo,&algorithm,
         &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);
          /*
           * path, filename, t_min_number_of_samples,
                    num_scouts_per_ensemble, number_trials_until_give_up,
                    thermo_jumps_every, max_iter_equilibrium,
                    n_points_per_decade, medium_beta, stops_at, includes_zero,
                    saving_itervals, random_jumps
           * */
          auto [min_fraction, n_points_per_decade_fraction, segments] =
              std::move(fraction_algo);
          
          
          auto [path, file_name, t_min_number_of_samples,
                num_scouts_per_ensemble, number_trials_until_give_up,
                thermo_jump_factor, max_iter_equilibrium, n_points_per_decade,
                 medium_beta, stops_at,
                includes_zero, saving_itervals, random_jumps] =
              std::move(algorithm);

          auto [adaptive_aproximation, recursive_approximation,
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

inline auto set_simulation_algorithm(bool includeN, std::size_t n_sub_dt) {
  return std::pair(includeN, n_sub_dt);
}

using simulation_algo_type = typename return_type<
    std::decay_t<decltype(&set_simulation_algorithm)>>::type;

inline std::string run_simulation(std::string filename_prefix,
                                  recording_type recording_file,
                                  experiment_type experiment,
                                  std::size_t myseed, std::string modelName,
                                  parameters_value_type parameter_files,
                                  simulation_algo_type sim_algo_type) {
  using namespace macrodr;

  auto [includeN, n_sub_dt] = sim_algo_type;
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
          [&filename_prefix, &experiment, &recording, &myseed, &parameter_files,
           n_sub_dt, includeN](auto model0ptr) -> std::string {
            auto &model0 = *model0ptr;
            myseed = calc_seed(myseed);
            mt_64i mt(myseed);
            using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
            auto Maybe_parameter_values = var::load_Parameters<MyModel>(
                parameter_files.first, parameter_files.second,
                model0.model_name(), model0.names());
            std::string ModelName = model0.model_name();

            std::string filename = filename_prefix + "_" + ModelName + "_" +
                                   time_now() + "_" + std::to_string(myseed);

            if (!Maybe_parameter_values)
              return Maybe_parameter_values.error()();
            else {
              auto param1 = std::move(Maybe_parameter_values.value());
              save_Parameter<Parameters<MyModel>> s(filename, 1);
              if (!includeN) {
                auto sim = Macro_DMR{}.sample(
                    mt, model0, param1, experiment,
                    Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)),
                    recording);
                if (!sim)
                  return "";
                else {
                    save_Recording(filename + "_simulation.csv",",", get<Recording>(sim.value()()));
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
                  save_Simulated_Recording(filename + "_N_simulation.csv",",", sim.value());
                  return filename + "_N_simulation.csv";
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
