#ifndef CLI_FUNCTION_TABLE_H
#define CLI_FUNCTION_TABLE_H

#include <cstddef>
#include <string>

#include "CLI_macro_dr_base.h"
//#include "cuevi.h"
#include "function_memoization.h"
#include "maybe_error.h"
#include "models_used.h"
#include "parallel_tempering.h"
#include "parameters.h"
//#include "parameters_derivative.h"
#include "qmodel.h"

namespace macrodr::cmd {

inline auto get_function_Table_maker_St(std::string filename, std::size_t sampling_interval,
                                        std::size_t max_number_of_values_per_iteration) {
    using namespace macrodr;
    return [filename, sampling_interval, max_number_of_values_per_iteration]() {
        return var::FuncMap_St(
            filename, sampling_interval, max_number_of_values_per_iteration,
            //           Time_it_st(F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{})),
            //           Time_it_st(F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{})),
            //           var::Time_it_st(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
            //                             cuevi::step_stretch_cuevi_mcmc_per_walker{})),
            var::Time_it_st(
                F(logLikelihood_f{},
                  [](auto&&... x) { return logLikelihood(std::forward<decltype(x)>(x)...); })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<true>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                                  uses_variance_aproximation<true>,
                                  uses_taylor_variance_correction_aproximation<false>>(
                              std::forward<decltype(x)>(x)...);
                  })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<false>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                                  uses_variance_aproximation<false>,
                                  uses_taylor_variance_correction_aproximation<false>>(
                              std::forward<decltype(x)>(x)...);
                  })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<false>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<false>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<false>,
                                  uses_averaging_aproximation<2>, uses_variance_aproximation<false>,
                                  uses_taylor_variance_correction_aproximation<false>>(
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
            //     var::Memoiza_all_values<Maybe_error<Qdt>, Agonist_step,
            //     double>{}, num_scouts_per_ensemble / 2),
            var::Single_Thread_Memoizer(
                var::F(Calc_Qdt_step{},
                       [](auto& f, auto& m, auto& t_step, double fs) {
                           if constexpr (false) {
                               auto test_der_t_Qdt = var::test_Derivative(
                                   [&t_step, &fs, &f](auto const& l_m) {
                                       auto ma = Macro_DMR{};
                                       return ma.calc_Qdt_agonist_step(f, l_m, t_step, fs);
                                   },
                                   1e-6, 1e-2, m);
                               if (true && !test_der_t_Qdt) {
                                   std::cerr << test_der_t_Qdt.error()();
                                   std::cerr << "\nt_step\n" << t_step;
                                   abort();
                               }
                           }
                           auto ma = Macro_DMR{};

                           return ma.calc_Qdt_agonist_step(std::forward<decltype(f)>(f), m, t_step, fs);
                       }),
                var::Memoiza_overload<
                    var::Memoiza_all_values<Maybe_error<Qdt>, Agonist_step, double>,
                    var::Memoiza_all_values<
                        Maybe_error<var::Derivative<Qdt, var::Parameters_transformed>>, Agonist_step,
                        double>>{}),
            // var::Parameters_transformed
            var::Single_Thread_Memoizer(
                var::F(Calc_Qdtm_step{},
                       [](auto& f, auto& m, auto& t_step, double fs) {
                           if constexpr (false) {
                               auto test_der_t_Qdtm = var::test_Derivative(
                                   [&t_step, &fs, &f](auto const& l_m) {
                                       auto ma = Macro_DMR{};
                                       return ma.calc_Qdtm_agonist_step(f, l_m, t_step, fs);
                                   },
                                   1e-6, 1e-2, m);
                               if (true && !test_der_t_Qdtm) {
                                   std::cerr << test_der_t_Qdtm.error()();
                                   std::cerr << "\nt_step\n" << t_step;
                                   std::cerr << "\nmodel\n" << m;

                                   abort();
                               }
                           }
                           auto ma = Macro_DMR{};

                           return ma.calc_Qdtm_agonist_step(std::forward<decltype(f)>(f), m, t_step,
                                                        fs);
                       }),
                var::Memoiza_overload<
                    var::Memoiza_all_values<Maybe_error<Qdtm>, Agonist_step, double>,
                    var::Memoiza_all_values<
                        Maybe_error<var::Derivative<Qdtm, var::Parameters_transformed>>, Agonist_step,
                        double>>{}),
            // var::Time_it(
            //     var::F(Calc_Qdt_step{},
            //            [](auto &&...x) {
            //              auto m = Macro_DMR{};
            //              return
            //              m.calc_Qdt_agonist_step(std::forward<decltype(x)>(x)...);
            //            })),

            var::F(Calc_Qdt{},
                   [](auto&&... x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                   }),
            F(Calc_Qx{},
              [](auto&&... x) {
                  auto m = Macro_DMR{};
                  return m.calc_Qx(std::forward<decltype(x)>(x)...);
              }),
            var::Single_Thread_Memoizer(
                F(Calc_eigen{},
                  [](auto& mo, auto& agonist, auto&&... x) {
                      auto m = Macro_DMR{};

                      auto out = m.calc_eigen(mo, agonist, std::forward<decltype(x)>(x)...);
                      if constexpr (false) {
                          if (out) {
                              auto test_der_eigen = var::test_Derivative(
                                  [&m, &agonist](auto l_m) { return m.calc_eigen(l_m, agonist); }, 1e-6,
                                  1e-2, mo);

                              if (!test_der_eigen) {
                                  std::cerr << test_der_eigen.error()();
                                  //                     return test_der_eigen.error();
                                  std::abort();
                              }
                          }
                      }
                      return out;
                  }),
                var::Memoiza_overload<
                    var::Memoiza_all_values<Maybe_error<Eigs>, Agonist_concentration>,
                    var::Memoiza_all_values<
                        Maybe_error<var::Derivative<Eigs, var::Parameters_transformed>>,
                        Agonist_concentration>>{})

            // var::Time_it(
            //     F(Calc_eigen{},
            //       [](auto &&...x) {
            //           auto m = Macro_DMR{};
            //           return
            //           m.calc_eigen(std::forward<decltype(x)>(x)...);
            //       }))

            )
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation<false>>,
            //                        ::V<uses_recursive_aproximation<true>>>,
            //         in_progress::S<::V<uses_averaging_aproximation<2>>>,
            //         in_progress::S<::V<uses_variance_aproximation<true>>>,
            //         in_progress::S<
            //             ::V<uses_taylor_variance_correction_aproximation<false>>>>{});
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation<false>>,
            //                        ::V<uses_recursive_aproximation<true>>>,
            //         in_progress::S<::V<uses_averaging_aproximation<0>>,
            //                        ::V<uses_averaging_aproximation<1>>,
            //                        ::V<uses_averaging_aproximation<2>>>,
            //         in_progress::S<::V<uses_variance_aproximation<false>>,
            //                        ::V<uses_variance_aproximation<true>>>,
            //         in_progress::S<
            //             ::V<uses_taylor_variance_correction_aproximation<false>>,
            //             ::V<uses_taylor_variance_correction_aproximation<true>>>>{})
            ;
    };
}

inline auto get_function_Table_maker_St_no_Qdt_memoization(
    std::string filename, std::size_t sampling_interval,
    std::size_t max_number_of_values_per_iteration) {
    using namespace macrodr;
    return [filename, sampling_interval, max_number_of_values_per_iteration]() {
        return var::FuncMap_St(
            filename, sampling_interval, max_number_of_values_per_iteration,
            //Time_it_st(F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{})),
            // Time_it_st(F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{})),
            //  var::Time_it_st(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
            //                    cuevi::step_stretch_cuevi_mcmc_per_walker{})),
            var::Time_it_st(
                F(logLikelihood_f{},
                  [](auto&&... x) { return logLikelihood(std::forward<decltype(x)>(x)...); })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<true>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                                  uses_variance_aproximation<true>,
                                  uses_taylor_variance_correction_aproximation<false>>(
                              std::forward<decltype(x)>(x)...);
                  })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<false>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<true>, uses_averaging_aproximation<2>,
                                  uses_variance_aproximation<false>,
                                  uses_taylor_variance_correction_aproximation<false>>(
                              std::forward<decltype(x)>(x)...);
                  })),
            var::Time_it_st(
                F(MacroR<uses_recursive_aproximation<false>, uses_averaging_aproximation<2>,
                         uses_variance_aproximation<false>>{},
                  [](auto&&... x) {
                      auto m = Macro_DMR{};
                      return m
                          .Macror<uses_recursive_aproximation<false>,
                                  uses_averaging_aproximation<2>, uses_variance_aproximation<false>,
                                  uses_taylor_variance_correction_aproximation<false>>(
                              std::forward<decltype(x)>(x)...);
                  })),
            var::F(Calc_Qdt_step{},
                   [](auto& f, auto& m, auto& t_step, double fs) {
                       if constexpr (false) {
                           auto test_der_t_Qdt = var::test_Derivative(
                               [&t_step, &fs, &f](auto const& l_m) {
                                   auto ma = Macro_DMR{};
                                   return ma.calc_Qdt_agonist_step(f, l_m, t_step, fs);
                               },
                               1e-6, 1e-2, m);
                           if (true && !test_der_t_Qdt) {
                               std::cerr << test_der_t_Qdt.error()();
                               std::cerr << "\nt_step\n" << t_step;
                               abort();
                           }
                       }
                       auto ma = Macro_DMR{};

                       return ma.calc_Qdt_agonist_step(std::forward<decltype(f)>(f), m, t_step, fs);
                   }),
            var::F(Calc_Qdtm_step{},
                   [](auto& f, auto& m, auto& t_step, double fs) {
                       if constexpr (false) {
                           auto test_der_t_Qdtm = var::test_Derivative(
                               [&t_step, &fs, &f](auto const& l_m) {
                                   auto ma = Macro_DMR{};
                                   return ma.calc_Qdtm_agonist_step(f, l_m, t_step, fs);
                               },
                               1e-6, 1e-2, m);
                           if (true && !test_der_t_Qdtm) {
                               std::cerr << test_der_t_Qdtm.error()();
                               std::cerr << "\nt_step\n" << t_step;
                               std::cerr << "\nmodel\n" << m;

                               abort();
                           }
                       }
                       auto ma = Macro_DMR{};

                       return ma.calc_Qdtm_agonist_step(std::forward<decltype(f)>(f), m, t_step, fs);
                   }),

            var::F(Calc_Qdt{},
                   [](auto&&... x) {
                       auto m = Macro_DMR{};
                       return m.calc_Qdt(std::forward<decltype(x)>(x)...);
                   }),
            F(Calc_Qx{},
              [](auto&&... x) {
                  auto m = Macro_DMR{};
                  return m.calc_Qx(std::forward<decltype(x)>(x)...);
              }),
            var::Single_Thread_Memoizer(
                F(Calc_eigen{},
                  [](auto& mo, auto& agonist, auto&&... x) {
                      auto m = Macro_DMR{};

                      auto out = m.calc_eigen(mo, agonist, std::forward<decltype(x)>(x)...);
                      if constexpr (false) {
                          if (out) {
                              auto test_der_eigen = var::test_Derivative(
                                  [&m, &agonist](auto l_m) { return m.calc_eigen(l_m, agonist); }, 1e-6,
                                  1e-2, mo);

                              if (!test_der_eigen) {
                                  std::cerr << test_der_eigen.error()();
                                  //                     return test_der_eigen.error();
                                  std::abort();
                              }
                          }
                      }
                      return out;
                  }),
                var::Memoiza_overload<
                    var::Memoiza_all_values<Maybe_error<Eigs>, Agonist_concentration>,
                    var::Memoiza_all_values<
                        Maybe_error<var::Derivative<Eigs, var::Parameters_transformed>>,
                        Agonist_concentration>>{})

            // var::Time_it(
            //     F(Calc_eigen{},
            //       [](auto &&...x) {
            //           auto m = Macro_DMR{};
            //           return
            //           m.calc_eigen(std::forward<decltype(x)>(x)...);
            //       }))

            )
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation<false>>,
            //                        ::V<uses_recursive_aproximation<true>>>,
            //         in_progress::S<::V<uses_averaging_aproximation<2>>>,
            //         in_progress::S<::V<uses_variance_aproximation<true>>>,
            //         in_progress::S<
            //             ::V<uses_taylor_variance_correction_aproximation<false>>>>{});
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation<false>>,
            //                        ::V<uses_recursive_aproximation<true>>>,
            //         in_progress::S<::V<uses_averaging_aproximation<0>>,
            //                        ::V<uses_averaging_aproximation<1>>,
            //                        ::V<uses_averaging_aproximation<2>>>,
            //         in_progress::S<::V<uses_variance_aproximation<false>>,
            //                        ::V<uses_variance_aproximation<true>>>,
            //         in_progress::S<
            //             ::V<uses_taylor_variance_correction_aproximation<false>>,
            //             ::V<uses_taylor_variance_correction_aproximation<true>>>>{})
            ;
    };
}

using tablefun_value_type =
    typename return_type<std::decay_t<decltype(&get_function_Table_maker_value)>>::type;

inline std::string run_simulation(std::string filename_prefix, recording_type recording_file,
                                  experiment_type experiment, std::size_t myseed,
                                  std::string modelName, parameters_value_type parameter_files,
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
                [&filename_prefix, &experiment, &recording, &myseed, &parameter_files, n_sub_dt,
                 includeN](auto model0ptr) -> std::string {
                    auto& model0 = *model0ptr;
                    myseed = calc_seed(myseed);
                    mt_64i mt(myseed);
                    using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                    auto Maybe_parameter_values =
                        var::load_Parameters(parameter_files.first, parameter_files.second,
                                             model0.model_name(), model0.names());
                    std::string ModelName = model0.model_name();

                    std::string filename = filename_prefix + "_" + ModelName + "_" + time_now() +
                                           "_" + std::to_string(myseed);

                    if (!Maybe_parameter_values)
                        return Maybe_parameter_values.error()();
                    else {
                        auto param1 = Maybe_parameter_values.value().standard_parameter();
                        save_Parameter<var::Parameters_transformed> s(filename, 1ul, 200ul);
                        if (!includeN) {
                            auto sim = Macro_DMR{}.sample(
                                mt, model0, param1, experiment,
                                Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)), recording);
                            if (!sim)
                                return "";
                            else {
                                save_Recording(filename + "_simulation.csv", ",",
                                               get<Recording>(sim.value()()));
                                return filename + "_simulation.csv";
                            }
                        } else {
                            auto sim = Macro_DMR{}.sample_N(
                                mt, model0, param1, experiment,
                                Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)), recording);
                            if (!sim)
                                return sim.error()();
                            else {
                                save_Simulated_Recording(filename + "_N_simulation.csv", ",",
                                                         sim.value());
                                return filename + "_N_simulation.csv";
                            }
                        }
                    }
                },
                model_v);
        }
    }
}

#ifdef ZOMBIE
namespace zombie {

inline auto set_CueviAlgorithm(
    std::size_t num_scouts_per_ensemble, std::size_t number_trials_until_give_up, double stops_at,
    double medium_beta, bool includes_zero, bool random_jumps, std::size_t max_iter_equilibrium,
    std::string path, double n_points_per_decade_beta_high, double n_points_per_decade_beta_low,
    bool average_the_agonist_evolution, std::string filename, std::size_t thermo_jumps_every,
    std::size_t sampling_interval, std::size_t max_number_of_values_per_iteration) {
    using namespace macrodr;

    return std::tuple(path, filename, average_the_agonist_evolution, num_scouts_per_ensemble,
                      number_trials_until_give_up, thermo_jumps_every, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
                      stops_at, includes_zero, random_jumps, sampling_interval,
                      max_number_of_values_per_iteration);
}

using cuevi_algo_type = typename return_type<std::decay_t<decltype(&set_CueviAlgorithm)>>::type;

inline void calc_likelihood_old(std::string outfilename, std::string model,
                                parameters_value_type par, likelihood_algo_type likelihood,
                                std::string recording, experiment_type experiment,
                                cuevi_algo_type algorithm, tablefun_value_type ft) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every, save_every)();

    auto Maybe_model_v = get_model(model);
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [outfilename, &par, &ftbl3, &experiment, &recording, &likelihood,
             &algorithm](auto model0ptr) {
                auto& model0 = *model0ptr;

                auto [path, filename, average_the_agonist_evolution, num_scouts_per_ensemble,
                      number_trials_until_give_up, thermo_jump_factor, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
                      stops_at, includes_zero, random_jumps, sampling_interval,
                      max_number_of_values_per_iteration] = std::move(algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction_approximation, variance_correction, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                auto Maybe_param1 = var::load_Parameters(par.first, par.second, model0.model_name(),
                                                         model0.names());
                Simulated_Recording<var::please_include<>> y;
                auto Maybe_y = load_simulation(recording, ",", y);
                if (!Maybe_param1.valid() || !Maybe_y.valid()) {
                    std::cerr << "---------ERROR_______________\n";
                    std::cerr << Maybe_param1.error()() << Maybe_y.error()();
                    std::cerr << "---------ERROR_______________\n";

                } else {
                    auto param1 = Maybe_param1.value().standard_parameter();

                    auto modelLikelihood = make_Likelihood_Model<
                        uses_adaptive_aproximation<true>, uses_recursive_aproximation<true>,
                        uses_averaging_aproximation<2>, uses_variance_aproximation<false>,
                        uses_taylor_variance_correction_aproximation<false>>(
                        model0, Simulation_n_sub_dt(n_sub_dt));
                    auto lik =
                        Macro_DMR{}
                            .log_Likelihood<
                                uses_adaptive_aproximation<true>, uses_recursive_aproximation<true>,
                                uses_averaging_aproximation<2>, uses_variance_aproximation<false>,
                                uses_taylor_variance_correction_aproximation<false>,
                                return_predictions<2>>(ftbl3, model0, param1, get<Recording>(y()),
                                                       experiment);
                    if (lik)
                        save_Likelihood_Predictions(outfilename, lik.value(), y, experiment);
                    else {
                        std::cerr << "---------ERROR_______________\n";
                        std::cerr << lik.error()();
                        std::cerr << "---------ERROR_______________\n";
                    }
                }
            },
            model_v);
    }
}

inline Maybe_error<std::string> calc_fraction_likelihood(
    std::string file_name, std::string model, parameters_value_type par,
    fractioned_simulation_type Maybe_frac_simulation, likelihood_algo_type likelihood_algo,
    tablefun_value_type ft) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    auto save_every = num_scouts_per_ensemble;

    auto ftbl3 = get_function_Table_maker_St(filename, save_every, save_every)();

    auto Maybe_model_v = get_model(model);
    if (!Maybe_model_v)
        return Maybe_model_v.error()();
    auto model_v = std::move(Maybe_model_v.value());
    if (!Maybe_frac_simulation)
        return Maybe_frac_simulation.error()();
    auto frac_simulation = std::move(Maybe_frac_simulation.value());
    return std::visit(
        [&ftbl3, &frac_simulation, &likelihood_algo, &par,
         &file_name](auto model0ptr) -> Maybe_error<std::string> {
            auto& model0 = *model0ptr;
            auto [experiment, simulation, fs, iniagonist] = frac_simulation;

            auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                  variance_correction, variance_correction_approximation, n_sub_dt] =
                likelihood_algo;

            using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

            std::string ModelName = model0.model_name();
            auto Maybe_param1 =
                var::load_Parameters(par.first, par.second, model0.model_name(), model0.names());

            std::vector<Experiment> xs;
            std::vector<Simulated_Recording<var::please_include<N_Ch_State_Evolution>>> ys;
            auto Maybe_e = load_fractioned_experiment(experiment, ",", fs, iniagonist, xs);

            auto Maybe_y = load_fractioned_simulation(simulation, ",", ys);

            if (xs.size() != ys.size())
                return error_message(
                    "number of fractions mismatch between "
                    "recordings and experiments");

            if (!(Maybe_e.valid() && Maybe_y.valid() && Maybe_param1.valid()))
                return error_message(Maybe_e.error()() + Maybe_y.error()() +
                                     Maybe_param1.error()());

            auto param1 = Maybe_param1.value().standard_parameter();

            std::string filename = file_name + "_" + ModelName + "_" + time_now();

            auto maybe_modelLikelihood =
                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              true>,
                    decltype(model0)>(
                    model0, Simulation_n_sub_dt(n_sub_dt),
                    uses_adaptive_aproximation_value(adaptive_aproximation),
                    uses_recursive_aproximation_value(recursive_approximation),
                    uses_averaging_aproximation_value(averaging_approximation),
                    uses_variance_aproximation_value(variance_correction),
                    uses_variance_correction_aproximation_value(variance_correction_approximation))
                    .get_variant();
            if (!maybe_modelLikelihood) {
                return maybe_modelLikelihood.error()();
            }
            auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

            return std::visit(
                [&ftbl3, &param1, &ys, &xs,
                 &filename](auto& modelLikelihood) -> Maybe_error<std::string> {
                    auto Maybe_lik =
                        fractioned_logLikelihoodPredictions(ftbl3, modelLikelihood, param1, ys, xs);
                    if (!Maybe_lik)
                        return Maybe_lik.error();
                    else {
                        save_fractioned_Likelihood_Predictions(filename + "frac_likelihood.csv",
                                                               ",", Maybe_lik.value(), ys, xs);
                        return filename;
                    }
                },
                modelLikelihood_v);
        },
        model_v);
}

inline void calc_fraction_evidence(std::string model, prior_value_type prior,
                                   likelihood_algo_type likelihood,
                                   fractioned_simulation_type Maybe_frac_experiment,
                                   cuevi_algo_type cuevi_algorithm, tablefun_value_type ft,
                                   std::size_t myseed) {
    using namespace macrodr;
    if (!Maybe_frac_experiment)
        return;
    auto frac_experiment = Maybe_frac_experiment.value();

    auto [filename, num_scouts_per_ensemble] = std::move(ft);

    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every, save_every)();

    auto Maybe_model_v = get_model(model);

    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &frac_experiment, &prior, &likelihood, &cuevi_algorithm,
             &myseed](auto model0ptr) {
                auto& model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                auto [fname_experiment, fname_simulation, fs, iniagonist] = frac_experiment;
                std::vector<Experiment> xs;
                std::vector<Recording> ys;
                auto Maybe_e = load_fractioned_experiment(fname_experiment, ",", fs, iniagonist, xs);

                auto Maybe_ys = load_fractioned_Recording(fname_simulation, ",", ys);

                auto [path, file_name, average_the_agonist_evolution, num_scouts_per_ensemble,
                      number_trials_until_give_up, thermo_jumps_every, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
                      stops_at, includes_zero, random_jumps, sampling_interval,
                      max_number_of_values_per_iteration] = std::move(cuevi_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior.first, prior.second, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    if (Maybe_ys.valid() && Maybe_e.valid()) {
                        std::string filename = file_name + "_" + ModelName + "_" + time_now() +
                                               "_" + std::to_string(myseed);

                        auto saving_itervals = Saving_intervals(Vector_Space(
                            Save_Evidence_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Likelihood_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Parameter_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_RateParameter_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Predictions_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration))));

                        auto cbc = new_cuevi_Model_already_fraction_by_iteration(
                            path, filename, num_scouts_per_ensemble, number_trials_until_give_up,
                            thermo_jumps_every, max_iter_equilibrium, n_points_per_decade_beta_high,
                            n_points_per_decade_beta_low, medium_beta, stops_at, includes_zero,
                            saving_itervals, random_jumps);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_taylor_variance_correction_aproximation, true>,
                                decltype(model0)>(
                                model0, Simulation_n_sub_dt(n_sub_dt),
                                uses_adaptive_aproximation_value(adaptive_aproximation),
                                uses_recursive_aproximation_value(recursive_approximation),
                                uses_averaging_aproximation_value(averaging_approximation),
                                uses_variance_aproximation_value(variance_correction),
                                uses_variance_correction_aproximation_value(
                                    variance_correction_approximation))
                                .get_variant();
                        if (!maybe_modelLikelihood)
                            std::cerr << maybe_modelLikelihood.error()();
                        else {
                            auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());
                            // using m2=typename
                            // decltype(modelLikelihood_v)::paseModelLikelihoodv;

                            std::visit(
                                [&ftbl3, &cbc, &param1_prior, &ys, &xs,
                                 &myseed](auto& modelLikelihood) {
                                    auto opt3 = cuevi::evidence_fraction(
                                        ftbl3, std::move(cbc), param1_prior, modelLikelihood, ys,
                                        xs, cuevi::Init_seed(myseed));
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

inline void calc_evidence(std::string model, prior_value_type prior,
                          likelihood_algo_type likelihood, std::string recording,
                          experiment_type experiment, fraction_algo_type fraction_algo,
                          cuevi_algo_type cuevi_algorithm, tablefun_value_type ft,
                          std::size_t myseed) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);

    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every, save_every)();

    auto Maybe_model_v = get_model(model);

    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &experiment, &recording, &prior, &likelihood, &fraction_algo, &cuevi_algorithm,
             &myseed](auto model0ptr) {
                auto& model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);

                auto [min_fraction, n_points_per_decade_fraction, segments] =
                    std::move(fraction_algo);
                auto [path, file_name, average_the_agonist_evolution, num_scouts_per_ensemble,
                      number_trials_until_give_up, thermo_jumps_every, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
                      stops_at, includes_zero, random_jumps, sampling_interval,
                      max_number_of_values_per_iteration] = std::move(cuevi_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior.first, prior.second, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        std::string filename = file_name + "_" + ModelName + "_" + time_now() +
                                               "_" + std::to_string(myseed);

                        auto Maybe_t_segments_used =
                            load_segments_length_for_fractioning(segments, ",");

                        if (Maybe_t_segments_used) {
                            auto t_segments_used = std::move(Maybe_t_segments_used.value());

                            auto saving_itervals = Saving_intervals(Vector_Space(
                                Save_Evidence_every(std::pair(sampling_interval,
                                                              max_number_of_values_per_iteration)),
                                Save_Likelihood_every(std::pair(
                                    sampling_interval, max_number_of_values_per_iteration)),
                                Save_Parameter_every(std::pair(sampling_interval,
                                                               max_number_of_values_per_iteration)),
                                Save_RateParameter_every(std::pair(
                                    sampling_interval, max_number_of_values_per_iteration)),
                                Save_Predictions_every(std::pair(
                                    sampling_interval, max_number_of_values_per_iteration))));

                            auto cbc = new_cuevi_Model_by_iteration(
                                path, filename, t_segments_used, average_the_agonist_evolution,
                                num_scouts_per_ensemble, number_trials_until_give_up, min_fraction,
                                thermo_jumps_every, max_iter_equilibrium,
                                n_points_per_decade_beta_low, n_points_per_decade_fraction,
                                medium_beta, stops_at, includes_zero, saving_itervals,
                                random_jumps);

                            auto maybe_modelLikelihood =
                                Likelihood_Model_regular<
                                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation,
                                                              true>,
                                    var::constexpr_Var_domain<bool, uses_recursive_aproximation,
                                                              true>,
                                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                    var::constexpr_Var_domain<bool, uses_variance_aproximation,
                                                              true>,
                                    var::constexpr_Var_domain<
                                        bool, uses_taylor_variance_correction_aproximation, true>,
                                    decltype(model0)>(
                                    model0, Simulation_n_sub_dt(n_sub_dt),
                                    uses_adaptive_aproximation_value(adaptive_aproximation),
                                    uses_recursive_aproximation_value(recursive_approximation),
                                    uses_averaging_aproximation_value(averaging_approximation),
                                    uses_variance_aproximation_value(variance_correction),
                                    uses_variance_correction_aproximation_value(
                                        variance_correction_approximation))
                                    .get_variant();
                            if (!maybe_modelLikelihood) {
                                std::cerr << maybe_modelLikelihood.error()();
                                return;
                            }
                            auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());
                            // using m2=typename
                            // decltype(modelLikelihood_v)::paseModelLikelihoodv;

                            std::visit(
                                [&ftbl3, &cbc, &param1_prior, &y, &experiment,
                                 &myseed](auto& modelLikelihood) {
                                    auto opt3 = cuevi::evidence(ftbl3, std::move(cbc), param1_prior,
                                                                modelLikelihood, y, experiment,
                                                                cuevi::Init_seed(myseed));
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

inline void calc_evidence_continuation(std::string model, prior_value_type prior,
                                       likelihood_algo_type likelihood, std::string recording,
                                       experiment_type experiment, fraction_algo_type fraction_algo,
                                       cuevi_algo_type algorithm, tablefun_value_type ft,
                                       std::size_t myseed) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every, save_every)();

    auto Maybe_model_v = get_model(model);
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &experiment, &recording, &prior, &likelihood, &fraction_algo, &algorithm,
             &myseed](auto model0ptr) {
                auto& model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                /*
* path, filename, average_the_agonist_evolution,
num_scouts_per_ensemble, number_trials_until_give_up,
thermo_jumps_every, max_iter_equilibrium,
n_points_per_decade, medium_beta, stops_at, includes_zero,
saving_itervals, random_jumps
* */
                auto [min_fraction, n_points_per_decade_fraction, segments] =
                    std::move(fraction_algo);

                auto [path, file_name, average_the_agonist_evolution, num_scouts_per_ensemble,
                      number_trials_until_give_up, thermo_jumps_every, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
                      stops_at, includes_zero, random_jumps, sampling_interval,
                      max_number_of_values_per_iteration] = std::move(algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction_approximation, variance_correction, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior.first, prior.second, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
                        auto number_of_traces = 7;
                        auto number_of_segments = t_segments.size();
                        t_segments.reserve(number_of_traces * t_segments.size());

                        for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
                            std::copy_n(t_segments.begin(), number_of_segments,
                                        std::back_inserter(t_segments));

                        std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};

                        bool average_the_agonist_evolution = true;

                        /**
               * @brief cbc cumulative evidence algorithm, ends using
               * convergence criteria
               */

                        std::string filename = file_name + "_" + ModelName + "_" + time_now() +
                                               "_" + std::to_string(myseed);

                        auto& t_segments_used = t_segments_7;

                        auto saving_itervals = Saving_intervals(Vector_Space(
                            Save_Evidence_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Likelihood_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Parameter_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_RateParameter_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration)),
                            Save_Predictions_every(
                                std::pair(sampling_interval, max_number_of_values_per_iteration))));

                        auto cbc = new_cuevi_Model_by_iteration(
                            path, filename, t_segments_used, average_the_agonist_evolution,
                            num_scouts_per_ensemble, number_trials_until_give_up, min_fraction,
                            thermo_jumps_every, max_iter_equilibrium, n_points_per_decade_beta_low,
                            n_points_per_decade_fraction, medium_beta, stops_at, includes_zero,
                            saving_itervals, random_jumps);

                        auto modelLikelihood = make_Likelihood_Model<
                            uses_adaptive_aproximation<true>, uses_recursive_aproximation<true>,
                            uses_averaging_aproximation<2>, uses_variance_aproximation<false>,
                            uses_taylor_variance_correction_aproximation<false>>(
                            model0, Simulation_n_sub_dt(n_sub_dt));
                        auto opt3 =
                            cuevi::evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood, y,
                                            experiment, cuevi::Init_seed(myseed));
                    }
                }
            },
            model_v);
    }
}

inline dsl::Compiler make_cuevi_compiler() {
    dsl::Compiler cm;
    cm.push_function(
        "set_CueviAlgorithm",
        dsl::to_typed_function<std::size_t, std::size_t, double, double, bool, bool, std::size_t,
                               std::string, double, double, std::size_t, std::string, std::size_t,
                               std::size_t, std::size_t>(
            &set_CueviAlgorithm, "num_scouts_per_ensemble", "number_trials_until_give_up",
            "stops_at", "medium_beta", "includes_zero", "random_jumps", "max_iter_equilibrium",
            "path", "n_points_per_decade_beta_high", "n_points_per_decade_beta_low",
            "average_the_agonist_evolution", "filename", "thermo_jumps_every", "sampling_interval",
            "max_number_of_values_per_iteration"));
    // Si hay comandos tipo evidence_fraction con cuevi, van aquí.
    return cm;
}

}  // namespace zombie

#endif
}  // namespace macrodr::cmd

#endif  // CLI_FUNCTION_TABLE_H
