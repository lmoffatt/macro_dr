#ifndef CLI_FUNCTION_TABLE_H
#define CLI_FUNCTION_TABLE_H

#include "CLI_macro_dr_base.h"
#include "cuevi.h"
#include "maybe_error.h"
#include "models_MoffattHume_allosteric.h"
#include "parameters.h"
#include "qmodel.h"
#include <cstddef>
#include <string>
namespace macrodr {
namespace cmd {

inline auto set_CueviAlgorithm(
    std::size_t num_scouts_per_ensemble,
    std::size_t number_trials_until_give_up, double stops_at,
    double medium_beta, bool includes_zero, bool random_jumps,
    std::size_t max_iter_equilibrium, std::string path,
    double n_points_per_decade_beta_high, double n_points_per_decade_beta_low,
    bool average_the_ATP_evolution, std::string filename,
    std::size_t thermo_jumps_every, std::size_t save_every_param_factor) {
    using namespace macrodr;
    
    return std::tuple(
        path, filename, average_the_ATP_evolution, num_scouts_per_ensemble,
        number_trials_until_give_up, thermo_jumps_every, max_iter_equilibrium,
        n_points_per_decade_beta_high, n_points_per_decade_beta_low, medium_beta,
        stops_at, includes_zero, random_jumps, save_every_param_factor);
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
            //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step,
            //     double>{}, num_scouts_per_ensemble / 2),
            var::Single_Thread_Memoizer(
                var::F(Calc_Qdt_step{},
                                               [](auto &&...x) {
                                                   auto m = Macro_DMR{};
                                            return m.calc_Qdt_ATP_step(
                                                       std::forward<decltype(x)>(x)...);
                                        }),
                                        var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{}),
            var::Single_Thread_Memoizer(
                var::F(Calc_Qdtm_step{},
                                               [](auto &&...x) {
                                                   auto m = Macro_DMR{};
                                            return m.calc_Qdtm_ATP_step(
                                                       std::forward<decltype(x)>(x)...);
                                        }),
                                        var::Memoiza_all_values<Maybe_error<Qdtm>, ATP_step, double>{}),
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
            //           return
            //           m.calc_eigen(std::forward<decltype(x)>(x)...);
            //       }))
            
            )
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation(false)>,
            //                        ::V<uses_recursive_aproximation(true)>>,
            //         in_progress::S<::V<uses_averaging_aproximation(2)>>,
            //         in_progress::S<::V<uses_variance_aproximation(true)>>,
            //         in_progress::S<
            //             ::V<uses_variance_correction_aproximation(false)>>>{});
            // .append_Fs<MacroR2>(
            //     in_progress::P<
            //         in_progress::S<::V<uses_recursive_aproximation(false)>,
            //                        ::V<uses_recursive_aproximation(true)>>,
            //         in_progress::S<::V<uses_averaging_aproximation(0)>,
            //                        ::V<uses_averaging_aproximation(1)>,
            //                        ::V<uses_averaging_aproximation(2)>>,
            //         in_progress::S<::V<uses_variance_aproximation(false)>,
            //                        ::V<uses_variance_aproximation(true)>>,
            //         in_progress::S<
            //             ::V<uses_variance_correction_aproximation(false)>,
            //             ::V<uses_variance_correction_aproximation(true)>>>{})
            ;
    };
}

using tablefun_value_type = typename return_type<
    std::decay_t<decltype(&get_function_Table_maker_value)>>::type;

using cuevi_algo_type =
    typename return_type<std::decay_t<decltype(&set_CueviAlgorithm)>>::type;

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
                        auto param1 = Maybe_parameter_values.value().standard_parameter();
                        save_Parameter<var::Parameters_transformed<MyModel>> s(filename,
                                                                               1);
                        if (!includeN) {
                            auto sim = Macro_DMR{}.sample(
                                mt, model0, param1, experiment,
                                Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)),
                                recording);
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
                                Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)),
                                recording);
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

inline void calc_likelihood(std::string outfilename, std::string model,
                            parameters_value_type par,
                            likelihood_algo_type likelihood,
                            recording_value_type recording,
                            experiment_type experiment,
                            cuevi_algo_type algorithm, tablefun_value_type ft) {
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
                
                auto [path, filename, average_the_ATP_evolution,
                      num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low,
                      medium_beta, stops_at, includes_zero, random_jumps,
                      save_every_param_factor] = std::move(algorithm);
                
                auto [adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction_approximation,
                      variance_correction, n_sub_dt] = likelihood;
                
                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                
                auto Maybe_param1 = var::load_Parameters<MyModel>(
                    par.first, par.second, model0.model_name(), model0.names());
                Simulated_Recording<includes_N_state_evolution(true)> y;
                auto Maybe_y = load_simulation(recording.first, recording.second, y);
                if (!Maybe_param1.valid() || !Maybe_y.valid()) {
                    std::cerr << "---------ERROR_______________\n";
                    std::cerr << Maybe_param1.error()() << Maybe_y.error()();
                    std::cerr << "---------ERROR_______________\n";
                    
                } else {
                    
                    auto param1 = Maybe_param1.value().standard_parameter();
                    
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
    fractioned_simulation_type Maybe_frac_simulation,
    likelihood_algo_type likelihood_algo, tablefun_value_type ft) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    auto save_every = num_scouts_per_ensemble;
    
    auto ftbl3 = get_function_Table_maker_St(filename, save_every)();
    
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
            auto &model0 = *model0ptr;
            auto [experiment, simulation, fs, iniATP] = frac_simulation;
            
            auto [adaptive_aproximation, recursive_approximation,
                  averaging_approximation, variance_correction,
                  variance_correction_approximation, n_sub_dt] = likelihood_algo;
            
            using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
            
            std::string ModelName = model0.model_name();
            auto Maybe_param1 = var::load_Parameters<MyModel>(
                par.first, par.second, model0.model_name(), model0.names());
            
            std::vector<Experiment> xs;
            std::vector<Simulated_Recording<includes_N_state_evolution(true)>> ys;
            auto Maybe_e =
                load_fractioned_experiment(experiment, ",", fs, iniATP, xs);
            
            auto Maybe_y = load_fractioned_simulation(simulation, ",", ys);
            
            if (xs.size() != ys.size())
                return error_message("number of fractions mismatch between "
                                     "recordings and experiments");
            
            if (!(Maybe_e.valid() && Maybe_y.valid() && Maybe_param1.valid()))
                return error_message(Maybe_e.error()() + Maybe_y.error()() +
                                     Maybe_param1.error()());
            
            auto param1 = Maybe_param1.value().standard_parameter();
            
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
                 &filename](auto &modelLikelihood) -> Maybe_error<std::string> {
                    auto Maybe_lik = fractioned_logLikelihoodPredictions(
                        ftbl3, modelLikelihood, param1, ys, xs);
                    if (!Maybe_lik)
                        return Maybe_lik.error();
                    else {
                        save_fractioned_Likelihood_Predictions(
                            filename + "frac_likelihood.csv", ",", Maybe_lik.value(),
                            ys, xs);
                        return filename;
                    }
                },
                modelLikelihood_v);
        },
        model_v);
}

inline void
calc_fraction_evidence(std::string model, prior_value_type prior,
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
    auto ftbl3 = get_function_Table_maker_St(filename, save_every)();
    
    auto Maybe_model_v = get_model(model);
    
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &frac_experiment, &prior, &likelihood, &cuevi_algorithm,
             &myseed](auto model0ptr) {
                auto &model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                auto [fname_experiment, fname_simulation, fs, iniATP] =
                    frac_experiment;
                std::vector<Experiment> xs;
                std::vector<Recording> ys;
                auto Maybe_e =
                    load_fractioned_experiment(fname_experiment, ",", fs, iniATP, xs);
                
                auto Maybe_ys = load_fractioned_Recording(fname_simulation, ",", ys);
                
                auto [path, file_name, average_the_ATP_evolution,
                      num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low,
                      medium_beta, stops_at, includes_zero, random_jumps,
                      save_every_param_factor] = std::move(cuevi_algorithm);
                
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
                    
                    if (Maybe_ys.valid() && Maybe_e.valid()) {
                        
                        std::string filename = file_name + "_" + ModelName + "_" +
                                               time_now() + "_" + std::to_string(myseed);
                        
                        auto saving_itervals = Saving_intervals(Vector_Space(
                            Save_Evidence_every(save_every_param_factor *
                                                param1_prior.size()),
                            Save_Likelihood_every(save_every_param_factor *
                                                  param1_prior.size()),
                            Save_Parameter_every(save_every_param_factor *
                                                 param1_prior.size()*4),
                            Save_Predictions_every(save_every_param_factor *
                                                   param1_prior.size() * 500)));
                        
                        auto cbc = new_cuevi_Model_already_fraction_by_iteration<MyModel>(
                            path, filename, num_scouts_per_ensemble,
                            number_trials_until_give_up, thermo_jumps_every,
                            max_iter_equilibrium, n_points_per_decade_beta_high,
                            n_points_per_decade_beta_low, medium_beta, stops_at,
                            includes_zero, saving_itervals, random_jumps);
                        
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
                        
                    } else
                        std::cerr << Maybe_ys.error()();
                }
            },
            model_v);
    }
}

inline void calc_evidence(std::string model, prior_value_type prior,
                          likelihood_algo_type likelihood,
                          std::string recording, experiment_type experiment,
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
            [&ftbl3, &experiment, &recording, &prior, &likelihood, &fraction_algo,
             &cuevi_algorithm, &myseed](auto model0ptr) {
                auto &model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                
                auto [min_fraction, n_points_per_decade_fraction, segments] =
                    std::move(fraction_algo);
                auto [path, file_name, average_the_ATP_evolution,
                      num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low,
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
                                Save_Predictions_every(num_scouts_per_ensemble * 500)));
                            
                            auto cbc = new_cuevi_Model_by_iteration<MyModel>(
                                path, filename, t_segments_used, average_the_ATP_evolution,
                                num_scouts_per_ensemble, number_trials_until_give_up,
                                min_fraction, thermo_jumps_every, max_iter_equilibrium,
                                n_points_per_decade_beta_low, n_points_per_decade_fraction,
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

inline void calc_evidence_continuation(
    std::string model, prior_value_type prior, likelihood_algo_type likelihood,
    std::string recording, experiment_type experiment,
    fraction_algo_type fraction_algo, cuevi_algo_type algorithm,
    tablefun_value_type ft, std::size_t myseed) {
    using namespace macrodr;
    auto [filename, num_scouts_per_ensemble] = std::move(ft);
    auto save_every = num_scouts_per_ensemble;
    auto ftbl3 = get_function_Table_maker_St(filename, save_every)();
    
    auto Maybe_model_v = get_model(model);
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&ftbl3, &experiment, &recording, &prior, &likelihood, &fraction_algo,
             &algorithm, &myseed](auto model0ptr) {
                auto &model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                /*
     * path, filename, average_the_ATP_evolution,
              num_scouts_per_ensemble, number_trials_until_give_up,
              thermo_jumps_every, max_iter_equilibrium,
              n_points_per_decade, medium_beta, stops_at, includes_zero,
              saving_itervals, random_jumps
     * */
                auto [min_fraction, n_points_per_decade_fraction, segments] =
                    std::move(fraction_algo);
                
                auto [path, file_name, average_the_ATP_evolution,
                      num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium,
                      n_points_per_decade_beta_high, n_points_per_decade_beta_low,
                      medium_beta, stops_at, includes_zero, random_jumps,
                      save_every_param_factor] = std::move(algorithm);
                
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
                        
                        bool average_the_ATP_evolution = true;
                        
                        /**
               * @brief cbc cumulative evidence algorithm, ends using
               * convergence criteria
               */
                        
                        std::string filename = file_name + "_" + ModelName + "_" +
                                               time_now() + "_" + std::to_string(myseed);
                        
                        auto &t_segments_used = t_segments_7;
                        
                        auto saving_itervals = Saving_intervals(Vector_Space(
                            Save_Evidence_every(param1_prior.size() *
                                                save_every_param_factor),
                            Save_Likelihood_every(param1_prior.size() *
                                                  save_every_param_factor),
                            Save_Parameter_every(param1_prior.size() *
                                                 save_every_param_factor),
                            Save_Predictions_every(param1_prior.size() *
                                                   save_every_param_factor * 50)));
                        
                        auto cbc = new_cuevi_Model_by_iteration<MyModel>(
                            path, filename, t_segments_used, average_the_ATP_evolution,
                            num_scouts_per_ensemble, number_trials_until_give_up,
                            min_fraction, thermo_jumps_every, max_iter_equilibrium,
                            n_points_per_decade_beta_low, n_points_per_decade_fraction,
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

} // namespace cmd
} // namespace macrodr

#endif // CLI_FUNCTION_TABLE_H
