#ifndef CLI_THERMO_EVIDENCE_DTS_H
#define CLI_THERMO_EVIDENCE_DTS_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include "CLI_function_table.h"
#include "experiment.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "qmodel.h"

namespace macrodr {
namespace cmd {

inline auto set_ThermoAlgorithm_dts(
    std::size_t num_scouts_per_ensemble, std::size_t number_trials_until_give_up,
    std::size_t max_iter_equilibrium, std::size_t beta_size, std::size_t thermo_jumps_every,
    std::size_t save_every_param_size_factor, std::size_t t_adapt_beta_every,
    std::string t_adapt_beta_equalizer, std::string t_adapt_beta_controler,
    std::string t_adapt_beta_variance, double t_adapt_beta_nu, double t_adapt_beta_t0,
    bool t_adjust_beta, double t_acceptance_upper_limit, double t_acceptance_lower_limit,
    double t_desired_acceptance) {
    using namespace macrodr;

    return std::tuple(num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance);
}

using thermo_algo_dts_type =
    typename return_type<std::decay_t<decltype(&set_ThermoAlgorithm_dts)>>::type;

inline void calc_thermo_evidence_dts(std::string id, std::string model, std::string prior,
                                     likelihood_algo_type likelihood, std::string recording,
                                     experiment_file_type experiment_file,
                                     thermo_algo_dts_type thermo_algorithm,
                                     std::size_t sampling_interval,
                                     std::size_t max_number_of_values_per_iteration,
                                     std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);

    if (true) {
        std::ofstream f("thermo_evidence_dts_" + id + ".txt");
        save_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, sampling_interval, max_number_of_values_per_iteration, myseed);
    }
    using namespace macrodr;
    auto experiment = get_Experiment(std::get<0>(experiment_file), std::get<1>(experiment_file),
                                     std::get<2>(experiment_file));

    auto ftbl3 = cmd::get_function_Table_maker_St(filename, sampling_interval,
                                                  max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (!Maybe_model_v) {
        std::cerr << Maybe_model_v.error()();
    } else {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(

            [&filename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &myseed, sampling_interval, max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;
                mt_64i mt(myseed);

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance] =
                    std::move(thermo_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        auto saving_intervals = Saving_intervals(Vector_Space(
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

                        auto tmi = new_thermo_Model_by_max_iter_dts(
                            "", filename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, saving_intervals, myseed,
                            t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                            t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0,
                            t_adjust_beta, t_acceptance_upper_limit, t_acceptance_lower_limit,
                            t_desired_acceptance);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_variance_correction_aproximation, true>,
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
                            [&ftbl3, &tmi, &param1_prior, &y, &experiment](auto& modelLikelihood) {
                                auto opt =
                                    thermo_evidence<true>(ftbl3, std::move(tmi), param1_prior,
                                                          modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}

inline void calc_thermo_evidence_dts_2(std::string id, std::string model, std::string prior,
                                       likelihood_algo_type likelihood, std::string recording,
                                       experiment_file_type experiment_file,
                                       thermo_algo_dts_type thermo_algorithm,
                                       std::size_t sampling_interval,
                                       std::size_t max_number_of_values_per_iteration,
                                       std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);

    if (true) {
        std::ofstream f("thermo_evidence_dts_" + id + ".txt");
        save_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, sampling_interval, max_number_of_values_per_iteration, myseed);
    }
    using namespace macrodr;
    auto experiment = get_Experiment(std::get<0>(experiment_file), std::get<1>(experiment_file),
                                     std::get<2>(experiment_file));

    auto ftbl3 = cmd::get_function_Table_maker_St_no_Qdt_memoization(
        filename, sampling_interval, max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (!Maybe_model_v) {
        std::cerr << Maybe_model_v.error()();
    } else {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(

            [&filename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &myseed, sampling_interval, max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;
                mt_64i mt(myseed);

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance] =
                    std::move(thermo_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        auto saving_intervals = Saving_intervals(Vector_Space(
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

                        auto tmi = new_thermo_Model_by_max_iter_dts(
                            "", filename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, saving_intervals, myseed,
                            t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                            t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0,
                            t_adjust_beta, t_acceptance_upper_limit, t_acceptance_lower_limit,
                            t_desired_acceptance);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_variance_correction_aproximation, true>,
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
                            [&ftbl3, &tmi, &param1_prior, &y, &experiment](auto& modelLikelihood) {
                                auto opt =
                                    thermo_evidence<true>(ftbl3, std::move(tmi), param1_prior,
                                                          modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}

inline void calc_thermo_evidence_dts_continuation(std::string id, std::size_t ith,
                                                  std::size_t myinit_seed) {
    std::string model;
    std::string prior;
    likelihood_algo_type likelihood;
    std::string recording;
    experiment_file_type experiment_file;
    thermo_algo_dts_type thermo_algorithm;
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;

    std::size_t myseed;
    std::string filename;
    if (true) {
        std::ifstream f("thermo_evidence_dts_" + id + ".txt");
        if (!f) {
            std::cerr << "Error!!\n thermo_evidence_dts_" + id + ".txt" + " not found";
            return;
        }
        load_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, sampling_interval, max_number_of_values_per_iteration, myseed);
    }
    auto experiment = get_Experiment(std::get<0>(experiment_file), std::get<1>(experiment_file),
                                     std::get<2>(experiment_file));

    std::string oldfilename = filename;
    if (ith > 1)
        oldfilename = filename + '_' + std::to_string(ith - 1);

    std::string newfilename = filename + '_' + std::to_string(ith);
    using namespace macrodr;

    auto ftbl3 = get_function_Table_maker_St(newfilename, sampling_interval,
                                             max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&oldfilename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &myinit_seed, &newfilename, sampling_interval,
             max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;
                auto myseed = calc_seed(myinit_seed);
                mt_64i mt(myseed);

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance] =
                    std::move(thermo_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        auto saving_intervals = Saving_intervals(Vector_Space(
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

                        auto tmi = new_thermo_Model_by_max_iter_dts(
                            "", newfilename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, saving_intervals, myseed,
                            t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                            t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0,
                            t_adjust_beta, t_acceptance_upper_limit, t_acceptance_lower_limit,
                            t_desired_acceptance);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_variance_correction_aproximation, true>,
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
                            [&oldfilename, &ftbl3, &tmi, &param1_prior, &y,
                             &experiment](auto& modelLikelihood) {
                                auto opt = thermo_evidence_continuation<true>(
                                    oldfilename, ftbl3, std::move(tmi), param1_prior,
                                    modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}

inline void calc_thermo_evidence_dts_continuation_2(std::string id, std::size_t ith,
                                                    std::size_t myinit_seed) {
    std::string model;
    std::string prior;
    likelihood_algo_type likelihood;
    std::string recording;
    experiment_file_type experiment_file;
    thermo_algo_dts_type thermo_algorithm;
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;

    std::size_t myseed;
    std::string filename;
    if (true) {
        std::ifstream f("thermo_evidence_dts_" + id + ".txt");
        if (!f) {
            std::cerr << "Error!!\n thermo_evidence_dts_" + id + ".txt" + " not found";
            return;
        }
        load_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, sampling_interval, max_number_of_values_per_iteration, myseed);
    }
    auto experiment = get_Experiment(std::get<0>(experiment_file), std::get<1>(experiment_file),
                                     std::get<2>(experiment_file));

    std::string oldfilename = filename;
    if (ith > 1)
        oldfilename = filename + '_' + std::to_string(ith - 1);

    std::string newfilename = filename + '_' + std::to_string(ith);
    using namespace macrodr;

    auto ftbl3 = get_function_Table_maker_St_no_Qdt_memoization(
        newfilename, sampling_interval, max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&oldfilename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &myinit_seed, &newfilename, sampling_interval,
             max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;
                auto myseed = calc_seed(myinit_seed);
                mt_64i mt(myseed);

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance] =
                    std::move(thermo_algorithm);

                auto [adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction, variance_correction_approximation, n_sub_dt] =
                    likelihood;

                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

                std::string ModelName = model0.model_name();

                auto Maybe_param1_prior =
                    var::load_Prior(prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());

                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        auto saving_intervals = Saving_intervals(Vector_Space(
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

                        auto tmi = new_thermo_Model_by_max_iter_dts(
                            "", newfilename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, saving_intervals, myseed,
                            t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                            t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0,
                            t_adjust_beta, t_acceptance_upper_limit, t_acceptance_lower_limit,
                            t_desired_acceptance);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, true>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_variance_correction_aproximation, true>,
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
                            [&oldfilename, &ftbl3, &tmi, &param1_prior, &y,
                             &experiment](auto& modelLikelihood) {
                                auto opt = thermo_evidence_continuation<true>(
                                    oldfilename, ftbl3, std::move(tmi), param1_prior,
                                    modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}
inline dsl::Compiler make_dts_compiler() {
    dsl::Compiler cm;
    cm.push_function(
        "set_ThermoAlgorithm_dts",
        dsl::to_typed_function<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                               std::size_t, std::size_t, std::string, std::string, std::string,
                               double, double, bool, double, double, double>(
            &set_ThermoAlgorithm_dts, "num_scouts_per_ensemble", "number_trials_until_give_up",
            "max_iter_equilibrium", "beta_size", "thermo_jumps_every",
            "save_every_param_size_factor", "adapt_beta_every", "adapt_beta_equalizer",
            "adapt_beta_controler", "adapt_beta_variance", "adapt_beta_nu", "adapt_beta_t0",
            "adjust_beta", "acceptance_upper_limit", "acceptance_lower_limit",
            "desired_acceptance"));

    cm.push_function(
        "thermo_evidence_dts",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, experiment_file_type, thermo_algo_dts_type, std::size_t,
                               std::size_t, std::size_t>(
            &calc_thermo_evidence_dts, "idname", "model", "prior", "likelihood_algorithm", "data",
            "experiment", "thermo_algorithm", "sampling_interval",
            "max_number_of_values_per_iteration", "init_seed"));
    cm.push_function(
        "thermo_evidence_dts_2",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, experiment_file_type, thermo_algo_dts_type, std::size_t,
                               std::size_t, std::size_t>(
            &calc_thermo_evidence_dts_2, "idname", "model", "prior", "likelihood_algorithm", "data",
            "experiment", "thermo_algorithm", "sampling_interval",
            "max_number_of_values_per_iteration", "init_seed"));
    cm.push_function(
        "thermo_evidence_dts_continuation",
        dsl::to_typed_function<std::string, std::size_t, std::size_t>(
            &calc_thermo_evidence_dts_continuation, "idname", "continuation_number", "init_seed"));
    cm.push_function("thermo_evidence_dts_continuation_2",
                     dsl::to_typed_function<std::string, std::size_t, std::size_t>(
                         &calc_thermo_evidence_dts_continuation_2, "idname", "continuation_number",
                         "init_seed"));
    return cm;
}

}  // namespace cmd
}  // namespace macrodr

#endif  // CLI_THERMO_EVIDENCE_DTS_H
