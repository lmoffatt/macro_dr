#ifndef CLI_THERMO_EVIDENCE_FRACTION_DTS_H
#define CLI_THERMO_EVIDENCE_FRACTION_DTS_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include "CLI_function_table.h"
#include "CLI_thermo_evidence_dts.h"
#include "experiment.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
#include "parallel_tempering_fraction.h"
#include "parallel_tempering_linear_regression.h"
#include "qmodel.h"
#include "report_thermo_fraction_evidence.h"

namespace macrodr {
namespace cmd {

inline auto set_ThermoAlgorithm_fraction_dts(
    std::size_t num_scouts_per_ensemble, std::size_t number_trials_until_give_up,
    std::size_t max_iter_equilibrium, std::size_t beta_size, std::size_t thermo_jumps_every,
    std::size_t save_every_param_size_factor, std::size_t t_adapt_beta_every,
    std::string t_adapt_beta_equalizer, std::string t_adapt_beta_controler,
    std::string t_adapt_beta_variance, double t_adapt_beta_nu, double t_adapt_beta_t0,
    double t_adapt_beta_threshold, bool t_adjust_beta, double t_acceptance_upper_limit,
    double t_acceptance_lower_limit, double t_desired_acceptance) {
    using namespace macrodr;

    return std::tuple(num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0,
                      t_adapt_beta_threshold, t_adjust_beta, t_acceptance_upper_limit,
                      t_acceptance_lower_limit, t_desired_acceptance);
}

using thermo_algo_fraction_dts_type =
    typename return_type<std::decay_t<decltype(&set_ThermoAlgorithm_fraction_dts)>>::type;

inline void calc_thermo_evidence_fraction_dts(std::string id, std::string model, std::string prior,
                                              likelihood_algo_type likelihood,
                                              fractioned_experiment_type Maybe_frac_experiment,
                                              thermo_algo_fraction_dts_type thermo_algorithm,
                                              std::size_t sampling_interval,
                                              std::size_t max_number_of_values_per_iteration,
                                              std::size_t myseed) {
    myseed = calc_seed(myseed);
    if (!Maybe_frac_experiment)
        return;
    auto frac_experiment = Maybe_frac_experiment.value();
    auto [fname_experiment, fname_simulation, fs, iniagonist] = frac_experiment;
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);

    if (true) {
        std::ofstream f("thermo_evidence_fraction_dts_" + id + ".txt");
        save_vars(f, filename, model, prior, likelihood, fname_simulation, fname_experiment, fs,
                  iniagonist, thermo_algorithm, sampling_interval, max_number_of_values_per_iteration,
                  myseed);
    }
    using namespace macrodr;
    std::vector<Experiment> xs;
    std::vector<Recording> ys;

    auto Maybe_e = load_fractioned_experiment(fname_experiment, ",", fs, iniagonist, xs);

    auto Maybe_ys = load_fractioned_Recording(fname_simulation, ",", ys);
    if (!Maybe_e.valid() || !Maybe_ys.valid()) {
        std::cerr << Maybe_e.error()() << Maybe_ys.error()();
        return;
    }

    auto ftbl3 = cmd::get_function_Table_maker_St(filename, sampling_interval,
                                                  max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (!Maybe_model_v) {
        std::cerr << Maybe_model_v.error()();
    } else {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(

            [&filename, &ftbl3, &xs, &ys, &prior, &likelihood, &thermo_algorithm, &myseed,
             sampling_interval, max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;
                mt_64i mt(myseed);

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0,
                      t_adapt_beta_threshold, t_adjust_beta, t_acceptance_upper_limit,
                      t_acceptance_lower_limit, t_desired_acceptance] = std::move(thermo_algorithm);

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
                        t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0,
                        t_adapt_beta_threshold, t_adjust_beta, t_acceptance_upper_limit,
                        t_acceptance_lower_limit, t_desired_acceptance);

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
                    if (!maybe_modelLikelihood) {
                        std::cerr << maybe_modelLikelihood.error()();
                        return
                    };
                    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());
                    // using m2=typename
                    // decltype(modelLikelihood_v)::paseModelLikelihoodv;

                    std::visit(
                        [&ftbl3, &tmi, &param1_prior, &ys, &xs](auto& modelLikelihood) {
                            auto opt = thermo_fraction_evidence<true>(
                                ftbl3, std::move(tmi), param1_prior, modelLikelihood, ys, xs);
                        },
                        modelLikelihood_v);
                }
            },
            model_v);
    }
}

inline void calc_thermo_evidence_fraction_dts_continuation(std::string id, std::size_t ith,
                                                           std::size_t myinit_seed) {
    std::string model;
    std::string prior;
    likelihood_algo_type likelihood;
    std::string recording;
    experiment_file_type experiment_file;
    thermo_algo_dts_type thermo_algorithm;
    std::size_t sampling_interval;
    std::size_t max_number_of_values_per_iteration;

    std::string fname_simulation;
    std::string fname_experiment;
    double fs;
    double iniagonist;

    std::size_t myseed;
    std::string filename;
    if (true) {
        std::ifstream f("thermo_evidence_fraction_dts_" + id + ".txt");
        if (!f) {
            std::cerr << "Error!!\n thermo_evidence_dts_" + id + ".txt" + " not found";
            return;
        }
        load_vars(f, filename, model, prior, likelihood, fname_simulation, fname_experiment, fs,
                  iniagonist, thermo_algorithm, sampling_interval, max_number_of_values_per_iteration,
                  myseed);
    }
    std::string oldfilename = filename;
    if (ith > 1)
        oldfilename = filename + '_' + std::to_string(ith - 1);

    std::string newfilename = filename + '_' + std::to_string(ith);

    using namespace macrodr;
    std::vector<Experiment> xs;
    std::vector<Recording> ys;

    auto Maybe_e = load_fractioned_experiment(fname_experiment, ",", fs, iniagonist, xs);

    auto Maybe_ys = load_fractioned_Recording(fname_simulation, ",", ys);
    if (!Maybe_e.valid() || !Maybe_ys.valid()) {
        std::cerr << Maybe_e.error()() << Maybe_ys.error()();
        return;
    }

    auto ftbl3 = cmd::get_function_Table_maker_St(filename, sampling_interval,
                                                  max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&oldfilename, &ftbl3, &xs, &ys, &prior, &likelihood, &thermo_algorithm, &myinit_seed,
             &newfilename, sampling_interval, max_number_of_values_per_iteration](auto model0ptr) {
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
                        t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0, t_adjust_beta,
                        t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance);

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
                    if (!maybe_modelLikelihood) {
                        std::cerr << maybe_modelLikelihood.error()();
                        return
                    };
                    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());
                    // using m2=typename
                    // decltype(modelLikelihood_v)::paseModelLikelihoodv;

                    std::visit(
                        [&oldfilename, &ftbl3, &tmi, &param1_prior, &ys,
                         &xs](auto& modelLikelihood) {
                            auto opt = thermo_fraction_evidence_continuation<true>(
                                oldfilename, ftbl3, std::move(tmi), param1_prior, modelLikelihood,
                                ys, xs);
                        },
                        modelLikelihood_v);
                }
            },
            model_v);
    }
}

}  // namespace cmd
}  // namespace macrodr

#endif  // CLI_THERMO_EVIDENCE_FRACTION_DTS_H
