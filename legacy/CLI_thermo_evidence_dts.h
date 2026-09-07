#ifndef CLI_THERMO_EVIDENCE_DTS_H
#define CLI_THERMO_EVIDENCE_DTS_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>

#include "CLI_function_table.h"
#include "experiment.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "qdtf_member.h"
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

// Drift-and-hold schedule for the beta ladder (ladder_schedule in
// parallel_tempering.h). Its own DSL object, like set_Acquisition_filter, so
// set_ThermoAlgorithm_dts and every script that calls it stay untouched: the
// thermo_evidence_dts overloads that take a ladder_schedule are its only
// consumers.
inline auto set_Ladder_schedule(std::size_t phase1_end, std::size_t drift, std::size_t hold,
                                std::size_t hold_burnin, std::size_t n_cycles,
                                double cycle_gain) {
    return std::tuple(phase1_end, drift, hold, hold_burnin, n_cycles, cycle_gain);
}
using ladder_schedule_type =
    typename return_type<std::decay_t<decltype(&set_Ladder_schedule)>>::type;

// No schedule: the ladder adapts all run long (the pre-2026-09-07 run).
inline ladder_schedule_type no_ladder_schedule() {
    return std::tuple(std::numeric_limits<std::size_t>::max(), std::size_t{0}, std::size_t{0},
                      std::size_t{0}, std::size_t{0}, 0.0);
}

inline ladder_schedule to_ladder_schedule(ladder_schedule_type const& t) {
    ladder_schedule s;
    std::tie(s.phase1_end, s.drift, s.hold, s.hold_burnin, s.n_cycles, s.cycle_gain) = t;
    return s;
}

// The schedule is read in two places: the tempering loop (when the ladder
// moves, when the statistics restart) and the evidence saver (when its
// telescopic windows may pool again after a move).
template <class Tmi>
void apply_ladder_schedule(Tmi& tmi, ladder_schedule_type const& t) {
    auto s = to_ladder_schedule(t);
    tmi.set_beta_schedule(s);
    std::get<save_Evidence>(tmi.reporter().m_m).ss_schedule = s;
}

// Core evidence run over an in-memory Experiment. Both entry points below
// (file-based, and inline as the eLife_2025 likelihood commands take it)
// funnel here; `filename` is the already-composed output prefix and `myseed`
// the already-resolved seed.
inline void run_thermo_evidence_dts(std::string filename, std::string model, std::string prior,
                                    likelihood_algo_type likelihood, std::string recording,
                                    const Experiment& experiment,
                                    thermo_algo_dts_type thermo_algorithm,
                                    ladder_schedule_type schedule, std::size_t sampling_interval,
                                    std::size_t max_number_of_values_per_iteration,
                                    std::size_t myseed) {
    using namespace macrodr;

    auto ftbl3 = cmd::get_function_Table_maker_St(filename, sampling_interval,
                                                  max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (!Maybe_model_v) {
        std::cerr << Maybe_model_v.error()();
    } else {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(

            [&filename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &schedule, &myseed, sampling_interval,
             max_number_of_values_per_iteration](auto model0ptr) {
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
                      variance_correction, taylor_variance_correction_approximation, n_sub_dt] =
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
                        apply_ladder_schedule(tmi, schedule);

                        auto maybe_modelLikelihood =
                            Likelihood_Model_regular<
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_taylor_variance_correction_aproximation, false>,
                                var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                                decltype(model0)>(
                                model0, Simulation_n_sub_dt(n_sub_dt),
                                uses_adaptive_aproximation_value(adaptive_aproximation),
                                uses_recursive_aproximation_value(recursive_approximation),
                                uses_averaging_aproximation_value(averaging_approximation),
                                uses_variance_aproximation_value(variance_correction),
                                uses_taylor_variance_correction_aproximation_value(
                                    taylor_variance_correction_approximation),
                                uses_family_aproximation_value(family_macro))
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

// Pre-schedule entry, same signature as before: the ladder adapts all run long.
inline void run_thermo_evidence_dts(std::string filename, std::string model, std::string prior,
                                    likelihood_algo_type likelihood, std::string recording,
                                    const Experiment& experiment,
                                    thermo_algo_dts_type thermo_algorithm,
                                    std::size_t sampling_interval,
                                    std::size_t max_number_of_values_per_iteration,
                                    std::size_t myseed) {
    run_thermo_evidence_dts(filename, model, prior, likelihood, recording, experiment,
                            thermo_algorithm, no_ladder_schedule(), sampling_interval,
                            max_number_of_values_per_iteration, myseed);
}

// File-based entry: loads the experiment from disk and writes the
// thermo_evidence_dts_<id>.txt restart file that the continuation commands
// reload. Unchanged behavior relative to the original command.
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
    auto experiment = get_Experiment(std::get<0>(experiment_file), std::get<1>(experiment_file),
                                     std::get<2>(experiment_file));
    run_thermo_evidence_dts(filename, model, prior, likelihood, recording, experiment,
                            thermo_algorithm, sampling_interval,
                            max_number_of_values_per_iteration, myseed);
}

// Inline-Experiment entry, same shape as the eLife_2025 likelihood commands
// (calc_dlikelihood_predictions and friends): the Experiment built with
// create_experiment in the .macroir flows straight into the run, no file
// round-trip, so the evidence lanes can sweep the experiment via
// dispatch-injected segments exactly like the figure lanes did. Stateless
// like those commands: it does NOT write the restart txt, so
// thermo_evidence_dts_continuation applies to file-based runs only.
inline void calc_thermo_evidence_dts(std::string id, std::string model, std::string prior,
                                     likelihood_algo_type likelihood, std::string recording,
                                     const Experiment& experiment,
                                     thermo_algo_dts_type thermo_algorithm,
                                     std::size_t sampling_interval,
                                     std::size_t max_number_of_values_per_iteration,
                                     std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);
    run_thermo_evidence_dts(filename, model, prior, likelihood, recording, experiment,
                            thermo_algorithm, sampling_interval,
                            max_number_of_values_per_iteration, myseed);
}

// Inline-Experiment entry with a drift-and-hold ladder schedule
// (set_Ladder_schedule). Same DSL name; the extra named argument selects it.
inline void calc_thermo_evidence_dts(std::string id, std::string model, std::string prior,
                                     likelihood_algo_type likelihood, std::string recording,
                                     const Experiment& experiment,
                                     thermo_algo_dts_type thermo_algorithm,
                                     ladder_schedule_type ladder_schedule,
                                     std::size_t sampling_interval,
                                     std::size_t max_number_of_values_per_iteration,
                                     std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);
    run_thermo_evidence_dts(filename, model, prior, likelihood, recording, experiment,
                            thermo_algorithm, ladder_schedule, sampling_interval,
                            max_number_of_values_per_iteration, myseed);
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
                      variance_correction, taylor_variance_correction_approximation, n_sub_dt] =
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
                                    bool, uses_taylor_variance_correction_aproximation, false>,
                                var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                                decltype(model0)>(
                                model0, Simulation_n_sub_dt(n_sub_dt),
                                uses_adaptive_aproximation_value(adaptive_aproximation),
                                uses_recursive_aproximation_value(recursive_approximation),
                                uses_averaging_aproximation_value(averaging_approximation),
                                uses_variance_aproximation_value(variance_correction),
                                uses_taylor_variance_correction_aproximation_value(
                                    taylor_variance_correction_approximation),
                                uses_family_aproximation_value(family_macro))
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
                      variance_correction, taylor_variance_correction_approximation, n_sub_dt] =
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
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_taylor_variance_correction_aproximation, false>,
                                var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                                decltype(model0)>(
                                model0, Simulation_n_sub_dt(n_sub_dt),
                                uses_adaptive_aproximation_value(adaptive_aproximation),
                                uses_recursive_aproximation_value(recursive_approximation),
                                uses_averaging_aproximation_value(averaging_approximation),
                                uses_variance_aproximation_value(variance_correction),
                                uses_taylor_variance_correction_aproximation_value(
                                    taylor_variance_correction_approximation),
                                uses_family_aproximation_value(family_macro))
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
                      variance_correction, taylor_variance_correction_approximation, n_sub_dt] =
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
                                var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                                var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                                var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                                var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                                var::constexpr_Var_domain<
                                    bool, uses_taylor_variance_correction_aproximation, false>,
                                var::constexpr_Var_domain<int, uses_family_aproximation, family_macro>,
                                decltype(model0)>(
                                model0, Simulation_n_sub_dt(n_sub_dt),
                                uses_adaptive_aproximation_value(adaptive_aproximation),
                                uses_recursive_aproximation_value(recursive_approximation),
                                uses_averaging_aproximation_value(averaging_approximation),
                                uses_variance_aproximation_value(variance_correction),
                                uses_taylor_variance_correction_aproximation_value(
                                    taylor_variance_correction_approximation),
                                uses_family_aproximation_value(family_macro))
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
// ── Bessel (qdtf) member: evidence with the acquisition filter explicit ─────
// The filter is rig metadata (pole count, −3 dB cutoff in Hz), never fitted:
// it enters the DSL as its own object so the SAME likelihood_algorithm can be
// passed to a box arm and a Bessel arm of a paired fit. n_poles = 0 selects
// the qdtf box configuration (the regression anchor against the av=2 member).
inline auto set_Acquisition_filter(std::size_t n_poles, double cutoff_hz) {
    return std::pair(n_poles, cutoff_hz);
}
using acquisition_filter_type =
    typename return_type<std::decay_t<decltype(&set_Acquisition_filter)>>::type;

// Reporter tuple for the qdtf member: save_Score and save_Predictions are
// EXCLUDED because the member is plain-likelihood only (no dlogLikelihood /
// logLikelihoodPredictions overloads exist for it; see the adapter note in
// qdtf_member.h). Everything else is the dts tuple unchanged: save_likelihood
// and save_Evidence read stored walker state, save_RateParameter only needs
// lik.m, so the telescopic/trapezoid evidence columns come out identical in
// form to the box runs.
inline auto new_thermo_Model_by_max_iter_dts_qdtf(
    std::string path, std::string filename, std::size_t num_scouts_per_ensemble,
    std::size_t thermo_jumps_every, std::size_t max_iter_equilibrium, std::size_t beta_size,
    Saving_intervals sint, std::size_t initseed, std::size_t t_adapt_beta_every,
    std::string t_adapt_beta_equalizer, std::string t_adapt_beta_constroler,
    std::string t_adapt_beta_variance, double t_adapt_beta_nu, double t_adapt_beta_t0,
    double t_adapt_beta_threshold, bool t_adjust_beta, double t_acceptance_upper_limit,
    double t_acceptance_lower_limit, double t_desired_acceptance) {
    return new_thermodynamic_integration(
        thermo_less_than_max_iteration(max_iter_equilibrium),
        save_mcmc<var::Parameters_transformed, save_Iter,
                  save_likelihood<var::Parameters_transformed>,
                  save_Parameter<var::Parameters_transformed>,
                  save_RateParameter<var::Parameters_transformed>, save_Evidence>(
            path, filename, std::pair(1ul, 1ul), get<Save_Likelihood_every>(sint())(),
            get<Save_Parameter_every>(sint())(), get<Save_RateParameter_every>(sint())(),
            get<Save_Evidence_every>(sint())()),
        num_scouts_per_ensemble, thermo_jumps_every, beta_size, initseed, t_adapt_beta_every,
        t_adapt_beta_equalizer, t_adapt_beta_constroler, t_adapt_beta_variance, t_adapt_beta_nu,
        t_adapt_beta_t0, t_adapt_beta_threshold, t_adjust_beta, t_acceptance_upper_limit,
        t_acceptance_lower_limit, t_desired_acceptance);
}

// Core evidence run for the qdtf member over an in-memory Experiment. Same
// skeleton as run_thermo_evidence_dts, but the likelihood is the concrete
// Qdtf_Likelihood_Model (no flag-domain variant: the member has no member
// flags, only the filter). Of likelihood_algorithm only n_sub_dt is used
// (for the concept-satisfying simulate overload); the av/variance/taylor
// flags do not apply to this member.
inline void run_thermo_evidence_dts_qdtf(std::string filename, std::string model,
                                         std::string prior, likelihood_algo_type likelihood,
                                         std::string recording, const Experiment& experiment,
                                         thermo_algo_dts_type thermo_algorithm,
                                         acquisition_filter_type acquisition_filter,
                                         ladder_schedule_type schedule,
                                         std::size_t sampling_interval,
                                         std::size_t max_number_of_values_per_iteration,
                                         std::size_t myseed) {
    using namespace macrodr;

    auto ftbl3 = cmd::get_function_Table_maker_St(filename, sampling_interval,
                                                  max_number_of_values_per_iteration)();

    auto Maybe_model_v = get_model(model);

    if (!Maybe_model_v) {
        std::cerr << Maybe_model_v.error()();
    } else {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&filename, &ftbl3, &experiment, &recording, &prior, &likelihood, &thermo_algorithm,
             &acquisition_filter, &schedule, &myseed, sampling_interval,
             max_number_of_values_per_iteration](auto model0ptr) {
                std::string sep = ",";
                auto& model0 = *model0ptr;

                auto [num_scouts_per_ensemble, number_trials_until_give_up, thermo_jumps_every,
                      max_iter_equilibrium, beta_size, save_every_param_size_factor,
                      t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                      t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, t_adjust_beta,
                      t_acceptance_upper_limit, t_acceptance_lower_limit, t_desired_acceptance] =
                    std::move(thermo_algorithm);

                auto n_sub_dt = std::get<5>(likelihood);

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

                        auto tmi = new_thermo_Model_by_max_iter_dts_qdtf(
                            "", filename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, saving_intervals, myseed,
                            t_adapt_beta_every, t_adapt_beta_equalizer, t_adapt_beta_controler,
                            t_adapt_beta_variance, t_adapt_beta_nu, t_adapt_beta_t0, 0.0,
                            t_adjust_beta, t_acceptance_upper_limit, t_acceptance_lower_limit,
                            t_desired_acceptance);
                        apply_ladder_schedule(tmi, schedule);

                        auto lik = Qdtf_Likelihood_Model<std::decay_t<decltype(model0)>>{
                            model0, Simulation_n_sub_dt(n_sub_dt),
                            static_cast<int>(acquisition_filter.first),
                            acquisition_filter.second};

                        auto opt = thermo_evidence<true>(ftbl3, std::move(tmi), param1_prior, lik,
                                                         y, experiment);
                    }
                }
            },
            model_v);
    }
}

// Pre-schedule qdtf entry, same signature as before.
inline void run_thermo_evidence_dts_qdtf(std::string filename, std::string model,
                                         std::string prior, likelihood_algo_type likelihood,
                                         std::string recording, const Experiment& experiment,
                                         thermo_algo_dts_type thermo_algorithm,
                                         acquisition_filter_type acquisition_filter,
                                         std::size_t sampling_interval,
                                         std::size_t max_number_of_values_per_iteration,
                                         std::size_t myseed) {
    run_thermo_evidence_dts_qdtf(filename, model, prior, likelihood, recording, experiment,
                                 thermo_algorithm, acquisition_filter, no_ladder_schedule(),
                                 sampling_interval, max_number_of_values_per_iteration, myseed);
}

// Inline-Experiment entry for the qdtf member, stateless like the inline box
// entry above (no restart txt; continuation applies to file-based runs only).
inline void calc_thermo_evidence_dts_qdtf(std::string id, std::string model, std::string prior,
                                          likelihood_algo_type likelihood, std::string recording,
                                          const Experiment& experiment,
                                          thermo_algo_dts_type thermo_algorithm,
                                          acquisition_filter_type acquisition_filter,
                                          std::size_t sampling_interval,
                                          std::size_t max_number_of_values_per_iteration,
                                          std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);
    run_thermo_evidence_dts_qdtf(filename, model, prior, likelihood, recording, experiment,
                                 thermo_algorithm, acquisition_filter, sampling_interval,
                                 max_number_of_values_per_iteration, myseed);
}

// Same, with a drift-and-hold ladder schedule (set_Ladder_schedule).
inline void calc_thermo_evidence_dts_qdtf(std::string id, std::string model, std::string prior,
                                          likelihood_algo_type likelihood, std::string recording,
                                          const Experiment& experiment,
                                          thermo_algo_dts_type thermo_algorithm,
                                          acquisition_filter_type acquisition_filter,
                                          ladder_schedule_type ladder_schedule,
                                          std::size_t sampling_interval,
                                          std::size_t max_number_of_values_per_iteration,
                                          std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename = id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);
    run_thermo_evidence_dts_qdtf(filename, model, prior, likelihood, recording, experiment,
                                 thermo_algorithm, acquisition_filter, ladder_schedule,
                                 sampling_interval, max_number_of_values_per_iteration, myseed);
}

// Bessel-filtered simulation entry (the plan's M1 truth generator), mirroring
// the legacy 7-arg simulate in filename and output conventions plus the
// acquisition filter argument. n_poles = 0 delegates to the unfiltered
// substep sampler bit-identically (the M1 gate). Substeps only; the
// include_N_states channel is not carried (plain Recording output).
inline Maybe_error<std::string> run_qdtf_simulation(
    std::string filename_prefix, std::string recording_file, const Experiment& experiment,
    std::size_t myseed, const std::string& modelName, parameters_value_type parameter_files,
    simulation_algo_type sim_algo_type, acquisition_filter_type acquisition_filter) {
    if (sim_algo_type.algorithm != "substeps")
        return error_message(
            "simulate with acquisition_filter supports number_of_substeps only, got algorithm \"" +
            sim_algo_type.algorithm + "\"");
    if (sim_algo_type.include_N_states)
        return error_message("simulate with acquisition_filter does not carry N-state evolution");
    Recording recording;
    auto Maybe_y = load_Recording_Data(recording_file, ",", recording);
    if (!Maybe_y)
        return Maybe_y.error();
    auto Maybe_model_v = get_model(modelName);
    if (!Maybe_model_v)
        return Maybe_model_v.error();
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&](auto model0ptr) -> Maybe_error<std::string> {
            auto& model0 = *model0ptr;
            auto seed = calc_seed(myseed);
            mt_64i mt(seed);
            auto Maybe_parameter_values = var::load_Parameters(
                parameter_files.first, parameter_files.second, model0.model_name(), model0.names());
            if (!Maybe_parameter_values)
                return Maybe_parameter_values.error();
            auto param1 = Maybe_parameter_values.value().standard_parameter();
            std::string filename = filename_prefix + "_" + model0.model_name() + "_" + time_now() +
                                   "_" + std::to_string(seed);
            auto sim = sample_bessel(mt, model0, param1, experiment,
                                     sim_algo_type.number_of_substeps,
                                     static_cast<int>(acquisition_filter.first),
                                     acquisition_filter.second, recording);
            if (!sim)
                return sim.error();
            save_Recording(filename + "_simulation.csv", ",", get<Recording>(sim.value()()));
            return filename + "_simulation.csv";
        },
        model_v);
}

inline dsl::Compiler<dsl::Lexer> make_dts_compiler() {
    dsl::Compiler<dsl::Lexer> cm;
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
            static_cast<void (*)(std::string, std::string, std::string, likelihood_algo_type,
                                 std::string, experiment_file_type, thermo_algo_dts_type,
                                 std::size_t, std::size_t, std::size_t)>(&calc_thermo_evidence_dts),
            "idname", "model", "prior", "likelihood_algorithm", "data", "experiment",
            "thermo_algorithm", "sampling_interval", "max_number_of_values_per_iteration",
            "init_seed"));
    cm.push_function(
        "thermo_evidence_dts",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, const Experiment&, thermo_algo_dts_type, std::size_t,
                               std::size_t, std::size_t>(
            static_cast<void (*)(std::string, std::string, std::string, likelihood_algo_type,
                                 std::string, const Experiment&, thermo_algo_dts_type, std::size_t,
                                 std::size_t, std::size_t)>(&calc_thermo_evidence_dts),
            "idname", "model", "prior", "likelihood_algorithm", "data", "experiment",
            "thermo_algorithm", "sampling_interval", "max_number_of_values_per_iteration",
            "init_seed"));
    cm.push_function(
        "thermo_evidence_dts_2",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, experiment_file_type, thermo_algo_dts_type, std::size_t,
                               std::size_t, std::size_t>(
            &calc_thermo_evidence_dts_2, "idname", "model", "prior", "likelihood_algorithm", "data",
            "experiment", "thermo_algorithm", "sampling_interval",
            "max_number_of_values_per_iteration", "init_seed"));
    cm.push_function("set_Acquisition_filter",
                     dsl::to_typed_function<std::size_t, double>(&set_Acquisition_filter,
                                                                 "n_poles", "cutoff_hz"));
    cm.push_function("set_Ladder_schedule",
                     dsl::to_typed_function<std::size_t, std::size_t, std::size_t, std::size_t,
                                            std::size_t, double>(
                         &set_Ladder_schedule, "phase1_end", "drift", "hold", "hold_burnin",
                         "n_cycles", "cycle_gain"));
    // Same DSL name, one extra named argument (ladder_schedule) selects the
    // drift-and-hold run; inline-Experiment form only (stateless).
    cm.push_function(
        "thermo_evidence_dts",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, const Experiment&, thermo_algo_dts_type,
                               ladder_schedule_type, std::size_t, std::size_t, std::size_t>(
            static_cast<void (*)(std::string, std::string, std::string, likelihood_algo_type,
                                 std::string, const Experiment&, thermo_algo_dts_type,
                                 ladder_schedule_type, std::size_t, std::size_t, std::size_t)>(
                &calc_thermo_evidence_dts),
            "idname", "model", "prior", "likelihood_algorithm", "data", "experiment",
            "thermo_algorithm", "ladder_schedule", "sampling_interval",
            "max_number_of_values_per_iteration", "init_seed"));
    // Same DSL name as the legacy 7-arg simulate; the extra named argument
    // (acquisition_filter) selects the Bessel-filtered truth generator.
    cm.push_function(
        "simulate",
        dsl::to_typed_function<std::string, std::string, const Experiment&, std::size_t,
                               const std::string&, parameters_value_type, simulation_algo_type,
                               acquisition_filter_type>(
            &run_qdtf_simulation, "filename_prefix", "recording", "experiment", "init_seed",
            "modelName", "parameter_values", "simulation_algorithm", "acquisition_filter"));
    // Same DSL name, one extra named argument (acquisition_filter) selects
    // the Bessel/qdtf member; inline-Experiment form only (stateless).
    cm.push_function(
        "thermo_evidence_dts",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, const Experiment&, thermo_algo_dts_type,
                               acquisition_filter_type, std::size_t, std::size_t, std::size_t>(
            static_cast<void (*)(std::string, std::string, std::string, likelihood_algo_type,
                                 std::string, const Experiment&, thermo_algo_dts_type,
                                 acquisition_filter_type, std::size_t, std::size_t, std::size_t)>(
                &calc_thermo_evidence_dts_qdtf),
            "idname", "model", "prior", "likelihood_algorithm", "data", "experiment",
            "thermo_algorithm", "acquisition_filter", "sampling_interval",
            "max_number_of_values_per_iteration", "init_seed"));
    // qdtf member with a drift-and-hold ladder schedule.
    cm.push_function(
        "thermo_evidence_dts",
        dsl::to_typed_function<std::string, std::string, std::string, likelihood_algo_type,
                               std::string, const Experiment&, thermo_algo_dts_type,
                               acquisition_filter_type, ladder_schedule_type, std::size_t,
                               std::size_t, std::size_t>(
            static_cast<void (*)(std::string, std::string, std::string, likelihood_algo_type,
                                 std::string, const Experiment&, thermo_algo_dts_type,
                                 acquisition_filter_type, ladder_schedule_type, std::size_t,
                                 std::size_t, std::size_t)>(&calc_thermo_evidence_dts_qdtf),
            "idname", "model", "prior", "likelihood_algorithm", "data", "experiment",
            "thermo_algorithm", "acquisition_filter", "ladder_schedule", "sampling_interval",
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
