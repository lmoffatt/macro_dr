#ifndef CLI_THERMO_LEVENBERG_EVIDENCE_H
#define CLI_THERMO_LEVENBERG_EVIDENCE_H
#include "CLI_function_table.h"

#include "experiment.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_levenberg_tempering.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

namespace macrodr {
namespace cmd {

inline auto set_ThermoLevenAlgorithm(
    std::size_t num_scouts_per_ensemble,
    std::size_t number_trials_until_give_up, double stops_at,
    double beta_upper_value, double beta_medium_value, bool includes_zero,
    std::size_t max_iter_equilibrium, std::size_t beta_size,
    std::size_t beta_upper_size, std::size_t beta_medium_size,
    std::size_t n_lambdas,
    std::size_t thermo_jumps_every, std::size_t save_every_param_size_factor) {
    using namespace macrodr;
    
    return std::tuple(num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jumps_every, max_iter_equilibrium, beta_size,
                      beta_upper_size, beta_medium_size, beta_upper_value,
                      beta_medium_value, n_lambdas,stops_at, includes_zero,
                      save_every_param_size_factor);
}

using thermo_leven_algo_type =
    typename return_type<std::decay_t<decltype(&set_ThermoLevenAlgorithm)>>::type;



inline void calc_thermo_levenberg_evidence(std::string id, std::string model,
                                 std::string prior,
                                 likelihood_algo_type likelihood,
                                 std::string recording,
                                 experiment_file_type experiment_file,
                                 thermo_leven_algo_type thermo_algorithm,
                                 std::size_t save_every, std::size_t myseed) {
    myseed = calc_seed(myseed);
    std::string filename =
        id + "_" + model + "_" + time_now() + "_" + std::to_string(myseed);
    
    if (true) {
        std::ofstream f("thermo_evidence_" + id + ".txt");
        save_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, save_every, myseed);
    }
    using namespace macrodr;
    auto experiment= get_Experiment(std::get<0>(experiment_file),std::get<1>(experiment_file),std::get<2>(experiment_file));
    
    auto ftbl3 = cmd::get_function_Table_maker_St(filename, save_every)();
    
    auto Maybe_model_v = get_model(model);
    
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&filename, &ftbl3, &experiment, &recording, &prior, &likelihood,
             &thermo_algorithm, &myseed](auto model0ptr) {
                std::string sep = ",";
                auto &model0 = *model0ptr;
                mt_64i mt(myseed);
                
                auto [num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium, beta_size,
                      beta_upper_size, beta_medium_size, beta_upper_value,
                      beta_medium_value, n_lambdas,stops_at, includes_zero,
                      save_every_param_size_factor] = std::move(thermo_algorithm);
                
                auto [adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction,
                      variance_correction_approximation, n_sub_dt] = likelihood;
                
                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                
                std::string ModelName = model0.model_name();
                
                auto Maybe_param1_prior = var::load_Prior(
                    prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());
                    std::size_t thermo_jumps_every =
                        param1_prior.size() *  thermo_jump_factor;
                    
                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        
                        auto saving_intervals = Saving_Levenberg_intervals(Vector_Space(
                            Save_Evidence_every(save_every_param_size_factor *
                                                param1_prior.size()),
                            Save_Likelihood_every(save_every_param_size_factor *
                                                  param1_prior.size()),
                            Save_Parameter_every(save_every_param_size_factor *
                                                 param1_prior.size()),
                            save_Levenberg_Lambdas_every(save_every_param_size_factor *
                                                 param1_prior.size()*5),
                            Save_Predictions_every(save_every_param_size_factor *
                                                   param1_prior.size() * 50)));
                        
                        auto tmi = thermo_levenberg_Model_by_max_iter(
                            "", filename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, beta_upper_size,
                            beta_medium_size, beta_upper_value, beta_medium_value,
                            n_lambdas,
                            
                            stops_at, includes_zero, saving_intervals, myseed);
                        
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
                            [&ftbl3, &tmi, &param1_prior, &y,
                             &experiment](auto &modelLikelihood) {
                                auto opt =
                                    thermo_levenberg_evidence(ftbl3, std::move(tmi), param1_prior,
                                                    modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}

inline void calc_thermo_levenberg_evidence_continuation(std::string id, std::size_t ith, std::size_t myinit_seed) {
    
    std::string model;
    std::string prior;
    likelihood_algo_type likelihood;
    std::string recording;
    experiment_file_type experiment_file;
    thermo_leven_algo_type thermo_algorithm;
    std::size_t save_every;
    std::size_t myseed;
    std::string filename;
    if (true) {
        std::ifstream f("thermo_evidence_" + id + ".txt");
        if (!f)
        {
            std::cerr<<"Error!!\n thermo_evidence_" + id + ".txt"+ " not found";
            return;
        }
        load_vars(f, filename, model, prior, likelihood, recording, experiment_file,
                  thermo_algorithm, save_every, myseed);
    }
    auto experiment= get_Experiment(std::get<0>(experiment_file),std::get<1>(experiment_file),std::get<2>(experiment_file));
    
    std::string oldfilename=filename;
    if (ith>1)
        oldfilename=filename+'_'+std::to_string(ith-1);
    
    std::string newfilename=filename+'_'+std::to_string(ith);
    using namespace macrodr;
    
    
    auto ftbl3 = get_function_Table_maker_St(newfilename, save_every)();
    
    auto Maybe_model_v = get_model(model);
    
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [&oldfilename, &ftbl3, &experiment, &recording, &prior, &likelihood,
             &thermo_algorithm, &myinit_seed, &newfilename](auto model0ptr) {
                std::string sep = ",";
                auto &model0 = *model0ptr;
                auto myseed = calc_seed(myinit_seed);
                mt_64i mt(myseed);
                
                auto [num_scouts_per_ensemble, number_trials_until_give_up,
                      thermo_jump_factor, max_iter_equilibrium, beta_size,
                      beta_upper_size, beta_medium_size, beta_upper_value,
                      beta_medium_value, n_lambda,stops_at, includes_zero,
                      save_every_param_size_factor] = std::move(thermo_algorithm);
                
                auto [adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction,
                      variance_correction_approximation, n_sub_dt] = likelihood;
                
                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                
                std::string ModelName = model0.model_name();
                
                auto Maybe_param1_prior = var::load_Prior(
                    prior, sep, model0.model_name(), model0.names());
                if (!Maybe_param1_prior) {
                    std::cerr << "\n-------------errror------------\n"
                              << Maybe_param1_prior.error()();
                } else {
                    auto param1_prior = std::move(Maybe_param1_prior.value());
                    std::size_t thermo_jumps_every =
                        param1_prior.size() * thermo_jump_factor;
                    
                    Recording y;
                    auto Maybe_y = load_Recording_Data(recording, ",", y);
                    if (Maybe_y) {
                        
                        auto saving_intervals = Saving_intervals(Vector_Space(
                            Save_Evidence_every(save_every_param_size_factor *
                                                param1_prior.size()),
                            Save_Likelihood_every(save_every_param_size_factor *
                                                  param1_prior.size()),
                            Save_Parameter_every(save_every_param_size_factor *
                                                 param1_prior.size()),
                            Save_Predictions_every(save_every_param_size_factor *
                                                   param1_prior.size() * 50)));
                        
                        auto tmi = new_thermo_Model_by_max_iter(
                            "", newfilename, num_scouts_per_ensemble, thermo_jumps_every,
                            max_iter_equilibrium, beta_size, beta_upper_size,
                            beta_medium_size, beta_upper_value, beta_medium_value,
                            
                            stops_at, includes_zero, saving_intervals, myseed);
                        
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
                            [&oldfilename,&ftbl3, &tmi, &param1_prior, &y,
                             &experiment](auto &modelLikelihood) {
                                auto opt =
                                    thermo_evidence_continuation(oldfilename,ftbl3, std::move(tmi), param1_prior,
                                                                 modelLikelihood, y, experiment);
                            },
                            modelLikelihood_v);
                    }
                }
            },
            model_v);
    }
}
}
}


#endif // CLI_THERMO_LEVENBERG_EVIDENCE_H
