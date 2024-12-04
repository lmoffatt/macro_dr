#ifndef CLI_LIKELIHOOD_H
#define CLI_LIKELIHOOD_H
#include "CLI_function_table.h"

#include "experiment.h"
#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "qmodel.h"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

namespace macrodr {
namespace cmd {

inline void calc_likelihood(std::string outfilename, std::string model,
                            parameters_value_type par,
                            likelihood_algo_type likelihood,
                            std::string recording, experiment_type experiment) {
    using namespace macrodr;
    auto ftbl3 = get_function_Table_maker_St(outfilename, 100, 100)();
    
    auto Maybe_model_v = get_model(model);
    if (Maybe_model_v) {
        auto model_v = std::move(Maybe_model_v.value());
        return std::visit(
            [outfilename, &par, &ftbl3, &experiment, &recording,
             &likelihood](auto model0ptr) {
                auto &model0 = *model0ptr;
                
                // auto [path, filename, average_the_ATP_evolution,
                //       num_scouts_per_ensemble, number_trials_until_give_up,
                //       thermo_jump_factor, max_iter_equilibrium,
                //       n_points_per_decade_beta_high, n_points_per_decade_beta_low,
                //       medium_beta, stops_at, includes_zero, random_jumps,
                //       sampling_interval,max_number_of_values_per_iteration] =
                //       std::move(algorithm);
                
                auto [adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction_approximation,
                      variance_correction, n_sub_dt] = likelihood;
                
                using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                
                auto Maybe_param1 = var::load_Parameters(
                    par.first, par.second, model0.model_name(), model0.names());
                Simulated_Recording<includes_N_state_evolution(false)> y;
                auto Maybe_y = load_simulation(recording, ",", y);
                if (!Maybe_param1.valid() || !Maybe_y.valid()) {
                    std::cerr << "---------ERROR_______________\n";
                    std::cerr << Maybe_param1.error()() << Maybe_y.error()();
                    std::cerr << "---------ERROR_______________\n";
                    
                } else {
                    
                    auto param1 = Maybe_param1.value().standard_parameter();
                    
                    if (recursive_approximation && averaging_approximation == 2) {
                        auto lik0 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(2),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(0)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        std::cerr << "<\nlog_Likelihood\n" << lik0 << "\n";
                        
                        auto lik1 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(2),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(1)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        //  std::cerr<<"\nlog_Likelihood 1\n"<<lik1<<"\n";
                        
                        auto lik2 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(2),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(2)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        if (lik2) {
                            save_Likelihood_Predictions(outfilename + ".csv", lik1.value(),
                                                        y, experiment);
                            save_Likelihood_Predictions(outfilename + "_cov.csv",
                                                        lik2.value(), y, experiment);
                        } else {
                            std::cerr << "---------ERROR_______________\n";
                            std::cerr << lik2.error()();
                            std::cerr << "---------ERROR_______________\n";
                        }
                    } else if (recursive_approximation &&
                               averaging_approximation == 1) {
                        auto lik0 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(0)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        std::cerr << "<\nlog_Likelihood\n" << lik0 << "\n";
                        
                        auto lik1 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(1)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        // std::cerr<<"\nlog_Likelihood 1\n"<<lik1<<"\n";
                        
                        auto lik2 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(2)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        if (lik2) {
                            save_Likelihood_Predictions(outfilename + ".csv", lik1.value(),
                                                        y, experiment);
                            save_Likelihood_Predictions(outfilename + "_cov.csv",
                                                        lik2.value(), y, experiment);
                        } else {
                            std::cerr << "---------ERROR_______________\n";
                            std::cerr << lik2.error()();
                            std::cerr << "---------ERROR_______________\n";
                        }
                    } else if (recursive_approximation &&
                               averaging_approximation == 0) {
                        auto lik0 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(0)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        std::cerr << "<\nlog_Likelihood\n" << lik0 << "\n";
                        
                        auto lik1 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(1)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        // std::cerr<<"\nlog_Likelihood 1\n"<<lik1<<"\n";
                        
                        auto lik2 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(true),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(2)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        if (lik2) {
                            save_Likelihood_Predictions(outfilename + ".csv", lik1.value(),
                                                        y, experiment);
                            save_Likelihood_Predictions(outfilename + "_cov.csv",
                                                        lik2.value(), y, experiment);
                        } else {
                            std::cerr << "---------ERROR_______________\n";
                            std::cerr << lik2.error()();
                            std::cerr << "---------ERROR_______________\n";
                        }
                    } else if (!recursive_approximation &&
                               averaging_approximation == 1) {
                        auto lik0 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(0)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        std::cerr << "<\nlog_Likelihood\n" << lik0 << "\n";
                        
                        auto lik1 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(1)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        // std::cerr<<"\nlog_Likelihood 1\n"<<lik1<<"\n";
                        
                        auto lik2 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(1),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(2)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        if (lik2) {
                            save_Likelihood_Predictions(outfilename + ".csv", lik1.value(),
                                                        y, experiment);
                            save_Likelihood_Predictions(outfilename + "_cov.csv",
                                                        lik2.value(), y, experiment);
                        } else {
                            std::cerr << "---------ERROR_______________\n";
                            std::cerr << lik2.error()();
                            std::cerr << "---------ERROR_______________\n";
                        }
                    } else if (!recursive_approximation &&
                               averaging_approximation == 0) {
                        
                        auto lik0 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(0)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        std::cerr << "<\nlog_Likelihood\n" << lik0 << "\n";
                        
                        auto lik1 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(1)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        
                        // std::cerr<<"\nlog_Likelihood 1\n"<<lik1<<"\n";
                        
                        auto lik2 = Macro_DMR{}
                                        .log_Likelihood<
                                            uses_adaptive_aproximation(true),
                                            uses_recursive_aproximation(false),
                                            uses_averaging_aproximation(0),
                                            uses_variance_aproximation(true),
                                            uses_variance_correction_aproximation(false),
                                            return_predictions(2)>(ftbl3, model0, param1,
                                                                   experiment,
                                                                   get<Recording>(y()));
                        if (lik2) {
                            save_Likelihood_Predictions(outfilename, lik1.value(), y,
                                                        experiment);
                            save_Likelihood_Predictions(outfilename, lik2.value(), y,
                                                        experiment);
                        } else {
                            std::cerr << "---------ERROR_______________\n";
                            std::cerr << lik2.error()();
                            std::cerr << "---------ERROR_______________\n";
                        }
                    }
                }
            },
            model_v);
    }
}

} // namespace cmd
} // namespace macrodr

#endif // CLI_LIKELIHOOD_H
