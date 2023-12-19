#include "CLI_macro_dr.h"
#include "distributions.h"
#include "experiment.h"
#include "function_measure_verification_and_optimization.h"
#include "function_memoization.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "maybe_error.h"
#include "parameters_derivative.h"
#include "parameters_distribution.h"
#include "variables.h"
#include "variables_derivative.h"
// #include "multivariate_normal_distribution.h"
#include "allosteric_models.h"
#include "cuevi.h"
#include "micror_stochastic.h"
#include "parallel_tempering.h"
#include "qmodel.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include "models.h"

using namespace macrodr;




int main(int argc, char **argv) {
    
    std::vector<std::string> arguments(argc);
    for (auto i=0; i<argc; ++i)
        arguments[i]=argv[i];
    constexpr bool test_dynamice_command_line_interprester = false;
    
    if constexpr (test_dynamice_command_line_interprester) {
        auto cm = dcli::Compiler{};
        
        cm.push_function("load_experiment",
                         dcli::to_typed_function<std::string, double, double>(
                             &macrodr::load_experiment, "filename",
                             "frequency_of_sampling", "initial_ATP"));
        
        auto filename = "../macro_dr/run_script.txt";
        std::ifstream f(filename);
        
        std::string s;
        while (f) {
            std::string line;
            std::getline(f, line);
            s += line + "\n";
        }
        std::cout << "\ntest file \n" << s << "\n";
        auto p = dcli::extract_program(s);
        
        std::cerr << p;
        
        if (p) {
            auto c = dcli::compile_program(cm, p.value());
            if (c) {
                auto exec = c.value().run();
            }
        }
        
        if (p) {
            auto ss = p.value().str();
            std::cerr << ss;
        } else
            std::cerr << p.error()();
    }
    
    
    constexpr bool new_cuevi_by_max_iter = true;
    if (new_cuevi_by_max_iter) {
        auto myseed = 0ul;
        
        myseed = calc_seed(myseed);
        std::cerr << "myseed =" << myseed << "\n";
        
        mt_64i mt(myseed);
        
        auto Efilename_all = "../macro_dr/Moffatt_Hume_2007_ATP_time.txt";
        
        auto [recording_conditions_all, recording_all] =
            macrodr::load_recording(Efilename_all);
        
        auto experiment_all = Experiment(
            std::move(recording_conditions_all), Frequency_of_Sampling(50e3),
            initial_ATP_concentration(ATP_concentration(0.0)));
        
        auto Efilename_7 = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt";
        auto Efilename_7_const_dt = "../macro_dr/Moffatt_Hume_2007_ATP_time_7_constant_dt.txt";
        
        auto [recording_conditions_7, recording_7] =
            macrodr::load_recording(Efilename_7);
        
        auto [recording_conditions_7_dt, recording_7_dt] =
            macrodr::load_recording(Efilename_7_const_dt);
        
        auto experiment_7 =
            Experiment(std::move(recording_conditions_7), Frequency_of_Sampling(50e3),
                       initial_ATP_concentration(ATP_concentration(0.0)));
        
        auto experiment_7_dt =
            Experiment(std::move(recording_conditions_7_dt), Frequency_of_Sampling(50e3),
                       initial_ATP_concentration(ATP_concentration(0.0)));
        
        auto &experiment = experiment_7;
        auto &recording = recording_7;
        
        auto &model0 = model6_Eff_no_inactivation;
        auto &param1Names = model0.names();
        auto &param1 = model0.parameters();
        std::string ModelName = "model6_Eff_no_inactivation";
        using MyModel = Allost1;
        
        
        /**
     * @brief myseed defines the random number seed so all runs are identical
     * for debugging purposes
     */
        //   auto myseed = 9762841416869310605ul;
        //    auto myseed = 2555984001541913735ul;
        
        /**
     * @brief num_scouts_per_ensemble number of scouts per ensemble in the
     * affine ensemble mcmc model
     */
        std::size_t num_scouts_per_ensemble = 16;
        
        /**
     * @brief max_num_simultaneous_temperatures when the number of parallel
     * temepratures reaches this number, it stops growing, the same scouts
     * drifts on temperature
     *
     */
        std::size_t max_num_simultaneous_temperatures = 1e5;
        
        /**
     * @brief stops_at minimum value of beta greater than zero
     */
        double stops_at = 1e-15;
        
        /**
     * @brief beta where the slope changes     */
        double medium_beta = 1e-2;
        
        
        /**
     * @brief includes_zero considers also beta equal zero
     */
        bool includes_zero = true;
        
        
        /**
     * @brief randomly tries thermodynamic jumps
     */
        bool random_jumps = true;
        
        /**
     * @brief max_iter maximum number of iterations on each warming step
     */
        std::size_t max_iter_warming = 50;
        
        /**
     * @brief max_iter maximum number of iterations on the equilibrium step
     */
        std::size_t max_iter_equilibrium = 50000;
        
        /**
     * @brief path directory for the output
     */
        std::string path = "";
        
        /**
     * @brief min_fraction fraction of the prior parameter size used as the
     * minimal sample used for the cumulative sequence
     */
        double min_fraction = 4;
        
        /**
     * @brief checks_derivative_every_model_size number of steps before every
     * check of the derivative against the beta thermo parameter for stopping
     */
        std::size_t checks_derivative_every_model_size = 10;
        
        /**
     * @brief max_ratio maximimum tolerated ratio for the beta derivative method
     */
        double max_ratio = 8000e16;
        
        /**
     * @brief n_points_per_decade number of points per 10 times increment in
     * beta thermodynamic parameter
     */
        
        double n_points_per_decade = 1;
        /**
     * @brief n_points_per_decade_fraction number of points per 10 times
     * increment in the number of samples
     */
        double n_points_per_decade_fraction = 6;
        
        /**
     * @brief thermo_jumps_every factor that multiplied by the model size it
     * produces the number of steps skipped until the next thermo jump
     */
        std::size_t thermo_jumps_every = param1().size() / 4;
        
        double prior_error = 2;
        
        auto param1_prior = var::prior_around(param1, prior_error);
        
        // auto& param1_prior = prior_model00_7;
        
        /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
        
        auto modelLikelihood = make_Likelihood_Model<
            uses_adaptive_aproximation(true), uses_recursive_aproximation(true),
            uses_averaging_aproximation(2), uses_variance_aproximation(false),                                     uses_variance_correction_aproximation(false)>(
            model0, Simulation_n_sub_dt(1000ul));
        
        auto sim = Macro_DMR{}.sample(
            mt, model0, param1, experiment,
            Simulation_Parameters(Simulation_n_sub_dt(1000ul)),
            recording);
        
        if (sim) {
            std::vector<std::size_t> t_segments = {73, 33, 22, 22, 1, 1, 1, 1};
            auto number_of_traces = 7;
            auto number_of_segments = t_segments.size();
            t_segments.reserve(number_of_traces * t_segments.size());
            
            for (std::size_t i = 0; i + 1 < number_of_traces; ++i)
                std::copy_n(t_segments.begin(), number_of_segments,
                            std::back_inserter(t_segments));
            std::cerr << "t_segments\n" << t_segments;
            std::cerr << "cum t_segments\n" << var::cumsum(t_segments);
            
            std::vector<std::size_t> t_segments_7 = {73, 33, 22, 22};
            
            std::size_t t_min_number_of_samples = 20;
            
            /**
       * @brief cbc cumulative evidence algorithm, ends using convergence
       * criteria
       */
            std::size_t bisection_count = 2ul;
            std::string filename_bisection =
                ModelName + "_bisection_" + std::to_string(bisection_count) + "_" +
                std::to_string(myseed) + "_" + time_now();
            
            bool all_at_once = true;
            
            std::string all_at_once_str =
                all_at_once ? "_all_at_once_" : "_progressive_";
            
            std::string n_points_per_decade_str =
                "_" + std::to_string(n_points_per_decade) + "_";
            
            std::string filename = ModelName + "_new_cuevi_sim_eig_4800ch_only_7_" +
                                   n_points_per_decade_str + time_now() + "_" +
                                   // std::to_string(bisection_count) + "_" +
                                   std::to_string(myseed);
            
            auto &t_segments_used = t_segments_7;
            
            auto saving_itervals=Saving_intervals(Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),Save_Likelihood_every(num_scouts_per_ensemble),Save_Parameter_every(num_scouts_per_ensemble),Save_Predictions_every(num_scouts_per_ensemble*20)));
            
            auto cbc = new_cuevi_Model_by_iteration<MyModel>(
                path, filename, t_segments_used, t_min_number_of_samples,
                num_scouts_per_ensemble, max_num_simultaneous_temperatures,
                min_fraction, thermo_jumps_every, max_iter_warming,
                max_iter_equilibrium, max_ratio, n_points_per_decade,
                n_points_per_decade_fraction, medium_beta, stops_at,includes_zero, myseed, saving_itervals, random_jumps);
            
            // auto opt3 = evidence(std::move(cbc), param1_prior, modelLikelihood,
            //                      sim.value()(), experiment);
            auto ftbl3 = FuncMap(
                path + filename,
                Time_it(F(cuevi::step_stretch_cuevi_mcmc{}, cuevi::step_stretch_cuevi_mcmc{}),
                        num_scouts_per_ensemble / 2),
                Time_it(F(cuevi::thermo_cuevi_jump_mcmc{}, cuevi::thermo_cuevi_jump_mcmc{}),
                        num_scouts_per_ensemble / 2),
                Time_it(F(thermo_cuevi_randomized_jump_mcmc{},
                          thermo_cuevi_randomized_jump_mcmc{}),
                        num_scouts_per_ensemble / 2),
                var::Time_it(F(cuevi::step_stretch_cuevi_mcmc_per_walker{},
                               cuevi::step_stretch_cuevi_mcmc_per_walker{}),
                             num_scouts_per_ensemble / 2),
                var::Time_it(F(logLikelihood_f{},
                               [](auto &&...x) {
                                   return logLikelihood(
                                       std::forward<decltype(x)>(x)...);
                               }),
                             num_scouts_per_ensemble / 2),
                var::Time_it(F(MacroR<uses_recursive_aproximation(true),
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
                var::Time_it(F(MacroR<uses_recursive_aproximation(true),
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
                var::Time_it(F(MacroR<uses_recursive_aproximation(false),
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
                //            [](auto &&...x) {
                //              auto m = Macro_DMR{};
                //                auto bisection_order=16ul;
                //              return m.calc_Qdt_bisection(
                //                  std::forward<decltype(x)>(x)...,bisection_order);
                //            }),
                //     var::Memoiza_all_values<Maybe_error<Qdt>, ATP_step, double>{},
                //     num_scouts_per_ensemble / 2),
                var::Thread_Memoizer(
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
                var::Thread_Memoizer(
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
            auto lik = Macro_DMR{}
                           .log_Likelihood<uses_adaptive_aproximation(false),
                                           uses_recursive_aproximation(true),
                                           uses_averaging_aproximation(2),
                                           uses_variance_aproximation(false),
                                           uses_variance_correction_aproximation(false),
                                           return_predictions(true)>(
                               ftbl3.fork(var::I_thread(0)), model0, param1,
                               experiment, sim.value()());
            report(filename+"_lik.csv",lik.value(),sim.value(), experiment);
            if (true)
            {
                auto opt3 = cuevi::evidence(ftbl3, std::move(cbc), param1_prior, modelLikelihood,
                                            sim.value()(), experiment, cuevi::Init_seed(myseed));
                //auto opt4= cuevi::continue_evidence(filename, 2* maxiter);
            }
        }
    }
    
    
}
