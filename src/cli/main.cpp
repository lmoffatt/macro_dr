#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr.h"
//#include "CLI_thermo_evidence.h"
#include "CLI_thermo_evidence_dts.h"

//#include "CLI_thermo_evidence_fraction_dts.h"
//#include "CLI_thermo_levenberg_evidence.h"
#include <macrodr/cmd/init_commands.h>
#include <macrodr/dsl/lexer_typed.h>

#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "experiment.h"
#include "lapack_headers.h"
#include "maybe_error.h"
#include "models_MoffattHume_allosteric.h"
#include "qmodel.h"
using namespace macrodr;

Maybe_error<std::string> read_from_input_files(std::vector<std::string> const args) {
    std::string s;
    for (std::size_t i = 1; i < args.size(); ++i) {
        auto filename = args[i];
        std::ifstream f(filename);
        if (!f)
            return error_message(filename + " does not exist or cannot be opened");
        else {
            while (f) {
                std::string line;
                std::getline(f, line);
                s += line + "\n";
            }
        }
    }
    return s;
}

Maybe_error<std::string> append_files_content(std::string&& s, std::string const& filename) {
    std::ifstream f(filename);
    if (!f)
        return error_message(filename + " does not exist or cannot be opened");
    else {
        while (f) {
            std::string line;
            std::getline(f, line);
            s += line + "\n";
        }
    }
    return std::move(s);
}

Maybe_error<std::string> append_command(std::string&& s, std::string const& line) {
    if (line[line.size() - 1] != '"')
        return error_message(
            "does not have a "
            " in" +
            line);
    s += line.substr(1, line.size() - 2) + "\n";
    return std::move(s);
}

Maybe_error<std::string> read_from_input_files_or_commands(std::vector<std::string> const args) {
    std::string s;
    for (std::size_t i = 1; i < args.size(); ++i) {
        auto filename_or_command = args[i];
        Maybe_error<std::string> Maybe_s;
        if (filename_or_command.starts_with("--")) {
            s += filename_or_command.substr(2) + "\n";
        } else {
            auto Maybe_s = append_files_content(std::move(s), filename_or_command);
            if (!Maybe_s)
                return Maybe_s.error();
            else
                s = std::move(Maybe_s.value());
        }
    }
    return s;
}

auto get_compiler() {
    auto cm = dsl::Compiler{};
    using namespace cmd;
    cm.push_function("get_random_Id",
                     dsl::to_typed_function<std::string>(&get_random_id, "prefix"));

    cm.push_function("get_number", dsl::to_typed_function<std::size_t>(&get_number, "n"));

    cm.push_function("write_text", dsl::to_typed_function<std::string, std::string>(
                                       &write_text, "filename", "text"));

    cm.push_function("write_script",
                     dsl::to_typed_function<std::string>(&write_script, "script_name"));

    cm.push_function("load_Parameter", dsl::to_typed_function<std::string, std::string>(
                                           &load_Parameter_value, "filename", "separator"));

    cm.push_function("load_Prior", dsl::to_typed_function<std::string, std::string>(
                                       &load_Prior_value, "filename", "separator"));

    cm.push_function("get_Observations",
                     dsl::to_typed_function<std::string>(&get_Observations, "filename"));

    cm.push_function("get_Experiment_file",
                     dsl::to_typed_function<std::string, double, double>(
                         &get_Experiment_file, "filename", "frequency_of_sampling", "initial_ATP"));

    cm.push_function("get_num_parameters",
                     dsl::to_typed_function<std::string>(&get_num_parameters, "model"));

    /**
*
* auto get_Prior(
double prior_error ,
const std::string& modelname*/

    cm.push_function("get_Prior", dsl::to_typed_function<double, std::string>(
                                      &get_Prior, "prior_error", "model"));

    /**

auto get_Likelihood(const std::string& model, bool adaptive_aproximation, bool
recursive_approximation, int averaging_approximation, bool
variance_correction_approximation,  std::size_t n_sub_dt)

*/

    cm.push_function(
        "set_Likelihood_algorithm",
        dsl::to_typed_function<bool, bool, int, bool, bool, std::size_t>(
            &set_Likelihood_algorithm, "adaptive_aproximation", "recursive_approximation",
            "averaging_approximation", "variance_correction_approximation",
            "variance_approximation", "n_sub_dt"));

    /**
inline auto set_Fraction_algorithm(
  double min_fraction ,
  double n_points_per_decade_fraction ,
  std::string segments)
* */

    // cm.push_function("set_fraction_algorithm", dsl::to_typed_function<double, double, std::string>(
    //                                                &set_Fraction_algorithm, "min_fraction",
    //                                                "n_points_per_decade_fraction", "segments"));

    /**
   * calc_experiment_fractions(std::string save_name,std::string recording,
   * experiment_type experiment, fraction_algo_type fraction_algo, std::string
   * model,std::size_t i_seed )
   * */

    // cm.push_function(
    //     "fraction_experiment",
    //     dsl::to_typed_function<std::string, std::string, experiment_type, fraction_algo_type,
    //                             Maybe_error<std::size_t>, std::size_t>(
    //         &calc_experiment_fractions, "save_name", "recording", "experiment",
    //         "fraction_algorithm", "number_of_parameters", "init_seed"));

    /**
   * inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_simulation_fractions(std::string save_name,std::string simulation,
experiment_type experiment,                                  fraction_algo_type
fraction_algo, std::string model,std::size_t i_seed )

   * */

    // if constexpr(false){
    // cm.push_function(
    //     "fraction_simulation",
    //     dsl::to_typed_function<std::string, std::string, cmd::experiment_type,
    //                             fraction_algo_type, Maybe_error<std::size_t>,
    //                             std::size_t>(
    //         &cmd::calc_simulation_fractions, "save_name", "simulation",
    //         "experiment", "fraction_algorithm", "number_of_parameters",
    //         "init_seed"));
    // }

    /**
   * calc_fraction_likelihood(const std::string file_name,
                                          std::string model,
                                   parameters_value_type par,
                                   fractioned_simulation_type
   Maybe_frac_simulation, likelihood_algo_type likelihood_algo,
                                   tablefun_value_type ft)
   * */

    // if constexpr(false){
    // cm.push_function(
    //     "fraction_likelihood",
    //     dsl::to_typed_function<std::string, std::string, parameters_value_type,
    //                             fractioned_simulation_type, likelihood_algo_type,
    //                             tablefun_value_type>(
    //         &calc_fraction_likelihood, "file_name", "model", "parameter",
    //         "fractioned_simulation", "likelihood_algorithm", "function_table"));
    // }

    /**
   * inline void calc_fraction_evidence(std::string model,
                        prior_value_type prior,
                        likelihood_algo_type likelihood,
                        fractioned_simulation_type Maybe_frac_experiment,
                        fraction_algo_type fraction_algo,
                        cuevi_algo_type cuevi_algorithm,
                        tablefun_value_type ft, std::size_t myseed)
   * */

    // if constexpr(false){
    // cm.push_function(
    //     "evidence_fraction",
    //     dsl::to_typed_function<std::string, prior_value_type,
    //                             likelihood_algo_type, fractioned_simulation_type,
    //                             cuevi_algo_type, tablefun_value_type,
    //                             std::size_t>(
    //         &calc_fraction_evidence, "model", "prior", "likelihood_algorithm",
    //         "fractioned_experiment", "cuevi_algorithm", "function_table",
    //         "init_seed"));
    // }

    /**    inline auto set_CueviAlgorithm(
  std::size_t num_scouts_per_ensemble = 16,
  std::size_t number_trials_until_give_up = 1e5,
  double stops_at = 1e-15,
  double medium_beta = 1e-2,
  bool includes_zero = true,
  bool random_jumps = true,
  std::size_t max_iter_equilibrium = 50000,
  std::string path = "",
  double n_points_per_decade = 1,
  bool average_the_ATP_evolution = true,
  std::string filename = "haha",
  std::size_t thermo_jumps_every = 10)   */
    /*
    cm.push_function(
        "set_CueviAlgorithm",
        dsl::to_typed_function<std::size_t, std::size_t, double, double, bool, bool, std::size_t,
                                std::string, double, double, std::size_t, std::string, std::size_t,
                                std::size_t, std::size_t>(
            &set_CueviAlgorithm, "num_scouts_per_ensemble", "number_trials_until_give_up",
            "stops_at", "medium_beta", "includes_zero", "random_jumps", "max_iter_equilibrium",
            "path", "n_points_per_decade_beta_high", "n_points_per_decade_beta_low",
            "average_the_ATP_evolution", "filename", "thermo_jumps_every", "sampling_interval",
            "max_number_of_values_per_iteration"));

    /**
   * set_ThermoAlgorithm(std::size_t num_scouts_per_ensemble,
                    std::size_t number_trials_until_give_up, double stops_at,
                    double beta_upper_value, double beta_medium_value,
                    bool includes_zero, std::size_t max_iter_equilibrium,
                    std::string path, std::size_t beta_size,
                    std::size_t beta_upper_size, std::size_t beta_medium_size,
                    std::string filename, std::size_t thermo_jumps_every,
                    std::size_t max_num_simultaneous_temperatures) {
 */

    //   if constexpr(false){
    //   cm.push_function(
    //       "set_ThermoAlgorithm",
    //       dsl::to_typed_function<std::size_t ,
    //                               std::size_t ,
    //                               double ,
    //                               double , double ,
    //                               bool , std::size_t ,
    //                               std::size_t ,
    //                               std::size_t , std::size_t ,
    //                               std::size_t ,
    //                               std::size_t >(
    //           &set_ThermoAlgorithm, "num_scouts_per_ensemble",
    //           "number_trials_until_give_up", "stops_at", "beta_upper_value",
    //           "beta_medium_value",
    //           "includes_zero",
    //           "max_iter_equilibrium", "beta_size", "beta_upper_size",
    //           "beta_medium_size",  "thermo_jumps_every",
    //           "save_every_param_size_factor"));
    // }

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

    // if constexpr (false){
    // cm.push_function(
    //     "set_ThermoAlgorithm_fraction_dts",
    //     dsl::to_typed_function<    std::size_t ,
    //                             std::size_t ,
    //                             std::size_t , std::size_t ,
    //                             std::size_t , std::size_t ,std::size_t ,
    //                             std::string,
    //                             std::string,
    //                             std::string,
    //                             double ,double, double,bool   ,
    //                             double ,
    //                             double ,
    //                             double  >(
    //         &set_ThermoAlgorithm_fraction_dts, "num_scouts_per_ensemble",
    //         "number_trials_until_give_up",
    //         "max_iter_equilibrium", "beta_size",  "thermo_jumps_every",
    //         "save_every_param_size_factor","adapt_beta_every","adapt_beta_equalizer","adapt_beta_controler","adapt_beta_variance","adapt_beta_nu","adapt_beta_t0","adapt_beta_threshold","adjust_beta","acceptance_upper_limit","acceptance_lower_limit","desired_acceptance"));

    // }

    /**
calc_thermo_evidence(std::string id,
                                              std::string model,
                                              std::string prior,
                                 likelihood_algo_type likelihood,
                                 std::string recording,
                                 experiment_type experiment,
                                 thermo_algo_type thermo_algorithm,
                                 std::string ft_filename,
                                 std::size_t sampling_interval,max_number_of_values_per_iteration,
                                 std::size_t myseed)
   * */

    //   if constexpr(false){
    //   cm.push_function(
    //       "thermo_evidence",
    //       dsl::to_typed_function<std::string,
    //                               std::string, std::string, likelihood_algo_type, std::string,
    //                               experiment_file_type, thermo_algo_type, std::size_t, std::size_t, std::size_t>(
    //           &calc_thermo_evidence, "idname","model", "prior", "likelihood_algorithm",
    //           "data", "experiment", "thermo_algorithm", "sampling_interval","max_number_of_values_per_iteration",
    //           "init_seed"));
    // }

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

    // if constexpr(false){

    // cm.push_function(
    //     "thermo_fraction_evidence_dts",
    //     dsl::to_typed_function<std::string,
    //                             std::string, std::string, likelihood_algo_type, fractioned_experiment_type, thermo_algo_fraction_dts_type, std::size_t,std::size_t, std::size_t>(
    //         &calc_thermo_evidence_fraction_dts, "idname","model", "prior", "likelihood_algorithm",
    //         "fractional_experiment",  "thermo_algorithm", "sampling_interval","max_number_of_values_per_iteration",
    //         "init_seed"));

    // }

    /**
   *  calc_thermo_evidence_continuation(std::string id, std::size_t ith)
   * */

    // if constexpr(false){
    // cm.push_function(
    //     "thermo_evidence_continuation",
    //     dsl::to_typed_function<std::string, std::size_t, std::size_t>(
    //         &calc_thermo_evidence_continuation, "idname","continuation_number",
    //         "init_seed"));

    // }

    cm.push_function(
        "thermo_evidence_dts_continuation",
        dsl::to_typed_function<std::string, std::size_t, std::size_t>(
            &calc_thermo_evidence_dts_continuation, "idname", "continuation_number", "init_seed"));

    cm.push_function("thermo_evidence_dts_continuation_2",
                     dsl::to_typed_function<std::string, std::size_t, std::size_t>(
                         &calc_thermo_evidence_dts_continuation_2, "idname", "continuation_number",
                         "init_seed"));

    /**
 * set_ThermoLevenAlgorithm(
    std::size_t num_scouts_per_ensemble,
    std::size_t number_trials_until_give_up, double stops_at,
    double beta_upper_value, double beta_medium_value, bool includes_zero,
    std::size_t max_iter_equilibrium, std::size_t beta_size,
    std::size_t beta_upper_size, std::size_t beta_medium_size,
    std::size_t thermo_jumps_every, std::size_t save_every_param_size_factor)
 * 
 * */

    // if constexpr(false){

    // cm.push_function(
    //     "set_ThermoLevenAlgorithm",
    //     dsl::to_typed_function<std::size_t ,
    //                             std::size_t ,
    //                             double ,
    //                             double , double ,
    //                             bool , std::size_t ,
    //                             std::size_t ,
    //                             std::size_t , std::size_t ,
    //                             std::size_t ,
    //                             std::string,
    //                             std::size_t,
    //                             std::size_t >(
    //         &set_ThermoLevenAlgorithm, "num_scouts_per_ensemble",
    //         "number_trials_until_give_up", "stops_at", "beta_upper_value",
    //         "beta_medium_value",
    //         "includes_zero",
    //         "max_iter_equilibrium", "beta_size", "beta_upper_size",
    //         "beta_medium_size",  "n_lambdas","lambda_adaptive_algorithm","thermo_jumps_every",
    //         "save_every_param_size_factor"));

    // }

    /**
   * inline void calc_thermo_levenberg_evidence(std::string id, std::string model,
                                 std::string prior,
                                 likelihood_algo_type likelihood,
                                 std::string recording,
                                 experiment_file_type experiment_file,
                                 thermo_leven_algo_type thermo_algorithm,
                                 std::size_t sampling_interval,max_number_of_values_per_iteration, std::size_t myseed) 
   * */

    // if constexpr (false){
    // cm.push_function(
    //     "thermo_levenberg_evidence",
    //     dsl::to_typed_function<std::string,
    //                             std::string, std::string, likelihood_algo_type, std::string,
    //                             experiment_file_type, thermo_leven_algo_type, std::size_t, std::size_t,std::size_t, double>(
    //         &calc_thermo_levenberg_evidence, "idname","model", "prior", "likelihood_algorithm",
    //         "data", "experiment", "thermo_levenberg_algorithm", "sampling_interval","max_number_of_values_per_iteration",
    //         "init_seed", "delta_par"));

    // }

    /**
   *  auto set_simulation_algorithm(bool includeN, std::size_t n_sub_dt)
   *
   * */

    cm.push_function("simulation_algorithm",
                     dsl::to_typed_function<bool, std::size_t>(
                         &set_simulation_algorithm, "include_N_states", "number_of_substeps"));

    /**
   * bool run_simulation(std::string filename,recording_type recording,
   experiment_type experiment, std::size_t myseed, std::string
   modelName,Matrix<double> parameter_values,bool includeN,std::size_t n_sub_dt)
   {

   * */

    cm.push_function(
        "simulate",
        dsl::to_typed_function<std::string, recording_type, experiment_type, std::size_t,
                               std::string, parameters_value_type, simulation_algo_type>(
            &run_simulation, "output", "recording", "experiment", "init_seed", "modelName",
            "parameter_values", "simulation_algorithm"));

    /**
   *inline void calc_likelihood(std::string outfilename,
                          std::string model,
                          parameters_value_type par,
                          likelihood_algo_type likelihood,
                          recording_value_type recording,
                          experiment_type experiment, cuevi_algo_type algorithm,
                          tablefun_value_type ft)
   */

    cm.push_function("likelihood",
                     dsl::to_typed_function<std::string, std::string, parameters_value_type,
                                            likelihood_algo_type, std::string, experiment_type>(
                         &calc_likelihood, "output", "model", "parameter_values",
                         "likelihood_algorithm", "recording", "experiment"));

    /**
   * std::string outfilename,
                                std::string model,
                                parameters_value_type par,
                                likelihood_algo_type likelihood,
                                std::string recording,
                                experiment_type experiment   * */
    // if constexpr(false){
    // cm.push_function("evidence",
    //                  dsl::to_typed_function<
    //                      std::string, prior_value_type, likelihood_algo_type,
    //                      std::string, experiment_type, fraction_algo_type,
    //                      cuevi_algo_type, tablefun_value_type, std::size_t>(
    //                      &calc_evidence, "model", "prior", "likelihood_algorithm",
    //                      "data", "experiment", "fraction_algorithm",
    //                      "cuevi_algorithm", "function_table", "init_seed"));
    // }

    return cm;
}
inline dsl::Compiler make_compiler() {
    dsl::Compiler cm;
    cm.merge(macrodr::cmd::make_utilities_compiler());
    cm.merge(macrodr::cmd::make_io_compiler());
    cm.merge(macrodr::cmd::make_experiment_compiler());
    cm.merge(macrodr::cmd::make_model_compiler());
    cm.merge(macrodr::cmd::make_simulation_compiler());
    // cm.merge(macrodr::cmd::make_frac_compiler());
    cm.merge(macrodr::cmd::make_dts_compiler());
    // cm.merge(macrodr::cmd::zombie::make_cuevi_compiler());
    //  cm.merge(macrodr::cmd::make_levenberg_compiler());
    return cm;
}

int main(int argc, char** argv) {
    print_model_Priors(2.0);
    std::vector<std::string> arguments(argc);
    for (auto i = 0; i < argc; ++i) arguments[i] = argv[i];

    auto Maybe_script = read_from_input_files_or_commands(arguments);
    if (!Maybe_script)
        std::cerr << "Error: \n" << Maybe_script.error()();
    else {
        auto cm = make_compiler();
        auto s = std::move(Maybe_script.value());
        std::cout << "\n read files " << arguments << "\n" << s << "\n";
        /// Horrible hack to force the script to write itself at the start
        cmd::write_text(cmd::temp_script_file, s);

        auto p = dsl::extract_program(s);

        std::cerr << p;
        
        dsl::Environment<dsl::Lexer,dsl::Compiler> env(cm);
        if (p) {
            auto c = dsl::compile_program(env, p.value());
            if (!c) {
                std::cerr << "\n --------------Error--------\n"
                          << c.error()() << "\n --------------Error--------\n";
            } else {
                auto exec = c.value().run(env);
            }
        }

        if (p) {
            auto ss = p.value().str();
            std::cerr << ss;
        } else
            std::cerr << p.error()();
    }
}
