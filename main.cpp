#include "CLI_function_table.h"
#include "CLI_macro_dr.h"
#include "CLI_thermo_evidence.h"
#include "experiment.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "maybe_error.h"
#include "models_MoffattHume_linear.h"

#include "qmodel.h"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace macrodr;

Maybe_error<std::string>
read_from_input_files(std::vector<std::string> const args) {
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

Maybe_error<std::string> append_files_content(std::string &&s,
                                              std::string const &filename) {
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

Maybe_error<std::string> append_command(std::string &&s,
                                        std::string const &line) {
  if (line[line.size() - 1] != '"')
    return error_message("does not have a "
                         " in" +
                         line);
  s += line.substr(1, line.size() - 2) + "\n";
  return std::move(s);
}

Maybe_error<std::string>
read_from_input_files_or_commands(std::vector<std::string> const args) {
  std::string s;
  for (std::size_t i = 1; i < args.size(); ++i) {
    auto filename_or_command = args[i];
    Maybe_error<std::string> Maybe_s;
    if (filename_or_command.find('"') != filename_or_command.npos) {
      s += filename_or_command + "\n";
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
  auto cm = dcli::Compiler{};
  using namespace cmd;
  cm.push_function("get_random_Id", dcli::to_typed_function<std::string>(
                                        &get_random_id, "prefix"));

  cm.push_function("write_text",
                   dcli::to_typed_function<std::string, std::string>(
                       &write_text, "filename", "text"));

  cm.push_function("write_script", dcli::to_typed_function<std::string>(
                                       &write_script, "script_name"));

  cm.push_function("load_Parameter",
                   dcli::to_typed_function<std::string, std::string>(
                       &load_Parameter_value, "filename", "separator"));

  cm.push_function("load_Prior",
                   dcli::to_typed_function<std::string, std::string>(
                       &load_Prior_value, "filename", "separator"));

  cm.push_function("get_Observations", dcli::to_typed_function<std::string>(
                                           &get_Observations, "filename"));

  cm.push_function("load_experiment",
                   dcli::to_typed_function<std::string, double, double>(
                       &macrodr::load_experiment, "filename",
                       "frequency_of_sampling", "initial_ATP"));

  cm.push_function("get_function_Table_maker",
                   dcli::to_typed_function<std::string, std::size_t>(
                       &get_function_Table_maker_value, "filename",
                       "num_scouts_per_ensemble"));

  /*
auto get_Experiment(
std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
double frequency_of_sampling = 50e3, double initial_ATP = 0)
*/
  cm.push_function(
      "get_Experiment",
      dcli::to_typed_function<std::string, double, double>(
          &get_Experiment, "filename", "frequency_of_sampling", "initial_ATP"));

  cm.push_function("get_Observations", dcli::to_typed_function<std::string>(
                                           &get_Observations, "filename"));

  cm.push_function("get_num_parameters", dcli::to_typed_function<std::string>(
                                             &get_num_parameters, "model"));

  /**
*
* auto get_Prior(
double prior_error ,
const std::string& modelname*/

  cm.push_function("get_Prior", dcli::to_typed_function<double, std::string>(
                                    &get_Prior, "prior_error", "model"));

  /**

auto get_Likelihood(const std::string& model, bool adaptive_aproximation, bool
recursive_approximation, int averaging_approximation, bool
variance_correction_approximation,  std::size_t n_sub_dt)

*/

  cm.push_function(
      "set_Likelihood_algorithm",
      dcli::to_typed_function<bool, bool, int, bool, bool, std::size_t>(
          &set_Likelihood_algorithm, "adaptive_aproximation",
          "recursive_approximation", "averaging_approximation",
          "variance_correction_approximation", "variance_approximation",
          "n_sub_dt"));

  /**
inline auto set_Fraction_algorithm(
  double min_fraction ,
  double n_points_per_decade_fraction ,
  std::string segments)
* */

  cm.push_function("set_fraction_algorithm",
                   dcli::to_typed_function<double, double, std::string>(
                       &set_Fraction_algorithm, "min_fraction",
                       "n_points_per_decade_fraction", "segments"));

  /**
   * calc_experiment_fractions(std::string save_name,std::string recording,
   * experiment_type experiment, fraction_algo_type fraction_algo, std::string
   * model,std::size_t i_seed )
   * */

  cm.push_function(
      "fraction_experiment",
      dcli::to_typed_function<std::string, std::string, experiment_type,
                              fraction_algo_type, Maybe_error<std::size_t>,
                              std::size_t>(
          &calc_experiment_fractions, "save_name", "recording", "experiment",
          "fraction_algorithm", "number_of_parameters", "init_seed"));

  /**
   * inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_simulation_fractions(std::string save_name,std::string simulation,
experiment_type experiment,                                  fraction_algo_type
fraction_algo, std::string model,std::size_t i_seed )

   * */

  cm.push_function(
      "fraction_simulation",
      dcli::to_typed_function<std::string, std::string, cmd::experiment_type,
                              fraction_algo_type, Maybe_error<std::size_t>,
                              std::size_t>(
          &cmd::calc_simulation_fractions, "save_name", "simulation",
          "experiment", "fraction_algorithm", "number_of_parameters",
          "init_seed"));

  /**
   * calc_fraction_likelihood(const std::string file_name,
                                          std::string model,
                                   parameters_value_type par,
                                   fractioned_simulation_type
   Maybe_frac_simulation, likelihood_algo_type likelihood_algo,
                                   tablefun_value_type ft)
   * */

  cm.push_function(
      "fraction_likelihood",
      dcli::to_typed_function<std::string, std::string, parameters_value_type,
                              fractioned_simulation_type, likelihood_algo_type,
                              tablefun_value_type>(
          &calc_fraction_likelihood, "file_name", "model", "parameter",
          "fractioned_simulation", "likelihood_algorithm", "function_table"));

  /**
   * inline void calc_fraction_evidence(std::string model,
                        prior_value_type prior,
                        likelihood_algo_type likelihood,
                        fractioned_simulation_type Maybe_frac_experiment,
                        fraction_algo_type fraction_algo,
                        cuevi_algo_type cuevi_algorithm,
                        tablefun_value_type ft, std::size_t myseed)
   * */

  cm.push_function(
      "evidence_fraction",
      dcli::to_typed_function<std::string, prior_value_type,
                              likelihood_algo_type, fractioned_simulation_type,
                              cuevi_algo_type, tablefun_value_type,
                              std::size_t>(
          &calc_fraction_evidence, "model", "prior", "likelihood_algorithm",
          "fractioned_experiment", "cuevi_algorithm", "function_table",
          "init_seed"));

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
  std::size_t t_min_number_of_samples = 20,
  std::string filename = "haha",
  std::size_t thermo_jumps_every = 10)   */

  cm.push_function(
      "set_CueviAlgorithm",
      dcli::to_typed_function<std::size_t, std::size_t, double, double, bool,
                              bool, std::size_t, std::string, double, double,
                              std::size_t, std::string, std::size_t>(
          &set_CueviAlgorithm, "num_scouts_per_ensemble",
          "number_trials_until_give_up", "stops_at", "medium_beta",
          "includes_zero", "random_jumps", "max_iter_equilibrium", "path",
          "n_points_per_decade_beta_high", "n_points_per_decade_beta_low",
          "t_min_number_of_samples", "filename", "thermo_jumps_every"));

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

  cm.push_function(
      "set_ThermoAlgorithm",
      dcli::to_typed_function<std::size_t ,
                              std::size_t ,
                              double ,
                              double , double ,
                              bool , std::size_t ,
                              std::string , std::size_t ,
                              std::size_t , std::size_t ,
                              std::string , std::size_t ,
                              std::size_t >(
          &set_ThermoAlgorithm, "num_scouts_per_ensemble",
          "number_trials_until_give_up", "stops_at", "beta_upper_value",
          "beta_medium_value",
          "includes_zero",
          "max_iter_equilibrium", "path", "beta_size", "beta_upper_size",
          "beta_medium_size", "filename", "thermo_jumps_every",
          "max_num_simultaneous_temperatures"));

  /**
   *calc_thermo_evidence(std::string model,
                               prior_value_type prior,
                               likelihood_algo_type likelihood,
                               std::string recording,
                               experiment_type experiment,
                               thermo_algo_type thermo_algorithm,
                               tablefun_value_type ft,
                               std::size_t myseed) {

   * */

  cm.push_function(
      "thermo_evidence",
      dcli::to_typed_function<
          std::string, prior_value_type, likelihood_algo_type, std::string,
          experiment_type, thermo_algo_type, tablefun_value_type, std::size_t>(
          &calc_thermo_evidence, "model", "prior", "likelihood_algorithm",
          "data", "experiment", "thermo_algorithm", "function_table",
          "init_seed"));

  /**
   *  auto set_simulation_algorithm(bool includeN, std::size_t n_sub_dt)
   *
   * */

  cm.push_function(
      "simulation_algorithm",
      dcli::to_typed_function<bool, std::size_t>(
          &set_simulation_algorithm, "include_N_states", "number_of_substeps"));

  /**
   * bool run_simulation(std::string filename,recording_type recording,
   experiment_type experiment, std::size_t myseed, std::string
   modelName,Matrix<double> parameter_values,bool includeN,std::size_t n_sub_dt)
   {

   * */

  cm.push_function(
      "simulate",
      dcli::to_typed_function<std::string, recording_type, experiment_type,
                              std::size_t, std::string, parameters_value_type,
                              simulation_algo_type>(
          &run_simulation, "output", "recording", "experiment", "init_seed",
          "modelName", "parameter_values", "simulation_algorithm"));

  /**
   *inline void calc_likelihood(std::string outfilename,
                          std::string model,
                          parameters_value_type par,
                          likelihood_algo_type likelihood,
                          recording_value_type recording,
                          experiment_type experiment, cuevi_algo_type algorithm,
                          tablefun_value_type ft)
   */

  cm.push_function(
      "likelihood",
      dcli::to_typed_function<std::string, std::string, parameters_value_type,
                              likelihood_algo_type, recording_value_type,
                              experiment_type, cuevi_algo_type,
                              tablefun_value_type>(
          &calc_likelihood, "output", "model", "parameter_values",
          "likelihood_algorithm", "recording", "experiment", "cuevi_algorithm",
          "function_table"));

  /**
   * std::string model,
                        prior_value_type prior,
                        likelihood_algo_type likelihood,
                        std::string recording, experiment_type experiment,
                        fraction_algo_type fraction_algo, cuevi_algo_type
   cuevi_algorithm, tablefun_value_type ft, std::size_t myseed     * */

  cm.push_function("evidence",
                   dcli::to_typed_function<
                       std::string, prior_value_type, likelihood_algo_type,
                       std::string, experiment_type, fraction_algo_type,
                       cuevi_algo_type, tablefun_value_type, std::size_t>(
                       &calc_evidence, "model", "prior", "likelihood_algorithm",
                       "data", "experiment", "fraction_algorithm",
                       "cuevi_algorithm", "function_table", "init_seed"));
  return cm;
}

int main(int argc, char **argv) {

  print_model_Priors(2.0);
  std::vector<std::string> arguments(argc);
  for (auto i = 0; i < argc; ++i)
    arguments[i] = argv[i];

  auto Maybe_script = read_from_input_files_or_commands(arguments);
  if (!Maybe_script)
    std::cerr << "Error: \n" << Maybe_script.error()();
  else {
    auto cm = get_compiler();

    auto s = std::move(Maybe_script.value());
    std::cout << "\n read files " << arguments << "\n" << s << "\n";
    /// Horrible hack to force the script to write itself at the start
    cmd::write_text(cmd::temp_script_file, s);

    auto p = dcli::extract_program(s);

    std::cerr << p;

    if (p) {
      auto c = dcli::compile_program(cm, p.value());
      if (!c) {
        std::cerr << "\n --------------Error--------\n"
                  << c.error()() << "\n --------------Error--------\n";
      } else {
        auto exec = c.value().run();
      }
    }

    if (p) {
      auto ss = p.value().str();
      std::cerr << ss;
    } else
      std::cerr << p.error()();
  }
}
