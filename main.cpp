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
#include "CLI_macrodr_commands.h"
#include "allosteric_models.h"
#include "cuevi.h"
#include "micror_stochastic.h"
#include "models.h"
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
using namespace macrodr;

int main(int argc, char **argv) {

  std::vector<std::string> arguments(argc);
  for (auto i = 0; i < argc; ++i)
    arguments[i] = argv[i];
  std::string filename = arguments[1];
  std::ifstream f(filename);

  std::string s;
  while (f) {
    std::string line;
    std::getline(f, line);
    s += line + "\n";
  }

  auto cm = dcli::Compiler{};

  cm.push_function("load_experiment",
                   dcli::to_typed_function<std::string, double, double>(
                       &macrodr::load_experiment, "filename",
                       "frequency_of_sampling", "initial_ATP"));

  cm.push_function(
      "get_function_Table_maker",
      dcli::to_typed_function<std::string, std::size_t>(
          &get_function_Table_maker, "filename", "num_scouts_per_ensemble"));

  /*
auto get_Experiment(
  std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
  double frequency_of_sampling = 50e3, double initial_ATP = 0)
*/
  cm.push_function(
      "get_Experiment",
      dcli::to_typed_function<std::string, double, double>(
          &get_Experiment, "filename", "frequency_of_sampling", "initial_ATP"));

  /**
   *
   * auto get_Prior(
    double prior_error ,
    const std::string& modelname*/

  cm.push_function("get_Prior", dcli::to_typed_function<double, std::string>(
                                    &get_Prior, "prior_error", "modelname"));

  /**

  auto get_Likelihood(const std::string& model, bool adaptive_aproximation, bool
  recursive_approximation, int averaging_approximation, bool
  variance_correction_approximation,  std::size_t n_sub_dt)

  */

  cm.push_function(
      "get_Likelihood",
      dcli::to_typed_function<std::string, bool, bool, int, bool, bool, std::size_t>(
          &get_Likelihood, "model", "adaptive_aproximation",
          "recursive_approximation", "averaging_approximation",
          "variance_correction_approximation", "variance_approximation","n_sub_dt"));

  /**    auto get_CueviAlgorithm(
          std::size_t num_scouts_per_ensemble = 16,
          std::size_t number_trials_until_give_up = 1e5, double stops_at =
     1e-15, double medium_beta = 1e-2, bool includes_zero = true, bool
     random_jumps = true, std::size_t max_iter_equilibrium = 50000, std::string
     path = "", double min_fraction = 4, double n_points_per_decade = 1, double
     n_points_per_decade_fraction = 6,

          std::size_t t_min_number_of_samples = 20, std::string filename =
     "haha", std::size_t thermo_jumps_every = 10);
     */
  cm.push_function(
      "get_CueviAlgorithm",
      dcli::to_typed_function<std::size_t, std::size_t, double, double, bool,
                              bool, std::size_t, std::string, double, double,
                              double, std::size_t, std::string, std::size_t>(
          &get_CueviAlgorithm, "num_scouts_per_ensemble",
          "number_trials_until_give_up", "stops_at", "medium_beta",
          "includes_zero", "random_jumps", "max_iter_equilibrium", "path",
          "min_fraction", "n_points_per_decade", "n_points_per_decade_fraction",
          "t_min_number_of_samples", "filename", "thermo_jumps_every"));
  
  cm.push_function(
      "evidence",
      dcli::to_typed_function<prior_type  , likelihood_type , recording_type ,experiment_type  ,
                              algo_type,tablefun_value_type, std::size_t >(
          &calc_evidence,  "prior", "likelihoodModel","data" , "experiment","algorithm" , "function_table" , "init_seed"));
  
  
  
  // (prior_type prior, likelihood_type likelihood,
  //  recording_type recording, experiment_type experiment,
  //  algo_type algorithm, tablefun_value_type ft,
  //  std::size_t myseed)
  
  std::cout << "\n read file " << filename << "\n" << s << "\n";
  auto p = dcli::extract_program(s);

  std::cerr << p;

  if (p) {
    auto c = dcli::compile_program(cm, p.value());
    if (!c) {
        std::cerr <<"\n --------------Error--------\n"<< c.error()()<<"\n --------------Error--------\n";
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
