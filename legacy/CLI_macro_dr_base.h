#ifndef CLI_MACRO_DR_BASE_H
#define CLI_MACRO_DR_BASE_H

#include <macrodr/dsl/lexer_typed.h>

#include <cstddef>
#include <fstream>
#include <string>

#include "CLI_macro_dr.h"
#include "experiment.h"
#include "mcmc.h"

namespace macrodr::cmd {

inline auto get_random_id(std::string prefix) {
    return prefix + "_" + std::to_string(calc_seed(0ul));
}

inline std::size_t get_number(std::size_t number) {
    return number;
}

inline void write_text(std::string filename, std::string s) {
    std::ofstream f(filename);
    f << s;
}

static constexpr inline auto temp_script_file = "temp_script_file.txt";
inline void write_script(std::string program_name) {
    rename(temp_script_file, program_name.c_str());
}

inline auto load_Parameter_value(const std::string filename, std::string separator) {
    return std::pair(filename, separator);
}

inline auto load_Prior_value(std::string filename, std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_Recording_value(std::string filename, std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_fractioned_Recording_value(std::string filename, std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_fractioned_Experiment(std::string filename, std::string separator = ",") {
    return std::pair(filename, separator);
}
inline auto get_function_Table_maker_value_St(std::string filename) {
    return filename;
}

using recording_value_type =
    typename return_type<std::decay_t<decltype(&load_Recording_value)>>::type;
using prior_value_type = typename return_type<std::decay_t<decltype(&load_Prior_value)>>::type;
using parameters_value_type =
    typename return_type<std::decay_t<decltype(&load_Parameter_value)>>::type;

inline auto get_Prior(double prior_error, const std::string modelname) {
    return std::pair(prior_error, modelname);
}

inline auto set_Likelihood_algorithm(bool adaptive_aproximation, bool recursive_approximation,
                                     int averaging_approximation,
                                     bool variance_correction_approximation,
                                     bool variance_correction, std::size_t n_sub_dt) {
    return std::tuple(adaptive_aproximation, recursive_approximation, averaging_approximation,
                      variance_correction_approximation, variance_correction, n_sub_dt);
}

using likelihood_algo_type =
    typename return_type<std::decay_t<decltype(&set_Likelihood_algorithm)>>::type;

inline auto get_function_Table_maker_value(std::string filename,
                                           std::size_t num_scouts_per_ensemble) {
    return std::pair(filename, num_scouts_per_ensemble);
}

inline auto set_simulation_algorithm(bool includeN, std::size_t n_sub_dt) {
    return std::pair(includeN, n_sub_dt);
}

using simulation_algo_type =
    typename return_type<std::decay_t<decltype(&set_simulation_algorithm)>>::type;

inline dsl::Compiler make_utilities_compiler() {
    dsl::Compiler cm;
    cm.push_function("get_random_Id",
                     dsl::to_typed_function<std::string>(&get_random_id, "prefix"));
    cm.push_function("get_number", dsl::to_typed_function<std::size_t>(&get_number, "n"));
    return cm;
}

inline auto make_io_compiler() {
    auto cm = dsl::Compiler{};
    using namespace cmd;
    cm.push_function("write_text", dsl::to_typed_function<std::string, std::string>(
                                       &write_text, "filename", "text"));

    cm.push_function("write_script",
                     dsl::to_typed_function<std::string>(&write_script, "script_name"));

    cm.push_function("load_Parameter", dsl::to_typed_function<std::string, std::string>(
                                           &load_Parameter_value, "filename", "separator"));

    cm.push_function("load_Prior", dsl::to_typed_function<std::string, std::string>(
                                       &load_Prior_value, "filename", "separator"));

    return cm;
}

inline dsl::Compiler make_experiment_compiler() {
    dsl::Compiler cm;
    cm.push_function("get_Observations",
                     dsl::to_typed_function<std::string>(&get_Observations, "filename"));
    cm.push_function("idealize_Experiment",
                     dsl::to_typed_function<std::string, std::string, std::string>(
                         &idealize_Experiment, "experiment_filename", "sep", "idealized_filename"));
    cm.push_function("get_function_Table_maker",
                     dsl::to_typed_function<std::string, std::size_t>(
                         &get_function_Table_maker_value, "filename", "num_scouts_per_ensemble"));
    cm.push_function("get_Experiment",
                     dsl::to_typed_function<std::string, double, double>(
                         &get_Experiment, "filename", "frequency_of_sampling", "initial_ATP"));
    cm.push_function("get_Experiment_file",
                     dsl::to_typed_function<std::string, double, double>(
                         &get_Experiment_file, "filename", "frequency_of_sampling", "initial_ATP"));
    return cm;
}

inline dsl::Compiler make_model_compiler() {
    dsl::Compiler cm;
    cm.push_function("get_Prior", dsl::to_typed_function<double, std::string>(
                                      &get_Prior, "prior_error", "model"));

    cm.push_function(
        "set_Likelihood_algorithm",
        dsl::to_typed_function<bool, bool, int, bool, bool, std::size_t>(
            &set_Likelihood_algorithm, "adaptive_aproximation", "recursive_approximation",
            "averaging_approximation", "variance_correction_approximation",
            "variance_approximation", "n_sub_dt"));
    // cm.push_function("set_fraction_algorithm", dsl::to_typed_function<double, double, std::string>(
    //                                                &set_Fraction_algorithm, "min_fraction",
    //                                                "n_points_per_decade_fraction", "segments"));

    cm.push_function("simulation_algorithm",
                     dsl::to_typed_function<bool, std::size_t>(
                         &set_simulation_algorithm, "include_N_states", "number_of_substeps"));
    return cm;
}
#ifdef ZOMBIE
namespace zombie {
inline dsl::Compiler make_frac_compiler() {
    dsl::Compiler cm;
    cm.push_function(
        "fraction_experiment",
        dsl::to_typed_function<std::string, std::string, experiment_type, fraction_algo_type,
                               Maybe_error<std::size_t>, std::size_t>(
            &calc_experiment_fractions, "save_name", "recording", "experiment",
            "fraction_algorithm", "number_of_parameters", "init_seed"));
    // Aquí irían fraction_simulation, fraction_likelihood, evidence_fraction si están activos
    return cm;
}
}  // namespace zombie
#endif

} // namespace macrodr::cmd


#endif  // CLI_MACRO_DR_BASE_H
