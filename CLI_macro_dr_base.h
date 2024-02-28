#ifndef CLI_MACRO_DR_BASE_H
#define CLI_MACRO_DR_BASE_H


#include "mcmc.h"
#include <cstddef>
#include <string>
#include <fstream>
namespace macrodr {
namespace cmd {

inline auto get_random_id(std::string prefix) {
    
    return prefix + "_" + std::to_string(calc_seed(0ul));
}

inline std::size_t get_number(std::size_t number){return number;}


inline void write_text(std::string filename, std::string s) {
    std::ofstream f(filename);
    f << s;
}

static constexpr inline auto temp_script_file = "temp_script_file.txt";
inline void write_script(std::string program_name) {
    rename(temp_script_file, program_name.c_str());
}

inline auto load_Parameter_value(const std::string filename,
                                 std::string separator) {
    return std::pair(filename, separator);
}

inline auto load_Prior_value(std::string filename,
                             std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_Recording_value(std::string filename,
                                 std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_fractioned_Recording_value(std::string filename,
                                            std::string separator = ",") {
    return std::pair(filename, separator);
}

inline auto load_fractioned_Experiment(std::string filename,
                                       std::string separator = ",") {
    return std::pair(filename, separator);
}
inline auto get_function_Table_maker_value_St(std::string filename) {
    return filename;
}

using recording_value_type =
    typename return_type<std::decay_t<decltype(&load_Recording_value)>>::type;
using prior_value_type =
    typename return_type<std::decay_t<decltype(&load_Prior_value)>>::type;
using parameters_value_type =
    typename return_type<std::decay_t<decltype(&load_Parameter_value)>>::type;


inline auto get_Prior(double prior_error, const std::string modelname) {
    
    return std::pair(prior_error, modelname);
}

inline auto set_Likelihood_algorithm(bool adaptive_aproximation,
                                     bool recursive_approximation,
                                     int averaging_approximation,
                                     bool variance_correction_approximation,
                                     bool variance_correction,
                                     std::size_t n_sub_dt) {
    
    return std::tuple(adaptive_aproximation, recursive_approximation,
                      averaging_approximation, variance_correction_approximation,
                      variance_correction, n_sub_dt);
}


using likelihood_algo_type = typename return_type<
    std::decay_t<decltype(&set_Likelihood_algorithm)>>::type;


inline auto
get_function_Table_maker_value(std::string filename,
                               std::size_t num_scouts_per_ensemble) {
    return std::pair(filename, num_scouts_per_ensemble);
}

inline auto set_simulation_algorithm(bool includeN, std::size_t n_sub_dt) {
    return std::pair(includeN, n_sub_dt);
}

using simulation_algo_type = typename return_type<
    std::decay_t<decltype(&set_simulation_algorithm)>>::type;


}
} // namespace macrodr



#endif // CLI_MACRO_DR_BASE_H
