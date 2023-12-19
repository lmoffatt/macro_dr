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
    
}
