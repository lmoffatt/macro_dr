
#include "../catch2/catch.hpp"
#include "../CLI_function_table.h"
#include "../lapack_headers.h"

TEST_CASE("calculation of the likelihood of  fractions of a simulation with N_states", "[likelihood of fraction]")
{
    using namespace  macrodr;
    using namespace cmd;
    
    std::string lik_name ="../macro_dr/test_macro_dr/examples/fractions_likelihood";

   std::string model ="scheme_4_inact";
  macrodr::cmd::parameters_value_type par=macrodr::cmd::load_Parameter_value("../macro_dr/models/scheme_4_inact_par.csv", ",");
    
    
    
    std::string save_name ="../macro_dr/test_macro_dr/examples/fractions_N_simulation";
    
    
    std::string fsimulation="../macro_dr/test_macro_dr/examples/N_simulation.csv";
    
    std::string exp_filename="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time.txt";
    std::string exp_segments="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time_segments.txt";
    
    experiment_type fexperiment=get_Experiment(exp_filename,50e3,  0);
    fraction_algo_type fraction_algo=set_Fraction_algorithm(5,6,exp_segments);
    
    std::string sep=",";
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(true)> sim;
    auto Maybe_frac=calc_simulation_fractions(save_name,fsimulation,  fexperiment,  fraction_algo,  8,1);
    REQUIRE(Maybe_frac);
    
    bool adaptive_aproximation=true;
    bool recursive_approximation=true;
    int averaging_approximation=2;
    bool variance_correction_approximation=true;
    bool variance_correction=false;
    std::size_t n_sub_dt=100;
    auto lik_algo=set_Likelihood_algorithm(adaptive_aproximation,recursive_approximation,
                                             averaging_approximation,variance_correction_approximation,variance_correction,n_sub_dt);
    std::string ft_name ="../macro_dr/test_macro_dr/examples/ft_likelihood";
    
   auto ft= get_function_Table_maker_value(ft_name,0);
    
    
    auto Maybe_lik=calc_fraction_likelihood(
        lik_name,  model,  par,
        Maybe_frac,
         lik_algo,  ft);

    
    REQUIRE(Maybe_lik);
}





TEST_CASE("simulate a model and compare the partial likelihood with the expected partial likelihood", "[likelihood of fraction]")
{
    using namespace  macrodr;
    using namespace cmd;
    
    // load model
    // perform simulation
    
    std::string model = "scheme_4_inact";
    std::string filename_prefix="../macro_dr/test_macro_dr/examples/simulation_";
    recording_type experiment_observations="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time_recording.txt";
    
    auto experiment_cond= get_Experiment( "../macro_dr/experiments/Moffatt_Hume_2007_ATP_time.txt",
                                     50e3,  0);
    
    parameters_value_type param_file=load_Parameter_value("../macro_dr/models/scheme_4_inact_par.csv",",");
    
    auto sim_algo=set_simulation_algorithm(true, 100);
    
    auto sim_file=run_simulation(filename_prefix,
                   experiment_observations,
                   experiment_cond,0ul, model,
                                   param_file, sim_algo);
    
    
    
    std::string lik_name ="../macro_dr/test_macro_dr/examples/sim_fractions_likelihood";
    
    
    std::string save_name ="../macro_dr/test_macro_dr/examples/sim_fractions_N_simulation";
    
    
    std::string fsimulation=sim_file;
    
    std::string exp_filename="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time.txt";
    std::string exp_segments="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time_segments.txt";
    
    experiment_type fexperiment=get_Experiment(exp_filename,50e3,  0);
    fraction_algo_type fraction_algo=set_Fraction_algorithm(5,6,exp_segments);
    
    std::string sep=",";
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(true)> sim;
    auto Maybe_frac=calc_simulation_fractions(save_name,fsimulation,  fexperiment,  fraction_algo,  8,1);
    REQUIRE(Maybe_frac);
    
    bool adaptive_aproximation=true;
    bool recursive_approximation=true;
    int averaging_approximation=2;
    bool variance_correction_approximation=false;
    bool variance_correction=false;
    std::size_t n_sub_dt=100;
    auto lik_algo=set_Likelihood_algorithm(adaptive_aproximation,recursive_approximation,
                                             averaging_approximation,variance_correction_approximation,variance_correction,n_sub_dt);
    std::string ft_name ="../macro_dr/test_macro_dr/examples/ft_likelihood";
    
    auto ft= get_function_Table_maker_value(ft_name,0);
    
    
    auto Maybe_lik=calc_fraction_likelihood(
        lik_name,  model,  param_file,
        Maybe_frac,
        lik_algo,  ft);
    
    
    REQUIRE(Maybe_lik);
}
