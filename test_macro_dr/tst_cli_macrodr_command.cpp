
#include "../general_output_operator.h"
#include "../catch2/catch.hpp"
#include "../experiment.h"
#include "../qmodel.h"
#include <cstddef>
#include <limits>


TEST_CASE("calculation of fractions of a simulation with N_states", "[save_fractioned_simulation/load_fractioned_simulation]")
{
    /**
     *inline Maybe_error<std::tuple<std::string, std::string, double, double>>
calc_simulation_fractions(std::string save_name,std::string simulation, experiment_type experiment,                                  fraction_algo_type fraction_algo, std::string model,std::size_t init_seed )
     */
    using namespace macrodr;
    using namespace cmd ;
    
    std::string save_name ="../macro_dr/test_macro_dr/examples/fractions_N_simulation";
    std::string save_name_2 ="../macro_dr/test_macro_dr/examples/fractions_N_simulation_2";
        
    
    std::string fsimulation="../macro_dr/test_macro_dr/examples/N_simulation.csv";
    
    std::string exp_filename="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time.txt";
    std::string exp_segments="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time_segments.txt";
    
    experiment_type fexperiment=get_Experiment(exp_filename,50e3,  0);
    fraction_algo_type fraction_algo=set_Fraction_algorithm(5,6,exp_segments);
    
    std::string model="scheme_4_inact";
    std::string sep=",";
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(true)> sim;
    auto Maybe_frac=calc_simulation_fractions(save_name,fsimulation,  fexperiment,  fraction_algo,  8,1);
    
    REQUIRE(Maybe_frac);
    auto [experiment, simulation, fs, iniATP] = Maybe_frac.value();
    
    std::vector<Experiment> xs;
    std::vector<Simulated_Recording<includes_N_state_evolution(true)>> ys;
    auto Maybe_e =
        load_fractioned_experiment(experiment, ",", fs, iniATP, xs);
    
    REQUIRE(Maybe_e);
    
    auto Maybe_y = load_fractioned_simulation(simulation, ",", ys);
    
    REQUIRE(Maybe_y);
    
    
    save_fractioned_experiment(save_name_2 + "_frac_experiment.csv", ",", xs);
    save_fractioned_simulation(save_name_2 + "_frac_simulation.csv", ",", ys);
    
    REQUIRE(compare_file_contents(simulation,save_name_2 + "_frac_simulation.csv"));
    REQUIRE(compare_file_contents(experiment,save_name_2 + "_frac_experiment.csv"));
    
    
    // now check with recordings
    
    std::vector<Recording> rs;
    
    auto Maybe_rs = load_fractioned_Recording(simulation, ",", rs);
    
    REQUIRE(Maybe_rs);
    
    for (std::size_t i_frac=0; i_frac<ys.size(); ++i_frac)
        REQUIRE(get<Recording>(ys[i_frac]())==rs[i_frac]);
    
    std::string rec_name ="../macro_dr/test_macro_dr/examples/fractions_recording.csv";
    save_fractioned_Recording(rec_name , ",", rs);
    
    std::vector<Recording> rs2;
    auto Maybe_rs2 = load_fractioned_Recording(rec_name, ",", rs2);
    
    REQUIRE(Maybe_rs2);
    REQUIRE(rs==rs2);
    
    
    
    
    
    
    
}





TEST_CASE("loading and saving experiments work together", "[save_experiment/load_experiment]")
{
    using namespace macrodr;
    using namespace cmd ;
    
    std::string save_name ="../macro_dr/test_macro_dr/examples/saved_experiment.txt";
    
    std::string exp_filename="../macro_dr/experiments/Moffatt_Hume_2007_ATP_time.txt";
    
    experiment_type fexperiment=get_Experiment(exp_filename,50e3,  0);
    std::string sep=",";
    save_experiment(save_name, sep,
                    fexperiment);
    
    Experiment exp_again;
    auto Maybe_exp=load_experiment(save_name,sep,50e3,0.0,exp_again);
    
    REQUIRE(Maybe_exp);
    
    CHECK(compare_contents(fexperiment, exp_again )==Maybe_error<bool>(true));
    
    std::string i_name ="../macro_dr/test_macro_dr/examples/idealized_experiment.txt";
    idealize_Experiment(exp_filename,sep,i_name);
    Experiment exp_ideal;
    auto Maybe_exp2=load_experiment(i_name,sep,50e3,0.0,exp_ideal);
    REQUIRE(Maybe_exp2);
    std::string i_name2 ="../macro_dr/test_macro_dr/examples/idealized_experiment_2.txt";
    
    save_experiment(i_name2, sep,
                    exp_ideal);
    
    CHECK(compare_file_contents(i_name,i_name2));
    
    
    
    
    
    
    
}
