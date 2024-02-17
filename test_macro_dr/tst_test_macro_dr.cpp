#include "../catch2/catch.hpp"
#include "../qmodel.h"

TEST_CASE("duality of saving and loading simulations", "[save_simulation/load_simulation]")
{
    using namespace macrodr;
    std::string fname="../macro_dr/test_macro_dr/examples/N_simulation.csv";
    std::string sep=",";
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(true)> sim;
    auto Maybe_sim=macrodr::load_simulation(fname,sep,sim);
    REQUIRE(Maybe_sim);
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(false)> y;
    auto Maybe_sim2= macrodr::load_simulation(fname, sep, y);
    
    REQUIRE(Maybe_sim2);
    
    REQUIRE(get<Recording>(sim())==get<Recording>(y()));
    
    
    std::string temp="../macro_dr/test_macro_dr/examples/N_simulation_temp.csv";
    std::string temp2="../macro_dr/test_macro_dr/examples/simulation_temp.csv";
        
    macrodr::save_Simulated_Recording(temp,sep,sim);
    macrodr::save_Simulated_Recording(temp2,sep,y);
    
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(true)> sim_again;
    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(false)> y_again;
    
    auto Maybe_sim_again=macrodr::load_simulation(temp,sep,sim_again);
    REQUIRE(Maybe_sim_again);
    auto Maybe_y_again=macrodr::load_simulation(temp,sep,y_again);
    REQUIRE(Maybe_y_again);
    
    REQUIRE(get<Recording>(sim())==get<Recording>(sim_again()));
    REQUIRE(get<Recording>(sim())==get<Recording>(y_again()));
    
    REQUIRE(get<N_Ch_State_Evolution>(sim())==get<N_Ch_State_Evolution>(sim_again()));
    
    
    std::ifstream f0(fname);
    std::stringstream b0;
    b0 << f0.rdbuf();
    
    std::ifstream f1(temp);
    std::stringstream b1;
    b1 << f1.rdbuf();
    REQUIRE(b0.str()==b1.str());    
}

