#include "experiment.h"
#include "lexer_untyped.h"
#include <fstream>
#include <iostream>

#include "lapack_headers.h"
#include "matrix.h"

#include "qmodel.h"
using namespace macrodr;
int main(int argc, char **argv) {
  //        auto filename=argv[1];
//  auto filename = "../macro_dr/test.txt";
  
  auto Efilename = "../macro_dr/Moffatt_Hume_2007_ATP_2.txt";
  
  auto recording=macrodr::load_experiment(Efilename);
  
  
 // std::cerr<<recording.value();
 // std::cerr<<std::numeric_limits<double>::quiet_NaN();
//  aram_0 = State_Model_Parameters ( values = { { "kon"  50 }  { "koff"   20 }  { "beta"  500 }  { "alfa"  400 } { "g"  16.59e-3 } {"Number_of_Channels"  100}  { "gaussian_noise"  1.0e-5 } } )
  
  auto experiment=Experiment(std::move(recording),Frequency_of_Sampling(50e3), initial_ATP_concentration(ATP_concentration(0.0)));
  
  Patch_Model model1(N_St(5),
                     Q0(Matrix<double>(5,5,0.0)), Qa(Matrix<double>(5,5,0.0)), g(Matrix<double>(5,1,0.0)), N_Ch_mean(100), N_Ch_std(2.0), curr_noise(100.0), min_P(1e-7), Probability_error_tolerance(1e-2), Conductance_variance_error_tolerance(1e-2));
  
  get<Qa>(model1)()(0,1)=18;
  get<Qa>(model1)()(1,2)=12;
  get<Qa>(model1)()(2,3)=6;
  get<Q0>(model1)()(1,0)=210;
  get<Q0>(model1)()(2,1)=420;
  get<Q0>(model1)()(3,2)=630;
  get<Q0>(model1)()(3,4)=1680;
  get<Q0>(model1)()(4,3)=54;
  
  get<g>(model1)()(4,0)=0.5;
  
  
  std::random_device rd;
  typename std::mt19937_64::result_type seed =rd();
  std::mt19937_64 mt(seed);
  auto sim=Macro_DMR{}.sample(mt,model1,experiment,Simulation_Parameters(Number_of_simulation_sub_steps(1)));
  
  
  std::cerr<<sim;
  
  auto lik=Macro_DMR{}.log_Likelihood(model1,sim.value()());
  
  std::cerr<<"likelihood!!! "<<lik<<"\n";
  
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
  if (p) {
    auto ss = p.value().str();
    std::cerr << ss;
  } else
    std::cerr << p.error()();
}
