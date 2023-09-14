#include "CLI_macro_dr.h"
#include "experiment.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "variables_derivative.h"
#include "parameters_derivative.h"
#include <fstream>
#include <iostream>
#include <map>


#include "qmodel.h"
using namespace macrodr;



int main(int argc, char **argv) {
  //        auto filename=argv[1];
  //  auto filename = "../macro_dr/test.txt";

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
  
  std::cerr<<p;
  
  if (p)
  { 
     auto c=dcli::compile_program(cm,p.value());
      if (c)
      {
         auto exec=c.value().run();
      }
  }

  auto Efilename = "../macro_dr/Moffatt_Hume_2007_ATP_2.txt";

  auto recording = macrodr::load_recording(Efilename);

  // std::cerr<<recording.value();
  // std::cerr<<std::numeric_limits<double>::quiet_NaN();
  //  aram_0 = State_Model_Parameters ( values = { { "kon"  50 }  { "koff"   20
  //  }  { "beta"  500 }  { "alfa"  400 } { "g"  16.59e-3 }
  //  {"Number_of_Channels"  100}  { "gaussian_noise"  1.0e-5 } } )

  auto experiment =
      Experiment(std::move(recording), Frequency_of_Sampling(50e3),
                 initial_ATP_concentration(ATP_concentration(0.0)));
  
  
  auto model1= Model1::Model([](const auto& p){
      return build<macrodr::Patch_Model>(N_St(5),
                                build<Q0>(var::build_<Matrix<double>>(5, 5, {{1,0}, {2,1},{3,2},{3,4},{4,3}},
                                                               {p()[5],p()[6],p()[7],p()[8],p()[9]})),
                                build<Qa>(var::build_<Matrix<double>>(5, 5, {{0, 1}, {1, 2},{2, 3}},{p()[2],p()[3],p()[4]})),
                                build<g>(var::build_<Matrix<double>>(5,1,{{4,0}},{p()[10]})),
                                build<N_Ch_mean>(p()[0]),
                                build<curr_noise>(p()[1]),
                                min_P(1e-7),
                                Probability_error_tolerance(1e-2),
                                Conductance_variance_error_tolerance(1e-2));
  });
  
  
  
  auto param1=Parameters<Model1>(Matrix<double>(1,11,std::vector<double>{100,0.05,18,12,6,210,420,630,1680,54,0.5}));
  
  
  auto dparam1=var::selfDerivative(param1);
  std::cerr<<"dparam1\n"<<dparam1;
  auto n= build<N_Ch_mean>(dparam1()[0]);
  
  auto dp0=dparam1()[0];
  
  auto qq=  var::build_<Matrix<double>>(5, 5, {{0, 1}, {1, 2},{2, 3}},{dparam1()[2],dparam1()[3],dparam1()[4]});
  
  
  auto dm=model1(dparam1);
  
  std::cerr<<"\ndm!!!\n"<<dm;
  auto fs = get<Frequency_of_Sampling>(experiment).value();
  auto ini = macrodr::Macro_DMR{}.init(dm, get<initial_ATP_concentration>(experiment));
  auto t_step=get<Recording>(experiment)()[0];
  auto t_Qx = macrodr::Macro_DMR{}.calc_eigen(dm, get<ATP_concentration>(t_step));
  
  auto t_Qdt = macrodr::Macro_DMR{}.calc_Qdt(dm, t_Qx.value(),
                        get<number_of_samples>(t_step).value() / fs);
  
  std::random_device rd;
  typename std::mt19937_64::result_type seed = rd();
  std::mt19937_64 mt(seed);
  auto sim = Macro_DMR{}.sample(
      mt, model1,param1, experiment,
      Simulation_Parameters(Number_of_simulation_sub_steps(10)));

  std::cerr << sim;
  
  
  auto lik = Macro_DMR{}.log_Likelihood(model1,dparam1, sim.value()());

  std::cerr << "likelihood!!! " << lik << "\n";

    if (p) {
    auto ss = p.value().str();
    std::cerr << ss;
  } else
    std::cerr << p.error()();
}
