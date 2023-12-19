#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "parallel_tempering.h"
#include "variables.h"
#include <cstddef>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
namespace macrodr {

using var::Var;
using var::Vector_Space;

class Time : public Var<Time, double> {};
class number_of_samples : public var::Constant<number_of_samples, double> {};

class ATP_concentration: public Var<ATP_concentration, double> {};

using ATP_step=var::Vector_Space<number_of_samples,ATP_concentration>;

using ATP_evoltype= std::variant<ATP_step,std::vector<ATP_step>>;

class ATP_evolution: public Var<ATP_evolution, std::variant<ATP_step,std::vector<ATP_step>>>{
    using Var<ATP_evolution, std::variant<ATP_step,std::vector<ATP_step>>>::Var;
};


inline std::ostream& put(std::ostream& f,std::string sep, std::size_t i_frac, std::size_t i_step, double time, std::vector<ATP_step> const& v)
{
    for (std::size_t i=0; i<v.size(); ++i)
        f<<i_frac<<sep<<i_step<<sep<<time<<sep<<i_step+(i+0.5)/v.size()<<sep<<get<number_of_samples>(v[i])<<sep
          <<get<ATP_concentration>(v[i])<<"\n";
    return f; 
}
inline std::ostream& put(std::ostream& f,std::string sep, std::size_t i_frac, std::size_t i_step, double time, ATP_step const& x)
{
        f<<i_frac<<sep<<i_step<<sep<<time<<sep<<i_step+0.5<<sep<<get<number_of_samples>(x)<<sep
          <<get<ATP_concentration>(x)<<"\n";
    return f; 
}

inline std::ostream& put(std::ostream& f,std::string sep, std::size_t i_frac, std::size_t i_step, double time, ATP_evolution const& v)
{
    return std::visit([&f,sep,i_frac,i_step,time](auto& e)->decltype(auto){return put(f,sep,i_frac,i_step,time,e);},v());
}


class Patch_current : public Var<Patch_current, double> {};




using Experiment_step=    var::Vector_Space<Time, ATP_evolution>;


class Recording : public Var<Recording,std::vector<Patch_current>>{};


class Recording_conditions : public Var<Recording_conditions,std::vector<Experiment_step>>{};

class Frequency_of_Sampling : public var::Var<Frequency_of_Sampling, double> {};

class initial_ATP_concentration
    : public Var<initial_ATP_concentration, ATP_concentration> {};



using Experiment=var::Vector_Space<Recording_conditions,Frequency_of_Sampling, initial_ATP_concentration>;
}

template<class Parameter>
class save_Parameter;
namespace macrodr {
template<class Parameter>
void report_model(save_Parameter<Parameter>& s, std::vector<Experiment> const & e)
{
    std::ofstream f(s.fname+"_experiment.csv");
    f<<std::setprecision(std::numeric_limits<double>::digits10 + 1);
    f<<"frequency_of_sampling\n"<<get<Frequency_of_Sampling>(e[0])()<<"\n";
    f<<"initial_ATP_concentration\n"<<get<initial_ATP_concentration>(e[0])()<<"\n";
    
    f<<"i_frac"<<s.sep<<"i_step"<<s.sep<<"time"<<s.sep<<"i_step_f"<<s.sep<<"number_of_samples"<<s.sep
      <<"ATP"<<"\n";
                                                                            
    for (auto i_frac=0ul; i_frac<e.size(); ++i_frac)
    {
        auto r=get<Recording_conditions>(e[i_frac]);
        for (auto i_step=0ul; i_step<r().size(); ++i_step)
        {
            auto sa=r()[i_step];
            auto t=get<Time>(sa)();
            auto ev=get<ATP_evolution>(sa);
            put(f,s.sep,i_frac,i_step,t,ev);
        }
        
    }
    
}


template<class Parameter>
void report_model(save_Parameter<Parameter>& s, std::vector<Recording> const & e)
{
    std::ofstream f(s.fname+"_recording.csv");
    f<<std::setprecision(std::numeric_limits<double>::digits10 + 1);
    
    f<<"i_frac"<<s.sep<<"i_step"<<s.sep<<"patch_current"<<"\n";
    
    for (auto i_frac=0ul; i_frac<e.size(); ++i_frac)
    {
        for (auto i_step=0ul; i_step<e[i_frac]().size(); ++i_step)
        {
            f<<i_frac<<s.sep<<i_step<<s.sep<<e[i_frac]()[i_step]()<<"\n";
        }
        
    }
    
}




inline auto &extract_double(std::istream &is, double &r) {
  std::string s;
  is >> s;
  if ((s == "nan") || (s == "NAN") || (s == "NaN"))
    r = std::numeric_limits<double>::quiet_NaN();
  else
    r = std::stod(s);
  return is;
}

inline std::tuple<Recording_conditions,Recording> load_recording(const std::string filename) {
  std::ifstream f(filename);
  std::vector<Experiment_step> out0;
  std::vector<Patch_current> out1;
  
  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    if (!line.empty()) {
      std::stringstream ss(line);
      double v_time, v_ns, v_ATP, v_current;
      extract_double(ss >> v_time >> v_ns >> v_ATP, v_current);
      ATP_evolution step;
      step()=ATP_step(number_of_samples(std::size_t(v_ns)), ATP_concentration(v_ATP));
      out0.emplace_back(Time(v_time/1000), std::move(step));
      out1.emplace_back(Patch_current(v_current));
    }
  }
  return std::tuple(Recording_conditions(out0), Recording(out1));
}

inline std::tuple<Experiment, Recording>
load_experiment(const std::string filename, double frequency_of_sampling, double initial_ATP) {
  auto [v_recording_conditions,v_recording]=load_recording(filename);
  
  return {Experiment(std::move(v_recording_conditions),Frequency_of_Sampling(frequency_of_sampling), initial_ATP_concentration(ATP_concentration(initial_ATP))),
          std::move(v_recording)};
  
}


} // namespace macrodr

#endif // EXPERIMENT_H
