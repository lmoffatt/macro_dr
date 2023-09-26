#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "variables.h"
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
namespace macrodr {

using var::Var;

class Time : public Var<Time, double> {};
class number_of_samples : public Var<Time, std::size_t> {};


class ATP_concentration
    : public Var<ATP_concentration, double> {};

class Patch_current : public Var<Patch_current, double> {};




using Experiment_step=    var::Vector_Space<Time, number_of_samples, ATP_concentration>;


class Recording : public Var<Recording,std::vector<Patch_current>>{};


class Recording_conditions : public Var<Recording,std::vector<Experiment_step>>{};

class Frequency_of_Sampling : public var::Var<Frequency_of_Sampling, double> {};

class initial_ATP_concentration
    : public Var<initial_ATP_concentration, ATP_concentration> {};



using Experiment=var::Vector_Space<Recording_conditions,Frequency_of_Sampling, initial_ATP_concentration>;






auto &extract_double(std::istream &is, double &r) {
  std::string s;
  is >> s;
  if ((s == "nan") || (s == "NAN") || (s == "NaN"))
    r = std::numeric_limits<double>::quiet_NaN();
  else
    r = std::stod(s);
  return is;
}

std::tuple<Recording_conditions,Recording> load_recording(const std::string filename) {
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
      out0.emplace_back(Time(v_time/1000), number_of_samples(std::size_t(v_ns)), ATP_concentration(v_ATP));
      out1.emplace_back(Patch_current(v_current));
    }
  }
  return std::tuple(Recording_conditions(out0), Recording(out1));
}

std::tuple<Experiment, Recording>
load_experiment(const std::string filename, double frequency_of_sampling, double initial_ATP) {
  auto [v_recording_conditions,v_recording]=load_recording(filename);
  
  return {Experiment(std::move(v_recording_conditions),Frequency_of_Sampling(frequency_of_sampling), initial_ATP_concentration(ATP_concentration(initial_ATP))),
          std::move(v_recording)};
  
}


} // namespace macrodr

#endif // EXPERIMENT_H
