#pragma once
#include <fstream>
#include "experiment.h"
#include "maybe_error.h"



namespace macrodr::cmd {

 Experiment load_experiment(const std::string& filename,
                           double frequency_of_sampling , double initial_agonist );

 Experiment create_experiment(std::vector<std::tuple<std::size_t,std::size_t,double>> repetitions,
    double fs,
    double agonist0, 
    double t0);

 Maybe_error<Recording> load_recording(const std::string& filename) ;

 Recording define_recording(std::vector<double> values);

 inline Maybe_error<std::string> write_csv(Experiment const& e, Recording const& r, std::string path)
 {
    path=path+".csv";
     std::ofstream f(path); 
     if (!f.is_open())
     {
      return error_message("cannot open ",path);
     }
     double fs=get<Frequency_of_Sampling>(e)();
     f<<"time"<<",";
     return path; 


 }


}
