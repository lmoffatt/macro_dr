#pragma once
#include "experiment.h"
#include "maybe_error.h"



namespace macrodr::cmd {

 Experiment load_experiment(const std::string& filename,
                           double frequency_of_sampling , double initial_ATP );

 Experiment create_experiment(std::vector<std::tuple<std::size_t,std::size_t,double>> repetitions,
    double fs,
    double atp0, 
    double t0);

 Maybe_error<Recording> load_recording(const std::string& filename) ;

}
