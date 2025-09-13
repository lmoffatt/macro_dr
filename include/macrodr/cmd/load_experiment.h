#pragma once
#include "experiment.h"
#include "maybe_error.h"



namespace macrodr::cmd {

 Experiment load_experiment(const std::string& filename,
                           double frequency_of_sampling , double initial_ATP );


 Maybe_error<Recording> load_recording(const std::string& filename) ;

}
