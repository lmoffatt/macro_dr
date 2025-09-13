#include <macrodr/cmd/load_experiment.h>

#include "experiment.h"
#include "lapack_headers.h"
#include "maybe_error.h"



namespace macrodr::cmd {

Experiment load_experiment(std::string const& filename, double frequency_of_sampling, double initial_ATP) {
   using namespace macrodr;
    Experiment e;
    auto Maybe_exp = load_experiment(filename, ",", frequency_of_sampling, initial_ATP, e);
    if (Maybe_exp) {
        return e;
    }          
    auto [recording_conditions, recording] = macrodr::load_recording(filename);

        return {std::move(recording_conditions),
                          Frequency_of_Sampling(frequency_of_sampling),
                          initial_ATP_concentration(ATP_concentration(initial_ATP))};
   
}

 Maybe_error<Recording> load_recording(const std::string& filename) {
    Recording e;
    auto Maybe_e= load_Recording_Data(filename, ",",e);
    if (!Maybe_e){
        return Maybe_e.error(); 

    }
    return e;
    }
} // namespace macrodr::cmd