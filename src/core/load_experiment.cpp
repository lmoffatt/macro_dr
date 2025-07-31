#include <macrodr/cmd/load_experiment.h>

#include "lapack_headers.h"

namespace macrodr {

namespace cmd {

Experiment load_experiment(std::string filename, double frequency_of_sampling, double initial_ATP) {
    Experiment e;
    auto Maybe_exp = load_experiment(filename, ",", frequency_of_sampling, initial_ATP, e);
    if (Maybe_exp) {
        return e;
    } else {
        auto [recording_conditions, recording] = macrodr::load_recording(filename);

        return Experiment(std::move(recording_conditions),
                          Frequency_of_Sampling(frequency_of_sampling),
                          initial_ATP_concentration(ATP_concentration(initial_ATP)));
    }
}
}  // namespace cmd
}  // namespace macrodr
