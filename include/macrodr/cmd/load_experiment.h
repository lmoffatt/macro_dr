#ifndef LOAD_EXPERIMENT_H
#define LOAD_EXPERIMENT_H
#include "experiment.h"

namespace macrodr {

namespace cmd {

Experiment load_experiment(std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
                           double frequency_of_sampling = 50e3, double initial_ATP = 0);

}
}  // namespace macrodr
#endif  // LOAD_EXPERIMENT_H
