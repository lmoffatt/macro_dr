#include <macrodr/cmd/load_experiment.h>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "experiment.h"
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

Experiment create_experiment(std::vector<std::tuple<std::size_t,std::size_t,double>> repetitions,
    double fs,
    double atp0, 
    double t0){
    Frequency_of_Sampling fs_out(fs);
    initial_ATP_concentration atp0_out{ATP_concentration{atp0}}; 
    Time t0_out(t0);

    std::vector<std::pair<std::size_t,std::vector<ATP_step>>> repetitions_out;
    repetitions_out.reserve(repetitions.size());
    for (auto& rep:repetitions){
       repetitions_out.emplace_back(std::get<0>(rep), 
       std::vector<ATP_step>{1,ATP_step{number_of_samples{get<1>(rep)}, ATP_concentration{get<2>(rep)}}}); 
    }
    return build_experiment(repetitions_out,fs_out,atp0_out,t0_out);
}

 Maybe_error<Recording> load_recording(const std::string& filename) {
    Recording e;
    auto Maybe_e= load_Recording_Data(filename, ",",e);
    if (!Maybe_e){
        return Maybe_e.error(); 

    }
    return e;
    }

 Recording define_recording(std::vector<double> values){
    Recording r;
    r().reserve(values.size());
    for (auto&& v:std::move(values))
    {
        r().emplace_back(v);
    }
    return r;

 }


} // namespace macrodr::cmd