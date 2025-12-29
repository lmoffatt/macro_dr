#include <macrodr/cmd/load_experiment.h>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "experiment.h"
#include "maybe_error.h"



namespace macrodr::cmd {

Experiment load_experiment(std::string const& filename, double frequency_of_sampling, double initial_agonist) {
   using namespace macrodr;
    Experiment e;
    auto Maybe_exp = load_experiment(filename, ",", frequency_of_sampling, initial_agonist, e);
    if (Maybe_exp) {
        return e;
    }          
    auto [recording_conditions, recording] = macrodr::load_recording(filename);

        return {std::move(recording_conditions),
                          Frequency_of_Sampling(frequency_of_sampling),
                          initial_agonist_concentration(Agonist_concentration(initial_agonist))};
   
}

Experiment create_experiment(std::vector<std::tuple<std::size_t,std::size_t,double>> repetitions,
    double fs,
    double agonist0, 
    double t0){
    Frequency_of_Sampling fs_out(fs);
    initial_agonist_concentration agonist0_out{Agonist_concentration{agonist0}}; 
    Time t0_out(t0);

    std::vector<std::pair<std::size_t,std::vector<Agonist_step>>> repetitions_out;
    repetitions_out.reserve(repetitions.size());
    for (auto& rep:repetitions){
       repetitions_out.emplace_back(std::get<0>(rep), 
       std::vector<Agonist_step>{1,Agonist_step{number_of_samples{get<1>(rep)}, Agonist_concentration{get<2>(rep)}}}); 
    }
    return build_experiment(repetitions_out,fs_out,agonist0_out,t0_out);
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
 Maybe_error<std::string> write_csv(Experiment const& e, Recording const& r, std::string  path)
 {
    auto path_=path+".csv";
     std::ofstream f(path_); 
     if (!f.is_open())
     {
      return error_message("cannot open ",path_);
     }
     auto& E=get<Recording_conditions>(e);

     auto n_samples=E().size(); 
     if (n_samples!= r().size()){
        return error_message("Experiment samples", n_samples," differ from Recording samples ",r().size());
     }
     
     double ns=1.0/get<Frequency_of_Sampling>(e)();
     f<<"n_step"<<","<<"step_start"<<","<<"step_end"<<","<<"step_middle"<<","<<"Agonist"<<","<<"patch_current"<<"\n";
      
     double step_start=0.0;
    
     for (std::size_t i=0; i<n_samples; ++i){
        auto const& ev= get<Agonist_evolution>(E()[i])();
        for (std::size_t j=0; j<ev.size(); ++j){
        auto n_step=static_cast<double>(i)+static_cast<double>(j)/static_cast<double>(ev.size());
        double step_end= step_start+ get<number_of_samples>(ev[j])()*ns;
        double step_middle=(step_start+step_end)/2;
        double Agonist = get<Agonist_concentration>(ev[j])();
        double Patch_current =r()[i]();
        f<<n_step<<","<<step_start<<","<<step_end<<","<<step_middle<<","<<Agonist<<","<<Patch_current<<"\n";
        step_start=step_end;    
        }
     }
     return path_; 
 }


 

} // namespace macrodr::cmd