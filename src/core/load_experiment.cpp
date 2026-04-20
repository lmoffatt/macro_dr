#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/detail/write_csv_common.h>
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

 Recording define_recording_with_gaps(std::size_t n_samples, double fill_value,
                                      std::vector<std::size_t> missing_indices) {
    Recording r;
    r().reserve(n_samples);
    for (std::size_t i = 0; i < n_samples; ++i) {
        r().emplace_back(fill_value);
    }
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t idx : missing_indices) {
        if (idx < n_samples) {
            r()[idx]() = nan_value;
        }
    }
    return r;
 }

 Recording define_recording_with_missing_range(std::size_t n_samples, double fill_value,
                                               std::size_t missing_start,
                                               std::size_t missing_end) {
    Recording r;
    r().reserve(n_samples);
    for (std::size_t i = 0; i < n_samples; ++i) {
        r().emplace_back(fill_value);
    }
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const std::size_t lo = std::min(missing_start, n_samples);
    const std::size_t hi = std::min(missing_end, n_samples);
    for (std::size_t i = lo; i < hi; ++i) {
        r()[i]() = nan_value;
    }
    return r;
 }
 Maybe_error<std::string> write_csv(Experiment const& e, Recording const& r, std::string  path)
 {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    const auto& conditions = get<Recording_conditions>(e);
    const auto n_samples = conditions().size();
    if (n_samples != r().size()) {
        return error_message("Experiment samples ", n_samples, " differ from Recording samples ",
                             r().size());
    }

    detail::CsvWriter writer(f, {});
    const double fs = get<Frequency_of_Sampling>(e)();
    for (std::size_t i = 0; i < n_samples; ++i) {
        double step_start = get<Time>(conditions()[i])();
        const auto& segments = get<Agonist_evolution>(conditions()[i])();
        for (std::size_t j = 0; j < segments.size(); ++j) {
            const double duration = get<number_of_samples>(segments[j])() / fs;
            const double step_end = step_start + duration;

            detail::CsvContext ctx;
            ctx.scope = "recording";
            ctx.sample_index = i;
            ctx.segment_index = j;
            ctx.n_step = static_cast<double>(i) +
                         static_cast<double>(j) / static_cast<double>(segments.size());
            ctx.time_start = step_start;
            ctx.time_end = step_end;
            ctx.time_middle = 0.5 * (step_start + step_end);
            ctx.agonist = get<Agonist_concentration>(segments[j])();
            ctx.patch_current = r()[i]();

            auto ok = detail::emit_named_component(writer, ctx, "patch_current", r()[i]);
            if (!ok || !ok.value()) {
                return ok.error()();
            }
            step_start = step_end;
        }
    }
    return path_;
 }


 

} // namespace macrodr::cmd
