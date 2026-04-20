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

 // Build a Recording of length `n_samples` filled with `fill_value`, then
 // overwrite each index in `missing_indices` with NaN. Indices ≥ n_samples
 // are silently ignored. NaN samples are detected downstream in the
 // likelihood path (qmodel.h: `isnan(y)` gating) as "no measurement": the
 // prediction is propagated forward without a Bayes update, and the
 // log-likelihood contribution at that step is zero. Useful for
 // pre-equilibration segments where the channel state should relax to
 // stationarity before any measurement is scored.
 Recording define_recording_with_gaps(std::size_t n_samples, double fill_value,
                                      std::vector<std::size_t> missing_indices);

 // Convenience for a single contiguous gap at the start, end, or interior
 // of the recording. Marks the half-open interval [missing_start, missing_end)
 // as NaN; everything else gets `fill_value`. Equivalent to
 // `define_recording_with_gaps(n_samples, fill_value, {missing_start, …,
 // missing_end-1})` but doesn't force the caller to list every index.
 Recording define_recording_with_missing_range(std::size_t n_samples, double fill_value,
                                               std::size_t missing_start,
                                               std::size_t missing_end);

Maybe_error<std::string> write_csv(Experiment const& e, Recording const& r, std::string path);


}
