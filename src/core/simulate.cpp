

#include <macrodr/cmd/simulate.h>
#include <macrodr/cmd/detail/write_csv_common.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <map>
#include <omp.h>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "CLI_macro_dr_base.h"
//#include "cuevi.h"
#include "experiment.h"
#include "maybe_error.h"
#include "mcmc.h"
#include "parallel_tempering.h"
#include "parameters.h"
//#include "parameters_derivative.h"
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/interface/IModel.h>
#include <variables.h>

#include "lapack_headers.h"
#include "models_used.h"
#include "qmodel.h"
#include "random_samplers.h"

namespace {

namespace csv_detail = macrodr::cmd::detail;
using namespace macrodr;

Maybe_error<Simulation_Parameters> simulation_parameters_from_algorithm(
    const std::string& simulation_algorithm) {
    if (simulation_algorithm == simulation_algorithm_uniformization_name) {
        return make_uniformization_simulation_parameters();
    }
    if (simulation_algorithm == simulation_algorithm_substeps_name) {
        return error_message(
            "simulation_algorithm=\"substeps\" is not supported by this overload; use "
            "number_of_substeps instead");
    }
    return error_message("unsupported simulation_algorithm \"", simulation_algorithm,
                         "\"; expected \"", simulation_algorithm_uniformization_name, "\"");
}

Maybe_error<Simulation_Parameters> simulation_parameters_for_substeps(std::size_t n_sub) {
    if (n_sub == 0) {
        return error_message("number_of_substeps must be greater than zero");
    }
    return make_substep_simulation_parameters(n_sub);
}

Maybe_error<Simulation_Parameters> simulation_parameters_from_legacy(
    const macrodr::cmd::simulation_algo_type& sim_algo_type) {
    if (sim_algo_type.algorithm == simulation_algorithm_uniformization_name) {
        return make_uniformization_simulation_parameters();
    }
    if (sim_algo_type.algorithm == simulation_algorithm_substeps_name) {
        return simulation_parameters_for_substeps(sim_algo_type.number_of_substeps);
    }
    return error_message("unsupported simulation algorithm \"", sim_algo_type.algorithm, "\"");
}

Maybe_error<bool> emit_simulation_rows(csv_detail::CsvWriter& writer, const Experiment& e,
                                       const Recording& recording,
                                       std::optional<std::size_t> simulation_index,
                                       const csv_detail::CsvContext& base_ctx = {}) {
    const auto& conditions = get<Recording_conditions>(e);
    const auto n_samples = conditions().size();
    if (n_samples != recording().size()) {
        return error_message("Experiment samples ", n_samples, " differ from Recording samples ",
                             recording().size());
    }

    const double fs = get<Frequency_of_Sampling>(e)();
    for (std::size_t i = 0; i < n_samples; ++i) {
        double step_start = get<Time>(conditions()[i])();
        const auto& segments = get<Agonist_evolution>(conditions()[i])();
        for (std::size_t j = 0; j < segments.size(); ++j) {
            const double duration = get<number_of_samples>(segments[j])() / fs;
            const double step_end = step_start + duration;

            auto ctx = base_ctx;
            ctx.scope = "simulation";
            ctx.simulation_index = simulation_index;
            ctx.sample_index = i;
            ctx.segment_index = j;
            ctx.n_step = static_cast<double>(i) +
                         static_cast<double>(j) / static_cast<double>(segments.size());
            ctx.time_start = step_start;
            ctx.time_end = step_end;
            ctx.time_middle = 0.5 * (step_start + step_end);
            ctx.agonist = get<Agonist_concentration>(segments[j])();
            ctx.patch_current = recording()[i]();

            auto ok =
                csv_detail::emit_named_component(writer, ctx, "patch_current", recording()[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
            step_start = step_end;
        }
    }
    return true;
}

template <class SubEvolution>
Maybe_error<bool> emit_simulation_rows_with_sub(csv_detail::CsvWriter& writer, const Experiment& e,
                                                const Recording& recording,
                                                const SubEvolution& sub_evolution,
                                                std::size_t n_sub,
                                                std::optional<std::size_t> simulation_index,
                                                const csv_detail::CsvContext& base_ctx = {}) {
    if (n_sub == 0) {
        return error_message("number_of_substates must be greater than zero");
    }

    const auto& conditions = get<Recording_conditions>(e);
    const auto n_samples = conditions().size();
    if (n_samples != recording().size()) {
        return error_message("Experiment samples ", n_samples, " differ from Recording samples ",
                             recording().size());
    }

    std::size_t expected_sub = 0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        expected_sub += get<Agonist_evolution>(conditions()[i])().size() * n_sub;
    }
    if (expected_sub != sub_evolution().size()) {
        return error_message("Sub evolution samples ", sub_evolution().size(),
                             " differ from expected ", expected_sub);
    }

    const double fs = get<Frequency_of_Sampling>(e)();
    std::size_t sub_i = 0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        double step_start = get<Time>(conditions()[i])();
        const auto& segments = get<Agonist_evolution>(conditions()[i])();
        for (std::size_t j = 0; j < segments.size(); ++j) {
            const double duration = get<number_of_samples>(segments[j])() / fs;
            const double step_end = step_start + duration;
            const double step_mid = 0.5 * (step_start + step_end);
            const double agonist = get<Agonist_concentration>(segments[j])();

            auto full_ctx = base_ctx;
            full_ctx.scope = "simulation";
            full_ctx.simulation_index = simulation_index;
            full_ctx.sample_index = i;
            full_ctx.segment_index = j;
            full_ctx.n_step = static_cast<double>(i) +
                              static_cast<double>(j) / static_cast<double>(segments.size());
            full_ctx.time_start = step_start;
            full_ctx.time_end = step_end;
            full_ctx.time_middle = step_mid;
            full_ctx.agonist = agonist;
            full_ctx.patch_current = recording()[i]();

            auto ok = csv_detail::emit_named_component(writer, full_ctx, "patch_current",
                                                       recording()[i]);
            if (!ok || !ok.value()) {
                return ok;
            }

            for (std::size_t k = 0; k < n_sub; ++k) {
                const double frac0 = static_cast<double>(k) / static_cast<double>(n_sub);
                const double frac1 = static_cast<double>(k + 1) / static_cast<double>(n_sub);
                const double sub_start = step_start + duration * frac0;
                const double sub_end = step_start + duration * frac1;

                auto sub_ctx = base_ctx;
                sub_ctx.scope = "simulation_sub";
                sub_ctx.simulation_index = simulation_index;
                sub_ctx.sample_index = i;
                sub_ctx.segment_index = j;
                sub_ctx.sub_index = k;
                sub_ctx.n_step =
                    static_cast<double>(i) +
                    (static_cast<double>(j) + static_cast<double>(k + 1) / static_cast<double>(n_sub)) /
                        static_cast<double>(segments.size());
                sub_ctx.time_start = sub_start;
                sub_ctx.time_end = sub_end;
                sub_ctx.time_middle = 0.5 * (sub_start + sub_end);
                sub_ctx.agonist = agonist;
                sub_ctx.patch_current = sub_evolution()[sub_i]();

                ok = csv_detail::emit_named_component(writer, sub_ctx, "patch_current",
                                                      sub_evolution()[sub_i]);
                if (!ok || !ok.value()) {
                    return ok;
                }
                ++sub_i;
            }
            step_start = step_end;
        }
    }
    return true;
}

template <class IndexedValue, class EmitRows>
Maybe_error<std::string> write_indexed_simulation_csv(const IndexedValue& indexed,
                                                      std::string path, EmitRows&& emit_rows) {
    auto valid = csv_detail::validate_indexed_write_csv_value(indexed, "simulation");
    if (!valid) {
        return valid.error();
    }

    return csv_detail::write_indexed_rows_csv(
        indexed.index_space(), {}, std::move(path),
        [&](csv_detail::CsvWriter& writer, const csv_detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_value = indexed.at(coord);
            if (!maybe_value) {
                return maybe_value.error();
            }
            return emit_rows(writer, base_ctx, maybe_value.value().get());
        });
}

}  // namespace

namespace macrodr::cmd {

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed) {
    myseed = calc_seed(myseed);

    mt_64i mt(myseed);
    auto maybe_sim = simulation_parameters_for_substeps(n_sub);
    if (!maybe_sim) {
        return maybe_sim.error();
    }
    auto sim = maybe_sim.value();
    return Macro_DMR{}.sample(mt, *model, par, e, sim, r);
   }

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed) {
    return run_simulations(model, par.to_value(), e, r, n_sub, myseed);
}

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_values& par, const Experiment& e,
    const Recording& r, std::string simulation_algorithm, std::size_t number_of_substeps, std::size_t myseed) {
    myseed = calc_seed(myseed);
    if (simulation_algorithm == simulation_algorithm_substeps_name) {
        return run_simulations(model, par, e, r, number_of_substeps, myseed);
    }


    auto maybe_sim = simulation_parameters_from_algorithm(simulation_algorithm);
    if (!maybe_sim) {
        return maybe_sim.error();
    }

    mt_64i mt(myseed);
    return Macro_DMR{}.sample(mt, *model, par, e, maybe_sim.value(), r);
}

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, std::string simulation_algorithm, std::size_t number_of_substeps, std::size_t myseed) {
    return run_simulations(model, par.to_value(), e, r, simulation_algorithm, number_of_substeps,  myseed);
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    std::size_t n_simulations, const Experiment& e, const Recording& r, std::size_t n_sub,
    std::size_t myseed) {
    myseed = calc_seed(myseed);

    // Sim params are constant across recordings — compute once.
    auto maybe_params = simulation_parameters_for_substeps(n_sub);
    if (!maybe_params) {
        return maybe_params.error();
    }
    auto sim = maybe_params.value();

    // One independent RNG per recording: seeds are pre-generated SERIALLY from the
    // master mt → reproducible regardless of thread schedule (NOT a per-thread pool,
    // whose state would depend on the runtime schedule). Recordings are independent
    // and sample() reads a const model + a per-recording mt (no shared state), so
    // the loop parallelizes. NB: not bit-identical to the old single-stream serial
    // sim — it is a different but equally valid, reproducible draw.
    mt_64i mt(myseed);
    std::vector<decltype(mt())> seeds(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) seeds[i] = mt();
    std::vector<std::optional<Simulated_Recording<var::please_include<>>>> slots(n_simulations);
    std::vector<std::string> errs(n_simulations);
#pragma omp parallel for schedule(dynamic) if(n_simulations >= 2 && !omp_in_parallel())
    for (std::size_t i = 0; i < n_simulations; ++i) {
        mt_64i mt_i(seeds[i]);
        auto maybe_recording = Macro_DMR{}.sample(mt_i, *model, par, e, sim, r);
        if (maybe_recording)
            slots[i].emplace(std::move(maybe_recording.value()));
        else
            errs[i] = maybe_recording.error()();
    }
    std::vector<Simulated_Recording<var::please_include<>>> result;
    result.reserve(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) {
        if (!errs[i].empty()) return error_message(errs[i]);
        result.push_back(std::move(slots[i].value()));
    }
    return result;
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    std::size_t n_simulations, const Experiment& e, const Recording& r, std::size_t n_sub,
    std::size_t myseed) {
    return run_n_simulations(model, par.to_value(), n_simulations, e, r, n_sub, myseed);
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_values& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::string simulation_algorithm,std::size_t number_of_substeps, 
    std::size_t myseed) {
    if (simulation_algorithm == simulation_algorithm_substeps_name) {
      return run_n_simulations(model, par, n_simulations, e, r, number_of_substeps, myseed);
    }
        
    myseed = calc_seed(myseed);

    auto maybe_sim_params = simulation_parameters_from_algorithm(simulation_algorithm);
    if (!maybe_sim_params) {
        return maybe_sim_params.error();
    }

    // One independent RNG per recording: seeds pre-generated serially from the
    // master mt → reproducible, schedule-independent. sample() reads a const model
    // + per-recording mt (no shared state) → safe to parallelize. Not bit-identical
    // to the old single-stream serial draw (a different, reproducible realization).
    mt_64i mt(myseed);
    std::vector<decltype(mt())> seeds(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) seeds[i] = mt();
    std::vector<std::optional<Simulated_Recording<var::please_include<>>>> slots(n_simulations);
    std::vector<std::string> errs(n_simulations);
#pragma omp parallel for schedule(dynamic) if(n_simulations >= 2 && !omp_in_parallel())
    for (std::size_t i = 0; i < n_simulations; ++i) {
        mt_64i mt_i(seeds[i]);
        auto maybe_sim = Macro_DMR{}.sample(mt_i, *model, par, e, maybe_sim_params.value(), r);
        if (maybe_sim)
            slots[i].emplace(std::move(maybe_sim.value()));
        else
            errs[i] = maybe_sim.error()();
    }
    std::vector<Simulated_Recording<var::please_include<>>> result;
    result.reserve(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) {
        if (!errs[i].empty()) return error_message(errs[i]);
        result.push_back(std::move(slots[i].value()));
    }
    return result;
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par, std::size_t n_simulations,
    const Experiment& e, const Recording& r, std::string simulation_algorithm, std::size_t number_of_substeps, std::size_t myseed) {
    return run_n_simulations(model, par.to_value(), n_simulations, e, r, simulation_algorithm, number_of_substeps, myseed);
}
 

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
    run_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_values& par,const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    myseed = calc_seed(myseed);
    mt_64i mt(myseed);
    auto maybe_sim = simulation_parameters_for_substeps(n_sub);
    if (!maybe_sim) {
        return maybe_sim.error();
    }
    auto sim = maybe_sim.value();
      return  Macro_DMR{}.sample_sub_y(mt, *model, par, e, sim, r);
  }

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
    run_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_transformed& par, const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    return run_simulations_with_sub_intervals(model, par.to_value(), e, r, n_sub,
                                                myseed);
}

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
run_simulations_with_sub_intervals(const ModelPtr&, const var::Parameters_values&,
                                   const Experiment&, const Recording&,
                                   std::string simulation_algorithm, std::size_t) {
    if (simulation_algorithm == simulation_algorithm_uniformization_name) {
        return error_message(
            "simulate_with_sub_intervals does not support "
            "simulation_algorithm=\"uniformization\"; use number_of_substeps instead");
    }
    if (simulation_algorithm == simulation_algorithm_substeps_name) {
        return error_message(
            "simulation_algorithm=\"substeps\" is not supported by this overload; use "
            "number_of_substeps instead");
    }
    return error_message("unsupported simulation_algorithm \"", simulation_algorithm,
                         "\" for simulate_with_sub_intervals");
}

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
run_simulations_with_sub_intervals(const ModelPtr& model, const var::Parameters_transformed& par,
                                   const Experiment& e, const Recording& r,
                                   std::string simulation_algorithm,
                                   std::size_t myseed) {
    return run_simulations_with_sub_intervals(model, par.to_value(), e, r, simulation_algorithm,
                                              myseed);
}


Maybe_error<std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>>
    run_n_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_values& par,
                                         std::size_t n_simulations, const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    myseed = calc_seed(myseed);
    // Sim params constant — compute once.
    auto maybe_params = simulation_parameters_for_substeps(n_sub);
    if (!maybe_params) {
        return maybe_params.error();
    }
    auto sim = maybe_params.value();

    // Independent RNG per recording (serial seed pre-gen → reproducible, schedule-
    // independent). sample_sub_y reads a const model + per-recording mt → safe to
    // parallelize. Not bit-identical to the old single-stream serial draw.
    mt_64i mt(myseed);
    std::vector<decltype(mt())> seeds(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) seeds[i] = mt();
    std::vector<std::optional<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>>
        slots(n_simulations);
    std::vector<std::string> errs(n_simulations);
#pragma omp parallel for schedule(dynamic) if(n_simulations >= 2 && !omp_in_parallel())
    for (std::size_t i = 0; i < n_simulations; ++i) {
        mt_64i mt_i(seeds[i]);
        auto maybe_recording = Macro_DMR{}.sample_sub_y(mt_i, *model, par, e, sim, r);
        if (maybe_recording)
            slots[i].emplace(std::move(maybe_recording.value()));
        else
            errs[i] = maybe_recording.error()();
    }
    std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> result;
    result.reserve(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) {
        if (!errs[i].empty()) return error_message(errs[i]);
        result.push_back(std::move(slots[i].value()));
    }
    return result;
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>>
    run_n_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_transformed& par,
                                         std::size_t n_simulations, const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    return run_n_simulations_with_sub_intervals(model, par.to_value(), n_simulations, e, r, n_sub,
                                                myseed);
}

Maybe_error<std::string> runsimulation(std::string filename_prefix, recording_type recording_file,
                                       experiment_type experiment, std::size_t myseed,
                                       const std::string& modelName,
                                       parameters_value_type parameter_files,
                                       simulation_algo_type sim_algo_type) {
    auto maybe_sim_params = simulation_parameters_from_legacy(sim_algo_type);
    if (!maybe_sim_params) {
        return maybe_sim_params.error();
    }
    auto Maybe_recording = load_Recording(std::move(recording_file), std::string(","));
    if (!Maybe_recording) {
        return Maybe_recording.error()();
    }
    auto recording = std::move(Maybe_recording.value());
    auto Maybe_model_v = get_model(modelName);
    if (!Maybe_model_v)
        return Maybe_model_v.error()();
    else {
        auto model_v = std::move(Maybe_model_v.value());

        return std::visit(
            [&filename_prefix, &experiment, &recording, &myseed, &parameter_files,
             includeN = sim_algo_type.include_N_states,
             sim_params = maybe_sim_params.value()](auto model0ptr) -> std::string {
                auto& model0 = *model0ptr;
                myseed = calc_seed(myseed);
                mt_64i mt(myseed);
                auto Maybe_parameter_values =
                    var::load_Parameters(parameter_files.first, parameter_files.second,
                                         model0.model_name(), model0.names());
                std::string ModelName = model0.model_name();

                std::string filename = filename_prefix + "_" + ModelName + "_" + time_now() + "_" +
                                       std::to_string(myseed);

                if (!Maybe_parameter_values)
                    return Maybe_parameter_values.error()();
                else {
                    auto param1 = Maybe_parameter_values.value().standard_parameter();
                    save_Parameter<var::Parameters_transformed> s(filename, 1ul, 200ul);
                    if (!includeN) {
                        auto sim = Macro_DMR{}.sample(mt, model0, param1, experiment, sim_params,
                                                      recording);
                        if (!sim)
                            return sim.error()();
                        else {
                            save_Recording(filename + "_simulation.csv", ",",
                                           get<Recording>(sim.value()()));
                            return filename + "_simulation.csv";
                        }
                    } else {
                        auto sim =
                            Macro_DMR{}.sample_N(mt, model0, param1, experiment, sim_params,
                                                 recording);
                        if (!sim)
                            return sim.error()();
                        else {
                            save_Simulated_Recording(filename + "_N_simulation.csv", ",",
                                                     sim.value());
                            return filename + "_N_simulation.csv";
                        }
                    }
                }
            },
            model_v);
    }
}

Maybe_error<std::string> write_csv(Experiment const& e,
                                   std::vector<Simulated_Recording<var::please_include<>>> const& simulation,
                                    std::string path) {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, {});
    for (std::size_t sim_index = 0; sim_index < simulation.size(); ++sim_index) {
        auto ok = emit_simulation_rows(writer, e, get<Recording>(simulation[sim_index]()),
                                       sim_index);
        if (!ok || !ok.value()) {
            return ok.error()();
        }
    }
    return path_;
}

Maybe_error<std::string> write_csv(
    Experiment const& e,
    const var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>& simulation,
    std::string path) {
    return write_indexed_simulation_csv(
        simulation, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const auto& recordings) -> Maybe_error<bool> {
            for (std::size_t sim_index = 0; sim_index < recordings.size(); ++sim_index) {
                auto ok = emit_simulation_rows(writer, e, get<Recording>(recordings[sim_index]()),
                                               sim_index, base_ctx);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        });
}


Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<var::please_include<>> const& simulation,
                                   std::string path) {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, {});
    auto ok = emit_simulation_rows(writer, e, get<Recording>(simulation()), std::nullopt);
    if (!ok || !ok.value()) {
        return ok.error()();
    }
    return path_;
}

Maybe_error<std::string> write_csv(
    Experiment const& e,
    const var::Indexed<Simulated_Recording<var::please_include<>>>& simulation, std::string path) {
    return write_indexed_simulation_csv(
        simulation, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const auto& item) -> Maybe_error<bool> {
            return emit_simulation_rows(writer, e, get<Recording>(item()), std::nullopt, base_ctx);
        });
}


Maybe_error<std::string> write_csv(
    Experiment const& e,
    Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const& simulation,
    std::size_t n_sub, std::string path) {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, {});
    auto ok = emit_simulation_rows_with_sub(writer, e, get<Recording>(simulation()),
                                            get<Only_Ch_Curent_Sub_Evolution>(simulation()), n_sub,
                                            std::nullopt);
    if (!ok || !ok.value()) {
        return ok.error()();
    }
    return path_;
}

Maybe_error<std::string> write_csv(
    Experiment const& e,
    const var::Indexed<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>&
        simulation,
    std::size_t n_sub, std::string path) {
    return write_indexed_simulation_csv(
        simulation, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const auto& item) -> Maybe_error<bool> {
            return emit_simulation_rows_with_sub(writer, e, get<Recording>(item()),
                                                 get<Only_Ch_Curent_Sub_Evolution>(item()), n_sub,
                                                 std::nullopt, base_ctx);
        });
}

Maybe_error<std::string> write_csv(
    Experiment const& e,
    std::vector<    Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const& simulation,
    std::size_t n_sub, std::string path) {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, {});
    for (std::size_t sim_index = 0; sim_index < simulation.size(); ++sim_index) {
        auto ok = emit_simulation_rows_with_sub(
            writer, e, get<Recording>(simulation[sim_index]()),
            get<Only_Ch_Curent_Sub_Evolution>(simulation[sim_index]()), n_sub, sim_index);
        if (!ok || !ok.value()) {
            return ok.error()();
        }
    }
    return path_;
}

Maybe_error<std::string> write_csv(
    Experiment const& e,
    const var::Indexed<
        std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>>&
        simulation,
    std::size_t n_sub, std::string path) {
    return write_indexed_simulation_csv(
        simulation, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const auto& recordings) -> Maybe_error<bool> {
            for (std::size_t sim_index = 0; sim_index < recordings.size(); ++sim_index) {
                auto ok = emit_simulation_rows_with_sub(
                    writer, e, get<Recording>(recordings[sim_index]()),
                    get<Only_Ch_Curent_Sub_Evolution>(recordings[sim_index]()), n_sub, sim_index,
                    base_ctx);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        });
}


// ===========================================================================
// load_simulations: INDEX-AWARE sibling of simulate(). Reads a CSV written by
// write_csv(experiment, simulations, path) and reconstructs the SAME type that
// simulate() produces through DSL lifting: var::Indexed<vector<Simulated_
// Recording>>, carrying the same axes. This keeps the axis chain alive so
// downstream (calc_dlikelihood / calc_*_fisher / write_csv) stays Indexed and
// the axis columns propagate into every output CSV — symmetric to the
// index-aware WRITE that put all combos into one CSV with axis columns.
//
// Axis detection: any header column NOT in the fixed CSV schema is treated as
// an axis (the schema's standard 24 columns are known; the writer appends axis
// columns beyond them). Distinct labels per axis are collected in first-
// appearance order (matching the writer's flat-coordinate iteration).
//
// Recording reconstruction: rows with scope=="simulation" and non-empty
// (simulation_index, sample_index, value) are grouped by (coordinate,
// simulation_index); the `value` column (= patch_current) builds each
// Recording in sample_index order.
//
// `replica_indices` filters which simulation_index values to load (empty =>
// all). Errors if any coordinate of the reconstructed index space has no data
// (a missing axis index — the CSV is incomplete).
//
// The original SeedNumber is not persisted in the CSV; loaded records get
// SeedNumber=0 (downstream reads only the Recording).
// ===========================================================================
Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>>
load_simulations(std::string filename,
                 std::vector<std::size_t> replica_indices) {
    using SimVec = std::vector<Simulated_Recording<var::please_include<>>>;
    std::ifstream f(filename);
    if (!f.is_open())
        return error_message("cannot open ", filename);

    std::string header_line;
    if (!std::getline(f, header_line))
        return error_message("empty file ", filename);

    std::vector<std::string> headers;
    {
        std::stringstream ss(header_line);
        std::string col;
        while (std::getline(ss, col, ',')) headers.push_back(col);
    }
    auto find_col = [&](const std::string& name) -> Maybe_error<std::size_t> {
        auto it = std::find(headers.begin(), headers.end(), name);
        if (it == headers.end())
            return error_message("column '", name, "' not found in header of ", filename);
        return static_cast<std::size_t>(std::distance(headers.begin(), it));
    };

    auto maybe_scope = find_col("scope");
    auto maybe_sim = find_col("simulation_index");
    auto maybe_sample = find_col("sample_index");
    auto maybe_value = find_col("value");
    if (!maybe_scope) return maybe_scope.error();
    if (!maybe_sim) return maybe_sim.error();
    if (!maybe_sample) return maybe_sample.error();
    if (!maybe_value) return maybe_value.error();
    const auto scope_col = maybe_scope.value();
    const auto sim_col = maybe_sim.value();
    const auto sample_col = maybe_sample.value();
    const auto value_col = maybe_value.value();

    // Fixed CSV schema (anything else in the header is an axis column).
    static const std::set<std::string> fixed_cols = {
        "scope", "simulation_index", "sample_index", "segment_index", "sub_index",
        "n_step", "step_start", "step_end", "step_middle", "agonist", "patch_current",
        "component_path", "variable", "operation", "value_row", "value_col", "probit",
        "calculus", "statistic", "quantile_level", "param_index", "param_col",
        "param_name", "value"};
    std::vector<std::size_t> axis_cols;
    std::vector<std::string> axis_names;
    for (std::size_t c = 0; c < headers.size(); ++c) {
        if (fixed_cols.find(headers[c]) == fixed_cols.end()) {
            axis_cols.push_back(c);
            axis_names.push_back(headers[c]);
        }
    }
    const std::size_t n_axes = axis_cols.size();

    // Per-axis label → index, in first-appearance order.
    std::vector<std::map<std::string, std::size_t>> label_index(n_axes);
    std::vector<std::vector<std::string>> labels(n_axes);

    std::size_t max_col = std::max({scope_col, sim_col, sample_col, value_col});
    for (auto c : axis_cols) max_col = std::max(max_col, c);

    const bool load_all = replica_indices.empty();
    std::set<std::size_t> requested(replica_indices.begin(), replica_indices.end());

    // Raw rows kept for a second pass once axes are known: (axis label idxs,
    // sim_idx, sample_idx, value).
    struct Row { std::vector<std::size_t> coord; std::size_t sim, sample; double value; };
    std::vector<Row> rows;

    std::string line;
    while (std::getline(f, line)) {
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string col;
        while (std::getline(ss, col, ',')) cols.push_back(col);
        if (cols.size() <= max_col) continue;
        if (cols[scope_col] != "simulation") continue;
        if (cols[sim_col].empty() || cols[sample_col].empty() || cols[value_col].empty())
            continue;

        std::size_t sim_idx, sample_idx;
        double value;
        try {
            sim_idx = std::stoull(cols[sim_col]);
            sample_idx = std::stoull(cols[sample_col]);
            value = std::stod(cols[value_col]);
        } catch (std::exception const& e) {
            return error_message("parse error in ", filename, ": ", e.what(), " (line: ", line, ")");
        }
        if (!load_all && requested.find(sim_idx) == requested.end()) continue;

        Row row;
        row.coord.resize(n_axes);
        for (std::size_t a = 0; a < n_axes; ++a) {
            const std::string& lbl = cols[axis_cols[a]];
            auto it = label_index[a].find(lbl);
            std::size_t idx;
            if (it == label_index[a].end()) {
                idx = labels[a].size();
                label_index[a].emplace(lbl, idx);
                labels[a].push_back(lbl);
            } else {
                idx = it->second;
            }
            row.coord[a] = idx;
        }
        row.sim = sim_idx;
        row.sample = sample_idx;
        row.value = value;
        rows.push_back(std::move(row));
    }

    // Build the index space and per-axis sizes (for flat indexing).
    std::vector<var::Axis> out_axes;
    out_axes.reserve(n_axes);
    std::vector<std::size_t> axis_size(n_axes);
    for (std::size_t a = 0; a < n_axes; ++a) {
        axis_size[a] = labels[a].size();
        out_axes.emplace_back(var::AxisId(axis_names[a]), labels[a]);
    }
    var::IndexSpace space(std::move(out_axes));
    const std::size_t n_coords = space.size();  // product of axis sizes (1 if no axes)

    auto flat_of = [&](const std::vector<std::size_t>& coord) -> std::size_t {
        std::size_t flat = 0, mult = 1;
        for (std::size_t a = 0; a < n_axes; ++a) {
            flat += coord[a] * mult;
            mult *= axis_size[a];
        }
        return flat;
    };

    // data[flat][sim_idx][sample_idx] = value.
    std::vector<std::map<std::size_t, std::map<std::size_t, double>>> data(n_coords);
    for (auto& row : rows) {
        data[flat_of(row.coord)][row.sim][row.sample] = row.value;
    }

    // Build the Indexed values: one vector<Sim> per coordinate.
    std::vector<SimVec> values(n_coords);
    for (std::size_t flat = 0; flat < n_coords; ++flat) {
        if (data[flat].empty())
            return error_message("load_simulations: index coordinate (flat ", flat,
                                 ") has no data in ", filename,
                                 " — the CSV is missing an axis index");
        SimVec sims;
        sims.reserve(data[flat].size());
        for (auto& [sim_idx, samples] : data[flat]) {
            Recording rec;
            rec().reserve(samples.size());
            for (auto& [sample_idx, v] : samples) rec().emplace_back(v);
            sims.push_back(Simulated_Recording<var::please_include<>>{
                {SeedNumber{0}, std::move(rec)}});
        }
        values[flat] = std::move(sims);
    }

    return var::Indexed<SimVec>(std::move(space), std::move(values));
}

Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>>
load_simulations(std::string filename) {
    return load_simulations(std::move(filename), std::vector<std::size_t>{});
}

}  // namespace macrodr::cmd
