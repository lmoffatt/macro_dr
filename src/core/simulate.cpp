

#include <macrodr/cmd/simulate.h>
#include <macrodr/cmd/detail/write_csv_common.h>

#include <cstddef>
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

Maybe_error<bool> emit_simulation_rows(csv_detail::CsvWriter& writer, const Experiment& e,
                                       const Recording& recording,
                                       std::optional<std::size_t> simulation_index) {
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

            csv_detail::CsvContext ctx;
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
                                                std::optional<std::size_t> simulation_index) {
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

            csv_detail::CsvContext full_ctx;
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

                csv_detail::CsvContext sub_ctx;
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

}  // namespace

namespace macrodr::cmd {

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed) {
    myseed = calc_seed(myseed);

    mt_64i mt(myseed);
    Simulation_Parameters sim = Simulation_Parameters(Simulation_n_sub_dt(n_sub));
    return Macro_DMR{}.sample(mt, *model, par, e, sim, r);
   }

Maybe_error<Simulated_Recording<var::please_include<>>> run_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    const Experiment& e, const Recording& r, std::size_t n_sub, std::size_t myseed) {
    return run_simulations(model, par.to_value(), e, r, n_sub, myseed);
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_values& par,
    std::size_t n_simulations, const Experiment& e, const Recording& r, std::size_t n_sub,
    std::size_t myseed) {
    myseed = calc_seed(myseed);

    mt_64i mt(myseed);
    std::vector<Simulated_Recording<var::please_include<>>> result;
    result.reserve(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) {
        Simulation_Parameters sim = Simulation_Parameters(Simulation_n_sub_dt(n_sub));
        auto maybe_sim = Macro_DMR{}.sample(mt, *model, par, e, sim, r);
        if (!maybe_sim) {
            return maybe_sim.error();
        }
        result.push_back(maybe_sim.value());
    }
    return result;
}

Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> run_n_simulations(
    const ModelPtr& model, const var::Parameters_transformed& par,
    std::size_t n_simulations, const Experiment& e, const Recording& r, std::size_t n_sub,
    std::size_t myseed) {
    return run_n_simulations(model, par.to_value(), n_simulations, e, r, n_sub, myseed);
}

Maybe_error<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>
    run_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_values& par,const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    myseed = calc_seed(myseed);
    mt_64i mt(myseed);
     Simulation_Parameters sim = Simulation_Parameters(Simulation_n_sub_dt(n_sub));
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


Maybe_error<std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>>>
    run_n_simulations_with_sub_intervals(const ModelPtr& model,
                                         const var::Parameters_values& par,
                                         std::size_t n_simulations, const Experiment& e,
                                         const Recording& r, std::size_t n_sub,
                                         std::size_t myseed) {
    myseed = calc_seed(myseed);
    mt_64i mt(myseed);
    std::vector<Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> result;
    result.reserve(n_simulations);
    for (std::size_t i = 0; i < n_simulations; ++i) {
        Simulation_Parameters sim = Simulation_Parameters(Simulation_n_sub_dt(n_sub));
        auto maybe_sim = Macro_DMR{}.sample_sub_y(mt, *model, par, e, sim, r);
        if (!maybe_sim) {
            return maybe_sim.error();
        }
        result.push_back(maybe_sim.value());
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
    auto [includeN, n_sub_dt] = sim_algo_type;
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
            [&filename_prefix, &experiment, &recording, &myseed, &parameter_files, n_sub_dt,
             includeN](auto model0ptr) -> std::string {
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
                        auto sim = Macro_DMR{}.sample(
                            mt, model0, param1, experiment,
                            Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)), recording);
                        if (!sim)
                            return "";
                        else {
                            save_Recording(filename + "_simulation.csv", ",",
                                           get<Recording>(sim.value()()));
                            return filename + "_simulation.csv";
                        }
                    } else {
                        auto sim = Macro_DMR{}.sample_N(
                            mt, model0, param1, experiment,
                            Simulation_Parameters(Simulation_n_sub_dt(n_sub_dt)), recording);
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


}  // namespace macrodr::cmd
