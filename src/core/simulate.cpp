

#include <macrodr/cmd/simulate.h>

#include "CLI_macro_dr_base.h"
//#include "cuevi.h"
#include "function_memoization.h"
#include "maybe_error.h"
#include "models_MoffattHume_allosteric.h"
#include "parallel_tempering.h"
#include "parameters.h"
//#include "parameters_derivative.h"
#include "qmodel.h"

namespace macrodr {
namespace cmd {

inline std::string run_simulation(std::string filename_prefix, recording_type recording_file,
                                  experiment_type experiment, std::size_t myseed,
                                  std::string modelName, parameters_value_type parameter_files,
                                  simulation_algo_type sim_algo_type) {
    auto [includeN, n_sub_dt] = sim_algo_type;
    auto Maybe_recording = load_Recording(recording_file, std::string(","));
    if (!Maybe_recording)
        return Maybe_recording.error()();
    else {
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
                    using MyModel = typename std::decay_t<decltype(model0)>::my_Id;
                    auto Maybe_parameter_values =
                        var::load_Parameters(parameter_files.first, parameter_files.second,
                                             model0.model_name(), model0.names());
                    std::string ModelName = model0.model_name();

                    std::string filename = filename_prefix + "_" + ModelName + "_" + time_now() +
                                           "_" + std::to_string(myseed);

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
}

}  // namespace cmd
}  // namespace macrodr
