#ifndef CLI_THERMO_EVIDENCE_H
#define CLI_THERMO_EVIDENCE_H
#include "CLI_function_table.h"

#include "parallel_tempering.h"
#include <cstddef>

namespace macrodr {
namespace cmd {

inline auto
set_ThermoAlgorithm(std::size_t num_scouts_per_ensemble,
                    std::size_t number_trials_until_give_up, double stops_at,
                    double beta_upper_value, double beta_medium_value,
                    bool includes_zero, std::size_t max_iter_equilibrium,
                    std::string path, std::size_t beta_size,
                    std::size_t beta_upper_size, std::size_t beta_medium_size,
                    std::string filename, std::size_t thermo_jumps_every,
                    std::size_t max_num_simultaneous_temperatures) {
  using namespace macrodr;

  auto saving_itervals = Saving_intervals(
      Vector_Space(Save_Evidence_every(num_scouts_per_ensemble),
                   Save_Likelihood_every(num_scouts_per_ensemble),
                   Save_Parameter_every(num_scouts_per_ensemble),
                   Save_Predictions_every(num_scouts_per_ensemble * 20)));

  return std::tuple(
      path, filename, num_scouts_per_ensemble, number_trials_until_give_up,
      thermo_jumps_every, max_iter_equilibrium, beta_size, beta_upper_size,
      beta_medium_size, beta_upper_value, beta_medium_value, stops_at,
      includes_zero, saving_itervals, max_num_simultaneous_temperatures);
}

using thermo_algo_type =
    typename return_type<std::decay_t<decltype(&set_ThermoAlgorithm)>>::type;

inline void calc_thermo_evidence(std::string model, prior_value_type prior,
                                 likelihood_algo_type likelihood,
                                 std::string recording,
                                 experiment_type experiment,
                                 thermo_algo_type thermo_algorithm,
                                 tablefun_value_type ft, std::size_t myseed) {
  using namespace macrodr;
  auto [filename, num_scouts_per_ensemble] = std::move(ft);

  auto save_every = num_scouts_per_ensemble;
  auto ftbl3 = get_function_Table_maker_St(filename, save_every)();

  auto Maybe_model_v = get_model(model);

  if (Maybe_model_v) {
    auto model_v = std::move(Maybe_model_v.value());
    return std::visit(
        [&ftbl3, &experiment, &recording, &prior, &likelihood,
         &thermo_algorithm, &myseed](auto model0ptr) {
          auto &model0 = *model0ptr;
          myseed = calc_seed(myseed);
          mt_64i mt(myseed);

          auto [path, file_name, num_scouts_per_ensemble,
                number_trials_until_give_up, thermo_jump_factor,
                max_iter_equilibrium, beta_size, beta_upper_size,
                beta_medium_size, beta_upper_value, beta_medium_value, stops_at,
                includes_zero, saving_intervals,
                max_num_simultaneous_temperatures] =
              std::move(thermo_algorithm);

          auto [adaptive_aproximation, recursive_approximation,
                averaging_approximation, variance_correction,
                variance_correction_approximation, n_sub_dt] = likelihood;

          using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

          std::string ModelName = model0.model_name();

          auto Maybe_param1_prior = var::load_Prior<MyModel>(
              prior.first, prior.second, model0.model_name(), model0.names());
          if (!Maybe_param1_prior) {
            std::cerr << Maybe_param1_prior.error()();
          } else {
            auto param1_prior = std::move(Maybe_param1_prior.value());
            std::size_t thermo_jumps_every =
                param1_prior.size() * thermo_jump_factor;

            Recording y;
            auto Maybe_y = load_Recording_Data(recording, ",", y);
            if (Maybe_y) {

              std::string filename = file_name + "_" + ModelName + "_" +
                                     time_now() + "_" + std::to_string(myseed);

              auto tmi = new_thermo_Model_by_max_iter<MyModel>(
                  path, filename, num_scouts_per_ensemble,
                  max_num_simultaneous_temperatures, thermo_jumps_every,
                  max_iter_equilibrium, beta_size, beta_upper_size,
                  beta_medium_size, beta_upper_value, beta_medium_value,

                  stops_at, includes_zero, saving_intervals, myseed);

              auto modelLikelihood_v = Likelihood_Model_v{}.bool_op(
                  uses_adaptive_aproximation(adaptive_aproximation),
                  uses_recursive_aproximation(recursive_approximation),
                  uses_averaging_aproximation(averaging_approximation),
                  uses_variance_aproximation(variance_correction),
                  uses_variance_correction_aproximation(
                      variance_correction_approximation),
                  model0, Simulation_n_sub_dt(n_sub_dt));
              // using m2=typename
              // decltype(modelLikelihood_v)::paseModelLikelihoodv;

              std::visit(
                  [&ftbl3, &tmi, &param1_prior, &y,
                   &experiment](auto &modelLikelihood) {
                    auto opt =
                        thermo_evidence(ftbl3, std::move(tmi), param1_prior,
                                        modelLikelihood, y, experiment);
                  },
                  modelLikelihood_v);
            }
          }
        },
        model_v);
  }
}

} // namespace cmd
} // namespace macrodr

#endif // CLI_THERMO_EVIDENCE_H
