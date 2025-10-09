#include <CLI_function_table.h>
#include <derivative_operator.h>
#include <macrodr/cmd/likelihood.h>
#include <parameters.h>

#include "macrodr/cmd/load_model.h"
#include "qmodel.h"

namespace macrodr::cmd {
auto calculate_likelihood(const interface::IModel<var::Parameters_values>& model0,
                          const var::Parameters_values& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto maybe_modelLikelihood =
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(model0)>(model0, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant();
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) {
            return logLikelihood(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood(const interface::IModel<var::Parameters_values>& model0,
                           const var::Parameters_values& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation) -> Maybe_error<dlogLs> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto par_transformed = par.to_transformed();

    auto dmodel = load_dmodel(model0.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood =
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant();
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par_transformed, &e, &r](auto& modelLikelihood) {
            return dlogLikelihood(ftbl3, modelLikelihood, par_transformed, r, e);
        },
        modelLikelihood_v);
}

auto calculate_diff_likelihood(const interface::IModel<var::Parameters_values>& model0,
                               const var::Parameters_values& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<dlogLs> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto par_transformed = par.to_transformed();

    auto dmodel = load_dmodel(model0.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood =
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant();
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par_transformed, &e, &r, delta_param](auto& modelLikelihood) {
            return diff_logLikelihood(ftbl3, modelLikelihood, par_transformed, r, e, delta_param);
        },
        modelLikelihood_v);
}

auto calculate_likelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_values& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Patch_State_Evolution> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto maybe_modelLikelihood =
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(model0)>(model0, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant();
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) {
            return logLikelihoodPredictions(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                       const var::Parameters_values& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto par_transformed = par.to_transformed();

    auto dmodel = load_dmodel(model0.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood =
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant();
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par_transformed, &e, &r](auto& modelLikelihood) {
            return dlogLikelihoodPredictions(ftbl3, modelLikelihood, par_transformed, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood_predictions_model(const std::string& model_name,
                                             const var::Parameters_values& par, const Experiment& e,
                                             const Recording& r, bool adaptive_approximation,
                                             bool recursive_approximation,
                                             int averaging_approximation,
                                             bool variance_approximation,
                                             bool taylor_variance_correction_approximation)
    -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto par_transformed = par.to_transformed();

    auto dmodel = get_model(model_name);
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = dmodel.value();
    return std::visit(
        [&](auto m_ptr)
            -> Maybe_error<var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
            auto& m = *m_ptr;
            auto maybe_modelLikelihood =

                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              false>,
                    decltype(m)>(m, nsub, uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
                    .get_variant();
            if (!maybe_modelLikelihood) {
                return maybe_modelLikelihood.error();
            }
            auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

            return std::visit(
                [&ftbl3, &par_transformed, &e, &r](auto& modelLikelihood)
                    -> Maybe_error<
                        var::Derivative<Patch_State_Evolution, var::Parameters_transformed>> {
                    return dlogLikelihoodPredictions(ftbl3, modelLikelihood, par_transformed, r, e);
                },
                modelLikelihood_v);
        },
        model0_d);
};

}  // namespace macrodr::cmd
