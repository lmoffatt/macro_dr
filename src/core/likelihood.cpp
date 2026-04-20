#include <CLI_function_table.h>
#include <bootstrap.h>
#include <derivative_fwd.h>
#include <derivative_operator.h>
#include <distributions.h>
#include <macrodr/cmd/likelihood.h>
#include <macrodr/dsl/type_name.h>
#include <matrix.h>
#include <moment_statistics.h>
#include <parameters.h>
#include <random_samplers.h>
#include <variables.h>

#include <concepts>
#include <cstddef>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "macrodr/cmd/load_model.h"
#include "qmodel.h"

namespace macrodr::cmd {

namespace {

template <class adaptive, class recursive, class averaging, class variance, class taylor,
          class Model, class FuncTable>
auto calculate_mdlikelihood_impl(
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor, Model>& lik,
    FuncTable& ftbl3, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r) -> Maybe_error<dMacro_State_Hessian_minimal> {
    auto dmodel = load_dmodel(lik.m.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto dlikelihood = Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor,
                                                  decltype(*model0_d)>(*model0_d, lik.n_sub_dt);
    return dlogLikelihood(ftbl3, dlikelihood, par, r, e);
}

template <class adaptive, class recursive, class averaging, class variance, class taylor,
          class Model, class FuncTable>
auto calculate_mdiff_likelihood_impl(
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor, Model>& lik,
    FuncTable& ftbl3, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, double delta_param) -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    auto dmodel = load_dmodel(lik.m.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto dlikelihood = Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor,
                                                  decltype(*model0_d)>(*model0_d, lik.n_sub_dt);
    return diff_logLikelihood(ftbl3, dlikelihood, par, r, e, delta_param);
}

template <class adaptive, class recursive, class averaging, class variance, class taylor,
          class Model, class FuncTable>
auto calculate_mdlikelihood_predictions_impl(
    const Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor, Model>& lik,
    FuncTable& ftbl3, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto dmodel = load_dmodel(lik.m.model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto dlikelihood = Likelihood_Model_constexpr<adaptive, recursive, averaging, variance, taylor,
                                                  decltype(*model0_d)>(*model0_d, lik.n_sub_dt);
    return dlogLikelihoodPredictions(ftbl3, dlikelihood, par, r, e);
}

}  // namespace

auto calculate_mlikelihood(const likelihood_algorithm_type& modelLikelihood_v,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r) -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto par_values = par.to_value();

    return std::visit(
        [&](const auto& modelLikelihood) {
            return logLikelihood(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_mdlikelihood(const likelihood_algorithm_type& modelLikelihood_v,
                            const var::Parameters_transformed& par, const Experiment& e,
                            const Recording& r) -> Maybe_error<dMacro_State_Hessian_minimal> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    return std::visit(
        [&](const auto& modelLikelihood) -> Maybe_error<dMacro_State_Hessian_minimal> {
            return calculate_mdlikelihood_impl(modelLikelihood, ftbl3, par, e, r);
        },
        modelLikelihood_v);
}

auto calculate_mdiff_likelihood(const likelihood_algorithm_type& modelLikelihood_v,
                                const var::Parameters_transformed& par, const Experiment& e,
                                const Recording& r, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    return std::visit(
        [&](const auto& modelLikelihood) -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
            return calculate_mdiff_likelihood_impl(modelLikelihood, ftbl3, par, e, r, delta_param);
        },
        modelLikelihood_v);
}

auto calculate_mlikelihood_predictions(const likelihood_algorithm_type& modelLikelihood_v,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_predictions> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto par_values = par.to_value();

    return std::visit(
        [&](const auto& modelLikelihood) {
            return logLikelihoodPredictions(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_mlikelihood_diagnostics(const likelihood_algorithm_type& modelLikelihood_v,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto par_values = par.to_value();

    return std::visit(
        [&](const auto& modelLikelihood) {
            return logLikelihoodDiagnostic(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_mdlikelihood_predictions(const likelihood_algorithm_type& modelLikelihood_v,
                                        const var::Parameters_transformed& par, const Experiment& e,
                                        const Recording& r)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    return std::visit(
        [&](const auto& modelLikelihood) -> Maybe_error<dMacro_State_Ev_gradient_all> {
            return calculate_mdlikelihood_predictions_impl(modelLikelihood, ftbl3, par, e, r);
        },
        modelLikelihood_v);
}

auto calculate_likelihood(const ModelPtr& model0,
                          const var::Parameters_transformed& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    const auto& model_ref = *model0;
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto maybe_modelLikelihood =
        merge_Maybe_variant(Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant(),
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation,1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant());
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihood(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood(const ModelPtr& model0,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto dmodel = load_dmodel(model0->model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood = merge_Maybe_variant(
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant(), 
               Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation,  true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant());
    ;
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) -> Maybe_error<dMacro_State_Hessian_minimal> {
            return dlogLikelihood(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_diff_likelihood(const ModelPtr& model0,
                               const var::Parameters_transformed& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto dmodel = load_dmodel(model0->model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood =merge_Maybe_variant(
         Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
             var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
             var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
             var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
             var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
             decltype(*model0_d)>(*model0_d, nsub,
                                  uses_adaptive_aproximation_value(adaptive_approximation),
                                  uses_recursive_aproximation_value(recursive_approximation),
                                  uses_averaging_aproximation_value(averaging_approximation),
                                  uses_variance_aproximation_value(variance_approximation),
                                  uses_taylor_variance_correction_aproximation_value(
                                      taylor_variance_correction_approximation))
            .get_variant(),
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation,  true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation,  true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,  true>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant());
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par, &e, &r, delta_param](auto& modelLikelihood) {
            return diff_logLikelihood(ftbl3, modelLikelihood, par, r, e, delta_param);
        },
        modelLikelihood_v);
}

auto calculate_likelihood_predictions(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions> {
    const auto& model_ref = *model0;
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto maybe_modelLikelihood =merge_Maybe_variant(
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant(),
            Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation,  1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation,  true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant());
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihoodPredictions(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_likelihood_diagnostics(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    const auto& model_ref = *model0;
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto maybe_modelLikelihood =merge_Maybe_variant(
        Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false >,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant(),
            Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation,  true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation,  1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation,  true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
            decltype(model_ref)>(model_ref, nsub,
                              uses_adaptive_aproximation_value(adaptive_approximation),
                              uses_recursive_aproximation_value(recursive_approximation),
                              uses_averaging_aproximation_value(averaging_approximation),
                              uses_variance_aproximation_value(variance_approximation),
                              uses_taylor_variance_correction_aproximation_value(
                                  taylor_variance_correction_approximation))
            .get_variant());
    if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihoodDiagnostic(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood_predictions(const ModelPtr& model0,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto dmodel = load_dmodel(model0->model_name());
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = std::move(dmodel.value());
    auto maybe_modelLikelihood =merge_Maybe_variant(
         Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
             var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
             var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
             var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
             var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
             decltype(*model0_d)>(*model0_d, nsub,
                                  uses_adaptive_aproximation_value(adaptive_approximation),
                                  uses_recursive_aproximation_value(recursive_approximation),
                                  uses_averaging_aproximation_value(averaging_approximation),
                                  uses_variance_aproximation_value(variance_approximation),
                                  uses_taylor_variance_correction_aproximation_value(
                                      taylor_variance_correction_approximation))
            .get_variant(), 
               Likelihood_Model_regular<
            var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
            var::constexpr_Var_domain<bool, uses_recursive_aproximation,  true>,
            var::constexpr_Var_domain<int, uses_averaging_aproximation, 1, 2>,
            var::constexpr_Var_domain<bool, uses_variance_aproximation,  true>,
            var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
            decltype(*model0_d)>(*model0_d, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
            .get_variant());
       if (!maybe_modelLikelihood) {
        return maybe_modelLikelihood.error();
    }
    auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

    return std::visit(
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) -> Maybe_error<dMacro_State_Ev_gradient_all> {
            return dlogLikelihoodPredictions(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood_predictions_model(
    const std::string& model_name, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto dmodel = get_model(model_name);
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = dmodel.value();
    return std::visit(
        [&](auto m_ptr) -> Maybe_error<dMacro_State_Ev_gradient_all> {
            auto& m = *m_ptr;
            auto maybe_modelLikelihood = merge_Maybe_variant(
                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, false, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 0, 1, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, false, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              false>,
                    decltype(m)>(m, nsub, uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
                    .get_variant(),
                Likelihood_Model_regular<
                    var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
                    var::constexpr_Var_domain<bool, uses_recursive_aproximation, true>,
                    var::constexpr_Var_domain<int, uses_averaging_aproximation, 1, 2>,
                    var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
                    var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation,
                                              true>,
                    decltype(m)>(m, nsub, uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
                    .get_variant());
            if (!maybe_modelLikelihood) {
                return maybe_modelLikelihood.error();
            }
            auto modelLikelihood_v = std::move(maybe_modelLikelihood.value());

            return std::visit(
                [&ftbl3, &par, &e,
                 &r](auto& modelLikelihood) -> Maybe_error<dMacro_State_Ev_gradient_all> {
                    return dlogLikelihoodPredictions(ftbl3, modelLikelihood, par, r, e);
                },
                modelLikelihood_v);
        },
        model0_d);
};

template <class VS>
struct vector_space_types;
template <class... Cs>
struct vector_space_types<var::Vector_Space<Cs...>> {
    using types = std::tuple<std::type_identity<Cs>...>;
};

template <class T>
struct state_component_types;

template <typename... Vars>
struct state_component_types<Macro_State<Vars...>> {
    using types = std::tuple<std::type_identity<logL>, std::type_identity<Patch_State>,
                             std::type_identity<Vars>...>;
};

template <typename... Vars>
struct state_component_types<dMacro_State<Vars...>> {
    using types =
        std::tuple<std::type_identity<var::Derivative<logL, var::Parameters_transformed>>,
                   std::type_identity<var::Derivative<Patch_State, var::Parameters_transformed>>,
                   std::type_identity<Moment_statistics<Hessian_minus_CovGrad>>,
                   std::type_identity<covariance<Grad>>, std::type_identity<Hessian>,
                   std::type_identity<Vars>...>;
};

template <class T>
double primitive_to_double(const T& x) {
    if constexpr (requires { x.primitive()(); }) {
        return x.primitive()();
    } else if constexpr (requires { x.primitive(); }) {
        return x.primitive();
    } else {
        return x();
    }
}

template <class Lik>
std::vector<std::string> get_param_names_if_any(const Lik& lik) {
    std::vector<std::string> param_names;
    // Prefer extracting dx() from a concrete derivative component (e.g. get<logL>(dMacro_State)),
    // since the container itself might not be a Derivative<> even when it stores derivatives.
    if constexpr (macrodr::has_var_c<Lik const&, logL> &&
                  requires { var::get_dx_of_dfdx(get<logL>(lik)).parameters().names(); }) {
        const auto& names = var::get_dx_of_dfdx(get<logL>(lik)).parameters().names();
        param_names.assign(names.begin(), names.end());
    } else if constexpr (requires { var::get_dx_of_dfdx(lik).parameters().names(); }) {
        const auto& names = var::get_dx_of_dfdx(lik).parameters().names();
        param_names.assign(names.begin(), names.end());
    } else if constexpr (macrodr::has_var_c<Lik const&, Grad> &&
                         requires { var::parameter_names(get<Grad>(lik)()); }) {
        param_names = var::parameter_names(get<Grad>(lik)());
    } else if constexpr (macrodr::has_var_c<Lik const&, FIM> &&
                         requires { var::parameter_names(get<FIM>(lik)()); }) {
        param_names = var::parameter_names(get<FIM>(lik)());
    } else if constexpr (requires { var::parameter_names(lik); }) {
        param_names = var::parameter_names(lik);
    }
    return param_names;
}

inline std::string append_index(std::string base, std::size_t i) {
    return base + "_" + std::to_string(i);
}

inline std::string append_index(std::string base, std::size_t i, std::size_t j) {
    return base + "_" + std::to_string(i) + "_" + std::to_string(j);
}

template <class Writer, class T>
Maybe_error<bool> emit_component_rows(Writer& w, std::string label, const T& x);

template <class Writer, class T>
Maybe_error<bool> emit_value_rows(Writer& w, std::string label, const T& x) {
    if constexpr (is_of_this_template_type_v<std::decay_t<T>, var::Vector_Space>) {
        using Tuple = typename vector_space_types<std::decay_t<T>>::types;
        return std::apply(
            [&](auto... type_tags) -> Maybe_error<bool> {
                Maybe_error<bool> ok = true;
                ((ok =
                      ok ? emit_component_rows(
                               w,
                               macrodr::dsl::type_name_no_namespace<
                                   var::untransformed_type_t<typename decltype(type_tags)::type>>(),
                               get<typename decltype(type_tags)::type>(x))
                         : ok),
                 ...);
                return ok;
            },
            Tuple{});
    } else if constexpr (requires { x(); }) {
        return emit_value_rows(w, std::move(label), x());
    } else if constexpr (std::convertible_to<T, double>) {
        w.write_value(std::move(label), static_cast<double>(x));
        return true;
    } else if constexpr (requires {
                             x.nrows();
                             x.ncols();
                             x(0ul, 0ul);
                         }) {
        if (x.nrows() == 1 || x.ncols() == 1) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                w.write_value(append_index(label, i), static_cast<double>(x[i]));
            }
            return true;
        }
        for (std::size_t r = 0; r < x.nrows(); ++r) {
            for (std::size_t c = 0; c < x.ncols(); ++c) {
                w.write_value(append_index(label, r, c), static_cast<double>(x(r, c)));
            }
        }
        return true;
    } else if constexpr (requires {
                             x.size();
                             x[0];
                         }) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            w.write_value(append_index(label, i), static_cast<double>(x[i]));
        }
        return true;
    } else {
        return error_message("write_csv_rows: unsupported value type ", label, " (",
                             macrodr::dsl::type_name<T>(), ")");
    }
}

template <class Writer, class Der>
Maybe_error<bool> emit_derivative_rows(Writer& w, std::string label, const Der& d) {
    if constexpr (requires { var::inside_out(d); }) {
        auto out = var::inside_out(d);
        if constexpr (requires {
                          out.nrows();
                          out.ncols();
                          out(0ul, 0ul);
                      }) {
            for (std::size_t r = 0; r < out.nrows(); ++r) {
                for (std::size_t c = 0; c < out.ncols(); ++c) {
                    auto maybe_ok = emit_derivative_rows(w, append_index(label, r, c), out(r, c));
                    if (!maybe_ok || !maybe_ok.value()) {
                        return maybe_ok;
                    }
                }
            }
            return true;
        } else if constexpr (requires {
                                 out.size();
                                 out[0];
                             }) {
            for (std::size_t i = 0; i < out.size(); ++i) {
                auto maybe_ok = emit_derivative_rows(w, append_index(label, i), out[i]);
                if (!maybe_ok || !maybe_ok.value()) {
                    return maybe_ok;
                }
            }
            return true;
        } else {
            return error_message("write_csv_rows: unsupported inside_out container for ", label,
                                 " (", macrodr::dsl::type_name<decltype(out)>(), ")");
        }
    } else if constexpr (requires { d(); }) {
        return emit_derivative_rows(w, std::move(label), d());
    } else if constexpr (var::is_derivative_v<Der> &&
                         is_of_this_template_type_v<var::untransformed_type_t<std::decay_t<Der>>,
                                                    var::Vector_Space>) {
        using VS = var::untransformed_type_t<std::decay_t<Der>>;
        using Tuple = typename vector_space_types<VS>::types;
        return std::apply(
            [&](auto... type_tags) -> Maybe_error<bool> {
                Maybe_error<bool> ok = true;
                ((ok =
                      ok ? emit_component_rows(
                               w,
                               macrodr::dsl::type_name_no_namespace<
                                   var::untransformed_type_t<typename decltype(type_tags)::type>>(),
                               get<typename decltype(type_tags)::type>(d))
                         : ok),
                 ...);
                return ok;
            },
            Tuple{});
    } else if constexpr (requires { d.derivative()(); }) {
        const double prim = primitive_to_double(d);
        // w.write_value(label, prim);

        const auto& m = d.derivative()();
        for (std::size_t r = 0; r < m.nrows(); ++r) {
            for (std::size_t c = 0; c < m.ncols(); ++c) {
                w.write_derivative(label, r, c, prim, static_cast<double>(m(r, c)));
            }
        }
        return true;
    } else {
        return error_message("write_csv_rows: unsupported derivative type ", label, " (",
                             macrodr::dsl::type_name<Der>(), ")");
    }
}

template <class Writer, class T>
Maybe_error<bool> emit_component_rows(Writer& w, std::string label, const T& x) {
    if constexpr (var::is_derivative_v<T>) {
        return emit_derivative_rows(w, std::move(label), x);
    } else {
        return emit_value_rows(w, std::move(label), x);
    }
}

template <class Writer, class Lik>
Maybe_error<bool> emit_state_rows(Writer& w, Lik const& lik) {
    using Types = typename state_component_types<std::decay_t<Lik>>::types;
    auto handle_component = [&](auto type_tag) -> Maybe_error<bool> {
        using Id = typename decltype(type_tag)::type;
        if constexpr (is_of_this_template_type_v<Id, Evolution_of>) {
            return true;
        } else if constexpr (var::gets_it_c<Lik const&, Id>) {
            return emit_component_rows(
                w, macrodr::dsl::type_name_no_namespace<var::untransformed_type_t<Id>>(),
                get<Id>(lik));
        } else {
            return true;
        }
    };

    Maybe_error<bool> ok = true;
    std::apply([&](auto... type_tags) { ((ok = ok ? handle_component(type_tags) : ok), ...); },
               Types{});
    return ok;
}

struct SampleWriter {
    std::ostream& f;
    std::vector<std::string> const& param_names;
    double n_step;
    std::size_t sub_step;
    double step_start;
    double step_end;
    double step_mid;
    double agonist;
    double patch;

    void write_value(std::string component, double primitive) {
        f << "value," << n_step << "," << sub_step << "," << step_start << "," << step_end << ","
          << step_mid << "," << agonist << "," << patch << "," << component << ",,,," << primitive
          << ",\n";
    }

    void write_derivative(const std::string& component, std::size_t r, std::size_t c, double prim,
                          double deriv_value) {
        const auto& pname = (r < param_names.size()) ? param_names[r] : std::string{};
        f << "derivative," << n_step << "," << sub_step << "," << step_start << "," << step_end
          << "," << step_mid << "," << agonist << "," << patch << "," << component << "," << r
          << "," << c << "," << pname << "," << prim << "," << deriv_value << "\n";
    }
};

struct SampleWriter_i {
    std::ostream& f;
    std::vector<std::string> const& param_names;
    std::size_t n_simulation;
    double n_step;
    std::size_t sub_step;
    double step_start;
    double step_end;
    double step_mid;
    double agonist;
    double patch;

    void write_value(std::string component, double primitive) {
        f << "value," << n_simulation << "," << n_step << "," << sub_step << "," << step_start
          << "," << step_end << "," << step_mid << "," << agonist << "," << patch << ","
          << component << ",,,," << primitive << ",\n";
    }

    void write_derivative(const std::string& component, std::size_t r, std::size_t c, double prim,
                          double deriv_value) {
        const auto& pname = (r < param_names.size()) ? param_names[r] : std::string{};
        f << "derivative," << n_simulation << "," << n_step << "," << sub_step << "," << step_start
          << "," << step_end << "," << step_mid << "," << agonist << "," << patch << ","
          << component << "," << r << "," << c << "," << pname << "," << prim << "," << deriv_value
          << "\n";
    }
};

struct StateWriter {
    std::ostream& f;
    std::vector<std::string> const& param_names;

    void write_value(std::string component, double primitive) {
        f << "value,,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << ",,,," << primitive << ",\n";
    }

    void write_derivative(const std::string& component, std::size_t r, std::size_t c, double prim,
                          double deriv_value) {
        const auto& pname = (r < param_names.size()) ? param_names[r] : std::string{};
        f << "derivative,,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << "," << r << "," << c << "," << pname << "," << prim << ","
          << deriv_value << "\n";
    }
};

struct StateWriter_scoped {
    std::ostream& f;
    std::vector<std::string> const& param_names;

    void write_value(std::string component, double primitive) {
        f << "state_value,,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << ",,,," << primitive << ",\n";
    }

    void write_derivative(const std::string& component, std::size_t r, std::size_t c, double prim,
                          double deriv_value) {
        const auto& pname = (r < param_names.size()) ? param_names[r] : std::string{};
        f << "state_derivative,,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << "," << r << "," << c << "," << pname << "," << prim << ","
          << deriv_value << "\n";
    }
};

struct StateWriter_i {
    std::ostream& f;
    std::vector<std::string> const& param_names;
    std::size_t n_simulation;

    void write_value(std::string component, double primitive) {
        f << "state_value," << n_simulation
          << ",,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << ",,,," << primitive << ",\n";
    }

    void write_derivative(const std::string& component, std::size_t r, std::size_t c, double prim,
                          double deriv_value) {
        const auto& pname = (r < param_names.size()) ? param_names[r] : std::string{};
        f << "state_derivative," << n_simulation
          << ",,,,,,,"  // n_step,sub_step,step_start,step_end,step_middle,agonist,patch
          << "," << component << "," << r << "," << c << "," << pname << "," << prim << ","
          << deriv_value << "\n";
    }
};

template <class Writer, class Lik>
Maybe_error<bool> emit_state_members(Writer& writer, detail::CsvContext ctx, const Lik& lik) {
    using Types = typename state_component_types<std::decay_t<Lik>>::types;
    Maybe_error<bool> ok = true;
    []<std::size_t... Is>(Writer& w, detail::CsvContext base_ctx, const Lik& value,
                          Maybe_error<bool>& ok_ref, std::index_sequence<Is...>) {
        (([&] {
            if (!ok_ref || !ok_ref.value()) {
                return;
            }
            using Id = typename std::tuple_element_t<Is, Types>::type;
            if constexpr (!is_of_this_template_type_v<Id, Evolution_of>) {
                if constexpr (var::gets_it_c<Lik const&, Id>) {
                    ok_ref = detail::emit_named_component(
                        w, base_ctx, detail::component_label<Id>(), get<Id>(value));
                }
            }
        }()),
         ...);
    }(writer, std::move(ctx), lik, ok, std::make_index_sequence<std::tuple_size_v<Types>>{});
    return ok;
}

template <class Writer, class Lik>
Maybe_error<bool> emit_state_rows_without_experiment(Writer& writer,
                                                     std::optional<std::size_t> simulation_index,
                                                     const Lik& lik,
                                                     const detail::CsvContext& base_ctx = {}) {
    auto state_ctx = base_ctx;
    state_ctx.scope = "state";
    state_ctx.simulation_index = simulation_index;

    auto ok = emit_state_members(writer, state_ctx, lik);
    if (!ok || !ok.value()) {
        return ok;
    }

    if constexpr (macrodr::has_var_c<Lik const&, Evolution>) {
        const auto& evolution = get<Evolution>(lik)();
        for (std::size_t i = 0; i < evolution.size(); ++i) {
            auto evo_ctx = base_ctx;
            evo_ctx.scope = "evolution";
            evo_ctx.simulation_index = simulation_index;
            evo_ctx.sample_index = i;
            ok = detail::emit_any(writer, std::move(evo_ctx), evolution[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
        }
    }

    return true;
}

template <class Writer, class Lik>
Maybe_error<bool> emit_state_rows_with_experiment(Writer& writer, const Experiment& e,
                                                  const Recording& recording,
                                                  std::optional<std::size_t> simulation_index,
                                                  const Lik& lik,
                                                  const detail::CsvContext& base_ctx = {}) {
    const auto& conditions = get<Recording_conditions>(e);
    if (conditions().size() != recording().size()) {
        return error_message("Experiment samples ", conditions().size(),
                             " differ from Recording samples ", recording().size());
    }

    auto ok = emit_state_rows_without_experiment(writer, simulation_index, lik, base_ctx);
    if (!ok || !ok.value()) {
        return ok;
    }

    const auto& evolution = get<Evolution>(lik)();
    if (evolution.size() != conditions().size()) {
        return error_message("Evolution samples ", evolution.size(),
                             " differ from Recording samples ", conditions().size());
    }

    const double fs = get<Frequency_of_Sampling>(e)();
    for (std::size_t i = 0; i < conditions().size(); ++i) {
        double step_start = get<Time>(conditions()[i])();
        const auto& segments = get<Agonist_evolution>(conditions()[i])();
        for (std::size_t j = 0; j < segments.size(); ++j) {
            const double duration = get<number_of_samples>(segments[j])() / fs;
            const double step_end = step_start + duration;

            auto evo_ctx = base_ctx;
            evo_ctx.scope = "evolution";
            evo_ctx.simulation_index = simulation_index;
            evo_ctx.sample_index = i;
            evo_ctx.segment_index = j;
            evo_ctx.n_step = static_cast<double>(i) +
                             static_cast<double>(j) / static_cast<double>(segments.size());
            evo_ctx.time_start = step_start;
            evo_ctx.time_end = step_end;
            evo_ctx.time_middle = 0.5 * (step_start + step_end);
            evo_ctx.agonist = get<Agonist_concentration>(segments[j])();
            evo_ctx.patch_current = recording()[i]();

            ok = detail::emit_any(writer, std::move(evo_ctx), evolution[i]);
            if (!ok || !ok.value()) {
                return ok;
            }
            step_start = step_end;
        }
    }

    return true;
}

template <class Lik>
std::vector<std::string> first_non_empty_param_names(const std::vector<Lik>& liks) {
    std::vector<std::string> param_names;
    for (const auto& lik : liks) {
        param_names = detail::get_param_names_if_any(lik);
        if (!param_names.empty()) {
            break;
        }
    }
    return param_names;
}

template <class T>
Maybe_error<std::reference_wrapper<const T>> indexed_front(const var::Indexed<T>& indexed) {
    auto maybe_coord = indexed.begin();
    if (!maybe_coord) {
        return maybe_coord.error();
    }
    auto maybe_value = indexed.at(maybe_coord.value());
    if (!maybe_value) {
        return maybe_value.error();
    }
    return maybe_value.value();
}

// (1) Experiment + Simulation + State with per-sample Evolution
template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path) {
    const auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, detail::get_param_names_if_any(lik));
    auto ok =
        emit_state_rows_with_experiment(writer, e, get<Recording>(simulation()), std::nullopt, lik);
    if (!ok || !ok.value()) {
        return ok.error()();
    }
    return path_;
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path) {
    auto valid = detail::validate_indexed_write_csv_value(simulation, "simulation");
    if (!valid) {
        return valid.error();
    }

    return detail::write_indexed_rows_csv(
        simulation.index_space(), detail::get_param_names_if_any(lik), std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_sim = simulation.at(coord);
            if (!maybe_sim) {
                return maybe_sim.error();
            }
            return emit_state_rows_with_experiment(writer, e,
                                                   get<Recording>(maybe_sim.value().get()()),
                                                   std::nullopt, lik, base_ctx);
        });
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path) {
    auto valid = detail::validate_indexed_write_csv_value(lik, "likelihood");
    if (!valid) {
        return valid.error();
    }
    auto maybe_first = indexed_front(lik);
    if (!maybe_first) {
        return maybe_first.error();
    }

    return detail::write_indexed_rows_csv(
        lik.index_space(), detail::get_param_names_if_any(maybe_first.value().get()),
        std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_lik = lik.at(coord);
            if (!maybe_lik) {
                return maybe_lik.error();
            }
            return emit_state_rows_with_experiment(writer, e, get<Recording>(simulation()),
                                                   std::nullopt, maybe_lik.value().get(),
                                                   base_ctx);
        });
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path) {
    auto sim_valid = detail::validate_indexed_write_csv_value(simulation, "simulation");
    if (!sim_valid) {
        return sim_valid.error();
    }
    auto lik_valid = detail::validate_indexed_write_csv_value(lik, "likelihood");
    if (!lik_valid) {
        return lik_valid.error();
    }
    auto maybe_space = detail::merge_indexed_write_csv_spaces(
        simulation.index_space(), "simulation", lik.index_space(), "likelihood");
    if (!maybe_space) {
        return maybe_space.error();
    }
    auto maybe_first = indexed_front(lik);
    if (!maybe_first) {
        return maybe_first.error();
    }

    return detail::write_indexed_rows_csv(
        maybe_space.value(), detail::get_param_names_if_any(maybe_first.value().get()),
        std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_sim = simulation.at(coord);
            if (!maybe_sim) {
                return maybe_sim.error();
            }
            auto maybe_lik = lik.at(coord);
            if (!maybe_lik) {
                return maybe_lik.error();
            }
            return emit_state_rows_with_experiment(writer, e,
                                                   get<Recording>(maybe_sim.value().get()()),
                                                   std::nullopt, maybe_lik.value().get(),
                                                   base_ctx);
        });
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   std::vector<Simulated_Recording<SimTag>> const& simulation,
                                   std::vector<TMacro_State<vVars...>> const& liks,
                                   std::string path) {
    auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }
    if (simulation.size() != liks.size()) {
        return error_message("number of Simulated_Recordings ", simulation.size(),
                             " differ from number of Likelihood states ", liks.size());
    }

    auto param_names = first_non_empty_param_names(liks);

    detail::CsvWriter writer(f, param_names);
    for (std::size_t sim_i = 0; sim_i < simulation.size(); ++sim_i) {
        auto ok = emit_state_rows_with_experiment(writer, e, get<Recording>(simulation[sim_i]()),
                                                  sim_i, liks[sim_i]);
        if (!ok || !ok.value()) {
            return ok.error()();
        }
    }
    return path_;
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    std::vector<TMacro_State<vVars...>> const& liks, std::string path) {
    auto valid = detail::validate_indexed_write_csv_value(simulation, "simulations");
    if (!valid) {
        return valid.error();
    }
    auto param_names = first_non_empty_param_names(liks);

    return detail::write_indexed_rows_csv(
        simulation.index_space(), param_names, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_sim = simulation.at(coord);
            if (!maybe_sim) {
                return maybe_sim.error();
            }
            const auto& batch = maybe_sim.value().get();
            if (batch.size() != liks.size()) {
                return error_message("number of Simulated_Recordings ", batch.size(),
                                     " differ from number of Likelihood states ", liks.size());
            }
            for (std::size_t sim_i = 0; sim_i < batch.size(); ++sim_i) {
                auto ok = emit_state_rows_with_experiment(writer, e, get<Recording>(batch[sim_i]()),
                                                          sim_i, liks[sim_i], base_ctx);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        });
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path) {
    auto valid = detail::validate_indexed_write_csv_value(liks, "likelihood");
    if (!valid) {
        return valid.error();
    }
    auto maybe_first = indexed_front(liks);
    if (!maybe_first) {
        return maybe_first.error();
    }
    auto param_names = first_non_empty_param_names(maybe_first.value().get());

    return detail::write_indexed_rows_csv(
        liks.index_space(), param_names, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_lik = liks.at(coord);
            if (!maybe_lik) {
                return maybe_lik.error();
            }
            const auto& batch = maybe_lik.value().get();
            if (simulation.size() != batch.size()) {
                return error_message("number of Simulated_Recordings ", simulation.size(),
                                     " differ from number of Likelihood states ", batch.size());
            }
            for (std::size_t sim_i = 0; sim_i < simulation.size(); ++sim_i) {
                auto ok =
                    emit_state_rows_with_experiment(writer, e, get<Recording>(simulation[sim_i]()),
                                                    sim_i, batch[sim_i], base_ctx);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        });
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path) {
    auto sim_valid = detail::validate_indexed_write_csv_value(simulation, "simulations");
    if (!sim_valid) {
        return sim_valid.error();
    }
    auto lik_valid = detail::validate_indexed_write_csv_value(liks, "likelihood");
    if (!lik_valid) {
        return lik_valid.error();
    }
    auto maybe_space = detail::merge_indexed_write_csv_spaces(
        simulation.index_space(), "simulations", liks.index_space(), "likelihood");
    if (!maybe_space) {
        return maybe_space.error();
    }
    auto maybe_first = indexed_front(liks);
    if (!maybe_first) {
        return maybe_first.error();
    }
    auto param_names = first_non_empty_param_names(maybe_first.value().get());

    return detail::write_indexed_rows_csv(
        maybe_space.value(), param_names, std::move(path),
        [&](detail::CsvWriter& writer, const detail::CsvContext& base_ctx,
            const var::Coordinate& coord) -> Maybe_error<bool> {
            auto maybe_sim = simulation.at(coord);
            if (!maybe_sim) {
                return maybe_sim.error();
            }
            auto maybe_lik = liks.at(coord);
            if (!maybe_lik) {
                return maybe_lik.error();
            }
            const auto& sims = maybe_sim.value().get();
            const auto& batch = maybe_lik.value().get();
            if (sims.size() != batch.size()) {
                return error_message("number of Simulated_Recordings ", sims.size(),
                                     " differ from number of Likelihood states ", batch.size());
            }
            for (std::size_t sim_i = 0; sim_i < sims.size(); ++sim_i) {
                auto ok = emit_state_rows_with_experiment(writer, e, get<Recording>(sims[sim_i]()),
                                                          sim_i, batch[sim_i], base_ctx);
                if (!ok || !ok.value()) {
                    return ok;
                }
            }
            return true;
        });
}

// (2) State without Experiment indexing
template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path) {
    const auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    detail::CsvWriter writer(f, detail::get_param_names_if_any(lik));
    auto ok = emit_state_rows_without_experiment(writer, std::nullopt, lik);
    if (!ok || !ok.value()) {
        return ok.error()();
    }
    return path_;
}

template <class... Ms, class... Ids2, class... Fs>
auto calculate_Likelihood_diagnostics_impl(const std::vector<dMacro_State_Ev_gradient_all>& dy,
                                           std::vector<std::size_t> indices,
                                           std::type_identity<Vector_Space<Ms...>> /*tag*/,
                                           std::type_identity<Vector_Space<Ids2...>> /*tag*/,
                                           Fs&&... fs) {
    static_assert(sizeof...(Ms) == sizeof...(Fs));
    return Vector_Space<Ms..., Ids2...>(Ms(dy, indices, std::forward<Fs>(fs))..., Ids2()...);
}

template <class... Ids, bool... include_covariance, class... Fs>
auto calculate_Likelihood_diagnostics_evolution_correlation_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> const& indices,
    std::type_identity<
        Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>> /*tag*/,
    Fs&&... fs) {
    return Vector_Space<Moment_statistics<Sum<Ids>, include_covariance>...>(
        Moment_statistics<Sum<Ids>, include_covariance>(
            dy, indices, [&f = fs](const dMacro_State_Ev_gradient_all& d) {
                return Sum<Ids>(get<Evolution>(d)(), f)();
            })...);
}

template <class... Ids, bool... include_covariance, class... Fs>
auto calculate_Likelihood_diagnostics_evolution_cross_correlation_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t max_lag,std::vector<std::size_t> const& indices,
    std::type_identity<
        Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>> /*tag*/,
    Fs&&... fs) {
    auto make_one = [&](auto id_tag, auto& g) {
        using Id = typename decltype(id_tag)::type;
        auto stats = Series_Moment_statistics<Id, true>(
            dy, max_lag,indices,
            [](const dMacro_State_Ev_gradient_all& d) { return get<Evolution>(d)(); },
            g);
        auto maybe = make_report_cross(stats);
        return maybe ? std::move(maybe.value()) : Report_cross<Id>{};
    };
    return Vector_Space<Report_cross<Ids>...>(
        make_one(std::type_identity<Ids>{}, fs)...);
}

template <class... Ids, bool... include_covariance, class... Fs>
auto build_series_reports_integral_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t max_lag,
    std::vector<std::size_t> const& indices,
    std::type_identity<Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>>
    /*tag*/,
    Fs&&... fs) {
    auto make_one = [&](auto id_tag, auto& g) {
        using Id = typename decltype(id_tag)::type;
        auto stats = Series_Moment_statistics<Id, true>(
            dy, max_lag, indices,
            [](const dMacro_State_Ev_gradient_all& d) { return get<Evolution>(d)(); }, g);
        auto maybe = make_report_integral(stats);
        return maybe ? std::move(maybe.value()) : Report_integral<Id>{};
    };
    return Vector_Space<Report_integral<Ids>...>(
        make_one(std::type_identity<Ids>{}, fs)...);
}

template <class... Ids, bool... include_covariance, class... Fs>
auto build_series_reports_local_var_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t max_lag,
    std::vector<std::size_t> const& indices,
    std::type_identity<Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>>
    /*tag*/,
    Fs&&... fs) {
    auto make_one = [&](auto id_tag, auto& g) {
        using Id = typename decltype(id_tag)::type;
        auto stats = Series_Moment_statistics<Id, true>(
            dy, max_lag, indices,
            [](const dMacro_State_Ev_gradient_all& d) { return get<Evolution>(d)(); }, g);
        auto maybe = make_report_local_var(stats);
        return maybe ? std::move(maybe.value()) : Report_local_var<Id>{};
    };
    return Vector_Space<Report_local_var<Ids>...>(
        make_one(std::type_identity<Ids>{}, fs)...);
}

template <class... Ids, bool... include_covariance, class... Fs>
auto build_series_reports_local_cov_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t max_lag,
    std::vector<std::size_t> const& indices,
    std::type_identity<Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>>
    /*tag*/,
    Fs&&... fs) {
    auto make_one = [&](auto id_tag, auto& g) {
        using Id = typename decltype(id_tag)::type;
        auto stats = Series_Moment_statistics<Id, true>(
            dy, max_lag, indices,
            [](const dMacro_State_Ev_gradient_all& d) { return get<Evolution>(d)(); }, g);
        auto maybe = make_report_local_cov(stats);
        return maybe ? std::move(maybe.value()) : Report_local_cov<Id>{};
    };
    return Vector_Space<Report_local_cov<Ids>...>(
        make_one(std::type_identity<Ids>{}, fs)...);
}




template <class... Ids, bool... include_covariance, class... Id2s, class... Fs>
auto calculate_Likelihood_diagnostics_evolution_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> const& indices,
    std::type_identity<
        Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>> /*tag*/,
    std::type_identity<Evolution_of<Vector_Space<Id2s...>>> /*tag2*/, Fs&&... fs) {
    Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>..., Id2s...>> out;
    auto n = get<Evolution>(dy[0])().size();
    for (std::size_t i = 0; i < n; ++i) {
        out().emplace_back(calculate_Likelihood_diagnostics_impl(
            dy, indices,
            std::type_identity<Vector_Space<Moment_statistics<Ids, include_covariance>...>>(),
            std::type_identity<Vector_Space<Id2s...>>() /*tag2*/,
            [i, &f = fs](const dMacro_State_Ev_gradient_all& d) {
                return f(get<Evolution>(d)()[i]);
            }...));
    }

    return out;
}

enum class Diagnostic_preset {
    basic,               // base + Report_integral for 5 headline V's
    series_var,          // base + Report_local_var for 7 V's + Per_sample_derived
    series_cov,          // base + Report_local_cov for 7 V's + GFI@local_var + Per_sample_derived
    series_kernel,       // base + Report_cross for 5 V's, no Per_sample_derived
    series_kernel_full   // base + Report_cross for 7 V's + GFI@local_var + Per_sample_derived
};

template <Diagnostic_preset preset>
auto calculate_Likelihood_diagnostics_preset_f(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> const& indices, std::size_t max_lag) {
    auto sum_moments = calculate_Likelihood_diagnostics_evolution_correlation_impl(
        dy, indices,
        std::type_identity<Evolution_of<
            Vector_Space<Moment_statistics<logL>, Moment_statistics<elogL>,
                         Moment_statistics<r_std>, Moment_statistics<r2_std>,
                         Moment_statistics<trust_coefficient>, Moment_statistics<dlogL, true>,
                         Moment_statistics<Gaussian_Fisher_Information, false>>>>{},
        // macroIR guards logL / elogL / dlogL in calculate_logL at source (y_is_nan
        // early-return → hard zero). But r_std, r2_std, trust_coefficient are
        // computed independently and inherit NaN from NaN-observation samples.
        // Wrap each scalar read with an isnan check so NaN samples contribute
        // zero to Σ_t and its moments — physically "no measurement → no signal".
        [](const auto& evo_i) { return primitive(get<logL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<elogL>(evo_i))(); },
        [](const auto& evo_i) {
            const double v = primitive(get<r_std>(evo_i))();
            return std::isnan(v) ? 0.0 : v;
        },
        [](const auto& evo_i) {
            const double v = primitive(get<r_std>(evo_i))();
            return std::isnan(v) ? 0.0 : sqr(v);
        },
        [](const auto& evo_i) {
            const double v = get<trust_coefficient>(evo_i)();
            return std::isnan(v) ? 0.0 : v;
        },
        [](const auto& evo_i) {
            return parameter_vector_payload(derivative(get<logL>(evo_i))(),
                                            var::get_dx_of_dfdx(get<logL>(evo_i)));
        },
        [](const auto& evo_i) {
            // Gaussian Fisher information per sample. We read prediction
            // derivatives (dy_mean/dθ, dy_var/dθ) rather than logL's derivative,
            // so the NaN-observation guard that zeros logL & dlogL in
            // calculate_logL does NOT propagate here. For samples over long
            // integration intervals (e.g. 2000-s equilibration steps), the
            // prediction derivatives can contain secular-growth terms that
            // produce non-finite entries — if we pass those into Σ_t GFI_t we
            // poison H and compute_psd_decomp rejects every bootstrap sample.
            // Guard: if the computed per-sample GFI contains any non-finite
            // entry, return a zero SPD. Physically correct — a sample whose
            // derivatives are not numerically usable (or whose observation was
            // NaN-masked) contributes zero Fisher information.
            auto gfi = sqr_X<true>(derivative(get<y_mean>(evo_i))()) *
                           (1.0 / primitive(get<y_var>(evo_i))()) +
                       sqr_X<true>(derivative(get<y_var>(evo_i))()) *
                           (1.0 / 2.0 / sqr(primitive(get<y_var>(evo_i))()));
            if (!lapack::matrix_has_only_finite(gfi))
                gfi = SymPosDefMatrix<double>(gfi.nrows(), gfi.ncols(), 0.0);
            return parameter_spd_payload(std::move(gfi),
                                         var::get_dx_of_dfdx(get<y_mean>(evo_i)));
        });

    auto evol_moments = calculate_Likelihood_diagnostics_evolution_impl(
        dy, indices,
        std::type_identity<Evolution_of<Vector_Space<
            Moment_statistics<logL>, Moment_statistics<elogL>, Moment_statistics<y_mean>,
            Moment_statistics<y_var>, Moment_statistics<r_std>,
            Moment_statistics<trust_coefficient>, Moment_statistics<dlogL, true>,
            Moment_statistics<Gaussian_Fisher_Information, false>>>>{},
        std::type_identity<
            Evolution_of<Vector_Space<Sample_Distortion_Matrix, Distortion_Induced_Bias>>>{},

        // Same NaN guards as in sum_moments: r_std and trust_coefficient inherit
        // NaN from NaN-observation samples (unlike logL/elogL/dlogL which are
        // hard-zeroed in calculate_logL). y_mean / y_var are predictions and
        // stay finite, so they don't need guarding.
        [](const auto& evo_i) { return primitive(get<logL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<elogL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<y_mean>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<y_var>(evo_i))(); },
        [](const auto& evo_i) {
            const double v = primitive(get<r_std>(evo_i))();
            return std::isnan(v) ? 0.0 : v;
        },
        [](const auto& evo_i) {
            const double v = get<trust_coefficient>(evo_i)();
            return std::isnan(v) ? 0.0 : v;
        },
        [](const auto& evo_i) {
            return parameter_vector_payload(derivative(get<logL>(evo_i))(),
                                            var::get_dx_of_dfdx(get<logL>(evo_i)));
        },
        [](const auto& evo_i) {
            // Gaussian Fisher information per sample. We read prediction
            // derivatives (dy_mean/dθ, dy_var/dθ) rather than logL's derivative,
            // so the NaN-observation guard that zeros logL & dlogL in
            // calculate_logL does NOT propagate here. For samples over long
            // integration intervals (e.g. 2000-s equilibration steps), the
            // prediction derivatives can contain secular-growth terms that
            // produce non-finite entries — if we pass those into Σ_t GFI_t we
            // poison H and compute_psd_decomp rejects every bootstrap sample.
            // Guard: if the computed per-sample GFI contains any non-finite
            // entry, return a zero SPD. Physically correct — a sample whose
            // derivatives are not numerically usable (or whose observation was
            // NaN-masked) contributes zero Fisher information.
            auto gfi = sqr_X<true>(derivative(get<y_mean>(evo_i))()) *
                           (1.0 / primitive(get<y_var>(evo_i))()) +
                       sqr_X<true>(derivative(get<y_var>(evo_i))()) *
                           (1.0 / 2.0 / sqr(primitive(get<y_var>(evo_i))()));
            if (!lapack::matrix_has_only_finite(gfi))
                gfi = SymPosDefMatrix<double>(gfi.nrows(), gfi.ncols(), 0.0);
            return parameter_spd_payload(std::move(gfi),
                                         var::get_dx_of_dfdx(get<y_mean>(evo_i)));
        });

    constexpr double k_psd_rtol = 1e-10;
    constexpr double k_psd_atol = 0.0;
    for (auto& m : evol_moments()) {
        auto H_t = get<mean<Gaussian_Fisher_Information>>(get<Gaussian_Fisher_Information>(m)())();
        auto g_t = get<mean<dlogL>>(get<dlogL>(m)())();

        // One eigendecomp of H_t reused by both SDM_t and DIB_t.
        auto maybe_decomp_H_t = lapack::compute_psd_decomp(H_t.value(), "H_t", k_psd_rtol, k_psd_atol);
        lapack::PSDDecomposition W_H_t;
        if (maybe_decomp_H_t) W_H_t = std::move(maybe_decomp_H_t.value());

        get<Sample_Distortion_Matrix>(m)() =
            parameter_spd_payload(
                lapack::apply_normalized_congruence(
                    W_H_t, get<covariance<dlogL>>(get<dlogL>(m)())().value(),
                    "sample distortion subspace matrix", k_psd_rtol, k_psd_atol)
                    .value_or(SymPosDefMatrix<double>{}),
                H_t.parameters_ptr());

        get<Distortion_Induced_Bias>(m)() =
            parameter_vector_payload(
                lapack::apply_inverse_vector(W_H_t, g_t.value(),
                                             "Distortion-induced bias")
                    .value_or(Matrix<double>{}),
                H_t.parameters_ptr());
    }

    auto sum_r_std = Sum<Moment_statistics<r_std>>(evol_moments(),
                                                   [](const auto& m) { return get<r_std>(m)(); });

    auto sum_dlogL = Sum<Moment_statistics<dlogL, true>>(
        evol_moments(), [](const auto& m) { return get<dlogL>(m)(); });

    auto sum_Gaussian_Fisher_Information =
        Sum<Moment_statistics<Gaussian_Fisher_Information, false>>(
            evol_moments(), [](const auto& m) { return get<Gaussian_Fisher_Information>(m)(); });

    auto H = get<mean<Sum<Gaussian_Fisher_Information>>>(
        get<Sum<Gaussian_Fisher_Information>>(sum_moments)());

    auto J = get<covariance<Sum<dlogL>>>(get<Sum<dlogL>>(sum_moments)());

    auto J_sample = get<covariance<dlogL>>(sum_dlogL());
    auto score_mean = get<mean<Sum<dlogL>>>(get<Sum<dlogL>>(sum_moments)());

    // Distortion quantities are evaluated on the retained informative subspace.
    // Eigendecompose each anchor once and reuse across all derived matrices
    // (IDM, SDM, DCC, FC, DIB all share H; CDM uses J_sample; IDM2 uses SDM).
    auto maybe_W_H = lapack::compute_psd_decomp(H().value(), "H", k_psd_rtol, k_psd_atol);
    lapack::PSDDecomposition W_H;
    if (maybe_W_H) W_H = std::move(maybe_W_H.value());

    auto idm = Information_Distortion_Matrix(
        lapack::apply_normalized_congruence(W_H, J().value(), "IDM subspace matrix", k_psd_rtol,
                                            k_psd_atol)
            .value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());
    auto log_det_idm = log_Det<Information_Distortion_Matrix>(idm);

    auto sdm = Sample_Distortion_Matrix(
        lapack::apply_normalized_congruence(W_H, J_sample().value(),
                                            "sample distortion subspace matrix", k_psd_rtol,
                                            k_psd_atol)
            .value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());
    auto log_det_sdm = log_Det<Sample_Distortion_Matrix>(sdm);

    auto maybe_W_Js = lapack::compute_psd_decomp(J_sample().value(), "J_sample", k_psd_rtol,
                                                  k_psd_atol);
    lapack::PSDDecomposition W_Js;
    if (maybe_W_Js) W_Js = std::move(maybe_W_Js.value());

    auto cdm = Correlation_Distortion_Matrix(
        lapack::apply_normalized_congruence(W_Js, J().value(),
                                            "correlation distortion subspace matrix", k_psd_rtol,
                                            k_psd_atol)
            .value_or(SymPosDefMatrix<double>{}),
        J_sample().parameters_ptr());
    auto log_det_cdm = log_Det<Correlation_Distortion_Matrix>(cdm);

    auto maybe_W_SDM = lapack::compute_psd_decomp(sdm().value(), "SDM", k_psd_rtol, k_psd_atol);
    lapack::PSDDecomposition W_SDM;
    if (maybe_W_SDM) W_SDM = std::move(maybe_W_SDM.value());

    auto idm2 = Information_Distortion_Reconstituted(
        lapack::apply_sqrt_congruence(W_SDM, cdm().value(),
                                      "reconstituted distortion matrix", k_psd_rtol, k_psd_atol)
            .value_or(SymPosDefMatrix<double>{}),
        sdm().parameters_ptr());

    // Spectral-form identifiability diagnostics rooted in H's decomposition.
    // DCC = H⁻¹ J H⁻¹ inherits H's null subspace (Hv=0 ⇒ vᵀ·DCC·v = 0), so the
    // null projector and effective rank are identical for FC and DCC. The FC
    // and H eigenspectra differ only by reciprocation; we emit H's sorted
    // spectrum for both, and downstream analysis recovers FC spectrum via 1/λ.
    auto spectrum_H_mat = lapack::eigenvalue_spectrum(W_H);
    auto eff_rank_H = lapack::effective_rank(W_H, k_psd_rtol, k_psd_atol);
    auto cond_H = lapack::spectrum_condition_number(W_H);
    auto null_proj_matrix = lapack::null_space_projector(W_H, k_psd_rtol, k_psd_atol);
    auto worst_proj_matrix = lapack::worst_subspace_projector(W_H);

    auto fc = Fisher_Covariance(lapack::apply_inverse_as_matrix(W_H), H().parameters_ptr());
    auto log_det_fc = log_Det<Fisher_Covariance>(fc);
    // Spectral-form correlations: bounded in [-1, 1] by construction, computed
    // directly from (V, λ) without reconstructing a p×p covariance.
    using fc_corr_payload = typename Correlation_Of<Fisher_Covariance>::payload_type;
    auto corr_fc = Correlation_Of<Fisher_Covariance>(
        fc_corr_payload(lapack::fc_correlation_from_decomp(W_H), H().parameters_ptr()));
    auto spectrum_fc = Eigenvalue_Spectrum<Fisher_Covariance>(spectrum_H_mat);
    auto eff_rank_fc = Effective_Rank<Fisher_Covariance>(eff_rank_H);
    auto cond_fc = Spectrum_Condition_Number<Fisher_Covariance>(cond_H);
    using fc_null_payload = typename Null_Space_Projector<Fisher_Covariance>::payload_type;
    auto null_proj_fc = Null_Space_Projector<Fisher_Covariance>(
        fc_null_payload(null_proj_matrix, H().parameters_ptr()));
    using fc_worst_payload = typename Worst_Subspace_Projector<Fisher_Covariance>::payload_type;
    auto worst_proj_fc = Worst_Subspace_Projector<Fisher_Covariance>(
        fc_worst_payload(worst_proj_matrix, H().parameters_ptr()));

    auto dcc = Distortion_Corrected_Covariance(
        lapack::apply_inverse_congruence(W_H, J().value(), "DCC subspace matrix", k_psd_rtol,
                                         k_psd_atol)
            .value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());
    auto log_det_dcc = log_Det<Distortion_Corrected_Covariance>(dcc);
    using dcc_corr_payload =
        typename Correlation_Of<Distortion_Corrected_Covariance>::payload_type;
    auto corr_dcc = Correlation_Of<Distortion_Corrected_Covariance>(
        dcc_corr_payload(lapack::dcc_correlation_from_decomp(W_H, J().value()),
                         H().parameters_ptr()));
    auto spectrum_dcc = Eigenvalue_Spectrum<Distortion_Corrected_Covariance>(spectrum_H_mat);
    auto eff_rank_dcc = Effective_Rank<Distortion_Corrected_Covariance>(eff_rank_H);
    auto cond_dcc = Spectrum_Condition_Number<Distortion_Corrected_Covariance>(cond_H);
    using dcc_null_payload =
        typename Null_Space_Projector<Distortion_Corrected_Covariance>::payload_type;
    auto null_proj_dcc = Null_Space_Projector<Distortion_Corrected_Covariance>(
        dcc_null_payload(null_proj_matrix, H().parameters_ptr()));
    using dcc_worst_payload =
        typename Worst_Subspace_Projector<Distortion_Corrected_Covariance>::payload_type;
    auto worst_proj_dcc = Worst_Subspace_Projector<Distortion_Corrected_Covariance>(
        dcc_worst_payload(worst_proj_matrix, H().parameters_ptr()));

    auto dib = Distortion_Induced_Bias(
        lapack::apply_inverse_vector(W_H, score_mean().value(), "Distortion-induced bias")
            .value_or(Matrix<double>{}),
        H().parameters_ptr());

    // Extract Per_sample_derived (SDM, DIB only) from the full evol_moments.
    // Built lazily below when the preset needs it.
    auto build_per_sample_derived = [&]() {
        using derived_vs = var::Vector_Space<Sample_Distortion_Matrix, Distortion_Induced_Bias>;
        Evolution_of<derived_vs> out;
        out().reserve(evol_moments().size());
        for (auto& m : evol_moments()) {
            out().emplace_back(derived_vs(get<Sample_Distortion_Matrix>(m),
                                          get<Distortion_Induced_Bias>(m)));
        }
        return out;
    };

    // Base-tier pack (same across all presets).
    auto base_pack = push_back_var(
        std::move(sum_moments), std::move(sum_r_std), std::move(sum_dlogL),
        std::move(sum_Gaussian_Fisher_Information),
        std::move(idm), std::move(log_det_idm), std::move(idm2),
        std::move(sdm), std::move(log_det_sdm),
        std::move(cdm), std::move(log_det_cdm),
        std::move(fc), std::move(log_det_fc),
        std::move(corr_fc), std::move(spectrum_fc),
        std::move(eff_rank_fc), std::move(cond_fc), std::move(null_proj_fc),
        std::move(worst_proj_fc),
        std::move(dcc), std::move(log_det_dcc),
        std::move(corr_dcc), std::move(spectrum_dcc),
        std::move(eff_rank_dcc), std::move(cond_dcc), std::move(null_proj_dcc),
        std::move(worst_proj_dcc),
        std::move(dib));

    // Lambda list for scalar series observables: logL / elogL / y_mean / y_var
    // / r_std / trust / dlogL. Used by series_var, series_cov,
    // series_kernel_full.
    auto tag_seven = std::type_identity<Evolution_of<Vector_Space<
        Moment_statistics<logL>, Moment_statistics<elogL>,
        Moment_statistics<macrodr::y_mean>, Moment_statistics<macrodr::y_var>,
        Moment_statistics<macrodr::r_std>,
        Moment_statistics<macrodr::trust_coefficient>,
        Moment_statistics<dlogL, true>>>>{};

    auto lam_logL = [](const auto& evo_i) { return primitive(get<logL>(evo_i))(); };
    auto lam_elogL = [](const auto& evo_i) { return primitive(get<elogL>(evo_i))(); };
    auto lam_ymean = [](const auto& evo_i) { return primitive(get<y_mean>(evo_i))(); };
    auto lam_yvar = [](const auto& evo_i) { return primitive(get<y_var>(evo_i))(); };
    auto lam_rstd = [](const auto& evo_i) { return primitive(get<r_std>(evo_i))(); };
    auto lam_trust = [](const auto& evo_i) { return get<trust_coefficient>(evo_i)(); };
    auto lam_dlogL = [](const auto& evo_i) {
        return parameter_vector_payload(derivative(get<logL>(evo_i))(),
                                        var::get_dx_of_dfdx(get<logL>(evo_i)));
    };
    auto lam_gfi = [](const auto& evo_i) {
        return parameter_spd_payload(sqr_X<true>(derivative(get<y_mean>(evo_i))()) *
                                             (1.0 / primitive(get<y_var>(evo_i))()) +
                                         sqr_X<true>(derivative(get<y_var>(evo_i))()) *
                                             (1.0 / 2.0 / sqr(primitive(get<y_var>(evo_i))())),
                                     var::get_dx_of_dfdx(get<y_mean>(evo_i)));
    };

    // Core-five tag (logL / y_mean / y_var / r_std / dlogL) — used by basic
    // and series_kernel presets that don't include elogL/trust auxiliaries.
    auto tag_five = std::type_identity<Evolution_of<Vector_Space<
        Moment_statistics<logL>, Moment_statistics<macrodr::y_mean>,
        Moment_statistics<macrodr::y_var>, Moment_statistics<macrodr::r_std>,
        Moment_statistics<dlogL, true>>>>{};

    // GFI tag for the standalone local_var emission in series_cov and
    // series_kernel_full.
    auto tag_gfi = std::type_identity<Evolution_of<Vector_Space<
        Moment_statistics<Gaussian_Fisher_Information, false>>>>{};

    if constexpr (preset == Diagnostic_preset::basic) {
        auto series = build_series_reports_integral_impl(dy, max_lag, indices, tag_five,
                                                         lam_logL, lam_ymean, lam_yvar, lam_rstd,
                                                         lam_dlogL);
        return concatenate(std::move(base_pack), std::move(series));
    } else if constexpr (preset == Diagnostic_preset::series_var) {
        auto series = build_series_reports_local_var_impl(dy, max_lag, indices, tag_seven,
                                                          lam_logL, lam_elogL, lam_ymean, lam_yvar,
                                                          lam_rstd, lam_trust, lam_dlogL);
        auto per_sample = build_per_sample_derived();
        auto combined = concatenate(std::move(base_pack), std::move(series));
        return push_back_var(std::move(combined), std::move(per_sample));
    } else if constexpr (preset == Diagnostic_preset::series_cov) {
        auto series = build_series_reports_local_cov_impl(dy, max_lag, indices, tag_seven,
                                                          lam_logL, lam_elogL, lam_ymean, lam_yvar,
                                                          lam_rstd, lam_trust, lam_dlogL);
        auto gfi = build_series_reports_local_var_impl(dy, max_lag, indices, tag_gfi, lam_gfi);
        auto per_sample = build_per_sample_derived();
        auto combined = concatenate(concatenate(std::move(base_pack), std::move(series)),
                                    std::move(gfi));
        return push_back_var(std::move(combined), std::move(per_sample));
    } else if constexpr (preset == Diagnostic_preset::series_kernel) {
        auto series = calculate_Likelihood_diagnostics_evolution_cross_correlation_impl(
            dy, max_lag, indices, tag_five, lam_logL, lam_ymean, lam_yvar, lam_rstd, lam_dlogL);
        return concatenate(std::move(base_pack), std::move(series));
    } else {  // series_kernel_full
        auto series = calculate_Likelihood_diagnostics_evolution_cross_correlation_impl(
            dy, max_lag, indices, tag_seven, lam_logL, lam_elogL, lam_ymean, lam_yvar, lam_rstd,
            lam_trust, lam_dlogL);
        auto gfi = build_series_reports_local_var_impl(dy, max_lag, indices, tag_gfi, lam_gfi);
        auto per_sample = build_per_sample_derived();
        auto combined = concatenate(concatenate(std::move(base_pack), std::move(series)),
                                    std::move(gfi));
        return push_back_var(std::move(combined), std::move(per_sample));
    }
}








auto calculate_Likelihood_derivative_basic_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_basic {
    auto mt = mt_64i(seed);
    return bootstrap_it_to_Probit(
        &calculate_Likelihood_diagnostics_preset_f<Diagnostic_preset::basic>, dy,
        n_boostrap_samples, cis, mt, max_lag);
}

auto calculate_Likelihood_derivative_series_var_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_var {
    auto mt = mt_64i(seed);
    return bootstrap_it_to_Probit(
        &calculate_Likelihood_diagnostics_preset_f<Diagnostic_preset::series_var>, dy,
        n_boostrap_samples, cis, mt, max_lag);
}

auto calculate_Likelihood_derivative_series_cov_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_cov {
    auto mt = mt_64i(seed);
    return bootstrap_it_to_Probit(
        &calculate_Likelihood_diagnostics_preset_f<Diagnostic_preset::series_cov>, dy,
        n_boostrap_samples, cis, mt, max_lag);
}

auto calculate_Likelihood_derivative_series_kernel_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel {
    auto mt = mt_64i(seed);
    return bootstrap_it_to_Probit(
        &calculate_Likelihood_diagnostics_preset_f<Diagnostic_preset::series_kernel>, dy,
        n_boostrap_samples, cis, mt, max_lag);
}

auto calculate_Likelihood_derivative_series_kernel_full_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel_full {
    auto mt = mt_64i(seed);
    return bootstrap_it_to_Probit(
        &calculate_Likelihood_diagnostics_preset_f<Diagnostic_preset::series_kernel_full>, dy,
        n_boostrap_samples, cis, mt, max_lag);
}

// Explicit instantiations for the CLI-registered overloads to avoid link errors.
template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
        Experiment const&, Simulated_Recording<var::please_include<>> const&,
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
        Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
        Experiment const&, Simulated_Recording<var::please_include<>> const&,
        var::Indexed<
            Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
        Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
        var::Indexed<
            Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Experiment const&, Simulated_Recording<var::please_include<>> const&,
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Experiment const&, Simulated_Recording<var::please_include<>> const&,
        var::Indexed<
            Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<>, Macro_State, elogL, vlogL,
              Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
        var::Indexed<
            Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>> const&,
        std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, Simulated_Recording<var::please_include<>> const&,
    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&, std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&, std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, Simulated_Recording<var::please_include<>> const&,
    var::Indexed<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, var::Indexed<Simulated_Recording<var::please_include<>>> const&,
    var::Indexed<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, std::vector<Simulated_Recording<var::please_include<>>> const&,
    std::vector<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>> const&,
    std::vector<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, std::vector<Simulated_Recording<var::please_include<>>> const&,
    var::Indexed<
        std::vector<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>> const&,
    var::Indexed<
        std::vector<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>>> const&,
    std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
              Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        Experiment const&,
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
        dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
              Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        Experiment const&,
        var::Indexed<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const&,
        dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
              Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        Experiment const&,
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
        var::Indexed<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
              Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        Experiment const&,
        var::Indexed<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> const&,
        var::Indexed<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<Macro_State, elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<Macro_State, elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
        std::string);

template Maybe_error<std::string>
    write_csv<dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
        std::string);

}  // namespace macrodr::cmd
