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
            .get_variant();
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
    auto maybe_modelLikelihood =
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
            .get_variant();
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
    auto maybe_modelLikelihood =
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
            .get_variant();
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

    auto maybe_modelLikelihood =
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
            .get_variant();
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

    auto maybe_modelLikelihood =
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
            .get_variant();
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
    auto maybe_modelLikelihood =
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
            .get_variant();
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
            auto maybe_modelLikelihood =

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
                    .get_variant();
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
                                                     const Lik& lik) {
    detail::CsvContext state_ctx;
    state_ctx.scope = "state";
    state_ctx.simulation_index = simulation_index;

    auto ok = emit_state_members(writer, state_ctx, lik);
    if (!ok || !ok.value()) {
        return ok;
    }

    if constexpr (macrodr::has_var_c<Lik const&, Evolution>) {
        const auto& evolution = get<Evolution>(lik)();
        for (std::size_t i = 0; i < evolution.size(); ++i) {
            detail::CsvContext evo_ctx;
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
                                                  const Lik& lik) {
    const auto& conditions = get<Recording_conditions>(e);
    if (conditions().size() != recording().size()) {
        return error_message("Experiment samples ", conditions().size(),
                             " differ from Recording samples ", recording().size());
    }

    auto ok = emit_state_rows_without_experiment(writer, simulation_index, lik);
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

            detail::CsvContext evo_ctx;
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

    std::vector<std::string> param_names;
    for (const auto& lik : liks) {
        param_names = detail::get_param_names_if_any(lik);
        if (!param_names.empty()) {
            break;
        }
    }

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
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> indices,
    std::type_identity<
        Evolution_of<Vector_Space<Moment_statistics<Ids, include_covariance>...>>> /*tag*/,
    Fs&&... fs) {
    return Vector_Space<Moment_statistics<Sum<Ids>, include_covariance>...>(
        Moment_statistics<Sum<Ids>, include_covariance>(
            dy, indices, [&f = fs](const dMacro_State_Ev_gradient_all& d) {
                return Sum<Ids>(get<Evolution>(d)(), f)();
            })...);
}

template <class... Ids, bool... include_covariance, class... Id2s, class... Fs>
auto calculate_Likelihood_diagnostics_evolution_impl(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> indices,
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

auto calculate_Likelihood_diagnostics_evolution_f(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::vector<std::size_t> const& indices) {
    auto sum_moments = calculate_Likelihood_diagnostics_evolution_correlation_impl(
        dy, indices,
        std::type_identity<Evolution_of<
            Vector_Space<Moment_statistics<logL>, Moment_statistics<elogL>,
                         Moment_statistics<r_std>, Moment_statistics<r2_std>,
                         Moment_statistics<trust_coefficient>, Moment_statistics<dlogL, true>,
                         Moment_statistics<Gaussian_Fisher_Information, false>>>>{},
        [](const auto& evo_i) { return primitive(get<logL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<elogL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<r_std>(evo_i))(); },
        [](const auto& evo_i) { return sqr(primitive(get<r_std>(evo_i))()); },
        [](const auto& evo_i) { return get<trust_coefficient>(evo_i)(); },
        [](const auto& evo_i) {
            return parameter_vector_payload(derivative(get<logL>(evo_i))(),
                                            var::get_dx_of_dfdx(get<logL>(evo_i)));
        },
        [](const auto& evo_i) {
            return parameter_spd_payload(
                sqr_X<true>(derivative(get<y_mean>(evo_i))()) *
                        (1.0 / primitive(get<y_var>(evo_i))()) +
                    sqr_X<true>(derivative(get<y_var>(evo_i))()) *
                        (1.0 / 2.0 / sqr(primitive(get<y_var>(evo_i))())),
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

        [](const auto& evo_i) { return primitive(get<logL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<elogL>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<y_mean>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<y_var>(evo_i))(); },
        [](const auto& evo_i) { return primitive(get<r_std>(evo_i))(); },
        [](const auto& evo_i) { return get<trust_coefficient>(evo_i)(); },
        [](const auto& evo_i) {
            return parameter_vector_payload(derivative(get<logL>(evo_i))(),
                                            var::get_dx_of_dfdx(get<logL>(evo_i)));
        },
        [](const auto& evo_i) {
            return parameter_spd_payload(
                sqr_X<true>(derivative(get<y_mean>(evo_i))()) *
                        (1.0 / primitive(get<y_var>(evo_i))()) +
                    sqr_X<true>(derivative(get<y_var>(evo_i))()) *
                        (1.0 / 2.0 / sqr(primitive(get<y_var>(evo_i))())),
                var::get_dx_of_dfdx(get<y_mean>(evo_i)));
        });

    for (auto& m : evol_moments()) {
        auto H_t = get<mean<Gaussian_Fisher_Information>>(get<Gaussian_Fisher_Information>(m)())();
        auto g_t = get<mean<dlogL>>(get<dlogL>(m)())();

        get<Sample_Distortion_Matrix>(m)() =
            parameter_spd_payload(
                sample_distortion_matrix_subspace(
                    H_t.value(), get<covariance<dlogL>>(get<dlogL>(m)())().value())
                    .value_or(SymPosDefMatrix<double>{}),
                H_t.parameters_ptr());

        get<Distortion_Induced_Bias>(m)() =
            parameter_vector_payload(
                distortion_induced_bias_subspace(H_t.value(), g_t.value())
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
    auto idm = Information_Distortion_Matrix(
        idm_matrix_subspace(H().value(), J().value()).value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());

    auto sdm = Sample_Distortion_Matrix(
        sample_distortion_matrix_subspace(H().value(), J_sample().value())
            .value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());

    auto cdm = Correlation_Distortion_Matrix(
        correlation_distortion_matrix_subspace(J_sample().value(), J().value())
            .value_or(SymPosDefMatrix<double>{}),
        J_sample().parameters_ptr());

    auto idm2 = Information_Distortion_Reconstituted(
        c_h_r_c_h_matrix_subspace(sdm().value(), cdm().value()).value_or(SymPosDefMatrix<double>{}),
        sdm().parameters_ptr());
    
    auto dcc = Distortion_Corrected_Covariance(
        dcc_matrix_subspace(H().value(), J().value()).value_or(SymPosDefMatrix<double>{}),
        H().parameters_ptr());

    auto dib = Distortion_Induced_Bias(
        distortion_induced_bias_subspace(H().value(), score_mean().value()).value_or(Matrix<double>{}),
        H().parameters_ptr());

    return push_back_var(std::move(sum_moments), std::move(sum_r_std), std::move(sum_dlogL),
                         std::move(sum_Gaussian_Fisher_Information), std::move(idm), std::move(idm2), std::move(sdm),
                         std::move(cdm), std::move(dcc), std::move(dib), std::move(evol_moments));
}

auto calculate_Likelihood_derivative_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed)
    -> Analisis_derivative_diagnostic {

    auto mt = mt_64i(seed);            
    return bootstrap_it_to_Probit(&calculate_Likelihood_diagnostics_evolution_f, dy,
                                  n_boostrap_samples, cis, mt);
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
              Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
        Experiment const&, Simulated_Recording<var::please_include<>> const&,
        Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
        std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, Simulated_Recording<var::please_include<>> const&,
    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&, std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
    Experiment const&, std::vector<Simulated_Recording<var::please_include<>>> const&,
    std::vector<dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>> const&,
    std::string);

template Maybe_error<std::string>
    write_csv<var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
              Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
        Experiment const&,
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
        dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
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
