#include <CLI_function_table.h>
#include <derivative_fwd.h>
#include <derivative_operator.h>
#include <macrodr/cmd/likelihood.h>
#include <macrodr/dsl/type_name.h>
#include <parameters.h>
#include <concepts>
#include <string_view>
#include <tuple>
#include <utility>

#include "macrodr/cmd/load_model.h"
#include "qmodel.h"

namespace macrodr::cmd {

namespace {

template <class adaptive, class recursive, class averaging, class variance, class taylor, class Model,
          class FuncTable>
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

template <class adaptive, class recursive, class averaging, class variance, class taylor, class Model,
          class FuncTable>
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

template <class adaptive, class recursive, class averaging, class variance, class taylor, class Model,
          class FuncTable>
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
                                       const Recording& r) -> Maybe_error<Macro_State_Ev_predictions> {
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
                                       const Recording& r) -> Maybe_error<Macro_State_Ev_diagnostic> {
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
                                        const Recording& r) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

   
    return std::visit(
        [&](const auto& modelLikelihood) -> Maybe_error<dMacro_State_Ev_gradient_all> {
            return calculate_mdlikelihood_predictions_impl(modelLikelihood, ftbl3, par, e, r);
        },
        modelLikelihood_v);
}



auto calculate_likelihood(const interface::IModel<var::Parameters_values>& model0,
                          const var::Parameters_transformed& par, const Experiment& e,
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

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihood(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood(const interface::IModel<var::Parameters_values>& model0,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Hessian_minimal> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

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
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) -> Maybe_error<dMacro_State_Hessian_minimal> {
            return dlogLikelihood(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_diff_likelihood(const interface::IModel<var::Parameters_values>& model0,
                               const var::Parameters_transformed& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

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
        [&ftbl3, &par, &e, &r, delta_param](auto& modelLikelihood) {
            return diff_logLikelihood(ftbl3, modelLikelihood, par, r, e, delta_param);
        },
        modelLikelihood_v);
}

auto calculate_likelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions> {
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

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihoodPredictions(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}

auto calculate_likelihood_diagnostics(const interface::IModel<var::Parameters_values>& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
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

    auto par_values = par.to_value();
    return std::visit(
        [&ftbl3, &par_values, &e, &r](auto& modelLikelihood) {
            return logLikelihoodDiagnostic(ftbl3, modelLikelihood, par_values, r, e);
        },
        modelLikelihood_v);
}




auto calculate_dlikelihood_predictions(const interface::IModel<var::Parameters_values>& model0,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

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
        [&ftbl3, &par, &e, &r](auto& modelLikelihood) -> Maybe_error<dMacro_State_Ev_gradient_all> {
            return dlogLikelihoodPredictions(ftbl3, modelLikelihood, par, r, e);
        },
        modelLikelihood_v);
}

auto calculate_dlikelihood_predictions_model(const std::string& model_name,
                                             const var::Parameters_transformed& par, const Experiment& e,
                                             const Recording& r, bool adaptive_approximation,
                                             bool recursive_approximation,
                                             int averaging_approximation,
                                             bool variance_approximation,
                                             bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    auto ftbl3 = get_function_Table_maker_St("dummy", 100, 100)();

    auto nsub = Simulation_n_sub_dt(100);

    auto dmodel = get_model(model_name);
    if (!dmodel) {
        return dmodel.error();
    }
    auto model0_d = dmodel.value();
    return std::visit(
        [&](auto m_ptr)
            -> Maybe_error<dMacro_State_Ev_gradient_all> {
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
                [&ftbl3, &par, &e, &r](auto& modelLikelihood)
                    -> Maybe_error<dMacro_State_Ev_gradient_all> {
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
                ((ok = ok ? emit_component_rows(w,  macrodr::dsl::type_name_no_namespace<
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
    } else if constexpr (requires { x.nrows(); x.ncols(); x(0ul, 0ul); }) {
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
    } else if constexpr (requires { x.size(); x[0]; }) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            w.write_value(append_index(label, i), static_cast<double>(x[i]));
        }
        return true;
    } else {
        return error_message("write_csv_rows: unsupported value type ", label,
                             " (", macrodr::dsl::type_name<T>(), ")");
    }
}

template <class Writer, class Der>
Maybe_error<bool> emit_derivative_rows(Writer& w, std::string label, const Der& d) {
    if constexpr (requires { var::inside_out(d); }) {
        auto out = var::inside_out(d);
        if constexpr (requires { out.nrows(); out.ncols(); out(0ul, 0ul); }) {
            for (std::size_t r = 0; r < out.nrows(); ++r) {
                for (std::size_t c = 0; c < out.ncols(); ++c) {
                    auto maybe_ok = emit_derivative_rows(w, append_index(label, r, c), out(r, c));
                    if (!maybe_ok || !maybe_ok.value()) {
                        return maybe_ok;
                    }
                }
            }
            return true;
        } else if constexpr (requires { out.size(); out[0]; }) {
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
    } else if constexpr (var::is_derivative_v<Der> &&
                         is_of_this_template_type_v<var::untransformed_type_t<std::decay_t<Der>>,
                                                   var::Vector_Space>) {
        using VS = var::untransformed_type_t<std::decay_t<Der>>;
        using Tuple = typename vector_space_types<VS>::types;
        return std::apply(
            [&](auto... type_tags) -> Maybe_error<bool> {
                Maybe_error<bool> ok = true;
                ((ok = ok ? emit_component_rows(w,  macrodr::dsl::type_name_no_namespace<
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
        return error_message("write_csv_rows: unsupported derivative type ", label,
                             " (", macrodr::dsl::type_name<Der>(), ")");
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

// (1) Experiment + Simulation + State with per-sample Evolution
template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e, Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path) {
    const auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    const auto& conds = get<Recording_conditions>(e);
    const auto& y = get<Recording>(simulation());
    if (conds().size() != y().size()) {
        return error_message("Experiment samples ", conds().size(),
                             " differ from Recording samples ", y().size());
    }

    const auto& evo = get<Evolution>(lik);
    const auto& evo_vec = evo();
    if (evo_vec.size() != conds().size()) {
        return error_message("Evolution samples ", evo_vec.size(),
                             " differ from Recording samples ", conds().size());
    }

    const auto param_names = get_param_names_if_any(lik);

    f << "row_kind,n_step,sub_step,step_start,step_end,step_middle,agonist,patch_current,"
         "component,param_index,param_col,param_name,primitive,derivative_value\n";

    using Element = typename std::decay_t<decltype(evo)>::element_type;
    using ComponentTuple = typename vector_space_types<Element>::types;

    const double fs = get<Frequency_of_Sampling>(e)();
    for (std::size_t i = 0; i < conds().size(); ++i) {
        double step_start = get<Time>(conds()[i])();
        const auto& ag = get<Agonist_evolution>(conds()[i])();
        const auto& el = evo_vec[i];
        for (std::size_t j = 0; j < ag.size(); ++j) {
            const double ns = get<number_of_samples>(ag[j])();
            const double duration = ns / fs;
            const double step_end = step_start + duration;
            const double step_mid = 0.5 * (step_start + step_end);
            const double agonist = get<Agonist_concentration>(ag[j])();
            const double n_step =
                static_cast<double>(i) + static_cast<double>(j) / static_cast<double>(ag.size());
            const double patch = y()[i]();

            SampleWriter w{f, param_names, n_step, j, step_start, step_end, step_mid, agonist,
                           patch};

            auto handle_component = [&](auto type_tag) -> Maybe_error<bool> {
                using Comp = typename decltype(type_tag)::type;
                return emit_component_rows(w, macrodr::dsl::type_name_no_namespace<var::untransformed_type_t<Comp>>(), get<Comp>(el));
            };

            Maybe_error<bool> ok = true;
            []<std::size_t... Is>(auto&& handle, Maybe_error<bool>& ok_ref,
                                  std::index_sequence<Is...>) {
                ((ok_ref = ok_ref ? handle(std::tuple_element_t<Is, ComponentTuple>{}) : ok_ref),
                 ...);
            }(handle_component, ok, std::make_index_sequence<std::tuple_size_v<ComponentTuple>>{});
            if (!ok || !ok.value()) {
                return ok.error()();
            }

            step_start = step_end;
        }
    }

    return path_;
}

// (2) State without Evolution (no Experiment/Simulation indexing)
template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path) {
    const auto path_ = path + ".csv";
    std::ofstream f(path_);
    if (!f.is_open()) {
        return error_message("cannot open ", path_);
    }

    const auto param_names = get_param_names_if_any(lik);
    f << "meta,n_params," << param_names.size() << "\n";
    for (std::size_t i = 0; i < param_names.size(); ++i) {
        f << "param," << i << "," << param_names[i] << "\n";
    }

    f << "row_kind,n_step,sub_step,step_start,step_end,step_middle,agonist,patch_current,"
         "component,param_index,param_col,param_name,primitive,derivative_value\n";

    StateWriter w{f, param_names};

	    auto handle_root = [&](auto type_tag) -> Maybe_error<bool> {
	        using Id = typename decltype(type_tag)::type;
	        if constexpr (is_of_this_template_type_v<Id, Evolution_of>) {
	            // Skip evolution containers in the state-only view.
	            return true;
	        } else if constexpr (var::gets_it_c<TMacro_State<vVars...> const&, Id>) {
	            return emit_component_rows(w, macrodr::dsl::type_name_no_namespace<var::untransformed_type_t<Id>>(), get<Id>(lik));
	        } else {
	            return true;
	        }
    };

    Maybe_error<bool> ok = true;
    ok = ok ? handle_root(std::type_identity<logL>{}) : ok;
    ok = ok ? handle_root(std::type_identity<Patch_State>{}) : ok;
    ((ok = ok ? handle_root(std::type_identity<vVars>{}) : ok), ...);
    if (!ok || !ok.value()) {
        return ok.error()();
    }

    return path_;
}

// Explicit instantiations for the CLI-registered overloads to avoid link errors.
template Maybe_error<std::string> write_csv<
    var::please_include<>, Macro_State, elogL, vlogL,
    Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
    Experiment const&, Simulated_Recording<var::please_include<>> const&,
    Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<
    var::please_include<>, Macro_State, elogL, vlogL,
    Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
    Experiment const&, Simulated_Recording<var::please_include<>> const&,
    Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
    std::string);

	template Maybe_error<std::string> write_csv<
	    var::please_include<>, dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
	    Experiment const&, Simulated_Recording<var::please_include<>> const&,
	    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
	    std::string);

	template Maybe_error<std::string> write_csv<
	    var::please_include<Only_Ch_Curent_Sub_Evolution>, dMacro_State,
	    Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
	    Experiment const&, Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
	    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
	    std::string);

template Maybe_error<std::string> write_csv<Macro_State, elogL, vlogL,
                                            Evolution_of<add_t<Vector_Space<>, predictions_element>>>(
    Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, predictions_element>>> const&,
    std::string);

template Maybe_error<std::string> write_csv<Macro_State, elogL, vlogL,
                                            Evolution_of<add_t<Vector_Space<>, diagnostic_element>>>(
    Macro_State<elogL, vlogL, Evolution_of<add_t<Vector_Space<>, diagnostic_element>>> const&,
    std::string);

	template Maybe_error<std::string>
	write_csv<dMacro_State, Evolution_of<add_t<Vector_Space<>, gradient_all_element>>>(
	    dMacro_State<Evolution_of<add_t<Vector_Space<>, gradient_all_element>>> const&,
	    std::string);

}  // namespace macrodr::cmd
