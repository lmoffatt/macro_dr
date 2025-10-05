#include <CLI_function_table.h>
#include <macrodr/cmd/likelihood.h>

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

}  // namespace macrodr::cmd
