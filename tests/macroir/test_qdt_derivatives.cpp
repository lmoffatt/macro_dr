// Qdt derivative tests extracted and isolated from legacy/qmodel.h
// Validates derivatives via var::test_derivative_clarke on calc_Qdtm.

#include <catch_amalgamated.hpp>

#include <derivative_test.h>
#include <experiment.h>
#include <function_measure_verification_and_optimization.h>
#include <maybe_error.h>
#include <parameters.h>
#include <parameters_derivative.h>
#include <qmodel.h>

#include <macrodr/interface/IModel.h>
#include <macrodr/cmd/load_model.h>

#include <utility>
#include <string>
#include <iostream>

// Bring commonly used types/names into scope for clarity
using macrodr::P;
using macrodr::gmean_i;
using macrodr::gtotal_ij;
using macrodr::gsqr_i;
using macrodr::gvar_i;

TEST_CASE("macroir: Qdt derivative (Clarke)", "[macroir][derivatives][qdt]") {
    // Finite-difference step for Clarke test
    const auto h = 1e-7;

    // Empty function registry (no memoized functions)
    auto f_no_memoi = var::create_empty_function_map();

    // Single ATP step: 10 samples at 1000 uM, sampling frequency 50 kHz
    auto t_step = macrodr::ATP_step(macrodr::number_of_samples(10), macrodr::ATP_concentration(1000));
    double fs = 50e3;

    // Load plain model to obtain parameter transformations
    auto maybe_m = macrodr::cmd::load_model("scheme_1");
    if (!maybe_m) { UNSCOPED_INFO(maybe_m.error()()); }
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    var::Parameters_Transformations par = model0->parameters_transformations();
    auto par_values = par.standard_parameter();
    auto par_transformed = par_values.to_transformed();

    // Load derivative-aware model
    auto dmodel = macrodr::cmd::load_dmodel("scheme_1");
    if (!dmodel) { UNSCOPED_INFO(dmodel.error()()); }
    REQUIRE(dmodel.valid());
    auto model0_d = std::move(dmodel.value());

    // Differentiate parameters and apply model
    auto dp = var::selfDerivative(par_transformed);
    auto dpp = dp.to_value();
    auto maybe_dm = (*model0_d)(dpp);
    if (!maybe_dm) { UNSCOPED_INFO(maybe_dm.error()()); }
    REQUIRE(maybe_dm.valid());
    auto dm = std::move(maybe_dm.value());

    // Run Clarke derivative test on calc_Qdtm with stabilizers disabled
    auto test_der_t_Qdtm = test_derivative_clarke(
        [&t_step, &fs, &f_no_memoi](auto const& l_m)
            -> Maybe_error<
                var::Transfer_Op_to<std::decay_t<decltype(l_m)>,
                                    var::Vector_Space<P, gmean_i, gtotal_ij, gsqr_i, gvar_i>>> {
            auto maybe_res = macrodr::Macro_DMR{}.calc_Qdtm<macrodr::StabilizerPolicyDisabled>(
                f_no_memoi, l_m, t_step, fs);
            if (!maybe_res.valid()) {
                return maybe_res.error();
            }
            return select<P, gmean_i, gtotal_ij, gsqr_i, gvar_i>(std::move(maybe_res.value()));
        },
        h, dm);

    if (!test_der_t_Qdtm) { UNSCOPED_INFO(test_der_t_Qdtm.error()()); }
    REQUIRE(test_der_t_Qdtm);
}
