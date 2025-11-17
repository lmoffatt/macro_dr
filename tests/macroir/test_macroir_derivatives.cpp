// Scaffolding for MacroIR (test_der_macroir) derivative tests, mirroring
// the pattern embedded in legacy/qmodel.h around lines ~4601–4629.
// These tests are set up to run in isolation once the required inputs
// (prior state, Qdtm, model parameters, channel count, etc.) are provided.

#include <catch_amalgamated.hpp>
#include <cstddef>

#include <derivative_test.h>   // var::test_derivative_clarke
#include <experiment.h>
#include <maybe_error.h>
#include <parameters.h>
#include <parameters_derivative.h>
#include <qmodel.h>

#include <macrodr/interface/IModel.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/patch_model.h>


TEST_CASE("macroir: MacroIR derivative (Clarke) — scaffolding", "[macroir][derivatives][macroir]") {
 


 
 
    // Finite-difference step for Clarke test
    const auto h = 1e-7;

    // Empty function registry (no memoized functions)
    auto f_no_memoi = var::create_empty_function_map();

    // Load plain model to obtain parameter transformations
    auto maybe_m = macrodr::cmd::load_model("scheme_1");
    if (!maybe_m) { UNSCOPED_INFO(maybe_m.error()()); }
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_parameter = macrodr::cmd::load_parameters(model0, "../data/scheme_1_PI_par.csv");
    if (!Maybe_parameter) { UNSCOPED_INFO(Maybe_parameter.error()()); }
    REQUIRE(Maybe_parameter.valid());
    auto par = std::move(Maybe_parameter.value());
const double fs=50e3;
    auto experiment = macrodr::cmd::load_experiment(
 "../data/Moffatt_Hume_2007_ATP_time_idealized_2.txt",fs,
0.0);

    auto maybe_observations = macrodr::cmd::load_recording("../data/Moffatt_Hume_2007_ATP_time_recording.txt");
    if (!maybe_observations) { UNSCOPED_INFO(maybe_observations.error()()); }
    REQUIRE(maybe_observations.valid());
    auto y = std::move(maybe_observations.value());

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
    using namespace macrodr; 
    auto Maybe_t_prior = Macro_DMR{}.init<return_predictions<0>>(dm, get<initial_agonist_concentration>(experiment));
    if (!Maybe_t_prior) { UNSCOPED_INFO(Maybe_t_prior.error()()); }
    REQUIRE(Maybe_t_prior.valid());
    auto t_prior = std::move(Maybe_t_prior.value());
    auto record = get<Recording_conditions>(experiment)();


    auto Nch = get<N_Ch_mean>(dm)()[0];
    
    Maybe_error<bool> result(true);
    for (std::size_t i_step = 0; i_step < record.size(); ++i_step)
    { 
    Agonist_evolution const& t_step =
                            get<Agonist_evolution>(get<Recording_conditions>(experiment)()[i_step]);


    auto Maybe_t_Qdtm = Macro_DMR{}.calc_Qdtm(f_no_memoi, dm, t_step, fs);
    if (!Maybe_t_Qdtm) { UNSCOPED_INFO(Maybe_t_Qdtm.error()()); }
    REQUIRE(Maybe_t_Qdtm.valid());
    auto t_Qdtm = std::move(Maybe_t_Qdtm.value());

    auto Maybe_prior_next=MacroR2<uses_recursive_aproximation<true>,
                                                uses_averaging_aproximation<2>,
                                                uses_variance_aproximation<true>,
                                                uses_taylor_variance_correction_aproximation<false>>{}(
                f_no_memoi, std::move(t_prior), t_Qdtm, dm, Nch, y()[i_step], fs);
    if (!Maybe_prior_next) { UNSCOPED_INFO(Maybe_prior_next.error()()); }
    REQUIRE(Maybe_prior_next.valid());
    auto t_prior = select<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL, eplogL, vplogL, macrodr::macror_algorithm>(std::move(Maybe_prior_next.value()));


    auto result_i = var::test_derivative_clarke<false>(
        [&](auto l_t_prior, auto const& l_Qdtm, auto const& l_m, auto const& l_Nch)
            -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(l_t_prior)>,
                var::Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL, eplogL, vplogL>>> {

                    
            auto Maybe_res = MacroR2<uses_recursive_aproximation<true>,
                                                uses_averaging_aproximation<2>,
                                                uses_variance_aproximation<true>,
                                                uses_taylor_variance_correction_aproximation<false>>{}(
                f_no_memoi, std::move(l_t_prior), l_Qdtm, l_m, l_Nch, y()[i_step], fs);
            if (!Maybe_res) return Maybe_res.error();
            return select<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL, eplogL, vplogL>(
                std::move(Maybe_res.value()));
        },
        h,
        t_prior, t_Qdtm, dm, Nch);
    result = result || result_i;          
                
            
    }

    REQUIRE(result);
    //
    // Once concrete types and builders are provided, adapt the snippet accordingly
    // and enable the actual derivative checks.
}
