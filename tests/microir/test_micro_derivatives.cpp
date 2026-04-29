// Per-step Clarke-bracket derivative test for MicroR2 (the micro-fold inner
// step used by dlog_Likelihood_micro). Mirrors tests/macroir/test_macroir_derivatives.cpp
// but exercises the discrete-mixture micro path on the simplest model
// (scheme_CO, k=2) with a small channel count (Nch=5) so the microstate
// enumeration stays cheap.
//
// What this verifies: for every step i of the experiment, given the current
// derivative-typed prior probability vector and the derivative-typed model,
// MicroR2's analytical d(logL_contribution)/d(parameter) matches finite
// differences along every parameter direction (Clarke first-order test).
//
// What it does NOT verify: the dispatch wiring in calculate_*_predictions, the
// Evolution_of bookkeeping, or the final project_Micro_to_Macro output. Those
// are end-to-end behaviors better covered by .macroir smoke scripts.

#include <catch_amalgamated.hpp>
#include <cstddef>

#include <derivative_test.h>
#include <experiment.h>
#include <maybe_error.h>
#include <parameters.h>
#include <parameters_derivative.h>
#include <qmodel.h>

#include <micro_full.h>

#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/patch_model.h>
#include <macrodr/interface/IModel.h>


TEST_CASE("microir: MicroR2 step derivative (Clarke) — averaging=2",
          "[microir][derivatives]") {
    using namespace macrodr;

    const double h = 1e-7;
    auto f_no_memoi = var::create_empty_function_map();

    // Plain model handle (used to load parameter transformations).
    auto maybe_m = cmd::load_model("scheme_CO");
    if (!maybe_m) UNSCOPED_INFO(maybe_m.error()());
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    // Small-Nch parameter file: scheme_CO with Num_ch_mean=5 keeps the
    // microstate enumeration to (5+1) = 6 states, so per-step MicroR2 +
    // derivative test is fast.
    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO(Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    const double fs = 50e3;
    auto experiment = cmd::load_experiment(
        "../data/Moffatt_Hume_2007_ATP_time_idealized_2.txt", fs, 0.0);

    auto maybe_obs = cmd::load_recording("../data/Moffatt_Hume_2007_ATP_time_recording.txt");
    if (!maybe_obs) UNSCOPED_INFO(maybe_obs.error()());
    REQUIRE(maybe_obs.valid());
    auto y = std::move(maybe_obs.value());

    auto par_values = par.standard_parameter();
    auto par_transformed = par_values.to_transformed();

    // Derivative-aware model.
    auto dmodel = cmd::load_dmodel("scheme_CO");
    if (!dmodel) UNSCOPED_INFO(dmodel.error()());
    REQUIRE(dmodel.valid());
    auto model0_d = std::move(dmodel.value());

    auto dp = var::selfDerivative(par_transformed);
    auto dpp = dp.to_value();
    auto maybe_dm = (*model0_d)(dpp);
    if (!maybe_dm) UNSCOPED_INFO(maybe_dm.error()());
    REQUIRE(maybe_dm.valid());
    auto dm = std::move(maybe_dm.value());

    Macro_DMR macro_dmr{};

    auto Maybe_init = macro_dmr.init(dm);
    if (!Maybe_init) UNSCOPED_INFO(Maybe_init.error()());
    REQUIRE(Maybe_init.valid());

    // Build the microstate enumeration once.
    std::size_t k_states = get<N_St>(dm)();
    auto Nchs = get<N_Ch_mean>(dm)();
    std::size_t N_ch = static_cast<std::size_t>(std::round(var::primitive(Nchs[std::size_t{0}])));
    auto full_states = create_Micro_state_Num_ch(N_ch, k_states);

    // Lift the init Patch_State to a Micro_Patch_State, extract the initial
    // derivative-typed probability vector.
    auto Maybe_mps = lift_Macro_to_Micro(full_states, Maybe_init.value());
    if (!Maybe_mps) UNSCOPED_INFO(Maybe_mps.error()());
    REQUIRE(Maybe_mps.valid());

    auto prior_probs = get<Micro_P_mean>(Maybe_mps.value()())();
    using T_probs = std::decay_t<decltype(prior_probs)>;

    auto record = get<Recording_conditions>(experiment)();
    // Cap the per-step loop. The recording was generated with Nch=5000, but
    // our test fixture uses Nch=5 for tractable microstate enumeration. After
    // ~150 steps the prior drifts into a regime where the discrete predicted
    // currents (which take only N+1 = 6 distinct values) are far enough from
    // the observed y that every microstate's Gaussian likelihood underflows
    // to zero and micro_full_step_avg2 fails with "no valid (n_s, Nij) with
    // finite likelihood" — a data/model-mismatch artifact, not a derivative
    // bug. Capping at 100 keeps the test in a regime where MicroR2 stays
    // numerically well-behaved while still exercising ~1000 Clarke-bracket
    // derivative checks (100 steps × ~5 parameters × 2 directions).
    const std::size_t step_cap = std::min<std::size_t>(record.size(), 100);
    Maybe_error<bool> result(true);

    using recursive = uses_recursive_aproximation<true>;
    using averaging = uses_averaging_aproximation<2>;
    using variance = uses_variance_aproximation<true>;

    // Walk the fold step-by-step. Per step: run test_derivative_clarke on a
    // closure that recomputes Qdt and runs MicroR2, returning the scalar
    // logL_contribution. Then advance prior_probs via the actual MicroR2
    // step (using outside_in to re-pack per-element new_probs into the
    // monolithic Derivative<Matrix> the next iteration expects).
    for (std::size_t i_step = 0; i_step < step_cap; ++i_step) {
        Agonist_evolution const& t_step =
            get<Agonist_evolution>(record[i_step]);
        double y_obs = y()[i_step].value();

        auto result_i = var::test_derivative_clarke<false>(
            [&macro_dmr, &t_step, fs, &full_states, y_obs, &f_no_memoi](
                auto const& l_dm, auto const& l_probs)
                -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(l_dm)>, logL>> {
                auto Maybe_l_Qdt = macro_dmr.calc_Qdt(f_no_memoi, l_dm, t_step, fs);
                if (!Maybe_l_Qdt) return Maybe_l_Qdt.error();
                auto const& l_Qdt = Maybe_l_Qdt.value();
                auto const& l_P_single = get<P>(l_Qdt)();
                auto const& l_gmean_ij = get<gmean_ij>(l_Qdt)();
                auto const& l_gvar_ij = get<gvar_ij>(l_Qdt)();
                double l_n =
                    static_cast<double>(get<number_of_samples>(l_Qdt).value());
                auto l_sigma2 = get<Current_Noise>(l_dm).value() * fs / l_n;
                auto Maybe_l_step =
                    MicroR2<recursive, averaging, variance>{}(
                        l_probs, l_P_single, l_gmean_ij, l_gvar_ij, l_sigma2,
                        full_states, y_obs);
                if (!Maybe_l_step) return Maybe_l_step.error();
                // test_clarke_brackets_init requires a Var-typed (or
                // Vector_Space) Derivative; wrap the raw Derivative<double>
                // logL_contribution into Derivative<logL>.
                return var::build<logL>(std::move(Maybe_l_step.value().logL_contribution));
            },
            h, dm, prior_probs);
        if (!result_i) UNSCOPED_INFO(result_i.error()());
        result = std::move(result) && std::move(result_i);

        // Advance the fold: compute the actual step output and update probs
        // (via outside_in to convert per-element new_probs to monolithic).
        auto Maybe_Qdt = macro_dmr.calc_Qdt(f_no_memoi, dm, t_step, fs);
        if (!Maybe_Qdt) UNSCOPED_INFO(Maybe_Qdt.error()());
        REQUIRE(Maybe_Qdt.valid());
        auto const& t_Qdt = Maybe_Qdt.value();
        auto const& P_single_m = get<P>(t_Qdt)();
        auto const& t_gmean_ij = get<gmean_ij>(t_Qdt)();
        auto const& t_gvar_ij = get<gvar_ij>(t_Qdt)();
        double n_samples =
            static_cast<double>(get<number_of_samples>(t_Qdt).value());
        auto sigma2_obs = get<Current_Noise>(dm).value() * fs / n_samples;

        auto Maybe_step = MicroR2<recursive, averaging, variance>{}(
            prior_probs, P_single_m, t_gmean_ij, t_gvar_ij, sigma2_obs,
            full_states, y_obs);
        if (!Maybe_step) UNSCOPED_INFO(Maybe_step.error()());
        REQUIRE(Maybe_step.valid());
        auto step = std::move(Maybe_step.value());

        // T_probs is always Derivative<Matrix<double>> in this test (dm is
        // built via selfDerivative). MicroR2 emits per-element new_probs;
        // re-pack into the monolithic Derivative<Matrix> the next iteration
        // expects.
        static_assert(var::is_derivative_v<T_probs>,
                      "test fixture must run with derivative-typed probs");
        prior_probs = var::outside_in(step.new_probs, prior_probs.dx());
    }

    if (!result) UNSCOPED_INFO(result.error()());
    REQUIRE(result.valid());
}


// Same fixture and shape as the MicroR2 (multinomial) test above, but exercises
// the lifted micro path: builds Patch_Model_micro via lift_Patch_Model_to_Micro,
// drives macro Calc_Qdt on it (yielding M×M P, gmean_ij, gvar_ij), then runs
// micro_full_step_avg2_lifted (plain Bayes per (s, s') pair, normal_pdf only,
// no per-cell logarithms) against test_derivative_clarke.
TEST_CASE("microir: micro_full_step_avg2_lifted derivative (Clarke) — averaging=2",
          "[microir][derivatives][lifted]") {
    using namespace macrodr;

    const double h = 1e-7;
    auto f_no_memoi = var::create_empty_function_map();

    auto maybe_m = cmd::load_model("scheme_CO");
    if (!maybe_m) UNSCOPED_INFO(maybe_m.error()());
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO(Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    const double fs = 50e3;
    auto experiment = cmd::load_experiment(
        "../data/Moffatt_Hume_2007_ATP_time_idealized_2.txt", fs, 0.0);

    auto maybe_obs = cmd::load_recording("../data/Moffatt_Hume_2007_ATP_time_recording.txt");
    if (!maybe_obs) UNSCOPED_INFO(maybe_obs.error()());
    REQUIRE(maybe_obs.valid());
    auto y = std::move(maybe_obs.value());

    auto par_values = par.standard_parameter();
    auto par_transformed = par_values.to_transformed();

    auto dmodel = cmd::load_dmodel("scheme_CO");
    if (!dmodel) UNSCOPED_INFO(dmodel.error()());
    REQUIRE(dmodel.valid());
    auto model0_d = std::move(dmodel.value());

    auto dp = var::selfDerivative(par_transformed);
    auto dpp = dp.to_value();
    auto maybe_dm = (*model0_d)(dpp);
    if (!maybe_dm) UNSCOPED_INFO(maybe_dm.error()());
    REQUIRE(maybe_dm.valid());
    auto dm_macro = std::move(maybe_dm.value());

    Macro_DMR macro_dmr{};

    auto Maybe_init = macro_dmr.init(dm_macro);
    if (!Maybe_init) UNSCOPED_INFO(Maybe_init.error()());
    REQUIRE(Maybe_init.valid());

    std::size_t k_states = get<N_St>(dm_macro)();
    auto Nchs = get<N_Ch_mean>(dm_macro)();
    std::size_t N_ch = static_cast<std::size_t>(std::round(var::primitive(Nchs[std::size_t{0}])));
    auto full_states = create_Micro_state_Num_ch(N_ch, k_states);

    auto Maybe_mps = lift_Macro_to_Micro(full_states, Maybe_init.value());
    if (!Maybe_mps) UNSCOPED_INFO(Maybe_mps.error()());
    REQUIRE(Maybe_mps.valid());

    // Lift the derivative-typed macro Patch_Model to its microstate twin.
    auto dm_micro = lift_Patch_Model_to_Micro(dm_macro, full_states);

    auto prior_probs = get<Micro_P_mean>(Maybe_mps.value()())();
    using T_probs = std::decay_t<decltype(prior_probs)>;

    auto record = get<Recording_conditions>(experiment)();
    const std::size_t step_cap = std::min<std::size_t>(record.size(), 100);
    Maybe_error<bool> result(true);

    for (std::size_t i_step = 0; i_step < step_cap; ++i_step) {
        Agonist_evolution const& t_step =
            get<Agonist_evolution>(record[i_step]);
        if (t_step().size() != 1) continue;  // lifted path is single-substep
        Agonist_step const& sub_step = t_step()[0];
        double y_obs = y()[i_step].value();

        auto result_i = var::test_derivative_clarke<false>(
            [&macro_dmr, &sub_step, fs, y_obs, &f_no_memoi, &full_states](
                auto const& l_dm_macro, auto const& l_probs)
                -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(l_dm_macro)>, logL>> {
                auto l_dm_micro = lift_Patch_Model_to_Micro(l_dm_macro, full_states);
                auto Maybe_l_Qdt = macro_dmr.calc_Qdt(f_no_memoi, l_dm_micro, sub_step, fs);
                if (!Maybe_l_Qdt) return Maybe_l_Qdt.error();
                auto const& l_Qdt = Maybe_l_Qdt.value();
                auto const& l_P = get<P>(l_Qdt)();
                auto const& l_gmean_ij = get<gmean_ij>(l_Qdt)();
                auto const& l_gvar_ij = get<gvar_ij>(l_Qdt)();
                double l_n =
                    static_cast<double>(get<number_of_samples>(l_Qdt).value());
                auto l_sigma2 = get<Current_Noise>(l_dm_micro).value() * fs / l_n;
                auto Maybe_l_step =
                    micro_full_step_avg2_lifted(l_probs, l_P, l_gmean_ij, l_gvar_ij,
                                                 l_sigma2, y_obs);
                if (!Maybe_l_step) return Maybe_l_step.error();
                return var::build<logL>(std::move(Maybe_l_step.value().logL_contribution));
            },
            h, dm_macro, prior_probs);
        if (!result_i) UNSCOPED_INFO(result_i.error()());
        result = std::move(result) && std::move(result_i);

        // Advance the fold using the actual lifted step.
        auto Maybe_Qdt = macro_dmr.calc_Qdt(f_no_memoi, dm_micro, sub_step, fs);
        if (!Maybe_Qdt) UNSCOPED_INFO(Maybe_Qdt.error()());
        REQUIRE(Maybe_Qdt.valid());
        auto const& t_Qdt = Maybe_Qdt.value();
        auto const& P_micro_m = get<P>(t_Qdt)();
        auto const& t_gmean_ij = get<gmean_ij>(t_Qdt)();
        auto const& t_gvar_ij = get<gvar_ij>(t_Qdt)();
        double n_samples =
            static_cast<double>(get<number_of_samples>(t_Qdt).value());
        auto sigma2_obs = get<Current_Noise>(dm_micro).value() * fs / n_samples;

        auto Maybe_step =
            micro_full_step_avg2_lifted(prior_probs, P_micro_m, t_gmean_ij, t_gvar_ij,
                                         sigma2_obs, y_obs);
        if (!Maybe_step) UNSCOPED_INFO(Maybe_step.error()());
        REQUIRE(Maybe_step.valid());
        auto step = std::move(Maybe_step.value());

        static_assert(var::is_derivative_v<T_probs>,
                      "test fixture must run with derivative-typed probs");
        prior_probs = var::outside_in(step.new_probs, prior_probs.dx());
    }

    if (!result) UNSCOPED_INFO(result.error()());
    REQUIRE(result.valid());
}
