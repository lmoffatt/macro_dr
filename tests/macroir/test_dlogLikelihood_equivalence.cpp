// Regression test: dlogLikelihood (returns dMacro_State_Hessian_minimal) and
// dlogLikelihoodPredictions (returns dMacro_State_Ev_gradient_all) must produce
// bit-identical logL value AND bit-identical score (= derivative of logL).
//
// Both functions ultimately call the SAME templated Macro_DMR::log_Likelihood<>,
// differing only by the State template argument they pass. The State controls
// what additional slots get accumulated (Evolution_of, evaluation_time, ...)
// but the per-step logL accumulation in update_macro_state is identical:
//   get<logL>(t_prior_all)() += t_logL();
// runs unchanged regardless of which other slots the State carries.
//
// This test is the precondition for an optimisation in
// calculate_mnumerical_fisher_information (src/core/likelihood.cpp:855):
// that function does central-difference F via 2p calls to
// calculate_mdlikelihood_predictions, extracting only the score
// (derivative(get<logL>(state))). The full Evolution_of payload is computed
// and discarded inside each call — pure overhead. If this test passes,
// swapping the inner call to calculate_mdlikelihood (returns the slimmer
// dMacro_State_Hessian_minimal) is safe and ~30-50% faster per F.

#include <catch_amalgamated.hpp>

#include <experiment.h>
#include <maybe_error.h>
#include <parameters.h>
#include <parameters_derivative.h>
#include <qmodel.h>

#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/simulate.h>

#include <cmath>
#include <string>
#include <vector>

namespace {

// Reuse figure_2 canonical parameter values, lifted Current_Baseline to 1 to
// escape the Log10(0) degeneracy (same as test_gauss_newton.cpp).
inline var::Parameters_transformed build_theta_sim(
    var::Parameters_Transformations const& par_T) {
    auto theta = par_T.standard_parameter().to_transformed();
    auto vec   = theta();
    vec(0, 0ul) = std::log10(10.0);    // on
    vec(1, 0ul) = std::log10(100.0);   // off
    vec(2, 0ul) = std::log10(1.0);     // unitary_current
    vec(3, 0ul) = std::log10(1e-4);    // Current_Noise (figure_2 canonical)
    vec(4, 0ul) = std::log10(1.0);     // Current_Baseline (lift from 0)
    vec(5, 0ul) = std::log10(100.0);   // Num_ch_mean
    return theta.create(std::move(vec));
}

struct TestSetup {
    macrodr::cmd::ModelPtr                   model;
    var::Parameters_Transformations          par_T;
    macrodr::Experiment                       experiment;
    macrodr::Recording                        recording;
    macrodr::cmd::likelihood_algorithm_type   lik_algorithm;
};

Maybe_error<TestSetup> make_setup() {
    auto maybe_model = macrodr::cmd::load_model("scheme_CO");
    if (!maybe_model) return maybe_model.error();
    auto model = std::move(maybe_model.value());

    auto par_T = model->parameters_transformations();
    auto theta_sim_local = build_theta_sim(par_T);

    // figure_2 3-segment experiment (pre / agonist / post), interval_in_tau index 0.
    auto experiment = macrodr::cmd::create_experiment(
        {{2, 500, 0.0}, {4, 500, 10.0}, {4, 500, 0.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording(
        std::vector<double>(5000, 0.0));

    auto maybe_sim = macrodr::cmd::run_simulations(
        model, theta_sim_local, experiment, observations,
        std::string("uniformization"), 1000, 0);
    if (!maybe_sim) return maybe_sim.error();
    auto recording = get<macrodr::Recording>(maybe_sim.value()());

    auto maybe_lik = macrodr::cmd::build_likelihood_function(
        model, /*adaptive=*/ false, /*recursive=*/ true, /*averaging=*/ 2,
        /*variance=*/ true, /*taylor_variance_correction=*/ false,
        /*micro=*/ false, /*taylor_qdt=*/ false);
    if (!maybe_lik) return maybe_lik.error();

    return TestSetup{std::move(model), std::move(par_T),
                     std::move(experiment), std::move(recording),
                     std::move(maybe_lik.value())};
}

}  // namespace

TEST_CASE("calculate_mdlikelihood vs calculate_mdlikelihood_predictions: "
          "logL and score are bit-identical",
          "[likelihood][regression][equivalence]") {
    auto maybe_bundle = make_setup();
    if (!maybe_bundle) {
        UNSCOPED_INFO(maybe_bundle.error()());
    }
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();

    auto theta = bundle.par_T.standard_parameter().to_transformed();
    auto vec   = theta();
    vec(0, 0ul) = std::log10(10.0);
    vec(1, 0ul) = std::log10(100.0);
    vec(2, 0ul) = std::log10(1.0);
    vec(3, 0ul) = std::log10(1e-4);
    vec(4, 0ul) = std::log10(1.0);
    vec(5, 0ul) = std::log10(100.0);
    theta = theta.create(std::move(vec));

    // Slim path: returns dMacro_State_Hessian_minimal.
    auto maybe_state_min = macrodr::cmd::calculate_mdlikelihood(
        bundle.lik_algorithm, theta, bundle.experiment, bundle.recording);
    if (!maybe_state_min) {
        UNSCOPED_INFO("calculate_mdlikelihood failed: " << maybe_state_min.error()());
    }
    REQUIRE(maybe_state_min);
    auto const& state_min = maybe_state_min.value();

    // Full path: returns dMacro_State_Ev_gradient_all.
    auto maybe_state_full = macrodr::cmd::calculate_mdlikelihood_predictions(
        bundle.lik_algorithm, theta, bundle.experiment, bundle.recording);
    if (!maybe_state_full) {
        UNSCOPED_INFO("calculate_mdlikelihood_predictions failed: "
                      << maybe_state_full.error()());
    }
    REQUIRE(maybe_state_full);
    auto const& state_full = maybe_state_full.value();

    SECTION("logL value is bit-identical between the two functions") {
        const double logL_min  = primitive(get<logL>(state_min))();
        const double logL_full = primitive(get<logL>(state_full))();
        INFO("logL_min  = " << logL_min);
        INFO("logL_full = " << logL_full);
        INFO("diff      = " << (logL_min - logL_full));
        // Same Macro_DMR.log_Likelihood<> instantiation up to the State template
        // arg, which only controls accumulated-but-unread slots. Expect exact
        // floating-point equality (no reorderings).
        REQUIRE(logL_min == logL_full);
    }

    SECTION("score (gradient of logL) is bit-identical between the two functions") {
        auto const& score_min  = derivative(get<logL>(state_min))();
        auto const& score_full = derivative(get<logL>(state_full))();
        REQUIRE(score_min.nrows() == 6);
        REQUIRE(score_full.nrows() == 6);
        REQUIRE(score_min.ncols() == score_full.ncols());
        for (std::size_t i = 0; i < 6; ++i) {
            INFO("score[" << i << "]: min = " << score_min(i, 0ul)
                          << ", full = " << score_full(i, 0ul)
                          << ", diff = " << (score_min(i, 0ul) - score_full(i, 0ul)));
            REQUIRE(score_min(i, 0ul) == score_full(i, 0ul));
        }
    }

    SECTION("Repeated calls are deterministic (same θ + recording → same result)") {
        auto maybe_state_min2 = macrodr::cmd::calculate_mdlikelihood(
            bundle.lik_algorithm, theta, bundle.experiment, bundle.recording);
        REQUIRE(maybe_state_min2);
        const double logL1 = primitive(get<logL>(state_min))();
        const double logL2 = primitive(get<logL>(maybe_state_min2.value()))();
        INFO("logL call 1 = " << logL1 << ", logL call 2 = " << logL2);
        REQUIRE(logL1 == logL2);
    }
}
