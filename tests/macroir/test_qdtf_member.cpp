// Regression anchor for the Bessel member (legacy/qdtf_member.h): in the box
// configuration (n_poles = 0) it computes the same likelihood as the
// existing av=2 recursive member, on the same simulated recording.
//
// Exact equality is not expected: the established member applies Bayesian
// shrinkage to gmean_ij (pseudo-count ε_mach·κ_V, qmodel.h:1694) and a
// simplex trust coefficient on the mean update, neither of which the qdtf
// engine carries. Both effects are tiny on a well-conditioned two-state
// problem, so the gate is a tight relative tolerance on the total logL,
// with the absolute difference reported for the record.
//
// The mathematical chain behind this gate: engine tables ≡ quadrature and
// ≡ the classic Ee/E3 route (tests/math/test_qdtf_engine.cpp), the
// recursion ≡ the boundary-state assembly and ≡ Monte Carlo of the exact
// filtered process (theory/macroir/notes/bessel_reference/).

#include <catch_amalgamated.hpp>

#include <experiment.h>
#include <maybe_error.h>
#include <parameters.h>
#include <qmodel.h>

#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/qdtf_likelihood.h>
#include <macrodr/cmd/simulate.h>

#include <cmath>
#include <string>
#include <vector>

namespace {

inline var::Parameters_transformed build_theta_sim(var::Parameters_Transformations const& par_T) {
    auto theta = par_T.standard_parameter().to_transformed();
    auto vec = theta();
    vec(0, 0ul) = std::log10(10.0);   // on
    vec(1, 0ul) = std::log10(100.0);  // off
    vec(2, 0ul) = std::log10(1.0);    // unitary_current
    vec(3, 0ul) = std::log10(1e-4);   // Current_Noise
    vec(4, 0ul) = std::log10(1.0);    // Current_Baseline
    vec(5, 0ul) = std::log10(100.0);  // Num_ch_mean
    return theta.create(std::move(vec));
}

// par_T lives here because Parameters_transformed keeps a raw pointer to its
// Parameters_Transformations and IModel::parameters_transformations() returns
// by value: theta must be built from bundle.par_T once the bundle is in place
// (same arrangement as test_dlogLikelihood_equivalence.cpp).
struct Setup {
    macrodr::cmd::ModelPtr model;
    var::Parameters_Transformations par_T;
    macrodr::Experiment experiment;
    macrodr::Recording recording;
};

Maybe_error<Setup> make_setup() {
    auto maybe_model = macrodr::cmd::load_model("scheme_CO");
    if (!maybe_model)
        return maybe_model.error();
    auto model = std::move(maybe_model.value());

    auto par_T = model->parameters_transformations();
    auto theta_sim_local = build_theta_sim(par_T);

    auto experiment = macrodr::cmd::create_experiment(
        {{2, 500, 0.0}, {4, 500, 10.0}, {4, 500, 0.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording(std::vector<double>(5000, 0.0));

    auto maybe_sim = macrodr::cmd::run_simulations(model, theta_sim_local, experiment,
                                                   observations, std::string("uniformization"),
                                                   1000, 0);
    if (!maybe_sim)
        return maybe_sim.error();
    auto recording = get<macrodr::Recording>(maybe_sim.value()());

    return Setup{std::move(model), std::move(par_T), std::move(experiment),
                 std::move(recording)};
}

}  // namespace

TEST_CASE("qdtf box configuration matches the established av=2 member",
          "[qdtf][member][regression]") {
    auto maybe_bundle = make_setup();
    if (!maybe_bundle) {
        UNSCOPED_INFO(maybe_bundle.error()());
    }
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();
    const auto theta = build_theta_sim(bundle.par_T);

    auto maybe_lik = macrodr::cmd::build_likelihood_function(
        bundle.model, /*adaptive=*/false, /*recursive=*/true, /*averaging=*/2,
        /*variance=*/true, /*taylor_variance_correction=*/false, /*micro=*/false,
        /*taylor_qdt=*/false);
    REQUIRE(maybe_lik);

    auto maybe_ref = macrodr::cmd::calculate_mlikelihood(maybe_lik.value(), theta,
                                                         bundle.experiment, bundle.recording);
    if (!maybe_ref) {
        UNSCOPED_INFO("reference member failed: " << maybe_ref.error()());
    }
    REQUIRE(maybe_ref);
    const double logL_ref = get<logL>(maybe_ref.value())();

    auto maybe_box = macrodr::cmd::calculate_qdtf_likelihood(
        bundle.model, theta, bundle.experiment, bundle.recording, 0, 0.0);
    if (!maybe_box) {
        UNSCOPED_INFO("qdtf box failed: " << maybe_box.error()());
    }
    REQUIRE(maybe_box);
    const double logL_box = get<logL>(maybe_box.value())();

    INFO("logL reference (av=2)   = " << logL_ref);
    INFO("logL qdtf box           = " << logL_box);
    INFO("difference              = " << (logL_box - logL_ref));
    CHECK(std::abs(logL_box - logL_ref) < 1e-4 * std::abs(logL_ref));
}

TEST_CASE("qdtf Bessel configuration runs and is deterministic", "[qdtf][member]") {
    auto maybe_bundle = make_setup();
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();
    const auto theta = build_theta_sim(bundle.par_T);

    auto maybe_b4 = macrodr::cmd::calculate_qdtf_likelihood(
        bundle.model, theta, bundle.experiment, bundle.recording, 4, 10e3);
    if (!maybe_b4) {
        UNSCOPED_INFO("qdtf bessel-4 failed: " << maybe_b4.error()());
    }
    REQUIRE(maybe_b4);
    const double l1 = get<logL>(maybe_b4.value())();
    CHECK(std::isfinite(l1));

    auto maybe_b4b = macrodr::cmd::calculate_qdtf_likelihood(
        bundle.model, theta, bundle.experiment, bundle.recording, 4, 10e3);
    REQUIRE(maybe_b4b);
    CHECK(get<logL>(maybe_b4b.value())() == l1);

    // 8-pole variant constructs and runs too
    auto maybe_b8 = macrodr::cmd::calculate_qdtf_likelihood(
        bundle.model, theta, bundle.experiment, bundle.recording, 8, 10e3);
    REQUIRE(maybe_b8);
    CHECK(std::isfinite(get<logL>(maybe_b8.value())()));
}
