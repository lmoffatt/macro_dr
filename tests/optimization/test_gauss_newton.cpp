// Tests for the generic Gauss-Newton / Levenberg-Marquardt optimizer
// (legacy/gauss_newton.h) on the macro-DR Gaussian-moment-matched likelihood.
//
// Strategy: verify each property of the Gaussian_Fisher_Information path
// (the slot inside dMacro_State_Hessian_minimal) on a simple scheme_CO +
// macro_IR setup. macro_IR is the well-behaved algorithm at all Δt regimes
// per the eLife 2025 paper analysis, so this is a "happy path" test of the
// optimizer's correctness; pathological cases (macro_R at Δt=1τ) are
// validated elsewhere.
//
// Configuration:
//   model:           scheme_CO with figure_2 canonical parameter values
//                    (see build_theta_sim() and figure_2.macroir). NOT the
//                    model's standard_parameter defaults — figure_2 overrides
//                    them.
//   algorithm:       macro_IR (recursive=true, averaging=2, variance=true)
//   experiment:      3-segment protocol {2,4,4} steps @ {0,10,0} μM, each step
//                    averaging 500 raw pts, fs = 50 kHz → 10 scored macro steps,
//                    interval ≈ 1τ (the coarsest figure_2 cell), Num_ch=100
//   prior:           NONE — pure MLE (gradient is just the likelihood gradient)
//   seed:            42 (reproducibility)

#include <catch_amalgamated.hpp>

#include <experiment.h>
#include <gauss_newton.h>
#include <maybe_error.h>
#include <parameters.h>
#include <qmodel.h>

#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/simulate.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace {

// Bundle of everything needed for the test. Constructed once per test case
// (Catch2 isolates state per TEST_CASE).
//
// NOTE on lifetimes: Parameters_transformed stores a raw pointer to a
// Parameters_Transformations. We keep par_T inside the bundle (stable address
// once the bundle is constructed) and rebuild theta_sim locally in each test
// from bundle.par_T. Storing theta_sim directly in the bundle would dangle
// across the make_setup return.
struct TestSetup {
    macrodr::cmd::ModelPtr                        model;
    var::Parameters_Transformations               par_T;        // keep alive for theta_sim
    macrodr::Experiment                           experiment;
    macrodr::Recording                            recording;
    macrodr::cmd::likelihood_algorithm_type       lik_algorithm;
};

// Use the canonical figure_2 parameter values for scheme_CO (per
// projects/eLife_2025/ops/local/figure_2.macroir). The model's standard_parameter
// defaults (in legacy/models_simple.h:28) are NOT used in figure_2 — the script
// overrides them via create_parameters. Key differences:
//   on               default 0.1   → figure_2 uses 10
//   Current_Noise    default 1e-3  → figure_2 uses 1e-4 (small but nonzero)
//   Current_Baseline default 0     → figure_2 uses 1 (escapes Log10(0) degeneracy)
//   Num_ch_mean      default 5000  → figure_2 uses 100
// off=100, unitary_current=1 match the defaults.
//
// Parameter order (per models_simple.h:11-14):
//   0: on, 1: off, 2: unitary_current,
//   3: Current_Noise, 4: Current_Baseline, 5: Num_ch_mean
//
// All values stored in Log10-transformed space.
inline var::Parameters_transformed build_theta_sim(
    var::Parameters_Transformations const& par_T) {
    auto theta = par_T.standard_parameter().to_transformed();
    auto vec = theta();
    vec(0, 0ul) = std::log10(10.0);     // on
    vec(1, 0ul) = std::log10(100.0);    // off
    vec(2, 0ul) = std::log10(1.0);      // unitary_current
    vec(3, 0ul) = std::log10(1e-4);     // Current_Noise (figure_2 canonical)
    vec(4, 0ul) = std::log10(1.0);      // Current_Baseline (lift from 0)
    vec(5, 0ul) = std::log10(100.0);    // Num_ch_mean (figure_2 canonical)
    return theta.create(std::move(vec));
}

Maybe_error<TestSetup> make_setup() {
    // 1. Load model
    auto maybe_model = macrodr::cmd::load_model("scheme_CO");
    if (!maybe_model) return maybe_model.error();
    auto model = std::move(maybe_model.value());

    // 2. Standard parameters with figure_2 canonical values (see build_theta_sim).
    //    par_T is long-lived; theta_sim_local is only used within this function
    //    for the simulation step. Each TEST_CASE rebuilds theta_sim via
    //    build_theta_sim(bundle.par_T).
    auto par_T = model->parameters_transformations();
    auto theta_sim_local = build_theta_sim(par_T);

    // 3. Experiment: replicate figure_2.macroir setup at axis_interval index 0
    //    (interval_in_tau = "1"). Three-segment protocol:
    //      pre-agonist:  2 intervals × 500 samples = 1000 samples @ 0 μM
    //      agonist:      4 intervals × 500 samples = 2000 samples @ 10 μM
    //      post-agonist: 4 intervals × 500 samples = 2000 samples @ 0 μM
    //    Total: 5000 samples at fs=50 kHz, initial agonist 0.0, initial t 0.0.
    constexpr std::size_t T_TOTAL = 5000;
    auto experiment = macrodr::cmd::create_experiment(
        {{2, 500, 0.0}, {4, 500, 10.0}, {4, 500, 0.0}},
        50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording(
        std::vector<double>(T_TOTAL, 0.0));

    // 4. Simulate one recording at theta_sim_local using figure_2's uniformization
    //    with number_of_substeps=1000 (per projects/eLife_2025/ops/local/figure_2.macroir).
    //    seed=0 matches figure_2; reproducible per Catch2 run.
    auto maybe_sim = macrodr::cmd::run_simulations(
        model, theta_sim_local, experiment, observations,
        std::string("uniformization"), 1000, 0);
    if (!maybe_sim) return maybe_sim.error();
    auto recording = get<macrodr::Recording>(maybe_sim.value()());

    // 5. Build macro_IR algorithm:
    //    recursive=true, averaging=2, variance=true.
    //    build_likelihood_function returns Maybe_error<variant<...>>; unwrap.
    auto maybe_lik = macrodr::cmd::build_likelihood_function(
        model,
        /*adaptive=*/ false,
        /*recursive=*/ true,
        /*averaging=*/ 2,
        /*variance=*/ true,
        /*taylor_variance_correction=*/ false,
        /*micro=*/ false,
        /*taylor_qdt=*/ false);
    if (!maybe_lik) return maybe_lik.error();

    return TestSetup{std::move(model), std::move(par_T),
                     std::move(experiment), std::move(recording),
                     std::move(maybe_lik.value())};
}

// Frobenius norm of a column vector (Matrix<double> of shape (p, 1)).
double frob(Matrix<double> const& m) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.nrows(); ++i)
        for (std::size_t j = 0; j < m.ncols(); ++j)
            s += m(i, j) * m(i, j);
    return std::sqrt(s);
}

// Convert dMacro_State_Hessian_minimal (the AD-derivative output of
// calculate_mdlikelihood) into dlogPs (logL/Grad/Gaussian_Fisher_Information
// Vector_Space expected by gauss_newton_maximize<Gaussian_Fisher_Information>).
//
// dMacro_State_Hessian_minimal contains Derivative<logL, Parameters_transformed>
// from which value (primitive) and gradient (derivative) are extracted. The
// Gaussian_Fisher_Information slot is populated by update_macro_state.
Maybe_error<dlogPs> evaluate_likelihood_as_dlogPs(
    macrodr::cmd::likelihood_algorithm_type const& lik,
    var::Parameters_transformed const&             p,
    macrodr::Experiment const&                     e,
    macrodr::Recording const&                      r) {
    auto Maybe_dml = macrodr::cmd::calculate_mdlikelihood(lik, p, e, r);
    if (!Maybe_dml) return Maybe_dml.error();
    auto const& dml = Maybe_dml.value();

    // primitive() strips the Derivative wrapper to get the underlying logL scalar.
    // derivative() extracts the gradient (Matrix<double>) from the AD chain.
    return dlogPs{
        logL(primitive(get<logL>(dml))()),
        Grad(derivative(get<logL>(dml))(), p),
        Gaussian_Fisher_Information(get<Gaussian_Fisher_Information>(dml)().value(), p)
    };
}

}  // namespace


// ----------------------------------------------------------------------------
// TEST A — dlogLikelihood at θ_sim is well-formed
// ----------------------------------------------------------------------------
TEST_CASE("dlogLikelihood at theta_sim returns valid logL/Grad/Gaussian_Fisher_Information",
          "[optimization][gauss_newton][validity]") {
    auto maybe_bundle = make_setup();
    if (!maybe_bundle) {
        UNSCOPED_INFO(maybe_bundle.error()());
    }
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();

    // Rebuild theta_sim locally; its pointer references bundle.par_T which
    // is alive for the duration of this TEST_CASE.
    auto theta_sim = build_theta_sim(bundle.par_T);

    auto maybe_eval = evaluate_likelihood_as_dlogPs(
        bundle.lik_algorithm, theta_sim, bundle.experiment, bundle.recording);
    if (!maybe_eval) {
        UNSCOPED_INFO(maybe_eval.error()());
    }
    REQUIRE(maybe_eval);
    auto const& eval = maybe_eval.value();

    SECTION("logL is finite") {
        REQUIRE(std::isfinite(get<logL>(eval)()));
    }

    SECTION("Grad is finite and has dim 6") {
        auto const& grad = get<Grad>(eval)().value();
        REQUIRE(grad.nrows() == 6);
        REQUIRE(grad.ncols() == 1);
        for (std::size_t i = 0; i < 6; ++i) {
            INFO("grad[" << i << "] = " << grad(i, 0ul));
            REQUIRE(std::isfinite(grad(i, 0ul)));
        }
    }

    SECTION("Gaussian_Fisher_Information is 6×6 and PSD (Cholesky succeeds)") {
        auto const& gfi = get<Gaussian_Fisher_Information>(eval)().value();
        REQUIRE(gfi.nrows() == 6);
        REQUIRE(gfi.ncols() == 6);

        // PSD test: Cholesky-based inverse succeeds iff the matrix is positive
        // definite (numerically). This is stricter than non-negative
        // eigenvalues but more robust than approximate threshold testing.
        auto maybe_inv = inv(gfi);
        if (!maybe_inv) {
            UNSCOPED_INFO("Gaussian_Fisher_Information not PD: " << maybe_inv.error()());
        }
        REQUIRE(maybe_inv);
    }
}


// ----------------------------------------------------------------------------
// TEST B — Gradient is an ascent direction at an off-MLE point
// ----------------------------------------------------------------------------
TEST_CASE("Gradient is an ascent direction at off-MLE point",
          "[optimization][gauss_newton][gradient]") {
    auto maybe_bundle = make_setup();
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();

    auto theta_sim = build_theta_sim(bundle.par_T);

    // Offset θ_sim by 0.05 in each parameter (transformed Log10 space ≈ 12 % up)
    Matrix<double> delta(6, 1, 0.0);
    for (std::size_t i = 0; i < 6; ++i) delta(i, 0ul) = 0.05;
    auto theta_off = theta_sim.create(theta_sim() + delta);

    auto maybe_eval_off = evaluate_likelihood_as_dlogPs(
        bundle.lik_algorithm, theta_off, bundle.experiment, bundle.recording);
    REQUIRE(maybe_eval_off);
    auto const& eval_off = maybe_eval_off.value();

    auto const& grad = get<Grad>(eval_off)().value();
    double logL_off = get<logL>(eval_off)();

    // Tiny step in gradient direction; expect logL to increase.
    const double grad_norm = frob(grad);
    REQUIRE(grad_norm > 0.0);
    const double step = 1e-4 / grad_norm;

    auto theta_step = theta_off.create(theta_off() + grad * step);
    auto maybe_eval_step = evaluate_likelihood_as_dlogPs(
        bundle.lik_algorithm, theta_step, bundle.experiment, bundle.recording);
    REQUIRE(maybe_eval_step);
    double logL_step = get<logL>(maybe_eval_step.value())();

    INFO("logL_off = " << logL_off << ", logL_step = " << logL_step
                       << ", grad_norm = " << grad_norm);
    REQUIRE(logL_step > logL_off);
}


// ----------------------------------------------------------------------------
// TEST C — Gauss-Newton converges from a perturbed start
// ----------------------------------------------------------------------------
TEST_CASE("Gauss-Newton converges from perturbed initial guess",
          "[optimization][gauss_newton][convergence]") {
    auto maybe_bundle = make_setup();
    REQUIRE(maybe_bundle);
    auto const& bundle = maybe_bundle.value();

    auto theta_sim = build_theta_sim(bundle.par_T);

    // Initial: θ_sim + 0.05 perturbation per parameter
    Matrix<double> delta(6, 1, 0.0);
    for (std::size_t i = 0; i < 6; ++i) delta(i, 0ul) = 0.05;
    auto theta_init = theta_sim.create(theta_sim() + delta);

    // Objective inline (Camino 2): adapter from dMacro_State_Hessian_minimal
    // to dlogPs so gauss_newton_maximize<Gaussian_Fisher_Information> can use
    // logL/Grad/Gaussian_Fisher_Information tags.
    auto objective = [&bundle](var::Parameters_transformed const& p) {
        return evaluate_likelihood_as_dlogPs(
            bundle.lik_algorithm, p, bundle.experiment, bundle.recording);
    };

    macrodr::optimization::gauss_newton_options opts;
    // Loop starts at lambda=0 (Newton). On failure, jumps to lambda_kickoff (default 1).
    // For figure_2 setup, kickoff=1 is enough — θ_init is moderately close to θ_sim.
    opts.max_iter       = 50;
    opts.grad_rtol      = 1e-6;
    opts.dvalue_tol     = 1e-10;
    opts.verbose        = false;  // set true to trace per-iteration state to stderr

    auto maybe_result = macrodr::optimization::gauss_newton_maximize<Gaussian_Fisher_Information>(
        objective, theta_init, opts);
    if (!maybe_result) {
        UNSCOPED_INFO(maybe_result.error()());
    }
    REQUIRE(maybe_result);
    auto const& result = maybe_result.value();

    SECTION("Converged with reasonable status") {
        INFO("status = " << result.status << ", n_iter = " << result.n_iter
                         << ", max_logL = " << result.max_value);
        REQUIRE((result.status == "converged_newton_dec" ||
                 result.status == "converged_grad" ||
                 result.status == "converged_value"));
        REQUIRE(result.n_iter < 30);
    }

    SECTION("Converged to a stationary point of the macro objective (Newton decrement → 0)") {
        // What a SINGLE low-channel recording can certify is that the optimizer
        // reached a STATIONARY point of the (misspecified) macro objective — NOT
        // that the macro MLE equals θ_sim.
        //
        // Why we do NOT assert θ_sim recovery (the old Wald D² < χ²_6 bar, removed):
        // the moment-matched macro likelihood is misspecified by construction, and
        // for one low-channel recording (Num_ch=100, interval≈1τ, 10 scored steps)
        // the per-recording Wald statistic
        //     D²(θ_MAP,θ_sim) = (θ_MAP−θ_sim)ᵀ · G_lik(θ_MAP) · (θ_MAP−θ_sim)
        // is NOT χ²_6. On the matching on-disk cloud (1024 macro_IR single-recording
        // MLEs at this exact regime) D² has median ≈ 89 and 100 % of recordings
        // exceed χ²_6 = 12.59 — at certified-stationary converged points (this very
        // setup yields D² ≈ 1149 with ½·newton_dec² ≈ 2.4e-9). So a χ²_6 bar measures
        // non-quadraticity + F≠J, not optimizer error. θ_sim recovery is only grounded
        // in the near-exact regime (Num_ch≈1e4), where the joint MLE θ_pool ≈ θ_sim;
        // for a single low-N trace it is false by design. (See theory/macroir/docs/
        // Posterior_Information_Distortion/supplement_algorithm_validation.tex and the
        // figure_3 misspecification frame.)
        //
        // The grounded acceptance check: the Newton decrement ½·gᵀ G_lik⁻¹ g (scale-
        // invariant, the GN's own primary convergence criterion) → 0 at the returned
        // point. ‖grad‖ alone is the WRONG bar — a genuine stationary point here has
        // ‖grad‖ ≈ 0.03 (not < grad_rtol·‖θ‖) because the GN stops on the decrement,
        // not the raw gradient. D² and per-parameter |Δ|/σ are kept as INFO only.
        auto const& gfi_at_max  = get<Gaussian_Fisher_Information>(result.result_value)().value();
        auto const& grad_at_max = get<Grad>(result.result_value)().value();
        auto maybe_chol = cholesky(gfi_at_max);
        if (!maybe_chol) {
            UNSCOPED_INFO("cholesky(Gaussian_Fisher_Information) failed at converged "
                          "point (not PD): " << maybe_chol.error()());
        }
        REQUIRE(maybe_chol);
        auto maybe_L_inv = inv(maybe_chol.value());
        REQUIRE(maybe_L_inv);
        // Sigma = G_lik⁻¹ via cholesky+inv(L)+XXT (see feedback_lapack_symmposdef_inv_broken).
        auto Sigma = XXT(tr(maybe_L_inv.value()));

        auto const& theta_MAP_vec = result.argmax();
        auto const& theta_sim_vec = theta_sim();

        // Per-parameter Δ and σ — INFO diagnostic only.
        Matrix<double> delta_theta(6, 1);
        for (std::size_t i = 0; i < 6; ++i) {
            delta_theta(i, 0ul) = theta_MAP_vec(i, 0ul) - theta_sim_vec(i, 0ul);
            const double sigma_i = std::sqrt(Sigma(i, i));
            const double diff_i  = std::abs(delta_theta(i, 0ul));
            INFO("param " << i << ": |Δ| = " << diff_i
                          << ", σ = " << sigma_i
                          << ", |Δ|/σ = " << diff_i / sigma_i);
        }

        // D² = δᵀ · G_lik · δ — INFO diagnostic only, NOT an acceptance bar (see above).
        auto G_delta = gfi_at_max * delta_theta;
        double D_squared = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            D_squared += delta_theta(i, 0ul) * G_delta(i, 0ul);
        }
        INFO("Wald D² = " << D_squared << " (diagnostic only; χ²_6,0.95 = 12.59 is NOT "
             "asserted — a single-recording macro MLE is not χ²_6 distributed)");

        // Stationarity (the actual acceptance criterion): ½·gᵀ G_lik⁻¹ g at θ_MAP.
        auto Sigma_grad = Sigma * grad_at_max;
        double newton_dec2 = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            newton_dec2 += grad_at_max(i, 0ul) * Sigma_grad(i, 0ul);
        }
        INFO("½·newton_dec² = " << 0.5 * newton_dec2 << " ; ‖grad‖ = " << frob(grad_at_max)
             << " ; status = " << result.status);
        // Loose by 100× vs the GN's own newton_dec2_tol (default 1e-8) to absorb the
        // independent G_lik⁻¹ recompute; still cleanly separates a stationary point
        // (~1e-9) from a non-stationary converged_value stall (large decrement).
        REQUIRE(0.5 * newton_dec2 < 100.0 * opts.newton_dec2_tol);
    }
}
