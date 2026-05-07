// Parity test: calc_Qdt_eig vs calc_Qdt_taylor and calc_Qdtg_eig vs
// calc_Qdtg_taylor. The eig path is treated as ground truth (it's what the
// macro pipeline uses by default and what micro_R has historically agreed
// with). The Taylor path is the scaling-and-squaring fallback used by the
// micro_full lifted code when taylor_qdt_approximation=true (figure_2_micro).
//
// Each component of the returned Vector_Space is compared independently so a
// failure points at the actual quantity (P, gmean_i, gmean_ij, gvar_i,
// gvar_ij, gtotal_ij, gtotal_sqr_ij, gtotal_var_ij). dt is varied to force
// the scaling factor n ∈ {0, 1, 4, 7, 9} so the seed (get_Qn) and the
// composition (sum_Qn) are exercised at multiple depths.

#include <catch_amalgamated.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include <general_algorithm_on_containers.h>
#include <maybe_error.h>
#include <parameters.h>
#include <qmodel.h>

#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/interface/IModel.h>

#include <parameters_derivative.h>

using macrodr::P;
using macrodr::P_half;
using macrodr::Qdt;
using macrodr::Qdtg;
using macrodr::g;
using macrodr::gmean_i;
using macrodr::gmean_ij;
using macrodr::gsqr_i;
using macrodr::gtotal_ij;
using macrodr::gtotal_sqr_ij;
using macrodr::gtotal_var_ij;
using macrodr::gvar_i;
using macrodr::gvar_ij;
using macrodr::number_of_samples;

namespace {

// Compute the scaling factor n that calc_Qdt_taylor / calc_Qdtg_taylor pick
// for a given Qx and dt. Mirrors the formulas at qmodel.h:3094-3098 and
// :3127-:3134. Useful for diagnostic output so a failure includes how many
// squarings the Taylor path is meant to perform.
template <class C_Qx>
int taylor_scaling_n(const C_Qx& t_Qx, double dt, double half_factor) {
    auto v_Qrun = t_Qx() * (dt * half_factor);
    double max = maxAbs(var::primitive(v_Qrun));
    constexpr double desired = 0.125 / 4.0;
    if (max <= 0.0)
        return 0;
    int k = static_cast<int>(std::ceil(std::log2(max / desired)));
    return std::max(0, k);
}

template <class Qdt_t>
void check_qdt_components(Qdt_t const& qdt_eig, Qdt_t const& qdt_tay,
                           double rel_tol, double abs_tol, std::string const& tag) {
    auto report = [&](std::string const& field, auto const& a, auto const& b) {
        auto cmp = var::compare_contents(a, b, rel_tol, abs_tol, 5);
        if (!cmp) {
            UNSCOPED_INFO(tag << " — " << field << " mismatch:\n" << cmp.error()());
        }
        CHECK(cmp.valid());
    };

    report("number_of_samples",
           get<number_of_samples>(qdt_eig)(), get<number_of_samples>(qdt_tay)());
    report("P", get<P>(qdt_eig)(), get<P>(qdt_tay)());
    report("gmean_i", get<gmean_i>(qdt_eig)(), get<gmean_i>(qdt_tay)());
    report("gtotal_ij", get<gtotal_ij>(qdt_eig)(), get<gtotal_ij>(qdt_tay)());
    report("gmean_ij", get<gmean_ij>(qdt_eig)(), get<gmean_ij>(qdt_tay)());
    report("gsqr_i", get<gsqr_i>(qdt_eig)(), get<gsqr_i>(qdt_tay)());
    report("gvar_i", get<gvar_i>(qdt_eig)(), get<gvar_i>(qdt_tay)());
    report("gtotal_sqr_ij", get<gtotal_sqr_ij>(qdt_eig)(), get<gtotal_sqr_ij>(qdt_tay)());
    report("gtotal_var_ij", get<gtotal_var_ij>(qdt_eig)(), get<gtotal_var_ij>(qdt_tay)());
    report("gvar_ij", get<gvar_ij>(qdt_eig)(), get<gvar_ij>(qdt_tay)());
}

}  // namespace

TEST_CASE("macroir: calc_Qdt eig vs taylor parity",
          "[macroir][qdt][taylor][parity]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO("load_parameters failed (CWD must be build/<preset>/tests so ../data resolves): "
                                   << Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    auto par_values = par.standard_parameter();
    auto Maybe_model = (*model0)(par_values);
    REQUIRE(Maybe_model.valid());
    auto m = std::move(Maybe_model.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    // Use agonist=10 so the off-diagonal entries are non-trivial. With
    // scheme_CO_small_par.csv: on=0.1, off=100 → at agonist=10, max|Q| ≈ 100.
    constexpr double agonist = 10.0;
    auto t_Qx = build<Qx>(macro_dmr.calc_Qx(m, Agonist_concentration(agonist)));
    auto Maybe_eig = macro_dmr.calc_eigen(t_Qx);
    REQUIRE(Maybe_eig.valid());
    auto t_eig = std::move(Maybe_eig.value());

    // Set of dt values chosen to force n ∈ {0, 1, 4, 7, 9}. With max|Q| ≈ 100
    // and desired = 0.03125: n = ⌈log2(max·dt/0.03125)⌉.
    struct Case {
        double dt;
        std::size_t ns;
        const char* label;
    };
    const Case cases[] = {
        {1e-4, 5,   "dt=1e-4 (n=0)"},
        {1e-3, 50,  "dt=1e-3 (n≈1)"},
        {1e-2, 500, "dt=1e-2 (n≈4)"},
        {5e-2, 2500, "dt=5e-2 (n≈7)"},
        {2e-1, 10000, "dt=2e-1 (n≈9)"},
    };

    // Tolerance: eig path uses spectral integrals, Taylor path uses an order-6
    // truncation followed by 2ⁿ squarings. P is tightly constrained (machine
    // precision after squaring), the moment integrals (gtotal_*) accumulate
    // truncation error per squaring, and gvar_* are derived as differences of
    // similar quantities so they're the loosest. Pick a single relaxed tolerance
    // that still flags any qualitative bug.
    const double rel_tol = 1e-6;
    const double abs_tol = 1e-8;

    for (auto const& c : cases) {
        auto Maybe_qdt_eig = macro_dmr.calc_Qdt_eig(
            f_no_memoi, m, t_eig, number_of_samples(c.ns), c.dt);
        REQUIRE(Maybe_qdt_eig.valid());

        auto Maybe_qdt_tay = macro_dmr.calc_Qdt_taylor(
            m, t_Qx, number_of_samples(c.ns), c.dt);
        REQUIRE(Maybe_qdt_tay.valid());

        int n_full = taylor_scaling_n(t_Qx, c.dt, 1.0);
        int n_half = taylor_scaling_n(t_Qx, c.dt, 0.5);
        std::cerr << "[Qdt parity] " << c.label
                  << " — n_full=" << n_full << " n_half=" << n_half << "\n";

        check_qdt_components(Maybe_qdt_eig.value(), Maybe_qdt_tay.value(),
                             rel_tol, abs_tol, std::string{c.label});
    }
}

TEST_CASE("macroir: calc_Qdtg eig vs taylor parity (P_half)",
          "[macroir][qdtg][taylor][parity]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO("load_parameters failed (CWD must be build/<preset>/tests so ../data resolves): "
                                   << Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    auto par_values = par.standard_parameter();
    auto Maybe_model = (*model0)(par_values);
    REQUIRE(Maybe_model.valid());
    auto m = std::move(Maybe_model.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    constexpr double agonist = 10.0;
    auto t_Qx = build<Qx>(macro_dmr.calc_Qx(m, Agonist_concentration(agonist)));
    auto Maybe_eig = macro_dmr.calc_eigen(t_Qx);
    REQUIRE(Maybe_eig.valid());
    auto t_eig = std::move(Maybe_eig.value());

    struct Case {
        double dt;
        std::size_t ns;
        const char* label;
    };
    const Case cases[] = {
        {1e-4, 5,    "dt=1e-4 (n=0)"},
        {1e-3, 50,   "dt=1e-3 (n≈0)"},
        {1e-2, 500,  "dt=1e-2 (n≈3)"},
        {5e-2, 2500, "dt=5e-2 (n≈6)"},
        {2e-1, 10000, "dt=2e-1 (n≈8)"},
    };

    const double rel_tol = 1e-6;
    const double abs_tol = 1e-8;

    for (auto const& c : cases) {
        auto Maybe_g_eig = macro_dmr.calc_Qdtg_eig(
            f_no_memoi, m, t_eig, number_of_samples(c.ns), c.dt);
        REQUIRE(Maybe_g_eig.valid());

        auto Maybe_g_tay = macro_dmr.calc_Qdtg_taylor(
            m, t_Qx, number_of_samples(c.ns), c.dt);
        REQUIRE(Maybe_g_tay.valid());

        int n_half = taylor_scaling_n(t_Qx, c.dt, 0.5);
        std::cerr << "[Qdtg parity] " << c.label
                  << " — n_half=" << n_half << "\n";

        auto const& g_eig = Maybe_g_eig.value();
        auto const& g_tay = Maybe_g_tay.value();

        auto cmp_ph = var::compare_contents(get<P_half>(g_eig)(), get<P_half>(g_tay)(),
                                             rel_tol, abs_tol, 5);
        if (!cmp_ph) {
            UNSCOPED_INFO(c.label << " — P_half mismatch:\n" << cmp_ph.error()());
        }
        CHECK(cmp_ph.valid());
    }
}

// Same parity check as the first case, but on a *Derivative-aware* model. The
// figure_2_micro pipeline goes through dlikelihood_predictions →
// dlog_Likelihood_micro → calc_Qdt_taylor(d_m, …) — the Derivative-typed
// instantiation. The user reports that flipping taylor_qdt_approximation off
// (eig path) gives correct micro_MR/IR logL while the Taylor path (with seed
// fixed) still produces wrong logL. If the plain-double parity passes but the
// Derivative .value() parity fails, the bug is in how the new seed propagates
// through Derivative arithmetic, not in the math.
TEST_CASE("macroir: calc_Qdt eig vs taylor parity (Derivative)",
          "[macroir][qdt][taylor][parity][derivative]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO("load_parameters failed: " << Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    auto par_values = par.standard_parameter();
    auto par_transformed = par_values.to_transformed();
    auto dp = var::selfDerivative(par_transformed);
    auto dpp = dp.to_value();

    auto dmodel = cmd::load_dmodel("scheme_CO");
    if (!dmodel) UNSCOPED_INFO("load_dmodel failed: " << dmodel.error()());
    REQUIRE(dmodel.valid());
    auto model0_d = std::move(dmodel.value());

    auto maybe_dm = (*model0_d)(dpp);
    REQUIRE(maybe_dm.valid());
    auto dm = std::move(maybe_dm.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    // Cover both agonist=10 (active segments) and agonist=0 (relaxation
    // segments) since figure_2_micro alternates the two and any sensitivity
    // to clustered eigenvalues / zero rates would show only at one of them.
    const double agonists[] = {10.0, 0.0};

    struct Case {
        double dt;
        std::size_t ns;
        const char* label;
    };
    const Case cases[] = {
        {1e-4, 5,    "dt=1e-4 (n=0)"},
        {1e-3, 50,   "dt=1e-3"},
        {1e-2, 500,  "dt=1e-2"},
        {5e-2, 2500, "dt=5e-2"},
        {2e-1, 10000, "dt=2e-1"},
    };

    const double rel_tol = 1e-5;
    const double abs_tol = 1e-7;

    for (double agonist : agonists) {
        auto t_Qx_d = build<Qx>(macro_dmr.calc_Qx(dm, Agonist_concentration(agonist)));
        auto Maybe_eig_d = macro_dmr.calc_eigen(t_Qx_d);
        if (!Maybe_eig_d) {
            std::cerr << "[Qdt-D parity] agonist=" << agonist
                      << " — eigen decomposition failed (skipping eig path): "
                      << Maybe_eig_d.error()() << "\n";
            continue;
        }
        auto t_eig_d = std::move(Maybe_eig_d.value());

        for (auto const& c : cases) {
            auto Maybe_eig = macro_dmr.calc_Qdt_eig(
                f_no_memoi, dm, t_eig_d, number_of_samples(c.ns), c.dt);
            REQUIRE(Maybe_eig.valid());

            auto Maybe_tay = macro_dmr.calc_Qdt_taylor(
                dm, t_Qx_d, number_of_samples(c.ns), c.dt);
            REQUIRE(Maybe_tay.valid());

            std::cerr << "[Qdt-D parity] agonist=" << agonist
                      << " " << c.label << "\n";

            // Compare the *value* component of each derivative-typed Qdt field.
            auto const& q_eig = Maybe_eig.value();
            auto const& q_tay = Maybe_tay.value();

            auto report = [&](std::string const& field,
                              auto const& a_d, auto const& b_d) {
                auto cmp = var::compare_contents(
                    var::primitive(a_d), var::primitive(b_d),
                    rel_tol, abs_tol, 5);
                if (!cmp) {
                    UNSCOPED_INFO("agonist=" << agonist << " " << c.label
                                  << " — " << field << " mismatch:\n"
                                  << cmp.error()());
                }
                CHECK(cmp.valid());
            };

            report("P", get<P>(q_eig)(), get<P>(q_tay)());
            report("gmean_i", get<gmean_i>(q_eig)(), get<gmean_i>(q_tay)());
            report("gtotal_ij", get<gtotal_ij>(q_eig)(), get<gtotal_ij>(q_tay)());
            report("gmean_ij", get<gmean_ij>(q_eig)(), get<gmean_ij>(q_tay)());
            report("gsqr_i", get<gsqr_i>(q_eig)(), get<gsqr_i>(q_tay)());
            report("gvar_i", get<gvar_i>(q_eig)(), get<gvar_i>(q_tay)());
            report("gtotal_sqr_ij", get<gtotal_sqr_ij>(q_eig)(), get<gtotal_sqr_ij>(q_tay)());
            report("gtotal_var_ij", get<gtotal_var_ij>(q_eig)(), get<gtotal_var_ij>(q_tay)());
            report("gvar_ij", get<gvar_ij>(q_eig)(), get<gvar_ij>(q_tay)());
        }
    }
}

// Parity test for the Schur+Parlett Qdt path. Plain-double instantiation
// only at step 7 — Derivative parity lands in step 8 alongside the Derivative
// implementation in calc_Qdt_schur.
TEST_CASE("macroir: calc_Qdt eig vs schur parity",
          "[macroir][qdt][schur][parity]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    if (!Maybe_par) UNSCOPED_INFO("load_parameters failed: " << Maybe_par.error()());
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    auto par_values = par.standard_parameter();
    auto Maybe_model = (*model0)(par_values);
    REQUIRE(Maybe_model.valid());
    auto m = std::move(Maybe_model.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    constexpr double agonist = 10.0;
    auto t_Qx = build<Qx>(macro_dmr.calc_Qx(m, Agonist_concentration(agonist)));
    auto Maybe_eig = macro_dmr.calc_eigen(t_Qx);
    REQUIRE(Maybe_eig.valid());
    auto t_eig = std::move(Maybe_eig.value());

    struct Case {
        double dt;
        std::size_t ns;
        const char* label;
    };
    const Case cases[] = {
        {1e-4, 5,    "dt=1e-4"},
        {1e-3, 50,   "dt=1e-3"},
        {1e-2, 500,  "dt=1e-2"},
        {5e-2, 2500, "dt=5e-2"},
        {2e-1, 10000, "dt=2e-1"},
    };

    const double rel_tol = 1e-7;
    const double abs_tol = 1e-9;

    for (auto const& c : cases) {
        auto Maybe_qdt_eig = macro_dmr.calc_Qdt_eig(
            f_no_memoi, m, t_eig, number_of_samples(c.ns), c.dt);
        REQUIRE(Maybe_qdt_eig.valid());

        auto Maybe_qdt_sch = macro_dmr.calc_Qdt_schur(
            m, t_Qx, number_of_samples(c.ns), c.dt);
        if (!Maybe_qdt_sch) UNSCOPED_INFO("calc_Qdt_schur failed at " << c.label
                                          << ": " << Maybe_qdt_sch.error()());
        REQUIRE(Maybe_qdt_sch.valid());

        std::cerr << "[Qdt schur parity] " << c.label << "\n";
        check_qdt_components(Maybe_qdt_eig.value(), Maybe_qdt_sch.value(),
                              rel_tol, abs_tol, std::string{c.label});
    }
}

// Derivative-aware parity for the Schur Qdt path. With the Derivative-case
// delegation to calc_Qdt_taylor (which itself was validated against
// calc_Qdt_eig in the [macroir][parity][derivative] test), this should pass
// at the same tolerance as the Derivative Taylor parity. Re-running it here
// catches any regression in the delegation wiring.
TEST_CASE("macroir: calc_Qdt eig vs schur parity (Derivative)",
          "[macroir][qdt][schur][parity][derivative]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());

    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());

    auto par_values = par.standard_parameter();
    auto par_transformed = par_values.to_transformed();
    auto dp = var::selfDerivative(par_transformed);
    auto dpp = dp.to_value();

    auto dmodel = cmd::load_dmodel("scheme_CO");
    REQUIRE(dmodel.valid());
    auto model0_d = std::move(dmodel.value());
    auto maybe_dm = (*model0_d)(dpp);
    REQUIRE(maybe_dm.valid());
    auto dm = std::move(maybe_dm.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    constexpr double agonist = 10.0;
    auto t_Qx_d = build<Qx>(macro_dmr.calc_Qx(dm, Agonist_concentration(agonist)));
    auto Maybe_eig_d = macro_dmr.calc_eigen(t_Qx_d);
    REQUIRE(Maybe_eig_d.valid());
    auto t_eig_d = std::move(Maybe_eig_d.value());

    const double dts[] = {1e-4, 1e-3, 1e-2, 5e-2, 2e-1};
    const double rel_tol = 1e-5;
    const double abs_tol = 1e-7;

    for (double dt : dts) {
        const std::size_t ns = static_cast<std::size_t>(dt * 5e4);
        auto Maybe_eig = macro_dmr.calc_Qdt_eig(
            f_no_memoi, dm, t_eig_d, number_of_samples(ns), dt);
        REQUIRE(Maybe_eig.valid());

        auto Maybe_sch = macro_dmr.calc_Qdt_schur(
            dm, t_Qx_d, number_of_samples(ns), dt);
        if (!Maybe_sch) UNSCOPED_INFO("calc_Qdt_schur (D) failed at dt=" << dt
                                        << ": " << Maybe_sch.error()());
        REQUIRE(Maybe_sch.valid());

        std::cerr << "[Qdt-D schur parity] dt=" << dt << "\n";
        auto const& q_eig = Maybe_eig.value();
        auto const& q_sch = Maybe_sch.value();

        auto report = [&](std::string const& field, auto const& a_d, auto const& b_d) {
            auto cmp = var::compare_contents(var::primitive(a_d), var::primitive(b_d),
                                              rel_tol, abs_tol, 5);
            if (!cmp) UNSCOPED_INFO("dt=" << dt << " " << field << " mismatch:\n"
                                          << cmp.error()());
            CHECK(cmp.valid());
        };

        report("P", get<P>(q_eig)(), get<P>(q_sch)());
        report("gmean_i", get<gmean_i>(q_eig)(), get<gmean_i>(q_sch)());
        report("gtotal_ij", get<gtotal_ij>(q_eig)(), get<gtotal_ij>(q_sch)());
        report("gmean_ij", get<gmean_ij>(q_eig)(), get<gmean_ij>(q_sch)());
        report("gsqr_i", get<gsqr_i>(q_eig)(), get<gsqr_i>(q_sch)());
        report("gvar_i", get<gvar_i>(q_eig)(), get<gvar_i>(q_sch)());
        report("gtotal_sqr_ij", get<gtotal_sqr_ij>(q_eig)(), get<gtotal_sqr_ij>(q_sch)());
        report("gtotal_var_ij", get<gtotal_var_ij>(q_eig)(), get<gtotal_var_ij>(q_sch)());
        report("gvar_ij", get<gvar_ij>(q_eig)(), get<gvar_ij>(q_sch)());
    }
}

TEST_CASE("macroir: calc_Qdtg eig vs schur parity (P_half)",
          "[macroir][qdtg][schur][parity]") {
    using namespace macrodr;

    auto maybe_m = cmd::load_model("scheme_CO");
    REQUIRE(maybe_m.valid());
    auto model0 = std::move(maybe_m.value());
    auto Maybe_par = cmd::load_parameters(model0, "../data/scheme_CO_small_par.csv");
    REQUIRE(Maybe_par.valid());
    auto par = std::move(Maybe_par.value());
    auto par_values = par.standard_parameter();
    auto Maybe_model = (*model0)(par_values);
    REQUIRE(Maybe_model.valid());
    auto m = std::move(Maybe_model.value());

    auto f_no_memoi = var::create_empty_function_map();
    Macro_DMR macro_dmr{};

    constexpr double agonist = 10.0;
    auto t_Qx = build<Qx>(macro_dmr.calc_Qx(m, Agonist_concentration(agonist)));
    auto Maybe_eig = macro_dmr.calc_eigen(t_Qx);
    REQUIRE(Maybe_eig.valid());
    auto t_eig = std::move(Maybe_eig.value());

    const double dts[] = {1e-4, 1e-3, 1e-2, 5e-2, 2e-1};
    for (double dt : dts) {
        const std::size_t ns = static_cast<std::size_t>(dt * 5e4);
        auto Maybe_eig_g = macro_dmr.calc_Qdtg_eig(
            f_no_memoi, m, t_eig, number_of_samples(ns), dt);
        REQUIRE(Maybe_eig_g.valid());
        auto Maybe_sch_g = macro_dmr.calc_Qdtg_schur(
            m, t_Qx, number_of_samples(ns), dt);
        if (!Maybe_sch_g) UNSCOPED_INFO("calc_Qdtg_schur failed at dt=" << dt
                                         << ": " << Maybe_sch_g.error()());
        REQUIRE(Maybe_sch_g.valid());

        auto cmp_ph = var::compare_contents(get<P_half>(Maybe_eig_g.value())(),
                                              get<P_half>(Maybe_sch_g.value())(),
                                              1e-9, 1e-11, 5);
        if (!cmp_ph) UNSCOPED_INFO("dt=" << dt << " P_half mismatch:\n" << cmp_ph.error()());
        CHECK(cmp_ph.valid());
    }
}
