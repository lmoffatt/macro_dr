// Scaffolding for MacroIR (test_der_macroir) derivative tests, mirroring
// the pattern embedded in legacy/qmodel.h around lines ~4601–4629.
// These tests are set up to run in isolation once the required inputs
// (prior state, Qdtm, model parameters, channel count, etc.) are provided.

#include <catch_amalgamated.hpp>

#include <derivative_test.h>   // var::test_derivative_clarke
#include <maybe_error.h>

#include <optional>

// Forward-looking context for wiring real MacroIR inputs.
// Fill this with concrete types (tt_prior, t_Qdtm, m, Nch, function set).
struct MacroIRTestContext {
    // TODO: add the concrete types for:
    //  - tt_prior (select<logL, ... , macror_algorithm>(t_prior))
    //  - t_Qdtm (result of calc_Qdtm / select<...>)
    //  - m (model parameters)
    //  - Nch (channel count)
    //  - f_no_memoi (function set without memoization)
};

static std::optional<MacroIRTestContext> build_macroir_test_context() {
    // Return a populated context once the appropriate variables are provided.
    return std::nullopt;
}

TEST_CASE("macroir: MacroIR derivative (Clarke) — scaffolding", "[macroir][derivatives][macroir]") {
    auto ctx = build_macroir_test_context();
    if (!ctx) {
        INFO("MacroIR derivative test is scaffolded; provide MacroIRTestContext to enable.");
        // Intentionally exit without assertions so CI remains green until wired.
        return;
    }

    // Example pattern, mirroring legacy/qmodel.h around lines 4609–4624:
    //
    // const auto h = 1e-7;
    // auto result = var::test_derivative_clarke<false>(
    //     [&](auto l_t_prior, auto const& l_Qdtm, auto const& l_m, auto const& l_Nch)
    //         -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(l_t_prior)>,
    //             var::Vector_Space<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL, eplogL, vplogL>>> {
    //         auto Maybe_res = MacroR2<recursive, averaging, variance, variance_correction>{}(
    //             f_no_memoi, std::move(l_t_prior), l_Qdtm, l_m, l_Nch, y()[i_step], fs);
    //         if (!Maybe_res) return Maybe_res.error();
    //         return select<logL, elogL, vlogL, P_mean, P_Cov, y_mean, y_var, plogL, eplogL, vplogL>(
    //             std::move(Maybe_res.value()));
    //     },
    //     h,
    //     tt_prior, t_Qdtm, m, Nch);
    // REQUIRE(result);
    //
    // Once concrete types and builders are provided, adapt the snippet accordingly
    // and enable the actual derivative checks.
}

