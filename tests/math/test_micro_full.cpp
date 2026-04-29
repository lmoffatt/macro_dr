#include "catch_amalgamated.hpp"

#include <cmath>
#include <cstddef>

#include "matrix.h"
#include "micro_full.h"
#include "qmodel.h"

using namespace macrodr;

namespace {

Patch_State make_patch_state(P_mean pm, P_Cov pc) {
    Patch_State ps;
    get<P_mean>(ps()) = std::move(pm);
    get<P_Cov>(ps()) = std::move(pc);
    return ps;
}

Micro_Patch_State make_micro_patch_state(Micro_P_mean pm) {
    Micro_Patch_State mps;
    get<Micro_P_mean>(mps()) = std::move(pm);
    return mps;
}

// Build a Patch_State corresponding to an i.i.d. multinomial with probability p0 for k=2.
// For k=2, P_Cov = (diag(p) - p p^T) / N.
Patch_State multinomial_patch(double p0, std::size_t N) {
    Matrix<double> mean(1, 2);
    mean[0] = p0;
    mean[1] = 1.0 - p0;

    SymmetricMatrix<double> cov(2, 2);
    double inv_N = 1.0 / static_cast<double>(N);
    cov.set(0, 0, (p0 - p0 * p0) * inv_N);
    cov.set(1, 0, (-p0 * (1.0 - p0)) * inv_N);
    cov.set(1, 1, ((1.0 - p0) - (1.0 - p0) * (1.0 - p0)) * inv_N);

    return make_patch_state(P_mean(std::move(mean)), P_Cov(std::move(cov)));
}

}  // namespace

TEST_CASE("Micro_state_Num_ch enumeration", "[micro_full][enumerate]") {
    SECTION("k=2, N=10 produces N+1 rows each summing to N") {
        auto micro = create_Micro_state_Num_ch(10, 2);
        REQUIRE(micro().nrows() == 11);
        REQUIRE(micro().ncols() == 2);

        for (std::size_t i = 0; i < micro().nrows(); ++i) {
            std::size_t row_sum = micro()(i, std::size_t{0}) + micro()(i, std::size_t{1});
            CHECK(row_sum == 10);
        }
    }

    SECTION("k=3, N=4 produces C(6,2)=15 rows each summing to N") {
        auto micro = create_Micro_state_Num_ch(4, 3);
        REQUIRE(micro().nrows() == 15);
        REQUIRE(micro().ncols() == 3);

        for (std::size_t i = 0; i < micro().nrows(); ++i) {
            std::size_t row_sum = 0;
            for (std::size_t j = 0; j < 3; ++j) {
                row_sum += micro()(i, j);
            }
            CHECK(row_sum == 4);
        }
    }

    SECTION("N_channels_of reads N from the enumeration") {
        auto micro = create_Micro_state_Num_ch(100, 4);
        CHECK(N_channels_of(micro) == 100);
    }

    SECTION("num_full_states_of matches C(N+k-1, k-1)") {
        CHECK(num_full_states_of(10, 2) == 11);
        CHECK(num_full_states_of(4, 3) == 15);
        CHECK(num_full_states_of(100, 5) == 4598126);
    }
}

TEST_CASE("lift_Macro_to_Micro produces a normalized distribution", "[micro_full][lift]") {
    std::size_t N = 50;
    auto micro = create_Micro_state_Num_ch(N, 2);
    auto patch = multinomial_patch(0.3, N);

    auto maybe_mps = lift_Macro_to_Micro(micro, patch);
    REQUIRE(maybe_mps.valid());
    auto mps = std::move(maybe_mps.value());

    auto const& prob = get<Micro_P_mean>(mps())();
    REQUIRE(prob.size() == N + 1);

    double sum = 0.0;
    for (std::size_t i = 0; i < prob.size(); ++i) {
        sum += prob[i];
    }
    CHECK(std::abs(sum - 1.0) < 1e-12);

    for (std::size_t i = 0; i < prob.size(); ++i) {
        CHECK(prob[i] >= 0.0);
    }
}

TEST_CASE("round-trip project(lift(·)) recovers 2-moment input within O(1/N)",
          "[micro_full][round-trip]") {
    SECTION("k=2, N=100, multinomial moments") {
        std::size_t N = 100;
        auto micro = create_Micro_state_Num_ch(N, 2);
        auto patch = multinomial_patch(0.3, N);

        auto maybe_mps = lift_Macro_to_Micro(micro, patch);
        REQUIRE(maybe_mps.valid());
        auto mps = std::move(maybe_mps.value());

        auto recovered = project_Micro_to_Macro(micro, mps);
        auto const& recovered_mean = get<P_mean>(recovered())();
        auto const& recovered_cov = get<P_Cov>(recovered())();
        auto const& original_cov = get<P_Cov>(patch())();

        CHECK(std::abs(recovered_mean[0] - 0.3) < 1e-2);
        CHECK(std::abs(recovered_mean[1] - 0.7) < 1e-2);

        CHECK(std::abs(recovered_cov(0, 0) - original_cov(0, 0)) < 1e-3);
        CHECK(std::abs(recovered_cov(1, 1) - original_cov(1, 1)) < 1e-3);
        CHECK(std::abs(recovered_cov(1, 0) - original_cov(1, 0)) < 1e-3);
    }

    SECTION("k=2, N=1000: tighter recovery") {
        std::size_t N = 1000;
        auto micro = create_Micro_state_Num_ch(N, 2);
        auto patch = multinomial_patch(0.3, N);

        auto maybe_mps = lift_Macro_to_Micro(micro, patch);
        REQUIRE(maybe_mps.valid());
        auto mps = std::move(maybe_mps.value());

        auto recovered = project_Micro_to_Macro(micro, mps);
        auto const& recovered_mean = get<P_mean>(recovered())();
        auto const& recovered_cov = get<P_Cov>(recovered())();
        auto const& original_cov = get<P_Cov>(patch())();

        CHECK(std::abs(recovered_mean[0] - 0.3) < 1e-3);
        CHECK(std::abs(recovered_mean[1] - 0.7) < 1e-3);

        CHECK(std::abs(recovered_cov(0, 0) - original_cov(0, 0)) < 1e-5);
        CHECK(std::abs(recovered_cov(1, 1) - original_cov(1, 1)) < 1e-5);
    }

    SECTION("k=2, N=50, non-multinomial P_Cov (narrower than multinomial)") {
        std::size_t N = 50;
        auto micro = create_Micro_state_Num_ch(N, 2);

        Matrix<double> mean(1, 2);
        mean[0] = 0.5;
        mean[1] = 0.5;

        SymmetricMatrix<double> cov(2, 2);
        double inv_N = 1.0 / static_cast<double>(N);
        double base_var = 0.25 * inv_N;
        cov.set(0, 0, 0.5 * base_var);
        cov.set(1, 0, -0.5 * base_var);
        cov.set(1, 1, 0.5 * base_var);

        auto patch = make_patch_state(P_mean(std::move(mean)), P_Cov(std::move(cov)));

        auto maybe_mps = lift_Macro_to_Micro(micro, patch);
        REQUIRE(maybe_mps.valid());
        auto mps = std::move(maybe_mps.value());

        auto recovered = project_Micro_to_Macro(micro, mps);
        auto const& recovered_mean = get<P_mean>(recovered())();
        auto const& recovered_cov = get<P_Cov>(recovered())();
        auto const& original_cov = get<P_Cov>(patch())();

        CHECK(std::abs(recovered_mean[0] - 0.5) < 1e-2);

        CHECK(std::abs(recovered_cov(0, 0) - original_cov(0, 0)) < 5e-3);
    }
}

TEST_CASE("project_Micro_to_Macro on a delta distribution", "[micro_full][project]") {
    std::size_t N = 20;
    auto micro = create_Micro_state_Num_ch(N, 2);

    std::size_t size = micro().nrows();
    Matrix<double> p_vals(1, size);
    for (std::size_t i = 0; i < size; ++i) {
        p_vals[i] = 0.0;
    }
    std::size_t delta_row = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (micro()(i, std::size_t{0}) == 6) {
            delta_row = i;
            break;
        }
    }
    p_vals[delta_row] = 1.0;

    auto mps = make_micro_patch_state(Micro_P_mean(std::move(p_vals)));
    auto result = project_Micro_to_Macro(micro, mps);
    auto const& mean = get<P_mean>(result())();
    auto const& cov = get<P_Cov>(result())();

    CHECK(std::abs(mean[0] - 0.3) < 1e-12);
    CHECK(std::abs(mean[1] - 0.7) < 1e-12);

    CHECK(std::abs(cov(0, 0)) < 1e-12);
    CHECK(std::abs(cov(1, 1)) < 1e-12);
    CHECK(std::abs(cov(1, 0)) < 1e-12);
}
