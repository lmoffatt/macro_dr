
#include "../catch2/catch.hpp"
#include "../qmodel.h"

TEST_CASE("all_probability_elements and all_covariance_elements", "[fancy]") {
    auto x = Matrix<double>(1, 5, {0.1, 0.1, 0.1, 0.4, 0.3});

    REQUIRE(macrodr::all_Probability_elements(x));

    x[1] = -0.1;
    REQUIRE(!macrodr::all_Probability_elements(x));
}
