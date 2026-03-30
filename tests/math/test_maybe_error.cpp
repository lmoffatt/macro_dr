#include <catch_amalgamated.hpp>

#include "maybe_error.h"

#include <memory>

TEST_CASE("apply_on_Maybe_error perfect-forwards plain rvalues", "[maybe_error]") {
    auto maybe = apply_on_Maybe_error(
        [](std::unique_ptr<int> p) { return *p + 1; },
        std::make_unique<int>(4));

    REQUIRE(maybe);
    CHECK(maybe.value() == 5);
}

TEST_CASE("apply_on_Maybe_error perfect-forwards Maybe_error rvalues", "[maybe_error]") {
    Maybe_error<std::unique_ptr<int>> arg(std::make_unique<int>(6));

    auto maybe = apply_on_Maybe_error(
        [](std::unique_ptr<int> p) { return *p + 2; },
        std::move(arg));

    REQUIRE(maybe);
    CHECK(maybe.value() == 8);
}
