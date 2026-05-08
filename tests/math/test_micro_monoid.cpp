#include "catch_amalgamated.hpp"

#include <cstddef>
#include <set>
#include <vector>

#include "micro_monoid.h"

using namespace macrodr;

namespace {

// Walk every idx in [0, count): check round-trip identity, length, sum, uniqueness,
// and strict lex-increasing order. count = num_full_states_of(N, k).
void exhaustive_round_trip(std::size_t N, std::size_t k) {
    std::size_t const count = num_full_states_of(N, k);
    std::set<std::vector<std::size_t>> seen;
    std::vector<std::size_t> prev;
    for (std::size_t idx = 0; idx < count; ++idx) {
        auto n = index_to_microstate(idx, N, k);

        REQUIRE(n.size() == k);

        std::size_t s = 0;
        for (auto x : n) s += x;
        CHECK(s == N);

        CHECK(microstate_to_index(n, N, k) == idx);

        CHECK(seen.insert(n).second);

        if (idx > 0) {
            CHECK(prev < n);  // std::vector lex comparison
        }
        prev = std::move(n);
    }
    CHECK(seen.size() == count);
}

}  // namespace

TEST_CASE("micro_monoid: num_full_states_of = C(N+k-1, k-1)", "[micro_monoid]") {
    CHECK(num_full_states_of(10, 2) == 11);
    CHECK(num_full_states_of(4, 3) == 15);
    CHECK(num_full_states_of(100, 5) == 4598126);
    CHECK(num_full_states_of(0, 3) == 1);  // single all-zero state
    CHECK(num_full_states_of(7, 1) == 1);  // single bin: must hold all N
    CHECK(num_full_states_of(0, 0) == 0);  // convention: empty alphabet, no states
}

TEST_CASE("micro_monoid: index <-> microstate round-trip", "[micro_monoid]") {
    SECTION("k=1, N=0..5") {
        for (std::size_t N = 0; N <= 5; ++N) {
            exhaustive_round_trip(N, 1);
        }
    }
    SECTION("k=2") {
        for (std::size_t N : {std::size_t{0}, std::size_t{1}, std::size_t{10}, std::size_t{50}}) {
            exhaustive_round_trip(N, 2);
        }
    }
    SECTION("k=3, N=0..8") {
        for (std::size_t N = 0; N <= 8; ++N) {
            exhaustive_round_trip(N, 3);
        }
    }
    SECTION("k=4, N=0..6") {
        for (std::size_t N = 0; N <= 6; ++N) {
            exhaustive_round_trip(N, 4);
        }
    }
    SECTION("k=5, N=4 (70 states)") {
        exhaustive_round_trip(4, 5);
    }
    SECTION("k=6, N=5 (252 states)") {
        exhaustive_round_trip(5, 6);
    }
}

TEST_CASE("micro_monoid: edge cases", "[micro_monoid]") {
    SECTION("k=0 yields empty vector") {
        CHECK(index_to_microstate(0, 0, 0).empty());
        CHECK(index_to_microstate(0, 7, 0).empty());
        CHECK(microstate_to_index({}, 0, 0) == 0);
    }
    SECTION("k=1 yields (N) at idx=0") {
        for (std::size_t N : {std::size_t{0}, std::size_t{1}, std::size_t{42}}) {
            auto v = index_to_microstate(0, N, 1);
            REQUIRE(v.size() == 1);
            CHECK(v[0] == N);
            CHECK(microstate_to_index(v, N, 1) == 0);
        }
    }
    SECTION("N=0 yields the all-zero tuple of length k at idx=0") {
        for (std::size_t k : {std::size_t{1}, std::size_t{3}, std::size_t{7}}) {
            auto v = index_to_microstate(0, 0, k);
            REQUIRE(v.size() == k);
            for (auto x : v) CHECK(x == 0);
            CHECK(microstate_to_index(v, 0, k) == 0);
        }
    }
}

TEST_CASE("micro_monoid: spot-check vs. hand-computed enumeration (N=4, k=3)",
          "[micro_monoid]") {
    // Lex order with leftmost coordinate most significant:
    //  0:(0,0,4) 1:(0,1,3) 2:(0,2,2) 3:(0,3,1) 4:(0,4,0)
    //  5:(1,0,3) 6:(1,1,2) 7:(1,2,1) 8:(1,3,0)
    //  9:(2,0,2) 10:(2,1,1) 11:(2,2,0)
    // 12:(3,0,1) 13:(3,1,0)
    // 14:(4,0,0)
    using V = std::vector<std::size_t>;
    std::vector<V> expected = {
        {0, 0, 4}, {0, 1, 3}, {0, 2, 2}, {0, 3, 1}, {0, 4, 0},
        {1, 0, 3}, {1, 1, 2}, {1, 2, 1}, {1, 3, 0},
        {2, 0, 2}, {2, 1, 1}, {2, 2, 0},
        {3, 0, 1}, {3, 1, 0},
        {4, 0, 0},
    };
    REQUIRE(num_full_states_of(4, 3) == expected.size());
    for (std::size_t idx = 0; idx < expected.size(); ++idx) {
        CHECK(index_to_microstate(idx, 4, 3) == expected[idx]);
        CHECK(microstate_to_index(expected[idx], 4, 3) == idx);
    }
}
