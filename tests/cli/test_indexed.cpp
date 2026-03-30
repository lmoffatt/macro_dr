#include <catch_amalgamated.hpp>

#include "indexed.h"

#include <string>
#include <vector>

namespace {

var::Axis axis(std::size_t id, std::size_t n) {
    return var::Axis{var::AxisId{id}, var::AxisSize{n}};
}

var::Coordinate coordinate(std::initializer_list<var::Axis> axes,
                           std::initializer_list<var::AxisIndex> indexes) {
    return var::Coordinate(std::vector<var::Axis>(axes), std::vector<var::AxisIndex>(indexes));
}

}  // namespace

TEST_CASE("Indexed apply_Index aligns shared axes and expands independent axes") {
    var::IndexSpace model_space{{axis(1, 2)}};
    var::IndexSpace algo_space{{axis(2, 3)}};

    var::Indexed<std::string> models(model_space, {"CCO", "M2"});
    var::Indexed<int> algos(algo_space, {10, 20, 30});

    auto maybe = var::apply_Index(
        [](const std::string& model, int algo) { return model + ":" + std::to_string(algo); },
        models, algos);

    REQUIRE(maybe);
    const auto& out = maybe.value();
    REQUIRE(out.index_space().m_axes.size() == 2);
    REQUIRE(out.size() == 6);
    CHECK(out[0] == "CCO:10");
    CHECK(out[1] == "M2:10");
    CHECK(out[2] == "CCO:20");
    CHECK(out[5] == "M2:30");

    auto maybe_first = out.at(coordinate({axis(1, 2), axis(2, 3)},
                                         {var::AxisIndex{0}, var::AxisIndex{0}}));
    REQUIRE(maybe_first);
    CHECK(maybe_first.value().get() == "CCO:10");
}

TEST_CASE("Indexed apply_Index rejects conflicting shared-axis sizes") {
    var::Indexed<int> left(var::IndexSpace{{axis(7, 2)}}, {1, 2});
    var::Indexed<int> right(var::IndexSpace{{axis(7, 3)}}, {10, 20, 30});

    auto maybe = var::apply_Index([](int x, int y) { return x + y; }, left, right);

    REQUIRE_FALSE(maybe);
    CHECK(maybe.error()().find("conflicting sizes") != std::string::npos);
}

TEST_CASE("Indexed coordinate validation checks bounds and axis-size agreement") {
    var::Indexed<std::string> models(var::IndexSpace{{axis(1, 2)}}, {"CCO", "M2"});

    auto bad_size =
        models.at(coordinate({axis(1, 3)}, {var::AxisIndex{1}}));
    REQUIRE_FALSE(bad_size);
    CHECK(bad_size.error()().find("size mismatch") != std::string::npos);

    auto bad_index =
        models.at(coordinate({axis(1, 2)}, {var::AxisIndex{3}}));
    REQUIRE_FALSE(bad_index);
    CHECK(bad_index.error()().find("out of range") != std::string::npos);

    auto duplicate_axis = coordinate({axis(1, 2), axis(1, 2)},
                                     {var::AxisIndex{0}, var::AxisIndex{1}});
    auto bad_duplicate = models.at(duplicate_axis);
    REQUIRE_FALSE(bad_duplicate);
    CHECK(bad_duplicate.error()().find("duplicate axis") != std::string::npos);
}

TEST_CASE("Indexed apply_Index returns an empty indexed payload for zero-sized axes") {
    var::Indexed<int> empty(var::IndexSpace{{axis(3, 0)}}, {});

    auto maybe = var::apply_Index([](int x) { return x + 1; }, empty);

    REQUIRE(maybe);
    CHECK(maybe.value().size() == 0);
    CHECK(maybe.value().index_space().size() == 0);
}

TEST_CASE("Indexed rejects payloads whose size does not match the index-space cardinality") {
    var::Indexed<int> bad(var::IndexSpace{{axis(9, 2)}}, {1, 2, 3});

    auto maybe = bad.at(coordinate({axis(9, 2)}, {var::AxisIndex{0}}));

    REQUIRE_FALSE(maybe);
    CHECK(maybe.error()().find("size mismatch") != std::string::npos);

    auto lifted = var::apply_Index([](int x) { return x + 1; }, bad);
    REQUIRE_FALSE(lifted);
    CHECK(lifted.error()().find("size mismatch") != std::string::npos);
}

TEST_CASE("Indexed rejects invalid standalone index spaces before iteration") {
    var::IndexSpace invalid{{axis(11, 2), axis(11, 2)}};

    auto maybe = var::apply_Index([](int x) { return x + 1; },
                                  var::Indexed<int>(invalid, {1, 2, 3, 4}));

    REQUIRE_FALSE(maybe);
    CHECK(maybe.error()().find("duplicate axis") != std::string::npos);
}

TEST_CASE("Indexed apply_Index handles const references, mutable references, and wrappers") {
    var::Indexed<int> xs(var::IndexSpace{{axis(13, 2)}}, {3, 5});

    auto by_const_ref = var::apply_Index([](const int& x) { return x + 1; }, xs);
    REQUIRE(by_const_ref);
    CHECK(by_const_ref.value()[0] == 4);
    CHECK(by_const_ref.value()[1] == 6);
    CHECK(xs[0] == 3);
    CHECK(xs[1] == 5);

    auto by_mutable_ref = var::apply_Index([](int& x) {
        x += 10;
        return x;
    }, xs);
    REQUIRE(by_mutable_ref);
    CHECK(by_mutable_ref.value()[0] == 13);
    CHECK(by_mutable_ref.value()[1] == 15);
    CHECK(xs[0] == 13);
    CHECK(xs[1] == 15);

    auto by_wrapper = var::apply_Index([](std::reference_wrapper<int> x) {
        x.get() *= 2;
        return x.get();
    }, xs);
    REQUIRE(by_wrapper);
    CHECK(by_wrapper.value()[0] == 26);
    CHECK(by_wrapper.value()[1] == 30);
    CHECK(xs[0] == 26);
    CHECK(xs[1] == 30);
}
