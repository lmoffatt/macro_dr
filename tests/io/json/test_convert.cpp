#include <catch_amalgamated.hpp>

#include <macrodr/io/json/convert.h>

#include <cstdint>
#include <macrodr/io/json/minijson.h>
#include <macrodr/dsl/type_name.h>
#include <matrix.h>
#include <variables.h>

namespace conv = macrodr::io::json::conv;
using macrodr::io::json::Json;
using Catch::Approx;


struct Foo : public var::Var<Foo, double> {
    using var::Var<Foo, double>::Var;
    friend std::string className(Foo) { return "Foo"; }
};

struct Bar : public var::Var<Bar, Matrix<double>> {
    using var::Var<Bar, Matrix<double>>::Var;
    friend std::string className(Bar) { return "Bar"; }
};

using TestSpace = var::Vector_Space<Foo, Bar>;

static Matrix<double> make_matrix(std::size_t rows, std::size_t cols) {
    Matrix<double> m(rows, cols, false);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            m(i, j) = static_cast<double>(i * cols + j + 1);
        }
    }
    return m;
}

TEST_CASE("Matrix round-trip") {
    auto m = make_matrix(2, 3);
    Json j = conv::to_json(m);
    Matrix<double> decoded;
    REQUIRE(conv::from_json(j, decoded));
    REQUIRE(decoded.nrows() == 2);
    REQUIRE(decoded.ncols() == 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t jcol = 0; jcol < 3; ++jcol)
            CHECK(decoded(i, jcol) == Approx(m(i, jcol)));
}

TEST_CASE("Vector and map round-trip") {
    std::vector<double> vec{1.0, 2.5, 3.75};
    Json jvec = conv::to_json(vec);
    std::vector<double> vec_out;
    REQUIRE(conv::from_json(jvec, vec_out));
    REQUIRE(vec_out == vec);

    std::map<std::string, int64_t> str_map{{"one", 1}, {"two", 2}};
    Json jmap = conv::to_json(str_map);
    std::map<std::string, int64_t> str_map_out;
    REQUIRE(conv::from_json(jmap, str_map_out));
    REQUIRE(str_map_out == str_map);

    std::map<std::pair<int64_t, int64_t>, double> pair_map{
        {{0, 1}, 1.5}, {{2, 3}, 4.0}};
    Json jpair_map = conv::to_json(pair_map);
    std::map<std::pair<int64_t, int64_t>, double> pair_map_out;
    REQUIRE(conv::from_json(jpair_map, pair_map_out));
    REQUIRE(pair_map_out == pair_map);
}

TEST_CASE("Tuple round-trip") {
    std::tuple<int64_t, double, std::string> tup{42, 2.5, "hello"};
    Json jtup = conv::to_json(tup);
    decltype(tup) tup_out;
    REQUIRE(conv::from_json(jtup, tup_out));
    REQUIRE(tup_out == tup);
}

TEST_CASE("Var round-trip") {
    Foo foo(3.14);
    Json jfoo = conv::to_json(foo);
    Foo foo_out;
    REQUIRE(conv::from_json(jfoo, foo_out));
    CHECK(foo_out() == Approx(foo()));

    Bar bar(make_matrix(2, 2));
    Json jbar = conv::to_json(bar);
    Bar bar_out;
    REQUIRE(conv::from_json(jbar, bar_out));
    auto& mat = bar();
    auto& decoded = bar_out();
    REQUIRE(decoded.nrows() == mat.nrows());
    REQUIRE(decoded.ncols() == mat.ncols());
    for (std::size_t i = 0; i < mat.nrows(); ++i)
        for (std::size_t j = 0; j < mat.ncols(); ++j)
            CHECK(decoded(i, j) == Approx(mat(i, j)));
}

TEST_CASE("Vector_Space round-trip") {
    TestSpace space;
    static_cast<Foo&>(space)() = 5.5;
    static_cast<Bar&>(space)() = make_matrix(2, 1);

    Json jspace = conv::to_json(space);
    TestSpace decoded;
    REQUIRE(conv::from_json(jspace, decoded));

    CHECK(static_cast<Foo&>(decoded)() == Approx(static_cast<Foo&>(space)()));
    auto& mat = static_cast<Bar&>(space)();
    auto& mat_decoded = static_cast<Bar&>(decoded)();
    REQUIRE(mat_decoded.nrows() == mat.nrows());
    REQUIRE(mat_decoded.ncols() == mat.ncols());
    for (std::size_t i = 0; i < mat.nrows(); ++i)
        for (std::size_t j = 0; j < mat.ncols(); ++j)
            CHECK(mat_decoded(i, j) == Approx(mat(i, j)));
}

