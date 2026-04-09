#include "catch_amalgamated.hpp"

#include <limits>
#include <map>
#include <set>
#include <vector>

#include "distributions.h"
#include "moment_statistics.h"

using ::Matrix;
using ::SymPosDefMatrix;

namespace {

Matrix<double> make_matrix(std::initializer_list<double> values) {
    Matrix<double> out(2, 2, false);
    auto it = values.begin();
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            out(i, j) = *it;
            ++it;
        }
    }
    return out;
}

SymPosDefMatrix<double> make_diag_spd(std::initializer_list<double> diag_vals) {
    const std::size_t n = diag_vals.size();
    SymPosDefMatrix<double> out(n, n, false);
    std::size_t i = 0;
    for (auto v : diag_vals) {
        out.set(i, i, v);
        ++i;
    }
    return out;
}

}  // namespace

TEST_CASE("matrix get_mean_Probits drops empty samples", "[probit][matrix]") {
    std::vector<Matrix<double>> samples{
        Matrix<double>{},
        make_matrix({1.0, 2.0, 3.0, 4.0}),
        make_matrix({5.0, 6.0, 7.0, 8.0})};

    auto [mean_matrix, probits,n] = get_mean_Probits(samples, std::set<double>{0.5}, std::identity{});

    REQUIRE(mean_matrix.nrows() == 2);
    REQUIRE(mean_matrix.ncols() == 2);
    CHECK(mean_matrix(std::size_t{0}, std::size_t{0}) == Catch::Approx(3.0));
    CHECK(mean_matrix(std::size_t{0}, std::size_t{1}) == Catch::Approx(4.0));
    CHECK(mean_matrix(std::size_t{1}, std::size_t{0}) == Catch::Approx(5.0));
    CHECK(mean_matrix(std::size_t{1}, std::size_t{1}) == Catch::Approx(6.0));

    REQUIRE(probits.size() == 1);
    const auto& q50 = probits.at(0.5);
    CHECK(q50(std::size_t{0}, std::size_t{0}) == Catch::Approx(3.0));
    CHECK(q50(std::size_t{0}, std::size_t{1}) == Catch::Approx(4.0));
    CHECK(q50(std::size_t{1}, std::size_t{0}) == Catch::Approx(5.0));
    CHECK(q50(std::size_t{1}, std::size_t{1}) == Catch::Approx(6.0));
    CHECK(n==2);
}

TEST_CASE("matrix get_mean_Probits uses first non-empty shape", "[probit][matrix]") {
    std::vector<Matrix<double>> samples{
        Matrix<double>{},
        make_matrix({2.0, 4.0, 6.0, 8.0})};

    auto [mean_matrix, probits, n] = 
        get_mean_Probits(samples, std::set<double>{0.5}, std::identity{});

    REQUIRE(mean_matrix.nrows() == 2);
    REQUIRE(mean_matrix.ncols() == 2);
    CHECK(mean_matrix(std::size_t{0}, std::size_t{0}) == Catch::Approx(2.0));
    CHECK(mean_matrix(std::size_t{0}, std::size_t{1}) == Catch::Approx(4.0));
    CHECK(mean_matrix(std::size_t{1}, std::size_t{0}) == Catch::Approx(6.0));
    CHECK(mean_matrix(std::size_t{1}, std::size_t{1}) == Catch::Approx(8.0));
    CHECK(n==1);

    REQUIRE(probits.size() == 1);
    const auto& q50 = probits.at(0.5);
    CHECK(q50(std::size_t{0}, std::size_t{0}) == Catch::Approx(2.0));
    CHECK(q50(std::size_t{0}, std::size_t{1}) == Catch::Approx(4.0));
    CHECK(q50(std::size_t{1}, std::size_t{0}) == Catch::Approx(6.0));
    CHECK(q50(std::size_t{1}, std::size_t{1}) == Catch::Approx(8.0));
}

TEST_CASE("matrix get_mean_Probits returns empty matrices when all samples are empty",
          "[probit][matrix]") {
    std::vector<Matrix<double>> samples{Matrix<double>{}, Matrix<double>{}};

    auto [mean_matrix, probits, n] = get_mean_Probits(samples, std::set<double>{0.25, 0.75}, std::identity{});

    CHECK(mean_matrix.size() == 0);
    CHECK(mean_matrix.nrows() == 0);
    CHECK(mean_matrix.ncols() == 0);
    CHECK(n==0);
    REQUIRE(probits.size() == 2);
    CHECK(probits.at(0.25).size() == 0);
    CHECK(probits.at(0.75).size() == 0);
}

TEST_CASE("Probit_statistics for distortion matrices inherits empty-sample filtering",
          "[probit][matrix]") {
    using Idm = Information_Distortion_Matrix;

    std::vector<Idm> samples{
        Idm(SymPosDefMatrix<double>{}),
        Idm(make_diag_spd({2.0, 6.0})),
        Idm(make_diag_spd({4.0, 10.0}))};

    auto stats = Probit_statistics<Idm>(samples, [](const auto& x) { return x(); },
                                        std::set<double>{0.5});

    const auto& mean_matrix = get<mean<Idm>>(stats())();
    CHECK(mean_matrix(std::size_t{0}, std::size_t{0}) == Catch::Approx(3.0));
    CHECK(mean_matrix(std::size_t{1}, std::size_t{1}) == Catch::Approx(8.0));

    const auto& q50 = get<Probits<Idm>>(stats())().at(0.5);
    CHECK(q50(std::size_t{0}, std::size_t{0}) == Catch::Approx(3.0));
    CHECK(q50(std::size_t{1}, std::size_t{1}) == Catch::Approx(8.0));
}

TEST_CASE("vector-space get_mean_Probits falls back to component means when a quantile is missing",
          "[probit][vector-space]") {
    using Sample = var::Vector_Space<logL, elogL>;
    const auto nan = std::numeric_limits<double>::quiet_NaN();

    std::vector<Sample> samples{Sample{logL(-1.0), elogL(nan)},
                                Sample{logL(-3.0), elogL(nan)}};

    auto [mean_value, probits, n] =
        get_mean_Probits(samples, std::set<double>{0.5}, std::identity{});

    CHECK(get<logL>(mean_value)() == Catch::Approx(-2.0));
    CHECK(std::isnan(get<elogL>(mean_value)()));
    CHECK(n == 0);

    REQUIRE(probits.size() == 1);
    const auto& q50 = probits.at(0.5);
    CHECK(get<logL>(q50)() == Catch::Approx(-2.0));
    CHECK(std::isnan(get<elogL>(q50)()));
}
