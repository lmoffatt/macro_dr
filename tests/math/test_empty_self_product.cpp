#include "catch_amalgamated.hpp"

#include "distributions.h"
#include "lapack_headers.h"
#include "matrix.h"
#include "moment_statistics.h"

TEST_CASE("vectorized self product accepts empty matrices", "[lapack][empty]") {
    Matrix<double> empty;

    auto cov = lapack::Lapack_Product_Self_Transpose_vectorized(empty);

    REQUIRE(cov.size() == 0);
    REQUIRE(cov.nrows() == 0);
    REQUIRE(cov.ncols() == 0);
}

TEST_CASE("moment statistics stays default-shaped for empty dlogL inputs", "[moment_statistics][empty]") {
    Moment_statistics<dlogL, true> from_empty_ids(std::vector<dlogL>{});

    CHECK(get<count<dlogL>>(from_empty_ids())() == 0);
    CHECK(get<mean<dlogL>>(from_empty_ids())().size() == 0);
    CHECK(get<covariance<dlogL>>(from_empty_ids())().size() == 0);

    std::vector<Matrix<double>> empty_values;
    std::vector<std::size_t> empty_indices;
    Moment_statistics<dlogL, true> from_empty_indices(
        empty_values, empty_indices, [](const Matrix<double>& value) { return value; });

    CHECK(get<count<dlogL>>(from_empty_indices())() == 0);
    CHECK(get<mean<dlogL>>(from_empty_indices())().size() == 0);
    CHECK(get<covariance<dlogL>>(from_empty_indices())().size() == 0);
}
