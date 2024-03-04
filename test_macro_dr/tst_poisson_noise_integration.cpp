
#include "../catch2/catch.hpp"
#include "../gsl_integrate.h"

TEST_CASE("Poisson_noise_normalization", "[fancy]")
{
    
    CHECK(std::get<0>(Poisson_noise_normalization_p(1,0)) == std::sqrt(2*std::numbers::pi));
    CHECK(std::get<0>(Poisson_noise_normalization_p(1,1)) == 5.272307);
    CHECK(std::get<0>(Poisson_noise_normalization_p(10,1)) == Poisson_noise_normalization(10,1));
    CHECK(std::get<0>(Poisson_noise_normalization_p(100,1)) == Poisson_noise_normalization(100,1));
    CHECK(std::get<0>(Poisson_noise_normalization_p(100,10000)) == Poisson_noise_normalization(100,10000));
//   CHECK(std::get<0>(Poisson_noise_normalization_p(100,100000)) == Poisson_noise_normalization(100,100000));
}
