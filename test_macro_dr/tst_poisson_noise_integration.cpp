
#include "../catch2/catch.hpp"
#include "../gsl_integrate.h"

TEST_CASE("Poisson_noise_normalization", "[fancy]")
{
    
    CHECK_THAT(std::get<0>(Poisson_noise_normalization_p(1,0)), Catch::Matchers::WithinRel(sqrt(2*std::numbers::pi),1e-5));
    CHECK_THAT(std::get<0>(Poisson_noise_normalization_p(1,1)), Catch::Matchers::WithinRel( 5.272307,1e-5));
    CHECK_THAT(std::get<0>(Poisson_noise_normalization_p(10,1)), Catch::Matchers::WithinRel( Poisson_noise_normalization(10,1),1e-5));
    CHECK_THAT(std::get<0>(Poisson_noise_normalization_p(100,1)), Catch::Matchers::WithinRel( Poisson_noise_normalization(100,1),1e-5));
    CHECK_THAT(std::get<0>(Poisson_noise_normalization_p(100,10000)), Catch::Matchers::WithinRel( Poisson_noise_normalization(100,10000),1e-5));
    CHECK_THAT(Poisson_noise_expected_logL(1,0),Catch::Matchers::WithinRel(-0.5 * std::log(2 * std::numbers::pi *1) - 0.5,1e-5));
    CHECK_THAT(Poisson_noise_expected_logL(1,1),Catch::Matchers::WithinRel(0,1e-5));
    
    
    CHECK_THAT(Poisson_noise_expected_logL(10,10),Catch::Matchers::WithinRel(Poisson_noise_expected_logL_2(10,10),1e-5));
    
//   CHECK(std::get<0>(Poisson_noise_normalization_p(100,100000)) == Poisson_noise_normalization(100,100000));
}
