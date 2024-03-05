
#include "../catch2/catch.hpp"
#include "../gsl_integrate.h"

TEST_CASE("Poisson_noise_normalization", "[fancy]")
{
    
    CHECK(std::get<0>(Poisson_noise_normalization_p(1,0)) == std::sqrt(2*std::numbers::pi));
    CHECK(std::get<0>(Poisson_noise_normalization_p(1,1)) == 5.272307);
    CHECK(std::get<0>(Poisson_noise_normalization_p(10,1)) == Poisson_noise_normalization(10,1));
    CHECK(std::get<0>(Poisson_noise_normalization_p(100,1)) == Poisson_noise_normalization(100,1));
    CHECK(std::get<0>(Poisson_noise_normalization_p(100,10000)) == Poisson_noise_normalization(100,10000));
    CHECK(Poisson_noise_expected_logL(1,0)==-0.5 * std::log(2 * std::numbers::pi *1) - 0.5);
    CHECK(Poisson_noise_expected_logL(1,1)==0);
    CHECK(Poisson_noise_expected_logL(100,100)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(1)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(10)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(100)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(1000)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(10000)==0);
    CHECK(Poisson_noise_expected_lik_logL_pr(100000)==0);
    
    
    CHECK(Poisson_noise_expected_logL(10,10)==Poisson_noise_expected_logL_2(10,10));
    
//   CHECK(std::get<0>(Poisson_noise_normalization_p(100,100000)) == Poisson_noise_normalization(100,100000));
}
