#ifndef GSL_INTEGRATE_H
#define GSL_INTEGRATE_H
#include <exponential_matrix.h>
#include <cmath>
#include <limits>
#include <tuple>

#include "derivative_fwd.h"
#include "gsl/gsl_integration.h"

namespace var{
template<class C_double>
requires (U<C_double,double>)
auto lik_Poisson_noise(C_double x, C_double variance, C_double Poisson_noise) {
    using std::exp;
    using std::abs;
    return exp(-1.0 / 2.0 * x * x / (variance + Poisson_noise * abs(x)));

}

template<class C_double>
requires (U<C_double,double>)
auto f_lik_Poisson_noise(C_double x, void* params) {
    auto variance = ((double*)params)[0];
    auto Poisson_noise = ((double*)params)[1];
    auto f = lik_Poisson_noise(x, variance, Poisson_noise);
    return f;
}

template<class C_double>
requires (U<C_double,double>)
auto f_lik_logL_Poisson_noise(C_double x, void* params) {
    auto variance = ((double*)params)[0];
    auto Poisson_noise = ((double*)params)[1];
    auto f = lik_Poisson_noise(x, variance, Poisson_noise);
    using std::log;
    if (f > 0)
        f = log(f) * f;
    return f;
}

template<class C_double>
requires (U<C_double,double>)
auto lik_Poisson_noise_f(C_double x, C_double Poisson_noise_ratio) {
    using std::exp;
    using std::abs;
    return exp(-1.0 / 2.0 * x * x / (1.0 + Poisson_noise_ratio * abs(x)));
}

template<class C_double>
requires (U<C_double,double>)
auto lik_logL_Poisson_noise_f(C_double x, C_double Poisson_noise_ratio) {
    using std::abs;
    using std::exp;
    auto r = -1.0 / 2.0 * x * x / (1.0 + Poisson_noise_ratio * abs(x));
    return r * exp(r);
}

template<class C_double>
requires (U<C_double,double>)
auto f_lik_Poisson_noise_f(C_double x, void* params) {
    auto Poisson_noise_ratio = ((double*)params)[0];
    auto f = lik_Poisson_noise_f(x, Poisson_noise_ratio);
    return f;
}

template<class C_double>
requires (U<C_double,double>)
auto f_lik_logL_Poisson_noise_f(C_double x, void* params) {
    auto Poisson_noise_ratio = ((double*)params)[0];
    auto f = lik_logL_Poisson_noise_f(x, Poisson_noise_ratio);
    return f;
}

template<class C_double>
requires (U<C_double,double>)
 auto Poisson_noise_normalization_p(C_double noise, C_double Poisson_noise) {
    auto eps = std::numeric_limits<double>::epsilon();
    //  auto epsabs=std::sqrt(eps);
    //  auto epsrel=std::sqrt(eps);
    using std::pow;
    auto epsabs = pow(eps, 0.25);
    auto epsrel = epsabs;
    double result;
    double abserr;
    auto limit = 30;
    auto workspace = gsl_integration_workspace_alloc(limit);
    gsl_function F;
    F.function = f_lik_Poisson_noise;  // Set integrand
    double params[] = {noise, Poisson_noise};
    F.params = params;
    gsl_integration_qagi(&F, epsabs, epsrel, limit, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    return std::tuple(result, abserr);
}



template<class C_double>
requires (U<C_double,double>)
auto Poisson_noise_expected_lik_logL_p(C_double noise, C_double Poisson_noise) {
    auto eps = std::numeric_limits<double>::epsilon();
    //  auto epsabs=std::sqrt(eps);
    //  auto epsrel=std::sqrt(eps);
    auto epsabs = std::pow(eps, 0.25);
    auto epsrel = epsabs;
    double result;
    double abserr;
    auto limit = 30;
    auto workspace = gsl_integration_workspace_alloc(limit);
    gsl_function F;
    F.function = f_lik_logL_Poisson_noise;  // Set integrand
    double params[] = {noise, Poisson_noise};
    F.params = params;
    gsl_integration_qagi(&F, epsabs, epsrel, limit, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    return std::tuple(result, abserr);
}

template<class C_double>
requires (U<C_double,double>)
 auto Poisson_noise_normalization_pr(C_double Poisson_noise_ratio) {
    if (Poisson_noise_ratio > 1e4)
        return 4 * Poisson_noise_ratio;
    else {
        auto eps = std::numeric_limits<double>::epsilon();
        //  auto epsabs=std::sqrt(eps);
        //  auto epsrel=std::sqrt(eps);
        auto epsabs = std::pow(eps, 0.25);
        auto epsrel = epsabs;
        double result;
        double abserr;
        auto limit = 30;
        auto workspace = gsl_integration_workspace_alloc(limit);
        gsl_function F;
        F.function = f_lik_Poisson_noise_f;  // Set integrand
        C_double params[] = {Poisson_noise_ratio};
        F.params = params;
        gsl_integration_qagi(&F, epsabs, epsrel, limit, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return primitive(Poisson_noise_ratio)/Poisson_noise_ratio* result;
    }
}

template<class C_double>
requires (U<C_double,double>)
auto Poisson_noise_expected_lik_logL_pr(C_double Poisson_noise_ratio) {
    if (Poisson_noise_ratio > 1e3)
        return 4 * Poisson_noise_ratio;
    else {
        auto eps = std::numeric_limits<double>::epsilon();
        //  auto epsabs=std::sqrt(eps);
        //  auto epsrel=std::sqrt(eps);
        auto epsabs = std::pow(eps, 0.25);
        auto epsrel = epsabs;
        double result;
        double abserr;
        auto limit = 30;
        auto workspace = gsl_integration_workspace_alloc(limit);
        gsl_function F;
        F.function = f_lik_logL_Poisson_noise_f;  // Set integrand
        double params[] = {primitive(Poisson_noise_ratio)};
        F.params = params;
        gsl_integration_qagi(&F, epsabs, epsrel, limit, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return result*primitive(Poisson_noise_ratio)/Poisson_noise_ratio;
    }
}


template<class C_double>
requires (U<C_double,double>)
 auto Poisson_noise_expected_lik_logL(C_double noise, C_double Poisson_noise) {
    using std::sqrt;
    auto s = sqrt(noise);

    auto Pn = Poisson_noise / s;
    return Poisson_noise_expected_lik_logL_pr(Pn) * s;
}

template<class C_double>
requires (U<C_double,double>)
auto Poisson_noise_normalization(C_double noise, C_double Poisson_noise) {
    using std::sqrt;
    auto s = sqrt(noise);
    auto Pn = Poisson_noise / var::primitive(s);
    return Poisson_noise_normalization_pr(Pn) * s;
}

template<class C_double>
requires (U<C_double,double>)
 auto Poisson_noise_expected_logL_2(C_double noise, C_double Poisson_noise) {
    auto z = Poisson_noise_normalization(noise, Poisson_noise);
    auto lik_logL = std::get<0>(Poisson_noise_expected_lik_logL_p(noise, Poisson_noise));
    return lik_logL / z - std::log(z);
}

template<class C_double>
requires (U<C_double,double>)
auto Poisson_noise_expected_logL(C_double noise, C_double Poisson_noise) {
    auto z = Poisson_noise_normalization(noise, Poisson_noise);
    auto lik_logL = Poisson_noise_expected_lik_logL(noise, Poisson_noise);
    using std::log;
    return lik_logL / z - log(z);
}
}
#endif  // GSL_INTEGRATE_H
