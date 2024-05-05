#ifndef GSL_INTEGRATE_H
#define GSL_INTEGRATE_H
#include <cmath>
#include "gsl/gsl_integration.h"
#include <limits>
#include <tuple>
#include "derivative_operator.h"

inline auto lik_Poisson_noise(double x, double variance, double Poisson_noise)
{
    return std::exp(-1.0/2.0*x*x/(variance+Poisson_noise*std::abs(x)));
   
}


inline double f_lik_Poisson_noise (double x, void * params) {
    double variance = ((double *) params)[0];
    double Poisson_noise = ((double *) params)[1];
    double f = lik_Poisson_noise(x,variance,Poisson_noise);
    return f;
}



inline double f_lik_logL_Poisson_noise (double x, void * params) {
    double variance = ((double *) params)[0];
    double Poisson_noise = ((double *) params)[1];
    double f = lik_Poisson_noise(x,variance,Poisson_noise);
    if (f>0)
        f= std::log(f)*f;
    return f;
}




inline auto lik_Poisson_noise_f(double x,  double Poisson_noise_ratio)
{
    return std::exp(-1.0/2.0*x*x/(1.0+Poisson_noise_ratio*std::abs(x)));
    
}

inline auto lik_logL_Poisson_noise_f(double x,  double Poisson_noise_ratio)
{
    auto r=-1.0/2.0*x*x/(1.0+Poisson_noise_ratio*std::abs(x));
    return r*std::exp(r);
    
}




inline double f_lik_Poisson_noise_f (double x, void * params) {
    double Poisson_noise_ratio = ((double *) params)[0];
    double f = lik_Poisson_noise_f(x,Poisson_noise_ratio);
    return f;
}

inline double f_lik_logL_Poisson_noise_f (double x, void * params) {
    double Poisson_noise_ratio = ((double *) params)[0];
    double f = lik_logL_Poisson_noise_f(x,Poisson_noise_ratio);
    return f;
}


inline auto Poisson_noise_normalization_p(double noise,double Poisson_noise){
    
    
    auto eps=std::numeric_limits<double>::epsilon();
  //  auto epsabs=std::sqrt(eps);
  //  auto epsrel=std::sqrt(eps);
    auto epsabs=std::pow(eps,0.25);
    auto epsrel=epsabs;
    double result;
    double abserr;
    auto limit=30;
    auto workspace=gsl_integration_workspace_alloc(limit);
    gsl_function F;
    F.function = f_lik_Poisson_noise;               // Set integrand
    double params[] = { noise, Poisson_noise };
    F.params = params;    
    gsl_integration_qagi(&F, epsabs, epsrel,  limit, workspace, &result, &abserr);
    gsl_integration_workspace_free (workspace);
    
    return std::tuple(result,abserr);
        
}

inline auto Poisson_noise_expected_lik_logL_p(double noise,double Poisson_noise){
    
    
    auto eps=std::numeric_limits<double>::epsilon();
    //  auto epsabs=std::sqrt(eps);
    //  auto epsrel=std::sqrt(eps);
    auto epsabs=std::pow(eps,0.25);
    auto epsrel=epsabs;
    double result;
    double abserr;
    auto limit=30;
    auto workspace=gsl_integration_workspace_alloc(limit);
    gsl_function F;
    F.function = f_lik_logL_Poisson_noise;               // Set integrand
    double params[] = { noise, Poisson_noise };
    F.params = params;    
    gsl_integration_qagi(&F, epsabs, epsrel,  limit, workspace, &result, &abserr);
    gsl_integration_workspace_free (workspace);
    
    return std::tuple(result,abserr);
    
}





inline auto Poisson_noise_normalization_pr(double Poisson_noise_ratio){
    
    if (Poisson_noise_ratio>1e4)
        return 4*Poisson_noise_ratio;
    else{
    
    auto eps=std::numeric_limits<double>::epsilon();
    //  auto epsabs=std::sqrt(eps);
    //  auto epsrel=std::sqrt(eps);
    auto epsabs=std::pow(eps,0.25);
    auto epsrel=epsabs;
    double result;
    double abserr;
    auto limit=30;
    auto workspace=gsl_integration_workspace_alloc(limit);
    gsl_function F;
    F.function = f_lik_Poisson_noise_f;               // Set integrand
    double params[] = {Poisson_noise_ratio };
    F.params = params;    
    gsl_integration_qagi(&F, epsabs, epsrel,  limit, workspace, &result, &abserr);
    gsl_integration_workspace_free (workspace);
    return result;
    }
}

inline auto Poisson_noise_expected_lik_logL_pr(double Poisson_noise_ratio){
    if (Poisson_noise_ratio>1e3)
        return 4*Poisson_noise_ratio;
    else{
    
        
        auto eps=std::numeric_limits<double>::epsilon();
        //  auto epsabs=std::sqrt(eps);
        //  auto epsrel=std::sqrt(eps);
        auto epsabs=std::pow(eps,0.25);
        auto epsrel=epsabs;
        double result;
        double abserr;
        auto limit=30;
        auto workspace=gsl_integration_workspace_alloc(limit);
        gsl_function F;
        F.function = f_lik_logL_Poisson_noise_f;               // Set integrand
        double params[] = {Poisson_noise_ratio };
        F.params = params;    
        gsl_integration_qagi(&F, epsabs, epsrel,  limit, workspace, &result, &abserr);
        gsl_integration_workspace_free (workspace);
        return result;
    }
}
inline auto Poisson_noise_expected_lik_logL(double noise,double Poisson_noise){
    auto s=std::sqrt(noise);
    auto Pn=Poisson_noise/s;
    return Poisson_noise_expected_lik_logL_pr(Pn)*s;
}


inline double Poisson_noise_normalization(double noise,double Poisson_noise){
    using std::sqrt;
    auto s=sqrt(noise);
    auto Pn=Poisson_noise/var::primitive(s);
    return Poisson_noise_normalization_pr(Pn)*s;
}

inline auto Poisson_noise_expected_logL_2(double noise,double Poisson_noise){
        auto z=Poisson_noise_normalization(noise,Poisson_noise);
        auto lik_logL=std::get<0>(Poisson_noise_expected_lik_logL_p(noise,Poisson_noise));
        return lik_logL/z-std::log(z);
}

inline auto Poisson_noise_expected_logL(double noise,double Poisson_noise){
    auto z=Poisson_noise_normalization(noise,Poisson_noise);
    auto lik_logL=Poisson_noise_expected_lik_logL(noise,Poisson_noise);
    return lik_logL/z-std::log(z);
}


#endif // GSL_INTEGRATE_H
