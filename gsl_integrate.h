#ifndef GSL_INTEGRATE_H
#define GSL_INTEGRATE_H
#include <cmath>
#include <gsl/gsl_integration.h>
#include <limits>
#include <tuple>


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

inline auto lik_Poisson_noise_f(double x,  double Poisson_noise_ratio)
{
    return std::exp(-1.0/2.0*x*x/(1.0+Poisson_noise_ratio*std::abs(x)));
    
}


inline double f_lik_Poisson_noise_f (double x, void * params) {
    double Poisson_noise_ratio = ((double *) params)[0];
    double f = lik_Poisson_noise_f(x,Poisson_noise_ratio);
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


inline auto Poisson_noise_normalization(double noise,double Poisson_noise){
    auto s=std::sqrt(noise);
    auto Pn=Poisson_noise/s;
    return Poisson_noise_normalization_pr(Pn)*s;
}




#endif // GSL_INTEGRATE_H
