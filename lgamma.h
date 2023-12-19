#ifndef LGAMMA_H
#define LGAMMA_H
#include <cmath>
#include <numbers>
namespace var {

inline double lgamma(double x)
{
    //Stirling approximation
    return x*std::log(x)-x+1.0/2*std::log(2* std::numbers::pi *x)+1.0/12.0/x-1.0/360.0/x/x/x;
}
} // namespace var

#endif // LGAMMA_H
