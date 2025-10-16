#include "catch_amalgamated.hpp"
#include "derivative_test.h"
#include "maybe_error.h"

using var::Derivative;
using var::test_taylor_convergence;

namespace {

struct LinearScalarFunctor {
    double slope{};

    Maybe_error<double> operator()(double x, double /*unused*/) const {
        return Maybe_error<double>(slope * x);
    }

    Maybe_error<Derivative<double, double>> operator()(const Derivative<double, double>& x,
                                                       const Derivative<double, double>&) const {
        return Maybe_error<Derivative<double, double>>(
            Derivative<double, double>(slope * x.primitive(), slope * x.derivative()(), x.dx()));
    }
};

}  // namespace

// TEST_CASE("test_derivative_clarke succeeds on linear scalar map", "[taylor][identity]") {
//     const double slope = 2.0;
//     double parameter_value = 1.2345;

//     LinearScalarFunctor functor{.slope = slope};
//     Derivative<double, double> parameter(parameter_value, 1.0, parameter_value);

//     auto result = var::test_derivative_clarke(functor, 1e-6, parameter, parameter);

//     REQUIRE(result);
// }
