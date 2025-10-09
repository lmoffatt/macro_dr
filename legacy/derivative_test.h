#ifndef DERIVATIVE_TEST_H
#define DERIVATIVE_TEST_H
#include "derivative_operator.h"
#include "maybe_error.h"
#include "parameters_derivative.h"

namespace var {

template <class F, class Parameters, class... Xs>
Maybe_error<bool> test_Derivative(F f, const Parameters x, double dx, double eps, Xs const&... xs) {
    auto MaybedY = f(xs...);
    auto MaybeY = f(primitive(xs)...);

    if (!(is_valid(MaybedY))) {
        return get_error(MaybedY);
    }
    if (!(is_valid(MaybeY))) {
        return get_error(MaybeY);
    }
    auto dY = std::move(MaybedY.value());
    auto Y = std::move(MaybeY.value());

    auto same_primitive = test_equality(primitive(dY), Y, eps);
    if (!same_primitive) {
        std::stringstream ss;
        ss << "\n-----different primitive parts!!! ---- \n";
        ss << "\n Test_Derivative on function ";
        ss << type_name<F>();
        ss << " with parameters " << type_name<Parameters>() << " equal to: ";
        ss << "eps =" << eps << "\n";
        ss << same_primitive.error()();
        ss << "\n-----end of primitive parts--------\n";
        return error_message(ss.str());
    }

    auto T0 = Taylor_first(dY, x, dx);
    auto MaybedY1 = std::invoke(f, Taylor_first(xs, x, dx)...);
    if (!is_valid(MaybedY1)) {
        return get_error(MaybedY1);
    }
    auto dY1 = std::move(MaybedY1.value());

    auto out = test_equality(T0, dY1, eps);

    //    auto dout = test_equality(get_value(T0) - primitive(get_value(dY)),
    //                              get_value(T1) - primitive(get_value(dY)), eps);
    if (!out) {
        std::stringstream ss;
        ss << "\n Test_Derivative on function ";
        ss << type_name<F>();
        ss << " with delta parameters " << type_name<Parameters>() << " equal to: ";
        ss << "\n x dx\n" << x() << "\ndx=" << dx << "\n";
        ss << "eps =" << eps << "\n";

        ss << "\n-----error---------------\n";
        ss << out.error()();

        ss << "\n--------------------\n";

        ss << "\n--------------------\n";

        return error_message(ss.str());
    }
    return true;
}

template <class F, class... Xs>
    requires((std::is_same_v<NoDerivative, decltype(get_dx_of_dfdx(std::declval<Xs>()...))>))

Maybe_error<bool> test_Derivative(F, double, double, Xs...) {
    auto out = Maybe_error<bool>(true);

    return out;
}

template <class F, class... Xs>
    requires(!(std::is_same_v<NoDerivative, decltype(get_dx_of_dfdx(std::declval<Xs>()...))>))
Maybe_error<bool> test_Derivative(F f, double dx, double eps, const Xs&... xs) {
    using Y = dx_of_dfdx_t<Xs...>;
    decltype(auto) x = get_dx_of_dfdx(xs...);

    auto out = Maybe_error<bool>(true);

    for (std::size_t i = 0; i < x().size(); ++i) {
        auto xi = x;
        xi() = xi() - xi();
        xi()[i] = 1.0;
        auto test_i = test_Derivative(f, xi, dx, eps, xs...);
        if (!test_i)
            out = error_message(out.error()() + "\n at " + std::to_string(i) +
                                "th parameter :" + test_i.error()());
    }
    return out;
}

}  // namespace var

#endif  // DERIVATIVE_TEST_H
