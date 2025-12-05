#ifndef DERIVATIVE_TEST_H
#define DERIVATIVE_TEST_H
#include <derivative_fwd.h>
#include <exponential_matrix.h>
#include <matrix.h>
#include <parameters.h>
#include <variables.h>

#include <algorithm>
#include <cassert>
#include <cstddef>

#include "macrodr/dsl/type_name.h"

// Toggleable asserts for derivative tests. Enable hard checks with
// -DMACRODR_STRICT_DX_ASSERT at compile time.
#ifndef MACRODR_STRICT_DX_ASSERT
#define MACRODR_TEST_ASSERT(cond) ((void)0)
#else
#define MACRODR_TEST_ASSERT(cond) assert(cond)
#endif
#include <cmath>
#include <functional>
#include <iomanip>
#include <ios>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "derivative_operator.h"
#include "maybe_error.h"
#include "parameters_derivative.h"

namespace var {
template <class T>
using underlying_t = std::decay_t<decltype(std::declval<T>()())>;

template <class T>
inline constexpr bool is_matrix_like_v =
is_Matrix_v<std::decay_t<T>> || is_Matrix_v<underlying_t<T>>;



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
    auto const& x = get_dx_of_dfdx(xs...);

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

// Generic norm helper with sane defaults
namespace detail_taylor {

struct MetricInfo {
    double diff{};
    std::string where;
};

template <class T>
struct is_tuple_like : std::false_type {};

template <class... Ts>
struct is_tuple_like<std::tuple<Ts...>> : std::true_type {};

template <class T>
inline constexpr bool is_tuple_like_v = is_tuple_like<std::decay_t<T>>::value;

template <class Metric, class F>
void for_each_metric(const Metric& metric, F&& f) {
    if constexpr (is_tuple_like_v<Metric>) {
        std::apply([&](auto const&... elems) { (for_each_metric(elems, f), ...); }, metric);
    } else {
        f(metric);
    }
}

template <class Metric>
std::size_t count_metrics(const Metric& metric) {
    std::size_t count = 0;
    for_each_metric(metric, [&](const MetricInfo&) { ++count; });
    return count;
}

template <class Metric, class F>
auto transform_metric(const Metric& metric, F&& f) {
    using Result = std::decay_t<decltype(f(std::declval<const MetricInfo&>()))>;
    std::vector<Result> out;
    out.reserve(count_metrics(metric));
    for_each_metric(metric, [&](const MetricInfo& info) { out.push_back(f(info)); });
    return out;
}

template <class Metric>
std::vector<MetricInfo> collect_metric_info(const Metric& metric) {
    return transform_metric(metric, [](const MetricInfo& info) { return info; });
}

template <class Metric>
double max_metric_diff(const Metric& metric) {
    double max_diff = 0.0;
    for_each_metric(metric,
                    [&](const MetricInfo& info) { max_diff = std::max(max_diff, info.diff); });
    return max_diff;
}

inline double default_nonsmooth_tolerance() {
    return 10.0 * std::sqrt(std::numeric_limits<double>::epsilon());
}

// Extract underlying value to allow safe scalar ops on leaf types (Var/Maybe_error)

template <class T>
concept HasValue = requires(const std::remove_reference_t<T>& t) { t.value(); };

template <class T>
decltype(auto) value_of(T&& y) {
    if constexpr (HasValue<T>)
        return std::forward<T>(y).value();  // preserves ref/value
    else
        return std::forward<T>(y);  // returns T&& or T& as appropriate
}

inline double norm_of(double x) {
    return std::abs(x);
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
inline double norm_of(const aMatrix<double>& M) {
    return norm_1(M);
}

template <class T>
inline double norm_of(const Maybe_error<T>& y) {
    return norm_of(get_value(y));
}

template <class... Ts>
inline double norm_of(const std::tuple<Ts...>& t) {
    double m = 0.0;
    std::apply([&m](auto const&... parts) { ((m = std::max(m, norm_of(parts))), ...); }, t);
    return m;
}

// Overloads for var::Var and var::Vector_Space
template <class Id, class T>
inline double norm_of(const var::Var<Id, T>& v) {
    return norm_of(v());
}

template <class... Vars>
inline double norm_of(const var::Vector_Space<Vars...>& vs) {
    double m = 0.0;
    ((m = std::max(m, norm_of(get<Vars>(vs)()))), ...);
    return m;
}

// Difference helper for tuples and matrices with location info
inline MetricInfo metric_default_info(const double& a, const double& b) {
    return {std::abs(a - b), std::string{"scalar"}};
}
inline double metric_default(const double& a, const double& b) {
    return metric_default_info(a, b).diff;
}

template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
inline MetricInfo metric_default_info(const aMatrix<double>& A, const aMatrix<double>& B) {
    double res = norm_1(A - B);
    double maxabs = 0.0;
    std::size_t mi = 0, mj = 0;
    for (std::size_t i = 0; i < A.nrows(); ++i)
        for (std::size_t j = 0; j < A.ncols(); ++j) {
            double d = std::abs(A(i, j) - B(i, j));
            if (d > maxabs) {
                maxabs = d;
                mi = i;
                mj = j;
            }
        }
    std::ostringstream where;
    where << "Matrix(" << mi << "," << mj << ")";
    return {res, where.str()};
}
template <template <class> class aMatrix>
    requires aMatrix<double>::is_Matrix
inline double metric_default(const aMatrix<double>& A, const aMatrix<double>& B) {
    return metric_default_info(A, B).diff;
}

template <class... Ts>
inline MetricInfo metric_default_info(const std::tuple<Ts...>& A, const std::tuple<Ts...>& B) {
    MetricInfo best{0.0, {}};
    std::size_t idx = 0;
    std::apply(
        [&](auto const&... a_parts) {
            (void)std::initializer_list<int>{(
                [&] {
                    auto child = metric_default_info(a_parts, std::get<idx>(B));
                    if (child.diff > best.diff) {
                        std::ostringstream where;
                        where << "tuple[" << idx << "]";
                        if (!child.where.empty())
                            where << " -> " << child.where;
                        best.diff = child.diff;
                        best.where = where.str();
                    }
                    ++idx;
                }(),
                0)...};
        },
        A);
    return best;
}
template <class... Ts>
inline double metric_default(const std::tuple<Ts...>& A, const std::tuple<Ts...>& B) {
    return metric_default_info(A, B).diff;
}

template <class Id, class T>
inline MetricInfo metric_default_info(const var::Var<Id, T>& a, const var::Var<Id, T>& b) {
    auto child = metric_default_info(a(), b());
    std::ostringstream where;
    where << type_name<Id>();
    if (!child.where.empty())
        where << " -> " << child.where;
    return {child.diff, where.str()};
}
template <class Id, class T>
inline double metric_default(const var::Var<Id, T>& a, const var::Var<Id, T>& b) {
    return metric_default_info(a, b).diff;
}

template <class T>
inline MetricInfo metric_default_info(const Maybe_error<T>& a, const Maybe_error<T>& b) {
    return metric_default_info(get_value(a), get_value(b));
}
template <class T>
inline double metric_default(const Maybe_error<T>& a, const Maybe_error<T>& b) {
    return metric_default_info(a, b).diff;
}

template <class Metric>
auto prefix_metric_where(Metric metric, const std::string& prefix) {
    if constexpr (is_tuple_like_v<Metric>) {
        return std::apply(
            [&](auto&&... elems) {
                return std::tuple{
                    prefix_metric_where(std::forward<decltype(elems)>(elems), prefix)...};
            },
            std::move(metric));
    } else {
        if (!prefix.empty()) {
            if (metric.where.empty())
                metric.where = prefix;
            else
                metric.where = prefix + " -> " + metric.where;
        }
        return metric;
    }
}

template <class Comp>
auto vector_component_metric_info(const Comp& a, const Comp& b) {
    auto prefix = type_name<Comp>();
    if constexpr (requires { metric_default_info(a, b); }) {
        auto result = metric_default_info(a, b);
        return prefix_metric_where(std::move(result), prefix);
    } else {
        auto result = metric_default_info(a(), b());
        return prefix_metric_where(std::move(result), prefix);
    }
}

template <class... Vars>
inline auto metric_default_info(const var::Vector_Space<Vars...>& A,
                                const var::Vector_Space<Vars...>& B) {
    return std::tuple(vector_component_metric_info(get<Vars>(A), get<Vars>(B))...);
}

template <class... Vars>
inline double metric_default(const var::Vector_Space<Vars...>& A,
                             const var::Vector_Space<Vars...>& B) {
    auto metrics = metric_default_info(A, B);
    return max_metric_diff(metrics);
}

}  // namespace detail_taylor

inline std::string print_delta_parameters(const var::Parameters_Transformations& x,
                                          const Matrix<double>& dx) {
    std::stringstream ss;

    for (std::size_t i = 0; i < dx.size(); ++i) {
        if (dx[i] != 0.0) {
            ss << "direction " << i << " (" << x.names()[i] << ") : dx=" << dx[i] << "\n";
        }
    }
    return ss.str();
}

inline Maybe_error<bool> concatenate_error(Maybe_error<bool> one, Maybe_error<bool> two) {
    if (one.valid() && two.valid()) {
        return {true};
    }
    return error_message(one.error()() + two.error()());
}

template <class... Id>
var::Vector_Space<var::Var<Id, Maybe_error<bool>>...> concatenate_error(
    var::Vector_Space<var::Var<Id, Maybe_error<bool>>...>&& one,
    var::Vector_Space<var::Var<Id, Maybe_error<bool>>...>&& two) {
    var::Vector_Space<var::Var<Id, Maybe_error<bool>>...> out;
    ((get<Id>(out) = concatenate_error(get<Id>(std::move(one))(), get<Id>(std::move(two))())), ...);
    return out;
}

template <class Id>
var::Var<Id, Maybe_error<bool>> concatenate_error(var::Var<Id, Maybe_error<bool>>&& one,
                                                  var::Var<Id, Maybe_error<bool>>&& two) {
    if (std::move(one)().valid() && std::move(two)().valid()) {
        return var::Var<Id, Maybe_error<bool>>(true);
    }

    return error_message(std::move(one)().error()() + std::move(two)().error()());
}
template <class Id>
Maybe_error<bool> condense_errors(var::Var<Id, Maybe_error<bool>> ok) {
    if (ok().valid()) {
        return {true};
    }
    return error_message(type_name<Id>(), ": ", std::move(ok)().error()());
}

template <class... Id>
Maybe_error<bool> condense_errors(var::Vector_Space<var::Var<Id, Maybe_error<bool>>...>&& ok) {
    if ((true && ... && get<Id>(ok)().valid())) {
        return {true};
    }
    return error_message(condense_errors(get<Id>(std::move(ok))).error()()...);
}

template <bool inspect = false, class T, class X>
    requires(!is_of_this_template_type_v<T, Vector_Space>)
var::Var<T, Maybe_error<bool>> test_clarke_brackets(T const&, T const&, const X&, double, const T&,
                                                    const T&) {
    return var::Var<T, Maybe_error<bool>>(true);
}




// --- helpers (adapt rows/cols/indexing to your T if needed) ------------------
template <class M>
inline std::pair<int, int> argmax_with_index(const M& A) {
    int bi = 0, bj = 0;
    auto best = A(0ul, 0ul);
    for (int i = 0; i < A.nrows(); ++i)
        for (int j = 0; j < A.ncols(); ++j)
            if (A(i, j) > best) {
                best = A(i, j);
                bi = i;
                bj = j;
            }
    return {bi, bj};
}
template <class M>
inline int count_positive(const M& A) {
    int c = 0;
    for (int i = 0; i < A.nrows(); ++i)
        for (int j = 0; j < A.ncols(); ++j)
            if (A(i, j) > 0)
                ++c;
    return c;
}

template <bool inspect = false, class T, class X>
requires(!is_of_this_template_type_v<T, Vector_Space> &&
is_matrix_like_v<T>)    
var::Var<T, Maybe_error<bool>> test_clarke_brackets(const Derivative<T, X>& Y0, const T& Y,
                                                    const X& u, double t, const T& Yp,
                                                    const T& Yn) {
    // ================================================================
    // Clarke brackets test
    // ---------------------------------------------------------------
    // For a differentiable function f : R^n → R and point x,
    // the directional derivative along a unit vector u satisfies:
    //
    //   D⁻f(x;u) ≤ f'(x)·u ≤ D⁺f(x;u)
    //
    // where
    //   D⁺f(x;u) ≈ (f(x + t u) - f(x)) / t      // forward difference
    //   D⁻f(x;u) ≈ (f(x) - f(x - t u)) / t      // backward difference
    //
    // The analytic gradient g is *consistent* if for small t:
    //
    //   D⁻ ≤ g·u ≤ D⁺
    //
    // Violations of this inequality indicate either:
    //   (a) f is not differentiable at x, or
    //   (b) the analytic gradient g is incorrect.
    // =========================================================

    MACRODR_TEST_ASSERT(Y0.has_dx() && "test_clarke_brackets: derivative missing dx pointer");
    auto g = derivative(Y0);

    auto Ydiff = Y() - Y0.primitive()();

    auto y_norm = maxAbs(Y());
    auto ydiff_norm = maxAbs(Ydiff);
    // Accept exact match or relative closeness; the zero-baseline case should pass.
    auto Ydiff_ok = (ydiff_norm <= eps) || (ydiff_norm <= y_norm * eps);

    auto h = t * u();
    auto J_u = g * u;  // directional derivative from analytic gradient

    auto D_plus = (Yp() - Y0.primitive()()) / t;
    auto D_minus = (Y0.primitive()() - Yn()) / t;

    auto mid = 0.5 * (D_plus + D_minus);
    auto width = zip(
        [](auto a, auto b) {
            using std::abs;
            return abs(a - b);
        },
        D_plus, D_minus);


    constexpr double M = 10.0;
    constexpr double rel_tol = 1e-6;
    auto tau_abs = (apply(
                        [](auto x) {
                            using std::abs;
                            return 1 + abs((double)x);
                        },
                        Yp()) +
                    apply(
                        [](auto x) {
                            using std::abs;
                            return 1 + abs((double)x);
                        },
                        Yn()) +
                    apply(
                        [](auto x) {
                            using std::abs;
                            return 1 + abs((double)x);
                        },
                        Y0.primitive()())) *
                   (4.0 * sqrt(std::numeric_limits<double>::epsilon()) / std::abs(t));

    auto band = M * (width * 0.5) + tau_abs +
                rel_tol * (apply(
                               [](auto x) {
                                   using std::abs;
                                   return abs((double)x);
                               },
                               J_u) +
                           apply(
                               [](auto x) {
                                   using std::abs;
                                   return abs((double)x);
                               },
                               mid));
    auto gap = zip(
                   [](auto a, auto b) {
                       using std::abs;
                       return abs(a - b);
                   },
                   J_u, mid) -
               band;

    bool ok = true;

    for (std::size_t i = 0; i < gap.size(); ++i) {
        // if the directional derivative change sign we accept the test
        if (gap[i] > 0 && (D_plus[i] * D_minus[i] > 0)) {
            ok = false;
        }
    }

    if (!ok || inspect || !Ydiff_ok) {
        // residual/gap: positive means violation

        // summarize failures
        int fail_count = count_positive(gap);
        auto [wi, wj] = argmax_with_index(gap);
        double worst = gap(wi, wj);

        // optional central check for context
        auto C = (Yp() - Yn()) / (2 * t);
        auto cen_err = zip(
            [](auto a, auto b) {
                using std::abs;
                return abs(a - b);
            },
            J_u, C);

        std::stringstream ss;
        ss << std::setprecision(10);
        ss << "-------------------------------------------------------------------------------"
              "-";
        ss << "-------------------------ERROR ";

        ss << "Y()\n" << Y() << "\n";
        ss << "Y0.primitive()()\n" << Y0.primitive()() << "\n";
        ss << "\nYdiff\n" << Ydiff;

        if (Y0.has_dx())
            ss << print_delta_parameters(Y0.dx().parameters(), h) << "-----------------------\n";
        else
            ss << "[missing dx pointer]" << "-----------------------\n";

        ss << "||h||=" << detail_taylor::norm_of(h) << "\n";

        // (1) closeness diagnostics
        ss << "D_minus=\n" << D_minus << "\n";
        ss << "J_u=\n" << J_u << "\n";
        ss << "D_plus=\n" << D_plus << "\n";

        // (2) failure location and local details
        ss << "fail_count=" << fail_count << " (of " << (gap.nrows() * gap.ncols()) << ")\n";
        ss << "worst_idx=(" << wi << "," << wj << "), gap=" << worst
           << "  // gap = |J_u - mid| - band\n";
        ss << "at worst idx:\n";
        ss << "  Ju=" << J_u(wi, wj) << "  D-=" << D_minus(wi, wj) << "  D+=" << D_plus(wi, wj)
           << "\n";
        ss << "  mid=" << mid(wi, wj) << "  width=" << width(wi, wj) << "  band=" << band(wi, wj)
           << "\n";

        // optional: central-diff info at the same spot
        ss << "  C=" << C(wi, wj) << "  |Ju-C|=" << cen_err(wi, wj) << "\n";

        ss << "------------------------------------------------------------------\n\n\n";
        return error_message(ss.str());
    }
    return ok;
}


template <bool inspect = false, class T, class X>
requires(!is_of_this_template_type_v<T, Vector_Space> &&
  !is_of_this_template_type_v<decltype(std::declval<T>()()), Vector_Space> && 
!is_matrix_like_v<T>)    
var::Var<T, Maybe_error<bool>> test_clarke_brackets(const Derivative<T, X>& Y0, const T& Y,
                                                    const X& u, double t, const T& Yp,
                                                    const T& Yn) {
    // ================================================================
    // Clarke brackets test
    // ---------------------------------------------------------------
    // For a differentiable function f : R^n → R and point x,
    // the directional derivative along a unit vector u satisfies:
    //
    //   D⁻f(x;u) ≤ f'(x)·u ≤ D⁺f(x;u)
    //
    // where
    //   D⁺f(x;u) ≈ (f(x + t u) - f(x)) / t      // forward difference
    //   D⁻f(x;u) ≈ (f(x) - f(x - t u)) / t      // backward difference
    //
    // The analytic gradient g is *consistent* if for small t:
    //
    //   D⁻ ≤ g·u ≤ D⁺
    //
    // Violations of this inequality indicate either:
    //   (a) f is not differentiable at x, or
    //   (b) the analytic gradient g is incorrect.
    // =========================================================

    MACRODR_TEST_ASSERT(Y0.has_dx() && "test_clarke_brackets: derivative missing dx pointer");
   // Analytic directional derivative J[u]
auto g = derivative(Y0);      // d_d_<double, X> for T = double
auto h = t * u();
auto J_u = g * u;           // dot with direction (see parameters_derivative.h)

// Forward/backward diffs
auto f0 = Y0.primitive()();
auto D_plus  = (Yp() - f0) / t;
auto D_minus = (f0 - Yn()) / t;

// Basic closeness guard (same as matrix flavor)
auto Ydiff = Y() - f0;
auto y_norm = std::abs(Y());
auto ydiff_norm = std::abs(Ydiff);
// Accept exact match or relative closeness; handle zero-baseline gracefully.
auto Ydiff_ok = (ydiff_norm <= eps) || (ydiff_norm <= y_norm * eps);

// Bracket band (scalar)
auto mid   = 0.5 * (D_plus + D_minus);
auto width = std::abs(D_plus - D_minus);

auto M = 10.0;
auto rel_tol = 1e-5;
auto abs_tol = 1e-5;

auto band = M * (0.5 * width) +  rel_tol * (std::abs(J_u) + std::abs(mid))+ abs_tol;
auto gap  = std::abs(J_u - mid) - band;

// Accept if inside band, or if the directional slope changes sign (non-smooth, acceptable)
bool ok = !(gap > 0.0 && (D_plus * D_minus > 0.0));
if ((!ok || inspect || !Ydiff_ok )&&!std::isnan(mid)) {
    // Optional central check for context
    auto C = (Yp() - Yn()) / (2.0 * t);
    auto cen_err = std::abs(J_u - C);

    std::stringstream ss;
    ss << std::setprecision(10);
    ss << "---------------- Clarke (scalar) ----------------\n";
    ss << "Y()=" << Y() << "  Y0=" << f0 << "  Ydiff=" << Ydiff << "\n";
    if (Y0.has_dx())
        ss << print_delta_parameters(Y0.dx().parameters(), h) << "-----------------------\n";
    else
        ss << "[missing dx pointer]-----------------------\n";
    ss << "||h||=" << detail_taylor::norm_of(h) << "\n";
    ss << "D-=" << D_minus << "  J[u]=" << J_u << "  D+=" << D_plus << "\n";
    ss << "mid=" << mid << "  width=" << width << "  band=" << band << "  gap=" << gap << "\n";
    ss << "C=" << C << "  |J-C|=" << cen_err << "\n";
    return error_message(ss.str());
}
return ok;
}



template <bool inspect = false, class T, class X>
    requires(!is_of_this_template_type_v<T, Vector_Space>)
var::Var<T, Maybe_error<bool>> test_clarke_brackets_old(const Derivative<T, X>& Y0, const X& u,
                                                        double t, const T& Yp, const T& Yn) {
    // ================================================================
    // Clarke brackets test
    // ---------------------------------------------------------------
    // For a differentiable function f : R^n → R and point x,
    // the directional derivative along a unit vector u satisfies:
    //
    //   D⁻f(x;u) ≤ f'(x)·u ≤ D⁺f(x;u)
    //
    // where
    //   D⁺f(x;u) ≈ (f(x + t u) - f(x)) / t      // forward difference
    //   D⁻f(x;u) ≈ (f(x) - f(x - t u)) / t      // backward difference
    //
    // The analytic gradient g is *consistent* if for small t:
    //
    //   D⁻ ≤ g·u ≤ D⁺
    //
    // Violations of this inequality indicate either:
    //   (a) f is not differentiable at x, or
    //   (b) the analytic gradient g is incorrect.
    // =========================================================

    MACRODR_TEST_ASSERT(Y0.has_dx() && "test_clarke_brackets: derivative missing dx pointer");
    auto g = derivative(Y0);
    auto h = t * u();
    auto J_u = g * u;  // directional derivative from analytic gradient

    auto D_plus = (Yp() - Y0.primitive()()) / t;
    auto D_minus = (Y0.primitive()() - Yn()) / t;
    auto lo = zip([](auto x, auto y) { return std::min(x, y); }, D_minus, D_plus);
    auto hi = zip([](auto x, auto y) { return std::max(x, y); }, D_minus, D_plus);

    auto C = (Yp() - Yn()) / (2 * t);  // central check
    constexpr const double abs_tol = 1e-10;
    constexpr const double rel_tol = 1e-6;

    auto mixed_tol = [&](auto a, auto b, auto c) {
        return abs_tol + rel_tol * (apply(abs, a) + apply(abs, b) + apply(abs, c));
    };

    // elementwise
    auto gap_low = J_u - lo;
    auto gap_high = hi - J_u;

    // pass if both sides within tolerance
    auto tol = mixed_tol(J_u, D_minus, D_plus);
    auto pass_bracket = (min(gap_low + tol) >= 0) && (min(gap_high + tol) >= 0);

    // central difference consistency
    auto cen_tol = mixed_tol(J_u, C, C);
    auto pass_central = (max(apply(abs, J_u - C) - cen_tol) <= 0);

    // final decision (single t)
    bool ok = pass_bracket || pass_central;

    if (!ok || inspect) {
        std::stringstream ss;
        ss << std::setprecision(10);
        ss << "--------------------------------------------------------------------------------";
        ss << "-------------------------ERROR ";
        if (Y0.has_dx())
            ss << print_delta_parameters(Y0.dx().parameters(), h) << "-----------------------\n";
        else
            ss << "[missing dx pointer]" << "-----------------------\n";
        ss << "||h||=" << detail_taylor::norm_of(h) << "\t";
        ss << "D_minus=" << D_minus << "\t";
        ss << "J_u=" << J_u << "\t";
        ss << "gap_low + tol " << gap_low + tol;
        ss << "gap_high + tol " << gap_high + tol;
        ss << "cen_tol" << cen_tol;

        ss << "------------------------------------------------------------------\n\n\n";
        return error_message(ss.str());
    }
    return ok;
}
template <bool inspect = false, class T, class X>
requires(is_of_this_template_type_v<decltype(std::declval<T>()()), Vector_Space>)
auto test_clarke_brackets(
    const Derivative<T, X>& Y0, const T& Y, const X& u, double t,
    const T& Yp, const T& Yn) {
    auto inner = test_clarke_brackets<inspect>(Y0(), Y(), u, t, Yp(), Yn());
    // Collapse nested Vector_Space error details into a single Maybe_error<bool>
    auto condensed = condense_errors(std::move(inner));
    return var::Var<T, Maybe_error<bool>>(std::move(condensed));
}



template <bool inspect = false, class... T, class X>
auto test_clarke_brackets(
    const Derivative<Vector_Space<T...>, X>& Y0, const Vector_Space<T...>& Y, const X& u, double t,
    const Vector_Space<T...>& Yp, const Vector_Space<T...>& Yn) {
    return var::Vector_Space(
        test_clarke_brackets<inspect>(get<T>(Y0), get<T>(Y), u, t, get<T>(Yp), get<T>(Yn))...);
}

template <class... T, class X>
auto test_clarke_brackets_init(
    const Derivative<Vector_Space<T...>, X>& x) {
    return var::Vector_Space(test_clarke_brackets_init(get<T>(x))...);
}

// Vars whose value() is a Vector_Space should still be treated as atomic
// variables for initialization so we don't try to build nested Vector_Space
// components (e.g., Patch_State wraps a Vector_Space<P_mean, P_Cov>).
template <class T, class X>
    requires(requires { std::decay_t<T>::is_variable; } &&
             is_of_this_template_type_v<decltype(std::declval<T>()()), Vector_Space>)
var::Var<T, Maybe_error<bool>> test_clarke_brackets_init(const Derivative<T, X>&) {
    return var::Var<T, Maybe_error<bool>>(true);
}

template <class T, class X>
requires(!is_of_this_template_type_v<T, Vector_Space>&& !is_of_this_template_type_v<decltype(std::declval<T>()()), Vector_Space>)
var::Var<T, Maybe_error<bool>> test_clarke_brackets_init(const Derivative<T, X>&) {
    return var::Var<T, Maybe_error<bool>>(true);
}

template <class T>
    requires(requires { std::decay_t<T>::is_variable; } &&
             !is_of_this_template_type_v<std::decay_t<T>, Vector_Space>)
var::Var<T, Maybe_error<bool>> test_clarke_brackets_init(const T&) {
    return var::Var<T, Maybe_error<bool>>(true);
}



template <class T, class X>
    requires(is_of_this_template_type_v<decltype(std::declval<T>()()), Vector_Space> &&
             !(requires { std::decay_t<T>::is_variable; }))
auto  test_clarke_brackets_init(const Derivative<T, X>& x)  {
    return test_clarke_brackets_init(x());
}



// Taylor convergence test using norm-based residuals and known bound on third derivative
// Accept if || f(x+h p) - (f(x) + J(x)[p] h) ||  <= max_third_deriv_norm * ||h||^3 /           6
// where J(x)[p] is the directional derivative in direction p
// Note: this is a much stronger test than the step-doubling test
// but requires knowledge of a bound on the third derivative

template <bool inspect = false, class F, class... Xs>
    requires(!std::is_same_v<NoDerivative, dx_of_dfdx_t<Xs...>> && std::is_invocable_v<F, Xs...>&& 
            is_derivative_v<decltype(get_value(std::invoke(std::declval<F>(), std::declval<Xs>()...)))>)
Maybe_error<bool> test_derivative_clarke(F f, double h, const Xs&... xs) {
    auto MaybedY = std::invoke(f, xs...);
    if (!is_valid(MaybedY)) {
        return get_error(MaybedY);
    }
    auto dY = std::move(get_value(MaybedY));
    static_assert(is_derivative_v<decltype(dY)>);
    auto MaybeY = std::invoke(f, primitive(xs)...);
    if (!is_valid(MaybeY)) {
        return get_error(MaybeY);
    }
    auto Y = std::move(get_value(MaybeY));

    // Direction space from the first derivative argument among xs...
    using Y_type = dx_of_dfdx_t<Xs...>;
    if constexpr (std::is_same_v<NoDerivative, Y_type>) {
        // Nothing to test if there is no derivative argument
        return true;
    }
    auto xdir = get_dx_of_dfdx(xs...);
    auto ok = test_clarke_brackets_init(dY);
    //using show_type = typename decltype(ok)::show_type;
    //auto s = show_type {}
    for (std::size_t i = 0; i < xdir().size(); ++i) {
        auto p = xdir;
        // Zero-initialize direction without relying on operator-
        {
            auto tmp = p();
            for (std::size_t k = 0; k < tmp.size(); ++k) {
                tmp[k] = 0.0;
            }
            p() = std::move(tmp);
        }
        p()[i] = 1.0;
        auto xs_orginal = std::tuple(xs...);
        auto xs_taylor_pos = std::tuple(Taylor_first(xs, p, h)...);
        auto xs_taylor_neg = std::tuple(Taylor_first(xs, p, -h)...);

        auto Maybe_Yp = std::invoke(f, Taylor_first(xs, p, h)...);
        auto Maybe_Yn = std::invoke(f, Taylor_first(xs, p, -h)...);

        auto Tp = Taylor_first(dY, p, h);
        auto Tn = Taylor_first(dY, p, -h);

        if (!is_valid(Maybe_Yp)) {
            return get_error(Maybe_Yp);
        }
        if (!is_valid(Maybe_Yn)) {
            return get_error(Maybe_Yn);
        }
        auto Y_p = std::move(get_value(Maybe_Yp));
        auto Y_n = std::move(get_value(Maybe_Yn));

        // Work with primitives for Clarke brackets so the vector-space
        // overloads see a consistent value type even if f returns
        // Derivative<...> in some instantiations.
        auto Y_p_prim = primitive(Y_p);
        auto Y_n_prim = primitive(Y_n);
        auto Y_prim = primitive(Y);

        auto diff_p = Y_p_prim - Tp;
        auto diff_n = Y_n_prim - Tn;
        auto new_ok =
            test_clarke_brackets<inspect>(dY, Y_prim, p, h, Y_p_prim, Y_n_prim);
        ok = concatenate_error(std::move(ok), std::move(new_ok));
    }
    return condense_errors(std::move(ok));
}

template <bool inspect = false, class F, class... Xs>
    requires((std::is_same_v<NoDerivative, decltype(get_dx_of_dfdx(std::declval<Xs>()...))>))
Maybe_error<bool> test_derivative_clarke(F f, double h, const Xs&... xs) {
    auto out = Maybe_error<bool>(true);

    return out;
}

template <class T, class X>
    requires(!is_of_this_template_type_v<T, Vector_Space>)
var::Var<T, Maybe_error<bool>> test_taylor_error(const Derivative<T, X>& Y0,
                                                 const Derivative<T, X>& Y1,
                                                 double max_third_deriv_norm = 10) {
    // Taylor error test for second-order Taylor expansion

    //  f(x+h) = f(x) + J(x) h + 1/2 h^T H(x) h + O(||h||^3)

    // now, we have J(x+h) = J(x) + H(x) h + O(||h||^2)
    // so we can estimate H(x) h as J(x+h) - J(x)
    // then we can estimate the second-order term as 1/2 h^T (J(x+h) - J(x))
    // and compare f(x+h) to f(x) + J(x) h + 1/2 h^T (J(x+h) - J(x))
    //  Error = f(x+h) - (f(x) + J(x) h + 1/2 h^T (J(x+h) - J(x))) = O(||h||^3)
    MACRODR_TEST_ASSERT(Y0.has_dx() && Y1.has_dx() &&
                        "test_taylor_error: derivative missing dx pointer");
    assert(Y0.has_dx() && Y1.has_dx());
    auto h = Y0.dx()() - Y1.dx()();
    auto J_h = Y1.derivative()() - Y0.derivative()();
    auto half_hT_Jh = 0.5 * TranspMult(h, J_h);
    auto Taylor_2nd = primitive(Y0)() + Y0.derivative()() * h + half_hT_Jh;
    auto Error = primitive(Y1)() - Taylor_2nd;
    auto err = detail_taylor::norm_of(Error);
    auto hnorm = detail_taylor::norm_of(h);
    double third_deriv_norm = err / (hnorm * hnorm * hnorm) / 6.0;
    bool ok = third_deriv_norm < max_third_deriv_norm;
    if (!ok) {
        std::stringstream ss;
        ss << std::setprecision(10);
        ss << print_delta_parameters(Y0.dx().parameters(), h) << "\t";
        ss << "||h||=" << hnorm << "\t";
        ss << "Error=" << err << "\t";
        ss << "Third Derivative norm, Error / ||h||^3/6 = " << third_deriv_norm << "\n";
        return error_message(ss.str());
    }
    return ok;
}
template <class... Vs, class X>
Vector_Space<var::Var<typename Vs::Id, Maybe_error<bool>>...> test_taylor_error(
    const Derivative<var::Vector_Space<Vs...>, X>& Y0,
    const Derivative<var::Vector_Space<Vs...>, X>& Y1, double max_third_deriv_norm = 10.0) {
    return Vector_Space(test_taylor_error(get<typename Vs::Id>(Y0), get<typename Vs::Id>(Y1),
                                          max_third_deriv_norm)...);
}

// Taylor convergence test using norm-based residuals and known bound on third derivative
// Accept if || f(x+h p) - (f(x) + J(x)[p] h) ||  <= max_third_deriv_norm * ||h||^3 /           6
// where J(x)[p] is the directional derivative in direction p
// Note: this is a much stronger test than the step-doubling test
// but requires knowledge of a bound on the third derivative

template <class F, class... Xs>
    requires(!std::is_same_v<NoDerivative, dx_of_dfdx_t<Xs...>> && std::is_invocable_v<F, Xs...>)
Maybe_error<bool> test_taylor_convergence_third(F f, double h, double max_third_deriv_norm,
                                                const Xs&... xs) {
    auto MaybedY = std::invoke(f, xs...);
    auto MaybeY0 = std::invoke(f, primitive(xs)...);
    if (!is_valid(MaybedY)) {
        return get_error(MaybedY);
    }
    if (!is_valid(MaybeY0)) {
        return get_error(MaybeY0);
    }
    auto dY = std::move(get_value(MaybedY));
    auto Y0 = std::move(get_value(MaybeY0));

    // Direction space from the first derivative argument among xs...
    using Y = dx_of_dfdx_t<Xs...>;
    if constexpr (std::is_same_v<NoDerivative, Y>) {
        // Nothing to test if there is no derivative argument
        return true;
    } else {
        auto xdir = get_dx_of_dfdx(xs...);
        auto ok = decltype(test_taylor_error(dY, dY, max_third_deriv_norm)){};
        //using show_type = typename decltype(ok)::show_type;
        //auto s = show_type {}
        for (std::size_t i = 0; i < xdir().size(); ++i) {
            auto p = xdir;
            // Zero-initialize direction without relying on operator-
            {
                auto tmp = p();
                for (std::size_t k = 0; k < tmp.size(); ++k) {
                    tmp[k] = 0.0;
                }
                p() = std::move(tmp);
            }
            p()[i] = 1.0;
            auto MaybedY_h = std::invoke(f, Taylor_first(xs, p, h)...);
            if (!is_valid(MaybedY_h)) {
                return get_error(MaybedY_h);
            }
            auto dY_h = std::move(get_value(MaybedY_h));
            auto new_ok = test_taylor_error(dY, dY_h, max_third_deriv_norm);
            ok = concatenate_error(std::move(ok), std::move(new_ok));
        }
        return condense_errors(std::move(ok));
    }
}

// Taylor convergence test using norm-based residuals and step-doubling.
// Accept if R(h/2) <= R(h)/ratio_target (default ~3.0 slack for ideal 4x reduction),
// where R(h) = || f(x+h p) - (f(x) + J(x)[p] h) || / (1 + ||f(x)||).

template <class F, class... Xs>
    requires(!std::is_same_v<NoDerivative, dx_of_dfdx_t<Xs...>> && std::is_invocable_v<F, Xs...>)
Maybe_error<bool> test_taylor_convergence_chat_gpt(F f, double h, double ratio_target,
                                                   double nonsmooth_tol, const Xs&... xs) {
    using detail_taylor::collect_metric_info;
    using detail_taylor::norm_of;
    using detail_taylor::transform_metric;
    using detail_taylor::value_of;

    auto MaybedY = std::invoke(f, xs...);
    auto MaybeY0 = std::invoke(f, primitive(xs)...);
    if (!is_valid(MaybedY))
        return get_error(MaybedY);
    if (!is_valid(MaybeY0))
        return get_error(MaybeY0);

    auto dY = std::move(get_value(MaybedY));
    auto Y0 = std::move(get_value(MaybeY0));

    // Baseline norm
    const double base = 1.0 + norm_of(Y0);

    // Direction space from the first derivative argument among xs...
    using Y = dx_of_dfdx_t<Xs...>;
    if constexpr (std::is_same_v<NoDerivative, Y>) {
        // Nothing to test if there is no derivative argument
        return true;
    } else {
        auto xdir = get_dx_of_dfdx(xs...);
        auto ok = Maybe_error<bool>(true);
        auto append_error = [&](const std::string& msg) {
            if (ok)
                ok = error_message(msg);
            else
                ok = error_message(ok.error()() + "\n" + msg);
        };

        for (std::size_t i = 0; i < xdir().size(); ++i) {
            auto p = xdir;
            // Zero-initialize direction without relying on operator-
            {
                auto tmp = p();
                for (std::size_t k = 0; k < tmp.size(); ++k) tmp[k] = 0.0;
                p() = std::move(tmp);
            }
            p()[i] = 1.0;

            // First-order predictions at h and h/2
            auto T_h = Taylor_first(dY, p, h);
            auto T_h2 = Taylor_first(dY, p, h * 0.5);

            // True values at perturbed points
            auto MaybeY_h = std::invoke(f, Taylor_first(xs, p, h)...);
            if (!is_valid(MaybeY_h))
                return get_error(MaybeY_h);
            auto Y_h = std::move(get_value(MaybeY_h));

            auto MaybeY_h2 = std::invoke(f, Taylor_first(xs, p, h * 0.5)...);
            if (!is_valid(MaybeY_h2))
                return get_error(MaybeY_h2);
            auto Y_h2 = std::move(get_value(MaybeY_h2));

            // Residuals (with diagnostic location)
            auto M1 = detail_taylor::metric_default_info(Y_h, T_h);
            auto M2 = detail_taylor::metric_default_info(Y_h2, T_h2);
            auto metrics1 = collect_metric_info(M1);
            auto metrics2 = collect_metric_info(M2);
            if (metrics1.size() != metrics2.size()) {
                return error_message("metric size mismatch in test_taylor_convergence");
            }
            auto R1 = transform_metric(
                M1, [base](const detail_taylor::MetricInfo& info) { return info.diff / base; });
            auto R2 = transform_metric(
                M2, [base](const detail_taylor::MetricInfo& info) { return info.diff / base; });

            // If we are at floating-point floor, accept (typical for linear maps).
            constexpr double eps = std::numeric_limits<double>::epsilon();
            const double floor_abs = 1000.0 * eps;  // generous floor on normalized residuals
            auto floor_ok =
                std::all_of(R1.begin(), R1.end(),
                            [floor_abs](double r) { return r <= floor_abs; }) &&
                std::all_of(R2.begin(), R2.end(), [floor_abs](double r) { return r <= floor_abs; });
            if (floor_ok)
                continue;

            // Otherwise, check reduction; allow slack for nonideal behavior.
            std::vector<std::size_t> failing_indices;
            failing_indices.reserve(R1.size());
            for (std::size_t idx = 0; idx < R1.size(); ++idx) {
                double rhs = R1[idx] / ratio_target;
                if (!(R2[idx] <= rhs))
                    failing_indices.push_back(idx);
            }
            if (failing_indices.empty())
                continue;

            // Gather slope comparisons to distinguish non-smooth behavior from derivative bugs
            auto MaybeY_hm = std::invoke(f, Taylor_first(xs, p, -h)...);
            if (!is_valid(MaybeY_hm))
                return get_error(MaybeY_hm);
            auto Y_hm = std::move(get_value(MaybeY_hm));

            auto ad_slope = (value_of(T_h) - value_of(Y0)) / h;
            auto forward_slope = (value_of(Y_h) - value_of(Y0)) / h;
            auto backward_slope = (value_of(Y0) - value_of(Y_hm)) / h;
            auto central_slope = (value_of(Y_h) - value_of(Y_hm)) / (2.0 * h);

            auto diff_forward = detail_taylor::metric_default_info(ad_slope, forward_slope);
            auto diff_backward = detail_taylor::metric_default_info(ad_slope, backward_slope);
            auto diff_central = detail_taylor::metric_default_info(ad_slope, central_slope);

            auto forward_diffs = collect_metric_info(diff_forward);
            auto backward_diffs = collect_metric_info(diff_backward);
            auto central_diffs = collect_metric_info(diff_central);

            if (forward_diffs.size() != metrics1.size() ||
                backward_diffs.size() != metrics1.size() ||
                central_diffs.size() != metrics1.size()) {
                return error_message("slope metric size mismatch in test_taylor_convergence");
            }

            double forward_norm = 1.0 + norm_of(forward_slope);
            double backward_norm = 1.0 + norm_of(backward_slope);
            double central_norm = 1.0 + norm_of(central_slope);

            std::string direction_name{};
            bool has_direction_name = false;
            if constexpr (requires {
                              xdir.parameters();
                              xdir.parameters().names();
                          }) {
                const auto& params = xdir.parameters();
                const auto& names = params.names();
                if (i < names.size()) {
                    direction_name = names[i];
                    has_direction_name = true;
                }
            }

            bool recorded_failure = false;
            for (std::size_t idx : failing_indices) {
                double forward_resid = forward_diffs[idx].diff / forward_norm;
                double backward_resid = backward_diffs[idx].diff / backward_norm;
                double central_resid = central_diffs[idx].diff / central_norm;

                if (forward_resid <= nonsmooth_tol || backward_resid <= nonsmooth_tol ||
                    central_resid <= nonsmooth_tol) {
                    continue;
                }

                std::ostringstream ss;
                ss << std::setprecision(15);
                ss << "Taylor residual did not reduce sufficiently in direction " << i;
                if (has_direction_name)
                    ss << " (" << direction_name << ")";
                ss << " at component " << metrics2[idx].where;
                ss << " : R(h)=" << R1[idx] << " R(h/2)=" << R2[idx] << " expected <= R(h)/"
                   << ratio_target;
                ss << "\n  slope diff forward (" << forward_diffs[idx].where
                   << ") = " << forward_resid;
                ss << "\n  slope diff backward (" << backward_diffs[idx].where
                   << ") = " << backward_resid;
                ss << "\n  slope diff central (" << central_diffs[idx].where
                   << ") = " << central_resid;
                append_error(ss.str());
                recorded_failure = true;
            }

            if (!recorded_failure)
                continue;
        }
        return ok;
    }
}

// Fallback when no derivative is present among arguments: skip test and succeed.
template <class F, class... Xs>
    requires(std::is_same_v<NoDerivative, dx_of_dfdx_t<Xs...>>)
Maybe_error<bool> test_taylor_convergence(F, double, double, double, const Xs&...) {
    return Maybe_error<bool>(true);
}

template <class F, class... Xs>
Maybe_error<bool> test_taylor_convergence(F f, double h, double ratio_target, const Xs&... xs) {
    return test_taylor_convergence(f, h, ratio_target, detail_taylor::default_nonsmooth_tolerance(),
                                   xs...);
}

}  // namespace var

#endif  // DERIVATIVE_TEST_H
