#ifndef MOMENT_STATISTICS_H
#define MOMENT_STATISTICS_H

#include <maybe_error.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <limits>

#include "matrix.h"
#include "lapack_headers.h"
#include "parameter_indexed.h"
#include "variables.h"

template<class T>
struct value_type_impl {
    using type = std::remove_cvref_t<T>;
};

template<class T>
requires requires (T const& v) { v(); }
struct value_type_impl<T> {
    using type = std::remove_cvref_t<decltype(std::declval<T const&>()())>;
};

template<class T>
using value_type_t = typename value_type_impl<T>::type;






template<typename T>
struct mean_value_type_impl {
    using type = std::conditional_t<std::is_arithmetic_v<std::remove_cvref_t<T>>, double,std::remove_cvref_t<T>>;
};
template<class T>
using mean_value_type_t = typename mean_value_type_impl<T>::type;



template<class T>
requires requires (T const& v) { v(); }
struct mean_value_type_impl<T> {
    using type = mean_value_type_t<value_type_t<T>>;
};


template<bool include_covariance, class T>
auto sqr_X( const std::vector<T>& x, std::size_t max_lag) {
using R= typename std::decay_t<decltype(prod_XY<include_covariance>(std::declval<T>(), std::declval<T>()))>;
    std::vector<std::vector<R>> out;
    out.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        out[i].reserve(std::min(x.size() - i, max_lag));
        for (std::size_t j = i; j < x.size() && j < i + max_lag ; ++j) {
            out[i].push_back(prod_XY<include_covariance>(x[i], x[j]));
        }
    }
    return out;
}

template<class T>
auto nested_add(const T& a, const T& b) {
    return a + b;
}

template<class T>
auto nested_sub(const T& a, const T& b) {
    return a - b;
}

template<class T>
auto nested_scale(const T& a, double s) {
    return a * s;
}

template<class T>
auto nested_div(const T& a, double s) {
    return a / s;
}

template<class T>
auto nested_add(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<std::decay_t<decltype(nested_add(a[0], b[0]))>> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        out.push_back(nested_add(a[i], b[i]));
    }
    return out;
}

template<class T>
auto nested_sub(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<std::decay_t<decltype(nested_sub(a[0], b[0]))>> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        out.push_back(nested_sub(a[i], b[i]));
    }
    return out;
}

template<class T>
auto nested_scale(const std::vector<T>& a, double s) {
    std::vector<std::decay_t<decltype(nested_scale(a[0], s))>> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        out.push_back(nested_scale(a[i], s));
    }
    return out;
}

template<class T>
auto nested_div(const std::vector<T>& a, double s) {
    std::vector<std::decay_t<decltype(nested_div(a[0], s))>> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        out.push_back(nested_div(a[i], s));
    }
    return out;
}


template <bool include_covariance, class R, class T>
    requires requires(T const& v) {
        { prod_XY<include_covariance>(v, v) } -> std::convertible_to<R>;
    }
auto sum_series_sqr_X(std::vector<std::vector<R>>&& sum, const std::vector<T>& x,
                      std::size_t max_lag) {
    if (sum.empty()) {
        sum.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            const std::size_t end = std::min(x.size(), i + max_lag);
            sum[i].reserve(end - i);
            for (std::size_t j = i; j < end; ++j) {
                sum[i].push_back(prod_XY<include_covariance>(x[i], x[j]));
            }
        }
    } else {
        for (std::size_t i = 0; i < x.size(); ++i) {
            const std::size_t end = std::min(x.size(), i + max_lag);
            for (std::size_t j = i; j < end; ++j) {
                sum[i][j - i] = nested_add(sum[i][j - i], prod_XY<include_covariance>(x[i], x[j]));
            }
        }
    }
    return sum;
}

template <bool include_covariance, typename R, class V, class G>
    requires requires(G&& g, V const& v) {
        {
            prod_XY<include_covariance>(std::forward<G>(g)(v), std::forward<G>(g)(v))
        } -> std::convertible_to<R>;
    }
auto sum_series_sqr_X(std::vector<std::vector<R>>&& sum, const std::vector<V>& x,
                      std::size_t max_lag, G&& g) {
    if (sum.empty()) {
        sum.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            const std::size_t end = std::min(x.size(), i + max_lag);
            sum[i].reserve(end - i);
            for (std::size_t j = i; j < end; ++j) {
                sum[i].push_back(
                    prod_XY<include_covariance>(std::forward<G>(g)(x[i]), std::forward<G>(g)(x[j])));
            }
        }
    } else {
        for (std::size_t i = 0; i < x.size(); ++i) {
            const std::size_t end = std::min(x.size(), i + max_lag);
            for (std::size_t j = i; j < end; ++j) {
                sum[i][j - i] = nested_add(
                    sum[i][j - i],
                    prod_XY<include_covariance>(std::forward<G>(g)(x[i]), std::forward<G>(g)(x[j])));
            }
        }
    }
    return sum;
}


template<typename T, class V, class G>
requires requires (G&& g, V const& v) { {std::forward<G>(g)(v)} -> std::convertible_to<T>; }
auto sum_series_X(std::vector<T>&& sum,const std::vector<V>& x, G&& g) {
    if (sum.empty()){
        sum.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            sum[i] = std::forward<G>(g)(x[i]);
        }
    } else {
        for (std::size_t i = 0; i < x.size(); ++i) {
            sum[i] = nested_add(sum[i], std::forward<G>(g)(x[i]));
        }
    }
    return sum;
}





template <class Va>
struct count : public var::Constant<count<Va>, std::size_t> {
    count() : var::Constant<count, std::size_t>{0} {
    }
    count(std::size_t n) : var::Constant<count<Va>, std::size_t>{n} {
    }
};


template <class Va>
struct series_count : public var::Constant<series_count<Va>, std::vector<std::size_t>> {
    series_count() : var::Constant<series_count<Va>, std::vector<std::size_t>>{std::vector<std::size_t>{}} {
    }
     series_count(std::vector<std::size_t>&& x) : var::Constant<series_count<Va>, std::vector<std::size_t>>{std::move(x)} {
    }
     series_count(std::vector<std::size_t> const& x) : var::Constant<series_count<Va>, std::vector<std::size_t>>{x} {
    }
};




template <class Va>
struct mean : public var::Var<mean<Va>, mean_value_type_t<Va>>{
   using mean_value_type= mean_value_type_t<Va>;
    using value_type= value_type_t<Va>;

    mean() : var::Var<mean<Va>, mean_value_type>{value_type_t<Va>{}} {
    }
    mean(Va&& x) : var::Var<mean<Va>, mean_value_type>{std::move(x)()} {
    }
    mean(Va const& x) : var::Var<mean<Va>, mean_value_type>{x()} {
    }

    template<typename T>
    requires std::convertible_to<T, mean_value_type>
    mean(T&& x) : var::Var<mean<Va>, mean_value_type>{std::forward<T>(x)} {
    }
 };


template <class Va>
struct series_mean : public var::Var<series_mean<Va>, std::vector<mean_value_type_t<Va>>>{
    using mean_value_type= mean_value_type_t<Va>;
    using value_type= value_type_t<Va>;
    using series_mean_value_type= std::vector<mean_value_type_t<Va>>;

    series_mean() : var::Var<series_mean<Va>, series_mean_value_type>{std::vector<mean_value_type_t<Va>>{}} {
    }
    series_mean(std::vector<mean_value_type_t<Va>>&& x) : var::Var<series_mean<Va>, series_mean_value_type>{std::move(x)} {
    }
    series_mean(Va const& x) : var::Var<series_mean<Va>, series_mean_value_type>{std::vector<mean_value_type_t<Va>>{x()}} {
    }

    template<typename T>
    requires std::convertible_to<T, series_mean_value_type>
    series_mean(T&& x) : var::Var<series_mean<Va>, series_mean_value_type>{std::forward<T>(x)} {
    }
 };





template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct median : public var::Var<median<Va>, double>{
   using mean_type= double;

    median() : var::Var<median<Va>, mean_type>{value_type_t<Va>{}} {
    }
    median(Va&& x) : var::Var<median<Va>, mean_type>{std::move(x)()} {
    }
    median(Va const& x) : var::Var<median<Va>, mean_type>{x()} {
    }
    median(value_type_t<Va>&& x) : var::Var<median<Va>, mean_type>{std::move(x)} {
    }
    median(value_type_t<Va> const& x) : var::Var<median<Va>, mean_type>{x} {
    }
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct stddev : public var::Var<stddev<Va>, double>{
   using mean_type= double;

    stddev() : var::Var<stddev<Va>, mean_type>{value_type_t<Va>{}} {
    }
    stddev(Va&& x) : var::Var<stddev<Va>, mean_type>{std::move(x)()} {
    }
    stddev(Va const& x) : var::Var<stddev<Va>, mean_type>{x()} {
    }
    stddev(value_type_t<Va>&& x) : var::Var<stddev<Va>, mean_type>{std::move(x)} {
    }
    stddev(value_type_t<Va> const& x) : var::Var<stddev<Va>, mean_type>{x} {
    }
};

template <class Va>
struct Probits: public var::Var<Probits<Va>, std::map<double,mean_value_type_t<Va>>>{
   using mean_value_type= mean_value_type_t<Va>;

    Probits() : var::Var<Probits<Va>, std::map<double,mean_value_type>>{std::map<double,mean_value_type>{}} {
    }
    Probits(std::map<double,mean_value_type>&& x) : var::Var<Probits<Va>, std::map<double,mean_value_type>>{std::move(x)} {
    }
    Probits(std::map<double,mean_value_type>const & x) : var::Var<Probits<Va>, std::map<double,mean_value_type>>{x} {
    }
};


template <class Va>
struct Probit_025 : public var::Var<Probit_025<Va>, value_type_t<Va>>{
   using value_type= value_type_t<Va>;

    Probit_025() : var::Var<Probit_025<Va>, value_type>{value_type_t<Va>{}} {
    }
    Probit_025(Va&& x) : var::Var<Probit_025<Va>, value_type>{std::move(x)()} {
    }
    Probit_025(Va const& x) : var::Var<Probit_025<Va>, value_type>{x()} {
    }
    Probit_025(value_type&& x) : var::Var<Probit_025<Va>, value_type>{std::move(x)} {
    }
    Probit_025(value_type const& x) : var::Var<Probit_025<Va>, value_type>{x} {
    }
};

template <class Va>
struct Probit_975 : public var::Var<Probit_975<Va>, value_type_t<Va>>{
   using value_type= value_type_t<Va>;

    Probit_975() : var::Var<Probit_975<Va>, value_type>{value_type_t<Va>{}} {
    }
    Probit_975(Va&& x) : var::Var<Probit_975<Va>, value_type>{std::move(x)()} {
    }
    Probit_975(Va const& x) : var::Var<Probit_975<Va>, value_type>{x()} {
    }
    Probit_975(value_type&& x) : var::Var<Probit_975<Va>, value_type>{std::move(x)} {
    }
    Probit_975(value_type const& x) : var::Var<Probit_975<Va>, value_type>{x} {
    }
};


template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_minus_2_sigma : public var::Var<Probit_minus_2_sigma<Va>, double>{
   using mean_type= double;

    Probit_minus_2_sigma() : var::Var<Probit_minus_2_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_minus_2_sigma(Va&& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_minus_2_sigma(Va const& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{x()} {
    }
    Probit_minus_2_sigma(mean_type&& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_minus_2_sigma(mean_type const& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{x} {
    }
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_minus_1_sigma : public var::Var<Probit_minus_1_sigma<Va>, double>{
   using mean_type= double;

    Probit_minus_1_sigma() : var::Var<Probit_minus_1_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_minus_1_sigma(Va&& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_minus_1_sigma(Va const& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{x()} {
    }
    Probit_minus_1_sigma(mean_type&& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_minus_1_sigma(mean_type const& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{x} {
    }
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_plus_1_sigma : public var::Var<Probit_plus_1_sigma<Va>, double>{
   using mean_type= double;

    Probit_plus_1_sigma() : var::Var<Probit_plus_1_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_plus_1_sigma(Va&& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_plus_1_sigma(Va const& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{x()} {
    }
    Probit_plus_1_sigma(mean_type&& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_plus_1_sigma(mean_type const& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{x} {
    }
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_plus_2_sigma : public var::Var<Probit_plus_2_sigma<Va>, double>{
   using mean_type= double;

    Probit_plus_2_sigma() : var::Var<Probit_plus_2_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_plus_2_sigma(Va&& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_plus_2_sigma(Va const& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{x()} {
    }
    Probit_plus_2_sigma(mean_type&& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_plus_2_sigma(mean_type const& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{x} {
    }
};




template <class Va>
struct Sum : public var::Var<Sum<Va>, value_type_t<Va>>{
   using sum_type= value_type_t<Va>;

    Sum() : var::Var<Sum<Va>, sum_type>{value_type_t<Va>{}} {
    }
    Sum(Va&& x) : var::Var<Sum<Va>, sum_type>{std::move(x)()} {
    }
    Sum(Va const& x) : var::Var<Sum<Va>, sum_type>{x()} {
    }
    template<class T>
    requires std::convertible_to<T, sum_type>
     Sum(T&& x) : var::Var<Sum<Va>, sum_type>{std::forward<T>(x)} {
    }

 template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { {f(v)} -> std::convertible_to<value_type_t<Va>>; }
    Sum(std::vector<VectorSpace>const& x,  F&& f)
        : var::Var<Sum<Va>, sum_type>{value_type_t<Va>{}} {
        for (const auto& v : x) {
         (*this)() = (*this)() + sum_type(std::invoke(f, v));
    }
       }


};



template <class Va>
struct variance : public var::Var<variance<Va>, std::decay_t<decltype(sqr_X<false>(std::declval<value_type_t<Va>>()))>> {
    using variance_type= std::decay_t<decltype(sqr_X<false>(std::declval<value_type_t<Va>>()))>;
    variance() : var::Var<variance<Va>, variance_type>{variance_type{}} {
    }
    variance(variance_type&& x) : var::Var<variance<Va>, variance_type>{std::move(x)} {
    }
    variance(variance_type const& x) : var::Var<variance<Va>, variance_type>{x} {
    }
    variance(Va const& x) : var::Var<variance<Va>, variance_type>{sqr_X<false>(x())} {
    }
};



template <class Va>
struct log_Det : public var::Var<log_Det<Va>, double> {
    using log_Det_type = double;

    log_Det() : var::Var<log_Det<Va>, log_Det_type>{log_Det_type{}} {}
    log_Det(log_Det_type&& x) : var::Var<log_Det<Va>, log_Det_type>{std::move(x)} {}
    log_Det(log_Det_type const& x) : var::Var<log_Det<Va>, log_Det_type>{x} {}

    // Va whose stored value is directly logdet-able (bare matrix wrapper).
    log_Det(Va const& x)
        requires requires(Va const& v) {
            { logdet(v()) } -> std::convertible_to<Maybe_error<double>>;
        }
        : var::Var<log_Det<Va>, log_Det_type>{
              logdet(x()).value_or(std::numeric_limits<double>::quiet_NaN())} {}

    // Va whose stored value is an Indexed/Parameters wrapper exposing the
    // matrix via .value() (e.g. Sample_Distortion_Matrix).
    log_Det(Va const& x)
        requires requires(Va const& v) {
            { logdet(v().value()) } -> std::convertible_to<Maybe_error<double>>;
        }
        : var::Var<log_Det<Va>, log_Det_type>{
              logdet(x().value()).value_or(std::numeric_limits<double>::quiet_NaN())} {}
};


// Pearson correlation matrix derived from an SPD covariance Va. Same payload
// shape as the source — a p×p SPD wrapped with parameter metadata — but with
// unit diagonal. Parameter pairs with |R[i,j]| near 1 flag collinear /
// non-identifiable directions.
template <class Va>
struct Correlation_Of
    : public var::Var<Correlation_Of<Va>,
                      var::ParameterIndexed<SymPosDefMatrix<double>,
                                            var::Parameters_transformed>> {
    using payload_type = var::ParameterIndexed<SymPosDefMatrix<double>,
                                               var::Parameters_transformed>;

    Correlation_Of() : var::Var<Correlation_Of<Va>, payload_type>{payload_type{}} {}
    Correlation_Of(payload_type&& x)
        : var::Var<Correlation_Of<Va>, payload_type>{std::move(x)} {}
    Correlation_Of(payload_type const& x) : var::Var<Correlation_Of<Va>, payload_type>{x} {}

    // From a Va whose stored payload exposes .value() / .parameters_ptr()
    // (e.g. Fisher_Covariance, Distortion_Corrected_Covariance).
    Correlation_Of(Va const& x)
        requires requires(Va const& v) {
            { v().value() } -> std::convertible_to<SymPosDefMatrix<double>>;
            { v().parameters_ptr() };
        }
        : var::Var<Correlation_Of<Va>, payload_type>{
              payload_type(correlation_matrix(x().value()), x().parameters_ptr())} {}
};

// Spectrum of the SPD source Va as a p×1 column vector, sorted descending.
// Dynamic range (κ = λ_max / λ_min) + location of any gap reveal sloppiness
// and effective rank.
template <class Va>
struct Eigenvalue_Spectrum : public var::Var<Eigenvalue_Spectrum<Va>, Matrix<double>> {
    Eigenvalue_Spectrum()
        : var::Var<Eigenvalue_Spectrum<Va>, Matrix<double>>{Matrix<double>{}} {}
    Eigenvalue_Spectrum(Matrix<double>&& x)
        : var::Var<Eigenvalue_Spectrum<Va>, Matrix<double>>{std::move(x)} {}
    Eigenvalue_Spectrum(Matrix<double> const& x)
        : var::Var<Eigenvalue_Spectrum<Va>, Matrix<double>>{x} {}
};

// Count of eigenvalues above the retention tolerance — the number of
// identifiable directions. Null defect = p − Effective_Rank.
template <class Va>
struct Effective_Rank : public var::Var<Effective_Rank<Va>, std::size_t> {
    Effective_Rank() : var::Var<Effective_Rank<Va>, std::size_t>{0} {}
    Effective_Rank(std::size_t x) : var::Var<Effective_Rank<Va>, std::size_t>{x} {}
};

// Condition number κ = λ_max / λ_min of the source Va. Complements
// Effective_Rank: large κ without rank deficit = sloppy-but-identified; κ = ∞
// = formally degenerate. Actionable thresholds: < 10⁴ healthy, 10⁴–10⁸ sloppy,
// 10⁸–10¹² practically non-identifiable, > 10¹² formally degenerate.
template <class Va>
struct Spectrum_Condition_Number
    : public var::Var<Spectrum_Condition_Number<Va>, double> {
    Spectrum_Condition_Number()
        : var::Var<Spectrum_Condition_Number<Va>, double>{1.0} {}
    Spectrum_Condition_Number(double x)
        : var::Var<Spectrum_Condition_Number<Va>, double>{x} {}
};

// Orthogonal projector onto the non-identifiable (null) subspace of the
// source Va. Π_{ii} ∈ [0,1] is parameter i's fraction of variance living in
// the null subspace; off-diagonals reveal entanglement. Unlike individual
// null eigenvectors, the projector is rotation-invariant within the null
// subspace — unambiguous when null defect ≥ 2.
template <class Va>
struct Null_Space_Projector
    : public var::Var<Null_Space_Projector<Va>,
                      var::ParameterIndexed<SymPosDefMatrix<double>,
                                            var::Parameters_transformed>> {
    using payload_type = var::ParameterIndexed<SymPosDefMatrix<double>,
                                               var::Parameters_transformed>;

    Null_Space_Projector()
        : var::Var<Null_Space_Projector<Va>, payload_type>{payload_type{}} {}
    Null_Space_Projector(payload_type&& x)
        : var::Var<Null_Space_Projector<Va>, payload_type>{std::move(x)} {}
    Null_Space_Projector(payload_type const& x)
        : var::Var<Null_Space_Projector<Va>, payload_type>{x} {}
};

// Orthogonal projector onto the *least-identified* subspace — always
// populated when the spectrum is sloppy (κ ≥ activation threshold), not
// only when a mode is formally null. Complements Null_Space_Projector:
// the latter asks "are any directions truly null?"; this one asks "where
// is my sloppiness?". For the sloppy-not-null case (κ ≫ 1 but λ_min still
// above machine noise), Null_Space_Projector is empty while this one
// correctly localizes the worst direction. See
// `lapack::worst_subspace_projector` for the activation and gap-growth
// algorithm.
template <class Va>
struct Worst_Subspace_Projector
    : public var::Var<Worst_Subspace_Projector<Va>,
                      var::ParameterIndexed<SymPosDefMatrix<double>,
                                            var::Parameters_transformed>> {
    using payload_type = var::ParameterIndexed<SymPosDefMatrix<double>,
                                               var::Parameters_transformed>;

    Worst_Subspace_Projector()
        : var::Var<Worst_Subspace_Projector<Va>, payload_type>{payload_type{}} {}
    Worst_Subspace_Projector(payload_type&& x)
        : var::Var<Worst_Subspace_Projector<Va>, payload_type>{std::move(x)} {}
    Worst_Subspace_Projector(payload_type const& x)
        : var::Var<Worst_Subspace_Projector<Va>, payload_type>{x} {}
};


template <class Va>
struct covariance : public var::Var<covariance<Va>, std::decay_t<decltype(sqr_X<true>(std::declval<value_type_t<Va>>()))>> {
    using covariance_type= std::decay_t<decltype(sqr_X<true>(std::declval<value_type_t<Va>>()))>;
    covariance() : var::Var<covariance<Va>, covariance_type>{covariance_type{}} {
    }
    covariance(covariance_type&& x) : var::Var<covariance<Va>, covariance_type>{std::move(x)} {
    }
    covariance(covariance_type const& x) : var::Var<covariance<Va>, covariance_type>{x} {
    }

};

template <class Va>
struct series_covariance_kernel :
public var::Var<series_covariance_kernel<Va>, std::decay_t<decltype(sqr_X<true>(
    std::vector<value_type_t<Va>>{}, std::size_t{}))>> {
     using variance_type= decltype(prod_XY<true>(std::declval<value_type_t<Va>>(),std::declval<value_type_t<Va>>()));
    using series_covariance_kernel_type= std::decay_t<decltype(sqr_X<true>(std::vector<value_type_t<Va>>{}, std::size_t{}))>;
    using kernel_type = series_covariance_kernel_type;

    series_covariance_kernel() : var::Var<series_covariance_kernel<Va>, series_covariance_kernel_type>{series_covariance_kernel_type{}} {
    }
    series_covariance_kernel(series_covariance_kernel_type&& x) : 
    var::Var<series_covariance_kernel<Va>, series_covariance_kernel_type>{std::move(x)} {
    }
    series_covariance_kernel(series_covariance_kernel_type const& x) :
     var::Var<series_covariance_kernel<Va>, series_covariance_kernel_type>{x} {
    }

};


template <class Va>
struct series_variance_kernel :
public var::Var<series_variance_kernel<Va>, std::decay_t<decltype(sqr_X<false>(
    std::vector<value_type_t<Va>>{}, std::size_t{}))>> {
     using variance_type= decltype(prod_XY<false>(std::declval<value_type_t<Va>>(),std::declval<value_type_t<Va>>()));
    using series_variance_kernel_type= std::decay_t<decltype(sqr_X<false>(std::vector<value_type_t<Va>>{}, std::size_t{}))>;
    using kernel_type = series_variance_kernel_type;

    series_variance_kernel() : var::Var<series_variance_kernel<Va>, series_variance_kernel_type>{series_variance_kernel_type{}} {
    }
    series_variance_kernel(series_variance_kernel_type&& x) : var::Var<series_variance_kernel<Va>, series_variance_kernel_type>{std::move(x)} {
    }
    series_variance_kernel(series_variance_kernel_type const& x) : var::Var<series_variance_kernel<Va>, series_variance_kernel_type>{x} {
    }

};



template <class Va, class Vb>
struct cross_covariance : public var::Var<cross_covariance<Va, Vb>, std::decay_t<decltype(prod_XY<true>(std::declval<value_type_t<Va>>(), std::declval<value_type_t<Vb>>()))>> {
    using cross_covariance_type= std::decay_t<decltype(prod_XY<true>(std::declval<value_type_t<Va>>(), std::declval<value_type_t<Vb>>()))>;
    cross_covariance() : var::Var<cross_covariance<Va, Vb>, cross_covariance_type>{cross_covariance_type{}} {
    }
    cross_covariance(cross_covariance_type&& x) : var::Var<cross_covariance<Va, Vb>, cross_covariance_type>{std::move(x)} {
    }
    cross_covariance(cross_covariance_type const& x) : var::Var<cross_covariance<Va, Vb>, cross_covariance_type>{x} {
    }

};






template <bool include_covariance>
struct variance_kind;

template <>
struct variance_kind<true> {
    template <class Va>
    using type = covariance<Va>;
};

template <>
struct variance_kind<false  > {
    template <class Va>
    using type = variance<Va>;
};

template<class Va,bool include_covariance>
using variance_t= typename variance_kind<include_covariance>::template type<Va>;


template <bool include_covariance>
struct variance_kernel_kind;

template <>
struct variance_kernel_kind<true> {
    template <class Va>
    using type = series_covariance_kernel<Va>;
};

template <>
struct variance_kernel_kind<false  > {
    template <class Va>
    using type = series_variance_kernel<Va>;
};

template<class Va,bool include_covariance>
using variance_kernel_t= typename variance_kernel_kind<include_covariance>::template type<Va>;

template <class Va>
using series_kernel_entry_t = std::decay_t<decltype(
    prod_XY<true>(std::declval<value_type_t<Va>>(), std::declval<value_type_t<Va>>()))>;

template <class Va>
struct series_diagonal_covariance_blocks
    : public var::Var<series_diagonal_covariance_blocks<Va>, std::vector<series_kernel_entry_t<Va>>> {
    using entry_type = series_kernel_entry_t<Va>;
    using diagonal_type = std::vector<entry_type>;

    series_diagonal_covariance_blocks()
        : var::Var<series_diagonal_covariance_blocks<Va>, diagonal_type>{diagonal_type{}} {}
    series_diagonal_covariance_blocks(diagonal_type&& x)
        : var::Var<series_diagonal_covariance_blocks<Va>, diagonal_type>{std::move(x)} {}
    series_diagonal_covariance_blocks(diagonal_type const& x)
        : var::Var<series_diagonal_covariance_blocks<Va>, diagonal_type>{x} {}
};

// Per-sample variance entries (just the diagonal elements, not the full block
// per sample). For scalar V this is the same data as
// series_diagonal_covariance_blocks; for vector V it stores only the p diagonal
// entries per sample, shrinking the storage by a factor of p compared to full
// per-sample blocks. Used by local_var-level reports where per-parameter
// correlation cross-terms per sample aren't required.
template <class Va>
struct series_diagonal_variance
    : public var::Var<series_diagonal_variance<Va>, std::vector<mean_value_type_t<Va>>> {
    using entry_type = mean_value_type_t<Va>;
    using diagonal_type = std::vector<entry_type>;

    series_diagonal_variance()
        : var::Var<series_diagonal_variance<Va>, diagonal_type>{diagonal_type{}} {}
    series_diagonal_variance(diagonal_type&& x)
        : var::Var<series_diagonal_variance<Va>, diagonal_type>{std::move(x)} {}
    series_diagonal_variance(diagonal_type const& x)
        : var::Var<series_diagonal_variance<Va>, diagonal_type>{x} {}
};

template <class Va>
struct series_cross_correlation_kernel
    : public var::Var<
          series_cross_correlation_kernel<Va>,
          std::decay_t<decltype(sqr_X<true>(std::vector<value_type_t<Va>>{}, std::size_t{}))>> {
    using entry_type = series_kernel_entry_t<Va>;
    using kernel_type = std::decay_t<decltype(sqr_X<true>(std::vector<value_type_t<Va>>{}, std::size_t{}))>;

    series_cross_correlation_kernel()
        : var::Var<series_cross_correlation_kernel<Va>, kernel_type>{kernel_type{}} {}
    series_cross_correlation_kernel(kernel_type&& x)
        : var::Var<series_cross_correlation_kernel<Va>, kernel_type>{std::move(x)} {}
    series_cross_correlation_kernel(kernel_type const& x)
        : var::Var<series_cross_correlation_kernel<Va>, kernel_type>{x} {}
};

// Per-sample forward sum of the normalized cross-correlation kernel:
//     cross_correlation_lag_forward[i] = sum_{offset=0..max_lag-1} kernel[i][offset]
// One matrix per sample, recording cumulative forward-lag correlation from that
// starting point. "Forward" because the kernel is stored only for non-negative
// lags (backward lags would require samples before i and are not representable
// at a fixed starting index).
template <class Va>
struct cross_correlation_lag_forward
    : public var::Var<cross_correlation_lag_forward<Va>,
                      std::vector<series_kernel_entry_t<Va>>> {
    using entry_type = series_kernel_entry_t<Va>;
    using forward_type = std::vector<entry_type>;

    cross_correlation_lag_forward()
        : var::Var<cross_correlation_lag_forward<Va>, forward_type>{forward_type{}} {}
    cross_correlation_lag_forward(forward_type&& x)
        : var::Var<cross_correlation_lag_forward<Va>, forward_type>{std::move(x)} {}
    cross_correlation_lag_forward(forward_type const& x)
        : var::Var<cross_correlation_lag_forward<Va>, forward_type>{x} {}
};

// Trace-global integral correlation lag:
//     integral_correlation_lag = C(0)^{-1/2} * Omega * C(0)^{-1/2}
// where
//     C(0) = sum_i kernel_raw[i][0]           (per-sample variance summed across samples)
//     Omega = sum_i sum_{offset} kernel_raw[i][offset]   (summed unnormalized kernel)
// Scalar for scalar Va (V / V²); matrix for vector Va (the truncated-kernel
// analog of the Correlation_Distortion_Matrix). Kish/Sokal integral
// autocorrelation time in a unified form.
template <class Va>
struct integral_correlation_lag
    : public var::Var<integral_correlation_lag<Va>, series_kernel_entry_t<Va>> {
    using entry_type = series_kernel_entry_t<Va>;

    integral_correlation_lag()
        : var::Var<integral_correlation_lag<Va>, entry_type>{entry_type{}} {}
    integral_correlation_lag(entry_type&& x)
        : var::Var<integral_correlation_lag<Va>, entry_type>{std::move(x)} {}
    integral_correlation_lag(entry_type const& x)
        : var::Var<integral_correlation_lag<Va>, entry_type>{x} {}
};

// Report family — four levels of per-variable series analysis, each a strict
// superset of the previous except for the local_var / local_cov pair which
// differ in storage of the per-sample diagonal (variance only vs full
// covariance block).
//
//   Report_integral<V>   : count + integral_correlation_lag
//   Report_local_var<V>  : integral + series_mean + series_diagonal_variance
//                          + cross_correlation_lag_forward
//   Report_local_cov<V>  : integral + series_mean + series_diagonal_covariance_blocks
//                          + cross_correlation_lag_forward
//   Report_cross<V>      : local_cov + series_cross_correlation_kernel
//
// For scalar V, local_var == local_cov (same storage).

template <class Id>
class Report_integral
    : public var::Var<Report_integral<Id>,
                      var::Vector_Space<count<Id>, integral_correlation_lag<Id>>> {
   public:
    using base_space_type = var::Vector_Space<count<Id>, integral_correlation_lag<Id>>;
    using base_type = var::Var<Report_integral<Id>, base_space_type>;

    Report_integral()
        : base_type{base_space_type{count<Id>{}, integral_correlation_lag<Id>{}}} {}
    Report_integral(base_space_type&& x) : base_type{std::move(x)} {}
    Report_integral(base_space_type const& x) : base_type{x} {}
    Report_integral(count<Id> n, integral_correlation_lag<Id> g)
        : base_type{base_space_type{std::move(n), std::move(g)}} {}

    constexpr auto& operator[](var::Var<Id>) { return *this; }
    constexpr auto& operator[](var::Var<Id>) const { return *this; }
};

template <class Id>
class Report_local_var
    : public var::Var<Report_local_var<Id>,
                      var::Vector_Space<count<Id>, series_mean<Id>,
                                        series_diagonal_variance<Id>,
                                        cross_correlation_lag_forward<Id>,
                                        integral_correlation_lag<Id>>> {
   public:
    using base_space_type = var::Vector_Space<count<Id>, series_mean<Id>,
                                              series_diagonal_variance<Id>,
                                              cross_correlation_lag_forward<Id>,
                                              integral_correlation_lag<Id>>;
    using base_type = var::Var<Report_local_var<Id>, base_space_type>;

    Report_local_var()
        : base_type{base_space_type{count<Id>{}, series_mean<Id>{},
                                    series_diagonal_variance<Id>{},
                                    cross_correlation_lag_forward<Id>{},
                                    integral_correlation_lag<Id>{}}} {}
    Report_local_var(base_space_type&& x) : base_type{std::move(x)} {}
    Report_local_var(base_space_type const& x) : base_type{x} {}
    Report_local_var(count<Id> n, series_mean<Id> m, series_diagonal_variance<Id> d,
                     cross_correlation_lag_forward<Id> f, integral_correlation_lag<Id> g)
        : base_type{base_space_type{std::move(n), std::move(m), std::move(d), std::move(f),
                                    std::move(g)}} {}

    constexpr auto& operator[](var::Var<Id>) { return *this; }
    constexpr auto& operator[](var::Var<Id>) const { return *this; }
};

template <class Id>
class Report_local_cov
    : public var::Var<Report_local_cov<Id>,
                      var::Vector_Space<count<Id>, series_mean<Id>,
                                        series_diagonal_covariance_blocks<Id>,
                                        cross_correlation_lag_forward<Id>,
                                        integral_correlation_lag<Id>>> {
   public:
    using base_space_type = var::Vector_Space<count<Id>, series_mean<Id>,
                                              series_diagonal_covariance_blocks<Id>,
                                              cross_correlation_lag_forward<Id>,
                                              integral_correlation_lag<Id>>;
    using base_type = var::Var<Report_local_cov<Id>, base_space_type>;

    Report_local_cov()
        : base_type{base_space_type{count<Id>{}, series_mean<Id>{},
                                    series_diagonal_covariance_blocks<Id>{},
                                    cross_correlation_lag_forward<Id>{},
                                    integral_correlation_lag<Id>{}}} {}
    Report_local_cov(base_space_type&& x) : base_type{std::move(x)} {}
    Report_local_cov(base_space_type const& x) : base_type{x} {}
    Report_local_cov(count<Id> n, series_mean<Id> m, series_diagonal_covariance_blocks<Id> d,
                     cross_correlation_lag_forward<Id> f, integral_correlation_lag<Id> g)
        : base_type{base_space_type{std::move(n), std::move(m), std::move(d), std::move(f),
                                    std::move(g)}} {}

    constexpr auto& operator[](var::Var<Id>) { return *this; }
    constexpr auto& operator[](var::Var<Id>) const { return *this; }
};

template <class Id>
class Report_cross
    : public var::Var<Report_cross<Id>,
                      var::Vector_Space<count<Id>, series_mean<Id>,
                                        series_diagonal_covariance_blocks<Id>,
                                        series_cross_correlation_kernel<Id>,
                                        cross_correlation_lag_forward<Id>,
                                        integral_correlation_lag<Id>>> {
   public:
    using base_space_type = var::Vector_Space<count<Id>, series_mean<Id>,
                                              series_diagonal_covariance_blocks<Id>,
                                              series_cross_correlation_kernel<Id>,
                                              cross_correlation_lag_forward<Id>,
                                              integral_correlation_lag<Id>>;
    using base_type = var::Var<Report_cross<Id>, base_space_type>;

    Report_cross()
        : base_type{base_space_type{count<Id>{}, series_mean<Id>{},
                                    series_diagonal_covariance_blocks<Id>{},
                                    series_cross_correlation_kernel<Id>{},
                                    cross_correlation_lag_forward<Id>{},
                                    integral_correlation_lag<Id>{}}} {}
    Report_cross(base_space_type&& x) : base_type{std::move(x)} {}
    Report_cross(base_space_type const& x) : base_type{x} {}
    Report_cross(count<Id> n, series_mean<Id> m, series_diagonal_covariance_blocks<Id> d,
                 series_cross_correlation_kernel<Id> r, cross_correlation_lag_forward<Id> f,
                 integral_correlation_lag<Id> g)
        : base_type{base_space_type{std::move(n), std::move(m), std::move(d), std::move(r),
                                    std::move(f), std::move(g)}} {}

    constexpr auto& operator[](var::Var<Id>) { return *this; }
    constexpr auto& operator[](var::Var<Id>) const { return *this; }
};

// Backward-compatibility alias. Existing write_csv specializations and traits
// that reference Series_Moment_report continue to work. Prefer Report_cross in
// new code; this alias can be removed once all callers are migrated.
template <class Id>
using Series_Moment_report = Report_cross<Id>;


template <class Id, bool include_covariance=false>
class Moment_statistics
    : public var::Var<Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, mean<Id>,
    variance_t<Id,include_covariance>>>
     {
   public:
    using variance_type = variance_t<Id, include_covariance>;
    using base_type =var::Var<Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, mean<Id>,
    variance_type>> ;

    using mean_type=typename mean<Id>::mean_value_type;

    Moment_statistics()
        : base_type{
              var::Vector_Space{count<Id>{0}, mean<Id>(Id{}()),
                                variance_type{}}} {
    }
    Moment_statistics(Id x)
        : base_type{var::Vector_Space{count<Id>{1}, mean<Id>(x()), variance_type(sqr_X<include_covariance>(x()-x()))}} {
    }
    Moment_statistics(value_type_t<Id> x)
        : base_type{var::Vector_Space{count<Id>{1}, mean<Id>(x), variance_type(sqr_X<include_covariance>(x-x))}} {
    }



    Moment_statistics(std::vector<Id>const& x)
        : base_type{
              var::Vector_Space{count<Id>{x.size()}, mean<Id>(Id{}()),
                                variance_type{}}} {
        if (x.empty()) {
            return;
        }
        auto sum=Id{}();
        auto sum2=variance_type{}();
        for (const auto& xi : x) {
            sum = sum + xi();
            sum2 = sum2 + sqr_X<include_covariance>(xi());
        }
        get<mean<Id>>((*this)())()= sum/x.size();
        auto var= x.size()>1 ? (sum2 - x.size() * sqr_X<include_covariance>(get<mean<Id>>((*this)())()))/(x.size()-1)
                             : variance_type{}();
        get<variance_type>((*this)())()= var;
    }

    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F&& f) { {f(v)} -> std::convertible_to<value_type_t<Id>>; }
    Moment_statistics(std::vector<VectorSpace>const& x, const std::vector<std::size_t>& indices, F&& f)
        : base_type{
              var::Vector_Space{count<Id>{indices.size()}, mean<Id>(Id{}()),
                                variance_type{}}} {
        if (indices.empty()) {
            return;
        }
        auto sum=Id{}();
        auto sum2=variance_type{}();
        for (const auto& i : indices) {
            auto xi = value_type_t<Id>(std::forward<F>(f)(x[i]));
            sum = sum + xi;
            sum2 = sum2 + sqr_X<include_covariance>(xi);
        }
        get<mean<Id>>((*this)())()= sum/indices.size();
        auto var= indices.size()>1 ? (sum2 - indices.size() * sqr_X<include_covariance>(get<mean<Id>>((*this)())()))/(indices.size()-1)
                                   : variance_type{}();
        get<variance_type>((*this)())()= var;
    }


    Moment_statistics(var::Vector_Space<count<Id>, mean<Id>, variance_type>&& vs) : base_type{std::move(vs)} {
    }

    Moment_statistics(var::Vector_Space<count<Id>, mean<Id>, variance_type>const& vs) : base_type{vs} {
    }


    Moment_statistics(count<Id> n,mean<Id> x, variance_type v) : base_type{var::Vector_Space{n, x, v}} {
    }

    friend Moment_statistics operator&(const Moment_statistics& one,
                                       const Moment_statistics& other) {
        double n0 = get<count<Id>>((one)())();
        double n1 = get<count<Id>>((other)())();

        auto m0 = get<mean<Id>>((one)())();
        auto m1 = get<mean<Id>>((other)())();

        auto v0 = get<variance_type>((one)())();
        auto v1 = get<variance_type>((other)())();


        double n = n0 + n1;

        auto m = m0 + n1 / n * (m1 - m0);

        auto malt = (m0 * n0 + m1 * n1) / n;

        //  auto v = n>1?v0 + (n1-1)/(n-1)* (v1-v0) + n0/(n-1) * sqr_X(m0-m) + n1/(n-1) * sqr_X(m1-m):0;
        auto v =
            n > 1
                ? ((n0 - 1) * v0 + n0 * sqr_X<include_covariance>(m0) + (n1 - 1) * v1 +
                n1 * sqr_X<include_covariance>(m1) - n * sqr_X<include_covariance>(m)) /
                      (n - 1)
                : sqr_X<include_covariance>(malt-malt);


        return Moment_statistics{count<Id>(n), mean<Id>(Id(m)), variance_type(v)};
    }




    Moment_statistics& operator&=(const Moment_statistics& other) {
        *this = *this & other;
        return *this;
    }

    Moment_statistics& operator&=(const Id& x) {
        auto n0 = get<count<Id>>((*this)())();

        auto m0 = get<mean<Id>>((*this)())();

        auto v0 = get<variance_type>((*this)())();
        auto n = n0 + 1;

        auto m = m0 + (x() - m0) / n;
        auto v = n > 1 ? v0 + (sqr_X<include_covariance>(x() - m) - v0) / n0 + sqr_X<include_covariance>(m0 - m) : Id{}();

        get<count<Id>>((*this)())() = n;
        get<variance_type>((*this)())()= std::move(v);
        get<mean<Id>>((*this)())() = std::move(m);
        return *this;
    }
    void reset() {
        get<count<Id>>((*this)())() = 0;
        get<mean<Id>>((*this)())() = Id{}();
        get<variance_type>((*this)())() = sqr_X<include_covariance>(Id{}());
    }

    friend Moment_statistics operator+(const Moment_statistics& one,
                                       const Moment_statistics& other) {
        auto n0 = get<count<Id> >((one)())();
        auto m0 = get<mean<Id>>((one)())();
        auto v0 = get<variance_type>((one)())();
        auto n1 = get<count<Id>>((other)())();
        auto m1 = get<mean<Id>>((other)())();
        auto v1 = get<variance_type>((other)())();
        return Moment_statistics(count<Id>(std::min(n0, n1)), mean<Id>(Id(m0 + m1)), variance_type(Id(v0 + v1)));

    }

    friend Moment_statistics operator*(double a, Moment_statistics const& one) {
        auto n0 = get<count<Id>>((one)())();
        auto m0 = get<mean<Id>>((one)())();
        auto v0 = get<variance_type>((one)())();
        return Moment_statistics(count<Id>(n0),mean<Id>(Id(m0 * a)), variance_type(Id(v0 * a * a)));
    }

    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};





template <class... Ids, bool include_covariance>
struct Moment_statistics<var::Vector_Space<Ids...>, include_covariance>
    : var::Vector_Space<Moment_statistics<Ids, include_covariance>...> {
    using base_type = var::Vector_Space<Moment_statistics<Ids, include_covariance>...>;

    Moment_statistics() : base_type{Moment_statistics<Ids, include_covariance>()...} {}
    explicit Moment_statistics(const var::Vector_Space<Ids...>& sample)
        : base_type{Moment_statistics<Ids, include_covariance>(get<Ids>(sample))...} {}

    template <class... Ms>
        requires(sizeof...(Ms) == sizeof...(Ids))
    explicit Moment_statistics(Ms&&... ms) : base_type{std::forward<Ms>(ms)...} {}

    void reset() { (get<Moment_statistics<Ids, include_covariance>>(*this).reset(), ...); }

    friend auto operator&(const Moment_statistics& a, const Moment_statistics& b) {
        return Moment_statistics((get<Moment_statistics<Ids, include_covariance>>(a) & get<Moment_statistics<Ids, include_covariance>>(b))...);
    }

    Moment_statistics& operator&=(const Moment_statistics& other) {
        ((get<Moment_statistics<Ids, include_covariance>>(*this) &= get<Moment_statistics<Ids, include_covariance>>(other)), ...);
        return *this;
    }

    Moment_statistics& operator&=(const var::Vector_Space<Ids...>& sample) {
        ((get<Moment_statistics<Ids, include_covariance>>(*this) &= get<Ids>(sample)), ...);
        return *this;
    }

    friend auto operator+(const Moment_statistics& a, const Moment_statistics& b) {
        return Moment_statistics((get<Moment_statistics<Ids, include_covariance>>(a) + get<Moment_statistics<Ids, include_covariance>>(b))...);
    }

    friend auto operator*(double k, const Moment_statistics& a) {
        return Moment_statistics((k * get<Moment_statistics<Ids, include_covariance>>(a))...);
    }

    template <class Id>
    friend auto& get(Moment_statistics& x) {
        return get<Moment_statistics<Id, include_covariance>>(static_cast<base_type&>(x));
    }
    template <class Id>
    friend auto const& get(Moment_statistics const& x) {
        return get<Moment_statistics<Id, include_covariance>>(static_cast<base_type const&>(x));
    }
};



template <class Id, bool include_covariance=false>
class Series_Moment_statistics
    : public var::Var<Series_Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, series_mean<Id>,
    variance_kernel_t<Id,include_covariance>>>
     {
   public:
    using series_variance_type = variance_kernel_t<Id, include_covariance>;
    using base_type =var::Var<Series_Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, series_mean<Id>,
    variance_kernel_t<Id,include_covariance>>>;
    using kernel_value_type = typename series_variance_type::kernel_type;

    using base_type::operator();
    using series_mean_type=typename series_mean<Id>::series_mean_value_type;

    Series_Moment_statistics()
        : base_type{
              var::Vector_Space{count<Id>{ 0UL},series_mean<Id>{},
                                series_variance_type{}}} {
    }
    Series_Moment_statistics(std::vector<value_type_t<Id>> x, std::size_t max_lag)
        : base_type{var::Vector_Space{
              count<Id>{1}, series_mean<Id>(x),
              series_variance_type(sqr_X<include_covariance>(nested_sub(x, x), max_lag))}} {
    }



    Series_Moment_statistics(std::vector<std::vector<value_type_t<Id>>> const& x,
                             std::size_t max_lag)
        : base_type{
              var::Vector_Space{count<Id>{x.size()}, series_mean<Id>{},
                                series_variance_type{}}} {
        if (x.empty()) {
            return;
        }
        auto sum = x[0];
        auto sum2 = sqr_X<include_covariance>(x[0], max_lag);
        for (std::size_t i = 1; i < x.size(); ++i) {
            sum = nested_add(sum, x[i]);
            sum2 = nested_add(sum2, sqr_X<include_covariance>(x[i], max_lag));
        }
        get<series_mean<Id>>((*this)())() = nested_div(sum, x.size());
        auto var = x.size() > 1
                       ? nested_div(nested_sub(sum2,
                                               nested_scale(sqr_X<include_covariance>(
                                                                get<series_mean<Id>>((*this)())(),
                                                                max_lag),
                                                            x.size())),
                                    x.size() - 1)
                       : kernel_value_type{};
        get<series_variance_type>((*this)())() = var;
    }

    template<class VectorSpace, class F, class G>
requires
    is_of_this_template_type_v<
        std::decay_t<std::invoke_result_t<F&, const VectorSpace&>>,
        std::vector
    >
    &&
    std::convertible_to<
        decltype(std::declval<G&>()(
            std::declval<typename std::decay_t<std::invoke_result_t<F&, const VectorSpace&>>::value_type const&>()
        )),
        value_type_t<Id>
    >
    Series_Moment_statistics(std::vector<VectorSpace>const& x, std::size_t max_lag,const std::vector<std::size_t>& indices, F&& f, G&& g)
        : base_type{
              var::Vector_Space{count<Id>{indices.size()}, series_mean<Id>{},
                                variance_kernel_t<Id,include_covariance>{}}} {
        if (indices.empty()) {
            return;
        }
        auto sum=std::vector<value_type_t<Id>>{};
        auto sum2=kernel_value_type{};
        for (auto idx : indices) {
            auto xi =std::forward<F>(f)(x[idx]);
            sum = sum_series_X(std::move(sum),xi, std::forward<G>(g));
            sum2 = sum_series_sqr_X<include_covariance>(std::move(sum2), xi, max_lag,
                                                       std::forward<G>(g));
        }
        get<series_mean<Id>>((*this)())() = nested_div(sum, indices.size());
        auto var= indices.size()>1 ? nested_div(
                                    nested_sub(sum2, nested_scale(sqr_X<include_covariance>(get<series_mean<Id>>((*this)())(), max_lag), indices.size())),
                                    indices.size()-1)
                             : kernel_value_type{};
        get<series_variance_type>((*this)())()= var;
    }


    Series_Moment_statistics(var::Vector_Space<count<Id>, series_mean<Id>, series_variance_type>&& vs) : base_type{std::move(vs)} {
    }

    Series_Moment_statistics(var::Vector_Space<count<Id>, series_mean<Id>, series_variance_type>const& vs) : base_type{vs} {
    }


    Series_Moment_statistics(count<Id> n,series_mean<Id> x, series_variance_type v) : base_type{var::Vector_Space{n, x, v}} {
    }

    friend Series_Moment_statistics operator&(const Series_Moment_statistics& one,
                                       const Series_Moment_statistics& other) {
        double n0 = get<count<Id>>((one)())();
        double n1 = get<count<Id>>((other)())();

        auto m0 = get<series_mean<Id>>((one)())();
        auto m1 = get<series_mean<Id>>((other)())();

        auto v0 = get<series_variance_type>((one)())();
        auto v1 = get<series_variance_type>((other)())();


        double n = n0 + n1;

        auto m = nested_add(m0, nested_scale(nested_sub(m1, m0), n1 / n));

        auto malt = nested_div(nested_add(nested_scale(m0, n0), nested_scale(m1, n1)), n);

        //  auto v = n>1?v0 + (n1-1)/(n-1)* (v1-v0) + n0/(n-1) * sqr_X(m0-m) + n1/(n-1) * sqr_X(m1-m):0;
        auto v =
            n > 1
                ? nested_div(
                    nested_sub(
                        nested_add(
                            nested_add(
                                nested_scale(v0, n0 - 1),
                                nested_scale(sqr_X<include_covariance>(m0), n0)),
                            nested_add(
                                nested_scale(v1, n1 - 1),
                                nested_scale(sqr_X<include_covariance>(m1), n1))),
                        nested_scale(sqr_X<include_covariance>(m), n)),
                    n - 1)
                : sqr_X<include_covariance>(malt-malt);


        return Series_Moment_statistics{count<Id>(n), series_mean<Id>(m), series_variance_type(v)};
    }




    Series_Moment_statistics& operator&=(const Series_Moment_statistics& other) {
        *this = *this & other;
        return *this;
    }

    Series_Moment_statistics& accumulate(const std::vector<value_type_t<Id>>& x,
                                         std::size_t max_lag) {
        *this &= Series_Moment_statistics(x, max_lag);
        return *this;
    }
    void reset() {
        get<count<Id>>((*this)())() = 0;
        get<series_mean<Id>>((*this)())() = {};
        get<series_variance_type>((*this)())() = {};
    }

    friend Series_Moment_statistics operator+(const Series_Moment_statistics& one,
                                       const Series_Moment_statistics& other) {
        auto n0 = get<count<Id> >((one)())();
        auto m0 = get<series_mean<Id>>((one)())();
        auto v0 = get<series_variance_type>((one)())();
        auto n1 = get<count<Id>>((other)())();
        auto m1 = get<series_mean<Id>>((other)())();
        auto v1 = get<series_variance_type>((other)())();
        return Series_Moment_statistics(count<Id>(std::min(n0, n1)), series_mean<Id>(nested_add(m0, m1)), series_variance_type(nested_add(v0, v1)));

    }

    friend Series_Moment_statistics operator*(double a, Series_Moment_statistics const& one) {
        auto n0 = get<count<Id>>((one)())();
        auto m0 = get<series_mean<Id>>((one)())();
        auto v0 = get<series_variance_type>((one)())();
        return Series_Moment_statistics(count<Id>(n0), series_mean<Id>(nested_scale(m0, a)), series_variance_type(nested_scale(v0, a * a)));
    }

    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};







inline Maybe_error<double> normalize_series_cross_correlation_block(
    double sigma_i, double block, double sigma_j, double /*rtol*/ = 1e-10,
    double /*atol*/ = 0.0, const std::string& /*left_name*/ = "left PSD block",
    const std::string& /*middle_name*/ = "cross block",
    const std::string& /*right_name*/ = "right PSD block") {
    if (!std::isfinite(sigma_i) || !std::isfinite(sigma_j) || !std::isfinite(block))
        return error_message("normalize_series_cross_correlation_block: non-finite scalar input");
    if (sigma_i <= 0.0 || sigma_j <= 0.0)
        return 0.0;
    return block / std::sqrt(sigma_i * sigma_j);
}

inline Maybe_error<Matrix<double>> normalize_series_cross_correlation_block(
    const Matrix<double>& sigma_i, const Matrix<double>& block, const Matrix<double>& sigma_j,
    double rtol = 1e-10, double atol = 0.0, const std::string& left_name = "left PSD block",
    const std::string& middle_name = "cross block",
    const std::string& right_name = "right PSD block") {
    return psd_whitened_cross_matrix_subspace(sigma_i, block, sigma_j, left_name, middle_name,
                                              right_name, rtol, atol);
}

template <class ValueT, class Params>
    requires var::ParameterMetadataLike<Params>
inline Maybe_error<var::ParameterIndexed<Matrix<double>, Params>>
normalize_series_cross_correlation_block(
    const var::ParameterIndexed<ValueT, Params>& sigma_i,
    const var::ParameterIndexed<ValueT, Params>& block,
    const var::ParameterIndexed<ValueT, Params>& sigma_j, double rtol = 1e-10,
    double atol = 0.0, const std::string& left_name = "left PSD block",
    const std::string& middle_name = "cross block",
    const std::string& right_name = "right PSD block") {
    return_error<var::ParameterIndexed<Matrix<double>, Params>> Error{
        "normalize_series_cross_correlation_block(ParameterIndexed): "};

    auto maybe_block = normalize_series_cross_correlation_block(
        sigma_i.value(), block.value(), sigma_j.value(), rtol, atol, left_name, middle_name,
        right_name);
    if (!maybe_block)
        return Error(maybe_block.error()());

    const Params* metadata = sigma_i.has_parameters()
                                 ? sigma_i.parameters_ptr()
                                 : (sigma_j.has_parameters() ? sigma_j.parameters_ptr()
                                                             : block.parameters_ptr());
    return var::ParameterIndexed<Matrix<double>, Params>(std::move(maybe_block.value()), metadata);
}

template <class Id>
auto series_diagonal_covariance(const Series_Moment_statistics<Id, true>& stats) {
    using diagonal_type = typename series_diagonal_covariance_blocks<Id>::diagonal_type;

    diagonal_type diagonal;
    const auto& kernel = get<variance_kernel_t<Id, true>>(stats())();
    diagonal.reserve(kernel.size());
    for (const auto& row : kernel) {
        if (row.empty())
            diagonal.emplace_back();
        else
            diagonal.push_back(row.front());
    }
    return series_diagonal_covariance_blocks<Id>(std::move(diagonal));
}

template <class Id>
Maybe_error<series_cross_correlation_kernel<Id>> series_cross_correlation(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return_error<series_cross_correlation_kernel<Id>> Error{"series_cross_correlation: "};

    using kernel_type = typename series_cross_correlation_kernel<Id>::kernel_type;
    kernel_type out;
    const auto diagonal = series_diagonal_covariance(stats);
    const auto& diagonal_blocks = diagonal();
    const auto& kernel = get<variance_kernel_t<Id, true>>(stats())();

    out.resize(kernel.size());
    if (kernel.empty())
        return Maybe_error<series_cross_correlation_kernel<Id>>(
            series_cross_correlation_kernel<Id>(std::move(out)));

    using block_entry_t = std::decay_t<decltype(diagonal_blocks.front())>;

    if constexpr (std::is_arithmetic_v<block_entry_t>) {
        // Scalar Id: amortize 1/sqrt(sigma_i) across all offsets. Halves the
        // sqrt count but the cost was always O(1) per block — minor win.
        std::vector<double> inv_sqrt_sigma(diagonal_blocks.size(), 0.0);
        for (std::size_t i = 0; i < diagonal_blocks.size(); ++i) {
            if (!std::isfinite(diagonal_blocks[i]))
                return Error("diagonal[" + std::to_string(i) + "]: non-finite");
            inv_sqrt_sigma[i] =
                (diagonal_blocks[i] > 0.0) ? 1.0 / std::sqrt(diagonal_blocks[i]) : 0.0;
        }
        for (std::size_t i = 0; i < kernel.size(); ++i) {
            out[i].reserve(kernel[i].size());
            for (std::size_t offset = 0; offset < kernel[i].size(); ++offset) {
                const std::size_t j = i + offset;
                const double v = kernel[i][offset];
                if (!std::isfinite(v))
                    return Error("kernel[" + std::to_string(i) + "][" + std::to_string(offset) +
                                 "]: non-finite");
                out[i].push_back(v * inv_sqrt_sigma[i] * inv_sqrt_sigma[j]);
            }
        }
    } else {
        // Matrix / ParameterIndexed: one eigendecomp per sample, reused across
        // all (i, i+offset) pairs for offset ∈ [0, max_lag). Dominant win when
        // max_lag is comparable to the number of parameters — typical figure-2
        // runs spend >90% of per-bootstrap-sample time in this loop's
        // eigendecompositions, which this refactor amortizes ~max_lag×.
        // Return the underlying matrix *preserving its static type* so that
        // SymPosDefMatrix inputs dispatch to the SymPosDefMatrix overload of
        // compute_psd_decomp (single retained eigendecomp) rather than the
        // dense Matrix<double> overload (which first validates via a full
        // eigendecomp, then calls the SymPosDef path — doubling the cost).
        auto get_matrix = [](const auto& x) -> decltype(auto) {
            if constexpr (requires { x.value(); })
                return (x.value());
            else
                return (x);
        };
        std::vector<lapack::PSDDecomposition> whitenings;
        whitenings.reserve(diagonal_blocks.size());
        for (std::size_t i = 0; i < diagonal_blocks.size(); ++i) {
            auto maybe_W = lapack::compute_psd_decomp(
                get_matrix(diagonal_blocks[i]),
                "series diagonal covariance[" + std::to_string(i) + "]", rtol, atol);
            if (!maybe_W)
                return Error(maybe_W.error()());
            whitenings.push_back(std::move(maybe_W.value()));
        }

        for (std::size_t i = 0; i < kernel.size(); ++i) {
            out[i].reserve(kernel[i].size());
            for (std::size_t offset = 0; offset < kernel[i].size(); ++offset) {
                const std::size_t j = i + offset;
                auto maybe_block = lapack::apply_psd_whitening(
                    whitenings[i], get_matrix(kernel[i][offset]), whitenings[j],
                    "series covariance kernel[" + std::to_string(i) + "," + std::to_string(j) +
                        "]");
                if (!maybe_block)
                    return Error(maybe_block.error()());
                if constexpr (requires { diagonal_blocks[i].parameters_ptr(); }) {
                    const auto* params = diagonal_blocks[i].has_parameters()
                                             ? diagonal_blocks[i].parameters_ptr()
                                             : (diagonal_blocks[j].has_parameters()
                                                    ? diagonal_blocks[j].parameters_ptr()
                                                    : nullptr);
                    out[i].push_back(
                        block_entry_t(std::move(maybe_block.value()), params));
                } else {
                    out[i].push_back(std::move(maybe_block.value()));
                }
            }
        }
    }

    return Maybe_error<series_cross_correlation_kernel<Id>>(
        series_cross_correlation_kernel<Id>(std::move(out)));
}

// Per-sample forward sum of the normalized cross-correlation kernel.
// For each starting sample i, sums kernel[i][offset] over offset = 0..max_lag-1.
// Result: one entry per sample, same shape as series_diagonal_covariance_blocks.
template <class Id>
cross_correlation_lag_forward<Id> series_cross_correlation_lag_forward(
    const series_cross_correlation_kernel<Id>& kernel) {
    using forward_type = typename cross_correlation_lag_forward<Id>::forward_type;
    using entry_type = typename cross_correlation_lag_forward<Id>::entry_type;
    forward_type forward;
    const auto& rows = kernel();
    forward.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.empty()) {
            forward.emplace_back();
            continue;
        }
        entry_type acc = row.front();
        for (std::size_t offset = 1; offset < row.size(); ++offset) {
            acc = acc + row[offset];
        }
        forward.push_back(std::move(acc));
    }
    return cross_correlation_lag_forward<Id>(std::move(forward));
}

// Trace-global integral correlation lag via the unified sandwich formula.
// Inputs: the unnormalized kernel (from Series_Moment_statistics) — entries are
// raw per-sample cross-products of the quantity V at (start, start+offset).
//
//     Omega = sum_i sum_{offset>=0} [ kernel[i][offset] + kernel[i][offset]^T ]
//             - sum_i kernel[i][0]    (avoid double-counting the lag-0 diagonal)
//     C(0)  = sum_i kernel[i][0]
//
// Then:
//     integral_correlation_lag = C(0)^{-1/2} * Omega * C(0)^{-1/2}
//
// For scalar V this reduces to Omega / C(0), i.e. the Kish/Sokal integral
// autocorrelation time (bidirectional: positive and negative lags counted
// symmetrically).
template <class Id>
Maybe_error<integral_correlation_lag<Id>> series_integral_correlation_lag(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10,
    double atol = 0.0) {
    return_error<integral_correlation_lag<Id>> Error{"series_integral_correlation_lag: "};

    using entry_type = typename integral_correlation_lag<Id>::entry_type;
    const auto& kernel = get<variance_kernel_t<Id, true>>(stats())();

    if (kernel.empty()) {
        return integral_correlation_lag<Id>(entry_type{});
    }

    // C(0) = sum_i kernel[i][0]   (sum of per-sample variances / lag-0 blocks)
    bool c0_initialised = false;
    entry_type c0{};
    for (std::size_t i = 0; i < kernel.size(); ++i) {
        if (kernel[i].empty()) continue;
        if (!c0_initialised) {
            c0 = kernel[i][0];
            c0_initialised = true;
        } else {
            c0 = nested_add(c0, kernel[i][0]);
        }
    }
    if (!c0_initialised) {
        return integral_correlation_lag<Id>(entry_type{});
    }

    // Omega = sum_i [ kernel[i][0] + 2 * sum_{offset>0} kernel[i][offset] ]
    // Accounts for positive AND negative lags via symmetry of the stationary
    // process: kernel[i][-offset] = kernel[i-offset][offset]^T, and summing over
    // all i folds the two directions into a factor of 2 on the positive side.
    bool omega_initialised = false;
    entry_type omega{};
    auto accumulate = [&](const entry_type& e) {
        if (!omega_initialised) {
            omega = e;
            omega_initialised = true;
        } else {
            omega = nested_add(omega, e);
        }
    };
    for (std::size_t i = 0; i < kernel.size(); ++i) {
        if (kernel[i].empty()) continue;
        accumulate(kernel[i][0]);
        for (std::size_t offset = 1; offset < kernel[i].size(); ++offset) {
            accumulate(kernel[i][offset]);
            accumulate(kernel[i][offset]);  // symmetric contribution from negative lag
        }
    }

    // Sandwich via the existing normalize_series_cross_correlation_block helper,
    // which dispatches on scalar/matrix/ParameterIndexed entry types.
    auto maybe_sandwich =
        normalize_series_cross_correlation_block(c0, omega, c0, rtol, atol,
                                                 "integral_correlation_lag C(0) left",
                                                 "integral_correlation_lag Omega",
                                                 "integral_correlation_lag C(0) right");
    if (!maybe_sandwich) return Error(maybe_sandwich.error()());
    return integral_correlation_lag<Id>(std::move(maybe_sandwich.value()));
}

template <class Id>
Maybe_error<Series_Moment_report<Id>> series_moment_report(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return_error<Series_Moment_report<Id>> Error{"series_moment_report: "};

    auto maybe_corr = series_cross_correlation(stats, rtol, atol);
    if (!maybe_corr)
        return Error(maybe_corr.error()());

    auto forward = series_cross_correlation_lag_forward<Id>(maybe_corr.value());

    auto maybe_integral = series_integral_correlation_lag<Id>(stats, rtol, atol);
    if (!maybe_integral) {
        return Error(maybe_integral.error()());
    }

    return Maybe_error<Series_Moment_report<Id>>(Series_Moment_report<Id>(
        get<count<Id>>(stats())(), get<series_mean<Id>>(stats())(),
        series_diagonal_covariance(stats), std::move(maybe_corr.value()), std::move(forward),
        std::move(maybe_integral.value())));
}

// Extract just the diagonal entries from the per-sample covariance blocks.
// For scalar Id, each block is a scalar and the diagonal IS the full block.
// For vector Id, each block is a p x p matrix and we extract the p diagonal
// entries, storing them as a p-vector per sample (type = mean_value_type_t<Id>).
template <class Id>
series_diagonal_variance<Id> extract_series_diagonal_variance(
    const Series_Moment_statistics<Id, true>& stats) {
    using out_entry = typename series_diagonal_variance<Id>::entry_type;
    using out_vec = typename series_diagonal_variance<Id>::diagonal_type;

    const auto& kernel = get<variance_kernel_t<Id, true>>(stats())();
    out_vec diag_out;
    diag_out.reserve(kernel.size());
    for (const auto& row : kernel) {
        if (row.empty()) {
            diag_out.emplace_back();
            continue;
        }
        const auto& block = row.front();  // lag-0 entry = per-sample covariance block
        if constexpr (std::is_arithmetic_v<std::decay_t<decltype(block)>>) {
            // Scalar Id: block is scalar, already the full "diagonal".
            diag_out.push_back(static_cast<out_entry>(block));
        } else if constexpr (requires { block.value(); block.parameters_ptr(); }) {
            // ParameterIndexed wrapper.
            if constexpr (std::is_constructible_v<out_entry, Matrix<double>&&,
                                                  decltype(block.parameters_ptr())>) {
                // entry_type's inner value is a plain Matrix<double> — extract
                // the covariance diagonal into a p-vector.
                const auto& mat = block.value();
                const std::size_t n = mat.nrows();
                Matrix<double> diag_vec(n, 1, false);
                for (std::size_t i = 0; i < n; ++i) {
                    diag_vec(i, 0) = mat(i, i);
                }
                diag_out.push_back(out_entry(std::move(diag_vec), block.parameters_ptr()));
            } else if constexpr (requires {
                                     typename out_entry::value_type(std::size_t{}, std::size_t{},
                                                                    double{});
                                     std::declval<typename out_entry::value_type&>().set(
                                         std::size_t{}, std::size_t{}, double{});
                                 }) {
                // Matrix-valued Id (e.g. GFI as SymPosDefMatrix<double>).
                // prod_XY<true>(Sym, Sym) packs via to_vector_half into an
                // m-vector (m = p(p+1)/2) then xxt → m×m block. Its diagonal
                // holds Var(GFI[i,j]) for each unique (i,j). Scatter those
                // back into a p×p symmetric matrix of per-entry variances —
                // compact (stores p(p+1)/2 values, matches entry_type shape).
                using Inner = typename out_entry::value_type;
                const auto& mat = block.value();
                const std::size_t m = mat.nrows();
                const std::size_t p =
                    (static_cast<std::size_t>(std::sqrt(8.0 * static_cast<double>(m) + 1.0)) - 1) /
                    2;
                Inner diag_mat(p, p, 0.0);
                std::size_t k = 0;
                for (std::size_t i = 0; i < p; ++i) {
                    for (std::size_t j = i; j < p; ++j) {
                        diag_mat.set(i, j, mat(k, k));
                        ++k;
                    }
                }
                diag_out.push_back(out_entry(std::move(diag_mat), block.parameters_ptr()));
            } else {
                diag_out.emplace_back();
            }
        } else if constexpr (requires { block.nrows(); block.ncols(); block(std::size_t{0}, std::size_t{0}); }) {
            if constexpr (std::is_constructible_v<out_entry, Matrix<double>&&>) {
                // Plain matrix: extract diagonal as a column vector.
                const std::size_t n = block.nrows();
                Matrix<double> diag_vec(n, 1, false);
                for (std::size_t i = 0; i < n; ++i) {
                    diag_vec(i, 0) = block(i, i);
                }
                diag_out.push_back(out_entry(std::move(diag_vec)));
            } else {
                // Matrix-valued Id whose entry_type isn't Matrix<double>
                // (e.g. SymPosDefMatrix): keep the block as-is.
                diag_out.push_back(out_entry(block));
            }
        } else {
            // Fallback: copy the block as-is. Safe for exotic entry types but
            // silently wrong if the block has non-trivial off-diagonals.
            diag_out.push_back(out_entry(block));
        }
    }
    return series_diagonal_variance<Id>(std::move(diag_out));
}

// Report builders — one per level.

template <class Id>
Maybe_error<Report_integral<Id>> make_report_integral(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return_error<Report_integral<Id>> Error{"make_report_integral: "};

    auto maybe_integral = series_integral_correlation_lag<Id>(stats, rtol, atol);
    if (!maybe_integral) {
        return Error(maybe_integral.error()());
    }
    return Maybe_error<Report_integral<Id>>(
        Report_integral<Id>(get<count<Id>>(stats())(), std::move(maybe_integral.value())));
}

template <class Id>
Maybe_error<Report_local_var<Id>> make_report_local_var(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return_error<Report_local_var<Id>> Error{"make_report_local_var: "};

    auto maybe_corr = series_cross_correlation(stats, rtol, atol);
    if (!maybe_corr) {
        return Error(maybe_corr.error()());
    }
    auto forward = series_cross_correlation_lag_forward<Id>(maybe_corr.value());
    auto maybe_integral = series_integral_correlation_lag<Id>(stats, rtol, atol);
    if (!maybe_integral) {
        return Error(maybe_integral.error()());
    }
    auto diag = extract_series_diagonal_variance<Id>(stats);

    return Maybe_error<Report_local_var<Id>>(Report_local_var<Id>(
        get<count<Id>>(stats())(), get<series_mean<Id>>(stats())(), std::move(diag),
        std::move(forward), std::move(maybe_integral.value())));
}

template <class Id>
Maybe_error<Report_local_cov<Id>> make_report_local_cov(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return_error<Report_local_cov<Id>> Error{"make_report_local_cov: "};

    auto maybe_corr = series_cross_correlation(stats, rtol, atol);
    if (!maybe_corr) {
        return Error(maybe_corr.error()());
    }
    auto forward = series_cross_correlation_lag_forward<Id>(maybe_corr.value());
    auto maybe_integral = series_integral_correlation_lag<Id>(stats, rtol, atol);
    if (!maybe_integral) {
        return Error(maybe_integral.error()());
    }

    return Maybe_error<Report_local_cov<Id>>(Report_local_cov<Id>(
        get<count<Id>>(stats())(), get<series_mean<Id>>(stats())(),
        series_diagonal_covariance(stats), std::move(forward),
        std::move(maybe_integral.value())));
}

template <class Id>
Maybe_error<Report_cross<Id>> make_report_cross(
    const Series_Moment_statistics<Id, true>& stats, double rtol = 1e-10, double atol = 0.0) {
    return series_moment_report<Id>(stats, rtol, atol);
}


template <typename T, class F = std::identity>
requires (std::is_arithmetic_v<std::decay_t<std::invoke_result_t<F, const T&>>>)
auto get_mean_Probits(std::vector<T> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) {
    assert(!bootstrap_estimates.empty());

    std::vector<double> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations

    double sum = 0.0;

    // 1. Evaluate, filter NaNs, and Accumulate (O(N))
    for (const auto& x : bootstrap_estimates) {
        double value = static_cast<double>(std::invoke(f, x));
        if (!std::isnan(value)) {
            valid_estimates.push_back(value);
            sum += value;
        }
    }

    std::size_t n = valid_estimates.size();

    // Guard against an array entirely filled with NaNs
    if (n == 0) {
        return std::make_tuple(std::numeric_limits<double>::quiet_NaN(), std::map<double, double>{}, 0UL);
    }

    double mean = sum / static_cast<double>(n);

    // 2. Direct Sort (Cache-friendly, no NaN-logic required)
    std::sort(valid_estimates.begin(), valid_estimates.end());

    // 3. Map Percentiles
    std::map<double, double> out;

    for (double ci : cis) {
        double real_idx = std::clamp(ci * (n + 1) - 1.0, 0.0, static_cast<double>(n - 1));

        std::size_t low = static_cast<std::size_t>(std::floor(real_idx));
        std::size_t high = static_cast<std::size_t>(std::ceil(real_idx));

        double weight = real_idx - static_cast<double>(low);

        double v_low  = valid_estimates[low];
        double v_high = valid_estimates[high];

        out[ci] = (1.0 - weight) * v_low + weight * v_high;
    }

    return std::make_tuple(mean, std::move(out), n);
}

template <typename T, class F = std::identity>
requires requires(T x, F f) { { std::invoke(f, x)() }; }
auto get_mean_Probits(std::vector<T> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) {
    assert(!bootstrap_estimates.empty());

    // Recursive call: Dispatch to scalar, Matrix, or Vector overloads
    // by unwrapping the Var/Functor layer.
    auto result = get_mean_Probits(bootstrap_estimates, cis,
        [&f](const auto& x) -> decltype(auto) {
            return std::invoke(f, x)();
        });

    // We can return the pair directly because the inner call already
    // produced the unwrapped value types (double, Matrix, etc.)
    return result;
}

template<class M, class F = std::identity>
 requires (is_Matrix_v<std::decay_t<std::invoke_result_t<F, const M&>>>)
auto get_mean_Probits(
    std::vector<M> const& bootstrap_estimates,
    const std::set<double>& cis, F f = {})
{
    assert(!bootstrap_estimates.empty());

    // Determine exactly what f returns
    using ProjResult = std::invoke_result_t<F, const M&>;
    using RawMatrix  = std::decay_t<ProjResult>;

    // 1. The Safe Storage Type:
    // If f returns an lvalue reference, store a reference_wrapper.
    // If f returns a value, store the raw matrix.
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawMatrix
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations

    // 2. Populate the valid estimates
    for (const auto& x : bootstrap_estimates) {
        decltype(auto) value = std::invoke(f, x);
        if (value.size() > 0) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                valid_estimates.emplace_back(std::ref(value));
            } else {
                valid_estimates.emplace_back(std::move(value));
            }
        }
    }

    if (valid_estimates.empty()) {
        return std::make_tuple(RawMatrix{}, std::map<double, RawMatrix>{}, 0UL);
    }

    // 3. Extract the first matrix safely
    const RawMatrix& first_matrix = [] (const StorageType& item) -> const RawMatrix& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get();
        else return item;
    }(valid_estimates.front());

    RawMatrix mean_matrix(first_matrix.nrows(), first_matrix.ncols());
    std::map<double, RawMatrix> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, RawMatrix(first_matrix.nrows(), first_matrix.ncols()));
    }

    std::size_t n = std::numeric_limits<std::size_t>::max();

    // 4. Delegate to the scalar 1D version
    for (std::size_t i = 0; i < mean_matrix.size(); ++i) {

        // Notice: `valid_estimates` already holds the matrices.
        // We just need a lambda that extracts the i-th element.
        auto prob = get_mean_Probits(valid_estimates, cis, [i](const StorageType& mat_ref) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                return mat_ref.get()[i];
            } else {
                return mat_ref[i];
            }
        });

        mean_matrix[i] = std::get<0>(prob);
        for (const auto& [ci, value] : std::get<1>(prob)) {
            probits_map[ci][i] = value;
        }

        if (std::get<2>(prob) < n) {
            n = std::get<2>(prob);
        }
    }

    if (n == std::numeric_limits<std::size_t>::max()) n = 0;

    return std::make_tuple(std::move(mean_matrix), std::move(probits_map), n);
}

template<class V, class F = std::identity> requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const V&>>, std::vector>)
auto get_mean_Probits(
    std::vector<V> const& bootstrap_estimates,
    const std::set<double>& cis, F f = {})
{
    assert(!bootstrap_estimates.empty());

    // Determine exactly what f returns
    using ProjResult = std::invoke_result_t<F, const V&>;
    using RawVector  = std::decay_t<ProjResult>;

    // Fallback: reference_wrapper for lvalues, raw vector for prvalues
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawVector
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations

    // 1. Evaluate, filter empty vectors, and store safely
    for (const auto& x : bootstrap_estimates) {
        decltype(auto) value = std::invoke(f, x);
        if (!value.empty()) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                valid_estimates.emplace_back(std::ref(value));
            } else {
                valid_estimates.emplace_back(std::move(value));
            }
        }
    }

    // Guard against an array entirely filled with empty vectors
    if (valid_estimates.empty()) {
        return std::make_tuple(RawVector{}, std::map<double, RawVector>{}, 0UL);
    }

    // 2. Extract the first vector safely to get sizes
    const RawVector& first_vector = [] (const StorageType& item) -> const RawVector& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get();
        else return item;
    }(valid_estimates.front());

    // Check size consistency only on the valid vectors
    assert(std::all_of(valid_estimates.begin(), valid_estimates.end(), [&](const auto& item) {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get().size() == first_vector.size();
        else return item.size() == first_vector.size();
    }));

    RawVector mean_vector(first_vector.size());
    std::map<double, RawVector> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, RawVector(first_vector.size()));
    }

    std::size_t n = std::numeric_limits<std::size_t>::max();

    // 3. Delegate to the scalar 1D version
    for (std::size_t i = 0; i < mean_vector.size(); ++i) {

        auto prob = get_mean_Probits(valid_estimates, cis, [i](const StorageType& vec_ref) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                return vec_ref.get()[i];
            } else {
                return vec_ref[i];
            }
        });

        mean_vector[i] = std::get<0>(prob);
        for (const auto& [ci, value] : std::get<1>(prob)) {
            probits_map[ci][i] = value;
        }

        if (std::get<2>(prob) < n) {
            n = std::get<2>(prob);
        }
    }

    if (n == std::numeric_limits<std::size_t>::max()) n = 0;

    return std::make_tuple(std::move(mean_vector), std::move(probits_map), n);
}

template <class ParamIndexed, class F = std::identity>
 requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const ParamIndexed&>>, var::ParameterIndexed>)
auto get_mean_Probits(std::vector<ParamIndexed> const& bootstrap_estimates,
                      const std::set<double>& cis, F f = {})
{
    assert(!bootstrap_estimates.empty());

    using ProjResult = std::invoke_result_t<F, const ParamIndexed&>;
    using RawParamIdx = std::decay_t<ProjResult>;
    using ValueT = typename RawParamIdx::value_type;
    using Params = typename RawParamIdx::params_type;

    // Fallback: reference_wrapper for lvalues, raw object for prvalues
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawParamIdx
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size());

    const Params* metadata = nullptr;

    // 1. Evaluate and store safely (Zero-Copy)
    for (const auto& estimate : bootstrap_estimates) {
        decltype(auto) projected = std::invoke(f, estimate);

        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            valid_estimates.emplace_back(std::ref(projected));
        } else {
            valid_estimates.emplace_back(std::move(projected));
        }

        // Grab metadata from the first item that has it
        if (metadata == nullptr) {
            const RawParamIdx& raw_item = [&]() -> const RawParamIdx& {
                if constexpr (std::is_lvalue_reference_v<ProjResult>) return valid_estimates.back().get();
                else return valid_estimates.back();
            }();

            if (raw_item.has_parameters()) {
                metadata = raw_item.parameters_ptr();
            }
        }
    }

    // 2. Delegate to the underlying ValueT implementation.
    // Instead of copying the values, we project them on the fly!
    auto prob = get_mean_Probits(valid_estimates, cis, [](const StorageType& item) -> const ValueT& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            return item.get().value();
        } else {
            return item.value();
        }
    });

    // 3. Re-wrap the results with the metadata
    var::ParameterIndexed<ValueT, Params> mean_out(std::move(std::get<0>(prob)), metadata);

    std::map<double, var::ParameterIndexed<ValueT, Params>> probits_out;
    for (auto& [level, val] : std::get<1>(prob)) {
        probits_out.emplace(level, var::ParameterIndexed<ValueT, Params>(std::move(val), metadata));
    }

    return std::make_tuple(std::move(mean_out), std::move(probits_out), std::get<2>(prob));
}

template<class F=std::identity, class VecSpace>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) ;


namespace old{
template <class Id>
class Probit_statistics
    : public var::Var<Probit_statistics<Id>,
    var::Vector_Space<mean<Id>,Probit_025<Id>, Probit_975<Id>>>
     {
   public:
    using base_type =var::Var<Probit_statistics<Id>,
    var::Vector_Space<mean<Id>,Probit_025<Id>, Probit_975<Id>>>
    ;

    using value_type= value_type_t<mean<Id>>;

    static constexpr double LOWER_PERCENTILE = 0.025;
    static constexpr double UPPER_PERCENTILE = 0.975;

    Probit_statistics()
        : base_type{var::Vector_Space{mean<Id>{}, Probit_025<Id>{}, Probit_975<Id>{}}} {}

    Probit_statistics(value_type&& m, value_type&& p025, value_type&& p975):
        base_type{var::Vector_Space{mean<Id>(std::move( m)), Probit_025<Id>(std::move(p025)), Probit_975<Id>(std::move(p975))}} {}

    Probit_statistics(std::vector<value_type>& values) {
        auto [mean_probit, probits] = get_mean_Probits(values, {LOWER_PERCENTILE, UPPER_PERCENTILE});
        *this = Probit_statistics(std::move(mean_probit), std::move(probits[0]), std::move(probits[1]));


    }
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { {f(v)} -> std::convertible_to<value_type_t<Id>>; }
    Probit_statistics(std::vector<VectorSpace>const& xs,  F&& f)
    {
        std::vector<value_type > values;
        values.reserve(xs.size());
        for (const auto& x : xs) {
            values.push_back(std::forward<F>(f)(x));
        }
        auto [mean_probit, probbits] = get_mean_Probits(values, {LOWER_PERCENTILE, UPPER_PERCENTILE});
        *this = Probit_statistics(std::move(mean_probit), std::move(probbits[0]), std::move(probbits[1]));
    }



    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};
}

template <class Id>
class Probit_statistics
    : public var::Var<Probit_statistics<Id>,
                      var::Vector_Space<mean<Id>, Probits<Id>, count<Id>>>
{
   public:
    using base_type = var::Var<Probit_statistics<Id>,
                               var::Vector_Space<mean<Id>, Probits<Id>, count<Id>>>;

    using value_type = value_type_t<mean<Id>>;
    using raw_value_type = value_type_t<Id>;

private:
    // Helper constructor to unpack the tuple directly into the main constructor
    template <typename TTuple>
    Probit_statistics(TTuple&& data)
        : Probit_statistics(std::move(std::get<0>(data)),
                            std::move(std::get<1>(data)),
                            std::get<2>(data)) {}

public:
    Probit_statistics()
        : base_type{var::Vector_Space{mean<Id>{}, Probits<Id>{}, count<Id>{}}} {}

    Probit_statistics(value_type&& m, std::map<double,value_type>&& pbits, std::size_t n)
        : base_type{var::Vector_Space{mean<Id>(std::move(m)), Probits<Id>(std::move(pbits)), count<Id>(n)}} {}

    template <class Dummy = void>
    requires (!std::is_same_v<raw_value_type, value_type> &&
              std::constructible_from<Id, raw_value_type>)
    Probit_statistics(raw_value_type&& m, std::map<double, raw_value_type>&& pbits, std::size_t n)
        : base_type{var::Vector_Space{
              mean<Id>(Id(std::move(m))),
              Probits<Id>([&pbits]() {
                  std::map<double, value_type> wrapped;
                  for (auto& [ci, value] : pbits) {
                      wrapped.emplace(ci, value_type(Id(std::move(value))));
                  }
                  return wrapped;
              }()),
              count<Id>(n)}} {}

    // Zero-copy delegation directly to the tuple constructor
    Probit_statistics(std::vector<value_type> const& values, const std::set<double>& cis)
        : Probit_statistics(get_mean_Probits(values, cis)) {}

    // Zero-copy functional pipeline! No intermediate vector allocations.
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { { std::invoke(f, v) } -> std::convertible_to<value_type_t<Id>>; }
    Probit_statistics(std::vector<VectorSpace> const& xs, F f, const std::set<double>& cis)
        : Probit_statistics(get_mean_Probits(xs, cis, std::move(f))) {}

    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};

template<class F, class VecSpace, class... Vs>
requires (std::is_same_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space<Vs...>>)
auto get_mean_Probits_impl(std::vector<VecSpace> const& bootstrap_estimates,
                           const std::set<double>& cis, F f,
                           std::type_identity<var::Vector_Space<Vs...>>)
{
    using ProjResult = std::invoke_result_t<F, const VecSpace&>;
    using RawVecSpace = std::decay_t<ProjResult>;

    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawVecSpace
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size());

    // Evaluate 'f' EXACTLY ONCE per estimate
    for (const auto& v : bootstrap_estimates) {
        decltype(auto) val = std::invoke(f, v);
        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            valid_estimates.emplace_back(std::ref(val));
        } else {
            valid_estimates.emplace_back(std::move(val));
        }
    }

    auto probit = var::Vector_Space<Probit_statistics<Vs>...>(
        Probit_statistics<Vs>(
            valid_estimates,
            [](const StorageType& st) -> decltype(auto) {
                if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                    return get<Vs>(st.get())();
                } else {
                    return get<Vs>(st)();
                }
            },
            cis)...
    );

    auto meanout = var::Vector_Space<Vs...>(
        Vs(std::move(get<mean<Vs>>(get<Probit_statistics<Vs>>(probit)())()))...
    );

    // FIX: Extract the minimum count as a std::size_t instead of a Vector_Space
    std::size_t min_n = valid_estimates.size();
    (([&]{
        std::size_t current_n = get<count<Vs>>(get<Probit_statistics<Vs>>(probit)())();
        if (current_n < min_n) {
            min_n = current_n;
        }
    }()), ...);

    std::map<double, var::Vector_Space<Vs...>> probitsout;
    for (double ci : cis) {
        probitsout.emplace(
            ci,
            var::Vector_Space<Vs...>(
                Vs([&]() {
                    const auto& component = get<Probit_statistics<Vs>>(probit)();
                    const auto& probits_map = get<Probits<Vs>>(component)();
                    if (auto it = probits_map.find(ci); it != probits_map.end()) {
                        return it->second;
                    }
                    return get<mean<Vs>>(component)();
                }())...
            )
        );
    }

    // Return min_n so it perfectly matches the Probit_statistics tuple constructor
    return std::make_tuple(std::move(meanout), std::move(probitsout), min_n);
}



template<class F, class VecSpace>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f) {
   using RVecSp = std::decay_t<std::invoke_result_t<F, const VecSpace&>>;
   return get_mean_Probits_impl(bootstrap_estimates, cis, std::move(f), std::type_identity<RVecSp>{});
}

#endif  // MOMENT_STATISTICS_H
