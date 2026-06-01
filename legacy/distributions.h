#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <random>
#include <sstream>
#include <type_traits>
#include <utility>

#include "lapack_headers.h"
#include "matrix.h"
#include "maybe_error.h"
#include "normal_distribution.h"
#include "parameter_indexed.h"
#include "random_samplers.h"
#include "variables.h"

template <class Distribution>
concept is_Distribution = requires(Distribution& m, Distribution const& m_const) {
    { m(std::declval<mt_64i&>()) };

    { m_const.logP(m(std::declval<mt_64i&>())) } -> std::convertible_to<Maybe_error<double>>;
};

template <class Distribution, class T>
concept is_Distribution_of = requires(Distribution& m, Distribution const& m_const) {
    { m(std::declval<mt_64i&>()) } -> std::convertible_to<T>;

    { m_const.logP(m(std::declval<mt_64i&>())) } -> std::convertible_to<Maybe_error<double>>;
};

template <class Sampler, class ProbabilityFunction>
class Custom_Distribution {
    std::size_t m_size;
    Sampler m_s;
    ProbabilityFunction m_p;

   public:
    using T = std::invoke_result_t<Sampler, mt_64i&>;
    Custom_Distribution(std::size_t t_size, Sampler&& s, ProbabilityFunction&& f)
        : m_size{t_size}, m_s{std::move(s)}, m_p{std::move(f)} {
    }

    auto operator()(mt_64i& mt) {
        return std::invoke(m_s, mt);
    }
    Maybe_error<double> logP(const T& x) const {
        return std::invoke(m_p, x);
    }

    auto size() const {
        return m_size;
    }
};

template <class Distribution>
concept is_vector_sampler = requires(Distribution& m) {
    { m(std::declval<mt_64i&>, 0) } -> index_accesible;
};

template <class T>
    requires(index_accesible<T>)
auto& get_at(T& v, std::size_t i) {
    return v[i];
}

template <class T>
auto get_at(T& v, std::size_t i) -> std::decay_t<decltype(v[i])> {
    return v[i];
}
inline double get_at(double const& x, std::size_t i) {
    assert(i == 0);
    return x;
}

inline double& get_at(double& x, std::size_t i) {
    assert(i == 0);
    return x;
}

template <class T>
    requires(has_size<T>)
std::size_t size(const T& x) {
    return x.size();
}

template <class T>
    requires(!has_size<T> && has_size<std::decay_t<decltype(std::declval<T>()())>>)
std::size_t size(const T& x) {
    return x().size();
}

constexpr std::size_t size(double) {
    return 1;
}

template <class F, class T>
auto operator|(F&& f, T&& x) -> std::invoke_result_t<F, T> {
    return std::invoke(std::forward<F>(f), std::forward<T>(x));
}

inline auto& append_to(Matrix<double>& m, std::size_t) {
    return m;
}

template <class X, class... Xs>
auto& append_to(Matrix<double>& m, std::size_t i, X&& x, Xs&&... xs) {
    for (std::size_t n = 0; n < size(x); ++n) m[i + n] = get_at(x, n);
    return append_to(m, i + size(x), std::forward<Xs>(xs)...);
}

template <class... Xs>
auto concatenate_to_columns(Xs&&... xs) {
    auto n = (size(xs) + ...);
    auto out = Matrix<double>(1, n);
    std::size_t i = 0;
    out = append_to(out, i, std::forward<Xs>(xs)...);
    return out;
}

template <class Dist>
concept Multivariate = requires(Dist& d) {
    { d(std::declval<mt_64i&>()) } -> std::convertible_to<Matrix<double>>;
};

inline double logP_impl(const Matrix<double>&, std::size_t, double partial_logP) {
    return partial_logP;
}
template <class Dist, class... Ds>
    requires(Multivariate<Dist>)
Maybe_error<double> logP_impl(const Matrix<double>& x, std::size_t ipos, double partial_logP,
                              const Dist& d, const Ds&... ds) {
    auto n = d.size();
    auto out = Matrix<double>(1, size(d));
    for (std::size_t i = 0; i < n; ++i) out[i] = x[ipos + i];
    auto logPi = d.logP(out);
    if (logPi)
        return logP_impl(x, ipos + n, partial_logP + logPi.value(), ds...);
    else
        return logPi.error() + "\n log_impl";
}
template <class Dist, class... Ds>
    requires(!Multivariate<Dist>)
Maybe_error<double> logP_impl(const Matrix<double>& x, std::size_t ipos, double partial_logP,
                              const Dist& d, const Ds&... ds) {
    auto logPi = d.logP(x[ipos]);
    if (logPi)
        return logP_impl(x, ipos + 1, partial_logP + logPi.value(), ds...);
    else
        return logPi.error() + "\n log_impl";
}

template <class D>
    requires(is_Distribution<D>)
auto sample(mt_64i& mt, D& d) {
    return d(mt);
}

template <class D>
    requires(is_Distribution<D>)
auto sample(mt_64i& mt, D const& d) {
    return d(mt);
}

template <class D>
    requires(is_Distribution<D>)
auto sampler(const D& d) {
    return D(d);
}

template <class D>
    requires(is_vector_sampler<D>)
auto sample(mt_64i& mt, D& d, std::size_t n) {
    return d(mt, n);
}

template <class D>
    requires(is_vector_sampler<D>)
auto sample(mt_64i& mt, D const& d, std::size_t n) {
    return d(mt, n);
}

template <class D>
    requires(is_Distribution<D>)
auto sample(mt_64i& mt, D&& d, std::size_t nrows, std::size_t ncols) {
    Matrix<double> out(nrows, ncols, false);
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = std::forward<D>(d)(mt);
    return out;
}

template <class C_double>
auto logit(const C_double& x) {
    using std::log;
    return log(x / (1.0 - x));
}

template <class C_double>
auto logit_inv(const C_double& x) {
    using std::exp;
    return exp(x) / (1.0 + exp(x));
}

template <class... ds>
    requires(is_Distribution<ds> && ...)
class distributions : public ds... {
   public:
    auto operator()(mt_64i& mt) {
        return concatenate_to_columns(ds::operator()(mt)...);
    }

    Maybe_error<double> logP(const Matrix<double>& x) const {
        return logP_impl(x, 0ul, 0.0, static_cast<ds const&>(*this)...);
    }

    explicit distributions(ds&&... d) : ds{std::move(d)}... {
    }
    explicit distributions(ds const&... d) : ds{d}... {
    }
};

class logEv : public var::Var<logEv, double> {
   public:
    friend std::string className(logEv) {
        return "logEvidence";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("logEvidence not finite: " + std::to_string((*this)()));
        return true;
    }
};

class logL : public var::Var<logL, double> {
    public:

    friend std::string className(logL) {
        return "logL";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("logL not finite: " + std::to_string((*this)()));
        return true;
    }
};

using parameter_vector_payload = var::ParameterIndexed<Matrix<double>, var::Parameters_transformed>;
using parameter_spd_payload =
    var::ParameterIndexed<SymPosDefMatrix<double>, var::Parameters_transformed>;
using parameter_symmetric_payload =
    var::ParameterIndexed<SymmetricMatrix<double>, var::Parameters_transformed>;

// Default relative tolerance for treating a parameter direction as
// structurally inactive (zero diagonal of an SPD payload). Matches the rtol
// used by lapack::compute_psd_decomp on the bigger picture (1e-10) but is set
// a couple of orders of magnitude tighter so we only drop parameters that are
// truly outside the gradient chain (e.g. Num_ch_mean under the lifted micro
// path), not ill-conditioned ones the PSD machinery should still handle.
inline constexpr double k_active_subspace_rtol = 1e-12;
inline constexpr double k_logdet_subspace_rtol = 1e-10;
inline constexpr double k_logdet_subspace_atol = 0.0;

// Indices of parameters with zero diagonal of M — i.e. structurally inactive
// directions in a Fisher-information / score-covariance / IDM payload. A
// parameter shows up here when *no* derivative chain in the model touches it
// (e.g. Num_ch_mean under the lifted micro path, where it sets the dimension
// of the microstate space rather than entering the gradient continuously).
inline std::vector<std::size_t> inactive_parameter_indices(
    parameter_spd_payload const& M, double rtol = k_active_subspace_rtol) {
    auto const& mat = M.value();
    std::size_t p = mat.nrows();
    double max_diag = 0.0;
    for (std::size_t i = 0; i < p; ++i) {
        max_diag = std::max(max_diag, mat(i, i));
    }
    double thr = rtol * std::max(max_diag, 1.0);
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < p; ++i) {
        if (mat(i, i) <= thr) {
            out.push_back(i);
        }
    }
    return out;
}

// Names of the inactive parameters detected by inactive_parameter_indices.
// Uses the payload's Parameters_transformed for the name lookup.
inline std::vector<std::string> inactive_parameter_names(
    parameter_spd_payload const& M, double rtol = k_active_subspace_rtol) {
    auto idx = inactive_parameter_indices(M, rtol);
    std::vector<std::string> names;
    names.reserve(idx.size());
    if (M.has_parameters()) {
        auto const& full = M.parameters().parameters().names();
        for (auto i : idx) {
            if (i < full.size()) {
                names.push_back(full[i]);
            }
        }
    }
    return names;
}

// log-determinant of the parameter-indexed SPD matrix restricted to its
// non-null eigenspace — i.e. excluding the zero-diagonal (structurally
// inactive) parameter directions. For full-rank payloads this is the regular
// log_det; for rank-deficient ones it returns a finite pseudo-log-determinant
// instead of −∞ / NaN. Use this in diagnostic chains where one or more model
// parameters are categorically absent from the gradient (e.g. Num_ch_mean for
// micro algorithms).
inline double logdet_active_subspace(parameter_spd_payload const& M,
                                      double rtol = k_logdet_subspace_rtol,
                                      double atol = k_logdet_subspace_atol,
                                      std::string const& name = "logdet_active_subspace") {
    return lapack::logdet_subspace(M.value(), name, rtol, atol);
}
class elogL : public var::Var<elogL, double> {
   public:
    friend std::string className(elogL) {
        return "elogL";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("elogL not finite: " + std::to_string((*this)()));
        return true;
    }
};
class vlogL : public var::Constant<vlogL, double> {
   public:
    friend std::string className(vlogL) {
        return "vlogL";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("vlogL not finite: " + std::to_string((*this)()));
        if ((*this)() <= 0)
            return error_message("vlogL negative or cero: " + std::to_string((*this)()));

        return true;
    }
};

// Wall-clock time (seconds) of one top-level likelihood evaluation. Captured
// at the dispatch sites in src/core/likelihood.cpp around dlog_Likelihood_micro
// / dlogLikelihoodPredictions / log_Likelihood_micro / logLikelihood. Stored
// per-recording in the returned MacroState so the diagnostic preset can
// aggregate Sum<evaluation_time> like any other scalar — gives per-algorithm
// mean ± std cost across the bootstrap of simulations.
class evaluation_time : public var::Var<evaluation_time, double> {
   public:
    friend std::string className(evaluation_time) { return "evaluation_time"; }
};

class dlogL : public var::Constant<dlogL, parameter_vector_payload> {
    using base_type = var::Constant<dlogL, parameter_vector_payload>;

   public:
    using base_type::base_type;
    dlogL() = default;
    dlogL(Matrix<double> value, var::Parameters_transformed const& params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}
    dlogL(Matrix<double> value, var::Parameters_transformed const* params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}

    friend std::string className(dlogL) {
        return "dlogL";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite(var::fullsum((*this)()))) {
            std::stringstream ss;
            ss << "dlogL  not finite: " << (*this)();
            return error_message(ss.str());
        }
        return true;
    }
};

class d2logL : public var::Constant<d2logL, SymPosDefMatrix<double>> {
   public:
    friend std::string className(d2logL) {
        return "d2logL";
    }
};  





class Residual_based_ESS : public var::Constant<Residual_based_ESS, double> {
   public:
    friend std::string className(Residual_based_ESS) {
        return "Residual_based_ESS";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("Residual_based_ESS not finite: " + std::to_string((*this)()));
        if ((*this)() <= 0)
            return error_message("Residual_based_ESS negative or cero: " + std::to_string((*this)()));

        return true;
    }
};

class Score_det_based_ESS : public var::Constant<Score_det_based_ESS, double> {
   public:
    friend std::string className(Score_det_based_ESS) {
        return "Score_det_based_ESS";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("Score_det_based_ESS not finite: " + std::to_string((*this)()));
        if ((*this)() <= 0)
            return error_message("Score_det_based_ESS negative or cero: " + std::to_string((*this)()));

        return true;
    }
};

class Score_trace_based_ESS : public var::Constant<Score_trace_based_ESS, double> {
   public:
    friend std::string className(Score_trace_based_ESS) {
        return "Score_trace_based_ESS";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite((*this)()))
            return error_message("Score_trace_based_ESS not finite: " + std::to_string((*this)()));
        if ((*this)() <= 0)
            return error_message("Score_trace_based_ESS negative or cero: " + std::to_string((*this)()));

        return true;
    }
};


class Gaussian_Fisher_Information
    : public var::Constant<Gaussian_Fisher_Information, parameter_spd_payload> {
    using base_type = var::Constant<Gaussian_Fisher_Information, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Gaussian_Fisher_Information() = default;
    Gaussian_Fisher_Information(SymPosDefMatrix<double> value,
                                var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Gaussian_Fisher_Information(SymPosDefMatrix<double> value,
                                var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Gaussian_Fisher_Information) {
        return "Gaussian_Fisher_Information";
    }
};

// Numerical (observed) Fisher information for the *whole* recording, computed
// by finite-differencing the analytic score:
//   F = -∂(score)/∂θ = -∂²logL/∂θ²    (stored as the positive-definite FIM-form)
// Built once per recording from 2·n_active extra calls to the score-producing
// likelihood pass at θ₀ ± h·e_i for each active parameter coordinate i, then
// symmetrized and negated. Truth for both macro and micro paths by construction
// — used as the H reference in Likelihood_Information_Distortion instead of the
// per-step Gaussian_Fisher_Information sum (which is exact for macro but only
// an approximation for the micro mixture likelihood).
class Likelihood_Numerical_Fisher_Information
    : public var::Constant<Likelihood_Numerical_Fisher_Information, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Numerical_Fisher_Information, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Numerical_Fisher_Information() = default;
    Likelihood_Numerical_Fisher_Information(SymPosDefMatrix<double> value,
                                 var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Numerical_Fisher_Information(SymPosDefMatrix<double> value,
                                 var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Numerical_Fisher_Information) {
        return "Likelihood_Numerical_Fisher_Information";
    }
};

class Score_Covariance_Matrix
    : public var::Constant<Score_Covariance_Matrix, parameter_spd_payload> {
    using base_type = var::Constant<Score_Covariance_Matrix, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Score_Covariance_Matrix() = default;
    Score_Covariance_Matrix(SymPosDefMatrix<double> value,
                            var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Score_Covariance_Matrix(SymPosDefMatrix<double> value,
                            var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Score_Covariance_Matrix) {
        return "Score_Covariance_Matrix";
    }
};  

class Score_Sample_Covariance_Matrix
    : public var::Constant<Score_Sample_Covariance_Matrix, parameter_spd_payload> {
    using base_type = var::Constant<Score_Sample_Covariance_Matrix, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Score_Sample_Covariance_Matrix() = default;
    Score_Sample_Covariance_Matrix(SymPosDefMatrix<double> value,
                                   var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Score_Sample_Covariance_Matrix(SymPosDefMatrix<double> value,
                                   var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Score_Sample_Covariance_Matrix) {
        return "Score_Sample_Covariance_Matrix";
    }
};  



class Likelihood_Information_Distortion
    : public var::Constant<Likelihood_Information_Distortion, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Information_Distortion, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Information_Distortion() = default;
    Likelihood_Information_Distortion(SymPosDefMatrix<double> value,
                                  var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Information_Distortion(SymPosDefMatrix<double> value,
                                  var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Information_Distortion) {
        return "Likelihood_Information_Distortion";
    }
};  

class Likelihood_Information_Distortion_Reconstituted
    : public var::Constant<Likelihood_Information_Distortion_Reconstituted, parameter_spd_payload> {
    using base_type =
        var::Constant<Likelihood_Information_Distortion_Reconstituted, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Information_Distortion_Reconstituted() = default;
    Likelihood_Information_Distortion_Reconstituted(SymPosDefMatrix<double> value,
                                         var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Information_Distortion_Reconstituted(SymPosDefMatrix<double> value,
                                         var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Information_Distortion_Reconstituted) {
        return "Likelihood_Information_Distortion_Reconstituted";
    }
};

// Quantifies how far the cheap analytic Gaussian_Fisher_Information G is from
// the numerical truth F = Likelihood_Numerical_Fisher_Information, in the same frame as
// the IDM:
//   Likelihood_Gaussian_Fisher_Distortion = inv(√G) · F · inv(√G)
// Equals the identity when the Gaussian-formula FIM matches the truth, so its
// eigenvalues' departure from 1 (and log-det) measures how misleading the
// single-Gaussian moment-match approximation is for the algorithm under test.
class Likelihood_Gaussian_Fisher_Distortion
    : public var::Constant<Likelihood_Gaussian_Fisher_Distortion, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Gaussian_Fisher_Distortion, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Gaussian_Fisher_Distortion() = default;
    Likelihood_Gaussian_Fisher_Distortion(SymPosDefMatrix<double> value,
                               var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Gaussian_Fisher_Distortion(SymPosDefMatrix<double> value,
                               var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Gaussian_Fisher_Distortion) {
        return "Likelihood_Gaussian_Fisher_Distortion";
    }
};

class Likelihood_Fisher_Covariance
    : public var::Constant<Likelihood_Fisher_Covariance, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Fisher_Covariance, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Fisher_Covariance() = default;
    Likelihood_Fisher_Covariance(SymPosDefMatrix<double> value,
                      var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Fisher_Covariance(SymPosDefMatrix<double> value,
                      var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Fisher_Covariance) {
        return "Likelihood_Fisher_Covariance";
    }
};  


class Likelihood_Distortion_Corrected_Covariance
    : public var::Constant<Likelihood_Distortion_Corrected_Covariance, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Distortion_Corrected_Covariance, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Distortion_Corrected_Covariance() = default;
    Likelihood_Distortion_Corrected_Covariance(SymPosDefMatrix<double> value,
                                    var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Distortion_Corrected_Covariance(SymPosDefMatrix<double> value,
                                    var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Distortion_Corrected_Covariance) {
        return "Likelihood_Distortion_Corrected_Covariance";
    }
};  

class Score_Mean : public var::Constant<Score_Mean, Matrix<double>> {
   public:
    friend std::string className(Score_Mean) {
        return "Score_Mean";
    }
};  

class Likelihood_Distortion_Induced_Bias : public var::Constant<Likelihood_Distortion_Induced_Bias, parameter_vector_payload> {
    using base_type = var::Constant<Likelihood_Distortion_Induced_Bias, parameter_vector_payload>;

   public:
    using base_type::base_type;
    Likelihood_Distortion_Induced_Bias() = default;
    Likelihood_Distortion_Induced_Bias(Matrix<double> value, var::Parameters_transformed const& params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}
    Likelihood_Distortion_Induced_Bias(Matrix<double> value, var::Parameters_transformed const* params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Distortion_Induced_Bias) {
        return "Likelihood_Distortion_Induced_Bias";
    }

    Maybe_error<bool> is_good() const {
        if (!std::isfinite(var::fullsum((*this)()))) {
            std::stringstream ss;
            ss << "Likelihood_Distortion_Induced_Bias not finite: " << (*this)();
            return error_message(ss.str());
        }
        return true;
    }
};



class Likelihood_Sample_Distortion
    : public var::Constant<Likelihood_Sample_Distortion, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Sample_Distortion, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Sample_Distortion() = default;
    Likelihood_Sample_Distortion(SymPosDefMatrix<double> value,
                             var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Sample_Distortion(SymPosDefMatrix<double> value,
                             var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Sample_Distortion) {
        return "Likelihood_Sample_Distortion";
    }
};  

class Likelihood_Correlation_Distortion
    : public var::Constant<Likelihood_Correlation_Distortion, parameter_spd_payload> {
    using base_type = var::Constant<Likelihood_Correlation_Distortion, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Likelihood_Correlation_Distortion() = default;
    Likelihood_Correlation_Distortion(SymPosDefMatrix<double> value,
                                  var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Likelihood_Correlation_Distortion(SymPosDefMatrix<double> value,
                                  var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Likelihood_Correlation_Distortion) {
        return "Likelihood_Correlation_Distortion";
    }
};  






class Grad : public var::Constant<Grad, parameter_vector_payload> {
    using base_type = var::Constant<Grad, parameter_vector_payload>;

   public:
    using base_type::base_type;
    Grad() = default;
    Grad(Matrix<double> value, var::Parameters_transformed const& params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}
    Grad(Matrix<double> value, var::Parameters_transformed const* params)
        : base_type(parameter_vector_payload(std::move(value), params)) {}

    friend std::string className(Grad) {
        return "Grad";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite(var::fullsum((*this)()))) {
            std::stringstream ss;
            ss << "Grad not finite: " << (*this)();
            return error_message(ss.str());
        }
        return true;
    }
};

class Hessian : public var::Constant<Hessian, parameter_spd_payload> {
    using base_type = var::Constant<Hessian, parameter_spd_payload>;

   public:
    using base_type::base_type;
    Hessian() = default;
    Hessian(SymPosDefMatrix<double> value, var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    Hessian(SymPosDefMatrix<double> value, var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(Hessian) {
        return "Hessian";
    }
};

class FIM : public var::Constant<FIM, parameter_spd_payload> {
    using base_type = var::Constant<FIM, parameter_spd_payload>;

   public:
    using base_type::base_type;
    FIM() = default;
    FIM(SymPosDefMatrix<double> value, var::Parameters_transformed const& params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}
    FIM(SymPosDefMatrix<double> value, var::Parameters_transformed const* params)
        : base_type(parameter_spd_payload(std::move(value), params)) {}

    friend std::string className(FIM) {
        return "FIM";
    }
    Maybe_error<bool> is_good() const {
        if (!std::isfinite(var::fullsum((*this)()))) {
            std::stringstream ss;
            ss << "FIM not finite: " << (*this)();
            return error_message(ss.str());
        }
        for (std::size_t i = 0; i < (*this)().nrows(); ++i)
            if ((*this)()(i, i) <= 0) {
                std::stringstream ss;
                ss << "FIM diagonal negative or cero: " << (*this)();
                return error_message(ss.str());
            }

        return true;
    }
};

class Hessian_minus_CovGrad
    : public var::Constant<Hessian_minus_CovGrad, parameter_symmetric_payload> {
    using base_type = var::Constant<Hessian_minus_CovGrad, parameter_symmetric_payload>;

   public:
    using base_type::base_type;
    Hessian_minus_CovGrad() = default;
    Hessian_minus_CovGrad(SymmetricMatrix<double> value,
                          var::Parameters_transformed const& params)
        : base_type(parameter_symmetric_payload(std::move(value), params)) {}
    Hessian_minus_CovGrad(SymmetricMatrix<double> value,
                          var::Parameters_transformed const* params)
        : base_type(parameter_symmetric_payload(std::move(value), params)) {}

    friend std::string className(Hessian_minus_CovGrad) {
        return "Hessian_minus_CovGrad";
    }
};



using logLs = var::Vector_Space<logL, elogL, vlogL>;

inline logLs logLs_0() {
    return logLs(logL(0.0), elogL(0.0), vlogL(0.0));
}

inline logLs nan_logL() {
    return logLs(logL(std::numeric_limits<double>::quiet_NaN()),
                 elogL(std::numeric_limits<double>::quiet_NaN()),
                 vlogL(std::numeric_limits<double>::quiet_NaN()));
}

using dMacro_State_Hessian_minimal = var::Vector_Space<logL, elogL, vlogL, Grad, FIM>;
using dlogPs = var::Vector_Space<logL, Grad, FIM>;

template <class... Vs>
inline Maybe_error<bool> is_good(const var::Vector_Space<Vs...>& x) {
    return (get<Vs>(x).is_good() && ...);
}

#endif  // DISTRIBUTIONS_H
