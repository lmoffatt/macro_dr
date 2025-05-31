#ifndef BAYESIAN_LINEAR_REGRESSION_H
#define BAYESIAN_LINEAR_REGRESSION_H

#include "distributions.h"
#include "matrix.h"
#include "matrix_random.h"
#include "multivariate_normal_distribution.h"
#include "variables.h"
//#include "parallel_tempering.h"
#include <cmath>

inline double digamma(double x) {
    if (x > 2)

        return std::log(x) - 0.5 / x - (1. / 12.) * std::pow(x, -2) +
               (1. / 120.) * std::pow(x, -4) - (1. / 252.) * std::pow(x, -6) +
               (1. / 240.) * std::pow(x, -8) - (5. / 660.) * std::pow(x, -10) +
               (691. / 32760.) * std::pow(x, -12) - (1. / 12.) * std::pow(x, -14);
    else {
        double eps = 1e-6;
        double xp = x * (1 + eps);
        return (var::lgamma(xp) - var::lgamma(x)) / (xp - x);
    }
}

struct linear_model {};

class covariance_model {
    std::size_t ndim_;
    Matrix<double> sigmas_;

   public:
    covariance_model(std::size_t ndim, Matrix<double> const& sigmas)
        : ndim_{ndim}, sigmas_{sigmas} {
    }
    auto operator()(mt_64i& mt) {
        auto X = sample(mt, normal_distribution{0, 1}, ndim_, ndim_);
        auto [Q, R] = qr(X);
        auto s = DiagPosDetMatrix(sigmas_);
        auto cov = AT_D_A(Q, s);
        return cov;
    }
};

class independent_variable_model {
    std::size_t npar_ = 40ul;
    double log10_std_par_ = 1.0;

    double mean_mean_par_ = 0.0;
    double std_mean_par_ = 10.0;

   public:
    independent_variable_model(std::size_t npar = 40ul, double std_log10_std_par = 1.0,

                               double mean_mean_par = 0.0, double std_mean_par = 10.0)
        : npar_{npar},
          log10_std_par_{std_log10_std_par},
          mean_mean_par_{mean_mean_par},
          std_mean_par_{std_mean_par} {
    }

    auto operator()(mt_64i& mt, std::size_t nsamples) {
        auto v_par_stds = sample(mt, normal_distribution(0, log10_std_par_), 1, npar_);
        auto par_stds = apply([](auto const& x) { return std::pow(10, x); }, v_par_stds);

        auto cov_par = covariance_model(npar_, par_stds)(mt);
        auto mean_par = sample(mt, normal_distribution(mean_mean_par_, std_mean_par_), 1, npar_);

        auto par_dist_ = make_multivariate_normal_distribution(mean_par, cov_par);

        auto par_dist = std::move(par_dist_.value());
        return std::tuple(par_dist(mt, nsamples), mean_par, cov_par);
    }
};

inline auto get_parameters(linear_model, const Matrix<double>& beta) {
    double logvar = beta[0];
    auto b = Matrix<double>(1, size(beta) - 1, false);
    for (std::size_t i = 0; i < b.ncols(); ++i) b[i] = beta[i + 1];
    return std::tuple(logvar, b);
}

template <class FunctionTable>
Maybe_error<var::Vector_Space<logL, elogL, vlogL>> logLikelihood(FunctionTable&& f, linear_model,
                                                                 const Matrix<double>& beta,
                                                                 const Matrix<double>& y,
                                                                 const Matrix<double>& X) {
    assert(beta.ncols() - 1 == X.ncols() && "beta has the right number");
    assert(y.nrows() == X.nrows() && "right number of rows in X");

    auto [logvar, b] = get_parameters(linear_model{}, beta);
    double var = std::exp(logvar);
    auto yfit = X * tr(b);
    auto ydiff = y - yfit;
    double n = y.nrows();
    double SS = xtx(ydiff);
    double chi2 = SS / var;
    double out = -0.5 * (n * log(2 * std::numbers::pi) + n * logvar + chi2);
    double t_elogL = -0.5 * (n * log(2 * std::numbers::pi) + n * logvar + 1);
    double t_vlogL = 0.5;
    if (std::isfinite(out))
        return var::Vector_Space(logL(out), elogL(t_elogL), vlogL(t_vlogL));
    else {
        std::cerr << std::string("likelihood error: ") + std::to_string(out) << "\n";
        return error_message(std::string("likelihood error: ") + std::to_string(out));
    }
}

inline auto simulate(mt_64i& mt, linear_model, const Matrix<double>& beta,
                     const Matrix<double>& X) {
    assert(beta.ncols() - 1 == X.ncols() && "beta has the right number");
    auto [logvar, b] = get_parameters(linear_model{}, beta);
    double var = std::exp(logvar);
    auto yfit = X * tr(b);
    auto noise = random_matrix_normal(mt, yfit.nrows(), yfit.ncols(), 0, std::sqrt(var));
    return yfit + noise;
}

template <class Cov>
class bayesian_linear_model : public linear_model,
                              public multivariate_gamma_normal_distribution<double, Cov> {
   public:
    using dist_type = multivariate_gamma_normal_distribution<double, Cov>;

    using dist_type::operator();
    using dist_type::logP;
    using dist_type::size;
    bayesian_linear_model(dist_type&& d) : linear_model{}, dist_type{std::move(d)} {
    }

    dist_type const& prior() const {
        return *this;
    }
    linear_model const& likelihood() const {
        return *this;
    }
};

template <class Cov>
    requires(Covariance<double, Cov>)
Maybe_error<bayesian_linear_model<Cov>> make_bayesian_linear_model(double prior_eps_df,
                                                                   double prior_eps_variance,
                                                                   Matrix<double>&& mean,
                                                                   Cov&& cov) {
    auto a = prior_eps_df / 2.0;
    auto b = prior_eps_df * prior_eps_variance / 2.0;
    auto prior = make_multivariate_normal_distribution(std::move(mean), std::move(cov));
    if (prior)
        return bayesian_linear_model<Cov>(multivariate_gamma_normal_distribution(
            log_inverse_gamma_distribution(a, b), std::move(prior.value())));
    else
        return error_message(prior.error()() + "\n in make_bayesian_linear_model");
}

inline auto make_bayesian_linear_model(double prior_eps_df, double prior_eps_variance,
                                       std::size_t npar, double mean_b, double std_b,
                                       double prior_error_factor) {
    return make_bayesian_linear_model(prior_eps_df, prior_eps_variance,
                                      Matrix<double>(1ul, npar, mean_b),
                                      IdM<double>(npar) * std_b * std_b * prior_error_factor);
}

struct conjugate {};

template <class Cova>
    requires Covariance<double, Cova>
auto posterior(conjugate, const multivariate_gamma_normal_distribution<double, Cova>& prior,
               linear_model, const Matrix<double>& y, const Matrix<double>& X) {
    auto a_0 = prior.alpha();
    ;
    auto prior_eps_df = 2.0 * a_0;
    auto b_0 = prior.beta();
    auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto beta_0 = prior.mean();
    auto L_n = L_0 + SSx;

    auto beta_n = tr(inv(L_n) * ((tr(X) * y) + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = xtx(ydiff.value());

    auto a_n = a_0 + n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

    auto posterior_Normal =
        make_multivariate_normal_distribution_from_precision(std::move(beta_n), std::move(L_n));

    return multivariate_gamma_normal_distribution<double, SymPosDefMatrix<double>>(
        log_inverse_gamma_distribution(a_n, b_n.value()), std::move(posterior_Normal.value()));
}

template <class Cova>
    requires Covariance<double, Cova>
auto evidence(conjugate, const multivariate_gamma_normal_distribution<double, Cova>& prior,
              const linear_model& li, const Matrix<double>& y, const Matrix<double>& X) {
    auto post = posterior(conjugate{}, prior, li, y, X);
    auto a_0 = prior.alpha();
    ;
    auto b_0 = prior.beta();
    auto L_0 = prior.Gamma();

    auto a_n = post.alpha();
    ;
    auto b_n = post.beta();
    auto L_n = post.Gamma();
    auto n = y.nrows();

    auto logE_n = -0.5 * n * std::log(2 * std::numbers::pi) + 0.5 * (logdet(L_0) - logdet(L_n)) +
                  a_0 * log(b_0) - a_n * log(b_n) + var::lgamma(a_n) - var::lgamma(a_0);
    return logE_n;
}

//template <class Cova>
//  requires Covariance<double, Cova>
//auto bayesian_linear_regression(
//    const multivariate_normal_distribution<double, Cova> &prior,
//    double prior_eps_df, double prior_eps_variance, const Matrix<double> &y,
//    const Matrix<double> &X) {
//  auto L_0 = prior.cov_inv() * prior_eps_variance;
//  auto SSx = XTX(X);
//  auto n = y.nrows();
//  auto a_0 = prior_eps_df / 2.0;
//  auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
//  auto beta_0 = prior.mean();
//  // auto beta_ml=inv(SSx)*tr(X)*y;
//  auto L_n = SSx + L_0;

//  auto beta_ = inv(SSx) * tr(X) * y;
//  auto beta_nn = tr(inv(L_n) * (SSx * beta_ + (L_0 * tr(prior.mean()))));
//  auto beta_n = tr(inv(L_n) * (tr(X) * y + (L_0 * tr(prior.mean()))));

//  std::cerr << "beta_n<<beta_nn<<beta_n-beta_nn" << beta_n << beta_nn
//            << beta_n - beta_nn;
//  auto yfit = X * tr(beta_n);
//  auto ydiff = y - yfit;
//  auto SS = xtx(ydiff.value());

//  auto posterior =
//      make_multivariate_normal_distribution_from_precision(beta_n, L_n);
//  auto a_n = a_0 + n / 2.0;
//  auto b_n = b_0 + 0.5 * (xtx(y) + xAxt(beta_0, L_0) - xAxt(beta_n, L_n));

//  auto b_n2 = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
//  // auto
//  // E_n=std::pow(2*std::numbers::pi,-n/2)*std::sqrt(det(L_0)/det(L_n))*std::pow(b_0,a_0)/std::pow(b_n,a_n)*std::tgamma(a_n)/std::tgamma(a_0);

//  std::cerr << "logdet(L_n)" << logdet(L_n);

//  auto logE_n = -0.5 * n * std::log(2 * std::numbers::pi) +
//                0.5 * (logdet(L_0) - logdet(L_n)) + a_0 * log(b_0) -
//                a_n * log(b_n) + var::lgamma(a_n) - var::lgamma(a_0);
//  return std::tuple(
//      logE_n, "beta_n", beta_n, "b_n", b_n, "b_n2", b_n2,

//      "xtx(y)", xtx(y),
//      "0.5*( logdet(L_0)) + a_0 * log(b_0) - var::lgamma(a_0)",
//      0.5 * (logdet(L_0)) + a_0 * log(b_0) - var::lgamma(a_0),
//      "0.5*(-logdet(L_n)) - a_n * log(b_n) + var::lgamma(a_n)",
//      0.5 * (-logdet(L_n)) - a_n * log(b_n) + var::lgamma(a_n),
//      "-0.5* n*std::log(2*std::numbers::pi)",
//      -0.5 * n * std::log(2 * std::numbers::pi),
//      "0.5*( logdet(L_0) -logdet(L_n))", 0.5 * (logdet(L_0) - logdet(L_n)),
//      "a_0 * log(b_0) - a_n * log(b_n)", a_0 * log(b_0) - a_n * log(b_n),
//      "var::lgamma(a_n) - var::lgamma(a_0)",
//      var::lgamma(a_n) - var::lgamma(a_0),

//      "sqrt(b_n / (a_n - 1))", sqrt(b_n / (a_n - 1)), "b_n / (a_n - 1)",
//      b_n / (a_n - 1));
//}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_Evidence(
    const multivariate_normal_distribution<double, Cova>& prior, const linear_model&,
    double prior_eps_df, double prior_eps_variance, const Matrix<double>& y,
    const Matrix<double>& X) {
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    auto L_n = SSx + L_0;
    auto beta_n = tr(inv(L_n) * (tr(X) * y + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = xtx(ydiff.value());
    auto a_n = a_0 + n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

    auto logE_n = -0.5 * n * std::log(2 * std::numbers::pi) + 0.5 * (logdet(L_0) - logdet(L_n)) +
                  a_0 * log(b_0) - a_n * log(b_n) + var::lgamma(a_n) - var::lgamma(a_0);
    return logE_n;
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_Evidence(
    const multivariate_gamma_normal_distribution<double, Cova>& prior, const linear_model&,
    const Matrix<double>& y, const Matrix<double>& X) {
    auto a_0 = prior.alpha();
    ;
    auto b_0 = prior.beta();
    auto L_0 = prior.Gamma();
    auto prior_eps_df = 2.0 * a_0;
    auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

    auto SSx = XTX(X);
    auto n = y.nrows();
    auto beta_0 = prior.mean();
    auto L_n = SSx + L_0;
    auto beta_n = tr(inv(L_n) * (tr(X) * y + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = xtx(ydiff.value());
    auto a_n = a_0 + n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

    auto logE_n = -0.5 * n * std::log(2 * std::numbers::pi) + 0.5 * (logdet(L_0) - logdet(L_n)) +
                  a_0 * log(b_0) - a_n * log(b_n) + var::lgamma(a_n) - var::lgamma(a_0);
    return logE_n;
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_posterior(
    const multivariate_normal_distribution<double, Cova>& prior, double prior_eps_df,
    double prior_eps_variance, const Matrix<double>& y, const Matrix<double>& X, double beta0) {
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    auto L_n = L_0 + beta0 * SSx;

    auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta0 * xtx(ydiff.value());

    auto a_n = a_0 + beta0 * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
    return std::tuple(std::log(b_n.value() / a_n), beta_n.value());
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_posterior(
    const multivariate_gamma_normal_distribution<double, Cova>& prior, const Matrix<double>& y,
    const Matrix<double>& X) {
    auto a_0 = prior.alpha();
    ;
    auto prior_eps_df = 2.0 * a_0;
    auto b_0 = prior.beta();
    auto prior_eps_variance = 2.0 * b_0 / prior_eps_df;

    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto beta_0 = prior.mean();
    auto L_n = L_0 + SSx;

    auto beta_n = tr(inv(L_n) * ((tr(X) * y) + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = xtx(ydiff.value());

    auto a_n = a_0 + n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

    auto posterior_Normal =
        make_multivariate_normal_distribution_from_precision(std::move(beta_n), std::move(L_n));

    return multivariate_gamma_normal_distribution<double, SymPosDefMatrix<double>>(
        log_inverse_gamma_distribution(a_n, b_n.value()), std::move(posterior_Normal.value()));
}

template <class Cova>
    requires Covariance<double, Cova>
auto mean_logLik(conjugate, const multivariate_gamma_normal_distribution<double, Cova>& prior,
                 const linear_model&, const Matrix<double>& y, const Matrix<double>& X,
                 double beta0) {
    auto a_0 = prior.alpha();
    ;
    auto b_0 = prior.beta();
    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto beta_0 = prior.mean();
    auto L_n = L_0 + beta0 * SSx;
    auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));
    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta0 * xtx(ydiff.value());

    auto a_n = a_0 + beta0 * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
    double d_a_n = 1.0 * n / 2.0;
    auto d_b_n = 0.5 * xtx(ydiff.value());
    auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                      a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;
    return mean_logLi;
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_gamma_normal_distribution<double, Cova>& prior, const linear_model&,
    const Matrix<double>& y, const Matrix<double>& X, double beta0) {
    auto a_0 = prior.alpha();
    ;
    auto b_0 = prior.beta();
    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto beta_0 = prior.mean();
    auto L_n = L_0 + beta0 * SSx;
    auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));
    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta0 * xtx(ydiff.value());

    auto a_n = a_0 + beta0 * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
    double d_a_n = 1.0 * n / 2.0;
    auto d_b_n = 0.5 * xtx(ydiff.value());
    auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                      a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;
    return mean_logLi;
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_var_logLik_diff(
    const multivariate_gamma_normal_distribution<double, Cova>& prior, const Matrix<double>& y,
    const Matrix<double>& X, double beta0, double eps) {
    double beta1 = std::max(eps, beta0 * (1 + eps));
    return (bayesian_linear_regression_calculate_mean_logLik(prior, y, X, beta1) -
            bayesian_linear_regression_calculate_mean_logLik(prior, y, X, beta0)) /
           (beta1 - beta0);
}

template <class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_normal_distribution<double, Cova>& prior, double prior_eps_df,
    double prior_eps_variance, const Matrix<double>& y, const Matrix<double>& X, double beta0) {
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    auto L_n = L_0 + beta0 * SSx;
    auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));
    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta0 * xtx(ydiff.value());

    auto a_n = a_0 + beta0 * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
    double d_a_n = 1.0 * n / 2.0;
    auto d_b_n = 0.5 * xtx(ydiff.value());
    auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                      a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;
    return mean_logLi;
}

template <class Parameters, class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik_posterior(
    const multivariate_normal_distribution<double, Cova>& prior, double prior_eps_df,
    double prior_eps_variance, const Matrix<double>& y, const Matrix<double>& X, double beta0) {
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    std::tuple<Maybe_error<double>, double, Parameters> mean_logLik;

    auto L_n = L_0 + beta0 * SSx;

    auto beta_n = tr(inv(L_n) * (beta0 * (tr(X) * y) + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta0 * xtx(ydiff.value());

    auto a_n = a_0 + beta0 * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);

    auto logE_n = -0.5 * beta0 * n * std::log(2 * std::numbers::pi) +
                  0.5 * (logdet(L_0) - logdet(L_n)) + a_0 * log(b_0) - a_n * log(b_n) +
                  var::lgamma(a_n) - var::lgamma(a_0);

    double d_a_n = 1.0 * n / 2.0;
    auto d_b_n = 0.5 * xtx(ydiff.value());

    auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) - 0.5 * Trace(inv(L_n) * SSx) -
                      a_n / b_n * d_b_n + (digamma(a_n) - log(b_n)) * d_a_n;

    return std::tuple(mean_logLi, std::log(b_n.value() / a_n), beta_n.value());
}

#endif  // BAYESIAN_LINEAR_REGRESSION_H
