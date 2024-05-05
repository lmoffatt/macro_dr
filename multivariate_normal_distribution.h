#ifndef MULTIVARIATE_NORMAL_DISTRIBUTION_H
#define MULTIVARIATE_NORMAL_DISTRIBUTION_H

#include "general_algorithm_on_containers.h"
#include "lgamma.h"
#include "matrix.h"
#include "maybe_error.h"
#include "random_samplers.h"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

template <typename T, typename C>
concept Covariance = std::convertible_to<C, DiagPosDetMatrix<T>> ||
                     std::convertible_to<C, SymPosDefMatrix<T>>;

class normal_distribution {
private:
  std::normal_distribution<double> n_;

public:
  normal_distribution(double mean, double stddev) : n_{mean, stddev} {}
  normal_distribution() : n_{} {}
  auto mean() const { return n_.mean(); }
  std::size_t size() const { return 1ul; }
  auto cov() const { return n_.stddev() * n_.stddev(); }
  auto stddev() const { return n_.stddev(); }

  double operator()(mt_64i &mt) { return n_(mt); }

  Maybe_error<double> logP(double x) const {
    double out = -0.5 * (log(2 * std::numbers::pi * cov()) +
                         std::pow(x - mean(), 2) / cov());
    if (std::isfinite(out))
      return out;
    else
      return error_message("likelihood not finite:" + std::to_string(out));
  }
};

template <class Parameters> class save_Parameter;
template <typename T, class Cova>
  requires Covariance<T, Cova>
class multivariate_normal_distribution {
private:
  std::normal_distribution<T> n_;
  Matrix<T> mean_;
  Cova cov_;
  using cholesky_type = std::decay_t<decltype(get_value(::cholesky(cov_)))>;
  cholesky_type cho_;
  Cova cov_inv_;
  double logdetCov_;

  static double calc_logdet(const cholesky_type &cho) {
    return std::visit([](auto const &m) { return 2 * logdet(m); }, cho);
  }

  
public:
  multivariate_normal_distribution(Matrix<T> &&mean, Cova &&cov,
                                   cholesky_type &&chol, Cova &&cov_inv,
                                   double logdetCov)
      : n_{}, mean_{std::move(mean)}, cov_{std::move(cov)},
      cho_{std::move(chol)}, cov_inv_{std::move(cov_inv)},
      logdetCov_{logdetCov} {}
  multivariate_normal_distribution() = default;
  template <class Mat, class Cov, typename Te,
            //= std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
            class Covae>
  // = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
    requires(Covariance<Te, Covae> && contains_value<Mat &&, Matrix<Te>> &&
             contains_value<Cov &&, Covae>)
  friend Maybe_error<multivariate_normal_distribution<Te, Covae>>
  make_multivariate_normal_distribution(Mat &&mean, Cov &&cov);

  template <
      class Mat, class Cov,
      typename Te, //= std::decay_t<decltype(get_value(std::declval<Mat>())(0,
                   // 0))>,
      class Covae> // = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
    requires(Covariance<Te, Covae> && contains_value<Mat &&, Matrix<Te>> &&
             contains_value<Cov &&, Covae>)
  friend Maybe_error<multivariate_normal_distribution<Te, Covae>>
  make_multivariate_normal_distribution_from_precision(Mat &&mean,
                                                       Cov &&cov_inv);

  auto &mean() const { return mean_; }

  std::size_t size() const { return mean().size(); }
  auto &cov() const { return cov_; }

  auto &cov_inv() const { return cov_inv_; }
  auto &cholesky() const { return cho_; }

  Matrix<double> operator()(mt_64i &mt) const {
    auto z =
        sample(mt, normal_distribution(0, 1), mean().nrows(), mean().ncols());
    if (mean().nrows() == 1)
      return z * tr(cholesky()) + mean();
    else
      return cholesky() * z + mean();
  }
  auto operator()(mt_64i &mt, std::size_t n) const {
    return sample(mt, normal_distribution(0.0, 1.0), n, mean().ncols()) *
           tr(cholesky());
  }

  auto logDetCov() const { return logdetCov_; }
  double chi2(const Matrix<T> &x) const {
    auto xdiff = x - mean();
    if (xdiff.nrows() == cov_inv().nrows())
      return xtAx(xdiff, cov_inv());
    else
      return xAxt(xdiff, cov_inv());
  }
  Maybe_error<double> logP(const Matrix<T> &x) const {
    assert(x.size() == mean().size());
    double out = -0.5 * (mean().size() * log(2 * std::numbers::pi) +
                         logDetCov() + chi2(x));
    if (std::isfinite(out))
      return out;
    else
      return error_message("likelihood not finite:" + std::to_string(out));
  }
  
  Matrix<T> score(const Matrix<T> &x) const {
      assert(x.size() == mean().size());
      auto xdiff =  mean()-x;
      if (xdiff.nrows() == cov_inv().nrows())
          return cov_inv()*xdiff;
      else
          return xdiff*cov_inv();
   }
  friend std::ostream &operator<<(std::ostream &os,
                                  const multivariate_normal_distribution &m) {

    os << "mean " << m.mean();
    os << "diag(cov) " << diag(m.cov());
    return os;
  }
};

template <class Mat, class Cov,
          typename T =
              std::decay_t<decltype(get_value(std::declval<Mat>())(0ul, 0ul))>,
          class Cova = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
  requires(Covariance<T, Cova> && contains_value<Mat &&, Matrix<T>> &&
           contains_value<Cov &&, Cova>)
Maybe_error<multivariate_normal_distribution<T, Cova>>
make_multivariate_normal_distribution(Mat &&mean, Cov &&cov) {
  return_error<multivariate_normal_distribution<T, Cova>> Error{
      "make_multivariate_normal_distribution"};
  if (!is_valid(mean))
    return Error(get_error(mean)());
  else if (!is_valid(cov))
    return Error(get_error(cov)());
  else {
    auto beta_cov = get_value(std::forward<Cov>(cov));
    auto chol = cholesky(beta_cov);

    if (chol) {
      auto inv = inv_from_chol(chol.value());
      if (inv) {
        auto meanbeta = get_value(std::forward<Mat>(mean));
        auto logDetCov = logdet(chol.value());
        if (logDetCov)
          return multivariate_normal_distribution<T, Cova>(
              std::move(meanbeta), std::move(beta_cov), std::move(chol.value()),
              std::move(inv.value()), logDetCov.value());
        else
          return Error(logDetCov.error()() + " log determinant fails");
      } else
        return Error(inv.error()() + " covariance cannot be inverted");
    } else
      return Error(chol.error()() +
                   " cholesky fails to build a normal distribution");
  }
}

template <
    class Mat, class Cov,
    typename T = std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
    class Cova = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
  requires(Covariance<T, Cova> && contains_value<Mat &&, Matrix<T>> &&
           contains_value<Cov &&, Cova>)
Maybe_error<multivariate_normal_distribution<T, Cova>>
make_multivariate_normal_distribution_from_precision(Mat &&mean,
                                                     Cov &&cov_inv) {
  return_error<multivariate_normal_distribution<T, Cova>> Error{
      "make_multivariate_normal_distribution"};
  if (!is_valid(mean))
    return Error(get_error(mean)());
  else if (!is_valid(cov_inv))
    return Error(get_error(cov_inv)());
  else {
    auto chol_inv = cholesky(cov_inv);

    if (chol_inv) {
      auto chol = inv(chol_inv.value());

      if (chol) {
        auto cov = XXT(tr(chol.value()));
        auto beta_mean = get_value(std::forward<Mat>(mean));
        auto beta_cov_inv = get_value(std::forward<Cov>(cov_inv));
        auto logDetCov = logdet(diag(chol.value()));
        if (logDetCov)
          return multivariate_normal_distribution<T, Cova>(
              std::move(beta_mean), std::move(cov), std::move(chol.value()),
              std::move(beta_cov_inv), logDetCov.value());
        else
          return Error(logDetCov.error()() + " log of determinant fails");

      } else
        return Error(chol.error()() + " cholesky cannot be inverted");
    } else
      return Error(chol_inv.error()() +
                   " cholesky fails to be built from precision");
  }
}

template <typename T, class Cova>
  requires Covariance<T, Cova>
class multivariate_normal_distribution_by_eig {
private:
  std::normal_distribution<T> n_;
  Matrix<T> mean_;
  Cova cov_;
  Matrix<T> m_Q;
  DiagonalMatrix<T> m_std;
  Cova cov_inv_;
  double logdetCov_;

  static double calc_logdet(const DiagonalMatrix<T> &t_var, double min_landa) {
    double out = 0;
    for (std::size_t i = 0; i < t_var.size(); ++i) {
      if (t_var[i] > min_landa)
        out = out + std::log(t_var[i]);
    }
    return out;
  }

  static double calc_std(const DiagonalMatrix<T> &t_var, double min_landa) {

    DiagonalMatrix<T> t_std(t_var.nrows(), t_var.ncols());
    for (std::size_t i = 0; i < t_var.size(); ++i) {
      if (t_var[i] > min_landa)
        t_std[i] = std::sqrt(t_var[i]);
    }
    return t_std;
  }

  static double calc_cov_inv(const Matrix<T> &Q, const DiagonalMatrix<T> &t_var,
                             double min_landa) {

    DiagonalMatrix<T> v_inv(t_var.nrows(), t_var.ncols());
    for (std::size_t i = 0; i < t_var.size(); ++i) {
      if (t_var[i] > min_landa)
        v_inv[i] = 1 / t_var[i];
    }
    return Q * v_inv * tr(Q);
  }

  multivariate_normal_distribution_by_eig(Matrix<T> &&t_mean, Cova &&t_cov,
                                          Matrix<T> &&Q,
                                          DiagonalMatrix<T> &&std,
                                          Cova &&cov_inv, double logdetCov)
      : n_{}, mean_{std::move(t_mean)}, cov_{std::move(t_cov)},
        m_Q{std::move(Q)}, m_std{std::move(std)}, cov_inv_{std::move(cov_inv)},
        logdetCov_{logdetCov} {}

public:
  template <class Mat, class Cov, typename Te,
            //= std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
            class Covae>
  // = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
    requires(Covariance<Te, Covae> && contains_value<Mat &&, Matrix<Te>> &&
             contains_value<Cov &&, Covae>)
  friend Maybe_error<multivariate_normal_distribution_by_eig<Te, Covae>>
  make_multivariate_normal_distribution_by_eig(Mat &&mean, Cov &&cov);

  auto &mean() const { return mean_; }

  std::size_t size() const { return mean().size(); }
  auto &cov() const { return cov_; }

  auto &cov_inv() const { return cov_inv_; }

  Matrix<double> operator()(mt_64i &mt) {
    auto z =
        sample(mt, normal_distribution(0, 1), mean().nrows(), mean().ncols());
    if (mean().nrows() == 1)
      return z * tr(m_Q * m_std) + mean();
    else
      return m_Q * m_std * z + mean();
  }
  auto operator()(mt_64i &mt, std::size_t n) {
    return sample(mt, normal_distribution(0.0, 1.0), n, mean().ncols()) *
           tr(m_Q * m_std);
  }

  auto logDetCov() const { return logdetCov_; }
  double chi2(const Matrix<T> &x) const {
    auto xdiff = x - mean();
    if (xdiff.nrows() == cov_inv().nrows())
      return xtAx(xdiff, cov_inv());
    else
      return xAxt(xdiff, cov_inv());
  }
  Maybe_error<double> logP(const Matrix<T> &x) const {
    assert(x.size() == mean().size());
    double out = -0.5 * (mean().size() * log(2 * std::numbers::pi) +
                         logDetCov() + chi2(x));
    if (std::isfinite(out))
      return out;
    else
      return error_message("likelihood not finite:" + std::to_string(out));
  }

  friend std::ostream &
  operator<<(std::ostream &os,
             const multivariate_normal_distribution_by_eig &m) {

    os << "mean " << m.mean();
    os << "diag(cov) " << diag(m.cov());
    return os;
  }
};

template <class Mat, class Cov,
          typename T =
              std::decay_t<decltype(get_value(std::declval<Mat>())(0ul, 0ul))>,
          class Cova = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
  requires(Covariance<T, Cova> && contains_value<Mat &&, Matrix<T>> &&
           contains_value<Cov &&, Cova>)
Maybe_error<multivariate_normal_distribution_by_eig<T, Cova>>
make_multivariate_normal_distribution_by_eig(Mat &&mean, Cov &&cov,
                                             double min_std) {
  if (!is_valid(mean))
    return get_error(mean);
  else if (!is_valid(cov))
    return get_error(cov);
  else {
    auto beta_cov = get_value(std::forward<Cov>(cov));
    auto Maybe_eig = eigs(beta_cov);

    if (!Maybe_eig)
      return Maybe_eig.error();
    else {
      auto [r_Q, r_var] = std::move(Maybe_eig.value());
      auto meanbeta = get_value(std::forward<Mat>(mean));
      auto logDetCov = calc_logdet(r_var);
      auto v_Cov_inv = calc_cov_inv(r_Q, r_var);
      auto v_std = calc_std(r_var, min_std);
      return multivariate_normal_distribution_by_eig<T, Cova>(
          std::move(meanbeta), std::move(beta_cov), std::move(r_Q),
          std::move(v_std), std::move(v_Cov_inv), logDetCov);
    }
  }
}

class multivariate_normal_distribution_of_probabilities
    : public multivariate_normal_distribution<double, SymPosDefMatrix<double>> {
private:
  std::size_t m_excluded_row;

  static Maybe_error<Matrix<double>> expand(const Matrix<double> &s,
                                            std::size_t excluded_row) {
    Matrix<double> out;
    if (s.nrows() == s.size()) {
      out = Matrix<double>(s.nrows() + 1, 1);
    } else {
      out = Matrix<double>(1, s.ncols() + 1);
    }
    for (std::size_t i = 0; i < out.size(); ++i) {
      if (i < excluded_row) {
        if ((s[i] < 0) || (s[i] > 1))
          return error_message("");
        else
          out[i] = s[i];
      } else if (i > excluded_row) {
        if ((s[i - 1] < 0) || (s[i - 1] > 1))
          return error_message("");
        else
          out[i] = s[i - 1];
      } else {
        if (var::sum(s) > 1)
          return error_message("");
        else
          out[i] = 1.0 - var::sum(s);
      }
    }
    return out;
  }
  static auto reduce(const Matrix<double> &s, std::size_t excluded_row) {
    Matrix<double> out;
    if (s.nrows() == s.size()) {
      out = Matrix<double>(s.nrows() - 1, 1);
    } else {
      out = Matrix<double>(1, s.ncols() - 1);
    }
    for (std::size_t i = 0; i < out.size(); ++i) {
      if (i < excluded_row) {
        out[i] = s[i];
      } else if (i > excluded_row) {
        out[i - 1] = s[i];
      }
    }
    return out;
  }

  static std::size_t reduce(std::size_t i, std::size_t excluded_row) {
    if (i < excluded_row)
      return i;
    else
      return i - 1;
  }
  static auto reduce_cov(const SymPosDefMatrix<double> s,
                         std::size_t excluded_row) {
    SymPosDefMatrix<double> out(s.nrows() - 1, s.ncols() - 1);
    for (std::size_t i = 0; i < s.nrows(); ++i) {
      if (i != excluded_row) {
        for (std::size_t j = 0; j < s.ncols(); ++j) {
          if (j != excluded_row) {
            out.set(reduce(i, excluded_row), reduce(j, excluded_row), s(i, j));
          }
        }
      }
    }
    return out;
  }

  using base_type =
      multivariate_normal_distribution<double, SymPosDefMatrix<double>>;

public:
  multivariate_normal_distribution_of_probabilities(base_type &&d,
                                                    std::size_t excluded_row)
      : base_type{std::move(d)}, m_excluded_row{excluded_row} {}
  multivariate_normal_distribution_of_probabilities() = delete;
  friend Maybe_error<multivariate_normal_distribution_of_probabilities>
  make_multivariate_normal_distribution_of_probabilities(
      Matrix<double> mean, SymPosDefMatrix<double> cov);

  friend Maybe_error<multivariate_normal_distribution_of_probabilities>
  make_multivariate_normal_distribution_of_probabilities(
      Matrix<double> const &mean, SymPosDefMatrix<double> const &cov);
  std::size_t size() const { return base_type::mean().size() + 1; }

  static auto determine_excluded_row(SymPosDefMatrix<double> cov) {
    std::size_t excluded_row = 0;
    for (std::size_t i = 0; i < cov.nrows(); ++i) {
      if (cov(i, i) > cov(excluded_row, excluded_row))
        excluded_row = i;
    }
    return excluded_row;
  }

  static auto logDet(SymPosDefMatrix<double> cov) {
    auto excluded_row = determine_excluded_row(cov);
    return logdet(reduce_cov(cov, excluded_row));
  }

  Matrix<double> operator()(mt_64i &mt) {
    Maybe_error<Matrix<double>> out(error_message(""));
    while (!out) {
      auto s = base_type::operator()(mt);
      out = expand(s, m_excluded_row);
    }
    return std::move(out.value());
  }

  auto logDetCov() const { return base_type::logDetCov(); }
  double chi2(const Matrix<double> &x) const {
    return base_type::chi2(reduce(x, m_excluded_row));
  }
  Maybe_error<double> logP(const Matrix<double> &x) const {
    assert(size() == x.size());
    return base_type::logP(reduce(x, m_excluded_row));
  }
};

inline Maybe_error<multivariate_normal_distribution_of_probabilities>
make_multivariate_normal_distribution_of_probabilities(
    Matrix<double> mean, SymPosDefMatrix<double> cov) {
  using T = multivariate_normal_distribution_of_probabilities;
  std::size_t excluded_row = T::determine_excluded_row(cov);
  auto Maybe_dist = make_multivariate_normal_distribution(
      T::reduce(std::move(mean), excluded_row),
      T::reduce_cov(std::move(cov), excluded_row));
  if (!Maybe_dist)
    return Maybe_dist.error();
  else
    return T(std::move(Maybe_dist.value()), excluded_row);
}

class inverse_gamma_distribution {
private:
  std::gamma_distribution<> g_;
  double cte_int_;

public:
  inverse_gamma_distribution(double _alpha, double _beta)
      : g_{_alpha, _beta},
        cte_int_{-var::lgamma(_alpha) + _alpha * std::log(_beta)} {}

  double operator()(mt_64i &mt) { return 1.0 / g_(mt); }

  double alpha() const { return g_.alpha(); }
  double beta() const { return g_.beta(); }
  Maybe_error<double> logP(double x) const {
    auto out = cte_int_ - (alpha() + 1.0) * std::log(x) - beta() / x;
    if (std::isfinite(out))
      return out;
    else
      return error_message("probability not finite:" + std::to_string(out));
  }
};

class log_inverse_gamma_distribution {
private:
  std::gamma_distribution<> g_;
  double cte_int_;

public:
  log_inverse_gamma_distribution(double _alpha, double _beta)
      : g_{_alpha, 1.0 / _beta},
        cte_int_{-var::lgamma(_alpha) + _alpha * std::log(_beta)} {}

  double operator()(mt_64i &mt) { return -std::log(g_(mt)); }

  double alpha() const { return g_.alpha(); }
  double beta() const { return 1.0 / g_.beta(); }
  Maybe_error<double> logP(double logx) const {
    auto out = cte_int_ - alpha() * logx - beta() * std::exp(-logx);
    if (std::isfinite(out))
      return out;
    else
      return error_message("probability not finite:" + std::to_string(out));
  }

  double expected_variance() const { return beta() / alpha(); }
};

template <typename T, class Cova>
  requires Covariance<T, Cova>
class multivariate_gamma_normal_distribution
    : private log_inverse_gamma_distribution,
      multivariate_normal_distribution<T, Cova> {

public:
  using m_normal = multivariate_normal_distribution<T, Cova>;
  using log_inverse_gamma_distribution::alpha;
  using log_inverse_gamma_distribution::beta;
  using m_normal::mean;

  using m_normal::chi2;

  std::size_t size() const { return 1 + m_normal::size(); }

  multivariate_gamma_normal_distribution(log_inverse_gamma_distribution &&g,
                                         m_normal &&n)
      : log_inverse_gamma_distribution{std::move(g)}, m_normal{std::move(n)} {}
  Matrix<double> operator()(mt_64i &mt) {
    auto k = m_normal::size();
    auto sample_logVar = log_inverse_gamma_distribution::operator()(mt);
    auto sample_b = m_normal::operator()(mt);
    auto out = Matrix<double>(1, k + 1, false);
    out[0] = sample_logVar;

    for (std::size_t i = 0; i < k; ++i)
      out[i + 1] = sample_b[i];
    return out;
  }

  Maybe_error<double> logP(const Matrix<T> &x) const {
    assert(x.size() == mean().size() + 1);
    double logvar = x[0];
    auto logPvar = log_inverse_gamma_distribution::logP(logvar);
    if (!logPvar)
      return error_message("likelihood of variance wrong " + logPvar.error()());
    else {
      double var = std::exp(logvar);
      auto k = mean().size();
      auto b = Matrix<double>(1, k, false);
      for (std::size_t i = 0; i < k; ++i)
        b[i] = x[i + 1];
      auto chi_2 = m_normal::chi2(b) / var;

      double out =
          logPvar.value() - 0.5 * (k * log(2 * std::numbers::pi) +
                                   m_normal::logDetCov() + k * logvar + chi_2);
      if (std::isfinite(out))
        return out;
      else
        return error_message("likelihood not finite:" + std::to_string(out));
    }
  }

  auto &Gamma() const { return m_normal::cov_inv(); }
};

class dirchlet_distribution {
private:
  Matrix<std::gamma_distribution<>> m_gammas;
  Matrix<double> m_alpha;
  double m_alpha0;
  double m_Ba;

public:
  dirchlet_distribution() = default;
  dirchlet_distribution(Matrix<double> alpha)
      : m_gammas{applyMap(
            [](double a) {
              return std::gamma_distribution<>{a, 1.0};
            },
            alpha)},
        m_alpha{std::move(alpha)}, m_alpha0{var::sum(m_alpha)} {}

  Matrix<double> operator()(mt_64i &mt) {
    Matrix<double> s = applyMap(
        [&mt](std::gamma_distribution<> &g) { return g(mt); }, m_gammas);
    return s / var::sum(s);
  }
};
inline std::size_t sample_N(mt_64i &mt, double N) {
  auto r = std::uniform_real_distribution{}(mt);
  std::size_t out = std::floor(N);
  auto p = N - out;
  if (r < p)
    return out + 1;
  else
    return out;
}

class multinomial_distribution {
private:
  Matrix<double> m_P;

public:
  auto size() const { return m_P.size(); }
  multinomial_distribution() = default;
  multinomial_distribution(Matrix<double> P) : m_P{P} {
    auto sumP = var::sum(P);
    assert(std::abs(var::sum(P) - 1.0) <
           std::sqrt(std::numeric_limits<double>::epsilon()));
  }

  Matrix<std::size_t> operator()(mt_64i &mt, double N) {

    Matrix<std::size_t> out(m_P.nrows(), m_P.ncols());
    std::size_t N_remaining = sample_N(mt, N);
    double p_remaining = 1;
    auto k = out.size();
    for (std::size_t i = 0; i + 1 < out.size(); ++i) {
      auto n = std::binomial_distribution<std::size_t>(
          N_remaining, m_P[i] / p_remaining)(mt);
      N_remaining -= n;
      p_remaining -= m_P[i];
      out[i] = n;
    }
    out[k - 1] = N_remaining;
    return out;
  }

  double logP(const Matrix<std::size_t> &Nij, std::size_t N) const {
    double out = var::lgamma(N + 1.0);
    for (std::size_t i = 0; i < Nij.size(); ++i) {
      if (Nij[i] > 0)
        out = out + Nij[i] * std::log(m_P[i]) - var::lgamma(Nij[i] + 1.0);
    }
    return out;
  }

  Maybe_error<double> logP(const Matrix<std::size_t> &Nij, double t_N) const {
    auto N = var::sum(Nij);
    auto r = N - t_N;

    if (std::abs(r) > 1)
      return error_message("count mismatch");
    else {
      auto p = 1 - std::abs(r);

      auto out = logP(Nij, N) + log(p);
      if (std::isfinite(out))
        return out;
      else
        return error_message("not finite");
    }
  }
};

class multinomial_transition_distribution {
private:
  Matrix<double> m_P;

public:
  multinomial_transition_distribution() = default;
  multinomial_transition_distribution(Matrix<double> t_P) : m_P{t_P} {}

  Matrix<std::size_t> operator()(mt_64i &mt, Matrix<std::size_t> Ni) {
    assert(Ni.ncols() == m_P.nrows());
    Matrix<std::size_t> out(m_P.nrows(), m_P.ncols(), 0ul);
    for (std::size_t i = 0; i < m_P.nrows(); ++i) {
      double p_remaining = 1;
      std::size_t N_remaining = Ni[i];
      auto k = m_P.ncols();
      for (std::size_t j = 0; j + 1 < k; ++j) {
        if (N_remaining > 0) {
          auto n = std::binomial_distribution<std::size_t>(
              N_remaining, m_P(i, j) / p_remaining)(mt);
          N_remaining -= n;
          p_remaining -= m_P(i, j);
          out(i, j) = n;
        }
      }
      out(i, k - 1) = N_remaining;
    }
    return out;
  }
  double logP(const Matrix<std::size_t> &Nij, std::size_t N,
              std::size_t i) const {
    double out = var::lgamma(N + 1.0);
    for (std::size_t j = 0; j < Nij.ncols(); ++j) {
      if (Nij(i, j) > 0)
        out = out + Nij(i, j) * std::log(m_P(i, j)) -
              var::lgamma(Nij(i, j) + 1.0);
    }
    return out;
  }

  Maybe_error<double> logP(const Matrix<std::size_t> &Nij,
                           Matrix<std::size_t> Ni) const {

    if (var::sum(Nij) != var::sum(Ni))
      return error_message("count mismatch");
    else {
      double out;
      for (std::size_t i = 0; i < Ni.size(); ++i)
        out = out + logP(Nij, Ni[i], i);
      if (std::isfinite(out))
        return out;
      else
        return error_message("not finite");
    }
  }
};

class dirchlet_transition_distribution {
private:
  Matrix<double> m_P;

public:
  dirchlet_transition_distribution() = default;
  dirchlet_transition_distribution(Matrix<double> t_P) : m_P{t_P} {}

  Matrix<std::size_t> operator()(mt_64i &mt, Matrix<double> Ni) {
    assert(Ni.nrows() == m_P.nrows());
    Matrix<std::size_t> out(m_P.nrows(), m_P.ncols(), 0ul);
    for (std::size_t i = 0; i < m_P.nrows(); ++i) {
      double p_remaining = 1;
      std::size_t N_remaining = Ni[i];
      auto k = m_P.ncols();
      for (std::size_t j = 0; j + 1 < k; ++j) {
        if (N_remaining > 0) {
          auto n = std::binomial_distribution<std::size_t>(
              N_remaining, m_P(i, j) / p_remaining)(mt);
          N_remaining -= n;
          p_remaining -= m_P[i];
          out(i, j) = n;
        }
      }
      out(i, k - 1) = N_remaining;
      return out;
    }
    return out;
  }
  double logP(const Matrix<std::size_t> &Nij, std::size_t N,
              std::size_t i) const {
    double out = var::lgamma(N + 1.0);
    for (std::size_t j = 0; j < Nij.ncols(); ++i) {
      if (Nij(i, j) > 0)
        out = out + Nij(i, j) * std::log(m_P(i, j)) -
              var::lgamma(Nij(i, j) + 1.0);
    }
    return out;
  }

  Maybe_error<double> logP(const Matrix<std::size_t> &Nij,
                           Matrix<std::size_t> Ni) const {

    if (var::sum(Nij) != var::sum(Ni))
      return error_message("count mismatch");
    else {
      double out;
      for (std::size_t i = 0; i < Ni.size(); ++i)
        out = out + logP(Nij, Ni[i], i);
      if (std::isfinite(out))
        return out;
      else
        return error_message("not finite");
    }
  }
};

#endif // MULTIVARIATE_NORMAL_DISTRIBUTION_H
