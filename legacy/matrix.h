#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "maybe_error.h"
template <class>
class Matrix;
template <class>
class SubMatrix;

template <class>
class SymmetricMatrix;
template <class>
class SymPosDefMatrix;
template <class>
class UpTrianMatrix;
template <class>
class DownTrianMatrix;
template <class>
class DiagonalMatrix;

template <class>
struct element_of_Matrix;

template <template <class> class Matrix, class T>
    requires Matrix<T>::is_Matrix
struct element_of_Matrix<Matrix<T>> {
    using type = T;
};

template <class T>
using element_of_Matrix_t = typename element_of_Matrix<T>::type;

template <class>
struct is_Matrix : public std::false_type {};

template <class Matrix>
    requires Matrix::is_Matrix
struct is_Matrix<Matrix> : public std::true_type {};

template <class Matrix>
constexpr bool is_Matrix_v = is_Matrix<Matrix>::value;

constexpr bool Matrix_uses_vector = true;

template <class T, bool>
constexpr auto empty_matrix_container() {
    if constexpr (Matrix_uses_vector)
        return std::vector<T>{};
    else
        return nullptr;
}

template <class T>
T sqr(T x) {
    return x * x;
}

template <class T, bool>
auto initialize_matrix_container(std::size_t t_size, bool initialize = false) {
    if constexpr (Matrix_uses_vector) {
        if (initialize)
            return std::vector<T>(t_size, T{});
        else
            return std::vector<T>(t_size);
    } else {
        if (initialize)
            return new T[t_size]();
        else
            return new T[t_size];
    }
}

template <std::size_t I, typename... Args>
void write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i) {
    if (i == I)
        os << std::get<I>(tu);
}

template <std::size_t... Is, typename... Args>
std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i,
                            std::index_sequence<Is...>) {
    (write_tuple_i<Is>(os, tu, i), ...);
    return os;
}
template <typename... Args>
std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i) {
    return write_tuple_i(os, tu, i, std::index_sequence_for<Args...>());
}

inline double elemDivSafe(double x, double y, double eps = std::numeric_limits<double>::epsilon()) {
    if (std::abs(y) > eps)
        return x / y;
    else
        return 0;
}

inline double elemDivSoftAbs(double x, double y,
                             double eps = std::numeric_limits<double>::epsilon()) {
    double denom = std::sqrt(y * y + eps * eps);
    if (denom == 0.0)
        return 0.0;
    return x / denom;
}

namespace lapack {
Matrix<double> Lapack_Full_Product(const Matrix<double>& x, const Matrix<double>& y,
                                   bool transpose_x, bool transpose_y, double alpha = 1.0,
                                   double beta = 0.0);
Matrix<double> Lapack_Sym_Product(const SymmetricMatrix<double>& x, const Matrix<double>& y,
                                  bool first_symmetric_c, double alpha = 1.0, double beta = 0.0);

Matrix<double> Lapack_Triang_Product(const Matrix<double>& a, const Matrix<double>& b,
                                     bool up_triangular_in_c, bool triangular_first_in_c,
                                     bool transpose_A_in_c, bool ones_in_diag, double alpha = 1);

SymPosDefMatrix<double> Lapack_Product_Self_Transpose(const Matrix<double>& a,
                                                      bool first_transposed_in_c,
                                                      char UPLO_in_c = 'U', double alpha = 1,
                                                      double beta = 0);

std::pair<Matrix<double>, Matrix<double>> Lapack_QR(const Matrix<double>& x);
Maybe_error<Matrix<double>> Lapack_Full_inv(const Matrix<double>& a);

Maybe_error<SymmetricMatrix<double>> Lapack_Symm_inv(const SymmetricMatrix<double>& a);

Maybe_error<SymPosDefMatrix<double>> Lapack_SymmPosDef_inv(const SymPosDefMatrix<double>& a);
Maybe_error<DownTrianMatrix<double>> Lapack_chol(const SymPosDefMatrix<double>& x);

template <class T>
Maybe_error<SymPosDefMatrix<T>> Lapack_UT_Cholesky_inv(const UpTrianMatrix<T>& x);
template <class T>
Maybe_error<SymPosDefMatrix<T>> Lapack_LT_Cholesky_inv(const DownTrianMatrix<T>& x);

inline Maybe_error<DownTrianMatrix<double>> Lapack_LT_inv(const DownTrianMatrix<double>& x,
                                                          bool ones_in_diag);

Maybe_error<UpTrianMatrix<double>> Lapack_UT_inv(const UpTrianMatrix<double>& x, bool ones_in_diag);

Maybe_error<std::tuple<DiagonalMatrix<double>, Matrix<double>, Matrix<double>>> Lapack_EigenSystem(
    const Matrix<double>& x, bool does_permutations = true, bool does_diagonal_scaling = true,
    bool computes_eigenvalues_condition_numbers = false,
    bool computes_eigenvectors_condition_numbers = false);

Maybe_error<std::tuple<DiagonalMatrix<double>, Matrix<double>, Matrix<double>>>
    Lapack_Symm_EigenSystem(const SymmetricMatrix<double>& x, std::string kind = "lower");

}  // namespace lapack

template <class T>
class DiagPosDetMatrix;

// template <class Matrix>
//   requires(Matrix::is_Matrix)
// auto operator*(const Matrix &a, double b) {
//   return apply([&b](auto x) { return x * b; }, a);
// }

template <template <class> class Matrix, typename S>
    requires(Matrix<S>::is_Matrix)
auto operator*(const Matrix<S>& a, double b) {
    return applyMap([&b](auto x) { return x * b; }, a);
}

template <class Matrix>
    requires(Matrix::is_Matrix)
auto operator-(const Matrix& a) {
    return apply([](auto x) { return -x; }, a);
}

template <class Matrix>
    requires(Matrix::is_Matrix)
auto operator-(Matrix&& a) {
    for (std::size_t i = 0; i < a.size; ++i) a[i] = -a[i];
    return a;
}

inline Maybe_error<double> inv(double a) {
    if (a == 0)
        return error_message("Division by zero");
    else
        return 1.0 / a;
}

template <class Matrix>
    requires Matrix::is_Matrix
auto operator/(const Matrix& a, double b) {
    return apply([&b](auto x) { return x / b; }, a);
}

template <class Matrix>
    requires Matrix::is_Matrix
auto operator*(double b, const Matrix& a) {
    return apply([&b](auto x) { return x * b; }, a);
}

template <class F, template <class> class T_Matrix>
    requires(!is_Maybe_error<std::invoke_result_t<F, T_Matrix<double>>>)
auto apply(F&& f, Matrix<T_Matrix<double>> const& a) {
    using S = std::decay_t<std::invoke_result_t<F, T_Matrix<double>>>;
    Matrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a[i]);
    return x;
}


template <class F, template <class> class T_Matrix>
    requires(is_Maybe_error<std::invoke_result_t<F, T_Matrix<double>>>)
auto apply_maybe(F&& f, Matrix<T_Matrix<double>> const& a) -> Maybe_error<
    Matrix<typename std::decay_t<std::invoke_result_t<F, T_Matrix<double>>>::value_type>> {
    using S = typename std::invoke_result_t<F, T_Matrix<double>>::value_type;
    // using test=typename S::egber;
    Matrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto Maybe_x = std::invoke(std::forward<F>(f), a[i]);
        if (!Maybe_x)
            return Maybe_x.error();
        x[i] = std::move(Maybe_x.value());
    }
    return x;
}

template <class F, template <class> class T_Matrix>
    requires(!is_Maybe_error<std::invoke_result_t<F, double>> && T_Matrix<double>::is_Matrix &&
             (T_Matrix<double>::is_Diagonal || !T_Matrix<double>::is_Symmetric))
auto apply(F&& f, T_Matrix<double> const& a) {
    using S = std::decay_t<std::invoke_result_t<F, double>>;
    T_Matrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a[i]);
    return x;
}

template <template <class> class Matrix, class T, class F>
auto applyMap_i(F&& f, Matrix<T> const& a) {
    using S = std::decay_t<std::invoke_result_t<F, T, std::size_t>>;
    Matrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::invoke(std::forward<F>(f), a[i], i);
    return x;
}

template <class T>
class Matrix {
   private:
    std::size_t size_ = 0;
    std::size_t nrows_ = 0;
    std::size_t ncols_ = 0;
    std::conditional_t<Matrix_uses_vector, std::vector<T>, T*> x_ =
        empty_matrix_container<T, Matrix_uses_vector>();

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = false;
    static constexpr bool is_Diagonal = false;

    auto ncols() const { return ncols_; }
    auto nrows() const { return nrows_; }
    auto size() const { return size_; }
    ~Matrix() {
        if constexpr (!Matrix_uses_vector)  //  std::cerr<<"release "<<size()<<"\n";
            delete[] x_;
    }
    static int cell_width() { return 12; }
    explicit Matrix() {}
    explicit Matrix(std::size_t _nrows, std::size_t _ncols, bool initialize = true)
        : size_{_nrows * _ncols},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_, initialize)} {}
    explicit Matrix(std::size_t _nrows, std::size_t _ncols, T value)
        : size_{_nrows * _ncols},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_)} {
        for (std::size_t i = 0; i < size_; ++i) x_[i] = value;
    }

    explicit Matrix(std::size_t _nrows, std::size_t _ncols, const std::vector<T>& value)
        : size_{_nrows * _ncols},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_)} {
        for (std::size_t i = 0; i < size_; ++i) x_[i] = value[i];
    }

    explicit Matrix(std::size_t _nrows, std::size_t _ncols,
                    std::vector<std::tuple<std::size_t, std::size_t, T>>&& value)
        : size_{_nrows * _ncols},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_)} {
        for (std::size_t k = 0; k < value.size(); ++k) {
            auto [i, j, x] = value[k];
            (*this)(i, j) = x;
        }
    }

    Matrix(const Matrix& x)
        : size_{x.size()},
          nrows_{x.nrows()},
          ncols_{x.ncols()},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(x.size())} {
        for (std::size_t i = 0; i < size_; ++i) x_[i] = x[i];
    }
    Matrix(Matrix&& x)
        : size_{x.size()},
          nrows_{x.nrows()},
          ncols_{x.ncols()},
          x_{std::exchange(x.x_, empty_matrix_container<T, Matrix_uses_vector>())} {}

    friend class DiagonalMatrix<T>;

    Matrix& operator=(const Matrix& x) {
        if (&x != this) {
            size_ = x.size();
            nrows_ = x.nrows();
            ncols_ = x.ncols();
            if constexpr (!Matrix_uses_vector) {
                if (size() != x.size()) {
                    delete[] x_;
                    x_ = new T[x.size()];
                }
                for (std::size_t i = 0; i < size_; ++i) x_[i] = x[i];
            } else {
                x_ = x.x_;
            }
        }
        return *this;
    }

    Matrix& operator=(Matrix&& x) {
        if (&x != this) {
            size_ = x.size();
            nrows_ = x.nrows();
            ncols_ = x.ncols();
            if constexpr (!Matrix_uses_vector) {
                delete[] x_;
                x_ = std::exchange(x.x_, nullptr);
            } else {
                x_ = std::move(x.x_);
            }
        }
        return *this;
    }

    auto& operator[](std::size_t i) {
        assert(i < size_);
        return x_[i];
    }
    auto& operator[](std::size_t i) const {
        assert(i < size_);
        return x_[i];
    }

    Matrix<double> operator[](std::pair<std::size_t, std::size_t> ij) const {
        if (size() == ncols()) {
            auto out = Matrix<double>(1ul, ij.second - ij.first + 1);
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = (*this)[ij.first + i];
            return out;
        } else {
            auto out = Matrix<double>(ij.second - ij.first + 1, 1ul);
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = (*this)[ij.first + i];
            return out;
        }
    }

    void set(std::pair<std::size_t, std::size_t> ij, Matrix<double> x) {
        for (std::size_t i = ij.first; i <= ij.second; ++i) (*this)[i] = x[i - ij.first];
    }

    auto& operator()(std::size_t i) {
        assert(i < size_);
        return x_[i];
    }
    auto& operator()(std::size_t i) const {
        assert(i < size_);
        return x_[i];
    }
    auto& operator()(std::size_t i, std::size_t j) {
        assert(i < nrows() && j < ncols());
        return x_[i * ncols_ + j];
    }

    void set(std::size_t i, std::size_t j, const T& x) { operator()(i, j) = x; }
    void set(std::size_t i, std::size_t j, T&& x) { operator()(i, j) = std::move(x); }
    void set(std::size_t i, T x) { operator[](i) = x; }

    auto operator()(std::size_t i, const char* ch) const {
        assert(i < nrows() && std::strcmp(ch, ":") == 0);
        Matrix out(1, ncols());
        for (std::size_t j = 0; j < ncols(); ++j) out(0, j) = (*this)(i, j);
        return out;
    }

    auto operator()(const char* ch, std::size_t j) const {
        assert(j < ncols() && std::strcmp(ch, ":") == 0);
        Matrix out(nrows(), 1);
        for (std::size_t i = 0; i < nrows(); ++i) out(i, 0) = (*this)(i, j);
        return out;
    }

    auto& operator()(std::size_t i, std::size_t j) const {
        assert((i < nrows()) && (j < ncols()));
        return x_[i * ncols_ + j];
    }

    template <class S>
        requires(std::is_same_v<S, T> && std::is_same_v<S, double>)
    friend auto operator*(const Matrix& a, const Matrix<S>& b) {
        if (a.size() == 0)
            return a;
        else if (b.size() == 0)
            return b;
        else
            return lapack::Lapack_Full_Product(a, b, false, false);
    }

    template <class S>
        requires(std::is_same_v<S, T> && !std::is_same_v<S, double>)
    friend auto operator*(const Matrix& a, const Matrix<S>& b) {
        if (a.size() == 0)
            return a;
        else if (b.size() == 0)
            return b;
        else {
            Matrix<T> out(a.nrows(), b.ncols(), T{});
            for (std::size_t i = 0; i < a.nrows(); ++i)
                for (std::size_t j = 0; j < b.ncols(); ++j)
                    for (std::size_t k = 0; k < a.ncols(); ++k) out(i, j) += a(i, k) * b(k, j);
            return out;
        }
    }

    friend auto operator+(const Matrix& a, const Matrix& b) {
        if (a.size() == 0)
            return b;
        else if (b.size() == 0)
            return a;
        return zip([](auto x, auto y) { return x + y; }, a, b);
    }

    template <class S>
        requires(S::is_Matrix)
    friend auto operator+(const Matrix& a, const Matrix<S>& b) {
        using R = std::decay_t<decltype(a[0] + b[0])>;
        if (a.size() == 0)
            return Matrix<R>(b);
        else if (b.size() == 0)
            return Matrix<R>(a);
        else
            return zip([](auto x, auto y) { return x + y; }, a, b);
    }

    template <class S>
        requires(std::is_arithmetic_v<S>)
    friend auto operator-(const Matrix& a, const Matrix<S>& b) {
        using R = std::decay_t<decltype(a[0] - b[0])>;
        if (a.size() == 0)
            return applyMap([](auto e) { return R(e); }, b);
        else if (b.size() == 0)
            return Matrix<R>(a);
        else
            return zip([](auto x, auto y) { return x - y; }, a, b);
    }

    friend auto operator-(const Matrix& a, const Matrix& b) {
        if (a.size() == 0)
            return a;  // avoid promotion issues for non-floating T
        else if (b.size() == 0)
            return a;
        else
            return zip([](auto x, auto y) { return x - y; }, a, b);
    }

    template <class F>
    friend auto zip(F&& f, const Matrix& x, const Matrix& y) {
        assert(same_dimensions(x, y) && "same size");
        Matrix out(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = f(x[i], y[i]);
        return out;
    }


    

    template <class F, class S>
        requires S::is_Matrix
    friend auto zip(F&& f, const Matrix& x, const Matrix<S>& y) {
        // assert(same_dimensions(x, y) && "same size");
        using R = std::invoke_result_t<F, T, S>;

        Matrix<R> out(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = f(x[i], y[i]);
        return out;
    }

    template <class F, class S>
        requires(std::is_arithmetic_v<S>)
    friend auto zip(F&& f, const Matrix& x, const Matrix<S>& y) {
        // assert(same_dimensions(x, y) && "same size");
        using R = std::invoke_result_t<F, T, S>;

        Matrix<R> out(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = f(x[i], y[i]);
        return out;
    }

    friend auto TranspMult(const Matrix& a, const Matrix& b) {
        if (a.size() == 0)
            return a;
        else if (b.size() == 0)
            return b;
        else
            return lapack::Lapack_Full_Product(a, b, true, false);
    }

    friend auto MultTransp(const Matrix& a, const Matrix& b) {
        return lapack::Lapack_Full_Product(a, b, true, false);
    }

    friend auto inv(const Matrix& a) { return lapack::Lapack_Full_inv(a); }

    friend auto tr(const Matrix& a) {
        Matrix out(a.ncols(), a.nrows(), false);
        for (std::size_t i = 0; i < out.nrows(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j) out(i, j) = a(j, i);
        return out;
    }

    template <class F>
        requires std::constructible_from<T, std::invoke_result_t<F, T>>
    friend auto apply(F&& f, Matrix&& x) {
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(x[i]);
        return x;
    }

    template <class F>
    friend auto applyMap(F&& f, Matrix const& a) {
        using S = std::decay_t<std::invoke_result_t<F, T>>;
        Matrix<S> x(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a[i]);
        return x;
    }

    template <class F>
    friend auto applyMap(F&& f, Matrix& a) {
        using S = std::decay_t<std::invoke_result_t<F, T&>>;
        Matrix<S> x(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a[i]);
        return x;
    }

    template <class F>
        requires(is_Maybe_error<std::invoke_result_t<F, T>>)
    friend Maybe_error<Matrix<T>> apply(F&& f, Matrix const& a) {
        Matrix x(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = 0; j < x.ncols(); ++j) {
                auto v = f(a(i, j));
                if (v)
                    x(i, j) = std::move(v.value());
                else {
                    return v.error()() + " at cell (" + std::to_string(i) + "," +
                           std::to_string(j) + ")";
                }
            }
        return x;
    }

    friend bool same_dimensions(const Matrix& x, const Matrix& y) {
        return x.size() == y.size() && x.nrows() == y.nrows() && x.ncols() == y.ncols();
    }
    template <class S>
        requires S::is_Matrix
    friend bool same_dimensions(const Matrix& x, const Matrix<S>& y) {
        return x.size() == y.size() && x.nrows() == y.nrows() && x.ncols() == y.ncols();
    }

    friend T xtx(const Matrix& x) {
        assert(x.ncols() == 1);
        auto out = T{};
        for (std::size_t i = 0; i < x.nrows(); ++i) out += x[i] * x[i];
        return out;
    }

    template <class F>
    friend auto reduce(F&& f, const Matrix& x) {
        auto cum = x[0];
        for (std::size_t i = 1; i < x.size(); ++i) cum = f(cum, x[i]);
        return cum;
    }

    template <class F, class X>
    friend auto reduce(F&& f, X init, const Matrix& x) {
        for (std::size_t i = 0; i < x.size(); ++i) init = f(init, x[i]);
        return init;
    }

    template <class F>
    friend auto reduce_ij(F&& f, const Matrix& x, T init) {
        auto cum = init;
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = 0; j < x.ncols(); ++j) cum = f(cum, i, j, x(i, j));
        return cum;
    }

    friend auto& operator<<(std::ostream& os, const Matrix& x) {
        os << "\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os << std::setw(cell_width()) << x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};

template <template <class> class T_Matrix, template <class> class S_Matrix>
    requires(T_Matrix<double>::is_Matrix && S_Matrix<double>::is_Matrix)
auto operator-(const Matrix<T_Matrix<double>>& a, const Matrix<S_Matrix<double>>& b) {
    if (a.size() == 0)
        return apply(
            [](S_Matrix<double> const& e) { return T_Matrix<double>(e.nrows(), e.ncols()) - e; },
            b);
    else if (b.size() == 0)
        return apply(
            [](T_Matrix<double> const& e) { return e - S_Matrix<double>(e.nrows(), e.ncols()); },
            a);
    else
        return zip([](T_Matrix<double> const& x, S_Matrix<double> const& y) { return x - y; }, a,
                   b);
}

template <class T>
class SymmetricMatrix : public Matrix<T> {
    using base_type = Matrix<T>;

   protected:
    explicit SymmetricMatrix(base_type&& x) : base_type{std::move(x)} {}
    explicit SymmetricMatrix(base_type const& x) : base_type{x} {}

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = true;
    static constexpr bool is_Diagonal = false;

    operator Matrix<T> const&() { return static_cast<Matrix<T> const&>(*this); }

    explicit SymmetricMatrix(std::size_t _nrows) : base_type(_nrows, _nrows, false) {}
    explicit SymmetricMatrix(std::size_t _nrows, std::size_t, T value)
        : base_type(_nrows, _nrows, value) {}
    explicit SymmetricMatrix(std::size_t _nrows, std::size_t) : base_type(_nrows, _nrows) {}

    SymmetricMatrix(const SymmetricMatrix& x) : base_type(x) {}
    SymmetricMatrix(SymmetricMatrix&& x) : base_type(std::move(x)) {}

    SymmetricMatrix() {}

    SymmetricMatrix& operator=(const SymmetricMatrix& x) {
        static_cast<base_type&>(*this) = static_cast<base_type const&>(x);
        return *this;
    }
    SymmetricMatrix& operator=(SymmetricMatrix&& x) {
        static_cast<base_type&>(*this) = std::move(static_cast<base_type&>(x));
        return *this;
    }

    void set(std::size_t i, std::size_t j, const T& x) {
        if (size() > 0) {
            base_type::operator()(i, j) = x;
            if (i != j)
                base_type::operator()(j, i) = x;
        }
    }
    auto operator()(std::size_t i, std::size_t j) const { return base_type::operator()(i, j); }
    auto& operator[](std::size_t i) const { return base_type::operator[](i); }
    auto& operator[](std::size_t i) { return base_type::operator[](i); }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~SymmetricMatrix() {}

    friend auto operator*(const SymmetricMatrix& a, const SymmetricMatrix& b) {
        return lapack::Lapack_Sym_Product(a, b, true);
    }
    friend auto operator*(const SymmetricMatrix& a, const Matrix<double>& b) {
        assert(a.ncols() == b.nrows());
        return lapack::Lapack_Sym_Product(a, b, true);
    }
    friend auto operator*(const Matrix<double>& b, const SymmetricMatrix& a) {
        assert(b.ncols() == a.nrows());

        return lapack::Lapack_Sym_Product(a, b, false);
    }

    friend SymmetricMatrix operator*(const SymmetricMatrix& a, T b) {
        return apply([&b](auto x) { return x * b; }, a);
    }

    friend SymmetricMatrix operator*(T b, const SymmetricMatrix& a) {
        return apply([&b](auto x) { return x * b; }, a);
    }

    template <class S>
        requires(std::is_same_v<decltype(T{} + S{}), T>)
    friend auto operator+(const SymmetricMatrix& a, const SymmetricMatrix<S>& b) {
        if (a.size() == 0)
            return b;
        else if (b.size() == 0)
            return a;
        else
            return zip([](auto x, auto y) { return x + y; }, a, b);
    }

    template <class S>
        requires(std::is_same_v<decltype(T{} - S{}), T>)
    friend auto operator-(const SymmetricMatrix& a, const SymmetricMatrix<S>& b) {
        if (a.size() == 0)
            return b * (-1.0);
        else if (b.size() == 0)
            return a;
        else
            return zip([](auto x, auto y) { return x - y; }, a, b);
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator+(const SymmetricMatrix& a, const DiagonalMatrix<S>& b) {
        auto out = a;
        if (out.size() == 0)
            out = SymmetricMatrix(b.nrows(), b.ncols());
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) + b(i, i));
        return out;
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator+(const DiagonalMatrix<S>& b, const SymmetricMatrix& a) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) + b(i, i));
        return out;
    }

    friend auto inv(const SymmetricMatrix& a) { return lapack::Lapack_Symm_inv(a); }

    friend auto& tr(const SymmetricMatrix& a) { return a; }

    template <class F>
    friend auto apply(F&& f, SymmetricMatrix&& x) {
        for (std::size_t i = 0; i < x.size(); ++i)
            for (std::size_t j = i; j < x.size(); ++j) x.set(i, j, f(x(i, j)));
        return x;
    }

    template <class F>
    friend auto zip(F&& f, const SymmetricMatrix& x, const SymmetricMatrix& y) {
        assert(x.nrows() == y.nrows() && "same size");

        SymmetricMatrix out(x.nrows(), false);
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = i; j < x.ncols(); ++j) out.set(i, j, f(x(i, j), y(i, j)));
        return out;
    }

    friend auto operator-(const SymmetricMatrix& x, const DiagonalMatrix<T>& y) {
        SymmetricMatrix out(x);
        for (std::size_t i = 0; i < x.nrows(); ++i) out.set(i, i, x(i, i) - y(i, i));
        return out;
    }
    friend auto operator-(const DiagonalMatrix<T>& x, const SymmetricMatrix& y) {
        SymmetricMatrix out(-y);
        for (std::size_t i = 0; i < x.nrows(); ++i) out.set(i, i, x(i, i) - y(i, i));
        return out;
    }
};

template <class T>
auto xtAx(const Matrix<T>& x, const SymmetricMatrix<T>& A) {
    assert(x.nrows() == A.nrows());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * x[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j) out += 2 * x[i] * A(i, j) * x[j];
    }
    return out;
}

template <class T>
auto xAxt(const Matrix<T>& x, const SymmetricMatrix<T>& A) {
    assert(x.ncols() == A.nrows());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * x[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j) out += 2 * x[i] * A(i, j) * x[j];
    }
    return out;
}

template <class T>
auto xAyt(const Matrix<T>& x, const SymmetricMatrix<T>& A, const Matrix<T>& y) {
    assert(x.ncols() == A.nrows());
    assert(y.ncols() == A.ncols());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * y[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j) out += 2 * x[i] * A(i, j) * y[j];
    }
    return out;
}

template <class T>
auto xtAy(const Matrix<T>& x, const SymmetricMatrix<T>& A, const Matrix<T>& y) {
    assert(x.nrows() == A.nrows());
    assert(y.nrows() == A.ncols());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * y[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j) out += 2 * x[i] * A(i, j) * y[j];
    }
    return out;
}

template <class F, class T>
auto apply(F&& f, const SymmetricMatrix<T>& x) {
    SymmetricMatrix<T> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i; j < x.ncols(); ++j) out.set(i, j, f(x(i, j)));
    return out;
}

inline auto elemMult(const Matrix<double>& x, const Matrix<double>& y) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = x[i] * y[i];
    return out;
}

template <class T, class S>
auto elemMult(const Matrix<S>& x, const Matrix<T>& y) -> Matrix<std::decay_t<decltype(T{} * S{})>> {
    Matrix<std::decay_t<decltype(T{} * S{})>> out(x.nrows(), x.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = x[i] * y[i];
    return out;
}

inline auto elemDiv(const Matrix<double>& x, const Matrix<double>& y) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = x[i] / y[i];
    return out;
}

inline auto elemDivSafe(const Matrix<double>& x, const Matrix<double>& y, double eps) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = elemDivSafe(x[i], y[i], eps);
    return out;
}

inline auto elemDivSafe(const Matrix<double>& x, double y, double eps) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = elemDivSafe(x[i], y, eps);
    return out;
}

inline auto elemDivSoftAbs(const Matrix<double>& x, const Matrix<double>& y, double eps) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = elemDivSoftAbs(x[i], y[i], eps);
    return out;
}

inline auto elemDivSoftAbs(const Matrix<double>& x, double y, double eps) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = elemDivSoftAbs(x[i], y, eps);
    return out;
}

template <class T>
void copy_UT_to_LT(Matrix<T>& x) {
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = 0; j < i; ++j) x(i, j) = x(j, i);
}

template <class T>
void copy_LT_to_UT(Matrix<T>& x) {
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i + 1; j < x.ncols(); ++j) x(i, j) = x(j, i);
}

template <class Ts>
UpTrianMatrix<Ts> fill_LT_zeros(SymPosDefMatrix<Ts>&& x);

template <class Ts>
DownTrianMatrix<Ts> fill_UT_zeros(SymPosDefMatrix<Ts>&& x);

template <class T>
class UpTrianMatrix : public Matrix<T> {
    using base_type = Matrix<T>;

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = false;
    static constexpr bool is_Diagonal = false;

    void fill_the_zeros() {
        for (std::size_t i = 0; i < nrows(); ++i)
            for (std::size_t j = 0; j < i; ++j) base_type::operator()(i, j) = 0;
    }
    operator Matrix<T> const&() { return static_cast<Matrix<T> const&>(*this); }
    explicit UpTrianMatrix(base_type&& x) : base_type{std::move(x)} {}

    explicit UpTrianMatrix(std::size_t _nrows, bool initialize = true)
        : base_type(_nrows, _nrows, initialize) {}
    explicit UpTrianMatrix(std::size_t _nrows, T value) : base_type(_nrows, _nrows, value) {
        fill_the_zeros();
    }

    UpTrianMatrix(const UpTrianMatrix& x) : base_type(x) {}
    UpTrianMatrix(UpTrianMatrix&& x) : base_type(std::move(x)) {}

    UpTrianMatrix& operator=(const UpTrianMatrix& x) {
        static_cast<base_type&>(*this) = static_cast<base_type const&>(x);
        return *this;
    }

    template <class Ts>
    friend UpTrianMatrix fillzeros(SymPosDefMatrix<Ts>&& x);

    void set(std::size_t i, std::size_t j, double x) {
        assert(j >= i);
        base_type::operator()(i, j) = x;
    }
    auto& operator()(std::size_t i, std::size_t j) const { return base_type::operator()(i, j); }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~UpTrianMatrix() {}

    friend auto inv(const UpTrianMatrix& a) { return lapack::Lapack_UT_inv(a, false); }

    friend auto operator*(const UpTrianMatrix& a, const Matrix<T>& b) {
        return lapack::Lapack_Triang_Product(a, b, true, true, false, false);
    }
    friend auto operator*(const Matrix<T>& a, const UpTrianMatrix& b) {
        return lapack::Lapack_Triang_Product(b, a, true, false, false, false);
    }

    friend auto& operator<<(std::ostream& os, const UpTrianMatrix& x) {
        os << "Upper triang part\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os << std::setw(base_type::cell_width()) << x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};

template <class T>
UpTrianMatrix<T> fill_LT_zeros(SymPosDefMatrix<T>&& x) {
    UpTrianMatrix<T> out(std::move(x));
    out.fill_the_zeros();
    return out;
}

template <class T>
class DownTrianMatrix : public Matrix<T> {
    using base_type = Matrix<T>;

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = false;
    static constexpr bool is_Diagonal = false;

    void fill_the_zeros() {
        for (std::size_t i = 0; i < nrows(); ++i)
            for (std::size_t j = i + 1; j < ncols(); ++j) base_type::operator()(i, j) = 0;
    }
    operator Matrix<T> const&() { return static_cast<Matrix<T> const&>(*this); }
    explicit DownTrianMatrix(base_type&& x) : base_type{std::move(x)} {}
    explicit DownTrianMatrix(base_type const& x) : base_type{x} {}

    explicit DownTrianMatrix(std::size_t _nrows, bool initialize = true)
        : base_type(_nrows, _nrows, initialize) {}
    explicit DownTrianMatrix(std::size_t _nrows, T value) : base_type(_nrows, _nrows, value) {
        fill_the_zeros();
    }

    DownTrianMatrix(const DownTrianMatrix& x) : base_type(x) {}
    DownTrianMatrix(DownTrianMatrix&& x) : base_type(std::move(x)) {}

    DownTrianMatrix& operator=(const DownTrianMatrix& x) {
        static_cast<base_type&>(*this) = static_cast<base_type const&>(x);
        return *this;
    }

    template <class Ts>
    friend DownTrianMatrix fillzeros(SymPosDefMatrix<Ts>&& x);

    void set(std::size_t i, std::size_t j, double x) {
        assert(j <= i);
        base_type::operator()(i, j) = x;
    }
    auto& operator()(std::size_t i, std::size_t j) const { return base_type::operator()(i, j); }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~DownTrianMatrix() {}

    friend auto inv(const DownTrianMatrix& a) { return lapack::Lapack_LT_inv(a, false); }

    friend auto operator*(const DownTrianMatrix& a, const Matrix<T>& b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }
    friend auto operator*(const DownTrianMatrix& a, const DownTrianMatrix& b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }
    friend auto operator*(const DownTrianMatrix& a, const UpTrianMatrix<T>& b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }

    friend auto operator*(const Matrix<T>& b, const DownTrianMatrix& a) {
        return lapack::Lapack_Triang_Product(a, b, false, false, false, false);
    }
    friend auto operator*(const UpTrianMatrix<T>& b, const DownTrianMatrix& a) {
        return lapack::Lapack_Triang_Product(a, b, false, false, false, false);
    }

    friend auto& operator<<(std::ostream& os, const DownTrianMatrix& x) {
        os << "Down triang part\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os << std::setw(base_type::cell_width()) << x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};

template <class T>
DownTrianMatrix<T> fill_UT_zeros(SymPosDefMatrix<T>&& x) {
    DownTrianMatrix<T> out(std::move(x));
    out.fill_the_zeros();
    return out;
}

template <class T>
auto tr(const UpTrianMatrix<T>& x) {
    return DownTrianMatrix<T>(tr(static_cast<Matrix<T> const&>(x)));
}

template <class T>
auto tr(const DownTrianMatrix<T>& x) {
    return UpTrianMatrix<T>(tr(static_cast<Matrix<T> const&>(x)));
}

template <class T>
Maybe_error<double> logdet(const UpTrianMatrix<T>& x) {
    double out = 0;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (x(i, i) > 0)
            out += std::log(x(i, i));
        else
            return error_message("logdet of negative value at (" + std::to_string(i) + "," +
                                 std::to_string(i) + ")");
    }
    return out;
}

template <class T>
Maybe_error<double> logdet(const DownTrianMatrix<T>& x) {
    double out = 0;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (x(i, i) > 0)
            out += std::log(x(i, i));
        else
            return error_message("logdet of negative value at (" + std::to_string(i) + "," +
                                 std::to_string(i) + ")");
    }
    return out;
}

template <class T>
class SymPosDefMatrix : public SymmetricMatrix<T> {
    using base_type = SymmetricMatrix<T>;

    explicit SymPosDefMatrix(base_type&& x) : base_type{std::move(x)} {}
    explicit SymPosDefMatrix(Matrix<T>&& x) : base_type{std::move(x)} {}

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = true;
    static constexpr bool is_Diagonal = false;

    template <class K>
    friend Maybe_error<SymPosDefMatrix<K>> lapack::Lapack_UT_Cholesky_inv(
        const UpTrianMatrix<K>& x);

    template <class K>
    friend Maybe_error<SymPosDefMatrix<K>> lapack::Lapack_LT_Cholesky_inv(
        const DownTrianMatrix<K>& x);
    operator Matrix<T> const&() { return static_cast<Matrix<T> const&>(*this); }
    operator SymmetricMatrix<T> const&() { return static_cast<SymmetricMatrix<T> const&>(*this); }

    explicit SymPosDefMatrix(std::size_t _nrows, std::size_t _ncols, bool initialize = true)
        : base_type(_nrows, _nrows, initialize) {}
    explicit SymPosDefMatrix(std::size_t _nrows, std::size_t, T value) : base_type(_nrows, value) {}

    explicit SymPosDefMatrix(std::size_t _nrows, T value) : base_type(_nrows, value) {}

    static SymPosDefMatrix I_sware_it_is_possitive(base_type x) {
        return SymPosDefMatrix(std::move(x));
    }
    SymPosDefMatrix(const SymPosDefMatrix& x) : base_type(x) {}
    SymPosDefMatrix(SymPosDefMatrix&& x) : base_type(std::move(static_cast<base_type&>(x))) {}

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    SymPosDefMatrix(const DiagPosDetMatrix<S>& b) : base_type{b.nrows(), b.ncols(), S{}} {
        for (std::size_t i = 0; i < b.nrows(); ++i) set(i, i, b(i, i));
    }

    friend auto operator/(SymPosDefMatrix&& x, double y) {
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = i; j < x.ncols(); ++j) x.set(i, j, x(i, j) / y);
        return x;
    }
    friend auto operator/(SymPosDefMatrix const& x, double y) {
        auto z = x;
        return std::move(z) / y;
    }

    SymPosDefMatrix& operator=(const SymPosDefMatrix& x) {
        base_type::operator=(static_cast<base_type const&>(x));
        return *this;
    }

    SymPosDefMatrix() {}

    void set(std::size_t i, std::size_t j, double x) { base_type::set(i, j, x); }
    auto operator()(std::size_t i, std::size_t j) const { return base_type::operator()(i, j); }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~SymPosDefMatrix() {}

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator+(const SymPosDefMatrix& a, const SymPosDefMatrix<S>& b) {
        if (a.size() == 0)
            return b;
        else if (b.size() == 0)
            return a;
        else
            return SymPosDefMatrix(zip([](auto x, auto y) { return x + y; }, a, b));
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator+(const SymPosDefMatrix& a, const DiagPosDetMatrix<S>& b) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) + b(i, i));
        return out;
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator+(const DiagPosDetMatrix<S>& b, const SymPosDefMatrix& a) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) + b(i, i));
        return out;
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator-(const SymPosDefMatrix& a, const DiagPosDetMatrix<S>& b) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) - b(i, i));
        return out;
    }
    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator-(const SymPosDefMatrix& a, const SymPosDefMatrix<S>& b) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i)
            for (std::size_t j = 0; j < b.ncols(); ++j) out.set(i, j, out(i, j) - b(i, j));
        return out;
    }

    template <class S>
        requires(std::is_same_v<T, std::decay_t<decltype(T{}, S{})>>)
    friend auto operator-(const DiagPosDetMatrix<S>& b, const SymPosDefMatrix& a) {
        auto out = a * (-1.0);
        for (std::size_t i = 0; i < b.nrows(); ++i) out.set(i, i, out(i, i) + b(i, i));
        return out;
    }
    friend SymPosDefMatrix operator*(const SymPosDefMatrix& a, T b) {
        return SymPosDefMatrix(static_cast<SymmetricMatrix<T> const&>(a) * b);
    }

    friend SymPosDefMatrix operator*(T b, const SymPosDefMatrix& a) {
        return SymPosDefMatrix(b * static_cast<SymmetricMatrix<T> const&>(a));
    }

    friend auto inv(const SymPosDefMatrix& a) { return lapack::Lapack_SymmPosDef_inv(a); }

    friend auto tr(const SymPosDefMatrix& a) { return a; }

    friend Maybe_error<DownTrianMatrix<T>> cholesky(const SymPosDefMatrix& x) {
        return lapack::Lapack_chol(x);
    }

    friend Maybe_error<double> logdet(const SymPosDefMatrix& x) {
        auto chol = cholesky(x);
        if (chol)
            return 2.0 * logdet(chol.value());
        else
            return Maybe_error<double>(chol.error()() + " in calculation of logdet");
    }
};
template <class T>
Maybe_error<DownTrianMatrix<T>> cholesky(const SymPosDefMatrix<T>& x);

template <class T>
Maybe_error<DiagPosDetMatrix<T>> cholesky(const DiagPosDetMatrix<T>& x);

template <class T>
class DiagonalMatrix {
   private:
    std::size_t size_ = 0;
    std::size_t nrows_ = 0;
    std::size_t ncols_ = 0;
    // T *x_ = nullptr;
    std::conditional_t<Matrix_uses_vector, std::vector<T>, T*> x_ =
        empty_matrix_container<T, Matrix_uses_vector>();

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = true;
    static constexpr bool is_Diagonal = true;

    explicit DiagonalMatrix(std::size_t _nrows, std::size_t _ncols, bool initialize = true)
        : size_{std::min(_nrows, _ncols)},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_, initialize)} {}
    explicit DiagonalMatrix(std::size_t _nrows, std::size_t _ncols, T value)
        : size_{std::min(_nrows, _ncols)},
          nrows_{_nrows},
          ncols_{_ncols},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_)} {
        for (std::size_t i = 0; i < size_; ++i) x_[i] = value;
    }
    //  explicit DiagonalMatrix(std::initializer_list<T> &&a)
    //      : size_{a.size()}, nrows_{a.size()}, ncols_{a.size()}, x_{new
    //      T[size_]} {
    //    std::copy(a.begin(), a.end(), x_);
    //  }

    explicit DiagonalMatrix(const Matrix<T>& a)
        : size_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : std::min(a.nrows(), a.ncols())},
          nrows_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.nrows()},
          ncols_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.ncols()},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(size_)} {
        if (a.ncols() == 1 || a.nrows() == 1) {
            for (std::size_t i = 0; i < size(); ++i) (*this)[i] = a[i];
        } else {
            for (std::size_t i = 0; i < size(); ++i) (*this)[i] = a(i, i);
        }
    }

    explicit DiagonalMatrix(const std::vector<T>& a)
        : DiagonalMatrix(Matrix<T>(a.size(), 1ul, a)) {}

    explicit DiagonalMatrix(Matrix<T>&& a)
        : size_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : std::min(a.nrows(), a.ncols())},
          nrows_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.nrows()},
          ncols_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.ncols()},
          x_{} {
        if (a.ncols() == 1 || a.nrows() == 1) {
            if constexpr (!Matrix_uses_vector)
                x_ = std::exchange(a.x_, nullptr);
            else
                x_ = std::move(a.x_);
        } else {
            x_ = initialize_matrix_container<T, Matrix_uses_vector>(size_);
            for (std::size_t i = 0; i < size(); ++i) (*this)[i] = a(i, i);
        }
    }
    DiagonalMatrix() {}
    DiagonalMatrix(const DiagonalMatrix& x)
        : size_{x.size()},
          nrows_{x.nrows()},
          ncols_{x.ncols()},
          x_{initialize_matrix_container<T, Matrix_uses_vector>(x.size())} {
        for (std::size_t i = 0; i < size_; ++i) x_[i] = x[i];
    }
    DiagonalMatrix(DiagonalMatrix&& x)
        : size_{x.size()},
          nrows_{x.nrows()},
          ncols_{x.ncols()},
          x_{std::exchange(x.x_, empty_matrix_container<T, Matrix_uses_vector>())} {}

    DiagonalMatrix& operator=(const DiagonalMatrix& x) {
        if (&x != this) {
            if (size() != x.size()) {
                size_ = x.size();
                nrows_ = x.nrows();
                ncols_ = x.ncols();
                if constexpr (!Matrix_uses_vector) {
                    delete[] x_;
                    x_ = new T[x.size()];
                    for (std::size_t i = 0; i < size_; ++i) x_[i] = x[i];
                } else {
                    x_ = x.x_;
                }
            }
        }
        return *this;
    }

    DiagonalMatrix& operator=(DiagonalMatrix&& x) {
        if (&x != this) {
            size_ = x.size();
            nrows_ = x.nrows();
            ncols_ = x.ncols();
            x_ = std::exchange(x.x_, empty_matrix_container<T, Matrix_uses_vector>());
        }
        return *this;
    }

    auto& operator[](std::size_t i) { return x_[i]; }
    auto& operator[](std::size_t i) const { return x_[i]; }
    auto operator()(std::size_t i, std::size_t j) const {
        if (i == j)
            return x_[i];
        else
            return T{};
    }

    auto ncols() const { return ncols_; }
    auto nrows() const { return nrows_; }
    auto size() const { return size_; }
    ~DiagonalMatrix() {
        if constexpr (!Matrix_uses_vector)
            delete[] x_;
    }
    static int cell_width() { return 12; }

    friend auto operator*(const DiagonalMatrix& a, const Matrix<double>& b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        Matrix<double> out(a.nrows(), b.ncols(), false);
        for (std::size_t i = 0; i < a.size(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j) out(i, j) = a[i] * b(i, j);
        for (std::size_t i = a.size(); i < out.nrows(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j) out(i, j) = 0;
        return out;
    }

    friend auto operator*(const Matrix<double>& a, const DiagonalMatrix& b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        Matrix<double> out(a.nrows(), b.ncols(), false);
        for (std::size_t i = 0; i < out.nrows(); ++i) {
            for (std::size_t j = 0; j < b.size(); ++j) out(i, j) = a(i, j) * b[j];
            for (std::size_t j = b.size(); j < out.ncols(); ++j) out(i, j) = 0;
        }
        return out;
    }

    template <class F>
    friend auto zip(F&& f, const DiagonalMatrix& x, const DiagonalMatrix& y) {
        assert(x.size() == y.size() && "same size");
        DiagonalMatrix out(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = f(x[i], y[i]);
        return out;
    }

    friend auto operator+(const DiagonalMatrix& a, const DiagonalMatrix& b) {
        if (a.size() == 0)
            return b;
        else if (b.size() == 0)
            return a;
        return zip([](auto x, auto y) { return x + y; }, a, b);
    }

    friend Matrix<T> operator+(const DiagonalMatrix& a, const Matrix<T>& b) {
        auto out = b;
        for (std::size_t i = 0; i < a.nrows(); ++i) out(i, i) = a(i, i) + out(i, i);
        return out;
    }
    friend Matrix<T> operator+(const Matrix<T>& b, const DiagonalMatrix& a) {
        auto out = b;
        for (std::size_t i = 0; i < a.nrows(); ++i) out(i, i) = out(i, i) + a(i, i);
        return out;
    }

    friend Matrix<T> operator-(const DiagonalMatrix& a, const Matrix<T>& b) {
        auto out = b * -1.0;
        for (std::size_t i = 0; i < a.nrows(); ++i) out(i, i) = a(i, i) + out(i, i);
        return out;
    }

    friend auto operator-(const DiagonalMatrix& a, const DiagonalMatrix& b) {
        if (a.size() == 0)
            return b;
        else if (b.size() == 0)
            return a;
        return zip([](auto x, auto y) { return x - y; }, a, b);
    }

    friend auto operator-(const Matrix<T>& a, const DiagonalMatrix& b) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i) out(i, i) -= b(i, i);
        return out;
    }

    friend auto operator*(const DiagonalMatrix& a, const DiagonalMatrix& b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        auto out = DiagonalMatrix(a.nrows(), b.ncols());
        for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] * b[i];
        return out;
    }
    friend auto inv(const DiagonalMatrix& a) {
        return apply([](auto const& a) { return inv(a); }, a);
    }

    friend auto tr(const DiagonalMatrix& a) { return a; }

    friend auto diag(const DiagonalMatrix& a) {
        Matrix out(a.size(), 1, false);
        for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i];
        return out;
    }

    friend Maybe_error<double> logdet(const DiagonalMatrix& x) {
        double out = 0;
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            if (x(i, i) > 0)
                out += std::log(x(i, i));
            else
                return error_message("logdet of negative value at (" + std::to_string(i) + "," +
                                     std::to_string(i) + ")");
        }
        return out;
    }

    void set(std::size_t i, T x) { (*this)[i] = x; }

    void set(std::size_t i, std::size_t j, T x) {
        assert(i == j);
        (*this)[i] = x;
    }

    template <class F>
    friend auto apply(F&& f, DiagonalMatrix&& x) {
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(x[i]);
        return x;
    }

    friend auto xtAx(const Matrix<T>& x, const DiagonalMatrix& A) {
        assert(x.nrows() == A.nrows());
        auto out = T{};
        for (std::size_t i = 0; i < A.nrows(); ++i) out += x[i] * A(i, i) * x[i];
        return out;
    }

    friend auto xAxt(const Matrix<T>& x, const DiagonalMatrix& A) {
        assert(x.ncols() == A.nrows());
        auto out = T{};
        for (std::size_t i = 0; i < A.nrows(); ++i) out += x[i] * A(i, i) * x[i];
        return out;
    }

    template <class F>
        requires(!is_Maybe_error<std::invoke_result_t<F, T>>)
    friend auto apply(F&& f, DiagonalMatrix const& x) {
        auto out = DiagonalMatrix(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = f(x[i]);
        return out;
    }

    template <class F>
        requires(is_Maybe_error<std::invoke_result_t<F, T>>)
    friend Maybe_error<DiagonalMatrix<T>> apply(F&& f, DiagonalMatrix const& x) {
        auto out = DiagonalMatrix(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) {
            auto v = f(x[i]);
            if (v)
                out[i] = std::move(v.value());
            else
                return error_message(v.error()() + " at the " + std::to_string(i) + "th row");
        }
        return out;
    }

    template <class F>
    friend auto reduce(F&& f, const DiagonalMatrix& x) {
        auto cum = x[0];
        for (std::size_t i = 1; i < x.size(); ++i) cum = f(cum, x[i]);
        return cum;
    }

    friend auto& operator<<(std::ostream& os, const DiagonalMatrix& x) {
        os << "Diagonal matrix "
           << "nrows: " << x.nrows() << " ncols: " << x.ncols() << "\n";
        for (std::size_t i = 0; i < x.size(); ++i) os << std::setw(cell_width()) << x[i] << " ";
        os << "\n";

        return os;
    }
};

inline auto elemMult(const DiagonalMatrix<double>& x, const DiagonalMatrix<double>& y) {
    auto out = x;
    for (std::size_t i = 0; i < x.size(); ++i) {
        out[i] = x[i] * y[i];
    }
    return out;
}

template <class Matrix>
    requires(Matrix::is_Matrix)
std::ostream& print(std::ostream& os, const Matrix& x) {
    return os << x;
}

template <class T>
DiagPosDetMatrix<T> XTX(const DiagonalMatrix<T>& a);

template <class T>
class DiagPosDetMatrix : public DiagonalMatrix<T> {
   private:
    DiagPosDetMatrix(DiagonalMatrix<T>&& x) : base_type{x} {}

   public:
    static constexpr bool is_Matrix = true;
    static constexpr bool is_Symmetric = true;
    static constexpr bool is_Diagonal = true;

    using base_type = DiagonalMatrix<T>;
    explicit DiagPosDetMatrix(std::size_t _nrows, std::size_t) : base_type{_nrows, _nrows, false} {}
    explicit DiagPosDetMatrix(std::size_t _nrows, std::size_t, T value)
        : base_type{_nrows, _nrows, value} {}
    explicit DiagPosDetMatrix(std::initializer_list<T>&& a) : base_type(std::move(a)) {}

    explicit DiagPosDetMatrix(const Matrix<T>& a) : base_type(a) {}

    explicit DiagPosDetMatrix(Matrix<T>&& a) : base_type(std::move(a)) {}

    explicit DiagPosDetMatrix(const std::vector<T>& a) : base_type(a) {}

    using base_type::operator[];
    using base_type::operator();
    DiagPosDetMatrix() {}
    ~DiagPosDetMatrix() {}

    template <class T0>
    friend DiagPosDetMatrix<T0> XXT(const DiagonalMatrix<T0>& a);

    friend DiagPosDetMatrix<T> XTX(const DiagonalMatrix<T>& a) {
        return DiagPosDetMatrix<T>(a * a);
    }

    friend auto operator*(const DiagPosDetMatrix& a, const DiagPosDetMatrix& b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        auto out = DiagPosDetMatrix(a.nrows(), b.ncols());
        for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] * b[i];
        return out;
    }
    friend auto operator*(const DiagPosDetMatrix& a, double b) {
        auto out = DiagPosDetMatrix(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] * b;
        return out;
    }

    friend auto& tr(const DiagPosDetMatrix& a) { return a; }

    friend Maybe_error<DiagPosDetMatrix> cholesky(const DiagPosDetMatrix& a) {
        DiagPosDetMatrix out(a.nrows(), a.nrows());
        for (std::size_t i = 0; i < out.nrows(); ++i)
            if (a(i, i) > 0)
                out[i] = std::sqrt(a(i, i));
            else if (a(i, i) == 0)
                return error_message(std::string("zero value at") + std::to_string(i));
            else
                return error_message("sqrt of negative value " + std::to_string(a(i, i)));
        return out;
    }
};

template <template <class> class T2_Matrix, template <class> class T3_Matrix>
    requires(is_Matrix_v<T2_Matrix<double>> && is_Matrix_v<T3_Matrix<double>>)
auto operator*(const Matrix<T2_Matrix<double>>& a, const T3_Matrix<double>& b) {
    using S = std::decay_t<decltype(a[0] * b)>;
    auto out = Matrix<S>(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] * b;
    return out;
}
template <template <class> class T2_Matrix, template <class> class T3_Matrix>
    requires(is_Matrix_v<T2_Matrix<double>> && is_Matrix_v<T3_Matrix<double>>)
auto operator*(const T3_Matrix<double>& b, const Matrix<T2_Matrix<double>>& a) {
    using S = std::decay_t<decltype(a[0] * b)>;
    auto out = Matrix<S>(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = b * a[i];
    return out;
}

template <class T>
void set(SymmetricMatrix<T>& m, std::size_t i, std::size_t j, const T& x) {
    m.set(i, j, x);
}
template <class T>
void set(Matrix<T>& m, std::size_t i, std::size_t j, const T& x) {
    m.set(i, j, x);
}

template <class F, class T, template <class> class aMatrix>
    requires(aMatrix<T>::is_Symmetric && !aMatrix<T>::is_Diagonal)
auto applyMap(F&& f, aMatrix<T> const& a) {
    using S = std::decay_t<std::invoke_result_t<F, T>>;
    aMatrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i; j < x.ncols(); ++j) set(x, i, j, f(a(i, j)));
    return x;
}

template <class F, class T, template <class> class aMatrix>
    requires(aMatrix<T>::is_Symmetric && !aMatrix<T>::is_Diagonal)
auto applyMap(F&& f, aMatrix<T>& a) {
    using S = std::decay_t<std::invoke_result_t<F, T&>>;
    aMatrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i; j < x.ncols(); ++j) set(x, i, j, f(a(i, j)));
    return x;
}

template <class F, class T, template <class> class aMatrix>
    requires(aMatrix<T>::is_Diagonal)
auto applyMap(F&& f, aMatrix<T> const& a) {
    using S = std::decay_t<std::invoke_result_t<F, T>>;
    aMatrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a(i, i));
    return x;
}

template <class F, class T, template <class> class aMatrix>
    requires(aMatrix<T>::is_Diagonal)
auto applyMap(F&& f, aMatrix<T>& a) {
    using S = std::decay_t<std::invoke_result_t<F, T&>>;
    aMatrix<S> x(a.nrows(), a.ncols());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = f(a(i, i));
    return x;
}

template <class F, class T, template <class> class aMatrix>
    requires(aMatrix<T>::is_Symmetric)
auto operator*(const aMatrix<T>& a, double b) {
    return applyMap([b](auto const& x) { return x * b; }, a);
}

template <class T>
auto IdM(std::size_t ndim) {
    return DiagPosDetMatrix<T>(ndim, ndim, T{1.0});
}

template <class T>
auto inv_from_chol(const DownTrianMatrix<T>& x) {
    return lapack::Lapack_LT_Cholesky_inv(x);
}

template <class T>
auto inv_from_chol(const UpTrianMatrix<T>& x) {
    return lapack::Lapack_UT_Cholesky_inv(x);
}

template <class T>
Maybe_error<DiagPosDetMatrix<T>> inv_from_chol(const DiagPosDetMatrix<T>& x) {
    DiagPosDetMatrix<T> out(x.nrows(), x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        if (x(i, i) > 0)
            out[i] = std::pow(x(i, i), -2);
        else
            return error_message("error inverting chol");
    return out;
}

template <class T>
Maybe_error<DiagPosDetMatrix<T>> inv(const DiagPosDetMatrix<T>& x) {
    DiagPosDetMatrix<T> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        if (out(i, i) > 0)
            out[i] = std::pow(x(i, i), -1);
        else
            return error_message("error inverting chol");
    return out;
}

template <template <class> class Matrix, class T>
    requires Matrix<T>::is_Matrix
double getvalue(const Matrix<T>& x) {
    assert(x.size() == 1);
    return x[0];
}

template <class T>
auto diag(const Matrix<T>& a) {
    return DiagonalMatrix<T>(a);
}
template <class T>
auto diag(Matrix<T>&& a) {
    return DiagonalMatrix<T>(std::move(a));
}

template <class T>
auto diag(std::initializer_list<T>&& a) {
    return DiagonalMatrix<T>(std::move(a));
}

template <class T>
auto diagpos(std::initializer_list<T>&& a) {
    return DiagPosDetMatrix<T>(std::move(a));
}

template <class T>
auto diagpos(const Matrix<T>& a) {
    return DiagPosDetMatrix<T>(std::move(a));
}

template <class T>
auto qr(const Matrix<T>& a) {
    return lapack::Lapack_QR(a);
}

inline auto eigs(const Matrix<double>& x, bool does_permutations = false,
                 bool does_diagonal_scaling = false,
                 bool computes_eigenvalues_condition_numbers = false,
                 bool computes_eigenvectors_condition_numbers = false) {
    return lapack::Lapack_EigenSystem(x, does_permutations, does_diagonal_scaling,
                                      computes_eigenvalues_condition_numbers,
                                      computes_eigenvectors_condition_numbers);
}

inline auto eigs(const SymmetricMatrix<double>& x) {
    return lapack::Lapack_Symm_EigenSystem(x);
}

// Enforce CTMC-generator conventions on eigendecomposition (primitive):
// - deterministically place the zero eigenvalue at index 0
// - apply the corresponding column permutation to VR and VL
// - biorthogonal rescaling of (VL, VR) so that v_i^T u_i = 1 for all i
inline Maybe_error<std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>>
    eig_enforce_q_mode(const Matrix<double>& /*Q*/,
                       const std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>& e,
                       double tol = 1e-12) {
    auto [VR, L, VL] = e;
    if (VR.ncols() == 0)
        return std::tuple(VR, L, VL);

    const std::size_t n = VR.ncols();

    // Identify (near-)zero eigenvalue cluster indices
    std::vector<std::size_t> zero_idx;
    zero_idx.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (std::abs(L[i]) <= tol)
            zero_idx.push_back(i);

    if (zero_idx.empty()) {
        // Fallback: pick smallest-magnitude eigenvalue as the zero-mode
        std::size_t k0 = 0;
        double amin = std::abs(L[0]);
        for (std::size_t i = 1; i < n; ++i)
            if (std::abs(L[i]) < amin) {
                amin = std::abs(L[i]);
                k0 = i;
            }
        zero_idx.push_back(k0);
    }

    // Build permutation bringing the entire zero cluster to the front, preserving relative order
    Matrix<double> Per(n, n, 0.0);
    std::vector<std::size_t> order(n);
    std::size_t w = 0;
    for (auto k : zero_idx) order[w++] = k;
    for (std::size_t i = 0; i < n; ++i) {
        bool in_zero = std::find(zero_idx.begin(), zero_idx.end(), i) != zero_idx.end();
        if (!in_zero)
            order[w++] = i;
    }
    for (std::size_t i = 0; i < n; ++i) Per(i, order[i]) = 1.0;
    auto PerT = tr(Per);

    // Apply permutation
    auto VRp = VR * PerT;
    auto VLp = VL * PerT;
    DiagonalMatrix<double> Lp(L.nrows(), L.ncols());
    for (std::size_t i = 0; i < n; ++i) Lp[i] = L[order[i]];

    // Biorthogonal diagonal rescaling: v_i^T u_i = 1
    for (std::size_t i = 0; i < n; ++i) {
        auto ui = VRp(":", i);
        auto vTi = tr(VLp(":", i));
        double s = getvalue(vTi * ui);
        if (std::abs(s) <= tol)
            return error_message("eig_enforce_q_mode: singular (v_i^T u_i == 0)");
        for (std::size_t r = 0; r < VLp.nrows(); ++r) VLp(r, i) = VLp(r, i) / s;
    }
    return std::tuple(VRp, Lp, VLp);
}

struct Block {
    std::size_t begin{};
    std::size_t end{};
};

inline std::vector<Block> cluster_blocks(const DiagonalMatrix<double>& L, double tol) {
    std::vector<Block> blocks;
    if (L.size() == 0)
        return blocks;

    std::size_t i = 0;
    const std::size_t n = L.size();
    while (i < n) {
        std::size_t j = i + 1;
        for (; j < n; ++j) {
            const double gap = std::abs(L[j] - L[j - 1]);
            const double scale = 1.0 + std::abs(L[j]) + std::abs(L[j - 1]);
            if (gap / scale > tol)
                break;
        }
        blocks.push_back(Block{i, j});
        i = j;
    }
    return blocks;
}

inline Matrix<double> select_columns(std::size_t ncols, const Block& b) {
    Matrix<double> S(ncols, b.end - b.begin, 0.0);
    for (std::size_t j = 0; j < b.end - b.begin; ++j) S(b.begin + j, j) = 1.0;
    return S;
}

inline Matrix<double> select_rows(std::size_t nrows, const Block& b) {
    Matrix<double> S(b.end - b.begin, nrows, 0.0);
    for (std::size_t i = 0; i < b.end - b.begin; ++i) S(i, b.begin + i) = 1.0;
    return S;
}

inline DiagonalMatrix<double> select_diag(std::size_t n, const Block& b) {
    DiagonalMatrix<double> D(n, n);
    for (std::size_t i = 0; i < n; ++i) D[i] = 0.0;
    for (std::size_t i = b.begin; i < b.end; ++i) D[i] = 1.0;
    return D;
}

template <class MatLike, class DiagLike>
struct SpectralProjectors {
    MatLike V;
    MatLike W;
    DiagLike L;
    std::vector<Block> blocks;
};

template <class MatLike, class DiagLike>
SpectralProjectors<MatLike, DiagLike> make_projectors(MatLike V, MatLike W, DiagLike L,
                                                      double tol) {
    SpectralProjectors<MatLike, DiagLike> sp{std::move(V), std::move(W), std::move(L), {}};
    auto Lprim = primitive(sp.L);
    sp.blocks = cluster_blocks(Lprim, tol);
    return sp;
}

inline auto XXT(const Matrix<double>& a) {
    return lapack::Lapack_Product_Self_Transpose(a, false);
}

template <class T>
auto XXT(const DiagonalMatrix<T>& a) {
    return DiagPosDetMatrix<T>(a * a);
}

inline auto XTX(const Matrix<double>& a) {
    if (a.size() == 0)
        return SymPosDefMatrix<double>{};
    else
        return lapack::Lapack_Product_Self_Transpose(a, true);
}
inline auto XTX(const Matrix<std::size_t>& a) {
    if (a.size() == 0)
        return SymPosDefMatrix<std::size_t>{};
    else {
        SymPosDefMatrix<std::size_t> out(a.ncols(), a.ncols(), 0ul);
        for (std::size_t i = 0; i < out.nrows(); ++i)
            for (std::size_t j = i; j < out.ncols(); ++j)
                for (std::size_t k = 0; k < a.nrows(); ++k)
                    out.set(i, j, out(i, j) + a(k, i) * a(k, j));
        return out;
    }
}

inline auto sqr_X(const Matrix<double>& a) {
    if (a.ncols() < a.nrows())
        return XXT(a);
    else
        return XTX(a);
}

inline auto sqr_X(double a) {
    return a * a;
}

template <template <class> class Matrix>
    requires(!Matrix<double>::is_Symmetric)
auto X_plus_XT(const Matrix<double>& x) {
    assert(x.ncols() == x.nrows());
    SymmetricMatrix<double> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); i++)
        for (std::size_t j = 0; j < i + 1; ++j) out.set(j, i, x(i, j) + x(j, i));
    return out;
}

template <template <class> class Matrix>
    requires(Matrix<double>::is_Symmetric)
auto X_plus_XT(const Matrix<double>& x) {
    return x * 2;
}

template <class T>
SymPosDefMatrix<T> AT_D_A(const Matrix<T>& A, const DiagPosDetMatrix<T>& D) {
    SymPosDefMatrix<T> out(A.ncols(), A.ncols(), T{});
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = i; j < out.ncols(); ++j)
            for (std::size_t k = 0; k < D.ncols(); ++k)
                out.set(i, j, out(i, j) + A(k, i) * D(k, k) * A(k, j));
    return out;
}

template <class T>
SymPosDefMatrix<T> AT_D_A_Safe(const Matrix<T>& A, const DiagPosDetMatrix<T>& D) {
    SymPosDefMatrix<T> out(A.ncols(), A.ncols(), T{});
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = i; j < out.ncols(); ++j)
            for (std::size_t k = 0; k < D.ncols(); ++k)
                if (std::isfinite(A(k, i) * D(k, k) * A(k, j)))
                    out.set(i, j, out(i, j) + A(k, i) * D(k, k) * A(k, j));
    return out;
}

template <class T>
SymPosDefMatrix<T> AT_B_A(const Matrix<T>& A, const SymmetricMatrix<T>& B) {
    SymPosDefMatrix<T> out(B.nrows(), 0.0);
    if (B.size() == 0)
        return out;
    else {
        auto AB = B * A;
        for (std::size_t i = 0; i < out.nrows(); ++i)
            for (std::size_t j = i; j < out.ncols(); ++j)
                for (std::size_t k = 0; k < A.nrows(); ++k)
                    out.set(i, j, out(i, j) + A(k, i) * AB(k, j));
        return out;
    }
}

template <class Matrix>
double Trace(Matrix const& x) {
    double out = 0;
    for (std::size_t i = 0; i < std::min(x.nrows(), x.ncols()); ++i) out += x(i, i);
    return out;
}

inline double norm_1(const double x) {
    return std::abs(x);
}

/**
    Maximum Value of the Sum of the absolute values in each row
   */
template <typename Matrix>
    requires Matrix::is_Matrix
double norm_1(const Matrix& x) {
    double n = 0;
    for (size_t i = 0; i < x.ncols(); ++i) {
        double sum = 0;
        for (size_t j = 0; j < x.nrows(); ++j) sum += norm_1(x(j, i));
        n = std::max(n, sum);
    }
    return n;
}

/**
    Maximum Value of the Sum of the absolute values in each row
   */
template <typename Matrix>
    requires Matrix::is_Matrix
double norm_inf(const Matrix& x) {
    double n(0);
    for (size_t i = 0; i < x.nrows(); ++i) {
        double sum(0);
        for (size_t j = 0; j < x.ncols(); ++j)
            if (x(i, j) > 0)
                sum += x(i, j);
            else
                sum -= x(i, j);

        n = std::max(n, sum);
    }
    return n;
}

inline Maybe_error<bool> test_equality(double x, double y, double eps) {
    if (std::abs(x - y) / (std::abs(x) + std::abs(y) + 1.0) > eps)
        return error_message(
            ToString(x) + " is not equal to " + ToString(y) +
            "|x-y|/(|x|+|y|) = " + ToString(std::abs(x - y) / (std::abs(x) + std::abs(y) + 1.0)) +
            " is greater than " + ToString(eps));
    else
        return true;
}

inline Maybe_error<bool> test_equality(std::integral auto x, std::integral auto y, double) {
    if (x != y)
        return error_message(ToString(x) + " is not equal to " + ToString(y));
    else
        return true;
}

template <template <class> class aMatrix>
Maybe_error<bool> test_equality(const aMatrix<double>& one, const aMatrix<double>& two,
                                double eps) {
    auto out = Maybe_error<bool>(true);
    for (std::size_t i = 0; i < one.nrows(); ++i)
        for (std::size_t j = 0; j < one.ncols(); ++j) {
            auto test = test_equality(one(i, j), two(i, j), eps);
            if (!test)
                out = error_message(out.error()() + "\n at (" + std::to_string(i) + "," +
                                    std::to_string(j) + "): " + test.error()());
        }
    if (!out) {
        std::stringstream ss;
        ss << "\nMatrix differs\n";
        ss << one << "\n differs from \n" << two;
        ss << out.error()();
        return error_message(ss.str());
    } else
        return out;
}

template <typename Matrix>
    requires Matrix::is_Matrix
bool operator<(const Matrix& x, const Matrix& y) {
    if (x.size() < y.size())
        return true;
    else if (y.size() < x.size())
        return false;
    else
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i] < y[i])
                return true;
            else if (y[i] < x[i])
                return false;
    return false;
}

template <typename Matrix>
    requires Matrix::is_Matrix
bool operator==(const Matrix& x, const Matrix& y) {
    if (x.size() != y.size())
        return false;
    else
        for (std::size_t i = 0; i < x.size(); ++i)
            if (!(x[i] == y[i]))
                return false;
    return true;
}

#endif  // MATRIX_H
