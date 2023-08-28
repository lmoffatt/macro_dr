#ifndef MATRIX_H
#define MATRIX_H

#include "maybe_error.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <iomanip>


template <class> class Matrix;
template <class> class SubMatrix;

template <class> class SymmetricMatrix;
template <class> class SymPosDefMatrix;
template <class> class UpTrianMatrix;
template <class> class DownTrianMatrix;
template <class> class DiagonalMatrix;



template<std::size_t I,typename ...Args >
void write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i)
{
    if (i==I)
        os<<std::get<I>(tu);
}

template<std::size_t...Is,typename ...Args>
std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i, std::index_sequence<Is...>)
{
    (write_tuple_i<Is>(os,tu, i),...);
    return os;
}
template<typename ...Args> std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i)
{
    return write_tuple_i(os,tu,i,std::index_sequence_for<Args...>());
}





namespace lapack {
Matrix<double> Lapack_Full_Product(const Matrix<double> &x,
                                   const Matrix<double> &y, bool transpose_x,
                                   bool transpose_y, double alpha = 1.0,
                                   double beta = 0.0);
Matrix<double> Lapack_Sym_Product(const SymmetricMatrix<double> &x,
                                  const Matrix<double> &y,
                                  bool first_symmetric_c, double alpha = 1.0,
                                  double beta = 0.0);

Matrix<double> Lapack_Triang_Product(const Matrix<double> &a,
                                     const Matrix<double> &b,
                                     bool up_triangular_in_c,
                                     bool triangular_first_in_c,
                                     bool transpose_A_in_c, bool ones_in_diag,
                                     double alpha = 1);

SymPosDefMatrix<double>
Lapack_Product_Self_Transpose(const Matrix<double> &a,
                              bool first_transposed_in_c, char UPLO_in_c = 'U',
                              double alpha = 1, double beta = 0);

std::pair<Matrix<double>, Matrix<double>> Lapack_QR(const Matrix<double> &x);
Maybe_error<Matrix<double>> Lapack_Full_inv(const Matrix<double> &a);

Maybe_error<SymmetricMatrix<double>>
Lapack_Symm_inv(const SymmetricMatrix<double> &a);

Maybe_error<SymPosDefMatrix<double>>
Lapack_SymmPosDef_inv(const SymPosDefMatrix<double> &a);
Maybe_error<DownTrianMatrix<double>>
Lapack_chol(const SymPosDefMatrix<double> &x);

template <class T>
Maybe_error<SymPosDefMatrix<T>>
Lapack_UT_Cholesky_inv(const UpTrianMatrix<T> &x);
template <class T>
Maybe_error<SymPosDefMatrix<T>>
Lapack_LT_Cholesky_inv(const DownTrianMatrix<T> &x);

Maybe_error<DownTrianMatrix<double>>
Lapack_LT_inv(const DownTrianMatrix<double> &x, bool ones_in_diag);

Maybe_error<UpTrianMatrix<double>> Lapack_UT_inv(const UpTrianMatrix<double> &x,
                                                 bool ones_in_diag);

Maybe_error<
    std::tuple<Matrix<double>, DiagonalMatrix<double>, Matrix<double>>>
Lapack_EigenSystem(const Matrix<double> &x, bool does_permutations=true,
                   bool does_diagonal_scaling=true,
                   bool computes_eigenvalues_condition_numbers=false,
                   bool computes_eigenvectors_condition_numbers=false);
}//namespace lapack

template <class T> class DiagPosDetMatrix;

template <class T> class Matrix {
private:
    std::size_t size_ = 0;
    std::size_t nrows_ = 0;
    std::size_t ncols_ = 0;
    T *x_ = nullptr;
    
public:
    static int cell_width(){return 12;}
    explicit Matrix(){}
    explicit Matrix(std::size_t _nrows, std::size_t _ncols,
                    bool initialize = true)
        : size_{_nrows * _ncols}, nrows_{_nrows}, ncols_{_ncols},
        x_{initialize ? new T[size_]() : new T[size_]} {}
    explicit Matrix(std::size_t _nrows, std::size_t _ncols, T value)
        : size_{_nrows * _ncols}, nrows_{_nrows}, ncols_{_ncols},
        x_{new T[size_]} {
        for (std::size_t i = 0; i < size_; ++i)
            x_[i] = value;
    }
    
    Matrix(const Matrix &x)
        : size_{x.size()}, nrows_{x.nrows()}, ncols_{x.ncols()},
        x_{new T[nrows_ * ncols_]} {
        for (std::size_t i = 0; i < size_; ++i)
            x_[i] = x[i];
    }
    Matrix(Matrix &&x)
        : size_{x.size()}, nrows_{x.nrows()}, ncols_{x.ncols()}, x_{x.x_} {
        x.x_ = nullptr;
    }
    
    friend class DiagonalMatrix<T>;
    
    Matrix &operator=(const Matrix &x) {
        if (&x != this) {
            if (size() != x.size()) {
                size_ = x.size();
                nrows_ = x.nrows();
                ncols_ = x.ncols();
                delete[] x_;
                x_ = new T[x.size()];
            }
            
            for (std::size_t i = 0; i < size_; ++i)
                x_[i] = x[i];
        }
        return *this;
    }
    
    Matrix &operator=(Matrix &&x) {
        if (&x != this) {
            size_ = x.size();
            nrows_ = x.nrows();
            ncols_ = x.ncols();
            std::swap(x_ , x.x_);
            
        }
        return *this;
    }
    
    
    auto &operator[](std::size_t i) { assert(i<size_); return x_[i]; }
    auto &operator[](std::size_t i) const {assert(i<size_);  return x_[i]; }
    auto &operator()(std::size_t i) { assert(i<size_); return x_[i]; }
    auto &operator()(std::size_t i) const {assert(i<size_);  return x_[i]; }
    auto &operator()(std::size_t i, std::size_t j) { assert(i<nrows()&& j<ncols()); return x_[i * ncols_ + j]; }
    auto &operator()(std::size_t i, std::size_t j) const {
        assert((i<nrows())&& (j<ncols()));
        return x_[i * ncols_ + j];
    }
    auto ncols() const { return ncols_; }
    auto nrows() const { return nrows_; }
    auto size() const { return size_; }
    ~Matrix() {
        //  std::cerr<<"release "<<size()<<"\n";
        delete[] x_; }
    
    friend auto operator*(const Matrix &a, const Matrix &b) {
        return lapack::Lapack_Full_Product(a, b, false, false);
    }
    
    friend auto operator*(const Matrix &a, double b) {
        return apply([&b](auto x) { return x * b; }, a);
    }
    
    friend auto operator/(const Matrix &a, double b) {
        return apply([&b](auto x) { return x / b; }, a);
    }
    
    
    friend auto operator*(double b, const Matrix &a) {
        return apply([&b](auto x) { return x * b; }, a);
    }
    
    friend auto operator+(const Matrix &a, const Matrix &b) {
        return zip([](auto x, auto y) { return x + y; }, a, b);
    }
    
    friend auto operator-(const Matrix &a, const Matrix &b) {
        return zip([](auto x, auto y) { return x - y; }, a, b);
    }
    
    friend auto inv(const Matrix &a) { return lapack::Lapack_Full_inv(a); }
    
    friend auto tr(const Matrix &a) {
        Matrix out(a.ncols(), a.nrows(), false);
        for (std::size_t i = 0; i < out.nrows(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j)
                out(i, j) = a(j, i);
        return out;
    }
    
    template <class F> friend auto apply(F &&f, Matrix &&x) {
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(x[i]);
        return x;
    }
    
    template <class F>
        requires(!is_Maybe_error<std::invoke_result_t<F, T>>)
    friend auto apply(F &&f, Matrix const &a) {
        Matrix x(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(a[i]);
        return x;
    }
    
    template <class F>
        requires(is_Maybe_error<std::invoke_result_t<F, T>>)
    friend Maybe_error<Matrix<T>> apply(F &&f, Matrix const &a) {
        Matrix x(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = 0; j < x.ncols(); ++j) {
                auto v = f(a(i, j));
                if (v)
                    x(i, j) = std::move(v.value());
                else {
                    return v.error() + " at cell (" + std::to_string(i) + "," +
                           std::to_string(j) + ")";
                }
            }
        return x;
    }
    
    friend bool same_dimensions(const Matrix &x, const Matrix& y) {
        return x.size() == y.size() && x.nrows() == y.nrows() &&
               x.ncols() == y.ncols();
    }
    
    friend T xtx(const Matrix &x) {
        assert(x.ncols() == 1);
        auto out = T{};
        for (std::size_t i = 0; i < x.nrows(); ++i)
            out += x[i] * x[i];
        return out;
    }
    
    template <class F> friend auto zip(F &&f, const Matrix &x, const Matrix &y) {
        assert(same_dimensions(x, y) && "same size");
        
        Matrix out(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i)
            out[i] = f(x[i], y[i]);
        return out;
    }
    
    template <class F> friend auto reduce(F &&f, const Matrix &x) {
        auto cum = x[0];
        for (std::size_t i = 1; i < x.size(); ++i)
            cum = f(cum, x[i]);
        return cum;
    }
    
    template <class F> friend auto reduce_ij(F &&f, const Matrix &x, T init) {
        auto cum = init;
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = 0; j < x.ncols(); ++j)
                cum = f(cum, i, j, x(i, j));
        return cum;
    }
    
    friend auto &operator<<(std::ostream &os, const Matrix &x) {
        os << "\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os << std::setw(cell_width())<<x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};





template <class T> class SymmetricMatrix : public Matrix<T> {
    using base_type = Matrix<T>;
    
protected:
    explicit SymmetricMatrix(base_type &&x) : base_type{std::move(x)} {}
    explicit SymmetricMatrix(base_type const &x) : base_type{x} {}
    
public:
    operator Matrix<T> const &() { return static_cast<Matrix<T> const &>(*this); }
    
    explicit SymmetricMatrix(std::size_t _nrows)
        : base_type(_nrows, _nrows, false) {}
    explicit SymmetricMatrix(std::size_t _nrows, T value)
        : base_type(_nrows, _nrows, value) {}
    
    SymmetricMatrix(const SymmetricMatrix &x) : base_type(x) {}
    SymmetricMatrix(SymmetricMatrix &&x) : base_type(std::move(x)) {}
    
    SymmetricMatrix &operator=(const SymmetricMatrix &x) {
        static_cast<base_type &>(*this) = static_cast<base_type const &>(x);
        return *this;
    }
    SymmetricMatrix &operator=(SymmetricMatrix &&x) {
        static_cast<base_type &>(*this) = std::move(static_cast<base_type &>(x));
        return *this;
    }
    
    void set(std::size_t i, std::size_t j, double x) {
        base_type::operator()(i, j) = x;
        if (i != j)
            base_type::operator()(j, i) = x;
    }
    auto &operator()(std::size_t i, std::size_t j) const {
        return base_type::operator()(i, j);
    }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~SymmetricMatrix() {}
    
    friend auto operator*(const SymmetricMatrix &a, const SymmetricMatrix &b) {
        return lapack::Lapack_Sym_Product(a, b, true);
    }
    friend auto operator*(const SymmetricMatrix &a, const Matrix<double> &b) {
        return lapack::Lapack_Sym_Product(a, b, true);
    }
    friend auto operator*(const Matrix<double> &b, const SymmetricMatrix &a) {
        return lapack::Lapack_Sym_Product(a, b, false);
    }
    
    friend SymmetricMatrix operator*(const SymmetricMatrix &a, double b) {
        return apply([&b](auto x) { return x * b; }, a);
    }
    
    friend SymmetricMatrix operator*(double b, const SymmetricMatrix &a) {
        return apply([&b](auto x) { return x * b; }, a);
    }
    
    friend auto operator+(const SymmetricMatrix &a, const SymmetricMatrix &b) {
        return zip([](auto x, auto y) { return x + y; }, a, b);
    }
    
    friend auto operator-(const SymmetricMatrix &a, const SymmetricMatrix &b) {
        return zip([](auto x, auto y) { return x - y; }, a, b);
    }
    
    friend auto inv(const SymmetricMatrix &a) {
        return lapack::Lapack_Symm_inv(a);
    }
    
    friend auto &tr(const SymmetricMatrix &a) { return a; }
    
    template <class F> friend auto apply(F &&f, SymmetricMatrix &&x) {
        for (std::size_t i = 0; i < x.size(); ++i)
            for (std::size_t j = i; j < x.size(); ++j)
                x.set(i, j, f(x(i, j)));
        return x;
    }
    
    template <class F>
    friend auto zip(F &&f, const SymmetricMatrix &x, const SymmetricMatrix &y) {
        assert(x.nrows() == y.nrows() && "same size");
        
        SymmetricMatrix out(x.nrows(), false);
        for (std::size_t i = 0; i < x.nrows(); ++i)
            for (std::size_t j = i; j < x.ncols(); ++j)
                out.set(i, j, f(x(i, j), y(i, j)));
        return out;
    }
};

template <class T> auto xtAx(const Matrix<T> &x, const SymmetricMatrix<T> &A) {
    assert(x.nrows() == A.nrows());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * x[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j)
            out += 2 * x[i] * A(i, j) * x[j];
    }
    return out;
}

template <class T> auto xAxt(const Matrix<T> &x, const SymmetricMatrix<T> &A) {
    assert(x.ncols() == A.nrows());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * x[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j)
            out += 2 * x[i] * A(i, j) * x[j];
    }
    return out;
}

template <class T> auto xAyt(const Matrix<T> &x, const SymmetricMatrix<T> &A, const Matrix<T>& y) {
    assert(x.ncols() == A.nrows());
    assert(y.ncols() == A.ncols());
    auto out = T{};
    for (std::size_t i = 0; i < A.nrows(); ++i) {
        out += x[i] * A(i, i) * x[i];
        for (std::size_t j = i + 1; j < A.ncols(); ++j)
            out += 2 * x[i] * A(i, j) * y[j];
    }
    return out;
}


template <class F, class T> auto apply(F &&f, const SymmetricMatrix<T> &x) {
    SymmetricMatrix<T> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i; j < x.ncols(); ++j)
            out.set(i, j, f(x(i, j)));
    return out;
}
template <class T> void copy_UT_to_LT(Matrix<T> &x) {
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = 0; j < i; ++j)
            x(i, j) = x(j, i);
}

template <class T> void copy_LT_to_UT(Matrix<T> &x) {
    for (std::size_t i = 0; i < x.nrows(); ++i)
        for (std::size_t j = i + 1; j < x.ncols(); ++j)
            x(i, j) = x(j, i);
}

template <class Ts> UpTrianMatrix<Ts> fill_LT_zeros(SymPosDefMatrix<Ts> &&x);

template <class Ts> DownTrianMatrix<Ts> fill_UT_zeros(SymPosDefMatrix<Ts> &&x);

template <class T> class UpTrianMatrix : public Matrix<T> {
    using base_type = Matrix<T>;
    
public:
    void fill_the_zeros() {
        for (std::size_t i = 0; i < nrows(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                base_type::operator()(i, j) = 0;
    }
    operator Matrix<T> const &() { return static_cast<Matrix<T> const &>(*this); }
    explicit UpTrianMatrix(base_type &&x) : base_type{std::move(x)} {}
    
    explicit UpTrianMatrix(std::size_t _nrows, bool initialize = true)
        : base_type(_nrows, _nrows, initialize) {}
    explicit UpTrianMatrix(std::size_t _nrows, T value)
        : base_type(_nrows, _nrows, value) {
        fill_the_zeros();
    }
    
    UpTrianMatrix(const UpTrianMatrix &x) : base_type(x) {}
    UpTrianMatrix(UpTrianMatrix &&x) : base_type(std::move(x)) {}
    
    UpTrianMatrix &operator=(const UpTrianMatrix &x) {
        static_cast<base_type &>(*this) = static_cast<base_type const &>(x);
        return *this;
    }
    
    template <class Ts> friend UpTrianMatrix fillzeros(SymPosDefMatrix<Ts> &&x);
    
    void set(std::size_t i, std::size_t j, double x) {
        assert(j >= i);
        base_type::operator()(i, j) = x;
    }
    auto &operator()(std::size_t i, std::size_t j) const {
        return base_type::operator()(i, j);
    }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~UpTrianMatrix() {}
    
    friend auto inv(const UpTrianMatrix &a) {
        return lapack::Lapack_UT_inv(a, false);
    }
    
    friend auto operator*(const UpTrianMatrix &a, const Matrix<T> &b) {
        return lapack::Lapack_Triang_Product(a, b, true, true, false, false);
    }
    friend auto operator*(const Matrix<T> &a, const UpTrianMatrix &b) {
        return lapack::Lapack_Triang_Product(b, a, true, false, false, false);
    }
    
    friend auto &operator<<(std::ostream &os, const UpTrianMatrix &x) {
        os << "Upper triang part\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os << std::setw(base_type::cell_width())<< x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};

template <class T> UpTrianMatrix<T> fill_LT_zeros(SymPosDefMatrix<T> &&x) {
    UpTrianMatrix<T> out(std::move(x));
    out.fill_the_zeros();
    return out;
}

template <class T> class DownTrianMatrix : public Matrix<T> {
    using base_type = Matrix<T>;
    
public:
    void fill_the_zeros() {
        for (std::size_t i = 0; i < nrows(); ++i)
            for (std::size_t j = i + 1; j < ncols(); ++j)
                base_type::operator()(i, j) = 0;
    }
    operator Matrix<T> const &() { return static_cast<Matrix<T> const &>(*this); }
    explicit DownTrianMatrix(base_type &&x) : base_type{std::move(x)} {}
    explicit DownTrianMatrix(base_type const &x) : base_type{x} {}
    
    explicit DownTrianMatrix(std::size_t _nrows, bool initialize = true)
        : base_type(_nrows, _nrows, initialize) {}
    explicit DownTrianMatrix(std::size_t _nrows, T value)
        : base_type(_nrows, _nrows, value) {
        fill_the_zeros();
    }
    
    DownTrianMatrix(const DownTrianMatrix &x) : base_type(x) {}
    DownTrianMatrix(DownTrianMatrix &&x) : base_type(std::move(x)) {}
    
    DownTrianMatrix &operator=(const DownTrianMatrix &x) {
        static_cast<base_type &>(*this) = static_cast<base_type const &>(x);
        return *this;
    }
    
    template <class Ts> friend DownTrianMatrix fillzeros(SymPosDefMatrix<Ts> &&x);
    
    void set(std::size_t i, std::size_t j, double x) {
        assert(j <= i);
        base_type::operator()(i, j) = x;
    }
    auto &operator()(std::size_t i, std::size_t j) const {
        return base_type::operator()(i, j);
    }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~DownTrianMatrix() {}
    
    friend auto inv(const DownTrianMatrix &a) {
        return lapack::Lapack_LT_inv(a, false);
    }
    
    friend auto operator*(const DownTrianMatrix &a, const Matrix<T> &b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }
    friend auto operator*(const DownTrianMatrix &a, const DownTrianMatrix &b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }
    friend auto operator*(const DownTrianMatrix &a, const UpTrianMatrix<T> &b) {
        return lapack::Lapack_Triang_Product(a, b, false, true, false, false);
    }
    
    friend auto operator*(const Matrix<T> &b, const DownTrianMatrix &a) {
        return lapack::Lapack_Triang_Product(a, b, false, false, false, false);
    }
    friend auto operator*(const UpTrianMatrix<T> &b, const DownTrianMatrix &a) {
        return lapack::Lapack_Triang_Product(a, b, false, false, false, false);
    }
    
    friend auto &operator<<(std::ostream &os, const DownTrianMatrix &x) {
        os << "Down triang part\n";
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            for (std::size_t j = 0; j < x.ncols(); ++j)
                os <<std::setw(base_type::cell_width())<< x(i, j) << " ";
            os << "\n";
        }
        return os;
    }
};

template <class T> DownTrianMatrix<T> fill_UT_zeros(SymPosDefMatrix<T> &&x) {
    DownTrianMatrix<T> out(std::move(x));
    out.fill_the_zeros();
    return out;
}

template <class T> auto tr(const UpTrianMatrix<T> &x) {
    return DownTrianMatrix<T>(tr(static_cast<Matrix<T> const &>(x)));
}

template <class T> auto tr(const DownTrianMatrix<T> &x) {
    return UpTrianMatrix<T>(tr(static_cast<Matrix<T> const &>(x)));
}

template <class T> Maybe_error<double> logdet(const UpTrianMatrix<T> &x) {
    double out = 0;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (x(i, i) > 0)
            out += std::log(x(i, i));
        else
            return error_message("logdet of negative value at (" + std::to_string(i) +
                               "," + std::to_string(i) + ")");
    }
    return out;
}

template <class T> Maybe_error<double> logdet(const DownTrianMatrix<T> &x) {
    double out = 0;
    for (std::size_t i = 0; i < x.nrows(); ++i) {
        if (x(i, i) > 0)
            out += std::log(x(i, i));
        else
            return error_message("logdet of negative value at (" + std::to_string(i) +
                               "," + std::to_string(i) + ")");
    }
    return out;
}

template <class T> class SymPosDefMatrix : public SymmetricMatrix<T> {
    using base_type = SymmetricMatrix<T>;
    
    explicit SymPosDefMatrix(base_type &&x) : base_type{std::move(x)} {}
    explicit SymPosDefMatrix(Matrix<T> &&x) : base_type{std::move(x)} {}
    
public:
    template <class K>
    friend Maybe_error<SymPosDefMatrix<K>>
    lapack::Lapack_UT_Cholesky_inv(const UpTrianMatrix<K> &x);
    
    template <class K>
    friend Maybe_error<SymPosDefMatrix<K>>
    lapack::Lapack_LT_Cholesky_inv(const DownTrianMatrix<K> &x);
    operator Matrix<T> const &() { return static_cast<Matrix<T> const &>(*this); }
    operator SymmetricMatrix<T> const &() {
        return static_cast<SymmetricMatrix<T> const &>(*this);
    }
    
    explicit SymPosDefMatrix(std::size_t _nrows, bool initialize = true)
        : base_type(_nrows, initialize) {}
    explicit SymPosDefMatrix(std::size_t _nrows, T value)
        : base_type(_nrows, value) {}
    
    static
        SymPosDefMatrix I_sware_it_is_possitive(base_type &&x) {return SymPosDefMatrix(std::move(x));}
    SymPosDefMatrix(const SymPosDefMatrix &x) : base_type(x) {}
    SymPosDefMatrix(SymPosDefMatrix &&x) : base_type(std::move(static_cast<base_type&>(x))) {}
    
    SymPosDefMatrix &operator=(const SymPosDefMatrix &x) {
        base_type::operator=(static_cast<base_type const &>(x));
        return *this;
    }
    
    void set(std::size_t i, std::size_t j, double x) {
        base_type::operator()(i, j) = x;
        if (i != j)
            base_type::operator()(j, i) = x;
    }
    auto &operator()(std::size_t i, std::size_t j) const {
        return base_type::operator()(i, j);
    }
    auto ncols() const { return base_type::ncols(); }
    auto nrows() const { return base_type::nrows(); }
    auto size() const { return base_type::size(); }
    ~SymPosDefMatrix() {}
    
    friend auto operator+(const SymPosDefMatrix &a, const SymPosDefMatrix &b) {
        return SymPosDefMatrix(zip([](auto x, auto y) { return x + y; }, a, b));
    }
    
    friend auto operator+(const SymPosDefMatrix &a,
                          const DiagPosDetMatrix<T> &b) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i)
            out(i, i) += b(i, i);
        return out;
    }
    
    friend auto operator+(const DiagPosDetMatrix<T> &b,const SymPosDefMatrix &a
                          ) {
        auto out = a;
        for (std::size_t i = 0; i < b.nrows(); ++i)
            out(i, i) += b(i, i);
        return out;
    }
    
    
    
    friend SymPosDefMatrix operator*(const SymPosDefMatrix &a, double b) {
        return SymPosDefMatrix(static_cast<SymmetricMatrix<T> const &>(a) * b);
    }
    
    friend SymPosDefMatrix operator*(double b, const SymPosDefMatrix &a) {
        return SymPosDefMatrix(b * static_cast<SymmetricMatrix<T> const &>(a));
    }
    
    friend auto inv(const SymPosDefMatrix &a) {
        return lapack::Lapack_SymmPosDef_inv(a);
    }
    
    friend auto tr(const SymPosDefMatrix &a) { return a; }
    
    friend Maybe_error<DownTrianMatrix<T>> cholesky(const SymPosDefMatrix &x) {
        return lapack::Lapack_chol(x);
    }
    
    friend Maybe_error<double> logdet(const SymPosDefMatrix &x) {
        auto chol = cholesky(x);
        if (chol)
            return 2.0*logdet(chol.value());
        else
            return Maybe_error<double>(chol.error() + " in calculation of logdet");
    }
};
template <class T>
Maybe_error<DownTrianMatrix<T>> cholesky(const SymPosDefMatrix<T> &x);

template <class T>
Maybe_error<DiagPosDetMatrix<T>> cholesky(const DiagPosDetMatrix<T> &x);

template <class T> class DiagonalMatrix {
private:
    std::size_t size_ = 0;
    std::size_t nrows_ = 0;
    std::size_t ncols_ = 0;
    T *x_ = nullptr;
    
public:
    explicit DiagonalMatrix(std::size_t _nrows, std::size_t _ncols,
                            bool initialize = true)
        : size_{std::min(_nrows, _ncols)}, nrows_{_nrows}, ncols_{_ncols},
        x_{initialize ? new T[size_]() : new T[size_]} {}
    explicit DiagonalMatrix(std::size_t _nrows, std::size_t _ncols, T value)
        : size_{std::min(_nrows, _ncols)}, nrows_{_nrows}, ncols_{_ncols},
        x_{new T[size_]} {
        for (std::size_t i = 0; i < size_; ++i)
            x_[i] = value;
    }
    //  explicit DiagonalMatrix(std::initializer_list<T> &&a)
    //      : size_{a.size()}, nrows_{a.size()}, ncols_{a.size()}, x_{new
    //      T[size_]} {
    //    std::copy(a.begin(), a.end(), x_);
    //  }
    
    explicit DiagonalMatrix(const Matrix<T> &a)
        : size_{(a.ncols() == 1 || a.nrows() == 1)
                    ? a.size()
                    : std::min(a.nrows(), a.ncols())},
        nrows_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.nrows()},
        ncols_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.ncols()},
        x_{new T[size_]} {
        if (a.ncols() == 1 || a.nrows() == 1) {
            for (std::size_t i = 0; i < size(); ++i)
                (*this)[i] = a[i];
        } else {
            for (std::size_t i = 0; i < size(); ++i)
                (*this)[i] = a(i, i);
        }
    }
    
    explicit DiagonalMatrix(Matrix<T> &&a)
        : size_{(a.ncols() == 1 || a.nrows() == 1)
                    ? a.size()
                    : std::min(a.nrows(), a.ncols())},
        nrows_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.nrows()},
        ncols_{(a.ncols() == 1 || a.nrows() == 1) ? a.size() : a.ncols()},
        x_{(a.ncols() == 1 || a.nrows() == 1) ? nullptr : new T[size_]} {
        if (a.ncols() == 1 || a.nrows() == 1) {
            std::swap(x_, a.x_);
            
        } else {
            for (std::size_t i = 0; i < size(); ++i)
                (*this)[i] = a(i, i);
        }
    }
    
    DiagonalMatrix(const DiagonalMatrix &x)
        : size_{x.size()}, nrows_{x.nrows()}, ncols_{x.ncols()},
        x_{new T[x.size()]} {
        for (std::size_t i = 0; i < size_; ++i)
            x_[i] = x[i];
    }
    DiagonalMatrix(DiagonalMatrix &&x)
        : size_{x.size()}, nrows_{x.nrows()}, ncols_{x.ncols()}, x_{x.x_} {
        x.x_ = nullptr;
    }
    
    DiagonalMatrix &operator=(const DiagonalMatrix &x) {
        if (size() != x.size()) {
            size_ = x.size();
            nrows_ = x.nrows();
            ncols_ = x.ncols();
            delete[] x_;
            x.x_ = new T[x.size()];
        }
        
        for (std::size_t i = 0; i < size_; ++i)
            x_[i] = x[i];
        return *this;
    }
    
    DiagonalMatrix &operator=(DiagonalMatrix &&x) {
        size_ = x.size();
        nrows_ = x.nrows();
        ncols_ = x.ncols();
        std::swap(x_, x.x_);
        return *this;
    }
    
    auto &operator[](std::size_t i) { return x_[i]; }
    auto &operator[](std::size_t i) const { return x_[i]; }
    auto operator()(std::size_t i, std::size_t j) const {
        if (i == j)
            return x_[i];
        else
            return T{};
    }
    
    auto ncols() const { return ncols_; }
    auto nrows() const { return nrows_; }
    auto size() const { return size_; }
    ~DiagonalMatrix() { delete[] x_; }
    static int cell_width(){return 12;}
    
    friend auto operator*(const DiagonalMatrix &a, const Matrix<double> &b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        Matrix<double> out(a.nrows(), b.ncols(), false);
        for (std::size_t i = 0; i < a.size(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j)
                out(i, j) = a[i] * b(i, j);
        for (std::size_t i = a.size(); i < out.nrows(); ++i)
            for (std::size_t j = 0; j < out.ncols(); ++j)
                out(i, j) = 0;
        return out;
    }
    
    friend auto operator*(const Matrix<double> &a, const DiagonalMatrix &b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        Matrix<double> out(a.nrows(), b.ncols(), false);
        for (std::size_t i = 0; i < out.nrows(); ++i) {
            for (std::size_t j = 0; j < b.size(); ++j)
                out(i, j) = a(i, j) * b[j];
            for (std::size_t j = b.size(); j < out.ncols(); ++j)
                out(i, j) = 0;
        }
        return out;
    }
    
    friend auto operator*(const DiagonalMatrix &a, const DiagonalMatrix &b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        auto out = DiagonalMatrix(a.nrows(), b.ncols());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = a[i] * b[i];
        return out;
    }
    
    friend auto tr(const DiagonalMatrix &a) { return a; }
    
    friend auto diag(const DiagonalMatrix &a) {
        Matrix out(a.size(), 1, false);
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = a[i];
        return out;
    }
    
    friend Maybe_error<double> logdet(const DiagonalMatrix &x) {
        double out = 0;
        for (std::size_t i = 0; i < x.nrows(); ++i) {
            if (x(i, i) > 0)
                out += std::log(x(i, i));
            else
                return error_message("logdet of negative value at (" + std::to_string(i) +
                                   "," + std::to_string(i) + ")");
        }
        return out;
    }
    
    template <class F> friend auto apply(F &&f, DiagonalMatrix &&x) {
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = f(x[i]);
        return x;
    }
    
    friend auto xtAx(const Matrix<T> &x, const DiagonalMatrix &A) {
        assert(x.nrows() == A.nrows());
        auto out = T{};
        for (std::size_t i = 0; i < A.nrows(); ++i)
            out += x[i] * A(i, i) * x[i];
        return out;
    }
    
    friend auto xAxt(const Matrix<T> &x, const DiagonalMatrix &A) {
        assert(x.ncols() == A.nrows());
        auto out = T{};
        for (std::size_t i = 0; i < A.nrows(); ++i)
            out += x[i] * A(i, i) * x[i];
        return out;
    }
    
    template <class F>
        requires(!is_Maybe_error<std::invoke_result_t<F, T>>)
    friend auto apply(F &&f, DiagonalMatrix const &x) {
        auto out = DiagonalMatrix(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i)
            out[i] = f(x[i]);
        return out;
    }
    
    template <class F>
        requires(is_Maybe_error<std::invoke_result_t<F, T>>)
    friend Maybe_error<DiagonalMatrix<T>> apply(F &&f, DiagonalMatrix const &x) {
        auto out = DiagonalMatrix(x.nrows(), x.ncols(), false);
        for (std::size_t i = 0; i < x.size(); ++i) {
            auto v = f(x[i]);
            if (v)
                out[i] = std::move(v.value());
            else
                return v.error() + " at the " + std::to_string(i) + "th row";
        }
        return out;
    }
    
    template <class F> friend auto reduce(F &&f, const DiagonalMatrix &x) {
        auto cum = x[0];
        for (std::size_t i = 1; i < x.size(); ++i)
            cum = f(cum, x[i]);
        return cum;
    }
    
    friend auto &operator<<(std::ostream &os, const DiagonalMatrix &x) {
        os << "Diagonal matrix "
           << "nrows: " << x.nrows() << " ncols: " << x.ncols() << "\n";
        for (std::size_t i = 0; i < x.size(); ++i)
            os << std::setw(cell_width())<<x[i] << " ";
        os << "\n";
        
        return os;
    }
};
template <class T> DiagPosDetMatrix<T> XTX(const DiagonalMatrix<T> &a);

template <class T> class DiagPosDetMatrix : public DiagonalMatrix<T> {
private:
    DiagPosDetMatrix(DiagonalMatrix<T> &&x) : base_type{x} {}
    
public:
    using base_type = DiagonalMatrix<T>;
    explicit DiagPosDetMatrix(std::size_t _nrows)
        : base_type{_nrows, _nrows, false} {}
    explicit DiagPosDetMatrix(std::size_t _nrows, T value)
        : base_type{_nrows, _nrows, value} {}
    explicit DiagPosDetMatrix(std::initializer_list<T> &&a)
        : base_type(std::move(a)) {}
    
    explicit DiagPosDetMatrix(const Matrix<T> &a) : base_type(a) {}
    
    explicit DiagPosDetMatrix(Matrix<T> &&a) : base_type(std::move(a)) {}
    
    using base_type::operator[];
    using base_type::operator();
    
    ~DiagPosDetMatrix() {}
    
    template <class T0>
    friend DiagPosDetMatrix<T0> XXT(const DiagonalMatrix<T0> &a);
    
    friend DiagPosDetMatrix<T> XTX(const DiagonalMatrix<T> &a) {
        return DiagPosDetMatrix<T>(a * a);
    }
    
    friend auto operator*(const DiagPosDetMatrix &a, const DiagPosDetMatrix &b) {
        assert(a.ncols() == b.nrows() && "matrix product dimensions mismatch");
        auto out = DiagPosDetMatrix(a.nrows(), b.ncols());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = a[i] * b[i];
        return out;
    }
    friend auto operator*(const DiagPosDetMatrix &a, double b) {
        auto out = DiagPosDetMatrix(a.nrows(), a.ncols());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = a[i] * b;
        return out;
    }
    
    friend auto &tr(const DiagPosDetMatrix &a) { return a; }
    
    friend Maybe_error<DiagPosDetMatrix> cholesky(const DiagPosDetMatrix &a) {
        DiagPosDetMatrix out(a.nrows());
        for (std::size_t i = 0; i < out.nrows(); ++i)
            if (a(i, i) > 0)
                out[i] = std::sqrt(a(i, i));
            else if (a(i, i) == 0)
                return std::string("zero value at") + std::to_string(i);
            else
                return "sqrt of negative value " + std::to_string(a(i, i));
        return out;
    }
};

template <class T> auto IdM(std::size_t ndim) {
    return DiagPosDetMatrix<T>(ndim, T{1.0});
}

template <class T> auto inv_from_chol(const DownTrianMatrix<T> &x) {
    return lapack::Lapack_LT_Cholesky_inv(x);
}

template <class T> auto inv_from_chol(const UpTrianMatrix<T> &x) {
    return lapack::Lapack_UT_Cholesky_inv(x);
}

template <class T>
Maybe_error<DiagPosDetMatrix<T>> inv_from_chol(const DiagPosDetMatrix<T> &x) {
    DiagPosDetMatrix<T> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        if (x(i, i) > 0)
            out[i] = std::pow(x(i, i), -2);
        else
            return std::string("error inverting chol");
    return out;
}

template <class T>
Maybe_error<DiagPosDetMatrix<T>> inv(const DiagPosDetMatrix<T> &x) {
    DiagPosDetMatrix<T> out(x.nrows());
    for (std::size_t i = 0; i < x.nrows(); ++i)
        if (out(i, i) > 0)
            out[i] = std::pow(x(i, i), -1);
        else
            return std::string("error inverting chol");
    return out;
}

template <class T> auto diag(const Matrix<T> &a) {
    return DiagonalMatrix<T>(a);
}
template <class T> auto diag(Matrix<T> &&a) {
    return DiagonalMatrix<T>(std::move(a));
}

template <class T> auto diag(std::initializer_list<T> &&a) {
    return DiagonalMatrix<T>(std::move(a));
}




template <class T> auto diagpos(std::initializer_list<T> &&a) {
    return DiagPosDetMatrix<T>(std::move(a));
}

template <class T> auto qr(const Matrix<T> &a) { return lapack::Lapack_QR(a); }

auto eigs(const Matrix<double> &x, bool does_permutations=true,
          bool does_diagonal_scaling=true,
          bool computes_eigenvalues_condition_numbers=false,
          bool computes_eigenvectors_condition_numbers=false) {
    return lapack::Lapack_EigenSystem(x, does_permutations,
                                      does_diagonal_scaling,
                                      computes_eigenvalues_condition_numbers,
                                      computes_eigenvectors_condition_numbers); }



auto XXT(const Matrix<double> &a) {
    return lapack::Lapack_Product_Self_Transpose(a, false);
}

template <class T> auto XXT(const DiagonalMatrix<T> &a) {
    return DiagPosDetMatrix<T>(a * a);
}

auto XTX(const Matrix<double> &a) {
    return lapack::Lapack_Product_Self_Transpose(a, true);
}

template <class T>
SymPosDefMatrix<T> AT_D_A(const Matrix<T> &A, const DiagPosDetMatrix<T> &D) {
    SymPosDefMatrix<T> out(A.nrows(), T{});
    for (std::size_t i = 0; i < out.nrows(); ++i)
        for (std::size_t j = i; j < out.ncols(); ++j)
            for (std::size_t k = 0; k < out.ncols(); ++k)
                out.set(i, j, out(i, j) + A(k, i) * D(k, k) * A(k, j));
    return out;
}



template<class Matrix>
double Trace(Matrix const& x)
{
    double out=0;
    for (std::size_t i=0; i<std::min(x.nrows(),x.ncols()); ++i)
        out+=x(i,i);
    return out;
}
#endif // MATRIX_H
