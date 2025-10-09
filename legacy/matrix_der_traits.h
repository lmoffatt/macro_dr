#pragma once

template <typename>
class SymmetricMatrix;
template <typename>
class DiagPosDetMatrix;

template <typename>
class SymPosDefMatrix;

template <typename>
class DiagonalMatrix;

template <class M>
struct M_der {
    using type = M;
};

template <>
struct M_der<SymPosDefMatrix<double>> {
    using type = SymmetricMatrix<double>;
};

template <>
struct M_der<DiagPosDetMatrix<double>> {
    using type = DiagonalMatrix<double>;
};

template <class M>
using M_der_t = typename M_der<M>::type;
