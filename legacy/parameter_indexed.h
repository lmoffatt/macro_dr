#ifndef PARAMETER_INDEXED_H
#define PARAMETER_INDEXED_H

#include <concepts>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "general_algorithm_on_containers.h"
#include "matrix.h"
#include "parameters.h"

namespace var {

template <class Params>
concept ParameterMetadataLike = requires(const Params& params) {
    { params.parameters().names() };
};

template <class T>
struct is_parameter_indexed : std::false_type {};

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
class ParameterIndexed {
    ValueT m_value;
    Params const* m_params = nullptr;

   public:
    using value_type = ValueT;
    using params_type = Params;

    constexpr ParameterIndexed() = default;
    constexpr ParameterIndexed(const ValueT& value) : m_value(value) {}
    constexpr ParameterIndexed(ValueT&& value) : m_value(std::move(value)) {}

    constexpr ParameterIndexed(const ValueT& value, const Params& params)
        : m_value(value), m_params(&params) {}
    constexpr ParameterIndexed(ValueT&& value, const Params& params)
        : m_value(std::move(value)), m_params(&params) {}

    constexpr ParameterIndexed(const ValueT& value, const Params* params)
        : m_value(value), m_params(params) {}
    constexpr ParameterIndexed(ValueT&& value, const Params* params)
        : m_value(std::move(value)), m_params(params) {}

    constexpr ValueT& value() { return m_value; }
    constexpr const ValueT& value() const { return m_value; }
    constexpr ValueT& matrix() { return m_value; }
    constexpr const ValueT& matrix() const { return m_value; }

    constexpr bool has_parameters() const { return m_params != nullptr; }
    constexpr const Params& parameters() const { return *m_params; }
    constexpr const Params* parameters_ptr() const { return m_params; }
    constexpr void set_parameters(const Params& params) { m_params = &params; }
    constexpr void set_parameters(const Params* params) { m_params = params; }

    constexpr auto size() const requires requires(const ValueT& x) { x.size(); } {
        return m_value.size();
    }

    constexpr auto nrows() const requires requires(const ValueT& x) { x.nrows(); } {
        return m_value.nrows();
    }

    constexpr auto ncols() const requires requires(const ValueT& x) { x.ncols(); } {
        return m_value.ncols();
    }

    constexpr decltype(auto) operator[](std::size_t i) requires requires(ValueT& x) { x[i]; } {
        return m_value[i];
    }

    constexpr decltype(auto) operator[](std::size_t i) const
        requires requires(const ValueT& x) { x[i]; }
    {
        return m_value[i];
    }

    constexpr decltype(auto) operator()(std::size_t r, std::size_t c)
        requires requires(ValueT& x) { x(r, c); }
    {
        return m_value(r, c);
    }

    constexpr decltype(auto) operator()(std::size_t r, std::size_t c) const
        requires requires(const ValueT& x) { x(r, c); }
    {
        return m_value(r, c);
    }

    friend double fullsum(const ParameterIndexed& x) { return fullsum(x.value()); }

    friend std::ostream& operator<<(std::ostream& os, const ParameterIndexed& x) {
        return os << x.value();
    }

    friend bool operator==(const ParameterIndexed& lhs, const ParameterIndexed& rhs) {
        return lhs.value() == rhs.value();
    }

    friend bool operator<(const ParameterIndexed& lhs, const ParameterIndexed& rhs) {
        return lhs.value() < rhs.value();
    }
};

template <class ValueT, class Params>
struct is_parameter_indexed<ParameterIndexed<ValueT, Params>> : std::true_type {};

template <class T>
inline constexpr bool is_parameter_indexed_v = is_parameter_indexed<std::remove_cvref_t<T>>::value;

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
const Params* first_parameter_metadata(const ParameterIndexed<ValueT, Params>& x) {
    return x.parameters_ptr();
}

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
const Params* first_parameter_metadata(const ParameterIndexed<ValueT, Params>& lhs,
                                       const ParameterIndexed<ValueT, Params>& rhs) {
    return lhs.has_parameters() ? lhs.parameters_ptr() : rhs.parameters_ptr();
}

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
auto operator+(const ParameterIndexed<ValueT, Params>& lhs,
               const ParameterIndexed<ValueT, Params>& rhs) {
    return ParameterIndexed<ValueT, Params>(lhs.value() + rhs.value(),
                                            first_parameter_metadata(lhs, rhs));
}

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
auto operator-(const ParameterIndexed<ValueT, Params>& lhs,
               const ParameterIndexed<ValueT, Params>& rhs) {
    return ParameterIndexed<ValueT, Params>(lhs.value() - rhs.value(),
                                            first_parameter_metadata(lhs, rhs));
}

template <class Scalar, class ValueT, class Params>
    requires ParameterMetadataLike<Params> && std::convertible_to<Scalar, double>
auto operator*(const ParameterIndexed<ValueT, Params>& x, Scalar a) {
    return ParameterIndexed<ValueT, Params>(x.value() * static_cast<double>(a), x.parameters_ptr());
}

template <class Scalar, class ValueT, class Params>
    requires ParameterMetadataLike<Params> && std::convertible_to<Scalar, double>
auto operator*(Scalar a, const ParameterIndexed<ValueT, Params>& x) {
    return ParameterIndexed<ValueT, Params>(static_cast<double>(a) * x.value(), x.parameters_ptr());
}

template <class Scalar, class ValueT, class Params>
    requires ParameterMetadataLike<Params> && std::convertible_to<Scalar, double>
auto operator/(const ParameterIndexed<ValueT, Params>& x, Scalar a) {
    return ParameterIndexed<ValueT, Params>(x.value() / static_cast<double>(a), x.parameters_ptr());
}

template <bool include_covariance, class ValueT, class Params>
    requires ParameterMetadataLike<Params>
auto sqr_X(const ParameterIndexed<ValueT, Params>& x) {
    return ParameterIndexed<std::decay_t<decltype(::sqr_X<include_covariance>(x.value()))>, Params>(
        ::sqr_X<include_covariance>(x.value()), x.parameters_ptr());
}

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
Maybe_error<bool> compare_contents(const ParameterIndexed<ValueT, Params>& s0,
                                   const ParameterIndexed<ValueT, Params>& s1, double rel_error,
                                   double abs_error, std::size_t max_errors) {
    return compare_contents(s0.value(), s1.value(), rel_error, abs_error, max_errors);
}

template <class ValueT, class Params>
    requires ParameterMetadataLike<Params>
std::vector<std::string> parameter_names(const ParameterIndexed<ValueT, Params>& x) {
    if (!x.has_parameters()) {
        return {};
    }
    const auto& names = x.parameters().parameters().names();
    return std::vector<std::string>(names.begin(), names.end());
}

}  // namespace var

#endif  // PARAMETER_INDEXED_H
