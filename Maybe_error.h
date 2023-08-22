#pragma once
#include <string>
#include <variant>

#ifndef MAYBE_ERROR_H
#define MAYBE_ERROR_H

namespace logic {
class error_message {
    std::string m_;

public:
    error_message(std::string error) : m_{error} {}

    auto operator()() const { return m_; }
};

template <class T> class Maybe_error : public std::variant<T, error_message> {
public:
    using std::variant<T, error_message>::variant;
    constexpr explicit operator bool() const noexcept {
        return this->index() == 0;
    }

    constexpr const T *operator->() const noexcept {
        return std::get_if<0>(*this);
    }

    constexpr T *operator->() noexcept { return std::get_if<0>(*this); }

    constexpr const T &operator*() const & noexcept { return std::get<0>(*this); }
    constexpr T &operator*() & noexcept { return std::get<0>(*this); }
    constexpr const T &&operator*() const && noexcept {
        return std::get<0>(std::move(*this));
    }
    constexpr T &&operator*() && noexcept {
        return std::get<0>(std::move(*this));
    }
    Maybe_error(const T& t) : std::variant<T, error_message>(t) {}
    Maybe_error(T&& t) : std::variant<T, error_message>(std::move(t)) {}
    Maybe_error(error_message m) : std::variant<T, error_message>(m) {}
    error_message error() const {
        if (*this)
            return error_message("");
        else
            return std::get<1>(*this);
    }
};

template<class T>
class build;


} // namespace logic



#endif // MAYBE_ERROR_H

