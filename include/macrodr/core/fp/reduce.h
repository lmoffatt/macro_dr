#pragma once

#include <functional>
namespace macrodr::core::fp {

template <typename Container, typename BinaryOp, class T>
T reduce(const Container& xs, BinaryOp&& f, T init) {
    for (const auto& x : xs) init = std::invoke(std::forward<BinaryOp>(f), init, x);
    return init;
}

template <typename Container, typename Func>
auto reduce(const Container& xs, Func f) {
    using T = typename Container::value_type;
    auto it = xs.begin();
    T acc = *it++;
    for (; it != xs.end(); ++it) acc = std::invoke(f, acc, *it);
    return acc;
}

}  // namespace macrodr::core::fp
