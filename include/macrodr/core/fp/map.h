#pragma once

#include <vector>
namespace macrodr::core::fp {
template <typename Container, typename Func>
Container map_to_self(const Container& xs, Func f) {
    Container out;
    for (const auto& x : xs) out.insert(out.end(), f(x));
    return out;
}

template <typename Container, typename Func>
auto map_to_vector(const Container& xs, Func f) {
    using Out = std::decay_t<decltype(f(*xs.begin()))>;
    std::vector<Out> out;
    out.reserve(xs.size());
    for (const auto& x : xs) out.push_back(f(x));
    return out;
}

}  // namespace macrodr::core::fp
