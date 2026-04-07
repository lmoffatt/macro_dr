#pragma once

#include <random_samplers.h>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector> // IWYU pragma: keep
#include <random>
#include <algorithm>
#include <cmath>
#include <moment_statistics.h>
#include <variables.h>

template<class T>
using bootstrap=std::vector<T>;


inline std::vector<std::vector<std::size_t>> generate_bootstrap_indices(std::size_t n, std::size_t num_replicates, mt_64i& gen) {
    // 1. Initialize the Random Number Generator
   
    // 2. Define the distribution [0, n-1]
    // This is the "Standard" way: every index has an equal 1/n chance 
    // of being picked for every slot in the resample.
    std::uniform_int_distribution<std::size_t> dist(0, n - 1);

    std::vector<std::vector<std::size_t>> all_indices;
    all_indices.reserve(num_replicates);

    for (std::size_t b = 0; b < num_replicates; ++b) {
        std::vector<std::size_t> replicate;
        replicate.reserve(n);
        
        for (std::size_t i = 0; i < n; ++i) {
            replicate.push_back(dist(gen));
        }
        
        all_indices.push_back(std::move(replicate));
    }

    return all_indices;
}

inline std::vector<std::size_t> generate_bootstrap_indices(std::size_t n, mt_64i& gen) {
    std::vector<std::size_t> all_indices; 
    if (n>0) {
    std::uniform_int_distribution<std::size_t> dist(0, n - 1);

    all_indices.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            all_indices.push_back(dist(gen));
        }
    }
    return all_indices;
}


template<class...Vs>
auto apply_Probit_statistics(const std::vector<var::Vector_Space<Vs...>>& samples, const std::set<double>& cis) {
     return var::Vector_Space<Probit_statistics<Vs>...>(
        Probit_statistics<Vs>(
            samples,
            [](const auto& v) -> decltype(auto) { return get<Vs>(v)(); },
            cis)...
    );  
}
   




template<class vectorSpace, class F>
requires requires(F&& f, const std::vector<vectorSpace> x, const std::vector<std::size_t> indices) 
{ { f(x, indices) }; }
auto bootstrap_it(F&& f,const std::vector<vectorSpace>& vs,  std::size_t n_replicates, mt_64i& gen) {
    using R=std::decay_t<decltype(f(vs,generate_bootstrap_indices(vs.size(), gen)))>;
    bootstrap<R> out;
    out.reserve(n_replicates);
    for (std::size_t i = 0; i < n_replicates; ++i) {
        out.push_back(std::forward<F>(f)(vs, generate_bootstrap_indices(vs.size(), gen)));
    }
    return out;
 }

 
template<class vectorSpace, class F>
requires requires(F&& f, const std::vector<vectorSpace> x, const std::vector<std::size_t> indices) 
{ { f(x, indices) }; }
auto bootstrap_it_to_Probit(F&& f,const std::vector<vectorSpace>& vs,  std::size_t n_replicates, std::set<double>const & cis, mt_64i& gen) {
    auto samples=bootstrap_it(f, vs, n_replicates, gen);
    return apply_Probit_statistics(samples, cis);
 }








