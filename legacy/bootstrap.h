#pragma once

#include <random_samplers.h>
#include <cstddef>
#include <format>
#include <type_traits>
#include <utility>
#include <vector> // IWYU pragma: keep
#include <optional>
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

template <class WrappedVecSpace, class... Vs>
auto apply_Probit_statistics_impl(const std::vector<WrappedVecSpace>& samples,
                                  const std::set<double>& cis,
                                  std::type_identity<var::Vector_Space<Vs...>>) {
    return var::Vector_Space<Probit_statistics<Vs>...>(
        Probit_statistics<Vs>(
            samples,
            [](const WrappedVecSpace& v) -> decltype(auto) { return get<Vs>(v())(); },
            cis)...);
}

template <class WrappedVecSpace>
    requires(is_of_this_template_type_v<value_type_t<WrappedVecSpace>, var::Vector_Space> &&
             !is_of_this_template_type_v<WrappedVecSpace, var::Vector_Space>)
auto apply_Probit_statistics(const std::vector<WrappedVecSpace>& samples, const std::set<double>& cis) {
    using raw_space_type = value_type_t<WrappedVecSpace>;
    return apply_Probit_statistics_impl(samples, cis, std::type_identity<raw_space_type>{});
}
   




template<class vectorSpace, class F, class ...T>
requires requires(F&& f, const std::vector<vectorSpace> x, const std::vector<std::size_t> indices, T&&...args ) 
{ { f(x, indices, std::forward<T>(args)...) }; }
auto bootstrap_it(F&& f,const std::vector<vectorSpace>& vs,  std::size_t n_replicates, mt_64i& gen, T&& ...args) {
    using R=std::decay_t<decltype(f(vs,generate_bootstrap_indices(vs.size(), gen), std::forward<T>(args)...))>;
    // 1) Serial pass: pre-generate every resample's indices, consuming `gen` in the
    //    same order as the old interleaved loop → bit-identical results (f never
    //    touches gen). This strips the only shared-mutable so f() can run parallel.
    std::vector<std::vector<std::size_t>> idxs(n_replicates);
    for (std::size_t i = 0; i < n_replicates; ++i)
        idxs[i] = generate_bootstrap_indices(vs.size(), gen);
    // 2) Parallel pass: replicates are independent; f must be pure/thread-safe.
    //    f/args are invoked n_replicates times, so pass as lvalues (no forward).
    std::vector<std::optional<R>> slots(n_replicates);
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n_replicates; ++i)
        slots[i].emplace(f(vs, idxs[i], args...));
    bootstrap<R> out;
    out.reserve(n_replicates);
    for (auto& s : slots) out.push_back(std::move(*s));
    return out;
 }

 
template<class vectorSpace, class F, class ...T>
requires requires(F&& f, const std::vector<vectorSpace> x, const std::vector<std::size_t> indices, T&&...args )
{ { f(x, indices, std::forward<T>(args)...) }; }
auto bootstrap_it_to_Probit(F&& f,const std::vector<vectorSpace>& vs,  std::size_t n_replicates, std::set<double>const & cis, mt_64i& gen, T&& ...args) {
    auto samples=bootstrap_it(f, vs, n_replicates, gen, std::forward<T>(args)...);
    return apply_Probit_statistics(samples, cis);
 }


// Two-input variant: draws independent bootstrap index sets for two parallel
// vectors per replicate. Used by the IDM diagnostic where dy (per-sample
// derivative evolution, size N) and F (per-decimated-sample numerical Fisher
// contributions, size N/decimate) are both bootstrap-resampled but with no
// pairing between their indices — IDM is a global property so the bootstrap
// is a comparison-of-estimators tool, not a re-recording-variance estimator
// (see derivation in eLife 2025 supplementary).
template<class V1, class V2, class F, class ...T>
requires requires(F&& f, const std::vector<V1> x1, const std::vector<std::size_t> i1,
                  const std::vector<V2> x2, const std::vector<std::size_t> i2, T&&...args)
{ { f(x1, i1, x2, i2, std::forward<T>(args)...) }; }
auto bootstrap_it_two(F&& f, const std::vector<V1>& vs1, const std::vector<V2>& vs2,
                       std::size_t n_replicates, mt_64i& gen, T&& ...args) {
    using R = std::decay_t<decltype(f(vs1, generate_bootstrap_indices(vs1.size(), gen),
                                       vs2, generate_bootstrap_indices(vs2.size(), gen),
                                       std::forward<T>(args)...))>;
    // 1) Serial RNG pass: pre-draw all replicate index vectors, in the same
    //    (idx1 then idx2, per replicate) order as the original loop, so results
    //    stay seed-reproducible independent of the parallel thread count below.
    std::vector<std::vector<std::size_t>> idx1s(n_replicates), idx2s(n_replicates);
    for (std::size_t i = 0; i < n_replicates; ++i) {
        idx1s[i] = generate_bootstrap_indices(vs1.size(), gen);
        idx2s[i] = generate_bootstrap_indices(vs2.size(), gen);
    }
    // 2) Parallel pass: replicates are independent; f must be pure. f/args are
    //    invoked n_replicates times, so pass as lvalues (no std::forward here).
    std::vector<std::optional<R>> slots(n_replicates);
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n_replicates; ++i)
        slots[i].emplace(f(vs1, idx1s[i], vs2, idx2s[i], args...));
    bootstrap<R> out;
    out.reserve(n_replicates);
    for (auto& s : slots) out.push_back(std::move(*s));
    return out;
}

template<class V1, class V2, class F, class ...T>
requires requires(F&& f, const std::vector<V1> x1, const std::vector<std::size_t> i1,
                  const std::vector<V2> x2, const std::vector<std::size_t> i2, T&&...args)
{ { f(x1, i1, x2, i2, std::forward<T>(args)...) }; }
auto bootstrap_it_two_to_Probit(F&& f, const std::vector<V1>& vs1, const std::vector<V2>& vs2,
                                 std::size_t n_replicates, std::set<double> const& cis,
                                 mt_64i& gen, T&& ...args) {
    auto samples = bootstrap_it_two(std::forward<F>(f), vs1, vs2, n_replicates, gen,
                                     std::forward<T>(args)...);
    return apply_Probit_statistics(samples, cis);
}


// Paired variant of bootstrap_it_two: vs1 and vs2 must have the same length;
// each replicate uses a SINGLE shared index vector applied to both. Removes the
// independent-bootstrap noise inflation that hurts ratio-of-quadratic-forms
// diagnostics (IDM, GFD) when J and F are correlated at the recording level
// under the information identity. Requires upstream F_per_recording to be
// computed without decimation (size == dy.size()).
template<class V1, class V2, class F, class ...T>
requires requires(F&& f, const std::vector<V1> x1, const std::vector<std::size_t> i1,
                  const std::vector<V2> x2, const std::vector<std::size_t> i2, T&&...args)
{ { f(x1, i1, x2, i2, std::forward<T>(args)...) }; }
auto bootstrap_it_two_paired(F&& f, const std::vector<V1>& vs1, const std::vector<V2>& vs2,
                              std::size_t n_replicates, mt_64i& gen, T&& ...args) {
    assert(vs1.size() == vs2.size() &&
           "bootstrap_it_two_paired requires matching sizes (use decimate=1)");
    // 1) Serial RNG pass: pre-draw one shared index vector per replicate, in the
    //    same order/count as the original (replicate 0 first), so results stay
    //    seed-reproducible regardless of the parallel thread count below.
    std::vector<std::vector<std::size_t>> idxs(n_replicates);
    for (std::size_t i = 0; i < n_replicates; ++i)
        idxs[i] = generate_bootstrap_indices(vs1.size(), gen);
    using R = std::decay_t<decltype(f(vs1, idxs[0], vs2, idxs[0],
                                       std::forward<T>(args)...))>;
    // 2) Parallel pass: replicates independent; f must be pure. f/args invoked
    //    n_replicates times → pass as lvalues (no std::forward here).
    std::vector<std::optional<R>> slots(n_replicates);
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n_replicates; ++i)
        slots[i].emplace(f(vs1, idxs[i], vs2, idxs[i], args...));
    bootstrap<R> out;
    out.reserve(n_replicates);
    for (auto& s : slots) out.push_back(std::move(*s));
    return out;
}

template<class V1, class V2, class F, class ...T>
requires requires(F&& f, const std::vector<V1> x1, const std::vector<std::size_t> i1,
                  const std::vector<V2> x2, const std::vector<std::size_t> i2, T&&...args)
{ { f(x1, i1, x2, i2, std::forward<T>(args)...) }; }
auto bootstrap_it_two_paired_to_Probit(F&& f, const std::vector<V1>& vs1,
                                        const std::vector<V2>& vs2,
                                        std::size_t n_replicates,
                                        std::set<double> const& cis, mt_64i& gen,
                                        T&& ...args) {
    auto samples = bootstrap_it_two_paired(std::forward<F>(f), vs1, vs2, n_replicates, gen,
                                            std::forward<T>(args)...);
    return apply_Probit_statistics(samples, cis);
}






