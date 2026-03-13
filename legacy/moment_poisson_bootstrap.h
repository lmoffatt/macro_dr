#ifndef MOMENT_POISSON_BOOTSTRAP_H
#define MOMENT_POISSON_BOOTSTRAP_H

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "moment_statistics.h"

template <class Id, bool include_covariance = false>
class Moment_poisson_bootstrap {
  public:
    using moment_statistics_type = Moment_statistics<Id, include_covariance>;

    explicit Moment_poisson_bootstrap(std::size_t n_replicates)
        : Moment_poisson_bootstrap(n_replicates, std::random_device{}()) {}

    Moment_poisson_bootstrap(std::size_t n_replicates, std::uint64_t seed)
        : replicates_(n_replicates),
          rng_(seed),
          poisson_(1.0) {}

    std::size_t size() const { return replicates_.size(); }

    void reset() {
        for (auto& r : replicates_) {
            r.reset();
        }
    }

    std::vector<moment_statistics_type> const& replicates() const { return replicates_; }
    std::vector<moment_statistics_type>& replicates() { return replicates_; }

    
    void add_sample(Id const& x) {
        for (auto& r : replicates_) {
            const auto k = poisson_(rng_);
            for (int i = 0; i < k; ++i) {
                r &= x;
            }
        }
    }

    Moment_poisson_bootstrap& operator&=(const Id& x) {
        add_sample(x);
        return *this;
    }

    template <class Rng>
    void add_sample(Id const& x, Rng& rng) {
        for (auto& r : replicates_) {
            const auto k = poisson_(rng);
            for (int i = 0; i < k; ++i) {
                r &= x;
            }
        }
    }

  private:
    std::vector<moment_statistics_type> replicates_;
    std::mt19937_64 rng_;
    std::poisson_distribution<int> poisson_;
};

#endif  // MOMENT_POISSON_BOOTSTRAP_H
