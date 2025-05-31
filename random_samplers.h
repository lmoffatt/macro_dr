#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <random>
#include <vector>
constexpr std::size_t large_prime = 654612646693;

class mt_64i {
    std::mt19937_64 m_s;
    unsigned long long m_i;

   public:
    using result_type = typename std::mt19937_64::result_type;

    explicit mt_64i(result_type value) : m_s{value}, m_i{0} {
    }

    template <class Sseq>
    explicit mt_64i(Sseq& s) : m_s{s}, m_i{0} {
    }

    static constexpr result_type min() {
        return std::mt19937_64::min();
    }
    static constexpr result_type max() {
        return std::mt19937_64::max();
    }

    void seed(typename std::mt19937_64::result_type value) {
        m_i = 0;
        m_s.seed(value);
    }
    template <class Sseq>
    void seed(Sseq& seq) {
        m_i = 0;
        m_s.seed(seq);
    }

    result_type operator()() {
        ++m_i;
        return m_s();
    }
    void discard(unsigned long long z) {
        m_i += z;
        m_s.discard(z);
    }

    friend std::ostream& operator<<(std::ostream& os, const mt_64i& e) {
        os << e.m_i << "\t " << e.m_s;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, mt_64i& e) {
        is >> e.m_i >> e.m_s;
        return is;
    }

    auto pos() const {
        return m_i;
    }
};

class normal_mt_64 : public mt_64i {};

template <class bidiiter>
bidiiter randomly_extract_n(mt_64i& mt, bidiiter begin, bidiiter end, std::size_t num_random) {
    std::size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, std::uniform_int_distribution<std::size_t>(0, left - 1)(mt));
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

#endif  // RANDOM_SAMPLERS_H
