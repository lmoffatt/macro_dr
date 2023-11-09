#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H
#include <random>
#include <vector>
#include <algorithm>
#include <cstdlib>


constexpr std::size_t large_prime=654612646693;

template<class bidiiter>
bidiiter randomly_extract_n(std::mt19937_64& mt,bidiiter begin, bidiiter end, std::size_t num_random ) {
    std::size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, std::uniform_int_distribution<std::size_t>(0, left-1)(mt));
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}




#endif // RANDOM_SAMPLERS_H
