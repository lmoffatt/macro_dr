#ifndef FUNCTION_MEMOIZATION_H
#define FUNCTION_MEMOIZATION_H
#include "distributions.h"
#include "function_measure_verification_and_optimization.h"
#include <concepts>
#include <cstddef>
#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace var {

template <class Y, class... Xs> class Memoiza_all_values {

  std::map<std::tuple<Xs...>, Y> m_values;

  template <class... Ts>
      requires(sizeof...(Ts)==sizeof...(Xs))
  static std::tuple<Xs...> match(Ts &&...ts) {
    return std::tuple<Xs...>(std::forward<Ts>(ts)...);
  }

  template <class T, class... Ts>
      requires(sizeof...(Ts)+1> sizeof...(Xs))
  static std::tuple<Xs...> match(T &&, Ts &&...ts) {
    return match(std::forward<Ts>(ts)...);
  }

public:
  constexpr Memoiza_all_values() {}
  void clear() {
      m_values.clear();
  }
  
  
  template <class F, class... Ts>
    requires(std::is_convertible_v<std::invoke_result_t<F, Ts...>, Y>)
  auto &get_or_calc(F &&f, Ts &&...xs) {
      auto x = match(std::forward<Ts>(xs)...);
    if (auto it = m_values.find(x); it != m_values.end())
      return it->second;
    else {
        m_values[x] = std::invoke(std::forward<F>(f), std::forward<Ts>(xs)...);
      return m_values[x];
    }
  }
};

template <class X, class Y, class... Ks> class Memoiza_closed_values {

  std::map<X, std::size_t> m_index;
  std::vector<Y> m_values;
  Y m_last_value;
  std::tuple<Ks...> m_other;

  static std::multimap<std::size_t, X>
  calc_reverse_count(std::map<X, std::size_t> x) {
    std::multimap<std::size_t, X> out;
    for (auto &e : x) {
      out.insert(std::pair(e.second, e.first));
    }
    return out;
  }

  static std::map<X, std::size_t> calc_count(const std::vector<X> x) {
    std::map<X, std::size_t> out;
    for (auto &e : x)
      ++out[e];
    return out;
  }

  static auto calc_min_count(std::map<X, std::size_t> map,
                             std::size_t max_size_buffer) {
    auto f_rev = calc_reverse_count(map);
    auto it = std::advance(f_rev.begin(), max_size_buffer + 1);
    return it->first - 1ul;
  }

  static std::map<X, std::size_t> build_index(const std::vector<X> x,
                                              std::size_t max_size_buffer) {

    auto counts = calc_count(x);
    auto min_count = calc_min_count(counts, max_size_buffer);
    std::map<X, std::size_t> out;
    auto i = 0;
    for (auto &e : counts)
      if (e.second >= min_count) {
        out.emplace(e.first, i);
        ++i;
      }
    return out;
  }

public:
  constexpr Memoiza_closed_values(std::vector<X> x_values,
                                  std::size_t size_buffer)
      : m_index{build_index(x_values, size_buffer)}, m_values{} {
    m_values.resize(m_index.size());
  }

  constexpr Memoiza_closed_values() {}

  template <class F>
    requires(std::is_convertible_v<std::invoke_result_t<F, X>, Y>)
  void calc_all(F &&f, Ks &&...ks) {
    m_other = std::make_tuple(ks...);
    for (auto &e : m_index)
      m_values[e.second] = std::apply(
          [&f, &e](auto &...ks) {
            return std::invoke(std::forward<F>(f), e.first, ks...);
          },
          m_other);
  }

  template <class F>
    requires(std::is_convertible_v<std::invoke_result_t<F, X, Ks...>, Y>)
  auto &get_this(F &&f, const X &x, Ks &&...ks) {
    assert(std::tuple(ks...) == m_other);
    if (auto it = m_index.find(x); it != m_index.end())
      return m_values[it->second];
    else {
      m_last_value = std::apply(
          [&f, &x](auto &...ks) {
            return std::invoke(std::forward<F>(f), x, ks...);
          },
          m_other);
      return m_last_value;
    }
  }
};

template <class, class, class> class Thread_Memoizer;

template <class Id, class... Fun, class Y, class... Xs,
          template <class...> class F,
          template <class, class...> class Memoizer>
//  requires(std::is_convertible_v<std::invoke_result_t<F<Id, Fun...>, Xs...>, Y>)
class Thread_Memoizer<Id, F<Id, Fun...>, Memoizer<Y, Xs...>> {
  F<Id, Fun...> m_f;
  std::vector<Memoizer<Y, Xs...>> m_memoiza;
  std::size_t m_n_threads;

public:
  static constexpr bool is_threadable= true; 
  
  auto &get_Fun() { return m_f.get_Fun(); }

  constexpr Thread_Memoizer(F<Id, Fun...> &&t_f, Memoizer<Y, Xs...>,
                            std::size_t n_threads)
      : m_f{std::move(t_f)}, m_memoiza{n_threads}, m_n_threads{n_threads} {}

  void clear(I_thread i) {
      m_memoiza[i.i % m_n_threads].clear();
  }
  constexpr Thread_Memoizer() {}
  auto &operator[](Id) { return *this; }
  auto &operator[](Id) const { return *this; }

  auto n_threads() const { return m_n_threads; }

  template <class... Ts>
  auto &operator()(I_thread i_thread, Ts&&... ts) {
      if constexpr( m_f.is_threadable)
          return m_memoiza[i_thread.i%m_n_threads].get_or_calc(m_f, i_thread,std::forward<Ts>(ts)...);
      else
          return m_memoiza[i_thread.i%m_n_threads].get_or_calc(m_f, std::forward<Ts>(ts)...);
  }
  
  template <class... Ts> friend auto apply_F(Thread_Memoizer& me, I_thread i_thread,Ts &&...ts) {
      
      return me(i_thread,std::forward<Ts>(ts)...);
  }
  
};

template <class Id, class... Fun, class Y, class... Xs,
          template <class...> class F,
          template <class, class...> class Memoizer>
Thread_Memoizer(F<Id, Fun...>, Memoizer<Y, Xs...>, std::size_t)
    -> Thread_Memoizer<Id, F<Id, Fun...>, Memoizer<Y, Xs...>>;

template <class Id, class... Fun, class Y, class... Xs,
         template <class...> class F,
         template <class, class...> class Memoizer>
Thread_Memoizer(F<Id, Fun...>, Memoizer<Y, Xs...>)
    -> Thread_Memoizer<Id, F<Id, Fun...>, Memoizer<Y, Xs...>>;


} // namespace var

#endif // FUNCTION_MEMOIZATION_H
