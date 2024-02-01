#ifndef FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
#define FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H

#include "maybe_error.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <functional>
#include <ostream>
#include <ratio>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
namespace var {

template <class Id, class Fun>
  requires(std::is_trivial_v<Id>)
class F {
  Fun m_f;

public:
  using myId=Id;
  static constexpr bool is_threadable = false;

  auto &get_Fun() { return m_f; }
  constexpr F(Id, Fun &&t_f) : m_f{std::move(t_f)} {}
  constexpr F(Fun &&t_f) : m_f{std::move(t_f)} {}
  template <class... Ts> constexpr auto operator()(Ts &&...ts) const {
      return m_f(std::forward<Ts>(ts)...);
      
  }

  template <class... Ts> friend auto apply_F(F const me, Ts &&...ts) {

    return me.m_f(std::forward<Ts>(ts)...);
  }
  constexpr F() {}
  constexpr auto &operator[](Id) { return *this; }
  constexpr auto &operator[](Id) const { return *this; }
  
  constexpr auto &operator+=(F) const { return *this; }
  
};

template <class, class> class Time_it;
template <class, class> class Time_it_st;

template <class, class...> class Test_it;

struct I_thread {
  std::size_t i;
};

template <typename T, typename I, typename = void>
struct has_clear_method : public std::false_type {};

template <typename T, typename I>
struct has_clear_method<
    T, I, std::void_t<decltype(std::declval<T>().clear(std::declval<I>()))>>
    : public std::true_type {
  using std::true_type::value;
};

template <typename T, typename I>
constexpr bool has_clear_method_v = has_clear_method<std::decay_t<T>, I>::value;

template <class T, class I> void clearit(T &&me, I i) {
  if constexpr (has_clear_method_v<T, I>)
    me.clear(i);
}
template <class T> void clearit(T &&me) {
  if constexpr (has_clear_method_v<T, void>)
    me.clear();
}

template <class Id, class... Fun, template <class...> class F>
    requires(std::is_trivial_v<Id>)
class Time_it_st<Id, F<Id, Fun...>> {
    F<Id, Fun...> m_f;
    std::chrono::nanoseconds m_sum;
    std::size_t m_count;
    
public:
    using myId=Id;
    void clear(I_thread i) { clearit(m_f, i); }
    
    auto &get_Fun() { return m_f.get_Fun(); }
    
    constexpr Time_it_st(F<Id, Fun...> &&t_f)
        : m_f{std::move(t_f)}, m_sum{std::chrono::nanoseconds::zero()},
        m_count{ 0ul} {}
    constexpr Time_it_st() = default;
    auto &operator[](Id) { return *this; }
    auto &operator[](Id) const { return *this; }
    
    auto naked_function() const { return m_f; }
    
    
    template <class... Ts> auto operator()(Ts &&...ts) {
        const auto start = std::chrono::high_resolution_clock::now();
             if constexpr (std::is_same_v<void, decltype(
                                                   m_f(std::forward<Ts>(ts)...))>) {
                std::invoke(m_f, std::forward<Ts>(ts)...);
                const auto end = std::chrono::high_resolution_clock::now();
                auto dur = end - start;
                m_sum += dur;
                ++m_count;
            } else {
                auto out = std::invoke(m_f, std::forward<Ts>(ts)...);
                const auto end = std::chrono::high_resolution_clock::now();
                auto dur = end - start;
                m_sum += dur;
                ++m_count;
                return out;
            }
    }
    template <class... Ts>
    friend auto apply_time(Time_it_st &me,  Ts &&...ts) {
        
        const auto start = std::chrono::high_resolution_clock::now();
    
            auto out = apply_F(me.m_f, std::forward<Ts>(ts)...);
            const auto end = std::chrono::high_resolution_clock::now();
            auto dur = end - start;
            me.m_sum += dur;
            ++me.m_count;
            return out;
        
    }
    
    auto& operator+=(Time_it_st const& other)
    {
        m_sum+=other.m_sum;
        m_count+=other.m_count;
    }
    
    
    
    auto mean_duration() const {
        m_sum/m_count;    
    }
    
    auto total_duration() const {
        return std::chrono::duration<double>(m_sum).count();
    }
    
    auto count() const { return m_count; }
    
    void reset() {
            m_sum= std::chrono::nanoseconds::zero();
            m_count = 0ul;
    }
    
    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    inline  friend std::ostream &report_title(std::ostream &os, const std::string &sep,
                                             Time_it_st const &) {
        os << ToString(Id{}) << "_sum_time" << sep << ToString(Id{}) << "_count"<< sep << ToString(Id{}) << "_average_time";
        return os;
       
    }
    
    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    friend std::ostream &report_point(std::ostream &os, const std::string &sep,
                                      Time_it_st &me) {
        
        os << me.total_duration() << sep << me.count()<< sep << me.mean_duration();
        me.reset();
        return os;
    }
};




template <class Id, class... Fun, template <class...> class F>
  requires(std::is_trivial_v<Id>)
class Time_it<Id, F<Id, Fun...>> {
  F<Id, Fun...> m_f;
  std::vector<std::chrono::nanoseconds> m_sum;
  std::vector<std::size_t> m_count;
  std::size_t m_n_threads{};

public:
  using myId=Id;
  
  static constexpr bool is_threadable = true;

  void clear(I_thread i) { clearit(m_f, i); }

  auto &get_Fun() { return m_f.get_Fun(); }

  constexpr Time_it(F<Id, Fun...> &&t_f, std::size_t n_threads)
      : m_f{std::move(t_f)}, m_sum{n_threads, std::chrono::nanoseconds::zero()},
        m_count{std::vector<std::size_t>(n_threads, 0ul)},
        m_n_threads{n_threads} {}
  constexpr Time_it() = default;
  auto &operator[](Id) { return *this; }
  auto &operator[](Id) const { return *this; }

  auto naked_function() const { return m_f; }

  auto n_threads() const { return m_n_threads; }

  template <class... Ts> auto operator()(I_thread i_thread, Ts &&...ts) {
    const auto start = std::chrono::high_resolution_clock::now();
    if constexpr (m_f.is_threadable) {
      if constexpr (std::is_same_v<void, decltype(std::invoke(
                                             m_f, i_thread,
                                             std::forward<Ts>(ts)...))>) {
        std::invoke(m_f, i_thread, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return void();
      } else {
        auto out = std::invoke(m_f, i_thread, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return out;
      }
    } else {
      if constexpr (std::is_same_v<void, decltype(
                                               m_f(std::forward<Ts>(ts)...))>) {
        std::invoke(m_f, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
      } else {
        auto out = std::invoke(m_f, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return out;
      }
    }
  }
  template <class... Ts>
  friend auto apply_time(Time_it &me, I_thread i_thread, Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    if constexpr (me.m_f.is_threadable) {
      auto out = apply_F(me.m_f, i_thread, std::forward<Ts>(ts)...);
      const auto end = std::chrono::high_resolution_clock::now();
      auto dur = end - start;
      me.m_sum[i_thread.i % me.n_threads()] += dur;
      ++me.m_count[i_thread.i % me.n_threads()];
      return out;
    } else {
      auto out = apply_F(me.m_f, std::forward<Ts>(ts)...);
      const auto end = std::chrono::high_resolution_clock::now();
      auto dur = end - start;
      me.m_sum[i_thread.i % me.n_threads()] += dur;
      ++me.m_count[i_thread.i % me.n_threads()];
      return out;
    }
  }

  auto mean_duration() const {
    std::vector<
        decltype(std::chrono::duration<float, std::chrono::nanoseconds::period>(
            m_sum[0] / m_count[0]))>
        out(m_sum.size());
    for (std::size_t i = 0; i < n_threads(); ++i)
      out[i] = std::chrono::duration<float, std::chrono::nanoseconds::period>(
          m_sum[i] / std::max(m_count[i], 1ul));
    return out;
  }

  auto total_duration() const {
    auto out = m_sum[0];
    for (std::size_t i = 1; i < n_threads(); ++i)
      out += m_sum[i];
    return std::chrono::duration<double>(out).count();
  }
  auto total_count() const {
    auto out = m_count[0];
    for (std::size_t i = 1; i < n_threads(); ++i)
      out += m_count[i];
    return out;
  }

  auto count() const { return m_count; }

  void reset() {
    for (auto &e : m_sum)
      e = std::chrono::nanoseconds::zero();
    for (auto &e : m_count)
      e = 0ul;
  }

  /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
 inline  friend std::ostream &report_title(std::ostream &os, const std::string &sep,
                                    Time_it const &me) {
    os << ToString(Id{}) << "_sum_time" << sep << ToString(Id{}) << "_count";
    for (std::size_t i = 0; i < me.m_count.size(); ++i)
      os << sep << ToString(Id{}) << "_" << i << "_sum_time" << sep
         << ToString(Id{}) << "_" << i << "_count" << sep << ToString(Id{})
         << "_" << i << "_average_time";
    return os;
  }

  /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
  friend std::ostream &report_point(std::ostream &os, const std::string &sep,
                                    Time_it &me) {

    os << me.total_duration() << sep << me.total_count();
    for (std::size_t i = 0; i < me.m_count.size(); ++i)
      os << sep << std::chrono::duration<double>(me.m_sum[i]).count() << sep
         << me.count()[i] << sep << me.mean_duration()[i];
    me.reset();
    return os;
  }
};

template <class Id, class... Fun, template <class...> class F>
Time_it(F<Id, Fun...>, std::size_t) -> Time_it<Id, F<Id, Fun...>>;

template <class Id, class... Fun, template <class...> class F>
Time_it(F<Id, Fun...>) -> Time_it<Id, F<Id, Fun...>>;

template <class Id, class... Fun, template <class...> class F,
          class Preconditions, class Postconditions>
  requires(std::is_trivial_v<Id>)
class Test_it<Id, F<Id, Fun...>, Preconditions, Postconditions> {
  F<Id, Fun...> m_f;
  std::vector<std::chrono::nanoseconds> m_sum;
  std::vector<std::size_t> m_count;
  std::size_t m_n_threads{};

public:
  using myId=Id;
  
  static constexpr bool is_threadable = true;

  void clear(I_thread i) { clearit(m_f, i); }

  auto &get_Fun() { return m_f.get_Fun(); }

  constexpr Test_it(F<Id, Fun...> &&t_f, std::size_t n_threads)
      : m_f{std::move(t_f)}, m_sum{n_threads, std::chrono::nanoseconds::zero()},
        m_count{std::vector<std::size_t>(n_threads, 0ul)},
        m_n_threads{n_threads} {}
  constexpr Test_it() = default;
  auto &operator[](Id) { return *this; }
  auto &operator[](Id) const { return *this; }

  auto naked_function() const { return m_f; }

  auto n_threads() const { return m_n_threads; }

  template <class... Ts> auto operator()(I_thread i_thread, Ts &&...ts) {
    const auto start = std::chrono::high_resolution_clock::now();
    if constexpr (m_f.is_threadable) {
      if constexpr (std::is_same_v<void, decltype(std::invoke(
                                             m_f, i_thread,
                                             std::forward<Ts>(ts)...))>) {
        std::invoke(m_f, i_thread, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return void();
      } else {
        auto out = std::invoke(m_f, i_thread, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return out;
      }
    } else {
      if constexpr (std::is_same_v<void, decltype(std::invoke(
                                             m_f, std::forward<Ts>(ts)...))>) {
        std::invoke(m_f, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
      } else {
        auto out = std::invoke(m_f, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        m_sum[i_thread.i % n_threads()] += dur;
        ++m_count[i_thread.i % n_threads()];
        return out;
      }
    }
  }
  template <class... Ts>
  friend auto apply_time(Test_it &me, I_thread i_thread, Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    if constexpr (me.m_f.is_threadable) {
      auto out = apply_F(me.m_f, i_thread, std::forward<Ts>(ts)...);
      const auto end = std::chrono::high_resolution_clock::now();
      auto dur = end - start;
      me.m_sum[i_thread.i % me.n_threads()] += dur;
      ++me.m_count[i_thread.i % me.n_threads()];
      return out;
    } else {
      auto out = apply_F(me.m_f, std::forward<Ts>(ts)...);
      const auto end = std::chrono::high_resolution_clock::now();
      auto dur = end - start;
      me.m_sum[i_thread.i % me.n_threads()] += dur;
      ++me.m_count[i_thread.i % me.n_threads()];
      return out;
    }
  }

  auto mean_duration() const {
    std::vector<
        decltype(std::chrono::duration<float, std::chrono::nanoseconds::period>(
            m_sum[0] / m_count[0]))>
        out(m_sum.size());
    for (std::size_t i = 0; i < n_threads(); ++i)
      out[i] = std::chrono::duration<float, std::chrono::nanoseconds::period>(
          m_sum[i] / std::max(m_count[i], 1ul));
    return out;
  }

  auto total_duration() const {
    auto out = m_sum[0];
    for (std::size_t i = 1; i < n_threads(); ++i)
      out += m_sum[i];
    return std::chrono::duration<double>(out).count();
  }
  auto total_count() const {
    auto out = m_count[0];
    for (std::size_t i = 1; i < n_threads(); ++i)
      out += m_count[i];
    return out;
  }

  auto count() const { return m_count; }

  void reset() {
    for (auto &e : m_sum)
      e = std::chrono::nanoseconds::zero();
    for (auto &e : m_count)
      e = 0ul;
  }

  /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
  inline friend std::ostream &report_title(std::ostream &os, const std::string &sep,
                                    Test_it const &me) {
    os << ToString(Id{}) << "_sum_time" << sep << ToString(Id{}) << "_count";
    for (std::size_t i = 0; i < me.m_count.size(); ++i)
      os << sep << ToString(Id{}) << "_" << i << "_sum_time" << sep
         << ToString(Id{}) << "_" << i << "_count" << sep << ToString(Id{})
         << "_" << i << "_average_time";
    return os;
  }

  /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
  friend std::ostream &report_point(std::ostream &os, const std::string &sep,
                                    Test_it &me) {

    os << me.total_duration() << sep << me.total_count();
    for (std::size_t i = 0; i < me.m_count.size(); ++i)
      os << sep << std::chrono::duration<double>(me.m_sum[i]).count() << sep
         << me.count()[i] << sep << me.mean_duration()[i];
    me.reset();
    return os;
  }
};

template <class F> class F_on_thread;
template <class... Fs> class FuncMap;

inline std::ostream &report_title(std::ostream &os, const std::string &, ...) {
  return os;
}

inline std::ostream &report_point(std::ostream &os, const std::string &, ...) {

  return os;
}

template <class F> class F_on_thread {
  F &m_f;
  I_thread i_thread;
  F_on_thread(F &f, I_thread i) : m_f{f}, i_thread{i} {}

public:
  explicit F_on_thread(F &f) : m_f{f}, i_thread{0ul} {}
  F_on_thread fork(I_thread i) { return F_on_thread(m_f, i); }
  template <class Id> decltype(auto) operator[](Id) { return m_f[Id{}]; }

  void clear() { m_f.clear(i_thread); }
  template <class Id, class... Ts> auto f(Id, Ts &&...ts) {
    auto &fun = (*this)[Id{}];
    return fun(i_thread, *this, std::forward<Ts>(ts)...);
  }

  template <class Id, class... Ts> auto fff(Id, Ts &&...ts) {
    auto &&fun = (*this)[Id{}];
    return apply_time(fun, i_thread, *this, std::forward<Ts>(ts)...);
  }

  template <class Id, class... Ts> auto fstop(Id, Ts &&...ts) {
    auto &fun = (*this)[Id{}];
    return std::invoke(fun, i_thread, std::forward<Ts>(ts)...);
  }
};

template <class... Fs> class FuncMap : public Fs... {
  std::ofstream m_file;
  std::string sep = ",";

public:
    
  
  using Fs::operator[]...;

  template <class Id> Nothing operator[](Id) const { return Nothing{}; }

  void clear(I_thread i) { (clearit(static_cast<Fs &>(*this), i), ...); }

  auto &file() { return m_file; }
  FuncMap(const std::string path, Fs &&...fs)
      : Fs{std::move(fs)}..., m_file{std::ofstream(path + "_time_it.cvs")} {m_file<<std::setprecision(std::numeric_limits<double>::digits10 + 1);}
  
  FuncMap( std::ofstream&& f, Fs &&...fs)
      : Fs{std::move(fs)}..., m_file{std::move(f)} {}
  
  FuncMap(const std::string path, Fs const&...fs)
      : Fs{fs}..., m_file{std::ofstream(path + "_time_it.cvs")} {m_file<<std::setprecision(std::numeric_limits<double>::digits10 + 1);}
  
  FuncMap( std::ofstream&& f, Fs const&...fs)
      : Fs{fs}..., m_file{std::move(f)} {m_file<<std::setprecision(std::numeric_limits<double>::digits10 + 1);}
  
  
  template <class... Context_data>
  friend auto &report_point(FuncMap<Fs...> &f, const Context_data &...s) {
    ((f.m_file << s << f.sep), ...);
    ((report_point(f.m_file, f.sep, static_cast<Fs &>(f)) << f.sep), ...);
    f.m_file << "\n";
    return f.m_file;
  }
  template <class... Context_string>
  friend auto &report_title(FuncMap<Fs...> &f, const Context_string &...s) {
    ((f.m_file << s << f.sep), ...);
    ((report_title(f.m_file, f.sep, static_cast<Fs &>(f)) << f.sep), ...);
    f.m_file << "\n";
    return f.m_file;
  }

  F_on_thread<FuncMap> fork(I_thread i) { return F_on_thread(*this).fork(i); }
};


template <class... Fs, class F, class G>
auto operator +(std::pair<FuncMap<Fs...>, F> && f, G const & g )
{
    if constexpr(std::is_same_v<typename G::myId, typename F::myId>)
        return std::pair(FuncMap<Fs...,F>(std::move(f.first.file()),std::move(f.first[typename Fs::myId{}])...,f.second),f.second);
    else
        return std::pair(FuncMap<Fs...,G>(std::move(f.first.file()),std::move(f.first[typename Fs::myId{}])...,g),f.second);
}

template <class G,class... Fs,  class F>
auto insert(const std::string & path,FuncMap<G,Fs...> const & fun,  F const & f )
{
   return (std::pair(FuncMap<G>(path,fun[typename G::myId{}]),f)+...+fun[typename Fs::myId{}]).first;
}


} // namespace var
#endif // FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
