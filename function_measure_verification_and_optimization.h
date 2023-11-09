#ifndef FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
#define FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <functional>
#include <ratio>
#include <type_traits>
#include <utility>
#include <vector>
namespace var {

template <class Id, class Fun>
  requires(std::is_trivial_v<Id>)
class F {
  Fun m_f;

public:
  constexpr F(Id, Fun &&t_f) : m_f{std::move(t_f)} {}
  constexpr F(Fun &&t_f) : m_f{std::move(t_f)} {}
  template <class... Ts> constexpr auto operator()(Ts &&...ts) const {
    return std::invoke(m_f, std::forward<Ts>(ts)...);
  }

  constexpr F() {}
  constexpr auto &f(Id) { return *this; }
  constexpr auto &f(Id) const { return *this; }
};

template <class> class Time_it;

template <class Id, class Fun>
  requires(std::is_trivial_v<Id>)
class Time_it<F<Id, Fun>> {
  F<Id, Fun> m_f;
  std::chrono::nanoseconds m_sum = std::chrono::nanoseconds::zero();
  std::size_t m_count = 0ul;
  std::string m_sep = ",";
  std::ofstream m_file;

public:
  constexpr Time_it(Id, Fun &&t_f, const std::string path)
      : m_f{std::move(t_f)},
        m_file{std::ofstream(path + "_" + ToString(Id{}) + "_time_it.cvs")} {}
  constexpr Time_it() {}
  constexpr auto &f(Id) { return *this; }
  constexpr auto &f(Id) const { return *this; }

  auto naked_function() const { return m_f; }

  template <class... Ts>
    requires(!std::is_void_v<std::invoke_result_t<F<Id, Fun>, Ts...>>)
  constexpr auto operator()(Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    auto out = std::invoke(m_f, std::forward<Ts>(ts)...);
    const auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - start;
    m_sum += dur;
    ++m_count;
    return out;
  }

  template <class... Ts>
    requires(std::is_void_v<std::invoke_result_t<F<Id, Fun>, Ts...>>)
  constexpr void operator()(Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    std::invoke(m_f, std::forward<Ts>(ts)...);
    const auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - start;
    m_sum += dur;
    ++m_count;
  }

  auto mean_duration() const {
    return std::chrono::duration<float, std::chrono::nanoseconds::period>(
          m_sum / std::max(m_count,1ul));
  }
  auto count() const { return m_count; }

  void reset() {
    m_sum = std::chrono::nanoseconds::zero();
    m_count = 0ul;
  }
  template <class... Context_string>
  friend void report_title(Time_it &me, Context_string const &...s) {
    ((me.m_file << s << me.m_sep), ...);
      me.m_file << "sum_time" << me.m_sep << "count"<<me.m_sep<<"average_time"
              << "\n";
  }

  template <class... Context_data>
  friend void report_point(Time_it &me, Context_data const &...d) {
    ((me.m_file << d << me.m_sep), ...);
      me.m_file << std::chrono::duration<double>(me.m_sum).count() << me.m_sep << me.count() <<me.m_sep<< me.mean_duration() << "\n";
    me.m_file.flush();
    me.reset();
  }
};
template <class Id, class Fun>
Time_it(Id, Fun, std::string) -> Time_it<F<Id, Fun>>;

template <class> class Time_it_parallel;

struct I_thread{
    std::size_t i; 
};

template <class Id, class Fun>
  requires(std::is_trivial_v<Id>)
class Time_it_parallel<F<Id, Fun>> {
  F<Id, Fun> m_f;
  std::vector<std::chrono::nanoseconds> m_sum;
  std::vector<std::size_t> m_count;
  std::string m_sep = ",";
  std::ofstream m_file;

public:
  constexpr Time_it_parallel(Id, Fun &&t_f, const std::string path,
                             std::size_t n_threads)
      : m_f{std::move(t_f)},
        m_sum{n_threads, std::chrono::nanoseconds::zero()},
        m_count{std::vector<std::size_t>(n_threads, 0ul)},
        m_file{std::ofstream(path + "_" + ToString(Id{}) + "_time_it.cvs")}
    {}
  constexpr Time_it_parallel() {}
  constexpr auto &f(Id) { return *this; }
  constexpr auto &f(Id) const { return *this; }

  auto naked_function() const { return m_f; }

  template <class... Ts>
    requires(!std::is_void_v<std::invoke_result_t<F<Id, Fun>, Ts...>>)
  constexpr auto operator()(I_thread i_thread, Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    auto out = std::invoke(m_f, std::forward<Ts>(ts)...);
    const auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - start;
    m_sum[i_thread.i] += dur;
    ++m_count[i_thread.i];
    return out;
  }

  template <class... Ts>
    requires(std::is_void_v<std::invoke_result_t<F<Id, Fun>, Ts...>>)
  constexpr void operator()(I_thread i_thread, Ts &&...ts) {

    const auto start = std::chrono::high_resolution_clock::now();
    std::invoke(m_f, std::forward<Ts>(ts)...);
    const auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - start;
    m_sum[i_thread.i] += dur;
    ++m_count[i_thread.i];
  }
  
  template <class T,class... Ts>
      requires(!std::is_same_v<std::decay_t<T>,I_thread>)
  constexpr auto operator()(T&& not_i_thread, Ts &&...ts)
  {
      return (*this)(I_thread(0ul),std::forward<T>(not_i_thread),std::forward<Ts>(ts)...);
  }
  
  auto mean_duration() const {
    std::vector<
        decltype(std::chrono::duration<float, std::chrono::nanoseconds::period>(
            m_sum[0] / m_count[0]))>
        out(m_sum.size());
    for (std::size_t i = 0; i < m_sum.size(); ++i)
      out[i] = std::chrono::duration<float, std::chrono::nanoseconds::period>(
            m_sum[i] / std::max(m_count[i],1ul));
    return out;
  }
  
  auto total_duration()const
  {
      auto out=m_sum[0];
      for (std::size_t i=1; i<m_sum.size(); ++i)
          out+=m_sum[i];
      return std::chrono::duration<double>(out).count(); 
  }
  auto total_count()const
  {
      auto out=m_count[0];
      for (std::size_t i=1; i<m_count.size(); ++i)
          out+=m_count[i];
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
  template <class... Context_string>
  friend void report_title(Time_it_parallel &me, Context_string const &...s) {
    ((me.m_file << s << me.m_sep), ...);
      me.m_file << "sum_time" << me.m_sep<<"count"<<me.m_sep;
      for (std::size_t i = 0; i < me.m_count.size(); ++i)
       me.m_file  << "sum_time_" << i << me.m_sep << "count_" << i<< me.m_sep << "average_time_" << i 
                << me.m_sep;
    me.m_file << "\n";
  }

  /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
  template <class... Context_data>
  friend void report_point(Time_it_parallel &me, Context_data const &...d) {
    
      ((me.m_file << d << me.m_sep), ...);
      me.m_file << me.total_duration()<< me.m_sep << me.total_count()<<me.m_sep;
      for (std::size_t i = 0; i < me.m_count.size(); ++i)
        me.m_file << std::chrono::duration<double>(me.m_sum[i]).count() << me.m_sep << me.count()[i]
                  << me.m_sep<< me.mean_duration()[i] << me.m_sep;
      me.m_file << "\n";
      me.reset();
  }
};
template <class Id, class Fun>
Time_it_parallel(Id, Fun, std::string, std::size_t)
    -> Time_it_parallel<F<Id, Fun>>;

template <class... Fs> class FuncMap : public Fs... {
public:
  using Fs::f...;
    
  FuncMap(Fs &&...fs) : Fs{std::move(fs)}... {}
};

template <class... Fs, class... Context_data>
void report_point(FuncMap<Fs...> &f, const Context_data &...s) {
  (report_point(static_cast<Fs &>(f), s...), ...);
}

template <class... Fs, class... Context_string>
void report_title(FuncMap<Fs...> &f, const Context_string &...s) {
  (report_title(static_cast<Fs &>(f), s...), ...);
}

template<class F>
class  F_on_thread
{
    F& m_f;
    I_thread i_thread;
public:
    F_on_thread(F& f, I_thread i):m_f{f},i_thread{i}{}
    
    template<class Id>
    auto f(Id)
    {
        return F_on_thread<std::decay_t<decltype(m_f.f(Id{}))>>(m_f.f(Id{}),i_thread);
    }
    template <class... Ts>
    constexpr auto operator()(Ts &&...ts) {
        return std::invoke(m_f,i_thread,std::forward<Ts>(ts)...);   
    }
};


} // namespace var
#endif // FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
