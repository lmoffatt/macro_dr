#ifndef FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
#define FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <functional>
#include <ostream>
#include <ratio>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "maybe_error.h"
#include "type_algebra.h"

// Forward declaration to avoid heavy includes where we only need the Id type
namespace var {

template <std::size_t N>
class Event_Timing {
    using Stype = std::decay_t<decltype(std::chrono::high_resolution_clock::now())>;
    Stype m_start;
    double last_dur;
    std::array<std::string, N> m_labels;
    std::array<double, N> m_dur;
    std::size_t m_i = 0;

   public:
    Event_Timing(const Stype& s) : m_start{s} {}
    void record(const std::string& s) {
        m_labels[m_i] = s;
        m_dur[m_i] =
            std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start)
                .count();
        ++m_i;
    }
    void record(const std::string& s, std::size_t ipos) {
        m_labels[m_i + ipos] = s + std::to_string(ipos);
        m_dur[m_i + ipos] =
            std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start)
                .count();
    }
    void advance(std::size_t ipos) { m_i += ipos; }
    auto& operator[](std::size_t i) const { return m_dur[i]; }

    void reset() { m_i = 0; }

    void report_title(std::ostream& os) {
        os << "Iter" << "," << "time_iter" << "," << "iter_duration";
        for (std::size_t i = 0; i < m_i; ++i) {
            os << "," << m_labels[i];
        }

        last_dur = m_dur[m_i - 1];
        os << "\n";
    }

    void report_iter(std::ostream& os, std::size_t iter) {
        os << iter << "," << m_dur[0] << "," << m_dur[m_i - 1] - last_dur << ","
           << m_dur[0] - last_dur;
        for (std::size_t i = 1; i < m_i; ++i) {
            os << "," << m_dur[i] - m_dur[i - 1];
        }
        last_dur = m_dur[m_i - 1];
        m_i = 0;

        os << "\n";
    }
};

template <class Id, class Fun>
    requires(std::is_trivial_v<Id>)
class F {
    Fun m_f;

   public:
    using myId = Id;
    static constexpr bool is_threadable = false;

    auto& get_Fun() { return m_f; }
    constexpr F(Id, Fun&& t_f) : m_f{std::move(t_f)} {}
    constexpr F(Fun&& t_f) : m_f{std::move(t_f)} {}
    constexpr F(Fun const& t_f) : m_f{t_f} {}
    template <class... Ts>
    requires std::is_invocable_v<Fun, Ts...>
    constexpr auto operator()(Ts&&... ts) const {
        return m_f(std::forward<Ts>(ts)...);
    }

    F create() const { return F(m_f); }

    F bare_function() const { return F(m_f); }

    template <class... Ts>
    friend auto apply_F(F const me, Ts&&... ts) {
        return me.m_f(std::forward<Ts>(ts)...);
    }
    constexpr F() {}
    constexpr auto& operator[](Id) { return *this; }
    constexpr auto& operator[](Id) const { return *this; }

    constexpr auto& operator+=(F) const { return *this; }

    void clear() const {}
};

template <template <auto...> class MacroR, auto... t>
auto make_F(MacroR<t...>) {
    return F(MacroR<t...>{},
             [](auto&&... x) { return MacroR<t...>{}(std::forward<decltype(x)>(x)...); });
}

template <template <class...> class MacroR, class... t>
auto make_F(MacroR<t...>) {
    return F(MacroR<t...>{},
             [](auto&&... x) { return MacroR<t...>{}(std::forward<decltype(x)>(x)...); });
}

template <template <class...> class MacroR, class... T>
auto make_F_P(in_progress::P<T...>) {
    //    using m00=typename Co<MacroR>::en_make_F_Ts;
    //   using m0=typename std::tuple<T...>::en_make_F_Ts;
    //   using m1=typename MacroR<T...>::en_make_F;

    return F(MacroR<T...>{},
             [](auto&&... x) { return MacroR<T...>{}(std::forward<decltype(x)>(x)...); });
}

template <class, class>
class Time_it_st;

struct I_thread {
    std::size_t i;
};

template <typename T, typename I, typename = void>
struct has_clear_method : public std::false_type {};

template <typename T, typename I>
struct has_clear_method<T, I, std::void_t<decltype(std::declval<T>().clear(std::declval<I>()))>>
    : public std::true_type {
    using std::true_type::value;
};

template <typename T, typename I>
constexpr bool has_clear_method_v = has_clear_method<std::decay_t<T>, I>::value;

template <class T, class I>
void clearit(T&& me, I i) {
    if constexpr (has_clear_method_v<T, I>)
        me.clear(i);
}
template <class T>
void clearit(T&& me) {
    me.clear();
}

template <class Id, class... Fun, template <class...> class F>
    requires(std::is_trivial_v<Id>)
class Time_it_st<Id, F<Id, Fun...>> {
    F<Id, Fun...> m_f;
    std::chrono::nanoseconds m_sum;
    std::size_t m_count;

   public:
    using myId = Id;

    auto& get_Fun() { return m_f.get_Fun(); }

    constexpr Time_it_st(F<Id, Fun...>&& t_f)
        : m_f{std::move(t_f)}, m_sum{std::chrono::nanoseconds::zero()}, m_count{0ul} {}
    constexpr Time_it_st() = default;
    auto& operator[](Id) { return *this; }
    auto& operator[](Id) const { return *this; }

    auto naked_function() const { return m_f; }
    auto create() const { return Time_it_st(m_f.create()); }

    F<Id, Fun...> bare_function() const { return m_f; }

    template <class... Ts>
    auto operator()(Ts&&... ts) {
        const auto start = std::chrono::high_resolution_clock::now();
        if constexpr (std::is_same_v<void, decltype(m_f(std::forward<Ts>(ts)...))>) {
            std::invoke(m_f, std::forward<Ts>(ts)...);
            const auto end = std::chrono::high_resolution_clock::now();
            auto dur = end - start;
            m_sum += dur;
            ++m_count;
        } else {
            auto out = m_f(std::forward<Ts>(ts)...);
            const auto end = std::chrono::high_resolution_clock::now();
            auto dur = end - start;
            m_sum += dur;
            ++m_count;
            return out;
        }
    }
    template <class... Ts>
    friend auto apply_time(Time_it_st& me, Ts&&... ts) {
        const auto start = std::chrono::high_resolution_clock::now();

        auto out = apply_F(me.m_f, std::forward<Ts>(ts)...);
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = end - start;
        me.m_sum += dur;
        ++me.m_count;
        return out;
    }

    auto& operator+=(Time_it_st const& other) {
        m_sum += other.m_sum;
        m_count += other.m_count;
        return *this;
    }

    auto mean_duration() const {
        return std::chrono::duration<double>(m_sum / std::max(m_count, 1ul)).count();
    }

    auto total_duration() const { return std::chrono::duration<double>(m_sum).count(); }

    auto count() const { return m_count; }

    void reset() {
        m_sum = std::chrono::nanoseconds::zero();
        m_count = 0ul;
    }

    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    inline friend std::ostream& report_title(std::ostream& os, const std::string& sep,
                                             Time_it_st const&) {
        os << sep << ToString(Id{}) << "_sum_time" << sep << ToString(Id{}) << "_count" << sep
           << ToString(Id{}) << "_average_time";
        return os;
    }

    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    friend std::ostream& report_point(std::ostream& os, const std::string& sep, Time_it_st& me) {
        os << sep << me.total_duration() << sep << me.count() << sep << me.mean_duration();
        me.reset();
        return os;
    }
};

template <class Id, class... Fun, template <class...> class F>
Time_it_st(F<Id, Fun...>) -> Time_it_st<Id, F<Id, Fun...>>;

inline std::ostream& report_title(std::ostream& os, const std::string&, ...) {
    return os;
}

inline std::ostream& report_point(std::ostream& os, const std::string&, ...) {
    return os;
}

template <class... Fs>
class FuncMap_St : public Fs... {
    std::string m_filename;
    std::string sep = ",";
    std::size_t m_sampling_interval;
    std::size_t m_max_number_of_values_per_iteration;

    std::ofstream m_file;

   public:
    using Fs::operator[]...;

    template <class Id>
    Nothing operator[](Id) const {
        return Nothing{};
    }

    void clear() { (clearit(static_cast<Fs&>(*this)), ...); }

    auto& file() const { return m_filename; }

    FuncMap_St(const std::string path, std::size_t sampling_interval,
               std::size_t max_number_of_values_per_iteration, Fs... fs)
        : Fs{std::move(fs)}...,
          m_filename{path},
          m_sampling_interval{sampling_interval},
          m_max_number_of_values_per_iteration{max_number_of_values_per_iteration} {}

    // FuncMap_St(const std::string path, std::size_t sampling_interval,
    //            std::size_t max_number_of_values_per_iteration, Fs const&... fs)
    //     : Fs{fs}...,
    //       m_filename{path},
    //       m_sampling_interval{sampling_interval},
    //       m_max_number_of_values_per_iteration{max_number_of_values_per_iteration} {}

    FuncMap_St create(const std::string& suffix) const {
        return FuncMap_St(file() + suffix, m_sampling_interval,
                          m_max_number_of_values_per_iteration,
                          static_cast<Fs const&>(*this).create()...);
    }

    auto to_bare_functions() const {
        return FuncMap_St<std::decay_t<decltype(static_cast<const Fs&>(*this).bare_function())>...>(
            m_filename, m_sampling_interval, m_max_number_of_values_per_iteration,
            static_cast<const Fs&>(*this).bare_function()...);
    }

    template <class Id, class... Ts>
   requires requires (FuncMap_St& me, Id id) {
        // 1. Check if the expression is valid
        { me[id] }; 
        
        // 2. Check the boolean condition (Nested Requirement)
        requires !std::is_same_v<std::remove_cvref_t<decltype(me[id])>, Nothing>;
    }
    auto f(Id, Ts&&... ts) {
        //     if constexpr(std::is_same_v<Nothing,decltype((*this)[Id{}])>)
        //         static_assert(false);
        auto& fun = (*this)[Id{}];
        return fun(*this, std::forward<Ts>(ts)...);
    }

    template <class Id, class... Ts>
    auto fstop(Id, Ts&&... ts) {
        auto& fun = (*this)[Id{}];
        return fun(std::forward<Ts>(ts)...);
    }

    template <class... Context_data>
    friend auto& report_point(FuncMap_St<Fs...>& f, std::size_t iter, const Context_data&... s) {
        if ((iter > 0) &&
            (iter % std::max(f.m_sampling_interval,
                             sizeof...(Fs) / f.m_max_number_of_values_per_iteration) ==
             0)) {
            if (!f.m_file.is_open()) {
                f.m_file.open(f.m_filename + "_funcmap.csv", std::ios::app);
                f.m_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
            }
            f.m_file << iter;
            ((f.m_file << f.sep << s), ...);
            (report_point(f.m_file, f.sep, static_cast<Fs&>(f)), ...);
            f.m_file << "\n";
            f.m_file.flush();
        }
        return f.m_file;
    }
    template <class... Context_string>
    friend auto& report_title(FuncMap_St<Fs...>& f, const std::string iter,
                              const Context_string&... s) {
        if (!f.m_file.is_open()) {
            f.m_file.open(f.m_filename + "_funcmap.csv");
            f.m_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        }
        f.m_file << iter;
        ((f.m_file << f.sep << s), ...);
        (report_title(f.m_file, f.sep, static_cast<Fs&>(f)), ...);
        f.m_file << "\n";
        f.m_file.flush();
        return f.m_file;
    }

    auto fork(std::size_t n) {
        std::vector<FuncMap_St> out;
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            out.emplace_back(this->create("_" + std::to_string(i)));
        }
        return out;
    }

    template <class... Fs2>
    auto append(Fs2 const&... fs2) const {
        return FuncMap_St<Fs..., Fs2...>(m_filename, m_sampling_interval,
                                         m_max_number_of_values_per_iteration,
                                         static_cast<Fs const&>(*this)..., fs2...);
    }

    auto& operator+=(const FuncMap_St& other) {
        ((static_cast<Fs&>(*this) += static_cast<Fs const&>(other)), ...);
        return *this;
    }
    auto& operator+=(const std::vector<FuncMap_St>& other) {
        for (auto& e : other) (*this) += e;
        return *this;
    }

    template <template <class...> class MacroR, class... Ps>
        requires(is_of_this_template_type_v<Ps, in_progress::P> && ...)
    auto append_Fs_S(in_progress::S<Ps...>) const {
        return append(make_F_P<MacroR>(Ps{})...);
    }

    template <template <class...> class MacroR, class... Ss>
        requires(is_of_this_template_type_v<Ss, in_progress::S> && ...)
    auto append_Fs(in_progress::P<Ss...>) const {
        auto ss = (... * Ss{});
        //using test=typename decltype(ss)::multi;
        return this->append_Fs_S<MacroR>(ss);
    }
};

template <class... Fs, class F, class G>
auto operator+(std::pair<FuncMap_St<Fs...>, F>&& f, G const& g) {
    if constexpr (std::is_same_v<typename G::myId, typename F::myId>)
        return std::pair(FuncMap_St<Fs..., F>(std::move(f.first.file()),
                                              std::move(f.first[typename Fs::myId{}])..., f.second),
                         f.second);
    else
        return std::pair(FuncMap_St<Fs..., G>(std::move(f.first.file()),
                                              std::move(f.first[typename Fs::myId{}])..., g),
                         f.second);
}

template <class G, class... Fs>
auto insert(const std::string& path, FuncMap_St<G, Fs...> const& fun) {
    return fun;
}

template <class G, class... Fs, class F, class... FFS>
auto insert(const std::string& path, FuncMap_St<G, Fs...> const& fun, F const& f,
            FFS const&... ff) {
    return insert(path,
                  (std::pair(FuncMap_St<G>(path, fun[typename G::myId{}]), f) + ... +
                   fun[typename Fs::myId{}])
                      .first,
                  ff...);
}

// -- Convenience helpers ----------------------------------------------------

// Factory: build an empty function map (no registered functions)
inline auto create_empty_function_map()
    -> FuncMap_St<>
{
    return FuncMap_St<>("filename", 1, 100);
}

// Compile-time detector: does the table provide F?
// Usage: if constexpr (var::has_it_defined<F>(f)) { ... }
template <class F,class FunctionTable>
constexpr bool has_it_defined() {
    return !std::is_same_v<Nothing,
                           decltype(std::declval<const FunctionTable&>()[F{}])>;
}

namespace partially_implemented {

template <class, class...>
class Test_it;

template <class Id, class... Fun, template <class...> class F, class Preconditions,
          class Postconditions>
    requires(std::is_trivial_v<Id>)
class Test_it<Id, F<Id, Fun...>, Preconditions, Postconditions> {
    F<Id, Fun...> m_f;
    std::vector<std::chrono::nanoseconds> m_sum;
    std::vector<std::size_t> m_count;
    std::size_t m_n_threads{};

   public:
    using myId = Id;

    static constexpr bool is_threadable = true;

    void clear(I_thread i) { clearit(m_f, i); }

    auto& get_Fun() { return m_f.get_Fun(); }

    constexpr Test_it(F<Id, Fun...>&& t_f, std::size_t n_threads)
        : m_f{std::move(t_f)},
          m_sum{n_threads, std::chrono::nanoseconds::zero()},
          m_count{std::vector<std::size_t>(n_threads, 0ul)},
          m_n_threads{n_threads} {}
    constexpr Test_it() = default;
    auto& operator[](Id) { return *this; }
    auto& operator[](Id) const { return *this; }

    auto naked_function() const { return m_f; }

    auto n_threads() const { return m_n_threads; }

    template <class... Ts>
    auto operator()(I_thread i_thread, Ts&&... ts) {
        const auto start = std::chrono::high_resolution_clock::now();
        if constexpr (m_f.is_threadable) {
            if constexpr (std::is_same_v<void, decltype(std::invoke(m_f, i_thread,
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
            if constexpr (std::is_same_v<void,
                                         decltype(std::invoke(m_f, std::forward<Ts>(ts)...))>) {
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
    friend auto apply_time(Test_it& me, I_thread i_thread, Ts&&... ts) {
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
        std::vector<decltype(std::chrono::duration<float, std::chrono::nanoseconds::period>(
            m_sum[0] / m_count[0]))>
            out(m_sum.size());
        for (std::size_t i = 0; i < n_threads(); ++i)
            out[i] = std::chrono::duration<float, std::chrono::nanoseconds::period>(
                m_sum[i] / std::max(m_count[i], 1ul));
        return out;
    }

    auto total_duration() const {
        auto out = m_sum[0];
        for (std::size_t i = 1; i < n_threads(); ++i) out += m_sum[i];
        return std::chrono::duration<double>(out).count();
    }
    auto total_count() const {
        auto out = m_count[0];
        for (std::size_t i = 1; i < n_threads(); ++i) out += m_count[i];
        return out;
    }

    auto count() const { return m_count; }

    void reset() {
        for (auto& e : m_sum) e = std::chrono::nanoseconds::zero();
        for (auto& e : m_count) e = 0ul;
    }

    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    inline friend std::ostream& report_title(std::ostream& os, const std::string& sep,
                                             Test_it const& me) {
        os << ToString(Id{}) << "_sum_time" << sep << ToString(Id{}) << "_count";
        for (std::size_t i = 0; i < me.m_count.size(); ++i)
            os << sep << ToString(Id{}) << "_" << i << "_sum_time" << sep << ToString(Id{}) << "_"
               << i << "_count" << sep << ToString(Id{}) << "_" << i << "_average_time";
        return os;
    }

    /**
   * @brief report_title
   * @param me
   * @param s
   *  Use outside parallel for
   */
    friend std::ostream& report_point(std::ostream& os, const std::string& sep, Test_it& me) {
        os << me.total_duration() << sep << me.total_count();
        for (std::size_t i = 0; i < me.m_count.size(); ++i)
            os << sep << std::chrono::duration<double>(me.m_sum[i]).count() << sep << me.count()[i]
               << sep << me.mean_duration()[i];
        me.reset();
        return os;
    }
};

}  // namespace partially_implemented

}  // namespace var
#endif  // FUNCTION_MEASURE_VERIFICATION_AND_OPTIMIZATION_H
