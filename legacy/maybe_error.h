#ifndef MAYBE_ERROR_H
#define MAYBE_ERROR_H

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

inline std::string git_Commit_Hash() {
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "0000000"  // 0000000 means uninitialized
#endif
    return GIT_COMMIT_HASH;
}

inline std::string leadingZero(int i) {
    if (i == 0)
        return "00";
    else if (i < 10)
        return "0" + std::to_string(i);
    else
        return std::to_string(i);
}

inline std::string time_now() {
    auto tc = std::chrono::system_clock::now();
    std::time_t rawtime = std::chrono::system_clock::to_time_t(tc);

    auto tcount = (tc.time_since_epoch().count() / 1000) % 1000000;

    struct std::tm* t;
    time(&rawtime);
    t = localtime(&rawtime);
    return git_Commit_Hash() + "_" + leadingZero(t->tm_hour) + leadingZero(t->tm_min) +
           leadingZero(t->tm_sec) + "s" + std::to_string(tcount);
}

constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();
constexpr auto eps = std::numeric_limits<double>::epsilon();

template <typename C, template <typename...> class V>
struct is_of_this_template_type : std::false_type {};
template <template <typename...> class V, typename... Ts>
struct is_of_this_template_type<V<Ts...>, V> : std::true_type {};
template <typename C, template <typename...> class V>
inline constexpr bool is_of_this_template_type_v =
    is_of_this_template_type<std::decay_t<C>, V>::value;

template <typename C, template <typename, auto...> class V>
struct is_of_this_template_value_type : std::false_type {};

template <template <typename, auto...> class V, typename T, auto Xs>
struct is_of_this_template_value_type<V<T, Xs>, V> : std::true_type {};

template <typename C, template <typename, auto...> class V>
inline constexpr bool is_of_this_template_value_type_v =
    is_of_this_template_value_type<C, V>::value;

template <class T>
struct return_type;

template <class R, class... A>
struct return_type<R (*)(A...)> {
    using type = R;
};

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
// explicit deduction guide (not needed as of C++20)
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

template <auto F>
std::string function_name();

struct Nothing {};

template <class T>
concept has_size = requires(const T& a) {
    { a.size() } -> std::convertible_to<std::size_t>;
};

template <class T>
concept index_accesible = requires(T a) {
    { a[std::size_t{}] } -> std::same_as<double&>;
};

template <class T>
    requires(!std::is_reference_v<T>)
class Maybe_error;

template <class T>
concept is_Maybe_error = is_of_this_template_type_v<std::decay_t<T>, Maybe_error>;

template <class T>
struct make_Maybe_error {
    using type = std::conditional_t<is_Maybe_error<T>, T, Maybe_error<T>>;
};

template <class T>
using Maybe_error_t = typename make_Maybe_error<T>::type;

template <class T>
using Maybe_unique = Maybe_error<std::unique_ptr<T>>;

template <class T>
constexpr bool is_valid(const Maybe_error<T>& v) {
    return v.valid();
}

template <class T>
    requires(!is_Maybe_error<T>)
constexpr bool is_valid(const T&) {
    return true;
}

template <class T>
    requires(!is_Maybe_error<T>)
auto& get_value(const T& x) {
    return x;
}

template <class T>
auto& get_value(const Maybe_error<T>& x) {
    return x.value();
}

template <class T>
T get_value(Maybe_error<T>&& x) {
    return std::move(x.value());
}

template <class C, class T>
concept contains_value = std::convertible_to<decltype(get_value(std::declval<C>())), T>;

template <class T>
auto get_error(const Maybe_error<T>& x) {
    return x.error();
}

class error_message {
    std::string m_;

   public:
    error_message() = default;
    error_message(std::string&& error) : m_{std::move(error)} {}

    template <class... T>
    error_message(T&&... args) {
        std::ostringstream os;
        (os << ... << std::forward<T>(args));
        m_ = os.str();
    }

    auto operator()() const { return m_; }
};
template <class T>
    requires(!is_Maybe_error<T>)
auto get_error(const T&) {
    return error_message("");
}

namespace err {

template <class T>
class Maybe_error : public std::variant<T, error_message> {
   public:
    using std::variant<T, error_message>::variant;
    constexpr explicit operator bool() const noexcept { return this->index() == 0; }

    constexpr const T* operator->() const noexcept { return std::get_if<0>(*this); }

    constexpr T* operator->() noexcept { return std::get_if<0>(*this); }

    constexpr const T& operator*() const& noexcept { return std::get<0>(*this); }
    constexpr T& operator*() & noexcept { return std::get<0>(*this); }
    constexpr const T&& operator*() const&& noexcept { return std::get<0>(std::move(*this)); }
    constexpr T&& operator*() && noexcept { return std::get<0>(std::move(*this)); }
    Maybe_error(const T& t) : std::variant<T, error_message>(t) {}
    Maybe_error(T&& t) : std::variant<T, error_message>(std::move(t)) {}
    Maybe_error(error_message m) : std::variant<T, error_message>(m) {}
    error_message error() const {
        if (*this)
            return error_message("");
        else
            return std::get<1>(*this);
    }
};
}  // namespace err

template <class T>
Maybe_unique<T> operator||(Maybe_unique<T>&& one, Maybe_unique<T>&& two) {
    if (one)
        return std::move(one);
    else if (two)
        return std::move(two);
    else {
        std::string error = one.error()() + "\n" + two.error()();
        return error_message(error);
    }
}

inline std::string ToString(const double& x) {
    std::stringstream ss;
    ss << std::setprecision(12);
    ss << x;
    return ss.str();
}

template <class T>
    requires(!std::is_reference_v<T>)
class Maybe_error : private std::variant<T, error_message> {
   public:
    using base_type = std::variant<T, error_message>;
    using value_type = T;
    template <typename U>
        requires std::constructible_from<T, U>
    Maybe_error(U&& u) : base_type(T(std::forward<U>(u))) {}

    template <typename U>
        requires std::is_same_v<T, std::reference_wrapper<std::remove_reference_t<U>>>
    Maybe_error(U& u) : base_type(std::ref(u)) {}

    Maybe_error() = default;
    Maybe_error(error_message&& x) : base_type{(std::move(x))} {}
    Maybe_error(const error_message& x) : base_type{x} {}
    constexpr explicit operator bool() const noexcept { return this->index() == 0; }
    [[nodiscard]] constexpr bool valid() const noexcept { return this->index() == 0; }

    constexpr const T& value() const& noexcept { return std::get<0>(*this); }
    constexpr T& value() & noexcept { return std::get<0>(*this); }
    constexpr const T&& value() const&& noexcept { return std::get<0>(std::move(*this)); }
    constexpr T&& value() && noexcept { return std::get<0>(std::move(*this)); }
    error_message error() const {
        if (*this)
            return error_message("");
        else
            return std::get<1>(*this);
    }

    // friend std::ostream &operator<<(std::ostream &os, const Maybe_error &x) {
    //   if (x)
    //     os << x.value();
    //   else
    //     os << x.error()();
    //   return os;
    // }
    // friend std::ostream &operator<<(std::ostream &os, Maybe_error &&x) {
    //   if (x)
    //     os << x.value();
    //   else
    //     os << x.error()();
    //   return os;
    // }
};

template <>
class Maybe_error<void> : private std::variant<std::monostate, error_message> {
    using variant_type = std::variant<std::monostate, error_message>;

   public:
    // success
    Maybe_error() : variant_type{std::in_place_index<0>} {}
    // failure
    Maybe_error(error_message err) : variant_type{std::in_place_index<1>, std::move(err)} {}

    constexpr explicit operator bool() const noexcept { return this->index() == 0; }
    [[nodiscard]] constexpr bool valid() const noexcept { return this->index() == 0; }

    constexpr void value() const& noexcept {}

    [[nodiscard]] error_message error() const {
        if (*this) {
            return {""};
        }
        return std::get<1>(*this);
    }
};

template <class T>
Maybe_error(T&&) -> Maybe_error<T>;

template <class T>
Maybe_error(T const&) -> Maybe_error<T>;

template <class T>
    requires std::is_object_v<typename T::base_type>
class Maybe_error<T*> : public Maybe_error<typename T::base_type*> {
    T* m_value;

   public:
    using base_type = Maybe_error<typename T::base_type*>;
    using base_type::operator bool;
    using base_type::error;
    using base_type::valid;

    Maybe_error() = default;
    Maybe_error(error_message&& x) : base_type{(std::move(x))} {}
    Maybe_error(const error_message& x) : base_type{x} {}
    Maybe_error(T* t) : base_type{t}, m_value{t} {}

    constexpr const T* value() const noexcept { return m_value; }
    constexpr T* value() { return m_value; }

    friend std::ostream& operator<<(std::ostream& os, const Maybe_error& x) {
        if (x)
            os << x.value();
        else
            os << x.error()();
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, Maybe_error&& x) {
        if (x)
            os << x.value();
        else
            os << x.error()();
        return os;
    }
};

// Your wrapper:
// template<class T> struct Maybe_error { using value_type = T; /* ... */ };

namespace detail {

// Base case: not a Maybe_error → underlying type is T itself
template <class T>
struct underlying_value_type_impl {
    using type = T;
};

// Specialization: unwrap Maybe_error<U> → underlying type is U
template <class U>
struct underlying_value_type_impl<Maybe_error<U>> {
    using type = U;
};

// (Optional) make it robust to cv/ref
template <class T>
using decay_maybe_t = std::remove_cvref_t<T>;

}  // namespace detail

// Public alias: T or U if T == Maybe_error<U>
template <class T>
using underlying_value_type_t =
    typename detail::underlying_value_type_impl<detail::decay_maybe_t<T>>::type;

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto operator*(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(get_value(x) * get_value(y));
    else
        return Maybe_error<std::decay_t<decltype(get_value(x) * get_value(y))>>(
            get_error(x)() + " multiplies " + get_error(y)());
}

inline Maybe_error<bool> operator>>(std::string&& context_message, Maybe_error<bool>&& x) {
    if (x)
        return x;
    else
        return error_message(context_message + x.error()());
}

inline Maybe_error<bool> operator==(Maybe_error<bool> const& one, Maybe_error<bool> const& two) {
    if ((one) && (two))
        return one.value() == two.value();
    else if ((!one) && (!two))
        return one.error()() == two.error()();
    else
        return false;
}

inline Maybe_error<bool> operator&&(Maybe_error<bool>&& one, Maybe_error<bool>&& two) {
    if ((one) && (two))
        return true;
    else
        return error_message(one.error()() + two.error()());
}

template <class T>
Maybe_error<std::vector<T>> promote_Maybe_error(std::vector<Maybe_error<T>> const& x) {
    std::vector<T> out(size(x));

    for (std::size_t i = 0; i < size(out); ++i)
        if (x[i])
            out[i] = x[i].value();
        else
            return error_message("error at " + std::to_string(i) + x[i].error()());
    return out;
}

inline Maybe_error<bool> consolidate(std::vector<Maybe_error<bool>> const& x) {
    Maybe_error<bool> err;
    for (std::size_t i = 0; i < size(x); ++i)
        if (!x[i]) {
            err.value() = false;
            err.error()() += x[i].error()();
        }
    return err;
}

/*
template <class... Ts>
Maybe_error<std::tuple<Ts...>>
promote_Maybe_error(std::tuple<Maybe_error<Ts>...> &&x) {
    return std::apply(
        [](auto &&...e) -> Maybe_error<std::tuple<Ts...>> {
            if ((e.valid() && ... && true))
                return std::tuple(std::move(e.value())...);
            else
                return error_message((e.error()() + ... + ""));
        },
        std::move(x));
}
*/
inline Maybe_error<std::tuple<>> promote_Maybe_error(std::tuple<>) {
    return {};
}

template <class... Ts>
Maybe_error<std::tuple<std::unique_ptr<Ts>...>> promote_Maybe_error(
    std::tuple<Maybe_unique<Ts>...>&& x) {
    return std::apply(
        [](auto&&... e) -> Maybe_error<std::tuple<std::unique_ptr<Ts>...>> {
            if ((e.valid() && ... && true))
                return std::tuple(std::move(e.value())...);
            else
                return error_message((e.error()() + ... + ""));
        },
        std::move(x));
}

template <class... Ts>
    requires(((!is_of_this_template_type_v<Ts, std::unique_ptr>)) && ...)
Maybe_error<std::tuple<Ts...>> promote_Maybe_error(std::tuple<Maybe_error<Ts>...>&& x) {
    return std::apply(
        [](auto&&... e) -> Maybe_error<std::tuple<Ts...>> {
            if ((e.valid() && ... && true)) {
                return std::tuple(std::move(e.value())...);
            }
            return error_message((e.error()() + ... + ""));
        },
        std::move(x));
}

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto operator+(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(get_value(x) + get_value(y));
    else
        return Maybe_error<std::decay_t<decltype(get_value(x) + get_value(y))>>(
            error_message(get_error(x)() + " multiplies " + get_error(y)()));
}

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto operator-(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(get_value(x) - get_value(y));
    else {
        if (!get_error(x)().empty() || !get_error(y)().empty())

            return Maybe_error<std::decay_t<decltype(get_value(x) - get_value(y))>>(
                get_error(x)() + " multiplies " + get_error(y)());
        else
            return Maybe_error<std::decay_t<decltype(get_value(x) - get_value(y))>>(
                error_message{});
    }
}

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto operator/(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(get_value(x) / get_value(y));
    else {
        if (!get_error(x)().empty() || !get_error(y)().empty())

            return Maybe_error<std::decay_t<decltype(get_value(x) / get_value(y))>>(
                get_error(x)() + " divides " + get_error(y)());
        else
            return Maybe_error<std::decay_t<decltype(get_value(x) / get_value(y))>>(
                error_message{});
    }
}

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto xtAx(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(xtAx(get_value(x), get_value(y)));
    else
        return Maybe_error<std::decay_t<decltype(xtAx(get_value(x), get_value(y)))>>(
            get_error(x) + " xtAx " + get_error(y));
}

template <class T, class S>
    requires(is_Maybe_error<T> || is_Maybe_error<S>)
auto xAxt(const T& x, const S& y) {
    if (is_valid(x) && is_valid(y))
        return Maybe_error(xAxt(get_value(x), get_value(y)));
    else
        return Maybe_error<std::decay_t<decltype(xtAx(get_value(x), get_value(y)))>>(
            get_error(x)() + " xAxt " + get_error(y)());
}

template <class T, class... Ts>
    requires(std::is_same_v<T, Ts> || ...)
Maybe_error<std::variant<Ts...>> make_Maybe_variant(Maybe_error<T>&& x) {
    if (x)
        return std::move(x);
    else
        return x.error();
}

template <class T>
    requires(is_Maybe_error<T>)
auto inv(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(inv(x.value()))>>;
    if (x)
        return return_type(inv(x.value()));
    else
        return return_type(x.error()() + "\n inverse");
}
template <class T>
    requires(is_Maybe_error<T>)
auto diag(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(diag(x.value()))>>;
    if (x)
        return return_type(diag(x.value()));
    else
        return return_type(x.error()() + "\n diag");
}
template <class T>
    requires(is_Maybe_error<T>)
auto XXT(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(XXT(x.value()))>>;
    if (x)
        return return_type(XXT(x.value()));
    else
        return return_type(x.error()() + "\n diag");
}

template <class T>
    requires(is_Maybe_error<T>)
auto operator-(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(-(x.value()))>>;
    if (x)
        return return_type(-(x.value()));
    else
        return return_type(x.error()() + "\n minus");
}

template <class T>
    requires(is_Maybe_error<T>)
auto logdet(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(logdet(x.value()))>>;
    if (x)
        return return_type(logdet(x.value()));
    else
        return return_type(x.error()() + "\n logdet");
}

using std::log;
template <class T>
    requires(is_Maybe_error<T>)
auto log(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(log(x.value()))>>;
    if (x)
        return return_type(log(x.value()));
    else
        return return_type(x.error()() + "\n log");
}

template <class T>
    requires(is_Maybe_error<T>)
auto cholesky(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(cholesky(x.value()))>>;
    if (x)
        return return_type(cholesky(x.value()));
    else
        return return_type(x.error()() + "\n log");
}
using std::sqrt;
template <class T>
    requires(is_Maybe_error<T>)
auto sqrt(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(sqrt(x.value()))>>;
    if (x)
        return return_type(sqrt(x.value()));
    else
        return return_type(x.error()() + "\n sqrt");
}

template <class T>
    requires(is_Maybe_error<T>)
auto tr(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(tr(x.value()))>>;
    if (x)
        return return_type(tr(x.value()));
    else
        return return_type(error_message(x.error()() + "\n transpose"));
}

template <class T>
    requires(is_Maybe_error<T>)
auto Trace(const T& x) {
    using return_type = Maybe_error_t<std::decay_t<decltype(Trace(x.value()))>>;
    if (x)
        return return_type(Trace(x.value()));
    else
        return return_type(x.error()() + "\n Trace");
}

template <class T, auto... F>
struct return_error {
    std::string function_label;
    Maybe_error<T> operator()(std::string&& error) {
        return Maybe_error<T>(function_label +
                              (function_name<F>() + ... + (": " + std::move(error))));
    }
};

template <class T>
    requires(!std::is_object_v<T>)
std::ostream& put(std::ostream& os, const T& t) {
    os << t << "\n";
    return os;
}

inline std::string file_contents(const std::string& x) {
    std::ifstream f0(x);
    std::stringstream b0;
    b0 << f0.rdbuf();

    return b0.str();
}

inline Maybe_error<bool> compare_file_contents(const std::string& one, const std::string& two,
                                               std::size_t max_errors = 10) {
    auto s0 = file_contents(one);
    auto s1 = file_contents(two);
    if (s0 == s1)
        return true;
    std::stringstream ss0(s0);
    std::stringstream ss1(s1);
    std::string line0;
    std::string line1;
    std::size_t i = 0;
    std::size_t n_errors = 0;
    std::string message;
    while ((n_errors < max_errors) && (std::getline(ss0, line0)) && (std::getline(ss1, line1))) {
        if (line0 != line1) {
            message += "\nline " + std::to_string(i) + ":\n\t" + line0 + "\n\t" + line1 + "\n";

            ++n_errors;
        }
        ++i;
    }
    return error_message(message);
}

#endif  // MAYBE_ERROR_H
