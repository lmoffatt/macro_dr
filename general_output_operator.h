#ifndef GENERAL_OUTPUT_OPERATOR_H
#define GENERAL_OUTPUT_OPERATOR_H
#include "maybe_error.h"
#include <istream>
#include <ostream>
#include <variant>

#include <map>
#include <optional>
#include <vector>
// #include <concepts>

template <class T>
  requires(std::is_arithmetic_v<T>)
std::ostream &print(std::ostream &os, T x) {
  os << x;
  return os;
}

inline std::ostream &print(std::ostream &os, const std::string &s) {
  return os << s;
}

template <class T>
std::ostream &operator<<(std::ostream &os, std::optional<T> const &x) {
  if (x)
    os << *x;
  return os;
}

template <class T>
std::ostream &print(std::ostream &os, std::optional<T> const &x) {
  if (x)
    print(os, *x);
  return os;
}

template <class... Ts>
std::ostream &print(std::ostream &os, std::variant<Ts...> const &x) {
  return std::visit(
      [&os](auto const &a) -> std::ostream & {
        print(os, a);
        return os;
      },
      x);
}

class septr : public std::string {
public:
  using std::string::string;
  //   separator(std::string s):std::string(std::move(s)){}
  septr(std::string x) : std::string(x) {}

  std::string operator()() const { return *this; }
  friend std::ostream &operator<<(std::ostream &os, const septr &se) {
    return os << se();
  }
  friend std::istream &operator>>(std::istream &is, const septr &se) {
    std::string ss = se();
    for (std::size_t i = 0; i < ss.size(); ++i) {
      is.get(ss[i]);
    }
    if (ss != se())
      is.setstate(std::ios::failbit);
    return is;
  }
};

class string_and_separator : public std::string {
  std::string m_sep;
   
public:
  using std::string::string;
  //   separator(std::string s):std::string(std::move(s)){}
  
  string_and_separator(const std::string& sep):m_sep{sep}{}
  
  std::string &operator()() { return *this; }
  std::string const &operator()() const { return *this; }
  friend std::ostream &operator<<(std::ostream &os,
                                  const string_and_separator &se) {
    return os << se() << se.m_sep;
  }
  friend std::istream &operator>>(std::istream &is, string_and_separator &se) {
    std::string string_candidate;
    std::string sep_candidate;
    char c;
    while ((sep_candidate != se.m_sep)&&(is.get(c)) ) {
      if (c == se.m_sep[sep_candidate.size()])
        sep_candidate.push_back(c);
      else
        string_candidate.push_back(c);
    }
    if (sep_candidate == se.m_sep) {
      se() = string_candidate;
    } else {
      is.setstate(std::ios::failbit);
    }
    return is;
  }
};

template <class K, class T>
std::ostream &operator<<(std::ostream &os, std::pair<K, T> const &x) {
  os << x.first << ", " << x.second << "\n";
  return os;
}

template <class K, class T>
std::istream &operator>>(std::istream &is, std::pair<K, T> &x) {
  is >> x.first >> septr(", ") >> x.second;
  return is;
}

template <class K, class T>
std::ostream &print(std::ostream &os, std::pair<K, T> const &x) {
  print(os, x.first);
  os << ", ";
  print(os, x.second);
  os << "\n";
  return os;
}

template <class K, class T>
std::ostream &operator<<(std::ostream &os, std::map<K, T> const &x) {
  for (auto it = x.begin(); it != x.end(); ++it)
    os << it->first << "--> " << it->second << "\n";
  return os;
}

template <class K, class T>
std::istream &operator>>(std::istream &is, std::map<K, T> &x) {
  K k;
  T e;
  while (is >> k >> septr("--> ") >> e >> septr("\n"))
    x.insert({k, e});

  if (!x.empty())
    is.setstate(std::ios::goodbit);
  return is;
}

template <class K, class T>
std::ostream &print(std::ostream &os, std::map<K, T> const &x) {
  for (auto it = x.begin(); it != x.end(); ++it) {
    print(os, it->first);
    os << "--> ";
    print(os, it->second);
    os << "\n";
  }
  os << "\n";
  return os;
}

template <class Ts>
std::istream &operator>>(std::istream &is, std::vector<Ts> &v) {
  Ts e;
  while (is >> e) {
    v.push_back(e);
    if constexpr (has_size<Ts>)
      is >> septr("\n ");
    else
      is >> septr(", ");
  }
  if (!v.empty())
    is.setstate(std::ios::goodbit);
  return is;
}

template <class T, class... Ts>
std::ostream &operator<<(std::ostream &os, const std::tuple<T, Ts...> &tu) {
  return std::apply(
      [&os](T const &x, Ts const &...xs) -> std::ostream & {
        ((os << x), ..., (os << ", " << xs));
        return os;
      },
      tu);
}

template <class T, class... Ts>
std::istream &operator>>(std::istream &is, std::tuple<T, Ts...> &tu) {
  return std::apply(
      [&is](T &x, Ts &...xs) -> std::istream & {
        ((is >> x), ..., (is >> septr(", ") >> xs));
        return is;
      },
      tu);
}

template <class... Ts>
std::ostream &print(std::ostream &os, const std::tuple<Ts...> &tu) {
  return std::apply(
      [&os](auto &...x) -> std::ostream & {
        ((print(os, x), os << "\n"), ...);
        os << "\n";
        return os;
      },
      tu);
}

template <class Ts>
  requires(!has_size<Ts>)
std::ostream &operator<<(std::ostream &os, const std::vector<Ts> &v) {
  for (std::size_t i = 0; i < v.size(); ++i)
    os << v[i] << ", ";
  //  os<<"\n";
  return os;
}

template <class Ts>
  requires(has_size<Ts>)
std::ostream &operator<<(std::ostream &os, const std::vector<Ts> &v) {
  for (std::size_t i = 0; i < v.size(); ++i)
    os << v[i] << "\n";
  return os;
}

template <class Ts>
std::ostream &print(std::ostream &os, const std::vector<Ts> &v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    print(os << i << ":\n", v[i]);
    os << "\n";
  }
  os << "\n";
  return os;
}

template <class... Ts>
std::ostream &operator<<(std::ostream &os, std::variant<Ts...> const &x) {
  return std::visit(
      [&os](auto const &a) -> std::ostream & {
        os << a;
        return os;
      },
      x);
}

template <class... Ts>
std::istream &operator>>(std::istream &is, std::variant<Ts...> &x) {
  return std::visit(
      [&is](auto &a) -> std::istream & {
        is >> a;
        return is;
      },
      x);
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Maybe_error<T> &x) {
  if (x)
    os << x.value();
  else
    os << x.error()();
  return os;
}

#endif // GENERAL_OUTPUT_OPERATOR_H
