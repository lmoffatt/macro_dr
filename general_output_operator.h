#ifndef GENERAL_OUTPUT_OPERATOR_H
#define GENERAL_OUTPUT_OPERATOR_H
#include <variant>
#include <ostream>
#include <map>
#include <optional>
#include <vector>
//#include <concepts>

template<class T>
    requires (std::is_arithmetic_v<T>)
std::ostream& print(std::ostream& os, T x)
{
    os<<x;
    return os;
}

std::ostream& print(std::ostream& os, const std::string& s)
{
    return os<<s;   
}




template<class T>
std::ostream& operator<<(std::ostream& os, std::optional<T>const & x)
{
    if (x)
        os<<*x;
    return os;
    
}

template<class T>
std::ostream& print(std::ostream& os, std::optional<T>const & x)
{
    if (x)
        print(os,*x);
    return os;
    
}



template<class ...Ts>
std::ostream& print(std::ostream& os, std::variant<Ts...>const & x)
{
    return std::visit([&os](auto const& a)->std::ostream&{
        print(os,a);
        return os;
    }, x);
}




template<class K, class T>
std::ostream& operator<<(std::ostream& os, std::pair<K,T>const & x)
{
    os<<x.first<<", "<<x.second<<"\n";
    return os;
    
}

template<class K, class T>
std::ostream& print(std::ostream& os, std::pair<K,T>const & x)
{
    print(os,x.first);
    os<<", ";
    print(os,x.second);
    os<<"\n";
    return os;
    
}



template<class K, class T>
std::ostream& operator<<(std::ostream& os, std::map<K,T>const & x)
{
    for (auto it=x.begin(); it!=x.end(); ++it)
        os<<it->first<<"--> "<<it->second<<"\n";
    return os;
    
}

template<class K, class T>
std::ostream& print(std::ostream& os, std::map<K,T>const & x)
{
    for (auto it=x.begin(); it!=x.end(); ++it)
    {    print(os,it->first);
    os<<"--> ";
    print(os,it->second);
    os<<"\n";
    }
    os<<"\n";
    return os;
    
}


template <class... Ts>
std::ostream &operator<<(std::ostream &os, const std::tuple<Ts...> &tu) {
    return std::apply(
        [&os](auto &...x) -> std::ostream & {
            (put(os, x), ...);
            return os;
        },
        tu);
}

template <class... Ts>
std::ostream &print(std::ostream &os, const std::tuple<Ts...> &tu) {
    return std::apply(
        [&os](auto &...x) -> std::ostream & {
            ((print(os, x), os<<"\n"),...);
            os<< "\n";
            return os;
        },
        tu);
}



template <class Ts>
std::ostream & operator<<(std::ostream &os, const std::vector<Ts> &v) {
    for (std::size_t i = 0; i < v.size(); ++i)
        os << v[i] << "\n";
    //  os<<"\n";
    return os;
}



template <class Ts>
std::ostream &print(std::ostream &os, const std::vector<Ts> &v) {
    for (std::size_t i = 0; i < v.size(); ++i)
    {print(os,v[i]);
        os<< "\n";
    }
    os<< "\n";
    return os;
}

template<class ...Ts>
std::ostream& operator<<(std::ostream& os, std::variant<Ts...>const & x)
{
    return std::visit([&os](auto const& a)->std::ostream&{
        os<<a;
        return os;
    }, x);
}


#endif // GENERAL_OUTPUT_OPERATOR_H
