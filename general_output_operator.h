#ifndef GENERAL_OUTPUT_OPERATOR_H
#define GENERAL_OUTPUT_OPERATOR_H
#include <variant>
#include <ostream>


template<class ...Ts>
std::ostream& operator<<(std::ostream& os, std::variant<Ts...>const & x)
{
    return std::visit([&os](auto const& a)->std::ostream&{
        os<<a;
        return os;
    }, x);
}





#endif // GENERAL_OUTPUT_OPERATOR_H
