#ifndef GRAMMAR_IDENTIFIER_H
#define GRAMMAR_IDENTIFIER_H

#include "Maybe_error.h"


#include <string>
namespace dcli {

using logic::Maybe_error;
using logic::error_message;

template<class Lexer>
class Identifier;
template<class Lexer>
Maybe_error<Identifier<Lexer>> to_Identifier(const std::string& t_id_candidate);
template<class Lexer>
class Identifier {
    std::string m_id;
    Identifier(std::string t_id) : m_id{t_id} {}

public:
    friend Maybe_error<Identifier>
    to_Identifier<Lexer>(const std::string& t_id_candidate);
    auto operator()() const { return m_id; }
};

template<class Lexer>
Maybe_error<Identifier<Lexer>> to_Identifier(const std::string& t_id_candidate) {
    if (t_id_candidate.empty())
        return error_message("empty identifier");
    else if (Lexer::alfa.find(t_id_candidate[0]) == std::string::npos)
        return error_message(t_id_candidate + " identifier starts with " +
                             t_id_candidate[0]);
    else {
        auto pos = t_id_candidate.find_first_not_of(Lexer::alfanum);
        if (pos == std::string::npos)
            return Identifier<Lexer>(t_id_candidate);
        else
            return error_message(t_id_candidate + " has a " + t_id_candidate[pos] +
                                 "at " + std::to_string(pos) + " position");
    }
}


} // namespace dcli


#endif // GRAMMAR_IDENTIFIER_H

#pragma once
