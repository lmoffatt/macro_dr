#ifndef GRAMMAR_IDENTIFIER_H
#define GRAMMAR_IDENTIFIER_H

#include <map>
#include <memory>
#include <string>
#include <string_view>

#include "maybe_error.h"

namespace macrodr::dsl {

template <class Lexer>
class Identifier;

template <class Lexer>
Maybe_error<Identifier<Lexer>> to_Identifier(std::string t_id_candidate);

template <class Lexer>
class Identifier {
    std::string m_id;
    Identifier(std::string&& t_id) : m_id{std::move(t_id)} {
    }

   public:
    friend Maybe_error<Identifier> to_Identifier<Lexer>(std::string t_id_candidate);
    auto operator()() const {
        return m_id;
    }

    friend bool operator<(const Identifier& id, const Identifier& id2) {
        return id() < id2();
    }
};

template <class Lexer>
Maybe_error<Identifier<Lexer>> to_Identifier(std::string t_id_candidate) {
    constexpr std::string_view kAlpha =
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";
    constexpr std::string_view kAlphanum =
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_";
    if (t_id_candidate.empty()) {
        return error_message("empty identifier");
    } if (kAlpha.find(t_id_candidate[0]) == std::string::npos)
        return error_message(t_id_candidate + " identifier starts with " + t_id_candidate[0]);
    else {
        auto pos = t_id_candidate.find_first_not_of(kAlphanum);
        if (pos == std::string::npos)
            return Identifier<Lexer>(std::move(t_id_candidate));
        else
            return error_message(t_id_candidate + " has a " + t_id_candidate[pos] + "at " +
                                 std::to_string(pos) + " position");
    }
}

template <class Abstract>
auto clone_vector(const std::vector<std::unique_ptr<Abstract>>& x) {
    std::vector<std::unique_ptr<Abstract>> out;
    out.reserve(x.size());
for (auto& ele : x) { out.push_back(std::unique_ptr<Abstract>(ele->clone()));
}
    return out;
}

template <class key, class Abstract>
auto clone_map(const std::map<key, std::unique_ptr<Abstract>>& x) {
    std::map<key, std::unique_ptr<Abstract>> out;
    for (auto& ele : x) { out.emplace(ele.first, ele.second->clone());
}
    return out;
}

}  // namespace macrodr::dsl

#endif  // GRAMMAR_IDENTIFIER_H

#pragma once
