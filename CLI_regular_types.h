#pragma once
#ifndef CLI_REGULAR_TYPES_H
#define CLI_REGULAR_TYPES_H

#include "grammar_Identifier.h"


namespace dcli {

template <typename T> struct T_s;

template <> struct T_s<int> {
    Identifier name() const { return *to_Identifier("integer"); }
};

template <class... Ts> struct Cs {};


} // namespace dcli


#endif // CLI_REGULAR_TYPES_H

