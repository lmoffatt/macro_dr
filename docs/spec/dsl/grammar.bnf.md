
---
# MacroIR DSL Grammar Specification (BNF)
# Status: DRAFT
# Last Updated: 2025-11-14
# Scope: untyped input syntax (prior to compilation)

```bnf

<program>         ::= <statement>*
# A program is a sequence of statements.
# Example:
#   simulate(...)
#   a = simulate(...)

<statement>       ::= <assignment>
                    | <function_call>
                    | <identifier>
                    | <literal>
                    | <string_literal>
                    | <argument_list>
                    | <vector_literal>
                    | <tuple_literal>
# Each line is parsed independently into one of the statement types.
# Vector and tuple literals are allowed as top-level statements
# for REPL-style exploration.
# Typing notes (informal, see type_rules.md for full details):
# - A `[e1, e2, ...]` literal has type `std::vector<T>` if all `ei` have type `T`.
# - A `{e1, e2, ...]` literal has type `std::tuple<T1, T2, ...>` where each `Ti`
#   is the type of `ei`.
# - Vector and tuple literals can appear anywhere an <expression> is allowed,
#   including as function arguments and in assignments.
# - Elements that are compiled as references (`T&`, `const T&`) must be
#   identifier expressions; literals and inline function calls are rejected
#   for such positions.

<assignment>      ::= <identifier> "=" <expression>
# Variable binding
# .str() → "a = simulate(...)"

<function_call>   ::= <identifier> "(" <argument_list>? ")"
# Function invocation with optional arguments
# .str() → "simulate(model, parameters)"

<argument_list>   ::= <argument> ("," <argument>)*
# Comma-separated arguments, which can be expressions or named assignments
# .str() → "model=lin, parameters=values"

<argument>        ::= <assignment>
                    | <expression>
# Accepts either named or positional argument

<expression>      ::= <literal>
                    | <string_literal>
                    | <identifier>
                    | <function_call>
                    | <argument_list>
                    | <vector_literal>
                    | <tuple_literal>
# General expression node (used in RHS). Now includes vector and
# tuple literals as first-class expressions.

<vector_literal>  ::= "[" <expression_list>? "]"
# Homogeneous collections (vectors/sets) are written with square brackets.
# Example:
#   [1, 2, 3]
#   [likelihood(m1), likelihood(m2)]

<tuple_literal>   ::= "{" <expression_list>? "}"
# Heterogeneous positional collections (tuples/pairs) are written with braces.
# Example:
#   {model1, param_set1}
#   {0.1, "fast"}

<expression_list> ::= <expression> ("," <expression>)*
# Generic comma-separated list of expressions used inside vectors and tuples.

<literal>         ::= <number>
<number>          ::= [0-9]+ ("." [0-9]+)?
# Numeric literals (int or float)
# .str() → "42" or "3.14"

<string_literal>  ::= "\"" .*? "\""
# Quoted string literal
# .str() → "\"data.csv\""

<identifier>      ::= [a-zA-Z_][a-zA-Z0-9_]*
# Function names and variable references
# .str() → "simulate", "Qdt", "beta1"

