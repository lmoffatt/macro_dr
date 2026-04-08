
---
# MacroIR DSL Grammar Specification (BNF)
# Status: WORKING SPEC
# Last Updated: 2026-03-24
# Scope: current untyped input syntax as accepted by the lexer/parser

```bnf

<program>         ::= <statement> ( <newline> <statement> )* <newline>*
# A program is newline-separated.
# Example:
#   simulate(...)
#   a = simulate(...)

<statement>       ::= <comment>
                    | <assignment>
                    | <function_call>
                    | <identifier>
                    | <literal>
                    | <string_literal>
                    | <argument_list>
                    | <vector_literal>
                    | <tuple_literal>
# Each line is parsed independently into one of the statement types.
# `// ...` comments are valid statements in the current lexer.
# Vector and tuple literals are allowed as top-level statements
# for REPL-style exploration.
# Typing notes (informal, see type_rules.md for full details):
# - A `[e1, e2, ...]` literal has type `std::vector<T>` if all `ei` have type `T`.
# - A `{e1, e2, ...}` literal has type `std::tuple<T1, T2, ...>` where each `Ti`
#   is the type of `ei`.
# - Vector and tuple literals can appear anywhere an <expression> is allowed,
#   including as function arguments and in assignments.
# - Elements that are compiled as references (`T&`, `const T&`) must be
#   identifier expressions; literals and inline function calls are rejected
#   for such positions.

<comment>         ::= "//" <comment_text>
# Line comment. Treated as a statement node by the current lexer.

<assignment>      ::= <identifier> <ws>* "=" <ws>* <expression>
# Variable binding
# .str() → "a = simulate(...)"

<function_call>   ::= <identifier> "(" <wsline>* <argument_list>? <wsline>* ")"
# Function invocation with optional arguments
# .str() → "simulate(model, parameters)"

<argument_list>   ::= <argument> ( <wsline>* "," <wsline>* <argument> )*
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

<vector_literal>  ::= "[" <wsline>* <expression_list>? <wsline>* "]"
# Homogeneous collections (vectors/sets) are written with square brackets.
# Example:
#   [1, 2, 3]
#   [likelihood(m1), likelihood(m2)]

<tuple_literal>   ::= "{" <wsline>* <expression_list>? <wsline>* "}"
# Heterogeneous positional collections (tuples/pairs) are written with braces.
# Example:
#   {model1, param_set1}
#   {0.1, "fast"}

<expression_list> ::= <expression> ( <wsline>* "," <wsline>* <expression> )*
# Generic comma-separated list of expressions used inside vectors and tuples.

<literal>         ::= <number>
<number>          ::= <sign>? <unsigned_number> <exponent_part>?
<unsigned_number> ::= <digits> ( "." <digits> )?
<exponent_part>   ::= ( "E" | "e" ) <sign>? <digits>
<sign>            ::= "+" | "-"
# Numeric literals may be signed and may use scientific notation.
# Examples:
#   42
#   -3.14
#   +2.0e-3

<string_literal>  ::= <double_quoted_string>
                    | <single_quoted_string>
<double_quoted_string> ::= "\"" <double_quoted_char>* "\""
<single_quoted_string> ::= "'" <single_quoted_char>* "'"
# The current lexer accepts both quote styles and stores the inner text.

<identifier>      ::= <identifier_start> <identifier_char>*
<identifier_start>::= [a-zA-Z_]
<identifier_char> ::= [a-zA-Z0-9_]
# Function names and variable references
# .str() → "simulate", "Qdt", "beta1"

<digits>          ::= [0-9]+
<newline>         ::= "\n"
<ws>              ::= " " | "\t"
<wsline>          ::= <ws> | <newline>
<comment_text>    ::= any characters up to end-of-line

```

## Notes

- This grammar is now meant to describe the current untyped surface accepted by
  `include/macrodr/dsl/lexer_untyped.h`, not a speculative future DSL.
- Operator expressions such as `a + b` are not documented here because they are
  not represented as first-class grammar productions in the current spec set.
- The typed meaning of vectors, tuples, references, and function overloads is
  still defined by the implementation and by `type_rules.md` / `run_rules.md`.
