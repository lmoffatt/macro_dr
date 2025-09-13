
---

# MacroIR DSL Grammar Specification (BNF)
# Status: DRAFT
# Last Updated: 2025-07-17
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
# Each line is parsed independently into one of the statement types.

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
# General expression node (used in RHS)

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
```

---

### ✅ Notes

* This grammar directly corresponds to the structures defined in:

  * `untyped_program`, `untyped_statement`, `untyped_function_evaluation`, `untyped_assignment`
* Function arguments are structurally parsed and support nesting (e.g., `simulate(model=init("lin"))`)
* Comments (`#`) are parsed and ignored outside of string literals

---

