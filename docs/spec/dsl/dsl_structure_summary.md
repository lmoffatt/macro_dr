# MacroDR DSL Architecture Summary

This document outlines the core components of the MacroDR embedded DSL and highlights how its single-pass compilation model impacts support for mutually recursive definitions, as well as the current symbol-registration mechanism.

---

## 1. Core Components

1. **Compiler (**``** )**

   - Maintains semantic context and symbol tables:
     - `m_func`: Registered functions (`std::map<Identifier, std::unique_ptr<base_function_compiler>>`)
     - `m_id`: Registered identifiers (`std::map<Identifier, std::unique_ptr<base_Identifier_compiler>>`)
   - Populated during DSL setup (`push_function`) and mutated per-statement via `` for assignments.

2. **Untyped Program (**``** )**

   - Parsed representation of raw text:
     - `std::vector<std::unique_ptr<untyped_statement>> m_statements`
   - Captures syntax without semantics.

3. **Typed Program (**``** )**

   - Semantically-checked AST:
     - Sequence of `typed_statement` nodes
   - Constructed via ``:
     - Iterates `untyped` statements, invoking `compile_statement(cm)`
     - Appends each typed node with ``
   - Exposes `run()` to execute against an environment.

4. **Environment (**``** )**

   - Runtime state mapping identifiers to values:
     - `std::map<Identifier, std::unique_ptr<base_typed_expression>> m_var`
   - Populated by `run_statement` calls during program execution.

---

## 2. Execution Pipeline

```
Raw Text → Lexer → Untyped Program (`p`)
           ↓
      compile_program(cm, p) → Typed Program (`c`)
                                  ↓
                         c.run() → Environment (`env`)
```

1. **Parsing**: Lexer tokenizes text; parser builds `untyped_program`.
2. **Compilation**: Single-pass walk over `p`:
   - `` registers new identifiers immediately in `cm`.
   - `` appends compiled statements to `c`.
3. **Execution**: `c.run()` executes typed statements, updating `env`.

---

## 3. Symbol Registration Mechanism

- **Compile-Time Registration**: `Compiler::push_back` populates `cm.m_id` with `base_Identifier_compiler` pointers for each symbol as soon as it's compiled, enabling name resolution during `compile_program`.
- **Typed Program Table**: `typed_program` defines an `insert(id, expr)` method to register `base_typed_expression` pointers in its own `m_identifier_table`, intended for compile-time lookup of compiled expressions.
- **Unused Table**: In practice, `compile_program` never calls `typed_program::insert`; instead, **symbol resolution is fully handled by `Compiler`**, making the typed-program-level table vestigial.
- **Rationale for Omission**: Since actual expression instances (`base_typed_expression`) are materialized during the **execution** (`run()`) phase within the `Environment`, maintaining a compile-time table of these objects is unnecessary.

------

## 4. Single-Pass Limitations

- **Immediate Binding**: Symbols must be defined before use; `compile_statement` registers identifiers on the fly.
- **Mutual Recursion**: Without forward declarations or a preliminary symbol-collection pass, references to later-defined identifiers cannot be resolved.

---

### Possible Remedies

- Introduce explicit **forward-declaration** syntax to pre-register names.
- Add a **first-pass** to collect all definitions before type-checking.
- Implement a **backpatching** stage to resolve unresolved identifiers post-compilation.

---

*Keep this document for reference when evolving the DSL's compilation strategy.*

