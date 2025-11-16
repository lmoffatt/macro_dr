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

---

## 5. Composite Values: Vectors and Tuples

### 5.1 Untyped Nodes

Two new untyped AST nodes represent composite literals:

- `untyped_vector_construction<Lexer,Compiler>` – parses `[e1, e2, …]`
- `untyped_tuple_construction<Lexer,Compiler>` – parses `{e1, e2, …}`

They can appear both as:
- Top-level statements (for REPL use), and
- Sub-expressions inside function calls and assignments.

### 5.2 Storage vs Exposed Types

The typed layer distinguishes:

- **Exposed type** – what the DSL user sees:
  - Vectors: `std::vector<T>`
  - Tuples: `std::tuple<Ts...>`
- **Storage type** – what the compiler actually stores for each element, via

```cpp
template<class Arg>
using function_argument_storage_t =
    detail::function_argument_storage_t<Arg>;
````

Rules (same as for function arguments):

* `T`, `const T` → stored as `T` (decayed)
* `T&`, `const T&` → stored as `std::reference_wrapper<T>` / `std::reference_wrapper<const T>`
* `const IModel<...>&` → stored as `std::unique_ptr<IModel<...>>`

This guarantees that vectors/tuples behave consistently with function parameters, including ownership and reference semantics.

### 5.3 element_compiler

`element_compiler<Lexer,Compiler,T>` is the unnamed analogue of `field_compiler`:

* Compiles **one positional element** of a vector/tuple to
  `typed_expression<Lexer,Compiler,T>` (or to the storage type).

* Handles:

  * identifiers (`untyped_identifier`)
  * literals (`untyped_literal`, when `literal_decodable<T>`)
  * function calls (`untyped_function_evaluation`)
  * nested `[ ... ]` and `{ ... }` when `T` is a `std::vector<...>` or `std::tuple<...>`.

* For reference elements (`T = std::reference_wrapper<U>`), only identifiers
  are accepted; literals and function calls are rejected by design.

### 5.4 vector_compiler and typed_vector_construction

For vectors:

```cpp
template<class Lexer, class Compiler, class T>
class vector_compiler {
    using storage_t = function_argument_storage_t<T>;
    // ...
    Maybe_unique<typed_expression<Lexer,Compiler,std::vector<T>>>
    compile_vector_construction(Environment<Lexer,Compiler> const&,
                                const untyped_argument_list<Lexer,Compiler>&);
};
```

* `vector_compiler` compiles each element using `element_compiler<Lexer,Compiler,storage_t>`.
* The result is a `typed_vector_construction<Lexer,Compiler,T>` that:

  * stores a `std::vector< unique_ptr<typed_expression<storage_t>> >`
  * on `run(env)`:

    * evaluates each element to `storage_t`
    * uses the same `adapt<T>(storage)` logic as function calls
    * returns `std::vector<T>`.

Thus vector literals `[e1, e2, …]` are first-class expressions with type `std::vector<T>`.

### 5.5 tuple_compiler and typed_tuple_construction

For tuples:

```cpp
template<class Lexer, class Compiler, class... Ts>
class tuple_compiler {
    using storage_t<U> = function_argument_storage_t<U>;
    // ...
    Maybe_unique<typed_expression<Lexer,Compiler,std::tuple<Ts...>>>
    compile_tuple_construction(Environment<Lexer,Compiler> const&,
                               const untyped_tuple_construction<Lexer,Compiler>&);
};
```

* `tuple_compiler` uses a `std::tuple` of `element_compiler<Lexer,Compiler,storage_t<Ts>>...`.

* It builds a `typed_tuple_construction<Lexer,Compiler,Ts...>` with:

  * internal storage: `std::tuple< unique_ptr<typed_expression<storage_t<Ts>>>... >`
  * exposed type: `std::tuple<Ts...>`.

* On `run(env)`, it:

  * evaluates each element into `storage_t<Tk>`
  * adapts them via the same `adapt<Tk>(storage)` logic
  * returns `std::tuple<Ts...>`.

### 5.6 Interaction with Single-Pass Compilation

* Vectors and tuples **do not change** the single-pass model:

  * They are compiled in place as ordinary expressions.
  * Elements are resolved using the same environment/Compiler lookup as scalar expressions.
* The same **“defined-before-use”** rule applies to identifiers inside `[ ... ]` and `{ ... ]`.
* Elements that compile to references (`std::reference_wrapper<...>`) preserve the current
  constraint: they must be identifiers bound in the environment at execution time.

````


*Keep this document for reference when evolving the DSL's compilation strategy.*

