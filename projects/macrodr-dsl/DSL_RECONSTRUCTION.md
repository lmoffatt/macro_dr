# MacroDR DSL Reconstruction

This document reconstructs the current MacroDR DSL from the existing
documentation under `docs/` and the current implementation entry points.

It is meant to be the working description of the DSL as it exists today,
without pretending that every older spec file is fully current.

## Sources used

### Documentation

- `docs/cli.md`
- `docs/spec/README.md`
- `docs/spec/dsl/dsl_structure_summary.md`
- `docs/spec/dsl/grammar.bnf.md`
- `docs/spec/dsl/type_rules.md`
- `docs/spec/dsl/run_rules.md`
- `docs/spec/dsl/functions.md`
- `docs/spec/cli/command_syntax.spec.md`

### Code anchors

- `include/macrodr/dsl/grammar_untyped.h`
- `include/macrodr/dsl/grammar_typed.h`
- `include/macrodr/dsl/lexer_typed.h`
- `include/macrodr/cli/command_manager.h`
- `src/cli/command_manager.cpp`
- `src/cli/main.cpp`
- `src/cli/app/main_flow.cpp`
- `src/cli/app/execute_program.cpp`

## Current DSL model

### 1. Process boundary vs DSL boundary

MacroDR has two distinct layers:

- the process-level CLI
- the embedded DSL executed by that CLI

The CLI is intentionally thin. The actual entrypoint is:

- `src/cli/main.cpp` -> `macrodr::cli::app::main_flow(...)`

`main_flow(...)` assembles the script, creates the compiler with
`macrodr::cli::make_compiler_new()`, and then passes control to
`execute_program(...)`.

That means the current canonical DSL registry surface is not `main.cpp`.
It is:

- `src/cli/command_manager.cpp`
- `include/macrodr/cli/command_manager.h`

### 2. Execution pipeline

The current pipeline is:

```text
CLI args / script files / --eval
  -> assembled script text
  -> extract_program(script)
  -> untyped_program
  -> compile_program(env, parsed)
  -> typed_program
  -> typed_program.run(env)
```

The relevant code anchors are:

- parsing and untyped AST: `include/macrodr/dsl/grammar_untyped.h`
- typed runtime objects and environment: `include/macrodr/dsl/grammar_typed.h`
- function registry and argument compilation: `include/macrodr/dsl/lexer_typed.h`
- process execution: `src/cli/app/execute_program.cpp`

### 3. Ownership of compile-time and runtime state

The current design separates:

- `Compiler`
  - registered functions
  - registered JSON-loadable types
  - overload sets for functions
- `Environment`
  - identifier compilers
  - runtime variable values
  - parameter-schema registry

This is important because some older DSL docs describe the compiler as if it
also owned identifier bindings. In the current code, identifier lookup during
compilation and variable lookup during execution live in `Environment`.

The key structures are:

- `Compiler` in `include/macrodr/dsl/lexer_typed.h`
  - `m_func`
  - `m_type_registry`
- `Environment` in `include/macrodr/dsl/grammar_typed.h`
  - `m_id`
  - `m_var`
  - parameter schema map

### 4. Untyped syntax model

The untyped layer currently supports:

- assignments
- identifiers
- numeric literals
- string literals
- function calls
- argument lists
- vector literals: `[e1, e2, ...]`
- tuple literals: `{e1, e2, ...}`

This is represented by:

- `untyped_program`
- `untyped_statement`
- `untyped_expression`
- `untyped_assignment`
- `untyped_identifier`
- `untyped_function_evaluation`
- `untyped_vector_construction`
- `untyped_tuple_construction`

The grammar files in:

- `docs/spec/dsl/grammar.bnf`
- `docs/spec/dsl/grammar.bnf.md`

have now been updated to track the current untyped surface more closely.
They should be treated as the working grammar description for the lexer-facing
DSL, while the implementation remains the final authority if a mismatch is
found.

### 5. Compilation model

Compilation is single-pass.

`compile_program(env, parsed)` walks each untyped statement in order, compiles
it immediately, and appends it to a `typed_program`.

Because identifier compilers are populated during compilation, the DSL still
has a defined-before-use constraint. In practice:

- later definitions are not visible earlier in the script
- forward declarations are not part of the current language
- mutually recursive or backpatched definitions are not supported

This part of the older spec is still correct in spirit.

### 6. Function registration model

The current registry model is:

- functions are registered into `dsl::Compiler`
- registration is done through `push_function(...)`
- registration helpers are `to_typed_function(...)` and related builders
- functions are now overloadable

The current `Compiler` stores:

- `std::map<Identifier, std::vector<std::unique_ptr<base_function_compiler>>>`

So one DSL name can map to multiple overloads.

The compiler also supports:

- `merge(const Compiler&)`
- `merge(Compiler&&)`
- `get_functions(...)`
- automatic type registration for JSON environment loading

This is more capable than the older docs imply.

### 7. Registry assembly in the current codebase

The canonical registry builder is `make_compiler_new()` in
`src/cli/command_manager.cpp`.

That builder currently merges:

- `macrodr::cmd::make_cli_meta_compiler()`
- `make_utilities_compiler()`
- `make_io_compiler()`
- `make_experiment_compiler()`
- `make_model_compiler()`
- `make_simulations_compiler()`
- `make_likelihood_compiler()`
- `make_dts_compiler()`

This means the DSL surface is currently a hybrid:

- some commands are already exposed through `include/macrodr/cmd/*`
- many registry builders still come from legacy headers

So the command surface is in transition, but `make_compiler_new()` is the
current authoritative assembly point.

### 8. Typed argument model

Function arguments are compiled through `field_compiler`.

The important rule is that the DSL distinguishes:

- exposed argument type
- internal storage type

The storage type is derived through
`detail::function_argument_storage_t<Arg>`.

Important cases:

- `T` -> stored as decayed `T`
- `const T&` / `T&` -> stored as `std::reference_wrapper<...>`
- `const IModel<...>&` -> stored as owning `std::unique_ptr<IModel<...>>`

This same pattern extends to composite values.

### 9. Composite values

The DSL currently supports:

- vectors
- sets through vector-like syntax when the target type is a `std::set`
- tuples

These are compiled through:

- `vector_compiler`
- `set_compiler`
- `tuple_compiler`

and represented in the typed layer through:

- `typed_vector_construction`
- `typed_set_construction`
- `typed_tuple_construction`

The important semantic rule is:

- the target type is context-sensitive
- composite literals are compiled against the type expected by the receiving
  parameter or assignment context

### 10. Reference semantics

Reference-like arguments are intentionally restricted.

For parameters compiled as:

- `std::reference_wrapper<T>`
- `std::reference_wrapper<const T>`

the DSL only accepts identifier expressions in those positions.

It rejects:

- raw literals
- inline function calls
- anonymous composite expressions

This is a real current rule in `field_compiler`, not just a draft design note.

### 11. Runtime execution

After compilation, `typed_program.run(env)` executes statements in order.

Assignments:

- evaluate the right-hand side
- insert the resulting typed expression into `env.m_var`

Identifiers at runtime:

- resolve from the environment
- re-run the bound typed expression as needed

This means the environment is the effective runtime store for script state.

### 12. Environment persistence

The current CLI supports environment save/load around the DSL.

The relevant features are described in `docs/cli.md` and implemented in
`src/cli/app/execute_program.cpp` together with
`include/macrodr/io/json/environment_io.h`.

Important consequences:

- the compiler has a type registry for JSON loading
- environment snapshots are now part of the runtime model
- only a subset of types is currently serializable in a useful way

This persistence layer is now part of the practical DSL story and should be
considered when redesigning diagnostics or script semantics.

## Current limitations

### 1. Single-pass visibility

The DSL still requires definitions before use.

### 2. Weak source-level diagnostics

The current compile/run path does not yet preserve enough script structure or
location context to make errors consistently easy to interpret in large active
workflows such as `projects/eLife_2025/ops/local/figure_1.macroir`.

### 3. Registry migration is incomplete

The public command surface is partly moving toward `include/macrodr/cmd/`,
but the actual registry still depends heavily on legacy builders.

### 4. Spec documents are unevenly current

Some spec files remain useful as conceptual drafts, but they do not all match
the current implementation.

## What should currently be treated as source of truth

For the DSL project, the current source of truth order should be:

1. implementation anchors in `include/macrodr/dsl/`, `src/cli/`, and
   `src/cli/app/`
2. command assembly in `src/cli/command_manager.cpp`
3. `docs/cli.md`
4. the better parts of `docs/spec/dsl/*.md`
5. older generated summaries such as `docs/spec/dsl/functions.md`

## Immediate next DSL work implied by this reconstruction

1. build a current command/function inventory from `make_compiler_new()`
   instead of relying on the stale generated catalog
2. add source-location and evaluation-context concepts so diagnostics can name
   the failing DSL statement more precisely
3. decide which parts of the old spec should be rewritten as current docs and
   which should be archived as draft/generated material
