# MacroDR DSL Documentation Changes

This document records what changed between the older DSL docs under `docs/`
and the reconstructed current DSL model in
`projects/macrodr-dsl/DSL_RECONSTRUCTION.md`.

## Reviewed documents

- `docs/cli.md`
- `docs/spec/README.md`
- `docs/spec/cli/command_syntax.spec.md`
- `docs/spec/dsl/dsl_structure_summary.md`
- `docs/spec/dsl/functions.md`
- `docs/spec/dsl/grammar.bnf.md`
- `docs/spec/dsl/type_rules.md`
- `docs/spec/dsl/run_rules.md`

## What still holds

These older docs remain broadly correct:

- the DSL has an untyped-to-typed compilation pipeline
- compilation is still effectively single-pass
- vectors and tuples are real parts of the language
- storage-vs-exposed type reasoning is a real part of the implementation
- runtime execution is statement-ordered through an environment

So the older spec is not useless. It is just uneven.

## What changed or needed correction

### 1. `main.cpp` is no longer the command-registration center

Older docs and generated notes still talk as if `main.cpp` or a generic
`make_compiler()` were the natural center of the DSL surface.

Current reality:

- `src/cli/main.cpp` is thin
- `src/cli/app/main_flow.cpp` owns the process flow
- `src/cli/command_manager.cpp` owns the canonical DSL compiler assembly via
  `make_compiler_new()`

This is a major architectural correction.

### 2. `Compiler` no longer matches the old summary

`docs/spec/dsl/dsl_structure_summary.md` describes a `Compiler` with:

- `m_func`
- `m_id`

That is not the current structure.

Current reality:

- `Compiler` owns function overloads and a type registry
- `Environment` owns identifier compilers and runtime variables

So identifier ownership moved conceptually out of the compiler description.

### 3. `typed_program` is simpler than the old summary says

The old summary claims `typed_program` has an identifier table and an `insert`
path that became vestigial.

Current reality:

- `typed_program` is just a sequence of typed statements
- it exposes `run(env)`
- there is no typed-program identifier table in the current code

So the old document is describing either an older design or a partially
remembered one, not the current implementation.

### 4. Function registration is now overload-aware

The old docs usually speak as if a function name maps to one compiler entry.

Current reality:

- `Compiler::m_func` stores a vector of function compilers per identifier
- `get_functions(...)` exists
- `merge(...)` preserves overload sets

This matters for reconstructing the language surface and for future
diagnostics, because ambiguity and overload selection are now real concerns.

### 5. Type registration is now part of the DSL runtime ecosystem

The older DSL summaries understate the role of JSON/environment persistence.

Current reality:

- `Compiler` owns `m_type_registry`
- built-in types are auto-registered
- JSON environment loading uses those registrations

So the DSL is no longer just "parse, compile, run". It also participates in
typed environment persistence.

### 6. `docs/spec/dsl/functions.md` is stale

This file still has value as a historical catalog, but it is not an accurate
current source of truth.

Problems:

- it says it was generated from `make_compiler()`
- it references `main.cpp` as a source anchor
- it foregrounds the older legacy registry surface
- it does not reflect the full current command-manager assembly

It should not be used as the canonical function inventory.

### 7. `docs/spec/README.md` is partly draft/generated

It still gives a useful overview of the spec area, but:

- it references files and schema areas that are not the current focus
- it ends with assistant-style prompt text
- it should be treated as a draft summary, not a fully curated canonical doc

### 8. `docs/spec/cli/command_syntax.spec.md` does not match the current CLI

This file describes a command-oriented shell interface like:

- `macroir run ...`
- `macroir eval ...`
- `macroir compile ...`
- `macroir describe ...`

Current reality in `docs/cli.md` and the code is:

- `macro_dr [options] <script...> [-e "<dsl>"]`
- `--check-syntax`
- `--env-save`, `--env-load`
- script assembly and run workspace persistence

So `command_syntax.spec.md` is a speculative or obsolete CLI direction, not
the current interface.

### 9. Grammar doc was upgraded into a working spec

The grammar files:

- `docs/spec/dsl/grammar.bnf`
- `docs/spec/dsl/grammar.bnf.md`

were updated as part of this reconstruction.

The main corrections were:

- add comment statements
- document signed and exponential numeric literals
- document both single-quoted and double-quoted strings
- make newline separation explicit
- fix the old tuple-brace typo

So grammar is no longer just a draft note. It is now a working untyped syntax
spec, still subordinate to the lexer implementation when there is a conflict.

### 10. Type and runtime rules are mostly conceptual reconstructions

`docs/spec/dsl/type_rules.md` and `docs/spec/dsl/run_rules.md` contain useful
high-level descriptions, but they are not direct executable truth.

Signals of draft/generated status:

- pseudo-C++ placeholders
- incomplete context-dependent types
- assistant-style trailing prompts

They remain useful explanatory documents, but they should not outrank the
current implementation files.

## Updated documentation rule for the DSL project

For current DSL work:

- use `projects/macrodr-dsl/DSL_RECONSTRUCTION.md` as the working current
  description
- use this file to understand what changed and what is stale
- use implementation files as the final authority when a doc and code disagree

## Recommended follow-up

1. replace the stale generated function catalog with a current command-surface
   inventory derived from `make_compiler_new()`
2. rewrite the most valuable parts of `docs/spec/dsl/` into a smaller curated
   canonical spec set
3. archive or clearly mark speculative/obsolete DSL and CLI spec files
