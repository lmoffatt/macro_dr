# DSL Dead Code Cleanup

This note records the main dead-code and near-dead-code candidates in the
current DSL implementation, with a conservative cleanup order.

It focuses on:

- definitely dead code
- stale legacy scaffolding
- weakly-live infrastructure that is still compiled but no longer part of the
  main execution path

## Current active path

The modern DSL function path is:

1. untyped parse tree
2. `function_compiler<F, Args...>`
3. `field_compiler<storage_t<Arg>>`
4. `typed_function_evaluation<..., storage_t<Args>...>`
5. runtime argument adaptation through `adapt_argument(...)`

Code anchors:

- [grammar_untyped.h](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_untyped.h)
- [lexer_typed.h](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h)
- [grammar_typed.h](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h)

The important implication is:

- the modern path works directly from `untyped_argument_list`
- it does not use the old typed argument-list machinery at runtime

## A. Definitely dead: commented legacy block in `grammar_typed.h`

There is a large commented block at:

- [grammar_typed.h:1243](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1243)

This includes old definitions for:

- `base_typed_argument_list`
- old `typed_assignment`
- duplicated `typed_identifier`
- old variadic `typed_argument_list<Lexer, Compiler, Ts...>`
- old `typed_function_evaluation`

These are not compiled at all. They are pure commented historical residue.

### Recommendation

Safe to delete immediately.

Reason:

- they are already inactive
- they duplicate concepts that now exist elsewhere
- they make `grammar_typed.h` harder to reason about

## B. Weakly-live: `typed_argument_list<Lexer, Compiler>`

The active generic `typed_argument_list<Lexer, Compiler>` is defined at:

- [grammar_typed.h:261](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L261)

It is still constructed by:

- [untyped_argument_list::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_untyped.h#L321)
- [untyped_argument_list::compile_expression](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_untyped.h#L311)

and populated through:

- [base_typed_assigment::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1225)
- [base_typed_expression::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1233)

### But the modern path does not consume it

The modern function system does not use `typed_argument_list<Lexer, Compiler>`
as a runtime structure.

Instead it consumes:

- `untyped_argument_list`

directly in:

- [base_function_compiler::compile_function_evaluation](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L215)
- [function_compiler::compile_function_evaluation_impl](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L1070)

### Assessment

`typed_argument_list<Lexer, Compiler>` is not obviously dead, but it is no
longer part of the main typed-function pipeline.

It is better described as:

- leftover generic argument-list infrastructure

### Recommendation

Do not delete immediately.

First verify whether any remaining DSL feature still depends on:

- `untyped_argument_list::compile_expression`
- `typed_argument_list<Lexer, Compiler>::arg_vector()`
- `typed_argument_list<Lexer, Compiler>::arg_map()`

If no meaningful consumer remains, this whole path can likely be removed.

## C. Likely dead: `typed_argument_typed_list`

Defined at:

- [grammar_typed.h:319](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L319)

I found no active references to:

- `typed_argument_typed_list<...>`

It looks like an abandoned typed tuple-of-arguments container that predates the
current `function_compiler` tuple path.

### Recommendation

Strong candidate for removal after one quick build verification.

## D. Legacy compile hooks tied to weakly-live argument-list machinery

These methods exist mainly to feed `typed_argument_list<Lexer, Compiler>`:

- [base_typed_assigment::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1225)
- [base_typed_expression::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1233)

If the generic typed argument-list path is removed, these virtual hooks can
likely disappear too.

### Recommendation

Treat them as second-wave cleanup, not first-wave cleanup.

They are still wired into:

- [base_typed_statement](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L157)
- [untyped_argument_list::compile_argument_list](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_untyped.h#L321)

## E. What is clearly live and should not be touched in this cleanup

Keep:

- `field_compiler`
- `element_compiler`
- `function_compiler`
- modern `typed_function_evaluation`
- `typed_vector_construction`
- `typed_tuple_construction`
- `typed_set_construction`
- `typed_literal`
- `typed_identifier`
- `typed_identifier_ref_const`
- `typed_identifier_ref`

These are part of the current DSL runtime path.

## Cleanup plan

### Phase 1: no-risk cleanup

Remove the commented legacy block at:

- [grammar_typed.h:1243](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L1243)

### Phase 2: weakly-live verification

Check whether any current behavior actually depends on:

- `typed_argument_list<Lexer, Compiler>`
- `untyped_argument_list::compile_expression`
- `base_typed_statement::compile_argument_list`

If not, remove them together.

### Phase 3: remove likely-dead typed container

Remove:

- `typed_argument_typed_list`

if the build confirms no hidden dependency remains.

## Main conclusion

The DSL has one obvious dead region and one likely-obsolete subsystem.

### Obvious dead region

- the commented legacy block at the bottom of `grammar_typed.h`

### Likely-obsolete subsystem

- generic `typed_argument_list` compilation infrastructure

The modern DSL no longer needs that subsystem for typed function calls, so it
is the next serious cleanup target after the dead commented block.
