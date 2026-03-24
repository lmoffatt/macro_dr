# ModelPtr Refactor Census

## Goal

Remove the DSL special case:

- `const IModel<...>&` -> stored as `std::unique_ptr<IModel<...>>`

and normalize model-like polymorphic objects around `ModelPtr`-shaped APIs.

Current intent:

- full internal sweep
- not only DSL registrations
- also `cmd/` and `src/core/` signatures where they currently expose
  `const interface::IModel<var::Parameters_values>&`

## Summary

This refactor touches four layers:

1. DSL storage/type adaptation
2. DSL registration surface
3. `cmd/` public APIs
4. `src/core/` implementations

There is also documentation fallout.

The main architectural effect is:

- the DSL should traffic in explicit owning model handles (`ModelPtr`)
- borrowed polymorphic access should be derived from those handles in wrappers or
  implementation code
- no DSL type-storage special case should be needed for `IModel`

## A. Must-change: DSL storage machinery

### 1. Remove the special storage rule

File:
- `include/macrodr/dsl/lexer_typed.h`

Current special case:
- `function_argument_storage<const macrodr::interface::IModel<ParamValues...>&>`

This is the special rule to remove.

### 2. Re-check generic adaptation logic

Files:
- `include/macrodr/dsl/lexer_typed.h`
- `include/macrodr/dsl/grammar_typed.h`

Reason:
- once `const IModel&` disappears from DSL-facing signatures, the generic rules
  should be enough:
  - `T`
  - `const T&` / `T&`
  - `std::unique_ptr<T>`

Important runtime sites:
- `function_compiler::adapt_argument(...)`
- `typed_vector_construction::adapt(...)`
- `typed_set_construction::adapt(...)`
- `typed_tuple_construction::adapt(...)`

These may not need logic changes, but they are part of the refactor surface and
must be revalidated.

### 3. Legacy helper mismatch

File:
- `legacy/CLI_macro_dr.h`

Problem:
- this helper still builds `field_compiler<Args>`
- the newer helper in `include/macrodr/dsl/function_builder.h` uses
  `field_compiler<detail::function_argument_storage_t<Args>>`

If `ModelPtr` reference forms are used through the old helper, behavior will
diverge from the new path.

Decision implied by this refactor:
- either update `legacy/CLI_macro_dr.h` to match the modern helper
- or stop using it for any DSL-facing registrations

## B. Must-change: DSL registration surface

File:
- `src/cli/command_manager.cpp`

This is the canonical registration assembly point.

### 1. Simulation registrations

Current registrations using `const IModel&`:
- `simulate` overload from values
- `simulate` overload from transformed values
- `simulate` returning `vector<Simulated_Recording<...>>` from values
- `simulate` returning `vector<Simulated_Recording<...>>` from transformed values
- `simulate_with_sub_intervals` overload from values
- `simulate_with_sub_intervals` overload from transformed values

These are currently in the `make_simulations_compiler()` block.

### 2. Likelihood / diagnostics registrations

Current registrations using `const IModel&`:
- `calc_likelihood`
- `calc_dlikelihood`
- `calc_diff_likelihood`
- `calc_likelihood_predictions`
- `calc_likelihood_diagnostic`
- `calc_dlikelihood_predictions`

including simulation/data overloads where applicable.

These are currently the registrations around:
- lines ~395 to ~509 in the current file

### 3. Registrations already using `ModelPtr`

These already align conceptually and should be kept consistent:
- `load_parameters`
- `patch_model` overloads
- `path_state` overload from model + values
- `build_likelihood_function`

These should be reviewed for signature consistency once the new rule is chosen:
- by-value `ModelPtr`
- `const ModelPtr&`
- or mixed ownership policy

## C. Must-change: public `cmd/` APIs

### 1. `include/macrodr/cmd/simulate.h`

Current `const IModel&` family:
- `run_simulations(...)` x2
- `run_simulations_with_sub_intervals(...)` x2
- `run_n_simulations(...)` x2
- `run_n_simulations_with_sub_intervals(...)`
- transformed `run_simulations_with_sub_intervals(...)` with `n_simulations`

These are all candidates to normalize to `ModelPtr`-shaped signatures.

### 2. `include/macrodr/cmd/likelihood.h`

Current `const IModel&` family:
- `calculate_likelihood(...)`
- `calculate_simulation_likelihood(...)`
- `calculate_dlikelihood(...)`
- `calculate_simulation_dlikelihood(...)`
- `calculate_diff_likelihood(...)`
- `calculate_simulation_diff_likelihood(...)`
- `calculate_likelihood_predictions(...)`
- `calculate_simulation_likelihood_predictions(...)`
- `calculate_likelihood_diagnostics(...)`
- `calculate_simulation_likelihood_diagnostics(...)`
- `calculate_dlikelihood_predictions(...)`
- `calculate_simulation_dlikelihood_predictions(...)`
- `calculate_simulation_sub_dlikelihood_predictions(...)`

Also review:
- helper wrappers in the same file that forward to the base functions

### 3. Files that already define the owning-handle vocabulary

These are not the main refactor targets, but they anchor the intended type:
- `include/macrodr/cmd/patch_model.h`
- `include/macrodr/cmd/load_model.h`

## D. Must-change: `src/core/` implementations

### 1. `src/core/simulate.cpp`

Current implementation family using `const IModel&`:
- `run_simulations(...)`
- `run_simulations_with_sub_intervals(...)`
- `run_n_simulations(...)`
- `run_n_simulations_with_sub_intervals(...)`

These must follow whatever new `cmd/` signature policy is chosen.

### 2. `src/core/likelihood.cpp`

Current implementation family using `const IModel&`:
- `calculate_likelihood(...)`
- `calculate_dlikelihood(...)`
- `calculate_diff_likelihood(...)`
- `calculate_likelihood_predictions(...)`
- `calculate_likelihood_diagnostics(...)`
- `calculate_dlikelihood_predictions(...)`

These must follow the same new `cmd/` policy.

## E. Documentation that must change

### 1. DSL project docs

Files:
- `projects/macrodr-dsl/DSL_RECONSTRUCTION.md`
- `projects/macrodr-dsl/DSL_CHANGES.md`

Reason:
- they currently describe the `const IModel& -> unique_ptr<IModel>` special case

### 2. Existing spec docs

Files:
- `docs/spec/dsl/dsl_structure_summary.md`
- `docs/spec/dsl/type_rules.md`

Reason:
- both currently document the `IModel` special storage case explicitly

### 3. Any future command inventory doc

If a current DSL function catalog is generated next, it must reflect the new
`ModelPtr`-shaped signatures rather than the old borrowed-polymorphic ones.

## F. Likely compatibility work

These are not necessarily separate files, but they are part of the refactor:

- wrapper functions may be needed temporarily if some call sites still prefer
  `const IModel&`
- overload resolution may change where `ModelPtr`, `const ModelPtr&`, and
  value-oriented parameters coexist
- tests or scripts relying on current implicit model unwrapping behavior will
  need review

## G. What does not appear to require direct change

These areas already fit the new direction or are orthogonal:

- `load_model(...)` returning an owning model handle
- `typed_literal<unique_ptr<T>>` support in `include/macrodr/dsl/grammar_typed.h`
- non-model scalar/reference storage rules
- JSON environment type registration, except for any model-related assumptions
  that may later be added

## H. Main design decision still required

Even under "full internal sweep", one detail still needs to be fixed explicitly:

Which model-handle shape should be standard in the new signatures?

Candidate shapes:
- `ModelPtr`
- `const ModelPtr&`
- `ModelPtr&&` for consuming APIs

This choice affects:
- how often ownership moves
- whether a loaded model can be reused across multiple DSL calls
- how much wrapper code is needed
