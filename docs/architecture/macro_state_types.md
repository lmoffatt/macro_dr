# Macro State Types

This note explains why `Macro_State`, `dMacro_State`, and `ddMacro_State` all exist and how they are used.
All three live in `legacy/qmodel.h`.

## Overview

The MacroDr likelihood pipeline needs three distinct notions of "state":

- value state (no derivatives),
- a vector space of first-derivative components, and
- a derivative wrapper around the whole state so the generic derivative machinery can reason about it.

Those are not interchangeable types in the current derivative framework, so each gets its own struct.

## Macro_State (value state)

`Macro_State<Vars...>` is a plain value container:

- Type: `Vector_Space<logL, Patch_State, Vars...>`
- Purpose: hold the accumulated log-likelihood, patch recursion state, and optional extras.
- Used in prediction/diagnostic paths and value-only log likelihood computations.
- The optional `Evolution_of<T>` captures per-time-step outputs when requested.

## dMacro_State (first-derivative state)

`dMacro_State<Vars...>` is a vector space whose components are derivatives:

- Type:
  `Vector_Space<Derivative<logL, Parameters_transformed>,
               Derivative<Patch_State, Parameters_transformed>,
               Vars...>`
- Purpose: carry first derivatives alongside the same structural layout as `Macro_State`.
- Special constructor seeds a shared `dx` into all derivative components. This is required to
  satisfy the `log_Likelihood` invariant that derivative-mode states must have `dx`.
- Keeps update logic symmetric with `Macro_State` while preserving derivative types per component.

## ddMacro_State (derivative-of-state wrapper)

`ddMacro_State<Vars...>` is not "Macro_State with derivative members".
It is a derivative wrapper around the *whole* vector space:

- Type:
  `Derivative<Vector_Space<logL, Patch_State, Vars...>, Parameters_transformed>`
- Purpose: let generic derivative operations (Transfer_Op_to, dx propagation, etc.)
  treat the entire state as a derivative object, not just its members.
- Specializations in `var` declare:
  - `transformation_type<ddMacro_State<...>>` as `Derivative_Op<Parameters_transformed>`
  - `is_derivative<ddMacro_State<...>> = true`
  - `dx_of_dfdx` for `ddMacro_State`

This is the key reason it cannot be trivially replaced by `Macro_State` or `dMacro_State`:
the derivative framework recognizes `Derivative<T, X>` as a derivative type, not arbitrary structs.

## Why keep all three

The three types encode different semantic roles:

- `Macro_State`: value-only accumulation and evolution tracking.
- `dMacro_State`: per-component derivatives (first-order), easy to update like a vector space.
- `ddMacro_State`: a derivative wrapper so the derivative infrastructure can propagate `dx`
  and transformation types across the whole state.

Merging them would either:

- lose the derivative-wrapper semantics needed by the generic machinery, or
- require new trait specializations and deeper refactors across the likelihood and CLI layers.

## Practical usage map

- `logLikelihood` / predictions: returns `Macro_State_*` variants.
- `dlogLikelihood` / derivative predictions: returns `dMacro_State_*` variants.
- Update overloads are specialized for each of the three types to keep recursion, accumulation,
  and `Evolution` handling consistent.

## If you revisit the design

Any attempt to unify these types should preserve:

- `dx` seeding behavior in derivative mode,
- `is_derivative` and `transformation_type` for the derivative-wrapper case, and
- the existing overload resolution in the update path.

Otherwise the derivative invariants enforced in `log_Likelihood` will break.
