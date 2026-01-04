# Moment Statistics

This note documents the `Moment_statistics` utilities extracted from `legacy/mcmc.h` into
`legacy/moment_statistics.h`.

## Purpose

Provide a small, reusable container for tracking running mean/variance with a count, independent
of MCMC logic. This is used for log-likelihood and evidence statistics but is also useful elsewhere.

## Types

- `count`: `var::Constant<count, std::size_t>` with default zero.
- `mean<Va>`: `var::Var<mean<Va>, Va>`.
- `variance<Va>`: `var::Var<variance<Va>, Va>` with explicit constructors.
- `Moment_statistics<Va>`: `var::Var<..., Vector_Space<mean<Va>, variance<Va>, count>>`.
- `Moment_statistics<var::Vector_Space<Ids...>>`: componentwise stats, stored as
  `var::Vector_Space<Moment_statistics<Ids>...>`. `get<Id>(...)` returns the corresponding
  `Moment_statistics<Id>` (via a custom `get`), and `operator[]` on the scalar form supports
  `Var`-style access.

`Va` is expected to behave like a `var::Var`-compatible scalar:
it is accessed via `x()` in updates, and must support basic arithmetic.

## Behavior

`Moment_statistics` supports:

- Default construction (mean/variance zero, count zero).
- Construction from a single value (mean = x, variance = 0, count = 1).
- Construction from explicit `(mean, variance, count)`.
- Merge of two statistics with `operator&` / `operator&=`.
- Online update with a single sample via `operator&=(const Va&)`.
- Reset to zeros.
- Linear operations: `operator+` and scaling by a scalar `operator*`.

The variance update uses `sqr_X` from `legacy/matrix.h`.

For `Moment_statistics<var::Vector_Space<Ids...>>`:

- Each field is tracked independently (componentwise mean/variance/count; no covariance).
- Merge/update/add/scale are applied componentwise using the scalar `Moment_statistics` logic.
- Updating with a sample `Vector_Space` uses the per-field sample values.

## Usage Notes

- `legacy/mcmc.h` still defines `logL_statistics` and `logEv_statistics` as thin wrappers around
  `Moment_statistics<logL>` and `Moment_statistics<logEv>`.
- Including `legacy/moment_statistics.h` avoids pulling in the full MCMC header.
- Because `Va` is treated as a variable-like type, plain `double` will not work without adaptation.

## If You Need Extensions

Consider adding:

- Accessors for mean/variance/count in a non-`var::get` style.
- A `merge` helper that documents the formula explicitly.
- Tests for numeric stability if the stats are used in new contexts.
- A covariance-aware (joint) variant when you need full cross-covariance rather than per-field
  variance.
