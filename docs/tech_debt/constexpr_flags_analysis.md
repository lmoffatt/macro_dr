# Compile-time vs runtime algorithm flags — tech-debt note

**Status:** open tech-debt. Architecture works correctly; analysis below
recommends migration when bandwidth allows.

## Current architecture

`Likelihood_Model_constexpr` carries algorithm choice as **template
parameters**: `uses_recursive_aproximation<bool>`,
`uses_averaging_aproximation<int>`, `uses_variance_aproximation<bool>`,
`uses_taylor_variance_correction_aproximation<bool>`,
`uses_micro_aproximation<bool>`, `uses_taylor_qdt<bool>`. Each combination is a
distinct type. `Likelihood_Model_regular::get_variant()` enumerates the
Cartesian product over allowed values, returning a
`std::variant<Likelihood_Model_constexpr<...>...>` with up to ~20 alternatives
under the current restricted domains. Dispatch happens via `std::visit` and
`if constexpr` inside each lambda.

## Theoretical benefits of compile-time flags

1. No runtime branch in the per-sample inner loop.
2. Inlining / dead-code elimination — false branches removed entirely.
3. Type-level safety: a `recursive=true` model can't accidentally call a
   `recursive=false`-only function.

## Why these benefits are negligible here

The inner loop iterates ~10⁴ samples per recording. Per sample:

| operation | cost (cycles) |
|---|---|
| `calc_Qdt` (eigen or Taylor) | ~10⁶ (matrix work, BLAS/LAPACK) |
| per-cell Bayes update | ~10⁴ (M² element ops) |
| half-dozen runtime `if (flag)` branches | ~25 |

Branch overhead is ~2.5e-5 of total runtime. **Inlining doesn't help much**
either — matrix routines are external library calls (BLAS/LAPACK), already
tight, can't be inlined further.

## Real costs of compile-time flags

1. **Variant explosion.** 6 dimensions → up to 96 variants in the worst case.
   Restricted domains in `build_likelihood_function` bring it to ~20
   actually-used variants. Each is a full template instantiation of every
   consumer function — `dlogLikelihoodPredictions`, `dlog_Likelihood_micro`,
   `calculate_mlikelihood`, write_csv overloads, etc.

2. **Compile time.** With ~20 variants × ~6 entry points × ~2 calls per
   site, ~240 instantiations per translation unit consuming the variant.
   Empirically ~30–60 s per file for `likelihood.cpp`. Adding the 7th flag
   (`taylor_qdt`) immediately doubled compile times until we restricted
   its domain to the micro branch only.

3. **Binary size.** Each variant emits its own object code; binary is
   5–10× larger than equivalent runtime-dispatched code. Slower link, worse
   instruction-cache behavior.

4. **Maintenance / extensibility.** Adding `taylor_qdt` took ~80 lines of
   plumbing across 6 files (constexpr Var class, `Likelihood_Model_constexpr`
   template, `Likelihood_Model_regular` cartesian/variant_impl/get_variant,
   `build_likelihood_function`, all `std::visit` lambdas, DSL bindings,
   figure_2 scripts).

5. **Phantom variants.** Adding a flag to *all* domains uniformly creates
   useless variants (e.g., `taylor_qdt=true` for the macro Kalman path,
   which doesn't go through the micro lifting where the flag matters).
   Mitigation: restrict the domain per branch, but this is brittle and the
   restriction logic must be re-evaluated each time domains change.

## Verdict

**Runtime flags would be a 10–100× compile-time win at <0.001× runtime cost.**

A future refactor would:

- Drop `constexpr_Var` / `constexpr_Var_domain` for algorithm flags.
- Make `Likelihood_Model` a single non-templated struct with plain
  `bool`/`int` fields.
- Keep one templated path per **major algorithm shape** (macro Kalman vs
  lifted micro — that's the genuine type-level distinction because the
  state representations differ). That's 2–3 functions, not 20+.
- Inside each path, `if (recursive)` / `if (taylor_qdt)` runtime branches.
- `std::visit` collapses to a normal function call.
- Type-level safety is preserved at the major-shape level (which is what
  actually matters); flag-combination errors become runtime errors with
  clear messages, which is acceptable since they're DSL-script-driven and
  always fail fast at the start of a run anyway.

## Why we're keeping the constexpr architecture for now

It works. The current diagnostic build (numerical Fisher info, GFD,
within_v, taylor_qdt, etc.) needs to ship. A refactor of this scope is
multi-day work touching every consumer of `likelihood_algorithm_type`.
Better to land the diagnostic infrastructure on the existing architecture,
then schedule the runtime-flags refactor as a separate cleanup pass.

## Trigger for revisiting

If any of these become true, prioritize the refactor:

- `likelihood.cpp` compile time exceeds ~2 minutes.
- We need to add 1–2 more algorithm flags (e.g. a 3rd Qdt method, or a
  per-cell weighting choice).
- The `merge_Maybe_variant` chain in `build_likelihood_function` grows
  beyond 5 branches.
- Build / link issues from binary size start affecting iteration.
