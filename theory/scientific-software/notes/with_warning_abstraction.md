# `With_warning<T>` — carrying soft diagnostics alongside results

## Motivation

The codebase uses `Maybe_error<T>` to encode a strict success/failure dichotomy.
On success it returns `T`; on failure it returns an `error_message` explaining
what went wrong, and the caller chains via short-circuit. This works well when
"something is wrong" implies "stop", but it has no slot for an intermediate
state: **"the result is usable but a soft invariant was missed."**

A concrete example is the simplex-canary family
(`to_Probability`, `to_Probability_displacement`, `to_Covariance_Probability`).
These functions check structural identities of the form
`|Σ x − target|/√N ≤ ε`. Real bugs sit at ratios near 1; pure FP-rounding sits
at machine precision; and there is a middle band — accumulated floating-point
drift from a slightly-imperfect upstream (e.g. `t_P`'s derivative payload not
being exactly row-stochastic) — where the canary fires but the result is still
mathematically usable. Currently the choices are:

- raise `ε` and lose all visibility into the middle band, or
- keep `ε` tight and reject computations that downstream consumers can absorb.

Neither matches the actual structure of the problem. The intent is to
**deliver the result and surface the warning** so the caller (or a logging
pipeline) can record it without crashing the run.

## Proposed abstraction

```cpp
template <class T>
struct With_warning {
    T value;
    std::string warning;   // empty ⇒ no warning attached
};

template <class T>
using Maybe_warning_error = Maybe_error<With_warning<T>>;
```

Semantics:

- `Maybe_warning_error<T>` distinguishes three states:
  - hard error: `Maybe_error::error(...)` — propagates via existing machinery.
  - soft warning: `Maybe_error::value(With_warning{result, msg})` with non-empty `warning`.
  - clean success: `Maybe_error::value(With_warning{result, ""})`.

- Callers that don't care about warnings can ignore the field — `.value().value`
  yields the underlying result.

- Callers that *do* care iterate, collect, log, or sink warnings as they see fit.

- Composition: a sequence of `Maybe_warning_error<T>`-returning calls can fold
  warnings together (concatenation, deduplication, or aggregation per-canary).

## Why the existing `Maybe_error<T>` won't do

`Maybe_error<T>` collapses warnings into one of two extremes:

- treat them as errors → run aborts on soft-FP drift → loses figures-of-merit
  the soft drift would otherwise allow us to compute.
- treat them as silent successes → the structural canary becomes purely
  decorative; we cannot tell whether a long run accumulated dangerous drift.

A side channel (`std::cerr`, a global counter, a per-thread sink) papers over
this in the short term — but it isn't composable, isn't typed, and gets lost
the moment the computation crosses a thread boundary or is serialized.

## Interim policy (until `With_warning<T>` lands)

The canaries adopt a two-tier scheme with the existing `Maybe_error<T>`:

- `ratio > eps_error` → hard error via `Maybe_error::error(...)` (existing path).
- `eps_warn < ratio ≤ eps_error` → log the violation to `std::cerr` once per
  call site, return success.
- `ratio ≤ eps_warn` → silent success.

This is **deliberately ugly**: it loses composability, can spam stderr, and is
not testable. The `With_warning<T>` migration is the proper replacement.

## Migration plan (sketch)

1. Land `With_warning<T>` and `Maybe_warning_error<T>` aliases.
2. Convert the simplex canaries first (`to_Probability`,
   `to_Probability_displacement`, `to_Covariance_Probability`) since they are
   the immediate motivation.
3. Add a fold helper `merge_warnings(Maybe_warning_error<T>...)` for callers
   that chain multiple canary checks.
4. Decide on a sink: per-thread vector of warnings drained at run boundaries,
   vs. immediate stderr emission, vs. structured log (CSV/JSON) for batch runs.
5. Extend gradually to other points in the pipeline where the
   success-with-caveat pattern shows up (e.g. `to_Transition_Probability`'s
   stabilizer policy).
