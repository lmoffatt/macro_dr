# Migration plan: inline derivative tests → catch2 tests

## Background

`test_derivative_clarke` (and `test_Derivative`) are runtime numerical checks
that compare an analytical derivative against finite differences along every
parameter direction. They were added inline inside the macro fold body
(`Macro_DMR::log_Likelihood`) during a development phase to catch derivative
regressions immediately. They do their job, but they leave ~330 lines of test
scaffolding embedded in the production hot path, gated by `constexpr bool`
toggles that can only be flipped by editing the source and rebuilding.

The cleaner pattern is to move these tests into dedicated `tests/`
catch2 files: production code stays clean, verification is always available
via `ctest`, and the tests run on small fixed inputs that exercise the same
code paths as production.

This was already done for the micro path:

- inline `if constexpr (test_micro_step_derivative) { ... }` blocks were
  drafted inside `dlog_Likelihood_micro` then immediately reverted.
- the verification lives at
  [`tests/microir/test_micro_derivatives.cpp`](../../tests/microir/test_micro_derivatives.cpp),
  using a small `scheme_CO` setup with `Num_ch_mean=5`
  ([`tests/data/scheme_CO_small_par.csv`](../../tests/data/scheme_CO_small_par.csv)).

This document is the matching plan for the macro path, which still has all
the inline blocks.

## Inventory: inline test sites in `Macro_DMR::log_Likelihood`

All inside the per-step lambda of the fold (`legacy/qmodel.h`).

| line | quantity tested                          | gate                                          | suggested destination                            |
| ---- | ---------------------------------------- | --------------------------------------------- | ------------------------------------------------ |
| 5620 | `constexpr bool test_derivative = true;` | —                                             | delete                                           |
| 5621 | `constexpr bool test_macroir_derivative = true;` | —                                     | delete                                           |
| 5659 | `Qx` analytical derivative                | `if constexpr (test_Qx)` (local true)         | `tests/macroir/test_qx_derivatives.cpp` (new)    |
| 5675 | eigen-decomposition derivative            | `if constexpr (test_eigen && false)` (already disabled) | `tests/macroir/test_qx_derivatives.cpp` (new) |
| 5705 | `Qdtm` (avg=0, stabilized) derivative     | `if constexpr (test_derivative)`              | `tests/macroir/test_qdt_derivatives.cpp` (exists) |
| 5726 | `Qdtm` non-stabilized fallback diagnostic | nested under above                            | same file as 5705                                 |
| 5732 | `Qdtm` non-stabilized fallback diagnostic (text variant) | nested under above             | same file as 5705                                 |
| 5762 | per-step `MacroR2` (the "macroir" inner step) | `if constexpr (test_macroir_derivative && averaging::value>0)` | `tests/macroir/test_macroir_derivatives.cpp` (already partially covers this — currently `[scaffolding]`-tagged and excluded from CI) |
| 5800 | `Qdt` (avg≥1) derivative                  | `if constexpr (test_derivative)`              | `tests/macroir/test_qdt_derivatives.cpp`          |
| 5813 | `Qdt` non-stabilized fallback diagnostic  | nested under above                            | same file as 5800                                 |
| 5853 | `Qdtm` (avg≥1) derivative                 | `if constexpr (test_derivative && averaging::value>0)` | `tests/macroir/test_qdt_derivatives.cpp` |
| 5864 | `Qdtm` non-stabilized fallback diagnostic | nested under above                            | same file as 5853                                 |
| 5936 | `Qdt` (final segment) `test_Derivative`   | `if constexpr (test_derivative)`              | `tests/macroir/test_qdt_derivatives.cpp`          |

11 active call sites (one already disabled), spanning ~330 lines.

## Migration plan

The micro precedent at
[`tests/microir/test_micro_derivatives.cpp`](../../tests/microir/test_micro_derivatives.cpp)
is the reference shape. Each new catch2 test should:

1. Load `scheme_CO` (or whichever model has the smallest viable state count
   for the quantity being tested).
2. Use a small-Nch parameter file in `tests/data/`.
3. Load the existing `Moffatt_Hume_2007_ATP_time_*` recording (already used
   by the macroir test) for step iteration.
4. Walk the fold step-by-step, calling `test_derivative_clarke` per step on
   the relevant `calc_*` / `MacroR2` lambda. `REQUIRE` the boolean result.

### Step 1 — Qdt / Qdtm tests

Extend [`tests/macroir/test_qdt_derivatives.cpp`](../../tests/macroir/test_qdt_derivatives.cpp)
with TEST_CASE blocks for each Qdt / Qdtm site (lines 5705, 5800, 5853, 5936 in
`qmodel.h`). The non-stabilized fallback diagnostics (5726, 5732, 5813, 5864)
are interesting but secondary — they fire only when the primary test fails,
to help localize whether the stabilizer is causing the regression. They can
either be folded into the same test as conditional follow-up checks, or
omitted on first migration and re-added later if a regression actually shows
up needing them.

### Step 2 — Qx / eigen test

New file [`tests/macroir/test_qx_derivatives.cpp`](../../tests/macroir/test_qx_derivatives.cpp)
covering the Qx site (5659) and optionally the currently-disabled eigen site
(5675 — gated by `test_eigen && false` so it has been off for a while; check
git history before re-enabling to understand why it was disabled).

### Step 3 — promote macroir-step test out of `[scaffolding]`

The existing `tests/macroir/test_macroir_derivatives.cpp` is tagged
`[macroir][derivatives][macroir]` with the comment "scaffolding". It is
excluded from CI by the `CATCH2_FILTER=~*scaffolding*` setting at
[`tests/CMakeLists.txt:55`](../../tests/CMakeLists.txt#L55). Once that test
is verified stable, drop the `scaffolding` tag (and the description suffix)
so it runs in CI alongside the new micro test.

### Step 4 — delete inline blocks from `qmodel.h`

After steps 1–3 are merged and the new tests pass reliably:

- delete the `constexpr bool test_derivative` and `test_macroir_derivative`
  toggles at qmodel.h:5620–5621.
- delete each `if constexpr (test_*) { ... }` block at the line numbers in
  the inventory above.
- the `Macro_DMR::log_Likelihood` body shrinks by ~330 lines and reads as
  pure algorithm.

## Why we are doing this

- Production code stays free of verification scaffolding (smaller, easier
  to read, faster to compile).
- Verification is always available via `ctest`; no rebuild, no env-var
  gating, no source edits to enable.
- Regressions are caught by CI rather than by happening to fire the inline
  blocks during a particular dev run.
- New tests added later (regimes that surfaced as bugs in production) are
  natural TEST_CASE additions to existing files — that's how
  `tests/macroir/test_macroir_derivatives.cpp` is intended to grow.

## What we are NOT doing

- The `MACRODR_DX_ASSERT` structural-correctness asserts (20 sites in
  `parameters_derivative.h` / `derivative_operator.h` / `qmodel.h:5628-5637`)
  are **not** numerical tests — they verify that a `Derivative<...>` carries
  its `dx()` pointer. They are cheap (single pointer null check), gated by
  the `MACRODR_DX_ASSERT_DEBUG` preprocessor symbol, compiled out unless
  explicitly enabled. They are fine where they are; do not migrate them.
- The actual `test_derivative_clarke` implementation in
  `legacy/derivative_test.h` stays put — only its callers move.
- The `MACRODR_PATH` env-var pattern stays as-is for path resolution.

## Order of operations (when we get to it)

1. Add `tests/macroir/test_qx_derivatives.cpp` covering site 5659.
2. Extend `tests/macroir/test_qdt_derivatives.cpp` with sites 5705, 5800,
   5853, 5936.
3. Run `ctest` on both, confirm they pass on `scheme_CO` (and ideally
   `scheme_1` and `scheme_10_d` as a multi-model spot-check — separate
   TEST_CASEs).
4. Delete inline blocks from `qmodel.h` in a single commit so the diff
   reads as "remove inline tests, keep inline algorithm" cleanly.
5. Drop the `[scaffolding]` tag from `tests/macroir/test_macroir_derivatives.cpp`
   in the same or a follow-up commit.
6. Update [`tests/CMakeLists.txt:55`](../../tests/CMakeLists.txt#L55) to
   remove the `~*scaffolding*` filter once nothing is tagged that way.
