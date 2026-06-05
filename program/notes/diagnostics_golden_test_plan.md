# Diagnostics numerical regression ("golden") test — plan

## Goal
A CI test that catches *numerical* regressions in the likelihood-derivative
diagnostics (IDM / DIB / Gaussian-Fisher distortion / Probit summaries), not just
structural presence. The current `macrodr_cli_likelihood_diagnostics_contains_dib`
is a **smoke test** — it only asserts that rows/components *exist*, says nothing
about whether the *values* are right. We want to convert "pipeline runs and emits
the expected shape" into "pipeline emits the expected numbers within tolerance."

## Why NOT now (preconditions — do this only once these hold)
- **The diagnostics output contract is still in flux.** Probit set, presets
  (basic / series_var / series_cov / series_kernel[_full]), IDM vs CDM, algorithm
  injection, and even component *labels* are actively changing (the stale
  `Distortion_Induced_Bias.Distortion_Induced_Bias` label just broke contains_dib).
  A golden pinned to a moving target is pure churn. Land this after the contract
  stabilizes.
- **The numerics-identity gate must hold first.** `OMP_NUM_THREADS=1` vs `=N`
  (combos serial) bit-identical, and bootstrap determinism (serial RNG pre-draw,
  per-index placement). Without that there is no stable value to pin.

## Hard design problems (the part that needs real thought)
1. **Platform FP variance — golden must be tolerance-based, never bit-exact.**
   CI (ubuntu, system BLAS), dirac (sequential MKL), Clementina (OpenBLAS, gcc15)
   differ in low-order bits; matrix-exp / eigendecomp amplify it. Measure the
   actual spread of each candidate quantity across the three toolchains and set a
   relative tolerance with margin — do NOT guess rtol.
2. **Ill-conditioning makes some quantities un-pinnable.** Per the κ(F)/√n
   analysis: soft-direction IDM / ratio-of-quadratic-forms quantities are
   heavy-tailed (possibly no finite variance) → pinning them invites flakiness.
   **Pin only well-conditioned, robust quantities.** Avoid mean±SD of ratios and
   tail quantiles.
3. **Bootstrap (seed, B) must be frozen.** Quantiles depend on both; pin them
   explicitly in the fixture. Use the *median*, not the mean, for any
   bootstrap-summarized quantity.
4. **Fixture must be tiny + fixed + fast** (CI seconds) yet numerically
   nontrivial: small Nch, small nsim, fixed seed, 1–2 intervals, a stable preset.

## Candidate quantities to pin (decide deliberately; favour stable + meaningful)
- **FIM eigenspectrum / condition number κ(F)** — deterministic, no bootstrap,
  scientifically central, the most stable. Strongest candidate.
- `log_Det` of `Fisher_Covariance` / `Sample_Distortion_Matrix` (scalar, robust).
- A **stiff-direction** distortion value (well-conditioned by construction).
- A **median** Probit quantile of a headline summary.
- NOT: 0.025/0.975 tail quantiles, soft-direction ratios, anything mean±SD on a
  ratio of quadratic forms.

## Candidate mechanism (sketch — refine later)
- `tests/data/diagnostics_golden.macroir`: fixed model (`scheme_CO_small`), fixed
  cell, fixed (seed, B), stable preset.
- A small expected table: `(component-key → value, rtol)` triples (committed),
  not a whole golden CSV (avoids re-baselining on every cosmetic format change).
- `tools/ctest_assert_csv_values_within.cmake`: for each expected key, find the
  row (keyed by the component/probit/parameter columns), compare
  `|golden − actual| / |golden| < rtol`.
- Before fixing rtol: run the same cell on ubuntu-CI / dirac / Clementina, take
  `max observed spread × margin`.

## Open questions to settle before writing it
- Has the output contract stabilized (presets, probit set, labels)?
- Which 3–5 quantities are both scientifically meaningful AND well-conditioned
  enough to pin?
- One golden cell, or two (a well-conditioned one for value asserts + a
  deliberately ill-conditioned one that asserts only κ(F), never the ratios)?

## Non-goals
- Not a replacement for the per-run scientific validation gate (numerics identity,
  eigenspectrum sanity, RSS plateau). This is a regression *guard* on a fixed tiny
  cell, not a correctness proof.

See also: [program/notes/figure2_parallelization_plan.md](figure2_parallelization_plan.md).
