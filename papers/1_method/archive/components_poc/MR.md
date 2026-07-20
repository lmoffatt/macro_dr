# MR — dimension

> type: `entity` (dimension) · group: `algorithms` · governed by `../_SPEC.md`
> **Canonical owner of MR's DEFINITION only.** Every fact/figure/datum/claim about MR is its OWN node (DRY);
> here = intrinsic coordinates + pointers. No cross-cutting facts here — that was `_LOG` failure-mode #2.

## Definition `[LOCKED nomenclature.md:28 · run-20260418-182133/script.macroir · 2026-07-14]`
`recursive = true` · `av = 1`.
- recursion axis: occupancy covariance IS propagated between intervals (a filter).
- conductance axis `M`: interval-mean conductance conditioned on the **start** state only (one endpoint).

## Pointers
- **neighbors:** `[[R]]` (av=0) · `[[IR]]` (av=2) · `[[MNR]]` (non-recursive, av=1) · `[[recursion-axis]]` · `[[window-axis]]`
- **code:** `legacy/qmodel.h` (av=1 branch) · `build_likelihood_function` flags (`recursive_approximation`, `averaging_approximation`)
- **data:** `[[data-1c2ae6f]]` — MR gaussian cells, noise 0.1, canonical N_ch
- **facts:**
  - `[[fact-MR-drops-boundary-crosscov]]` — MR omits `N·γᵀΣγ` that IR keeps (relational MR↔IR). src `nomenclature.md:52`
  - `[[fact-MR-overconfident]]` — emp/Fisher 1.5–2.1. owner `[[D-4-ranking]]`
  - `[[fact-window-nonmonotone]]` — MR worse than R. owner `[[window-axis]]`
- **claims:** `[[MR-sign]]` — observable-variance vs parameter-covariance. owner `decisions/D-4`
- **figures:** `[[fig-cond-Nch]]` · `[[fig-ranking]]` (pending)
- **decisions:** `decisions/D-4` (ranking; owns MR's verdict)
- **not-to-confuse:** `[[MRT]]` — Taylor variance-correction variant, cut from paper
- **open:** observable `y_var` direction MR vs IR (settle by `awk` over Fig 3/4 dumps) · `MR-sign` reconciliation (D-4)
