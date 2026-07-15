# sandwich — operation

> type: `operation` (PROV activity: inputs → output) · group: `diagnostics` · governed by `../_SPEC.md`
> **Defined ONCE here.** Consumers carry `computed-via: [[sandwich]]`; the applications list is a DERIVED
> backlink (`_SPEC` Rule 5/7) — `grep -rl '\[\[sandwich\]\]' components/`. Never hand-maintained.

## Signature `[LOCKED 03_metrics_diagnostics.md · decision_log "Method" · 2026-07-15]`
- **inputs:** `H` = Gaussian Fisher (the *bread*) · `J` = covariance of the score (the *meat*)
- **output:** `C` = information-distortion matrix; sandwich covariance = `H⁻¹ J H⁻¹`
- robust/QMLE estimator (White). `C` anchored on the **Gaussian** Fisher (not the numerical FD Fisher).

## Exact algebra — NOT owned here `[OPEN → theory/.../supplement_information_distortion_main.tex]`
Owner = the supplement. Two facts already locked elsewhere, pointed to (not restated):
- IDM reconstruction = `K·CDM·Kᵀ`, `K = H^(−1/2) J_s^(1/2)` — **NOT** `SDM^½·CDM·SDM^½` (that printed form is false).
  `[LOCKED correction_idm_reconstruction.md — E-2]`
- decomposition: `[[correlation-distortion]]` (R, cross-interval score corr = Milescu 2005's error) `+` `[[sample-distortion]]` (3rd+4th cumulants).

## Pointers
- **inputs (dimensions/facts):** `[[H-gaussian-fisher]]` · `[[J-score-covariance]]`
- **outputs:** `[[C-distortion-matrix]]`
- **downstream operations:** `[[correlation-distortion]]` · `[[sample-distortion]]` (consume C)
- **code:** `src/core/likelihood.cpp` (IDM call site ~3501; E-2 fix pending)
- **theory:** `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex` · `03_metrics_diagnostics.md`
- **applications:** DERIVED — do not list here. `grep -rl '\[\[sandwich\]\]' components/`
