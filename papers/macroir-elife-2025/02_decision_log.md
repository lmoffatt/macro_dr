# Paper 2 (MacroIR / eLife 2025) — Decision Log

Ledger of **settled** decisions, so we can "rewind" via git history. Operational plan: `00_master_plan_v2.md`. Framing: `From molecular mechanisms to data back and forth PROGRAM.md`. Open (unsettled) decisions live in `00_master_plan_v2.md` §8 (D-1…D-7); they move here once decided.

## Confirmed (current, MacroIR 13, 2026-07)

### Framing
- The paper is **bridge 2** of the research program: the robustness study of the likelihood algorithm. Comm Biol 2025 (P2X2) is the demonstration (phase 1); this characterizes the components (phase 2).
- **Thesis:** not "present MacroIR" but characterize where the macroscopic Gaussian likelihood approximations hold and where they break; MacroIR is the only member calibrated across the practical regime.

### Scope
- Minimal **two-state** model (`scheme_CO`), single K_on/K_off. **Non-stationary** protocol only.
- Five macroscopic algorithms: `NR`, `NMR`, `R`, `MR`, `IR`. Standardize naming on **NMR** (scripts have used MNR). Published-name bridge: IR = MacroIR, NMR = MacroINR.
- **Likelihood-only** analysis. The posterior information-distortion framework and full model-comparison results are OUT (later program components). The **likelihood-side evidence correction** (volume ½ log det C, effective-sample α⋆ = p/tr C) stays IN as **motivation**; derivation deferred to the bridge-3 study.
- MLE / Gauss-Newton local max kept, only to obtain the empirical parameter covariance.
- Out of scope = later program components, not exclusions: MicroIR, >2 states, stationary regime, experimental data, Taylor variants (IRT/MRT).

### Method
- **Gaussian Fisher** is the distortion anchor; the numerical finite-difference Fisher only gauges how good the Gaussian one is.
- Diagnostics: residual mean/var/whiteness; score bias; Var[score] vs Gaussian Fisher → distortion matrix C (sandwich), decomposed into correlation + sample/geometric distortion; plus the direct empirical-vs-sandwich covariance test.
- Distortion/bias evaluated at the optimum / θ_pool; θ_sim exposes the bias.

### Ranking (the verdict)
- **IR** sole survivor (distortion → ~1 in the Gaussian regime, ≤~1.3 in the corners).
- **R** ~factor-2 residual variance distortion.
- **MR** strawman, overestimates variance (drops the boundary cross-cov N·γᵀΣγ that IR keeps).
- **NR** biased, variance inflation ∝ N.
- **NMR** only other unbiased method, underestimates variance; possible speed niche.

### Figures / data
- Regime maps over **N_ch × noise**, with the acquisition **interval** as a third axis carried *inside* every battery cell (7 values of `interval_in_tau`). Matches the abstract (channels, interval, noise). Definitive figures anchored on the Gaussian Fisher.
- **D-0 settled (2026-07-15):** freeze at **`1c2ae6f`**; multi-commit provenance accepted (each CSV self-stamps its engine hash). `433ed13` kept as the numerical-Fisher equivalence demo, not re-run; `87889e6` (micro) out of scope. The only fill is the four non-IR algorithms {NR, NMR, R, MR} × canonical N_ch {10, 100, 1000, 10000} × noise {0.1, 1, 10} onto the Gaussian anchor — 36 cells, ~6,500 CPU-h, running on Dirac; complete when that batch finishes. E-1…E-5 decoupled to `main` as code hygiene (seeds were never logged; existing ensembles are statistically equivalent, not bit-reproducible, and cannot be fixed retroactively). Detail: `decisions/D-0_freeze_and_rerun_scope.md`.

### Repo / reproducibility
- **One paper = one repo.** This paper carves out to a dedicated repo (venue-agnostic name, e.g. `macroir-validity`) at **code freeze** (when the Gaussian-Fisher figures are final); `macro_dr` (public) referenced by pinned tag/hash. See `09_carve_plan.md`. Until freeze, keep developing in the monorepo and running on Dirac from it.
- Head manuscript is `elife_paper.tex`; superseded drafts stay as history.
- Audio sources and transcripts both tracked (user preference).

## Working plan (not final): Comm Biol erratum (gvar_i)
- Intent: disclose, done properly. **Decouple** the erratum (correctness: re-run at the same fidelity to isolate the gvar_i fix's effect on the Bayes factors) from any **Bessel-filter high-fidelity redo** (a separate study; bundling it confounds attribution). **Triage first**, cheaply, via the distortion correction and/or reweighting the deposited MCMC samples, to learn whether the fix moves the ranking before committing to a re-run. Low external urgency (small field), but load-bearing for this paper's use of Comm Biol as the demonstration.

## Superseded (kept for rewind)
- Validity-map axes N_ch / Δτ / Noise → briefly recorded as **N_ch × K_off** → corrected 2026-07-15 back to **N_ch × noise** (+ interval internal). K_off was never swept (fixed at 100 on disk); the "N_ch × K_off" framing never had data behind it. See D-0.
- LID↔Evidence recorded as a standalone "finding" (Δlog Z = ½ log det C) → now the bridge-3 **motivation** (½ log det C volume + α⋆ peak), likelihood-side, derivation deferred.
- "Keep MicroIR out to reduce attack surface" → reframed as a named later program component.
- `elife-macroir-merged.tex` as manuscript source of truth → now `elife_paper.tex`.
