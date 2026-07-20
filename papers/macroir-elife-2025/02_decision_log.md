# Paper 2 (MacroIR / eLife 2025) — Decision Log

> Updated: 2026-07-20 (usage-map reframe: LSE as method zero, noise axis through the gating crossover).

Ledger of **settled** decisions, so we can "rewind" via git history. Operational plan: `00_master_plan_v2.md`. Framing: `From molecular mechanisms to data back and forth PROGRAM.md`. Open (unsettled) decisions live in `00_master_plan_v2.md` §8 (D-0…D-7); they move here once decided.

## Confirmed (current, MacroIR 13, 2026-07)

### Framing
- The paper is **bridge 2** of the research program: the robustness study of the likelihood algorithm. Comm Biol 2025 (P2X2) is the demonstration (phase 1); this characterizes the components (phase 2).
- **Thesis (revised 2026-07-20):** the paper delivers a **usage map for the methods available to fit macroscopic currents**, from least squares on the mean through the five Gaussian likelihood approximations, saying per regime which is the cheapest method that still reports its uncertainty honestly. Not "present MacroIR", and no longer "MacroIR is the only calibrated member" either: that sentence was true only of the band of the noise axis that had been sampled. Rationale and the band structure: `00_master_plan_v2.md` §0 and §1a.

### Scope
- Minimal **two-state** model (`scheme_CO`), single K_on/K_off. **Non-stationary** protocol only.
- **Six methods on two levels** (settled 2026-07-20). Method zero, off the lattice: classical nonlinear least squares on the mean current, data key `nonlinearsqr`, display label `LSE`, engine flag `family_approximation = 2`. Then the five macroscopic likelihood algorithms `NR`, `NMR`, `R`, `MR`, `IR`. Standardize naming on **NMR** (scripts have used MNR). Published-name bridge: IR = MacroIR, NMR = MacroINR.
- **LSE is not a rung of the family.** In the dispatcher it carries the same two knob settings as NMR (`recursive=false, averaging=1`) and is distinguished only by the third flag. The "one object with two knobs" framing is retired; the structure is a root question (model the gating fluctuations or not) with the five-member ladder hanging from it.
- **The noise axis is swept through the gating-noise crossover** (settled 2026-07-20). Previously every cell sat at or below the single-channel noise scale, the band that favours the gating-aware likelihoods by construction. Band definitions in `00_master_plan_v2.md` §1a; units per `decisions/D-2_parameter_units.md`.
- **Likelihood-only** analysis. The posterior information-distortion framework and full model-comparison results are OUT (later program components). The **likelihood-side evidence correction** (volume ½ log det C, effective-sample α⋆ = p/tr C) stays IN as **motivation**; derivation deferred to the bridge-3 study.
- MLE / Gauss-Newton local max kept, only to obtain the empirical parameter covariance.
- Out of scope = later program components, not exclusions: MicroIR, >2 states, stationary regime, experimental data, Taylor variants (IRT/MRT).

### Method
- **Gaussian Fisher** is the distortion anchor; the numerical finite-difference Fisher only gauges how good the Gaussian one is.
- Diagnostics: residual mean/var/whiteness; score bias; Var[score] vs Gaussian Fisher → distortion matrix C (sandwich), decomposed into correlation + sample/geometric distortion; plus the direct empirical-vs-sandwich covariance test.
- Distortion/bias evaluated at the optimum / θ_pool; θ_sim exposes the bias.

### Ranking → the usage map (revised 2026-07-20)

**There is no longer a single ranking to settle.** The verdict is per band, and only band A has been measured. The former ranking is retained in `00_master_plan_v2.md` §4 as the band-A column, with its two contested cells still owned by `decisions/D-4_ranking_verdict.md`. Nothing about bands B and C is settled, and nothing about LSE is settled anywhere.

Retired phrasings, so they are not re-copied out of the older documents: "IR sole survivor", "only MacroIR stays calibrated across the practical regime", "MR strawman". The first two overstate a result measured in one band; the third is a ranking word with no meaning on a map, where a method has a domain, possibly empty.

### Figures / data
- Regime maps over **N_ch × noise**, with the acquisition **interval** as a third axis carried *inside* every battery cell (7 values of `interval_in_tau`). Matches the abstract (channels, interval, noise). Definitive figures anchored on the Gaussian Fisher.
- **Band-B/C fill dispatched 2026-07-20.** `dispatch_figure_3_LSE.sh` for `nonlinearsqr` and `dispatch_figure_3_G.sh` for {IR, R, MR, NMR}, both at `NCHS="10 100 1000 10000"`, `N_SIMS=1000`, `GROUP_SIZE="10 100"`, and `N_NOISE="100 1000 10000 100000"`. Since `N_NOISE` is paired index-wise with `NCHS`, that is `label = 10·N_ch`, which is the B/C boundary at Δ = τ; the interval axis then carries each cell from the boundary into band C. **NR is not in this dispatch** (see the new D-0).
- **Do not pool cells across n_sims.** The band-B/C fill is `n_sims = 1000`, the canonical band-A cells are 10000, and `433ed13` also holds 200. Distortion scalars carry a Jensen bias in n_sims, so mixing them inside one panel manufactures an apparent regime effect. Hold n_sims fixed per panel or use the debiased quadratic.
- **D-0 settled (2026-07-15):** freeze at **`1c2ae6f`**; multi-commit provenance accepted (each CSV self-stamps its engine hash). `433ed13` kept as the numerical-Fisher equivalence demo, not re-run; `87889e6` (micro) out of scope. The only fill is the four non-IR algorithms {NR, NMR, R, MR} × canonical N_ch {10, 100, 1000, 10000} × noise {0.1, 1, 10} onto the Gaussian anchor — 36 cells, ~6,500 CPU-h, running on Dirac; complete when that batch finishes. E-1…E-5 decoupled to `main` as code hygiene (seeds were never logged; existing ensembles are statistically equivalent, not bit-reproducible, and cannot be fixed retroactively). Detail: `decisions/D-0_freeze_and_rerun_scope.md`.

### Repo / reproducibility
- **One paper = one repo.** This paper carves out to a dedicated repo (venue-agnostic name, e.g. `macroir-validity`) at **code freeze** (when the Gaussian-Fisher figures are final); `macro_dr` (public) referenced by pinned tag/hash. See `09_carve_plan.md`. Until freeze, keep developing in the monorepo and running on Dirac from it.
- Head manuscript is `elife_paper.tex`; superseded drafts stay as history.
- Audio sources and transcripts both tracked (user preference).

## Working plan (not final): Comm Biol erratum (gvar_i)
- Intent: disclose, done properly. **Decouple** the erratum (correctness: re-run at the same fidelity to isolate the gvar_i fix's effect on the Bayes factors) from any **Bessel-filter high-fidelity redo** (a separate study; bundling it confounds attribution). **Triage first**, cheaply, via the distortion correction and/or reweighting the deposited MCMC samples, to learn whether the fix moves the ranking before committing to a re-run. Low external urgency (small field), but load-bearing for this paper's use of Comm Biol as the demonstration.

## Superseded (kept for rewind)
- **"Five algorithms" as the closed roster** → six methods on two levels (2026-07-20). LSE was previously present only as cited background in the Introduction and abstract, describing what the field does; it is now a measured arm.
- **The ranking as the deliverable** → the usage map (2026-07-20). The ranking survives as one band's column.
- **The noise axis as a nuisance third dimension** → the map's organizing axis (2026-07-20), because it carries the two crossovers that decide which method is needed.
- Validity-map axes N_ch / Δτ / Noise → briefly recorded as **N_ch × K_off** → corrected 2026-07-15 back to **N_ch × noise** (+ interval internal). K_off was never swept (fixed at 100 on disk); the "N_ch × K_off" framing never had data behind it. See D-0.
- LID↔Evidence recorded as a standalone "finding" (Δlog Z = ½ log det C) → now the bridge-3 **motivation** (½ log det C volume + α⋆ peak), likelihood-side, derivation deferred.
- "Keep MicroIR out to reduce attack surface" → reframed as a named later program component.
- `elife-macroir-merged.tex` as manuscript source of truth → now `elife_paper.tex`.
