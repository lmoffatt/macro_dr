# D-3 — Novelty: concept + evidence map

> **BANNER 2026-07-31, updated 2026-08-01.** `MNR`/`NMR` here is the **defective build** (missing the
> `N·ms` interval-variance term); the corrected member is `INR` (`../../_program/nomenclature.md`). The
> two quantitative lines that pool it with `NR` — the `N_ch²` conditioning scaling and the ×10–16
> overconfidence — were measured on that build and must be rechecked against the re-run before they
> enter the paper.
> **The re-run has landed** (`figures/data/1f7138b/`) and closes the bias and information-distortion
> half only: see the first-moment table under "Result magnitudes" and
> `recompute/d3_interval_vs_recursion_2x2.py`. **Both starred lines are still owed**, because each is a
> different statistic (a Fisher-spectrum scaling and a standard-error ratio) and neither can be read
> off that script. Expect real movement, not rounding: `INR`'s distortion is 22 where `NMR` pooled with
> `NR` at 79.

> Substrate, not prose. Concepts and their evidence, telegraphic. The paper text gets BUILT from this;
> nothing here is copy-paste. Rebuilt 2026-07-15 with Luciano (replaces the earlier paste-ready draft,
> which over-conceded "no new likelihood" and opened on a negation).
> Prior-art authority: `docs/bibliography/MacroIR_prior_art_map.md` (non-destructive; killed claims archived there).

---

## The question
What does THIS paper claim as new — distinct from MacroIR (Comm Biol 2025) and from the imported machinery.

---

## Claims by status

### DEAD — prior art kills; do not claim
- **"No macroscopic likelihood integrates the acquisition window"**
  - against: Qin, Auerbach & Sachs 2000 (multi-channel FIR-HMM, filter inside likelihood); Fredkin & Rice 1992; Michalek 2000; **self** (MacroIR does it → self-refuting)
  - verdict: dead as written. Survives only **scoped to macroscopic-N via scaling** (Qin = handful of channels, cost kᴸ + combinatorial in N, never reaches N=10³–10⁴)
- **"First/only exact CTMC integrated-observation filter"**
  - against: Kilic 2021 (exact within-window integral, any K, single molecule, Poisson, MCMC); Bäuerle (PDE exact up to numerics k>2)
  - verdict: dead. Narrow to **many-channel ensemble** (Kilic single-molecule; Bäuerle/Blackwell single-chain)
- **"exact-CTMC more accurate than their LNA/Gaussian"**
  - against: coincide to 1e-8 for linear channel dynamics (`verify_IR_vs_augmentation.py`)
  - verdict: dead as accuracy. It is a realization/efficiency/channel-native point
- **"the device / Kronecker Q⊕Q / k²↔(k+1) equivalence is ours"**
  - against: Zadrozny 88, Harvey 89 (IMKF); Albertsen & Hansen 1994 (Kronecker already in this field); k²-equiv is matrix-analytic folklore
  - verdict: concede
- **"the temporal correlation of macroscopic currents is unused / MacroR is the classical method for it"** (ADDED 2026-07-28)
  - against: **Katz & Miledi 1970/1972** and **Anderson & Stevens 1973** (closing rate α from the Lorentzian POWER SPECTRUM of the macroscopic current; the spectrum is the Fourier transform of the autocovariance, so kinetics have come out of temporal correlation since 1973); **Conti, Neumcke, Nonner & Stämpfli 1980** and **Sigworth 1981** (the exact nonstationary two-time covariance `C(t1,t2) = N i² p(t1)[p11(t2|t1) − p(t2)]`, the same object MacroR propagates); **Celentano & Hawkes 2004** (a covariance likelihood, three years BEFORE MacroR: *"CVF fits both the magnitude of the recorded current and the strength of the correlations between different time points"*)
  - verdict: **dead, and it is a one-line kill from any electrophysiologist over fifty.** Survives only as **validity**: the correlation's kinetic content has been known since the 1970s, a likelihood for it exists, the field still fits the deterministic mean, and **no criterion exists for deciding when a correlation-exploiting method is working**. That is the paper.
  - what still distinguishes us: Sigworth 1981 estimates the covariance **empirically from ensembles of 256-504 repeated sweeps**, uses it to **discriminate** schemes, never fits rates by likelihood, and states no uncertainty. The spectral route needs stationarity. Ours is a single non-stationary record with a calibrated interval.
  - sources: `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`; §A.10 of `docs/bibliography/MacroIR_prior_art_map.md`; PDFs and bib entries filed
- **"nobody has modelled correlated residuals in ion-channel fitting"** (ADDED 2026-07-28)
  - against: **Lei, Ghosh, Mirams et al. 2020** (Phil Trans R Soc A 378:20190349): ARMA(2,2) plus two Gaussian-process discrepancy models on hERG current residuals
  - verdict: dead as absence, **alive as mechanism**. Their model is phenomenological: it correlates residuals in time without deriving the correlation from finite-channel gating, so it widens posteriors without adding information, and they report that *"in some cases, the best predictions were still made by ignoring discrepancy"*. Cite it as the demonstration that the phenomenological route was tried in this domain and did not pay. **Say that their target is MODEL discrepancy, not gating noise, before a referee does.**

### CONCEDE — true; cite not claim; place LATE and embedded (never open on it)
- MacroIR ≈ integrated-measurement Kalman filter (1e-8). BUT own derivation CTMC-native, exact, O(k³)/interval indep of N — that is **Comm Biol's**, cite it, do not disown it.
- Sandwich / information-matrix-equality = White 1982 / Huber / Godambe. Imported apparatus.
- Both α⋆ and ½ log det C corrections = published objects (Pauli; Lv & Liu).

### LIVE — this paper's actual novelty
- ~~**New likelihoods: MR, MNR**~~ **RETRACTED for MNR on priority, and ONLY on priority (revised 2026-08-01).** It cannot be claimed as new by this paper because it is **already published, by us**: `MNR`/`NMR` is `MacroINR`, the control of Comm Biol 2025. **MR survives** as the cautionary intermediate and now lives in a supplement.
  - **The second reason given on 2026-07-29 is dead and must not be re-copied**: "the 'speed niche' does not exist, on the freeze it is numerically indistinguishable from NR, so it is prior work of ours that buys nothing a body member does not". That indistinguishability was `NMR`, the build missing the `N·ms` interval-variance term. The corrected `INR` separates from `NR` in **both** moments: median |bias| in `N_ch` 0.003 against 0.099 log10, and `k_off` information distortion at noise 0.1 / N_ch 10⁴ of 22.2 against 78.9, with the sample component at 1.00 against 47.2 (`recompute/d3_interval_vs_recursion_2x2.py`, on `figures/data/1f7138b/`).
  - **And it buys something no body member does**: it is the only member that isolates the interval-window margin, so it is what turns "the window restores the mean, the recursion restores the variance" from an assertion into a measurement. Being prior work of ours bars a **novelty** claim; it does not bar a **body column**. Those were conflated in the retracted sentence. Roster consequence is Q-5 in `../decisions.md`.
- **The Comm Biol control is this paper's real-data demonstration** (ADDED 2026-07-29). Comm Biol 2025 already ranked **nine** kinetic schemes by Bayesian evidence with MacroIR and repeated the ranking with MacroINR as a control: *"Model ranking was sensitive to likelihood approximation: the control method (MacroINR) systematically underestimated evidence for schemes with conformational intermediates"*, and *"leading to systematically different evidence values (Supplementary Table S1)"*. So **"an invalid likelihood corrupts model comparison" is not a promise in this paper, it is a published, real-data observation**, and what eLife adds is why, by how much, and where. Two consequences to handle deliberately: the winning margin over the top non-conformational alternative is only **6.4** (5.5-6.7), which is the scale a systematic likelihood distortion can move; and those evidences were computed with the version carrying the `gvar_i` defect, which is undisclosed and whose correction is an open question. See `../decisions.md`.
- **Positioning against Münch et al. 2022** (ADDED 2026-07-28, M-6). It is a published Bayesian Kalman filter for ion channels **in the target journal**, so the topic demonstrably clears the desk and the novelty bar is correspondingly higher. The abstract must offer **the test of whether such filters tell the truth about their own uncertainty, plus the map of where they fail** — not another filter. Münch also already published diagnostic (i), residual whiteness, plus an N_ch rule of thumb; what they never compute is the score, the Fisher information, the score covariance or any sandwich. **That is the delta, and it should be stated as the delta.**
- **Measurement, not test** (domain-first; map Part III #1)
  - ev: process exactly simulable → H analytic, J = Monte-Carlo over replicates; no null hypothesis; the literature's object is a test statistic with a notorious finite-sample defect (White/Godambe), ours is a measurement
  - claim with finesse — NOT "it's just White's". The epistemic move (exact simulator as ground truth) is the enabler nobody had
- **Distortion has a physical identity** (map Part III #2)
  - ev: C_sample = 3rd + 4th cumulants of the interval current; R = cross-interval score correlation = **exactly Milescu 2005's named error** ("local time correlation of the current")
  - payoff: says WHICH physical feature each approximation discards and in WHICH parameter direction
- **New empirical facts — the conceptual payload** (this is what makes it eLife, not a validation note)
  - Fisher → 0 on relaxation (where the information lives)
  - conditioning ∝ N_ch² for NR/NMR vs ∝ N_ch for IR (variance-inflation visible in the Fisher spectrum; measured 2026-07-15)
  - window axis non-monotone: MR (×1.5–2.1) WORSE than R (×1.3) — naive averaging backfires; only IR closes it
  - the validity map (which approximation, which regime)
- **Bounded license** (the deliverable, not an apology)
  - IR stays calibrated across the tested regime; frontier explicit: 2-state, p=0.5 plateau, simulation, no experimental data

---

## Result magnitudes — evidence grounding the recursion framing
- recursion = the big fix: NR/NMR **×10–16** overconfident → R **×1.3** (factor ~10)
- interval finishes it: R ×1.3 → IR **×1**; but MR ×1.5–2.1 (non-monotone)
- ⇒ recursion NECESSARY, not SUFFICIENT; IR (recursion + interval, done right) the only calibrated one
- ⇒ Luciano's "recursive methods are under-used and unproven" is the DOMINANT empirical effect, not a hunch

**The first-moment half, added 2026-08-01, and it is a second axis rather than a correction.** Everything above is the second moment. On the first moment the ordering is different and the window, not the recursion, is what acts. Median |bias| over every cell and interval on disk, log10, in `N_ch` / `i`:

| member | window | recursion | N_ch | i | k_off |
|---|---|---|---|---|---|
| `LSE` | — | no | 0.001 | Fixed | 0.001 |
| `NR` | no | no | 0.099 | 0.085 | 0.002 |
| `INR` | yes | no | **0.003** | **0.003** | 0.000 |
| `R` | no | yes | 0.110 | 0.101 | 0.002 |
| `MR` | one end | yes | 0.120 | 0.103 | 0.001 |
| `VR` | one end | yes | 0.058 | 0.066 | 0.001 |
| `IR` | both ends | yes | **0.002** | **0.001** | 0.000 |

- the bias is an **amplitude-pair** effect, ~25–30 %, and it is flat in N_ch; `k_off` is unbiased for everyone
- **recursion alone does not remove it and slightly worsens it** (NR 0.099 → R 0.110 → MR 0.120); the two members carrying the interval window are the two without it
- so the two axes are not two projections of one story: **window → first moment, recursion → second moment**, and IR is the only member that closes both. That is a stronger and more legible claim than the ladder, and it is the 2026-08-01 audio's decomposition
- source: `recompute/d3_interval_vs_recursion_2x2.py`, which also carries the LSE parameter-index trap (LSE renumbers its free vector, so the figure-4 digest mislabels its amplitude rows)

---

## Construction notes (how the text gets built from this)
- **Lead with what we did** (positive). Concessions go late, embedded, one graceful clause. Never open on a negation or a disclaimer.
- **Do not over-disclaim** — referees credit even already-published results; giving novelty away preemptively is a loss, not honesty.
- **MacroIR:** present clean, cite as prior, do not re-explain the mechanics; but drop the terror of "re-announcing" — that was too defensive.
- **Hook = temporal correlation** (Luciano's spine): field discards it (deterministic mean-fits) → recursive likelihoods use it but under-adopted, on cost (Del Core & Mirams 2025) and unproven validity (self, Moffatt 2007) → this paper measures exactly that correlation (Milescu's own term) and maps which approximation captures it. This unifies the recursion framing + measurement + physical-identity.
- Two axes, two projections of one 2D story: recursion (the reader's decision, the hook) × window (where the prior-art gap lived). Family = a coordinate system, not a zoo of strawmen.

---

## Open evidence gaps — close before writing the text
- **Qin 2000**: read (how many channels? what cost?) — the retraction leans on it
- **Michalek**: year is **2000** not 1999 (Crossref-verified, T-3)
- `[VERIFY]` cites: Fredkin & Rice 1992, Fatehi & Huang 2017 (venue/pages)
- ~~**MR sign**~~ **CLOSED (D-4 §2.1): over-confident, ~1.6 at the headline cell, recomputed on the freeze. The prose copies that say "overestimates variance" have it backwards. Do NOT assert an observable-variance direction for MR; the production code points the other way.**
- retracted-claim sweep still open. **Note the targets are stale**: `introduction_plan.md`, `abstract_draft.md` and `theory_plan.md` no longer exist. The live files are `../SPINE.md` (Introduction), `../SPINE.md` (Abstract), `../SPINE.md` (Theory).
- **Every `% src:` in this folder points at deleted files** (`results_plan.md`, `discussion_plan.md`, `abstract_draft.md`, `methods_plan.md`, `00_master_plan_v2.md`, `02_decision_log.md`, `diagnostics_plan.md`, `theory_plan.md`, `introduction_plan.md`). LINT-SRC (`../01_writing_plan.md` §1) requires every number to trace to the file that computed it, so the chain is currently broken for the whole novelty map. D-4 was repaired 2026-07-28 by recomputing and **committing** its scripts under `recompute/`; do the same here rather than repointing at prose.
- **Lei 2020 author list** was taken from the arXiv version; verify against the published Phil Trans version before it enters the manuscript.
