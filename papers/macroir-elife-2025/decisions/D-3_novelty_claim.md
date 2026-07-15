# D-3 — Novelty: concept + evidence map

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

### CONCEDE — true; cite not claim; place LATE and embedded (never open on it)
- MacroIR ≈ integrated-measurement Kalman filter (1e-8). BUT own derivation CTMC-native, exact, O(k³)/interval indep of N — that is **Comm Biol's**, cite it, do not disown it.
- Sandwich / information-matrix-equality = White 1982 / Huber / Godambe. Imported apparatus.
- Both α⋆ and ½ log det C corrections = published objects (Pauli; Lv & Liu).

### LIVE — this paper's actual novelty
- **New likelihoods: MR, MNR** (grid completion). Modest but REAL — B-3's "no new likelihood" was false. MNR should have a speed niche; MR is the cautionary intermediate.
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
- **MR sign** (→ D-4): docs say "overestimates variance"; data say over-confident. One is backwards.
- retracted-claim sweep still open: introduction_plan.md:54,:116; abstract_draft.md:30; theory_plan.md:23
