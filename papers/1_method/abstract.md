# Abstract — what it must do, and a draft

> Updated: 2026-07-20. Paper 1 (method). Reframed from the single-paper usage-map draft: **the
> least-squares opener moved to paper 2.** Paper 1 is the within-family method paper, in the
> gating-dominated regime; its abstract opens on the recursive filter that already exists and the
> question it left unanswered, not on what the field does wrong.
>
> Title lives in `title.md`; the pivot means the three live title candidates need re-checking against
> this scope. Owns the abstract and the Impact Statement while drafting.

## The job the abstract has to do

Answer, early: **given that a recursive likelihood for macroscopic currents already exists (MacroIR,
Comm Biol 2025), what does this paper add?** Not the method — that is published. What is new is (1)
the validation machinery that tests any such likelihood against the process it approximates, possible
because the forward process is exactly simulable while its likelihood is not; (2) the finding, read
directionally, that the members below IR misreport the information the data carry; and (3) the
mechanism, made concrete by the `R → MR → VR → IR` ladder, of *which part* of the interval structure
earns the calibration — the variance form (MR→VR) and the boundary gain (VR→IR) isolated one at a time.

**What is no longer this abstract's job:** the least-squares premise. That the field fits the mean
current by least squares, and what that costs, is **paper 2's** opener (`../_program/paper-2.md`). The
verified literature behind it (Clerx 2019, IonBench, Milescu 2005, Stepanyuk 2014, Del Core & Mirams
2025) is preserved in the shared bibliography and git history; do not re-import it as this paper's
frame.

## The literature this paper's abstract does stand on

`R` ≈ Moffatt 2007 and Münch 2022 (a published Bayesian Kalman filter); `IR` = MacroIR (Comm Biol
2025). The recursive filter is published and in use. What was never done is characterize *which
conditioning* it needs to stay calibrated, and where even it degrades. Two open questions in the
field's own words land here (full quotes in `docs/bibliography/` at the repo root, and `introduction.md`):

- **Milescu 2005**, on recursion: whether estimates are "intrinsically biased if obtained by any method
  that is **not** a Bayesian filtering algorithm", and his own dominant error source, "the local time
  correlation of the current" — exactly the correlation term of the distortion decomposition.
- **Moffatt 2007**, on averaging: a formulation that "can deal with time-averaged signals is therefore
  fundamental". MacroIR lifted that barrier; this paper measures what the intermediate members lose by
  not lifting it fully.

## Hard constraints

Accuracy constraints, not taste. Carried from the title review; still binding.

1. **Distortion, not loss.** Bidirectional: some approximations over-state the information the data
   carry, others under-state it. Never "information loss".
2. **The approximation distorts, not the averaging.** Time averaging is the physical reality; IR handles
   it correctly. Blaming the averaging inverts the causation.
3. **Continuous, not a threshold.** Avoid "breaks down at", "fails beyond".
4. **Scope is stated, not hidden**, and **in band terms.** Two-state model, non-stationary protocol,
   simulated macroscopic currents, **the gating-dominated regime**. Name the companion papers for the
   rest (`../_program/program.md` §7). The stated scope is what keeps the paper from reading as though
   it lived only where IR wins by construction.
5. **No method promotion.** No "exact", "unbiased", "resolves a longstanding problem". That is Comm
   Biol's claim.
6. **Register.** Plain scientific prose, the register of the 2007 and 2025 papers. No aphorisms, no
   triads, no em-dashes.
7. **Do not claim the test.** The diagnostic is classical (Huber 1967, White 1982; generalized form
   Golden, Henley, White & Kashner): the score has mean zero at the truth, its covariance equals the
   Fisher information under correct specification, the gap is the sandwich, and testing it is White's
   information matrix test. "There was no way to test whether a likelihood is faithful" is false and a
   statistically trained reviewer will catch it. The claim that holds: this test had never been applied
   to macroscopic-current likelihoods; applying it is possible because the process is exactly simulable
   while its likelihood is not; what is new is the matrix read directionally, the mechanism ladder, and
   the map. Cite the classical work in Introduction and Methods (`introduction.md`): importing a known
   diagnostic into a field that had not used it is a strength.

## What stays out, and why

- **The least-squares premise and the whole-noise-axis map.** Paper 2.
- **The Fisher-to-zero result.** Your call (11 July audio: *"un poco demasiado"*). Results and
  Discussion, not the abstract.
- **The evidence-correction derivation.** One clause at most; derived in a later component.
- **Micro, more than two states, the stationary regime, experimental data.** Out of scope; the
  Discussion's open-doors sentence, and `micro_IR` appears only as paper 1's one attribution anchor.
- **The acronyms.** Describe functionally in the abstract; name in the paper.
- **MacroIR as a new method.** Named once, parenthetically, for searchability. Not presented.

## Draft (paper 1)

> Macroscopic ion-channel currents carry information about the gating kinetics in their fluctuations,
> and using it requires a likelihood for a signal that each sample averages over a finite acquisition
> window. A recursive filter that conditions on the acquisition interval was recently introduced
> (MacroIR), but which part of that conditioning is load-bearing, and where the filter itself
> degrades, was never characterized. Here we test a graded family of recursive Gaussian likelihoods
> against the process they approximate, which is possible because the forward process can be simulated
> exactly although its likelihood cannot. For each member we ask whether the standardized residuals
> are white with unit variance, whether the score vanishes at the true parameters, and whether the
> covariance of the score matches the Fisher information the likelihood reports; the mismatch between
> the last two defines an information distortion matrix. In a minimal two-state channel model, in the
> regime where gating dominates the instrumental noise, we find that conditioning the conductance on
> the interval mean is not enough: a likelihood that uses the interval mean but the wrong variance
> misreports the parameter uncertainty, and only conditioning on both interval endpoints, which keeps
> the boundary cross-covariance in the filter gain, restores calibration. The same matrix supplies
> corrected estimates of parameter bias and variance. The result is a mechanism, not a ranking: it
> says which piece of the interval structure a macroscopic-current likelihood must keep, and what is
> lost by dropping each piece.

~220 words, over the eLife budget; cut per the order below. The closing sentence is the "so what": a
mechanism a reader can act on, not a leaderboard.

## Length

**150 to 200 words** (the live eLife guide; the LaTeX class still says ≤150). Write to 150, treat 200
as the ceiling. eLife rule that binds us: *"if the biological system is not in the title, it must be in
the abstract"* — name the system (macroscopic ion-channel currents) in the first sentence and state the
band scope in the middle. Cut order, if over: the corrected-bias clause; the residual-whiteness test
(keep the two score tests the matrix is built from); the "not a ranking" restatement if space is tight,
since the mechanism sentence already carries it.

## Register: the slot structure eLife abstracts have

| Slot | Does | Words |
|---|---|---|
| 1 | the object, and how it is measured | ~25 |
| 2 | the limitation, one line (here: which conditioning is load-bearing was never characterized) | ~15 |
| 3 | "Here we test X against Y" — the whole paper in one sentence | ~40 |
| 4 | validation of the approach | ~20 |
| 5–7 | results, one carrying the mechanism (MR→VR→IR) | ~45 |
| 8 | "The result is a mechanism …" | ~25 |

Our draft spends too long on the epistemic setup; that is the first cut.

## Impact Statement (eLife submission-form field)

One sentence, **15–30 words**, **third person**, complements the title, no unfamiliar acronyms, no "We
show…". Drafts **[Q]**:

1. *A simulation-based test shows which part of a macroscopic-current likelihood's interval
   conditioning keeps its reported uncertainty honest, and what is lost by dropping each part.* (26)
2. *For recursive likelihoods of time-averaged ion-channel currents, conditioning on both interval
   endpoints is what restores calibration; using the interval mean alone is not enough.* (26)
3. *Testing kinetic likelihoods against an exactly simulable process reveals that only boundary
   conditioning keeps their reported parameter uncertainty faithful in the gating-dominated regime.* (24)

Draft 1 leads with the mechanism and the payoff and avoids the title's wording; recommended, pending the
title pivot.

## Verify before submission

- **The MR → VR → IR mechanism sentence** against the runs. VR's sign (does the wrong-variance member
  come out over-confident?) is a prediction until measured; the abstract must not state it as found
  before `decisions.md`'s VR run exists.
- **The overconfidence numbers**, if any are quoted, at a fixed n_sims (`../_program/machinery.md` §6.1).
- **"Only boundary conditioning restores calibration"** — the band-A result; state it scoped to the
  gating-dominated regime, not globally.
