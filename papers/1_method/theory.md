# Theory section — what it must do, and a build plan

> Updated: 2026-07-20. Paper 1 (method). Covers the manuscript section *"The macroscopic interval
> likelihood and its approximations"* (`docs/manuscript-drafts/elife_paper.tex` §2), which eLife
> permits between Introduction and Results ("A Methods/Model section may appear after the Introduction
> where sensible", `../_program/elife-author-instructions.md`).
> Naming lives in `../_program/nomenclature.md` and is assumed here.

## The job

This section makes paper 1's four methods **one object graded by how much of the interval structure the
likelihood uses**, so the Results read as a traverse of a coordinate system rather than a bake-off
among unrelated codes. If the reader leaves able to say "there are two Gaussian approximations, and
these methods differ only in how much of the acquisition interval the hidden state is conditioned on,
and in whether the interval variance is handled correctly", the section has done its work.

**Scope, stated once.** Paper 1's roster is `R`, `MR`, `VR`, `IR`: the recursive, gating-aware family
in the gating-dominated regime. The prior question — whether to model the gating fluctuations at all,
which is where least squares and the non-recursive members live — is the *root* of the two-level
structure and belongs to paper 2 (`../_program/program.md` §1). This section is the lower level, the
ladder. Say that in one sentence and move on; do not re-derive the root question here.

Two pages, for an audience that mostly does not read filtering theory.

## The one-sentence spine

> A macroscopic current is the sum over a population of channels, each a continuous-time Markov chain;
> the exact likelihood of a time-averaged recording would require the distribution of the
> interval-averaged conductance of the whole population, which is intractable, so every practical
> likelihood replaces it with a Gaussian, and these methods differ in *what they condition that
> Gaussian on* and *how they account its variance*.

Everything below is that sentence, unpacked.

## Structure

### T1 — The observable, and why it is not a sample of the process

The physical statement first, because it is the one thing the field has systematically not modelled.
The recorded value at sample *t* is not the current at an instant; it is the average of the current
over the acquisition window Δ, after the anti-aliasing filter. What the likelihood must describe is the
distribution of

  ȳ_t = (1/Δ) ∫ over the interval of the population current, plus instrumental noise,

conditioned on the past. The distinction is invisible when Δ is much shorter than the fastest
relaxation and dominant when it is not.

**Register warning (`title.md`).** The averaging is the reality; it is not a degradation. Write "the
observable is an interval average" and never "the averaging degrades the signal".

### T2 — The two Gaussian approximations, named once and used everywhere

The conceptual core (`00_plan.md` §1), displayed as a boxed pair, because the Results refer back to it.

1. **The macro (occupancy) approximation.** The joint distribution of channel occupancies is exactly
   multinomial when the channels are independent (the maximum-entropy closure given the mean
   occupancies). It is replaced by a multivariate Gaussian, valid for large N_ch by the central limit
   theorem. It degrades toward few channels, the **multinomial regime** — **which is paper 3's
   subject.** Paper 1 meets it only at its N_ch = 10 floor, where one micro anchor cell attributes
   IR's degradation there to this closure rather than the next one (`decisions.md`).
2. **The interval-likelihood approximation.** The distribution of the interval-averaged conductance is
   a complicated stochastic-telegraph object. It is replaced by a Gaussian, valid when the interval
   contains many transitions or when instrumental noise dominates. It degrades toward intervals much
   shorter than the relaxation time, the **telegraphic regime** — **this is paper 1's subject.**

And the sentence that ties theory to diagnostics:

> The higher moments discarded by these closures do not disappear. They reappear in the data as
> temporal correlation the likelihood does not predict, which is exactly what the correlation term of
> the distortion decomposition measures.

Three regimes fall out, named here so the Results have vocabulary: **multinomial** (few channels,
paper 3), **telegraphic** (very short intervals, paper 1's failure edge), **Gaussian** (many channels,
moderate intervals, enough instrumental noise). Say plainly that IR, being an approximation, *must*
fail somewhere, and that paper 1's job is to show its failure is the predicted degradation of the
interval closure and not a bug. That sentence is the paper's intellectual honesty in one line.

### T3 — The ladder: the four methods as one graded object

The heart of paper 1. The roster is a **monotone progression in how much of the interval structure the
likelihood uses**, and `VR` is the rung that makes the progression fine enough to isolate a mechanism.

| Rung | Method | Conductance conditioned on | Interval variance | recursion / av |
|---|---|---|---|---|
| 0 | `R` | no endpoints (instantaneous conductance; the averaging is ignored) | — | yes / 0 |
| 1 | `MR` | one endpoint (the interval-mean given the start state) | **total** per start state, `gsqr_i − gmean_i²` | yes / 1 |
| 2 | `VR` | one endpoint (same as MR) | **residual**, `Σⱼ P_ij·gvar_ij` | yes / 1 |
| 3 | `IR` | two endpoints (the boundary; interior marginalized) | residual, and the boundary cross-covariance kept in the **gain** | yes / 2 |
| — | (exact) | the full trajectory | — | intractable; the simulation supplies it as ground truth |

The three steps up the ladder are the paper's argument, each isolating one thing:

- **R → MR:** start conditioning the conductance on the interval (the averaging enters the mean).
- **MR → VR:** *only the variance changes.* MR carries the total per-start-state variance, which
  double-counts the spread across end states; VR uses the residual (boundary-conditioned) variance.
  Nothing else moves.
- **VR → IR:** *only the gain changes.* IR conditions on both endpoints and keeps the boundary
  cross-covariance N·γᵀΣγ in the gain; VR does not.

This is the structural payoff, and it is new. The claim "MR's problem is the gain, not the variance"
was previously an algebraic assertion. The ladder **measures** it: if VR (residual variance, no
boundary gain) comes out *more* over-confident than MR rather than calibrated, then removing the
variance without gaining the boundary information made things worse, and the missing piece is the gain,
exactly as the algebra says. VR is the control that turns the assertion into a result.

Two further things the table buys:

- **IR is the top rung below intractability.** "IR is the one that stays calibrated" stops being an
  empirical surprise and becomes what the ladder predicts. The Results *confirm* a structure.
- **The non-recursive siblings sit off this ladder, in paper 2.** `NR` and `NMR` share the conductance
  column with `R` and `MR` but drop the recursion; they belong to the "is a likelihood needed at all"
  question and are measured there. Say so in one sentence so a reader does not think the roster is
  arbitrary.

**Literature on the ladder.** `R` ≈ Moffatt 2007 and Münch 2022 (a published Bayesian Kalman filter);
`IR` = MacroIR (Comm Biol 2025). Say it here as well as in the Introduction: the recursive filter is
published, and what paper 1 adds is the account of *which part of the interval structure* earns the
calibration.

### T4 — The boundary state and the two conductance objects, defined once

Definition, in the paper's words (`../_program/nomenclature.md`; the scoping sentence is mandatory):

> A **boundary state** is the pair (i₀, i_t) of a channel's states at the two ends of an acquisition
> interval. Chaining intervals, the end state of one is the start state of the next, so there is a
> single state variable per junction: the filter conditions on the states at the interval boundaries
> and leaves the trajectory between them free, marginalizing it analytically. This is the same device
> as static condensation, or the spatial Markov property, applied on the time axis.

And immediately: it is **not** a transition state in the mechanistic sense. Avoid the word *transition*
anywhere near it (`project_boundary_state_naming`).

The two conductance objects the ladder needs:

- `M`: the interval-averaged conductance conditioned on the state at the interval's start, a K-vector,
  (γ̄₀)_{i₀} = Σ_{i_t} P_{i₀→i_t}(Δ) · Γ̄_{i₀→i_t}.
- `I`: the same average conditioned on both endpoints, so the object stays K×K instead of collapsing.

And the variance objects that separate MR from VR, which the ladder now needs displayed because `VR`
depends on it (`theory/macroir/notes/gvar_i_overcount_audit.md`, mechanism only):

- **total** per start state: `gsqr_i − gmean_i²` (what MR uses);
- **residual** (expected boundary-conditioned): `Σⱼ P_ij · gvar_ij` (what VR and IR use).

They differ by `Var_j[gmean_ij | i]`, the spread of the interval-mean across end states. **That term
is the whole MR→VR step**, and it should be an equation.

**The single algebraic difference that VR still lacks and IR has:** the boundary cross-covariance
N·γᵀΣγ, kept in IR's gain, dropped by VR. Display it; it is the mechanism behind the last rung.

### T5 — The filter step, and the emission variance

The Kalman-like update: predictive occupancy mean and covariance, the predicted observable mean and
variance, the gain, the correction. One display of each. Salvage from the archived derivations
(`docs/manuscript-drafts/archives/elife-macroir-merged.tex`, `revised3.tex`: the interval average, the
E₂ kernel, the Q⊕Q factorization, the meta-state box, the E₃ variance kernel).

**The emission-variance decomposition is load-bearing.** The verified code-level decomposition
(`project_emission_variance_decomposition`) is

  y_var = e + N·gΣg + N·ms

instrumental noise, plus the population gating variance through the occupancy covariance, plus the
within-interval (per-channel) term. Reconcile the manuscript skeleton's fill-hint
(`eps^2/t + white + pink + N_ch[gating]`) against this and write the version the code computes.

**Correction to the earlier draft (D-2, `../_program/axes.md`).** An earlier version of this note
called the swept `Current_Noise` a "misnomer" the Methods must apologize for. That was overturned: the
figure axis is a **dimensionless** noise, ν = `Current_Noise`·k_off/g² (label = 10·ν), a legitimate
natural unit, and the two crossovers of the noise axis are what make paper 1's map two-dimensional.
State the unit; do not call it a misnomer.

### T6 — What the exact likelihood would be, and why the simulation is ground truth

Close by naming the object at the top of the ladder: the exact likelihood conditions on the full
trajectory and is intractable, which is precisely why the simulator, which *can* realize the full
trajectory exactly, is the only available ground truth. This hands the Results their licence, the same
asymmetry stated in the Introduction, now with the algebra behind it. Two sentences.

## What stays out

- **Any k-state generality beyond what the two-state model needs.** Carry K symbolically where it costs
  nothing (the ladder table, the M/I definitions); do not develop the general case. **[Q]** Is there a
  result that needs K > 2? If not, the general derivation goes to the Supplement.
- **The non-recursive members and least squares.** Paper 2. Name them once (T3) as the off-ladder
  siblings; do not measure them here.
- **The Taylor variance-correction variants (IRT, MRT).** Cut. **And keep them clear of `VR`:** the
  engine flag is `taylor_variance_correction`, `VR` is not a Taylor variant, and Methods must say the
  two are unrelated or a reader will conflate them (`../_program/program.md` §9).
- **The PSD trust coefficients (α_μ, α_σ).** Methods at most (`project_psd_trust_redundant`).
- **The Kalman prior-art connection.** Discussion, not here (`discussion.md`). It reads defensively in
  a Theory section.
- **The distortion machinery.** `../_program/machinery.md`; the section that presents it is
  `diagnostics.md`.

## Open questions

- **[Q] Theory section, or Methods?** Recommendation: keep a compact Theory section (T1–T4, two pages,
  a few displays), push the filter algebra (T5) into Methods with a pointer, keep T6 as Theory's last
  paragraph.
- **[Q] How much algebra survives in the main text?** Proposal: four displayed equations. (i) the
  observable as an interval average; (ii) the boundary-conditioned conductance; (iii) the total-vs-
  residual variance difference `Var_j[gmean_ij]` (the MR→VR step); (iv) the boundary cross-covariance
  N·γᵀΣγ (the VR→IR step). The last two are the mechanism the paper measures. Test: could an
  electrophysiologist read the Results with only those four?
- **[Q] Is the "maximum entropy closure" framing kept?** Recommendation: keep, as a parenthesis; it
  explains why the discarded moments reappear as correlation.

## Sources to lift from

- `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md`, `macroir_derivation.tex` — the algorithm.
- `docs/manuscript-drafts/archives/elife-macroir-merged.tex`, `elife-macroir-revised3.tex` — the
  boundary-state derivation, the E₂/E₃ kernels, the Q⊕Q factorization, the emission-variance
  decomposition.
- `theory/macroir/docs/Macro_MRT/macromrt_macromrt_paper_section.md`,
  `theory/macroir/notes/gvar_i_overcount_audit.md` — the MR/VR variance mechanism (use the mechanism,
  **not** the inverted May-2026 verdict).
- `../_program/nomenclature.md` — the ladder, the boundary-state definition, the scoping sentence.

## Verify before submission

- **The emission-variance decomposition**, against the code, not the archived drafts (T5). The drafts
  predate the `gvar_i` fix and the SymmetricMatrix storage fix.
- **The total-vs-residual variance forms and the N·γᵀΣγ term** — the two mechanism equations, re-derived
  from the current code before printing. `VR`'s very definition is the residual form, so if the code's
  `gvar_i` is not what T4 claims, the roster is wrong, not just the prose
  (`feedback_verify_dont_assume`).
