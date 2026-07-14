# Theory section — what it must do, and a build plan

> Working doc, same genre as `abstract_draft.md`. Opened 2026-07-14. Comment and correct inline.
> Covers the section the manuscript skeleton calls *"The macroscopic interval likelihood and its approximations"* (`docs/manuscript-drafts/elife_paper.tex` §2), which eLife explicitly permits between Introduction and Results ("A Methods/Model section may appear after the Introduction where sensible", `docs/elife-author-instructions.md`).
> Naming decisions live in `nomenclature.md` and are assumed here.

## The job

This section exists to make the five algorithms **one object with two knobs**, so that the Results read as a traverse of a coordinate system rather than a bake-off among five unrelated codes. If the reader leaves this section able to say "there are two Gaussian approximations, and the family differs only in how much of the acquisition interval the hidden state is conditioned on", the section has done its work, and every figure afterwards is interpretable.

It has to do that in roughly two pages, for an audience that mostly does not read filtering theory.

## The one-sentence spine

> A macroscopic current is the sum over a population of channels, each a continuous-time Markov chain; the exact likelihood of a time-averaged recording would require the distribution of the interval-averaged conductance of the whole population, which is intractable, so every practical likelihood replaces it with a Gaussian, and the members of the family differ in *what they condition that Gaussian on*.

Everything below is that sentence, unpacked.

## Structure

### T1 — The observable, and why it is not a sample of the process

The physical statement first, because it is the one thing the field has systematically not modelled. The recorded value at sample *t* is not the current at an instant; it is the average of the current over the acquisition window Δ, after the anti-aliasing filter. What the likelihood must therefore describe is the distribution of

  ȳ_t = (1/Δ) ∫ over the interval of the population current, plus instrumental noise,

conditioned on the past. The distinction is invisible when Δ is much shorter than the fastest relaxation and dominant when it is not.

**Register warning (rule 2 of `title_options.md`).** The averaging is the reality; it is not a degradation. Write "the observable is an interval average" and never "the averaging degrades the signal".

### T2 — The two Gaussian approximations, named once and used everywhere

This is the conceptual core (from `00_master_plan_v2.md` §1) and it should be displayed, probably as a boxed pair, because the entire Results section refers back to it.

1. **The macro (occupancy) approximation.** The joint distribution of channel occupancies is exactly multinomial when the channels are independent, which is also the maximum-entropy closure given the mean occupancies. It is replaced by a multivariate Gaussian, valid for large N_ch by the central limit theorem. It degrades toward few channels (the **multinomial regime**).
2. **The interval-likelihood approximation.** The distribution of the interval-averaged conductance is a complicated stochastic-telegraph object. It is replaced by a Gaussian, valid when the interval contains many transitions or when instrumental noise dominates. It degrades toward intervals much shorter than the relaxation time (the **telegraphic regime**).

And the sentence that ties the theory to the diagnostics, which is worth its space:

> The higher moments discarded by these closures do not disappear. They reappear in the data as temporal correlation the likelihood does not predict, which is exactly what the correlation term of the distortion decomposition measures.

Three regimes fall out and should be named here so the maps in the Results have vocabulary: **multinomial** (few channels), **telegraphic** (very short intervals), **Gaussian** (many channels, moderate intervals, enough instrumental noise). Say plainly that MacroIR, being an approximation, *must* fail somewhere, and that the paper's job is to show its failure is the predicted degradation of these two closures and not a bug. That sentence is the paper's intellectual honesty in one line and it should not be cut for space.

### T3 — The endpoint ladder: the family as a coordinate system

Use the ladder from `nomenclature.md` verbatim; it is the best presentation we have and it makes the verdict feel structural.

| Conditioned on | Members | Conductance model |
|---|---|---|
| no endpoints | `NR`, `R` | instantaneous conductance; the averaging done by the acquisition is ignored |
| one endpoint (the start) | `NMR`, `MR` | interval-mean conductance given the initial state |
| two endpoints (the boundary) | `IR` | interval-mean conductance given both boundary states; interior marginalized |
| the full trajectory | (exact) | intractable; supplied by the stochastic simulation as ground truth |

Crossed with the second axis, recursion (is the occupancy covariance propagated between intervals, conditioned on the data, or not), the five members are exactly the implemented grid: `NR` (no, av=0), `R` (yes, av=0), `NMR` (no, av=1), `MR` (yes, av=1), `IR` (yes, av=2).

Two things this table buys, both of them worth more than the space:

- **MacroIR is the top rung below intractability.** Stated this way, "IR is the one that stays calibrated" stops being an empirical surprise and becomes what the ladder predicts. The Results then *confirm* a structure instead of reporting a winner.
- **The literature sits on the ladder.** `NR` ≈ Milescu 2005, `R` ≈ Moffatt 2007 and Münch 2022, `IR` = MacroIR. Say it here as well as in the Introduction; a reader who skipped the Introduction must not think we built four strawmen.

### T4 — The boundary state, defined once and scoped

Definition, in the paper's words (from `nomenclature.md`, and the scoping sentence is mandatory):

> A **boundary state** is the pair (i₀, i_t) of a channel's states at the two ends of an acquisition interval. Chaining intervals, the end state of one is the start state of the next, so there is a single state variable per junction: the filter conditions on the states at the interval boundaries and leaves the trajectory between them free, marginalizing it analytically. This is the same device as static condensation, or the spatial Markov property, applied on the time axis.

And immediately: it is **not** a transition state in the mechanistic sense. Avoid the word *transition* anywhere near it. (The term "boundary state" is confirmed; `project_boundary_state_naming`.)

Then the two conductance objects the ladder needs, defined properly (per `nomenclature.md`):

- `M`: the interval-averaged conductance conditioned on the state at the interval's start, a K-vector, (γ̄₀)_{i₀} = Σ_{i_t} P_{i₀→i_t}(Δ) · Γ̄_{i₀→i_t}.
- `I`: the same average conditioned on both endpoints, so the object stays K×K instead of collapsing to a vector.

**The single concrete algebraic difference in the whole family**, and it should be displayed as an equation because it is the mechanism behind the MR verdict: `MR` drops the boundary cross-covariance term N·γᵀΣγ that `IR` keeps.

### T5 — The filter step, and the emission variance

The Kalman-like update: predictive occupancy mean and covariance, the predicted observable mean and variance, the gain, the correction. One display of each. This is where the archived derivations get salvaged (`docs/manuscript-drafts/archives/elife-macroir-merged.tex` and `revised3.tex`: the interval average, the E₂ kernel, the Q⊕Q factorization, the meta-state box, the E₃ variance kernel).

**The emission-variance decomposition is load-bearing and its current statement in the repo is inconsistent.** The manuscript skeleton's fill-hint says `eps^2/t + white + pink + N_ch[gating]`; the verified code-level decomposition (`project_emission_variance_decomposition`) is

  y_var = e + N·gΣg + N·ms

that is, instrumental noise, plus the population gating variance through the occupancy covariance, plus the within-interval (per-channel) term. **These are not obviously the same statement.** Reconcile them against the code before this subsection is written, and write the version the code actually computes. (Also verified there: the swept "noise" parameter is `Current_Noise`, which scales as N⁰ and as 1/Δ, and carries no per-channel conductance and no time constant. It is a misnomer and the Methods must say what it is, or the noise axis of every map is uninterpretable.)

### T6 — What the exact likelihood would be, and why the simulation is the ground truth

Close the section by naming the object at the top of the ladder: the exact likelihood conditions on the full trajectory and is intractable, which is precisely why the simulator, which *can* realize the full trajectory exactly, is the only available ground truth. This paragraph is what hands the Results their licence, and it is the same asymmetry stated in the Introduction, now with the algebra behind it. Two sentences.

## What stays out

- **Any k-state generality beyond what the two-state model needs.** The archived drafts derive the family for arbitrary K. That was written when the paper was a method paper. Now the model is two states; carry K symbolically where it costs nothing (the ladder table, the definitions of M and I) and do not develop the general case. **[Q]** Confirm: is there a result in the paper that needs K > 2 anywhere? If not, the general derivation goes to the Supplement, not the main text.
- **The Taylor variance-correction variants (IRT, MRT).** Cut (`00_master_plan_v2.md` §2). They are not in the family, they are not in the figures, and mentioning them opens a line of questioning about a PSD failure mode we chose not to fix (`project_macromrt_irt_psd_downdate`).
- **The PSD trust coefficients (α_μ, α_σ).** Implementation-level safeguards. They belong in Methods at most, and the α_σ story in particular is a piece of internal history that no reader needs (`project_psd_trust_redundant`).
- **The Kalman prior-art connection.** Discussion, not here (see `discussion_plan.md` D3). It interrupts the derivation and it reads defensively in a Theory section.
- **The distortion machinery.** Separate section; see `diagnostics_plan.md`.

## Open questions

- **[Q] Theory section, or Methods?** eLife allows a Model section after the Introduction, and the Results are unreadable without T2 and T3. Recommendation: keep a compact Theory section (T1-T4, roughly two pages, two displays), push the filter algebra (T5) into Methods with a pointer, and keep T6 as the last paragraph of Theory. That gives the Results a conceptual spine without a five-page derivation in front of them.
- **[Q] How much algebra survives in the main text?** Proposal: three displayed equations only. (i) the observable as an interval average; (ii) the definition of the boundary-conditioned conductance; (iii) the MR-vs-IR difference term N·γᵀΣγ. Everything else to Methods and Supplement. The test: could an electrophysiologist read the Results with only those three? If yes, that is the right cut.
- **[Q] Is the "maximum entropy closure" framing kept?** It is elegant (the multinomial is the max-ent distribution given the mean occupancies, and the Gaussian is its large-N limit) and it explains why the discarded moments show up later as correlation. It also costs a sentence that some readers will bounce off. Recommendation: keep, as a parenthesis.

## Sources to lift from

Pending the theory-source inventory (in progress), the ones named in `00_master_plan_v2.md` §7:

- `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md` and `macroir_derivation.tex` — the algorithm.
- `docs/manuscript-drafts/archives/elife-macroir-merged.tex`, `elife-macroir-revised3.tex` — the boundary-state derivation, the E₂/E₃ kernels, the Q⊕Q factorization, the meta-state box, the emission-variance decomposition.
- `theory/macroir/docs/Macro_MRT/macromrt_macromrt_paper_section.md` and `theory/macroir/notes/gvar_i_overcount_audit.md` — the MR mechanism (use the mechanism, **not** the inverted May-2026 verdict).
- `nomenclature.md` — the ladder, the boundary-state definition, the scoping sentence.

## Verify before submission

- **The emission-variance decomposition**, against the code, not against the archived drafts (see T5). The drafts predate the `gvar_i` fix and the SymmetricMatrix storage fix, and at least one of them states the variance in a form the code does not compute.
- **The N·γᵀΣγ term** as the sole algebraic difference between MR and IR. It is asserted in `nomenclature.md` and in the master plan; it is the mechanism sentence in both the Theory section and the Discussion table, and it should be re-derived from the current code once before it is printed.
</content>
