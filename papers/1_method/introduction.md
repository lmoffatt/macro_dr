# Introduction — what it must do, and a paragraph plan

> Working doc, same genre as `abstract.md` and `title.md`. Opened 2026-07-14.
> The literature block it rests on is the one verified in `abstract.md` (three independent checks, 2026-07-14); it is not repeated here except where the Introduction has to quote it.
> Comment and correct inline.

## The job

The abstract has three sentences to say what is new. The Introduction has about eight hundred words to make the reader *feel the gap before it is filled*. If the reader arrives at "we provide that test" already convinced that no such test existed and that the field is currently choosing likelihoods by habit, the paper is sold. If they arrive there thinking "this is the MacroIR paper again", it is not.

So the Introduction has one hard job and one soft one:

- **Hard:** establish that the epistemic gap is real, is the field's own, and is open right now. The evidence for this is quotations, not assertions. The field has written the gap down in its own words twice (Milescu 2005 on recursion, Moffatt 2007 on averaging) and is actively arguing about it in 2025 (Del Core and Mirams on whether filtering is worth its cost).
- **Soft:** make the *method of attack* feel inevitable. The forward process can be simulated exactly and its likelihood cannot be computed exactly. That asymmetry is the whole licence for the paper. Once the reader sees it, "use the simulator as ground truth and test the likelihood against it" is the obvious move, and the reader half-invents the paper before reading the Results.

## Hard constraints

The five accuracy rules from the title review (`title.md`) and the six from `abstract.md` apply verbatim here. The two that bite hardest in an Introduction:

1. **Do not re-announce MacroIR.** The bridge is published (Comm Biol 2025). The Introduction cites it as prior work, in the same breath as Milescu 2005 and Münch 2022, and does not present it. The moment the Introduction starts explaining how MacroIR works, the paper has become a method paper and the reviewer says "you published this".
2. **The approximation distorts, not the time averaging.** Time averaging is the physical reality. Any sentence of the form "time averaging degrades the information" inverts the causation.
3. Also: distortion is bidirectional (never "loss"), continuous (never "breaks down at"), and the scope (two states, simulation) is stated, not hidden.

## The reader

Two readers, and the Introduction has to hold both.

- **The electrophysiologist** who fits kinetic schemes to macroscopic currents. They fit the mean by least squares, probably in a deterministic ODE framework, and they have never computed a score. For them the Introduction must explain, plainly, what the fluctuations carry and why a likelihood is what it takes to use it. They will not know what a sandwich estimator is and must not need to.
- **The inference-methods reader** (the Mirams / Münch / Del Core axis, and the statisticians they work with). For them, the Introduction must say precisely what statistical object is being tested and against what, because they will otherwise assume this is another benchmark-by-RMSE exercise. The words that buy their attention are *misspecification*, *score*, and *information matrix equality*. One sentence each, no more.

## Paragraph plan

Six paragraphs. The order is chosen so the gap is fully open before the paper is mentioned.

### P1 — What the fluctuations carry, and what the standard practice throws away

Open on the data, not on the method. A macroscopic current is the sum over a channel population; its mean carries the kinetics only through the deterministic relaxation, while its fluctuations carry, in addition, the number of channels, the unitary conductance and, through the temporal structure of the noise, the rates themselves. The standard practice is least squares on the mean current: Clerx et al. 2019 (*Four Ways to Fit an Ion Channel Model*), where all four ways are deterministic ODE fits with an RMSE objective; IonBench (Owen and Mirams 2025), the community benchmark, which states outright that it does not include optimisation approaches for stochastic models.

The cost of this is not a matter of taste, and the field has said so: Del Core and Mirams 2025, *"standard methods for deterministic models do not distinguish between stochastic channel gating and measurement error noise, resulting in biased estimates."* Quote it. It is the strongest possible opening evidence and it is the field auditing itself.

**Do not** claim here that everyone should switch to likelihoods. The paper does not demonstrate that. Claim only that the fluctuations carry information and that the mean-fit discards it.

**Naming trap to avoid:** non-stationary fluctuation analysis *is* widely used and *does* use the variance. It must be named and distinguished, or an electrophysiologist reader will think the gap is already filled. It returns the unitary current, the channel number and the peak open probability; it is a two-moment regression, not a likelihood, and it does not identify a kinetic scheme. Stepanyuk et al. 2014 say it themselves: *"the unitary current is virtually the only parameter that can be reliably obtained from this type of analysis"*, and *"To the best of our knowledge, kinetic rates have never been estimated for any synaptic receptors in their intrinsic environment."* One sentence plus that second quotation ends the objection.

### P2 — Using the fluctuations requires a likelihood, and the ones that exist are approximations along two axes

The small literature, named as a lineage rather than a list: Celentano and Hawkes 2004 (full covariance, 12 to 36 hours per fit), Milescu, Akk and Sachs 2005 (diagonal likelihood, in QuB), Moffatt 2007 (recursive filter), Stepanyuk et al. 2011 and 2014 (semiseparable covariance), Münch et al. 2022 (Bayesian Kalman), Del Core and Mirams 2025. Milescu's own framing of the niche: *"exponential fitting remains the most widely used analytical method ... few procedures for direct estimation of rate constants are available."*

Then the structural point, which is the paper's organising idea and belongs here rather than in the Theory section: **every one of these makes the same two Gaussian approximations, and they differ along exactly two axes.**

- **Axis 1, recursion.** Whether the hidden occupancy distribution is conditioned on the data as it arrives (a filter) or propagated open-loop from the initial condition.
- **Axis 2, the acquisition window.** Whether each sample is treated as an instantaneous observation of the process, or as what it physically is, an average of the process over a finite window.

On axis 2 the literature is unanimous and it is unanimous in the wrong place: **no published likelihood for macroscopic currents integrates the observable over the acquisition window.** The three quotations that establish it (all in `abstract.md`) are Milescu 2005 (*"the method does not handle the effects of decreased variance due to low-pass filtering"*), Celentano and Hawkes 2004 (who mention the window only to dismiss it), and Münch 2022, who file finite integration time under *"misspecifications of the likelihood"* and prescribe an experimental workaround: *"the sampling should be faster than the fastest eigenvalues to avoid biased results."*

This paragraph is also where the family gets its honest attribution: **the algorithm family is the literature, not an invention.** `NR` is essentially Milescu 2005; `R` is Moffatt 2007 and Münch 2022; `IR` is MacroIR (Comm Biol 2025). `MR` and `NMR` are the natural intermediates that fill the grid. Say this plainly and early. It converts what could look like five strawmen into a map of the actual field, and it is the single sentence that most protects the paper from the "you invented four bad methods to beat" objection.

### P3 — The field has stated the gap twice, in its own words

Short paragraph, two quotations, no commentary beyond framing.

- **On recursion.** Milescu 2005: *"It remains an open question whether all estimates should be intrinsically biased if obtained by any method that is not a Bayesian filtering algorithm."* He also names his own dominant error source as *"the local time correlation of the current"*, which is, verbatim, the correlation term of our distortion decomposition. That is not a coincidence worth hiding; it is the strongest single link between the 2005 diagnosis and the 2026 measurement, and it should be said here in half a sentence and paid off in the Discussion.
- **On averaging.** Moffatt 2007: *"the measurements are considered to be instantaneous, where electrophysiological measurements always result from some time-averaging ... A formulation of the algorithms that can deal with time-averaged signals is therefore fundamental to the widespread use of the algorithms presented here."*

Self-citation caution: the second quote is the author's own. Quote it without ceremony, in the same list-of-the-field register as Milescu's. It is evidence, and it happens to be ours.

### P4 — Why they stayed open: nobody could test whether an approximate likelihood is faithful

This is the pivot, and it is the paragraph the whole Introduction exists to reach.

The two questions above are not hard because the algebra is hard. They are hard because **an approximate likelihood is a misspecified model, and there was no accepted way to ask, of a given approximation, how badly it misrepresents the process.** Goodness of fit does not answer it: an approximation can fit the mean current beautifully and still report an information matrix that is off by a factor of twenty. Comparing two approximations to each other does not answer it either, because both may be wrong.

Then the licence, stated as the asymmetry:

> The forward process can be simulated exactly. Its likelihood cannot be computed exactly. That asymmetry is what makes the test possible: simulation supplies the ground truth that the likelihood is an approximation *of*, so any candidate likelihood can be interrogated against data drawn from the very process it claims to describe.

And the statistical machinery that the asymmetry unlocks, in one sentence each so the methods reader knows exactly what is coming: the score has zero expectation at the truth for a correct likelihood, and the covariance of the score equals the information the likelihood reports (the information matrix equality). Under misspecification both fail, and *how* they fail is measurable. Name it: this is the classical misspecification apparatus (Huber, White), which is standard in statistics and, as far as we have found, has never been used to audit a likelihood approximation in this field.

**[Q]** Do we cite White 1982 / Huber 1967 here or in Methods? Recommendation: one citation here, because it tells the methods reader that we are not inventing a diagnostic, we are importing a known one into a field that had not used it. That is a strength, not a concession.

### P5 — What this paper does, with the scope stated

Three things, in the order the paper delivers them:

1. **A test.** Three checks, at the likelihood level: standardized residuals white with unit variance; the score vanishing at the true parameters; the covariance of the score matching the Fisher information the likelihood reports.
2. **A measurement.** The mismatch between the last two defines an **information distortion matrix**: the object that says, in units the reader already understands (parameter variance), by how much and in which directions the approximation misreports what the data know.
3. **A map.** Across channel number, acquisition interval and instrumental noise, for the five members of the family, in a minimal two-state channel model, in simulation, with no experimental data.

The scope sentence is not a hedge and it should not be written apologetically. The minimal model is a *choice*: it isolates the error the algorithm itself commits, with no confounding from scheme misspecification, no identifiability pathology of a rich model, and a regime where the exact simulation is cheap enough to run thousands of replicates. Say that in one sentence and the reviewer who was reaching for "only two states" puts the pen down.

### P6 — What it buys

The closing paragraph is the "so what", and it must land on the live argument rather than on a general appeal to rigour.

Del Core and Mirams 2025 argue, in print and against filtering, on cost: *"approaches that are based on filtering techniques are computationally expensive because an integration step of the differential moment equations (DMEs), and the corresponding updates in the correction step, are computed between every consecutive time points where the measurements are collected."* The field is deciding, right now, whether the filter is worth paying for, and it has nothing with which to decide. A distortion map is exactly what that decision needs: it says which approximation misreports the uncertainty, by how much, and in which corner of the regime.

One clause, not more, on the further payoff: the same matrix is the correction that keeps Bayesian model comparison valid under an approximate likelihood, which is where this work goes next. Forward-pointing, derivation deferred (`00_plan.md` §3).

## What stays out of the Introduction

- **The mechanics of MacroIR.** No boundary state, no Kalman gain, no interval kernel. That is the Theory section.
- **The Fisher-to-zero result.** It is a Results and Discussion item. Foreshadowing it here weakens the Results.
- **The research-program framing** (`../_program/research_program.md`, three bridges). It is genuinely the reason this paper exists, but an Introduction that opens with a research program reads as a grant proposal and *invites* the "this is a component of a bigger thing, come back when it's done" review. The field's own two open questions are a stronger and more self-contained motivation. **The program belongs in the Discussion**, in one paragraph, after the evidence has been delivered.
- **The Comm Biol P2X2 biology.** One citation, no case study.
- **Any claim about experimental data.** There is none in this paper.

## Open questions

- **[Q] How much statistics vocabulary in P4?** The tension is real: *score*, *information matrix equality* and *misspecification* are what convince the methods reader and are what lose the electrophysiologist. Current recommendation: use all three words, and immediately gloss each in the same sentence in physical terms ("the score, the gradient of the log-likelihood, which must average to zero at the true parameters if the likelihood is right"). Cost: about forty words. Worth it.
- **[Q] Where does the two-axis framing (recursion × window) first appear?** It is proposed above for P2, which means the Introduction, not the Theory section, is where the family is organised. That is deliberate: it makes the five algorithms look like a coordinate system rather than a zoo. If it feels heavy there, the fallback is a one-line version in P2 and the full ladder in Theory (the endpoint ladder in `../_program/nomenclature.md` is the cleanest presentation and could carry it).
- **[Q] Do we name the algorithms in the Introduction at all?** Recommendation: yes, once, in P2, tied to their sources (NR ≈ Milescu, R ≈ Moffatt/Münch, IR = MacroIR). The alternative (describe functionally, name only in Methods) reads coy given that the acronyms then run through every figure.
- **[Q] Length.** eLife's guidance is under 5,000 words for the main text excluding Methods, references and legends. Six paragraphs at ~130 words is ~800, which is right for a paper whose Results carry six figures. Do not let the literature review swell: the quotations are load-bearing, the paraphrase around them is not.

## Verify before submission

- **The "no published likelihood integrates the observable over the acquisition window" claim.** It is the Introduction's strongest and most falsifiable sentence, and one counterexample from an adjacent field (a Kalman filter with an integrated observation in cardiac or synaptic modelling) turns it into a liability. It survived three independent checks in July 2026 (see `abstract.md`) and the prior-art map (`docs/bibliography/MacroIR_prior_art_map.md`) already concedes the *device* is known outside the field (integrated-measurement augmented Kalman). Phrase it as scoped to macroscopic ion-channel likelihoods, and cite the prior-art map's concession in the Discussion, not here.
- **Del Core and Mirams 2025 and Owen and Mirams 2025** are both 2025 and both central to the opening. Re-read them at submission; if either has been revised, the quotations must be re-checked against the version of record.
</content>
</invoke>
