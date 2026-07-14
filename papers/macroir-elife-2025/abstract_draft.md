# Abstract — what it must do, and a draft

> Working doc. Title is settled (`title_options.md`): *Information distortion in likelihood approximations for macroscopic ion-channel currents*. The title is plain and functional, so the abstract carries the hook. Comment and correct inline.

## The job the abstract has to do

Answer, in the first three sentences, the question raised in the very first audio of this paper (19 May): **what does this say that Communications Biology did not?** The bridge (MacroIR) is already published. If the abstract re-announces the method, a reviewer from the field will say so.

The new thing is three-fold and the abstract must land all three:
1. **bringing the misspecification test to this problem**, which is possible because the forward process is exactly simulable while its likelihood is not (the test itself is classical, see constraint 7),
2. the **information distortion matrix** that quantifies the failure, and
3. the **map**: which approximations distort, by how much, and where.

---

## The literature, verified (2026-07-14)

This changed the premise. An earlier draft opened with "several approximations are in routine use", which is **false**. Three independent checks (two agents reading the PDFs, one literature search) converged:

**The field fits the deterministic mean current by least squares.**
- Clerx, Beattie, Gavaghan & Mirams 2019, *Four Ways to Fit an Ion Channel Model* (Biophys J 117:2420): all four ways are deterministic ODE fits with an RMSE objective. Gating stochasticity is never mentioned.
- Owen & Mirams 2025, *IonBench* (PLoS Comput Biol 21:e1013319), the community benchmark: sum of squared errors on deterministic simulations, and explicitly *"We do not include optimisation approaches for stochastic models."*
- Del Core & Mirams 2025 (Phil Trans A 383:20240224): *"Despite many methods being available to estimate the parameters of deterministic models of ion channel gating, only a few have been developed for the stochastic case"*, and *"standard methods for deterministic models do not distinguish between stochastic channel gating and measurement error noise, resulting in biased estimates."*
- Wang et al. 2012 (PLoS ONE 7:e35208): a macroscopic-fitting paper that uses least squares on the mean, cites Milescu 2005, and **declines the likelihood on computational grounds**.

**Non-stationary fluctuation analysis is common but does a different job**: it returns the unitary current, the channel number and peak open probability. It is a two-moment regression, not a likelihood, and it does not identify a kinetic scheme. Stepanyuk et al. 2014: *"the unitary current is virtually the only parameter that can be reliably obtained from this type of analysis"*, and **"To the best of our knowledge, kinetic rates have never been estimated for any synaptic receptors in their intrinsic environment."**

**Likelihood-based stochastic fitting of macroscopic currents is a niche**, four or five groups over twenty years: Celentano & Hawkes 2004 (full covariance, 12 to 36 hours per fit), Milescu, Akk & Sachs 2005 (diagonal likelihood, in QuB), Moffatt 2007 (recursive filter), Stepanyuk et al. 2011/2014 (semiseparable covariance), Münch et al. 2022 (Bayesian Kalman), Moffatt & Pierdominici-Sottile 2025 (recursive interval filter, MacroIR), Del Core & Mirams 2025. Milescu 2005: *"exponential fitting remains the most widely used analytical method ... few procedures for direct estimation of rate constants are available."*

**All of them but one treat the observation as an instantaneous sample.** The exception is ours: MacroIR (Moffatt & Pierdominici-Sottile, Comm Biol 2025) integrates the observable over the acquisition window, and it is the only published likelihood for macroscopic currents that does. That is exactly why the acquisition window cannot be this paper's premise; the window is already handled, and opening on it re-announces the published method (see "The barrier this paper lifts" below). The rest of the literature either declares the window out of scope or works around it:
- Milescu 2005: *"Currently, the method does not handle the effects of decreased variance due to low-pass filtering."*
- Celentano & Hawkes 2004 mention the window only to dismiss it: *"Because the rise time of the filter (0.17 ms) is shorter that the sampling interval, no correction of the covariance matrix is made for the effects of the filter."*
- Münch 2022 files finite integration time under *"misspecifications of the likelihood"* to be robust to, measures the bias it causes, and prescribes an **experimental** workaround: *"the sampling should be faster than the fastest eigenvalues to avoid biased results."*

**Two questions have been open for twenty years, and both are in the field's own words.**
- **Milescu 2005**, on recursion: *"It remains an open question whether all estimates should be intrinsically biased if obtained by any method that is **not** a Bayesian filtering algorithm."* He also identifies his own dominant error source as *"the local time correlation of the current"*, which is exactly the correlation term of our distortion decomposition.
- **Moffatt 2007**, on averaging: *"the measurements are considered to be instantaneous, where electrophysiological measurements always result from some time-averaging ... **A formulation of the algorithms that can deal with time-averaged signals is therefore fundamental to the widespread use of the algorithms presented here.**"*

Neither was answered, because nobody had tested these likelihoods against the process they approximate. The test for it is classical statistics (Huber 1967, White 1982): under a correct likelihood the score vanishes at the true parameters and the covariance of the score equals the Fisher information the likelihood reports, and the information matrix test asks exactly that. It had simply never been brought to macroscopic-current likelihoods. What lets us bring it is the asymmetry of this problem, that the forward process can be simulated exactly while its likelihood cannot. That is the gap this paper fills.

**The algorithm family is the literature, not an invention.** `NR` (non-recursive, instantaneous) is essentially Milescu 2005; `R` (recursive, instantaneous) is Moffatt 2007 and Münch 2022; `IR` is MacroIR (Moffatt & Pierdominici-Sottile 2025). `MR` and `NMR` are the natural intermediates, and the paper should say so plainly. Celentano-Hawkes and Stepanyuk are the full-covariance non-recursive branch.

---

## The barrier this paper lifts, and the one it does not

Two barriers, and only one of them belongs to this paper. Getting this wrong is how the abstract ends up re-announcing Comm Biol.

**The time-averaging barrier is already down, and not by this paper.** Moffatt 2007 named it: *"A formulation of the algorithms that can deal with time-averaged signals is therefore fundamental to the widespread use of the algorithms presented here."* MacroIR lifted it, in Comm Biol 2025. Announcing that in the title, or in the abstract's opening, re-announces the published method, and a reviewer from the field will say so. This was the fatal flaw of the rejected "bridge" title, and it is the same trap.

**The barrier this paper lifts is epistemic.** Nobody had checked whether any of these likelihoods was faithful to the process, so nobody could know. With nothing to tell them apart, there was no rational basis for choosing among them, and none for paying the cost of the expensive one. The map supplies that basis. The barrier was never the absence of a test in statistics; it was that the test had not been run here, and running it needs a process you can simulate exactly and a likelihood you cannot compute exactly, which is what this problem gives you.

**And the argument is live, right now.** Del Core & Mirams 2025 argue explicitly against filtering, on cost: *"approaches that are based on filtering techniques are computationally expensive because an integration step of the differential moment equations (DMEs), and the corresponding updates in the correction step, are computed between every consecutive time points where the measurements are collected."* The field is deciding whether the filter is worth paying for, and it has nothing with which to decide. This paper is that evidence, and it lands directly on an argument published in 2025.

**Where this belongs: the last sentence of the abstract, and the Discussion. Not the title.** "The barrier is lifted" is a claim about adoption and impact that the paper does not demonstrate (it demonstrates a validity map on a two-state model in simulation), it re-centers attention on the method, and it is the promotional register that rule 5 below rules out. As the abstract's closing sentence it is honest and strong, because there it follows the evidence that earns it.

## Hard constraints

From the title review (four AI reviews plus a nine-agent panel). Accuracy constraints, not taste.

1. **Distortion, not loss.** Bidirectional: some approximations over-state the information the data carry, others under-state it. Never "information loss" or "preserving information".
2. **The approximation distorts, not the averaging.** Time averaging is the physical reality; the recursive interval likelihood handles it correctly. Blaming time-averaging inverts the causation.
3. **Continuous, not a threshold.** Avoid "breaks down at", "fails beyond".
4. **Scope is stated, not hidden.** Two-state model, non-stationary protocol, simulated macroscopic currents. Say it.
5. **No method promotion.** No "exact", "unbiased", "resolves a longstanding problem". That is Comm Biol's claim.
6. **Register.** Plain scientific prose, the register of the 2007 and 2025 papers. No aphorisms, no colon drama, no triads, no em-dashes.
7. **Do not claim the test.** The diagnostic is classical and well known in statistics: the score has mean zero at the truth, the covariance of the score equals the Fisher information under a correct specification, the gap between them is the sandwich, and testing that gap is White's information matrix test (Huber 1967, White 1982; the modern generalized form is Golden, Henley, White & Kashner). "There was no way to test whether a likelihood is faithful" is false as written and a reviewer with any statistical training will catch it. The claim that holds is narrower and still strong: this test had never been applied to macroscopic-current likelihoods, applying it here is possible because the forward process is exactly simulable while its likelihood is not, and what is new is the distortion matrix read directionally plus the map it produces. Cite the classical work in the Introduction and Methods, per `introduction_plan.md`: importing a known diagnostic into a field that had not used it is a strength.

## What stays out, and why

- **The Fisher-to-zero result** (no further information about the original channel number once the open population relaxes). Your own call in the 11 July audio: *"no creo que vaya el abstract porque es un poco demasiado."* It belongs in Results and Discussion. (Note: Del Core & Mirams 2025 declare exactly this problem still open, so it lands hard where it does go.)
- **The derivation of the evidence correction.** One clause here; derived in the later mechanism-discrimination study.
- **Micro, more than two states, the stationary regime, experimental data.** Out of scope; they belong in the Discussion's open-doors sentence.
- **The algorithm acronyms.** Describe them functionally in the abstract; name them in the paper.
- **MacroIR as a new method.** Named once, parenthetically, for searchability. Not presented.

## Draft

> Kinetic schemes are usually fitted to macroscopic ion-channel currents by least squares on the mean current, which discards the gating fluctuations that carry the information about the rates, the number of channels and the unitary conductance. Using those fluctuations requires a likelihood. The few that exist differ along two axes: whether they condition the hidden state on the data as it arrives, and whether they take each sample to be an instantaneous observation of the process or the average over the finite acquisition window that an electrophysiological measurement actually is. What either choice costs has never been measured, although the means to measure it are standard: a correct likelihood has a score that vanishes at the true parameters and a score covariance equal to the Fisher information it reports, and any gap between the two is misspecification. Testing that here is possible because the forward process can be simulated exactly while its likelihood cannot. In simulation we ask, for each approximation, whether the standardized residuals are white with unit variance, whether the score vanishes at the true parameters, and whether the covariance of the score matches the Fisher information the likelihood reports. The mismatch between the last two defines an information distortion matrix. In a minimal two-state channel model, and with no experimental data, we map this distortion across channel number, acquisition interval and instrumental noise. It runs in both directions and grows continuously toward the corners of the regime. The non-recursive likelihood inflates the information it reports, by a factor that grows with the number of samples per interval, so its confidence intervals are far too tight; the mid-interval recursion deflates it; the recursive interval likelihood (MacroIR) stays near calibration across the practical regime, degrading gradually where few channels or very short intervals carry the data away from the Gaussian limit. The same matrix supplies corrected estimates of parameter bias and variance, and propagates into Bayesian model comparison. The map turns the choice of likelihood from a matter of habit into a matter of evidence: it says when the cheap approximation will misreport the uncertainty, and when the cost of the recursive interval filter is worth paying.

~360 words. The closing sentence is the impact statement and it answers, directly, the 2025 argument that filtering is too expensive to be worth it.

### Length

**Settled, from the completed instructions (`docs/elife-author-instructions.md`, corrected 2026-07-14).** The live guide says **150 to 200 words**; the LaTeX class (`elife.cls` v1.11) still says no more than 150. **Write to 150 and treat 200 as the ceiling.** Subheadings are discouraged, and the register eLife asks for is "clear, measured, and concise", which is constraint 6 in different words.

**One rule that binds us specifically:** *"If the biological system is not in the title, it must be in the abstract."* Our title names macroscopic ion-channel currents but no organism, receptor or preparation, and the paper deliberately has none (a two-state model, in simulation). The abstract satisfies the rule by naming the system in the first sentence and stating the scope in the middle. Do not cut both.

The 360-word draft is over budget by a factor of two. Cut in this order:

1. The Bayesian-model-comparison clause. Keep the closing impact sentence; it is the "so what" and it answers a live 2025 argument. Cut the evidence clause before you cut the impact.
2. The justification clause (*"The forward process can be simulated exactly while its likelihood cannot"*), which the Introduction can carry.
2b. The "although the means to measure it are standard" clause. It can go from the abstract if space is tight, because "has never been measured" is true on its own and does not claim the test. What it must not do is come back as "there was no test" (constraint 7). If it is cut here, the Introduction has to carry the Huber/White attribution explicitly, and it does.
3. The final subordinate clause (*"degrading gradually where few channels…"*) — the scope is already stated.
4. Compress the premise to two sentences by merging the first two.

### Register: the structure eLife abstracts actually have

A current eLife abstract (single-molecule FRET on the ABC transporter TmrAB, 167 words) runs like this, and it is worth copying slot for slot. It is the summary-paragraph structure Nature prescribes:

| Slot | What it does | Words in the model |
|---|---|---|
| 1 | The object, and how it has been studied so far | ~25 |
| 2 | The limitation of that, in one line | ~15 |
| 3 | **"Here, we use X to do Y"** — the whole paper in one sentence | ~40 |
| 4 | Validation of the approach | ~25 |
| 5–7 | Results, one of them carrying a hard number | ~45 |
| 8 | **"These results provide the first … and establish …"** | ~30 |

What the model abstract does **not** contain, and what our 360-word draft is full of: literature critique, epistemology, and justification of why the approach is possible. No sentence in it argues that the field could not previously know something. The gap is one clause in slot 2, and then the paper is stated. Our draft spends three sentences on the epistemic setup, and that is the first thing to cut. The model also claims priority plainly ("the first single-molecule characterization"), which is licensed as long as the claim is the narrow one (constraint 7).

### Version A — at the ceiling (200 words)

> Kinetic schemes for ion channels are largely fitted to the mean macroscopic current, discarding the gating fluctuations that carry much of the information about the rates, the channel number and the unitary conductance. The likelihoods that use those fluctuations differ above all in whether they treat each sample as an instantaneous observation or as an average over the acquisition window, and what that choice costs has never been measured. Here, we measure it, because macroscopic currents can be simulated exactly although their likelihood cannot. For each approximation we test whether the score vanishes at the true parameters and whether its covariance matches the Fisher information the likelihood reports. The mismatch defines an information distortion matrix, which we map in a two-state channel model across channel number, acquisition interval and instrumental noise. The instantaneous approximations overstate the information the data carry by an order of magnitude and more, so their confidence intervals are far too tight, while the recursive interval likelihood (MacroIR) stays near calibration and errs conservatively where it degrades. The same matrix supplies corrected estimates of parameter bias and variance. These results assess macroscopic-current likelihoods against the process they approximate, and give a quantitative criterion for choosing among them.

### Version B — closer to the template's 150 (182 words)

> Kinetic schemes for ion channels are largely fitted to the mean macroscopic current, discarding the gating fluctuations that carry much of the information about the rates and the channel number. The likelihoods that use those fluctuations differ above all in whether they treat each sample as an instantaneous observation or as an average over the acquisition window, and what that choice costs has never been measured. Here, we measure it, because macroscopic currents can be simulated exactly although their likelihood cannot. For each approximation we test whether the score vanishes at the true parameters and whether its covariance matches the Fisher information reported by the likelihood. The mismatch defines an information distortion matrix, which we map in a two-state channel model across channel number, acquisition interval and instrumental noise. The instantaneous approximations overstate the information the data carry by an order of magnitude and more, so their confidence intervals are far too tight; the recursive interval likelihood (MacroIR) stays near calibration. The same matrix corrects the reported bias and variance, and shows when the cost of the interval filter is worth paying.

**What both versions give up, deliberately:** the residual whiteness test (slot 4 keeps only the two score tests, which are the ones the distortion matrix is built from), the Huber and White attribution (the Introduction must carry it explicitly, per constraint 7), the two-axis description of the likelihood family, and the evidence-propagation clause. Version B also drops the unitary conductance from the opening list and the closing "errs conservatively" clause.

**Two things to settle before either is final.**

- **[Q] The hard number.** The model abstract's slot 5 carries "∼300 ms". Ours says "by an order of magnitude and more", which is vaguer than it needs to be. `Figure_S3_caption.md` gives 10 to 16 for the non-recursive pair; `figure_2_pub` gives 14 to 21 for the same comparison on other cells. Settle which cells the abstract is quoting and put the number in, because that is the sentence a reader remembers.
- **[Q] "The first assessment".** Version A claims priority in the narrow form the model abstract uses. It is defensible (constraint 7 allows the narrow claim), but check it against the Introduction's phrasing so the two do not drift.

The directionality clause in both follows the corrected carrier in the Verify section (window-ignoring approximations inflate, the interval likelihood errs conservatively where it degrades), and not the earlier "the mid-interval recursion deflates it", which the figure inventory does not support.

---

## Impact Statement — a required deliverable we did not know about

eLife still requires one, and the eLife assessment did **not** replace it. It is a **submission-form field**, not part of the LaTeX file.

**Rules:** one sentence, **15–30 words**, **third person** (no "we", no "our"), states the single most important finding, **complements rather than repeats the title**, no unfamiliar acronyms. Do not open with "We show…", "This study…", "Our work…".

**Drafts [Q] — pick one:**

1. *Approximations that treat each sample as instantaneous misreport how tightly macroscopic currents constrain channel kinetics; a validity map shows where the error matters and where the recursive filter earns its cost.* (30 words)
2. *A simulation-based test reveals when likelihood approximations misreport parameter uncertainty in macroscopic currents, and when the cost of a recursive interval filter is worth paying.* (24 words)
3. *The information a macroscopic current carries about channel kinetics is systematically misreported by the likelihoods in use, in both directions, by amounts that a validity map now quantifies.* (28 words)

4. *Likelihood approximations in routine use for macroscopic ion-channel currents overstate the information the data carry by more than tenfold, and the resulting overconfidence can be measured and corrected.* (28 words)

Draft 1 states the finding and the payoff and stays clear of the title's wording. Draft 2 is the safest. Draft 3 foregrounds the bidirectionality, which depends on the unresolved MR sign below. **Draft 4 is the recommendation**, because it is the only one carrying a number, and a number is what an impact statement is read for. It depends on settling the tenfold figure (the same [Q] as the abstract's hard number), and on "in routine use" being defensible, which the literature section supports for the non-recursive family.

## Decisions taken (2026-07)

1. **Name MacroIR once, in parentheses.** Continuity and searchability with Comm Biol, without presenting the method.
2. **Directionality explicit.** Inflate for some, deflate for others, rather than the smoother "misstate".
3. **Say the work is in simulation, outright.** "In simulation", "with no experimental data".
4. **Article type: open.** eLife's *Tools and Resources* requires that "new methods must be benchmarked against existing methods", which is literally what this paper does. Research Article is the default because the contribution is a characterization. Decide at submission.

## Verify before submission

> **Update 2026-07-14, from the figure inventory in `results_plan.md`. Two corrections, and the second one changes the draft.**
>
> 1. **The overconfidence factor is 10 to 16, not 14 to 21.** `Figure_S3_caption.md` (the authoritative source): empirical/Fisher variance ratio ≈ 10 to 16 for NR and NMR, 1.2 to 1.4 for R, 1.5 to 2.1 for MR, and the sandwich-corrected ratio ≈ 1 in every cell. Fix the number below and anywhere else it appears.
> 2. **The bidirectionality claim has a different, and better, carrier.** Every non-IR member appears to *inflate*. What *deflates* is IR itself: `figure_5_distortion_algo_grid.Rmd` reports IR's own departures running in the under-confident (conservative) direction, down to ~0.5 in the low-N_ch corner. So the honest sentence is that the approximations which ignore the acquisition window over-state what the data know, and the one that respects it, where it fails, under-states it. That is bidirectional, it is safe (the conservative direction does not mislead a user), and it does not depend on the disputed MR sign. Rewrite the draft's directionality clause on this.
> 3. The apparent MR sign discrepancy below is probably a **category error** that dissolves: MR overestimates the *predicted observable* variance (it drops the boundary cross-covariance) while still *under*-reporting the *parameter* covariance (over-confident, 1.5 to 2.1). Different variances; both can be true. Verify, then fix the master plan's ambiguous "Variance distortion" column header.

**The deflation side of the directionality claim is not yet confirmed.** The inflation side is solid: the non-recursive approximations are over-confident by roughly 14 to 21 fold, growing with samples per interval. The deflation side rests on MacroMR **overestimating** the variance (which would make it under-confident), and there is an unresolved sign discrepancy in the record: the ranking in `00_master_plan_v2.md` says MR overestimates variance, while the `Figure_S3` caption reports MR as **over**-confident (empirical/Fisher 1.5 to 2.1). One of the two has the sign backwards.

**Check the MR sign against the Gaussian-Fisher data before this sentence is submitted.** If MR is over-confident like the non-recursive ones, the bidirectionality claim needs a different carrier, or the honest statement becomes the distortion matrix's own anisotropy (a single algorithm inflating in some parameter directions and deflating in others).
