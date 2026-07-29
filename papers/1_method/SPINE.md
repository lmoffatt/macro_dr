# SPINE — one block per manuscript section

> **Created 2026-07-29**, consolidating the eight section plans (`abstract.md`, `introduction.md`,
> `theory.md`, `diagnostics.md`, `results.md`, `discussion.md`, `methods.md`, `title.md`), which are
> archived under `archive/section-plans-20260729/`. This is the file `elife_paper.tex:48` already
> cites as "SPINE".
>
> **What this file holds, and nothing else:** per section, the *job*, the *constraints* that govern
> the writing, what is *open or blocked*, and the *verify-before-submission* list. The prose lives in
> `docs/manuscript-drafts/sections/*.tex`, and the provenance of every number lives beside it as a
> `% src:` comment (144 of them at last count). **If a fact is in the `.tex`, this file cites it and
> does not restate it.** The one thing here that looks like content is the list of claims the data do
> not support, and that is the negative space, which a manuscript cannot contain by construction.
>
> **Written from the manuscript, not from the old plans.** The `.tex` was rewritten 2026-07-28 in the
> merged frame and is ahead of everything else in `papers/`; the archived plans are pre-merge and
> several of their scope sentences were false when harvested. Where a harvested item has been overtaken
> it is marked **[SETTLED SINCE]** rather than deleted.
>
> Roster, figure set and every cross-section decision: `decisions.md`. Cross-paper: `../_program/`.

---

## Rules that govern every section

Accuracy constraints, not taste. They came out of a four-review plus nine-judge convergence on the
title and they bind the whole manuscript.

1. **Distortion, not loss.** The effect is bidirectional: some approximations over-state the
   information the data carry, others under-state it, and `IR` itself runs to both sides in the
   few-channel corner. Any phrasing built on "information loss" or "preserving information" is
   factually wrong.
2. **The approximation distorts, not the time averaging.** Time averaging is the physical reality that
   the boundary-conditioned likelihood handles correctly. "Time averaging degrades the information"
   inverts the causation and contradicts the thesis.
3. **Continuous, not a phase transition.** Avoid "breaks down at", "fails beyond". **[SETTLED SINCE]**
   This once forbade the region map outright; it does not. The reconciliation is the map's own
   sentence, and it should be quoted rather than re-derived: *the map is a concept map, not a phase
   diagram; its boundaries are level sets of continuous diagnostics, so a looser criterion moves each
   of them by up to a decade in noise without changing the layout.*
4. **No "only", no universal claims.** The study is two-state. "Only X stays calibrated" invites a
   scope objection that the data cannot answer.
5. **No method promotion.** The paper characterises; it does not re-present MacroIR. That bridge is
   published (Comm Biol 2025). Every sentence beginning "MacroIR is…" is a candidate for deletion.
6. **Do not claim the test.** The diagnostic is classical: Huber 1967, White 1982, and the generalized
   form of Golden, Henley, White & Kashner. "There was no way to test whether a likelihood is faithful"
   is false and a statistically literate referee will say so. What is new is the *measurement* of it
   for this class of likelihood, enabled by a process that can be simulated exactly.
7. **The gap is validity, never absence.** The temporal correlation of macroscopic currents has carried
   kinetics since 1973. Do not write that it was unused. See `decisions.md` and
   `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`.
8. **Register.** Plain scientific prose, the register of the 2007 and 2025 papers. No aphorisms, no
   rule-of-three triads, no em-dashes.

**Two readers, and every section has to hold both.** The *electrophysiologist* who fits schemes to
macroscopic currents, has never computed a score, and must not need to know what a sandwich estimator
is; and the *inference-methods reader* (the Mirams / Münch / Del Core axis), for whom the words that
buy attention are misspecification, score, and information matrix equality. One sentence each, no more,
and gloss each in physical terms in the same sentence.

---

## Title

**Chosen, for now:** *Information distortion in likelihood approximations for macroscopic ion-channel
currents.* "for", not "of": what is approximated is the likelihood, not the currents.

**Currently in the manuscript** (`elife_paper.tex:32`): *A validity criterion and a usage map for
likelihoods of macroscopic ion-channel currents.* **[OPEN]** The two disagree and the `.tex` one is
post-merge, so it is probably the survivor; decide and record which, because rule 5 above and the
abstract both key off the title.

eLife rule that binds: if the biological system is not in the title, it must be in the abstract. It is
currently in both.

---

## Abstract — `elife_paper.tex:52-58`

**The job.** Make two words easy for an editor to assign. Under eLife's model there is no private
accept decision; there is a permanent public assessment written with a controlled vocabulary, and the
abstract is where it comes from. Significance runs landmark > fundamental > important > valuable >
useful; strength of evidence runs exceptional > compelling > convincing > solid > incomplete >
inadequate. Target: **compelling** on evidence, earned by the ground truth being an exact simulation of
the process the likelihood approximates; and **valuable to important** on significance, earned by
saying who can act on the result.

**Constraints.**
- 200 to 220 words. **[SETTLED SINCE]** The old rule was "write to 150, treat 200 as the ceiling".
  Measured against 41 recent eLife research articles: median 197, 39% over 200, and the
  computational/structural/biophysics subset has a median of **210**. Do not spend words getting to 200.
- Name the biological system in the first sentence.
- **The paraphrase test.** Every technical term must survive being restated by an editor who is not a
  statistician, because that restatement is published and permanent. If an editor cannot write the
  assessment without the word "score", gloss it.
- **[SETTLED SINCE]** The old slot table was Nature's funnel. Measured, eLife does not use it: only
  half of abstracts contain "Here we", and among those, method verbs beat result verbs 13 to 8; only
  27% carry an explicit comparison to prior belief. What *is* the eLife norm, and what this abstract
  already has, is the broader-perspective closer (68%).

**Impact Statement** (submission-form field, not part of the PDF): one sentence, 15 to 30 words, third
person, complements the title, no "We show". **[OPEN]** three drafts are in the archived `abstract.md`
and none is post-merge.

**Open.**
- The significance sentence. The current close, "this supplies a test… and a map…", describes the
  deliverable rather than what changes. The formulation to work from (Luciano, 2026-07-29): *a valid
  likelihood is what enables separating inter-experiment variability from stochastic error, and the
  evaluation of Bayes factors.* Order them by evidence, and make the Bayes half concrete via the
  distortion matrix rather than generic. See `decisions.md`.
- **The NSFA collision.** The map's largest region is sold on the unitary current and the channel
  count, and `01_introduction.tex` concedes both to non-stationary fluctuation analysis while arguing
  NSFA returns no kinetic scheme. Resolve before drafting: what is new is the *combination*, a
  calibrated joint estimate of rates and amplitudes from a single non-stationary record, where the
  1980-81 covariance methods needed ensembles of 256 to 504 repeated sweeps and stated no uncertainty.

**Verify.** Every quoted number at a fixed n_sims. The direction of the non-recursive result, which
inherits the unverified sign convention (`% TODO-SIGN` in the `.tex`). No channel-number threshold
unless D-J is closed — and note the manuscript already prefers the measured reachability floors
(N_ch below 17 for the channel count, below 52 for the opening rate) over "about 100".

---

## Introduction — `sections/01_introduction.tex`

**The job.** Take a reader who fits the mean by least squares and show them, in their own vocabulary,
what the fluctuations carry, why using them takes a likelihood, and why no one can currently tell
whether such a likelihood is working.

**The five moves** (2026-07-28 voice notes): independence as the foundational assumption of least
squares, and what violating it costs; Markov chains as the way to model dependence while still yielding
record-universal constants tied to biophysical structure; the ladder of methods; why least squares still
reigns (you can see, point by point, whether it predicts, which is visually convincing, while a
recursive filter follows the data closely and holds no surprises, so a working method and a bug look
alike); and why autoregressive alternatives do not solve it.

**Constraints.**
- **Do not re-announce MacroIR.** Cite it as prior work in the same breath as Milescu 2005 and
  Münch 2022. The moment the Introduction explains how it works, the referee says "you published this".
- **The NSFA naming trap.** Non-stationary fluctuation analysis *is* widely used and *does* use the
  variance. Name it and distinguish it or an electrophysiologist will think the gap is already filled.
  Stepanyuk 2014's own words end the objection: *"the unitary current is virtually the only parameter
  that can be reliably obtained from this type of analysis"*, and *"kinetic rates have never been
  estimated for any synaptic receptors in their intrinsic environment"*.
- **On ARIMA**, four points, all derivable from machinery already in the paper: an ARMA error model is
  stationary by construction while the gating covariance tracks the mean current and restarts at every
  jump; it contains no channel count and no unitary current; its timescales are free where the Markov
  model ties them to the rates that generate the mean; and whiteness is not falsification, since enough
  ARMA terms whiten anything.

**Stays out.** The mechanics of MacroIR (that is Theory). The Fisher-to-zero result. The
research-program framing, which reads as a grant proposal and belongs in the Discussion. The Comm Biol
biology beyond one citation. Any claim about experimental data.

**Open.** How much statistics vocabulary in the "why nobody could test this" move: recommendation is
all three terms, each glossed in physical terms in the same sentence, at a cost of about forty words.

**Verify.** Del Core & Mirams 2025 and Owen & Mirams 2025 quotations against the version of record.
**[SETTLED SINCE]** The old top verify item was the "no published likelihood integrates the acquisition
window" claim; it is dead as written and survives only scoped to macroscopic-N by a scaling argument
(`decisions/D-3_novelty_claim.md`).

---

## Theory — `sections/02_theory.tex`

**The job.** Make the methods **one object graded by how much of the interval structure the likelihood
uses**, so the Results read as a traverse of a coordinate system rather than a bake-off among unrelated
codes.

**The one-sentence spine.** A macroscopic current is the sum over a population of channels, each a
continuous-time Markov chain; the exact likelihood of a time-averaged recording would need the
distribution of the interval-averaged conductance of the whole population, which is intractable, so
every practical likelihood replaces it with a Gaussian, and the methods differ in *what they condition
that Gaussian on* and *how they account its variance*.

**[SETTLED SINCE]** The archived plan scoped this section to `R, MR, VR, IR` and handed least squares
and the non-recursive members to "paper 2". There is no paper 2. The section covers the root question
and the ladder together (`decisions.md`).

**Stays out.** k-state generality beyond what two states need. The Taylor variants, and keep them clear
of `VR`: the engine flag is `taylor_variance_correction`, `VR` is not a Taylor variant, and Methods
must say so or a reader will conflate them. The PSD trust coefficients. The Kalman prior-art
connection, which reads defensively here and belongs in the Discussion. The distortion machinery
itself, which is `../_program/machinery.md` and is presented in Diagnostics.

**Open.** How much algebra survives in the main text. Proposal: four displayed equations, the
observable as an interval average, the boundary-conditioned conductance, the total-versus-residual
variance difference (the MR→VR step) and the boundary cross-covariance (the VR→IR step). Test: could an
electrophysiologist read the Results with only those four?

**Verify.** The emission-variance decomposition against the code, not the archived drafts, which
predate the `gvar_i` fix. The total-versus-residual variance forms and the boundary cross-covariance
term, re-derived from current code before printing: `VR`'s definition *is* the residual form, so if the
code's `gvar_i` is not what Theory claims, the roster is wrong and not just the prose.

---

## Diagnostics — `sections/03_diagnostics.tex`

**The job.** Convert "we tested whether the likelihood is faithful" from a claim into a procedure a
sceptic could run on their own likelihood tomorrow. It is the part a reader might reuse, and it has to
survive a statistician, which means being explicit about which identities are exact algebra, which are
first-order, and which are conventions.

**The three tests, in the wording that wins.** `../_program/machinery.md` §2 owns them and its wording
beats the voice-note paraphrase on two points that matter: test 1 includes **residual whiteness**, which
is the sharpest discriminator on record (score ACF: IR ≈ 0.005 against 0.78 for the non-recursive
members), and test 3 anchors on **the model's own Gaussian Fisher**, not on "the Fisher as a proxy for
the Hessian". Do not overwrite machinery with the looser phrasing.

**The anchor, and the objection to pre-empt.** The anchor is the model's own Gaussian Fisher, computed
from its predictive mean and variance and their derivatives at no cost beyond the score. For a macro
algorithm the likelihood is Gaussian by construction, so this is the model's own Fisher and the
sandwich is exactly the information-matrix-equality test in the model's own frame. It is PSD by
construction, so the diagnostic never has to be rescued from an indefinite anchor; and the numerical
finite-difference Fisher is computed separately only to gauge how faithful the Gaussian one is, not as
the anchor, because in this regime it is widely indefinite per replicate. **State that reason.** A
referee who knows the sandwich literature will otherwise ask why the Hessian was not used.

**The circularity answer, which must be in the paper and not just in our heads.** The Gaussian Fisher
is the information the likelihood *claims*; the score covariance is what it *delivers* under data from
the true process. They are different objects and their mismatch is what misspecification means. The
anchor being internal is not a circularity, it is the definition.

**Stays out.** The evidence correction beyond one motivating paragraph. The entire posterior-distortion
family. The SAFE/MARGINAL/UNRELIABLE/BIASED categorization, which gates on the likelihood Hessian's
spectrum and belongs to the posterior framework. The implementation trust region beyond a Methods
paragraph.

**Verify.** The `C = C_s^{1/2} R C_s^{1/2}` identity: prove it under a stated condition or demote it to
the determinant version, which is exact. Do not print it as written. n_sims uniformity across every
cell of every heatmap. **The direction convention, once, against the producer** — see the standing
blocker below.

---

## Results — `sections/04_results.tex`

**The job.** Walk the reader from one filter step to the usage map without asking them to take a claim
on trust. Each figure answers the objection the previous one raises.

**The through-line, six body figures** (`decisions.md`, "The figure set"): the filter step along the
cost ladder; recovery clouds at one cell; the calibration cascade in time; the design space split by
moment; the information budget per parameter; and the usage map. **[SETTLED SINCE]** The archived plan
had five figures and an `R`-versus-`IR` arc.

**Framing the maps force.** The three named regimes (multinomial, telegraphic, Gaussian) are the two
approximations' **asymptotic corners, not territories with borders**. The map shows a gradient; the
corners explain its direction. Write that into both Theory and Results so they do not read as two
different claims.

**The threshold is 1.15, read as a variance ratio**, which is a 7% error on the reported standard
deviation. **[SETTLED SINCE]** The repo carried 1.05, 1.1, ±15% and 1.5 in four places; the built map
chose 1.15 and glossed it, and `../_program/machinery.md` §8 must now cite that rather than propose
1.1.

**Claims the data do NOT currently support. Say none of these until the gap is closed.**
- *"Sample × correlation = total distortion."* Asserted, never plotted.
- *The N\* ∝ noise law as a validated prediction.* N\* is measured but never overlaid on its predicted
  law, and the fitted exponent is 0.80 with R² = 0.59, which is neither a clean 1/2 nor 1. Report the
  measurement; do not present the law as confirmed.
- *The micro-versus-macro comparison.* Blocked on data and out of scope; it must not sneak in through a
  figure.
- *Anything from the flagged-wrong scripts*: the argmax ridge (it is noise), the `relCI_by_Nch` panel
  that contradicts its own sibling, and the schematic that hard-codes a boundary contradicting the
  paper's own N\* result.
- **[ADDED 2026-07-29]** *Separating inter-experiment variability from stochastic error.* This is the
  paper's best significance claim and it is currently a capability statement, not a measurement. It is
  cheap to demonstrate from data already on disk: every run is at fixed θ, so Var(θ_true) = 0 is known
  by construction and the between-experiment spread is entirely stochastic; inject a known variability
  onto the per-group MLEs and ask each method how much real variation it infers. Until that is run, the
  claim is an implication and must be worded as one.

**Verify.** The ±15% self-consistency fractions on whatever anchor the final figures use. The
overconfidence gap: confirm it is the recursive/non-recursive split and not a recursive/averaging split,
because the sentence reporting it will be read as a claim about *which* approximation matters.

---

## Discussion — `sections/05_discussion.tex`

**The job.** Say what it means for someone about to fit a macroscopic recording next week, and close the
loops the Introduction opened. Three deliverables, in order of importance: **a decision rule** (not
"MacroIR is best" but the conditions, which is what gets cited); **the two answers**, to Milescu 2005 on
whether non-filtering estimates are intrinsically biased and to Moffatt 2007 on the need for a
time-averaged formulation, both named explicitly; and **the honest perimeter**.

**The trap.** This is where a paper drifts back into being about its method. The subject of this
Discussion is *the distortion*; MacroIR is where the distortion happens to be small. Keep the grammar
that way and the paper stays what it claims to be.

**The three points** (2026-07-28 voice notes): MacroIR is calibrated over almost the whole measured
plane, so one could simply always use it; where it fails, which is few channels and telegraph
(non-Gaussian) noise; and why the intermediates fail, which is that the double conditioning at both
interval ends is what is load-bearing and conditioning on the start alone contributes nothing.

**Open, and both are strategic rather than editorial.**
- **The Comm Biol coupling.** The published P2X2 evidence ranking was run with the buggy `gvar_i`. This
  paper's second significance claim leans on that ranking as its only real-data demonstration, so the
  two are linked in print. Options: say nothing and handle the erratum separately; one neutral sentence
  that the numbers are being re-checked; or the strong version, in which this paper is the machinery
  that would have caught it and shows it working on the author's own published result. Decide before
  drafting, not after.
- **The 6.4 margin.** Comm Biol's winner beats its symmetric counterpart by a Bayes factor above 5000,
  which is immune to almost anything, but beats the top non-conformational alternative by only 6.4
  (5.5-6.7). That is the scale a systematic likelihood distortion can move, and this paper now has the
  map that says how much distortion that regime carried. Computing it is a strong result either way,
  and better found by us than by a referee.
- Does the research-program framing get a paragraph? Recommendation: two or three sentences, not a
  section. The paper must stand without it.

**Verify.** Every number in the decision-rule table against the freeze recompute
(`decisions/D-4_ranking_verdict.md`), the direction of the MR result in particular.
**[SETTLED SINCE]** The old verify list asked whether NMR is unbiased, on which a whole paragraph
rested; NMR is now a supplement member and measurably indistinguishable from NR, so that paragraph
needs rewriting rather than checking.

---

## Methods — `sections/06_methods.tex`

**The job.** Two readers: a referee checking the numbers could have come out this way, and someone in
2029 trying to re-run it. Write to the second. **This paper's claims are about algorithms, so the
Methods *are* the result**; a vague Methods section does not inconvenience the reader, it makes the
paper unfalsifiable.

**A naming correction the manuscript must not inherit.** The intervals in the CSVs come from a
**percentile bootstrap over groups**, not from a probit transform, despite every column being named
`probit_*`. Write: nonparametric bootstrap over groups, with percentile intervals at 2.5, 16, 50, 84 and
97.5%. Counts: 100 replicates for the diagnostic battery with `max_lag = 10`; 200 for the empirical
Gaussian distortion; `min_groups_for_bootstrap = 10`, below which interval slots are NaN-filled, which
is one source of grey cells in the maps. `Model_Parameters_Hat` is **not** bootstrapped by the program;
the raw MLE cloud is what is on disk and any bootstrap of it happens downstream in R.

**The reproducibility problem, and it is real.** There is no run manifest for either production run.
The dispatcher defaults do not match the data: defaults say n_sims 256 and group sizes {1, 10, 100},
the data says 10,000 and {10, 100}. The actual overrides were reconstructed from filenames and CSV
columns and are recorded nowhere. **Write the manifest from that reconstruction and commit it beside
the data.** Then state plainly in Data Availability that output directories are keyed by the engine's
baked git hash, the engine is pinned by tag, and each figure script names its data folder — that part
of the pipeline is genuinely good and should be described. Environment for the record: SLURM, 32 CPUs
per task, 48 GB (96 for the figure-2 dispatch), BLAS threads pinned to 1 with OpenMP taking the
simulation and bootstrap loops, and `MACRODR_AXIS_SERIAL=1`, which is load-bearing because without it
memory grows with concurrent combinations and the jobs OOM.

**The anchors, which Methods has to own.** Map panels come from the Gaussian-Fisher runs; the
time-resolved figures are a separate analytic per-step computation on the dumps, not the numerical
battery. State the anchor per figure. **[SETTLED SINCE]** `433ed13` is the numerical-Fisher demo and
**no paper number is quoted from it**; the verdict was recomputed on `1c2ae6f` + `87889e6` + `0ffbda7`
and survives, which is exactly what `433ed13` was run to show. Say that once, here.

**Numerical safeguards.** One paragraph, because the code runs it and because a reader implementing the
filter hits the simplex problem on day one. The paper needs the *what*, not the how-we-got-here.

**Owed.** The recording-condition table for the usage map: roughly 200 lines of commented derivation
for excised patch, whole cell and oocyte currently live inside `figure_6.Rmd`, backed by 88 sourced
files. Lift it into a Methods table, one source per row. It is the part a referee checks first.

---

## Back matter — `sections/07_backmatter.tex`

**The R and Python packages.** Already written, with tolerances, and honestly scoped: the core was
checked against the engine on eight cells and agrees on the total log-likelihood to between 10⁻¹⁶ and
6×10⁻¹⁵ relative; the two language packages are bindings over that one core and were run against it
**once, on a single cell**. That last sentence is the weakest in the submission — a usability claim
inside a paper about when methods are valid, backed by one cross-check. Run the check across the design
grid and record the tolerance per cell, or keep the sentence exactly as honest as it currently is.

---

## Standing blockers, which belong to no single section

1. **The direction convention.** `../_program/machinery.md` §4 asserts it three times with two different
   signs and it has never been checked against the producer. Every verdict and every boundary of the
   usage map inherits it, and the manuscript carries eight `% TODO-SIGN` markers waiting on it. It is
   about an hour of reading code and it should be done before any result sentence is finalised.
2. **`Figure_2.pdf` is not rendered in `paper_both`.** The notebook and the caption are there and the
   caption was updated 2026-07-23, but the PDF exists only in `paper_1` from before the merge. It is the
   only one of the manuscript's 23 referenced graphics that is missing, and the `.tex` already tracks it
   as `PENDING-FIG2`. One R render.
3. **LSE at n_sims 1000** against the likelihood arm at 10⁴. Mixing them in a panel manufactures a
   regime effect through the Jensen bias. Re-run, or restrict every LSE panel to a uniform n_sims and
   say so in the caption.
4. **Telegraph noise.** The Discussion asserts MacroIR fails there. Theory names the regime, but no
   simulator capability and no run exists on any freeze commit. Either measure it or demote it to an
   argument with a citation, and say which.
5. **`check.sh` tests `N_CAP >= 5`.** The figure set is six. Until that is raised the done-oracle stays
   green on a manuscript short one body figure.
