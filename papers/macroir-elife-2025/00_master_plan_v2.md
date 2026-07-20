# Paper 2 (MacroIR / eLife 2025): Master Plan (v2)

> Updated: 2026-07-20. Draft rewritten from the MacroIR 13 audio consensus (2026-07-11), then **reframed to the usage map** on 2026-07-20 (audios of 10:54–10:59: LSE as method zero, noise axis extended past the gating crossover). Supersedes `00_master_plan.md` (now retired to a stub; prior content in git history). Comment and correct inline; open items for the author are marked **[Q]**.
>
> **Propagation owed by this edit** (`00_master_list.md` §6): `02_decision_log.md` (thesis + roster + ranking), `05_experiment_grid.md` (axes), `nomenclature.md` and `theory_plan.md` (the two-knob framing is now false), `methods_plan.md` (three flags, not two), `results_plan.md` (the arc), `introduction_plan.md` and `abstract_draft.md` (LSE moves from cited background to measured arm), and every caption that says "five approximations".
>
> This paper is **bridge 2** of the research program in `From molecular mechanisms to data back and forth PROGRAM.md`: the robustness study of the likelihood algorithm. That file carries the framing (theory → model → likelihood → evidence → ranking); this file is the operational plan for this paper. Repo split plan: `09_carve_plan.md`.

## The argument

If we had an exact method, we could just demonstrate its validity mathematically.

On the other hand, if we have an approximation, a method we know is not exact, what do we do?

That's the question I was wondering myself.

The response I arrived at: find the situations where the method is very precise and the situations where it is not, and quantify both.

So the idea is to prove that the method works where it is supposed to work and breaks where it is supposed to break.

So, what is the method trying to solve?

The problem consists in approximating the likelihood function of time-averaged observations of a given Markovian process. We can simulate the unfolding of this Markovian process and its observations with very good accuracy; the inverse, going from an observed experiment to a likelihood, is the hard direction, and that asymmetry is exactly what lets us use the simulation as ground truth to judge any proposed approximation.

We have alternative approximations of the likelihood function and we want to see if they are good. After some introspection and reading the literature, we conclude that we can test whether a likelihood function faithfully represents a simulation process by a three-step analysis:

1. the normalized residuals should have mean zero, variance one, and no temporal autocorrelation;
2. the expected score should be zero, evaluated at the true parameters;
3. the variance of the score should equal the negative Hessian of the log-likelihood, evaluated at the optimum where the score is zero.

For a misspecified approximation, tests 2 and 3 land at different points, and that gap is the whole story of the bias.

Then we accept that we are dealing not with true likelihoods but with approximations, and we ask whether we can find a procedure to improve the approximation. We found such a procedure: an expression for the corrected bias of the MLE parameters and one for their corrected variance. Near the maximum, the same distortion propagates to the Bayesian evidence through two channels, a volume correction ½ log det C and an effective-sample rescaling α⋆ = p / tr C, so the correction is what keeps model comparison valid under an approximate likelihood. (Stated here as motivation; the derivation is deferred to the program's bridge-3 study.)

Now, what is the nature of the main algorithm we are presenting? The algorithm represents the predictive (prior) state of the channel population, the density that assigns a probability to each combination of channel states, by a multivariate Gaussian. That density is exactly multinomial when channels are independent (equivalently, the maximum-entropy closure given only the mean occupancies), and the Gaussian is its large-N limit. The higher moments the closure discards do not vanish; they reappear later as the correlation distortion the diagnostics measure.

## 0) Goal in one sentence

Produce a **usage map for the methods actually available to fit macroscopic ion-channel currents**, from least squares on the mean current through the five Gaussian likelihood approximations, saying for each regime of channel number, acquisition interval and instrumental noise which is the cheapest method that still reports its own uncertainty honestly, in a minimal two-state channel model, in simulation.

The reframe from v2 (2026-07-20): the deliverable is a **usage map**, not a ranking. Two changes force it.

1. **Least squares enters as method zero.** Classical nonlinear least squares on the mean current (`nonlinearsqr`, displayed `LSE`; Levenberg-Marquardt, the method of Moffatt & Hume 2007 JGP) is what the field overwhelmingly does. It was previously cited as background and never measured. Measuring it is the point: the reader of this paper is someone who uses it.
2. **The noise axis is extended past the gating-noise crossover.** Every cell run before 2026-07-20 sat at or below the single-channel noise scale, which is the regime that favours the gating-aware likelihoods by construction. Showing IR wins there is close to tautological. The map now spans the regime where instrumental noise exceeds the gating noise, where the cheap methods are expected to be adequate, and the paper says so.

Together these turn "MacroIR is the only survivor" into "here is when you need to pay for MacroIR and here is when you do not", which is the claim that gets cited and the one that answers, on its own terms, the 2025 argument that filtering costs too much to be worth it.

The v1→v2 reframe still holds underneath: this is not "here is MacroIR", it is a validation machinery plus the map it produces, MacroIR's own failures included.

## 1) Core conceptual frame

**The structure is a root question and a ladder hanging from it, not a single ladder.** Earlier versions of this plan described the methods as "one object with two knobs", which was true of the five likelihoods and is false once least squares is in the comparison: in the dispatcher LSE carries the same two knob settings as NMR (`recursive=false, averaging=1`) and is distinguished only by a third flag, `family_approximation=2`. It has no rung. So:

- **Root question: do you model the gating fluctuations at all?** LSE answers no. It fits the deterministic mean current and absorbs everything else into one fitted noise scale. Every other method answers yes.
- **Given yes, what is the Gaussian conditioned on?** That is the two-knob lattice below, whose five members are ordered by the endpoint ladder, and where MacroIR is the top rung below intractability.

This preserves what the ladder was doing for the argument. The verdict stays structural rather than empirical, because both levels are predicted before any run: the root question is decided by whether gating fluctuations are visible above the instrumental noise (§1a), and the rung is decided by the ladder. The Results confirm a structure; they do not report a winner.

Every algorithm **below the root**, that is the five likelihoods, makes two Gaussian approximations:

1. **Macro approximation.** The multinomial distribution of channel occupancies is replaced by a multivariate Gaussian. Valid for large N_ch by the central limit theorem; degrades for few channels (multinomial regime).
2. **Interval-likelihood approximation.** The conductance distribution over one measurement interval is replaced by a Gaussian, when the exact object is a complex stochastic-telegraph / Poisson-like distribution. Degrades for intervals much shorter than the relaxation time (telegraphic regime).

This yields **three regimes**: multinomial (few channels), telegraphic/Poissonian (very short intervals), and Gaussian (many channels, moderate intervals, enough instrumental noise), where MacroIR is ideal. Because it is an approximation, MacroIR must fail somewhere; the paper locates that failure and shows it is the predicted degradation of these two approximations, not a coding artifact.

### 1a) The noise axis has two reference scales, and they cut it into three bands

The instrumental noise is not one axis with one landmark. It has **two natural reference scales**, and which one you compare it against depends on the method:

- the **single-channel** scale, the gating variance contributed by one channel, which does not depend on N_ch;
- the **gating** scale, the total gating variance of the population, which is N_ch times the first.

The two landmarks are therefore separated by exactly a factor N_ch, and that separation is what makes the map two-dimensional rather than one.

Units, following the settled convention in `decisions/D-2_parameter_units.md` (dimensionless ν = `Current_Noise`·k_off/g², figure label = 10·ν, so the single-channel crossover ν = 1 sits at **label 10**). Across the production interval range (Δ ≤ τ, so the within-interval averaging does not yet suppress the gating term) the two boundaries are:

| Boundary | In label units | Meaning |
|---|---|---|
| A/B | `label = 10 · interval` | instrumental noise reaches the single-channel gating scale |
| B/C | `label = 10 · N_ch · interval` | instrumental noise reaches the total gating noise |

which gives three bands and a predicted answer in each:

- **Band A**, `label < 10·interval`. Gating fluctuations dominate and are resolved finely. Root question: yes. Rung: the top of the ladder, IR.
- **Band B**, between the two. Gating still carries the information but is coarsely resolved. Root question: yes. Rung: R may suffice.
- **Band C**, `label > 10·N_ch·interval`. Instrumental noise swamps the gating noise; the residual scatter is white and the fitted noise scale absorbs it. Root question: **no**. LSE (and the non-recursive members) should be adequate, and paying for the interval filter buys nothing.

**These are predictions, made before the runs, from the two reference scales alone.** That is the strongest position the paper can be in, and it is why the noise axis is the map's organizing axis rather than a third nuisance dimension.

**Coverage, honestly stated.** Every cell produced before 2026-07-20 used labels {0.05 … 10}. At the headline cell (N_ch 100, Δ = 0.1 τ) that is band A throughout; at short intervals the old grid does reach band B, and the single corner N_ch 10 / Δ = 0.01 τ / label 10 grazes the B/C boundary. What was systematically missing is band C at the canonical channel numbers and intervals. The 2026-07-20 dispatch fills exactly that: `N_NOISE` is paired index-wise with `NCHS` as (10, 100), (100, 1000), (1000, 10000), (10000, 100000), which is `label = 10·N_ch`, the B/C boundary at Δ = τ. The interval axis inside each cell then carries it from the boundary down into band C.

## 2) Scope

### In scope

- Minimal **two-state** channel model (closed-open, `scheme_CO`), single K_on / K_off.
- **Macroscopic** interval-averaged currents.
- **Six methods on two levels** (§1):
  - **method zero**, off the lattice: classical nonlinear least squares on the mean current, data key `nonlinearsqr`, display label `LSE`, engine flag `family_approximation = 2`, noise scale marginalized (Jeffreys). Design and seams: `theory/macroir/notes/nonlinearsqr_lse_plan.md`.
  - **the five Gaussian likelihood approximations**: `NR`, `NMR`, `R`, `MR`, `IR` (script naming has used `MNR`; standardize on `NMR`).
- **Non-stationary** protocol (single concentration jump).
- Three control variables: N_ch, interval Δ relative to the kinetic relaxation τ, and instrumental noise **swept through the gating-noise crossover** (§1a), which is the change that makes the map a map.
- **Likelihood-only** analysis. MLE / Gauss-Newton local maximum kept, only to obtain the empirical parameter covariance.
- **Gaussian Fisher** as the distortion anchor; numerical finite-difference Fisher kept only to gauge how good the Gaussian one is.
- The **evidence-correction payoff** (how the likelihood distortion shifts the Bayesian evidence) stated as motivation; derivation deferred (see §3).

### Out of scope (later components of the program, not exclusions; see PROGRAM.md)

- MicroIR / microscopic-recursive (the multinomial-regime solver).
- More than two states / allosteric schemes.
- Stationary regime.
- Experimental data (biophysical relevance already established in the P2X2 Comm Biol paper).
- The **posterior information-distortion framework** and full **model-comparison results** (the likelihood-side evidence correction stays in scope as motivation, above and §3).
- Taylor variance-correction variants (IRT, MRT).

## 3) Diagnostics (the validation machinery)

Four indicators, all at the likelihood level, evaluated at the global MLE / optimum:

1. **Standardized residual** mean and variance (R²): should be 0 and 1, with no temporal autocorrelation.
2. **Score bias**: E[score] at the truth, projected through the parameter covariance into a bias vector.
3. **Information distortion matrix**: Var[score] against the Gaussian Fisher H, as the symmetric sandwich C = H^(-1/2) J H^(-1/2), decomposed into **correlation distortion** (temporal, missing higher moments appearing as ghost state-correlation) and **sample/geometric distortion** (per-sample non-Gaussianity).
4. **Direct covariance test**: empirical covariance of the MLE cloud against the sandwich-predicted covariance.

Indicators 1 and 2 should vanish; 3 and 4 quantify how much a misspecified likelihood distorts parameter uncertainty. The anchor is the model's own Gaussian Fisher; the numerical Fisher enters only as the F-vs-G bridge.

The distortion matrix is not only a diagnostic of the error bars: near the maximum it is the correction to the Bayesian evidence (volume ½ log det C, effective-sample α⋆ = p / tr C), which is the paper's "so what" and the bridge to the program's model-comparison component. Stated as motivation; derivation deferred to the bridge-3 study.

## 4) The usage map (what replaces the ranking)

**The deliverable is no longer a ranking, it is a recommendation per regime.** A ranking answers "which is best", which is the wrong question once band C exists, because in band C the best method is the cheapest one. The map answers "which is the cheapest method that still reports its uncertainty honestly, here".

The ranking table below is retained as the **band-A column** of the map, which is the only band measured so far. It is no longer the verdict; it is one column of it.

| Algo | Bias | Variance distortion (band A) | Standing in band A |
|---|---|---|---|
| **IR** | none | to ~1 in Gaussian regime; conservative (~0.5) in the few-channel corner | Calibrated; the default when gating noise dominates |
| **R** | none | residual distortion, size disputed (see D-4) | Usable, correlation distortion at short Δ |
| **MR** | none | over-confident 1.5–2.1 (drops the boundary cross-covariance N·γᵀΣγ that IR keeps) | Cautionary intermediate; the window axis is non-monotone, MR is worse than R |
| **NR** | biased | inflation growing with N_ch | Fails; often does not even reach the MLE |
| **NMR** | none | over-confident ~10–16 | Unbiased but over-confident; possible speed niche (non-recursive, parallelizable) |
| **LSE** | expected biased | not yet measured | **Unmeasured.** The whole point of adding it |

Mechanism of MR's over-confidence: see `theory/macroir/docs/Macro_MRT/macromrt_macromrt_paper_section.md` and `theory/macroir/notes/gvar_i_overcount_audit.md` (use the mechanism, not the inverted May-2026 verdict). Note the audited finding that MR's error is genuine misspecification and not a code defect: it dumps the end-state spread into the observation variance instead of resolving it through the gain.

**Three cells of the map are open and only the runs can close them.**

1. **What LSE does, anywhere.** Expected to be biased and over-confident in bands A and B, and adequate in band C. Both halves are predictions.
2. **Where each method's boundary actually falls** against the predicted A/B and B/C lines of §1a. The interest is in whether the predicted boundary is where the measured one is.
3. **Whether NR survives the comparison at all.** Once LSE occupies the cheap niche, NR may be dominated everywhere: biased where gating matters, and no simpler than LSE where it does not. If so that is a finding, not a gap, and it should be stated. Note that NR is absent from the 2026-07-20 dispatch (`N_ALGO="macro_IR macro_R macro_MR macro_NMR"`), so the new bands will carry no NR column unless it is re-run.

**Do not pre-write any of these three.** The previous version of this table asserted verdicts that the data then contradicted in two cells (D-4), which is how the four-copy drift started.

## 5) Figures (new arc)

Numbering is in flux; this is the intended set.

> **The arc is owned by `results_plan.md`** (`00_master_list.md` C-9 rules it wins over this list). What follows is the intent this plan commits to; the register of what exists lives there.

- **Fig 1: Mechanism.** One filter step across the methods (prior, observation, update, logL), showing where each approximation enters (recursion, averaging, boundary state). Anchor-independent. LSE joins as a sixth column: `ops/local/figure_1_plus_lse.macroir` already builds it, and its diagnostic panel is blocked only until the guard is routed for `family==2`.
- **Fig 2: Recovery clouds.** MLE clouds with three ellipses (empirical, Gaussian-Fisher, sandwich-corrected), per method, showing over/under-confidence. LSE included; its ellipse comes from the routed numerical Fisher (see the LSE plan's fig2 note).
- **Fig 3: Time-resolved calibration.** logL gap, R² residual, score bias, per-interval and accumulated J_t/F_t, score autocorrelation. LSE included, with the two caveats the LSE plan flags: the shared (n/SSE) factor makes its per-interval scores non-martingale, and r̄²_std ≡ 1 is a tautology for LSE, so that row reads "calibrated" by construction and must be annotated.
- **Fig 4: Fisher profiles.** Per-step Fisher; the information about N_ch, k_on and i falls to zero once the open-channel count stops rising and relaxes, while k_off stays informative through the decay.
- **Fig 5: The usage map (the deliverable).** Over **N_ch × noise**, with interval carried inside each cell, spanning bands A through C, with the predicted A/B and B/C boundaries of §1a drawn on it and the measured boundary next to them. This is the figure the paper is for. The old **N_ch × K_off** axes are dead: no K_off axis exists in the scripts (`figure_3_mle_G.macroir` hard-fixes `off = 100`) and building one was costed at ~35,000 CPU-hours.
- **Fig 6: Distortion decomposition.** Correlation vs sample distortion, explaining where and why the approximation fails.
- **Supplement.** Two-path reconstruction of the distortion matrix (`IDM = K·CDM·Kᵀ`, per `correction_idm_reconstruction.md`); band schematic; S1 score-mean; S3 corner; S4 bias and autocorrelation.

**The figure count is an open decision and adding LSE does not settle it.** LSE enters as a column in Figs 1–5, which changes no count. What would change the count is a dedicated LSE-versus-family figure, and nothing yet argues for one. `01_writing_plan.md` §0 item 6 hard-gates at six while D-2 below leaves the number open; that contradiction is now the older of the two and should be resolved toward D-2.

## 6) Data and reproducibility anchor

- Definitive figures anchor on the **Gaussian Fisher**; the all-algorithm Gaussian rerun is in progress and becomes the reference basis.
- **Write against now** (already Gaussian-anchored): Fig 3, Fig 4, Fig 4-alt, S4-acf, and `projects/eLife_2025/figures/paper/sample_correlation_distortion_analysis.md`.
- **Numbers pending the rerun**: Fig 2, S3, figure_5 panel B, distortion_algo_grid.
- Grid: **N_ch × noise**, with interval carried inside every cell, spanning bands A to C (§1a); `scheme_CO` only (drop the 3-state `scheme_CCO`); drop the Taylor `V` variants. The K_off axis is dead, see §5.
- Repo split at code freeze (engine by pinned version): see `09_carve_plan.md`.
- **The n_sims-uniformity hazard, now sharper.** The band-C dispatch of 2026-07-20 runs `n_sims = 1000`; the canonical band-A cells on disk are `n_sims = 10000`, and `433ed13` also mixes 200. Every scalar summary of the distortion matrix is biased in n_sims by Jensen's inequality (`project_distortion_measure_N_dependence`), so a panel that puts a new band-C cell beside an old band-A cell is comparing two different biases and will read as a regime effect. Either hold n_sims fixed within any panel, or use the debiased quadratic (Wald T² − dof). This is the likeliest way for a wrong result to reach print and it now sits on the paper's headline figure.

## 7) Manuscript target

Build structurally on `papers/macroir-elife-2025/docs/manuscript-drafts/elife_paper.tex`. Rewrite both abstracts and the Results spine to the where-it-fails frame; drop the P2X case study; narrow the k-state theory to two states; restage the MacroIR ≈ time-augmented integrated Kalman link as a Discussion hypothesis (verified to ~1e-8 in the prior-art map) rather than a Theory fact; add a Methods section.

Theory sources ready to feed Methods/Supplement:
- Algorithm: `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md`, `macroir_derivation.tex`.
- Distortion machinery: `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex`, `theory/macroir/docs/Gaussian_Fisher_Distortion_Family.md`, `sample_correlation_distortion_analysis.md`.
- Results seed: `papers/macroir-elife-2025/analysis_figure_S1_score_mean.md`.

## 8) Open decisions

- **D-0 [Q] (new, 2026-07-20)** Does **NR** stay in the roster? It is absent from the band-B/C dispatch, and once LSE occupies the cheap niche NR may be dominated everywhere. Three options: re-run it in the new bands; keep it in band A only and say so; or retire it to a sentence explaining why it is dominated. Cheapest is the second.
- **D-1 [Q]** NMR and MR in main text or supplement? Proposal: MR in main text as the strawman (part of the argument for IR); NMR to supplement unless the speed niche is shown. **Reopened by the reframe:** "strawman" was a ranking word. In a usage map a method is not a strawman, it is a method with a domain, possibly an empty one. Re-decide on the map.
- **D-2 [Q]** Number of main figures for eLife (5 or 6); whether Fig 5 (maps) and Fig 6 (decomposition) collapse into one two-row figure.
- **D-3 [Q]** Placement of the key result (Fisher to 0 once the open count relaxes): main Fig 4 as drawn, or supplement. Author reaction in the audio: "me voló la cabeza" but "does not go in the abstract".
- **D-4 [Q]** "Valid" thresholds on the regime map. Candidate: distortion < 1.1 (10% error), bootstrap error < 1.1.
- **D-5 [Q]** Anchor point for bias/distortion. Consensus: optimum / θ_pool; θ_sim used to expose the bias.
- **D-6 [Q]** eLife vs a more specialized venue. Leave the door open, address in the cover letter.
- **D-7 [Q]** Language of these planning docs: keep English (current) or switch to Spanish.

## 9) Next actions

Moved out, 2026-07-14. This section is a pointer, not a list. It was a task ledger inside a document whose own header says "freeze this plan", and a frozen document cannot host a ledger that mutates weekly: within a day it had forked from the one in `00_master_list.md`, with zero overlap between them.

- **Writing order, section gates, definition of done:** `01_writing_plan.md`.
- **Engine work that must land in the frozen commit:** `09_carve_plan.md`, Freeze preconditions.
- **Who owns which topic, and where the documents contradict each other:** `00_master_list.md`.

Framing and program context live in `From molecular mechanisms to data back and forth PROGRAM.md`; the repo carve plan in `09_carve_plan.md`.
