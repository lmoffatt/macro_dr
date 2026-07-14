# Paper 2 (MacroIR / eLife 2025): Master Plan (v2)

> Draft rewritten from the MacroIR 13 audio consensus (2026-07-11). Supersedes `00_master_plan.md` (now retired to a stub; prior content in git history). Comment and correct inline; open items for the author are marked **[Q]**.
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

Characterize where the Gaussian macroscopic likelihood approximations for time-integrated Markov currents stay reliable and where they break, using a minimal two-state channel model, and show that among the family (NR, NMR, R, MR, IR) only **MacroIR stays calibrated across the practical regime**, while mapping the number of channels, interval length and instrumental noise at which even MacroIR begins to distort inference.

The reframe from v1: the paper is no longer "here is MacroIR". It is "here is a validation machinery for time-averaged likelihoods, and here is exactly where each approximation, MacroIR included, succeeds and fails". Finding the failure, and showing it is coherent with the approximations made, is the result.

## 1) Core conceptual frame

Every algorithm in the family makes two Gaussian approximations:

1. **Macro approximation.** The multinomial distribution of channel occupancies is replaced by a multivariate Gaussian. Valid for large N_ch by the central limit theorem; degrades for few channels (multinomial regime).
2. **Interval-likelihood approximation.** The conductance distribution over one measurement interval is replaced by a Gaussian, when the exact object is a complex stochastic-telegraph / Poisson-like distribution. Degrades for intervals much shorter than the relaxation time (telegraphic regime).

This yields **three regimes**: multinomial (few channels), telegraphic/Poissonian (very short intervals), and Gaussian (many channels, moderate intervals, enough instrumental noise), where MacroIR is ideal. Because it is an approximation, MacroIR must fail somewhere; the paper locates that failure and shows it is the predicted degradation of these two approximations, not a coding artifact.

## 2) Scope

### In scope

- Minimal **two-state** channel model (closed-open, `scheme_CO`), single K_on / K_off.
- **Macroscopic** interval-averaged currents.
- Five Gaussian-approximation algorithms: `NR`, `NMR`, `R`, `MR`, `IR` (script naming has used `MNR`; standardize on `NMR`).
- **Non-stationary** protocol (single concentration jump).
- Three control variables: N_ch, interval Δ relative to the kinetic relaxation τ, and instrumental noise.
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

## 4) Algorithm ranking (the verdict)

| Algo | Bias | Variance distortion | Verdict |
|---|---|---|---|
| **IR** | none | to ~1 in Gaussian regime; up to ~1.3 in the corners | Sole survivor; recommended default |
| **R** | none | ~factor 2 residual | Usable, correlation distortion at short Δ |
| **MR** | none | overestimates (double-counts between-final-state variance; drops the boundary cross-covariance N·γᵀΣγ that IR keeps) | Strawman; worse than R |
| **NR** | biased | inflation proportional to N (hundreds) | Fails; often does not even reach the MLE |
| **NMR** | none | underestimates | Only other unbiased method; possible niche where speed matters (non-recursive, parallelizable) |

Mechanism of the MacroMR strawman and the total-variance double-count: see `theory/macroir/docs/Macro_MRT/macromrt_macromrt_paper_section.md` and `theory/macroir/notes/gvar_i_overcount_audit.md` (use the mechanism, not the inverted May-2026 verdict).

## 5) Figures (new arc)

Numbering is in flux; this is the intended set.

- **Fig 1: Mechanism.** One filter step across the five approximations (prior, observation, update, logL), showing where each approximation enters (recursion, averaging, boundary state). Anchor-independent.
- **Fig 2: Recovery clouds.** MLE clouds with three ellipses (empirical, Gaussian-Fisher, sandwich-corrected), per algorithm, showing over/under-confidence. *Regenerate on the Gaussian rerun.*
- **Fig 3: Time-resolved calibration.** logL gap, R² residual, score bias, per-interval and accumulated J_t/F_t, score autocorrelation. *Already Gaussian-anchored.*
- **Fig 4: Fisher profiles, key result.** Per-step Fisher; the information about N_ch, k_on and i falls to zero once the open-channel count stops rising and relaxes, while k_off stays informative through the decay. *Already Gaussian-anchored.*
- **Fig 5: Regime maps.** Bias and distortion heatmaps over **N_ch × K_off**, all algorithms. *Regenerate on the Gaussian rerun.*
- **Fig 6: Distortion decomposition.** Correlation vs sample distortion, explaining where and why the approximation fails.
- **Supplement.** Two-path reconstruction of the distortion matrix (~1 but not exactly, flagged open question); three-regime domain schematic; S1 score-mean; S3 corner; S4 bias and autocorrelation.

## 6) Data and reproducibility anchor

- Definitive figures anchor on the **Gaussian Fisher**; the all-algorithm Gaussian rerun is in progress and becomes the reference basis.
- **Write against now** (already Gaussian-anchored): Fig 3, Fig 4, Fig 4-alt, S4-acf, and `projects/eLife_2025/figures/paper/sample_correlation_distortion_analysis.md`.
- **Numbers pending the rerun**: Fig 2, S3, figure_5 panel B, distortion_algo_grid.
- Grid: N_ch × K_off, plus interval and noise sweeps; `scheme_CO` only (drop the 3-state `scheme_CCO`); drop the Taylor `V` variants.
- Repo split at code freeze (engine by pinned version): see `09_carve_plan.md`.

## 7) Manuscript target

Build structurally on `papers/macroir-elife-2025/docs/manuscript-drafts/elife_paper.tex`. Rewrite both abstracts and the Results spine to the where-it-fails frame; drop the P2X case study; narrow the k-state theory to two states; restage the MacroIR ≈ time-augmented integrated Kalman link as a Discussion hypothesis (verified to ~1e-8 in the prior-art map) rather than a Theory fact; add a Methods section.

Theory sources ready to feed Methods/Supplement:
- Algorithm: `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md`, `macroir_derivation.tex`.
- Distortion machinery: `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex`, `theory/macroir/docs/Gaussian_Fisher_Distortion_Family.md`, `sample_correlation_distortion_analysis.md`.
- Results seed: `papers/macroir-elife-2025/analysis_figure_S1_score_mean.md`.

## 8) Open decisions

- **D-1 [Q]** NMR and MR in main text or supplement? Proposal: MR in main text as the strawman (part of the argument for IR); NMR to supplement unless the speed niche is shown.
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
