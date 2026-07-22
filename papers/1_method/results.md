# Results — what it must do, the claim-by-claim spine, and what the data actually support

> Updated: 2026-07-20 (A-strict, paper 1). Written against a full inventory of `projects/eLife_2025/figures/`.
> Every number below is quoted with its source file. Numbers with no source are not in the paper.
>
> **A-strict reframe.** Paper 1's roster is `R`, `MR`, `VR`, `IR` (all recursive). The non-recursive
> members `NR`, `NMR` are named once in Theory and **measured in no figure**; their quantitative
> results (the 87-nat logL gap, the 10–16× overconfidence) are **relocated to paper 2** and must not
> anchor a paper-1 claim (`decisions.md`). The figures currently on disk are five-algorithm and need
> re-rendering to R/MR/VR/IR once `VR` is run.

## The job

The Results walk the reader from "here is one filter step" to "here is the within-family validity map"
without making them take a claim on trust. Each figure answers the objection the previous one raises.
That is what turns the figures into an argument.

The through-line, one sentence per figure (A-strict):

1. The likelihood is a recursive filter over intervals, and its two poles — no interval conditioning (R) and full boundary conditioning (IR) — are visibly different objects. (**Fig 1**)
2. That difference shows up as miscalibrated parameter uncertainty. (**Fig 2**)
3. The miscalibration is not a fitting artifact: it is present in the likelihood itself, interval by interval. (**Fig 3**)
4. And it comes from correlation across time, not from any per-sample error. (**Fig 4**)
5. Here is how large it is across the gating-dominated regime, for the recursive family. (**Fig 5**)
6. And here is *which piece of the interval structure* produces it: the `MR → VR → IR` mechanism, variance step then gain step. (**Fig 6**)

The spine's centre of gravity moved from Fig 5 (the map) to Fig 6 (the mechanism): with the map's
dramatic cross-family contrasts handed to paper 2, paper 1's novel result is the mechanism, and Fig 6
is where it is proved.

## Three corrections to the record, before anything is written

The figure inventory contradicts three things currently written down elsewhere. Fix the sources, do not propagate them.

### 1. The overconfidence factor is 10 to 16, not 14 to 21

`abstract.md:102` says "over-confident by roughly 14 to 21 fold". The figure says otherwise. `Figure_S3_caption.md`: empirical/Fisher variance ratio ≈ **10 to 16** for NR and NMR, **1.2 to 1.4** for R, **1.5 to 2.1** for MR, and the sandwich-corrected ratio ≈ **1** in every cell of all four. **The caption is the authoritative number.** Correct the abstract.

### 2. The "MR sign discrepancy" is probably a category error, and it dissolves

`abstract.md` and `discussion.md` both flag an apparent contradiction: the master plan says MR *overestimates variance* while Figure S3 reports MR as *over-confident* (ratio 1.5 to 2.1, meaning the empirical parameter variance exceeds what the Fisher predicts).

These are two different variances and both statements can be true at once:

- MR overestimates the **predicted observable variance** at each interval, because it drops the boundary cross-covariance N·γᵀΣγ that IR keeps.
- MR nevertheless **under-reports the parameter covariance** (over-confident), because the information it claims is not delivered once the cross-time correlation is accounted for.

**Check this against the data before writing either sentence, then fix the master plan's column header**, which says "Variance distortion" without saying which variance. If the two statements really are about different objects, they are both usable, and the ambiguity is the only bug.

### 3. The bidirectionality claim has a different carrier than the abstract thinks

The abstract needs some algorithm to *deflate* the information, and it currently hangs that on MR. But every non-IR member appears to inflate (S3 ratios all > 1). The real carrier is elsewhere and it is better:

`figure_5_distortion_algo_grid.Rmd` reports that **IR's own departures from calibration run in the under-confident (conservative) direction, down to ~0.5 in the low-N_ch corner**, while NR, NMR, R and MR all inflate. So the distortion is bidirectional in exactly the way that matters: the approximations that ignore the acquisition window over-state what the data know, and the one that respects it, when it fails, under-states it. **That is a better sentence than the one in the abstract, and it is safe (it is the direction that does not mislead a user).** Verify the direction convention once against the code, then rewrite the abstract's directionality clause on this.

## Figure-by-figure spine

> **A-strict scoping banner.** The sections below were written for the five-algorithm arc. Read every
> NR/NMR-specific quantity as **paper 2's**, not paper 1's: the 87-nat recursive-vs-non-recursive logL
> gap (Fig 3), the 10–16× overconfidence (Fig 2/S3), and the "NR/NMR blow up" framing (Fig 6). The
> "recursive trio / naive pair" partition is gone; paper 1's members are all recursive. Where a section
> quotes an NR/NMR number, substitute the recursive-family analog (R vs the averaged members) for
> paper 1, and leave the dramatic number for paper 2. The `VR` column is new and unrun; any figure
> that needs it is blocked on the implementation.

### Fig 1 — One filter step along the ladder, R → MR → VR → IR (mechanism, no statistics yet)

**BUILT 2026-07-21.** `figures/paper/figure_1.Rmd` (roster) + `figure_1_panels.R` (panels, shared) →
`Figure_1.pdf`, supplement `Figure_S1.pdf`, caption rewritten. The six-algorithm version lives beside
it as `figure_1_all.Rmd` → `Figure_1_all.pdf`, for orientation only.

**Four columns, the recursive ladder (2026-07-21, Luciano).** Recursion is held fixed and the columns
vary only what the interval-averaged conductance is conditioned on: R (instantaneous), MR (interval
average conditioned on the start state, total variance), VR (same mean and update, residual variance),
IR (both endpoints, boundary cross-covariance in the gain). Rows: prior P_open, observation and
innovation, posterior, cumulative logL. Window 6 to 10 ms. Cell: N_ch = 20, f_s = 50 kHz, 100 samples
per interval.

**This supersedes the A-strict two-column call of 2026-07-20** (`decisions.md`), whose reasoning was
that MR and VR would be near-duplicate columns separated only by a variance band. The rendered figure
does not bear that out: VR's band is visibly the narrower one wherever the channels gate, and because
the predictive variance divides the gain, VR's posterior and propagated mean also part from MR's
within the window. The four columns show the *path*; Fig 6 still owns the quantitative version of the
same step.

**What it claims:** the likelihood is a recursive filter, and the reader can *see* the two ends of the
conductance axis — instantaneous vs boundary-conditioned — and where the conditioning enters the
update. It is the visual companion of the Theory ladder (`theory.md` T3).

**What it must not claim:** anything about accuracy. No error bars, no verdict.

### Fig 2 — Recovery clouds: the miscalibration, in the units the reader cares about

`figures/paper/figure_2.Rmd` → `Figure_2.pdf`. Caption rewritten. **Rebuilt 2026-07-22 with the
A-strict roster on the GAUSSIAN anchor**: R, MR, VR, IR, with VR read from `0ffbda7` and the rest from
`1c2ae6f`, all nsim 10000. Cell: N_ch = 100, interval 0.1τ, group_size = 10, 10,000 recordings.

Part A: MLE clouds, four algorithms × two parameter pairs (kinetic: k_on, k_off; amplitude: N_ch, i), each with three ellipses (empirical, Fisher-predicted, sandwich-corrected), plus area ratios and bias markers. Part B: the full four-parameter corner, IR alone.

**The claim, and it is the paper's headline claim:** the ellipse the likelihood *reports* and the
ellipse the estimates *actually occupy* are different objects, and the sandwich correction closes the
gap to ≈ 1. Measured empirical-over-Fisher area ratios at this cell: **R 1.32, MR 1.97, VR 2.18,
IR 1.02** on the kinetic pair and **1.09, 1.53, 1.77, 1.00** on the amplitude pair; corrected, all of
them within a tenth of one. Being able to *fix* it, not just measure it, is what makes the distortion
matrix a tool rather than a complaint.

**And this panel is where the mechanism claim is decided.** VR is the most over-confident rung, worse
than MR, which is the falsifier's prediction: taking the variance out without the boundary gain makes
the reported uncertainty worse, so the gain is what recovers calibration. Write the confirming branch
(`decisions.md`). The 10-to-16 factor that used to be quoted here is an NR/NMR number and belongs to
paper 2.

**Anchor note, now resolvable:** the two anchors agree at this cell (Gaussian R 1.32 / MR 1.97 /
IR 1.02 against the numeric run's R 1.18-1.42 / MR 1.52-2.07 / IR ≈ 1, `decisions/D-4_ranking_verdict.md`),
so moving the main text onto the Gaussian anchor costs no result and the open two-anchor decision at
the Fig 5 section can be closed the same way.

**Watch the group_size trap.** `_mle_cloud_runs.csv` mixes group_size 10 and 100 and averaging across them is meaningless (`project_cloud_group_size_column`). Figure 2 uses group_size = 10; if any panel of the Results quotes a cloud-mean bias, it must state which.

### Fig 3 — The calibration cascade over time (the likelihood itself, not the fit)

`figures/paper/figure_3.Rmd` → `Figure_3.pdf`. Caption rewritten. **Rebuilt 2026-07-22** on the
A-strict roster over regenerated dumps (engine `0ffbda7`, `seed = 20260722`, 1000 recordings,
N_ch = 100, interval 0.1τ). Every number below is printed by the notebook's caption-numbers chunk on
each knit; re-read it after a re-render rather than copying these forward.

Six rows × four algorithms: (A) output and logL, (B) standardized residual r²_std, (C) score bias,
(D) per-interval J_t / F_t, (E) accumulated J_T / F_T, (F) score autocorrelation.

**The claims:**
- Total logL: R −147.2, MR −153.9, VR −148.9, IR −143.4, each ± ~0.2. The ladder spans **10.4 nats**
  and IR ends highest. **The old ~87-nat gap is gone from this figure and belongs to paper 2**: it was
  the recursive-versus-naive contrast, and the naive pair is no longer a column.
- **Row E separates them**: R 1.07 (CI 0.99–1.17), MR 1.19 (1.10–1.29), VR 1.20 (1.09–1.30),
  IR 1.08 (1.00–1.17). MR and VR are the pair whose interval clears one.
- **CORRECTION, 2026-07-22, and it sharpens the hinge.** The long-standing claim that "row D ≈ 1 for
  all algorithms, so per interval everyone reports its information correctly" is **not what the data
  say**, and Figure 4 now measures it per parameter: median per-step log10(J_t/F_t) is about −0.18 for
  R, −0.25 for MR, −0.15 for VR and **+0.01 for IR**, and only IR's bootstrap interval covers zero at
  most steps (64–91% of them against 0–28% for the rest). Per step the score varies *less* than the
  information predicts, by up to ~1.5× for MR.
  **The hinge is better stated as a sign flip:** per step J_t < F_t, yet accumulated J_T/F_T > 1. A
  discrepancy that reverses sign between the per-step and the accumulated statistic cannot be a
  per-sample modelling error; it is exactly what positive cross-time score correlation does. Say it
  that way in the text and in Fig 4's caption; the old wording was a compression that the numbers do
  not support.
- **VR is no better than MR here either**, which is the time-domain companion of Figure 2's parameter
  ratios and the same verdict: removing the variance without the boundary gain does not help.

**Row F is where the correlation distortion becomes visible**, and it is the empirical form of
Milescu's own diagnosis ("the local time correlation of the current"). At lag one: IR −0.004 ± 0.004,
indistinguishable from zero, against R 0.191, VR 0.256, MR 0.357. **Only the fully conditioned filter
produces a white score**, and that is a cleaner statement than the old one because it now separates
members that share the recursion. Say it here in one clause and pay it off in the Discussion.

**Do not equate row E with Figure 2's overconfidence factor.** Same direction, same over-confident
pair, different scalar: R sits at 1.07 here against 1.32 for its kinetic-pair ellipse in Figure 2,
because E is a single-parameter (k_off) accumulated ratio and Figure 2's is a two-parameter ellipse
area. The previous caption asserted the identification; it does not survive the numbers.

**Note on the dumps:** every evolution row is duplicated ×2 in these files (blank-segment plus segment=0 copies; `project_figure3_dlik_per_interval_score`). The scripts dedup. Anyone re-deriving a number from the raw CSVs must too, or they will find a spurious factor of 2.

### Fig 4 — Where the information lives, and where the overconfidence comes from

`figures/paper/figure_4.Rmd` → `Figure_4.pdf`. Two captions written (main and alt). **Finished.** Same data as Fig 3.

Per-step Fisher F_t against per-step score variance J_t, four parameters, five algorithms.

**Two claims, and they are the two best results in the paper.**

1. **The identity J_t = F_t holds per step for all five algorithms** (the ratio row sits on zero within bootstrap CI). Therefore the non-recursive overconfidence is **entirely** cross-time score correlation, and not a per-sample modelling error. This is a mechanism, it is clean, and it is exactly the answer to Milescu's open question.
2. **The information about k_on, the unitary current i, and N_ch falls to zero once the agonist is removed and the open population relaxes; k_off stays informative through the decay.** It is a property of the macroscopic observable, so it survives every approximation in the family. Del Core and Mirams 2025 declare this problem open.

**Placement decision (Q-3 in `00_plan.md` §8, labelled D-3 before 2026-07-21) is now easy: this is main text.** The author's own reaction ("me voló la cabeza") plus a 2025 open problem plus algorithm-independence is three reasons. It stays out of the abstract; it does not stay out of the Results.

### Fig 5 — Where IR itself stops being faithful (BUILT 2026-07-22)

`figures/paper/figure_5.Rmd` → `Figure_5.pdf`, caption written. **Subject changed**: this is no longer
the within-family validity map, it is the figure `00_plan.md` §1 promised and never staffed, the one
that locates MacroIR's own failure. `macro_IR` alone; MR and VR have done their work by Fig 4.

Three blocks: **A** grouped MLE distributions at three design points (distorted / biased / faithful),
**B** the distortion-induced bias (first moment, anchored at θ_sim), **C** the information distortion
(second moment, anchored at θ_pool). Columns N_ch {10, 20, 50, 100, 1000, 10000}; inside each panel
x = Δ·k_off over the seven free intervals, y = the six noise levels. Every colour is CI-shrunk, so
what is coloured is significant at 95%.

**The claims, all measured** (numbers printed by the notebook on each knit, `figures_build_plan.md`
§3d):

- **IR fails, and where the theory said it would.** The telegraphic edge, worst at few channels and
  low noise, plus a second narrower edge at the coarsest interval. The dependence on Δ is a **U**,
  minimum at Δ·k_off = 0.5, so the acquisition interval has an optimum rather than a direction of
  improvement.
- **The two moments fail at opposite ends of that axis.** Bias worst at coarse Δ, distortion worst at
  short Δ. **No single interval optimises both**, and that is the decision rule this figure hands the
  Discussion.
- **The failure is signed and direction-dependent**: k_off over-confident where N_ch is conservative,
  by comparable amounts, so any scalar averaged over parameters cancels them and reports a faithful
  algorithm. State the map per direction, never as one number.
- **Attribution** (measured, not plotted here): the short-Δ failure is correlation-led about eight to
  one, i.e. the interval closure fails to decorrelate successive intervals. Where this goes is open,
  see below.
- **The cheap diagnostics are blind to all of it**: r̄²_std off by ≤5.5% and integrated autocorrelation
  ≤4.8% over the same corner. That is the argument for why the machinery is needed, and it belongs in
  the text.

**Two caveats the text must carry, not bury.** The ellipse of panel A is a marginalised sandwich and
is systematically milder than the map; at the current point 1 its N_ch direction is unresolved by 100
fits and the claim rests on block C. And the MLE distribution is **not normal at few channels**
(skewness 4.13 at N_ch 10 / group 10), which limits what any ellipse or covariance can mean there.
That is a limitation of the *estimator's* finite-sample behaviour and **not a defect of the
likelihood** — a correct likelihood can produce a skewed MLE — so it deserves its own sentence rather
than being folded into the distortion story. Empirical coverage of the nominal 95% interval is
0.80–0.93 at N_ch 20 against 0.94–0.98 at N_ch 200.

**[Q] Open:** where the sample/correlation decomposition goes now that it is out of Fig 5 (candidate:
Fig 6 beside the R-vs-IR map, or supplement), and whether the paper stays gated at six figures.

### Fig 6 — Why: the mechanism (MR → VR → IR), and the distortion decomposition

**Not finished, and now paper 1's headline figure** (the map's cross-family drama went to paper 2, so
the mechanism is what paper 1 is for). Two jobs.

**Job 1, the new one: the `MR → VR → IR` mechanism, each step isolating one thing.** This is paper 1's
central novel result and it needs the `VR` column, which is not yet run (`decisions.md`;
`theory/macroir/notes/vr_variance_form_plan.md`).
- **MR → VR (the variance step).** MR carries the total per-start-state variance and double-counts the
  spread across end states `Var_j[gmean_ij|i]`; VR uses the residual (boundary-conditioned) variance.
  The panel: MR's predicted y_var sits above VR's above the realized residual variance, and VR matches
  the realized. **The prediction is that VR is *over*-confident in parameter space** (it removed real
  observation variance without gaining the boundary information), and if it is not, the mechanism thesis
  is wrong (write both branches).
- **VR → IR (the gain step).** IR keeps the boundary cross-covariance N·γᵀΣγ in the gain; VR does not.
  The panel isolates the calibration recovered by the gain alone.
This turns "MR's problem is the gain, not the variance" from an algebraic assertion into a measurement,
which is exactly what a method paper should do.

**Job 2, carried over: the sample/correlation decomposition for the recursive family.**
- **The two failure modes are separable and live in different corners.** The correlation term dominates
  at short intervals, the sample term at coarse intervals. (The dramatic version of this, NR/NMR
  blowing up, is **paper 2**; paper 1 shows it on the recursive family, where it is milder and is the
  stringent test of the machinery.)
- **The sample distortion is non-monotonic in channel number**, peaking at N\* = 2σ²/G, and the peak marches from N ≈ 10 to N ≈ 1000 as the noise rises (fitted N\* ∝ noise^0.80, R² = 0.59). Instrumental noise *helps* the Gaussian approximation, and there is a worst channel number for a given noise level. This is genuinely novel, it is counter-intuitive in the right way, and it is the single most interesting mechanistic result in the paper.
- **A correctness point that belongs in the text, not the caption**: at short intervals the large correlation distortion is **not a model failure**. It is the correct information-redundancy cost of oversampling. Without that sentence the map's left edge reads as a bug.

**Do not claim the multiplicative decomposition as demonstrated.** `figure_5_PLAN.md:99`: *"none shows sample × correlation = total J/F on the same axes. The 'multiplicative' claim is asserted, never demonstrated."* Either produce the panel that demonstrates it (the additivity of the log-determinants is exact for IR, R and MR per `figure_5_master_STATUS.md`, so the panel is cheap) or state it as a definition rather than a finding.

### Fig 7 family — IR-only contours (currently seven competing layouts)

`figures/paper/figure_7_*.Rmd`, seven scripts, all IR-only, all on `1c2ae6f`, **no captions, no selection made**. They carry the noise axis and the r̄_std / r̄²_std calibration rows.

**Recommendation: this is not a main figure.** It is either (a) the supplement that supports Fig 6's mechanism claims, or (b) the source of one or two panels that get absorbed into Fig 5. Seven layouts of the same object is a sign the question they answer has not been decided. **[Q] Decide what question Fig 7 answers that Figs 5 and 6 do not.** If there is no such question, drop it to supplement and move on; the paper does not need it.

One finding from that family that *is* worth main-text space, wherever it lands: **at few channels IR over-predicts the observable variance by about 6%** (r̄²_std against the θ_sim anchor; it must be `battery_sim_G`, because the pool anchor averages the effect back to ≈ 1 and hides it; `project_rstd_calibration_sim_not_pool`). That is IR's own honest failure and the paper is better for showing it.

## Supplements, as they stand

| Figure | Status | What it carries |
|---|---|---|
| S1 (`figure_1.Rmd`, whole recording) | finished | the filter step over the full protocol |
| S2 | finished | Spearman ρ(logL, Mahalanobis outlier distance) ≈ 0 in every cell: **outliers are not driving the clouds** |
| S3 | finished | the corner plots for NR, NMR, R, MR, and the authoritative 10-16 / 1.2-1.4 / 1.5-2.1 ratios |
| S4, S5 | finished | bias and autocorrelation at residual and score level |
| score-mean supplement | drafted (`analysis_figure_S1_score_mean.md`, in Spanish) | significant score-bias steps out of 40: NR 6/7/7/6, NMR 6/5/8/6, R 24/16/38/24, MR 26/16/38/27, **IR 2/5/2/4 (chance level is 2)**; max abs mean score NR up to 10.03, IR 0.30 to 0.71 |

**Naming collision:** `analysis_figure_S1_score_mean.md` calls its figure S1, and `figure_1.Rmd` already emits `Figure_S1.pdf`. Renumber before anything cites either.

**The score-mean supplement contains a puzzle worth resolving before it is published:** R and MR show *more* significant score-bias steps (24-38 of 40) than NR and NMR (5-8 of 40), while IR sits at chance. If that survives scrutiny, the "unbiased score" story is not monotone along the ladder, and the Discussion's clean "recursion fixes the uncertainty, not the bias" line needs re-checking. **[Q] Is this a power effect** (the recursive likelihoods have smaller score variance, so the same absolute bias becomes significant), rather than a larger bias? Almost certainly. If so, report the standardized bias, not the count of significant steps, and the puzzle dissolves. If not, it is a finding.

## Results claims that the data do NOT currently support

Say none of these until the gap is closed.

- **"Definitive figures are anchored on the Gaussian Fisher."** Not true for anything cross-algorithm (see Fig 5). NMR has no Gaussian-Fisher run at all.
- **"Sample × correlation = total distortion."** Asserted, never plotted (`figure_5_PLAN.md:99`).
- **The N\* ∝ noise law as a validated prediction.** N\* is measured but never overlaid on its own predicted law N\* = 2σ²/G (`figure_5_PLAN.md:100`), and the fitted exponent is 0.80 with R² = 0.59, which is not a clean 1/2 or 1. Report the measurement; do not present the law as confirmed.
- **The micro-versus-macro comparison** (`figures/in_progress/figure_6_micro_macro_linear.Rmd`, the newest script in the tree). Blocked on data: micro_IR exists at N_ch = 5 only, micro_R has no run at N_ch ≥ 1000. It is out of scope for this paper anyway (`00_plan.md` §2) and should not sneak in through a figure.
- **Anything from the flagged-wrong scripts.** `figure_5_PLAN.md` names them: `IR_distortion_ridge` (the argmax ridge is noise), the `relCI_by_Nch` panel of `distortion_gauss_vs_numeric_lines` (median-over-intervals artifact that contradicts its own sibling), `domains_schematic` (hard-codes an N_ch = 30 boundary that contradicts the paper's own N\* result), and three scripts with stale headers (`decomp_gaussian`, `IR_covariance_design`, `IR_interval_loss`). **Do not lift a number from any of these.**

## Framing decision the maps force

`figures/in_progress/figure_7_validity_map.Rmd` already states it: *"Validity map, NOT a regime map ... no hard regime boundaries: the distortion is a continuous gradient in (N, Δ, noise), not domains."*

That is right, and it is rule 3 of `title.md` (continuous, not a threshold). But it sits awkwardly with the Theory section's three named regimes (multinomial, telegraphic, Gaussian; `theory.md` T2). The resolution, and it should be written into both sections: **the three regimes are named as the two approximations' asymptotic corners, not as territories with borders.** The map shows a gradient; the corners explain its direction. Any figure that draws a hard line on the plane (the `domains_schematic`) is off-message and is already flagged as method-wrong.

The trustworthiness contour, if one is drawn, is at distortion = **1.05** in `figure_7_validity_map.Rmd` and at **±15%** in `figure_5_distortion_algo_grid.Rmd`. **Pick one and use it everywhere.** (Q-4 in `00_plan.md` §8, labelled D-4 before 2026-07-21, proposed 1.1. Three different thresholds in three places is exactly the kind of thing a reviewer notices.)

## Verify before submission

- The direction convention of the distortion ratio (which way is over-confident), once, against the code. Every verdict in the paper hangs on it and it is currently asserted in three places with two different signs.
- The ±15% self-consistency fractions, recomputed on whatever anchor the final Fig 5 uses.
- The 87-nat logL gap: confirm it is the recursive/non-recursive split and not a recursive/averaging split, because the sentence that reports it will be read as a claim about *which* approximation matters.
</content>
