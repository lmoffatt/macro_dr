# Results — what it must do, the claim-by-claim spine, and what the data actually support

> Working doc, same genre as `abstract.md`. Opened 2026-07-14, written against a full inventory of `projects/eLife_2025/figures/` (scripts, captions, STATUS/PLAN docs, and what is on disk under `figures/data/`), **not** against `00_plan.md` §5, which is now behind the figures.
> Every number below is quoted with its source file. Numbers with no source are not in the paper.

## The job

The Results have to walk the reader from "here is one filter step" to "here is the map" without ever making them take a claim on trust. The order below is chosen so that each figure answers the objection raised by the previous one. That is the only ordering principle that matters, and it is what turns six figures into an argument.

The through-line, in one sentence per figure:

1. The five algorithms differ, visibly, in a single filter step. (**Fig 1**)
2. That difference shows up as miscalibrated parameter uncertainty. (**Fig 2**)
3. The miscalibration is not a fitting artifact: it is present in the likelihood itself, interval by interval. (**Fig 3**)
4. And it comes from one place: correlation across time, not from any per-sample error. (**Fig 4**)
5. Here is how large it is, everywhere in the regime, for every algorithm. (**Fig 5**)
6. And here is why it is large where it is large. (**Fig 6**)

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

### Fig 1 — One filter step, five algorithms (mechanism, no statistics yet)

`figures/paper/figure_1.Rmd` → `Figure_1.pdf`; whole-recording version `Figure_S1.pdf`. Caption written. **Finished.**

Columns = the five algorithms; rows = prior P_open, observation and innovation, posterior, cumulative logL. Window 6 to 10 ms. Cell: N_ch = 20, f_s = 50 kHz, 100 samples per interval.

**What it claims:** the family is a coordinate system, not a zoo, and the reader can *see* where each approximation enters (does the posterior get corrected; is the conductance instantaneous, start-conditioned, or boundary-conditioned).

**What it must not claim:** anything about accuracy. No error bars, no verdict. It is the visual companion of the Theory section's endpoint ladder (`theory.md` T3), and its whole job is to make the ladder concrete.

### Fig 2 — Recovery clouds: the miscalibration, in the units the reader cares about

`figures/paper/figure_2.Rmd` → `Figure_2.pdf`. Caption written. **Finished.** Data: `433ed13` (`_mle_cloud_runs`, `_battery_pool`, `_battery_sim`). Cell: N_ch = 100, interval 0.1τ, group_size = 10, 10,000 recordings.

Part A: MLE clouds, five algorithms × two parameter pairs (kinetic: k_on, k_off; amplitude: N_ch, i), each with three ellipses (empirical, Fisher-predicted, sandwich-corrected), plus area ratios and bias markers. Part B: the full four-parameter corner, IR alone.

**The claim, and it is the paper's headline claim:** the ellipse the likelihood *reports* and the ellipse the estimates *actually occupy* are different objects, by a factor of 10 to 16 in variance for the non-recursive members (`Figure_S3_caption.md`), and the sandwich correction closes that gap to ≈ 1. Being able to *fix* it, not just measure it, is what makes the distortion matrix a tool rather than a complaint.

**Watch the group_size trap.** `_mle_cloud_runs.csv` mixes group_size 10 and 100 and averaging across them is meaningless (`project_cloud_group_size_column`). Figure 2 uses group_size = 10; if any panel of the Results quotes a cloud-mean bias, it must state which.

### Fig 3 — The calibration cascade over time (the likelihood itself, not the fit)

`figures/paper/figure_3.Rmd` → `Figure_3.pdf`. Caption written. **Finished.** Data: the `figure_3_time_dlik_*` dumps (1000 simulated recordings, N_ch = 100, interval 0.1τ).

Six rows × five algorithms: (A) output and logL, (B) standardized residual r²_std, (C) score bias, (D) per-interval J_t / F_t, (E) accumulated J_T / F_T, (F) score autocorrelation.

**The claims:**
- The recursive trio (R, MR, IR) reaches total logL ≈ **−145**; the non-recursive pair (NR, NMR) ≈ **−230**. A gap of **~87 nats**: under the recursive likelihoods the same data are ≈ e^87 times more probable. Quote it as a log-likelihood gap and let the reader do the exponentiating; the "e^87 times more probable" line belongs in the caption at most.
- **Row D is ≈ 1 for all five algorithms**, and **row E diverges only for the non-recursive pair.** This is the hinge of the entire paper and it deserves to be stated in the text as a result, not left to the figure: *per interval*, every member of the family reports its information correctly. The failure is entirely in how the information accumulates.

**The residual-whiteness row (F) is where the correlation distortion becomes visible**, and it is the empirical form of Milescu's own diagnosis ("the local time correlation of the current"). Say so here, in one clause, and pay it off in the Discussion.

**Note on the dumps:** every evolution row is duplicated ×2 in these files (blank-segment plus segment=0 copies; `project_figure3_dlik_per_interval_score`). The scripts dedup. Anyone re-deriving a number from the raw CSVs must too, or they will find a spurious factor of 2.

### Fig 4 — Where the information lives, and where the overconfidence comes from

`figures/paper/figure_4.Rmd` → `Figure_4.pdf`. Two captions written (main and alt). **Finished.** Same data as Fig 3.

Per-step Fisher F_t against per-step score variance J_t, four parameters, five algorithms.

**Two claims, and they are the two best results in the paper.**

1. **The identity J_t = F_t holds per step for all five algorithms** (the ratio row sits on zero within bootstrap CI). Therefore the non-recursive overconfidence is **entirely** cross-time score correlation, and not a per-sample modelling error. This is a mechanism, it is clean, and it is exactly the answer to Milescu's open question.
2. **The information about k_on, the unitary current i, and N_ch falls to zero once the agonist is removed and the open population relaxes; k_off stays informative through the decay.** It is a property of the macroscopic observable, so it survives every approximation in the family. Del Core and Mirams 2025 declare this problem open.

**Placement decision (D-3 in the master plan) is now easy: this is main text.** The author's own reaction ("me voló la cabeza") plus a 2025 open problem plus algorithm-independence is three reasons. It stays out of the abstract; it does not stay out of the Results.

### Fig 5 — The validity map (the paper's deliverable)

**Not finished. This is the critical-path figure.** Composition is decided (`figures/in_progress/figure_5_PLAN.md`, 2026-07-08, which sorts 44 candidates into 5 MAIN / ~19 SUPP / 16 DROP), but no caption exists and the panel selection has not been executed into a final PDF.

The PLAN's main-text arc:

- **`figure_5_distortion_algo_grid.Rmd`** (433ed13, five algorithms × two parameters × three noise levels). The headline number: the fraction of the (Δ·k_off × N_ch) plane where the likelihood is self-consistent within ±15%, as noise goes 0.1 → 1 → 10: **IR 0.88 → 0.99 → 1.00; R 0.27 → 0.50 → 0.75; MR 0.05 → 0.34; NR and NMR ≈ 0.** This single row of numbers is the paper's verdict, and it is quantitative, continuous, and it names the conditions. It is what the abstract's last sentence promises.
- **`figure_5_bias_empirical_grid.Rmd`** (433ed13): IR's bias stays inside the tolerance over 98% of the plane; R, MR and NR fail over roughly half.
- **`figure_5_IR_channel_pooled_loss.Rmd`** (1c2ae6f): concentrating channels into fewer patches is ≈ free for k_off, costs ≈ **4×** for k_on, and up to ≈ **500×** for the unitary current i. This is the one panel that gives the experimentalist something to do differently on Monday, and it should not be buried.

**The anchor problem, and it is real.** Cross-algorithm panels can only come from `433ed13`, which is the **numeric-Fisher** run and the only commit with all five algorithms. The IR mechanism panels come from `1c2ae6f`, the **Gaussian-Fisher** run, which has no NMR at all. So the main text is not all on one anchor. `figure_5_PLAN.md:5` states this and says a reader must be told. **Two options and they need a decision:**
  (a) tell the reader, in Methods and in the caption, and justify why it does not matter (Fig 4's per-step J_t = F_t identity is arguably the justification: it says the two anchors agree per step);
  (b) run the missing Gaussian-Fisher cells (at minimum NMR, which is entirely absent from 1c2ae6f) and put the whole main text on one anchor.
  This contradicts `00_plan.md` §6, which asserts the definitive figures anchor on the Gaussian Fisher. **They currently cannot.** Decide before drafting, because option (b) is a compute request, not a writing task.

### Fig 6 — Why: the distortion decomposition

**Not finished.** PLAN's proposal: `figure_5_distortion_decomp_grid.Rmd` + `figure_5_IR_sample_peak_in_Nch.Rmd`.

**The claims:**
- **The two failure modes are separable and they live in different corners.** NR and NMR blow up in the **correlation** term at short intervals and in the **sample** term at coarse intervals. Two Gaussian approximations, two failure modes, each visible where the theory says it should be. This is the section that makes the whole thing feel explained rather than measured.
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

The trustworthiness contour, if one is drawn, is at distortion = **1.05** in `figure_7_validity_map.Rmd` and at **±15%** in `figure_5_distortion_algo_grid.Rmd`. **Pick one and use it everywhere.** (D-4 in the master plan proposed 1.1. Three different thresholds in three places is exactly the kind of thing a reviewer notices.)

## Verify before submission

- The direction convention of the distortion ratio (which way is over-confident), once, against the code. Every verdict in the paper hangs on it and it is currently asserted in three places with two different signs.
- The ±15% self-consistency fractions, recomputed on whatever anchor the final Fig 5 uses.
- The 87-nat logL gap: confirm it is the recursive/non-recursive split and not a recursive/averaging split, because the sentence that reports it will be read as a claim about *which* approximation matters.
</content>
