# The macro paper — decision log

> Updated: 2026-07-28. **This file governs the MERGED macro paper** (former papers 1 and 2, fused
> 2026-07-23; `../_program/decisions.md` §1). Its **own** settled decisions; everything that binds
> the micro paper too moved to `../_program/decisions.md`, which this file cites and does not restate.
>
> Open decisions live in `00_plan.md` §8.

## What this paper is

- **The question, in the reader's terms: what is the best I can do with a macroscopic recording, and
  how would I know?** The paper answers it by measuring, over the design space, which methods report
  their own uncertainty honestly and which do not.
- **The comparison anchor is least squares** (2026-07-28). MacroR has essentially no uptake, so a
  paper pitting one unused algorithm against another is unsellable. The paper compares **the method
  everyone uses** against the new one and shows what more is available and where.
- **Body roster: `LSE`, `NR`, `R`, `IR`**, faceted by N_ch. As a ladder of cost it reads: fit the mean
  and discard the fluctuations → use the gating variance without a filter → filter on instantaneous
  samples → condition on the interval. Reading LSE against IR answers "do you need a likelihood?";
  reading R against IR answers "what does conditioning buy?".
- **The four macro members are also a 2×2, and since 2026-08-01 that is the reading that carries the
  mechanism.** Interval window × recursion: `NR` neither, `INR` window only, `R` recursion only, `IR`
  both. **Each margin moves a different moment**, which the ladder ordering hides because it presents
  one axis. Measured on the corrected `INR` run at `1f7138b`; the table and its provenance are
  `decisions/recompute/d3_interval_vs_recursion_2x2.py` (Luciano stated the decomposition from the
  2026-08-01 audio; the numbers below are the recompute).
  - **The window fixes the first moment. Recursion does not, and alone it makes it worse.** Median
    |bias| over every cell and interval on disk, log10: in `N_ch` / `i`, `NR` 0.099 / 0.085, `R`
    0.110 / 0.101, `MR` 0.120 / 0.103, against `INR` 0.003 / 0.003 and `IR` 0.002 / 0.001. `VR` sits
    between at 0.058 / 0.066. `k_off` is unbiased for every member (≤ 0.002). So the bias is an
    **amplitude-pair** effect of a quarter to a third, it does not shrink with N_ch, and the only two
    members free of it are the two that carry the interval window.
  - **Recursion fixes the second moment, and the window does not — beyond the per-sample part.**
    Information distortion at noise
    0.1, N_ch 10⁴, on `k_off` (total / sample / correlation): `NR` 78.9 / 47.2 / 21.2, `INR` 22.2 /
    1.00 / 22.1, `R` 1.46 / 1.00 / 1.47, `IR` 1.01 / 1.01 / 1.00. **`INR`'s distortion is
    all correlation** — total and correlation agree to the third digit, and they do so across the
    whole noise fan, which is Luciano's "toda la distorsión se da a nivel de la correlación" as a
    measurement. The window does buy the *sample* part (`NR` 47.2 → `INR` 1.00 at N_ch 10⁴), so
    "the window buys nothing without recursion" was false in both moments.
  - **One end is worse than none, which is the sharp form of the both-ends claim.** `MR` against `R`
    on `k_off` at N_ch 10⁴: total 1.94 against 1.46, correlation 2.00 against 1.47; and in `N_ch`
    bias 0.120 against 0.110. Conditioning the interval mean on the start state alone degrades both
    moments relative to not conditioning at all, and only conditioning on both ends closes it
    (`IR` 1.01, bias 0.002). This is the measured version of "si condicionás a uno solo, todo da mal".
  - **`LSE` is unbiased too and fails the same way `INR` does**, which is why the first moment cannot
    carry the paper's argument: 0.001 in both parameters it fits, and a distortion of 13.3 on `k_off`
    at N_ch 10⁴ that is 14.1 correlation and 0.95 sample. LSE and INR are the same story about the
    second moment told at different cost.
- **Supplement roster: `MR`, `VR`, `INR`** (2026-07-29; `NMR` → `INR` 2026-07-31). `MR` and `VR` split
  the R → IR step and carry the mechanism. `INR` is the published `MacroINR`, the control of Comm Biol
  2025.
  - **`INR`'s job is now known, and it is a body job, not a supplement job — REOPENED, Luciano's
    call** (2026-08-01). The old reading (it sits beside `NR` because the two are numerically
    indistinguishable, stated as *the interval-mean conductance buys nothing without recursion*) is
    **dead**: that measurement was made on `NMR`, the build missing the `N·ms` interval-variance term,
    so the near-identity had a mechanism and was never evidence about the method. The corrected `INR`
    separates from `NR` in both moments (numbers above), and it is the only member that isolates the
    window margin. **The body roster question Q-1/Q-2 is therefore reopened**: `LSE, NR, R, IR` is
    three of the four cells of the factorial, and adding `INR` completes it at the cost of one column.
    Do not restate the old conclusion anywhere.
  - **`VR` keeps its name** (2026-07-28, `../_program/nomenclature.md`), displayed as "Variance
    Recursive". It is the control that turns "MR's problem is the gain, not the variance" from algebra
    into measurement, and it fired. **Do not describe it as "MR→VR changes only the variance, VR→IR
    only the gain":** the predictive variance divides the gain, so the variance step moves the update
    too, and **at the same prior** IR's *total* predicted variance is algebraically equal to MR's
    (`figures_build_plan.md` §F1-2). Along a recording it is not: the gains drive the priors apart and
    MR then runs 69% to 81% above IR (figure-1 dumps, re-verified 2026-08-06). VR is below both from
    any state.
  - **`INR` is a supplement member, not dropped** (2026-07-29). The old "no literature attribution"
    reason was false: it is `MacroINR` in print. See `../_program/decisions.md` §2.
  - **`NMR` is kept as a named object, second class but real** (Luciano, 2026-07-31): the defective
    implementation that ran between `a3e0a89` and `1f7138b` and produced every `macro_NMR` file in the
    freeze. It is not a legacy spelling of `INR` and must not be swept out of the data or the `.Rmd`
    readers. Whether it also earns a *column* — `NMR` beside `INR` is a direct read-out of what the
    interval-variance term does — is an open figure call, not settled here.
  - **The supplement panels carry the FULL LATTICE, all six** (Luciano, 2026-07-29): `NR`, `INR`, `R`,
    `MR`, `VR`, `IR`. Not the three demoted members alone. They are *defined* by their position
    relative to the body members, so a panel without `R` and `IR` cannot state "MR is worse than R,
    VR worse still, IR closes it", and a panel without `NR` cannot place `INR` against the rung below.
    Three consequences, all favourable: `figure_1_all.Rmd` already renders exactly those six columns,
    so this is columns added to existing scripts and not a new figure; the columns shared with the body
    are a free consistency check, because a body and a supplement that disagree on `NR` or `IR` means a
    bug; and **`LSE` stays out**, which is both conceptually right (it is not on the lattice, and the
    supplement is about what conditioning buys within the lattice) and keeps **n_sims uniform at 10⁴**,
    since the LSE arm is at 1000 and mixing them is the Jensen hazard.
- **Literature positioning.** `R` carries the recursive lineage (Moffatt 2007; Münch 2022, a published
  Bayesian Kalman filter **in the target journal**, so the abstract must position against it: what is
  offered is the test of whether such filters tell the truth about their own uncertainty, not another
  filter).
- **The gap claim is about validity, never about absence** (2026-07-28). The temporal correlation of
  macroscopic currents has carried kinetics since 1973 (Lorentzian spectra: Katz & Miledi 1970/1972,
  Anderson & Stevens 1973), the exact nonstationary two-time covariance was derived and used in
  1980-1981 (Conti et al.; Sigworth), and a covariance likelihood predates MacroR by three years
  (Celentano & Hawkes 2004). Writing "nobody used the correlation" is a one-line kill from any referee
  over fifty. Write instead that **nobody characterised when a method that uses it is valid**.
  Sources, the replacement paragraph and the ARMA-on-residuals precedent (Lei et al. 2020):
  `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`, and §A.10 of
  `docs/bibliography/MacroIR_prior_art_map.md`.

## Scope, and the declaration paper 1 must carry

- Minimal two-state `scheme_CO`, non-stationary protocol, macroscopic currents, likelihood-only. All
  cross-paper; owned in `../_program/decisions.md` §2–3.
- **Paper 1 characterizes the gating-dominated regime** (bands A and B of `../_program/axes.md`), and
  **must say so in band terms**, naming the companion papers for the rest. This is not a hedge; it is
  what closes the flank that adding LSE was meant to close, now closed by a stated scope instead
  (`../_program/program.md` §7). Without it the paper reads as living only where IR wins by
  construction.
- **N_ch 10 … 10⁴.** The floor of 10 is where the microscopic boundary begins, which is why paper 1
  needs the micro attribution anchor below. The three papers' N_ch ranges are chosen jointly
  (`../_program/program.md` §2, still open).

## A-strict is DEAD (reversed 2026-07-23)

**Superseded, kept for rewind.** A-strict said: restrict to the recursive roster R/MR/VR/IR, name the
non-recursive members once in Theory, measure them in no figure. The merge reverses it. **`NR` is
measured and it is in the body**, and the numbers it carries (the overconfidence factor, the
information gap) come back into this paper because there is no paper 2 to hold them.

The original rationale is preserved because half of it survives: the validation machinery is indeed
tested harder on the subtle within-recursive distortions than on NR's gross one, and that is why
`MR` and `VR` stay in a supplement rather than disappearing. What does not survive is the conclusion.
Showing NR is now the paper's job, because the ladder is what makes the cost argument legible to a
reader who is currently using least squares.

**Consequence to sweep:** every "belongs to paper 2" routing of the 87-nat gap and the 10-16×
overconfidence is dead. **Trap while sweeping:** the caption was already refreshed from 87 nats to
10.4 once NR and NMR left the panel (`figures_build_plan.md` §321, §349), so restoring the routing
without checking reinstates a superseded number. And the overconfidence factor is now **10 to 15**,
not 10 to 16, once the non-recursive interval member leaves the panel and the anchor moves to the
freeze — and **recheck it after the 2026-07-31 re-run**, since it was computed against `NMR`, the
build missing the `N·ms` term
(`decisions/D-4_ranking_verdict.md` §3).
**Status 2026-08-01: the re-run has landed** (`figures/data/1f7138b/`, `macro_INR` at 17 cells, plus
`macro_VR`) **and the recheck is still owed.** What the re-run settles is the information distortion
and the bias, recomputed in `decisions/recompute/d3_interval_vs_recursion_2x2.py`. The overconfidence
factor is a **different statistic** — the empirical-over-Fisher standard-error ratio — so nothing above
may be updated from that script. Recheck it with D-4's own recipe against `macro_INR`, and note that
the answer will not be a small correction: `INR`'s distortion is 22 where `NMR`'s pooled with `NR`.

## Fig 1 = the window × recursion lattice: NR, INR, R, IR

**Decided 2026-08-12 (Luciano), built and rendered the same day.** `paper_both/figure_1.Rmd` reads
`COLS_TARGET <- c("NR", "INR", "R", "IR")`: the 2×2 of the two axes, which is the reading that
carries the mechanism (block above, 2026-08-01). The two least-squares columns come out, and the
reason is measured rather than aesthetic. On the figure's own recording `LSE` predicts the same mean
current as `NR` to 8.9e-16 pA and `ILSE` the same as `INR` to 2.2e-16, so rows A and B of each
least-squares column were a copy of its neighbour's and row C was blank; the whole content of the
pair was one constant predictive sd (3.536 pA for `LSE`, 3.559 for `ILSE`, against `NR`'s 0.224 to
2.010), and in the body's two-interval crop "constant" does not even read, only "wider". What is
bought is width: panels go from ~1.0 in to ~1.55 in, and `IR`'s two occupancy rows, their dashed
covariance tie and the two-tone boundary disc are what needed it.

**The fact those two columns carried is not lost, it moved into words**, in two places, and both are
measured: the Theory paragraph that introduces the two arms now says each is the open-loop member at
its own window setting with the gating variance replaced by one constant fitted over the record, and
the Fig 1 caption says the least-squares arms are the open-loop columns with a band of constant
height. The caption was NOT allowed to carry the full sentence: at 326 words the float overflows the
page by ~12 pt, and at 287 it fits with ~14 pt of slack.

**The caption was rewritten the same day, on four corrections from Luciano, and three of them are
about what the figure is FOR.** It said "recorded" where the record is simulated. The reason the
noise is the lowest cell of the design plane was nowhere: it is there so the reader can see that a
window holds several gating transitions (four and five in the two drawn, measured off
`figure_1_simulation.csv`) and that the recording keeps one number for each. The four undrawn
members are now covered by naming what each drawn column REPRESENTS, which is stronger than saying
they are not drawn: least squares is the open-loop column of its own window setting with a constant
band, and `MR`/`VR` are the recursive column carrying `INR`'s prediction and innovation, which is
the structural claim Table 1's caption already makes. And the x axis stopped being the clock: ticks
name the acquisition interval, dotted verticals mark its boundaries in every panel, the duration is
a caption fact. The caption is now 331 words and the float takes a page of its own with no slack
left, verified by rendering page 9 of the build.

**Rejected the same day: R vs IR alone** (Luciano floated it). Four reasons. The 2×2 is what the
Results measure and each margin moves a different moment, so a figure showing one margin under-serves
them; `NR` and `INR` are the two members with literature attribution (Milescu's independent-interval
likelihood and the published `MacroINR`), and dropping them leaves the only mechanistic picture
showing house members alone; the blank row-C cells are the one place the paper draws what open loop
means; and the width no longer buys anything, the `IR` mechanism being fully legible at four columns.

**Files.** `figure_1.Rmd` (roster) over the shared `figure_1_panels.R`, which is unchanged and still
defines all eight columns, so the roster is one vector. No supplement (see the supplement table).
The six-column predecessor is `archive/figure_1_superseded_20260812.{Rmd,pdf}` and its stale caption
file, which still described the 2026-07-21 R/MR/VR/IR roster, is
`archive/Figure_1_caption_superseded_20260812.md`; the caption now lives only in `02_framework.tex`.
The all-algorithm version is `archive/figure_1_all.Rmd`, which writes its own `Figure_1_all.pdf` and
carries no measured claim.

**Superseded, kept for rewind (2026-07-23, the merge):** Fig 1 is the four-column cost ladder
`LSE, NR, R, IR`, the non-recursive members fitting the same grammar, LSE adding a flat global σ̂²
band and a blank row C. Extended to six on 2026-08-05 (`LSE, ILSE, NR, INR, R, IR`, the
instantaneous/interval pair at each of the three levels of the gating description), which is the
version that was superseded on 2026-08-12.

**Superseded, kept for rewind (2026-07-21):** Fig 1 shows the recursive ladder **R, MR, VR, IR**, with
recursion held fixed so the columns vary only the conductance conditioning. That was the A-strict
roster; the merge replaced it with the cost ladder. MR and VR keep a column in the supplement version.

**Superseded, kept for rewind (2026-07-20):** Fig 1 shows *R vs IR only*, on the ground that MR and VR
would be visually identical in everything a single filter step shows except the predicted observable
variance band, making two near-duplicate columns. The rendered four-column figure does not bear that
out. VR's band is visibly narrower wherever the channels gate, and the same predictive variance
divides the gain, so VR's posterior parts from MR's at the second gated interval and its propagated
mean one step later. The distinction is still *also* a quantitative claim and still belongs in the
mechanism/decomposition figure; Fig 1 now shows the path, Fig 2's clouds and the Fig 4 supplements measure it.

**One phrase from the old call has to be retired wherever it was copied:** MR→VR is *not* purely a
variance step and VR→IR is *not* purely a gain step. Changing the variance changes the gain, because
the predictive variance divides it. See `figures_build_plan.md` §F1-2 for the measurement.

## VR must run before the figures, and the arc branches on its sign

- **Three of paper 1's figures depend on the VR column existing.** VR's engine work has landed and
  it **ran for Figure 1 on 2026-07-21**; what is still missing is the *grid* run that Figs 5 and 6
  need (spec: `theory/macroir/notes/vr_variance_form_plan.md`; order: `figures_build_plan.md` §3).
  Which three figures they are has never been enumerated, and `docs/manuscript-drafts/sections/04_results.md` implies more.
- **The branch is resolved: VR came out over-confident, and MORE than MR (2026-07-22).** Measured on
  the Gaussian anchor at N_ch 100 / noise 0.1 / Delta = 0.1 tau, empirical-over-Fisher ellipse-area
  ratios: R 1.32, MR 1.97, **VR 2.18**, IR 1.02 on the kinetic pair, and R 1.09, MR 1.53, **VR 1.77**,
  IR 1.00 on the amplitude pair; the sandwich correction returns all of them to within a tenth of one.
  This is exactly the prediction the falsifier stated (`theory/macroir/notes/vr_variance_form_plan.md`:
  C_ii > 1 for VR, more than MR), so **paper 1's central claim stands**: removing the variance without
  the boundary gain makes the reported uncertainty worse, and the step that recovers calibration is the
  gain. The Results arc no longer needs two branches; write the confirming one.
  One cell so far. The claim is scoped to the regime, so it needs the N_ch x noise grid before it can
  be stated at that scope (`figures_build_plan.md` section 4).

## The band-A result (the former ranking)

- **The ranking is now the band-A column of paper 1's usage map**, not a global verdict. The table is
  in `00_plan.md` §4; its two contested cells (MR's variance direction, IR's corner bound) stay owned
  by `decisions/D-4_ranking_verdict.md`, rescoped to band A rather than rewritten.
- Retired phrasings, so they are not re-copied from the older drafts: *"IR sole survivor"*, *"only
  MacroIR stays calibrated across the practical regime"*, *"MR strawman"* (`../_program/decisions.md`
  §6).

## The micro attribution anchor (what paper 3 owes paper 1)

- Paper 1 runs to 10 channels and reports IR's own degradation there but **cannot attribute it** — all
  four of its methods share the macro closure. `micro_IR` at N_ch 10 / n_sims 10⁴ / noise 0.1, already
  on disk (`../_program/provenance.md`), decides whether that degradation is the interval closure
  (paper 1's subject) or the occupancy closure (paper 3's).
- **One or two annotated cells, not a column** — a full micro column re-opens the roster question the
  three-way split just closed. Paper 1 forward-cites paper 3 as the companion for the full boundary.
- Caution: the anchor exists at **one noise level only**; if the few-channel claim is made across the
  noise fan, more cells are needed. And only that one cell pairs cleanly on n_sims.

## Manuscript

- Head manuscript is `docs/manuscript-drafts/elife_paper.tex`; it is the vessel, every claim owned
  upstream. Superseded drafts stay as history.
- **Article type: Research Article. Settled long before 2026-08-26 and recorded here that day,
  because it was settled nowhere.** It lived only in the first line of `cover_letter.md` ("for
  consideration as a Research Article") and in the author's memory, so `check.sh` item 8 and
  `01_writing_plan.md` §3 both kept reporting it as the open decision `D-1`, and it resurfaced in
  review at least twice. This is rule 3 of `../_program/00_index.md` doing exactly what it exists to
  prevent. Tools and Resources was the alternative and is dead: it obliges benchmarking against
  existing methods as new science and it buys nothing here. What the choice binds: eLife sets **no
  hard word limit** for a Research Article (the 5,000-word main-text figure is an advisory that 55%
  of published articles exceed, `../_program/elife_main_text_length_survey.md`), and the abstract
  norm is the research-article one, median 204 words, not the shorter Tools and Resources median of
  157 (`../_program/elife_abstract_length_survey.md`). Both surveys already assumed this answer.

## The figure set (2026-07-22, revised same day)

**REVISED 2026-08-10: SEVEN body figures, and the last three renumbered.** The τ_int decision panel
entered the body as **Figure 5**, so the corner figure became **Figure 6** and the usage map
**Figure 7**. The live roster is: 1 the filter step on the cost ladder (in Theory), 2 the recovery
clouds at one cell, 3 the calibration cascade in time, 4 the design plane in both moments, 5 the
distortion against the memory left in the least-squares residual, 6 where the calibrated member
departs, 7 the usage map. The numbered list below is the 2026-07-28 state and is kept as history;
where the two disagree, this block is the live one.

*Why a seventh, since it is the one figure of the set with no abstract sentence behind it.* Checked
against the hard rule in `../_program/figure_supplements_decision_guide.md` §6: every other body
figure carries an abstract claim, and this one carries the Introduction's closing paragraph instead,
which makes the offer (measure τ_int on your own record, no ensemble, no known truth) that had no
display item anywhere in the paper. The count is not the constraint: that file's own census of 1,079
eLife Version-of-Record articles puts the median at six body figures with 5–7 the modal band
(n = 631), and volume correlates with the assessment's strength-of-evidence term at ρ = −0.028.

*Two design choices settled with the promotion, and the producer's header owns both.* FOUR members
(`ILSE`, `INR`, `R`, `IR`), the factorial of the two axes, because `R` is the only
recursion-without-averaging member and without it the reader attributes the repair to the rest of the
ladder in a lump. ONE shared y, because the vertical distance between members at one x is the
quantity the figure delivers, and because `IR`'s internal structure on a free axis is the estimation
floor of the statistic (1.043 ± 0.010), so a per-panel scale would draw the bootstrap at full panel
height. The free-scale, own-τ_int reading is Figure 5—figure supplement 2 (supplement 1 until the 2026-08-27 renumber).

*File renames, 2026-08-10, in `figures/paper_both/`.* `figure_6B_four.Rmd` → `figure_5.Rmd`
(`Figure_5.pdf`); `figure_6BB.Rmd` → `figure_5_S1.Rmd` (`Figure_5_S1.pdf`); the old `figure_5.Rmd` →
`figure_6.Rmd` and `figure_6.Rmd` → `figure_7.Rmd`, with their `.html`, `.pdf` and `_caption.md`
siblings; `figure_6B.Rmd` → `figure_5_variant_2members.Rmd` and `figure_6B_all.Rmd` →
`figure_5_variant_8members.Rmd`, neither of them a display item. *Superseded in part 2026-08-27:*
`figure_5_S1` and `figure_5_S2` swapped names and outputs when the two supplements swapped slots,
so the notebook this paragraph calls `figure_5_S1.Rmd` (the own-τ_int companion, ex `figure_6BB.Rmd`)
is now `figure_5_S2.Rmd`. Same for `figure_4_S2` and `figure_4_S3`. NOT renamed, on purpose: the parked
`figure_5_budget*` and `figure_5_IR_*` families, whose 5 is the pre-2026-07-22 numbering and is
already historical, and `figures/archive/figure_6_superseded_20260731/`, whose name is dated.

*Owed on the new Figure 5, both in its producer's header.* The ill-conditioned corner (N_ch 10⁴,
Δ·k_off 1) is not excluded although Figure 4 greys it by a rule it owns (`figure_4_common.R:216-224`,
off-scale **and** κ > 3e4), and those cells alone set the floor of the shared axis: `R` has six below
magnitude 0.9 and `LSE` twelve, every one of them in that corner, against none for `ILSE`, `INR`,
`NR` and `IR`. And the x tick labels collide at four panels across 7 in; `c(1, 2, 5, 10, 20)` fits.

**Five figures in the body, with supplements attached to parents.** The earlier eight-figure list
(the one still described in `figures_build_plan.md` §6, now stale) was superseded within the day (see
next block); the renumbering merged the IR-only map into the R-vs-IR map and demoted three figures to
supplements. eLife imposes no limit ("No limit on display items",
`../_program/elife-author-instructions.md:89`; the "up to 8 figures" at :135 is the Review Article
rule), so this is a density-of-argument decision, not a space one, and the six-figure gate at
`01_writing_plan.md:35` is moot.

**Revised 2026-07-28: SIX body figures.** The region map is the paper's sell and cannot be a
supplement. `check.sh` check 6 tests `N_CAP >= 5` and must be raised, or the done-oracle stays green
on a paper short one body figure.

The body:

1. the filter step along the cost ladder (**LSE, NR, R, IR**)
2. recovery clouds, the ladder at one cell
3. the calibration cascade in time
4. **the design space, split by moment** (`Figure_4_bias.pdf` and `Figure_4_distortion.pdf`, both
   built), plus the grouped-MLE clouds that give the maps a referent. Splitting by moment rather than
   by family is deliberate: splitting by family separates LSE from IR and kills the headline contrast
5. **the information budget per parameter**, across methods. Rates live in the deterministic
   transient's shape, so LSE gets them and recursion buys little; the amplitude pair needs the
   fluctuations, and **LSE holds `unitary_current` Fixed so it cannot estimate `i` at all**. This
   shields the thesis: where the information is nearly equal, the only variable left is the honesty of
   the error bar. Must use the distortion-corrected covariance. **Scope caveat obligatory:** all of
   this is the non-stationary regime; under stationarity the mean is flat and the balance inverts
6. **the usage map** (`figure_6.Rmd`, built 2026-07-26/27; `Figure_6.pdf`, `Figure_6_regions.pdf`,
   `Figure_6_frontiers.pdf`, `Figure_6_caption.md`). One plane, N_ch × instrumental noise, partitioned
   by four boundaries into the five regions named in `../_program/nomenclature.md`, with the three
   real recording configurations (excised patch, whole cell, oocyte) placed on it. It is the only
   figure in the set that is not a heatmap: the cells are the input and the **lines are the output**.
   Read the notebook header before writing about it; it is more current than any planning file.

**Owed on Figure 6, and none of it is redrawing:** lift the ~200 lines of commented derivation for the
three preparations into a citable Methods table, one source per row (the sourcing is in
`docs/bibliography/recording_configurations/`, 88 files, **currently untracked**); write the figure
into `04_results.tex`, which has no region map at all; and decide whether the lower vertex, currently
extrapolated below N_ch 10, is pinned by the ten cells the notebook's own open item asks for
(N_ch 2 and 5 at noise 0.1 to 10). That last one is `../_program/program.md` §2.

**THE SUPPLEMENT SET, SETTLED 2026-08-11/12 (Luciano), one parent at a time.** Eleven, declared
inside each parent's `figure` environment and every one of them cited in the main text in ascending
order per parent, which was checked file by file and is the one hygiene item the eLife census
measures (4.3% of reviews complain about uncited or out-of-order supplements). The list below
supersedes everything above it in this block.

| Parent | # | What it carries |
|---|---|---|
| 1 filter step | — | none. The declared one drew samples 0–4 where the body draws 3–4, so it CONTAINED the body figure: a second crop, not a supplement, and never cited. A parent at zero is the eLife norm. |
| 2 clouds | 1 | Mahalanobis Q-Q, the sandwich against the empirical multivariate distribution. MOVED here from Figure 4: its citing sentence is in this subsection, it runs on one cloud cell, and it uses no part of `figure_4_layout.R`. |
| 3 cascade in time | 3 | per-step information vs score variance on k_on/k_off, the same on i/N_ch, and the checks the parent has no room for. The first two are ONE figure split by page height (8 rows = 11.2 in against a 9.2 in text block), not two answers. |
| 4 plane | 3 | the distortion factored into per-sample and correlation; its size against its shape; and the same plane for the four rungs of the recursive family. |
| 5 memory | 3 | each member against its OWN residual memory; the collapse through r = S̃/N_ch; and τ_int over the design plane. |
| 6 IR's corner | — | none, settled 2026-07-31 and unchanged. |
| 7 usage map | 1 | the distortion-corrected standard error over the plane: the continuous field under the map's solid boundary. |

**What left, and on what test.** Every candidate was asked for the one sentence in the body that
cites it. The line version of the sample/correlation split went out because it broke this parent's
grammar (lines against N_ch where every other supplement here is the plane) and drew two members
where its replacement draws six. `lag_kappa` was SPLIT: τ_int survives as Figure 5's third
supplement and the condition number was dropped outright, since its only job in the paper is the
masking rule κ > 3e4, which Methods states in prose. Out with no sentence to write: `se_kappa`
(its side A duplicates the standard error and it drops k_on, one of the four directions the usage
map's boundary is built on), `coverage` and `qq_grid` (both ask what Figure 2's supplement already
answers, and only they would add the plane, which no sentence claims), `other_parameters` (two
pages of completeness), `eigendirections`. `Figure_4_supplement_residuals.pdf` turned out to be a
dead file whose producer was replaced by `lag_kappa` on 2026-07-31.

**The interval-pairs page is OUT, decided 2026-08-12, and the reason is a measurement rather than
a preference.** It was the only candidate that put the acquisition interval on an axis of its own,
and it was built to show that each pair closes as the window shrinks: at a hundredth of a time
constant the interval average and the instantaneous sample differ by about two per cent, so the two
members of a pair must coincide, a limit the implementation has to satisfy and one that does not
come from the code under test. Measured at the shortest swept window, noise label 0.1, on the
per-sample distortion, the premise fails and it fails asymmetrically. `R`/`IR` close to $0.07\%$ on
the closing rate at ten channels and their gap shrinks with the window, which is the expected
behaviour. `NR`/`INR` close to $8.3\%$ there and their gap does not go to zero. `LSE`/`ILSE` sit
flat at about one per cent with no trend, which is what two members that differ only in the mean
model look like.

So the page cannot be read as a validation: a test that does not separate "both right" from "both
wrong" certifies nothing, and at ten channels the `R`/`IR` pair closes partly because `IR` is itself
departing there. Two things were learnt on the way out and they are worth more than the page was.
The catastrophe of `NR` is exactly the omitted `N*ms`: the term is added only under
`if constexpr (variance::value && averaging::value > 0)` (legacy/qmodel.h:4585-4618), `ms` is the
within-window variance of ONE channel and is identically zero for an instantaneous observation, and
`NR`'s per-sample distortion in the channel-number direction runs 1.12, 1.15, 22, 569, 3020,
1.5e4, 3.5e4 across the seven windows at ten thousand channels while `INR` stays at 1.00. That
finding survives, in the Figure 4—figure supplement 1 caption. And what sorts the closure is not the
parameter, since the two pairs close equally in the channel-number direction at few channels (2.2%
against 2.0%) and differ by 120-fold on the rate; it appears to be the recursion, which absorbs a
variance misspecification through the update where an open-loop member passes it straight into the
score. That last is a conjecture, measured at one noise level, and is not in the manuscript.

**Still open, blocked on Luciano.**  The figure for the expected
Fisher not being positive definite, from the 2026-08-09 audio, which nothing on hand answers:
`figure_4_gaussian_vs_numeric_fisher.Rmd` asks a different question and its own header says it
cannot reach the two members that matter.

**One roster question left inside a declared supplement.** Figure 4's supplement 2 carries `ILSE`,
`NR`, `INR`, `R`, `MR`, `IR`: six members, but not the body figure's six. `MR` is deliberate, the
paragraph it serves names it. `LSE` being absent is not settled, and the abstract's "no rescaling
repairs it" covers that member.

## Figure 4 is R against IR, and the fusion (2026-07-22)

**Decided (Luciano, 2026-07-22).** The former Figure 5 (IR alone: where IR stops being faithful) and
the planned R-vs-IR decision map were **merged into one body figure**, not kept as figure plus
supplement. The reason is that the R-vs-IR map *contains* the IR-only map: the IR-only figure's bias
and distortion blocks were a strict subset of the merged figure's IR rows, so keeping both showed the
same panels twice. IR's own failure is still located, in its own rows, with the same numbers; the
promise of `00_plan.md` §1 (paper 1 locates IR's own failure) is kept by the paper locating it, not
by a figure dedicated to it. The old `figure_5*` are archived under
`figures/archive/paper_superseded_20260722/`.

**Layout.** Columns are N_ch (the reader enters by their own channel count, which is given; Δ is a
design choice and noise is instrumental). Rows are parameter (outer) then algorithm (inner), so each
parameter's R and IR panels are adjacent for a vertical glance. Two versions are written: `long`
(both k_off and N_ch in both blocks, the body figure, `Figure_4.pdf`) and `short_variant` (one
direction per block: bias on N_ch, distortion on k_off, the two directions that carry each moment's
failure). Block A keeps the grouped-MLE clouds because they are the referent a distortion map needs;
the MLE-clouds are NOT demoted (Luciano: they carry the figure's intuitive force). The design points
moved to noise 0.1, N_ch 10/10/10000, because R exists at neither of the old figure's noise levels
(0.05, 0.5) and N_ch 20 is no longer a column.

**What the merged figure measures (noise 0.1, medians over the 7 intervals):**

- **More channels do not rescue R.** Its k_off distortion is a flat floor, 1.37 / 1.44 / 1.44 / 1.42
  across N_ch 10 / 100 / 1000 / 10⁴, while IR converges to one (1.32 / 1.10 / 1.00 / 1.00). Its N_ch
  direction crosses from conservative to over-confident, 0.87 / 1.10 / 1.23 / 1.26, and keeps
  growing, while IR converges (0.85 / 0.94 / 0.99 / 1.00). The two agree at N_ch 10 and diverge from
  there, so the 10⁴ column is the most informative, not the emptiest.
- **New, and first-moment: R carries a large, flat bias in N_ch**, ~0.13–0.15 log10 (≈ 35–40%) across
  the four decades, not shrinking with channels; IR's is ~0. Both algorithms' k_off bias is ~0. This
  was not in the old ranking.
- **The correlation is the mechanism** (supplements 3 and 4). The distortion splits total ≈ sample ×
  correlation; R's sample part goes to one with N_ch exactly as IR's does (the per-sample Gaussian
  approximation becomes exact for both), while R's correlation part stays elevated and *rises* with
  N_ch (k_off 1.35 → 1.47). So R and IR do not differ in per-sample fidelity; they differ in temporal
  decorrelation, which is what conditioning on the interval buys. More channels push sample down and
  correlation up by the same amount, leaving the total flat.
- **The distortion is diffuse for R, confined for IR** (supplements 1 and 2). R departs by ≥15% in
  52 / 68 / 27 / 38% of cells for k_on / k_off / i / N_ch; IR in 2 / 11 / 0 / 4%. Only the
  instrumental-noise direction is clean for both. Conditioning localises the distortion, it does not
  merely shrink it.
- **Noise cures R only at few channels** (supplement 4). R's correlation distortion starts lower at
  higher instrumental noise but re-emerges to ~1.3–1.5 by N_ch 10⁴ at every noise level, because the
  gating variance it distorts grows with N_ch. This is the quantitative form of the band-B/band-C
  boundary the paper must state.

**The noise axis stays RAW, not noise/N_ch (Luciano).** Dividing by N_ch would flatten the staircase
by construction and assume the very N_ch-scaling the figure can measure; the slope of the calibration
frontier in log-log *is* the exponent. The grid is therefore declared as (N_ch, noise) pairs, not a
product, and no-data panels render grey (a CI-shrunk calibrated cell is white, so grey ≠ white is
what distinguishes "measured and fine" from "never run").

**The pathological cell, carried openly.** macro_R, N_ch 10⁴, Δ·k_off = 1: bias in N_ch reads −2.78
log10 and distortion 0.46, both off scale, reproduced across noise 0.1/1/10 and across nsim 1000 and
10000. Not a broken cell and not R-specific ill-conditioning (κ of the Gaussian Fisher covariance is
75272 for R and 99459 for IR there, so IR is worse and its bias is ~0). The reading: the bias is H⁻¹
times a miscalibration, IR's numerator is ~0, R's is not and gets amplified along the near-null
direction. Sign and existence real, magnitude not a physical bias. **Open decision, now spanning
Figure 4 and its four supplements: grey those cells as unidentified, or plot bias/SE.**

**The R coverage correction.** `figures_build_plan.md` §2 said R and MR exist at one noise level only
on the Gaussian anchor, so the roster map "cannot be built". That is wrong: both have noise 0.1/1/10
across all four decades on the G anchor, plus an N_ch-scaled level (noise = 10·N_ch, constant
relative noise r = 10). What is exclusive to IR is the finer N_ch grid and the sub-0.1 noise levels.
So Figure 4 was not blocked, only incomplete.

**The pending R/IR grid, retargeted by measurement.** To fill the calibration-frontier diagonal at
constant relative noise r = noise/N_ch, the cells to run are **r = 0.01, 0.1, 1, NOT 0.1, 1, 10**:
r = 10 is already answered and past the transition (the `82b956f` probe measures R at 1.02–1.06 there
in all four decades, and it is not underpowered — its bootstrap halfwidth is 0.04 log10 against the
0.158 a 1.44 distortion would need). The crossover sits between r = 0.01 and r = 0.1. Six cells per
algorithm: N_ch 100 at noise 100; N_ch 1000 at 100 and 1000; N_ch 10⁴ at 100, 1000, 10⁴. Dispatch
must pass `GROUP_SIZE="10 100"` (the dispatcher default is `1 10 100`, and no existing run has group
1); the same fix applies to the VR grid commands in `figures_build_plan.md` §4, which omit it.

**Q-3, and the orphan block deleted 2026-07-28.** Until today this file said in two places that the
Fisher-to-zero result **stays in the body as Figure 4**, and in one place that it is demoted to
`Fig 3—figure supplement 1`. The demotion is the later of the two and `00_plan.md:220` records the
sequence ("RESOLVED 2026-07-22, then revised the same day"), so **the supplement reading wins** and
the body-figure paragraph has been removed. It also referred to "Figure 4", a slot the current set
assigns to the design space.

What survives from the deleted paragraph, because it is the real argument and it survives demotion:
the figure carries two readings of one measurement, the information level F_t (*when* each parameter
is measured) and the ratio J_t/F_t (whether that information is honest), and the second cannot be
posed without the first, so the figure does not split. The claim is a **shape** and a shape is shown
or lost; a supplement is still a shown, legended figure.

**Carry the caveat with it:** that result holds for all four rungs and therefore **discriminates
nothing** within paper 1's roster. It is a property of the macroscopic observable, not of the
closure, and it must be presented as scope rather than as a finding about conditioning. It answers an
open question of Del Core and Mirams 2025, which is why it earns the space.

**A constraint on any future demotion, not a style choice.** eLife has no independently numbered
"Supplementary Figure". A supplementary *figure* is always `Figure N—figure supplement M`: attached
to a parent, numbering restarting at 1 under each parent, and no limit on how many (eLife's own
`elife-template.tex:209`, implemented in `elife.cls:316-412`, and visible in any published article,
e.g. elifesciences.org/articles/78075/figures, which renders
"Figure 2 with 9 supplements"). The only independently numbered supplementary items are
`Supplementary file 1, 2, …`, which are files and not figures: no figure legend, not in the figure
list. So "send it to supplementary" always means "make it a supplement *of* some figure", and that
parent relationship becomes part of the argument. **Do not re-derive this**: the transcription in
`../_program/elife-author-instructions.md` does not mention figure supplements at all, and that
absence is not evidence — the mechanics live in the template and the class.

**What this does NOT cost, corrected 2026-07-22.** An earlier version of this paragraph said the
supplement must be **named in the parent figure's legend**, and a demotion was priced as if it meant
hand-writing a sentence into each parent legend. That is wrong, twice over. First, `:209` is the body
text of eLife's own demo `\figsupp` call and says supplements *"should be referred to"* in the parent
legend: an editorial recommendation, not a requirement. Second, the mechanical part is not the
author's to write. `\figsupp[short]{legend}{graphic}` sits **inside the parent's `figure` environment**
(that is what `\thefigure` numbers off) and **auto-emits** the `Figure N—figure supplement M.` line
under the parent's caption; `elife.cls:400` gives an explicit opt-out, `\figsupp[none]{...}{...}`,
which suppresses it entirely. And what appears on elifesciences.org is eLife *production's*
typesetting of the accepted manuscript, not the author's LaTeX output, so a published parent caption
that names none of its supplements (checked against a real article) is consistent with all of the
above. **Demotion therefore costs no legend prose in the parent.** What it does cost is the parent
relationship itself, which is the part that belongs in the argument.

## The abstract, introduction and discussion skeletons (2026-07-28)

Recorded here because they are decisions about what the paper argues; the prose belongs to
`docs/manuscript-drafts/sections/00_abstract.md`, `docs/manuscript-drafts/sections/01_introduction.md` and `docs/manuscript-drafts/sections/05_discussion.md`.

**Abstract.** Macroscopic currents are generally analysed by least squares, which ignores the
information in the temporal correlation. Methods that use that correlation have little uptake,
probably because their validity was never characterised. This paper characterises it. The only method
without bias or distortion above **[threshold pending, D-J]** is MacroIR; and where the records carry
no autocorrelation, least squares is just as good.

**Amended 2026-08-01, and the amendment is a division of labour, not a ranking** (Luciano, audio
2026-08-01 08:47). "The only method without bias **or** distortion" fuses two failures that the
measurement separates and that a referee will separate anyway. The three arms, each with what it can
and cannot deliver:
- **least squares** carries no meaningful bias in what it fits (0.001 log10), but it holds the unitary
  current Fixed, so the amplitude pair reaches it only through the product, and its own error bar is
  distorted by an order of magnitude (13.3 on `k_off`, essentially all correlation);
- **`INR`** returns the amplitudes unbiased (0.003) at a wider error than the recursive members, and
  its reported error is distorted by the same order (22.2, all correlation). So it buys the estimate
  and not the interval;
- **`IR`** is the only member whose reported error can be believed (1.01).
Two different failures — *cannot estimate it* and *cannot tell you how well it estimated it* — and the
abstract should carry both, in that order. This does not change what was swept or measured; it changes
which sentence the results are packed into. The live prose is
`docs/manuscript-drafts/sections/00_abstract.tex` and it has its own rules file; do not rewrite it from
this paragraph without reading them.

**Introduction, five moves.** (1) Independence of residuals is the foundational assumption of least
squares, and violating it costs you the degrees-of-freedom count. (2) Markov chains model temporal
dependence while still yielding constants universal to the whole record, and they tie biophysical
structure to observables, so the parameters mean something. (3) The ladder of methods, as in the
roster above. (4) Why least squares still reigns: it lets you see, point by point, whether the
prediction is good, which is visually convincing; recursive methods use the previous datum to predict
the next, so they follow the data closely, hold no surprises, and give you no way to tell a working
method from a bug. **There is no easy validity criterion, and that is what this paper supplies.**
(5) Autoregressive alternatives do not solve it: an ARMA error model is stationary by construction
while the gating covariance tracks the mean current and restarts at every jump, it contains no channel
count and no unitary current, and its timescales are free where the Markov model ties them to the
same rates that generate the mean. See
`docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md` for the sourced version and for
Lei et al. 2020, who tried exactly this on ion-channel data.

**Discussion, four points.** (1) MacroIR is calibrated over almost the whole measured plane, so one
could simply always use it. (2) Where it fails: few channels, and telegraph (non-Gaussian) noise.
**Flag: "telegraph" appears nowhere under `papers/`, no simulator capability exists on any freeze
commit, and no cell has been run.** This is the only new claim in the set with zero data behind it;
either run it or demote it to an argument with a citation. (3) Why the intermediates fail — **rewritten
2026-08-01, the old version was measurably false.** It used to read "you need the double conditioning
at both interval ends, conditioning on the start alone contributes nothing". Conditioning on the start
alone contributes the **first moment**: it removes the amplitude bias entirely (`INR` 0.003 log10
against `NR` 0.099 in N_ch) and drives the per-sample distortion to one. What it does not touch is the
**second moment**, and that is the whole point: `INR`'s distortion is pure correlation, ~22 on `k_off`,
flat in N_ch and unmoved by the window. Recursion is what removes it, and only recursion conditioned
on **both** ends removes it completely — one end (`MR`, 1.94) is worse than none (`R`, 1.46).
So the sentence to write is *the window restores the mean, the recursion restores the variance, and
the recursion has to see both ends of the interval because that is where the information crosses from
one interval to the next.* Numbers: the roster block above, and
`decisions/recompute/d3_interval_vs_recursion_2x2.py`.
(4) **The bootstrap concession, and it has to be made in print** (Luciano, audio 2026-07-29). An
experimenter who does not need the conductance can fit by least squares and get an error bar by
resampling, and there is nothing wrong with that: an empirical error incorporates sources of
variability an internal one does not, so in that sense it is better. This concedes the "trustworthy
error bar" half of the ordering result to anyone willing to resample, and it must be conceded before a
referee does it — the route is in this literature already (Moffatt & Hume 2007 bootstrap; Stepanyuk
2011 bootstraps a filter it calls Kalman). What survives the concession, and is the argument to make
there: a resampled interval cannot separate **variability between recordings** from the **stochastic
error of the method on one recording**, so it cannot tell you whether two currents differ because the
rates differ or because the channel counts do. That separation needs an internal, calibrated error, it
is the second of the three advantages in the 2026-07-29 audio, and it is a capability statement — see
the closer's own scope note in `docs/manuscript-drafts/sections/00_abstract.tex` §1(g).

**The risk to fix before the abstract is written, not after review.** The region map's largest region
is sold on the unitary conductance and the channel count, and `introduction.md:43` already concedes
both to non-stationary fluctuation analysis while arguing that NSFA does not return a kinetic scheme.
The headline therefore leans on the parameter the existing defence gave away, and the 1980-1981
literature makes it worse. What is actually new is the **combination**: a calibrated joint estimate of
rates *and* amplitudes with honest intervals, from a **single non-stationary record**, where the
covariance methods needed ensembles of 256 to 504 repeated sweeps and produced no uncertainty
statement at all, and where NSFA's own practitioners report the unitary current is close to the only
parameter it recovers reliably. "More conductance information" on its own is not defensible.

## The R and Python ports (2026-07-28)

An R and a Python implementation of MacroIR were extracted from the C++ and cross-checked against it
once, via `tools/cross_language_check.py`. This is an availability claim and it belongs in Code
Availability (`../_program/carve_plan.md` owns the topic). **It is currently the weakest sentence in
the submission**: a usability claim inside a paper about when methods are valid, backed by a single
cross-check. Run the cross-check across the design grid and record the tolerance per cell before it
goes in.

## Open decisions

Live in `00_plan.md` §8, and they are labelled **`Q-n`** since 2026-07-21 — the `D-n` labels belong to
the manuscript-production briefs in `decisions/`, and the two registers used to collide. In brief:
**Q-1** (does NR stay — for paper 1 the answer is no, it moved, but confirm); **Q-2** (MR and VR main
text or supplement — "strawman" is retired, so re-decide on the map footing); **Q-3** (resolved: the Fisher-to-zero result is Fig 3—figure supplement 1); the figure count is
settled at five body figures plus supplements (the six-figure gate is moot); the title (three live
versions, `docs/manuscript-drafts/sections/README.md`).

**Q-5, opened 2026-08-01: does `INR` join the body roster?** The corrected run makes it the only
member that isolates the interval-window margin, and `LSE, NR, R, IR` is three of the four cells of a
2×2 whose fourth cell is `INR` (see "What this paper is" above). Adding it costs one column in Figs 1,
2 and 4 and turns the ladder reading into the factorial reading; leaving it in the supplement keeps the
cost-ladder narrative and forces the mechanism argument to be made in prose against a supplement panel.
Interacts with Q-1 and Q-2, which is why all three should be decided in one pass. **Luciano's call.**
