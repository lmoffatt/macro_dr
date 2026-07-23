# Paper 1 (method) — decision log

> Updated: 2026-07-20. Paper 1's **own** settled decisions. Everything that binds more than one paper
> moved to `../_program/decisions.md`; this file cites it and does not restate it. If a decision here
> turns out to bind papers 2 or 3 too, it belongs in `_program/`, not here (`../_program/00_index.md`
> rule 2).
>
> Open decisions live in `00_plan.md` §8.

## What paper 1 is

- **The method paper.** Given that you will compute a likelihood, **what must it condition on?** This
  is the interval-likelihood closure (`../_program/program.md` §1). The other closure is paper 3; the
  prior question of whether you need a likelihood at all is paper 2.
- **Roster: `R`, `MR`, `VR`, `IR`.** A monotone progression in how much of the interval structure the
  likelihood uses. Cross-paper roster and why each method sits where it does: `../_program/decisions.md`
  §2 and `../_program/program.md` §1.
  - **`NR` and `NMR` are not in paper 1.** NR moved to paper 2 (a cheap non-recursive likelihood is
    that paper's subject, and its link to Milescu 2005 / QuB); NMR was dropped from the program (no
    literature attribution, no mechanistic role).
  - **`VR` is in**, displayed as "Variance Recursive" (2026-07-21, `figure_1.Rmd`; overrule there if
    the paper prefers another wording). It is the control that turns the "MR's problem is the gain,
    not the variance" claim from algebra into measurement. **Do not describe it as "MR→VR changes only
    the variance, VR→IR only the gain":** the predictive variance divides the gain, so the variance
    step moves the update too, and IR's *total* predicted variance is algebraically equal to MR's
    (`figures_build_plan.md` §F1-2).
- **Literature anchor after the exclusions: `R`.** Moffatt 2007; Münch 2022 (a published Bayesian
  Kalman filter). The frame is "the published recursive filters treat each sample as instantaneous;
  here is what conditioning on the interval buys, and which half of the mechanism does it." Paper 1
  does not need NR to anchor itself.

## Scope, and the declaration paper 1 must carry

- Minimal two-state `scheme_CO`, non-stationary protocol, macroscopic currents, likelihood-only. All
  cross-paper; owned in `../_program/decisions.md` §2–3.
- **Paper 1 characterizes the gating-dominated regime** (bands A and B of `../_program/axes.md`), and
  **must say so in band terms**, naming the companion papers for the rest. This is not a hedge; it is
  what closes the flank that adding LSE was meant to close, now closed by a stated scope instead
  (`../_program/program.md` §7). Without it the paper reads as living only where IR wins by
  construction.
- **N_ch 10 … 10⁴.** The floor of 10 is where the multinomial boundary begins, which is why paper 1
  needs the micro attribution anchor below. The three papers' N_ch ranges are chosen jointly
  (`../_program/program.md` §2, still open).

## The non-recursive members: named, not measured (A-strict, 2026-07-20)

**Decided (Luciano).** Paper 1 restricts to the recursive roster R/MR/VR/IR. The non-recursive members
`NR`, `NMR` are **named once in Theory** — to locate paper 1's ladder as the conductance axis at
recursion=on — and **measured in no figure**. The recursion axis is established prior work (Moffatt
2007, Comm Biol 2025); paper 1 holds it fixed and varies the conductance conditioning.

Rationale is concept-first, not a numbers pitch: that NR is dramatically worse is old and roughly
obvious, so re-showing it is not paper 1's job. The validation machinery is tested *harder* on the
subtle within-recursive distortions (catching MR's ~1.5× miscalibration and attributing it) than on
NR's gross one, so restricting to the recursive family is where the method proves it works. The 87-nat
gap and the 10–16× overconfidence are NR/NMR results and belong to paper 2.

## Fig 1 = the four-column ladder R, MR, VR, IR

**Decided (Luciano, 2026-07-21), and built the same day.** Fig 1 (the filter step, no statistics)
shows the recursive ladder in order: **R, MR, VR, IR**. Recursion is held fixed across the four, so
the columns vary only what the interval-averaged conductance is conditioned on. Files:
`figure_1.Rmd` (roster) over the shared `figure_1_panels.R`; the six-algorithm version is
`figure_1_all.Rmd`, which writes its own `Figure_1_all.pdf` and carries no measured claim.

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
  Which three figures they are has never been enumerated, and `results.md` implies more.
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

## The figure set (2026-07-22, revised same day)

**Five figures in the body, with supplements attached to parents.** The earlier eight-figure list
(the one still described in `figures_build_plan.md` §6, now stale) was superseded within the day (see
next block); the renumbering merged the IR-only map into the R-vs-IR map and demoted three figures to
supplements. eLife imposes no limit ("No limit on display items",
`../_program/elife-author-instructions.md:89`; the "up to 8 figures" at :135 is the Review Article
rule), so this is a density-of-argument decision, not a space one, and the six-figure gate at
`01_writing_plan.md:35` is moot.

The body:

1. the filter step along the ladder (R, MR, VR, IR)
2. recovery clouds, the ladder at one cell
3. the calibration cascade in time
4. **R against IR over the design space** (bias and information distortion maps, plus the grouped-MLE
   clouds that give the maps a referent)
5. the design trade-off (the covariance figure, ex-8)

Attached supplements, each declared inside its parent's `figure` environment (eLife mechanics:
`reference_elife_figure_guidelines` memory; naming the supplement in the parent legend is a SHOULD,
auto-emitted, not a must):

- **Fig 3—figure supplement 1**: per-step Fisher profiles across the ladder (ex-Figure 4, the
  Fisher-to-zero result). Demoted because it discriminates nothing within the roster: the
  when-each-parameter-is-measured shape holds for all four rungs, so it is a property of the
  macroscopic observable, not of the closure. Q-3's "stays in the body" is thereby reversed; the
  reason Q-3 gave (a shape is shown or lost) survives demotion, since a supplement is still a shown,
  legended figure. Built.
- **Fig 4—figure supplement 1**: distortion-induced bias, all five parameters, R vs IR. Built.
- **Fig 4—figure supplement 2**: information distortion, all five parameters, R vs IR. Built.
- **Fig 4—figure supplement 3**: the distortion decomposed into sample and correlation, on k_off,
  R vs IR (ex-Figure 6, the decomposition). Built.
- **Fig 4—figure supplement 4**: sample and correlation vs N_ch, k_off and N_ch, resolved by noise
  (lines, not maps). The mechanism as a shape. Built.
- **Pending, and probably a child of Fig 2 not Fig 4**: MR and VR added to the design-space map. It
  generalises Fig 2's one-cell ladder measurement, so the parentage is on that claim, not on the
  visual grammar. Blocked on the VR grid run.

The old S1–S5 still carry the independent-numbering that eLife does not have and must be re-parented
in the same pass: S1 → Fig 1—figure supplement 1 (from `figure_1.Rmd`), S3 (corner) → a Fig 2
supplement, S4/S5 (dlik-fed) → Fig 3 supplements. Not yet done.

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

**Q-3 is resolved: the Fisher-to-zero result stays in the body**, as Figure 4, and the reason is not
that it is striking. The figure carries two readings of one measurement — the information level F_t,
which says *when* each parameter is measured, and the ratio J_t/F_t, which asks whether that
information is honest — and the second cannot be posed without the first, so the figure does not
split. The deciding test was whether the text could carry the result without the figure: it cannot,
because the claim is a *shape* (the information about k_on, i and N_ch dies the moment the agonist is
removed while k_off stays informative through the decay), and a shape is shown or lost.

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

## Paper 1's open decisions

Live in `00_plan.md` §8, and they are labelled **`Q-n`** since 2026-07-21 — the `D-n` labels belong to
the manuscript-production briefs in `decisions/`, and the two registers used to collide. In brief:
**Q-1** (does NR stay — for paper 1 the answer is no, it moved, but confirm); **Q-2** (MR and VR main
text or supplement — "strawman" is retired, so re-decide on the map footing); **Q-3** (resolved: the Fisher-to-zero result is Fig 3—figure supplement 1); the figure count is
settled at five body figures plus supplements (the six-figure gate is moot); the title (three live
versions, `title.md`).
