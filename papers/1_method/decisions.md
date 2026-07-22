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
mechanism/decomposition figure; Fig 1 now shows the path, Fig 6 measures it.

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

## Paper 1's open decisions

Live in `00_plan.md` §8, and they are labelled **`Q-n`** since 2026-07-21 — the `D-n` labels belong to
the manuscript-production briefs in `decisions/`, and the two registers used to collide. In brief:
**Q-1** (does NR stay — for paper 1 the answer is no, it moved, but confirm); **Q-2** (MR and VR main
text or supplement — "strawman" is retired, so re-decide on the map footing); **Q-3** (where the
Fisher-to-zero result goes); the figure count (`01_writing_plan.md` gates at six while the count is
open); the title (three live versions, `title.md`).
