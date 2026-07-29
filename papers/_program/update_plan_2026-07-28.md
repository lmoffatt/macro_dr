# Update plan, 2026-07-28

> **This is a work order, not a fact owner.** It settles nothing. Every decision it points at is
> settled in the file named beside it. Delete this document when the last batch lands.
>
> Built from: the voice notes of 2026-07-28
> (`program/source-notes/audios/audios/Chat de WhatsApp con MacroIR 13/`, four files, 15.08 to 15.33),
> the prior-art work of the same day (`docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`),
> and a twelve-agent audit of all 36 planning documents (workflow `wf_9cd9f5ab-91e`), whose line
> numbers were then re-verified against the working tree.

> ## STATUS 2026-07-29 — most of this is done, and two batches are superseded
>
> **Done:** the five questions of §1 are answered (`../1_method/decisions.md`, `decisions.md`).
> Batch 1 (programme spine) and Batch 3 (paper hub) landed in `6026587`, with the D-4 verdict
> recomputed on the freeze and its scripts committed. The `433ed13` anchor call, the NMR reopening and
> its resolution to the supplement landed in `b861ad8` and `109a030`. The bibliography migration and
> the `\graphicspath` repoint (both listed below as build blockers) are done.
>
> **Superseded:** batches 6a and 6b. The eight section plans were not edited into the new frame; they
> were **consolidated into `../1_method/SPINE.md`** and archived under
> `../1_method/archive/section-plans-20260729/`. The reason is in §0: the manuscript is ahead of the
> plans, so the right move was to write one spine from the `.tex` rather than update eight pre-merge
> files. References to `SPINE.md (Section)` below are the automatic repoint and read a little oddly;
> they mean the corresponding block of that file.
>
> **Still open:** batch 2 (`machinery.md`, `axes.md`), batch 5 (figure reality), batch 7 (routing),
> and every item in §3 and §4 that is not struck above.

## 0. The headline: the documents are behind the work, not ahead of it

The audit's implicit premise was that the plans lead the code. Checked, and it is backwards
everywhere:

| Thing the plans call unbuilt or absent | What is actually on disk |
|---|---|
| the region map (D-B), described as a figure to design | `figures/paper_both/figure_6.Rmd`, 53 KB, plus `Figure_6.pdf`, `Figure_6_regions.pdf`, `Figure_6_frontiers.pdf`, `Figure_6_caption.md`. Encoding decisions dated 2026-07-26 in its header, scaling ledger dated 2026-07-27 |
| the recording-condition overlay, not mentioned anywhere | ~200 lines of commented derivation in `figure_6.Rmd` chunk `regions`, backed by `docs/bibliography/recording_configurations/` (88 files, `SOURCES.md`, `EXPRESSION_LEVELS.md`), **all untracked** |
| roster R / MR / VR / IR | `figure_1.Rmd:37` and `figure_3.Rmd:55` both read `c("LSE","NR","R","IR")`; `figure_4_common.R:27-28` maps the four |
| Figure 4 as one figure | `figure_4_bias.Rmd` and `figure_4_distortion.Rmd` exist: the by-moment split M-4 asks for |
| four Fig-4 supplements | eleven exist |
| Figure 5 unbuilt (`results.md:247`, `provenance.md:191`, `CONTINUE_HERE.md:30`) | `Figure_5.pdf`, `figure_5.Rmd`, `Figure_5_caption.md`, `Figure_5_supplement_1.pdf`, plus four IR-affine variants |

**Consequence for how to do this work.** Most of the update pass is transcription from the notebooks
into the documents, not decision-making. Read `figure_6.Rmd`'s header and `Figure_6_caption.md` before
writing anything about the region map: they are more current than every planning file.

**The second headline.** The manuscript is finished, 11,747 words across seven `sections/*.tex`, and it
is finished in the retired frame. `check.sh` reports 7 pass / 1 fail / 1 warn, and the one fail is a
dead graphics path, so nothing in the project can currently detect that the deliverable and the
decisions have diverged.

## 1. Five questions to answer before any file is edited

These gate batches 2, 3 and 6b and they are one conversation.

1. **Does `R` sit in the body or the supplement?** M-2 as recorded puts it in both.
2. **Which anchor carries the verdict**, `433ed13` or the `1c2ae6f` Gaussian anchor? `D-4_ranking_verdict.md`
   computes on the first, the body figures run on the second, and D-J's number must come off one dataset.
3. **Which distortion cutoff, and in which units?** 1.1, 1.15 and 1.5 all live in the repo. `machinery.md:231`
   admits it. The built Figure 6 already uses **1.15** and glosses it as "a 7% error on the reported
   standard deviation" (a variance ratio); `machinery.md` §8 calls 1.1 "about 10% error", which reads it
   the other way. Three of the five regions are defined by this crossing, so it stops being housekeeping.
4. **Q-3, body or supplement?** `1_method/decisions.md:143-148` demotes the Fisher-to-zero result and
   says so; `:238-241` in the same file still asserts it stays in the body. One is an orphan.
5. **Six body figures or five?** D-B adds the region map. Nine documents plus `check.sh` say five.

## 2. Edit order

Batches 1 and 5 are independent and could run in either order. Batch 5 is fact-gathering rather than
judgement, so it is the one to delegate.

**Batch 1, programme spine.** 4 files, ~300 lines. `_program/decisions.md`, `_program/program.md`,
`_program/nomenclature.md`; archive `_program/paper-2.md` with a tombstone and a migration ledger.
Settles: two papers not three, the merge and its reason, LSE as the anchor, the roster and its display
names, the diagonal run policy (M-5), and D-I as a retired-phrasing entry. Carry `paper-2.md`'s four
still-live items across rather than losing them (the `family==2` guard status, the dispatch-versus-roster
reconciliation, the duplicate `dispatch_figure_3_G.sh` submissions, the two LSE diagnostic caveats).
While here: `nomenclature.md:22` states in bold that LSE "has no rung and no gloss", which is the
sentence D-D's ladder reverses.

**Batch 2, the measurement layer.** 2 files, ~150 lines. `_program/machinery.md`, `_program/axes.md`.
The threshold and its units, the sign convention, the five regions expressed in grid coordinates, the
grey-is-not-a-region rule, and the LSE n_sims conflict (M-8 puts 1000-sim LSE cells beside 10⁴-sim
band-A cells, which the corpus's own Jensen rule forbids). Partly blocked, see §3.

**Batch 3, the paper hub.** 4 files, ~450 lines. `1_method/decisions.md` first, then
`decisions/D-3_novelty_claim.md`, `D-4_ranking_verdict.md`, `D-0_freeze_and_rerun_scope.md`. The
heaviest sitting and the one with the most downstream leverage. D-A, D-B, M-2, M-3, M-4 and M-8 all
land in `1_method/decisions.md`; D-I lands in D-3 as a new DEAD class plus new CONCEDE entries, which
is that file's native shape; D-J lands in D-4, the only file with the `% src:` recompute discipline.

**Batch 4, plans and rules.** 3 files, ~300 lines. `00_plan.md`, `1_method/docs/manuscript-drafts/sections/README.md`, `01_writing_plan.md`.
Do the bibliography migration in the same sitting: `elife_paper.tex:47` is `\bibliography{biblio}` and
`biblio.bib` has 15 entries, so the first D-I citation breaks the build. What needs copying across is
the pre-1985 lineage (`katz1970membrane`, `katz1972statistical`, `anderson1973voltage`,
`conti1980conductance`, `sigworth1981covariance`) and `lei2020considering`. `celentano2004use` is
already there. Note that `01_writing_plan.md:31-32` and `:155` cite "106 entries, 9 cited, 97 uncited",
which describes neither file.

**Batch 5, figure reality.** 5 files plus 1 archive, ~500 lines. Do an inventory pass over
`figures/paper_both/` and `figures/data/` first, then `provenance.md`, `figures_build_plan.md`,
`CONTINUE_HERE.md`, `figures_system.md`, `1_method/docs/manuscript-drafts/sections/04_results.md`. Archive `04_figures_storyboard.md` and repoint
`figures_system.md:118` and `:149` in the same commit, because they currently route the reader to it as
step 1 of the figure workflow while its own line 3 declares it retired. `1_method/docs/manuscript-drafts/sections/04_results.md` alone is about half
this batch.

**Batch 6a, the conceptual sections.** 5 files, ~450 lines. `1_method/docs/manuscript-drafts/sections/01_introduction.md` (delete the literal
`</content>` and `</invoke>` at `:118-119` first, that is free), `1_method/docs/manuscript-drafts/sections/02_theory.md`, `1_method/docs/manuscript-drafts/sections/03_diagnostics.md`,
`analysis_figure_S1_score_mean.md`, `1_method/docs/manuscript-drafts/sections/06_methods.md`.

**Batch 6b, the claim sections.** 2 files, ~250 lines. `1_method/docs/manuscript-drafts/sections/00_abstract.md`, `1_method/docs/manuscript-drafts/sections/05_discussion.md`. Both carry the
D-J number, so they run last among the prose.

**Batch 7, routing and closure.** 7 small files, ~120 lines. `_program/00_index.md` **goes last**,
because its status column is a claim about the state of everything else. Then `sources.md`,
`research_program.md`, `paper-3.md`, `carve_plan.md`, `1_method/README.md`, and `papers/README.md` cut
to a three-line forwarder.

### The expensive edits: one fact in many places

Sweep these together or the update will be half-done.

- **The dead `figures/paper/` path (M-7).** 14 lines in 6 live files: `sources.md:21`,
  `carve_plan.md:19`, `provenance.md:13,92,203`, `figures_build_plan.md:50`,
  `results.md:71,97,125,168,189,291`, `01_writing_plan.md:159,182`. Plus
  `elife_paper.tex:9`, which is why the paper does not compile. Two are load-bearing:
  `01_writing_plan.md:182` is task W-3's Inputs list, and `carve_plan.md:19` is the carve boundary, so
  the carve currently ships the wrong notebooks.
- **The five-body-figure gate.** 18 locations across 9 files, plus `check.sh:164` and its
  `N_CAP >= 5` test at `:171`. The gate is check 6 of the done-oracle, so a half-done edit leaves the
  oracle green on a paper short one body figure.
- **The 87-nat gap and the 10-16x overconfidence, both routed to "paper 2".** 17 lines across 4 files.
  M-1 removes the destination and M-2 puts NR back in the body, so all of it re-enters the paper.
  **Trap:** `figures_build_plan.md:321,349` record that the caption was already refreshed from 87 nats
  to 10.4 once NR and NMR left the panel. An update that fixes only the routing reinstates a superseded
  number.
- **Retired phrasing that never propagated.** The retirement is logged twice
  (`1_method/decisions.md:100-102`, `nomenclature.md:100-103`) and the phrasing survives at
  `D-4:40` ("Sole survivor; calibrated across the practical regime"), `D-4:42` ("Strawman
  intermediate"), `sources.md:33`, `components/_MAP.md:16`, `research_program.md:30`. D-F turns this
  from a style rule into a substantive error.
- **The paper count.** Owning files are `00_index.md`, `program.md`, `_program/decisions.md`. Fix those
  three and everything else is a sweepable downstream mention.

## 3. Blocked edits, and what unblocks each

**D-J, the ~100 threshold.** Blocks the punchline of `1_method/docs/manuscript-drafts/sections/00_abstract.md`, `1_method/docs/manuscript-drafts/sections/05_discussion.md` D1, the new
threshold section in D-4, `1_method/docs/manuscript-drafts/sections/README.md`'s rule 4, and the verify lists in `1_method/docs/manuscript-drafts/sections/04_results.md`, `machinery.md` §11
and `1_method/docs/manuscript-drafts/sections/06_methods.md`. Three separate things unblock it.

1. Choose the anchor (question 2 of §1).
2. State the criterion. The data do not give a clean 100 under any unstated rule: IR's `k_off`
   distortion at noise 0.1 runs 1.32 / 1.10 / 1.00 / 1.00 at N_ch 10 / 100 / 1000 / 10⁴, so at exactly
   100 channels IR is still about 10% off, and Figure 6's closure boundary sits at **N_ch ~ 70-150** and
   moves with noise.
3. Channels or openings **cannot be settled by measurement**: P_open is fixed at 0.5 in every run and
   never swept, so the two differ by a constant factor of two and the design cannot distinguish them.
   This is a definitional call. If the answer is openings, `D-2_parameter_units.md` needs a
   derived-quantity row.

Until all three land, write the sentence with a placeholder rather than writing 100 and hoping.

**Threshold value and units.** See question 3 of §1. Blocks the region boundaries in `axes.md`, the
verdict cells in D-4, and Figure 6's caption. 1.15 is the incumbent because the built figure uses it.

**Distortion sign convention.** `machinery.md` §4 asserts it three times with two different signs and
admits it was never verified against the producer. Every region boundary and every verdict cell
inherits it. This is a code read, about an hour, not a document edit.

**LSE at body n_sims.** M-8 puts LSE in figures 1 to 4, but the LSE fill is n_sims 1000 while the
band-A cells are 10⁴. Blocks every body panel containing LSE and therefore the headline gap numbers.
Unblocked by re-running LSE at 10⁴, or by restricting those panels to a uniform n_sims and saying so in
the caption.

**The `family == nonlinearsqr` guard in `src/core/likelihood.cpp`.** Verify before writing either
state: part of this has clearly landed, because `figure_3_time_dlik_LSE.csv` exists at 1.1 GB and
`figure_3.Rmd:55` already reads the four-algorithm roster.

**Telegraph, non-Gaussian noise (D-F, D-G point 2).** The word appears nowhere under `papers/`, no
simulator capability exists on any freeze commit, and no cell has been run. This is the only new claim
in the set with literally zero data behind it. Either add the capability and run a few cells, or demote
it to an argument with a citation and say so.

**Physical units on the region-map overlay.** Needs the current unit and the time unit, which are open
questions 1 and 2 of `D-2_parameter_units.md`, and they are the author's to answer. Do not promote the
pA guess to a stated unit. If the map must be drawn first, keep the dimensionless label and mark the
overlay provisional.

## 4. Work that is not a document edit

1. **Repoint `\graphicspath` to `figures/paper_both/` and rebuild.** The paper does not compile.
2. **Bibliography migration** (batch 4, above).
3. **Track `docs/bibliography/recording_configurations/`.** 88 files, untracked. A body-figure input
   currently sits outside version control.
4. **Lift the recording-condition derivation into a Methods table**, one source per row. It is done but
   buried in `figure_6.Rmd` comments, and it is the part a referee checks first.
5. **Reconcile the region vocabulary.** The built map labels regions 0 to 4 bottom-up in noise, from
   "Gaussian closure fails" to "nothing estimable". The voice note numbers them the other way. Two
   numbering schemes for one figure is the same collision §8 of `00_plan.md` just spent a page
   untangling for `D-n`. Settle it in `nomenclature.md`.
6. **Pin Figure 6's lower vertex.** Its own source comment says the vertex where the map closes is
   extrapolated to N_ch 3-9 and proposes N_ch 2 and 5 at noise 0.1-10, ten cells. Separately, the
   `macro_R` battery exists only at N_ch 100 / 1000 / 10000.
7. **Write Figure 6 into `04_results.tex`.** It has five body figures and no region map.
8. **Test the R and Python ports across the grid** (D-H). Cross-checked once via
   `tools/cross_language_check.py`. A usability claim inside a paper about when methods are valid,
   backed by a single cross-check, is the weakest sentence in the submission.
9. **The two engine fixes that gate the freeze** (carry-over): the IDM reconstruction still uses the
   symmetric square root against the corrected `K = H^{-1/2} J_s^{1/2}` (`likelihood.cpp` ~3575), and
   `emit_state_rows_with_experiment` still writes every evolution row twice (~2484 and ~2495). Both
   must land before the tag, because the binary stamps its own git hash into every output file.

## 5. The risk nobody has written down

**D-B's headline payoff is what non-stationary fluctuation analysis has delivered since 1980, and no
document holds that thought and the prior-art thought at the same time.**

The region map's largest region is labelled "i + N_ch + honest CI", and the voice note's message is
that more information is available in almost every recording configuration, "especially conductance".
But the Introduction's own defence against NSFA, at `introduction.md:43` and in `01_introduction.tex`,
**concedes exactly those two quantities**: it grants that NSFA "returns the unitary current, the
channel number and the peak open probability" and argues that it does not return a kinetic scheme.
D-B then moves the headline onto the parameter that defence gave away. D-I, collected the same day,
makes it worse by establishing that reading conductance out of macroscopic fluctuations is old and
standard.

The claim survives, but it has to be restated before the abstract is written rather than after review.
What is actually new is the **combination**: a calibrated joint estimate of rates *and* amplitudes with
honest confidence intervals, from a **single non-stationary record**, where the 1980-1981 covariance
methods needed ensembles of 256 to 504 repeated sweeps and produced no uncertainty statement at all,
and where NSFA's own practitioners report that the unitary current is close to the only parameter it
recovers reliably (Stepanyuk 2014, already cited). "More conductance information" on its own is not
defensible, and it is what D-B currently asks the paper to lead with.

## 6. Leave alone, and archive rather than edit

**Leave completely alone.** `_program/notation_map.md` (current). `1_method/00_master_plan.md` (a
tombstone whose only job is the redirect `00_plan.md:12` depends on). `decisions/D-2_parameter_units.md`,
except to note that its four open questions are now blocking rather than cosmetic; the answers are the
author's and must not be filled in silently.

**Blocks inside edited files that must survive untouched.** `machinery.md` §5 (the K-factor identity)
and §2 (the three tests: **machinery's wording wins over the voice note's**, see below); `axes.md` §1
(the label definition and the K_off decision with its 35,000 CPU-hour costing); `provenance.md` §2, §6
and §9; `1_method/docs/manuscript-drafts/sections/06_methods.md` M1-M11 apart from the four flagged statements; `1_method/docs/manuscript-drafts/sections/README.md`'s two evidence surveys;
D-3's DEAD and CONCEDE lists; D-4 §2.1 and §7; `1_method/docs/manuscript-drafts/sections/05_discussion.md`'s anchor-problem section; `1_method/docs/manuscript-drafts/sections/04_results.md`'s
"claims the data do not support" section; and `00_plan.md`'s "The argument", which already states the
three tests and the force-it-to-fail principle almost verbatim.

**On the three validity tests.** The voice note's phrasing drops "and no autocorrelation" from test 1
and says "the Fisher approximation as a proxy for the Hessian" for test 3. `machinery.md` §2 keeps the
whiteness clause, which is the sharpest IR discriminator on record (score ACF: IR ~0.005 against
NR/NMR ~0.78, `results.md:154`), and anchors test 3 on the model's own Gaussian Fisher, which is what
`03_diagnostics.tex` was written to defend against the circularity objection. **Do not overwrite
`machinery.md` with the voice-note wording.**

**Archive rather than edit.** `_program/paper-2.md`: obsolete under M-1, but it is the ancestor of both
the current body roster and the region map, so archive it with a supersession header and a migration
ledger. `1_method/04_figures_storyboard.md`: already self-declared retired, still routed to as live.
`papers/README.md`: every path in it is dead and it sits outside the index's completeness guarantee, so
cut it to a forwarder rather than rewriting it into a second router.

**Label hazard.** `D-n` now means three things: the manuscript-production briefs in `decisions/`, the
pre-2026-07-21 label for what are now `Q-n`, and the new `D-A..D-J`. Note also that D-3 (novelty) and
D-I (prior art) are the same subject under two registers, as are D-4 (ranking) and D-J (threshold).
Say once, in `1_method/decisions.md`, where each register lives. **Do not renumber the existing briefs.**

## 7. Findings that were checked and are NOT work

Recorded so nobody re-does them.

- **theory.md's "misnomer" note is already fixed.** `00_index.md:95` says `theory.md:85` still calls the
  noise parameter a misnomer. `:85` is a row of the ladder table; the misnomer text is at `:169-172` and
  already records the overturn. Delete the index note and fix the line reference. (A separate live
  question does remain: whether the *name* `noise_in_conductance_tau` is right given that the swept
  quantity is `Current_Noise` with no tau in it. Flag it, do not resolve it here.)
- **The `1_method/docs/manuscript-drafts/sections/03_diagnostics.md` vs `analysis_figure_S1` conflict is a verify, not a fix.** They measure
  different objects: `diagnostics.md:22` states test 2 over the whole record, the score-mean note
  reports max per-step magnitudes and significance counts, and the note's own reading is that NR's bias
  concentrates at the transitions and returns near zero on the plateau, so contributions may cancel in
  the aggregate. Compute the signed record-level score before editing either sentence.
- **`abstract.md:30` needs scoping, not deletion.** "The recursive filter is published and in use" is
  defensible for the originating groups. What D-A kills is uptake *by others*. Do not overshoot into the
  sentence D-I kills: the temporal correlation itself has been used since 1973.
- **Three agents proposed building a region map that is already built.** What is actually missing is
  narrower: the overlay is not in the caption, the region numbering collides, and the lower vertex needs
  ten cells.

## 8. Pre-existing contradictions the audit surfaced, none of them caused by the new decisions

These predate 2026-07-28 and are cheap to close while the files are open.

- **Is NMR in the programme?** `_program/decisions.md:30-33` lists it as live and `:37-38` drops it,
  eight lines apart in the file that owns cross-paper decisions. Both readings have propagated:
  `axes.md:42` (in), `1_method/decisions.md:18-20` (out), `D-4:44` (in, with a full verdict row),
  `D-3:37` (in, as a novelty contribution).
- **Does the distortion fall monotonically along the ladder?** Yes at `D-4:49` and `discussion.md:38`,
  no at `D-3:47` and `00_plan.md:133`. D-4 contradicts itself two lines later at `:50-51`.
- **Is MR's variance direction resolved?** Resolved at `results.md:43-48` and D-4 §2.1 (a category
  error: the predicted observable variance and the reported parameter covariance are different
  objects). Still carried as open at `discussion.md:36` and `D-3:75`, and `title.md:31` states the
  losing side as settled fact inside rule 1, the rule that governs the abstract.
- **Every number in the approved ranking verdict traces to a file that no longer exists.** The
  `% src:` chain in `decisions/` points at nine deleted files (`results_plan.md`, `discussion_plan.md`,
  `abstract_draft.md`, and six more), cited from `D-4:49,67,82,92,103,104,119,132,154,162,163,214`,
  `D-3:76`, `D-2:27,35`. LINT-SRC in `01_writing_plan.md` §1 requires every number to trace to the file
  that computed it.
- **Three topic-index rows route to files that do not exist**: `model_and_sim.md`, `submission.md`,
  `grid.md`. So D-2's units, D-B's grid cells and D-H's code-availability statement each have nowhere
  to land.
- **The hard-boundary ban.** Four documents forbid drawing a hard line on the design plane
  (`title.md:33`, `abstract.md:49`, `introduction.md:22`, `results.md:337-339`, the last being the
  sharpest: "already flagged as method-wrong"). The figure they forbid is built and shipped. Its own
  caption carries the reconciliation: *"The map is a concept map, not a phase diagram. Its boundaries
  are level sets of continuous diagnostics, so a looser criterion moves each of them by up to a decade
  in noise without changing the layout."* Adopt that sentence in all four places rather than deriving
  four versions of it.
