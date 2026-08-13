# Carving Theory and Diagnostics down to what the Results need

Written 2026-08-11. Successor to `08_length_plan.md`, which took the manuscript from 22,727 to
13,765 counted words and stopped with Theory at its Phase-1 floor. This plan is a different cut,
asked for by Luciano in the same terms twice: **the pre-Results sections keep only what a reader
needs to understand the Results; every derivation and every justification goes to Methods or to an
appendix.**

Instrument: `../_program/wordcount.py` (eLife's definition of main text: comments, captions,
equations, Methods, back matter and appendices all excluded). Every number below was measured with
it, not estimated.

**Line ranges: `03_diagnostics.tex` md5 20666e1d, unchanged since this plan was written, so that
table is still good. `02_theory_full.tex` moved the same evening (48d89887 → d6af8bab) and its
ranges are already approximate; re-derive them by heading before executing that table.** The ranges
drift the moment editing starts, in this file as in any other. See the drift log at the bottom: the
manuscript was being rewritten while this plan was being written, and the baseline below is the
re-measured one.

## Decisions already taken (Luciano, 2026-08-11)

1. **Theory and Diagnostics merge into one section** before the Results.
2. **`eq:macror` stays in the body.** It is the one algorithmic anchor a reader gets before the
   Results, and it costs 141 words.
3. **Figure 1 stays in the body.** It is cited from Theory and Methods and from nowhere in the
   Results, and it is the only display item in the paper that measures nothing; it is kept on
   purpose, because it fixes the member vocabulary in pictures at no cost in counted words.
4. Open, under discussion, **not** to be settled by whoever executes this plan: whether the "what is
   measured" block stays inside the merged pre-Results section or opens the Results. See the last
   section of this file.

## The state this plan starts from

Re-measured 2026-08-11 at 19:30, after the Discussion rewrite of that evening.

| section | counted | published median | note |
|---|---|---|---|
| Introduction | 1,310 | 873 | not touched here |
| Theory | 3,356 | — | no counterpart in the survey |
| Diagnostics | 1,260 | — | idem |
| Results | 4,747 | 3,066 | not touched here |
| Discussion | 1,860 | 1,184 | rewritten 2026-08-11, was 3,122 |
| **main text** | **12,533** | 5,209 | target 9,200 |
| Methods (uncounted) | **5,868** | 2,249 | **p95 is 5,046, the plan's own cap is 6,000: 132 words of headroom** |
| Appendix 1 + 2 | 5,348 + 1,294 | — | no stated limit |

## The evidence that sets the criterion

Three facts, measured, not argued.

**The Results cite no equations.** Full reference map over `sections/`: `04_results.tex` cites six
figures and `app:members` once. Zero `\ref{eq:...}`, zero tables. `05_discussion.tex` cites one
equation, `eq:mu-post`. Every numbered equation in Theory and Diagnostics exists for its own
section, for Methods, or for an appendix.

**Two body display items are cited from nowhere outside their own section.** `tab:members` (Table 1,
Theory) and `fig:ladder` (Figure 1). The Table-1 defect was already logged in `08_length_plan.md`
and never repaired.

**Table 1's columns are pointers into the appendix.** It indexes `eq:ypred`, `eq:pred-var`,
`eq:mu-prop`, `eq:sig-prop`, `eq:mu-post`, `eq:sig-post`, all six of which now live in Appendix 1.
A body table on page ~10 whose columns name equations printed twenty pages later.

### The criterion

A block stays in the body if the Results **use** it: it names a member, names a region, names a
plotted quantity, or carries a caveat the Results invoke. Everything else is derivation or
justification and goes down.

### The destination rule

- **Methods**: only what is needed to *run* the thing (flags, values, safeguards, anchors, how a
  quantity was computed). **It receives nothing.** The 150 words it was to receive from Theory
  arrived on 2026-08-11 in the Discussion pass, and at 5,868 it is 132 words from the plan's own cap
  and 822 past the published p95. It cannot be the destination for derivation, and it has to shed.
- **Appendix 1** (`08_appendix_derivation.tex`): every derivation and every justification of a
  modelling choice.
- **A new appendix, "The diagnostics, derived"** (~550 words): the diagnostics algebra, which today
  has no owner. Part of it sits inside the derivation appendix, part of it points at a
  Supplementary Information document that does not exist.
- **Discussion**: the positioning against other people's work that currently sits in Theory.
- **Dies**: anything already stated in Methods or in an appendix. Nothing in the three protected
  categories of `08_length_plan.md` (scope caveats, limitations, concessions) is cut; they are
  compressed or relocated, never dropped.

## Theory: 3,303 → ~1,710

| block | lines | now | keeps | where the rest goes |
|---|---|---|---|---|
| opening: model scope, the two obstructions | 28-65 | 247 | 110 | Ap.1 (the cost of the microscopic filter) |
| the observable is an interval average | 67-84 | 155 | 150 | stays whole; it is the paper's central object |
| uniform window, shared by every member | 85-92 | 117 | 40 | Methods (which already states it at :73-78) |
| the Bessel bound | 94-111 | 259 → **~295** | 60 | Ap.1, **with the deficit table from `decisions/recompute/bessel_bound.py`**, which today lives only in a comment. **Grew on 2026-08-11** (the Colquhoun rise-time citation, the first-order/second-order sentence, the tenfold-cost clause from the Discussion). Protected on the way down: the first-order-for-white, second-order-for-correlated sentence, which is what makes the kernel concession and the noise-model concession one concession; the canonical BESSEL-BOUND note says so and it was the last thing printed of it |
| two closures, prose | 234-274 | 180 | 60 | folds into the box |
| the closure BOX | 275-308 | 298 | 200 | compress; the two *Validity* entries are what predicts the diagonal boundary of Fig. 7 and must survive |
| telegraph-like, no closed form | 310-352 | 169 | 35 | **dies**: Ap.1:24-41 already carries it at more length |
| the three regimes | 354-362 | 122 | 90 | stays; the Results name all three |
| least squares is not a member | 372-393 | 119 | 80 | stays; the copy in Methods goes instead |
| the two axes | 395-412 | 196 | 120 | stays |
| the naming rule and the lattice | 414-442 | 323 | 200 | Ap.2 (why VR's variance form is a choice) |
| the closures from the method side | 477-492 | 132 | 40 | folds into the box |
| IR at the top, prior art, the exact top rung | 494-505 | 154 | 110 | stays |
| the boundary state | 512-524 | 152 | 70 | Ap.1 (static condensation, the K² claim) |
| Eq. macror and the gain | 540-566 | 141 | **141** | **stays whole, decision 2** |
| per-channel moments, the N_ch factors | 583-593 | 76 | 20 | **dies**: Ap.1:42-51 repeats it verbatim |
| sufficiency, and Münch | 595-615 | 189 | 55 | Discussion (the Münch quotations, as relation-to-other-work) |
| the three-step bridge | 629-651 | 169 | 60 | Ap.1 |
| the posterior as the next prior is itself a closure | 653-660 | 105 | 70 | stays; it licenses "model-based information" |
| | | **3,303** | **1,711** | |

## Diagnostics: 1,253 → ~705

| block | lines | now | keeps | where the rest goes |
|---|---|---|---|---|
| the two anchors | 3 | 116 | 45 | Methods is already the owner (:668-677) |
| the standardized residual | 5 | 180 | 90 | stays; it is the horizontal axis of Fig. 5 |
| the score and the bias b | 30 | 100 | 70 | stays; Fig. 4A |
| C and the sign convention | 32-44 | 170 | 110 | the Bartlett-identities sentence → Discussion |
| the Gaussian-Fisher anchor | 46-59 | 76 | 30 | **the equation dies**: identical display at Methods:707-710, and `eq:gaussian-fisher` is cited only from inside this section |
| the least-squares qualification | 61-65 | 170 | 40 | the Results restate it whole at :130-136; single owner = Results |
| the sandwich | 67-71 | 76 | 60 | the Taylor expansion → new appendix (today it points at an SI that does not exist) |
| the split into C_s and R | 73-77 | 147 | 100 | stays; Fig. 4-supp 1 and Fig. 6 are read through it |
| κ and T_eff | 79 | 79 | 60 | stays |
| m, a, D_AI, error-bar units | 81-105 | 139 | 100 | D_AI → new appendix |
| | | **1,253** | **705** | |

## The merged section

One section before the Results, about 2,400 counted words, in this order:

1. the observable, with `eq:obs-avg` and the uniform window as a scope statement;
2. the two Gaussian approximations, as one box, and the three regimes;
3. least squares, the two axes, the naming rule, **Table 1**;
4. the boundary state, the interval update with `eq:macror`, the posterior-as-prior closure,
   **Figure 1**;
5. what is measured: the three checks, C with its sign, the split, κ and T_eff, m and a, and the
   conversion to error-bar units.

Title to settle when it is on the page. It has to name both halves without promising a derivation:
"The likelihood family and how it is tested" was the working one. SETTLED 2026-08-12 as **"The
likelihood family and its diagnostics"**, in the pass that rewrote every body heading against a
verbatim survey of six eLife articles (the rule and the survey are in the comment block at the
first subsection of 04_results.tex). Same two halves, one noun instead of the appositive clause.
Still open in that pass and left to Luciano: all six surveyed articles carry exactly four top-level
sections, so a fifth may be queried; the alternative is demoting this section to the first Results
subsection, which moves decision 1 and not just a heading.

**Table 1 is rebuilt, not moved.** Body version: member, what the states are conditioned on, one
line of plain English. No equation columns. The equation columns move into
`tab:appendix-members`, which today is cited from nowhere and becomes the equation-level
specification. Methods keeps `tab:family` (flags and data keys) unchanged. Three roster tables
become two with disjoint jobs, and the forward-reference defect disappears with them.

## What Methods sheds

Revised 2026-08-11 19:30. Methods is at 5,868 and receives **nothing**: the uniform-window sentence
this plan was going to send it arrived that evening on its own, in the Discussion pass. The three
rows below are no longer optional, and even with all three it lands at ~5,520, still past the
published p95. **If a block has nowhere obvious to go, it goes to an appendix, not to Methods.**

| out of Methods | lines | words | to |
|---|---|---|---|
| the least-squares-is-not-a-member argument, a verbatim-in-substance copy of Theory's | 152-158 | ~90 | a pointer |
| the White 1982 restatement (Diagnostics owns the definition) | 742-748 | ~80 | a pointer |
| the finite-difference Fisher construction and the superseded 433ed13 battery | 717-738 | ~180 | Supplementary File 1 |
| the limitation half of "That license has a boundary" | 137-153 | ~115 of 155 | **back to the Discussion**, see below |

### The block that has to go back

The Discussion rewrite of 2026-08-11 moved three blocks into Methods. Two of them belong there and
stay: the open-probability paragraph (:36-45, 66 words: during the step $k_{on}[A]=k_{off}$, so the
sweep works at $p=0.5$, where the gating variance is largest and the skewness vanishes) and the
perimeter declaration (:125-136, 80 words: simulated by design and not by convenience, likelihood
only, no prior, no evidence). Both are specification of what was run, the Discussion now points at
them rather than restating them, and the transfer rule $4Np(1-p)$ stays on the Discussion side where
a reader with another receptor needs it.

The third does not. "That license has a boundary" (:137-153, 155 words) is the Discussion's old
*third boundary*: the diagnostic compares two independent implementations, so everything downstream
of where they separate is tested and nothing upstream of it is, because they share one
specification. That is a limitation of the study, and eLife prints Methods after the Discussion, so
a limitation that lives only there is one most readers never reach. `08_length_plan.md`'s rule is
that limitations are compressed and never cut; relocating one out of the body is a demotion by
another route. Its last third is not method under any reading: "the general remedy is different in
kind, verifying each object against the theorems it must satisfy ... a piece of work with its own
design, not undertaken here" is future work.

The split: Methods keeps the operative sentence, about 40 words, that the two implementations share
one specification upstream and that this is why the section writes it out symbol by symbol. The
limitation returns to "What this paper does not settle", compressed to about 50. The general remedy
travels with it or dies. The Discussion's closing clause at :131 then stops being an orphan pointer
and becomes the tail of a stated limitation.

## Appendix reorganisation

**Appendix 1 receives** the Bessel bound as a section of its own (with the deficit table against
f_c·Δ, the first-order/second-order argument, the pole-count independence and the group-delay
caveat, all of which exist today only as a comment block at `02_theory_full.tex:161-205`), the
opening's cost argument, the boundary state's static-condensation reading, and the three-step
bridge. Net about +700.

**Appendix 1 gives away** the subsection "How the two parts of the distortion compose". It is
diagnostics algebra, not derivation of the likelihood, and it belongs with the rest of it (~200).

**New appendix, "The diagnostics, derived"** (~550): the Taylor expansion behind the sandwich and
the bias vector b, which today is promised as "(see Supplementary Information)" and exists nowhere;
C = K R Kᵀ with the non-symmetric-square-root caveat and the log-det additivity, arriving from
Appendix 1; D_AI; the finite-sample floor of 1.08 on the anisotropy statistic, which the Results
quote (`04_results.tex:205`) and which has no owner in the manuscript, only a `% src` to
`decisions/recompute/fig2_shape_and_floor.R`; and the constant-variance premise of the
least-squares branch.

Input order: derivation, members, diagnostics, repairs. The first two keep their numbers, so no
cross-reference renumbers.

## Where the paper lands

| | before | after |
|---|---|---|
| main text | 12,533 | ~10,400 |
| of which the pre-Results section | 4,616 | ~2,400 |
| Methods | 5,868 | ~5,400 |
| appendices | 6,642 | ~7,700, plus the new one |

About 1,200 over the 9,200 target, not the 2,400 this plan first estimated: the Discussion rewrite
of 2026-08-11 took 1,262 out on its own. **The remainder does not come from here.** It comes from
the Results (4,747, mostly numeric series a figure already draws, and 122 words measured as verbatim
duplication with its own captions). Say this out loud before cutting the merged section to the bone,
so that nobody discovers it afterwards and starts cutting caveats.

## Execution order

**Relocate first, with no sentence edited, and compile after each move.** Then compress what is
left. The reverse order compresses text that was about to move, which is how a cutting pass runs out
of energy before it reaches the sections that needed it. This is `08_length_plan.md`'s principle and
it held.

Phase A, relocation: the six "dies" and "to Ap.1 / to the new appendix" rows above, plus the three
Methods rows. Phase B, the table rebuild. Phase C, compression of what remains, block by block
against the "keeps" column. Phase D, re-measure and re-read the seams.

### The reference checks this owes

1. `tab:members` must acquire a citation in the Results opening. Pre-existing defect; the rebuild
   exposes it rather than creating it.
2. `tab:appendix-members` must be cited from Appendix 2's own text once it carries the equation
   columns.
3. The equations that change from arabic to A-numbering must have their citing sentences read, not
   just recompiled: a Methods sentence that now points forward into an appendix may need rewording.
4. The dangling "(see Supplementary Information)" at `03_diagnostics.tex:69` must resolve to the new
   appendix. It is the only pointer in the manuscript that resolves to nothing.
5. Run `_program/check.sh` and `_program/wordcount.py --by-subsection` after each phase and log the
   number here.

## The open item: does "what is measured" open the Results?

Not settled. What is at stake, so the discussion does not have to rediscover it:

- At ~705 words the diagnostics block is not a section, which is why the merge happened. The
  question is whether it is the *last part of the pre-Results section* or the *first part of the
  Results*.
- Word budget: inside the merged section those 705 words sit in a section that has no counterpart in
  the eLife survey and is therefore judged on its own terms. Moved into the Results they count
  against a section that is already at 4,753 against a median of 3,066.
- Reading order: the Results' first subsection (retitled 2026-08-12; it is the one-design-cell
  recovery subsection, now "Every member recovers the rates while ...") uses the sandwich, the
  distortion and the two anchors in its second paragraph, so the definitions have to arrive before
  it either way. The difference is whether the reader meets them as theory or as method.
- Precedent in the niche: worth checking against Münch e62714 before deciding, since it is the one
  published article in this problem class that carries a comparable apparatus.

## Drift log

The manuscript was being rewritten while this plan was being written, which is the normal state of
this file and the reason every table above carries a stamp.

**2026-08-11, evening.** The Discussion pass. Discussion 3,122 → 1,860; Methods 5,627 → 5,868;
Theory 3,320 → 3,356; Results 4,753 → 4,747; main text 13,765 → 12,533. Three blocks went from the
Discussion into Methods (above, "The block that has to go back"), the Bessel paragraph in Theory
grew rather than shrank, and `02_theory_full.tex` changed md5, so this file's Theory line ranges
became approximate within the hour.

What this changes in the plan: Methods stops being a destination for anything at all, its three
optional sheds become mandatory and a fourth is added, and the gap left after the carve falls from
~2,400 to ~1,200.

## Log

| date | phase | main text | note |
|---|---|---|---|
| 2026-08-11 | baseline | 13,765 | this plan written |
| 2026-08-11 19:30 | re-baseline | 12,533 | after the Discussion rewrite; Methods at 5,868 |
| 2026-08-12 | A+B+C executed | 11,415 | merged section 3,006; Methods 5,584; appendices 7,907 |

## What the execution of 2026-08-12 actually did

Phases A (relocation), B (the tables) and C (compression) were run in one pass, on Luciano's
instruction. **Nothing was compiled**, by standing instruction; the build is his to run, and the
checks below were static.

Done: `02_framework.tex` created and both originals retired to `archives/` with a dated suffix;
`10_appendix_diagnostics.tex` created as Appendix 3 and repairs renamed to `11_`; the composition
identity moved out of Appendix 1 into it; the Bessel bound printed as an Appendix 1 subsection with
its table and its four properties, which existed only as a comment; the static-condensation reading
moved to Appendix 1; Table 1 rebuilt without equation pointers and `tab:appendix-members` cited from
its own appendix; Methods shed the least-squares duplicate, the White restatement and the numerical
Fisher battery (to Supplementary File 1, new S1.5); the T2 limitation returned to the Discussion and
the Münch quotations moved there from the body; `elife_paper.tex` and `_program/wordcount.py`
repointed; every `(Theory)`/`(Diagnostics)` pointer in the other sections rewritten, and the stale
line numbers in the `[BESSEL-BOUND]` pointer notes stripped and re-aimed at the appendix.

Static checks run: the label/reference map over all sections shows **no undefined reference and no
duplicate label**; `tab:members` is now cited from the Results opening, which was reference check 1.

**Where it landed against the plan.** The merged section is 3,006 counted words against the ~2,400
this plan projected. The gap is real and it is diffuse: after two compression rounds every remaining
block is at or near its "keeps" number except three, the closure box and its lead-in (594 against
435), the naming and axes block (598 against 555) and "What is measured" (864 against 705). Getting
to 2,400 from here means removing content rather than prose, so it stops being this plan's call.
Two candidates, both Luciano's: making the closure box a display item with a caption (about 300
words, and captions are not counted, but that is only honest if the box is genuinely reference
material the reader returns to), and cutting the three-regime vocabulary paragraph, which the
Results use twice.

Main text is 11,415 against the 9,200 target. As predicted, the remainder is in the Results, which
grew to 5,048 in the author's own editing of 2026-08-11/12 and is now the largest section in the
paper.
