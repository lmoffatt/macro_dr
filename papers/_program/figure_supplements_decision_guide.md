# Deciding the figure supplements, one parent at a time

> Updated: 2026-07-31.
>
> **Owns:** the *rule* for deciding how many supplements a figure carries and which vehicle each piece
> of overflow material takes. **Does not own:** the figure set itself (`1_method/decisions.md`, "The
> figure set"), the visual system (`figures_system.md`), eLife's mechanical requirements
> (`elife-author-instructions.md`, `figures_system.md` §1), or build order
> (`1_method/figures_build_plan.md`).
>
> Written to be used **per figure**: §4 is a worksheet to paste into a figure's own thread, §5 gives
> that figure's current candidates so the thread starts from something concrete.

## 1. The rule, and the measurement behind it

**Count questions, not figures.** A figure carries as many supplements as it has distinct objections a
referee can raise against it. That number is not uniform across the set, and forcing it to be uniform
is the one thing the evidence says is off-norm.

Measured 2026-07-31 from `api.elifesciences.org`: 1,079 eLife Version-of-Record research articles,
counting assets labelled `Figure N—figure supplement M`, joined to the public reviews, decision
letters and the reviewer "Recommendations" blocks nested inside `authorResponse`.

| What | Value |
|---|---|
| Median body figures / median supplements | 6 / **7** (p25 3, p75 11, p90 15, max 56) |
| Same, articles with 5-7 body figures (n=631) | 6 / **7** (p25 4, p75 10, p90 14) |
| Articles with zero supplements | 15.7% |
| Pooled over 7,169 parent figures | **46.6% carry none**; 85.9% carry ≤2; 92.2% ≤3; 97.5% ≤5 |
| Structural Biology & Mol. Biophysics | median 9.5 (p90 19) |
| Ion-channel niche / statistical-method niche | median 7 / 7.5 |
| Trend 2012→2026 | flat, no inflation |

The four facts that decide the policy:

1. **Volume is not penalised.** Spearman ρ between supplement count and the eLife Assessment
   strength-of-evidence term is **−0.028** (n=662); with significance +0.019; with max-per-parent
   −0.019. Zero in every direction. eLife states verbatim: *"There is no limit on the number of figure
   supplements for any one primary figure."* The Version-of-Record checklist has no count item.
2. **Volume complaints exist below 1%** and are real: `e07367` (*"There are too many supplementary
   figures, making it hard to read this paper"* — the 2015 policy it cites is dead), `e59371` (editor:
   *"condense … in one (or maximum two) per main figure"*), `e76211` (*"too many supplementary figures
   which aren't all helpful … I can't find a reference for Supplementary Figure 1"*). **All three
   settled at ≤3 per parent.** Note the third is really a citation complaint.
3. **Review is a net adder.** 26.9% of reviews explicitly ask for more analysis; 25-37% of author
   responses announce a newly created figure supplement. Every supplement cut is an answer not
   available in round one.
4. **Figure supplements do get read**, unlike other supplement types. Price 2018 (BMJ Open
   8:e021753): extra figures and tables are the only supplementary category ~60% of reviewers read in
   full; everything else is under 36%.

What is **not** measured anywhere, at eLife or elsewhere: any effect of supplement *count* on
acceptance, assessment, citations or readership. Flanagin 2018 (JAMA 319:410) is presence-vs-absence,
points the helpful way, and its own authors decline the causal reading. Treat any confident claim that
"16 hurts relative to 8" as invention. PDFs in `docs/bibliography/publishing_norms/`.

Two weak signals, both with overlapping Wilson intervals, reported as weak:

- Readability complaints run 14.4% overall and are **U-shaped** in supplement count: 18.3% at zero,
  10.9% at 1-4, 12.0% at 5-9, 16.5% at 10-14, 17.5% at 15+. The minimum is 1-9.
- The same complaint is **flat** in max-per-parent (14.2 / 13.9 / 15.2 / 14.8 across 0-1, 2-3, 4-5,
  6+). So there is no empirical case for a per-parent cap. The cap below is a convention argument.

One signal that is neither weak nor about volume: **uncited or out-of-order supplements** draw
complaints at 4.3%, flat in count. That is hygiene, and it is what `e76211` was actually about.

## 2. The procedure, per parent

For each candidate attached to this figure, answer in order.

**Q1. Write the one sentence in the body that cites it.** Not a topic, the actual sentence.
- Cannot write it → the candidate is not a supplement. Repo only, or cut.
- The sentence is a **conclusion sentence, or a claim that appears in the abstract** → it belongs in
  the **body**, not in a supplement. This is the hard rule; see §6.
- Otherwise continue.

**Q2. Which job does it do?**

| Job | Example in this paper | Vehicle |
|---|---|---|
| A. Kills a named alternative explanation | conditioning (κ) is not why IR wins | figure supplement |
| B. Shows the body claim holds for the parameters/conditions not shown | all five parameters, not two | figure supplement |
| C. Carries the mechanism behind a body claim | sample × correlation decomposition | figure supplement |
| D. A standard check the field expects | multivariate calibration (Mahalanobis Q-Q) | figure supplement |
| E. Mathematical development, definitions, derivations | the affine-invariant metric; the factorisation | **appendix figure** |
| F. Build variant, exploratory, superseded | seed sweeps, old validity maps | repo only, cited in Data Availability |

**Q3. Is it the *same question* as another candidate on this parent?** If two candidates answer one
objection with two scalars, merge them into one multi-row supplement. This is the most common case in
the Figure 4 stack.

**Q4. Is the parent right?** Parentage follows the *claim*, not the visual grammar. A candidate whose
sentence lives in the Figure 5 paragraph is a Figure 5 supplement even if it is drawn on Figure 4's
plane.

**Q5. Count.** Target 0-3 per parent. Above 3, re-run Q3 and Q4 before accepting it: 92.2% of eLife
parents sit at ≤3, and all three documented complaint cases settled at ≤3. Nothing in the data
penalises 4+, so 4 is allowed when Q3 and Q4 leave four genuinely distinct objections standing. A
parent at **zero is normal** (46.6% of eLife parents) and is not a hole.

## 3. Whole-paper sanity check

Run once, after the six threads are done, not during.

- Total supplements: 7 is the median for this shape, 14 is p90, ~16 is p92. Anywhere in there is
  unremarkable. Do not pad up to a total and do not cut down to one.
- Distribution should be **clustered, not uniform**: most supplements on the 2-4 figures carrying the
  contested claims, others at 0-1. Uniform 2-3 across six parents is the off-norm shape.
- Two architectures coexist in this niche and this paper is both. Theory-led papers use appendices and
  zero supplements (Münch `e62714`: 12 body figures, **0 supplements**, 9 appendices, 13 appendix
  figures). Benchmark-led tool papers use many supplements and no appendices (Tapqir 9/10, miniML
  9/11, DISC 5/11, EPI 5/18, DeepFRET 4/19). Expect the final shape to be hybrid: supplements on the
  calibration figures, an appendix for the mathematical development.

## 4. Per-figure worksheet

Paste into the figure's thread and fill it there.

```
FIGURE N — supplement decision            (guide: papers/_program/figure_supplements_decision_guide.md)

Body claim of this figure (one sentence):
Claims of this figure that appear in the abstract:

Objections a referee can raise against THIS figure, one line each:
  O1.
  O2.
  O3.

Candidates on hand (from §5), one row each:
  file | citing sentence (Q1) | job A-F (Q2) | duplicates? (Q3) | right parent? (Q4) | verdict

Verdict per candidate: BODY | SUPP | APPENDIX | REPO | CUT

Count: ___ supplements, ___ appendix figures.  If >3 supplements, state which Q3/Q4 pass was run.
Every SUPP has a citing sentence, and they are cited in ascending order in the text: yes/no
```

## 5. Candidates on hand, by parent

Inventory of what is already built, 2026-07-31. **~31 candidates**, not the 16 the planning documents
have been working with; the extra came from `paper_1` and from the `figure_5_IR_*` family. Titles are
the `.Rmd` front matter. Paths relative to `projects/eLife_2025/figures/`.

**Figure 1 — the filter step over the window × recursion lattice** (roster `NR, INR, R, IR` since
2026-08-12; the body figure has no supplement)
- `archive/figure_1_all.Rmd` — one filter step across the whole roster (NR, MNR, R, MR, VR, IR)
- `archive/figure_1_superseded_20260812.Rmd` — the six-column version, `LSE ILSE NR INR R IR`
- `paper_1/Figure_S1*.pdf` + `paper_both/Figure_1_seed_{4,16,47,48}.pdf` — the same step on other seeds

**Figure 2 — recovery clouds**
- `paper_1/figure_S2.Rmd` — cloud colour-coded by per-recording maximised logL
- `paper_1/figure_S3.Rmd` — full correlation corners, NR / MNR / R / MR

**Figure 3 — the calibration cascade in time** (2 supplements, settled 2026-07-31 after the parent
itself moved). **The parent absorbed its own first supplement.** Figure 3 had walked four members
and supplement 1 carried all seven; once `macro_INR` ran and every column turned out to carry a
claim the paper already makes in prose, holding two figures apart for the sake of two columns
stopped making sense, and the seven-column figure was promoted whole. The four-column predecessor is
archived at `figures/archive/figure_3_4col_superseded_20260731/`, still runnable as the only
notebook that builds the figure straight from the dumps.
- `supplement_1` — **per-step information against score variance, all four identified parameters**:
  all seven members across the columns, the parameters down the page as bands of two rows. Revived
  from the archived `figure_3_supplement_2.Rmd`. Cited at `03_diagnostics.tex:47` and twice more.
- `supplement_2` — **the checks the body figure has no room for**: the mean standardized residual,
  and the score bias and autocorrelation for the three parameters the body does not colour (k_on,
  the unitary current, the noise level). Job D, completeness. Cited in Results after supplement 1.

**Three rules this parent produced, and they generalise.**

1. **PARAMETER IS A ROW, ALGORITHM IS A COLUMN**, and a panel carries at most two or three series
   when one of them is a different KIND of object. A first attempt at supplement 1 drew four
   parameters as four colours inside the parent's panels and was rejected on sight: four curves and
   four ribbons in a 1.3 in panel is unreadable, and the two amplitude parameters overlie during the
   pulse so one is invisible. Rows are the cheap axis, since eLife sets no height limit while panel
   width is fixed.
2. **THE CONTAINMENT TEST.** Before building a supplement, list its rows against the parent's and
   count how many are new *and* informative. A residual triptych (mean, variance, autocorrelation)
   was built and scrapped on this test: two of three rows were already in the parent, and of the two
   new things one separates nobody. One informative new row is not a supplement, it is a row that
   belongs in the parent — which is where the residual autocorrelation went, into row G beside the
   score, where the gap between them can be read off one panel.
   The same test run the other way is what promoted supplement 1 into the body: when the parent
   is contained in its own supplement, the supplement is the figure.
3. **A COLUMN COSTS PROSE.** A body column obliges the text to say something about it. That was the
   argument against promotion here and it was wrong, but only because MR and VR turned out to carry
   the non-obvious half of the claim (one interval end is not enough; MR is *worse* than R rather
   than between R and IR). Check what a column is worth before deciding it is decoration.

**Two rendering defects found while building this set**, both fixed in every file:
the family's agonist band used `ymin = -1e10, ymax = 1e10`, which rasterises to *nothing* on a
linear axis whose window is of order one, so every LINEAR row rendered unshaded while the log rows
were fine. `-Inf`/`Inf` is the idiom that works, and is NaN on a log axis, hence a tiny positive
lower bound there. And in a grid where the LEAD column carries the band title, row letter and y
axis, replacing that one cell with a `theme_void()` placeholder silently strips all four from the
whole row; an excluded cell has to be an ordinary panel with no data layers.

**Open on this parent:** two sentences in Results and Discussion name the unitary current among the
parameters whose information dies at washout. Measured, only k_on and N_ch reach the numerical
floor; the unitary current falls by a factor of thirty-six and levels at k_off's order. Flagged in
both `.tex` files as `TODO(Luciano)` and not edited: it changes a stated result.

**Figure 4 — the design plane** (11 built; `supplement_3` has no PDF yet)
- `1`, `2` — bias and distortion, all five parameters, R vs IR
- `3`, `4` — the distortion split into sample and correlation (map, then lines vs N_ch)
- `5` — Mahalanobis Q-Q, sandwich vs the empirical multivariate MLE distribution
- `6` — κ of the Gaussian Fisher over the same plane (the conditioning alternative)
- `7` — distortion-corrected standard error over the plane (what is actually achieved)
- `8`, `9` — residual and score autocorrelation over the plane
- `10`, `11` — whole-matrix affine distortion, and its sample/correlation decomposition

**Figure 5 — the information budget** (5 built)
- `Figure_5_budget_supplement_1.pdf`
- `figure_5_IR_acceptable_distortion`, `figure_5_IR_affine_collapse`,
  `figure_5_IR_affine_validity`, `figure_5_IR_affine_validity_byNch` — the IR-only validity family;
  check against `1_method/decisions.md` before reusing, the IR-only map was merged into Figure 4

**Figure 6 — the usage map**: **zero, settled 2026-07-31**, and the thread's own question answered
itself. The referee's obvious attack was the threshold (1.15, factor 2), and the choice was between
a supplement and drawing the sensitivity on the body figure. The body figure now draws it: the
ribbons are the boundary swept over acquisition intervals, and the caption's closing sentence says
the map is a concept map whose boundaries move by up to a decade in noise under a looser criterion.
The two variants that had been standing in for that answer, `regions` and `frontiers`, were
superseded when the canonical map became `Figure_6` plain; `figure_6.Rmd` writes them to
`figures/archive/figure_6_superseded_20260731/` and their `\figsupp` declarations were removed,
which is what had been failing the compile. Zero is the eLife norm (46.6% of parents) and is not a
hole. What no display item carries any more: the five named regions and the three recording
footprints. Bringing those back is a new figure, not a restored reference.

**Unparented, decide the parent before the verdict**
- `paper_1/figure_7_{sim_by_Nch, sim_by_interval, sim_by_noise, contour, contour_by_Nch,
  contour_by_interval, lines}.Rmd` — r̄²_std calibration and residual distortion over the plane
- `paper_1/figure_S4_S5.Rmd` — bias and autocorrelation diagnostics, data and inference level

**Not built, and each answers an objection nothing on hand answers.** Ranked by how likely a referee
is to raise it: (a) does the result depend on the two-state scheme; (b) sensitivity of Figure 6's
boundaries to the thresholds; (c) computational cost per algorithm. Building these outranks trimming
the list above.

## 6. Hard rules

1. **A claim in the abstract, or in a conclusion sentence, is a body figure.** Münch's Appendix 9
   states that scheme *topology*, not state count, sets the scale of over-confidence — a conclusion
   that appears nowhere in his body ("topolog" occurs 0 times there). That is the failure to avoid,
   and it costs more than any count decision.
2. **Cite every supplement in the main text, in ascending order.** 4.3% of reviews complain about
   uncited or out-of-order supplementary figures, and that rate does not fall with fewer supplements.
3. **One page, one panel-set, no composite, no tables.** eLife's guide: *"We also encourage authors to
   avoid composite figure supplements wherever possible."* Composites are the format rule staff check.
4. **Count and main-text real estate are independent.** `\figsupp[none]{...}{...}` suppresses the
   auto-emitted mention line, so a supplement can exist without spending a sentence in the body.
   Deciding to carry one costs no prose.
5. **Fallback venues change the vehicle, not the content.** Biophysical Journal has no cap but wants
   one flat-numbered PDF (S1…Sn) with supporting references merged: renumber-and-merge, cheap. JGP
   discourages supplementary material as policy: that transfer is a rebuild, and is a reason to keep
   the mathematical development in an appendix rather than in supplements.

## 7. Open

- **To confirm, not yet acted on:** the requirement to name each supplement in its parent's legend
  appears to have been **deleted from eLife's live guide at the 2023 reviewed-preprint launch**,
  surviving only in the stale `elife.cls` v1.11. `figures_system.md` and the
  `reference_elife_figure_guidelines` memory still carry it as a SHOULD. Verify against
  `docs/bibliography/publishing_norms/eLife_2026_AuthorGuide_*` before editing either.
- The `figure_5_IR_*` family predates the Figure 4/5 fusion. Some of it is superseded, not
  supplementary. Resolve in the Figure 5 thread.
- `1_method/decisions.md` declares 5 supplements. This file lists ~31 candidates. The declaration is
  the owner and must be updated once the six threads return their verdicts.

**Sources.** Census and reviewer-text analysis: own measurement, `api.elifesciences.org`, 2026-07-31,
n=1,079 articles; scripts were session-temporary, the numbers in §1 are the record. Journal policy and
meta-research: `docs/bibliography/publishing_norms/` (37 files, added 2026-07-31, **not yet in
`biblio.bib`**).
