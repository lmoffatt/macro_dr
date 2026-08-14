# Bringing paper 1 to length

Written 2026-08-09. Target and evidence: `../_program/elife_main_text_length_survey.md` (413
published eLife articles, the two subject areas this programme samples). Instrument:
`../_program/wordcount.py`, which counts the main text eLife's way (comments, captions, equations,
Methods, back matter and appendices all excluded).

## The target

**About 9,200 counted words of main text.** Revised upward from 8,000 on 2026-08-09, after Phase 1
made the retained Theory measurable rather than estimated. Luciano's call, with the trimming beyond
it left to him.

The original target was 8,000, the 90th percentile of published articles, on a budget that put
Theory at 1,700. That 1,700 was wrong: it counted the two-axis block and the closure box at their
compressed targets and forgot everything else in the section. Measured against what Phase 1 actually
left, the floor is about 2,650: apertura and notation 250, the observable 255, the two closures 490,
least squares and the two axes 850 which are the vocabulary the Results are written in and cannot
move, the boundary state 143, the update with the sufficiency concession and the posterior-as-prior
closure about 400, the top of the ladder 71, and 200 for the new paragraph on how the algorithms
were derived. The budget below therefore sums to about 9,150.

9,200 is the 95th percentile. Two things make it defensible where 11,000 would not be. Münch
`e62714`, the article eLife already published in this problem class, runs 11,515 and sits at the
99.3rd; and this paper carries a Theory section, which almost no article in the survey does, so the
comparison that matters is against the theory-led shape rather than against the median.

The alternative, taken and rejected: Results to 3,000 and the Discussion to 1,300, both near the
published median, which lands 8,850. It was rejected because the Discussion is where the paper's
concessions live, and the rule below says the target moves before a concession does.

Baseline, measured today: **22,727**. That is 2.83× the p90 and 1.97× Münch e62714, the longest
comparable article in the niche.

| section | baseline | after phase 1 | target | still to go |
|---|---|---|---|---|
| Introduction | 2,026 | 1,853 | its floor | done |
| Theory | 8,883 | 3,533 | its floor | done |
| Diagnostics | 2,392 | 1,336 | its floor | done |
| Results (prose) | 5,371 | 4,197 | its floor | done |
| Discussion | 4,055 | 3,129 | its floor | done |
| **total** | **22,727** | **14,048** | — | **the pass is done; 8,679 words came out** |

Published medians for scale: Introduction 873, Results 3,066, Discussion 1,184. Every target above
sits between the median and the p75 of what eLife prints, except Theory and Diagnostics, which have
no counterpart in the survey because almost no article carries them.

Not counted, but not unbounded either. Methods is at 4,986 against a published median of 2,249 and a
p95 of 5,046, so it can absorb roughly a thousand words and no more before it is remarkable in its
own right. Appendix prose has no stated limit; 7,500 words would put the paper in the top 2% of the
413, with precedent to 13,791 and with Münch at 5,889 across nine appendices.

## The principle that sets the order

**Relocate before compressing.** Compressing a paragraph that is about to move is wasted work, and
it is the way a cutting pass runs out of energy before it reaches the sections that needed it. Phase
1 moves whole blocks with no prose edited at all and takes the count from 22,727 to about 9,900 by
itself. Everything after that is small, and by then the budget is visible.

**Nothing is deleted to save words.** Three categories are compressed and never cut: scope caveats,
limitations, and concessions to other methods. If a target cannot be met without touching one, the
target moves, not the caveat.

## Phase 0. The instrument (done)

`../_program/wordcount.py`, with `--by-subsection`. Run it after every phase and record the number
in the log at the bottom of this file. Without it the plan is a feeling.

## Phase 1. Relocation, no rewriting (about one day)

Create `sections/09_appendix_derivation.tex` as a second `appendixbox`, input after
`08_appendix_members`. Move whole blocks, cut and paste, no sentence edited. Compile after each move
and confirm the reference count is unchanged.

Out of Theory, 3,896 words. Every row measured with the instrument over the line range given, on
2026-08-09; an earlier version of this table was built from a triage that counted captions and
different spans, and it overstated two rows by about 40%.

| block | lines in 02_theory_full.tex | words | to |
|---|---|---|---|
| The boundary construction: three steps, Eqs. mu-bnd through sig-post | 824–1016 | 1,096 | Appendix 2 |
| The linear-filtering frame: Eqs. kf, kf-mean, kf-var, kf-cross, mr-kf | 1328–1439 | 1,003 | Appendix 2 |
| MR against IR by total variance and covariance: Eqs. ltv, ltc, cdiff, and the D²-against-D asymmetry | 1087–1231 | 974 | Appendix 2 |
| The start-conditioned members stated forward | 1044–1085 | 374 | Appendix 2 |
| The four contractions and the per-state variance generalization | 751–794 | 275 | Appendix 2 |
| What was verified, and to what tolerance | 1519–1537 | 174 | Appendix 2 |

Out of Theory into Methods, 942 words (moment machinery 498–687: 661; what is enough to reproduce
this 1478–1516: 213; the symbol table 1441–1476: 68), all of it formula-to-evaluate rather than
derivation, and all of it already pointed at from Methods today:

- the tilted-semigroup derivatives and the two conductance objects, Eqs. tilted-derivatives and
  gammabar-vbar;
- the spectral evaluation with the divided differences and the coincidence limits, which belongs
  beside the pseudo-count regularizer already in Methods;
- the endpoint limit Eq. endpoint-limit;
- "What is enough to reproduce this": Eq. loglik and the initial condition Eq. init;
- the symbol-to-implementation table.

If Methods runs past 6,000 after this, the symbol table and the spectral recipe go to Appendix 2
instead. Methods must land under 6,000.

Out of Diagnostics: the composition identity C = K R Kᵀ with the non-symmetric-square-root caveat
and the log-det additivity, and the sandwich Taylor expansion, to Appendix 2 (about 300 words). The
near-singular anchor rule with its 1e-10 tolerance and the grey-cell convention, to Methods (149).

Out of Results: the `fig:limit` subsection (was "Where the calibrated member departs"; retitled
2026-08-12 to "In the few-channel corner ...") with Figure 5 and its caption,
to Appendix 2. 306 words of prose, a 231-word caption and a full-width float. `fig:limit` is cited
from nowhere outside Results, which was checked, so nothing breaks.

Out of the Introduction: the ARMA and generalized-least-squares half of the two-alternatives
paragraph (240) to Appendix 2 or to the Discussion's relation-to-other-work paragraph.

Out of Methods to a deposited supplementary file: the grid coverage manifest and the five reduction
stages (491 down to 80), the SLURM and environment paragraph (70), and the optimiser damping
schedule (40).

Expected after Phase 1: main text about **16,900**, with Theory at about 4,050. Relocation moves
5,800 of the 14,500 words that have to go, which is 40% of the job and not the bulk of it. The
ordering principle stands, the claim that everything after Phase 1 is small does not: Phase 3 carries
about 8,700 words of compression, of which about 2,350 are inside what is left of Theory.

### The mechanics Phase 1 owes, three of them verified against a probe build

**A `figure` float cannot go inside `appendixbox`.** It is a hard LaTeX error, "Not in outer par
mode", the caption is dropped and the label goes undefined, because elife.cls wraps appendixbox
content in `mdframed`. The manuscript already knows this: 08_appendix_members.tex:148 carries the
note and the workaround. Relocate Figure 5 as `center` plus `\captionof{figure}{...}` with its
label, matching what `tab:appendix-members` already does. And drop `fullwidth` for it: inside an
appendixbox the fullwidth `adjustwidth` escapes the panel and overprints the left rule and the
margin, with no overfull warning, so the build looks clean and the page is wrong. Re-render Figure 5
at the appendix's inner text width rather than scaling it.

**The appendix order should be reversed.** 08_appendix_members.tex cites twelve Theory labels, eight
of which are cited from nowhere else, and its opening sentence says "The Theory section gives the
cycle for the boundary-conditioned member", which stops being true. Input the derivation appendix
*before* the member spec so the derivation is Appendix 1, and rewrite the two sentences at
08_appendix_members.tex:19-21 and :47 to point at it.

**Five hard-coded "Appendix~1" strings will silently become wrong**, at 02_theory_full.tex:1123 and
:1319, 04_results.tex:147, and 06_methods.tex:238 and :317. None is a `\ref`, and no appendixbox
carries a label. Add `\label{app:derivation}` and `\label{app:members}` after each
`\begin{appendixbox}` and convert the five literals before touching the order.

### The three reference checks Phase 1 owes

1. Theory equations that move into the appendix acquire A-numbers. Four of them are cited from
   outside: mu-post (from Methods and the Discussion), sig-post, tilde and vector-tilde (from
   Methods). The other two of the six on the earlier version of this list, gammabar-vbar and
   endpoint-limit, go to Methods and stay arabic, so they need no attention. LaTeX renumbers the
   four, but read those sentences afterwards: a Methods sentence that now points forward into an
   appendix may need its wording adjusted.
2. `tab:members` (Table 1) is currently cited from nowhere outside Theory, although the Results use
   its eight member codes throughout. Add the citation in the Results opening. This is a defect the
   move exposes rather than creates.
3. `tab:appendix-members` is cited from nowhere at all. Cite it in Appendix 1 or drop it.

## Phase 2. The two summaries (about one day)

The only genuinely new writing in the plan.

**Theory, about 1,700 words, which is its floor and is derived rather than chosen.** What the
Results cannot be read without: the interval-average observable with Eq. obs-avg and the uniform
window as a scope statement (250); the two Validity entries of the closure box, which predict before
any measurement that the failing boundary runs diagonally across the channel-count-by-noise plane
(260); the three regimes, least squares, the two axes and the naming rule (777); the update stated
in words with Eq. macror (180); MR against IR in words, with Eq. mr-vr and Eq. cdiff and no
derivation (150); the top of the ladder (76). Table 1 and Figure 1 stay and their captions are not
counted.

Plus, and this is the paragraph that replaces what left, **about 200 words on how the algorithms
were derived**, not on what they deliver: run the Gaussian update on the boundary pairs instead of
the K states, charge to the noise term whatever the conductance still varies by once the pair is
fixed, marginalize the start index the next window has no use for; and the pair-conditioned moments
are two derivatives of the tilted semigroup, in closed form, with no object larger than K×K formed
anywhere. Close with the pointer to Appendix 2.

**Diagnostics, about 700 words**, folded to a short block: the three checks named, C with its sign
convention, m and a, κ and T_eff, and the conversion to error-bar units (1.15 is 7%, 2 is 41%, 4
turns a nominal 95% interval into about 68%). Whether it stays a top-level section at 700 words or
becomes the opening of the Results is a presentation call to take after seeing it on the page.

## Phase 3. Compression (three to four days)

Block by block against the triage of 2026-08-09, which named the target and the sentences for every
block. Order by return.

- **Discussion, 4,055 to 1,400.** The largest remaining gap and the biggest relative outlier in the
  paper: 3.4× the published median and within 800 words of the longest Discussion in the survey.
  The map paragraph goes from 295 to about 45 because every claim in it is already in the Results or
  verbatim in the Figure 6 caption. The bootstrap-alternatives paragraph (216) goes to Appendix 2.
  "What this paper does not settle" goes from 869 to about 500 by compression only.
- **Results prose, 5,371 to 3,200.** Mostly the numeric series that a figure already draws, and the
  122 words measured as verbatim duplication with this file's own captions.
- **Introduction, 2,026 to 1,200.** Eleven paragraphs is roughly twice what an eLife introduction
  runs. After the ARMA block leaves in Phase 1, the design-plane specification moves to Methods and
  the error-bar conversion table moves to Diagnostics, where it belongs anyway.

## Phase 4. Captions and floats (half a day)

Captions do not count toward the target and do cost pages. Figure 2 from 379 to 250, Figure 4 from
400 to 260, both against a published median of 150 to 260. Then re-measure the float geometry: the
7 pt type floor is measured from the non-panel space, and this file already records a case where a
40-word caption growth pushed a legend onto the footer.

## Phase 5. The defects, fixed while the text is open

Not length work, but every one of these lives in text that Phases 1 to 3 are rewriting, so fixing
them later means opening the same paragraphs twice.

1. Methods states that the Gaussian and finite-difference Fisher constructions agree; the repo's own
   measurement says they do not for the least-squares arm. One of the two has to change. LUCIANO'S
   CALL, and the ground was checked on 2026-08-10 so that the call is informed rather than open.
   What is NOT wrong, verified against the producer: the anchor claim. `src/core/likelihood.cpp:3462`
   builds `G_b = mean<Sum<Gaussian_Fisher_Information>>`, the ensemble mean of the per-recording sum,
   and `figure_4_common.R:103` reads `Probit_statistics_Likelihood_Gaussian_Information_Distortion`
   off the `battery_pool_G` and `battery_sim_G` families, so every number in the figures really is on
   the Gaussian anchor. Note in passing that the codebase carries a second object,
   `Likelihood_Information_Distortion`, anchored on the numerical Fisher instead ("F_b as the H
   reference, numerical truth, not the cheap Gaussian-formula approximation", same file), and the
   manuscript uses none of it. What remains wrong is only the sentence in Methods that offers battery
   `433ed13` as the demonstration that the two constructions target the same information, against a
   measurement that says they do not for the least-squares arm. Two honest exits: delete the claim
   and say the battery is retained without being used, or report the disagreement, which makes it a
   result and touches the least-squares distortion numbers.
2. "Per interval every member reports its information correctly" is contradicted by its own bound:
   ILSE sits at +0.184 in log₁₀, a factor of 1.53, and MR at −0.229. Scope the sentence to the
   gating-aware members.
3. H is defined with an expectation in Diagnostics and without one in Methods. For the recursive
   members those are different objects. The Phase 1 deduplication forces a single statement; make it
   the right one.
4. The 93 / 70 / 31 percentages are computed over different cell sets (560, 238 and 224 cells).
   Recount on the common sub-grid. No runs needed.
5. The least-squares distortion is called non-commensurable with a gating-aware one in Diagnostics
   and then plotted as the same level set on Figure 6. Distinguish the two readings in one sentence:
   commensurable as a calibration statement, not as a mechanism.
6. The enumeration collision at the top of Diagnostics, where "the first / the second" runs twice in
   consecutive paragraphs for different things.
7. The stale notes: Figures 2 and 3 already carry their eight columns, and the bibliography compiles
   clean with all 26 keys, so the `NEEDS RE-RENDER` and `PENDING-BIB` markers are spent.

## Phase 6. Verification

Recount with the instrument. Confirm no undefined references and no orphan labels. Confirm that no
claim in the abstract rests on relocated material, which is the programme's own hard rule. Confirm
Methods is under 6,000 and the appendix total is stated in the log. Rebuild and check the page
count of the article proper, which was 53 pages of the 61.

The abstract is not touched until this phase. It is at 211 words and settled, and it should be read
last, against the paper that then exists.

## Log

| date | phase | main text | note |
|---|---|---|---|
| 2026-08-09 | baseline | 22,727 | commit 5e9f738. Theory 8,883, Discussion 4,055, Results 5,371, Diagnostics 2,392, Introduction 2,026 |
| 2026-08-09 | 1, plumbing + move 1 | 21,724 | single master (`elife_paper.tex` at `02_theory_full`), `02_theory.tex` and `elife_paper_full.tex` retired to archives; Appendix 1 created as the derivation and the member spec renumbered to Appendix 2, both labelled, five hard-coded "Appendix~1" strings converted; the linear-filtering frame moved (−1,003). Not yet built. |
| 2026-08-09 | 1, complete | 17,638 | Diagnostics: the composition identity to Appendix 1, the near-singular anchor convention to Methods (2,392 → 2,141). Methods: the cell manifest, the five reduction stages, the cluster settings and the optimiser schedule to Supplementary File 1, four pointers left behind (5,412 → 5,236). Two items of the plan did NOT execute. The Introduction's ARMA block is reassigned to Phase 3: the Discussion is itself 3× over its own budget so it is the wrong host, and a derivation appendix is the wrong genre for a comparison with the statistics literature. And Figure 5, see the entry below. Appendices hold 6,064 words. |
| 2026-08-10 | 3, Introduction | 16,953 | 2,026 → 1,853. Five restorations from the loss ledger; the fixed-spectral-density clause moved to Methods rather than cut. Over its 1,250 target: eleven paragraphs against a published median of 873, six of them flagged untouchable. |
| 2026-08-10 | 3, Diagnostics | 16,148 | 2,141 → 1,336. Eleven restorations, five from the ledger, including the name of the two Bartlett identities, the rule-6 concession that in the classical setting the identities can only be tested, and the limitation that an approximation can pass the residual and score checks and still violate the information one. The double "the first / the second" enumeration is gone. |
| 2026-08-10 | 3, Results | 14,970 | 5,371 → 4,193 of prose, captions untouched. Three ledger restorations and four accuracy repairs, among them the two orange numbers that had been dropped from the sentence that then referred to them. |
| 2026-08-10 | 3, Discussion | 14,042 | 4,055 → 3,127. Far over its 1,500 target and the arithmetic says it must be: about 2,270 words of the section are caveat, limitation and concession before a sentence of argument, and the rule is that the target moves first. What went is duplication with the Results and the Figure 6 caption. Six restorations, and one promotion out of a comment line vetoed. |
| 2026-08-10 | 4 and 5, partial | 14,048 | Captions measured against the survey rather than trimmed: the six body captions run 230 to 412 words against a published band of 150 to 260 and an observed maximum of 447, so all are inside the range and the trim is not worth the risk to legends that must stand alone. The float geometry still needs re-measuring on the next render, which is the one live item. Defects fixed: the per-interval identity claim scoped to the gating-aware members, the H definition in Methods given its expectation, and the four stale markers closed (two PENDING-BIB, one PENDING-FIG2, one NEEDS RE-RENDER), each verified before closing. |
| 2026-08-09 | 2, the Theory core | 17,126 | Theory 4,045 → 3,533. Opening and notation box 606 → 369, the observable 365 → 285, the two closures 1,343 → 908, and the new 218-word paragraph on how the algorithms were derived, inserted at the end of the interval update. Every block came in over its target and the three judges were right that the floors are real: the closures block at 490 had dropped the lumpability limitation, the non-factorizing posterior with its citation, the sentence that averaging over the window is itself what brings the law closer to Gaussian (which grep finds nowhere else in the manuscript), the instantaneous member's misspecification statement, and $\Gamma=\mathrm{diag}(\gamma)$, which is defined once and used eight times in Appendix 1. All five restored; the block sits at 908. The tilted semigroup and $\Gamma$ now live in the new paragraph, which is what the appendix's two "above" pointers resolve to. |
| 2026-08-13 | outside read, pass 1 | 12,668 | A fourth outside read (GPT, on the built PDF of that morning, 69 pages) returned a 54-item cutting list. **Its arithmetic does not close and that is the first thing to know about it:** it promises 2,250-3,350 words, but 500-800 of those are Methods, ~200 are captions and ~50 are appendices, none of which this instrument counts, so the countable part is 1,550-2,300 against a gap of 3,832. It is also calibrated on "10-15% shorter" and on the judgement that the manuscript is not bloated, which is a different target from this file's 9,200. WHAT WAS TAKEN, 364 counted words and 234 caption words: the Introduction's "They are also harder to watch" paragraph deleted and three compressed; the Kullback-Leibler reading of row (A) out of the Figure 3 lead, which is the same call Luciano made on 2026-08-12 on the block below it; the "Conditioning localises the error" paragraph deleted as the per-parameter twin of the 93/70/31 paragraph, with its corner clause folded in; six numeric series moved from the line to their `% src` under this file's own test (replace the number by "measurably" and see whether the sentence dies); four captions trimmed, Figure 1 331 -> 275, Figure 4 408 -> 323, Figure 5 441 -> 363, Figure 7 283 -> 268, against a published band of 150-260. Backmatter: the availability sentence now says "regenerable ... up to the simulation seed, which was not recorded", which removes a real contradiction with Methods. WHAT WAS REFUSED, each recorded at its own site: rewriting Methods to say the seeds WERE recorded (false of every ensemble behind every figure, and only a full re-run makes it true); compressing the conserved-quantity conjecture to "the origin is unknown", which is the wording approach.md forbids; moving the three numerical safeguards and the least-squares marginal out of Methods, where the 2026-08-04 replication audit records four independent readers stopping on exactly those two blocks; deleting the finite-difference Fisher paragraph, whose battery is in the deposit; deleting the "no probit transform" sentence, which is what stops a deposit reader misreading the column names; deleting the automatic-differentiation paragraph, which is in on Luciano's instruction; cutting the Munch conjugacy quotations, moved there as a concession on 2026-08-12; cutting the Anderson 1973 concession; and moving the dimensionless-units block to Methods, which prints after the Results. **The remaining 3,468 is still where the carve plan said it was, in the Results.** |
| 2026-08-09 | 1, reverted | | **Figure 5 stays in the body, and the rule that says so is worth writing down.** The move to an appendix was made and undone the same day. It was scored on words and it should have been scored on the programme's own hard rule, that a claim in the abstract is body material: the abstract says the filter departs "toward few channels and low noise, where the Gaussian moment closure on the occupancy is misspecified", and Figure 5 is the figure behind that claim. Two further costs the word count does not see: elife.cls renames a figure inside an appendixbox, so Figure 5 would have been cited as "Appendix 3—figure 1", and a float cannot live in an appendixbox at all, so the move also cost the float. What it bought was 222 counted words, 1.5% of the 14,500 that have to go. The subsection is reassigned to Phase 3, where compressing it to about 150 words saves most of the same and breaks no rule. GENERAL LESSON for the rest of the pass: run every relocation past the abstract before running it past the budget. |
| 2026-08-09 | 1, the Theory relocation | 17,889 | eight more blocks out of Theory, 3,835 words: the boundary-conditioned moments, the two generalization paragraphs, the three-step construction, the start-conditioned members, what separates them, both members against the exact law, and the verification tolerances to Appendix 1; the scoring equation, the initial condition and the symbol table to Methods. Theory 8,883 → 4,045, Methods 4,986 → 5,281, Appendix 1 at 4,654. Conservation checked: 3,835 out, 3,835 in plus 111 words of new connective text. Deviation from the plan, deliberate: the moment machinery went to the appendix rather than to Methods, which is where it belongs by the plan's own test (it derives γ̄ and v̄ rather than telling anyone what to evaluate) and which keeps Methods at p96 instead of p99. Not yet built. |
