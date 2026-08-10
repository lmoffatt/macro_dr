# Bringing paper 1 to length

Written 2026-08-09. Target and evidence: `../_program/elife_main_text_length_survey.md` (413
published eLife articles, the two subject areas this programme samples). Instrument:
`../_program/wordcount.py`, which counts the main text eLife's way (comments, captions, equations,
Methods, back matter and appendices all excluded).

## The target

**8,000 counted words of main text, with 8,500 as the line above which the draft is not accepted.**

8,000 is the 90th percentile of published articles. 7,000 is the 82nd and buys nothing an editor
distinguishes, while forcing the Discussion below the published median of 1,184, which is where the
paper's concessions live. Above 8,500 the paper is in the top 8% and needs an argument for it that
nothing in the manuscript currently makes.

Baseline, measured today: **22,727**. That is 2.83× the p90 and 1.97× Münch e62714, the longest
comparable article in the niche.

| section | now | target | delta |
|---|---|---|---|
| Introduction | 2,026 | 1,200 | −826 |
| Theory | 8,883 | 1,700 | −7,183 |
| Diagnostics | 2,392 | 700 | −1,692 |
| Results (prose) | 5,371 | 3,200 | −2,171 |
| Discussion | 4,055 | 1,400 | −2,655 |
| **total** | **22,727** | **8,200** | **−14,527** |

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

Out of Results: the subsection "Where the calibrated member departs" with Figure 5 and its caption,
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
   measurement says they do not for the least-squares arm. One of the two has to change.
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
| 2026-08-09 | 1, the Theory relocation | 17,889 | eight more blocks out of Theory, 3,835 words: the boundary-conditioned moments, the two generalization paragraphs, the three-step construction, the start-conditioned members, what separates them, both members against the exact law, and the verification tolerances to Appendix 1; the scoring equation, the initial condition and the symbol table to Methods. Theory 8,883 → 4,045, Methods 4,986 → 5,281, Appendix 1 at 4,654. Conservation checked: 3,835 out, 3,835 in plus 111 words of new connective text. Deviation from the plan, deliberate: the moment machinery went to the appendix rather than to Methods, which is where it belongs by the plan's own test (it derives γ̄ and v̄ rather than telling anyone what to evaluate) and which keeps Methods at p96 instead of p99. Not yet built. |
