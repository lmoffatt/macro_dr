# Discussion — the brief for `05_discussion.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `05_discussion.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Say what it means for someone about to fit a macroscopic recording next week, and close the
loops the Introduction opened. Three deliverables, in order of importance: **a decision rule** (not
"MacroIR is best" but the conditions, which is what gets cited); **the two answers**, to Milescu 2005 on
whether non-filtering estimates are intrinsically biased and to Moffatt 2007 on the need for a
time-averaged formulation, both named explicitly; and **the honest perimeter**.

**The trap.** This is where a paper drifts back into being about its method. The subject of this
Discussion is *the distortion*; MacroIR is where the distortion happens to be small. Keep the grammar
that way and the paper stays what it claims to be.

**The three points** (2026-07-28 voice notes): MacroIR is calibrated over almost the whole measured
plane, so one could simply always use it; where it fails, which is few channels and telegraph
(non-Gaussian) noise; and why the intermediates fail, which is that the double conditioning at both
interval ends is what is load-bearing and conditioning on the start alone contributes nothing.

**Open, and both are strategic rather than editorial.**
- **The Comm Biol coupling.** The published P2X2 evidence ranking was run with the buggy `gvar_i`. This
  paper's second significance claim leans on that ranking as its only real-data demonstration, so the
  two are linked in print. Options: say nothing and handle the erratum separately; one neutral sentence
  that the numbers are being re-checked; or the strong version, in which this paper is the machinery
  that would have caught it and shows it working on the author's own published result. Decide before
  drafting, not after.
- **The 6.4 margin.** Comm Biol's winner beats its symmetric counterpart by a Bayes factor above 5000,
  which is immune to almost anything, but beats the top non-conformational alternative by only 6.4
  (5.5-6.7). That is the scale a systematic likelihood distortion can move, and this paper now has the
  map that says how much distortion that regime carried. Computing it is a strong result either way,
  and better found by us than by a referee.
- Does the research-program framing get a paragraph? Recommendation: two or three sentences, not a
  section. The paper must stand without it.

**Verify.** Every number in the decision-rule table against the freeze recompute
(`decisions/D-4_ranking_verdict.md`), the direction of the MR result in particular.
**[SETTLED SINCE]** The old verify list asked whether NMR is unbiased, on which a whole paragraph
rested; NMR is now a supplement member and measurably indistinguishable from NR, so that paragraph
needs rewriting rather than checking.

---
