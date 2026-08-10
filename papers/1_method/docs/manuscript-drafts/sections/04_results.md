# Results — the brief for `04_results.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `04_results.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Walk the reader from one filter step to the usage map without asking them to take a claim
on trust. Each figure answers the objection the previous one raises.

**The through-line, six body figures** (`../../../decisions.md`, "The figure set"): the filter step along the
cost ladder; recovery clouds at one cell; the calibration cascade in time; the design space split by
moment; the information budget per parameter; and the usage map. **[SETTLED SINCE]** The archived plan
had five figures and an `R`-versus-`IR` arc.

**Framing the maps force.** The three named regimes (microscopic, telegraphic, Gaussian) are the two
approximations' **asymptotic corners, not territories with borders**. The map shows a gradient; the
corners explain its direction. Write that into both Theory and Results so they do not read as two
different claims.

**The threshold is 1.15, read as a variance ratio**, which is a 7% error on the reported standard
deviation. **[SETTLED SINCE]** The repo carried 1.05, 1.1, ±15% and 1.5 in four places; the built map
chose 1.15 and glossed it, and `../../../../_program/machinery.md` §8 must now cite that rather than propose
1.1.

**Claims the data do NOT currently support. Say none of these until the gap is closed.**
- *"Sample × correlation = total distortion."* Asserted, never plotted.
- *The N\* ∝ noise law as a validated prediction.* N\* is measured but never overlaid on its predicted
  law, and the fitted exponent is 0.80 with R² = 0.59, which is neither a clean 1/2 nor 1. Report the
  measurement; do not present the law as confirmed.
- *The micro-versus-macro comparison.* Blocked on data and out of scope; it must not sneak in through a
  figure.
- *Anything from the flagged-wrong scripts*: the argmax ridge (it is noise), the `relCI_by_Nch` panel
  that contradicts its own sibling, and the schematic that hard-codes a boundary contradicting the
  paper's own N\* result.
- **[ADDED 2026-07-29]** *Separating inter-experiment variability from stochastic error.* This is the
  paper's best significance claim and it is currently a capability statement, not a measurement. It is
  cheap to demonstrate from data already on disk: every run is at fixed θ, so Var(θ_true) = 0 is known
  by construction and the between-experiment spread is entirely stochastic; inject a known variability
  onto the per-group MLEs and ask each method how much real variation it infers. Until that is run, the
  claim is an implication and must be worded as one.

**Verify.** The ±15% self-consistency fractions on whatever anchor the final figures use. The
overconfidence gap: confirm it is the recursive/non-recursive split and not a recursive/averaging split,
because the sentence reporting it will be read as a claim about *which* approximation matters.

---
