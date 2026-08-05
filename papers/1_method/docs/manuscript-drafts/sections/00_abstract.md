# Abstract — the brief for `00_abstract.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `00_abstract.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Make two words easy for an editor to assign. Under eLife's model there is no private
accept decision; there is a permanent public assessment written with a controlled vocabulary, and the
abstract is where it comes from. Significance runs landmark > fundamental > important > valuable >
useful; strength of evidence runs exceptional > compelling > convincing > solid > incomplete >
inadequate. Target: **compelling** on evidence, earned by the ground truth being an exact simulation of
the process the likelihood approximates; and **valuable to important** on significance, earned by
saying who can act on the result.

**Constraints.**
- 200 to 220 words. **[SETTLED SINCE]** The old rule was "write to 150, treat 200 as the ceiling".
  Measured against 41 recent eLife research articles: median 197, 39% over 200, and the
  computational/structural/biophysics subset has a median of **210**. Do not spend words getting to 200.
- Name the biological system in the first sentence.
- **The paraphrase test.** Every technical term must survive being restated by an editor who is not a
  statistician, because that restatement is published and permanent. If an editor cannot write the
  assessment without the word "score", gloss it.
- **[SETTLED SINCE]** The old slot table was Nature's funnel. Measured, eLife does not use it: only
  half of abstracts contain "Here we", and among those, method verbs beat result verbs 13 to 8; only
  27% carry an explicit comparison to prior belief. What *is* the eLife norm, and what this abstract
  already has, is the broader-perspective closer (68%).

**Impact Statement** (submission-form field, not part of the PDF): one sentence, 15 to 30 words, third
person, complements the title, no "We show". **[OPEN]** three drafts are in the archived `abstract.md`
and none is post-merge.

**Open.**
- The significance sentence. The current close, "this supplies a test… and a map…", describes the
  deliverable rather than what changes. The formulation to work from (Luciano, 2026-07-29): *a valid
  likelihood is what enables separating inter-experiment variability from stochastic error, and the
  evaluation of Bayes factors.* Order them by evidence, and make the Bayes half concrete via the
  distortion matrix rather than generic. See `../../../decisions.md`.
- **The NSFA collision.** The map's largest region is sold on the unitary current and the channel
  count, and `01_introduction.tex` concedes both to non-stationary fluctuation analysis while arguing
  NSFA returns no kinetic scheme. Resolve before drafting: what is new is the *combination*, a
  calibrated joint estimate of rates and amplitudes from a single non-stationary record, where the
  1980-81 covariance methods needed ensembles of 256 to 504 repeated sweeps and stated no uncertainty.
- **[OPEN 2026-08-05, Luciano] "The test is classical:" reads wrong at the seam.** Raised against the
  live sentence 2-to-3 transition: *"... a filter tracks the data even when the uncertainty it reports
  is wrong. The test is classical: at the true parameters the score, the log-likelihood's derivative,
  averages to zero ..."*. Two things in one. (i) The transition jars: nothing before it has said a test
  is wanted, so "The test" arrives with no antecedent and the definite article asserts a shared referent
  the reader does not have. (ii) The claim as heard is *"the test for likelihoods is classical?"*, which
  is a stronger and more doubtful statement than the one intended. What is classical is the information
  matrix equality and its use as a misspecification diagnostic (Huber 1967, White 1982, Golden et al.);
  that there is a standing classical test of whether an approximate likelihood for **this** class is
  faithful is not established, and rule 6 only asks that the identities not be claimed as new, not that
  a test be asserted as pre-existing.
  **Do not simply re-flip it.** The phrase has already oscillated twice and the history is in
  `00_abstract.tex:253-258` and `:332-336`: it was cut for word budget, then restored *because* two
  prior review rounds produced the rule-6 objection. Deleting it re-opens rule 6; keeping it produces
  this objection. A third form is needed, one that hands rule 6 its disclaimer without the bare definite
  article, e.g. attaching "classical" to the identity rather than to a test. Whatever is chosen, the
  gloss "the log-likelihood's derivative" stays (the brief's paraphrase test) and "variability should
  match" stays out (not the identity that was measured).
  Checked while registering: "at the true parameters" is **not** part of this problem. The abstract
  sentence states the correctly-specified identity, so the true parameters are the right anchor there;
  the pseudo-true caveat that `01_introduction.tex` header note (d) guards applies to the approximation's
  own score and is handled in `03_diagnostics.tex:3`.

**Verify.** Every quoted number at a fixed n_sims. The direction of the non-recursive result, which
inherits the unverified sign convention (`% TODO-SIGN` in the `.tex`). No channel-number threshold
unless D-J is closed — and note the manuscript already prefers the measured reachability floors
(N_ch below 17 for the channel count, below 52 for the opening rate) over "about 100".


**

---
