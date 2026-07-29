# Introduction — the brief for `01_introduction.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `01_introduction.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Take a reader who fits the mean by least squares and show them, in their own vocabulary,
what the fluctuations carry, why using them takes a likelihood, and why no one can currently tell
whether such a likelihood is working.

**The five moves** (2026-07-28 voice notes): independence as the foundational assumption of least
squares, and what violating it costs; Markov chains as the way to model dependence while still yielding
record-universal constants tied to biophysical structure; the ladder of methods; why least squares still
reigns (you can see, point by point, whether it predicts, which is visually convincing, while a
recursive filter follows the data closely and holds no surprises, so a working method and a bug look
alike); and why autoregressive alternatives do not solve it.

**Constraints.**
- **Do not re-announce MacroIR.** Cite it as prior work in the same breath as Milescu 2005 and
  Münch 2022. The moment the Introduction explains how it works, the referee says "you published this".
- **The NSFA naming trap.** Non-stationary fluctuation analysis *is* widely used and *does* use the
  variance. Name it and distinguish it or an electrophysiologist will think the gap is already filled.
  Stepanyuk 2014's own words end the objection: *"the unitary current is virtually the only parameter
  that can be reliably obtained from this type of analysis"*, and *"kinetic rates have never been
  estimated for any synaptic receptors in their intrinsic environment"*.
- **On ARIMA**, four points, all derivable from machinery already in the paper: an ARMA error model is
  stationary by construction while the gating covariance tracks the mean current and restarts at every
  jump; it contains no channel count and no unitary current; its timescales are free where the Markov
  model ties them to the rates that generate the mean; and whiteness is not falsification, since enough
  ARMA terms whiten anything.

**Stays out.** The mechanics of MacroIR (that is Theory). The Fisher-to-zero result. The
research-program framing, which reads as a grant proposal and belongs in the Discussion. The Comm Biol
biology beyond one citation. Any claim about experimental data.

**Open.** How much statistics vocabulary in the "why nobody could test this" move: recommendation is
all three terms, each glossed in physical terms in the same sentence, at a cost of about forty words.

**Verify.** Del Core & Mirams 2025 and Owen & Mirams 2025 quotations against the version of record.
**[SETTLED SINCE]** The old top verify item was the "no published likelihood integrates the acquisition
window" claim; it is dead as written and survives only scoped to macroscopic-N by a scaling argument
(`decisions/D-3_novelty_claim.md`).

---
