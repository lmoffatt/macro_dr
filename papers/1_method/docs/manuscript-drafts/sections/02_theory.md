# Theory — the brief for `02_theory.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `02_theory.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Make the methods **one object graded by how much of the interval structure the likelihood
uses**, so the Results read as a traverse of a coordinate system rather than a bake-off among unrelated
codes.

**The one-sentence spine.** A macroscopic current is the sum over a population of channels, each a
continuous-time Markov chain; the exact likelihood of a time-averaged recording would need the
distribution of the interval-averaged conductance of the whole population, which is intractable, so
every practical likelihood replaces it with a Gaussian, and the methods differ in *what they condition
that Gaussian on* and *how they account its variance*.

**[SETTLED SINCE]** The archived plan scoped this section to `R, MR, VR, IR` and handed least squares
and the non-recursive members to "paper 2". There is no paper 2. The section covers the root question
and the ladder together (`../../../decisions.md`).

**Stays out.** k-state generality beyond what two states need. The Taylor variants, and keep them clear
of `VR`: the engine flag is `taylor_variance_correction`, `VR` is not a Taylor variant, and Methods
must say so or a reader will conflate them. The PSD trust coefficients. The Kalman prior-art
connection, which reads defensively here and belongs in the Discussion. The distortion machinery
itself, which is `../../../../_program/machinery.md` and is presented in Diagnostics.

**Open.** How much algebra survives in the main text. Proposal: four displayed equations, the
observable as an interval average, the boundary-conditioned conductance, the total-versus-residual
variance difference (the MR→VR step) and the boundary cross-covariance (the VR→IR step). Test: could an
electrophysiologist read the Results with only those four?

**Verify.** The emission-variance decomposition against the code, not the archived drafts, which
predate the `gvar_i` fix. The total-versus-residual variance forms and the boundary cross-covariance
term, re-derived from current code before printing: `VR`'s definition *is* the residual form, so if the
code's `gvar_i` is not what Theory claims, the roster is wrong and not just the prose.

---
