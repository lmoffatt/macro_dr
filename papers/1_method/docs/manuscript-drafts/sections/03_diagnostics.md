# Diagnostics — the brief for `03_diagnostics.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `03_diagnostics.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Convert "we tested whether the likelihood is faithful" from a claim into a procedure a
sceptic could run on their own likelihood tomorrow. It is the part a reader might reuse, and it has to
survive a statistician, which means being explicit about which identities are exact algebra, which are
first-order, and which are conventions.

**The three tests, in the wording that wins.** `../../../../_program/machinery.md` §2 owns them and its wording
beats the voice-note paraphrase on two points that matter: test 1 includes **residual whiteness**, which
is the sharpest discriminator on record (score ACF: IR ≈ 0.005 against 0.78 for the non-recursive
members), and test 3 anchors on **the model's own Gaussian Fisher**, not on "the Fisher as a proxy for
the Hessian". Do not overwrite machinery with the looser phrasing.

**The anchor, and the objection to pre-empt.** The anchor is the model's own Gaussian Fisher, computed
from its predictive mean and variance and their derivatives at no cost beyond the score. For a macro
algorithm the likelihood is Gaussian by construction, so this is the model's own Fisher and the
sandwich is exactly the information-matrix-equality test in the model's own frame. It is PSD by
construction, so the diagnostic never has to be rescued from an indefinite anchor; and the numerical
finite-difference Fisher is computed separately only to gauge how faithful the Gaussian one is, not as
the anchor, because in this regime it is widely indefinite per replicate. **State that reason.** A
referee who knows the sandwich literature will otherwise ask why the Hessian was not used.

**The circularity answer, which must be in the paper and not just in our heads.** The Gaussian Fisher
is the information the likelihood *claims*; the score covariance is what it *delivers* under data from
the true process. They are different objects and their mismatch is what misspecification means. The
anchor being internal is not a circularity, it is the definition.

**Stays out.** The evidence correction beyond one motivating paragraph. The entire posterior-distortion
family. The SAFE/MARGINAL/UNRELIABLE/BIASED categorization, which gates on the likelihood Hessian's
spectrum and belongs to the posterior framework. The implementation trust region beyond a Methods
paragraph.

**Verify.** The `C = C_s^{1/2} R C_s^{1/2}` identity: prove it under a stated condition or demote it to
the determinant version, which is exact. Do not print it as written. n_sims uniformity across every
cell of every heatmap. **The direction convention, once, against the producer** — see the standing
blocker below.

---
