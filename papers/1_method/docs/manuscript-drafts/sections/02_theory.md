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

**[OPEN 2026-08-05, Luciano] The paper prints the derivation route for `gmean_ij` / `gvar_ij`; the code
runs a different one.** Both compute the same object in exact arithmetic, so this is not an error in the
mathematics. It is a provenance gap: a replicator following `02_theory.tex` writes code that is not the
code the results came from, and two of the differences are estimator choices, not implementation detail.

*What the paper prints* (`02_theory.tex:204-249`): the tilted-generator route, augmented generators
`Q_2` (2K x 2K) and `Q_3` (3K x 3K), `F_1` and `F_2` read off their off-diagonal blocks, `\citet{vanloan1978computing}`,
then `Gammabar = F_1 / (Delta P)` and `Vbar = 2 F_2 / (Delta^2 P) - Gammabar^2`, with the endpoint average
as the exact limit where `P -> 0`. Faithful to `theory/macroir/notes/Gmean_ij_gvarij/gmean_ij_gvar_by_blocks.tex`
**section 5**, "Block-matrix representation", which opens by saying the derivatives can be had *without
explicit diagonalization of Q*. The manuscript repeats that clause at `02_theory.tex:209`.

*What the code runs* (`legacy/qmodel.h`): it **diagonalizes Q**. `V`, `W` and `lambda*dt`, then the same
`U` and `W` of the note's section 4 built from **divided differences of the exponential over the spectrum**,
`E2` at `:780` and `E3` at `:844`, contracted with `WgV = W Gamma V`: `gtotal_ij = V (WgV o E2) W` at `:1628`,
`gtotal_sqr_ij = 2 V (WgV WgV o E3) W` at `:1640`. Then, instead of the quotient by `P_ij`, a
**conjugate-prior shrinkage**, `calc_g_ij_bayes` at `:1189` with its derivation in the comment at `:1102-1144`:
`gmean_ij = (gtotal_ij + prior*eps)/(P + eps)`, `prior_gmean_ij = (g_i + g_j)/2`, `prior_gsqr_ij = (g_i^2 + g_j^2)/2`,
pseudo-count `eps = eps_mach * kappa_F(V)` on the eigen paths and `10 sqrt(N) eps_mach` on taylor/schur.
Source: `theory/macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md`.

*Three specific consequences.*
1. `02_theory.tex:209` says these objects are available "without diagonalizing **Q**". True of the printed
   route, **false of the implementation**. As written it tells a replicator that the reference results avoided
   a step they in fact took, and the conditioning of `V` is exactly what the pseudo-count is scaled by.
2. The paper gives the endpoint average as an exact **limit** at `P = 0`. The code runs a **smooth blend**
   toward that same value with weight `eps`. They agree at `P = 0` exactly (the prior *is* the paper's limit),
   so nothing in the pre-pulse segment is affected; they differ for small-but-nonzero `P`, which the paper
   does not mention at all.
3. The comment at `02_theory.tex:258-260` argues against the `1/sqrt(P^2+eps^2)` smoothing. That is **not the
   smoothing the code uses**. The rejected form and the implemented form are different estimators, so the
   comment reads as a settled question that was never actually asked of the implementation.

*Why no audit caught it.* The two replication audits of 2026-08-04 (`../../../audits/raw_replication_20260804.json`,
`..._run2_20260804.json`) asked whether a student can **build** the paper from the paper. Run 2's `members`
stage reports that the van Loan construction worked and gave `Gammabar`, `Vbar`, `gbar0`, `G` and `H` directly;
it then stopped on the trust coefficient and on `NR`. Nothing in either audit asked whether the printed route
is the route that produced the figures. That question has never been put to the manuscript.

*Fix, cheap.* Keep the block construction as the definition, it is the clearer derivation. Add one Methods
sentence naming the implemented route and the shrinkage with its pseudo-count, and strike or qualify
"without diagonalizing Q". The full derivation of both routes and of the shrinkage is already written in the
two theory notes and is SI material, which connects to the standing promise below.

**[OPEN 2026-08-05] The promised Supplementary Information has no file.** `02_theory.tex:261-263` states that
the derivation of why `F_1` and `F_2` are those derivatives, and the argument for the endpoint average, "stays
out of the body and belongs in the Supplementary Information". There is no SI document anywhere in
`manuscript-drafts/`: `elife_paper.tex` inputs eight section files and `08_appendix_members.tex` as Appendix 1,
and that is all. Either the SI gets created (the material exists, in `gmean_ij_gvar_by_blocks.tex` and
`bayesian_prior_regularization_of_Qdt.md`) or the comment stops promising it.

**Verify.** The emission-variance decomposition against the code, not the archived drafts, which
predate the `gvar_i` fix. The total-versus-residual variance forms and the boundary cross-covariance
term, re-derived from current code before printing: `VR`'s definition *is* the residual form, so if the
code's `gvar_i` is not what Theory claims, the roster is wrong and not just the prose.

---
