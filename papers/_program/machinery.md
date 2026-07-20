# The validation machinery

> Updated: 2026-07-20. Merged from `03_metrics_diagnostics.md`, `correction_idm_reconstruction.md`
> and the definitional half of `diagnostics_plan.md`.
>
> **Owns:** the diagnostics, their definitions, sign conventions, thresholds, what composes with what,
> the numerical conventions, and the known hazards. All three papers cite this file; none restates it.
> **Does not own:** how paper 1's Diagnostics *section* is written (`1_method/diagnostics.md`), nor
> the model and simulator (`model_and_sim.md`), nor which CSV holds what (`provenance.md`).
>
> Primary sources: `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex`
> and `theory/macroir/docs/Gaussian_Fisher_Distortion_Family.md`.

## 1. Premise and licence

An approximate likelihood is a **misspecified model**, and misspecification has a classical signature
(Huber 1967, White 1982): the score no longer has zero mean at the truth, and the information matrix
equality fails, so the covariance of the score and the information the model reports stop agreeing.

Both are testable **only if you can draw data from the true process**. Here we can: the forward
process is an exact continuous-time Markov chain, simulated exactly by uniformization, so there is
unlimited data from the process the likelihood claims to describe, with the true parameters known.
That asymmetry — the process is exactly simulable while its likelihood is not — is the whole licence
for this program.

**The test is classical and is not ours to claim.** What is new is applying it to macroscopic-current
likelihoods, reading the resulting matrix directionally, and the maps that follow. Any sentence of
the form "there was no way to test whether a likelihood is faithful" is false and a statistically
trained reviewer will catch it.

## 2. The three tests, in increasing strength

**T1 — Residuals.** r_t = (y_t − ŷ_t)/σ̂_t must have mean 0, variance 1, and no autocorrelation.
Tests the *predictive* distribution one step at a time. Cheapest, and **the only one an
experimentalist can run on real data**, because it needs no known truth. Say so in every paper: it
tells the reader which part of the machinery transfers to their own recordings.

**T2 — Score.** E[∇ℓ(θ*)] = 0 at the true parameters. Tests the first moment of the gradient; it is
what governs bias in the estimate. Projected through the parameter covariance it becomes a bias
vector in parameter units, **DIB**, b = H⁻¹ E[S] — the reader-facing "how wrong is the estimate".

**T3 — Information.** Cov[∇ℓ] = H, the information the likelihood reports. The information matrix
equality, and the strongest of the three: **an approximation can pass T1 and T2 and still fail this
one by an order of magnitude.** T2 and T3 need a known truth, therefore a simulator.

## 3. The anchor

H is the model's **own Gaussian Fisher**:

  H = Σ_t E[ (∇μ_t)(∇μ_t)ᵀ/σ_t² + (∇σ_t²)(∇σ_t²)ᵀ/(2σ_t⁴) ]

computed from the likelihood's own predictive mean and variance and their derivatives, at no cost
beyond the score. For a macro algorithm the likelihood is Gaussian by construction, so this is the
model's own Fisher and H^(−1/2) J H^(−1/2) is exactly the information-matrix-equality test in the
model's own frame.

Three properties, all of which belong in the papers:

- **PSD by construction**, so the diagnostic never has to be rescued from an indefinite anchor.
- Costs nothing extra.
- The **numerical finite-difference Fisher is not the anchor.** It is computed only to gauge how
  faithful the Gaussian one is (the F-vs-G bridge). The reason is empirical and must be stated: in
  this regime the per-replicate FD Fisher is widely indefinite; pooling and eigenvalue-scaling
  established that the *true* bias is only in NR and that the rest is finite-difference noise.
  **If this is not stated, a reviewer who knows the sandwich literature will ask why the Hessian was
  not used.**

**The circularity objection, and its answer.** Comparing the approximation to its own Fisher is not
circular. H is the information the likelihood *claims*; J is the information it *delivers* against
data from the true process. Different objects, computed from different things, and their disagreement
is the definition of misspecification. Write that sentence.

## 4. The distortion matrix, and the direction convention

  **IDM = C = H^(−1/2) J H^(−1/2)**,  J = J_total = Cov(S), including cross-interval terms.

C = I under correct specification. The eigen-directions say *which combinations* of parameters are
misreported, which no scalar summary can.

> **C_ii > 1 means the likelihood under-reports the uncertainty in parameter i, i.e. over-confident.**
>
> **Every verdict in every paper depends on this convention.** It is currently asserted in three
> places with two different signs. **Verify it once against the code and fix the other copies.**

**It is a tool, not a complaint.** The same matrix repairs the error bars: the sandwich
Σ = H⁻¹ J H⁻¹ = H^(−1/2) C H^(−1/2) is the corrected parameter covariance, and the corrected ellipse
lands on the empirical cloud (ratio ≈ 1 in every measured cell of all four non-IR algorithms). A
diagnostic that also fixes what it diagnoses is a much stronger result than one that only scores
algorithms.

## 5. The decomposition, and what actually composes

Two failure modes, two matrices, with H the anchor (F_b or G_b), J_s = J_sample, J = J_total:

  IDM = H^(−1/2) J H^(−1/2)   SDM = H^(−1/2) J_s H^(−1/2)   CDM = J_s^(−1/2) J J_s^(−1/2)

- **SDM, sample / geometric distortion** — per-sample non-Gaussianity, in the H frame.
- **CDM, correlation distortion** — temporal correlation, in the J_s frame. Missing higher moments
  appearing as ghost state-correlation.

### The identity, corrected

`IDM = SDM^(1/2) · CDM · SDM^(1/2)` **is false.** The exact identity is

  **IDM = K · CDM · Kᵀ,  K := H^(−1/2) J_s^(1/2)**   (exact, always; verified to 1e-16)

K satisfies K Kᵀ = SDM, so it is a legitimate factor of SDM in the congruence sense, but it is **not**
the symmetric principal square root: K² ≠ SDM and K is not symmetric (a product of two SPD matrices is
symmetric only if they commute). The two differ by an orthogonal factor, K = SDM^(1/2) Q, with Q = I
**iff H and J_sample commute**. The printed identity silently sets Q = I.

Being SPD is not the issue; every matrix in the chain is SPD. **Non-commutativity is.**

### What composes and what does not

**Composes, unconditionally:** the determinant, hence

  **log det IDM = log det SDM + log det CDM**

with no assumption. Verified on data (IR, N_ch = 10⁴, noise 0.1, pool): 0.0661 vs 0.0662. This is the
version the maps use, and the sentence the figures support is that the distortion of the information
**volume** splits exactly into a per-sample part and a temporal part.

**Does not compose:** trace, eigenvalues, Frobenius norm, and the per-parameter diagonal. In a 2×2
test with H anisotropic (100:1), J_s rotated 45°, J rotated 30°: diag(C) = [16.5, 0.455] against
[14.16, 0.103] for the printed product, a factor 4.4 in the second parameter; trace 16.96 vs 14.26;
leading eigenvalue 16.88 vs 14.18; relative Frobenius error 22%. With near-isotropic random SPD
matrices the same error is 0.6%. **The size of the error is exactly the H-to-J_s frame misalignment,
which is the anisotropy the papers want to report.**

### Empirical size (433ed13, battery_pool, noise 0.1, bootstrap median)

| algo | N_ch | ‖IDM − I‖_F | ‖REC − IDM‖/‖IDM‖ | ‖Q − I‖_F |
|---|---|---|---|---|
| IR | 10 | 1.00 | 3.1% | 0.30 |
| IR | 10000 | 0.13 | 0.06% | 0.03 |
| R | 10000 | 0.73 | 0.34% | 0.03 |
| MR | 10000 | 1.76 | 0.64% | 0.04 |
| NMR | 10000 | 262.6 | 3.7% | 0.13 |
| NR | 10000 | 224.4 | 9.1% | 0.20 |

The reconstruction error tracks ‖Q − I‖ row by row, as predicted. It is **not** numerical noise and
more bootstrap replicates will not remove it.

### A by-product worth keeping

**‖Q − I‖ is a meaningful observable.** Q = SDM^(−1/2)·K is the rotation carrying the J_sample frame
into the Fisher frame; it is the identity exactly when the per-sample distortion is aligned with the
model's information geometry. It orders the algorithms the same way the distortion does. If the
two-path reconstruction panel is kept rather than dropped, this is what it should plot.

### Pending code fix (C1)

`src/core/likelihood.cpp:3501` implements the wrong version via
`apply_sqrt_congruence(W_SDM, cdm)`. Replace with a congruence by K, whose factors are already in
scope as `PSDDecomposition`s: `W_F_b` (`inv_sqrt` → H^(−1/2)) and `W_Js` (`sqrt_vals` → J_s^(1/2)).
New helper `apply_factor_congruence` beside the existing congruence family in
`legacy/lapack_headers.h` (`apply_normalized_congruence` :2844, `apply_inverse_congruence` :2865,
`apply_sqrt_congruence` :2886). K mixes two retained bases, so the result must be embedded back
through the H basis, and if either subspace is empty the result stays **undefined**, not zero-filled.
After the fix `idm_reconstituted == idm` to machine precision, turning an unexplained near-miss into a
genuine self-check.

**Also wrong, off the active path:** `Lapack_C_h_R_C_h` (`legacy/lapack_headers.h:2634`) uses the
Cholesky factor of SDM. L Lᵀ = SDM too, and L is also not K, so it is wrong the same way by a
different orthogonal. Fix or delete; do not leave a third variant of one identity in the tree.

**No published or plotted number changes.** IDM, SDM and CDM are each computed directly from
(H, J, J_s), never through one another (`likelihood.cpp:3398`, `:3474`, `:3488`; Gaussian twins
`:3575`, `:3610`). No re-run, no re-plot.

**Do not present log-det additivity as a finding.** It is an identity. It is a fine numerical guard —
it correctly flags the indefinite NR/NMR cells where the log-dets are ill-defined — but a
"residual = 0.000" line is a tautology, not evidence that the decomposition is meaningful.

## 6. Hazards

### 6.1 Every scalar summary of C is biased, and the bias depends on n_sims

C is estimated from finitely many simulated recordings. Every scalar computed from Ĉ (log det,
affine-invariant distance, Frobenius norm, the diagonal) is a **nonlinear** function of it, so by
Jensen it is biased, with the bias scaling as 1/n_sims. This is the classical problem with White's
information matrix test.

- **A distortion magnitude that shrinks as you run more simulations is measuring the estimator, not
  the algorithm.** Any statement "the distortion is X" must name its n_sims, or be debiased.
- **No map may mix cells with different n_sims.** The program's grid is ragged: band-A cells at 10⁴,
  the 2026-07-20 fill at 1000, `433ed13` also holding 200, micro cells at 100/1000/10⁴. A heatmap
  drawing across them shows partly the Jensen bias as a gradient. **This is the single most likely way
  for a wrong result to reach print**, and after the three-paper split it sits on more than one
  headline figure.
- **The principled fix is identified:** a debiased quadratic form (Wald T² minus its degrees of
  freedom, analytic floor = dof). `Wald_T2` is already wired in the code; it needs only the p-value
  and the (T² − p)/N_groups normalization. The signed log-det (Log-Det GIMT) is canonical but
  directional and anisotropy-blind, so it is a complement, not a magnitude.
- **[Q]** In scope, or report the raw statistic at fixed n_sims and state the bias as a limitation?
  Recommendation unchanged: **check first whether the final maps are uniform in n_sims.** If they are,
  one Methods sentence covers it. If not, it must be fixed, and it is a re-plot, not a re-run.

### 6.2 Two diagnostics do not mean for LSE what they mean for a likelihood

Wherever an LSE cell sits beside a likelihood cell (paper 2), annotate both:

- **r̄²_std ≡ 1 is a tautology for LSE.** Σ r_std² = n identically once σ̂² = SSE/n, so any
  residual-variance panel reads "calibrated" for LSE by construction and carries no information.
- **F = Var(score) requires homoscedasticity**, which LSE assumes and the exact simulator violates
  (gating variance N·gSg peaks in the transient and vanishes at the plateaus). The identity therefore
  **breaks** for LSE, and the size of the break is a **result**: the classical method's overconfidence,
  measured. Related: the shared (n/SSE) factor makes LSE's per-interval scores non-martingale through
  a rank-1 common mode.

### 6.3 The singular-Fisher convention

H can be near-singular where a parameter is unidentifiable — and that happens **by design** once the
open population relaxes. The convention (supplement §8) is to work in the retained eigenspace above a
numerical cutoff and to leave the quantity **undefined**, not finite, when the retained eigenspace is
empty. Undefined-rather-than-wrong is correct, and it is why the maps have grey cells.
**Every caption must say grey means undefined, not zero.**

## 7. Effective sample size

κ = Var(ΣX)/ΣVar(X), so T_eff = T/κ. The reader-friendly translation of the correlation distortion: a
non-recursive likelihood believes it has T independent observations when it effectively has T/κ.

**If one number from this program gets quoted in a talk, it is this one.** **[Q]** It is not clear κ
is reported in any current figure. If not, it is a cheap addition and worth more to a reader than
another matrix norm.

## 8. Thresholds (open)

- Candidate "valid": distortion < 1.1 per parameter (≈10% error); bootstrap error < 1.1. The repo
  currently carries 1.1, 1.15 and 1.5 in three places.
- Residual placeholders until empirical variability is seen: |mean(r)| < 0.1; 0.8 < var(r) < 1.2;
  max ACF (lag>0) < 0.1; max |mean(score)/sd(score)| < 0.2 per parameter.
- Start rank-based (best→worst per cell); add absolute cutoffs once the anchor grid is final.

## 9. Optional

**PIT:** u_t = F_t(y_t) ≈ Uniform(0,1) for a correct predictive. A good extra check, not one of the
main indicators.

## 10. Explicitly out of scope for all three papers

- **The evidence correction** (Laplace peak + volume, α-calibration). Motivation only, one paragraph,
  derivation deferred to a later component. Near the maximum the distortion propagates to the evidence
  through two scalar summaries of eig(C): volume Δ log Z = ½ log det C (geometric mean, O(1) in T,
  decisive for Bayes factors) and peak / effective samples α⋆ = p / tr C (arithmetic mean, O(T)).
- **The posterior information-distortion family** (five substantial supplements). A later component,
  and well-written enough to be a paper of its own.
- **The safety categorization** (SAFE / MARGINAL / UNRELIABLE / BIASED). It gates on the spectrum of
  the likelihood Hessian and belongs to the posterior framework. It appears in no figure here and
  must appear in no text here.
- **The implementation-level trust region** (α_μ, the simplex projection, the PSD guard). Methods at
  most, one paragraph, and only because the code runs it.

## 11. Verify before any submission

- **The direction convention** (§4), once, against the code.
- **n_sims uniformity** across every cell of every heatmap (§6.1). Highest-risk item in the program.
- **The corrected identity** (§5) reaches the supplement and the code, and the false one is printed
  nowhere.
- **Whether `Wald_T2` is computed on the final runs**, and whether its debiased form changes any map's
  ordering. If it does not, say so in one sentence and keep the simpler statistic.
