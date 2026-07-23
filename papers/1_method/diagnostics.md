# Diagnostics section — the validation machinery, what it claims, and two things in it that are not what they look like

> Working doc, same genre as `abstract.md`. Opened 2026-07-14.
> Covers the manuscript's third section, *"Diagnostics: testing a likelihood against its simulation"* (`docs/manuscript-drafts/elife_paper.tex`), which is the paper's genuinely new contribution and therefore the section most exposed to a methods reviewer.
> Source of the definitions: `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex` (373 lines, PDF built) and `theory/macroir/docs/Gaussian_Fisher_Distortion_Family.md`. Operational summary: `../_program/machinery.md`.

## The job

This section has to convert "we tested whether the likelihood is faithful" from a claim into a procedure a sceptic could run on their own likelihood tomorrow. It is the part of the paper that a reader might actually reuse, and it is the reason *Tools and Resources* is a live article-type option (`abstract.md`, decision 4).

It also has to survive a statistician. That means being precise about which identities are exact algebra, which are first-order approximations, and which are conventions. The supplement is honest about this internally; the manuscript must be too.

## The logic, in the order it should be presented

**Premise.** An approximate likelihood is a misspecified model. Misspecification has a classical signature (Huber, White): the score no longer has zero mean at the truth, and the information matrix equality fails, so the covariance of the score and the negative expected Hessian stop agreeing. Both are *testable* if you can draw data from the true process.

**The licence.** We can. The forward process is an exact continuous-time Markov chain and it is simulated exactly (uniformization). So we have unlimited data from the process the likelihood claims to describe, and the true parameters are known.

**Three tests, in increasing strength.**

1. **Residuals.** Standardized residuals r_t = (y_t − ŷ_t)/σ̂_t must have mean 0, variance 1, and no autocorrelation. This tests the *predictive* distribution one step at a time. It is the cheapest test and it is the one an experimentalist can run.
2. **Score.** E[∇ℓ(θ*)] = 0 at the true parameters. This tests the *first moment* of the likelihood's gradient and it is what governs bias in the estimate. Projected through the parameter covariance it becomes a bias vector in parameter units (**DIB**, b = H⁻¹ E[S]), which is the reader-facing "how wrong is the estimate".
3. **Information.** Cov[∇ℓ] = H, the information the likelihood reports. This is the information matrix equality, it is the strongest of the three, and it is the one the whole paper turns on, because an approximation can pass tests 1 and 2 and still fail this one by an order of magnitude. **NR and NMR are exactly that case**, and saying so is how the section earns its existence.

**The measurement.** The failure of test 3 is not a p-value, it is a matrix:

  **C = H^(−1/2) J_total H^(−1/2)**,  J_total = Cov(S), H = the Gaussian Fisher.

C = I under correct specification. C_ii > 1 means the likelihood under-reports the uncertainty in parameter i (over-confident). The eigen-directions say *which combinations* of parameters are misreported, which a scalar summary cannot.

**Why it is a tool and not a complaint.** The same matrix repairs the error bars: the sandwich Σ = H⁻¹ J H⁻¹ = H^(−1/2) C H^(−1/2) is the corrected parameter covariance, and Figure 2 shows the corrected ellipse landing on the empirical cloud (ratio ≈ 1 in every cell of all four non-IR algorithms, `Figure_S3_caption.md`). **A diagnostic that also fixes what it diagnoses is a much stronger paper than one that only scores algorithms.** Lead with that.

## The anchor: why the Gaussian Fisher, and the objection to pre-empt

The anchor H is the model's **own** Gaussian Fisher,

  H = Σ_t E[ (∇μ_t)(∇μ_t)ᵀ/σ_t² + (∇σ_t²)(∇σ_t²)ᵀ/(2σ_t⁴) ],

computed from the likelihood's own predictive mean and variance and their derivatives, at no extra cost beyond the score.

The justification, which `Gaussian_Fisher_Distortion_Family.md` states in one sentence worth quoting nearly verbatim in Methods: **for a macro algorithm the likelihood is Gaussian by construction, so this is the model's own Fisher, and H^(−1/2) J H^(−1/2) is exactly the information-matrix-equality test in the model's own frame.**

Three consequences, all of which belong in the paper:

- It is **PSD by construction**, so the diagnostic never has to be rescued from a numerically indefinite anchor.
- It costs nothing extra (derivable from the same derivative pass).
- The **numerical finite-difference Fisher** is computed separately, and only to gauge how faithful the Gaussian one is (the F-vs-G bridge). It is *not* the anchor, and the reason is empirical: in this regime the finite-difference Fisher is widely indefinite per replicate (`project_numerical_fisher_negative_eig_figure2`). Pooling and eigenvalue-scaling arguments established that the *true* bias is only in NR; the rest of the indefiniteness is finite-difference noise. **State this. It is the honest reason for the anchor switch, and if it is not stated, a reviewer who knows the sandwich literature will ask why the Hessian was not used.**

**The circularity objection, and its answer** (also in `discussion.md`): comparing the approximation to its own Fisher is not circular. H is the information the likelihood *claims*; J is the information it *delivers* against data from the true process. They are different objects, computed from different things, and their disagreement is the definition of misspecification. Write that sentence in the paper.

## Two things in this section are not what they look like

These are the section's technical soft spots. Both are fixable in a paragraph, and both will be found by a competent reviewer if they are not.

### 1. The "multiplicative decomposition", as written, is off by a rotation (and the fix is free)

The supplement (§10) defines the **sample distortion** C_s = H^(−1/2) J_s H^(−1/2) (per-sample non-Gaussianity, in the H frame) and the **correlation distortion** R = J_s^(−1/2) J J_s^(−1/2) (temporal correlation, in the J_s frame), and writes

  C = C_s^(1/2) · R · C_s^(1/2).   (as printed)

**That matrix identity is false in general**, and the reason is worth stating precisely because it is a genuinely easy mistake. Substituting the definitions, the identity requires C_s^(1/2) = H^(−1/2) J_s^(1/2). Now, the matrix

  **K := H^(−1/2) J_s^(1/2)**

does satisfy **K Kᵀ = C_s**, so it is a legitimate square-root *factor* of C_s in the congruence sense. But it is **not the principal square root**: K² ≠ C_s, and K is not even symmetric. The two notions of "square root" coincide only in one dimension, or when H and J_s commute. That is the whole bug.

**The exact identity does exist, and it is with that factor:**

  **C = K R Kᵀ,  K = H^(−1/2) J_s^(1/2)**   (exact, always; verified numerically to 1e-16)

Since any two factors with K Kᵀ = C_s differ by an orthogonal matrix, K = C_s^(1/2) Q, the true statement is C = C_s^(1/2) (Q R Qᵀ) C_s^(1/2), and the printed version drops Q, the rotation that carries the J_s frame into the H frame.

**What survives regardless: the determinant.** det C = det J / det H, det C_s = det J_s / det H, det R = det J / det J_s, so

  **log det C = log det C_s + log det R**

holds identically, with no assumption at all.

**What does not survive: trace, eigenvalues, Frobenius norm, and the per-parameter diagonal.** In a 2×2 test with H anisotropic (100:1), J_s rotated 45° and J rotated 30°: diag(C) = [16.5, 0.455] against [14.16, 0.103] for the printed product, a factor of 4.4 in the second parameter; trace 16.96 against 14.26; leading eigenvalue 16.88 against 14.18; relative Frobenius error 22%. With random (near-isotropic) SPD matrices the same error was 0.6%. **The size of the error is the misalignment between the H and J_s frames, and frame misalignment is exactly the anisotropy the paper wants to report.**

**The fix: use K, and keep every definition.** The chain composes exactly with K, so nothing has to be redefined:

  **C = K R Kᵀ,  K = H^(−1/2) J_s^(1/2),  K Kᵀ = C_s.**

This is strictly better than the alternative of redefining the correlation factor in the H frame (R_H := C_s^(−1/2) C C_s^(−1/2), which also makes the product exact), because K leaves R **frame-independent**, and that is a property the code deliberately exploits: one CDM serves both the numerical-Fisher and the Gaussian-Fisher anchors (`src/core/likelihood.cpp:3570`).

**This is not hypothetical: the code implements the wrong version, and the symptom was already on the record.** `src/core/likelihood.cpp:3501` computes `Likelihood_Information_Distortion_Reconstituted` as `apply_sqrt_congruence(W_SDM, cdm)` = SDM^(1/2) · CDM · SDM^(1/2), and `00_plan.md` §5 records the two-path reconstruction as landing "~1 but not exactly, flagged open question". That near-miss is the missing rotation. Measured on `433ed13` (battery_pool, noise 0.1): the reconstruction error is 0.06% for IR at N_ch = 10⁴, 3.1% for IR at N_ch = 10, and **9.1% for NR at N_ch = 10⁴**, tracking ‖Q − I‖ row by row.

**Full write-up and the change list: `../_program/machinery.md`.** No figure or reported number changes; IDM, SDM and CDM are each computed directly from (H, J, J_s) and never through one another.

Two consequences for the figures, either way:

- **`figure_5_master_STATUS.md`'s "additivity exact for IR/R/MR, residual 0.000" is a tautology**, not a finding. log-det additivity holds by construction. It is a fine numerical sanity check (it correctly flags the indefinite NR/NMR cells where the log-dets are ill-defined) but it is not evidence that the decomposition is meaningful, and it must not be presented as if it were.
- **`figure_5_PLAN.md:99` ("the multiplicative claim is asserted, never demonstrated") is right, and the resolution is a definition, not a panel.** Once R_H is the definition, the claim is an identity and there is nothing left to demonstrate; what the figures then show is the *attribution* (how much of the volume distortion is per-sample and how much is temporal), which is the interesting part.

### 2. Every scalar summary of C is biased upward, and the bias depends on the number of simulations

C is estimated from a finite number of simulated recordings, so Ĉ is a noisy estimate of C. Every scalar we compute from it (log det, affine-invariant distance, Frobenius norm, the diagonal) is a **nonlinear** function of Ĉ, so by Jensen it is biased, and the bias scales as 1/n_sims (`project_distortion_measure_N_dependence`). This is the classical problem with White's information matrix test.

Three consequences the paper cannot avoid:

- **A distortion "magnitude" that shrinks as you run more simulations is measuring the estimator, not the algorithm.** Any statement of the form "the distortion is X" must say at what n_sims, or be debiased.
- **The maps must not mix cells with different n_sims.** `433ed13` contains cells at n_sims = 10000, 1000 and 200. If a single heatmap draws from cells with different n_sims, the gradient across the map is partly the Jensen bias. **Check this on the final Fig 4 and its supplements (the maps and the decomposition) before they are drawn.** This is the single most likely way for a wrong result to reach print.
- The principled fix is already identified: a **debiased quadratic form** (Wald T² minus its degrees of freedom, analytic floor = dof), and `Wald_T2` is already wired in the code. It needs only the p-value and the (T² − p)/N_groups normalization. The signed log-det (log det C, the Log-Det GIMT statistic) is canonical but directional and anisotropy-blind, so it is a complement, not a magnitude.

**[Q]** Is this fix in scope for this paper, or is the honest move to report the raw statistic at a fixed n_sims (10,000 everywhere in the production runs, which is uniform) and state the bias as a limitation? Recommendation: **check first whether the final maps are uniform in n_sims.** If they are, the problem is cosmetic and one sentence in Methods covers it. If they are not, it must be fixed, and it is a re-plot, not a re-run.

## What else must be said in this section, and briefly

- **The residual test is the one the experimentalist can run on real data**, without knowing the true parameters. The score and information tests need a known truth, and therefore a simulator. Say this explicitly: it tells the reader which parts of the machinery transfer to their own recordings and which do not. It is a generous thing to say and it costs one sentence.
- **The singular-Fisher convention.** H can be near-singular in cells where a parameter is unidentifiable (and, by Fig 3—figure supplement 1, that happens *by design* once the open population relaxes). The supplement's §8 convention is to work in the retained eigenspace above a numerical cutoff, and to leave the quantity **undefined**, rather than finite, when the retained eigenspace is empty. Undefined-rather-than-wrong is the right choice, and it is also why the maps have grey cells. **The caption must say grey means undefined, not zero.**
- **Effective sample size.** κ = Var(ΣX)/ΣVar(X), so T_eff = T/κ. This is the reader-friendly translation of the correlation distortion: a non-recursive likelihood believes it has T independent observations when it effectively has T/κ. If one number from this paper is going to be quoted in a talk, it is this one. **[Q] Is κ reported anywhere in the current figures?** If not, it is a cheap addition to Fig 3 or Fig 4 and it is worth more to the reader than another matrix norm.

## What stays out

- **The evidence correction** (Laplace peak + volume, α-calibration; `supplement_evidence_correction.tex`, 651 lines). Motivation only, one paragraph in the Discussion, derivation deferred (`../_program/machinery.md` §10).
- **The entire posterior-distortion family** (`docs/Posterior_Information_Distortion/`, five substantial supplements). Cut. It is a later component and it is well-written enough to be a paper of its own.
- **The safety categorization** (SAFE / MARGINAL / UNRELIABLE / BIASED). It gates on the spectrum of the likelihood Hessian and it belongs to the posterior framework. It does not appear in this paper's figures and it should not appear in its text.
- **The implementation-level trust region** (α_μ, the simplex projection, the PSD guard). Methods, at most a paragraph, and only because the code runs it. Not here.

## Sources, with their readiness

| Source | Lift-readiness |
|---|---|
| `Likelihood_Information_Distortion/supplement_information_distortion_main.tex` | **Supplement-ready.** §1-§5 and §9-§10 carry the definitions; §6-§7 are first-order (Taylor plus expected-Hessian) and must be *labelled* as such; §8 is the singular-space convention. |
| `Likelihood_Information_Distortion/small section Information Distortion.tex` (44 lines) | **The best candidate for the main-text paragraph.** It already contains the sentence the section needs: if the Gaussian predictive model fully decorrelates the measurement process, the score behaves as an innovation process and C = I. |
| `Likelihood_Information_Distortion/section Information Geometry.tex` (69 lines) | Terser alternative for the same three objects. |
| `Gaussian_Fisher_Distortion_Family.md` | Definitions Methods-ready; the file is written as a maintainer note (file lists, CSV column-order warnings) and must be re-prosed, not lifted. Its two load-bearing assumptions (group scale is mandatory; the Gaussian Fisher's minimum eigenvalue is not an indefiniteness flag) both belong in Methods. |
| `implementation_subspace_information_distortion.md` | The numerical-methods paragraph (how singular H and J are handled in code). |

## Verify before submission

- **The C = C_s^(1/2) R C_s^(1/2) identity.** Either prove it under a stated condition, or demote it to the determinant version. Do not print it as written.
- **n_sims uniformity across every cell of every heatmap.** (See above. This is the highest-risk item in the paper.)
- **The direction convention**, once, in code: C_ii > 1 means over-confident. Every verdict depends on it and the repo currently states it with two different signs in three places (`results.md`).
- **Whether `Wald_T2` is computed on the final runs** and whether its debiased form changes any map's ordering. If it does not, say so in one sentence and keep the simpler statistic.
</content>
