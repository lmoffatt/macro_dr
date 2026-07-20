> **RETIRED 2026-07-20. Superseded by `../../_program/machinery.md`.**
> Its operational definitions were merged there verbatim; the CSV-field pointer went to _program/provenance.md.
> Do not read as current and do not edit. Kept for its content; the live owner is the successor.

---

# Paper 2 (MacroIR / eLife 2025) — Metrics & Diagnostics

Detail layer beneath `00_master_plan_v2.md` §3 (the authoritative summary); this file holds the operational definitions. Anchor throughout: the model's own **Gaussian Fisher** H; the numerical finite-difference Fisher enters only to gauge how faithful the Gaussian one is. All tests evaluated at the global MLE / optimum (for a misspecified likelihood, θ_sim and the optimum differ, and that gap is the bias).

## 1) Must-pass metrics (main text)

### M1 — Score bias
- **Claim:** at the true parameters θ*, E[∇ log L(θ*)] = 0.
- **Test:** simulate many datasets at θ*, compute the score each time; 0 should lie in the CI per parameter, CI width shrinking ~ 1/√(n_simulations). Project the residual E[score] through the parameter covariance → a bias vector (DIB), the reader-facing "how wrong is the estimate".
- **Why:** a biased score biases the MLE.

### M2 — Fisher consistency (Gaussian Fisher vs score covariance)
- Correct specification ⇒ Var(score) = H, the Gaussian Fisher.
- **Anchor:** H is the model's analytic Gaussian Fisher (PSD by construction), NOT a numerical −E[Hessian]. The numerical FD-Fisher is computed only to measure how faithful H is (the F-vs-G bridge).
- Compare Cov(per-dataset score) against H; summarize via the distortion matrix (§2).

### M3 — Normalized residual sanity
- r_t = (y_obs,t − y_pred,t) / σ_pred,t.
- Good predictive: mean(r) ≈ 0, var(r) ≈ 1, and **no temporal autocorrelation** (ACF ≈ 0 beyond lag 0). Whiteness is not optional: residual autocorrelation is the correlation-distortion signal.

## 2) Distortion matrix (the core quantitative object)

Symmetric sandwich, anchored on the Gaussian Fisher:
- **C = H^(−1/2) J H^(−1/2)**, with J = Cov(score) (J_total, incl. cross-interval terms).
- C = I under correct specification; deviation quantifies how much the approximation distorts parameter uncertainty.
- **Decomposition:**
  - **correlation distortion** — temporal, missing higher moments appearing as ghost state-correlation (short Δ, few channels);
  - **sample / geometric distortion** — per-sample non-Gaussianity (J_sample-anchored).
- **Direct check:** empirical covariance of the MLE cloud vs the sandwich-predicted covariance.

### Evidence-correction payoff (motivation only)
Near the maximum the distortion propagates to the Bayesian evidence through two scalar summaries of eig(C):
- volume: Δ log Z = ½ log det C (geometric mean, O(1) in T, decisive for Bayes factors);
- peak / effective samples: α⋆ = p / tr C (arithmetic mean, O(T)).
Derivation deferred to the bridge-3 study.

## 3) Optional predictive check
- **PIT:** u_t = F_t(y_obs,t) ≈ Uniform(0,1) for a correct predictive. Good extra check; not one of the four main indicators.

## 4) Thresholds (open — D-4 in v2 §8)
- Candidate "valid": distortion < 1.1 per parameter (≈10% error); bootstrap error < 1.1.
- Residual placeholders until empirical variability is seen: |mean(r)| < 0.1; 0.8 < var(r) < 1.2; max ACF (lag>0) < 0.1; max |mean(score)/sd(score)| < 0.2 per parameter.
- Start rank-based (best→worst per cell), add absolute cutoffs after the Gaussian rerun.

## 5) Repo outputs feeding these
- Diagnostic CSVs expose logL / elogL / vlogL, y_mean / y_var, per-interval prior/posterior state summaries, and the score + Gaussian-Fisher fields — enough for residuals, score summaries, and the distortion matrix. Exact field names: the current `figures/paper/*.Rmd` are the source of truth once the Gaussian rerun is final.
