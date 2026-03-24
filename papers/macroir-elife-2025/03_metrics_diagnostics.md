# Paper 2 (MacroIR / eLife 2025) — Metrics & Diagnostics

This file defines *what we measure* to claim an approximation is “good”.

## 1) Minimal “must-pass” metrics (main text)

### M1 — Score mean (unbiasedness at true parameters)

- **Claim:** At the true parameters θ\*, `E[∇ log L(θ\*)] = 0`.
- **Operational test:** simulate many independent datasets at θ\*, compute the score each time, and show:
  - 0 lies inside the CI for each parameter, and
  - CI width shrinks ~ `1/sqrt(n_simulations)`.
- **Why it matters:** if the score is biased, the MLE/MAP/posterior will be biased.

### M2 — FIM consistency (two estimators)

For a correctly specified likelihood, under regularity conditions:

- `FIM = Cov(score) = -E[Hessian(log L)]`

Operationally, compare:

- **Estimator A:** covariance of per-dataset score (from simulations)
- **Estimator B:** expected negative Hessian (or a Hessian estimator available from the diagnostic output)

Summaries:

- elementwise relative error, matrix norm ratio, determinant ratio, etc.

### M3 — Normalized residual sanity

Define residual per interval:

- `r_t = (y_obs,t - y_pred,t) / sigma_pred,t`

Expected behavior for a good predictive distribution:

- mean(r) ≈ 0
- var(r) ≈ 1
- weak temporal correlation (autocorrelation near 0 beyond lag 0)

This is a simple, robust “first check”.

## 2) Strong diagnostics (recommended)

### D1 — Score correlation diagnostic (hidden dependence)

Compute:

- `Cov( Σ_t score_t )`  (covariance of the total score)
- `Σ_t Cov(score_t)`    (sum of per-interval score covariances)

Their difference reflects cross-time correlation of score contributions.

Interpretation (from audio notes): large discrepancy implies missing “memory” or temporal dependence not captured by the approximation; MacroIR should reduce this relative to MacroMR.

### D2 — Variance inflation / sandwich factor (impact on inference)

Define two information-like matrices:

- `J = -E[Hessian(log L)]`
- `K = Cov(score)`

A standard misspecification diagnostic is the “sandwich”:

- `J^{-1} K J^{-1}`

We can summarize inflation relative to the ideal case (`J = K`) using:

- eigenvalues of `J^{-1}K` (inflation/deflation factors)
- trace / determinant ratios

This provides a direct “reader-facing” interpretation: **how much parameter uncertainty is distorted** by using a wrong likelihood approximation (and implications for evidence).

### D3 — PIT (optional; good for predictive checks)

If the predictive CDF per interval is `F_t`, compute:

- `u_t = F_t(y_obs,t)`

For a correct predictive model, `u_t` should be approximately Uniform(0,1).

## 3) How these connect to existing repo outputs

Existing CSV outputs from `projects/eLife_2025/ops/local/*` include fields like:

- `logL`, `elogL`, `vlogL`
- `y_mean`, `y_var`
- per-interval prior/posterior state probability summaries (`P_mean_*`, `P_Cov_*`)

These are sufficient to build:

- normalized residuals (`y_mean`, `y_var`)
- score summaries (from derivative outputs in `figure_2` datasets)

What may still be needed (plan-only for now):

- explicit extraction of Hessian or a consistent FIM estimator from diagnostics/predictions
- scripts to aggregate metrics into one table per algorithm × regime cell

## 4) Thresholds (to be decided)

We need to decide practical cutoffs to color the validity map.

Suggested start (easy, robust):

- Use **rank-based** (best to worst) per cell for the first iteration.
- Then add absolute thresholds, e.g.:
  - `|mean(r)| < 0.1`
  - `0.8 < var(r) < 1.2`
  - maximum residual autocorrelation at lag>0 below 0.1
  - max |score mean / score sd| below 0.2 (per parameter)

These are placeholders until we see empirical variability.

