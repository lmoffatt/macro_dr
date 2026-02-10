# Paper 2 (MacroIR / eLife 2025) — Validity Map Grid Spec

## 1) Purpose

Define a 2D regime grid that is:

- scientifically meaningful
- easy to compute / reproduce
- aligned with the paper’s core question (when to use each approximation)

## 2) Axes (selected)

### x-axis: Δ / τ_min

- **Δ**: interval duration (seconds). For our experiments it is determined by:
  - sampling frequency (`fs`) and number of samples per interval (`n_samples`)
  - roughly, `Δ = n_samples / fs`
- **τ_min**: fastest characteristic timescale of the kinetics (seconds).
  - Working definition: `τ_min = 1 / max_k |Re(λ_k)|` where `λ_k` are eigenvalues of the generator `Q` at a specified agonist condition (e.g. maximal agonist).

We should standardize how τ_min is computed and recorded in metadata for the grid runs.

### y-axis: N_ch

- effective number of channels in the patch (ensemble size).

## 3) Proposed grid (initial)

This is a starter; adjust after the first pilot run.

### N_ch values

- 10, 30, 100, 300, 1000, 3000

### Δ / τ_min values

- 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2

## 4) Models / schemes

Default selection:

- Main: `scheme_CO`
- Secondary: `scheme_CCO` (robustness/supplement)

If time allows, add one more scheme with more complexity for stress testing.

## 5) Algorithms to compare

Include all presets used in the current scripts:

- `NR`, `R`, `MNR`, `MR`, `MNRV`, `MRV`, `IR`, `IRV`

## 6) Per-cell experimental protocol (high level)

For each grid cell (Δ/τ_min, N_ch), and for each algorithm:

1. Choose a standard experiment structure (pulse + washout) with fixed shape.
2. Simulate `n_simulations` datasets at θ\* (true parameters).
3. Compute:
   - score (and Hessian/FIM estimators if available)
   - per-interval predictive mean/variance (residuals)
4. Reduce to a small set of scalar metrics (Section 7).

## 7) Scalar metrics to store per cell (minimal)

Per algorithm × cell:

- `score_bias_norm` (e.g., max |mean(score_i)/sd(score_i)| across parameters)
- `fim_mismatch` (norm of difference between two FIM estimators)
- `residual_mean`, `residual_var`
- `residual_autocorr_max_lag` (summary)
- `inflation_eigs_summary` (e.g., max eigenvalue of J^{-1}K)

Store raw data as needed, but always also store this reduced table so the validity map is easy to plot.

