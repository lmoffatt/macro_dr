# Paper 2 (MacroIR / eLife 2025) — Experiment Grid Spec

Detail layer beneath `00_master_plan_v2.md` §6. Defines the regime sweep for the validity / distortion maps.

## 1) Control variables
From the two-state `scheme_CO`:
- **N_ch** — effective number of channels.
- **K_off** — the kinetic axis (paired with N_ch in the regime maps).
- **interval Δ / τ** — Δ = n_samples / fs; τ = 1 / max_k|Re(λ_k)| of the generator Q at the agonist condition. Standardize how τ is recorded in run metadata.
- **instrumental noise** — swept separately.

## 2) Regime maps
- Primary: **N_ch × K_off** heatmaps (bias, distortion), per algorithm. (Supersedes the old Δ/τ_min × N_ch axes.)
- Additional sweeps: interval and noise, for the MacroIR-focused line plots and the three-domain picture.
- Critical cube from the July audios: N_ch 10–1000, noise 0.05–1, across the interval range — where the distortion is measurable and structured.

## 3) Model
- `scheme_CO` only (two states). `scheme_CCO` (3-state) is out of scope (>2 states = later program component).

## 4) Algorithms
- `NR`, `NMR`, `R`, `MR`, `IR`. Drop the Taylor variance-correction variants (MNRV/MRV/IRV — cut, taylor=false). Naming: MNR → NMR.

## 5) Per-cell protocol
For each cell and algorithm:
1. Fixed non-stationary experiment (single concentration jump / pulse + washout).
2. Simulate n_simulations datasets at θ* (reuse the same samples across intervals/algorithms where possible).
3. Compute the score, the Gaussian Fisher H, and per-interval predictive mean/variance (residuals). Numerical FD-Fisher only to gauge H.
4. MLE / Gauss-Newton local max → empirical parameter covariance.
5. Reduce to the scalar metrics below.

## 6) Scalar metrics per cell (per algorithm)
- score bias norm (max |mean(score_i)/sd(score_i)|) and the projected DIB;
- distortion matrix C = H^(−1/2) J H^(−1/2): det (→ ½ log det C) and its correlation vs sample/geometric split;
- residual mean, var, max ACF (lag>0);
- empirical-vs-sandwich covariance ratio.

Always store the reduced table (so the maps plot cheaply); raw data stays out of git (see `09_carve_plan.md`).

## 7) Exact grid values
The definitive values live in the `ops/` dispatch configs (`figure_2*.macroir`, the `figure_3*` Gaussian-Fisher dispatch) — the source of truth once the Gaussian rerun is final. Do not hardcode a second copy here.
