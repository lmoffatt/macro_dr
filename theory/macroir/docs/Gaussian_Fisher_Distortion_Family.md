# Gaussian-Fisher-anchored distortion family

Implementation note for maintainers / other threads. Status: **written 2026-07-04/05, not yet compiled or run** (builds are done by the user).

## What and why

Every distortion diagnostic in the codebase is anchored on a **Fisher matrix** `H`: it whitens the empirical score covariance `J` (or the empirical parameter covariance `Cov_emp`) by `H` and measures how far the result is from the identity. Until now `H` was always the **numerical Fisher** `F_b` — a finite-difference Hessian of `logL`, computed by a separate expensive command (`calc_numerical_fisher_information`, ~`2·n_params` extra likelihood passes per recording).

This change adds a **parallel family anchored on the Gaussian Fisher `G_b`** — the cheap analytic moment-matched Fisher `G_b = Σ_t GFI_t`, where per sample `GFI_t = (∂μ)ᵀ(∂μ)/σ² + (∂σ²)ᵀ(∂σ²)/(2σ⁴)`. `G_b` is derivable from the dlikelihood **alone**, so the whole numerical-Fisher stage can be dropped.

Two payoffs:
- **Frame comparison.** For a *macro* algorithm the likelihood is Gaussian by construction, so `G_b` is the model's *own* Fisher and `G_b⁻¹ᐟ²·J·G_b⁻¹ᐟ²` is the classic information-matrix-equality (misspecification) test in the model's frame. Emitting both the F- and G-anchored versions lets the R analysis pick the better estimate.
- **No numerical Fisher.** A run can compute only the G-family and skip `calc_numerical_fisher_information` entirely.

The mode is decided by **whether the numerical Fisher is supplied**, not by a flag.

## Phase 1 — the diagnostic batteries (`likelihood_derivative_*_diagnostics`)

New G-anchored twins in the battery output, computed inside `calculate_Likelihood_diagnostics_preset_f` (`src/core/likelihood.cpp`):

| G tag | twin of | formula |
|---|---|---|
| `Likelihood_Gaussian_Information_Distortion` (gidm) | `Likelihood_Information_Distortion` | `G_b⁻¹ᐟ²·J·G_b⁻¹ᐟ²` |
| `Gaussian_Fisher_Covariance` (gfc) | `Likelihood_Fisher_Covariance` | `G_b⁻¹` |
| `Gaussian_Distortion_Corrected_Covariance` (gdcc) | `Likelihood_Distortion_Corrected_Covariance` | `G_b⁻¹·J·G_b⁻¹` |
| `Gaussian_Distortion_Induced_Bias` (gdib) | `Likelihood_Distortion_Induced_Bias` | `G_b⁻¹·score_mean` |
| `Gaussian_Sample_Distortion` (gsdm) | `Likelihood_Sample_Distortion` | `G_b⁻¹ᐟ²·J_sample·G_b⁻¹ᐟ²` |

Deliberately **not** twinned: `Likelihood_Correlation_Distortion` (CDM) has no Fisher in it (`J_sample⁻¹ᐟ²·J·J_sample⁻¹ᐟ²`), so it is frame-independent and shared; `Likelihood_Gaussian_Fisher_Distortion` (GFD, `G_b⁻¹ᐟ²·F_b·G_b⁻¹ᐟ²`) is the F-vs-G bridge and stays F-gated.

**Mode = presence of the fim argument, resolved by DSL arity.** Each non-paired command name is registered twice:
- 6-arg (`…, numerical_fisher_information, …`) → **both** families + GFD.
- 5-arg (no fim) → **gaussian only**; the F-family + GFD are left empty and drop out of the CSV.

So a config that already omits the fim (e.g. `figure_3.macroir`) now resolves to gaussian-only instead of erroring. The `_paired` command has **no** gaussian overload (pairing asserts equal vector sizes), so gaussian-only always goes through the non-paired command.

The G-family in gaussian-only mode is safe with the empty F-family: `get_mean_Probits` skips empty (`size()==0`) matrices and returns a count-0 empty result, so the writer emits zero rows for the absent F-columns. (Two F-family scalars, `fisher_covariance_effective_rank` / `dcc_effective_rank`, read `0` instead of `NaN` because `effective_rank` is `size_t`-backed; `0` is itself a valid "no Fisher" sentinel — accepted, not fixed.)

## Phase 2 — the figure_3 empirical-distortion capstone

`calc_empirical_distortion` composes `Cov_emp` (empirical θ̂ covariance) with the numerical Fisher into `ECD_Fisher`, `ECD_Corrected`, `Optimum`. New **separate** command (a trimmed clone, not an overload):

```
calc_empirical_distortion_gaussian(cloud, dlik_sim, dlik_bar, n_bootstrap, seed, probit_cis)
```

anchors on the **group-scale Gaussian Fisher** instead. New tags in `legacy/distributions.h`:

| G tag | twin of | formula |
|---|---|---|
| `Empirical_Covariance_Gaussian_Distortion` | `Empirical_Covariance_Fisher_Distortion` | `Ḡ_bar^{1/2}·Cov_emp·Ḡ_bar^{1/2}` |
| `Empirical_Covariance_Gaussian_Corrected_Distortion` | `Empirical_Covariance_Corrected_Distortion` | whiten ECD_G by `GIDM = Ḡ_bar⁻¹ᐟ²·J·Ḡ_bar⁻¹ᐟ²` |
| `Optimum_Gaussian_Distortion` | `Optimum_Fisher_Distortion` | `Ḡ_bar⁻¹ᐟ²·Ḡ_sim·Ḡ_bar⁻¹ᐟ²` (sim frame vs pool frame) |

Key points:
- `dlik_sim` is a **new input** (the F version has none): it supplies `Ḡ_sim`, the Optimum numerator = the θ_sim-vs-θ_pool frame comparison. `dlik_bar` supplies the anchor `Ḡ_bar` and `J`.
- `J` is **unchanged and shared** — it is the score covariance, not a Fisher.
- **Group scale is mandatory.** Per-recording `G_r = Σ_t GFI_t` is fed through the *existing* `sum_fisher` (mean over recordings × `gsize`), exactly like the numerical fim, so `Ḡ` sits at `Cov_emp`'s per-group scale. Using the battery's recording-scale `G_b` (no `×gsize`) would leave `ECD_G ≈ (1/gsize)·I` instead of `≈ I`.
- `Optimum_Gaussian`'s `Min_Eigenvalue` is **not** an indefiniteness flag (the numerical version's was): the Gaussian Fisher is PSD by construction, so it can't go negative. It reads the sim-vs-pool frame gap.
- Per-recording accessor: `Sum<Gaussian_Fisher_Information>(get<Evolution>(dlik[r])(), gfi_lambda)()`.

## Files touched

- `legacy/distributions.h` — 8 new tag types (5 battery + 3 capstone).
- `src/core/likelihood.cpp` — battery G-blocks in `calculate_Likelihood_diagnostics_preset_f`; new `calc_empirical_distortion_gaussian` (+ explicit instantiation for `dMacro_State_Hessian_minimal_param`).
- `include/macrodr/cmd/likelihood.h` — battery typedef additions; `Empirical_Distortion_Gaussian_Analysis`/`_Bootstrap` typedefs; decl; `write_csv` overloads.
- `src/cli/command_manager.cpp` — 5-arg gaussian battery overloads; `calc_empirical_distortion_gaussian` DSL wrapper + registration + `write_csv` registrations.
- `projects/eLife_2025/ops/local/figure_3_mle_G.macroir` — mle_G config: no stage-3 numerical Fisher; batteries via the non-paired gaussian overload; stage 6 via `calc_empirical_distortion_gaussian`; outputs `_empirical_G`, `_battery_*_G`.
- `projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh` — dedicated dispatcher (clone of `dispatch_figure_3.sh`, no `axis_h_fim`/`h_rel_value`, targets `figure_3_mle_G.macroir`, prefix `figure_3_G_`).

## Invariants for future edits

- **Pack order == typedef order.** Every `var::Vector_Space` here is positional. The constructor call in the `.cpp` (`base_pack`, and the capstone return pack) must list fields in the exact order of the typedef in `cmd/likelihood.h`. A mismatch silently mislabels every CSV column, it does not error. Insert G-members twin-adjacent in both, in lockstep.
- **G-family always uses `nan_spd()` (p×p NaN) on per-replicate failure; F-family uses empty 0×0 when there is no fim** (so it drops from the CSV). Do not conflate the two fallbacks.
- **Anchor `Ḡ` must go through `sum_fisher` (×gsize)** in the capstone; the battery's recording-scale `G_b` is the wrong scale there.
- CDM stays shared (no Gaussian twin) — twinning it would double-count the correlation axis and break the SDM∘CDM ≈ IDM reconstruction.
