# Materials and Methods — what it must contain, with the actual run values

> Working doc, same genre as `abstract.md`. Opened 2026-07-14.
> Every value below was read out of the run scripts and the production data, not out of a plan document. Path given for each. Where a value could not be found, it says so, and those are the items that must be settled before submission.
> eLife excludes Methods from the ~5,000-word main-text target (`../_program/elife-author-instructions.md`), so this section can afford to be complete. It should be.

## The job

Two readers. A referee checking that the numbers could have come out the way we say, and someone in 2029 trying to re-run this. The second one is the harder standard and it is the one to write to.

There is also a specific hazard here that most Methods sections do not have: **this paper's claims are about algorithms, so the Methods *are* the result.** A vague Methods section does not merely inconvenience the reader; it makes the paper unfalsifiable.

## M1 — Model and parameters

Two-state scheme (`scheme_CO`, `legacy/models_simple.h:10-66`): closed C and open O.

- Agonist-dependent association C → O with rate `kon`; agonist-independent dissociation O → C with rate `off`.
- Only the open state conducts (`v_g_formula()[1] = "unitary_current"`); the code applies a sign flip, so the current is **inward**.
- Initial condition: all channels closed, P_initial = [1, 0].
- Noise model: **white instrumental noise only.** Pink noise and proportional noise are both set to zero. This matters and it should be said, because "noise" in this paper means one specific thing (below) and readers will assume 1/f noise is in there.

Six parameters, all estimated in **log10 coordinates** (from `projects/eLife_2025/ops/local/figure_3_mle_G.macroir`, and identical in `figure_1/2/3*`):

| Parameter | Value used to simulate |
|---|---|
| `on` (k_on) | 10 |
| `off` (k_off) | 100 |
| `unitary_current` | 1 |
| `Current_Baseline` | 1 |
| `Current_Noise` | swept (see M3) |
| `Num_ch_mean` | swept (see M4) |

**Units are not written anywhere in the repo.** Not in the model header, not in any `.macroir` script, not in any dispatcher. The agonist concentration `10.0` is called "10 µM" in a single code comment. **This has to be fixed by the author; nobody else can.** Without units the Methods are not reproducible and the parameter values are not interpretable. It is the single most concrete gap this document found.

Note also that `projects/eLife_2025/data/models/scheme_CO_par.csv` (on = 10, off = 100, unitary = 1, noise = 0.001, baseline = 1, Num_ch = 50) is **not** what the figures use, and neither are the model header's defaults. Only the `.macroir` scripts are authoritative. Say so, or delete the misleading files before the repo is carved (`../_program/carve_plan.md`).

## M2 — Protocol

A single concentration jump, three segments: **pre-pulse at 0, agonist at 10, post-pulse at 0** (`figure_2.macroir`, `figure_3_mle{,_G}.macroir`). Sampling at **50 kHz** throughout.

Each segment is specified as (number of measurement steps, number of raw 50 kHz samples averaged per step, agonist concentration). The averaging of raw samples into one measurement step **is** the acquisition window the whole paper is about, and the Methods should say that in one sentence: the interval length is not a filter setting, it is how many raw samples are averaged into each recorded value.

## M3 — The acquisition-interval sweep, and the τ problem

The interval is swept through seven levels (identical injection in every dispatcher):

| `interval_in_tau` | 1 | 0.5 | 0.2 | 0.1 | 0.05 | 0.02 | 0.01 |
|---|---|---|---|---|---|---|---|
| raw samples averaged per step | 500 | 250 | 100 | 50 | 25 | 10 | 5 |
| step duration | 10 ms | 5 ms | 2 ms | 1 ms | 0.5 ms | 0.2 ms | 0.1 ms |
| total steps in the record | 10 | 20 | 50 | 100 | 200 | 500 | 1000 |

**Total record duration is held fixed at 100 ms across the sweep.** That is a good design (it varies the interval without varying the amount of data) and it should be stated as a design choice, not left for the reader to derive.

**The τ inconsistency, which must be resolved before Methods are written.** The label `interval_in_tau` implies the interval is expressed in units of a relaxation time τ. The numbers are consistent with **τ = 10 ms = 1/k_off**. But `../_program/axes.md` defines τ differently, as 1/max_k |Re λ_k| of the generator at the agonist condition, and it also notes "standardize how τ is recorded in run metadata". At the agonist condition the relaxation rate of a two-state channel is k_on·[A] + k_off, not k_off, so **the two definitions do not agree and the axis label of every map depends on which one we mean.** Pick one, state it, and check that the figure axes are labelled accordingly.

## M4 — Channel number and the noise sweep

**N_ch** (`Num_ch_mean`) is swept over {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000} in the Gaussian-Fisher run and {5, 10, 20, 50, 100, 200, 500, 1000, 10000} in the numerical-Fisher run.

**The noise axis is a naming trap and the Methods must defuse it.** The figure labels ("noise 0.1", "noise 1") are **not** the noise parameter. The dispatcher (`ops/slurm/dispatch_figure_3_G.sh`) maps label → `Current_Noise` = label / 1000:

| label | 0.05 | 0.1 | 0.2 | 0.5 | 1 | 10 |
|---|---|---|---|---|---|---|
| `Current_Noise` | 5e-5 | 1e-4 | 2e-4 | 5e-4 | 1e-3 | 1e-2 |

Those six are the levels present in the production data. And what `Current_Noise` *is*, verified at the code level (`project_emission_variance_decomposition`): it enters the emission variance as a term independent of channel number (scaling as N⁰) and inversely proportional to the interval length. It carries no per-channel conductance and no time constant. **It is white instrumental noise expressed in a unit the name does not reveal, and every map in the paper has it on an axis.** One sentence in Methods, and the same sentence in the first figure caption that uses the axis.

## M5 — Simulation (the ground truth)

`simulate(..., simulation_algorithm = "uniformization", number_of_substeps = 1000, seed = 0)`.

The forward model is simulated **exactly** in the sense that matters: uniformization of the continuous-time Markov chain, with 1000 substeps per measurement step, so the interval average of the conductance is realized from the actual trajectory rather than from any moment closure. This is what licenses using it as ground truth, and the Methods should say why 1000 substeps is enough (or show a convergence check; **[Q] does one exist?**).

Figure 1 uses `simulate_with_sub_intervals(..., number_of_substeps = 250)`, which is a different call, for the within-interval visualization. State the difference.

**n_simulations = 10,000 recordings per cell** in both production runs, despite the dispatchers' default of 256. (The defaults in the scripts do not match what was run. See M9.)

## M6 — The five algorithms, as flags

The family is not five codebases; it is one code with two flags (`ops/slurm/dispatch_*.sh`, all five agree):

| Label | `recursive_approximation` | `averaging_approximation` |
|---|---|---|
| `macro_NR` | false | 0 |
| `macro_NMR` | false | 1 |
| `macro_R` | true | 0 |
| `macro_MR` | true | 1 |
| `macro_IR` | true | 2 |

with `taylor_variance_correction = false`, `micro_approximation = false`, `adaptive_approximation = false`, `variance_approximation = true` throughout. **The `averaging_approximation` flag literally counts conditioned endpoints** (0 = instantaneous conductance, 1 = conditioned on the interval's start, 2 = conditioned on both boundary states), which is the endpoint ladder of `../_program/nomenclature.md` and is worth pointing out: the implementation is self-documenting and the family is a genuine grid, not a selection.

**`taylor_variance_correction = false` means the MR that runs is MacroMR, not MacroMRT.** The theory documents in `theory/macroir/docs/Macro_MRT/` describe MRT. The structural sentence transfers (MacroMR lacks the boundary cross-correlation term N_ch·(γᵀΣγ)~ that MacroIR keeps; that is the central difference between the M and I families) but the Taylor variance correction is *not* in the paper and must not be described as if it were.

Useful reduction to quote in Theory or Methods: setting the boundary-conditioned conductance variance to zero recovers the Moffatt 2007 update (`macromrt_macromrt_paper_section.md`). It nails `R` to its published ancestor.

## M7 — Estimation, and the two anchors

Two-stage, per cell (`figure_3_mle_G.macroir`):

1. **The cloud.** `calc_MLE_per_group_of_replicates` with `group_size` = 10 and 100, warm-started at the true parameters, Gauss-Newton with Levenberg damping (`lambda_kickoff = 10`, `lambda_factor = 3`, `lambda_max = 1e10`, `max_iter = 100`, `grad_rtol = 1e-6`, `dvalue_tol = 1e-10`). With 10,000 recordings this gives **1,000 groups of 10** and **100 groups of 100**.
2. **The pooled fit.** The same call with `group_size = n_simulations`, giving a single joint MLE **θ_pool**.

**θ_sim** (the truth) and **θ_pool** (the pooled fit) are different points for a misspecified likelihood, and the gap between them is the bias. The diagnostics battery is run at **both** anchors (`_battery_sim` and `_battery_pool` CSVs), and which one a figure uses changes what it can show: at θ_pool the score is ≈ 0 by construction, so bias washes out; at θ_sim the bias is visible but the Fisher can be indefinite off the maximum. **Every figure must state its anchor.** (`project_gaussian_fisher_distortion_family`.)

**Warm-starting the MLE at the truth** is a defensible choice for a study of the likelihood (it isolates the likelihood's error from the optimizer's) and an indefensible one for a study of an estimator. Say which this is, in one sentence, before a reviewer says it for us.

## M8 — Uncertainty quantification, and a naming correction the paper must not inherit

The intervals in the CSVs come from a **percentile bootstrap over groups**, not from a probit transform, despite every column being named `probit_*` (`project_probit_statistics_is_percentile_bootstrap`; `get_mean_Probits` applies no normal transform, and its output was verified equal to R's `quantile(type = 6)`, difference 0).

**Do not carry the word "probit" into the manuscript.** Write: nonparametric bootstrap over groups, with percentile intervals at 2.5, 16, 50, 84 and 97.5%.

Bootstrap counts, from the scripts:
- diagnostic battery: **100** bootstrap replicates, `max_lag = 10` for the autocorrelation;
- `calc_empirical_distortion_gaussian`: **200** bootstrap replicates over groups;
- `min_groups_for_bootstrap = 10` (below that the interval slots are NaN-filled, which is one source of grey cells in the maps).

Also note: `Model_Parameters_Hat` is not bootstrapped by the program. The raw MLE cloud is what is in `_mle_cloud_runs.csv`, and any bootstrap of it is done downstream in R.

## M9 — The reproducibility problem, and it is real

**There is no run manifest for either production run.** The dispatchers' defaults (`NCHS`, `N_SIMS`, `N_NOISE`, `GROUP_SIZE`) do **not** match what is in the data: the defaults say n_sims = 256 and group sizes {1, 10, 100}; the data says n_sims = 10,000 and group sizes {10, 100}. The actual environment overrides used to launch `1c2ae6f` and `433ed13` were reconstructed from output filenames and CSV columns, and are recorded nowhere in the repo.

Two things follow, and the first is not optional:

1. **Write the manifest now**, from the reconstruction, and commit it next to the data. A Methods section that cannot state the exact invocation is a Methods section that will be challenged.
2. State it plainly in Data Availability: output directories are keyed by the engine's baked git hash (which is why the data folders are `1c2ae6f`, `433ed13`), the engine is pinned by tag, and the figure scripts name their data folder. That part of the pipeline is genuinely good and it should be described.

Environment, for the record: SLURM, 32 CPUs per task, 48 GB (96 GB for the figure-2 dispatch), BLAS threads pinned to 1 with OpenMP taking the simulation and bootstrap loops, and `MACRODR_AXIS_SERIAL=1`, which is load-bearing (without it, memory grows with concurrent combinations and the jobs OOM).

## M10 — The two-anchor split, which Methods has to own

The main text draws its map panels (Fig 2 clouds, Fig 4 maps) from the **Gaussian-Fisher** runs (`1c2ae6f` plus `0ffbda7` for VR); the time-resolved figures (Fig 3 and Fig 3—figure supplement 1) are a separate analytic per-step computation on the `0ffbda7` dumps, not the numerical battery. `00_plan.md` §6 wants the definitive figures anchored on the Gaussian Fisher, and as of 2026-07-22 **they are** (`results.md`; `figures_build_plan.md` §3b).

The earlier open question — justify a numerical/Gaussian split, or run the missing Gaussian cells — is **closed**: Fig 2 was rebuilt on the Gaussian anchor and the two anchors agree at its cell (Gaussian R 1.32 / MR 1.97 / IR 1.02 against the numeric run, `figures_build_plan.md` §3b), so the split is gone for the body. Methods should state the anchor per figure and cite that agreement.

## M11 — Numerical safeguards (short, honest, and no more)

The implementation constrains the occupancy mean to the simplex and the covariance to stay PSD, via a trust-region shrink on the Kalman mean update (`theory/macroir/docs/Macro_IR/macro_ir_variance_inflation_correction.tex`, plus `Binomial_magical_number = 5`, `min_P = 1e-12`, and the two error tolerances in the model header).

One paragraph. It exists because the code runs it, and because a reader implementing the filter will hit the same simplex problem on day one. Do not turn it into a story: the covariance down-date is PSD by construction in the standard branch, the diagonal σ-side trust coefficient turned out to be redundant, and the whole α_σ episode is internal history (`project_psd_trust_redundant`, `project_fisher_instability_trust_argmin`). The paper needs the *what*, not the *how we got here*.

## What to lift, and from where

| Manuscript piece | Source | Readiness |
|---|---|---|
| The filter, in Methods prose | `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md` (50 lines) | **Lift nearly as-is.** It is already written as a Methods block. |
| The full derivation, Supplement | `Macro_IR/macroir_derivation.tex` (1305 lines, compiles) | Supplement-ready. §8 (the measurement update) is labelled a sketch and is the least finished part; finish it or cite the compact version. |
| The compact supplement | `Macro_IR/macroir_macroir_supplement.tex` (140 lines) | Ready. The better choice if the long derivation feels disproportionate for a two-state paper. |
| The diagnostics | `Likelihood_Information_Distortion/supplement_information_distortion_main.tex` | Supplement-ready, with the caveats in `diagnostics.md`. |
| The Gaussian-Fisher anchor | `Gaussian_Fisher_Distortion_Family.md` | Definitions ready; the file is a maintainer note and must be re-prosed. |
| Numerics | `Macro_IR/macro_ir_variance_inflation_correction.tex` | Ready, cut to one paragraph. |

**Do not use** (both carry self-flagged banners): `Macro_IR/Macro_IR_deepseek.md` ("mixes raw-count and per-channel normalizations, contains sign and weighting inconsistencies") and `Macro_IR/MacroIR_tilde.md` ("predates the corrected normalization used by the implementation"). The LLM-merge files in the same folder are provenance, not sources.

## Verify before submission

1. **Units.** All six parameters and the agonist concentration. Author-only. Blocking.
2. **The definition of τ**, and that every axis labelled in τ uses it.
3. **The convergence of the simulator** at 1000 substeps: is there a check, and if not, run one. It is the ground truth; it is worth one supplementary panel.
4. **The run manifest** for `1c2ae6f` and `433ed13`, reconstructed and committed.
5. **Whether `figure_3_time.macroir`'s stale header comments** ("300 simulations", "noise = 1 → Current_Noise = 0.001") have propagated into any caption. The code says 1000 simulations and 1e-4. Trust the code, and fix the comments (`feedback_verify_captions_from_sources`).
</content>
