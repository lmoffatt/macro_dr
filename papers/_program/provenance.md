# How the figure results were obtained

> Opened 2026-07-14. Updated 2026-07-20. Shared across the three papers; the run manifest they all
> cite. Companion to each paper's Methods (`1_method/methods.md`; the model half migrates to
> `model_and_sim.md` when created), not a replacement for it.
>
> `1_method/methods.md` (M1 to M11) says what the Methods section must *claim*: the model, the protocol, the sweeps, the five algorithms, the estimator, the uncertainty quantification. This document says how the numbers on disk were actually *produced*: which script, which invocation, which data directory, which notebook, which PDF, and what part of that chain can and cannot be reproduced today. It is the run manifest that `1_method/methods.md` M9 asks for, plus the five things M9 does not cover (the R environment, the render step, the seed, data availability, and which cluster ran what).
>
> Every value below was read from the run scripts, the dispatchers, the C++ source or the production CSVs. Path and line are given. Where a code comment and the code disagree, the code wins and the disagreement is flagged. Where something could not be established, it says so, and those items are collected in §9.

## 1. The shape of the pipeline

There are four production runs behind the figures, not one. Each is a single `.macroir` script executed by the `macrodr_cli` binary. Two are launched by hand on a laptop, two are dispatched to a cluster. They write comma-separated value (CSV) files into `projects/eLife_2025/figures/data/`, and a set of R Markdown (`.Rmd`) notebooks under `projects/eLife_2025/figures/paper/` turn those into the figure PDFs.

| Run | Script | Dispatcher | Output | Feeds |
|---|---|---|---|---|
| 1 | `ops/local/figure_1.macroir` | none, run directly | loose `figures/data/figure_1_*.csv` | Figure 1, Fig 1—figure supplement 1 (ex-S1) |
| 2 | `ops/local/figure_3_mle.macroir` | `ops/slurm/dispatch_figure_3.sh` | `figures/data/433ed13/` (487 CSVs) | Figure 2 + its supplements (ex-S2, ex-S3) — numerical anchor, superseded for the body |
| 3 | `ops/local/figure_3_time.macroir` | none, run directly | `figures/data/figure_3_time_dlik_{NR,NMR,R,MR,IR}.csv` (5 files, about 1 GB each) | Figure 3, Fig 3—figure supplement 1 (ex-Fig 4), Fig 3 supplements (ex-S4, ex-S5) |

> **Anchor/seed note (2026-07-22).** The body Fig 2 and Fig 4 were moved to the **Gaussian anchor**
> (`1c2ae6f`/`0ffbda7`), so run 2's `433ed13` above now feeds the *superseded* numerical versions and
> the numerical supplements; and the run-3 dumps were **regenerated at `0ffbda7`, seed 20260722** (the
> `433ed13-dirty` provenance in §6 describes the old-seed dumps). Both are pending a full provenance
> refresh here. Figure renumbering (old Fig 4 → Fig 3—figure supplement 1, old Fig 5 → merged into
> Fig 4, old Fig 6 → Fig 4—figure supplements 3–4, old Fig 8 → Fig 5) follows `1_method/decisions.md`
> "The figure set".
| 4 | `ops/local/figure_3_mle_G.macroir` | `ops/slurm/dispatch_figure_3_G.sh` | `figures/data/1c2ae6f/` (355 CSVs) | exploratory `figure_7_*` notebooks only |

Paths are relative to `projects/eLife_2025/`.

Three naming traps are built into this table and every one of them has already misled a reader of the repository. They must not reach the manuscript.

**The notebook called `figure_2.Rmd` is fed by the script called `figure_3_mle.macroir`.** Every file in `433ed13` is named `figure_3_nch_*`, because that is the script that wrote them, and `figure_2.Rmd:59-60` builds its filenames with exactly that prefix. A Methods paragraph describing Figure 2 must describe run 2. The `figure_2.macroir` and `dispatch_figure_2.sh` lane is a different, older diagnostic battery that never calls the maximum likelihood estimation (MLE) stage at all; its outputs sit loose in `figures/data/` and are read only by archived notebooks.

**The notebooks called `figure_3.Rmd` and `figure_4.Rmd` are not fed by the `figure_3` MLE battery either.** They read the five time-resolved dumps of run 3, which is a separate script with no dispatcher.

**The axis labelled `noise_in_conductance_tau` is not the noise parameter.** Every dispatcher maps the label to `Current_Noise = label / 1000` (`ops/slurm/dispatch_figure_3_G.sh:143`, and the same case block in `dispatch_figure_3.sh:162-168`). The production "noise 0.1" cell is `Current_Noise = 1e-4`. This is stated in `1_method/methods.md` M4 and is repeated here because the label, not the value, is what appears on every figure axis and in every filename.

## 2. The engine, the script language, and how a grid gets swept

A `.macroir` file is a plain-text script in a small typed domain-specific language (DSL), one statement per line, which `macrodr_cli` parses, compiles and runs in a single process. There is no subcommand: `main()` forwards `argv` straight to `main_flow` (`src/cli/main.cpp:3`), which treats every positional argument as a script file and every unrecognised `--` token as an inline DSL line, then concatenates files and inline lines in `argv` order into one program (`src/cli/cli_parser.cpp:140-150`, `src/cli/script_loader.cpp:27-41`). So the invocation is simply:

```
build/gcc-release/macrodr_cli ops/local/figure_3_time.macroir
```

**Parameter injection.** Because inline `--name = value` tokens become DSL statements at their `argv` position, a shell dispatcher can treat a script as a function of its arguments. The dispatchers place every injection *before* the script path, so the injected assignments are the first statements of the assembled program, and they define names that the script deliberately leaves undefined (the injectable names are present in each script but commented out). Assignment is last-writer-wins (`include/macrodr/dsl/grammar_typed.h:115` uses `insert_or_assign`), which is why each dispatcher carries a "FILE CONTRACT" comment in its header listing the names the script must not define itself (for example `ops/slurm/dispatch_figure_3.sh:9-13`). If the script were to assign one of those names, it would silently clobber the injection.

**Axis broadcasting.** The grid is swept inside the process, not by a shell loop. `axis(name, labels)` declares a named axis and `indexed_double_by`, `indexed_size_by`, `indexed_int_by`, `indexed_bool_by` and `indexed_string_by` attach one value per label. Whenever any argument of a command is an indexed value, the whole call is lifted: its index space is the union of the axes of all its arguments, and the underlying C++ function is invoked once per coordinate of the Cartesian product, returning an indexed result (`include/macrodr/dsl/grammar_typed.h:589-631, 790-842`; `legacy/indexed.h:269-319, 540-566`). Vector and tuple constructions lift the same way, which is how one `create_experiment(...)` call fans out over the seven acquisition intervals. One SLURM job is therefore one `.macroir` script covering a whole sub-grid, and the shell layer stays deliberately thin.

The axes do not appear in the output filename. The file stem comes from the injected `filepath` string, and each axis becomes a CSV *column* carrying that coordinate's label (`include/macrodr/cmd/detail/write_csv_common.h:666-673, 750-789`).

**`MACRODR_AXIS_SERIAL=1` is load-bearing.** The axis-combination loop is an OpenMP parallel-for by default. This environment variable flips it to serial (`include/macrodr/dsl/grammar_typed.h:807-811`), so the inner per-simulation loop becomes the active parallel level. Without it, each concurrent combination holds its own per-simulation state and memory grows with the number of combinations, which is how these jobs run out of memory. The SLURM dispatchers set it in the `sbatch --export` list (`dispatch_figure_3.sh:205`); it is not an `argv` token, and the compiled-in default is off.

**Provenance stamping.** The CSV writers used by the analysis commands begin every file with the build's short git commit hash on a line of its own, before the column header (`write_provenance_row`, `include/macrodr/cmd/detail/write_csv_common.h:37-39`). That is why the R notebooks read with `skip = 1`, and it is why the data directories are named `433ed13` and `1c2ae6f`: the dispatchers ask the binary for its own hash (`BIN --commit`) and use it as the output folder name, so two code versions can never write into each other's results (`dispatch_figure_3.sh:65-70`). Checking row 1 of all 970 files across the six commit-named directories gives zero mismatches with the directory name.

The stamp is not universal, and one figure depends on that. The simulation writer (`src/core/simulate.cpp:549`) does *not* stamp, so `figures/data/figure_1_simulation.csv` starts directly with its header, and `figure_1.Rmd:49` correctly reads it with no `skip`. A Methods sentence saying "all CSVs carry the hash and are read with `skip = 1`" would be false for Figure 1.

**The run ledger.** On every invocation the binary snapshots the program it actually compiled: `runs/run-<YYYYMMDD>-<HHMMSS>/script.macroir` holds the assembled text with the injections included, and `meta.json` holds the working directory, a Unix timestamp and the full `argv` (`src/cli/app/workspace_persistence.cpp:22-53`). There are 737 such directories under `projects/eLife_2025/runs/`. This is the strongest provenance in the pipeline and two of the four production runs can be pinned to a specific ledger entry (§4 and §6).

**Build and dispatch.** `ops/build_cluster.sh <cluster> [tag]` sources `ops/clusters/<cluster>.sh` for modules, BLAS, scratch and account, resubmits itself as a batch job if it is not already inside an allocation (login-node builds get killed for memory), and configures into `build/<cluster>-<tag>/` with `-DMACRODR_GIT_COMMIT_OVERRIDE=<tag>`, where `tag` defaults to `git rev-parse --short HEAD`. That override is what makes the build directory name, the binary's `--commit` output and the CSV stamp agree. The SLURM dispatchers submit `ops/slurm/run_macroir.sh` with one job per (algorithm, channel-number) pair, requesting 32 CPUs, 48 GB (96 GB for the figure-2 lane) and a two-day limit, and pinning the BLAS to one thread so OpenMP owns the simulation and bootstrap loops.

## 3. The scientific inputs, in one paragraph

All four runs share a front end, and `1_method/methods.md` M1 to M6 is the authority on it. In brief: the compiled-in two-state scheme `scheme_CO` (`legacy/models_simple.h:10-98`), closed to open at rate `on` times agonist concentration and open to closed at rate `off`, only the open state conducting, all channels closed at t = 0, white instrumental noise only. Six parameters, all fitted in base-10 logarithmic coordinates, simulated at `on` = 10, `off` = 100, `unitary_current` = 1, `Current_Baseline` = 1, with `Current_Noise` and `Num_ch_mean` swept. The scripts call `create_parameters(...)` inline and never load `scheme_CO_par.csv` or `scheme_CO_prior.csv`, so those files are not what the figures use. There is no prior anywhere in the figure pipeline: every figure is likelihood-based, none is Bayesian. Ground truth is generated by exact simulation of the continuous-time Markov chain using uniformization with 1000 sub-steps per measurement interval, which is what licenses treating it as the reference against which the Gaussian approximations are judged.

The five algorithms are one code path under two flags, `recursive_approximation` and `averaging_approximation`, where the averaging flag counts conditioned interval endpoints (0, 1 or 2):

| Label in data | Label in figures | `recursive` | `averaging` |
|---|---|---|---|
| `macro_NR` | NR | false | 0 |
| `macro_NMR` | MNR | false | 1 |
| `macro_R` | R | true | 0 |
| `macro_MR` | MR | true | 1 |
| `macro_IR` | IR | true | 2 |

Note the relabelling: the data key is `NMR` and the displayed label is `MNR`. Both appear in the repository and the notebooks map between them explicitly (`figure_2.Rmd:43`).

## 4. Run 1: Figure 1 and Figure S1

**Producer.** `ops/local/figure_1.macroir`, run with no dispatcher and no injections. The ledger entry is `runs/run-20260711-174017/`.

```
build/gcc-release/macrodr_cli projects/eLife_2025/ops/local/figure_1.macroir
```

**What it computes.** A single recording, no grid, no ensemble, no MLE, no bootstrap. Twenty channels, six measurement intervals at 50 kHz with 100 raw samples averaged into each, giving 2 ms per interval and a 12 ms recording; agonist is 0 during interval 0 and 10 for intervals 1 to 5. One trajectory is simulated with 250 sub-steps per interval. That one recording is then passed through the likelihood algorithms via `calc_likelihood_diagnostic`, writing one long-format CSV per algorithm plus the sub-resolved trajectory (`figure_1_simulation.csv`). Two further IR-only files are written and no notebook reads them. **Since 2026-07-21 the roster is six**: `macro_VR` was added between MR and IR (`figure_1_likelihood_diagnostic_{NR,R,MNR,MR,VR,IR}.csv`). VR is the one build that uses `build_likelihood_function_with_family`, because the bool builder does not expose the variance form; a sibling script `figure_1_plus_lse.macroir` carries the same VR block plus an LSE arm and **writes to these same paths**, so whichever of the two runs last owns the data.

**One inconsistency to be aware of.** This script builds NR and R with `variance_approximation = 0` while MNR, MR, VR and IR get `variance_approximation = 1`. Runs 2, 3 and 4 use `variance_approximation = true` for all five. Figure 1 is illustrative rather than quantitative, so this does not propagate into any reported number, but it should not be described as "the same five builds as the rest of the paper".

**Notebook.** Since 2026-07-21 there are **two** notebooks over one shared body, `figures/paper/figure_1_panels.R`: `figure_1.Rmd` renders the paper's four recursive columns (R, MR, VR, IR) to `Figure_1.pdf` / `Figure_S1.pdf`, and `figure_1_all.Rmd` renders all six to `Figure_1_all.pdf` / `Figure_S1_all.pdf`. They differ by one vector, `COLS`, and cannot overwrite each other. Both read those CSVs (seven since VR, with no `skip`, see §2; the algorithm label is carried by the FILENAME and attached POSITIONALLY from two parallel vectors, since these CSVs have no `algorithm` column), and assemble a four-row panel: prior open probability, predicted current with its spread and the innovation, posterior (NR and MNR, when present, are marked "no update (open loop)"), and cumulative log-likelihood. Each calls `build_figure` twice, producing the zoom from samples 3 and 4 (6 to 10 ms) and the supplement from samples 0 to 4. Both are saved at 7.0 by 7.5 inches.

**Caption correction needed.** `Figure_1_caption.md:5` calls Figure S1 the whole recording. It is not: `sel` filters `sample_index %in% 0:4` (`figure_1.Rmd:508`) and the recording has samples 0 to 5, so sample 5 is dropped and S1 shows 0 to 10 ms of a 12 ms record. Figure S1 also has no caption file of its own.

## 5. Run 2: the numerical-Fisher MLE battery (`433ed13`), feeding Figure 2 and its supplements (ex-S2, ex-S3) — numerical anchor, superseded for the body by the Gaussian rebuild

**Producer.** `ops/local/figure_3_mle.macroir`, dispatched one SLURM job per (algorithm, channel-number) by `ops/slurm/dispatch_figure_3.sh`. A sequential local equivalent with the same injection contract is `ops/local/dispatch_figure_3_local.sh`.

**What one job does.** It simulates `n_simulations` recordings over the injected seven-point acquisition-interval axis, then runs six stages: a per-group Gauss-Newton MLE refit at each `group_size` (an axis, so one job sweeps several group sizes) giving the estimate cloud; a joint MLE over all replicates giving the pooled estimate θ_pool; a central-difference numerical Fisher matrix at both θ_sim and θ_pool; the analytic score at both anchors; the diagnostic battery paired with the Fisher matrix at each anchor (`_battery_sim`, `_battery_pool`); and the empirical-versus-theoretical distortion capstone anchored at θ_pool (`_empirical`).

The numerical Fisher matrix is a central difference *of the analytic score* (the score itself comes from automatic differentiation, so this is one numerical derivative, not two), with a per-coordinate step h_i = h_rel · max(|θ_i|, 1) and h_rel injected as 1e-5 through the `axis_h_fim` axis, symmetrized as (F + Fᵀ)/2 (`src/core/likelihood.cpp:1921-2007`). Non-finite columns are filled with NaN, not zeroed, despite the comment above the function saying otherwise (`likelihood.cpp:1995`).

θ_sim is the simulation truth. θ_pool is the single stationary point of the pooled approximate likelihood. They differ because the Gaussian macroscopic likelihood is misspecified with respect to the exact simulator, and θ_pool minus θ_sim is that misspecification bias. The mean of the per-group cloud minus θ_pool is a separate, finite-sample effect. Which anchor a figure uses changes what it can show, and every figure must state it.

**The grid actually on disk** (99 cells, 487 CSVs), which is *not* the dispatcher's default grid:

- all five algorithms, N_ch in {10, 100, 1000, 10000}, noise labels {0.1, 1, 10}, `n_simulations` = 10000;
- `macro_IR` additionally at N_ch in {20, 50, 200, 500} at all three noise labels;
- `macro_IR` at N_ch = 5, complete only at noise label 10 (at labels 0.1 and 1 only the estimate cloud was written);
- `macro_IR` at `n_simulations` = 200 and 1000, noise label 0.1, as smaller exploratory runs;
- `group_size` = {10, 100} for the 10000- and 1000-recording runs, and {1, 10, 100} for the 200-recording runs.

Five CSV families are written per cell: `_mle_cloud_runs`, `_pool_runs`, `_battery_sim`, `_battery_pool`, `_empirical`.

**Reconstructed invocation.** The dispatcher's committed defaults (`N_SIMS=256`, `N_ALGO=macro_IR`, `GROUP_SIZE=(1 10 100)`) do not match the data, so the production run used environment overrides that were never recorded. The invocation consistent with what is on disk is:

```
NCHS="10 100 1000 10000" N_SIMS="10000 10000 10000 10000" \
N_NOISE="0.1 0.1 0.1 0.1" N_ALGO="macro_NR macro_NMR macro_R macro_MR macro_IR" \
GROUP_SIZE="10 100" \
projects/eLife_2025/ops/slurm/dispatch_figure_3.sh <cluster>
```

repeated for noise labels 1 and 10, and again with `N_ALGO="macro_IR"` and the extended channel list. Note that the dispatcher pairs `N_NOISE` with `NCHS` *by index* inside the same loop, so it does not form a cross-product over noise: each noise level is a separate invocation.

**Notebooks.** `figure_2.Rmd` slices exactly one cell of this grid: N_ch = 100, 10000 recordings, noise label 0.1, group size 10, acquisition interval 0.1 τ (`figure_2.Rmd:35-40`). It plots the cloud of per-group MLEs for each algorithm together with three ellipses (the empirical covariance of the cloud, the Fisher prediction F_b⁻¹, and the sandwich correction F_b⁻¹ J F_b⁻¹) and a predicted-bias marker F_b⁻¹ times the mean score, all rescaled by 1/group_size = 1/10. It applies no convergence filter to the cloud. Saved at 7.0 by 9.5 inches. `figure_S2.Rmd` and `figure_S3.Rmd` read the same directory.

`figure_5_master*.Rmd` and `figure_6_precision.Rmd` also lived in `paper/` and read `433ed13`, but they were **never** the paper's figures; both are **archived** under `figures/archive/paper_superseded_20260722/` (2026-07-22). The first was self-described as a correctness audit (`figure_5_master_STATUS.md:14`) and the second carried the header `*** DATA SOURCE IS A DRAFT ***`. The body's design figure is the new **Fig 5** (ex-8). See §8.

## 6. Run 3: the time-resolved dumps, feeding Figure 3, Fig 3—figure supplement 1 (ex-Fig 4) and the Fig 3 supplements (ex-S4, ex-S5)

**Producer.** `ops/local/figure_3_time.macroir`, run directly with no dispatcher, so every value is script-defined and nothing is injected. The ledger entry is `runs/run-20260630-202545/`, whose `meta.json` records the exact invocation:

```
cd projects/eLife_2025
../../build/gcc-release/macrodr_cli ops/local/figure_3_time.macroir
```

**What it computes.** One shared ensemble of 1000 recordings of a 100-channel patch at `Current_Noise` = 1e-4, with 100 measurement intervals of 1 ms each (50 raw samples at 50 kHz), agonist applied during intervals 20 to 59, for a 100 ms recording. That single ensemble is then run through all five likelihood builds and dumped by `calc_dlikelihood_predictions` into five long-format CSVs of roughly 1 GB each.

**Header drift, already checked.** The script's own header comment (`figure_3_time.macroir:7`) says "150-sample protocol" and a noise of 0.001. The code says 100 intervals (20 + 40 + 40, line 33) and `Current_Noise` = 1e-4 (line 30), and the CSVs agree with the code. Neither stale number has reached a caption, and the comment should be fixed.

**Three traps in the dump, all confirmed against the writer and the data.**

1. *Every evolution row is written twice.* `emit_state_rows_with_experiment` calls `emit_state_rows_without_experiment` as its first act (`src/core/likelihood.cpp:2421`), so each per-interval record appears once with a blank segment index and once annotated. The notebooks deduplicate by averaging over (`simulation_index`, `sample_index`); a naive sum double-counts everything.
2. *There is no per-interval Fisher column.* Per-interval Fisher information is rebuilt in R from the dumped sensitivities as I_t = (∂y_mean/∂θ)²/y_var + ½ (∂y_var/∂θ)²/y_var². The 36 entries of the whole-recording Gaussian Fisher matrix are present, but only in the `state`-scope rows.
3. *Derivatives are zero before the agonist step*, which makes the first 20 intervals look empty rather than broken.

The per-interval log-likelihood and score are incremental, not cumulative. This was verified numerically: for the first IR recording the 100 deduplicated interval values sum to -147.3591333890, against a whole-recording total of -147.3591333890.

**Notebooks.** `figure_3.Rmd` renders a six-row time-resolved calibration cascade (predicted output and log-likelihood, standardized squared residual, score bias, per-interval score variance against per-interval Fisher, the same accumulated, and score autocorrelation), with 400 bootstrap replicates over the 1000 recordings, saved at 7.0 by 11.0 inches. `figure_4.Rmd` forms the per-step pair F_t (mean Fisher over recordings) and J_t (variance of the score over recordings) and plots log₁₀(J_t/F_t) with 200 bootstrap replicates and a 90% percentile interval. `figure_S4_S5.Rmd` reads the same five files and emits both `Figure_S4.pdf` (bias) and `Figure_S5.pdf` (autocorrelation). Reading these files takes about 5 GB of memory, which the notebook header warns about.

**Provenance caveat.** Row 1 of all five dumps reads `433ed13-dirty`. The trailing `dirty` means the binary was built from a working tree with uncommitted changes, so the hash does not uniquely identify the source that produced Figures 3 and 4. The ledger snapshot in `runs/run-20260630-202545/script.macroir` pins the *script*, but not the engine.

## 7. Run 4: the Gaussian-Fisher family (`1c2ae6f`)

**Producer.** `ops/local/figure_3_mle_G.macroir` via `ops/slurm/dispatch_figure_3_G.sh`. This is run 2 with the finite-difference Fisher stage deleted and every distortion diagnostic re-anchored on the analytic Gaussian Fisher matrix G_b = Σ_t GFI_t, with GFI_t = (∂μ_t/∂θ)(∂μ_t/∂θ)ᵀ/σ²_t + (∂σ²_t/∂θ)(∂σ²_t/∂θ)ᵀ/(2σ⁴_t) (`legacy/qmodel.h:6096-6098`). Because G_b follows from the analytic derivative alone, the numerical Fisher routine, which costs 2·n_params extra likelihood passes per recording, is never called, and there is no `axis_h_fim` column in the output.

The mode is selected by call *arity*, not by a flag: the five-argument form of `likelihood_derivative_basic_diagnostics` resolves to the Gaussian overload, which passes an empty Fisher vector (`src/cli/command_manager.cpp:127-133`). This is worth knowing before anyone tries to switch modes by editing a boolean.

**Coverage, and the gap that matters.**

| Algorithm | N_ch covered | Noise labels |
|---|---|---|
| `macro_IR` | 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000 | 0.05, 0.1, 0.2, 0.5, 1, 10 |
| `macro_R`, `macro_MR` | 10, 100, 1000, 10000 | 0.1 only |
| `macro_NR` | 10, 100, 1000 | 0.1 only |
| `macro_NMR` | none | none |

All at 10000 recordings. **MacroNMR is entirely absent from the Gaussian family, and MacroNR is missing its 10000-channel cell.** This is the single fact that governs what the definitive figures can currently be anchored on, and it is the operational content of the anchor split described in `1_method/methods.md` M10: any cross-algorithm statement has to come from `433ed13` (numerical Fisher, all five algorithms), and only the IR mechanism panels can come from `1c2ae6f`.

**Consumers.** Only the `figure_7_*` notebooks read this directory, four of them anchored on θ_sim (`battery_sim_G`) and three on θ_pool (`battery_pool_G`). None of them is a numbered figure in the manuscript. A separate directory, `87889e6`, holds `micro_R` and `micro_IR` runs of the same pipeline and feeds no paper figure.

## 8. Which figures exist, and which do not

The canonical arc is `1_method/results.md` (updated 2026-07-22), and the figure **set** (numbering,
body-vs-supplement) is owned by `1_method/decisions.md` "The figure set"; both supersede the figure list
in `1_method/00_plan.md` §5.

| Figure (new) | Notebook | PDF | Data | Status |
|---|---|---|---|---|
| 1 | `figure_1.Rmd` | `Figure_1.pdf` | run 1 | finished, captioned |
| 2 | `figure_2.Rmd` | `Figure_2.pdf` | Gaussian anchor (`1c2ae6f` + `0ffbda7`) | rebuilt 2026-07-22 on the Gaussian anchor, captioned |
| 3 | `figure_3.Rmd` | `Figure_3.pdf` | run 3, regenerated at `0ffbda7` seed 20260722 | finished, captioned |
| 3—fig. supp. 1 (ex-Fig 4) | `figure_3_supplement_1.Rmd` | `Figure_3_supplement_1.pdf` | run 3 | built (per-step Fisher / Fisher-to-zero) |
| 4 (R vs IR) | `figure_4.Rmd` | `Figure_4.pdf`, `Figure_4_short_variant.pdf` | Gaussian anchor | built 2026-07-22; merges the old IR-only map |
| 4—fig. supp. 1–2 | `figure_4_supplement_{1,2}.Rmd` | `Figure_4_supplement_{1,2}.pdf` | Gaussian anchor | built (bias, distortion, all five params, R vs IR) |
| 4—fig. supp. 3 (ex-Fig 6) | `figure_4_supplement_3.Rmd` | `Figure_4_supplement_3.pdf` | Gaussian anchor | built (sample/correlation decomposition, k_off) |
| 4—fig. supp. 4 | `figure_4_supplement_4.Rmd` | `Figure_4_supplement_4.pdf` | Gaussian anchor | built (sample/correlation vs N_ch) |
| 5 (design trade-off, ex-8) | none selected | none | intended `1c2ae6f`, corrected covariance | **not built.** Candidates in `figures/in_progress/figure_5_IR_*` |
| 1—fig. supp. 1 (ex-S1) | `figure_1.Rmd` | `Figure_S1.pdf` | run 1 | built, no caption file, caption text wrong (§4) |
| 2—fig. supp. 1 (ex-S2) | `figure_S2.Rmd` | `Figure_S2.pdf` | `433ed13` | finished |
| 2—fig. supp. 2 (ex-S3) | `figure_S3.Rmd` | `Figure_S3.pdf` | `433ed13` | finished |
| 3—fig. supp. 2–3 (ex-S4, S5) | `figure_S4_S5.Rmd` | `Figure_S4.pdf`, `Figure_S5.pdf` | run 3 | finished; `Figure_S4_acf_caption.md` describes the file saved as `Figure_S5.pdf` |

The manuscript itself is not yet assembled: `docs/manuscript-drafts/elife_paper.tex` loads `graphicx` but contains zero `\includegraphics` and zero figure environments. The figures are named only in `%` comments, two of which (Figure 2 and Figure 5) end with "[regenerate on Gaussian rerun]". Nothing has been staged into `papers/1_method/figures/`, which holds only `instructions.md`.

## 9. The R layer, and what is not reproducible

**Packages.** The notebooks load `tidyverse`, `patchwork`, `data.table` and `scales`, and call into `MASS`, `dplyr`, `tibble`, `ggplot2` and `knitr`. There is no `renv.lock`, no `.Rprofile`, no `sessionInfo()` dump and no R project file anywhere in the repository, so no R version is pinned to the rendered PDFs. The PDFs were rendered between 2 and 11 July 2026; the machine currently carries R 4.6.1, ggplot2 4.0.3, patchwork 1.3.2, data.table 1.18.4 and rmarkdown 2.30, but nothing ties those versions to the files.

**There is no render step in version control.** No Makefile, no `Rscript` call, no `_output.yml`, nothing in `ops/`. The notebooks are knit by hand. The reconstructed command, from `figures/paper/`, is `Rscript -e 'rmarkdown::render("figure_3.Rmd")'`.

**Random seeds, and a correction to `1_method/methods.md`.** M5 currently records `seed = 0` as though it fixed the random stream. It does the opposite. `calc_seed(0)` draws from `std::random_device` (`legacy/mcmc.h:37-45`), so zero is the sentinel meaning "seed randomly". Every production script passes `seed = 0` to `simulate(...)`, and the resolved seed is written neither to `meta.json` nor to the CSVs. **The simulated ensembles behind every figure are therefore not reproducible.** A rerun of the same script with the same code gives a statistically equivalent, numerically different ensemble. The same literal zero passed to the diagnostic, MLE and bootstrap commands goes straight into `mt_64i(0)` and *is* deterministic, so the resampling is reproducible given the data. On the R side, `set.seed(1)` appears only in `figure_3.Rmd:128` and `figure_4.Rmd:80`; the bootstraps in the other notebooks are unseeded.

This is fixable going forward (log the resolved seed in `meta.json`), and it is worth doing before any rerun, but it cannot be fixed retroactively for `433ed13` or `1c2ae6f`.

**The data are not in version control.** `.gitignore:52` is `*.csv`. None of the 487 files in `433ed13`, the 355 in `1c2ae6f`, the five 1 GB dumps or the loose Figure 1 inputs is tracked. The figure PDFs are tracked; their inputs are not, and no repository, archive or digital object identifier is named anywhere for them. Whatever is said in Data Availability has to be made true first.

**The environment is under-recorded.** The cluster profiles do pin modules (`clusters/dirac.sh` uses gnu14 with Intel MKL 2019.5.281 and CMake 3.30.8; `clusters/clementina.sh` uses GCC 15.1.0 with OpenBLAS 0.3.29), but **which cluster produced `433ed13` and which produced `1c2ae6f` is written nowhere**. The local builds behind runs 1 and 3 used an unrecorded compiler; `CMakePresets.json:21` pins the path `/usr/bin/g++` with no version, and the standard is C++20. The copy from cluster scratch back into the repository is not scripted anywhere either.

**Cost is unrecorded.** The SLURM requests are known (32 CPUs, 48 GB or 96 GB, two-day limit) but no job output was kept, so actual wall-clock and CPU consumption cannot be reported.

## 10. What to fix before submission

Ordered by how much a referee would care.

1. **Correct the seed claim** in `1_method/methods.md` M5, and decide what the paper says about it. The honest sentence is that the ensembles are large enough (1000 to 10000 recordings) that the reported statistics are stable, but the exact datasets cannot be regenerated bit for bit. Add the resolved seed to `meta.json` before any rerun.
2. **Deposit the data.** Two directories and five large dumps. Until they are archived and cited, the provenance story in §2, which is genuinely good, is a story about files on one laptop.
3. **Commit the run manifest.** §5 and §7 reconstruct the invocations from filenames and CSV columns. Put them next to the data, since the dispatcher defaults actively contradict what was run.
4. ~~**Decide the anchor split**~~ **Resolved for the body (2026-07-22):** Fig 2 and Fig 4 are on the Gaussian anchor and the two anchors agree at the Fig-2 cell (`1_method/figures_build_plan.md` §3b); NMR is out of paper 1, so the MacroNMR gap no longer forces a split. State the anchor per figure and cite the agreement.
5. ~~**Build Figures 5 and 6**~~ **Renumbered/done (2026-07-22):** the R-vs-IR map is body **Fig 4** (built) and the decomposition is **Fig 4—figure supplements 3 and 4** (built); the design trade-off is **Fig 5** (ex-8), still `figures/in_progress/` scripts pending promotion. See `1_method/decisions.md`.
6. **Pin the R environment** (an `renv.lock` and a one-line render script would close the last unscripted hop in the chain).
7. **Fix the stale comments and captions** identified here: the Figure S1 caption, the missing Figure S1 caption file, the S4/S5 caption filename mismatch, and the header comments in `figure_3_time.macroir` and `figure_3_mle.macroir`.
