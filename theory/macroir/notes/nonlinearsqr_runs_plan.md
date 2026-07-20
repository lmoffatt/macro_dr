# nonlinearsqr (LSE) — runs plan (data generation for the figures)

Companion to `nonlinearsqr_lse_plan.md` (implementation) + `nonlinearsqr_cpp_spec.md`. This is the OPS plan: which runs produce the data the LSE figures need. Telegraphic.

## Goal / completion criteria
LSE (`nonlinearsqr`) present in the figure-3 validity/bias story (θ̂ cloud + pool) and, via the same `_mle_cloud`, the fig-2 Fisher-ellipse-vs-cloud story, across the SAME grid as the macro production so LSE is comparable. DONE = `figure_3_LSE_*_mle_cloud_runs.csv` + `_pool_runs.csv` exist for the full grid at a NEW committed data hash, and the R figures render the `nonlinearsqr` series.

**SCOPE REVERSED 2026-07-19:** LSE is a FULL member — fig 2 (cloud/ellipse), fig 3 (bias/validity), **fig 4 (cumulative J_T/F_T)** and **fig 5 (distortion)**. The earlier "fig 4/5 excluded, no per-step GFI" was wrong (see `nonlinearsqr_cpp_spec.md` §H: per-interval y_var ≡ σ̂²=SSE/n with zero derivative makes every existing formula produce the LSE quantity; ONE guard flips at likelihood.cpp:928, nothing downstream changes). So the runs now cover the FULL figure_3_mle_G pipeline (Stages 1-6), not just the MLE:
- Stages 1-2 → `_mle_cloud` + `_pool` (fig 3 bias, fig 2 cloud; the per-θ̂ `Gaussian_Fisher_Information` also gives fig 2's ellipse).
- Stage 4 `calc_dlikelihood_predictions` → the per-interval Evolution → **`_time_dlik` (fig 4)**; also feeds Stages 5-6.
- Stages 5-6 gaussian battery + `calc_empirical_distortion_gaussian` → `_battery_*_G` + `_empirical_G` (**fig 5**). NOTE the provenance seam: LSE can only produce the **gaussian `_battery_*_G`** variety (the `_paired` overload needs a numerical Fisher); figure_5_master.Rmd currently reads the numeric-Fisher `_battery_*` — reconcile before rendering fig 5.
- fig 4 additionally needs a `figure_3_time.macroir`-style run (its own `_time_dlik` block) switched to `build_likelihood_function_with_family(family=2)`.

## Prereqs (state)
- C++ compiles; driver validated (θ̂ recovers truth on the smoke). [done]
- naming-under-Fixed fix (`parameter_indexed.h`) + 2 warnings. [done, needs the recompile = dirac job 117190]
- Dedicated scripts: `ops/local/figure_3_mle_LSE.macroir` (Stages 1-2, Fixed unitary_current+Current_Noise, `_with_family` family=2), `ops/slurm/dispatch_figure_3_LSE.sh` (case-arm `nonlinearsqr`, family injection), `ops/local/run_figure_3_mle_LSE_local.sh` (one-point local). [done]
- PENDING before any grid: (1) 117190 build finishes → new dirac binary; (2) local validation (`run_figure_3_mle_LSE_local.sh`) confirms param_name aligned + θ̂ near truth.

## The grid (matched to the macro production, inferred from figures/data/*)
Two pieces, same as macro:
1. **Main grid** — the N_ch axis at the canonical noise:
   - N_ch ∈ {10, 100, 1000, 10000}, noise = 0.1, group_size ∈ {1, 10, 100}, 7 intervals (swept in-job via axis_interval [1,0.5,0.2,0.1,0.05,0.02,0.01]), nsim = 1024, algo = `nonlinearsqr`.
2. **Noise sweep** — noise-dependence, concentrated at N_ch=10 @ high nsim:
   - N_ch = 10, noise ∈ {0.05, 0.1, 0.2, 0.5, 1, 10}, group_size ∈ {1, 10, 100}, 7 intervals, nsim = 10000, algo = `nonlinearsqr`.

Optional 2nd variant (DEFERRED, decide later): av=0 (instantaneous mean) as a second algo label `nonlinearsqr_av0` — doubles the runs. Fix the av=0 driver bug first ([[project_nonlinearsqr_lse_algorithm]]). MVP = av=1 only.

## Dispatch commands (dirac; after the build + local validation)
Main grid (one job per (algo, N_ch); group_size + intervals sweep inside each job):
```
NCHS="10 100 1000 10000" N_SIMS="1024 1024 1024 1024" N_NOISE="0.1 0.1 0.1 0.1" \
  N_ALGO="nonlinearsqr" GROUP_SIZE="1 10 100" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_LSE.sh dirac
```
Noise sweep (N_ch=10, high nsim):
```
NCHS="10 10 10 10 10 10" N_SIMS="10000 10000 10000 10000 10000 10000" \
  N_NOISE="0.05 0.1 0.2 0.5 1 10" N_ALGO="nonlinearsqr" GROUP_SIZE="1 10 100" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_LSE.sh dirac
```
Chain both after the build: `DEPEND=117190 … dispatch_figure_3_LSE.sh dirac` (afterok on the build job).

## CPU cost (small vs macro)
Per job, LSE is MUCH cheaper than macro: mean-only (NO P_Cov propagation, NO Kalman down-date, NO per-step Gaussian Fisher), and Stages 4-6 are dropped. calc_Qdt (shared) is the only heavy kinetics and is reused. And it is ONE algo vs 5. So the LSE grid ≈ a small fraction of one macro-algo's grid. Expect well within the remaining dirac budget; confirm from the first job's wall-time × job count before launching the full sweep. (10 jobs main + 6 noise = 16 jobs, each sweeping 7 intervals × 3 group_size.)

## Sequence
1. 117190 build finishes → new dirac binary (`-current` repointed). [in flight]
2. Local: `run_figure_3_mle_LSE_local.sh` → confirm param_name + θ̂ (Phase 2). Pegar el `_mle_cloud_runs.csv`.
3. Dispatch main grid (16 jobs total with the noise sweep). Harvest to `$WORKDIR/figures/data/figure_3_LSE_*`.
4. Copy the CSVs into the repo under a NEW `figures/data/<hash>/`, commit → the frozen hash.
5. R (Phase 4): add `nonlinearsqr` to ALGOS + ALGO_LAB + filters in the fig2/fig3 target Rmds, repoint DATA_DIR to the new hash, render. Confirm the `nonlinearsqr` series appears (guard against the silent NA-drop: the CSV `algorithm` cell = `nonlinearsqr` must match the R ALGOS key exactly).

## Open decisions (gate the FULL grid)
- **paper vs internal**: full grid above (paper) vs a few points (internal). If paper: also the manuscript reframe ("five → six approximations"), [[project_nonlinearsqr_lse_algorithm]].
- **av=0**: include as a 2nd variant (doubles runs, needs the av=0 fix) or av=1 only.
- **nsim**: 1024 main / 10000 noise-sweep matches macro; drop if a coarser cloud suffices.

## Completion checklist
- [ ] build 117190 done, new binary current
- [ ] local validation green (param_name + θ̂)
- [ ] main grid dispatched + finished (N_ch {10,100,1000,10000} @ noise 0.1)
- [ ] noise sweep dispatched + finished (N_ch=10 @ 6 noise levels)
- [ ] CSVs copied to figures/data/<newhash>/ + committed
- [ ] R figures fig2/fig3 render the nonlinearsqr series (no NA-drop)
