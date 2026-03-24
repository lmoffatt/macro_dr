# Paper 2 (MacroIR / eLife 2025) — Repro Pipeline Notes

This describes what already exists in the repo, and what we will standardize.

## Current “prototype” pipeline

### MacroIR script entrypoints (data generation)

- `projects/eLife_2025/ops/local/figure_1.macroir`
- `projects/eLife_2025/ops/local/figure_2.macroir`

These scripts:

- load a model (e.g., `scheme_CO`)
- load parameters from CSV
- define an experiment structure
- simulate data
- compute likelihood diagnostics and/or derivative predictions
- write CSV outputs

### Figure generation (R)

- `projects/eLife_2025/figures/figure_1.Rmd`
- `projects/eLife_2025/figures/figure_2.Rmd`

These read from:

- `projects/eLife_2025/figures/data/*.csv`

## What we standardize for the paper

### Output naming and locations

For each figure, define:

- one MacroIR script that writes into `projects/eLife_2025/figures/data/`
- one `.Rmd` that reads only from that directory
- one “meta” JSON per run (even locally) containing:
  - git commit hash
  - grid cell parameters (Δ, τ_min, N_ch, scheme)
  - algorithm choice
  - random seed(s)

### Minimal “rebuild” instructions (to be finalized)

We will add concrete commands once we confirm:

- the built binary path/name in your local build
- any environment assumptions (R packages, etc.)

For now, the principle is:

1. Run the MacroIR scripts to regenerate CSVs.
2. Knit the Rmd to regenerate plots.
3. Commit both plan doc changes and any figure outputs you want tracked.

