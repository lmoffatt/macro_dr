# Data note: what a reader of the deposit needs and the paper does not carry

Opened 2026-08-26, executing the `repository documentation` MOVEs of
`papers/1_method/docs/shortening_plan.md`. Everything here was **moved out of the manuscript, not
deleted**: it documents the deposited scripts, files and directories for someone reading them, which
is what a data note is for and what Methods is not.

Companion destinations: statistical specification that left Methods went to
`papers/1_method/docs/manuscript-drafts/supplementary_file_1.tex` (Supplementary File 1), not here.

## Parameter values: read the run scripts, not the tables

The parameter values behind every reported figure are read from the production run scripts
(`projects/eLife_2025/ops/local/*.macroir`), **not** from the parameter tables carried in the
repository: those tables hold stale values the figures never used.

## `number_of_substeps` is inert on this branch

The production scripts do pass a `number_of_substeps` argument, and on this branch it has no effect:
the simulation parameters are built from the algorithm name alone, and the substep count they carry
is zero.

- `src/core/simulate.cpp:41-52, 325-345` (`simulation_parameters_from_algorithm` takes the name only)
- `legacy/qmodel_types.h:1624-1627` (uniformization sets `Simulation_n_sub_dt(0)`)
- `legacy/qmodel.h:8079-8095` (the uniformization branch calls `uniformization_sample(mt, ..., fs)`
  with no substep count)

The one place a substep count is real is Figure 1's single 20-channel illustrative trace, which is
simulated at 250 substeps per interval; that run is described in Methods and does not enter any
measurement.

## The three data directories the figures read

The reported results are read from three directories, each named by the git commit hash of the
build that wrote it: `1c2ae6f`, `87889e6` and `0ffbda7`, **searched in that order with the first
hit winning**. Every cell entering a body figure is at `n_sim = 1e4`, which the figure code
enforces by matching `nsim_10000` in the filename. The later differenced-Fisher battery is in
`a202e03` (see Supplementary File 1, S1.5).

The figure code auto-detects which cells exist on disk at render time and prints the detected grid
into the render log, so the cell set behind each figure is recoverable from that log.

## Parameter indices in the least-squares output files

The four-parameter configuration is `ops/local/figure_4_LSE.macroir`, which produces the grid
cells, together with `ops/local/figure_3_mle_LSE.macroir`, the same model run with the per-group
estimate cloud retained for the few cells that need one. Its parameter indices are
`0 = k_on`, `1 = k_off`, `2 = current baseline`, `3 = mean channel number`. **There is no
unitary-current index in these files, and index 2 is the baseline.**

The six-parameter configuration is `ops/local/figure_3_time.macroir`. Keeping six parameters free
there preserves the parameter-index alignment with the macroscopic dumps that the figure code joins
against; fixing two would emit four indices against the macroscopic files' six and misalign them
without any visible error.

## Bootstrap internals not stated in the paper

Despite the internal column names, no probit transform is applied: the reported intervals are
ordinary empirical percentiles, verified equal to the R quantile of type 6. The program does run a
group bootstrap of the estimate cloud, with a floor of 10 groups below which it falls back to a
single unresampled pass, but that result is held in memory and never written; the cloud's only
output file is the per-group run log, one row per group and component.

## Run-ledger provenance

On every invocation the binary snapshots the assembled script and its full argument list into a run
ledger, so provenance is documented per file across multiple commits, beyond the commit-hash stamp
that every output file carries in its first row.

## Cross-language verification

The R and Python packages of `macroir` are bindings over one C++ core; they were run against it on
a single cell by `tools/cross_language_check.py`, which is the evidence for the statement in Data
availability that the three interfaces return the same numbers.
