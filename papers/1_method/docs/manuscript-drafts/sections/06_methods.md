# Methods — the brief for `06_methods.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `06_methods.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Two readers: a referee checking the numbers could have come out this way, and someone in
2029 trying to re-run it. Write to the second. **This paper's claims are about algorithms, so the
Methods *are* the result**; a vague Methods section does not inconvenience the reader, it makes the
paper unfalsifiable.

**A naming correction the manuscript must not inherit.** The intervals in the CSVs come from a
**percentile bootstrap over groups**, not from a probit transform, despite every column being named
`probit_*`. Write: nonparametric bootstrap over groups, with percentile intervals at 2.5, 16, 50, 84 and
97.5%. Counts: 100 replicates for the diagnostic battery with `max_lag = 10`; 200 for the empirical
Gaussian distortion; `min_groups_for_bootstrap = 10`, below which interval slots are NaN-filled, which
is one source of grey cells in the maps. `Model_Parameters_Hat` is **not** bootstrapped by the program;
the raw MLE cloud is what is on disk and any bootstrap of it happens downstream in R.

**The reproducibility problem, and it is real.** There is no run manifest for either production run.
The dispatcher defaults do not match the data: defaults say n_sims 256 and group sizes {1, 10, 100},
the data says 10,000 and {10, 100}. The actual overrides were reconstructed from filenames and CSV
columns and are recorded nowhere. **Write the manifest from that reconstruction and commit it beside
the data.** Then state plainly in Data Availability that output directories are keyed by the engine's
baked git hash, the engine is pinned by tag, and each figure script names its data folder — that part
of the pipeline is genuinely good and should be described. Environment for the record: SLURM, 32 CPUs
per task, 48 GB (96 for the figure-2 dispatch), BLAS threads pinned to 1 with OpenMP taking the
simulation and bootstrap loops, and `MACRODR_AXIS_SERIAL=1`, which is load-bearing because without it
memory grows with concurrent combinations and the jobs OOM.

**The anchors, which Methods has to own.** Map panels come from the Gaussian-Fisher runs; the
time-resolved figures are a separate analytic per-step computation on the dumps, not the numerical
battery. State the anchor per figure. **[SETTLED SINCE]** `433ed13` is the numerical-Fisher demo and
**no paper number is quoted from it**; the verdict was recomputed on `1c2ae6f` + `87889e6` + `0ffbda7`
and survives, which is exactly what `433ed13` was run to show. Say that once, here.

**Numerical safeguards.** One paragraph, because the code runs it and because a reader implementing the
filter hits the simplex problem on day one. The paper needs the *what*, not the how-we-got-here.

**Owed.** The recording-condition table for the usage map: roughly 200 lines of commented derivation
for excised patch, whole cell and oocyte currently live inside `figure_6.Rmd`, backed by 88 sourced
files. Lift it into a Methods table, one source per row. It is the part a referee checks first.

---
