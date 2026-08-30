# How the figure results were obtained

> **History rewrite, 2026-08-29 — read before resolving any hash in this document.** Repository
> history was rewritten on 2026-08-29. Every commit hash used as a data-directory name below is a
> *pre-rewrite* hash and no longer resolves: `git cat-file -t 433ed13` fails. The directory names on
> disk and the stamps inside the CSVs were **not** rewritten, and must not be, because the stamp is
> the provenance record. To resolve one, look it up in `.git/filter-repo/commit-map` (1084 lines,
> `old<TAB>new`). The mappings for every directory on disk, verified 2026-08-29:
>
> | on disk | rewritten to | date | subject |
> |---|---|---|---|
> | `433ed13` | `41eb69f` | 2026-06-17 | fig3: adaptive Gauss-Newton tolerance |
> | `1c2ae6f` | `26dae12` | 2026-07-05 | gaussian also for corrected covariance |
> | `87889e6` | `f30c8c5` | 2026-07-11 | claude script, minor .gitignore change |
> | `0ffbda7` | `d4a4534` | 2026-07-21 | legacy call sites for the new `variance_form` template parameter |
> | `1f7138b` | `a202289` | 2026-07-31 | macro: restore the interval conductance variance in the non-recursive path; NMR → INR |
> | `a202e03` | `d6402be` | 2026-08-05 | correction in LSE |
> | `ccd26f9` | `2e94384` | 2026-08-05 | many things (the run-3 dumps) |
> | `82b956f` | `6ff3567` | 2026-07-19 | finishing nonlinearsqrfit |
> | `5d9b43b` | `fa90eef` | 2026-06-17 | dispatch for 2 days |
> | `8fc274d` | `6149496` | 2026-06-16 | ECD_Corrected via ECD_Fisher and IDM |
> | `dfa842d` | `c52c1a9` | 2026-06-16 | fix(ops): closing quote in dispatch_figure_3.sh |
> | `a3e0a89` | `313420e` | 2025-12-02 | theoretical results (the NMR regression point) |

> **`NMR` vs `INR`, 2026-07-31 — read before using any `macro_NMR` file.** This document is a record of
> what ran, so its `NMR` references stay. They denote a **defective build**: the non-recursive
> interval member with the `N·ms` interval-variance term missing, which is what the engine computed
> between `a3e0a89` (2025-12-02) and `1f7138b` (2026-07-31). The corrected member is `INR`, which is
> the published `MacroINR` (`nomenclature.md`). So **every `macro_NMR` file on disk is that defective
> build**, and `INR` was re-dispatched from scratch on 2026-07-31; its data carries `macro_INR` and
> a later commit hash. The two are not interchangeable and the `.Rmd` readers must keep telling them
> apart. Nothing else is affected: the fix is confined to the non-recursive path with `av > 0`, so
> `NR`, `R`, `MR`, `VR`, `IR` and `LSE` are unchanged.

> Opened 2026-07-14. Updated 2026-07-20, 2026-07-31, and **audited against disk 2026-08-29** (this
> pass). Shared across the papers; the run manifest they all cite. Companion to the Methods
> (`1_method/docs/manuscript-drafts/sections/06_methods.md` and its `.tex`), not a replacement for it.
>
> The Methods plan (M1 to M11) says what the Methods section must *claim*: the model, the protocol,
> the sweeps, the algorithms, the estimator, the uncertainty quantification. This document says how
> the numbers on disk were actually *produced*: which script, which invocation, which data directory,
> which notebook, which PDF, and what part of that chain can and cannot be reproduced today. It is
> the run manifest M9 asks for, plus the five things M9 does not cover (the R environment, the render
> step, the seed, data availability, and which cluster ran what).
>
> Every value below was read from the run scripts, the dispatchers, the C++ source, `git`, or the
> production CSVs themselves. Path and line are given. Where a code comment and the code disagree,
> the code wins and the disagreement is flagged. Where something could not be established, it says
> so, and those items are collected in §9. Do not cite this file by line number: cite it by section
> anchor (`provenance.md#5`), because the line numbers move on every audit.

## 1. The shape of the pipeline

There are **seven** producing lanes with data on disk, not one and not four. Each is one or more
`.macroir` scripts executed by the `macrodr_cli` binary. Three are launched by hand on a laptop, four
are dispatched to a cluster. They write comma-separated value (CSV) files into
`projects/eLife_2025/figures/data/`, and a set of R Markdown (`.Rmd`) notebooks and plain R scripts
under `projects/eLife_2025/figures/paper_both/` turn those into the figure PDFs.

| Run | Script | Dispatcher | Output prefix | Directories on disk |
|---|---|---|---|---|
| 1 | `ops/local/figure_1_plus_lse.macroir` (and its smaller sibling `figure_1.macroir`) | none, run directly | loose `figure_1_*.csv` | `figures/data/` (14 loose files) |
| 2 | `ops/local/figure_3_mle.macroir` | `ops/slurm/dispatch_figure_3.sh` | `figure_3_` | `433ed13` (487), and the older `8fc274d` (38), `dfa842d` (10), `5d9b43b` (25) |
| 3 | `ops/local/figure_3_time.macroir` and `figure_3_time_LSE_g.macroir` | none, run directly | `figure_3_time_dlik_*` | 10 loose dumps, 998 MB to 1.16 GB each |
| 4 | `ops/local/figure_3_mle_G.macroir` | `ops/slurm/dispatch_figure_3_G.sh` | `figure_3_G_` | see the attribution note below |
| 5 | `ops/local/figure_4.macroir`, `figure_4_LSE.macroir` | `ops/slurm/dispatch_figure_4.sh` | `figure_3_G_` (macro/micro), `figure_3_LSE_` (LSE) | see the attribution note below |
| 6 | `ops/local/figure_3_mle_LSE.macroir` | `ops/slurm/dispatch_figure_3_LSE.sh` | `figure_3_LSE_` | see the attribution note below |
| 7 | `ops/local/figure_3_mle_fisher_only.macroir`, `figure_3_mle_LSE_fisher_only.macroir` | `ops/slurm/dispatch_figure_3_fisher_only.sh` | `figure_3_fim_` | `a202e03` (111 of its 207) |
| 8 | `ops/local/figure_6_traces.macroir`, `figure_6_traces_D0.01.macroir` | none, run directly | `trace_*` | `figures/data/figure_6_traces/` (19 CSVs, 43 MB) |

`ops/slurm/dispatch_figure_3_LSE_numfim.sh` exists and writes the prefix `figure_3_LSE_numfim_`
(`:189`). No file with that prefix is on disk, so that lane has not produced anything that reaches a
figure.

Paths are relative to `projects/eLife_2025/`.

**Attribution note, and it is a real gap.** Three dispatchers write the prefix `figure_3_G_`
(`dispatch_figure_3_G.sh:229` and `dispatch_figure_4.sh:263` via `prefix="figure_3_G_"` at `:223`)
and two write `figure_3_LSE_` (`dispatch_figure_3_LSE.sh:170`, `dispatch_figure_4.sh:222`). The
filename therefore does **not** identify the lane. Only `figure_3_` (run 2) and `figure_3_fim_`
(run 7) are unambiguous. Which of runs 4, 5 and 6 produced `1c2ae6f`, `87889e6`, `0ffbda7`,
`1f7138b` and the LSE half of `a202e03` **is not established anywhere in this repository**,
because the cluster jobs wrote their run ledgers into cluster scratch and those were never copied
back (§2, "The run ledger"). It could be recovered from the SLURM job scripts if any were kept; none
were. `82b956f` is the one that was reconstructed from evidence rather than from a ledger; see §7.1.

### What is actually in `figures/data/`, 2026-08-29

Ten commit-named directories, 2020 CSVs, plus ten loose dumps, plus the loose Figure 1 inputs, plus
about a hundred scratch directories from the exploratory phase. 86 GB in total.

| Directory | CSVs | Size | Row-1 stamp | Prefixes |
|---|---|---|---|---|
| `433ed13` | 487 | 8.1 G | `433ed13`, all | `figure_3_` |
| `1c2ae6f` | 371 | 6.2 G | `1c2ae6f` ×360, **`0ffbda7` ×11** | `figure_3_G_` |
| `0ffbda7` | 322 | 2.0 G | `0ffbda7`, all | `figure_3_G_` ×223, `figure_3_LSE_` ×99 |
| `87889e6` | 261 | 3.9 G | `87889e6`, all | `figure_3_G_` |
| `1f7138b` | 211 | 1.1 G | `1f7138b`, all | `figure_3_G_` |
| `a202e03` | 207 | 1.1 G | `a202e03`, all | `figure_3_fim_` ×111, `figure_3_LSE_` ×96 |
| `82b956f` | 88 | 402 M | `82b956f`, all | `figure_3_G_` ×80, `figure_3_LSE_` ×8 |
| `8fc274d` | 38 | 769 M | `8fc274d`, all | `figure_3_` |
| `5d9b43b` | 25 | 340 M | `5d9b43b`, all | `figure_3_` |
| `dfa842d` | 10 | 182 M | `dfa842d`, all | `figure_3_` |
| `digest` | 9 `.rds` | 103 M | not a CSV | `figure_3_digest_*` |
| `figure_6_traces` | 19 + 2 logs | 43 M | **unstamped** | `trace_*` |

**The eleven mis-stamped files are a live defect.** `figures/data/1c2ae6f/` contains eleven
`macro_R` files whose row 1 reads `0ffbda7`:
`figure_3_G_nch_{100,1000}_nsim_10000_macro_R_noise_100_{battery_pool_G,battery_sim_G,empirical_G,mle_cloud_runs,pool_runs}.csv`
and `figure_3_G_nch_100_nsim_10000_macro_R_noise_1000_mle_cloud_runs.csv`. They were produced by a
later build and copied into an older directory. The file's own stamp is the record, so nothing is
wrong with the *data*, but any statement of the form "directory name equals engine version" is false
for this directory and must not be made. The previous version of this document claimed zero
mismatches across all directories; that claim was wrong when written or has been overtaken.

Three naming traps are built into this pipeline and every one of them has already misled a reader of
the repository. They must not reach the manuscript.

**The notebook called `figure_2.Rmd` is fed by scripts called `figure_3_*`.** Every battery file
anywhere is named `figure_3_...`, because those are the scripts that write them. A Methods paragraph
describing Figure 2 must describe the Gaussian lanes, not a "figure 2" lane. The `figure_2.macroir`
and `dispatch_figure_2.sh` lane is a different, older diagnostic battery that never calls the maximum
likelihood estimation (MLE) stage at all; its outputs sit in the scratch directories and are read only
by archived notebooks.

**The notebooks called `figure_3*.Rmd` are not fed by any MLE battery.** They read digests built from
the time-resolved dumps of run 3, which is a separate script with no dispatcher.

**The axis labelled `noise_in_conductance_tau` is neither the noise parameter nor the dimensionless
noise the axis is named after.** Every dispatcher maps the label to `Current_Noise = label / 1000`
(`dispatch_figure_3.sh:162-172`, `dispatch_figure_3_G.sh:209`, `dispatch_figure_4.sh:242`, and the
same case block in three more scripts). The production "noise 0.1" cell is `Current_Noise = 1e-4`.
Separately, the dimensionless noise the reader should be shown is
`S̃ = Current_Noise · k_off / i²`, and **the stored label is ten times `S̃`**: the dispatcher divides
by 1000 where the definition asks for `k_off/i² = 100`, folding the reference sampling interval into
a quantity that is supposed to be independent of it. The single definition lives in
`figures/paper_both/noise_units.R:8-27`, the full account in `projects/eLife_2025/NOISE_AXIS_UNITS.md`.
Join and filter on the raw label; convert on display only.

## 2. The engine, the script language, and how a grid gets swept

A `.macroir` file is a plain-text script in a small typed domain-specific language (DSL), one
statement per line, which `macrodr_cli` parses, compiles and runs in a single process. There is no
subcommand: `main()` forwards `argv` straight to `main_flow` (`src/cli/main.cpp:3-4`), which treats
every positional argument as a script file and every unrecognised `--` token as an inline DSL line
(`src/cli/cli_parser.cpp:138-151`), then concatenates files and inline lines in `argv` order into one
program (`src/cli/script_loader.cpp:27-41`). So the invocation is simply:

```
build/gcc-release/macrodr_cli ops/local/figure_3_time.macroir
```

**Parameter injection.** Because inline `--name = value` tokens become DSL statements at their `argv`
position, a shell dispatcher can treat a script as a function of its arguments. The dispatchers place
every injection *before* the script path, so the injected assignments are the first statements of the
assembled program, and they define names that the script deliberately leaves undefined (the injectable
names are present in each script but commented out). Assignment is last-writer-wins
(`include/macrodr/dsl/grammar_typed.h:115` uses `insert_or_assign`), which is why each dispatcher
carries a "FILE CONTRACT" comment in its header listing the names the script must not define itself
(for example `dispatch_figure_3.sh:9-13`). If the script were to assign one of those names, it would
silently clobber the injection.

**Axis broadcasting.** The grid is swept inside the process, not by a shell loop. `axis(name, labels)`
declares a named axis and `indexed_double_by`, `indexed_size_by`, `indexed_int_by`, `indexed_bool_by`
and `indexed_string_by` attach one value per label. Whenever any argument of a command is an indexed
value, the whole call is lifted: its index space is the union of the axes of all its arguments, and the
underlying C++ function is invoked once per coordinate of the Cartesian product, returning an indexed
result (`include/macrodr/dsl/grammar_typed.h`, the lifting machinery around the three combo loops at
803, 1184 and 1474; `legacy/indexed.h`). Vector and tuple constructions lift the same way, which is how
one `create_experiment(...)` call fans out over the acquisition intervals. One SLURM job is therefore
one `.macroir` script covering a whole sub-grid, and the shell layer stays deliberately thin.

The axes do not appear in the output filename. The file stem comes from the injected `filepath` string,
and each axis becomes a CSV *column* carrying that coordinate's label
(`include/macrodr/cmd/detail/write_csv_common.h`, the `axis_names_`/`axis_values` machinery from 506).

**`MACRODR_AXIS_SERIAL=1` is load-bearing.** The axis-combination loop is an OpenMP parallel-for by
default. This environment variable flips it to serial (`grammar_typed.h:803-811`, and the two
identical guards at 1184-1192 and 1474-1482), so the inner per-simulation loop becomes the active
parallel level. Without it, each concurrent combination holds its own per-simulation state and memory
grows with the number of combinations, which is how these jobs run out of memory. The SLURM
dispatchers set it in the `sbatch --export` list (`dispatch_figure_3.sh:217`); it is not an `argv`
token, and the compiled-in default is off.

**Provenance stamping, and its three exemptions.** The CSV writers used by the *analysis* commands
begin every file with the build's short git commit hash on a line of its own, before the column header
(`write_provenance_row`, `include/macrodr/cmd/detail/write_csv_common.h:37-39`, called from
`src/core/likelihood.cpp:3036`, `src/core/load_experiment.cpp:107` and three sites inside
`write_csv_common.h` itself). That is why the R notebooks read those files with `skip = 1`, and it is
why the data directories are named after commits: the dispatchers ask the binary for its own hash
(`BIN --commit`) and use it as the output folder name, so two code versions do not normally write into
each other's results (`dispatch_figure_3.sh:65-70`).

The stamp is **not** universal, and three parts of the figure set depend on that:

1. The simulation writer (`src/core/simulate.cpp:545`) does not stamp, so
   `figures/data/figure_1_simulation.csv` starts directly with its header.
2. Neither does the likelihood-diagnostic writer used by run 1. **None** of the fourteen loose
   `figure_1_*.csv` files carries a stamp, checked by `head -1` on each. `figure_1_panels.R:50` reads
   them correctly with `read.csv(f, check.names = FALSE)` and no `skip`.
3. Neither do the nineteen `figures/data/figure_6_traces/trace_*.csv` of run 8, which come through the
   same simulation writer.

A Methods sentence saying "all CSVs carry the hash and are read with `skip = 1`" would be false for
Figures 1 and 6.

**The run ledger, and what it does not cover.** On every invocation the binary snapshots the program
it actually compiled: `runs/run-<YYYYMMDD>-<HHMMSS>/script.macroir` holds the assembled text with the
injections included, and `meta.json` holds the working directory, a Unix timestamp and the full `argv`
(`src/cli/app/workspace_persistence.cpp:38-43`). There are **763** such directories under
`projects/eLife_2025/runs/`.

They are all **local** runs. Grepping every `meta.json` for the script path returns 733 hits across
fourteen `ops/local/*.macroir` scripts and **zero** for `figure_3_mle_G.macroir`, `figure_4.macroir`,
`figure_4_LSE.macroir` or any `fisher_only` script. The cluster jobs wrote their ledgers into cluster
scratch and those were never copied back. So the ledger pins runs 1, 3 and 8 exactly (§4, §6, §8) and
pins nothing at all about runs 2 and 4 through 7.

**Build and dispatch.** `ops/build_cluster.sh <cluster> [tag]` sources `ops/clusters/<cluster>.sh` for
modules, BLAS, scratch and account, resubmits itself as a batch job if it is not already inside an
allocation (login-node builds get killed for memory), and configures into `build/<cluster>-<tag>/`
with `-DMACRODR_GIT_COMMIT_OVERRIDE=<tag>` (`build_cluster.sh:87`), where `tag` defaults to
`git rev-parse --short HEAD` (`:31`). That override is what makes the build directory name, the
binary's `--commit` output and the CSV stamp agree. The SLURM dispatchers submit
`ops/slurm/run_macroir.sh` with one job per (algorithm, channel-number, noise) cell, requesting 32
CPUs, 48 GB (96 GB for the figure-2 lane) and a two-day limit (`dispatch_figure_3.sh:212-214`), and
pinning the BLAS to one thread so OpenMP owns the simulation and bootstrap loops.

## 3. The scientific inputs, in one paragraph

All the MLE lanes share a front end, and the Methods M1 to M6 is the authority on it. In brief: the
compiled-in two-state scheme `scheme_CO` (`legacy/models_simple.h:10-98`), closed to open at rate `on`
times agonist concentration and open to closed at rate `off`, only the open state conducting, all
channels closed at t = 0, white instrumental noise only. Six parameters, all fitted in base-10
logarithmic coordinates, simulated at `on` = 10, `off` = 100, `unitary_current` = 1,
`Current_Baseline` = 1, with `Current_Noise` and `Num_ch_mean` swept. The scripts call
`create_parameters(...)` inline and never load `scheme_CO_par.csv` or `scheme_CO_prior.csv`, so those
files are not what the figures use. There is no prior anywhere in the figure pipeline: every figure is
likelihood-based, none is Bayesian. Ground truth is generated by exact simulation of the
continuous-time Markov chain using uniformization, with 1000 sub-steps per measurement interval in the
ensembles and 250 in the Figure 1 illustration, which is what licenses treating it as the reference
against which the Gaussian approximations are judged.

The macroscopic members are one code path under two flags, `recursive_approximation` and
`averaging_approximation`, where the averaging flag counts conditioned interval endpoints (0, 1 or 2);
`variance_approximation` selects the variance form, which is what separates `VR`. The least-squares
family is `algo_family_approximation = 2` and takes the same averaging flag.

| Label in data | Label in figures | `recursive` | `averaging` | note |
|---|---|---|---|---|
| `nonlinearsqr_g` | LSE | n/a | 0 | family 2, instantaneous |
| `nonlinearsqr` | ILSE | n/a | 1 | family 2, interval-averaged |
| `macro_NR` | NR | false | 0 | |
| `macro_NMR` | NMR (defective build) | false | 1 | do not use, see the banner |
| `macro_INR` | INR | false | 1 | the corrected member |
| `macro_R` | R | true | 0 | |
| `macro_MR` | MR | true | 1 | |
| `macro_VR` | VR | true | 1 | variance form |
| `macro_IR` | IR | true | 2 | |

The `.pre()` helper in `figure_4_common.R:52` routes the two families to their two file prefixes:
anything starting `nonlinearsqr` reads `figure_3_LSE_`, everything else `figure_3_G_`. The crossing of
the two least-squares tokens is a trap held in one place, `figure_1_panels.R` and `figure_2.Rmd:78-81`:
on disk `..._LSE` is averaging = 1, therefore **ILSE**, and `..._LSE_av0` is averaging = 0, therefore
**LSE**.

**The `macro_NMR` / `macro_INR` split is not a relabelling, it is two algorithms** (2026-07-31, see the
banner at the top). Same flags, different code: `NMR` omits the `N·ms` interval-variance term and `INR`
carries it. Until 2026-07-29 the repository also used `NMR` and `MNR` as two spellings of one thing;
that mapping now has to resolve to `INR` for new data and stay on `NMR` for the frozen files.

## 4. Run 1: Figure 1

**Producer.** `ops/local/figure_1_plus_lse.macroir`, run with no dispatcher and no injections. It and
its smaller sibling `figure_1.macroir` **write to the same paths**, so whichever ran last owns the
data. The data currently on disk is stamped 2026-08-05 22:54 and includes `..._LSE_av0.csv`, which only
the `_plus_lse` sibling writes (`:226`), so the sibling is the owner.

**Ledger entry.** `runs/run-20260805-225419/`, whose `meta.json` records
`["../../build/gcc-release/macrodr_cli", "ops/local/figure_1_plus_lse.macroir"]` at timestamp
1785981259 (22:54:19), matching the output mtimes.

```
cd projects/eLife_2025
../../build/gcc-release/macrodr_cli ops/local/figure_1_plus_lse.macroir
```

**What it computes.** A single recording, no grid, no ensemble, no MLE, no bootstrap. **Twenty
channels** (`figure_1_plus_lse.macroir:60`), six measurement intervals at 50 kHz with 100 raw samples
averaged into each (`:62`), giving 2 ms per interval and a 12 ms recording; agonist is 0 during
interval 0 and 10 for intervals 1 to 5. One trajectory is simulated with **250 sub-steps** per interval
(`:75`). That one recording is then passed through the likelihood algorithms via
`calc_likelihood_diagnostic`, writing one long-format CSV per algorithm plus the sub-resolved
trajectory (`figure_1_simulation.csv`). Two further files, `figure_1_likelihood_predictions.csv` and
`figure_1_dlikelihood_predictions.csv`, are written and no notebook reads them.

**The roster on disk is eight files, not six.** `figure_1_plus_lse.macroir:181-226` writes
`figure_1_likelihood_diagnostic_{NR,R,INR,MR,VR,IR,LSE,LSE_av0}.csv`. The **`INR`** file replaced the
old `MNR` one on 2026-08-05; `figure_1_likelihood_diagnostic_MNR.csv` is still on disk with mtime
2026-07-23 and is read by nothing. Three further leftovers, `..._IRV.csv`, `..._MNRV.csv` and
`..._MRV.csv` (2026-06-23), and a bare `figure_1_likelihood_diagnostic.csv` (2025-12-18), are also
dead. VR is the one build that uses `build_likelihood_function_with_family`, because the bool builder
does not expose the variance form.

**One inconsistency to be aware of.** This script builds NR and R with `variance_approximation = 0`
(`:87`, `:93`) while INR, MR, VR and IR get `variance_approximation = 1` (`:99`, `:105`, `:118`,
`:125`); the two least-squares arms are at `:144` and `:177`. The MLE lanes use
`variance_approximation = 1` throughout. Figure 1 is illustrative rather than quantitative, so this
does not propagate into any reported number, but it should not be described as "the same builds as the
rest of the paper". The script's own comment at `:84` records the same thing.

**Notebook.** The notebooks live in `figures/paper_both/` (not `figures/paper/`, renamed with the 1+2
merge; `figures/paper_1/` holds the pre-merge single-paper set) over one shared body,
`figure_1_panels.R`, which defines all eight columns and is untouched by roster changes.
`figure_1.Rmd` renders the paper's Figure 1 and nothing else: **four columns, `NR`, `INR`, `R`, `IR`**
(`figure_1.Rmd:58`), the window × recursion lattice, one `build_figure` call on samples 3 and 4
(`figure_1.Rmd:73`) at 7.0 by 6.05 inches, output `Figure_1.pdf`. A column auto-drops if its dump is
missing (`figure_1.Rmd:63`). `figure_1_all.Rmd`, now in `figures/archive/`, renders every algorithm to
its own `Figure_1_all.pdf` and cannot overwrite the paper's file. The reader carries the algorithm
label in the FILENAME and attaches it POSITIONALLY from two parallel vectors
(`figure_1_panels.R:20-27`), since these CSVs have no `algorithm` column. Rows: prior open probability,
predicted current with its spread and the innovation, posterior (the open-loop columns are marked "no
update (open loop)"), cumulative log-likelihood.

**The two least-squares dumps stay on disk though they are not drawn.** The columns came out on
2026-08-12; Methods and Theory both quote a mean-model identity measured from them (max|LSE − NR| =
8.882e-16 pA, max|ILSE − INR| = 2.220e-16 pA, recorded in the caption source note at
`02_framework.tex:477`), and the roster map in `figure_1_panels.R` still resolves their crossed file
tokens.

**Retired 2026-08-12: the supplement and its caption file.** `Figure_S1.pdf` was withdrawn because it
contained the body figure rather than supplementing it, so the `build_figure(0:4, ...)` call is gone
and the file is deleted; **Figure 1 has no supplement**, and the `\includegraphics` census in §8
confirms none is included. `Figure_1_caption.md` went to
`figures/archive/Figure_1_caption_superseded_20260812.md`. The live caption is the one in
`02_framework.tex`, and there is no caption file for this figure any more.

**This run is not reproducible.** `figure_1_plus_lse.macroir:75` passes `seed = 0`, the
`std::random_device` sentinel (§9). The script's own header says so at `:12`. No number is quoted from
this figure and none should be.

## 5. Run 2: the numerical-Fisher MLE battery (`433ed13`), superseded for the body

**Producer.** `ops/local/figure_3_mle.macroir`, dispatched one SLURM job per cell by
`ops/slurm/dispatch_figure_3.sh` (file prefix `figure_3_`, `:186`). A sequential local equivalent with
the same injection contract is `ops/local/dispatch_figure_3_local.sh`.

**What one job does.** It simulates `n_simulations` recordings over the injected seven-point
acquisition-interval axis, then runs six stages: a per-group Gauss-Newton MLE refit at each
`group_size` (an axis, so one job sweeps several group sizes) giving the estimate cloud; a joint MLE
over all replicates giving the pooled estimate θ_pool; a central-difference numerical Fisher matrix at
both θ_sim and θ_pool; the analytic score at both anchors; the diagnostic battery paired with the
Fisher matrix at each anchor (`_battery_sim`, `_battery_pool`); and the empirical-versus-theoretical
distortion capstone anchored at θ_pool (`_empirical`).

The numerical Fisher matrix is a central difference *of the analytic score* (the score itself comes
from automatic differentiation, so this is one numerical derivative, not two), with a per-coordinate
step h_i = h_rel · max(|θ_i|, 1) (`src/core/likelihood.cpp:1117`) and h_rel injected as 1e-5 through
the `axis_h_fim` axis, symmetrized as (F + Fᵀ)/2 (`likelihood.cpp:1162-1167`). Non-finite columns are
filled with NaN, not zeroed, and the comment above that branch (`likelihood.cpp:1145-1157`) now agrees
with the code; the disagreement flagged in earlier versions of this document has been fixed.

θ_sim is the simulation truth. θ_pool is the single stationary point of the pooled approximate
likelihood. They differ because the Gaussian macroscopic likelihood is misspecified with respect to the
exact simulator, and θ_pool minus θ_sim is that misspecification bias. The mean of the per-group cloud
minus θ_pool is a separate, finite-sample effect. Which anchor a figure uses changes what it can show,
and every figure must state it.

**The grid actually on disk** (99 cells, 487 CSVs), which is *not* the dispatcher's default grid:

- all five algorithms of that era (`NR`, `NMR`, `R`, `MR`, `IR`), N_ch in {10, 100, 1000, 10000}, noise
  labels {0.1, 1, 10}, `n_simulations` = 10000;
- **all five again at `n_simulations` = 1000**, the same four channel counts, noise label 0.1 only
  (20 cells, 100 files). This tier is not IR-only, as an earlier version of this document said;
- `macro_IR` additionally at N_ch in {20, 50, 200, 500} at all three noise labels;
- `macro_IR` at N_ch = 5, complete only at noise label 10; at labels 0.1 and 1 only the estimate cloud
  was written;
- `macro_IR` alone at `n_simulations` = 200, noise label 0.1, four channel counts, as the smallest
  exploratory tier;
- `group_size` = {10, 100} for the 10000- and 1000-recording runs, and {1, 10, 100} for the
  200-recording runs. The cloud CSV carries `group_size` as a column, so **a reader who does not filter
  on it averages two different group sizes together**.

Five CSV families are written per cell: `_mle_cloud_runs`, `_pool_runs`, `_battery_sim`,
`_battery_pool`, `_empirical`. Two cells have only the cloud, which is why the family counts are 99,
97, 97, 97, 97.

**Reconstructed invocation.** The dispatcher's committed defaults (`N_SIMS=256` at `:93`,
`N_ALGO=macro_IR` at `:95`, `GROUP_SIZE=(1 10 100)` at `:104`) do not match the data, so the production
run used environment overrides that were never recorded and no ledger entry exists (§2). The invocation
consistent with what is on disk is:

```
NCHS="10 100 1000 10000" N_SIMS="10000 10000 10000 10000" \
N_NOISE="0.1 0.1 0.1 0.1" N_ALGO="macro_NR macro_NMR macro_R macro_MR macro_IR" \
GROUP_SIZE="10 100" \
projects/eLife_2025/ops/slurm/dispatch_figure_3.sh <cluster>
```

repeated for noise labels 1 and 10, again at `N_SIMS="1000 ..."`, and again with `N_ALGO="macro_IR"`
and the extended channel list. The dispatcher pairs `N_NOISE` with `NCHS` *by index* inside the same
loop (`:145`), so it does not form a cross-product over noise: each noise level is a separate
invocation.

**Who reads it now: nobody in the live set.** `figure_2.Rmd` reaches this directory only through the
dead branch `ANCHOR == "numeric"` (`figure_2.Rmd:59-60`). The live setting is `ANCHOR <- "gaussian"`
(`:46`), which reads the five Gaussian directories (§8). The two anchors were compared and found to
agree at the Figure 2 cell (`1_method/figures_build_plan.md` §3b); that agreement is itself a reported
result, which is why the numeric branch is kept rather than deleted. Three older siblings of this lane,
`8fc274d`, `dfa842d` and `5d9b43b`, hold 73 CSVs between them at `nsim` 200, 1000 and 1024 and feed
nothing.

`figure_5_master*.Rmd` and `figure_6_precision.Rmd` also read `433ed13`, but they were **never** the
paper's figures; both are **archived** under `figures/archive/paper_superseded_20260722/` (verified
present). The first was self-described as a correctness audit
(`figures/in_progress/figure_5_master_STATUS.md:14`) and the second carried the header
`*** DATA SOURCE IS A DRAFT ***`.

## 6. Run 3: the time-resolved dumps, feeding Figure 3 and its three supplements

**Producer.** `ops/local/figure_3_time.macroir`, run directly with no dispatcher, so every value is
script-defined and nothing is injected. A sibling, `figure_3_time_LSE_g.macroir`, adds the
instantaneous least-squares arm.

**Ledger entries, exact.** `runs/run-20260805-231348/meta.json` records
`["../../build/gcc-release/macrodr_cli", "ops/local/figure_3_time.macroir"]` at timestamp 1785982428
(23:13:48); the nine dumps it wrote carry mtimes 23:14 to 23:19.
`runs/run-20260805-231234/` ran `figure_3_time_LSE_g.macroir` and matches the 23:13 mtime of
`figure_3_time_dlik_LSE_g.csv`.

```
cd projects/eLife_2025
../../build/gcc-release/macrodr_cli ops/local/figure_3_time.macroir
```

**What it computes.** One shared ensemble of **1000** recordings (`figure_3_time.macroir:44`) of a
100-channel patch (`:30`) at `Current_Noise` = 1e-4 (`:28`), with 100 measurement intervals of 1 ms
each (50 raw samples at 50 kHz), agonist applied during intervals 20 to 59 (`:33`), for a 100 ms
recording, at 1000 sub-steps per interval (`:47`). That single ensemble is then run through all the
likelihood builds and dumped by `calc_dlikelihood_predictions` into long-format CSVs.

**Ten dumps, not five.** `figure_3_time_dlik_{NR,R,INR,MR,VR,IR,LSE,LSE_av0,LSE_g,NMR}.csv`, 998 MB to
1.16 GB each.

**Provenance, and it is better than it used to be.** Row 1 of nine of the ten reads `ccd26f9-dirty`
(`ccd26f9` → `2e94384`, 2026-08-05). The trailing `dirty` means the binary was built from a working
tree with uncommitted changes, so the hash does not uniquely identify the source. The tenth,
`figure_3_time_dlik_NMR.csv`, still reads `0ffbda7` and dates from 2026-07-23: it is the defective
build (see the banner) and is kept only as a record.

**This run *is* reproducible.** `figure_3_time.macroir:53` is `seed = 20260722`, fixed on 2026-07-22
with a comment at `:49-52` recording exactly why. Repeating the seed regenerates the ensemble. This is
the one production lane where that is true (§9).

**Header drift, three stale numbers, all still there.** The script's own header says a noise of 0.001
(`:7`), a "150-sample protocol" (`:8`) and "300 simulations" (`:10`). The code says `Current_Noise` =
1e-4 (`:28`), 100 intervals of 50 samples (`:33`), and `n_simulations = 1000` (`:44`), and the CSVs
agree with the code. None of the three has reached a caption. Fix the comment.

**Three traps in the dump, all confirmed against the writer and the data.**

1. *Every evolution row is written twice.* `emit_state_rows_with_experiment` (`likelihood.cpp:2484`)
   calls `emit_state_rows_without_experiment` as its first act (`likelihood.cpp:2495`), so each
   per-interval record appears once with a blank segment index and once annotated. Deduplicate by
   averaging over (`simulation_index`, `sample_index`) or by filtering `segment_index == 0`; a naive
   sum double-counts everything, which is where the spurious "score ×2" came from.
2. *There is no per-interval Fisher column.* Per-interval Fisher information is rebuilt in R from the
   dumped sensitivities as I_t = (∂y_mean/∂θ)²/y_var + ½ (∂y_var/∂θ)²/y_var². The 36 entries of the
   whole-recording Gaussian Fisher matrix are present, but only in the `state`-scope rows, and they are
   global per recording, not per step.
3. *Derivatives are zero before the agonist step*, which makes the first 20 intervals look empty rather
   than broken. Mask on FIM < 1e-6, never on a percentage of the peak.

The per-interval log-likelihood and score are incremental, not cumulative. This was verified
numerically: for the first IR recording the 100 deduplicated interval values sum to -147.3591333890,
against a whole-recording total of -147.3591333890.

**Consumers, through a digest layer.** No live notebook reads the 1 GB CSVs. `figure_3_digest.R`
reduces them into `figures/data/digest/figure_3_digest_*.rds`, nine files totalling 103 MB, and
`figure_3.Rmd`, `figure_3_S1_S2.Rmd` and `figure_3_S3.Rmd` all read those
(`../data/digest/figure_3_digest_`). `figure_3.Rmd` renders the body figure's time-resolved calibration
cascade, seeded at `figure_3.Rmd:200` (`set.seed(20260731)`) and `:346` (`set.seed(20260812)`);
`figure_3_S1_S2.Rmd` writes `Figure_3_S1.pdf` and `Figure_3_S2.pdf` (per-step information against score
variance, split across parameters), seeded at `:154`; `figure_3_S3.Rmd` writes `Figure_3_S3.pdf`,
seeded at `:117`. The digest files are `.rds` and therefore gitignored (`.gitignore:161`).

**The finding this lane exists for.** Var(s_t)/E[I_t] ≈ 1 per step for every algorithm, so the global
F ≠ J separation (IR ≈ 1 against NR/INR ≈ 35 to 40) lives entirely in the *temporal correlation* of the
score, not in its per-step magnitude. The clean discriminator is the score autocorrelation: IR ≈ 0.005,
white; the recursive members ≈ 0.78.

## 7. Run 4 and its successors: the Gaussian-Fisher family

**Producers.** `ops/local/figure_3_mle_G.macroir` via `ops/slurm/dispatch_figure_3_G.sh`, and later
`ops/local/figure_4.macroir` / `figure_4_LSE.macroir` via `ops/slurm/dispatch_figure_4.sh`. All three
write the prefix `figure_3_G_` (macro/micro) or `figure_3_LSE_` (least squares), so **the file does not
say which dispatcher produced it** (§1). This is run 2 with the finite-difference Fisher stage deleted
and every distortion diagnostic re-anchored on the analytic Gaussian Fisher matrix
G_b = Σ_t GFI_t, with GFI_t = (∂μ_t/∂θ)(∂μ_t/∂θ)ᵀ/σ²_t + (∂σ²_t/∂θ)(∂σ²_t/∂θ)ᵀ/(2σ⁴_t)
(`legacy/qmodel.h:6148-6162`, the formula itself at 6152-6153). Because G_b follows from the analytic
derivative alone, the numerical Fisher routine, which costs 2·n_params extra likelihood passes per
recording, is never called, and there is no `axis_h_fim` column in the output.

The mode is selected by call *arity*, not by a flag: the Gaussian overload of
`likelihood_derivative_basic_diagnostics` is registered separately
(`src/cli/command_manager.cpp:127` and its registration at `:1533-1536`) and passes an empty Fisher
vector. This is worth knowing before anyone tries to switch modes by editing a boolean.

**Coverage, by directory, read off the `_battery_pool_G` files at `nsim` 10000.**

`1c2ae6f` (371 CSVs):

| Algorithm | N_ch | Noise labels |
|---|---|---|
| `macro_IR` | 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000 | 0.05, 0.1, 0.2, 0.5, 1, 10 |
| `macro_NR`, `macro_R`, `macro_MR` | 10, 100, 1000, 10000 | 0.1 |
| `macro_R` (mis-stamped, §1) | 100, 1000 | 100 |
| `macro_NMR` | none | none |

**`macro_NR` at 10000 channels is present**, contrary to what earlier versions of this document said;
`figure_3_G_nch_10000_nsim_10000_macro_NR_noise_0.1_*` is a complete five-family set. What is genuinely
absent from `1c2ae6f` is `macro_NMR`, and that is correct and permanent: the Gaussian family was later
filled for `NMR` in `87889e6`, and that fill is the defective build; the corrected member `INR` was run
into `1f7138b`.

`0ffbda7` (322 CSVs) carries the noise sweep out to label 10⁷ for `macro_IR`, `macro_NR`, `macro_R` and
`nonlinearsqr`, adds `macro_VR` at labels 0.1 and 1, and adds N_ch = 5 cells.

`1f7138b` (211 CSVs) is the post-fix run: `macro_INR` over N_ch {10, 100, 1000, 10000} at labels 0.1
through 10⁷, plus `macro_MR` and `macro_VR` over the same high-noise range.

`87889e6` (261 CSVs) holds `macro_NMR` at all four N_ch × labels {0.1, 1, 10}, `macro_NR`, `macro_R`
and `macro_MR` at labels {1, 10}, `macro_R` and `macro_IR` at small N_ch, and 13 `micro_R` / `micro_IR`
files.

`a202e03` (207 CSVs) is the `fisher_only` lane (111 files, prefix `figure_3_fim_`) plus 96
`figure_3_LSE_` files. It is **the only directory carrying the un-averaged least-squares arm**
(`nonlinearsqr_g`). A smoke run at `nsim` 8 (`macro_VR`, N_ch 10) is mixed in; do not glob the
directory blindly.

**Consumers: the whole Figure 4, 5, 6 and 7 set.** `figure_4_common.R:44-45` defaults to
`c("../data/1c2ae6f", "../data/87889e6", "../data/0ffbda7")` as a search path, first hit wins, and
seven notebooks override it with the five-directory form (§8). `figure_7.Rmd:91` uses the default
three. The grid is auto-detected from disk at render time: `figure_4_common.R` scans `DATA_DIRS` and
includes every (N_ch, noise) cell that has **both** `battery_pool_G` and `battery_sim_G` at `nsim`
10000, so a new run appears in the figure by itself.

**`a202e03` belongs on every `FIG4_DATA_DIRS`, drawn or not.** The reason is not the roster: the shared
source CSVs under `figures/figure_4_source_data/` are stamped by the count of input files
(`figure_4_data.R:246-252`, `fig4_stamp(paths)`), so a notebook that cannot see that directory rebuilds
them WITHOUT the un-averaged arm, and the next figure to render then draws one column fewer with
nothing to say so. That is a silent regression across the whole set from one missing path in one
notebook.

### 7.1 `figures/data/82b956f` — the constant-relative-noise probe at 1000 recordings

**What it is.** 88 CSV files, 421 MB, produced by `ops/local/figure_3_mle_G.macroir` (80 files) and `ops/local/figure_3_mle_LSE.macroir` (8 files). Same schema and same five output kinds as the other Gaussian-anchor directories, but a different experiment: it is a four-cell **diagonal** at constant relative noise, not a rectangle.

**Code.** Row 1 of all 88 files reads `82b956f`, **clean** — no `-dirty`, which the build would have appended for a modified tracked tree (`cmake/GenerateGitCommit.cmake:26-34`). After the 2026-08-29 history rewrite that hash resolves through `.git/filter-repo/commit-map`:

```
82b956fe1eb1bfc05dbe16f308f6c543d7c29d8b -> 6ff3567035ac14b18888a46a479b362eadaf98d1
6ff35670  2026-07-19 13:55:31 -0300  "finishing nonlinearsqrfit"
```

The **dispatch scripts were newer than the binary**. At `6ff3567` the noise `case` in `dispatch_figure_3_G.sh` accepted only labels up to 100 (1000–100000 arrived in `7d9e6618`, 2026-07-20 00:43) and `dispatch_figure_3_LSE.sh` did not yet exist (added in `7efa5d8c`, 2026-07-20 00:17). So the shell layer came from a tree at or after 2026-07-20 00:43 while the binary was a clean `82b956f` build. This is the ordinary BIN-pinning path: the folder name is the binary's own `--commit`, never the tree's.

**When.** 2026-07-20 13:50:53 to 2026-07-21 08:38:20, 18.8 h, one campaign. mtimes cluster into 20 job groups, one per (algorithm, N_ch): the four `nonlinearsqr` jobs first (13:50–14:17), then `macro_IR`, `macro_R`, `macro_MR`, `macro_NMR`. The groups overlap in time, so at least two ran concurrently; that excludes the serial local driver `ops/local/run_figure_3_G_local.sh`.

**Coverage.**

- Members: `macro_R`, `macro_MR`, `macro_NMR`, `macro_IR` (20 files each: 4 cells × `mle_cloud_runs`, `pool_runs`, `battery_sim_G`, `battery_pool_G`, `empirical_G`) and `nonlinearsqr` (8 files: 4 cells × `mle_cloud_runs`, `pool_runs`).
- Grid, four cells, `noise label = 10 · N_ch` throughout (instrumental-to-gating ratio `r = 10` constant): `(N_ch, label)` = (10, 100), (100, 1000), (1000, 10000), (10000, 100000). The label maps to `current_noise = label/1000`, and to the dimensionless `S̃ = label/10`, so `S̃ = N_ch` in every cell.
- `n_simulations = 1000` in every file. This is the only Gaussian-anchor directory at 1000; the others are at 10000.
- Interval axis: the full seven values 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01 in units of τ = 1/k_off. Six parameters. `group_size` ∈ {10, 100}, with no group size 1.

**Two caveats on the contents.**

1. The `nonlinearsqr` cells carry only 2 of the 5 kinds because at `82b956f` the derivative path refused that family: `calculate_mdlikelihood_predictions_visit` (`src/core/likelihood.cpp`) returns `error_message("unsupported for family==nonlinearsqr")`, and the script calls `calc_dlikelihood_predictions` right after the `_pool` write, so each LSE job aborted there. The two files it did write are complete (all seven intervals, both group sizes).
2. `82b956f` sits inside the `macro_NMR` defect window — `a3e0a89` (2025-12-02) is its ancestor and it is an ancestor of `1f7138b` (2026-07-31), both verified with `git merge-base`. Its 20 `macro_NMR` files are therefore the defective non-recursive build described at the head of this document, and are not `INR`.

**What reads it.** `figures/paper_1/figure_4_diag_nsim1000.Rmd` (and its archived twin) reads the directory by literal path; that notebook is explicitly a diagnostic and is not on any figure's path. Thirteen scripts under `figures/in_progress/ir_reliability_20260828/` list it last in `DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")` and glob `macro_IR.*battery_sim_G.csv` from it — but **none of its rows survives their `group_by(nch,z) |> slice(1)`**: the 86 matching files reduce to 80 unique cells, won by `1c2ae6f` (60) and `0ffbda7` (20), and all four `82b956f` cells already exist in `0ffbda7` at nsim 10000, which precedes it in `DIRS`. It contributes nothing to the digests behind Figure 6—figure supplement 1. Outside the paper path it is genuinely used by `tmp/refute_occ_clt/extract.py` and `tmp/fitgrade/rlog_defect.Rmd`, both of which keep `nsim` as a column. `paper_both/figure_4_logL_data.R` and `figure_4_logL.Rmd` name it only to exclude it.

**What could not be established.**

- **The launching invocation.** No run ledger entry exists. `projects/eLife_2025/runs/` (763 entries) has nothing between 2026-07-19 00:35:24 and 2026-07-21 17:33:46, which brackets the whole window, and no `meta.json` mentions `figure_3_mle_G`, `figure_3_G` or `82b956f`. This is structural rather than a lost record: `persist_run_workspace` (`src/cli/app/workspace_persistence.cpp:23`) writes the ledger to `runs/` relative to the process cwd, and `run_macroir.sh` cds into `$WORKDIR`, which the dispatcher sets to `${SCRATCH_MACRO}/eLife_2025/<commit>`. A cluster run's ledger stays on cluster scratch. None of `1c2ae6f`, `0ffbda7`, `87889e6` or `1f7138b` has a repo ledger entry either. The argv, and with it the `NCHS` / `N_SIMS` / `N_NOISE` / `N_ALGO` / `GROUP_SIZE` actually passed, is gone; the grid above is read off the filenames and confirmed against the data columns.
- **Which machine ran it.** No SLURM log, no `slurm-*.out`, and nothing else in the repo from the window except unrelated source edits and audio files. `dirac` is the only cluster profile in the tree, but no evidence connects this run to it. All that is established is that two or more jobs ran concurrently.
- **The exact dispatcher revision.** Only a lower bound: at or after `7d9e6618` (2026-07-20 00:43).
- **Bit reproducibility.** `figure_3_mle_G.macroir` passes `seed = 0`, the `std::random_device` sentinel. The drawn seed is written to no file. These numbers cannot be regenerated.
- **Why it was run at 1000 recordings.** The header of `figure_4_diag_nsim1000.Rmd` gives a purpose for the data (a constant-`r` probe before committing roughly 3300 CPU-hours), but that notebook was written afterwards against data already on disk; it is a reader's account, not a record of the launch.

## 8. Which figures exist, and which do not

The canonical arc is `1_method/docs/manuscript-drafts/sections/04_results.tex` and the figure **set**
(numbering, body-vs-supplement) is owned by `1_method/decisions.md` "The figure set"; both supersede
the figure list in `1_method/00_plan.md` §5.

**Nineteen display items.** The list below is the complete `\includegraphics` census of
`papers/1_method/docs/manuscript-drafts/sections/*.tex` (one in `02_framework.tex`, eighteen in
`04_results.tex`), each matched to the producer whose `ggsave` writes that exact name.

| Display item | Producer (`figures/paper_both/`) | Graphics file | Data dirs the producer declares |
|---|---|---|---|
| Figure 1 | `figure_1.Rmd` | `Figure_1.pdf` | none declared (loose `figure_1_*` via `figure_1_panels.R`) |
| 2 | `figure_2.Rmd` | `Figure_2.pdf` | a202e03 1f7138b 1c2ae6f 0ffbda7 87889e6 (`:56`) |
| 2—fig. supp. 1 | `figure_2_S1.Rmd` | `Figure_2_S1.pdf` | 0ffbda7 and 1c2ae6f, two cells named in the file |
| 3 | `figure_3.Rmd` | `Figure_3.pdf` | `../data/digest` |
| 3—fig. supp. 1 | `figure_3_S1_S2.Rmd` | `Figure_3_S1.pdf` | `../data/digest` |
| 3—fig. supp. 2 | `figure_3_S1_S2.Rmd` | `Figure_3_S2.pdf` | `../data/digest` |
| 3—fig. supp. 3 | `figure_3_S3.Rmd` | `Figure_3_S3.pdf` | `../data/digest` |
| 4 | `figure_4.Rmd` | `Figure_4.pdf` | 1c2ae6f 87889e6 0ffbda7 1f7138b a202e03 (`:331`) |
| 4—fig. supp. 1 | `figure_4_S1.Rmd` | `Figure_4_S1.pdf` | same five (`:114`) |
| 4—fig. supp. 2 | `figure_4_S2.Rmd` | `Figure_4_S2.pdf` | same five (`:100`) |
| 4—fig. supp. 3 | `figure_4_S3.Rmd` | `Figure_4_S3.pdf` | same five (`:167`) |
| 4—fig. supp. 4 | `figure_4_S4.Rmd` | `Figure_4_S4.pdf` | same five (`:138`) |
| 5 | `figure_5.Rmd` | `Figure_5.pdf` | `../figure_4_source_data` (`:44`) |
| 5—fig. supp. 1 | `figure_5_S1.Rmd` | `Figure_5_S1.pdf` | `../figure_4_source_data` (`:95`) |
| 5—fig. supp. 2 | `figure_5_S2.Rmd` | `Figure_5_S2.pdf` | `../figure_4_source_data` (`:44`) |
| 6 | `figure_6.Rmd` | `Figure_6.pdf` | none declared; parses `figure_4_common.R` (`:50`) |
| 6—fig. supp. 1 | `figure_6_S1.R` | `Figure_6_S1.pdf` | `../in_progress/ir_reliability_20260828/digest` |
| 6—fig. supp. 2 | `figure_6_traces.R` | `Figure_6_S2.pdf` | `../data/figure_6_traces` + the same digest |
| 7 | `figure_7.Rmd` | **`Figure_7_hatch.pdf`** | 1c2ae6f 87889e6 0ffbda7 (`:91`) |

**Four corrections this pass made to this table.** (i) `Figure_4_S4`, `Figure_5_S2`, `Figure_6_S1` and
`Figure_6_S2` are display items and were missing. (ii) `Figure_5_S3.pdf` was listed and is **not**
included by any `.tex`; its producer `figure_5_S3.Rmd` still exists and still declares the five-directory
path (`:131`), so it renders, but nothing consumes it. (iii) The Figure 7 row said `Figure_7.pdf`. The
manuscript includes `Figure_7_hatch` (`04_results.tex:1663`), and `figure_7.Rmd:1324` selects between
the two names on a `HATCH_BAD` flag with the hatch sibling as the body figure since 2026-08-26. The
plain `Figure_7.pdf` is on disk and unused. (iv) The row "7—fig. supp. 1 | `figure_7_supplement_standard_error.Rmd`"
named a **file that does not exist**; only its stale `.html` remains. There is no Figure 7 supplement.

Also on disk and not display items: `Figure_4_logL.pdf` (from `figure_4_logL.Rmd`), `Figure_6_traces.pdf`
and `Figure_6_traces_D0.01.pdf` (working output of `figure_6_traces.R`), `Figure_7.pdf`,
`Figure_7_hatch.png`, `Figure_7.png`, and the three `Figure_6_{frontiers,lines,regions}` pairs that
`figure_7.Rmd` also writes.

**Do not hand-edit this table.** The display-item column is every `\includegraphics` in
`sections/*.tex`; the producer column is the file whose `ggsave` writes that exact name. Regenerate it
by running both greps, or it drifts, and the drift changes decisions: the July table sent
`figure_3_supplement_1.Rmd` to a dead PDF no producer writes.

**Naming rule, 2026-08-12.** Every display item is `Figure_<parent>[_S<n>]` and its producer
`figure_<parent>[_S<n>].Rmd`, with `figure_3_S1_S2.Rmd` carrying both of the pages it writes in its
name. Rename producer and output together and never one alone: renaming only the file on disk recreates
the defect the hand-copied `Figure_4_S3.pdf` had, a display item no `ggsave` writes, which a re-render
silently fails to update. Superseded files are suffixed `_superseded_20260812` in the archive.

**Every graphic in `figures/paper_both/` that is not a display item was moved on 2026-08-12 to
`figures/archive/not_display_items_20260812/`**, **118 files**. The rule, from Luciano: a graphic is
either a body figure, or a figure supplement, or it lives in the archive. Producers were NOT moved,
because their `../data/...` paths are relative and moving them breaks the read.

## 9. The R layer, and what is not reproducible

**Packages.** The notebooks load `tidyverse`, `patchwork`, `data.table` and `scales`, and call into
`MASS`, `dplyr`, `tibble`, `ggplot2` and `knitr`. There is no `renv.lock`, no `.Rprofile`, no
`sessionInfo()` dump and no R project file anywhere in the repository, so no R version is pinned to the
rendered PDFs. The PDFs were rendered between **2026-07-31 and 2026-08-29**; the machine currently
carries R 4.6.1 (2026-06-24), ggplot2 4.0.3, patchwork 1.3.2, data.table 1.18.4, rmarkdown 2.30, scales
1.4.0 and tidyverse 2.0.0, but nothing ties those versions to the files. Figures 6—supp. 1 and 2 are
drawn with `device = pdf` and Helvetica; check the PDF, not the on-screen device.

**There is a render step in version control, and it is stale.**
`figures/paper_both/render_figure_4_set.sh` is tracked and drives
`Rscript -e "rmarkdown::render('<nb>.Rmd', quiet = TRUE)"` with a two-layer staleness check (the script
decides whether to start R; `cached()` in `figure_4_common.R` decides whether to re-parse the CSVs into
`figures/.cache`). It runs one notebook per R process on purpose, because the roster globals persist for
the session. But its `JOBS` list names ten notebooks of which **five no longer exist**
(`figure_4_supplement_standard_error`, `_lag_kappa`, `_se_kappa`, `_coverage`, `_qq_grid`,
`_other_parameters`) and **omits `figure_4_S4`**, which writes a real display item. Everything outside
the Figure 4 set is still knit by hand; the reconstructed command, from `figures/paper_both/`, is
`Rscript -e 'rmarkdown::render("figure_3.Rmd")'`.

**Random seeds: three regimes, not one.** `calc_seed(0)` draws from `std::random_device`
(`legacy/mcmc.h:37-45`), so zero is the sentinel meaning "seed randomly", and the resolved seed is
written neither to `meta.json` nor to the CSVs.

- **Not reproducible.** Run 1 (`figure_1_plus_lse.macroir:75`, `seed = 0`) and the MLE batteries
  (`figure_3_mle.macroir:74`, `figure_3_mle_G.macroir:74`, `figure_4.macroir:74`,
  `figure_4_LSE.macroir:81`, all `seed = 0`). A rerun of the same script with the same code gives a
  statistically equivalent, numerically different ensemble. This cannot be fixed retroactively for
  `433ed13`, `1c2ae6f`, `0ffbda7`, `1f7138b` or `87889e6`.
- **Reproducible.** Run 3 (`figure_3_time.macroir:53`, `seed = 20260722`, fixed 2026-07-22) and run 8
  (`figure_6_traces.macroir:14-39` and the `_D0.01` sibling, `seed = 20260829`).
- **Injected.** The `fisher_only` lane takes `sim_seed` as a dispatcher injection
  (`figure_3_mle_fisher_only.macroir:125`), which also pairs the ensemble across algorithms. That is the
  pattern the other lanes should adopt; its header at `:52-56` explains why.

The same literal zero passed to the diagnostic, MLE and bootstrap commands goes straight into `mt_64i(0)`
and *is* deterministic, so the resampling is reproducible given the data. On the R side there is **no
`set.seed(1)` anywhere**: the Figure 3 family is seeded with dated values (`figure_3.Rmd:200` and `:346`,
`figure_3_S1_S2.Rmd:154`, `figure_3_S3.Rmd:117`), and the Figure 4 and Figure 5 families are unseeded.

**The data are not in version control.** `.gitignore:52` is `*.csv` and `.gitignore:161` is `*.rds`.
None of the 2020 CSVs in the ten commit-named directories, the ten dumps, the nine digest `.rds` files,
the `figure_6_traces` CSVs or the loose Figure 1 inputs is tracked: 86 GB in total. The figure PDFs are
tracked (25 in `paper_both/`); their inputs are not, and no repository, archive or digital object
identifier is named anywhere for them. Whatever is said in Data Availability has to be made true first.

**The environment is under-recorded.** Five cluster profiles exist (`ops/clusters/{capitan,clementina,dirac,serafin,tupac}.sh`)
and two are known to have been used; the profiles pin modules (`dirac.sh:15-19` uses gnu14 with Intel
MKL 2019.5.281 and CMake 3.30.8; `clementina.sh:26-28` uses GCC 15.1.0 with OpenBLAS 0.3.29 and CMake
3.31.11), but **which cluster produced which data directory is written nowhere**, and neither is which
dispatcher did (§1, attribution note). The local builds behind runs 1, 3 and 8 used an unrecorded
compiler; `CMakePresets.json:20` pins the path `/usr/bin/g++` with no version, and the standard is
C++20. The copy from cluster scratch back into the repository is not scripted anywhere either, which is
also how eleven files ended up in the wrong directory (§1).

**Cost is unrecorded.** The SLURM requests are known (`dispatch_figure_3.sh:212-214`: 32 CPUs, 48 GB,
two-day limit; 96 GB for the figure-2 lane) but no job output was kept, so actual wall-clock and CPU
consumption cannot be reported. This is not recoverable.

## 10. What to fix before submission

Ordered by how much a referee would care.

1. **Fix the two wrong citations this document put into the manuscript.** `06_methods.tex:664` quotes
   `legacy/qmodel.h:6096-6098` for the Gaussian Fisher definition; the correct range is 6148-6162, the
   formula at 6152-6153. `06_methods.tex:590` cites `provenance.md:138` by line number; cite
   `provenance.md#5` instead, since this file is re-audited and its lines move.
2. **State the seed regime per figure, not globally.** Figures 3 and 6 are reproducible from a fixed
   seed; Figures 1, 2, 4, 5 and 7 are not, because their ensembles were drawn from `random_device` and
   the resolved value was never logged. The honest sentence is that the ensembles are large enough
   (1000 to 10000 recordings) that the reported statistics are stable, but those datasets cannot be
   regenerated bit for bit. The `sim_seed` injection of the `fisher_only` lane is the fix for any rerun.
3. **Deposit the data.** 86 GB across ten commit-named directories plus ten dumps. Until they are
   archived and cited, the provenance story in §2, which is genuinely good, is a story about files on
   one laptop.
4. **Resolve the eleven mis-stamped files.** Either move the eleven `0ffbda7`-stamped `macro_R` files
   out of `1c2ae6f/` into `0ffbda7/`, or state in the deposit that the directory name is a scratch
   grouping and the file's own row 1 is the provenance record. `figure_2.Rmd:50-55` already takes the
   second position explicitly; make it the documented one.
5. **Commit the run manifest, and start logging the cluster ledgers.** §5 and §7 reconstruct
   invocations from filenames because the cluster `runs/` directories were never copied back and the
   dispatcher defaults actively contradict what was run. Copying `runs/` back with the data costs
   nothing and closes the largest gap in this document.
6. **Repair `render_figure_4_set.sh`.** Five of its ten jobs name notebooks that no longer exist and
   `figure_4_S4` is missing, so the one scripted render step in the repository does not render the
   current set.
7. **Pin the R environment.** An `renv.lock` and a one-line render script for the notebooks outside the
   Figure 4 set would close the last unscripted hop in the chain.
8. **Fix the stale comments** identified here: the three wrong numbers in `figure_3_time.macroir:7-10`,
   and the header of `figure_3_mle.macroir`.
9. **Decide what to do with the orphans.** `Figure_5_S3.pdf` renders and nothing includes it;
   `figure_7_supplement_standard_error.html` has no producer; `Figure_4_logL.pdf` is not a display item.
   Either promote them or archive them, per the rule in §8.
