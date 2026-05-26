# eLife_2025 ops — building and running macro_dr on clusters

How we build `macrodr_cli` on a cluster and run `.macroir` scripts under SLURM.

## Philosophy

**The `.macroir` file is the unit of work.** A modern `.macroir` (e.g. `local/figure_2.macroir`)
declares its own axes — algorithm, scheme, noise, Nchannels — and `macrodr_cli` dispatches
that grid internally across the task's CPUs via OpenMP. So the default is **one `.macroir`
= one single-node SLURM job**; there is no bash-level scheme/experiment/continuation fan-out
like the old `projects/p2x2/ops/multi_task/` scripts.

Fan out into multiple jobs only when you actually need to:
- the full grid won't fit one node's cores / memory / walltime,
- you want fault isolation (one grid point crashing shouldn't kill the figure), or
- variants are *structural* (different scheme/experiment shape) rather than a single axis sweep.

## Layout

```
projects/eLife_2025/ops/
├── clusters/<name>.sh     # per-cluster profile: modules, PATH_MACRO, CLUSTER, PARTITION, …
├── build_cluster.sh       # build macrodr_cli into build/<cluster>-<git-hash>/
├── slurm/run_macroir.sh   # SBATCH wrapper: runs a .macroir (+ optional --key=value injections)
├── local/*.macroir        # the scripts themselves (figure_2.macroir, …)
└── README.md              # this file
```

(`projects/p2x2/ops/` holds the prior project's copies; it is historical — ignore it.)

## One-time per cluster: the profile

`clusters/<name>.sh` is sourced before building or submitting. It must:
`cd` to the repo, `module load` the toolchain, and `export PATH_MACRO`, `CLUSTER`,
a default `PARTITION`, and `SCRATCH_MACRO`. To add a cluster, copy `dirac.sh` and adjust.

```bash
source projects/eLife_2025/ops/clusters/dirac.sh
which gcc g++ ninja cmake     # sanity-check the toolchain resolved
```

## Build

`macrodr_cli` is heavy: the worst translation units (`command_manager.cpp`, `qmodel.h`
users) need **~16 GB each**. **Never build on a login node** — it gets OOM-killed.
Build inside an allocation:

```bash
salloc --partition=batch --cpus-per-task=2 --mem=64G --time=08:00:00
source projects/eLife_2025/ops/clusters/dirac.sh
projects/eLife_2025/ops/build_cluster.sh dirac
```

- Output: `build/dirac-<git-hash>/macrodr_cli`, plus a `build/macrodr_cli-dirac-current` symlink.
- Tag defaults to the git short hash → **a new commit means a new build dir = full recompile**
  (not incremental). Old build dirs stay intact, so in-flight jobs keep their binary.
- Parallelism: `-j` = `SLURM_CPUS_PER_TASK` (so `--cpus-per-task=2` → `-j 2`, which is safe;
  `-j 4` OOMs even at 64 GB). Override without re-allocating: `BUILD_JOBS=1 build_cluster.sh dirac`.
- Custom tag: `build_cluster.sh dirac nyquist` → `build/dirac-nyquist/`.

## Run a single .macroir

`run_macroir.sh` is the SBATCH script. It reads three env vars and takes the binary's
argument vector as positionals.

| env       | meaning                                                              |
|-----------|----------------------------------------------------------------------|
| `CLUSTER` | selects `clusters/<CLUSTER>.sh`                                      |
| `BIN`     | absolute path to the `macrodr_cli` to pin this job to                |
| `WORKDIR` | dir to `cd` into first (default cwd); outputs land here → use scratch |

```bash
cd ~/Code/macro_dr/macro_dr
sbatch --partition=batch --cpus-per-task=32 --mem=32G --time=1-00:00:00 --job-name=fig2 \
  --export=ALL,CLUSTER=dirac,BIN=$PWD/build/dirac-<hash>/macrodr_cli,WORKDIR=/scratch/$USER/macro_dr/eLife_2025 \
  projects/eLife_2025/ops/slurm/run_macroir.sh \
  projects/eLife_2025/ops/local/figure_2.macroir
```

Monitor: `squeue -u $USER`, then `tail -f slurm-<jobid>.out`. You can log out; the job survives.
Pinning `BIN` to an explicit build dir means a later rebuild can't disturb a running job.

## Parameterize: inject `--key=value`, treat the file as a function

Positional args are the binary's argv **in order**. Each is either a **script file**
(absolutized) or an **inline injection** (`--key=value`, the `--` is stripped and the rest
is concatenated into the program at that position). The file behaves like a function whose
parameters are supplied by injections.

**Ordering rule (critical):** assignments use last-writer-wins, so an injection only takes
effect on uses that come *after* it. Therefore the parameters a `.macroir` exposes must be
**used but not re-assigned** in the file — inject them *before* the file, and they become the
definition. Inject **scalars** (an injected value can't reference identifiers defined later
in the file):

```
# in the .macroir:  Num_ch = indexed_double_by(axis= axis_Nchanels, values=[Num_ch_value])
run_macroir.sh '--Num_ch_value = 100' figure_2.macroir
```

Defaults live in the **caller**, layered as injections (later wins, all before the file):

```bash
run_macroir.sh \
  '--Num_ch_value = 100' '--noise_value = 0.0001' \   # baseline defaults
  '--Num_ch_value = 1000' \                            # this job's override
  figure_2.macroir
```

Quoting: single-quote a literal injection (`'--Num_ch_value = 100'`) — no backslashes.
Only string-valued injections need inner DSL quotes; prefer to avoid them (see fan-out below).
A `concat(a, b)` DSL function exists for building strings (e.g. output paths) from parts.

## Fan out across jobs

Loop in bash, one `sbatch` per grid point. Put per-job **uniqueness in `WORKDIR`**, not in a
string injection — that keeps every injection numeric and quote-free, and each job writes its
`figures/data/…` into its own directory:

```bash
BIN=$PWD/build/dirac-<hash>/macrodr_cli
for nch in 1 10 100 1000 10000; do
  sbatch --partition=batch --cpus-per-task=32 --mem=32G --time=1-00:00:00 --job-name=fig2_$nch \
    --export=ALL,CLUSTER=dirac,BIN=$BIN,WORKDIR=/scratch/$USER/macro_dr/eLife_2025/nch_$nch \
    projects/eLife_2025/ops/slurm/run_macroir.sh \
    "--Num_ch_value = $nch" \
    projects/eLife_2025/ops/local/figure_2.macroir
done
```

For many *string-shaped* knobs (scheme names, algorithm labels), prefer **template rendering**
— a `*_template.macroir` with `@PLACEHOLDER@`s, `sed`-substituted into a concrete file per job
(quotes live in the template, no shell escaping; the rendered file is a provenance record).

## dirac quick reference

- Login: `ssh <user>@dirac.df.uba.ar` · repo at `~/Code/macro_dr/macro_dr` · docs `dirac.df.uba.ar`
- Toolchain (in `clusters/dirac.sh`): `gnu14` + `mkl/2019.5.281` + `gsl` + `ninja` + `cmake/3.30.8`,
  `BLA_VENDOR=Intel10_64lp_seq`. Mostly-Intel cluster, so MKL is the BLAS.
- Storage: home **162 GB** (source) · `/scratch/$USER` up to **8 TB** (outputs; files >30 days
  auto-deleted). Always write outputs under `SCRATCH_MACRO=/scratch/$USER/macro_dr`.
- Partitions (max CPU/job · jobs/user · max walltime):

  | partition | CPU | jobs | walltime          |
  |-----------|-----|------|-------------------|
  | batch     | 128 | 5    | 4 d (default 2 d) |
  | debug     | 4   | 1    | 2 h               |
  | long      | 8   | 2    | 10 d (default 3 d)|
  | gpu       | 64  | —    | 2 d               |
  | gpu_debug | —   | 1    | 1 h               |

## Troubleshooting

- **`cc1plus … señal Terminado (killed)`** during build — OOM. You're on a login node, or
  `-j` too high. Use `salloc --cpus-per-task=2 --mem=64G`; if a single TU still dies, `BUILD_JOBS=1`.
- **`Job … exceeded its time limit`** — `--time` too short, or you ran in a too-short `salloc`.
  batch allows up to 4 days; raise `--time`. For unattended work prefer `sbatch` over interactive
  `salloc` (survives disconnects).
- **`COLORTERM: variable sin asignar`** then silent exit — fixed; scripts use `set -eo pipefail`
  (no `-u`) because OpenHPC's `/etc/profile.d/*` references unset vars.
- **`No existe el fichero … clusters/dirac.sh`** after `salloc` — you're in `~`, not the repo;
  `salloc` drops you in home. `cd ~/Code/macro_dr/macro_dr` first.
- **A `--key=value` override seems ignored** — it was placed *after* the file (too late) or the
  file re-assigns that variable. Inject *before* the file, and don't assign the param inside it.
