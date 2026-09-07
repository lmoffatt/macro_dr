# Project - MacroIR Next Execution

## Identity

- Project id: `macroir_next`
- Program: `macroir`
- Type: execution / validation / paper-support
- Status: active (created 2026-09-04)
- Publication relevance: direct (next paper after the 2026 macro validation
  paper, which after the eLife desk rejection of 2026-09-03 goes to the
  Journal of General Physiology)

## Purpose

Own the operational path for the post-eLife program: Bayesian evidence
computation validated end to end, the CCO/COC equilibrium null test, the
Bessel acquisition filter, and from there conformational models with
inter-subunit interaction and, later, a CUDA port.

Phase sequence (from the 2026-09-04 weekly review audio):

1. Smoke test of `thermo_evidence_dts` after the 2026-09-02 retarget
   (canonical IR member: adaptive=0, recursive=1, averaging=2,
   variance_approximation=1, variance_correction_approximation=0).
2. CCO vs COC evidence null test on stationary data, WITHOUT Bessel.
   Expected result: the log-evidence difference stays bounded (order 1
   nat) and does not grow with record length; twin parameters in
   `data/models/scheme_COC_twin_par.csv` (COC exactly equivalent to
   scheme_CCO at 10 uM). A growing difference indicts the likelihood
   approximation, not the equivalence theorem.
3. Same with the Bessel filter, after validating the filter against the
   channel time constants.
4. Conformational models version 2.0 (inter-subunit interaction) with
   calibrated evidence; then the CUDA port.

## Canonical theory notes (pointers, dated 2026-09-04)

- Evidence machinery audit, retarget, telescopic (stepping-stone)
  estimator, save_Score: `theory/macroir/notes/evidence_dts_audit_20260901.md`
- CCO/COC equivalence theorem, twin map, 3-state classification, the two
  experiments: `theory/macroir/notes/cco_coc_discrimination_plan.md`
- Bessel plan: `theory/macroir/notes/macroir_bessel_plan.md`

## Layout and conventions

Same architecture as `projects/eLife_2025`:

- `data/experiments/`, `data/models/`: inputs referenced by RELATIVE path
  from the project directory. Scripts run with cwd = this directory.
- `ops/local/`: `.macroir` configs + local dispatch scripts (env-knob
  injection contract as in eLife_2025/ops/local).
- `ops/slurm/`, `ops/clusters/`: empty until cluster runs; reuse
  `projects/eLife_2025/ops/build_cluster.sh` and its cluster profiles
  rather than duplicating them.
- `runs/`: created by the binary itself when invoked with `--env-save`
  from this directory (one `run-YYYYMMDD-HHMMSS/` per invocation with
  script.macroir, environment.json, meta.json, snapshots).
- `figures/data/<git-hash>/`: digested outputs; `figures/in_progress/`:
  analysis in development (the R readers for the telescopic bracket and
  the score tests land here).

HYGIENE RULE (learned from eLife_2025, which needed a cleanup plan):
no outputs in the project top level, ever. Everything lands under
`runs/`, `logs/`, or `figures/data/<hash>/`.

## Seeds (2026-09-04)

`data/models/` was seeded from the REPO-ROOT copies of the scheme CSVs,
NOT from eLife_2025/data/models: the eLife copies of scheme_COC predate
the parameter renaming (they still say kon/koff/gating_on/gating_off,
while the current model in legacy/models_simple.h names them
on/off/inactivating_on/inactivating_off). The root copies match the
current names.

- `scheme_CCO_par.csv`, `scheme_CCO_prior.csv`: kon 6.73, koff 166,
  gating_on 157, gating_off 45.3 (log10-normal priors, variance 2).
  gating_on was 743 until 2026-09-07 (the e743fd5 campaign): that puts
  P_open at 0.826 at 10 uM and saturates the occupancy; 157 gives 0.50, the
  design value of the eLife figures. Both files and the twin pair below
  are written by `make_cco_coc_twin.py` in this directory.
- `scheme_COC_par.csv`, `scheme_COC_prior.csv`: placeholder values
  (copied from CCO historically), kept for reference.
- `scheme_COC_twin_par.csv`, `scheme_COC_twin_prior.csv`: the COC exactly
  equivalent to scheme_CCO at 10 uM agonist (on 2.926624, off 27.859054,
  inactivating_on 17.440946, inactivating_off 361.033762 for gating_on
  157; for 743 it was on 5.423005, off 9.347195, inactivating_on
  35.952805, inactivating_off 922.069949), from the closed-form twin map;
  equilibrium equivalence checked by the generator before writing
  (P_open and K_out coincide; the full check to 1e-15 is
  tmp/twin_coc_from_cco.py, 2026-09-01). Prior centered on the twin,
  variance 2, mirroring the CCO prior shape.

Open detail to settle at the smoke test: Current_Baseline is 0 in the
root/seeded CSVs (transformed_mean = -inf under Log10) while the eLife
copies used 1; check how the prior machinery treats the -inf row before
trusting evidence numbers.

## Experiments: inline vs file (checked 2026-09-04)

Experiments are normally DESIGNED INLINE in the .macroir via
`create_experiment(experiment_structure = [{n_step, n_samp, agonist},...],
frequency_of_sampling, initial_agonist, ...)`, with the segments injected
by the dispatch script (that is how every eLife_2025 lane works;
data/experiments/ there only holds the real recording used by figure 1).

RESOLVED 2026-09-04: `thermo_evidence_dts` historically took the
experiment only as a FILE (experiment_file_type). It now has a second
overload (same DSL name, dispatched by argument type) that takes the
inline `Experiment` object directly, same shape as the eLife likelihood
commands, so the evidence lanes can sweep the experiment via
dispatch-injected segments exactly like the figure lanes did. The core
was factored into `run_thermo_evidence_dts(..., const Experiment&, ...)`
and both entries call it (legacy/CLI_thermo_evidence_dts.h). The inline
entry is stateless: it does not write the restart txt, so
`thermo_evidence_dts_continuation` applies to file-based runs only.
`data/experiments/` stays only for real recordings, if any; the
stationary protocol lives in the .macroir via create_experiment.

## Figure 1 lane (written 2026-09-04, pending review and build)

`ops/local/`: `figure_1_sim.macroir` + `figure_1_evidence.macroir` +
`dispatch_figure_1_local.sh`. `ops/slurm/`: `dispatch_figure_1.sh <cluster>`
+ payload `run_figure_1.sh` (the three stages inside one allocation), same
grid and knobs; cluster profiles referenced from eLife_2025/ops/clusters.
The sim seed is baked into every label and filename in BOTH dispatchers,
so a local and a dirac run with the same BASE_SEED pair up file by file;
byte-identical results additionally require equal thread counts (one RNG
stream per OMP thread: local THREADS = cluster CPUS) and a numerically
identical BLAS. Design: evidence confusion matrix,
{truth CCO, COC-twin} x {episodic, stationary} x 10 replicas, both models
fitted to each replica (paired Delta logZ). Stationary realizes the
unmeasured equilibration as a nan-masked segment of the recording (the
likelihood's isnan gating propagates without update, zero logL). Injection
contract identical to eLife's dispatch_figure_3_local.sh; outputs under
`runs/<commit>/` (hygiene rule). Smoke line:
`REPLICAS=1 TRUTHS="CCO" PROTOCOLS="stationary" MAX_ITER=200
 ops/local/dispatch_figure_1_local.sh`.
Everything waits on Luciano's build of the 2026-09-02/04 code batch
(retarget + telescopic columns + save_Score full-ladder + inline
Experiment overload + ladder-reset tolerance).
