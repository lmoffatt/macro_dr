# ops

Same contract as projects/eLife_2025/ops:

- `local/`: .macroir configs plus dispatch scripts that run the binary
  directly. The dispatch scripts INJECT the swept names via environment
  knobs; the .macroir file must not define the injected names itself
  (see the FILE CONTRACT header in
  projects/eLife_2025/ops/local/dispatch_figure_3_local.sh).
- `slurm/`, `clusters/`: empty until cluster runs. Use
  projects/eLife_2025/ops/build_cluster.sh and the profiles in
  projects/eLife_2025/ops/clusters rather than duplicating them here.

Run from the project directory so that data/ paths resolve and the
binary writes runs/run-*/ snapshots (invoke with --env-save).
