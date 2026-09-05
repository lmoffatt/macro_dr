#!/bin/bash
# Local runner for the figure_0 parallelism characterization / probe: fills
# the environment the payload needs and runs it in-process (no SLURM).
#
# Usage (from anywhere):
#   projects/macroir_next/ops/local/dispatch_figure_0_local.sh full    # once
#   projects/macroir_next/ops/local/dispatch_figure_0_local.sh         # probe
#
# Outputs land in projects/macroir_next/runs/<commit>/figure_0/:
#   figure_0_table.txt  the measured table (config, threads, s/iter, node rate)
#   tuning.env          CPUS_RECOMMENDED, JOBS_PER_NODE, SECONDS_PER_ITER,
#                       BETA_EFF — dispatch_figure_1_local.sh sources it
#   baseline.env        written by full mode; probes compare against it
#
# Tunables via env: BIN, WORKDIR, TOTAL_CPUS, SLICE, DISC_SLICE, THREADS_LIST,
# PACK_LIST, and the reference-job knobs (NCH, TAU_MS, INTERVAL_IN_TAU, ...).

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../macroir_next/ops/local
BASE="$(readlink -f "$HERE/../../../..")"           # repo base
PROJ="$(readlink -f "$HERE/../..")"                 # .../projects/macroir_next
MODE="${1:-probe}"

BIN="${BIN:-$(readlink -f "$BASE/build/gcc-release/macrodr_cli")}"
[ -x "$BIN" ] || {
    echo "[fig0-local] binary not found: $BIN" >&2
    echo "             build first:  cmake --build --preset gcc-release" >&2
    exit 1
}
commit="$("$BIN" --commit)" || { echo "[fig0-local] '$BIN --commit' failed" >&2; exit 1; }

WORKDIR="${WORKDIR:-$PROJ/runs/$commit}"
mkdir -p "$WORKDIR"

MODE="$MODE" BIN="$BIN" PROJ="$PROJ" WORKDIR="$(readlink -f "$WORKDIR")" \
    SIM_SCRIPT="$HERE/figure_1_sim.macroir" \
    EVI0_SCRIPT="$HERE/figure_0_evidence.macroir" \
    bash "$HERE/../slurm/run_figure_0.sh"
