#!/bin/bash
# Submit the figure_0 parallelism characterization (or its per-campaign probe)
# on a cluster. One exclusive-node job; the payload run_figure_0.sh does all
# configurations inside the allocation (scaling sequentially, packing
# concurrently), reads iterations/second from the binary's own __i_iter.csv,
# and writes tuning.env (+ baseline.env in full mode) under
# <scratch>/macroir_next/<commit>/figure_0/. dispatch_figure_1.sh sources that
# tuning.env for its CPUS default.
#
# Usage (from the repo base):
#   projects/macroir_next/ops/slurm/dispatch_figure_0.sh <cluster> full   # once
#   projects/macroir_next/ops/slurm/dispatch_figure_0.sh <cluster>        # probe
# Local run (no SLURM), same payload:
#   MODE=full BIN=build/gcc-release/macrodr_cli WORKDIR=projects/macroir_next/runs/local \
#     PROJ=projects/macroir_next SIM_SCRIPT=... EVI0_SCRIPT=... \
#     projects/macroir_next/ops/slurm/run_figure_0.sh
# Tunables via env: TOTAL_CPUS, SLICE, THREADS_LIST, PACK_LIST, TOL, PARTITION,
# ACCOUNT, TIME, BIN, RUN_DIR, plus the reference-job knobs (NCH, TAU_MS, ...).

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"
PROJ="$(readlink -f "$HERE/../..")"
PAYLOAD="$HERE/run_figure_0.sh"
SIM_SCRIPT="$(readlink -f "$HERE/../local/figure_1_sim.macroir")"
EVI0_SCRIPT="$(readlink -f "$HERE/../local/figure_0_evidence.macroir")"

CLUSTER="${1:?Usage: $0 <cluster> [full|probe]}"
MODE="${2:-probe}"

PROFILE="$(readlink -f "$HERE/../../../eLife_2025/ops/clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[fig0] no such cluster profile: $PROFILE" >&2; exit 1; }
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$PROFILE"
set -e
export MACRODR_PROFILE="$PROFILE"

BIN="${BIN:-$(readlink -f "build/macrodr_cli-${CLUSTER}-current")}"
[ -x "$BIN" ] || {
    echo "[fig0] binary not found: $BIN" >&2
    echo "       build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
    exit 1
}

if ! commit="$("$BIN" --commit)"; then
    echo "[fig0] could not query commit hash" >&2; exit 1
fi
run="${RUN_DIR:-$commit}"
case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/macroir_next/$run" ;;
esac
mkdir -p "$WORKDIR/figure_0"

TOTAL_CPUS="${TOTAL_CPUS:-32}"

jobid=$(sbatch --parsable \
    --partition="${PARTITION:-batch}" \
    ${ACCOUNT:+--account="$ACCOUNT"} \
    --exclusive \
    --cpus-per-task="$TOTAL_CPUS" \
    --time="${TIME:-1:30:00}" \
    --job-name="fig0_${MODE}" \
    --output="$WORKDIR/figure_0/slurm-%j.out" \
    --export=ALL,MODE="$MODE",BIN="$BIN",WORKDIR="$WORKDIR",PROJ="$PROJ",SIM_SCRIPT="$SIM_SCRIPT",EVI0_SCRIPT="$EVI0_SCRIPT",MACRODR_PROFILE="$PROFILE",TOTAL_CPUS="$TOTAL_CPUS",CLUSTER="$CLUSTER" \
    "$PAYLOAD")
echo "[fig0] submitted $MODE as job $jobid; results in $WORKDIR/figure_0/"
echo "[fig0] chain the campaign with: DEPEND=$jobid projects/macroir_next/ops/slurm/dispatch_figure_1.sh $CLUSTER"
