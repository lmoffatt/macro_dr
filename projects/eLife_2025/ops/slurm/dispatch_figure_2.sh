#!/bin/bash
# Dispatch figure_2 across Nchannels values — one SLURM job per Nch.
#
# Each job injects the figure's named parameters before the .macroir (the file
# behaves as a function of them). Injection order matters: axis_Nchanels must
# precede Num_ch, which references it.
#
# FILE CONTRACT — figure_2.macroir must NOT define these itself (else its own
# assignment clobbers the injection):
#   axis_Nchanels, Num_ch, n_simulations, filepath
# In particular, comment out the in-file `axis_Nchanels = axis(...)` line.
#
# Prereq: a built binary for the cluster (run once, when code changes):
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/slurm/dispatch_figure_2.sh <cluster>     # e.g. dirac
# Tunables via env: NCHS, N_SIMULATIONS, CPUS, MEM, TIME, PARTITION, BIN.

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
SCRIPT="$HERE/../local/figure_2.macroir"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

# Load the cluster profile for CLUSTER/PARTITION/SCRATCH_MACRO and the repo cd.
# (module loads here are harmless — jobs get modules via run_macroir.sh on the node.)
[ -f /etc/profile ] && source /etc/profile
PROFILE="$HERE/../clusters/${CLUSTER}.sh"
[ -f "$PROFILE" ] || { echo "[dispatch] no such cluster profile: $PROFILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$PROFILE"

# Default BIN to this cluster's latest build; override (export BIN=…) to pin one.
BIN="${BIN:-$(readlink -f "build/macrodr_cli-${CLUSTER}-current")}"
[ -x "$BIN" ] || {
    echo "[dispatch] binary not found: $BIN" >&2
    echo "           build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
    exit 1
}

# This run's grid
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMULATIONS="${N_SIMULATIONS:-1024}"

# Shared output dir on scratch; jobs write nch-distinct filenames into it.
WORKDIR="${SCRATCH_MACRO:-/scratch/$USER/macro_dr}/eLife_2025"
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"

for nch in "${NCHS[@]}"; do
    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = %s' "$N_SIMULATIONS")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_2_nch_%s"' "$nch")

    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-32G}" \
        --time="${TIME:-1-00:00:00}" \
        --job-name="fig2_nch_${nch}" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR" \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$SCRIPT")

    echo "submitted fig2_nch_${nch}  job=${jobid}  n_sim=${N_SIMULATIONS}  -> ${WORKDIR}/figures/data/figure_2_nch_${nch}_*"
done

echo
echo "watch: squeue -u \$USER   logs: tail -f ${WORKDIR}/logs/slurm-<jobid>.out"
