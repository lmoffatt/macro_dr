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
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), CPUS, MEM,
# TIME, PARTITION, BIN. Example:
#   NCHS="10 100" N_SIMS="2048 1024" projects/.../dispatch_figure_2.sh dirac

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
SCRIPT="$HERE/../local/figure_2.macroir"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

# Load the cluster profile for CLUSTER/PARTITION/SCRATCH_MACRO and the repo cd.
# (module loads here are harmless — jobs get modules via run_macroir.sh on the node.)
[ -f /etc/profile ] && source /etc/profile
PROFILE="$(readlink -f "$HERE/../clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[dispatch] no such cluster profile: $PROFILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$PROFILE"
# Passed to run_macroir.sh so it finds the profile from the sbatch spool copy
# (where its own $0 no longer points at the repo).
export MACRODR_PROFILE="$PROFILE"

# Default BIN to this cluster's latest build; override (export BIN=…) to pin one.
BIN="${BIN:-$(readlink -f "build/macrodr_cli-${CLUSTER}-current")}"
[ -x "$BIN" ] || {
    echo "[dispatch] binary not found: $BIN" >&2
    echo "           build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
    exit 1
}

# This run's grid. NCHS and N_SIMS are parallel arrays paired by index:
# job i runs NCHS[i] channels with N_SIMS[i] simulations.
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMS=(${N_SIMS:-1024 1024 1024 1024})

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}

# Shared output dir on scratch; jobs write nch-distinct filenames into it.
WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/eLife_2025"
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"

for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    # get_number(n=...) → size_t; a bare literal would be a double and
    # simulate's n_simulations expects unsigned long.
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_2__nch_%s_nsim_%s"' "$nch" "$nsim")

    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-32G}" \
        --time="${TIME:-1-00:00:00}" \
        --job-name="fig2_nch_${nch}" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE" \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$SCRIPT")

    echo "submitted fig2_nch_${nch}  job=${jobid}  n_sim=${nsim}  -> ${WORKDIR}/figures/data/figure_2__nch_${nch}_nsim_${nsim}_*"
done

echo
echo "watch: squeue -u \$(whoami)   logs: tail -f ${WORKDIR}/logs/slurm-<jobid>.out"
