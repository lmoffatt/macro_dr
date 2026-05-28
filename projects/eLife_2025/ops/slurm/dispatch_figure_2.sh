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
N_SIMS=(${N_SIMS:-4096 4096 4096 4096})
N_NOISE=(${N_NOISE:-0.1 0.1 0.1 0.1})
N_ALGO=(${N_ALGO:-macro_R })
#N_ALGO=(${N_ALGO:-macro_IR macro_R macro_MR macro_NR macro_NMR})


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
    nnoise="${N_NOISE[$i]}"
 for j in "${!N_ALGO[@]}"; do

    algo="${N_ALGO[$j]}"
    
    # Map the algorithm label to its (recursive, averaging) settings. Mirrors the
    # commented-out indexed_bool_by/indexed_int_by rows in figure_2.macroir.
    case "$algo" in
        macro_NR)  recursive=false; averaging=0 ;;
        macro_R)   recursive=true;  averaging=0 ;;
        macro_NMR) recursive=false; averaging=1 ;;
        macro_MR)  recursive=true;  averaging=1 ;;
        macro_IR)  recursive=true;  averaging=2 ;;
        *) echo "[dispatch] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,IR})" >&2; exit 1 ;;
    esac

   case "$nnoise" in
        0.1)  vnoise=0.0001;;
        1)  vnoise=0.001;;
        10)  vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[dispatch] unknown noise level '$nnoise' (want 0.1, 1, 10, 100)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    # get_number(n=...) → size_t; a bare literal would be a double and
    # simulate's n_simulations expects unsigned long.
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_2_nch_%s_nsim_%s_%s_noise_%s"' "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    # algorithm — injected the same way as noise: the axis, plus the (recursive,
    # averaging) settings the label maps to. algorithm_axis must come first; the
    # other two reference it.
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")


    # MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop so the
    # per-simulation / bootstrap loops become the active OpenMP level. Without
    # it the combo loop stays parallel and memory ∝ concurrent combos → OOM
    # (cf. jobs 110042–110045, killed at ~129 GB). Load-bearing, not optional.
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-96G}" \
        --time="${TIME:-2-00:00:00}" \
        --job-name="f_${nch}c_${nsim}s" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",MACRODR_AXIS_SERIAL=1 \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" \
        "$SCRIPT")

    echo "submitted fig2_nch_${nch}  job=${jobid}  n_sim=${nsim}  -> ${WORKDIR}/figures/data/figure_2_nch_${nch}_nsim_${nsim}_*"
done
done

echo
echo "watch: squeue -u \$(whoami)   logs: tail -f ${WORKDIR}/logs/slurm-<jobid>.out"
