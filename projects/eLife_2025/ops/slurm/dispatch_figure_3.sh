#!/bin/bash
# Dispatch figure_3_mle (the MLE validation pipeline) across the grid — one
# SLURM job per (algorithm, Nchannels). Clone of dispatch_figure_2.sh: same
# cluster-profile sourcing, same run_macroir.sh wrapper, same printf injection
# idiom (axis before the indexed_*_by that references it), same load-bearing
# MACRODR_AXIS_SERIAL=1. Adds the per-group-MLE knobs (group_size, bootstrap
# counts, GN max_iter) as injected size_t literals.
#
# FILE CONTRACT — figure_3_mle.macroir must NOT define the injected names itself
# (they must stay commented/undefined): axis_Nchanels, Num_ch, axis_noise,
# current_noise, n_simulations, filepath, algorithm_axis, algo_*_approximation,
# axis_interval, exp_n_step_/exp_n_samp_, axis_h_fim, h_rel_value, group_size_axis,
# group_size, n_bootstrap_samples, min_groups_for_bootstrap, gn_max_iter.
#
# Prereq: a built binary for the cluster:
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/slurm/dispatch_figure_3.sh <cluster>     # e.g. dirac
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), N_NOISE,
# N_ALGO, H_RELS, GROUP_SIZE (a DSL axis — broadcast within one job; cloud CSV
# gains a group_size column), N_BOOT, MIN_GROUPS, GN_MAX_ITER, CPUS, MEM, TIME,
# PARTITION, BIN. Example:
#   NCHS="10000" N_SIMS="256" N_ALGO="macro_IR" GROUP_SIZE="1 10 100" \
#     projects/eLife_2025/ops/slurm/dispatch_figure_3.sh dirac

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
SCRIPT="$HERE/../local/figure_3_mle.macroir"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

PROFILE="$(readlink -f "$HERE/../clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[dispatch] no such cluster profile: $PROFILE" >&2; exit 1; }
# /etc/profile.d scripts (e.g. debuginfod.sh runs `cat /dev/null /etc/debuginfod/*.urls`
# on an empty glob → non-zero) and Lmod return non-zero; under `set -eo pipefail`
# that silently aborts before any job is submitted. Guard the env setup, restore after.
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$PROFILE"
set -e
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

# Per-commit output isolation. Ask the binary for its baked git hash (the same
# string it stamps as row 1 of every CSV) so two concurrently-running commit
# versions of macrodr_cli write to DISJOINT folders. RUN_DIR overrides the folder:
#   RUN_DIR=egesij56                  -> continue an old run under that commit's
#                                        scratch folder with a NEWER binary
#                                        (continuation over old data)
#   RUN_DIR=/abs/path/to/run          -> use that directory verbatim
if ! commit="$("$BIN" --commit)"; then
    echo "[dispatch] could not query commit hash: '$BIN --commit' failed" >&2
    exit 1
fi
[ -n "$commit" ] || { echo "[dispatch] '$BIN --commit' returned empty" >&2; exit 1; }
run="${RUN_DIR:-$commit}"

# This run's grid. NCHS and N_SIMS are parallel arrays paired by index:
# job i runs NCHS[i] channels with N_SIMS[i] simulations. figure_3 defaults to
# FEWER simulations than figure_2 — the per-group MLE refit is the expensive part
# (N_groups = n_simulations/group_size GN refits), so keep n_simulations modest.
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMS=(${N_SIMS:-256 256 256 256})
N_NOISE=(${N_NOISE:-0.1 0.1 0.1 0.1})
N_ALGO=(${N_ALGO:-macro_IR})
# h_rel: relative step for the central-difference numerical Fisher (injected as a
# single-value axis_h_fim so the output matches figure_2). Sweep with H_RELS="1e-5 1e-6".
H_RELS=(${H_RELS:-1e-5})

# Per-group MLE knobs (the cost levers — injected as size_t literals).
# GROUP_SIZE is a DSL axis (broadcast within one job): 1 = per-replicate; >1 pools
# recordings per refit. Multi-value sweeps the cloud over group size in one run,
# e.g. GROUP_SIZE="1 10 100".
GROUP_SIZE=(${GROUP_SIZE:- 1 10 100})
N_BOOT="${N_BOOT:-100}"                # bootstrap replicates over groups
MIN_GROUPS="${MIN_GROUPS:-10}"         # below this, probit slots NaN-filled
GN_MAX_ITER="${GN_MAX_ITER:-100}"      # GN per-group refit iteration cap

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}

# Shared output dir on scratch, under the per-commit (or RUN_DIR) folder; jobs
# write nch/algo-distinct filenames into it. A bare RUN_DIR/commit nests under
# scratch; an absolute RUN_DIR is used as-is.
case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/eLife_2025/$run" ;;
esac
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"
echo "[dispatch] commit=${commit}  run=${run}  WORKDIR=${WORKDIR}"

join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }

# Loop-invariant injections (built once). axis_h_fim must precede h_rel_value.
axis_h_arg=$(printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= [%s])' "$(join_qcsv "${H_RELS[@]}")")
h_rel_arg=$( printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$(join_csv "${H_RELS[@]}")")
# Per-group MLE knobs: get_number(n=...) → size_t (a bare literal would be a double).
# group_size is a DSL AXIS (group_size_axis must precede the indexed_size_by that
# references it). calc_MLE broadcasts over it, so one job sweeps all group sizes and
# the cloud CSV gains a group_size column. Multi-value via GROUP_SIZE="1 10 100".
group_size_axis_arg=$(printf -- '--group_size_axis = axis(name= "group_size", labels= [%s])' "$(join_qcsv "${GROUP_SIZE[@]}")")
gs_arg=$( printf -- '--group_size = indexed_size_by(axis= group_size_axis, values=[%s])' "$(join_csv "${GROUP_SIZE[@]}")")
nboot_arg=$( printf -- '--n_bootstrap_samples = get_number(n=%s)' "$N_BOOT")
mingrp_arg=$(printf -- '--min_groups_for_bootstrap = get_number(n=%s)' "$MIN_GROUPS")
gnmaxit_arg=$(printf -- '--gn_max_iter = get_number(n=%s)' "$GN_MAX_ITER")

for j in "${!N_ALGO[@]}"; do

for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"

    algo="${N_ALGO[$j]}"

    # Map the algorithm label to its (recursive, averaging, taylor, micro) flags.
    case "$algo" in
        macro_NR)  recursive=false; averaging=0 ; taylor=false ; micro=false ;;
        macro_R)   recursive=true;  averaging=0 ; taylor=false ; micro=false ;;
        macro_NMR) recursive=false; averaging=1 ; taylor=false ; micro=false ;;
        macro_MR)  recursive=true;  averaging=1 ; taylor=false ; micro=false ;;
        macro_IR)  recursive=true;  averaging=2 ; taylor=false ; micro=false ;;
        macro_IRT) recursive=true;  averaging=2 ; taylor=true  ; micro=false ;;
        micro_R)   recursive=true;  averaging=0 ; taylor=false ; micro=true ;;
        micro_MR)  recursive=true;  averaging=1 ; taylor=false ; micro=true ;;
        micro_IR)  recursive=true;  averaging=2 ; taylor=false ; micro=true ;;
        *) echo "[dispatch] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,IR,IRT})" >&2; exit 1 ;;
    esac

    case "$nnoise" in
        0.1)  vnoise=0.0001;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[dispatch] unknown noise level '$nnoise' (want 0.1, 1, 10, 100)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_3_nch_%s_nsim_%s_%s_noise_%s"' "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    micro_arg=$( printf -- '--algo_micro_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$micro")

    # interval_in_tau grid (identical to figure_2). axis_interval precedes exp_n_*.
    axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= ["1","0.5","0.2","0.1","0.05","0.02","0.01"])')
    exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2,4,10,20,40,100,200])')
    exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')

    # MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop so the
    # per-simulation / bootstrap loops become the active OpenMP level (load-
    # bearing; without it memory ∝ concurrent combos → OOM). figure_3's MLE
    # refit has a smaller sim count than figure_2, so mem/time defaults are lower.
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-48G}" \
        --time="${TIME:-12:00:00}" \
        --job-name="f3_${nch}c_${nsim}s" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",MACRODR_AXIS_SERIAL=1 \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$micro_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$axis_h_arg" "$h_rel_arg" \
        "$group_size_axis_arg" "$gs_arg" "$nboot_arg" "$mingrp_arg" "$gnmaxit_arg" \
        "$SCRIPT")

    echo "submitted fig3_nch_${nch}  job=${jobid}  n_sim=${nsim}  algo=${algo}  gs_axis=[${GROUP_SIZE[*]}]  -> ${WORKDIR}/figures/data/figure_3_nch_${nch}_nsim_${nsim}_${algo}_noise_${nnoise}_*"
done
done
