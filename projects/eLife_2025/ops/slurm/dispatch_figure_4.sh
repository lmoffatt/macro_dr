#!/bin/bash
# Dispatch figure_3_mle_G (the GAUSSIAN-Fisher MLE validation pipeline) across the
# grid — one SLURM job per (algorithm, Nchannels). Dedicated clone of
# dispatch_figure_3.sh: same cluster-profile sourcing, run_macroir.sh wrapper,
# printf injection idiom, and MACRODR_AXIS_SERIAL=1. The ONLY difference from the
# numerical-Fisher dispatcher is that it targets figure_3_mle_G.macroir and does
# NOT inject axis_h_fim / h_rel_value (there is no finite-difference Fisher here).
#
# FILE CONTRACT — figure_3_mle_G.macroir must NOT define the injected names itself
# (they must stay commented/undefined): axis_Nchanels, Num_ch, axis_noise,
# current_noise, n_simulations, filepath, algorithm_axis, algo_*_approximation,
# axis_interval, exp_n_step_/exp_n_samp_, group_size_axis, group_size,
# n_bootstrap_samples, min_groups_for_bootstrap, gn_max_iter.
# (No axis_h_fim / h_rel_value — no numerical Fisher.)
#
# Prereq: a built binary for the cluster:
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh <cluster>     # e.g. dirac
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), N_NOISE,
# N_ALGO, GROUP_SIZE (a DSL axis — broadcast within one job; cloud CSV gains a
# group_size column), INTERVALS (the acquisition-interval axis, in units of tau;
# any multiple of 1/500, default the historical seven-value grid),
# N_BOOT, MIN_GROUPS, GN_MAX_ITER, CPUS, MEM, TIME, PARTITION,
# BIN, DEPEND (job id to wait for — unset = no dependency). Example:
#   NCHS="10000" N_SIMS="256" N_ALGO="macro_IR" GROUP_SIZE="1 10 100" \
#     projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh dirac
#   # chain after a build job:  DEPEND=<jobid> projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh dirac

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
# LEAN figure_4 producers (only mle_cloud + battery_sim_G + battery_pool_G; NO Stage-6 empirical
# capstone). One dispatch, all algorithms: family==2 (LSE) needs i/Current_Noise Fixed, so it routes
# to figure_4_LSE.macroir; macro/micro use figure_4.macroir. The two differ only in those transforms.
SCRIPT_MACRO="$HERE/../local/figure_4.macroir"
SCRIPT_LSE="$HERE/../local/figure_4_LSE.macroir"

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
#   RUN_DIR=/abs/path/to/run          -> use that directory verbatim
if ! commit="$("$BIN" --commit)"; then
    echo "[dispatch] could not query commit hash: '$BIN --commit' failed" >&2
    exit 1
fi
[ -n "$commit" ] || { echo "[dispatch] '$BIN --commit' returned empty" >&2; exit 1; }
run="${RUN_DIR:-$commit}"

# Optional SLURM job dependency: hold every dispatched job until another job
# finishes — e.g. chain figure_3_G after a compile/build job. Set DEPEND to either
# a bare job id (→ afterok) or a full SLURM dependency expression (passed verbatim).
DEP_SPEC=""
if [ -n "${DEPEND:-}" ]; then
    case "$DEPEND" in
        *[!0-9]*) DEP_SPEC="$DEPEND" ;;         # has a non-digit → full SLURM expr
        *)        DEP_SPEC="afterok:$DEPEND" ;; # bare job id → afterok
    esac
    echo "[dispatch] job dependency: --dependency=$DEP_SPEC"
fi

# This run's grid. NCHS and N_SIMS are parallel arrays paired by index.
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMS=(${N_SIMS:-256 256 256 256})
N_NOISE=(${N_NOISE:-0.1 0.1 0.1 0.1})
N_ALGO=(${N_ALGO:-macro_IR})

# Per-group MLE knobs (the cost levers — injected as size_t literals).
GROUP_SIZE=(${GROUP_SIZE:- 1 10 100})
N_BOOT="${N_BOOT:-0}"                  # Stage-1 cloud bootstrap OFF by default (no figure reads the cloud Probit CIs)
MIN_GROUPS="${MIN_GROUPS:-10}"         # below this, probit slots NaN-filled
GN_MAX_ITER="${GN_MAX_ITER:-100}"      # GN per-group refit iteration cap

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}

# Shared output dir on scratch, under the per-commit (or RUN_DIR) folder.
case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/eLife_2025/$run" ;;
esac
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"
echo "[dispatch] commit=${commit}  run=${run}  WORKDIR=${WORKDIR}  (gaussian, no numerical Fisher)"

join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }

# Per-group MLE knobs: get_number(n=...) → size_t (a bare literal would be a double).
# group_size is a DSL AXIS (group_size_axis must precede the indexed_size_by that
# references it). calc_MLE broadcasts over it, so one job sweeps all group sizes.
group_size_axis_arg=$(printf -- '--group_size_axis = axis(name= "group_size", labels= [%s])' "$(join_qcsv "${GROUP_SIZE[@]}")")
gs_arg=$( printf -- '--group_size = indexed_size_by(axis= group_size_axis, values=[%s])' "$(join_csv "${GROUP_SIZE[@]}")")

# ---- interval axis, DERIVED rather than tabulated ------------------------------------------------
# The seven-value grid that used to be hardcoded here is not a table, it is a formula. At the
# configured sampling frequency (50 kHz) and tau = 1/k_off = 0.01 s, one acquisition interval of
# n_samp raw samples is n_samp/50000 s, i.e. exactly n_samp/500 in units of tau. Hence
#
#     n_samp   = 500 * delta            (raw samples averaged into one interval)
#     n_step_1 =   2 / delta            (intervals in the pre-agonist segment)
#     n_step_2 =   4 / delta = n_step_3 (intervals in each agonist segment)
#
# so the recording holds a fixed number of RAW samples and only their grouping changes. The old
# hardcoded arrays are reproduced exactly by these three expressions; verified for all seven values.
#
# CONSEQUENCE, and it is a hard limit rather than a rounding convenience: only multiples of 1/500 are
# reachable, because n_samp is a count. delta = 0.002 is one raw sample per interval and is the floor;
# anything finer needs a higher acquisition frequency in the .macroir, not a change here.
#
#   INTERVALS="1 0.05"      # run two intervals instead of the full grid
#   INTERVALS="1 0.5 0.2 0.1 0.05 0.02 0.01"   # the default, i.e. the historical grid
#
# COST NOTE: the likelihood runs once per INTERVAL, so cost scales with sum(1/delta), not with the
# number of intervals. On the default grid delta=0.01 alone is 53% of the work and delta=1 is 0.5%.
INTERVALS=(${INTERVALS:-1 0.5 0.2 0.1 0.05 0.02 0.01})

iv_lab=(); iv_samp=(); iv_step1=(); iv_step2=()
for d in "${INTERVALS[@]}"; do
    read -r nsamp st1 st2 eff <<<"$(awk -v d="$d" 'BEGIN{
        if (d <= 0) { print "0 0 0 0"; exit }
        ns = d*500; nsr = int(ns + 0.5)
        printf "%d %d %d %.10g", nsr, int(2/d + 0.5), int(4/d + 0.5), nsr/500 }')"
    [ "${nsamp:-0}" -ge 1 ] || {
        echo "[dispatch] interval '$d' would need n_samp = 500*$d = ${nsamp:-0} raw samples per interval;" >&2
        echo "           the floor at this acquisition frequency is 0.002 (one sample)." >&2
        exit 1; }
    if [ "$(awk -v a="$d" -v b="$eff" 'BEGIN{print (a==b) ? 1 : 0}')" != 1 ]; then
        echo "[dispatch] interval '$d' is not a multiple of 1/500; running $eff (n_samp=$nsamp) instead." >&2
    fi
    # Two requested intervals can snap to the same reachable value (0.002 and 0.001 both give one raw
    # sample). That would put two entries with the same label on the axis, and the second would
    # silently overwrite the first in the output. Refuse instead.
    for prev in ${iv_lab[@]+"${iv_lab[@]}"}; do
        [ "$prev" = "$eff" ] && {
            echo "[dispatch] intervals '$d' and an earlier one both resolve to $eff (n_samp=$nsamp);" >&2
            echo "           duplicate axis labels would collide in the output. Drop one." >&2
            exit 1; }
    done
    iv_lab+=("$eff"); iv_samp+=("$nsamp"); iv_step1+=("$st1"); iv_step2+=("$st2")
done
echo "[dispatch] intervals: ${iv_lab[*]}  (n_samp: ${iv_samp[*]})"

axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= [%s])' "$(join_qcsv "${iv_lab[@]}")")
exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_step1[@]}")")
exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_samp[@]}")")
exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_step2[@]}")")
exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_samp[@]}")")
exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_step2[@]}")")
exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${iv_samp[@]}")")
nboot_arg=$( printf -- '--n_bootstrap_samples = get_number(n=%s)' "$N_BOOT")
mingrp_arg=$(printf -- '--min_groups_for_bootstrap = get_number(n=%s)' "$MIN_GROUPS")
gnmaxit_arg=$(printf -- '--gn_max_iter = get_number(n=%s)' "$GN_MAX_ITER")

for j in "${!N_ALGO[@]}"; do

for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"
    algo="${N_ALGO[$j]}"

    # Map the algorithm label to its (recursive, averaging, taylor, family,
    # variance_form) flags. family: 0=macro, 1=micro, 2=nonlinearsqr (LSE, which
    # has its own dispatcher). variance_form: 0=total, 1=residual.
    #
    # macro_VR is macro_MR's mean and gain (averaging=1, start-conditioned, no
    # boundary cross-covariance) with macro_IR's residual interval variance. It is
    # the control that splits the MR->IR gap into a variance step and a gain step.
    case "$algo" in
        macro_NR)  recursive=false; averaging=0 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_R)   recursive=true;  averaging=0 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_NMR) recursive=false; averaging=1 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_MR)  recursive=true;  averaging=1 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_VR)  recursive=true;  averaging=1 ; taylor=false ; family=0 ; variance_form=1 ;;
        macro_IR)  recursive=true;  averaging=2 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_IRT) recursive=true;  averaging=2 ; taylor=true  ; family=0 ; variance_form=0 ;;
        micro_R)   recursive=true;  averaging=0 ; taylor=false ; family=1 ; variance_form=0 ;;
        micro_MR)  recursive=true;  averaging=1 ; taylor=false ; family=1 ; variance_form=0 ;;
        micro_IR)  recursive=true;  averaging=2 ; taylor=false ; family=1 ; variance_form=0 ;;
        nonlinearsqr) recursive=false; averaging=1 ; taylor=false ; family=2 ; variance_form=0 ;;
        *) echo "[dispatch] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,VR,IR,IRT}, micro_{R,MR,IR}, or nonlinearsqr)" >&2; exit 1 ;;
    esac

    # Producer + file prefix by family. LSE (family 2) fixes i/Current_Noise (figure_4_LSE.macroir)
    # and keeps its own figure_3_LSE_ prefix, so it reuses the figure_2 LSE run and figure_2/4's
    # per-algo prefix logic; macro/micro use figure_4.macroir and the figure_3_G_ prefix.
    if [ "$family" = 2 ]; then SCRIPT="$SCRIPT_LSE";   prefix="figure_3_LSE_"
    else                       SCRIPT="$SCRIPT_MACRO"; prefix="figure_3_G_"; fi

    case "$nnoise" in                    # label -> current_noise (vnoise = label / 1000)
        0.05) vnoise=0.00005;;
        0.1)  vnoise=0.0001;;
        0.2)  vnoise=0.0002;;
        0.3)  vnoise=0.0003;;
        0.5)  vnoise=0.0005;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        1000)  vnoise=1;;
        10000) vnoise=10;;
        100000) vnoise=100;;
        *) echo "[dispatch] unknown noise level '$nnoise' (want 0.05, 0.1, 0.2, 0.3, 0.5, 1, 10, 100, 1000, 10000, 100000; vnoise = label/1000)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/%snch_%s_nsim_%s_%s_noise_%s"' "$prefix" "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    family_arg=$( printf -- '--algo_family_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$family")
    variance_form_arg=$( printf -- '--algo_variance_form_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$variance_form")

    # interval axis: built once above from INTERVALS, loop-invariant.

    # MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop so the
    # per-simulation / bootstrap loops become the active OpenMP level (load-bearing).
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        ${DEPEND:+--dependency="$DEP_SPEC"} \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-48G}" \
        --time="${TIME:-2-00:00:00}" \
        --job-name="f3G_${nch}c_${nsim}s" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",MACRODR_AXIS_SERIAL=1 \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$family_arg" "$variance_form_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$group_size_axis_arg" "$gs_arg" "$nboot_arg" "$mingrp_arg" "$gnmaxit_arg" \
        "$SCRIPT")

    echo "submitted fig4_nch_${nch}  job=${jobid}  n_sim=${nsim}  algo=${algo}  gs_axis=[${GROUP_SIZE[*]}]  -> ${WORKDIR}/figures/data/${prefix}nch_${nch}_nsim_${nsim}_${algo}_noise_${nnoise}_*"
done
done
