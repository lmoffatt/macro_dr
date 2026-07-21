#!/bin/bash
# Local (no-SLURM) twin of ops/slurm/dispatch_figure_3_G.sh — runs figure_3_mle_G
# directly on your machine so you can smoke-test the pipeline (in particular the
# newly enabled micro_R / micro_MR / micro_IR MLE path) without a cluster.
#
# Same injection set and algorithm->flags mapping as the SLURM dispatcher, but:
#   - no sbatch / cluster profile: the binary is exec'd directly;
#   - small, fast DEFAULTS (Nch=10, 16 sims, micro_R, group_size=1) — a smoke,
#     not a production grid;
#   - output lands in a local run dir you can inspect, not scratch.
#
# Prereq: a local build. Point BIN at it or let it auto-detect
#   build/gcc-release/macrodr_cli (falls back to gcc-debug).
#
# Usage (from the repo root):
#   projects/eLife_2025/ops/local/run_figure_3_G_local.sh
#   N_ALGO="micro_R macro_R" projects/eLife_2025/ops/local/run_figure_3_G_local.sh
#   BIN=build/gcc-debug/macrodr_cli N_SIMS=8 projects/eLife_2025/ops/local/run_figure_3_G_local.sh
#
# Tunables via env (same names as the dispatcher where they overlap):
#   BIN, WORKDIR, NCHS, N_SIMS (parallel arrays), N_NOISE, N_ALGO, GROUP_SIZE,
#   N_BOOT, MIN_GROUPS, GN_MAX_ITER, OMP_NUM_THREADS.

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/local
REPO_ROOT="$(readlink -f "$HERE/../../../..")"    # repo base
SCRIPT="$HERE/figure_3_mle_G.macroir"
[ -f "$SCRIPT" ] || { echo "[local] script not found: $SCRIPT" >&2; exit 1; }

# ---- Binary: explicit BIN wins; else auto-detect a local build ---------------
if [ -z "${BIN:-}" ]; then
    for cand in "$REPO_ROOT/build/gcc-release/macrodr_cli" \
                "$REPO_ROOT/build/gcc-release/bin/macrodr_cli" \
                "$REPO_ROOT/build/gcc-debug/macrodr_cli" \
                "$REPO_ROOT/build/gcc-debug/bin/macrodr_cli"; do
        [ -x "$cand" ] && { BIN="$cand"; break; }
    done
fi
BIN="$(readlink -f "${BIN:-/nonexistent}")"
[ -x "$BIN" ] || {
    echo "[local] macrodr_cli not found/executable: ${BIN}" >&2
    echo "        build it first, or pass BIN=/path/to/macrodr_cli" >&2
    exit 1
}

# ---- Output dir (binary writes figures/data/... relative to cwd) --------------
WORKDIR="$(readlink -f "${WORKDIR:-$HERE/run_local}")"
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"

# ---- Threading: mirror the cluster model (per-sim/bootstrap OMP; serial BLAS) -
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${BLAS_THREADS:-1}"
export BLIS_NUM_THREADS="${BLAS_THREADS:-1}"
export MACRODR_AXIS_SERIAL=1

# ---- The (small by default) grid. NCHS and N_SIMS are parallel arrays ---------
NCHS=(${NCHS:-10})
N_SIMS=(${N_SIMS:-16})
N_NOISE=(${N_NOISE:-0.1})
N_ALGO=(${N_ALGO:-micro_R})

GROUP_SIZE=(${GROUP_SIZE:-1})       # keep <= min(N_SIMS): N_groups = nsim/gsize
N_BOOT="${N_BOOT:-20}"
MIN_GROUPS="${MIN_GROUPS:-10}"
GN_MAX_ITER="${GN_MAX_ITER:-100}"

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[local] NCHS (${#NCHS[@]}) and N_SIMS (${#N_SIMS[@]}) must be same length" >&2
    exit 1
}

# Stage-1 groups the recordings by each group_size in the axis: N_groups = nsim/gsize.
# A gsize larger than a job's nsim yields 0 groups and the binary aborts that combo.
# Catch it here (across all jobs) with a clear message instead of a per-combo failure.
min_nsim=2147483647
for s in "${N_SIMS[@]}"; do (( s < min_nsim )) && min_nsim=$s; done
max_gsize=0
for gs in "${GROUP_SIZE[@]}"; do (( gs > max_gsize )) && max_gsize=$gs; done
if (( max_gsize > min_nsim )); then
    echo "[local] group_size=${max_gsize} exceeds the smallest N_SIMS=${min_nsim}: no group forms" >&2
    echo "        (N_groups = nsim/gsize). Use GROUP_SIZE values <= your smallest N_SIMS" >&2
    echo "        (e.g. GROUP_SIZE=\"1\" for the smoke), or raise N_SIMS." >&2
    exit 1
fi

echo "[local] bin=${BIN}"
echo "[local] workdir=${WORKDIR}   threads=${OMP_NUM_THREADS}"
echo "[local] grid: algos=[${N_ALGO[*]}]  nch=[${NCHS[*]}]  nsim=[${N_SIMS[*]}]  noise=[${N_NOISE[*]}]  gsize=[${GROUP_SIZE[*]}]"

join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }

# Per-group MLE knobs (size_t literals). group_size is a DSL axis.
group_size_axis_arg=$(printf -- '--group_size_axis = axis(name= "group_size", labels= [%s])' "$(join_qcsv "${GROUP_SIZE[@]}")")
gs_arg=$( printf -- '--group_size = indexed_size_by(axis= group_size_axis, values=[%s])' "$(join_csv "${GROUP_SIZE[@]}")")
nboot_arg=$( printf -- '--n_bootstrap_samples = get_number(n=%s)' "$N_BOOT")
mingrp_arg=$(printf -- '--min_groups_for_bootstrap = get_number(n=%s)' "$MIN_GROUPS")
gnmaxit_arg=$(printf -- '--gn_max_iter = get_number(n=%s)' "$GN_MAX_ITER")

cd "$WORKDIR"

for j in "${!N_ALGO[@]}"; do
for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"
    algo="${N_ALGO[$j]}"

    # Algorithm label -> (recursive, averaging, taylor, family, variance_form) flags.
    # family: 0=macro, 1=micro. variance_form: 0=total, 1=residual.
    # macro_VR = MR's mean and gain (averaging=1) with IR's residual variance.
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
        *) echo "[local] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,VR,IR,IRT} or micro_{R,MR,IR})" >&2; exit 1 ;;
    esac

    case "$nnoise" in                    # label -> current_noise (vnoise = label / 1000)
        0.05) vnoise=0.00005;;
        0.1)  vnoise=0.0001;;
        0.2)  vnoise=0.0002;;
        0.3)  vnoise=0.0003;;
        0.5)  vnoise=0.0005;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[local] unknown noise level '$nnoise'" >&2; exit 1 ;;
    esac

    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_3_G_nch_%s_nsim_%s_%s_noise_%s"' "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    family_arg=$( printf -- '--algo_family_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$family")
    variance_form_arg=$( printf -- '--algo_variance_form_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$variance_form")

    axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= ["1","0.5","0.2","0.1","0.05","0.02","0.01"])')
    exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2,4,10,20,40,100,200])')
    exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')

    log="$WORKDIR/logs/fig3G_${algo}_nch_${nch}_nsim_${nsim}_noise_${nnoise}.log"
    # CHECK_ONLY=1 -> parse the script+injections and exit (seconds, no compute).
    check_flag=(); [ -n "${CHECK_ONLY:-}" ] && check_flag=(-n)
    if [ -n "${CHECK_ONLY:-}" ]; then
        echo "[local] CHECK-SYNTAX only: algo=${algo} nch=${nch}"
    else
        echo "[local] running algo=${algo} nch=${nch} nsim=${nsim} noise=${nnoise}  -> ${log}"
    fi

    set +e
    "$BIN" "${check_flag[@]}" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$family_arg" "$variance_form_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$group_size_axis_arg" "$gs_arg" "$nboot_arg" "$mingrp_arg" "$gnmaxit_arg" \
        "$SCRIPT" 2>&1 | tee "$log"
    rc="${PIPESTATUS[0]}"
    set -e

    if [ "$rc" -ne 0 ]; then
        echo "[local] FAILED (exit $rc) algo=${algo} nch=${nch} — see ${log}" >&2
        exit "$rc"
    fi
    echo "[local] OK algo=${algo} nch=${nch}; CSVs:"
    ls -1 "$WORKDIR"/figures/data/figure_3_G_nch_"${nch}"_nsim_"${nsim}"_"${algo}"_noise_"${nnoise}"* 2>/dev/null | sed 's/^/    /' || true
done
done

echo "[local] done. outputs under ${WORKDIR}/figures/data/"
