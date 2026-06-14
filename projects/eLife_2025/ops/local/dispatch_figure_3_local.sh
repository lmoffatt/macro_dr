#!/bin/bash
# Local dispatcher for figure_3_mle — same injection contract as the dirac
# dispatcher (ops/slurm/dispatch_figure_3.sh), but runs the binary DIRECTLY on
# this machine: NO sbatch, NO cluster profile, NO run_macroir.sh wrapper. Jobs
# run SEQUENTIALLY in a loop. Defaults to a SMALL grid for end-to-end validation;
# scale up via the same env knobs as the dirac one.
#
# FILE CONTRACT — figure_3_mle.macroir must NOT define the injected names itself
# (same as the dirac dispatcher): axis_Nchanels, Num_ch, axis_noise, current_noise,
# n_simulations, filepath, algorithm_axis, algo_*_approximation, axis_interval,
# exp_n_step_/exp_n_samp_, axis_h_fim, h_rel_value, group_size,
# n_bootstrap_samples, min_groups_for_bootstrap, gn_max_iter.
#
# Prereq: a local build (the user compiles): build/gcc-release/macrodr_cli.
#
# Usage (from anywhere; paths are resolved relative to the repo base):
#   projects/eLife_2025/ops/local/dispatch_figure_3_local.sh
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), N_NOISE,
# N_ALGO, H_RELS, GROUP_SIZE, N_BOOT, MIN_GROUPS, GN_MAX_ITER, THREADS, BIN, WORKDIR.
# Example (reproduce a dirac cell locally):
#   NCHS="1000" N_SIMS="64" N_ALGO="macro_IR" \
#     projects/eLife_2025/ops/local/dispatch_figure_3_local.sh

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../eLife_2025/ops/local
BASE="$(readlink -f "$HERE/../../../..")"           # repo base (macro_dr)
SCRIPT="$(readlink -f "$HERE/figure_3_mle.macroir")"

# Local binary (the gcc-release build the user compiles). Override with BIN=…
BIN="${BIN:-$(readlink -f "$BASE/build/gcc-release/macrodr_cli")}"
[ -x "$BIN" ] || {
    echo "[local] binary not found: $BIN" >&2
    echo "        build first:  cmake --build --preset gcc-release" >&2
    exit 1
}

# SMALL validation grid by default (the joint MLE + two batteries per cell are the
# cost). NCHS and N_SIMS are parallel arrays paired by index.
NCHS=(${NCHS:-1000})
N_SIMS=(${N_SIMS:-32})
N_NOISE=(${N_NOISE:-0.1})
N_ALGO=(${N_ALGO:-macro_IR})
H_RELS=(${H_RELS:-1e-5})

# Per-group MLE knobs (same as dirac; injected as size_t literals).
GROUP_SIZE="${GROUP_SIZE:-1}"          # 1 = per-replicate; >1 pools recordings per refit
N_BOOT="${N_BOOT:-100}"                # bootstrap replicates over groups
MIN_GROUPS="${MIN_GROUPS:-10}"         # below this, probit slots NaN-filled
GN_MAX_ITER="${GN_MAX_ITER:-100}"      # GN per-group refit iteration cap

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[local] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}

# Output dir: the repo's eLife_2025 tree, so figures/data/ lands next to the
# dirac downloads and the R notebooks read them. Override with WORKDIR=…
WORKDIR="${WORKDIR:-$BASE/projects/eLife_2025}"
WORKDIR="$(readlink -f "$WORKDIR")"
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"

# Threads: OpenMP across the per-simulation / bootstrap loop. BLAS stays single-
# threaded (parallelism lives at the OMP level; threaded BLAS would oversubscribe).
# MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop (load-bearing —
# without it memory grows with concurrent combos -> OOM).
export OMP_NUM_THREADS="${THREADS:-$(nproc)}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${BLAS_THREADS:-1}"
export BLIS_NUM_THREADS="${BLAS_THREADS:-1}"
export MACRODR_AXIS_SERIAL=1

join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }

# Loop-invariant injections (built once). axis_h_fim must precede h_rel_value.
axis_h_arg=$(printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= [%s])' "$(join_qcsv "${H_RELS[@]}")")
h_rel_arg=$( printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$(join_csv "${H_RELS[@]}")")
gs_arg=$(    printf -- '--group_size = get_number(n=%s)' "$GROUP_SIZE")
nboot_arg=$( printf -- '--n_bootstrap_samples = get_number(n=%s)' "$N_BOOT")
mingrp_arg=$(printf -- '--min_groups_for_bootstrap = get_number(n=%s)' "$MIN_GROUPS")
gnmaxit_arg=$(printf -- '--gn_max_iter = get_number(n=%s)' "$GN_MAX_ITER")

cd "$WORKDIR"
echo "[local] base=$BASE"
echo "[local] bin=$BIN"
echo "[local] cwd=$WORKDIR  threads=$OMP_NUM_THREADS"

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
        *) echo "[local] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,IR,IRT})" >&2; exit 1 ;;
    esac

    case "$nnoise" in
        0.1)  vnoise=0.0001;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[local] unknown noise level '$nnoise' (want 0.1, 1, 10, 100)" >&2; exit 1 ;;
    esac

    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_3_nch_%s_nsim_%s_%s_noise_%s_gs_%s"' "$nch" "$nsim" "$algo" "$nnoise" "$GROUP_SIZE")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    micro_arg=$( printf -- '--algo_micro_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$micro")

    # interval_in_tau grid (identical to figure_2/dirac). axis_interval precedes exp_n_*.
    axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= ["1","0.5","0.2","0.1","0.05","0.02","0.01"])')
    exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2,4,10,20,40,100,200])')
    exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')

    fp_label="figure_3_nch_${nch}_nsim_${nsim}_${algo}_noise_${nnoise}_gs_${GROUP_SIZE}"
    log="$WORKDIR/logs/${fp_label}.log"
    echo "[local] running ${fp_label}  (n_sim=${nsim}, algo=${algo}, gs=${GROUP_SIZE}) -> ${log}"

    "$BIN" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$micro_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$axis_h_arg" "$h_rel_arg" \
        "$gs_arg" "$nboot_arg" "$mingrp_arg" "$gnmaxit_arg" \
        "$SCRIPT" 2>&1 | tee "$log"

    echo "[local] done ${fp_label} -> ${WORKDIR}/figures/data/${fp_label}_*"
done
done

echo "[local] all cells finished."
