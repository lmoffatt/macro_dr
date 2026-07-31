#!/bin/bash
# LOCAL (no-SLURM) version of dispatch_figure_2.sh — runs sequentially on this
# machine with the locally-built binary. Mirrors the SLURM dispatcher's
# injections (axis_Nchanels/Num_ch/n_simulations/filepath/axis_noise/
# current_noise/algorithm_axis/algo_*_approximation) plus an extra h_rel sweep
# axis for the central-diff Fisher step-size scan.
#
# Each (algo, nch, nsim, noise, h_rel) combination produces its own CSV under
# WORKDIR/figures/data with the filename encoding all five tags.
#
# Prereq: a local build that includes the new write_csv overload taking
# (dlikelihood_predictions, numerical_fisher_information, path):
#   cmake --preset gcc-release && cmake --build --preset gcc-release
#
# Usage (from anywhere — paths are resolved absolutely):
#   projects/eLife_2025/ops/local/dispatch_figure_2_local.sh
# Tunables via env:
#   NCHS       parallel array, job i = NCHS[i] channels        (default: 10000)
#   N_SIMS     parallel array with NCHS                         (default: 16384)
#   N_NOISE    parallel array with NCHS                         (default: 0.1)
#   N_ALGO     algorithm list (Cartesian product against NCHS)  (default: macro_IR)
#   H_RELS     h_rel step-size sweep (Cartesian against the above)
#                                                               (default: 1e-5)
#   BIN        path to macrodr_cli  (default: build/gcc-release/macrodr_cli)
#   WORKDIR    output dir           (default: /tmp/macro_dr/eLife_2025)
#   OMP_NUM_THREADS                 (default: nproc)
# Example (step-size scan for the bug hunt at 16k samples, IR, noise 0.1):
#   H_RELS="1e-5 3e-6 3e-5" NCHS=10000 N_SIMS=16384 N_NOISE=0.1 N_ALGO=macro_IR \
#     projects/eLife_2025/ops/local/dispatch_figure_2_local.sh

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/local
SCRIPT="$HERE/figure_2.macroir"
REPO="$(cd "$HERE/../../../.." && pwd)"           # repo base

BIN="${BIN:-$REPO/build/gcc-release/macrodr_cli}"
[ -x "$BIN" ] || {
    echo "[local] binary not found: $BIN" >&2
    echo "        build: cmake --preset gcc-release && cmake --build --preset gcc-release" >&2
    echo "        or set BIN=/path/to/macrodr_cli" >&2
    exit 1
}

NCHS=(${NCHS:-10000})
N_SIMS=(${N_SIMS:-4096})
N_NOISE=(${N_NOISE:-0.1})
N_ALGO=(${N_ALGO:-macro_IR})
H_RELS=(${H_RELS:-1e-4 1e-5 1e-6})

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[local] NCHS (${#NCHS[@]}) and N_SIMS (${#N_SIMS[@]}) must be same length" >&2
    exit 1
}
[ "${#NCHS[@]}" -eq "${#N_NOISE[@]}" ] || {
    echo "[local] NCHS (${#NCHS[@]}) and N_NOISE (${#N_NOISE[@]}) must be same length" >&2
    exit 1
}

WORKDIR="${WORKDIR:-$REPO/projects/eLife_2025}"
mkdir -p "$WORKDIR/figures/data"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
# Serialize the axis-combo loop so the inner per-simulation parallel-for owns
# the OpenMP threads. Without this the combo loop runs in parallel, holding
# per-combo state concurrently — memory ∝ #combos × per-sim footprint, which
# OOMs at production nsim (cf. SLURM jobs 110042-110045 killed at ~129 GB).
# Default ON; override with MACRODR_AXIS_SERIAL=0 to benchmark combo-parallel.
export MACRODR_AXIS_SERIAL="${MACRODR_AXIS_SERIAL:-1}"

echo "[local] bin=$BIN"
echo "[local] workdir=$WORKDIR  threads=$OMP_NUM_THREADS  axis_serial=$MACRODR_AXIS_SERIAL"
echo "[local] grid: nch=${NCHS[*]}  nsim=${N_SIMS[*]}  noise=${N_NOISE[*]}  algo=${N_ALGO[*]}  h_rel=${H_RELS[*]}  intervals=${INTERVALS[*]:-0.01}"
cd "$WORKDIR"

# h_rel goes as a multi-value axis INSIDE each (algo, nch) job — all H_RELS in
# one CSV per job. Run the script multiple times with different H_RELS for
# per-h_rel files. Axis args don't depend on the loop variables, so build once.
join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }
join_tag()  { local IFS=_; echo "$*"; }
axis_h_arg=$( printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= [%s])' "$(join_qcsv "${H_RELS[@]}")")
h_rel_arg=$(  printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$(join_csv "${H_RELS[@]}")")
h_tag=$(join_tag "${H_RELS[@]}")

# ---- interval_in_tau injection -------------------------------------------
# figure_2.macroir now has axis_interval + the experiment step/sample sizes
# injection-ready (commented), so the dispatcher controls which intervals run.
# INTERVALS selects a subset of the 7 production intervals; per-interval step /
# sample arrays are sliced consistently. Default = "0.01" only (the interval
# where the numeric-FIM bug shows up). For the full production grid pass
# INTERVALS="1 0.5 0.2 0.1 0.05 0.02 0.01".
INTERVALS=(${INTERVALS:-0.01})

# Full 7-interval reference arrays (parallel to ALL_INTERVAL_LABELS). These are
# the exact values that used to be hardcoded in figure_2.macroir.
ALL_INTERVAL_LABELS=(1 0.5 0.2 0.1 0.05 0.02 0.01)
ALL_STEP_1=(2 4 10 20 40 100 200)      ; ALL_SAMP_1=(500 250 100 50 25 10 5)
ALL_STEP_2=(4 8 20 40 80 200 400)      ; ALL_SAMP_2=(500 250 100 50 25 10 5)
ALL_STEP_3=(4 8 20 40 80 200 400)      ; ALL_SAMP_3=(500 250 100 50 25 10 5)

interval_index() {
    local target="$1" k
    for k in "${!ALL_INTERVAL_LABELS[@]}"; do
        [ "${ALL_INTERVAL_LABELS[$k]}" = "$target" ] && { echo "$k"; return 0; }
    done
    echo "[local] unknown interval '$target' (want one of: ${ALL_INTERVAL_LABELS[*]})" >&2
    return 1
}

sel_labels=(); sel_step_1=(); sel_samp_1=(); sel_step_2=(); sel_samp_2=(); sel_step_3=(); sel_samp_3=()
for lbl in "${INTERVALS[@]}"; do
    idx=$(interval_index "$lbl") || exit 1
    sel_labels+=("$lbl")
    sel_step_1+=("${ALL_STEP_1[$idx]}"); sel_samp_1+=("${ALL_SAMP_1[$idx]}")
    sel_step_2+=("${ALL_STEP_2[$idx]}"); sel_samp_2+=("${ALL_SAMP_2[$idx]}")
    sel_step_3+=("${ALL_STEP_3[$idx]}"); sel_samp_3+=("${ALL_SAMP_3[$idx]}")
done

# axis_interval MUST be injected before the exp_n_* (which reference it).
axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= [%s])' "$(join_qcsv "${sel_labels[@]}")")
exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_step_1[@]}")")
exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_samp_1[@]}")")
exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_step_2[@]}")")
exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_samp_2[@]}")")
exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_step_3[@]}")")
exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[%s])' "$(join_csv "${sel_samp_3[@]}")")
# --------------------------------------------------------------------------

OUTS=()  # base paths produced by this run (for the final listing)
for j in "${!N_ALGO[@]}"; do
for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"
    algo="${N_ALGO[$j]}"

    # Map the algorithm label to its (recursive, averaging, taylor, micro) flags.
    # Mirrors the SLURM dispatcher's case statement.
    case "$algo" in
        macro_NR)  recursive=false; averaging=0; taylor=false; micro=false ;;
        macro_R)   recursive=true;  averaging=0; taylor=false; micro=false ;;
        macro_INR) recursive=false; averaging=1; taylor=false; micro=false ;;
        macro_MR)  recursive=true;  averaging=1; taylor=false; micro=false ;;
        macro_IR)  recursive=true;  averaging=2; taylor=false; micro=false ;;
        macro_IRT) recursive=true;  averaging=2; taylor=true;  micro=false ;;
        micro_R)   recursive=true;  averaging=0; taylor=false; micro=true ;;
        micro_MR)  recursive=true;  averaging=1; taylor=false; micro=true ;;
        micro_IR)  recursive=true;  averaging=2; taylor=false; micro=true ;;
        *) echo "[local] unknown algorithm '$algo' (want macro_{NR,R,INR,MR,IR,IRT} or micro_{R,MR,IR})" >&2; exit 1 ;;
    esac

    case "$nnoise" in
        0.1)  vnoise=0.0001;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[local] unknown noise level '$nnoise' (want 0.1, 1, 10, 100)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_base="figures/data/figure_2_nch_${nch}_nsim_${nsim}_${algo}_noise_${nnoise}_h_${h_tag}"
    OUTS+=("$fp_base")
    fp_arg=$(  printf -- '--filepath = "%s"' "$fp_base")
    axis_noise_arg=$(   printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$(    printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    micro_arg=$(     printf -- '--algo_micro_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$micro")

    echo
    echo "=== nch=${nch}  nsim=${nsim}  algo=${algo}  noise=${nnoise}  h_rel=${H_RELS[*]} ==="
    "$BIN" "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
           "$axis_noise_arg" "$current_noise_arg" \
           "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" "$micro_arg" \
           "$axis_interval_arg" \
           "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
           "$exp_step_3_arg" "$exp_samp_3_arg" \
           "$axis_h_arg" "$h_rel_arg" \
           "$SCRIPT"
done
done

echo
echo "[local] done. outputs (this run only):"
for base in "${OUTS[@]}"; do
    ls -1 "$WORKDIR/$base"* 2>/dev/null || true
done
