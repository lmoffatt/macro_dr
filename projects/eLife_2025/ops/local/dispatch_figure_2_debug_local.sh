#!/bin/bash
# LOCAL DEBUG variant — runs figure_2_debug.macroir instead of figure_2.macroir.
# Same dispatcher injections as dispatch_figure_2_local.sh (axes for nch /
# nsim / noise / algo / h_rel), but the .macroir target adds the per-sample
# numerical Fisher calculation + write_csv on top of the standard pipeline.
# Used for the bug hunt on FD instabilities in score derivatives — produces
# an extra _per_sample_F CSV per (algo, nch) combo.
#
# Keep nsim small (8-64) — _per_sample_F CSV scales as
# n_replicas × n_timesteps × p² and grows fast.
#
# Prereq: same local build as dispatch_figure_2_local.sh + the
# calc_per_sample_numerical_fisher_information overload.
#
# Usage (from anywhere — paths are resolved absolutely):
#   projects/eLife_2025/ops/local/dispatch_figure_2_debug_local.sh
# Tunables: same as dispatch_figure_2_local.sh (NCHS, N_SIMS, N_NOISE,
# N_ALGO, H_RELS, BIN, WORKDIR, OMP_NUM_THREADS, MACRODR_AXIS_SERIAL).
# Example (smoke-run at nsim=8 with h_rel=1e-5, macro_IR, nch=10000):
#   N_SIMS=8 H_RELS="1e-5" NCHS=10000 N_ALGO=macro_IR \
#     projects/eLife_2025/ops/local/dispatch_figure_2_debug_local.sh

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/local
SCRIPT="$HERE/figure_2_debug.macroir"
REPO="$(cd "$HERE/../../../.." && pwd)"           # repo base

BIN="${BIN:-$REPO/build/gcc-release/macrodr_cli}"
[ -x "$BIN" ] || {
    echo "[local] binary not found: $BIN" >&2
    echo "        build: cmake --preset gcc-release && cmake --build --preset gcc-release" >&2
    echo "        or set BIN=/path/to/macrodr_cli" >&2
    exit 1
}

NCHS=(${NCHS:-10000})
N_SIMS=(${N_SIMS:-64})
N_NOISE=(${N_NOISE:-0.1})
N_ALGO=(${N_ALGO:-macro_IR})
H_RELS=(${H_RELS:-1e-6})

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
echo "[local] grid: nch=${NCHS[*]}  nsim=${N_SIMS[*]}  noise=${N_NOISE[*]}  algo=${N_ALGO[*]}  h_rel=${H_RELS[*]}"
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
        macro_NMR) recursive=false; averaging=1; taylor=false; micro=false ;;
        macro_MR)  recursive=true;  averaging=1; taylor=false; micro=false ;;
        macro_IR)  recursive=true;  averaging=2; taylor=false; micro=false ;;
        macro_IRT) recursive=true;  averaging=2; taylor=true;  micro=false ;;
        micro_R)   recursive=true;  averaging=0; taylor=false; micro=true ;;
        micro_MR)  recursive=true;  averaging=1; taylor=false; micro=true ;;
        micro_IR)  recursive=true;  averaging=2; taylor=false; micro=true ;;
        *) echo "[local] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,IR,IRT} or micro_{R,MR,IR})" >&2; exit 1 ;;
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
           "$axis_h_arg" "$h_rel_arg" \
           "$SCRIPT"
done
done

echo
echo "[local] done. outputs (this run only):"
for base in "${OUTS[@]}"; do
    ls -1 "$WORKDIR/$base"* 2>/dev/null || true
done
