#!/bin/bash
# Focused dispatcher for the per-sample numerical Fisher Information bug hunt.
#
# Runs figure_2_per_sample_F.macroir on a SUBSET of replicas loaded from an
# existing _simulation.csv (saves minutes by skipping simulate()). Computes
# per-timestep F contribution for each loaded replica and writes a CSV with
# axes columns properly populated.
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/local/dispatch_per_sample_F.sh
# Tunables via env:
#   SIM_CSV     path to _simulation.csv to load replicas from
#                 (default: existing 64-replica run at nch=10000, macro_IR,
#                  noise=0.1, h_1e-4_1e-5_1e-6)
#   REPLICAS    space-separated list of replica indices to load
#                 (default: 0 1 5 10 44 — replica 44 is the sick one
#                  identified from the first replicates-CSV analysis)
#   H_REL       step-size for the central-diff Fisher  (default: 1e-5)
#   NCH         channels                                (default: 10000)
#   NOISE       noise_in_conductance_tau label          (default: 0.1)
#   ALGO        algorithm label                         (default: macro_IR)
#   BIN         macrodr_cli path
#   WORKDIR     output dir                              (default: repo's eLife_2025)
#   OMP_NUM_THREADS                                      (default: nproc)
#   MACRODR_AXIS_SERIAL                                   (default: 1)

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"
SCRIPT="$HERE/figure_2_per_sample_F.macroir"
REPO="$(cd "$HERE/../../../.." && pwd)"

BIN="${BIN:-$REPO/build/gcc-release/macrodr_cli}"
[ -x "$BIN" ] || {
    echo "[per_sample_F] binary not found: $BIN" >&2
    echo "        build: cmake --build --preset gcc-release -j2" >&2
    exit 1
}

# Defaults: point at the existing run with the sick replica 44.
SIM_CSV="${SIM_CSV:-$REPO/projects/eLife_2025/figures/data/figure_2_nch_10000_nsim_64_macro_IR_noise_0.1_h_1e-4_1e-5_1e-6_simulation.csv}"
REPLICAS=(${REPLICAS:-0 1 5 10 44})
H_REL="${H_REL:-1e-5}"
NCH="${NCH:-10000}"
NOISE="${NOISE:-0.1}"
ALGO="${ALGO:-macro_IR}"

[ -f "$SIM_CSV" ] || {
    echo "[per_sample_F] simulation CSV not found: $SIM_CSV" >&2
    echo "        run figure_2 first or override SIM_CSV=..." >&2
    exit 1
}

# Map algo label → (recursive, averaging, taylor, micro) — same case as the
# main figure_2 dispatcher.
case "$ALGO" in
    macro_NR)  recursive=false; averaging=0; taylor=false; micro=false ;;
    macro_R)   recursive=true;  averaging=0; taylor=false; micro=false ;;
    macro_NMR) recursive=false; averaging=1; taylor=false; micro=false ;;
    macro_MR)  recursive=true;  averaging=1; taylor=false; micro=false ;;
    macro_IR)  recursive=true;  averaging=2; taylor=false; micro=false ;;
    macro_IRT) recursive=true;  averaging=2; taylor=true;  micro=false ;;
    micro_R)   recursive=true;  averaging=0; taylor=false; micro=true ;;
    micro_MR)  recursive=true;  averaging=1; taylor=false; micro=true ;;
    micro_IR)  recursive=true;  averaging=2; taylor=false; micro=true ;;
    *) echo "[per_sample_F] unknown algo '$ALGO'" >&2; exit 1 ;;
esac

case "$NOISE" in
    0.1)  vnoise=0.0001 ;;
    1)    vnoise=0.001  ;;
    10)   vnoise=0.01   ;;
    100)  vnoise=0.1    ;;
    *) echo "[per_sample_F] unknown noise '$NOISE'" >&2; exit 1 ;;
esac

WORKDIR="${WORKDIR:-$REPO/projects/eLife_2025}"
mkdir -p "$WORKDIR/figures/data"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export MACRODR_AXIS_SERIAL="${MACRODR_AXIS_SERIAL:-1}"

# Build the replica_indices DSL expression: [0,1,5,10,44]
REPLICAS_DSL=$(IFS=,; echo "[${REPLICAS[*]}]")

# Build the injection arguments.
axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$NCH")
num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$NCH")
fp_arg=$(  printf -- '--filepath = "figures/data/figure_2_per_sample_F_nch_%s_%s_noise_%s_h_%s_replicas_%s"' \
            "$NCH" "$ALGO" "$NOISE" "$H_REL" "$(IFS=_; echo "${REPLICAS[*]}")")
axis_noise_arg=$(   printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$NOISE")
current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
axis_algo_arg=$(printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$ALGO")
recursive_arg=$(printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
averaging_arg=$(printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
taylor_arg=$(   printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
micro_arg=$(    printf -- '--algo_micro_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$micro")
axis_h_arg=$(   printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= ["%s"])' "$H_REL")
h_rel_arg=$(    printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$H_REL")
sim_csv_arg=$(  printf -- '--simulation_csv_path = "%s"' "$SIM_CSV")
replicas_arg=$( printf -- '--replica_indices = %s' "$REPLICAS_DSL")

cd "$WORKDIR"
echo "[per_sample_F] bin=$BIN"
echo "[per_sample_F] cwd=$WORKDIR  threads=$OMP_NUM_THREADS  axis_serial=$MACRODR_AXIS_SERIAL"
echo "[per_sample_F] sim_csv=$SIM_CSV"
echo "[per_sample_F] replicas=${REPLICAS[*]}  h_rel=$H_REL  algo=$ALGO  noise=$NOISE  nch=$NCH"
echo

"$BIN" "$axis_arg" "$num_arg" "$fp_arg" \
       "$axis_noise_arg" "$current_noise_arg" \
       "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" "$micro_arg" \
       "$axis_h_arg" "$h_rel_arg" \
       "$sim_csv_arg" "$replicas_arg" \
       "$SCRIPT"

echo
echo "[per_sample_F] done. output:"
ls -1 "$WORKDIR"/figures/data/figure_2_per_sample_F_* 2>/dev/null || true
