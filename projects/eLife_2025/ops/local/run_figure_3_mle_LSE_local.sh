#!/bin/bash
# Local one-point validation of figure_3_mle_LSE.macroir (no SLURM). Mirrors the
# injection idiom of dispatch_figure_3_LSE.sh but runs the release binary directly
# for a SINGLE grid cell, so we can eyeball the θ̂ cloud + the (post-fix) param_name
# labelling before committing the full grid on dirac.
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/local/run_figure_3_mle_LSE_local.sh
# Override the cell via env: NCH, NSIM, NOISE, VNOISE, INTERVAL_*, GROUP_SIZE, BIN.
set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"
SCRIPT="$HERE/figure_3_mle_LSE.macroir"
BIN="${BIN:-$(readlink -f build/gcc-release/macrodr_cli)}"
[ -x "$BIN" ] || { echo "[run] binary not found: $BIN (build gcc-release first)" >&2; exit 1; }

# One grid cell (small, for a quick local check).
NCH="${NCH:-100}"
NSIM="${NSIM:-64}"
NOISE="${NOISE:-0.1}"          # label
VNOISE="${VNOISE:-0.0001}"     # current_noise = label/1000
GROUP_SIZE="${GROUP_SIZE:-1}"
N_BOOT="${N_BOOT:-100}"
MIN_GROUPS="${MIN_GROUPS:-10}"
GN_MAX_ITER="${GN_MAX_ITER:-100}"
# Single interval_in_tau = 1 (exp_n_* below match the dispatch's first column).
FP="${FP:-figures/data/figure_3_LSE_local_nch_${NCH}_${NOISE}}"

mkdir -p figures/data

# Injections in dependency order (an axis must precede the indexed_*_by that uses it).
"$BIN" \
  "--axis_Nchanels = axis(name= \"Num_ch\", labels= [\"$NCH\"])" \
  "--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[$NCH])" \
  "--n_simulations = get_number(n=$NSIM)" \
  "--filepath = \"$FP\"" \
  "--axis_noise = axis(name= \"noise_in_conductance_tau\", labels= [\"$NOISE\"])" \
  "--current_noise = indexed_double_by(axis= axis_noise, values=[$VNOISE])" \
  "--algorithm_axis = axis(name= \"algorithm\", labels= [\"nonlinearsqr\"])" \
  "--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[false])" \
  "--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[1])" \
  "--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[false])" \
  "--algo_family_approximation = indexed_int_by(axis= algorithm_axis, values=[2])" \
  "--axis_interval = axis(name= \"interval_in_tau\", labels= [\"1\"])" \
  "--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2])" \
  "--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500])" \
  "--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4])" \
  "--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500])" \
  "--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4])" \
  "--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500])" \
  "--group_size_axis = axis(name= \"group_size\", labels= [\"$GROUP_SIZE\"])" \
  "--group_size = indexed_size_by(axis= group_size_axis, values=[$GROUP_SIZE])" \
  "--n_bootstrap_samples = get_number(n=$N_BOOT)" \
  "--min_groups_for_bootstrap = get_number(n=$MIN_GROUPS)" \
  "--gn_max_iter = get_number(n=$GN_MAX_ITER)" \
  "$SCRIPT"

echo "[run] done -> ${FP}_mle_cloud_runs.csv  (+ _pool_runs.csv)"
