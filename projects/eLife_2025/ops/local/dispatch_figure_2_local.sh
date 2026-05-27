#!/bin/bash
# LOCAL (no-SLURM) version of dispatch_figure_2.sh — runs the Nch grid
# sequentially on this machine with the locally-built binary. Meant for testing
# the plumbing (injection → concat → "# binary=" stamp → CSV output) before
# spending cluster time; use tiny NCHS / N_SIMS so it finishes in seconds.
#
# Prereq: a local build that includes concat + the CSV hash stamp, e.g.
#   cmake --preset gcc-release && cmake --build --preset gcc-release
#
# Usage (from anywhere — paths are resolved absolutely):
#   projects/eLife_2025/ops/local/dispatch_figure_2_local.sh
# Tunables via env:
#   NCHS, N_SIMS   parallel arrays, same length (job i = NCHS[i] ch, N_SIMS[i] sims)
#   BIN            path to macrodr_cli  (default: build/gcc-release/macrodr_cli)
#   WORKDIR        output dir           (default: /tmp/macro_dr/eLife_2025)
#   OMP_NUM_THREADS                     (default: nproc)
# Example:
#   NCHS="10 100 1000" N_SIMS="32 16 8" dispatch_figure_2_local.sh

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/local
SCRIPT="$HERE/figure_2.macroir"
REPO="$(cd "$HERE/../../../.." && pwd)"           # repo base

BIN="${BIN:-$REPO/build/gcc-release/macrodr_cli}"
[ -x "$BIN" ] || {
    echo "[local] binary not found: $BIN" >&2
    echo "        build it: cmake --preset gcc-release && cmake --build --preset gcc-release" >&2
    echo "        or set BIN=/path/to/macrodr_cli" >&2
    exit 1
}

NCHS=(${NCHS:-10 100})
N_SIMS=(${N_SIMS:-16 16})
[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[local] NCHS (${#NCHS[@]}) and N_SIMS (${#N_SIMS[@]}) must be the same length" >&2
    exit 1
}

WORKDIR="${WORKDIR:-/tmp/macro_dr/eLife_2025}"
mkdir -p "$WORKDIR/figures/data"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"

echo "[local] bin=$BIN"
echo "[local] workdir=$WORKDIR  threads=$OMP_NUM_THREADS"
cd "$WORKDIR"

for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = %s' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_2_nch_%s_nsim_%s"' "$nch" "$nsim")

    echo
    echo "=== nch=${nch}  nsim=${nsim} ==="
    "$BIN" "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" "$SCRIPT"
done

echo
echo "[local] done. outputs:"
ls -1 "$WORKDIR"/figures/data/figure_2_nch_* 2>/dev/null || true
