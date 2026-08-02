#!/usr/bin/env bash
# Run figure_1.macroir once per seed, each in its own workdir, so that the trace used for
# Figure 1 can be chosen from a set instead of from whatever random_device happened to give.
#
# Why per-seed workdirs and not a seed axis: `seed` is already an argument of
# simulate_with_sub_intervals (figure_1.macroir:30), and the diagnostic filenames inside the
# script are literals (figures/data/figure_1_likelihood_diagnostic_<ALGO>), which the binary
# writes relative to cwd. So one cwd per seed is enough and nothing permanent has to change.
#
# Why not seed = 0: calc_seed (legacy/mcmc.h:37-45) treats 0 as "draw from random_device" and
# run_simulations_with_sub_intervals (src/core/simulate.cpp:377) never returns or records the
# resolved value, unlike run_simulation which stamps it into the filename. A trace produced
# with seed = 0 cannot be recovered.
#
# Usage:
#   projects/eLife_2025/ops/local/sweep_figure_1_seeds.sh
#   SEEDS="1 2 3" projects/eLife_2025/ops/local/sweep_figure_1_seeds.sh
#   BIN=build/gcc-debug/macrodr_cli projects/eLife_2025/ops/local/sweep_figure_1_seeds.sh
#
# Env: BIN, SEEDS (default 1..20), OUTROOT, OMP_NUM_THREADS, BLAS_THREADS.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../../../.." && pwd)"
SCRIPT="$HERE/figure_1.macroir"

[ -f "$SCRIPT" ] || { echo "[seeds] missing $SCRIPT" >&2; exit 1; }

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
    echo "[seeds] macrodr_cli not found/executable: ${BIN}" >&2
    echo "        build it first, or pass BIN=/path/to/macrodr_cli" >&2
    exit 1
}

OUTROOT="$(readlink -f "${OUTROOT:-$HERE/figure_1_seeds}")"
SEEDS="${SEEDS:-$(seq 1 20)}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${BLAS_THREADS:-1}"
export BLIS_NUM_THREADS="${BLAS_THREADS:-1}"

# The script has exactly one `seed = 0`, on the simulate_with_sub_intervals call. Assert it,
# so a future edit that adds another one fails here instead of silently seeding half the run.
n_seed_lines="$(grep -c 'seed = 0' "$SCRIPT" || true)"
[ "$n_seed_lines" = "1" ] || {
    echo "[seeds] expected exactly one 'seed = 0' in $SCRIPT, found $n_seed_lines" >&2
    exit 1
}

mkdir -p "$OUTROOT"
echo "[seeds] binary : $BIN"
echo "[seeds] outroot: $OUTROOT"
echo "[seeds] seeds  : $(echo $SEEDS | tr '\n' ' ')"

for s in $SEEDS; do
    wd="$OUTROOT/seed_$s"
    if [ -s "$wd/figures/data/figure_1_likelihood_diagnostic_IR.csv" ]; then
        echo "[seeds] seed=$s already done, skipping"
        continue
    fi
    mkdir -p "$wd/figures/data" "$wd/logs"
    sed "s/seed = 0/seed = $s/" "$SCRIPT" > "$wd/figure_1_seed_$s.macroir"

    echo "[seeds] running seed=$s -> $wd"
    set +e
    ( cd "$wd" && "$BIN" "figure_1_seed_$s.macroir" ) 2>&1 | tee "$wd/logs/run.log"
    rc="${PIPESTATUS[0]}"
    set -e
    [ "$rc" = "0" ] || echo "[seeds] seed=$s FAILED rc=$rc (see $wd/logs/run.log)" >&2
done

echo "[seeds] done. Rank them with:"
echo "  Rscript projects/eLife_2025/figures/paper_both/rank_figure_1_seeds.R $OUTROOT"
