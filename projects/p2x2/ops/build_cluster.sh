#!/bin/bash
# Build macrodr_cli on a named cluster into a git-tagged directory.
#
# Each build lands in build/<cluster>-<tag>/ so a long-running job's binary
# is never overwritten by a later compile — the running process holds its
# file open, and future jobs are submitted against a fresh path.
#
# Usage:
#   projects/p2x2/ops/build_cluster.sh <cluster> [tag]
#
# Examples:
#   projects/p2x2/ops/build_cluster.sh dirac           # tag = git short hash (default)
#   projects/p2x2/ops/build_cluster.sh dirac nyquist   # tag = "nyquist" (custom)
#   projects/p2x2/ops/build_cluster.sh tupac
#
# Output:
#   build/<cluster>-<tag>/macrodr_cli
#   build/macrodr_cli-<cluster>-current   (symlink to the latest build for that cluster)

set -euo pipefail

OPS_DIR="$(dirname "$(readlink -f "$0")")"

CLUSTER="${1:?Usage: $0 <cluster> [tag]    (cluster profile must exist at clusters/<cluster>.sh)}"
TAG="${2:-$(git rev-parse --short HEAD)}"

PROFILE="${OPS_DIR}/clusters/${CLUSTER}.sh"
[ -f "$PROFILE" ] || { echo "[build_cluster] No such cluster profile: $PROFILE" >&2; exit 1; }

# Load module command (compute/login nodes that don't auto-init it)
[ -f /etc/profile ] && source /etc/profile

# shellcheck source=/dev/null
source "$PROFILE"   # loads modules, sets PATH_MACRO + CLUSTER, cd's to repo root

BUILD_DIR="build/${CLUSTER}-${TAG}"

cmake -S . -B "$BUILD_DIR" -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cmake --build "$BUILD_DIR" -j

BIN="$(readlink -f "$BUILD_DIR/macrodr_cli")"
mkdir -p build
ln -sfn "${CLUSTER}-${TAG}/macrodr_cli" "build/macrodr_cli-${CLUSTER}-current"

echo
echo "[build_cluster] cluster=${CLUSTER} tag=${TAG}"
echo "[build_cluster] binary:  ${BIN}"
echo "[build_cluster] current: build/macrodr_cli-${CLUSTER}-current"
echo
echo "Submit example:"
echo "  sbatch --export=ALL,BIN=${BIN} \\"
echo "         projects/eLife_2025/ops/slurm/run_macroir.sh \\"
echo "         projects/eLife_2025/ops/local/figure_2.macroir \\"
echo "         projects/eLife_2025"
