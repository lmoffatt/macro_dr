#!/bin/bash
# Build macrodr_cli on a named cluster into a git-tagged directory.
#
# Each build lands in build/<cluster>-<tag>/ so a long-running job's binary
# is never overwritten by a later compile — the running process holds its
# file open, and future jobs are submitted against a fresh path.
#
# Usage:
#   projects/eLife_2025/ops/build_cluster.sh <cluster> [tag]
#
# Examples:
#   projects/eLife_2025/ops/build_cluster.sh dirac           # tag = git short hash (default)
#   projects/eLife_2025/ops/build_cluster.sh dirac nyquist   # tag = "nyquist" (custom)
#   projects/eLife_2025/ops/build_cluster.sh tupac
#
# Output:
#   build/<cluster>-<tag>/macrodr_cli
#   build/macrodr_cli-<cluster>-current   (symlink to the latest build for that cluster)

# -e: abort on first error; pipefail: pipelines too.
# Intentionally no -u: OpenHPC's /etc/profile.d/*.sh references unset env vars
# (e.g. COLORTERM) and we must source /etc/profile to get the module function.
set -eo pipefail

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
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# cc1plus needs up to ~16 GB on the heaviest TUs (command_manager.cpp,
# legacy/qmodel.h users — their template instantiations are very wide).
# Recommended allocation: --cpus-per-task=2 --mem=64G  (32 GB per parallel
# job). -j defaults to the SLURM allocation; outside SLURM falls back to -j 1
# to avoid login-node OOM. Override with BUILD_JOBS=N to retune without
# re-allocating, e.g. BUILD_JOBS=1 projects/.../build_cluster.sh dirac.
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-1}}"
echo "[build_cluster] compiling with -j ${BUILD_JOBS}"
cmake --build "$BUILD_DIR" -j "${BUILD_JOBS}"

BIN="$(readlink -f "$BUILD_DIR/macrodr_cli")"
mkdir -p build
ln -sfn "${CLUSTER}-${TAG}/macrodr_cli" "build/macrodr_cli-${CLUSTER}-current"

echo
echo "[build_cluster] cluster=${CLUSTER} tag=${TAG}"
echo "[build_cluster] binary:  ${BIN}"
echo "[build_cluster] current: build/macrodr_cli-${CLUSTER}-current"
echo
echo "Submit example:"
echo "  sbatch --export=ALL,CLUSTER=${CLUSTER},BIN=${BIN},WORKDIR=/scratch/\$USER/macro_dr/eLife_2025 \\"
echo "         projects/eLife_2025/ops/slurm/run_macroir.sh \\"
echo "         projects/eLife_2025/ops/local/figure_2.macroir"
echo
echo "  # with overrides (head + injection + body):"
echo "  #   ... run_macroir.sh fig2_head.macroir \"--Num_ch = ...\" fig2_body.macroir"
