#!/bin/bash
# Run a .macroir script via SLURM, single-node, OMP-parallel within the binary.
#
# The .macroir file is the unit of work: its axes (algorithm, scheme, noise,
# Nchannels, ...) are dispatched internally by macrodr_cli via OMP across the
# task's CPUs, so this wrapper is intentionally thin — no SLURM_LOCALID arrays,
# no continuation chain, no scheme/lik/experiment env-var fan-out.
#
# Required env (set in your interactive shell before sbatch; SLURM passes them
# in via --export=ALL by default):
#   CLUSTER   : cluster name (selects projects/p2x2/ops/clusters/<CLUSTER>.sh).
#               Set by sourcing that profile in your shell first.
#   PATH_MACRO: parent of the macro_dr repo (set by the cluster profile).
#   BIN       : absolute path to the macrodr_cli to pin this job to,
#               typically $PWD/build/<cluster>-<git-tag>/macrodr_cli.
#               Pinning lets later rebuilds proceed without disturbing
#               in-flight jobs (each holds its binary's inode open).
#
# Usage:
#   source projects/p2x2/ops/clusters/dirac.sh
#   projects/p2x2/ops/build_cluster.sh dirac     # produces build/dirac-<hash>/macrodr_cli
#
#   sbatch --export=ALL,BIN=$PWD/build/dirac-<hash>/macrodr_cli \
#       [--cpus-per-task=N] [--time=…] [--partition=…] \
#       projects/eLife_2025/ops/slurm/run_macroir.sh \
#       <input.macroir> [workdir]

#SBATCH --job-name=macroir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

INPUT="${1:?Usage: sbatch $0 <input.macroir> [workdir]}"
WORKDIR="${2:-$(pwd)}"
INPUT_ABS="$(readlink -f "$INPUT")"
WORKDIR_ABS="$(readlink -f "$WORKDIR")"

: "${CLUSTER:?Set CLUSTER (source projects/p2x2/ops/clusters/<name>.sh in your shell first)}"
: "${PATH_MACRO:?PATH_MACRO not set — cluster profile incomplete?}"
: "${BIN:?Set BIN=<absolute path to macrodr_cli>, e.g. \$PWD/build/${CLUSTER}-<hash>/macrodr_cli}"

# Compute-node shells may not have the module command initialised
[ -f /etc/profile ] && source /etc/profile

# Re-source the cluster profile so module loads happen in this shell too
# shellcheck source=/dev/null
source "${PATH_MACRO}/macro_dr/projects/p2x2/ops/clusters/${CLUSTER}.sh"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export BLIS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "$WORKDIR_ABS"
echo "[run_macroir] cluster=${CLUSTER}"
echo "[run_macroir] cwd=${PWD}"
echo "[run_macroir] bin=${BIN}"
echo "[run_macroir] input=${INPUT_ABS}"
echo "[run_macroir] threads=${OMP_NUM_THREADS}"

srun "$BIN" "$INPUT_ABS"
