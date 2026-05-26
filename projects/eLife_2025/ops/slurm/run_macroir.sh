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
#   CLUSTER   : cluster name (selects projects/eLife_2025/ops/clusters/<CLUSTER>.sh).
#   BIN       : absolute path to the macrodr_cli to pin this job to,
#               typically $PWD/build/<cluster>-<git-tag>/macrodr_cli.
#               Pinning lets later rebuilds proceed without disturbing
#               in-flight jobs (each holds its binary's inode open).
# Optional env:
#   WORKDIR   : directory to cd into before running (default: cwd). The binary
#               writes its relative output paths (figures/data/...) from here,
#               so point it at scratch.
#
# Positional args = the binary's argv, in order. Each is either:
#   - a script file (any arg NOT starting with --): resolved to an absolute
#     path before the cd, so relative paths work; or
#   - an inline DSL injection (starts with --): forwarded verbatim. The binary
#     strips the -- and concatenates it into the program at THIS position, so
#     an injection placed between a head and body file overrides a value before
#     its consumer runs (a trailing injection lands too late — see the split-
#     template pattern in the ops docs).
#
# Usage:
#   source projects/eLife_2025/ops/clusters/dirac.sh
#   projects/eLife_2025/ops/build_cluster.sh dirac     # produces build/dirac-<hash>/macrodr_cli
#
#   sbatch --export=ALL,CLUSTER=dirac,BIN=$PWD/build/dirac-<hash>/macrodr_cli,WORKDIR=/scratch/$USER/macro_dr/eLife_2025 \
#       [--cpus-per-task=N] [--time=…] [--partition=…] \
#       projects/eLife_2025/ops/slurm/run_macroir.sh \
#       fig2_head.macroir "--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[100])" fig2_body.macroir

#SBATCH --job-name=macroir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%j.out

# -e: abort on first error; pipefail: pipelines too.
# Intentionally no -u: OpenHPC's /etc/profile.d/*.sh references unset env vars
# (e.g. COLORTERM) and we must source /etc/profile to get the module function.
set -eo pipefail

# CLUSTER selects the profile; BIN pins the binary. PATH_MACRO is NOT required
# up front — the cluster profile sets it when sourced below.
: "${CLUSTER:?Set CLUSTER=<name> (must match a profile in ../clusters/)}"
: "${BIN:?Set BIN=<absolute path to macrodr_cli>, e.g. \$PWD/build/${CLUSTER}-<hash>/macrodr_cli}"

[ "$#" -ge 1 ] || {
    echo "usage: [WORKDIR=dir] $0 <file.macroir | --key=value> ..." >&2
    exit 1
}

WORKDIR="${WORKDIR:-$(pwd)}"
WORKDIR_ABS="$(readlink -f "$WORKDIR")"

# Build the binary's argv: absolutize file args (before the cd below makes
# relative paths stale); forward --key=value injections verbatim, in order.
prog_args=()
for a in "$@"; do
    case "$a" in
        --*) prog_args+=("$a") ;;
        *)   prog_args+=("$(readlink -f "$a")") ;;
    esac
done

# Compute-node shells may not have the module command initialised
[ -f /etc/profile ] && source /etc/profile

# Re-source the cluster profile so module loads happen in this shell too.
# Resolved relative to this script so the wrapper survives future moves.
WRAPPER_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck source=/dev/null
source "${WRAPPER_DIR}/../clusters/${CLUSTER}.sh"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export BLIS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "$WORKDIR_ABS"
echo "[run_macroir] cluster=${CLUSTER}"
echo "[run_macroir] cwd=${PWD}"
echo "[run_macroir] bin=${BIN}"
echo "[run_macroir] threads=${OMP_NUM_THREADS}"
echo "[run_macroir] argv:"
for a in "${prog_args[@]}"; do echo "    ${a}"; done

srun "$BIN" "${prog_args[@]}"
