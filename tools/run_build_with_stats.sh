#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 [--preset <name>] [--skip-tests] [--dry-run]

  --preset <name>   CMake preset to build (default: gcc-release)
  --skip-tests      Build only; skip ctest
  --dry-run         Run the aggregator in dry-run mode (no CSV updates)
  -h, --help        Show this help
USAGE
  exit 0
}

PRESET="gcc-release"
RUN_TESTS=true
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --preset)
      shift
      [[ $# -gt 0 ]] || usage
      PRESET="$1"
      ;;
    --preset=*)
      PRESET="${1#*=}"
      ;;
    --skip-tests)
      RUN_TESTS=false
      ;;
    --dry-run)
      DRY_RUN=true
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      ;;
  esac
  shift
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build/${PRESET}"
STATS_DIR="${ROOT_DIR}/stats/${PRESET}"

# Configure (only if needed)
if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
  cmake --preset "${PRESET}"
fi

# Build (clean)
cmake --build --preset "${PRESET}" --clean-first

# Tests
if "${RUN_TESTS}"; then
  ctest --test-dir "${BUILD_DIR}" -V
fi

# Collect artifacts for stats
mkdir -p "${STATS_DIR}"
cp "${BUILD_DIR}/.ninja_log" "${STATS_DIR}/ninja_log"
if [[ -f "${BUILD_DIR}/compile_commands.json" ]]; then
  cp "${BUILD_DIR}/compile_commands.json" "${STATS_DIR}/compile_commands.json"
fi
lscpu > "${STATS_DIR}/cpu.txt"

# Run aggregator (append to docs/perf/compile/*.csv)
AGG_CMD=(
  python "${ROOT_DIR}/tools/stats/aggregate_ninja_log.py"
  --ninja-log "${STATS_DIR}/ninja_log"
  --preset "${PRESET}"
  --files-csv "${ROOT_DIR}/docs/perf/compile/files.csv"
  --summary-csv "${ROOT_DIR}/docs/perf/compile/summary.csv"
  --cpu-info "${STATS_DIR}/cpu.txt"
)
if "${DRY_RUN}"; then
  AGG_CMD+=(--dry-run)
fi
"${AGG_CMD[@]}"
