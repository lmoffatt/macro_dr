#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: tools/clang-tidy-complexity.sh [-p PATH] [-j N] [-o LOG] [-- EXTRA_ARGS]

Runs clang-tidy's readability-function-cognitive-complexity across the repo
using the detected compile_commands.json and prints a Top Offenders list and a
histogram.

Options:
  -p PATH   Build dir or directory containing compile_commands.json (auto by default)
  -j N      Parallel jobs (default: number of CPUs)
  -o LOG    Output log file (default: mktemp under /tmp)
  -h        Show this help

Examples:
  tools/clang-tidy-complexity.sh
  tools/clang-tidy-complexity.sh -p build/gcc-debug -j 8 -o /tmp/ctidy.log
  tools/clang-tidy-complexity.sh -- src include              # restrict files
USAGE
}

auto_detect_build_path() {
  # Prefer repo-root symlink if present
  if [[ -f ./compile_commands.json ]]; then
    echo "."
    return 0
  fi
  # Common local build locations
  for d in \
    build/gcc-debug \
    build/gcc-release \
    stats/gcc-release \
    build \
    .
  do
    if [[ -f "$d/compile_commands.json" ]]; then
      echo "$d"
      return 0
    fi
  done
  return 1
}

JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 4)"
BUILD_PATH=""
LOG_FILE=""

while getopts ":p:j:o:h" opt; do
  case "$opt" in
    p) BUILD_PATH="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    o) LOG_FILE="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -$OPTARG" >&2; usage; exit 2 ;;
    \?) echo "Unknown option -$OPTARG" >&2; usage; exit 2 ;;
  esac
done
shift $((OPTIND-1))

if [[ -z "${BUILD_PATH}" ]]; then
  if ! BUILD_PATH="$(auto_detect_build_path)"; then
    echo "Could not locate compile_commands.json. Build with a CMake preset first." >&2
    echo "Hint: cmake --preset gcc-debug && cmake --build --preset gcc-debug" >&2
    exit 1
  fi
fi

if ! command -v run-clang-tidy >/dev/null 2>&1; then
  echo "run-clang-tidy not found. Install clang-tools (e.g., apt-get install clang-tools) or ensure it's on PATH." >&2
  exit 1
fi

if [[ -z "${LOG_FILE}" ]]; then
  LOG_FILE="$(mktemp /tmp/clang-tidy-complexity.XXXXXX.log)"
fi

echo "Using compile commands from: $BUILD_PATH" >&2
echo "Writing full clang-tidy output to: $LOG_FILE" >&2

# Allow callers to pass file filters after a literal --
EXTRA_ARGS=("$@")

# Only run the single check we care about here.
run-clang-tidy -p "$BUILD_PATH" -j "$JOBS" \
  -checks='-*,readability-function-cognitive-complexity' \
  "${EXTRA_ARGS[@]}" 2>&1 | tee "$LOG_FILE"

echo
echo "Top offenders (score, function, file):"
grep -E "cognitive complexity of [0-9]+" "$LOG_FILE" \
  | sed -E 's#^(.*):[0-9]+:[0-9]+: warning: Function '\''([^'\'']+)\'\'' has cognitive complexity of ([0-9]+).*#\3\t\2\t\1#g' \
  | sort -k1,1nr | head -n 30

echo
echo "Histogram (bins of 25):"
grep -Eo "cognitive complexity of [0-9]+" "$LOG_FILE" \
  | awk '{n=$4; b=int(n/25)*25; h[b]++} END{for(k in h) printf "%3d–%3d: %d\n",k,k+24,h[k]}' \
  | sort -n

echo
echo "Done. Full output at: $LOG_FILE"

