#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 [--preset <name>] [--clean]

--preset <name>   Build preset to use (default: gcc-release)
--clean           Run a clean build (invokes --clean-first)
USAGE
  exit 1
}

PRESET="gcc-release"
CLEAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --preset)
      shift || usage
      [[ $# -gt 0 ]] || usage
      PRESET="$1"
      shift
      ;;
    --preset=*)
      PRESET="${1#*=}"
      shift
      ;;
    --clean)
      CLEAN=true
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      ;;
  esac
done

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
BUILD_DIR="${ROOT_DIR}/build/${PRESET}"
TIMING_DIR="${BUILD_DIR}/timing"
CSV="${BUILD_DIR}/build_stats.csv"

mkdir -p "${BUILD_DIR}" "${TIMING_DIR}"

BUILD_CMD=(cmake --build --preset "${PRESET}")
if [[ "${CLEAN}" == true ]]; then
  BUILD_CMD+=(--clean-first)
fi

START=$(date +%s)
/usr/bin/time -f "total_time_seconds,%e" -o "${TIMING_DIR}/total_time.txt" "${BUILD_CMD[@]}" 2>&1 | tee "${TIMING_DIR}/build.log"
END=$(date +%s)
TOTAL=$((END - START))

{
  echo "run_at,$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  cat "${TIMING_DIR}/total_time.txt"
  echo "total_wall_seconds,${TOTAL}"
} >> "${CSV}"

find "${BUILD_DIR}" -name "*.o-*.ftime" -o -name "*.o-*.time" 2>/dev/null | while read -r rpt; do
  file=$(basename "$rpt")
  user=$(grep -m1 "Total user time" "$rpt" | awk '{print $4}')
  wall=$(grep -m1 "Total wall time" "$rpt" | awk '{print $4}')
  heap=$(grep -m1 "Total Heap" "$rpt" | awk '{print $4}')
  if [[ -n "$user" ]]; then
    echo "tu,$file,user,$user,wall,$wall,heap,$heap" >> "${CSV}"
  fi
done
