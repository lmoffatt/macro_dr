#!/usr/bin/env bash
set -euo pipefail

# Always operate from the repo root (script may be called from anywhere)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

tool=${1:-gcc}   # gcc|clang
cfg=${2:-debug}  # debug|release
preset="${tool}-${cfg}"

echo ">> configure preset: ${preset}"
cmake --preset "${preset}"

echo ">> build preset: ${preset}"
cmake --build --preset "${preset}" -j

echo ">> point clangd to active preset"
ln -sf "build/${preset}/compile_commands.json" compile_commands.json
echo "OK: active preset = ${preset} (repo: ${repo_root})"
