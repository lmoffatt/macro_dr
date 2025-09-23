#!/usr/bin/env python3
"""Aggregate compilation timings from a Ninja log.

This script extracts per-translation-unit and per-build summaries from a
`.ninja_log` file (produced by Ninja, used by CMake builds).  It appends rows
into two CSV files so the repository can track compile-time trends over time.

Usage (example):
    python tools/stats/aggregate_ninja_log.py \
        --ninja-log stats/gcc-release/ninja_log \
        --preset gcc-release \
        --files-csv docs/perf/compile/files.csv \
        --summary-csv docs/perf/compile/summary.csv \
        --cpu-info stats/gcc-release/cpu.txt

The script reads environment variables that are typically available in CI:
    * GITHUB_SHA (commit hash)
    * GITHUB_REF (ref name)
    * MACRODR_VERSION (optional override for project version)
If they are not present the script falls back to local Git commands.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional


@dataclass
class LogEntry:
    start_ms: int
    end_ms: int
    output: str
    command_hash: str

    @property
    def duration_seconds(self) -> float:
        duration = max(0, self.end_ms - self.start_ms)
        return duration / 1000.0


TARGET_RE = re.compile(r"^CMakeFiles/(.+?)\\.dir/(.+)\\.o$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ninja-log", required=True, type=Path,
                        help="Path to the .ninja_log file")
    parser.add_argument("--preset", required=True,
                        help="Name of the build preset that produced the log")
    parser.add_argument("--files-csv", required=True, type=Path,
                        help="CSV file where per-file rows will be appended")
    parser.add_argument("--summary-csv", required=True, type=Path,
                        help="CSV file where per-build rows will be appended")
    parser.add_argument("--cpu-info", type=Path,
                        help="Optional path to a file containing CPU information")
    parser.add_argument("--timestamp", help="Override timestamp in ISO-8601 (UTC)")
    parser.add_argument("--commit", help="Override commit hash")
    parser.add_argument("--version", help="Override version string")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print rows instead of writing to CSV")
    return parser.parse_args()


CPU_MODEL_PATTERNS = (
    "model name",
    "nombre del modelo",
)


def _parse_cpu_info_text(content: str) -> Optional[str]:
    for raw_line in content.splitlines():
        line = raw_line.strip()
        lower = line.lower()
        for pattern in CPU_MODEL_PATTERNS:
            if pattern in lower:
                return line.split(":", 1)[-1].strip()
    return None


def read_cpu_model(cpu_info: Optional[Path]) -> str:
    if cpu_info and cpu_info.exists():
        content = cpu_info.read_text(errors="ignore")
        parsed = _parse_cpu_info_text(content)
        if parsed:
            return parsed
    # Fall back to /proc/cpuinfo (locale independent)
    try:
        proc_info = Path("/proc/cpuinfo")
        if proc_info.exists():
            parsed = _parse_cpu_info_text(proc_info.read_text(errors="ignore"))
            if parsed:
                return parsed
    except OSError:
        pass
    # Final fallback: lscpu in current locale
    try:
        out = subprocess.check_output(["lscpu"], text=True)
        parsed = _parse_cpu_info_text(out)
        if parsed:
            return parsed
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    return "unknown"


def detect_commit(commit_override: Optional[str]) -> str:
    if commit_override:
        return commit_override
    env_commit = os.environ.get("GITHUB_SHA")
    if env_commit:
        return env_commit
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unknown"


def detect_version(version_override: Optional[str]) -> str:
    if version_override:
        return version_override
    env_version = os.environ.get("MACRODR_VERSION")
    if env_version:
        return env_version
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--always", "--dirty"], text=True
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unknown"


def detect_timestamp(ts_override: Optional[str]) -> str:
    if ts_override:
        return ts_override
    return dt.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def parse_ninja_log(path: Path) -> Dict[str, LogEntry]:
    if not path.exists():
        raise FileNotFoundError(f"Ninja log not found: {path}")

    entries: Dict[str, LogEntry] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split('\t')
            if len(parts) != 5:
                # unexpected format; skip the line to be resilient
                continue
            start_str, end_str, _mtime, output, cmd_hash = parts
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue
            entries[output] = LogEntry(start, end, output, cmd_hash)
    return entries


def object_to_source(object_path: str) -> tuple[str, str]:
    """Infer the target and source path from a Ninja object file."""
    match = TARGET_RE.match(object_path)
    if match:
        target = match.group(1)
        source = match.group(2)
        source = source.replace('\\', '/')
        return target, source
    # fallback: treat entire object path as target with unknown source
    # Example: external objects or custom rules
    target = Path(object_path).name
    return target, "unknown"


def ensure_csv_header(path: Path, header: Iterable[str]) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(list(header))


def append_rows(path: Path, header: Iterable[str], rows: Iterable[Iterable[str]]) -> None:
    ensure_csv_header(path, header)
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(list(row))


def main() -> None:
    args = parse_args()
    entries = parse_ninja_log(args.ninja_log)

    if not entries:
        print(f"No entries parsed from {args.ninja_log}")
        return

    timestamp = detect_timestamp(args.timestamp)
    commit = detect_commit(args.commit)
    version = detect_version(args.version)
    cpu_model = read_cpu_model(args.cpu_info)
    preset = args.preset

    per_file_header = [
        "timestamp_utc",
        "commit",
        "version",
        "preset",
        "target",
        "source",
        "object",
        "compile_seconds",
        "cpu_model",
    ]
    summary_header = [
        "timestamp_utc",
        "commit",
        "version",
        "preset",
        "total_files",
        "compile_seconds_total",
        "link_seconds_total",
        "cpu_model",
    ]

    per_file_rows = []
    total_compile = 0.0
    total_files = 0
    total_link = 0.0

    ignored_patterns = ("/_deps/", "/third_party/")

    for output, entry in entries.items():
        duration = entry.duration_seconds
        if output.endswith(".o"):
            target, source = object_to_source(output)
            per_file_rows.append([
                timestamp,
                commit,
                version,
                preset,
                target,
                source,
                output,
                f"{duration:.6f}",
                cpu_model,
            ])
            if not any(p in output or p in source for p in ignored_patterns):
                total_compile += duration
                total_files += 1
        else:
            if not any(p in output for p in ignored_patterns):
                total_link += duration

    summary_rows = [[
        timestamp,
        commit,
        version,
        preset,
        str(total_files),
        f"{total_compile:.6f}",
        f"{total_link:.6f}",
        cpu_model,
    ]]

    if args.dry_run:
        print("Per-file rows:")
        for row in per_file_rows:
            print(",".join(row))
        print("Summary rows:")
        for row in summary_rows:
            print(",".join(row))
        return

    append_rows(Path(args.files_csv), per_file_header, per_file_rows)
    append_rows(Path(args.summary_csv), summary_header, summary_rows)
    print(f"Appended {len(per_file_rows)} per-file rows and 1 summary row.")


if __name__ == "__main__":
    main()
