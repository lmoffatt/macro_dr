# Compile-Time Statistics

CI appends two CSV files here after every build:

- `files.csv` — one row per compiled translation unit.
- `summary.csv` — one row per build preset summarizing the run.

Both CSVs include timestamps, commit information, preset name, version, and
CPU model so we can track regressions over time. The summary totals exclude
compilation units living under third-party directories (e.g. `_deps/`,
`third_party/`) so external dependencies do not dominate the charts.  Please
keep these files append-only.
