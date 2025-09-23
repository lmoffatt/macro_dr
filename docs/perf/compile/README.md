# Compile-Time Statistics

CI appends two CSV files here after every build:

- `files.csv` — one row per compiled translation unit.
- `summary.csv` — one row per build preset summarizing the run.

Both CSVs include timestamps, commit information, preset name, version, and
CPU model so we can track regressions over time.  Please keep these files
append-only.
