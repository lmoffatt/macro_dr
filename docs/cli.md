# MacroDR CLI Overview

This short guide documents the refreshed command-line interface that now
wraps MacroDR. The intent is to keep process-oriented concerns in the CLI
while all modelling semantics remain in the DSL.

## Basic Usage

```
macro_dr [options] <script1> [script2 ...] [-e "<dsl>" ...]
```

Scripts and `--eval` lines are processed strictly in the order they appear.
Each referenced file is concatenated into the assembled script; every eval
argument adds a single DSL line.

### Common Options

- `-e, --eval <dsl>` – append an inline DSL statement (repeatable).
- `-f, --file <path>` – enqueue a script file explicitly (repeatable).
- `-n, --check-syntax` – parse and compile, but do not execute.
- `-v, --verbose` – increase logging verbosity (repeatable).
- `-C, --chdir <dir>` – change the working directory before any IO.
- `--path <dir>` – add a search path for relative assets (repeatable).
- `-h, --help` / `--version` – print metadata and exit immediately.
- `--` – end of CLI options; subsequent arguments are treated as script
  files even if they start with `-`.

Legacy inline arguments of the form `--some_command(...)` are still accepted
for now and are internally rewritten to `--eval "some_command(...)"`. A
deprecation warning is emitted so existing ops scripts keep running while we
migrate them to the explicit `-e/--eval` form.

## Run Workspace Layout

Every invocation creates a time-stamped folder under `runs/`:

```
runs/run-YYYYMMDD-HHMMSS/
  script.macroir    # assembled script for reproducibility
  meta.json         # minimal metadata (cwd, timestamp, argv)
  ...               # future checkpoints, logs, artifacts
```

When `-v` is supplied the CLI prints the chosen run directory and the
assembled script to aid debugging.

## Path Resolution

Relative paths mentioned in scripts and DSL helpers are resolved using:

1. The current working directory (or `--chdir` target).
2. Any directories passed via `--path` (in order).
3. Directories listed in the `MACRODR_PATH` environment variable.

All matched paths are canonicalised before being opened so the DSL receives
absolute paths.

## Help and Version Inside DSL

Two convenience commands are now registered in the DSL:

- `help()` – prints the same usage summary as `macro_dr --help`.
- `version()` – prints and returns the version string.

They can be used within scripts or with `--eval` blocks.

## Next Steps

The CLI remains intentionally thin. Future work will build on the same
structure to add features such as continuations (`--run-id`), script export,
and richer verbosity controls without breaking compatibility with existing
DSL scripts.
