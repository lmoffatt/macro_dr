# Temporary Error Captures

This folder is for transient captured errors that are useful during active
debugging but do not belong to the curated program structure.

Typical examples:

- current compilation errors
- current runtime errors
- captured command output used for diagnosis

## Intended use

Store short-lived files such as:

- `current_build_error.txt`
- `current_run_error.txt`

These files are useful because they provide a stable artifact that can be read,
compared, and analyzed by humans or agents.

## What does not belong here

- canonical documentation
- project deliverables
- releases
- long-term archival debugging history

## Placement rule

Use `tmp/errors/` for temporary cross-cutting debugging captures.

If an error log becomes important as part of a specific project history or
result set, move it into that project's own area, for example:

- `projects/<name>/results/errors/`
- `projects/<name>/runs/errors/`
- `projects/<name>/archive/`

## Status

This folder is operational scratch space. It is outside the curated program
architecture.
