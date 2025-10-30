# Testing and Debugging – CLI Recipes

This guide shows how to build, run, filter, and debug the test suite entirely from the command line, without relying on editor integrations. It complements `docs/testing.md`.

## Presets and Build

- Configure Debug:
  - `cmake --preset gcc-debug`
- Build tests target:
  - `cmake --build --preset gcc-debug --target macrodr_tests`

Equivalent commands apply to Release via `gcc-release`.

## Running Tests with CTest

- Run all tests (Debug):
  - `ctest --test-dir build/gcc-debug -V --output-on-failure`

- Run only the Catch2 suite (skipping CLI smoke tests):
  - `ctest --test-dir build/gcc-debug -R ^macrodr_tests$ -V`

### Single test with runtime filter

There is a single CTest entry `macrodr_tests`. You can filter what it runs at runtime using an environment variable.

- MacroIR tests:
  - `CATCH2_FILTER="[macroir]" ctest --test-dir build/gcc-debug -R macrodr_tests -V`
- QDT derivatives only:
  - `CATCH2_FILTER="[macroir][qdt]" ctest --test-dir build/gcc-debug -R macrodr_tests -V`
- MacroIR scaffolding only (exclude QDT):
  - `CATCH2_FILTER="[macroir]~[qdt]" ctest --test-dir build/gcc-debug -R macrodr_tests -V`
- List available tests via CTest wrapper:
  - `CATCH2_LIST_TESTS=1 ctest --test-dir build/gcc-debug -R macrodr_tests -V`
- List tags:
  - `CATCH2_LIST_TAGS=1 ctest --test-dir build/gcc-debug -R macrodr_tests -V`

The CTest entry runs with working directory set to `build/<preset>/tests` and `OMP_NUM_THREADS=1` for determinism, so your relative data paths resolve.

## Running the Test Binary Directly (Recommended for Debugging)

The test binary is Catch2 v3 and supports tag expressions. Running it directly gives you full control and a clean debugging context.

- Change to the tests build dir so relative data paths resolve:
  - `cd build/gcc-debug/tests`
- List tests and tags:
  - `./macrodr_tests --list-tests`
  - `./macrodr_tests --list-tags`
- Run MacroIR tests by tag:
  - `./macrodr_tests "[macroir]"`
- Run only QDT derivatives:
  - `./macrodr_tests "[macroir][qdt]"`
- Exclude QDT (MacroIR scaffolding only):
  - `./macrodr_tests "[macroir]~[qdt]"`

## Debugging with GDB

Set the working directory and pass a tag filter. The build sets up test data under `build/<preset>/data`, and tests expect `../data/...` relative to `build/<preset>/tests`.

- `cd build/gcc-debug/tests`
- `OMP_NUM_THREADS=1 gdb --args ./macrodr_tests "[macroir]"`

Useful GDB one-liners:

- Start in the correct directory:
  - `gdb -ex "cd $(pwd)" --args ./macrodr_tests "[macroir]"`
- Break on first assertion failure message path (example):
  - `break tests/macroir/test_macroir_derivatives.cpp:1` (adjust as needed)

## CLI Smoke Tests

CLI tests (help/version/check-syntax) are also registered in CTest and can be invoked via:

- `ctest --test-dir build/gcc-debug -R macrodr_cli_ -V`

## Troubleshooting

- “File not found” for test data:
  - Ensure the working directory is `build/<preset>/tests` (CTest entries already set this). Data is copied to `build/<preset>/data` during configuration.
- Filters run “too many” or “too few” tests:
  - Use `--list-tests` and `--list-tags` to refine the Catch2 tag expression.
- Non-deterministic behavior while debugging:
  - Keep `OMP_NUM_THREADS=1` (CTest entries set this). Increase only when needed for performance.
