# Testing MacroDR

## Prerequisites

- CMake 3.18+
- Ninja or Make
- GCC/Clang with C++20 support
- BLAS/LAPACK and GSL development packages (`libblas-dev liblapack-dev libgsl-dev`)

## Unit Tests (Catch2)

The project vendors Catch2 v3.5.2 under `third_party/catch2/`. Tests live in
`tests/` and are built when `BUILD_TESTING` is ON (default when using
CMakePresets).

### Running

```
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --test-dir build/gcc-release -V
```

This runs both the Catch2 unit suite (`macrodr_tests`) and CLI smoke tests
(`macrodr_cli_help`, `macrodr_cli_version`, `macrodr_cli_check_syntax`).

## DSL Regression Tests

### Dependency notes

### Measuring build time

Run `tools/measure_build.sh` (optionally with a preset and `--clean`) to build and append timing information to `build/<preset>/build_stats.csv`. The script wraps `cmake --build`, captures `/usr/bin/time`, and scrapes GCC `-ftime-report` data so you can track per-translation-unit cost.

- Default behaviour measures an incremental build.
- Pass `--clean` to measure a clean rebuild (CI uses this option).

CI runs `tools/measure_build.sh --clean --preset gcc-release` and uploads `build/gcc-release/build_stats.csv` as an artifact for trend analysis.

If BLAS/LAPACK or GSL are installed in a non-standard prefix you can set
`MACRODR_BLAS_DIR` and `MACRODR_GSL_DIR` to point at the root containing `lib/`
and `include/`. Otherwise install the packages listed above.


The CLI smoke tests exercise the binary with `--help`, `--version`, and a
`--check-syntax` run using an inline DSL command. Additional regression scripts
can be wrapped using `add_test` in `tests/CMakeLists.txt`.

### Adding a new CLI regression

```
add_test(NAME macrodr_cli_mycase
         COMMAND $<TARGET_FILE:macrodr_cli> --check-syntax my_script.macroir)
set_tests_properties(macrodr_cli_mycase PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/path/to/scripts)
```

## GitHub Actions

A workflow can run the same steps (`cmake --build`, `ctest`) on every push.
Create `.github/workflows/ci.yml` with a standard Ubuntu runner to reuse the
commands above.

## Derivative Tests and Eigen‑Based Flows

When a pipeline depends on an eigendecomposition of a nonsymmetric matrix, take care in how you test
derivatives:

- Avoid asserting derivatives of raw eigenvectors `V`/`U`. Eigenvectors are gauge‑dependent; their
  derivatives depend on the normalization/phase convention and can differ from finite‑difference
  estimates based on LAPACK’s gauge.
- Prefer gauge‑invariant objects: full reconstruction `A ≈ V Λ U^T`, spectral projectors
  `P_i = v_i u_i^T`, or the small‑time kernel `Qdt = V e^{Λ dt} U^T`.
- Disable non‑smooth stabilizers for derivative tests (e.g., use `StabilizerPolicyDisabled`) to avoid
  clamping/masking steps that break differentiability.
- If you must compare eigenvectors across parameter perturbations, align columns by maximum
  correlation with a reference and flip signs to ensure positive correlation before differencing.

See `docs/adr/ADR-003-eigenvector-derivatives-and-testing.md` and `docs/math/eigensystem.md` for
background and formulas.
