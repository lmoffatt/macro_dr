# Testing MacroDR

## Prerequisites

- CMake 3.18+
- Ninja or Make
- GCC/Clang with C++20 support

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

