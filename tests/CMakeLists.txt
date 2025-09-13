cmake_minimum_required(VERSION 3.18)

# Toggle legacy (vendored Catch2 v2) test suite
option(MACRODR_ENABLE_LEGACY_TESTS "Build legacy Catch2 v2 tests" OFF)

if(MACRODR_ENABLE_LEGACY_TESTS)
    add_subdirectory(legacy_v2/test_macro_dr)
endif()

