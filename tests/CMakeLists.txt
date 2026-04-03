set(TEST_SOURCES
    cmd/test_write_csv.cpp
    cli/test_cli_helpers.cpp
    cli/test_dsl_core_refactor.cpp
    cli/test_env_persistence.cpp
    cli/test_indexed.cpp
    io/json/test_convert.cpp
    math/test_eigs.cpp
    math/test_eigensystem_contract.cpp
    math/test_c_h_r_c_h.cpp
    math/test_dcc_matrix.cpp
    math/test_empty_self_product.cpp
    math/test_idm_eigs.cpp
    math/test_idm_matrix.cpp
    math/test_maybe_error.cpp
    math/test_subspace_distortion.cpp
    math/test_probit_statistics.cpp
    math/test_taylor_convergence.cpp
    macroir/test_qdt_derivatives.cpp
    macroir/test_macroir_derivatives.cpp
    ${PROJECT_SOURCE_DIR}/third_party/catch2/catch_amalgamated.cpp
)

add_executable(macrodr_tests ${TEST_SOURCES})

target_include_directories(macrodr_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/third_party/catch2)

target_link_libraries(macrodr_tests PRIVATE
    macrodr_core
    macrodr_cli_helpers)

# Make test data available at runtime.
# Tests refer to files via relative paths like "../data/..." assuming the
# working directory is the tests binary dir. Copy the source data directory
# into the top-level build dir so ../data resolves to ${CMAKE_BINARY_DIR}/data.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/
     DESTINATION ${CMAKE_BINARY_DIR}/data)

# Single CTest entry that can be filtered at runtime via env var.
# Usage:
#   CATCH2_FILTER="[macroir]~[qdt]" ctest -R macrodr_tests -V --test-dir build/<preset>
add_test(NAME macrodr_tests
         COMMAND ${CMAKE_COMMAND}
                 -DPROGRAM=$<TARGET_FILE:macrodr_tests>
                 -DWORKING_DIR=${CMAKE_CURRENT_BINARY_DIR}
                 -P ${PROJECT_SOURCE_DIR}/tools/ctest_catch2_wrapper.cmake)
set_tests_properties(macrodr_tests PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # Keep unstable MacroIR Clarke scaffolding as opt-in diagnostic:
    #   CATCH2_FILTER="[macroir][derivatives][macroir]" ctest -R macrodr_tests -V --test-dir build/<preset>
    ENVIRONMENT "OMP_NUM_THREADS=1;CATCH2_FILTER=~*scaffolding*")

add_test(NAME macrodr_cli_help COMMAND $<TARGET_FILE:macrodr_cli> --help)
set_tests_properties(macrodr_cli_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage")

add_test(NAME macrodr_cli_version COMMAND $<TARGET_FILE:macrodr_cli> --version)
set_tests_properties(macrodr_cli_version PROPERTIES PASS_REGULAR_EXPRESSION "macro_dr")

add_test(NAME macrodr_cli_check_syntax COMMAND $<TARGET_FILE:macrodr_cli> --check-syntax -e "version()")

add_test(NAME macrodr_cli_env_flags
         COMMAND $<TARGET_FILE:macrodr_cli>
                 --env-save end
                 --env-save-path ${CMAKE_CURRENT_BINARY_DIR}/env_flags
                 --check-syntax
                 -e "alpha = 1.0")

add_test(NAME macrodr_cli_likelihood_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/likelihood_smoke.macroir --env-save step)
set_tests_properties(macrodr_cli_likelihood_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME macrodr_cli_write_csv_vector_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/write_csv_vector_smoke.macroir)
set_tests_properties(macrodr_cli_write_csv_vector_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME macrodr_cli_likelihood_diagnostics_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/likelihood_diagnostics_smoke.macroir)
set_tests_properties(macrodr_cli_likelihood_diagnostics_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME macrodr_cli_likelihood_diagnostics_contains_dib
         COMMAND ${CMAKE_COMMAND}
                 -DPROGRAM=$<TARGET_FILE:macrodr_cli>
                 -DWORKING_DIR=${CMAKE_CURRENT_BINARY_DIR}
                 -DSCRIPT=../data/likelihood_diagnostics_smoke.macroir
                 -DOUTPUT_FILE=likelihood_diagnostics_smoke.csv
                 -DSUMMARY_COMPONENT=Probit_statistics_Distortion_Induced_Bias
                 -DEVOLUTION_COMPONENT=Distortion_Induced_Bias.Distortion_Induced_Bias
                 -P ${PROJECT_SOURCE_DIR}/tools/ctest_run_and_assert_file_contains.cmake)
set_tests_properties(macrodr_cli_likelihood_diagnostics_contains_dib PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME macrodr_cli_write_csv_likelihood_overload_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/write_csv_likelihood_overload_smoke.macroir)
set_tests_properties(macrodr_cli_write_csv_likelihood_overload_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME macrodr_cli_write_csv_likelihood_indexed_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/write_csv_likelihood_indexed_smoke.macroir)
set_tests_properties(macrodr_cli_write_csv_likelihood_indexed_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME macrodr_cli_figure_2_indexed_smoke
         COMMAND $<TARGET_FILE:macrodr_cli> ../data/figure_2_indexed.macroir)
set_tests_properties(macrodr_cli_figure_2_indexed_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1")
