set(TEST_SOURCES
    cli/test_cli_helpers.cpp
    ${PROJECT_SOURCE_DIR}/third_party/catch2/catch_amalgamated.cpp
)

add_executable(macrodr_tests ${TEST_SOURCES})

target_include_directories(macrodr_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/third_party/catch2)

target_link_libraries(macrodr_tests PRIVATE
    macrodr_core
    macrodr_cli_helpers)

add_test(NAME macrodr_tests COMMAND macrodr_tests)

add_test(NAME macrodr_cli_help COMMAND $<TARGET_FILE:macrodr_cli> --help)
set_tests_properties(macrodr_cli_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage")

add_test(NAME macrodr_cli_version COMMAND $<TARGET_FILE:macrodr_cli> --version)
set_tests_properties(macrodr_cli_version PROPERTIES PASS_REGULAR_EXPRESSION "macro_dr")

add_test(NAME macrodr_cli_check_syntax COMMAND $<TARGET_FILE:macrodr_cli> --check-syntax -e "version()")
