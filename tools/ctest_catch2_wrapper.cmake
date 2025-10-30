# A minimal CTest wrapper for Catch2 binaries that allows runtime filtering
# via environment variables, avoiding multiple add_test() entries.

if(NOT DEFINED PROGRAM)
  message(FATAL_ERROR "PROGRAM not set")
endif()

# Prefer CATCH2_FILTER, fallback to CATCH_FILTER
set(FILTER "")
if(DEFINED ENV{CATCH2_FILTER} AND NOT "$ENV{CATCH2_FILTER}" STREQUAL "")
  set(FILTER "$ENV{CATCH2_FILTER}")
elseif(DEFINED ENV{CATCH_FILTER} AND NOT "$ENV{CATCH_FILTER}" STREQUAL "")
  set(FILTER "$ENV{CATCH_FILTER}")
endif()

# Optional listing flags: set CATCH2_LIST_TESTS=1 or CATCH2_LIST_TAGS=1
set(LIST_FLAG "")
if(DEFINED ENV{CATCH2_LIST_TESTS} AND NOT "$ENV{CATCH2_LIST_TESTS}" STREQUAL "")
  set(LIST_FLAG "--list-tests")
elseif(DEFINED ENV{CATCH2_LIST_TAGS} AND NOT "$ENV{CATCH2_LIST_TAGS}" STREQUAL "")
  set(LIST_FLAG "--list-tags")
endif()

set(cmd "${PROGRAM}")
if(NOT "${LIST_FLAG}" STREQUAL "")
  list(APPEND cmd "${LIST_FLAG}")
endif()
if(NOT "${FILTER}" STREQUAL "")
  list(APPEND cmd "${FILTER}")
endif()

execute_process(
  COMMAND ${cmd}
  WORKING_DIRECTORY "${WORKING_DIR}"
  RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Test failed with exit code ${result}")
endif()

