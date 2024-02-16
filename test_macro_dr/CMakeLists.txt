cmake_minimum_required(VERSION 3.5)

project(test_macro_dr LANGUAGES CXX)

enable_testing()

SET(CMAKE_CXX_STANDARD 20)


add_executable(test_macro_dr tst_test_macro_dr.cpp main.cpp
    examples/N_simulation.csv)
add_test(NAME test_macro_dr COMMAND test_macro_dr)


if (DEFINED ENV{CATCH_INCLUDE_DIR})
    set(CATCH_INCLUDE_DIR $ENV{CATCH_INCLUDE_DIR})
else ()
    set(CATCH_INCLUDE_DIR "/home/lmoffatt/Code/macro_dr/macro_dr/catch2") # set by Qt Creator wizard
endif ()
if (CATCH_INCLUDE_DIR STREQUAL "")
    message("CATCH_INCLUDE_DIR is not set, assuming Catch2 can be found automatically in your system")
elseif (EXISTS ${CATCH_INCLUDE_DIR})
    include_directories(${CATCH_INCLUDE_DIR})
endif ()
