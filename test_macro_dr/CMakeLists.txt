cmake_minimum_required(VERSION 3.5)

project(test_macro_dr LANGUAGES CXX)

enable_testing()

SET(CMAKE_CXX_STANDARD 20)
find_package(BLAS )
find_package(LAPACK )
find_package(OpenMP REQUIRED)

set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=20 -Werror=return-type -fconcepts-diagnostics-depth=4 -fmax-errors=5" )

if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    set(lapackblas_libraries /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so)
endif()



add_executable(test_macro_dr tst_test_macro_dr.cpp main.cpp
    examples/N_simulation.csv
    tst_cli_macrodr_command.cpp
    tst_cli_likelihood.cpp
    tst_macror_validation.cpp)
add_test(NAME test_macro_dr COMMAND test_macro_dr)

target_link_libraries(test_macro_dr ${lapackblas_libraries}  OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS test_macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

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
