cmake_minimum_required(VERSION 3.5)

project(test_macro_dr LANGUAGES CXX)

enable_testing()

SET(CMAKE_CXX_STANDARD 20)
find_package(GSL )
find_package(BLAS )
find_package(LAPACK )
find_package(OpenMP REQUIRED)

set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=20 -Werror=return-type  -fmax-errors=5 -ftime-report" )
														    
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    set(lapackblas_libraries /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so)
endif()


add_executable(test_macro_dr tst_test_macro_dr.cpp main.cpp
    examples/N_simulation.csv
    tst_cli_macrodr_command.cpp
    tst_cli_likelihood.cpp
    tst_macror_validation.cpp
    tst_fractioner.cpp
    tst_poisson_noise_integration.cpp
    tst_qdtm.cpp
)
add_test(NAME test_macro_dr COMMAND test_macro_dr)

target_link_libraries(test_macro_dr ${lapackblas_libraries} ${GSL_LIBRARY} OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS test_macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if (DEFINED ENV{CATCH_INCLUDE_DIR})
    set(CATCH_INCLUDE_DIR $ENV{CATCH_INCLUDE_DIR})
else ()
    # Use the vendored legacy single-header Catch2 colocated with these tests
    set(CATCH_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../catch2")
endif ()

if (EXISTS ${CATCH_INCLUDE_DIR})
    # Limit scope to this target instead of globally
    target_include_directories(test_macro_dr PRIVATE ${CATCH_INCLUDE_DIR})
else ()
    message(WARNING "CATCH_INCLUDE_DIR not found; relying on relative includes to ../catch2/catch.hpp")
endif ()
