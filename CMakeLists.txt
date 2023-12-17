cmake_minimum_required(VERSION 3.5)

project(macro_dr LANGUAGES CXX)
find_package(OpenMP REQUIRED)

# git commit hash macro
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")




set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=0 -Werror=return-type -fconcepts-diagnostics-depth=4" )
#find_package(Catch2 REQUIRED)
#target_link_libraries(tests Catch2::Catch2)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

add_executable(macro_dr main.cpp
   # CLI_base.h
   maybe_error.h
   CLI_regular_types.h
   qmodel.h
   run_script.txt
   CLI_macro_dr.h
   CLI_grammar.h
   grammar_untyped.h
    grammar_Identifier.h
    lexer_untyped.h
    test.txt
    matrix.h
    lapack_headers.h
    experiment.h
    variables.h
    fold.h
    derivative_operator.h
    parameters.h
    grammar_typed.h
    lexer_typed.h
    indexed.h
    matrix_derivative.h
    parameters_derivative.h
    variables_derivative.h
    derivative_test.h
    multivariate_normal_distribution.h
    parallel_tempering.h
    mcmc.h
    distributions.h
    bayesian_linear_regression.h
    parallel_tempering_linear_regression.h
    cuevi.h
    matrix_random.h
    random_samplers.h
    parameters_distribution.h
    general_algorithm_on_containers.h
    general_output_operator.h
    allosteric_models.h
    exponential_matrix.h
    function_measure_verification_and_optimization.h
    function_memoization.h
    micror_stochastic.h
    lgamma.h
    s_first
    Moffatt_Hume_2007_ATP_time.txt
    Moffatt_Hume_2007_ATP_time_7.txt
    Moffatt_Hume_2007_ATP_time_7_constant_dt.txt
    s_debug.sh
    s_long.sh
    continuation.h
    
)
target_link_libraries(macro_dr ${lapackblas_libraries}  OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
