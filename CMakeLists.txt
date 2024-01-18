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
   scripts/run_script.txt
   CLI_macro_dr.h
   CLI_grammar.h
   grammar_untyped.h
    grammar_Identifier.h
    lexer_untyped.h
    scripts/test.txt
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
    slurm/s_first
    experiments/Moffatt_Hume_2007_ATP_time.txt
    experiments/Moffatt_Hume_2007_ATP_time_7.txt
    experiments/Moffatt_Hume_2007_ATP_time_7_constant_dt.txt
    slurm/s_debug.sh
    slurm/s_long.sh
    continuation.h
    models.h
    old_code_in_main.cpp
    CLI_macrodr_commands.h
    scripts/run_evidence.txt
    scripts/prior_model4.txt
    scripts/likelihood_model4.txt
    scripts/simulation_mode4_7.txt
    experiments/experiment_7.txt
    scripts/algorithm_cuevi_standard.txt
    scripts/function_table_standard.txt
    models/model6_prior.csv models/model6_parameters.csv
    scripts/mode6_Eff_no_inactivation.txt
    scripts/experiment_7.txt
    scripts/simulation.txt
    scripts/evidence_simulation.txt
    scripts/evidence_data.txt
    scripts/model.txt
    models/model9_prior.csv models/model9_par.csv models/model8_prior.csv models/model8_par.csv models/model7_prior.csv models/model7_par.csv models/model6_Eff_no_inactivation_prior.csv models/model6_Eff_no_inactivation_par.csv models/model6_no_inactivation_prior.csv models/model6_no_inactivation_par.csv models/model6_prior.csv models/model6_par.csv models/model4_g_lin_prior.csv models/model4_g_lin_par.csv models/model4_prior.csv models/model4_par.csv models/model01_prior.csv models/model01_par.csv models/model00_7_prior.csv models/model00_7_par.csv models/model00_prior.csv models/model00_par.csv
    multi_task/multi_task_slurm.sh
    slurm/s_multi_debug.sh
    multi_task/simulation_4.sh
    models/model00_7.txt
    models/model6.txt models/model9.txt models/model8.txt models/model7.txt models/model4.txt models/model4_g_lin.txt models/model01.txt models/model00.txt models/model6_no_inactivation.txt models/model6_Eff_no_inactivation.txt models/model00_7.txt
    models_MoffattHume_linear.h
    models_MoffattHume_allosteric.h
    
)
target_link_libraries(macro_dr ${lapackblas_libraries}  OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
