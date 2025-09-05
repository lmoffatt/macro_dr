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
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_BUILD_TYPE "Release")

#set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=20 -Werror=return-type -fconcepts-diagnostics-depth=4 -fmax-errors=5 -ftime-report" )

set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=20 -Werror=return-type   -fconcepts-diagnostics-depth=4 -fmax-errors=5 -ftime-report" )
#set_property(TARGET macro_dr PROPERTY CXX_STANDARD 20)


#find_package(Catch2 REQUIRED)
#target_link_libraries(tests Catch2::Catch2)
#set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
#set(THREADS_PREFER_PTHREAD_FLAG TRUE)
#find_package(Threads REQUIRED)

find_package(BLAS )
find_package(LAPACK )
find_package(GSL )


if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    set(lapackblas_libraries /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so)
endif()



if(GSL_FOUND)
  set(gsl_library ${GSL_LIBRARY})
else()
    set(gsl_library /nfs/home/lmoffatt/.nix-profile/lib/libgsl.so )
endif()

if(GSL_FOUND)
  set(gsl_dir ${GSL_INCLUDE_DIR})
else()
    set(gsl_dir /nix/store/snnl5y7hbk7b49cvfi8mhvgzrknjwz1p-gsl-2.7.1-dev/include/ )
endif()


#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)


#if(LAPACK_FOUND AND BLAS_FOUND)
#  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#endif()


add_executable(macro_dr main.cpp
   # CLI_base.h
   maybe_error.h
   CLI_regular_types.h
   qmodel.h
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
    
    experiments/Moffatt_Hume_2007_ATP_time.txt
    experiments/Moffatt_Hume_2007_ATP_time_7.txt
    experiments/Moffatt_Hume_2007_ATP_time_7_constant_dt.txt
    continuation.h
    scripts/run_evidence.txt
    scripts/experiment_7.txt
    scripts/simulation.txt
    scripts/evidence_simulation.txt
    scripts/evidence_data.txt
    models_MoffattHume_linear.h
    models_MoffattHume_allosteric.h
    models/scheme_1.txt
    models/scheme_2.txt
    models/scheme_3.txt
    models/scheme_4.txt
    models/scheme_4_prior.csv models/scheme_4_par.csv models/scheme_4_inact_prior.csv models/scheme_4_inact_par.csv models/scheme_3_prior.csv models/scheme_3_par.csv models/scheme_3_inact_prior.csv models/scheme_3_inact_par.csv models/scheme_2_prior.csv models/scheme_2_par.csv models/scheme_2_inact_prior.csv models/scheme_2_inact_par.csv models/scheme_1_prior.csv models/scheme_1_par.csv models/scheme_1_inact_prior.csv models/scheme_1_inact_par.csv
    experiments/Moffatt_Hume_2007_ATP_time_7_recording.txt
    scripts/experiment.txt
    experiments/Moffatt_Hume_2007_ATP_time_7_segments.txt
    experiments/Moffatt_Hume_2007_ATP_time_segments.txt
    
    
    
    
    
    models/scheme_2_inact.txt models/scheme_3_inact.txt models/scheme_4_inact.txt
    
    
    
    experiments/Moffatt_Hume_2007_ATP_time_segments.txt experiments/Moffatt_Hume_2007_ATP_time_recording.txt
    
    scripts/simulation_seed.txt
    type_algebra.h
    scripts/likelihood.txt
    CLI_function_table.h
    CLI_macro_dr_base.h
    scripts/evidence_thermo_simulation.txt
    scripts/evidence_thermo_data.txt
    CLI_thermo_evidence.h
    scripts/evidence_thermo_data.txt
    scripts/evidence_thermo_simulation.txt
    scripts/evidence_fraction_data.txt scripts/evidence_fraction_simulation.txt
    
    
    
    
    
    scripts/likelihood_NR.txt
    models/scheme_4_inact_NP_par.csv models/scheme_4_inact_NP_prior.csv models/scheme_4_NP_par.csv models/scheme_4_NP_prior.csv
    scripts/load_simulation.txt
    experiments/Sim_scheme_4_NP_N_simulation.csv
    
    models/scheme_4_NP.txt
    models/scheme_4_inact_NP.txt
     
    
    gsl_integrate.h
    
    
    scripts/load_simulation_7.txt
    experiments/Sim_6465237192204259400_scheme_4_inact_NP_N_simulation.csv
    
     
     
    
     
     
    
    
    
    scripts/idealize_experiment.txt
    scripts/idealize_experiment_7.txt
    scripts/run_idealize_experiment.txt
    
    
    experiments/Moffatt_Hume_2007_ATP_time_idealized.txt
    
    
    
    
    scripts/evidence_thermo_continuation.txt
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    models/scheme_1_inact_NP_par.csv
    models/scheme_1_inact_NP_prior.csv
    
    models/scheme_1_inact_NP.txt
    models/scheme_2_inact_NP_par.csv
    models/scheme_2_inact_NP_prior.csv
    models/scheme_2_inact_NP.txt
    models/scheme_3_inact_NP_par.csv
    models/scheme_3_inact_NP_prior.csv
    models/scheme_3_inact_NP.txt
    models/scheme_1_inact.txt
    models/scheme_6.txt
    models/scheme_6_inact.txt
    models/scheme_6_PI.txt
    models/scheme_6_inact_PI.txt
    models/scheme_7.txt
    models/scheme_8.txt
    models/scheme_9.txt
    models/scheme_7_inact.txt
    models/scheme_7_inact_PI.txt
    models/scheme_8_inact.txt
    models/scheme_9_inact.txt
    models/scheme_7_inact_par.csv 
    models/scheme_7_inact_prior.csv 
    models/scheme_7_inact_PI_par.csv 
    models/scheme_7_inact_PI_prior.csv 
    models/scheme_7_par.csv 
    models/scheme_7_prior.csv 
    models/scheme_8_inact_par.csv 
    models/scheme_8_inact_prior.csv 
    models/scheme_8_par.csv 
    models/scheme_8_prior.csv
    models/scheme_9_inact_par.csv
    models/scheme_9_inact_prior.csv
    models/scheme_9_par.csv
    models/scheme_9_prior.csv
    models/scheme_6_inact_par.csv
    models/scheme_6_inact_prior.csv
    models/scheme_6_par.csv
    models/scheme_6_prior.csv
    models/scheme_6_inact_PI_par.csv
    models/scheme_6_inact_PI_prior.csv
    models/scheme_6_PI_par.csv
    models/scheme_6_PI_prior.csv
    models/scheme_9_inact_PI_par.csv
    models/scheme_9_inact_PI_prior.csv
    models/scheme_9_PI_par.csv
    models/scheme_9_PI_prior.csv
    models/scheme_9_inact_PI.txt
    models/scheme_9_PI.txt
    models/scheme_8_PI.txt
    models/scheme_8_PI_par.csv
    models/scheme_8_PI_prior.csv
    models/scheme_8_inact_PI.txt
    models/scheme_8_inact_PI_prior.csv
    models/scheme_8_inact_PI_par.csv
    models/scheme_7_PI.txt
    models/scheme_7_PI_prior.csv
    models/scheme_7_PI_par.csv
    models/scheme_9_NP.txt
    models/scheme_7_inact_PI.txt
    models/scheme_6_inact_NP.txt models/scheme_6_inact_NP_par.csv models/scheme_6_inact_NP_prior.csv models/scheme_6_NP.txt models/scheme_6_NP_par.csv models/scheme_6_NP_prior.csv models/scheme_7_inact_NP.txt models/scheme_7_inact_NP_par.csv models/scheme_7_inact_NP_prior.csv models/scheme_7_NP.txt models/scheme_7_NP_par.csv models/scheme_7_NP_prior.csv models/scheme_8_inact_NP.txt models/scheme_8_inact_NP_par.csv models/scheme_8_inact_NP_prior.csv models/scheme_8_NP.txt models/scheme_8_NP_par.csv models/scheme_8_NP_prior.csv models/scheme_9_inact_NP.txt models/scheme_9_inact_NP_par.csv models/scheme_9_inact_NP_prior.csv models/scheme_9_NP.txt models/scheme_9_NP_par.csv models/scheme_9_NP_prior.csv models/scheme_1_inact_PI.txt models/scheme_1_inact_PI_par.csv models/scheme_1_inact_PI_prior.csv models/scheme_2_inact_PI.txt models/scheme_2_inact_PI_par.csv models/scheme_2_inact_PI_prior.csv models/scheme_3_inact_PI.txt models/scheme_3_inact_PI_par.csv models/scheme_3_inact_PI_prior.csv models/scheme_4_inact_PI.txt models/scheme_4_inact_PI_par.csv models/scheme_4_inact_PI_prior.csv models/scheme_4_PI.txt models/scheme_4_PI_par.csv models/scheme_4_PI_prior.csv
    
    
    
    
    
    
    
    
    
    
    
    models/scheme_5_inact_par.csv models/scheme_5_inact_prior.csv models/scheme_5_par.csv models/scheme_5_prior.csv
    models/scheme_5.txt
    models/scheme_5_inact.txt
    models/scheme_5_PI.txt
    models/scheme_5_inact_PI.txt
    models/scheme_5_inact_PI_par.csv
    models/scheme_5_inact_PI_prior.csv
    models/scheme_5_PI_par.csv
    models/scheme_5_PI_prior.csv
    models/scheme_5_inact_NP_par.csv
    models/scheme_5_inact_NP.txt
    models/scheme_5_inact_NP_prior.csv
    models/scheme_5_NP.txt
    models/scheme_5_NP_prior.csv
    models/scheme_5_NP_par.csv
    models/scheme_10_inact_par.csv models/scheme_10_inact_prior.csv models/scheme_10_par.csv models/scheme_10_prior.csv
    models/scheme_10.txt
    models/scheme_10_inact.txt
    models/scheme_10_inact_NP.txt models/scheme_10_inact_PI.txt models/scheme_10_inact.txt models/scheme_10_NP.txt models/scheme_10.txt models/scheme_10_PI.txt models/scheme_10_inact_NP_par.csv models/scheme_10_inact_NP_prior.csv models/scheme_10_NP_par.csv models/scheme_10_NP_prior.csv models/scheme_10_inact_PI_prior.csv models/scheme_10_inact_PI_par.csv models/scheme_10_PI_prior.csv models/scheme_10_PI_par.csv models/scheme_10_inact_par.csv models/scheme_10_inact_prior.csv models/scheme_10_par.csv models/scheme_10_prior.csv
    scripts/beta_32.txt
    scripts/beta_64.txt
    multi_task/M_scheme_N.sh
    
    
    models_Ag/scheme_5_NP.txt models_Ag/scheme_5_NP_par.csv models_Ag/scheme_5_PI_prior.csv models_Ag/scheme_5_NP_prior.csv models_Ag/scheme_5_inact_NP_prior.csv models_Ag/scheme_5_inact_NP_par.csv models_Ag/scheme_5_inact_NP.txt models_Ag/scheme_5_PI_par.csv models_Ag/scheme_5_inact_PI_prior.csv models_Ag/scheme_5_inact_PI_par.csv models_Ag/scheme_5_inact_PI.txt models_Ag/scheme_5_PI.txt models_Ag/scheme_5_inact.txt models_Ag/scheme_5.txt models_Ag/scheme_5_inact_par.csv models_Ag/scheme_5_inact_prior.csv models_Ag/scheme_5_par.csv models_Ag/scheme_5_prior.csv models_Ag/scheme_6_inact_NP_par.csv models_Ag/scheme_7_NP_par.csv models_Ag/scheme_6_NP_par.csv models_Ag/scheme_7_inact_NP_par.csv models_Ag/scheme_6_NP.txt models_Ag/scheme_7_inact_NP.txt models_Ag/scheme_6_inact_NP.txt models_Ag/scheme_7_NP.txt models_Ag/scheme_4_PI.txt models_Ag/scheme_2_inact_PI.txt models_Ag/scheme_1_inact_PI.txt models_Ag/scheme_4_inact_PI.txt models_Ag/scheme_3_inact_PI.txt models_Ag/scheme_6_NP_prior.csv models_Ag/scheme_6_inact_NP_prior.csv models_Ag/scheme_7_inact_NP_prior.csv models_Ag/scheme_7_NP_prior.csv models_Ag/scheme_2_inact_PI_par.csv models_Ag/scheme_3_inact_PI_prior.csv models_Ag/scheme_1_inact_PI_prior.csv models_Ag/scheme_2_inact_PI_prior.csv models_Ag/scheme_4_PI_prior.csv models_Ag/scheme_4_inact_PI_prior.csv models_Ag/scheme_3_inact_PI_par.csv models_Ag/scheme_4_inact_PI_par.csv models_Ag/scheme_4_PI_par.csv models_Ag/scheme_1_inact_PI_par.csv models_Ag/scheme_7_inact_PI.txt models_Ag/scheme_6_PI.txt models_Ag/scheme_6_inact_PI.txt models_Ag/scheme_6_inact_PI_par.csv models_Ag/scheme_6_inact_PI_prior.csv models_Ag/scheme_6_PI_par.csv models_Ag/scheme_6_PI_prior.csv models_Ag/scheme_7_inact_PI_par.csv models_Ag/scheme_7_inact_PI_prior.csv models_Ag/scheme_7_PI_par.csv models_Ag/scheme_7_PI_prior.csv models_Ag/scheme_7_PI.txt models_Ag/scheme_7_inact_par.csv models_Ag/scheme_7_inact_prior.csv models_Ag/scheme_7_par.csv models_Ag/scheme_7_prior.csv models_Ag/scheme_7_inact.txt models_Ag/scheme_7.txt models_Ag/scheme_6_inact_par.csv models_Ag/scheme_6_inact_prior.csv models_Ag/scheme_6_inact.txt models_Ag/scheme_6_par.csv models_Ag/scheme_6_prior.csv models_Ag/scheme_6.txt models_Ag/scheme_3.txt models_Ag/scheme_2_inact_NP.txt models_Ag/scheme_3_inact_NP.txt models_Ag/scheme_1_inact.txt models_Ag/scheme_1_inact_prior.csv models_Ag/scheme_3_inact_par.csv models_Ag/scheme_1_par.csv models_Ag/scheme_2_inact_par.csv models_Ag/scheme_2_inact_prior.csv models_Ag/scheme_1_prior.csv models_Ag/scheme_3_inact_prior.csv models_Ag/scheme_1_inact_par.csv models_Ag/scheme_4_inact_NP_prior.csv models_Ag/scheme_4_inact_NP_par.csv models_Ag/scheme_4_NP_prior.csv models_Ag/scheme_4_NP_par.csv models_Ag/scheme_4_inact_prior.csv models_Ag/scheme_4_inact_par.csv models_Ag/scheme_1_inact_NP.txt models_Ag/scheme_3_inact_NP_prior.csv models_Ag/scheme_3_inact_NP_par.csv models_Ag/scheme_2_inact_NP_prior.csv models_Ag/scheme_2_inact_NP_par.csv models_Ag/scheme_1_par_copy.csv models_Ag/scheme_1_inact_NP_prior.csv models_Ag/scheme_1_inact_NP_par.csv models_Ag/scheme_1.txt models_Ag/scheme_4_inact_NP.txt models_Ag/scheme_4_inact.txt models_Ag/scheme_4_NP.txt models_Ag/scheme_3_inact.txt models_Ag/scheme_2_inact.txt models_Ag/scheme_2.txt models_Ag/scheme_4.txt models_Ag/scheme_4_prior.csv models_Ag/scheme_4_par.csv models_Ag/scheme_2_prior.csv models_Ag/scheme_2_par.csv models_Ag/scheme_3_prior.csv models_Ag/scheme_3_par.csv
    models_Ag/scheme_8.txt models_Ag/scheme_8_inact.txt models_Ag/scheme_8_inact_NP.txt models_Ag/scheme_8_inact_NP_par.csv models_Ag/scheme_8_inact_NP_prior.csv models_Ag/scheme_8_inact_par.csv models_Ag/scheme_8_inact_PI.txt models_Ag/scheme_8_inact_PI_par.csv models_Ag/scheme_8_inact_PI_prior.csv models_Ag/scheme_8_inact_prior.csv models_Ag/scheme_8_NP.txt models_Ag/scheme_8_NP_par.csv models_Ag/scheme_8_NP_prior.csv models_Ag/scheme_8_par.csv models_Ag/scheme_8_PI.txt models_Ag/scheme_8_PI_par.csv models_Ag/scheme_8_PI_prior.csv models_Ag/scheme_8_prior.csv models_Ag/scheme_9.txt models_Ag/scheme_9_inact.txt models_Ag/scheme_9_inact_NP.txt models_Ag/scheme_9_inact_NP_par.csv models_Ag/scheme_9_inact_NP_prior.csv models_Ag/scheme_9_inact_par.csv models_Ag/scheme_9_inact_PI.txt models_Ag/scheme_9_inact_PI_par.csv models_Ag/scheme_9_inact_PI_prior.csv models_Ag/scheme_9_inact_prior.csv models_Ag/scheme_9_NP.txt models_Ag/scheme_9_NP_par.csv models_Ag/scheme_9_NP_prior.csv models_Ag/scheme_9_par.csv models_Ag/scheme_9_PI.txt models_Ag/scheme_9_PI_par.csv models_Ag/scheme_9_PI_prior.csv models_Ag/scheme_9_prior.csv models_Ag/scheme_10.txt models_Ag/scheme_10_inact.txt models_Ag/scheme_10_inact_NP.txt models_Ag/scheme_10_inact_NP_par.csv models_Ag/scheme_10_inact_NP_prior.csv models_Ag/scheme_10_inact_par.csv models_Ag/scheme_10_inact_PI.txt models_Ag/scheme_10_inact_PI_par.csv models_Ag/scheme_10_inact_PI_prior.csv models_Ag/scheme_10_inact_prior.csv models_Ag/scheme_10_NP.txt models_Ag/scheme_10_NP_par.csv models_Ag/scheme_10_NP_prior.csv models_Ag/scheme_10_par.csv models_Ag/scheme_10_PI.txt models_Ag/scheme_10_PI_par.csv models_Ag/scheme_10_PI_prior.csv models_Ag/scheme_10_prior.csv
    models/scheme_11_inact.txt models/scheme_11_inact_NP.txt models/scheme_11_inact_NP_par.csv models/scheme_11_inact_NP_prior.csv models/scheme_11_inact_par.csv models/scheme_11_inact_PI.txt models/scheme_11_inact_PI_par.csv models/scheme_11_inact_PI_prior.csv models/scheme_11_inact_prior.csv models/scheme_11_NP.txt models/scheme_11_NP_par.csv models/scheme_11_NP_prior.csv models/scheme_11_par.csv models/scheme_11_PI.txt models/scheme_11_PI_par.csv models/scheme_11_PI_prior.csv models/scheme_11_prior.csv
    models_Ag/scheme_11.txt models_Ag/scheme_11_inact.txt models_Ag/scheme_11_inact_NP.txt models_Ag/scheme_11_inact_NP_par.csv models_Ag/scheme_11_inact_NP_prior.csv models_Ag/scheme_11_inact_par.csv models_Ag/scheme_11_inact_PI.txt models_Ag/scheme_11_inact_PI_par.csv models_Ag/scheme_11_inact_PI_prior.csv models_Ag/scheme_11_inact_prior.csv models_Ag/scheme_11_NP.txt models_Ag/scheme_11_NP_par.csv models_Ag/scheme_11_NP_prior.csv models_Ag/scheme_11_par.csv models_Ag/scheme_11_PI.txt models_Ag/scheme_11_PI_par.csv models_Ag/scheme_11_PI_prior.csv models_Ag/scheme_11_prior.csv
    models/scheme_15.txt models/scheme_15_inact.txt models/scheme_15_inact_par.csv models/scheme_15_inact_prior.csv models/scheme_15_par.csv models/scheme_15_prior.csv
    models/scheme_15.txt models/scheme_15_inact.txt models/scheme_15_inact_NP.txt models/scheme_15_inact_NP_par.csv models/scheme_15_inact_NP_prior.csv models/scheme_15_inact_par.csv models/scheme_15_inact_PI.txt models/scheme_15_inact_PI_par.csv models/scheme_15_inact_PI_prior.csv models/scheme_15_inact_prior.csv models/scheme_15_NP.txt models/scheme_15_NP_par.csv models/scheme_15_NP_prior.csv models/scheme_15_par.csv models/scheme_15_PI.txt models/scheme_15_PI_par.csv models/scheme_15_PI_prior.csv models/scheme_15_prior.csv
    models_Ag/scheme_15.txt models_Ag/scheme_15_inact.txt models_Ag/scheme_15_inact_NP.txt models_Ag/scheme_15_inact_NP_par.csv models_Ag/scheme_15_inact_NP_prior.csv models_Ag/scheme_15_inact_par.csv models_Ag/scheme_15_inact_PI.txt models_Ag/scheme_15_inact_PI_par.csv models_Ag/scheme_15_inact_PI_prior.csv models_Ag/scheme_15_inact_prior.csv models_Ag/scheme_15_NP.txt models_Ag/scheme_15_NP_par.csv models_Ag/scheme_15_NP_prior.csv models_Ag/scheme_15_par.csv models_Ag/scheme_15_PI.txt models_Ag/scheme_15_PI_par.csv models_Ag/scheme_15_PI_prior.csv models_Ag/scheme_15_prior.csv
    
    
    
    
    models/scheme_12.txt models/scheme_12_inact.txt models/scheme_12_inact_NP.txt models/scheme_12_inact_NP_par.csv models/scheme_12_inact_NP_prior.csv models/scheme_12_inact_par.csv models/scheme_12_inact_PI.txt models/scheme_12_inact_PI_par.csv models/scheme_12_inact_PI_prior.csv models/scheme_12_inact_prior.csv models/scheme_12_NP.txt models/scheme_12_NP_par.csv models/scheme_12_NP_prior.csv models/scheme_12_par.csv models/scheme_12_PI.txt models/scheme_12_PI_par.csv models/scheme_12_PI_prior.csv models/scheme_12_prior.csv models/scheme_13.txt models/scheme_13_inact.txt models/scheme_13_inact_NP.txt models/scheme_13_inact_NP_par.csv models/scheme_13_inact_NP_prior.csv models/scheme_13_inact_par.csv models/scheme_13_inact_PI.txt models/scheme_13_inact_PI_par.csv models/scheme_13_inact_PI_prior.csv models/scheme_13_inact_prior.csv models/scheme_13_NP.txt models/scheme_13_NP_par.csv models/scheme_13_NP_prior.csv models/scheme_13_par.csv models/scheme_13_PI.txt models/scheme_13_PI_par.csv models/scheme_13_PI_prior.csv models/scheme_13_prior.csv models/scheme_14.txt models/scheme_14_inact.txt models/scheme_14_inact_NP.txt models/scheme_14_inact_NP_par.csv models/scheme_14_inact_NP_prior.csv models/scheme_14_inact_par.csv models/scheme_14_inact_PI.txt models/scheme_14_inact_PI_par.csv models/scheme_14_inact_PI_prior.csv models/scheme_14_inact_prior.csv models/scheme_14_NP.txt models/scheme_14_NP_par.csv models/scheme_14_NP_prior.csv models/scheme_14_par.csv models/scheme_14_PI.txt models/scheme_14_PI_par.csv models/scheme_14_PI_prior.csv models/scheme_14_prior.csv
    models_Ag/scheme_12.txt models_Ag/scheme_12_inact.txt models_Ag/scheme_12_inact_NP.txt models_Ag/scheme_12_inact_NP_par.csv models_Ag/scheme_12_inact_NP_prior.csv models_Ag/scheme_12_inact_par.csv models_Ag/scheme_12_inact_PI.txt models_Ag/scheme_12_inact_PI_par.csv models_Ag/scheme_12_inact_PI_prior.csv models_Ag/scheme_12_inact_prior.csv models_Ag/scheme_12_NP.txt models_Ag/scheme_12_NP_par.csv models_Ag/scheme_12_NP_prior.csv models_Ag/scheme_12_par.csv models_Ag/scheme_12_PI.txt models_Ag/scheme_12_PI_par.csv models_Ag/scheme_12_PI_prior.csv models_Ag/scheme_12_prior.csv models_Ag/scheme_13.txt models_Ag/scheme_13_inact.txt models_Ag/scheme_13_inact_NP.txt models_Ag/scheme_13_inact_NP_par.csv models_Ag/scheme_13_inact_NP_prior.csv models_Ag/scheme_13_inact_par.csv models_Ag/scheme_13_inact_PI.txt models_Ag/scheme_13_inact_PI_par.csv models_Ag/scheme_13_inact_PI_prior.csv models_Ag/scheme_13_inact_prior.csv models_Ag/scheme_13_NP.txt models_Ag/scheme_13_NP_par.csv models_Ag/scheme_13_NP_prior.csv models_Ag/scheme_13_par.csv models_Ag/scheme_13_PI.txt models_Ag/scheme_13_PI_par.csv models_Ag/scheme_13_PI_prior.csv models_Ag/scheme_13_prior.csv models_Ag/scheme_14.txt models_Ag/scheme_14_inact.txt models_Ag/scheme_14_inact_NP.txt models_Ag/scheme_14_inact_NP_par.csv models_Ag/scheme_14_inact_NP_prior.csv models_Ag/scheme_14_inact_par.csv models_Ag/scheme_14_inact_PI.txt models_Ag/scheme_14_inact_PI_par.csv models_Ag/scheme_14_inact_PI_prior.csv models_Ag/scheme_14_inact_prior.csv models_Ag/scheme_14_NP.txt models_Ag/scheme_14_NP_par.csv models_Ag/scheme_14_NP_prior.csv models_Ag/scheme_14_par.csv models_Ag/scheme_14_PI.txt models_Ag/scheme_14_PI_par.csv models_Ag/scheme_14_PI_prior.csv models_Ag/scheme_14_prior.csv
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    parallel_levenberg_tempering.h
    CLI_thermo_levenberg_evidence.h
    scripts/evidence_levenberg_data.txt
    experiments/Moffatt_Hume_2007_ATP_time_idealized_2.txt
    scripts/idealize_experiment_2.txt
    scripts/beta_16.txt
    models/scheme_1_inact_PI_prior_0001.csv
    models/scheme_1_inact_PI_0001.txt
    
    
    
    
    
    slurm/M_scheme_N_tasks.sh
    scripts/calculate_likelihood.txt
    
    
    
    
    
    
    experiments/Moffatt_Hume_2007_ATP_time_idealized_old.txt
    
    CLI_thermo_evidence_dts.h
    scripts/evidence_thermo_dts_data.txt
    scripts/beta_4.txt
    scripts/evidence_thermo_dts_continuation.txt
    
    
    
    
    
    
    
    
    
    
    
    
    
    slurm/run_CPU_32_schemes_11.sh slurm/run_CPU_32_schemes_9_10.sh slurm/run_CPU_32_schemes_7_8.sh slurm/run_CPU_32_schemes_5_6.sh slurm/run_CPU_32_schemes_3_4.sh slurm/run_CPU_16_schemes_9_to_11.sh slurm/run_CPU_16_schemes_5_to_8.sh slurm/run_CPU_16_schemes_1_to_4.sh slurm/run_CPU_8_schemes_9_to_11.sh slurm/run_CPU_8_schemes_1_to_8.sh slurm/run_CPU_4_schemes_1_to_11.sh slurm/run_CPU_2_schemes_1_to_11.sh slurm/run_CPU_1_schemes_1_to_11.sh
    scripts/likelihood_R.txt
    scripts/likelihood_DR.txt
    scripts/likelihood_SNR.txt
    scripts/likelihood_AR.txt
    scripts/likelihood_ADR.txt
    slurm/run_Lik_schemes_1_2_3_5.sh
    slurm/run_Lik_schemes_4_6_8_9.sh
    slurm/run_Lik_schemes_10_11.sh
    slurm/run_Lik_schemes_7.sh
    slurm/set_default_variables.sh
    slurm/run_32_CPU_4_schemes_1_to_8.sh
    slurm/run_32_CPU_4_schemes_1_to_8.sh slurm/run_32_CPU_4_schemes_9_to_11.sh slurm/run_32_CPU_8_schemes_1_to_4.sh slurm/run_32_CPU_8_schemes_5_to_8.sh slurm/run_32_CPU_8_schemes_9_to_11.sh slurm/run_32_CPU_16_schemes_3_4.sh slurm/run_32_CPU_16_schemes_5_6.sh slurm/run_32_CPU_16_schemes_7_8.sh slurm/run_32_CPU_16_schemes_9_10.sh slurm/run_32_CPU_16_schemes_11.sh
    slurm/run_32_CPU_32_scheme_x.sh
    slurm/run_CPU_64_scheme_x.sh
    slurm/run_32_scheme_6_9.sh
    slurm/run_32_schemes_8_10.sh
    slurm/run_32_schemes_11.sh
    slurm/run_32_schemes_4.sh
    slurm/run_64_scheme_7.sh
    slurm/run_Continuation_1.sh slurm/run_Continuation_2.sh
    clusters/capitan.sh clusters/dirac.sh clusters/serafin.sh clusters/tupac.sh
    slurm/run_32_scheme_7.sh
    multi_task/M_scheme_N_tasks.sh
    slurm/run_32_scheme_10.sh
    models_posterior/scheme_10_inact_PI_posterior_par.csv models_posterior/scheme_10_inact_PI_prior.csv models_posterior/scheme_10_inact_PI.txt
    experiments/Sim_scheme_10_simulation.txt
    scripts/simulated_experiment.txt
    slurm/run_32_scheme_10_sim.sh
    scripts/full_experiment.txt
    experiments/Moffatt_Hume_2007_ATP_full_experiment.txt experiments/Moffatt_Hume_2007_ATP_full_recording.txt
    scripts/simulated_experiment_bassel.txt
    experiments/Sim_scheme_10_bessel.txt
    slurm/run_32_scheme_10_bassel.sh
    experiments/Moffatt_Hume_2007_ATP_Nyquist_experiment.txt experiments/Moffatt_Hume_2007_ATP_Nyquist_recording.txt
    scripts/nyquist_experiment.txt
    slurm/run_32_scheme_10_nyquist.sh
    slurm/run_32_scheme_10_NR.sh
    slurm/run_32_scheme_10_SNR.sh
    slurm/run_32_scheme_10_DR.sh
    slurm/run_32_scheme_10_ADR
    scripts/likelihood_SR.txt
    scripts/likelihood_ASR.txt
    slurm/run_32_scheme_10_AR.sh
    slurm/run_32_scheme_10_ASR.sh
    command_lines/run_Evidence.txt
    command_lines/test_likelihood.txt
    CLI_likelihood.h
    slurm/run_32_scheme_10_no_inact_ADR
    models_different_inact/scheme_10_inact_PI_prior.csv
    models_different_inact/scheme_10_inact_PI_par.csv models_different_inact/scheme_10_inact_PI.txt
    slurm/run_32_scheme_10_ADR_inact_prior.sh
    models_Ag_log_baseline/scheme_1.txt models_Ag_log_baseline/scheme_10.txt models_Ag_log_baseline/scheme_10_inact.txt models_Ag_log_baseline/scheme_10_inact_NP.txt models_Ag_log_baseline/scheme_10_inact_NP_par.csv models_Ag_log_baseline/scheme_10_inact_NP_prior.csv models_Ag_log_baseline/scheme_10_inact_par.csv models_Ag_log_baseline/scheme_10_inact_PI.txt models_Ag_log_baseline/scheme_10_inact_PI_par.csv models_Ag_log_baseline/scheme_10_inact_PI_prior.csv models_Ag_log_baseline/scheme_10_inact_prior.csv models_Ag_log_baseline/scheme_10_NP.txt models_Ag_log_baseline/scheme_10_NP_par.csv models_Ag_log_baseline/scheme_10_NP_prior.csv models_Ag_log_baseline/scheme_10_par.csv models_Ag_log_baseline/scheme_10_PI.txt models_Ag_log_baseline/scheme_10_PI_par.csv models_Ag_log_baseline/scheme_10_PI_prior.csv models_Ag_log_baseline/scheme_10_prior.csv models_Ag_log_baseline/scheme_11.txt models_Ag_log_baseline/scheme_11_inact.txt models_Ag_log_baseline/scheme_11_inact_NP.txt models_Ag_log_baseline/scheme_11_inact_NP_par.csv models_Ag_log_baseline/scheme_11_inact_NP_prior.csv models_Ag_log_baseline/scheme_11_inact_par.csv models_Ag_log_baseline/scheme_11_inact_PI.txt models_Ag_log_baseline/scheme_11_inact_PI_par.csv models_Ag_log_baseline/scheme_11_inact_PI_prior.csv models_Ag_log_baseline/scheme_11_inact_prior.csv models_Ag_log_baseline/scheme_11_NP.txt models_Ag_log_baseline/scheme_11_NP_par.csv models_Ag_log_baseline/scheme_11_NP_prior.csv models_Ag_log_baseline/scheme_11_par.csv models_Ag_log_baseline/scheme_11_PI.txt models_Ag_log_baseline/scheme_11_PI_par.csv models_Ag_log_baseline/scheme_11_PI_prior.csv models_Ag_log_baseline/scheme_11_prior.csv models_Ag_log_baseline/scheme_12.txt models_Ag_log_baseline/scheme_12_inact.txt models_Ag_log_baseline/scheme_12_inact_NP.txt models_Ag_log_baseline/scheme_12_inact_NP_par.csv models_Ag_log_baseline/scheme_12_inact_NP_prior.csv models_Ag_log_baseline/scheme_12_inact_par.csv models_Ag_log_baseline/scheme_12_inact_PI.txt models_Ag_log_baseline/scheme_12_inact_PI_par.csv models_Ag_log_baseline/scheme_12_inact_PI_prior.csv models_Ag_log_baseline/scheme_12_inact_prior.csv models_Ag_log_baseline/scheme_12_NP.txt models_Ag_log_baseline/scheme_12_NP_par.csv models_Ag_log_baseline/scheme_12_NP_prior.csv models_Ag_log_baseline/scheme_12_par.csv models_Ag_log_baseline/scheme_12_PI.txt models_Ag_log_baseline/scheme_12_PI_par.csv models_Ag_log_baseline/scheme_12_PI_prior.csv models_Ag_log_baseline/scheme_12_prior.csv models_Ag_log_baseline/scheme_13.txt models_Ag_log_baseline/scheme_13_inact.txt models_Ag_log_baseline/scheme_13_inact_NP.txt models_Ag_log_baseline/scheme_13_inact_NP_par.csv models_Ag_log_baseline/scheme_13_inact_NP_prior.csv models_Ag_log_baseline/scheme_13_inact_par.csv models_Ag_log_baseline/scheme_13_inact_PI.txt models_Ag_log_baseline/scheme_13_inact_PI_par.csv models_Ag_log_baseline/scheme_13_inact_PI_prior.csv models_Ag_log_baseline/scheme_13_inact_prior.csv models_Ag_log_baseline/scheme_13_NP.txt models_Ag_log_baseline/scheme_13_NP_par.csv models_Ag_log_baseline/scheme_13_NP_prior.csv models_Ag_log_baseline/scheme_13_par.csv models_Ag_log_baseline/scheme_13_PI.txt models_Ag_log_baseline/scheme_13_PI_par.csv models_Ag_log_baseline/scheme_13_PI_prior.csv models_Ag_log_baseline/scheme_13_prior.csv models_Ag_log_baseline/scheme_14.txt models_Ag_log_baseline/scheme_14_inact.txt models_Ag_log_baseline/scheme_14_inact_NP.txt models_Ag_log_baseline/scheme_14_inact_NP_par.csv models_Ag_log_baseline/scheme_14_inact_NP_prior.csv models_Ag_log_baseline/scheme_14_inact_par.csv models_Ag_log_baseline/scheme_14_inact_PI.txt models_Ag_log_baseline/scheme_14_inact_PI_par.csv models_Ag_log_baseline/scheme_14_inact_PI_prior.csv models_Ag_log_baseline/scheme_14_inact_prior.csv models_Ag_log_baseline/scheme_14_NP.txt models_Ag_log_baseline/scheme_14_NP_par.csv models_Ag_log_baseline/scheme_14_NP_prior.csv models_Ag_log_baseline/scheme_14_par.csv models_Ag_log_baseline/scheme_14_PI.txt models_Ag_log_baseline/scheme_14_PI_par.csv models_Ag_log_baseline/scheme_14_PI_prior.csv models_Ag_log_baseline/scheme_14_prior.csv models_Ag_log_baseline/scheme_15.txt models_Ag_log_baseline/scheme_15_inact.txt models_Ag_log_baseline/scheme_15_inact_NP.txt models_Ag_log_baseline/scheme_15_inact_NP_par.csv models_Ag_log_baseline/scheme_15_inact_NP_prior.csv models_Ag_log_baseline/scheme_15_inact_par.csv models_Ag_log_baseline/scheme_15_inact_PI.txt models_Ag_log_baseline/scheme_15_inact_PI_par.csv models_Ag_log_baseline/scheme_15_inact_PI_prior.csv models_Ag_log_baseline/scheme_15_inact_prior.csv models_Ag_log_baseline/scheme_15_NP.txt models_Ag_log_baseline/scheme_15_NP_par.csv models_Ag_log_baseline/scheme_15_NP_prior.csv models_Ag_log_baseline/scheme_15_par.csv models_Ag_log_baseline/scheme_15_PI.txt models_Ag_log_baseline/scheme_15_PI_par.csv models_Ag_log_baseline/scheme_15_PI_prior.csv models_Ag_log_baseline/scheme_15_prior.csv models_Ag_log_baseline/scheme_1_inact.txt models_Ag_log_baseline/scheme_1_inact_NP.txt models_Ag_log_baseline/scheme_1_inact_NP_par.csv models_Ag_log_baseline/scheme_1_inact_NP_prior.csv models_Ag_log_baseline/scheme_1_inact_par.csv models_Ag_log_baseline/scheme_1_inact_PI.txt models_Ag_log_baseline/scheme_1_inact_PI_par.csv models_Ag_log_baseline/scheme_1_inact_PI_prior.csv models_Ag_log_baseline/scheme_1_inact_prior.csv models_Ag_log_baseline/scheme_1_par.csv  models_Ag_log_baseline/scheme_1_prior.csv models_Ag_log_baseline/scheme_2.txt models_Ag_log_baseline/scheme_2_inact.txt models_Ag_log_baseline/scheme_2_inact_NP.txt models_Ag_log_baseline/scheme_2_inact_NP_par.csv models_Ag_log_baseline/scheme_2_inact_NP_prior.csv models_Ag_log_baseline/scheme_2_inact_par.csv models_Ag_log_baseline/scheme_2_inact_PI.txt models_Ag_log_baseline/scheme_2_inact_PI_par.csv models_Ag_log_baseline/scheme_2_inact_PI_prior.csv models_Ag_log_baseline/scheme_2_inact_prior.csv models_Ag_log_baseline/scheme_2_par.csv models_Ag_log_baseline/scheme_2_prior.csv models_Ag_log_baseline/scheme_3.txt models_Ag_log_baseline/scheme_3_inact.txt models_Ag_log_baseline/scheme_3_inact_NP.txt models_Ag_log_baseline/scheme_3_inact_NP_par.csv models_Ag_log_baseline/scheme_3_inact_NP_prior.csv models_Ag_log_baseline/scheme_3_inact_par.csv models_Ag_log_baseline/scheme_3_inact_PI.txt models_Ag_log_baseline/scheme_3_inact_PI_par.csv models_Ag_log_baseline/scheme_3_inact_PI_prior.csv models_Ag_log_baseline/scheme_3_inact_prior.csv models_Ag_log_baseline/scheme_3_par.csv models_Ag_log_baseline/scheme_3_prior.csv models_Ag_log_baseline/scheme_4.txt models_Ag_log_baseline/scheme_4_inact.txt models_Ag_log_baseline/scheme_4_inact_NP.txt models_Ag_log_baseline/scheme_4_inact_NP_par.csv models_Ag_log_baseline/scheme_4_inact_NP_prior.csv models_Ag_log_baseline/scheme_4_inact_par.csv models_Ag_log_baseline/scheme_4_inact_PI.txt models_Ag_log_baseline/scheme_4_inact_PI_par.csv models_Ag_log_baseline/scheme_4_inact_PI_prior.csv models_Ag_log_baseline/scheme_4_inact_prior.csv models_Ag_log_baseline/scheme_4_NP.txt models_Ag_log_baseline/scheme_4_NP_par.csv models_Ag_log_baseline/scheme_4_NP_prior.csv models_Ag_log_baseline/scheme_4_par.csv models_Ag_log_baseline/scheme_4_PI.txt models_Ag_log_baseline/scheme_4_PI_par.csv models_Ag_log_baseline/scheme_4_PI_prior.csv models_Ag_log_baseline/scheme_4_prior.csv models_Ag_log_baseline/scheme_5.txt models_Ag_log_baseline/scheme_5_inact.txt models_Ag_log_baseline/scheme_5_inact_NP.txt models_Ag_log_baseline/scheme_5_inact_NP_par.csv models_Ag_log_baseline/scheme_5_inact_NP_prior.csv models_Ag_log_baseline/scheme_5_inact_par.csv models_Ag_log_baseline/scheme_5_inact_PI.txt models_Ag_log_baseline/scheme_5_inact_PI_par.csv models_Ag_log_baseline/scheme_5_inact_PI_prior.csv models_Ag_log_baseline/scheme_5_inact_prior.csv models_Ag_log_baseline/scheme_5_NP.txt models_Ag_log_baseline/scheme_5_NP_par.csv models_Ag_log_baseline/scheme_5_NP_prior.csv models_Ag_log_baseline/scheme_5_par.csv models_Ag_log_baseline/scheme_5_PI.txt models_Ag_log_baseline/scheme_5_PI_par.csv models_Ag_log_baseline/scheme_5_PI_prior.csv models_Ag_log_baseline/scheme_5_prior.csv models_Ag_log_baseline/scheme_6.txt models_Ag_log_baseline/scheme_6_inact.txt models_Ag_log_baseline/scheme_6_inact_NP.txt models_Ag_log_baseline/scheme_6_inact_NP_par.csv models_Ag_log_baseline/scheme_6_inact_NP_prior.csv models_Ag_log_baseline/scheme_6_inact_par.csv models_Ag_log_baseline/scheme_6_inact_PI.txt models_Ag_log_baseline/scheme_6_inact_PI_par.csv models_Ag_log_baseline/scheme_6_inact_PI_prior.csv models_Ag_log_baseline/scheme_6_inact_prior.csv models_Ag_log_baseline/scheme_6_NP.txt models_Ag_log_baseline/scheme_6_NP_par.csv models_Ag_log_baseline/scheme_6_NP_prior.csv models_Ag_log_baseline/scheme_6_par.csv models_Ag_log_baseline/scheme_6_PI.txt models_Ag_log_baseline/scheme_6_PI_par.csv models_Ag_log_baseline/scheme_6_PI_prior.csv models_Ag_log_baseline/scheme_6_prior.csv models_Ag_log_baseline/scheme_7.txt models_Ag_log_baseline/scheme_7_inact.txt models_Ag_log_baseline/scheme_7_inact_NP.txt models_Ag_log_baseline/scheme_7_inact_NP_par.csv models_Ag_log_baseline/scheme_7_inact_NP_prior.csv models_Ag_log_baseline/scheme_7_inact_par.csv models_Ag_log_baseline/scheme_7_inact_PI.txt models_Ag_log_baseline/scheme_7_inact_PI_par.csv models_Ag_log_baseline/scheme_7_inact_PI_prior.csv models_Ag_log_baseline/scheme_7_inact_prior.csv models_Ag_log_baseline/scheme_7_NP.txt models_Ag_log_baseline/scheme_7_NP_par.csv models_Ag_log_baseline/scheme_7_NP_prior.csv models_Ag_log_baseline/scheme_7_par.csv models_Ag_log_baseline/scheme_7_PI.txt models_Ag_log_baseline/scheme_7_PI_par.csv models_Ag_log_baseline/scheme_7_PI_prior.csv models_Ag_log_baseline/scheme_7_prior.csv models_Ag_log_baseline/scheme_8.txt models_Ag_log_baseline/scheme_8_inact.txt models_Ag_log_baseline/scheme_8_inact_NP.txt models_Ag_log_baseline/scheme_8_inact_NP_par.csv models_Ag_log_baseline/scheme_8_inact_NP_prior.csv models_Ag_log_baseline/scheme_8_inact_par.csv models_Ag_log_baseline/scheme_8_inact_PI.txt models_Ag_log_baseline/scheme_8_inact_PI_par.csv models_Ag_log_baseline/scheme_8_inact_PI_prior.csv models_Ag_log_baseline/scheme_8_inact_prior.csv models_Ag_log_baseline/scheme_8_NP.txt models_Ag_log_baseline/scheme_8_NP_par.csv models_Ag_log_baseline/scheme_8_NP_prior.csv models_Ag_log_baseline/scheme_8_par.csv models_Ag_log_baseline/scheme_8_PI.txt models_Ag_log_baseline/scheme_8_PI_par.csv models_Ag_log_baseline/scheme_8_PI_prior.csv models_Ag_log_baseline/scheme_8_prior.csv models_Ag_log_baseline/scheme_9.txt models_Ag_log_baseline/scheme_9_inact.txt models_Ag_log_baseline/scheme_9_inact_NP.txt models_Ag_log_baseline/scheme_9_inact_NP_par.csv models_Ag_log_baseline/scheme_9_inact_NP_prior.csv models_Ag_log_baseline/scheme_9_inact_par.csv models_Ag_log_baseline/scheme_9_inact_PI.txt models_Ag_log_baseline/scheme_9_inact_PI_par.csv models_Ag_log_baseline/scheme_9_inact_PI_prior.csv models_Ag_log_baseline/scheme_9_inact_prior.csv models_Ag_log_baseline/scheme_9_NP.txt models_Ag_log_baseline/scheme_9_NP_par.csv models_Ag_log_baseline/scheme_9_NP_prior.csv models_Ag_log_baseline/scheme_9_par.csv models_Ag_log_baseline/scheme_9_PI.txt models_Ag_log_baseline/scheme_9_PI_par.csv models_Ag_log_baseline/scheme_9_PI_prior.csv models_Ag_log_baseline/scheme_9_prior.csv
    slurm/run_32_scheme_11_DR.sh
    slurm/run_32_scheme_9_DR.sh
    slurm/run_32_scheme_4_DR.sh
    
    slurm/run_64_CPU_8_schemes_1_2.sh
    models_Ag_log_baseline/scheme_2_PI_prior.csv models_Ag_log_baseline/scheme_3_PI_prior.csv models_Ag_log_baseline/scheme_1_PI_par.csv models_Ag_log_baseline/scheme_3_PI.txt models_Ag_log_baseline/scheme_2_PI_par.csv models_Ag_log_baseline/scheme_3_PI_par.csv models_Ag_log_baseline/scheme_1_PI_prior.csv models_Ag_log_baseline/scheme_1_PI.txt models_Ag_log_baseline/scheme_2_PI.txt
    slurm/run_64_CPU_8_schemes_3_5.sh
    slurm/run_64_CPU_16_schemes_6.sh
    slurm/run_64_CPU_8_schemes_8_9.sh
    
    slurm/run_64_CPU_16_schemes_10.sh
    slurm/run_64_CPU_32_schemes_11_DR.sh
    slurm/run_64_CPU_16_schemes_4.sh
    slurm/run_64_CPU_32_schemes_11_SNR.sh
    slurm/run_64_CPU_64_schemes_7_DR_no_inact.sh
    slurm/run_64_CPU_64_schemes_7_DR_inact.sh
    command_lines/run_cuevi.txt
    experiments/Moffatt_Hume_2007_ATP_time_idealized_segments.txt
    CLI_thermo_evidence_fraction_dts.h
    parallel_tempering_fraction.h
    report_thermo_fraction_evidence.h
    scripts/evidence_fraction_thermo_dts_data.txt
    scripts/evidence_fraction_thermo_dts_continuation.txt
    command_lines/run_fraction_dts.txt
    slurm/run_64_CPU_8_schemes_8_9_fraction.sh
    slurm/run_32_scheme_9_DR_fraction.sh
    scripts/evidence_thermo_dts_data_2.txt
    scripts/evidence_thermo_dts_continuation_2.txt
    scripts/idealize_experiment_3.txt
    slurm/run_32_scheme_10_11_DR_continuation.sh
    slurm/run_32_scheme_10_11_SNR.sh
    slurm/run_32_scheme_4_9_SNR.sh
    slurm/run_32_scheme_6_DR_SNR.sh
    slurm/run_32_scheme_7_DR_SNR.sh
    slurm/run_32_scheme_4_9_DR_continuation.sh
    slurm/run_64_CPU_16_schemes_1_2_DR_SNR.sh
    slurm/run_64_CPU_16_schemes_3_5_DR_SNR.sh
    slurm/run_64_CPU_32_schemes_4_DR_SNR_continuation.sh
    slurm/run_64_CPU_32_schemes_9_DR_SNR_continuation.sh
    slurm/run_64_CPU_32_schemes_8_DR_SNR_continuation.sh
    slurm/run_64_CPU_32_schemes_7_DR_SNR_semi_continuation.sh
    slurm/run_64_CPU_32_schemes_10_DR.sh
    command_lines/run_Simulation.txt
    command_lines/run_Simulation_full
    models_posterior/scheme_10_inact_PI_posterior_full_par.csv
    models_posterior/scheme_10_inact_PI_full.txt
    slurm/run_64_CPU_32_schemes_10_DR_sim_bessel.sh
    slurm/run_64_CPU_32_schemes_10_DR_rep_nyquist.sh
    models_posterior/scheme_9_inact_PI.txt models_posterior/scheme_9_inact_PI_posterior_par.csv models_posterior/scheme_9_inact_PI_prior.csv
    experiments/Sim_16288650277403656846_scheme_10_inact_0f15336_072300s883831_7319076573725717296_simulation.csv
    experiments/Sim_scheme_9_simulation.txt
    slurm/run_32_scheme_10_DR_nyquist.sh
    slurm/run_32_scheme_9_DR_nyquist.sh
    slurm/run_32_scheme_11_DR_nyquist.sh
    slurm/run_32_scheme_4_DR_nyquist.sh
    command_lines/run_Continuation.txt
    slurm/run_32_scheme_6_DR_nyquist.sh
    slurm/run_32_scheme_16_DR_nyquist.sh
    models_Ag_log_baseline/scheme_16_inact_PI_par.csv
    models_Ag_log_baseline/scheme_16_inact_PI_prior.csv
    models_Ag_log_baseline/scheme_16_inact_PI.txt
    models_Ag_log_baseline/scheme_16.txt
    slurm/run_EP_16_scheme_10_DR_nyquist.sh
    slurm/run_EP_16_scheme_16_DR_nyquist.sh
    models_Ag_log_baseline/scheme_16.txt
    models_Ag_log_baseline/scheme_17_inact_PI.txt
    models_Ag_log_baseline/scheme_17_inact_PI_prior.csv
    models_Ag_log_baseline/scheme_17_inact_PI_par.csv
    models_Ag_log_baseline/scheme_18_inact_PI_par.csv
    models_Ag_log_baseline/scheme_18_inact_PI_prior.csv
    models_Ag_log_baseline/scheme_18_inact_PI.txt
    slurm/run_32_scheme_17_DR_nyquist.sh
    slurm/run_32_scheme_18_DR_nyquist.sh
    slurm/run_32_scheme_8_DR_nyquist.sh
    
    
    
    
    
    
    
    
    
    
    
    
    
    
)
include_directories( ${gsl_dir})

target_link_libraries(macro_dr ${lapackblas_libraries} ${gsl_library} OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
