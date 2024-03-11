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
#set(CMAKE_BUILD_TYPE "Release")

set (CMAKE_CXX_FLAGS "-ftemplate-backtrace-limit=20 -Werror=return-type -fconcepts-diagnostics-depth=4 -fmax-errors=5 -ftime-report" )


#find_package(Catch2 REQUIRED)
#target_link_libraries(tests Catch2::Catch2)

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
    slurm_Serafin/s_first
    experiments/Moffatt_Hume_2007_ATP_time.txt
    experiments/Moffatt_Hume_2007_ATP_time_7.txt
    experiments/Moffatt_Hume_2007_ATP_time_7_constant_dt.txt
    slurm_Serafin/s_debug.sh
    slurm_Serafin/s_long.sh
    continuation.h
    models.h
   # old_code_in_main.cpp
    CLI_macrodr_commands.h
    scripts/run_evidence.txt
    models/model6_prior.csv models/model6_parameters.csv
    
    scripts/experiment_7.txt
    scripts/simulation.txt
    scripts/evidence_simulation.txt
    scripts/evidence_data.txt
    
    models/model9_prior.csv models/model9_par.csv models/model8_prior.csv models/model8_par.csv models/model7_prior.csv models/model7_par.csv models/model6_Eff_no_inactivation_prior.csv models/model6_Eff_no_inactivation_par.csv models/model6_no_inactivation_prior.csv models/model6_no_inactivation_par.csv models/model6_prior.csv models/model6_par.csv models/model4_g_lin_prior.csv models/model4_g_lin_par.csv models/model4_prior.csv models/model4_par.csv models/model01_prior.csv models/model01_par.csv models/model00_7_prior.csv models/model00_7_par.csv models/model00_prior.csv models/model00_par.csv
    multi_task/multi_task_slurm.sh
    slurm_Serafin/s_multi_debug.sh
    multi_task/simulation_4.sh
    models/model00_7.txt
    models/model6.txt models/model9.txt models/model8.txt models/model7.txt models/model4.txt models/model4_g_lin.txt models/model01.txt models/model00.txt models/model6_no_inactivation.txt models/model6_Eff_no_inactivation.txt models/model00_7.txt
    models_MoffattHume_linear.h
    models_MoffattHume_allosteric.h
    multi_task/MH_Li7.sh
    models/scheme_1.txt
    models/scheme_2.txt
    models/scheme_3.txt
    models/scheme_4.txt
    models/scheme_4_prior.csv models/scheme_4_par.csv models/scheme_4_inact_prior.csv models/scheme_4_inact_par.csv models/scheme_3_prior.csv models/scheme_3_par.csv models/scheme_3_inact_prior.csv models/scheme_3_inact_par.csv models/scheme_2_prior.csv models/scheme_2_par.csv models/scheme_2_inact_prior.csv models/scheme_2_inact_par.csv models/scheme_1_prior.csv models/scheme_1_par.csv models/scheme_1_inact_prior.csv models/scheme_1_inact_par.csv
    experiments/Moffatt_Hume_2007_ATP_time_7_recording.txt
    scripts/experiment.txt
    experiments/Moffatt_Hume_2007_ATP_time_7_segments.txt
    experiments/Moffatt_Hume_2007_ATP_time_segments.txt
    multi_task/MH_Li7_sim.sh
    slurm_Serafin/s_MH_Li7_sim.sh
    slurm_Serafin/s_MH_Li7.sh
    multi_task/MH_Li_sim.sh
    multi_task/MH_Li.sh
    models/scheme_2_inact.txt models/scheme_3_inact.txt models/scheme_4_inact.txt
    slurm_Serafin/s_MH_Li.sh
    slurm_Serafin/s_MH_Li_sim.sh
    slurm_Serafin/s_MH_Li_short.sh
    experiments/Moffatt_Hume_2007_ATP_time_segments.txt experiments/Moffatt_Hume_2007_ATP_time_recording.txt
    slurm_tupac/s_MH_Li7.sh slurm_tupac/s_MH_Li7_sim.sh slurm_tupac/s_MH_Li_sim.sh slurm_tupac/s_MH_Li.sh slurm_tupac/s_MH_Li_short.sh
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
    single_task/MH_Li_4.sh
    slurm_tupac/s_MH_Li_4.sh
    slurm_tupac/s_MH_Li_4_sim.sh
    slurm_tupac/s_MH_Li_2.sh
    slurm_tupac/s_MH_Li_2_sim.sh
    scripts/likelihood_NR.txt
    models/scheme_4_inact_NP_par.csv models/scheme_4_inact_NP_prior.csv models/scheme_4_NP_par.csv models/scheme_4_NP_prior.csv
    scripts/load_simulation.txt
    experiments/Sim_scheme_4_NP_N_simulation.csv
    slurm_tupac/s_v10_MacroNR.sh
    models/scheme_4_NP.txt
    models/scheme_4_inact_NP.txt
    slurm_tupac/s_v10_MacroDRPD.sh slurm_tupac/s_v10_MacroDRP.sh slurm_tupac/s_v10_MacroNRPD.sh 
    slurm_tupac/s_v10_MacroNRP.sh  slurm_tupac/s_v10_MacroDR.sh
    gsl_integrate.h
    slurm_tupac/s_v12_MacroDRPD.sh
    slurm_tupac/s_v12_MacroNRPD.sh
    scripts/load_simulation_7.txt
    experiments/Sim_6465237192204259400_scheme_4_inact_NP_N_simulation.csv
    slurm_tupac/s_v13_4_sim_DR.sh
    slurm_tupac/s_v13_4_DR.sh 
    slurm_tupac/s_v13_4_NR.sh 
    slurm_tupac/s_v13_4_sim_NR.sh
    slurm_tupac/s_v13_4_DRPD.sh 
    slurm_tupac/s_v13_4_sim_NRPD.sh 
    slurm_tupac/s_v13_4_sim_DRPD.sh
    slurm_tupac/s_v13_4_sim_NRPD.sh
    slurm_tupac/s_v13_4_NRPD.sh
    scripts/idealize_experiment.txt
    scripts/idealize_experiment_7.txt
    scripts/run_idealize_experiment.txt
    slurm_tupac/s_v15_4_DR.sh
    slurm_tupac/s_v15_4_sim_DR.sh
    experiments/Moffatt_Hume_2007_ATP_time_idealized.txt 
    
    
    
)
include_directories( ${gsl_dir})

target_link_libraries(macro_dr ${lapackblas_libraries} ${gsl_library} OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
