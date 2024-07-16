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
    
   # old_code_in_main.cpp
   
    scripts/run_evidence.txt
    
    scripts/experiment_7.txt
    scripts/simulation.txt
    scripts/evidence_simulation.txt
    scripts/evidence_data.txt
    
    multi_task/multi_task_slurm.sh
    slurm_Serafin/s_multi_debug.sh
    multi_task/simulation_4.sh
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
    multi_task/MH_Li_old.sh
    models/scheme_2_inact.txt models/scheme_3_inact.txt models/scheme_4_inact.txt
    slurm_Serafin/s_MH_Li.sh
    slurm_Serafin/s_MH_Li_sim.sh
    slurm_Serafin/s_MH_Li_short.sh
    experiments/Moffatt_Hume_2007_ATP_time_segments.txt experiments/Moffatt_Hume_2007_ATP_time_recording.txt
    slurm_tupac_old/s_MH_Li7.sh slurm_tupac_old/s_MH_Li7_sim.sh slurm_tupac_old/s_MH_Li_sim.sh slurm_tupac_old/s_MH_Li.sh slurm_tupac_old/s_MH_Li_short.sh
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
    slurm_tupac_old/s_MH_Li_4.sh
    slurm_tupac_old/s_MH_Li_4_sim.sh
    slurm_tupac_old/s_MH_Li_2.sh
    slurm_tupac_old/s_MH_Li_2_sim.sh
    scripts/likelihood_NR.txt
    models/scheme_4_inact_NP_par.csv models/scheme_4_inact_NP_prior.csv models/scheme_4_NP_par.csv models/scheme_4_NP_prior.csv
    scripts/load_simulation.txt
    experiments/Sim_scheme_4_NP_N_simulation.csv
    slurm_tupac_old/s_v10_MacroNR.sh
    models/scheme_4_NP.txt
    models/scheme_4_inact_NP.txt
    slurm_tupac_old/s_v10_MacroDRPD.sh slurm_tupac_old/s_v10_MacroDRP.sh slurm_tupac_old/s_v10_MacroNRPD.sh 
    slurm_tupac_old/s_v10_MacroNRP.sh  slurm_tupac_old/s_v10_MacroDR.sh
    gsl_integrate.h
    slurm_tupac_old/s_v12_MacroDRPD.sh
    slurm_tupac_old/s_v12_MacroNRPD.sh
    scripts/load_simulation_7.txt
    experiments/Sim_6465237192204259400_scheme_4_inact_NP_N_simulation.csv
    slurm_tupac_old/s_v13_4_sim_DR.sh
    slurm_tupac_old/s_v13_4_DR.sh 
    slurm_tupac_old/s_v13_4_NR.sh 
    slurm_tupac_old/s_v13_4_sim_NR.sh
    slurm_tupac_old/s_v13_4_DRPD.sh 
    slurm_tupac_old/s_v13_4_sim_NRPD.sh 
    slurm_tupac_old/s_v13_4_sim_DRPD.sh
    slurm_tupac_old/s_v13_4_sim_NRPD.sh
    slurm_tupac_old/s_v13_4_NRPD.sh
    scripts/idealize_experiment.txt
    scripts/idealize_experiment_7.txt
    scripts/run_idealize_experiment.txt
    slurm_tupac_old/s_v15_4_DR.sh
    slurm_tupac_old/s_v15_4_sim_DR.sh
    experiments/Moffatt_Hume_2007_ATP_time_idealized.txt
    slurm_tupac_old/s_v16_DR_no_cores.sh
    slurm_tupac_old/number_cores/s_1core.sh
    slurm_tupac_old/number_cores/s_64core.sh slurm_tupac_old/number_cores/s_16core.sh slurm_tupac_old/number_cores/s_32core.sh slurm_tupac_old/number_cores/s_2core.sh slurm_tupac_old/number_cores/s_4core.sh slurm_tupac_old/number_cores/s_8core.sh
    slurm_tupac_old/number_cores/run_cores.sh
    scripts/evidence_thermo_continuation.txt
    slurm_tupac_old/cross_model/simulated_cross_model_comparison.sh
    slurm_tupac/MHLi.sh
    slurm/slurm.sh
    
    multi_task/MH_Li.sh
    multi_task/MH_Li_continuation_1.sh
    multi_task/run_MH_Li.sh
    
    slurm_tupac/run_MHLi.sh
    slurm_tupac/MHLi_continuation.sh
    multi_task/MH_Li_IE.sh
    multi_task/MH_Li_IE_continuation.sh
    slurm_tupac/MHLi_IE_continuation.sh
    slurm_tupac/MHLi_IE.sh
    slurm_tupac/run_MHLi_IE.sh
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
    multi_task/MH_Al.sh
    multi_task/MH_Al_continuation_1.sh
    multi_task/MH_Al_IE.sh
    multi_task/MH_Al_IE_continuation.sh
    slurm_tupac/MHAl.sh
    slurm_tupac/MHAl_IE.sh
    slurm_tupac/MHAl_continuation.sh
    slurm_tupac/MHAl_IE_continuation.sh
    slurm_tupac/run_MHAl.sh
    slurm_tupac/run_MHAl_IE.sh
    multi_task/MH_Li34.sh
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
    multi_task/M_s2-s10.sh
    scripts/beta_32.txt
    scripts/beta_64.txt
    multi_task/M_s2-s10_continuation.sh
    slurm_tupac/M_s2-s10.sh
    slurm_tupac/M_s2-s10_continuation.sh
    slurm_tupac/run_M_s2-s10.sh
    slurm_tupac/run_M_s2-s10_continuation.sh
    multi_task/M_scheme_N.sh
    slurm_tupac/M_scheme_N.sh
    slurm_tupac/run_M_scheme_N.sh
    slurm_Serafin/run_M_scheme_N.sh
    slurm_Serafin/M_scheme_N.sh
    models_Ag/scheme_5_NP.txt models_Ag/scheme_5_NP_par.csv models_Ag/scheme_5_PI_prior.csv models_Ag/scheme_5_NP_prior.csv models_Ag/scheme_5_inact_NP_prior.csv models_Ag/scheme_5_inact_NP_par.csv models_Ag/scheme_5_inact_NP.txt models_Ag/scheme_5_PI_par.csv models_Ag/scheme_5_inact_PI_prior.csv models_Ag/scheme_5_inact_PI_par.csv models_Ag/scheme_5_inact_PI.txt models_Ag/scheme_5_PI.txt models_Ag/scheme_5_inact.txt models_Ag/scheme_5.txt models_Ag/scheme_5_inact_par.csv models_Ag/scheme_5_inact_prior.csv models_Ag/scheme_5_par.csv models_Ag/scheme_5_prior.csv models_Ag/scheme_6_inact_NP_par.csv models_Ag/scheme_7_NP_par.csv models_Ag/scheme_6_NP_par.csv models_Ag/scheme_7_inact_NP_par.csv models_Ag/scheme_6_NP.txt models_Ag/scheme_7_inact_NP.txt models_Ag/scheme_6_inact_NP.txt models_Ag/scheme_7_NP.txt models_Ag/scheme_4_PI.txt models_Ag/scheme_2_inact_PI.txt models_Ag/scheme_1_inact_PI.txt models_Ag/scheme_4_inact_PI.txt models_Ag/scheme_3_inact_PI.txt models_Ag/scheme_6_NP_prior.csv models_Ag/scheme_6_inact_NP_prior.csv models_Ag/scheme_7_inact_NP_prior.csv models_Ag/scheme_7_NP_prior.csv models_Ag/scheme_2_inact_PI_par.csv models_Ag/scheme_3_inact_PI_prior.csv models_Ag/scheme_1_inact_PI_prior.csv models_Ag/scheme_2_inact_PI_prior.csv models_Ag/scheme_4_PI_prior.csv models_Ag/scheme_4_inact_PI_prior.csv models_Ag/scheme_3_inact_PI_par.csv models_Ag/scheme_4_inact_PI_par.csv models_Ag/scheme_4_PI_par.csv models_Ag/scheme_1_inact_PI_par.csv models_Ag/scheme_7_inact_PI.txt models_Ag/scheme_6_PI.txt models_Ag/scheme_6_inact_PI.txt models_Ag/scheme_6_inact_PI_par.csv models_Ag/scheme_6_inact_PI_prior.csv models_Ag/scheme_6_PI_par.csv models_Ag/scheme_6_PI_prior.csv models_Ag/scheme_7_inact_PI_par.csv models_Ag/scheme_7_inact_PI_prior.csv models_Ag/scheme_7_PI_par.csv models_Ag/scheme_7_PI_prior.csv models_Ag/scheme_7_PI.txt models_Ag/scheme_7_inact_par.csv models_Ag/scheme_7_inact_prior.csv models_Ag/scheme_7_par.csv models_Ag/scheme_7_prior.csv models_Ag/scheme_7_inact.txt models_Ag/scheme_7.txt models_Ag/scheme_6_inact_par.csv models_Ag/scheme_6_inact_prior.csv models_Ag/scheme_6_inact.txt models_Ag/scheme_6_par.csv models_Ag/scheme_6_prior.csv models_Ag/scheme_6.txt models_Ag/scheme_3.txt models_Ag/scheme_2_inact_NP.txt models_Ag/scheme_3_inact_NP.txt models_Ag/scheme_1_inact.txt models_Ag/scheme_1_inact_prior.csv models_Ag/scheme_3_inact_par.csv models_Ag/scheme_1_par.csv models_Ag/scheme_2_inact_par.csv models_Ag/scheme_2_inact_prior.csv models_Ag/scheme_1_prior.csv models_Ag/scheme_3_inact_prior.csv models_Ag/scheme_1_inact_par.csv models_Ag/scheme_4_inact_NP_prior.csv models_Ag/scheme_4_inact_NP_par.csv models_Ag/scheme_4_NP_prior.csv models_Ag/scheme_4_NP_par.csv models_Ag/scheme_4_inact_prior.csv models_Ag/scheme_4_inact_par.csv models_Ag/scheme_1_inact_NP.txt models_Ag/scheme_3_inact_NP_prior.csv models_Ag/scheme_3_inact_NP_par.csv models_Ag/scheme_2_inact_NP_prior.csv models_Ag/scheme_2_inact_NP_par.csv models_Ag/scheme_1_par_copy.csv models_Ag/scheme_1_inact_NP_prior.csv models_Ag/scheme_1_inact_NP_par.csv models_Ag/scheme_1.txt models_Ag/scheme_4_inact_NP.txt models_Ag/scheme_4_inact.txt models_Ag/scheme_4_NP.txt models_Ag/scheme_3_inact.txt models_Ag/scheme_2_inact.txt models_Ag/scheme_2.txt models_Ag/scheme_4.txt models_Ag/scheme_4_prior.csv models_Ag/scheme_4_par.csv models_Ag/scheme_2_prior.csv models_Ag/scheme_2_par.csv models_Ag/scheme_3_prior.csv models_Ag/scheme_3_par.csv
    models_Ag/scheme_8.txt models_Ag/scheme_8_inact.txt models_Ag/scheme_8_inact_NP.txt models_Ag/scheme_8_inact_NP_par.csv models_Ag/scheme_8_inact_NP_prior.csv models_Ag/scheme_8_inact_par.csv models_Ag/scheme_8_inact_PI.txt models_Ag/scheme_8_inact_PI_par.csv models_Ag/scheme_8_inact_PI_prior.csv models_Ag/scheme_8_inact_prior.csv models_Ag/scheme_8_NP.txt models_Ag/scheme_8_NP_par.csv models_Ag/scheme_8_NP_prior.csv models_Ag/scheme_8_par.csv models_Ag/scheme_8_PI.txt models_Ag/scheme_8_PI_par.csv models_Ag/scheme_8_PI_prior.csv models_Ag/scheme_8_prior.csv models_Ag/scheme_9.txt models_Ag/scheme_9_inact.txt models_Ag/scheme_9_inact_NP.txt models_Ag/scheme_9_inact_NP_par.csv models_Ag/scheme_9_inact_NP_prior.csv models_Ag/scheme_9_inact_par.csv models_Ag/scheme_9_inact_PI.txt models_Ag/scheme_9_inact_PI_par.csv models_Ag/scheme_9_inact_PI_prior.csv models_Ag/scheme_9_inact_prior.csv models_Ag/scheme_9_NP.txt models_Ag/scheme_9_NP_par.csv models_Ag/scheme_9_NP_prior.csv models_Ag/scheme_9_par.csv models_Ag/scheme_9_PI.txt models_Ag/scheme_9_PI_par.csv models_Ag/scheme_9_PI_prior.csv models_Ag/scheme_9_prior.csv models_Ag/scheme_10.txt models_Ag/scheme_10_inact.txt models_Ag/scheme_10_inact_NP.txt models_Ag/scheme_10_inact_NP_par.csv models_Ag/scheme_10_inact_NP_prior.csv models_Ag/scheme_10_inact_par.csv models_Ag/scheme_10_inact_PI.txt models_Ag/scheme_10_inact_PI_par.csv models_Ag/scheme_10_inact_PI_prior.csv models_Ag/scheme_10_inact_prior.csv models_Ag/scheme_10_NP.txt models_Ag/scheme_10_NP_par.csv models_Ag/scheme_10_NP_prior.csv models_Ag/scheme_10_par.csv models_Ag/scheme_10_PI.txt models_Ag/scheme_10_PI_par.csv models_Ag/scheme_10_PI_prior.csv models_Ag/scheme_10_prior.csv
    models/scheme_11_inact.txt models/scheme_11_inact_NP.txt models/scheme_11_inact_NP_par.csv models/scheme_11_inact_NP_prior.csv models/scheme_11_inact_par.csv models/scheme_11_inact_PI.txt models/scheme_11_inact_PI_par.csv models/scheme_11_inact_PI_prior.csv models/scheme_11_inact_prior.csv models/scheme_11_NP.txt models/scheme_11_NP_par.csv models/scheme_11_NP_prior.csv models/scheme_11_par.csv models/scheme_11_PI.txt models/scheme_11_PI_par.csv models/scheme_11_PI_prior.csv models/scheme_11_prior.csv
    models_Ag/scheme_11.txt models_Ag/scheme_11_inact.txt models_Ag/scheme_11_inact_NP.txt models_Ag/scheme_11_inact_NP_par.csv models_Ag/scheme_11_inact_NP_prior.csv models_Ag/scheme_11_inact_par.csv models_Ag/scheme_11_inact_PI.txt models_Ag/scheme_11_inact_PI_par.csv models_Ag/scheme_11_inact_PI_prior.csv models_Ag/scheme_11_inact_prior.csv models_Ag/scheme_11_NP.txt models_Ag/scheme_11_NP_par.csv models_Ag/scheme_11_NP_prior.csv models_Ag/scheme_11_par.csv models_Ag/scheme_11_PI.txt models_Ag/scheme_11_PI_par.csv models_Ag/scheme_11_PI_prior.csv models_Ag/scheme_11_prior.csv
    models/scheme_15.txt models/scheme_15_inact.txt models/scheme_15_inact_par.csv models/scheme_15_inact_prior.csv models/scheme_15_par.csv models/scheme_15_prior.csv
    models/scheme_15.txt models/scheme_15_inact.txt models/scheme_15_inact_NP.txt models/scheme_15_inact_NP_par.csv models/scheme_15_inact_NP_prior.csv models/scheme_15_inact_par.csv models/scheme_15_inact_PI.txt models/scheme_15_inact_PI_par.csv models/scheme_15_inact_PI_prior.csv models/scheme_15_inact_prior.csv models/scheme_15_NP.txt models/scheme_15_NP_par.csv models/scheme_15_NP_prior.csv models/scheme_15_par.csv models/scheme_15_PI.txt models/scheme_15_PI_par.csv models/scheme_15_PI_prior.csv models/scheme_15_prior.csv
    models_Ag/scheme_15.txt models_Ag/scheme_15_inact.txt models_Ag/scheme_15_inact_NP.txt models_Ag/scheme_15_inact_NP_par.csv models_Ag/scheme_15_inact_NP_prior.csv models_Ag/scheme_15_inact_par.csv models_Ag/scheme_15_inact_PI.txt models_Ag/scheme_15_inact_PI_par.csv models_Ag/scheme_15_inact_PI_prior.csv models_Ag/scheme_15_inact_prior.csv models_Ag/scheme_15_NP.txt models_Ag/scheme_15_NP_par.csv models_Ag/scheme_15_NP_prior.csv models_Ag/scheme_15_par.csv models_Ag/scheme_15_PI.txt models_Ag/scheme_15_PI_par.csv models_Ag/scheme_15_PI_prior.csv models_Ag/scheme_15_prior.csv
    slurm_tupac/run_M_scheme_11.sh
    slurm_Serafin/run_M_scheme_11.sh
    slurm_Serafin/run_M_scheme_15.sh
    slurm_tupac/run_M_scheme_15.sh
    models/scheme_12.txt models/scheme_12_inact.txt models/scheme_12_inact_NP.txt models/scheme_12_inact_NP_par.csv models/scheme_12_inact_NP_prior.csv models/scheme_12_inact_par.csv models/scheme_12_inact_PI.txt models/scheme_12_inact_PI_par.csv models/scheme_12_inact_PI_prior.csv models/scheme_12_inact_prior.csv models/scheme_12_NP.txt models/scheme_12_NP_par.csv models/scheme_12_NP_prior.csv models/scheme_12_par.csv models/scheme_12_PI.txt models/scheme_12_PI_par.csv models/scheme_12_PI_prior.csv models/scheme_12_prior.csv models/scheme_13.txt models/scheme_13_inact.txt models/scheme_13_inact_NP.txt models/scheme_13_inact_NP_par.csv models/scheme_13_inact_NP_prior.csv models/scheme_13_inact_par.csv models/scheme_13_inact_PI.txt models/scheme_13_inact_PI_par.csv models/scheme_13_inact_PI_prior.csv models/scheme_13_inact_prior.csv models/scheme_13_NP.txt models/scheme_13_NP_par.csv models/scheme_13_NP_prior.csv models/scheme_13_par.csv models/scheme_13_PI.txt models/scheme_13_PI_par.csv models/scheme_13_PI_prior.csv models/scheme_13_prior.csv models/scheme_14.txt models/scheme_14_inact.txt models/scheme_14_inact_NP.txt models/scheme_14_inact_NP_par.csv models/scheme_14_inact_NP_prior.csv models/scheme_14_inact_par.csv models/scheme_14_inact_PI.txt models/scheme_14_inact_PI_par.csv models/scheme_14_inact_PI_prior.csv models/scheme_14_inact_prior.csv models/scheme_14_NP.txt models/scheme_14_NP_par.csv models/scheme_14_NP_prior.csv models/scheme_14_par.csv models/scheme_14_PI.txt models/scheme_14_PI_par.csv models/scheme_14_PI_prior.csv models/scheme_14_prior.csv
    models_Ag/scheme_12.txt models_Ag/scheme_12_inact.txt models_Ag/scheme_12_inact_NP.txt models_Ag/scheme_12_inact_NP_par.csv models_Ag/scheme_12_inact_NP_prior.csv models_Ag/scheme_12_inact_par.csv models_Ag/scheme_12_inact_PI.txt models_Ag/scheme_12_inact_PI_par.csv models_Ag/scheme_12_inact_PI_prior.csv models_Ag/scheme_12_inact_prior.csv models_Ag/scheme_12_NP.txt models_Ag/scheme_12_NP_par.csv models_Ag/scheme_12_NP_prior.csv models_Ag/scheme_12_par.csv models_Ag/scheme_12_PI.txt models_Ag/scheme_12_PI_par.csv models_Ag/scheme_12_PI_prior.csv models_Ag/scheme_12_prior.csv models_Ag/scheme_13.txt models_Ag/scheme_13_inact.txt models_Ag/scheme_13_inact_NP.txt models_Ag/scheme_13_inact_NP_par.csv models_Ag/scheme_13_inact_NP_prior.csv models_Ag/scheme_13_inact_par.csv models_Ag/scheme_13_inact_PI.txt models_Ag/scheme_13_inact_PI_par.csv models_Ag/scheme_13_inact_PI_prior.csv models_Ag/scheme_13_inact_prior.csv models_Ag/scheme_13_NP.txt models_Ag/scheme_13_NP_par.csv models_Ag/scheme_13_NP_prior.csv models_Ag/scheme_13_par.csv models_Ag/scheme_13_PI.txt models_Ag/scheme_13_PI_par.csv models_Ag/scheme_13_PI_prior.csv models_Ag/scheme_13_prior.csv models_Ag/scheme_14.txt models_Ag/scheme_14_inact.txt models_Ag/scheme_14_inact_NP.txt models_Ag/scheme_14_inact_NP_par.csv models_Ag/scheme_14_inact_NP_prior.csv models_Ag/scheme_14_inact_par.csv models_Ag/scheme_14_inact_PI.txt models_Ag/scheme_14_inact_PI_par.csv models_Ag/scheme_14_inact_PI_prior.csv models_Ag/scheme_14_inact_prior.csv models_Ag/scheme_14_NP.txt models_Ag/scheme_14_NP_par.csv models_Ag/scheme_14_NP_prior.csv models_Ag/scheme_14_par.csv models_Ag/scheme_14_PI.txt models_Ag/scheme_14_PI_par.csv models_Ag/scheme_14_PI_prior.csv models_Ag/scheme_14_prior.csv
    slurm_Serafin/run_M_scheme_12sh
    slurm_Serafin/run_M_scheme_13.sh
    slurm_Serafin/run_M_scheme_14.sh
    slurm_tupac/run_M_scheme_12.sh
    slurm_tupac/run_M_scheme_13.sh
    slurm_tupac/run_M_scheme_14.sh
    multi_task/M_scheme_N_2_tasks.sh
    multi_task/M_scheme_N_4_tasks.sh
    multi_task/M_scheme_N_8_tasks.sh
    multi_task/M_scheme_N_16_tasks.sh
    slurm_tupac/M_scheme_N_tasks.sh
    slurm_tupac/run_M_scheme_N_2_tasks.sh
    slurm_tupac/run_M_scheme_N_4_tasks.sh
    slurm_tupac/run_M_scheme_N_8_tasks.sh
    slurm_tupac/run_M_scheme_N_16_tasks.sh
    parallel_levenberg_tempering.h
    CLI_thermo_levenberg_evidence.h
    scripts/evidence_levenberg_data.txt
    experiments/Moffatt_Hume_2007_ATP_time_idealized_2.txt
    scripts/idealize_experiment_2.txt
    scripts/beta_16.txt
    models/scheme_1_inact_PI_prior_0001.csv
    models/scheme_1_inact_PI_0001.txt
    slurm_tupac/run_M_schemes_1_2.sh
    slurm_tupac/run_M_schemes_3_4_5.sh
    slurm_tupac/run_M_schemes_6_8.sh
    slurm_tupac/run_M_schemes_9_10.sh
    slurm_tupac/run_M_scheme_7.sh
    slurm_Serafin/M_scheme_N_tasks.sh
    scripts/calculate_likelihood.txt
    slurm_tupac/run_M_schemes_11_12.sh
    slurm_tupac/run_M_schemes_13_14.sh
    slurm_tupac/run_M_schemes_7_15.sh
    slurm_Serafin/run_M_scheme_11_12.sh
    slurm_Serafin/run_M_scheme_13_14.sh
    slurm_Serafin/run_M_scheme_7_15.sh
    experiments/Moffatt_Hume_2007_ATP_time_idealized_old.txt
    slurm_Serafin/run_M_scheme_7.sh
    CLI_thermo_evidence_dts.h
    
    
    
    
    
    
    
    
    
    
    
    
    
)
include_directories( ${gsl_dir})

target_link_libraries(macro_dr ${lapackblas_libraries} ${gsl_library} OpenMP::OpenMP_CXX)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
