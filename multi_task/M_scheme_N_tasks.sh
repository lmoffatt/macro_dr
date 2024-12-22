#!/bin/bash


SCHEME_DIRS=(${SCHEME_DIR_0}  ${SCHEME_DIR_1}  ${SCHEME_DIR_2}  ${SCHEME_DIR_3} ${SCHEME_DIR_4}  ${SCHEME_DIR_5}  ${SCHEME_DIR_6}  ${SCHEME_DIR_7} ${SCHEME_DIR_8}  ${SCHEME_DIR_9}  ${SCHEME_DIR_10}  ${SCHEME_DIR_11} ${SCHEME_DIR_12}  ${SCHEME_DIR_13}  ${SCHEME_DIR_14}  ${SCHEME_DIR_15} )
SCHEME_FILES=( ${SCHEME_0}  ${SCHEME_1} ${SCHEME_2}  ${SCHEME_3} ${SCHEME_4}  ${SCHEME_5} ${SCHEME_6}  ${SCHEME_7} ${SCHEME_8}  ${SCHEME_9} ${SCHEME_10}  ${SCHEME_11} ${SCHEME_12}  ${SCHEME_13} ${SCHEME_14}  ${SCHEME_15}  )

PATH_MACRO_DRS=( ${PATH_MACRO_DR_0} ${PATH_MACRO_DR_1}  ${PATH_MACRO_DR_2} ${PATH_MACRO_DR_3} ${PATH_MACRO_DR_4} ${PATH_MACRO_DR_5}  ${PATH_MACRO_DR_6} ${PATH_MACRO_DR_7} ${PATH_MACRO_DR_8} ${PATH_MACRO_DR_9}  ${PATH_MACRO_DR_10} ${PATH_MACRO_DR_11} ${PATH_MACRO_DR_12} ${PATH_MACRO_DR_13}  ${PATH_MACRO_DR_14} ${PATH_MACRO_DR_15}  )

EXPERIMENTS=(  ${EXPERIMENT_0} ${EXPERIMENT_1}  ${EXPERIMENT_2} ${EXPERIMENT_3}  ${EXPERIMENT_4} ${EXPERIMENT_5}  ${EXPERIMENT_6} ${EXPERIMENT_7} ${EXPERIMENT_8} ${EXPERIMENT_9}  ${EXPERIMENT_10} ${EXPERIMENT_11}  ${EXPERIMENT_12} ${EXPERIMENT_13}  ${EXPERIMENT_14} ${EXPERIMENT_15}  )


if [ -n "${LIK_0}" ];   then                                         
    LIKS=( ${LIK_0} ${LIK_1}  ${LIK_2} ${LIK_3}  ${LIK_4} ${LIK_5}  ${LIK_6} ${LIK_7} ${LIK_8} ${LIK_9}  ${LIK_10} ${LIK_11}  ${LIK_12} ${LIK_13}  ${LIK_14} ${LIK_15}  )
fi

if [ -n "${IDNAME_0}" ]; then                                           
    IDNAMES=( ${IDNAME_0} ${IDNAME_1}  ${IDNAME_2} ${IDNAME_3}  ${IDNAME_4} ${IDNAME_5}  ${IDNAME_6} ${IDNAME_7} ${IDNAME_8} ${IDNAME_9}  ${IDNAME_10} ${IDNAME_11}  ${IDNAME_12} ${IDNAME_13}  ${IDNAME_14} ${IDNAME_15}  )
fi

if [ -n "${CONTINUATION_NUMBER_0}" ]; then                                           
    CONTINUATION_NUMBERS=( ${CONTINUATION_NUMBER_0CONTINUATION_NUMBER_0} ${CONTINUATION_NUMBER_1}  ${CONTINUATION_NUMBER_2} ${CONTINUATION_NUMBER_3}  ${CONTINUATION_NUMBER_4} )
fi



SCHEME_DIR=${SCHEME_DIRS[$SLURM_LOCALID]}

SCHEME=${SCHEME_FILES[$SLURM_LOCALID]}

PATH_MACRO_DR=${PATH_MACRO_DRS[$SLURM_LOCALID]}

EXPERIMENT=${EXPERIMENTS[$SLURM_LOCALID]}


if [ -n "${LIKS}" ]; then                                           
   LIK=${LIKS[$SLURM_LOCALID]}
fi  

if [ -z "${CP}" ]; then                                           
  CP=$SLURM_CPUS_PER_TASK
fi  



LIKELIHOOD=$([ -z "${LIK}" ] && echo "likelihood" || echo "likelihood_"${LIK} )

LIK_AB=$([ -z "${LIK}" ] && echo "" || echo "_"${LIK} )




EXPER_ABR=$([ "$EXPERIMENT" = "idealize_experiment_2" ] && echo "_IE" || echo "")

SCH_ABR=$([ "$SCHEME_DIR" = "models_Ag" ] && echo "_Ag" ||[ "$SCHEME_DIR" = "models_Ag_log_baseline" ] && echo "_logbaseline" || echo "")

LOCAL_ID=$([ "$USE_LOCAL_ID" = 1 ] && echo ${SLURM_LOCALID} || echo "")


if [ -n "${IDNAMES}" ];  then                                          
   IDNAME=${IDNAMES[$SLURM_LOCALID]}
else
   IDNAME=${PATH_MACRO_DR}${EXPER_ABR}${LIK_AB}_${SLURM_CPUS_PER_TASK}c_${N_SCOUTS}s_${N_BETA}b_${SCHEME}${SCH_ABR}_${LOCAL_ID}_0   
fi   

if [ -n "${CONTINUATION_NUMBERS}" ];  then                                          
   CONTINUATION_NUMBER=${CONTINUATION_NUMBERS[$SLURM_LOCALID]}
fi   




cd ${PATH_MACRO}/${WORKING_DIRECTORY}
export BLIS_NUM_THREADS=$SLURM_CPUS_PER_TASK


if [ "$CONTINUATION_NUMBER" = 0 ]; then
    ${PATH_MACRO}/${PATH_MACRO_DRX}/macro_dr ${PATH_MACRO}/macro_dr/${SCHEME_DIR}/${SCHEME}.txt ${PATH_MACRO}/macro_dr/scripts/${EXPERIMENT}.txt  ${PATH_MACRO}/macro_dr/scripts/simulation.txt "--runIdName= \"${IDNAME}\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ${PATH_MACRO}/macro_dr/scripts/${LIKELIHOOD}.txt ${PATH_MACRO}/macro_dr/scripts/beta_${N_BETA}.txt  ${PATH_MACRO}/macro_dr/scripts/evidence_${EVIDENCE_ALGORITHM}_data.txt 
else
${PATH_MACRO}/${PATH_MACRO_DRX}/macro_dr ${PATH_MACRO}/macro_dr/${SCHEME_DIR}/${SCHEME}.txt ${PATH_MACRO}/macro_dr/scripts/${EXPERIMENT}.txt  ${PATH_MACRO}/macro_dr/scripts/simulation.txt "--runIdName= \"${IDNAME}\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ${PATH_MACRO}/macro_dr/scripts/${LIKELIHOOD}.txt ${PATH_MACRO}/macro_dr/scripts/beta_${N_BETA}.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})" ${PATH_MACRO}/macro_dr/scripts/evidence_${EVIDENCE_ALGORITHM}_continuation.txt 
fi





