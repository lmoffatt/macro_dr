#!/bin/bash



SCHEME_DIRS=(${SCHEME_DIR_0}  ${SCHEME_DIR_1} )
SCHEME_FILES=( ${SCHEME_0}  ${SCHEME_1}  )

PATH_MACRO_DRS=( ${PATH_MACRO_DR_0} ${PATH_MACRO_DR_1}  )

EXPERIMENTS=(  ${EXPERIMENT_0} ${EXPERIMENT_1} )


SCHEME_DIR=${SCHEME_DIRS[$SLURM_LOCALID]}

SCHEME=${SCHEME_FILES[$SLURM_LOCALID]}

PATH_MACRO_DR=${PATH_MACRO_DRS[$SLURM_LOCALID]}

EXPERIMENT=${EXPERIMENTS[$SLURM_LOCALID]}

EXPER_ABR=$([ "$EXPERIMENT" = "idealize_experiment" ] && echo "_IE" || echo "")

SCH_ABR=$([ "$SCHEME_DIR" = "models_Ag" ] && echo "_Ag" || echo "")

LOCAL_ID=$([ "$USE_LOCAL_ID" = 1 ] && echo _${SLURM_LOCALID} || echo "")

CP=$SLURM_CPUS_PER_TASK



if [ "$CONTINUATION_NUMBER" = 0 ]; then
    ${PATH_MACRO}/${PATH_MACRO_DR}/macro_dr ../macro_dr/${SCHEME_DIR}/${SCHEME}.txt ../macro_dr/scripts/${EXPERIMENT}.txt  ../macro_dr/   scripts/simulation.txt "--runIdName= \"${PATH_MACRO_DR}${EXPER_ABR}_${CP}c_${N_SCOUTS}s_${N_BETA}b_${SCHEME}${SCH_ABR}${LOCAL_ID}_0\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/beta_${N_BETA}.txt  ../macro_dr/scripts/evidence_thermo_data.txt 
else
${PATH_MACRO}/${PATH_MACRO_DR}/macro_dr ../macro_dr/${SCHEME_DIR}/${SCHEME}.txt ../macro_dr/scripts/${EXPERIMENT}.txt  ../macro_dr/   scripts/simulation.txt "--runIdName= \"${PATH_MACRO_DR}${EXPER_ABR}_${CP}c_${N_SCOUTS}s_${N_BETA}b_${SCHEME}${SCH_ABR}_${LOCAL_ID}_0\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/beta_${N_BETA}.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})" ../macro_dr/scripts/evidence_thermo_continuation.txt 
fi





