#!/bin/bash



SCHEME_DIRS=("models" "models_Ag")
SCHEME_FILES=( ${SCHEME_0}  )


EXPERIMENTS=(  "idealize_experiment" )

NUM_BETAS=( 32  )
NUM_SCOUTS=( 32 )

MAX_ITERS=( 1000000 ) 

IND=$(( SLURM_LOCALID % ${#SCHEME_DIRS[@]}))
SCHEME_DIR=${SCHEME_DIRS[$IND]}


IND=$(( SLURM_LOCALID % ${#SCHEME_FILES[@]}))

SCHEME=${SCHEME_FILES[$IND]}


IND=$(( SLURM_LOCALID % ${#EXPERIMENTS[@]}))
EXPERIMENT=${EXPERIMENTS[$IND]}

IND=$(( SLURM_LOCALID % ${#NUM_BETAS[@]}))
N_BETA=${NUM_BETAS[$IND]}

IND=$(( SLURM_LOCALID % ${#NUM_SCOUTS[@]}))
N_SCOUTS=${NUM_SCOUTS[$IND]}

IND=$(( SLURM_LOCALID % ${#MAX_ITERS[@]}))
MAX_ITER=${MAX_ITERS[$IND]}

EXPER_ABR=$([ "$EXPERIMENT" = "idealize_experiment" ] && echo "_IE" || echo "")

SCH_ABR=$([ "$SCHEME_DIR" = "models_Ag" ] && echo "_Ag" || echo "")

LOCAL_ID=$([ "$USE_LOCAL_ID" = 1 ] && echo ${SLURM_LOCALID} || echo "")

CP=$SLURM_CPUS_PER_TASK



if [ "$CONTINUATION_NUMBER" = 0 ]; then
    ${PATH_MACRO}/${PATH_MACRO_DRX}/macro_dr ../macro_dr/${SCHEME_DIR}/${SCHEME}.txt ../macro_dr/scripts/${EXPERIMENT}.txt  ../macro_dr/   scripts/simulation.txt "--runIdName= \"${PATH_MACRO_DR}${EXPER_ABR}_${CP}c_${N_SCOUTS}s_${N_BETA}b_${SCHEME}${SCH_ABR}_${LOCAL_ID}_0\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/beta_${N_BETA}.txt  ../macro_dr/scripts/evidence_thermo_data.txt 
else
${PATH_MACRO}/${PATH_MACRO_DRX}/macro_dr ../macro_dr/${SCHEME_DIR}/${SCHEME}.txt ../macro_dr/scripts/${EXPERIMENT}.txt  ../macro_dr/   scripts/simulation.txt "--runIdName= \"${PATH_MACRO_DR}${EXPER_ABR}_${CP}c_${N_SCOUTS}s_${N_BETA}b_${SCHEME}${SCH_ABR}_${LOCAL_ID}_0\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/beta_${N_BETA}.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})" ../macro_dr/scripts/evidence_thermo_continuation.txt 
fi





