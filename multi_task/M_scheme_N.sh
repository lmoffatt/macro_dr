#!/bin/bash

ECHO $SCHEME_FILES

IND=$(( SLURM_LOCALID % ${#SCHEME_FILES[@]}))

SCHEME=${SCHEME_FILES[$IND]}

IND=$(( SLURM_LOCALID % ${#PATH_MACRO_DRS[@]}))
PATH_MACRO_DR=${PATH_MACRO_DRS[$IND]}

IND=$(( SLURM_LOCALID % ${#EXPERIMENTS[@]}))
EXPERIMENT=${EXPERIMENTS[$IND]}

IND=$(( SLURM_LOCALID % ${#NUM_BETAS[@]}))
N_BETA=${NUM_BETAS[$IND]}

IND=$(( SLURM_LOCALID % ${#NUM_SCOUTS[@]}))
N_SCOUTS=${NUM_SCOUTS[$IND]}

IND=$(( SLURM_LOCALID % ${#MAX_ITERS[@]}))
MAX_ITER=${MAX_ITERS[$IND]}

EXPER_ABR=$([ "$EXPERIMENT" = "idealize_experiment" ] && echo "_IE" || echo "")



$PATH_MACRO_DR/macro_dr ../macro_dr/models/${SCHEME}.txt ../macro_dr/scripts/${EXPERIMENT}.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21${EXPER_ABR}_${N_SCOUTS}s_${N_BETA}b_${SCHEME}_0\""  "--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})" "--max_iter_equilibrium = get_number(n=${MAX_ITER})" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/beta_${N_BETA}.txt  ../macro_dr/scripts/evidence_thermo_data.txt



