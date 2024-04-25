#!/bin/bash

cd ~/Code/macro_dr/macro_dr



RUNTIME="2-00:00"

export PATH_MACRO=/nfs/home/lmoffatt/Code/macro_dr/



NTASKS=2
CPUSPERTASK=32

N_SCH=7

export SCHEME_DIR_0=models
export SCHEME_DIR_1=models_Ag

export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=$SCHEME_0


export PATH_MACRO_DR_0=v23
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0


export EXPERIMENT_0=idealize_experiment
export EXPERIMENT_1=$EXPERIMENT_0


#export CONTINUATION_NUMBER=0

#JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

JOBID1=8704

export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 


#NTASKS=8
#CPUSPERTASK=8

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 

#NTASKS=4
#CPUSPERTASK=16

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 


#NTASKS=2
#CPUSPERTASK=32

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 
