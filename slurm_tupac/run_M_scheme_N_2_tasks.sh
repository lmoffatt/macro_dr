#!/bin/bash

cd ~/Code/macro_dr/macro_dr



RUNTIME="2-00:00"

export PATH_MACRO=/nfs/home/lmoffatt/Code/macro_dr/



export NTASKS=2
CPUSPERTASK=32

N_SCH=11


export USE_LOCAL_ID=1

export SCHEME_DIR_0=models
export SCHEME_DIR_1=models_Ag

export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=$SCHEME_0


export PATH_MACRO_DR_0=v28
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0

export PATH_MACRO_DRX=v28


export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0

export N_BETA=32

export N_SCOUTS=32

export MAX_ITER=1000000

export EVIDENCE_ALGORITHM=levenberg

export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

#JOBID1=9015

#export CONTINUATION_NUMBER=1
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

#export CONTINUATION_NUMBER=2
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

#export CONTINUATION_NUMBER=3
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

#export CONTINUATION_NUMBER=4
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

#export CONTINUATION_NUMBER=5
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 


#NTASKS=8
#CPUSPERTASK=8

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 

#NTASKS=4
#CPUSPERTASK=16

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 


#NTASKS=2
#CPUSPERTASK=32

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 
