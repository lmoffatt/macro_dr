#!/bin/bash

cd ~/Code/macro_dr/macro_dr



RUNTIME="2-00:00"

export PATH_MACRO=/nfs/home/lmoffatt/Code/macro_dr/



export NTASKS=2
CPUSPERTASK=32


export SCHEME_0=scheme_11_inact_PI
export EVIDENCE_ALGORITHM=levenberg


export PATH_MACRO_DR=v31
export PATH_MACRO_DRX=v31

SCM_N=7
export CONTINUATION_NUMBER=0

export N_BETA=16

export N_SCOUTS=32

export MAX_ITER=1000000

export SCHEME_DIR_0=models
export SCHEME_DIR_2=models_Ag


export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=$SCHEME_0

export PATH_MACRO_DR_0=v31
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0

export PATH_MACRO_DRX=v31

JOBID1=$(sbatch --parsable --job-name=R${SCM_N}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

#export CONTINUATION_NUMBER=1
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

#export CONTINUATION_NUMBER=2
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

#export CONTINUATION_NUMBER=3
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

#export CONTINUATION_NUMBER=4
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

#export CONTINUATION_NUMBER=5
#JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 



