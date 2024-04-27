#!/bin/bash

cd ~/macro_dr/macro_dr



RUNTIME="2-00:00"

export PATH_MACRO=/home/lmoffatt/macro_dr/

# Cargar los módulos para la tarea
module load cmake
module load gcc

module load amdblis
module load amdlibflame
module load gsl


NTASKS=2
CPUSPERTASK=32
export USE_LOCAL_ID=1

export CONTINUATION_NUMBER=0


export SCHEME_0=scheme_14_inact_PI

export PATH_MACRO_DR=v26
export PATH_MACRO_DRX=v27

SCM_N=14

#JOBID1=$(sbatch --parsable --job-name=R${SCM_N}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi  slurm_Serafin/M_scheme_N.sh) 

JOBID1=193348

export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi  slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${SCM_N}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

