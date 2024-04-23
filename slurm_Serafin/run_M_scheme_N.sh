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


NTASKS=16
CPUSPERTASK=4


export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable --job-name=R1_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C1_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C1_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C1_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C1_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C1_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_tupac/M_scheme_N.sh) 


#NTASKS=8
#CPUSPERTASK=8

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 

#NTASKS=4
#CPUSPERTASK=16

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 


#NTASKS=2
#CPUSPERTASK=32

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 
