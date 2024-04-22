#!/bin/bash

cd ~/Code/macro_dr/macro_dr



SCHEME_FILES=( scheme_7_inact_PI  )

PATH_MACRO_DRS=("/nfs/home/lmoffatt/Code/macro_dr/v22/" "/nfs/home/lmoffatt/Code/macro_dr/v22/" )

EXPERIMENTS=( "idealize_experiment" )

N_BETA=( 32 )
N_SCOUTS=( 32 )

MAX_ITERS=( 20 ) 


NTASKS=16
CPUSPERTASK=4

JOBID1=$(sbatch --parsable --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} slurm_tupac/M_scheme_N.sh )


NTASKS=8
CPUSPERTASK=8

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} slurm_tupac/M_scheme_N.sh ) 

NTASKS=4
CPUSPERTASK=16

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} slurm_tupac/M_scheme_N.sh ) 

NTASKS=2
CPUSPERTASK=32

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} slurm_tupac/M_scheme_N.sh ) 

NTASKS=1
CPUSPERTASK=64

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} slurm_tupac/M_scheme_N.sh ) 
