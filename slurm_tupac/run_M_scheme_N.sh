#!/bin/bash

cd ~/Code/macro_dr/macro_dr



RUNTIME="0-01:00"

export PATH_MACRO=/nfs/home/lmoffatt/Code/macro_dr/



NTASKS=2
CPUSPERTASK=32


export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable --job-name=R_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --job-name=R_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --job-name=R_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 


#NTASKS=8
#CPUSPERTASK=8

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 

#NTASKS=4
#CPUSPERTASK=16

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 


#NTASKS=2
#CPUSPERTASK=32

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 
