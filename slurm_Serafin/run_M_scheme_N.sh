#!/bin/bash

cd ~/macro_dr/macro_dr



RUNTIME="0-01:00"

export PATH_MACRO=/home/lmoffatt/macro_dr/

# Cargar los módulos para la tarea
module load amdblis
module load amdlibflame
module load gsl

NTASKS=16
CPUSPERTASK=4

sbatch --parsable --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_Serafin/M_scheme_N.sh 


NTASKS=8
CPUSPERTASK=8

sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_Serafin/M_scheme_N.sh 

NTASKS=4
CPUSPERTASK=16

sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_Serafin/M_scheme_N.sh 

NTASKS=2
CPUSPERTASK=32

sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_Serafin/M_scheme_N.sh 

#NTASKS=1
#CPUSPERTASK=64

#sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_Serafin/M_scheme_N.sh
