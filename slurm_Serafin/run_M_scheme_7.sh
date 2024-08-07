#!/bin/bash

export PATH_MACRO=/home/lmoffatt/macro_dr/


# Cargar los módulos para la tarea
module load cmake
module load gcc

module load amdblis
module load amdlibflame
module load gsl

RUNTIME="2-00:00"

export NTASKS=2
CPUSPERTASK=32
export CP=CPUSPERTASK



export USE_LOCAL_ID=0

export SCHEME_DIR_0=models
export SCHEME_DIR_1=models_Ag



export PATH_MACRO_DR_0=v23
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0

export PATH_MACRO_DRX=v32


export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0

export N_BETA=32

export N_SCOUTS=32

export MAX_ITER=1000000

export EVIDENCE_ALGORITHM=thermo

export CONTINUATION_NUMBER=0

export N_SCH=7
export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=$SCHEME_0

JOBID1=8710


#JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N_tasks.sh) 

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N_tasks.sh) 

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=8
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=9
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 

export CONTINUATION_NUMBER=10
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi slurm_Serafin/M_scheme_N.sh) 


