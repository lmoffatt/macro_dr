#!/bin/bash

cd ~/macro_dr/macro_dr


export WORKING_DIRECTORY=data_w3

RUNTIME="0-06:00"

export PATH_MACRO=/home/lmoffatt/macro_dr

# Cargar los módulos para la tarea
module load cmake
module load gcc

module load amdblis
module load amdlibflame
module load gsl



export EVIDENCE_ALGORITHM=thermo_dts


export NTASKS=8
CPUSPERTASK=8
export CP=$CPUSPERTASK
export USE_LOCAL_ID=1

export N_SCH=10
N_SCH2=11

export SCHEME_DIR_0=models_Ag
export SCHEME_DIR_1=models_Ag
export SCHEME_DIR_2=models_Ag
export SCHEME_DIR_3=models_Ag
export SCHEME_DIR_4=models_Ag
export SCHEME_DIR_5=models_Ag
export SCHEME_DIR_6=models_Ag
export SCHEME_DIR_7=models_Ag
export SCHEME_DIR_8=models_Ag
export SCHEME_DIR_9=models_Ag
export SCHEME_DIR_10=models_Ag
export SCHEME_DIR_11=models_Ag
export SCHEME_DIR_12=models_Ag
export SCHEME_DIR_13=models_Ag
export SCHEME_DIR_14=models_Ag
export SCHEME_DIR_15=models_Ag


export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=scheme_${N_SCH2}_inact_PI
export SCHEME_2=SCHEME_0
export SCHEME_3=SCHEME_1
export SCHEME_4=$SCHEME_0
export SCHEME_5=$SCHEME_1
export SCHEME_6=$SCHEME_0
export SCHEME_7=$SCHEME_1
export SCHEME_8=$SCHEME_0
export SCHEME_9=$SCHEME_1
export SCHEME_10=$SCHEME_0
export SCHEME_11=$SCHEME_1
export SCHEME_12=$SCHEME_0
export SCHEME_13=$SCHEME_1
export SCHEME_14=$SCHEME_0
export SCHEME_15=$SCHEME_1


export PATH_MACRO_DR_0=w3_
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0
export PATH_MACRO_DR_2=$PATH_MACRO_DR_0
export PATH_MACRO_DR_3=$PATH_MACRO_DR_0
export PATH_MACRO_DR_4=$PATH_MACRO_DR_0
export PATH_MACRO_DR_5=$PATH_MACRO_DR_0
export PATH_MACRO_DR_6=$PATH_MACRO_DR_0
export PATH_MACRO_DR_7=$PATH_MACRO_DR_0
export PATH_MACRO_DR_8=$PATH_MACRO_DR_0
export PATH_MACRO_DR_9=$PATH_MACRO_DR_0
export PATH_MACRO_DR_10=$PATH_MACRO_DR_0
export PATH_MACRO_DR_11=$PATH_MACRO_DR_0
export PATH_MACRO_DR_12=$PATH_MACRO_DR_0
export PATH_MACRO_DR_13=$PATH_MACRO_DR_0
export PATH_MACRO_DR_14=$PATH_MACRO_DR_0
export PATH_MACRO_DR_14=$PATH_MACRO_DR_0


export PATH_MACRO_DRX=w3


export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0
export EXPERIMENT_2=$EXPERIMENT_0
export EXPERIMENT_3=$EXPERIMENT_0
export EXPERIMENT_4=$EXPERIMENT_0
export EXPERIMENT_5=$EXPERIMENT_0
export EXPERIMENT_6=$EXPERIMENT_0
export EXPERIMENT_7=$EXPERIMENT_0
export EXPERIMENT_8=$EXPERIMENT_0
export EXPERIMENT_9=$EXPERIMENT_0
export EXPERIMENT_10=$EXPERIMENT_0
export EXPERIMENT_11=$EXPERIMENT_0
export EXPERIMENT_12=$EXPERIMENT_0
export EXPERIMENT_13=$EXPERIMENT_0
export EXPERIMENT_14=$EXPERIMENT_0
export EXPERIMENT_15=$EXPERIMENT_0


export LIK_0=ADR
export LIK_1=$LIK_0
export LIK_2=DR
export LIK_3=$LIK_2
export LIK_4=R
export LIK_5=$LIK_4
export LIK_6=NR
export LIK_7=$LIK_6


export N_BETA=4

export N_SCOUTS=32

export MAX_ITER=1000000



JOBID1=223948


export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi  ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh) 
 


export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=8
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=9
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)

export CONTINUATION_NUMBER=10
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} --partition=multi ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh)


