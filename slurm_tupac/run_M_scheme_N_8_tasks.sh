#!/bin/bash

cd ~/Code/macro_dr/macro_dr



RUNTIME="2-00:00"

export PATH_MACRO=/nfs/home/lmoffatt/Code/macro_dr/

export EVIDENCE_ALGORITHM=thermo


export NTASKS=8
CPUSPERTASK=8

N_SCH=3
N_SCH2=4
N_SCH3=4
N_SCH4=5

export USE_LOCAL_ID=1

export SCHEME_DIR_0=models
export SCHEME_DIR_1=$SCHEME_DIR_0
export SCHEME_DIR_2=$SCHEME_DIR_0
export SCHEME_DIR_3=$SCHEME_DIR_0
export SCHEME_DIR_4=models_Ag
export SCHEME_DIR_5=$SCHEME_DIR_4
export SCHEME_DIR_6=$SCHEME_DIR_4
export SCHEME_DIR_7=$SCHEME_DIR_4



export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=scheme_${N_SCH2}_inact_PI
export SCHEME_2=scheme_${N_SCH3}_inact_PI
export SCHEME_3=scheme_${N_SCH4}_inact_PI
export SCHEME_4=scheme_${N_SCH}_inact_PI
export SCHEME_5=scheme_${N_SCH2}_inact_PI
export SCHEME_6=scheme_${N_SCH3}_inact_PI
export SCHEME_7=scheme_${N_SCH4}_inact_PI


export PATH_MACRO_DR_0=v31
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0
export PATH_MACRO_DR_2=$PATH_MACRO_DR_0
export PATH_MACRO_DR_3=$PATH_MACRO_DR_0
export PATH_MACRO_DR_4=$PATH_MACRO_DR_0
export PATH_MACRO_DR_5=$PATH_MACRO_DR_0
export PATH_MACRO_DR_6=$PATH_MACRO_DR_0
export PATH_MACRO_DR_7=$PATH_MACRO_DR_0

export PATH_MACRO_DR=v31

export PATH_MACRO_DRX=v32


export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0
export EXPERIMENT_2=$EXPERIMENT_0
export EXPERIMENT_3=$EXPERIMENT_0
export EXPERIMENT_4=$EXPERIMENT_0
export EXPERIMENT_5=$EXPERIMENT_0
export EXPERIMENT_6=$EXPERIMENT_0
export EXPERIMENT_7=$EXPERIMENT_0


export N_BETA=32

export N_SCOUTS=32

export MAX_ITER=1000000



#JOBID1=8710


export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable  --job-name=R${N_SCH}_${CPUSPERTASK}  --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_tupac/M_scheme_N.sh) 

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

##NTASKS=8
##CPUSPERTASK=8

##sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK} --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 

##NTASKS=4
##CPUSPERTASK=16

##sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 


##NTASKS=2
##CPUSPERTASK=32

##sbatch --parsable  --job-name=T_${CPUSPERTASK} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME} slurm_tupac/M_scheme_N.sh 
