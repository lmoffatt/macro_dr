#!/bin/bash

export NTASKS=8
CPUSPERTASK=8

export CP=$CPUSPERTASK

export USE_LOCAL_ID=0

export N_SCH=1
N_SCH2=2


MODEL_TYPE0=inact_PI
MODEL_TYPE1=PI

ALG_TYPE0=DR
ALG_TYPE1=SNR



export SCHEME_0=scheme_${N_SCH}_${MODEL_TYPE0}
export SCHEME_1=scheme_${N_SCH}_${MODEL_TYPE1}
export SCHEME_2=scheme_${N_SCH}_${MODEL_TYPE0}
export SCHEME_3=scheme_${N_SCH}_${MODEL_TYPE1}
export SCHEME_4=scheme_${N_SCH2}_${MODEL_TYPE0}
export SCHEME_5=scheme_${N_SCH2}_${MODEL_TYPE1}
export SCHEME_6=scheme_${N_SCH2}_${MODEL_TYPE0}
export SCHEME_7=scheme_${N_SCH2}_${MODEL_TYPE1}



export SCHEME_DIR_0=models_Ag_log_baseline
export SCHEME_DIR_1=${SCHEME_DIR_0}
export SCHEME_DIR_2=${SCHEME_DIR_0}
export SCHEME_DIR_3=${SCHEME_DIR_0}
export SCHEME_DIR_4=${SCHEME_DIR_0}
export SCHEME_DIR_5=${SCHEME_DIR_0}
export SCHEME_DIR_6=${SCHEME_DIR_0}
export SCHEME_DIR_7=${SCHEME_DIR_0}

export PATH_MACRO_DR_0=w9
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0
export PATH_MACRO_DR_2=$PATH_MACRO_DR_0
export PATH_MACRO_DR_3=$PATH_MACRO_DR_0
export PATH_MACRO_DR_4=$PATH_MACRO_DR_0
export PATH_MACRO_DR_5=$PATH_MACRO_DR_0
export PATH_MACRO_DR_6=$PATH_MACRO_DR_0
export PATH_MACRO_DR_7=$PATH_MACRO_DR_0

export PATH_MACRO_DRX=w9




export LIK_0=${ALG_TYPE0}
export LIK_1=${ALG_TYPE0}
export LIK_2=${ALG_TYPE1}
export LIK_3=${ALG_TYPE1}
export LIK_4=${ALG_TYPE0}
export LIK_5=${ALG_TYPE0}
export LIK_6=${ALG_TYPE1}
export LIK_7=${ALG_TYPE1}





export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0
export EXPERIMENT_2=$EXPERIMENT_0
export EXPERIMENT_3=$EXPERIMENT_0
export EXPERIMENT_4=$EXPERIMENT_0
export EXPERIMENT_5=$EXPERIMENT_0
export EXPERIMENT_6=$EXPERIMENT_0
export EXPERIMENT_7=$EXPERIMENT_0


export N_BETA=4

export N_SCOUTS=32

export MAX_ITER=1000000




export CONTINUATION_NUMBER=0

JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK}  --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh) 


for i in $(seq 1 8);
do
    export CONTINUATION_NUMBER=$i
    JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER}   --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh) 
done



