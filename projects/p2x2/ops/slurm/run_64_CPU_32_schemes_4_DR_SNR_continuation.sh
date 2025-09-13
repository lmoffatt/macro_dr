#!/bin/bash

export NTASKS=2
CPUSPERTASK=32

export CP=$CPUSPERTASK

export USE_LOCAL_ID=0

export IDNAME_0=w9_IE_DR_8c_32s_4b_scheme_4_inact_PI_logbaseline__0
export IDNAME_1=w9_IE_SNR_8c_32s_4b_scheme_4_inact_PI_logbaseline__0

export N_SCH=4
N_SCH2=4


MODEL_TYPE0=inact_PI
MODEL_TYPE1=$MODEL_TYPE0

ALG_TYPE0=DR
ALG_TYPE1=SNR


export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=scheme_${N_SCH2}_inact_PI

export LIK_0=${ALG_TYPE0}
export LIK_1=${ALG_TYPE1}


export SCHEME_DIR_0=models_Ag_log_baseline
export SCHEME_DIR_1=${SCHEME_DIR_0}

export PATH_MACRO_DR_0=w9
export PATH_MACRO_DR_1=$PATH_MACRO_DR_0

export PATH_MACRO_DRX=w9






export EXPERIMENT_0=idealize_experiment_2
export EXPERIMENT_1=$EXPERIMENT_0

export N_BETA=4

export N_SCOUTS=32

export MAX_ITER=1000000




export CONTINUATION_NUMBER=5

JOBID_4DR=21343
JOBID_4SNR=21343


JOBID1=$(sbatch --parsable  --dependency=afterany:${JOBID_4DR}:${JOBID_4SNR} --job-name=R${N_SCH}_${CPUSPERTASK}  --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh) 


for i in $(seq 6 15);
do
    export CONTINUATION_NUMBER=$i
    JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER}   --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  ${PATH_MACRO}/macro_dr/slurm/M_scheme_N_tasks.sh) 
done



