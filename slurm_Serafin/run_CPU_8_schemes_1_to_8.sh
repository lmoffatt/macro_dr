#!/bin/bash

cd ~/Code/macro_dr/macro_dr


export WORKING_DIRECTORY=data_CPU8

RUNTIME="0-01:00"

export PATH_MACRO=/home/lmoffatt/macro_dr/

# Cargar los módulos para la tarea
module load cmake
module load gcc

module load amdblis
module load amdlibflame
module load gsl


export PARTITION=short

export EVIDENCE_ALGORITHM=thermo_dts

export NTASKS=8
CPUSPERTASK=8

export CP=CPUSPERTASK

export USE_LOCAL_ID=1

export N_SCH=1
N_SCH2=2
N_SCH3=3
N_SCH4=4
N_SCH5=5
N_SCH6=6
N_SCH7=7
N_SCH8=8
N_SCH9=9
N_SCH10=10
N_SCH11=11

export SCHEME_0=scheme_${N_SCH}_inact_PI
export SCHEME_1=scheme_${N_SCH2}_inact_PI
export SCHEME_2=scheme_${N_SCH3}_inact_PI
export SCHEME_3=scheme_${N_SCH4}_inact_PI
export SCHEME_4=scheme_${N_SCH5}_inact_PI
export SCHEME_5=scheme_${N_SCH6}_inact_PI
export SCHEME_6=scheme_${N_SCH7}_inact_PI
export SCHEME_7=scheme_${N_SCH8}_inact_PI
export SCHEME_8=scheme_${N_SCH9}_inact_PI
export SCHEME_9=scheme_${N_SCH10}_inact_PI
export SCHEME_10=scheme_${N_SCH11}_inact_PI
export SCHEME_11=scheme_${N_SCH11}_inact_PI
export SCHEME_12=scheme_${N_SCH10}_inact_PI
export SCHEME_13=scheme_${N_SCH9}_inact_PI
export SCHEME_14=scheme_${N_SCH8}_inact_PI
export SCHEME_15=scheme_${N_SCH7}_inact_PI


export IDNAME_0=w1_IE_8c_32s_4b_scheme_1_inact_PI_Ag_0_0
export IDNAME_1=w1_IE_8c_32s_4b_scheme_2_inact_PI_Ag_1_0
export IDNAME_2=w1_IE_8c_32s_4b_scheme_3_inact_PI_Ag_2_0
export IDNAME_3=w1_IE_8c_32s_4b_scheme_4_inact_PI_Ag_3_0
export IDNAME_4=w1_IE_8c_32s_4b_scheme_5_inact_PI_Ag_4_0
export IDNAME_5=w1_IE_8c_32s_4b_scheme_6_inact_PI_Ag_5_0
export IDNAME_6=w1_IE_8c_32s_4b_scheme_7_inact_PI_Ag_6_0
export IDNAME_7=w1_IE_8c_32s_4b_scheme_8_inact_PI_Ag_7_0
export IDNAME_8=w1_IE_8c_32s_4b_scheme_9_inact_PI_Ag_0_0
export IDNAME_9=w1_IE_8c_32s_4b_scheme_10_inact_PI_Ag_1_0

export IDNAME_10=w1_IE_8c_32s_4b_scheme_11_inact_PI_Ag_2_0
export IDNAME_11=w1_IE_8c_32s_4b_scheme_11_inact_PI_Ag_3_0
export IDNAME_12=w1_IE_8c_32s_4b_scheme_10_inact_PI_Ag_4_0

export IDNAME_13=w1_IE_8c_32s_4b_scheme_9_inact_PI_Ag_5_0
export IDNAME_14=w1_IE_8c_32s_4b_scheme_8_inact_PI_Ag_6_0
export IDNAME_15=w1_IE_8c_32s_4b_scheme_7_inact_PI_Ag_7_0


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

export PATH_MACRO_DR_0=w4_
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


export PATH_MACRO_DRX=w4


export LIK_0=DR
export LIK_1=$LIK_0
export LIK_2=$LIK_0
export LIK_3=$LIK_0
export LIK_4=DR
export LIK_5=$LIK_4
export LIK_6=$LIK_4
export LIK_7=$LIK_4
export LIK_8=DR
export LIK_9=$LIK_8
export LIK_10=$LIK_8
export LIK_11=$LIK_8
export LIK_12=DR
export LIK_13=$LIK_12
export LIK_14=$LIK_12
export LIK_15=$LIK_12





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


export N_BETA=4

export N_SCOUTS=32

export MAX_ITER=1000000



JOBID1=12707 


export CONTINUATION_NUMBER=0

#JOBID1=$(sbatch --parsable --job-name=R${N_SCH}_${CPUSPERTASK}  --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_Serafin/M_scheme_N_tasks.sh) 

export CONTINUATION_NUMBER=1
JOBID1=$(sbatch --parsable --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER}   --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_Serafin/M_scheme_N_tasks.sh) 



for i in $(seq 1 0);
do
    export CONTINUATION_NUMBER=$i
    JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=C${N_SCH}_${CPUSPERTASK}_${CONTINUATION_NUMBER}   --partition=${PARTITION} --ntasks-per-node=${NTASKS} --cpus-per-task=${CPUSPERTASK}  --time=${RUNTIME}  slurm_Serafin/M_scheme_N_tasks.sh) 
done



