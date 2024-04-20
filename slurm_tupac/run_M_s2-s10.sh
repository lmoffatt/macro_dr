#!/bin/bash

cd ~/Code/macro_dr/macro_dr

export NUM_SCOUTS=32
export NUM_BETA=32

JOBID1=$(sbatch --parsable --job-name=M0_33 slurm_tupac/M_s2-s10.sh )

export CONTINUATION_NUMBER=1

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M1_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M2_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M3_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M4_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M5_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M6_33 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M7_33 slurm_tupac/M_s2-s10_continuation.sh)

export NUM_SCOUTS=32
export NUM_BETA=64

JOBID1=$(sbatch --parsable --job-name=M0_36 slurm_tupac/M_s2-s10.sh )

export CONTINUATION_NUMBER=1

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M1_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M2_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M3_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M4_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M5_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M6_36 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M7_36 slurm_tupac/M_s2-s10_continuation.sh)


export NUM_SCOUTS=64
export NUM_BETA=64

JOBID1=$(sbatch --parsable --job-name=M0_66 slurm_tupac/M_s2-s10.sh )

export CONTINUATION_NUMBER=1

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M1_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M2_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M3_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M4_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M5_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M6_66 slurm_tupac/M_s2-s10_continuation.sh)

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 --job-name=M7_66 slurm_tupac/M_s2-s10_continuation.sh)

