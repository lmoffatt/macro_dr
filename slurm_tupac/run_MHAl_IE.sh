#!/bin/bash

cd ~/Code/macro_dr/macro_dr

JOBID1=$(sbatch --parsable slurm_tupac/MHAl_IE.sh )

export CONTINUATION_NUMBER=1

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=4
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=5
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=6
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)

export CONTINUATION_NUMBER=7
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHAl_IE_continuation.sh)
