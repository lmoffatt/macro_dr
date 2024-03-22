#!/bin/bash

cd ~/Code/macro_dr/macro_dr

JOBID1=$(sbatch --parsable slurm_tupac/MHLi.sh )

export CONTINUATION_NUMBER=1

JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHLi_continuation.sh)

export CONTINUATION_NUMBER=2
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHLi_continuation.sh)

export CONTINUATION_NUMBER=3
JOBID1=$(sbatch --parsable --dependency=afterany:$JOBID1 slurm_tupac/MHLi_continuation.sh)


