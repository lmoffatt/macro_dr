#!/bin/bash

cd ~/Code/macro_dr/macro_dr

JOBID1=$(sbatch --parsable slurm_tupac/MHLi_IE.sh )


sbatch --dependency=afterany:$JOBID1 slurm_tupac/MHLi_IE_continuation.sh
