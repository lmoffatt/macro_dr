#!/bin/bash

cd ~/Code/macro_dr/macro_dr

JOBID1=$(sbatch --parsable slurm_tupac/MHLi.sh )

sbatch --dependency=afterany:$JOBID1 slurm_tupac/MHLi_continuation.sh
