#!/bin/bash

#case $(($SLURM_TASK_ID % 4))
case $(($$ % 4)) in
0)
$PATH_MACRO_DR/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

1)
$PATH_MACRO_DR/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

2)
$PATH_MACRO_DR/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

3)
$PATH_MACRO_DR/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

esac

