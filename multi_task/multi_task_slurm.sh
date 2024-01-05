#!/bin/bash

#case $(($SLURM_TASK_ID % 4))
case $(($$ % 4)) in
0)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

1)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

2)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

3)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/scripts/model.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_simulation.txt
;;

esac

