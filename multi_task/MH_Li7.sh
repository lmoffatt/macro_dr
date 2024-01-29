#!/bin/bash

case $(($SLURM_TASK_ID % 4)) in
#case $(($$ % 4)) in
0)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/models/scheme_1.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"EviData_7\"" ../macro_dr/scripts/evidence_data.txt
;;

1)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/models/scheme_2.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"EviData_7\"" ../macro_dr/scripts/evidence_data.txt
;;

2)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/models/scheme_3.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"EviData_7\"" ../macro_dr/scripts/evidence_data.txt
;;

3)
/home/lmoffatt/Code/macro_dr/build-macro_dr-gcc-Release/macro_dr ../macro_dr/models/scheme_4.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"EviData_7\"" ../macro_dr/scripts/evidence_data.txt
;;

esac

