#!/bin/bash

case $(($SLURM_TASK_ID % 4)) in
0)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model00_7.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

1)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model00.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

2)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model01.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

3)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model4.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

4)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model6_Eff_no_inactivation.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
5)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model7.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
6)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model8.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
7)
/home/lmoffatt/macro_dr/macro_dr/macro_dr ../macro_dr/models/model9.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

esac

