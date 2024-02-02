#!/bin/bash

case $(($SLURM_TASK_ID % 8)) in
0)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model00_7.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

1)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model00.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

2)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model01.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

3)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model4.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

4)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model6_Eff_no_inactivation.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
5)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model7.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
6)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model8.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;
7)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/model9.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt ../macro_dr/scripts/evidence_data.txt
;;

esac

