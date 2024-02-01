#!/bin/bash

case $SLURM_LOCALID in
0)
/home/lmoffatt/macro_dr/macro_dr/macro_dr  ../macro_dr/models/scheme_1.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v2\"" ../macro_dr/scripts/evidence_data.txt
;;

1)
/home/lmoffatt/macro_dr/macro_dr/macro_dr  ../macro_dr/models/scheme_2.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v2\"" ../macro_dr/scripts/evidence_data.txt
;;

2)
/home/lmoffatt/macro_dr/macro_dr/macro_dr  ../macro_dr/models/scheme_3.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v2\"" ../macro_dr/scripts/evidence_data.txt
;;

3)
/home/lmoffatt/macro_dr/macro_dr/macro_dr  ../macro_dr/models/scheme_4.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v2\"" ../macro_dr/scripts/evidence_data.txt
;;

esac

