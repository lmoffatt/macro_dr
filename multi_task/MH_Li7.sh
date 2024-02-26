#!/bin/bash

case $SLURM_LOCALID in
0)
$PATH_MACRO_DR/macro_dr  ../macro_dr/models/scheme_1.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v8a\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_fraction_data.txt
;;

1)
$PATH_MACRO_DR/macro_dr  ../macro_dr/models/scheme_2.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v8a\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_fraction_data.txt
;;

2)
$PATH_MACRO_DR/macro_dr  ../macro_dr/models/scheme_3.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v8a\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_fraction_data.txt
;;

3)
$PATH_MACRO_DR/macro_dr  ../macro_dr/models/scheme_4.txt ../macro_dr/scripts/experiment_7.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li7_v8a\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_fraction_data.txt
;;

esac

