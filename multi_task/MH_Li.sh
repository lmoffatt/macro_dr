#!/bin/bash


case $SLURM_LOCALID in
0)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_1_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li_inact_v7\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

1)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_2_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li_inact_v7\"" ../macro_dr/scripts/likelihood.txt  ../macro_dr/scripts/evidence_thermo_data.txt
;;

2)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li_inact_v7\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

3)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li_inact_v7\"" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

esac

