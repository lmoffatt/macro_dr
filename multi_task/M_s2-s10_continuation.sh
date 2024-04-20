#!/bin/bash


case $SLURM_LOCALID in
0)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_3_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

1)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_4_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

2)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_5_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_5_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

3)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_6_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_6_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

4)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_7_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_7_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

5)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_8_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_8_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

6)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_9_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_9_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

7)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_10_inact_PI.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_10_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

8)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_3_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

9)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_4_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

10)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_5_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_5_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

11)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_6_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_6_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

12)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_7_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_7_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

13)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_8_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_8_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

14)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_9_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_9_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;

15)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_10_inact_PI.txt ../macro_dr/scripts/idealize_experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v21_IE_${NUM_SCOUTS}s_${NUM_BETA}b_scheme_10_PI_0\""  "--num_scouts_per_ensemble = get_number(n=${NUM_SCOUTS})" ../macro_dr/scripts/likelihood.txt "--continuation_number=get_number(n=${CONTINUATION_NUMBER})"  ../macro_dr/scripts/beta_${NUM_BETA}.txt  ../macro_dr/scripts/evidence_thermo_continuation.txt
;;


esac

