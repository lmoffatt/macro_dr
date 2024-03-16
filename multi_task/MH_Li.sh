#!/bin/bash


case $SLURM_LOCALID in
0)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_1_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme1_NP_0\""  "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

1)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_2_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme2_NP_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt  ../macro_dr/scripts/evidence_thermo_data.txt
;;

2)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme3_NP_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

3)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme4_NP_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

4)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_1_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme1_NP_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

5)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_2_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme2_NP_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt  ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

6)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact_NP.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme3_NP_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

7)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme4_NP_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

8)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_1_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme1_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

9)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_2_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme2_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt  ../macro_dr/scripts/evidence_thermo_data.txt
;;

10)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme3_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

11)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme4_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
;;

12)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_1_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme1_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

13)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_2_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme2_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt  ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

14)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_3_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme3_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

15)
$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "--runIdName= \"v19_scheme4_sim_0\"" "--num_scouts_per_ensemble = get_number(n=32)" ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_simulation.txt
;;

esac

