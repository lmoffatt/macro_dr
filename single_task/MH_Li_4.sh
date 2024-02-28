#!/bin/bash


$PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact.txt ../macro_dr/scripts/experiment.txt  ../macro_dr/scripts/simulation.txt "runIdNamePre= \"MH_Li_inact_9a\""  "num_scouts_per_ensemble = get_number(n=32)"  ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt

