#!/bin/bash

cd ~/macro_dr/macro_dr

module load cmake
module load gcc
module load amdblis
module load amdlibflame
module load gsl
export PATH_MACRO=${HOME}/macro_dr
export CLUSTER=serafin
export RUNTIME=2-00:00
export EVIDENCE_ALGORITHM=thermo_dts
export MAX_ITER=1000000

if [ -z "${PARTITION}" ]; then                                           
   export PARTITION=multi
fi  
