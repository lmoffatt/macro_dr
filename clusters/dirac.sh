#!/bin/bash

cd ~/Code/macro_dr/macro_dr


module load gsl
module load gnu14
#module load openblas/0.3.28
export PATH_MACRO=${HOME}/Code/macro_dr
export CLUSTER=dirac


if [ -z "${PARTITION}" ]; then                                           
   export PARTITION=batch
fi  


