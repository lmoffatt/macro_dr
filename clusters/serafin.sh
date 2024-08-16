#!/bin/bash

cd ~/macro_dr/macro_dr

module load cmake
module load gcc
module load amdblis
module load amdlibflame
module load gsl
export PATH_MACRO=${HOME}/macro_dr
export CLUSTER=serafin

if [ -z "${PARTITION}" ]; then                                           
   export PARTITION=multi
fi  
