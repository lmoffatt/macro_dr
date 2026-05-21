#!/bin/bash

cd ~/Code/macro_dr/macro_dr

module use /opt/ohpc/pub/apps/modules/all/
module load GSL/2.7-GCC-12.2.0
module load intel/2022.0.2

export PATH_MACRO=${HOME}/Code/macro_dr
export CLUSTER=capitan

if [ -z "${PARTITION}" ]; then                                           
   export PARTITION=eth_hi
fi  
