#!/bin/bash

cd ~/Code/macro_dr/macro_dr

module use /opt/ohpc/pub/apps/modules/all/
module load cmake
module load intel/2022.0.2
module load GSL/2.7-GCC-12.2.0module load GSL

export PATH_MACRO=${HOME}/Code/macro_dr/


export PARTITION=eth_hi


