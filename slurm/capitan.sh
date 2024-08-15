#!/bin/bash

cd ~/Code/macro_dr/macro_dr

module use /opt/ohpc/pub/apps/modules/all/
module load cmake
export CC=$HOME/local/gcc-14.2.0/bin/gcc
export CXX=$HOME/local/gcc-14.2.0/bin/g++
module load GSL
module load OpenBLAS

export PATH_MACRO=${HOME}/Code/macro_dr/


export PARTITION=eth_hi


