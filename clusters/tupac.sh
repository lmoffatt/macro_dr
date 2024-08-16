#!/bin/bash
cd ~/Code/macro_dr/macro_dr

export PATH_MACRO=${HOME}/Code/macro_dr/
export CLUSTER=tupac
if [ -n "${PARTITION}" ]; then                                           
   export PARTITION=free-rider
fi  

