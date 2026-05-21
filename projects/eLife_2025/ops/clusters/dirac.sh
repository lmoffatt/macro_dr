#!/bin/bash
# Cluster profile: dirac.df.uba.ar (UBA Departamento de Física, OpenHPC stack).
#
# Hardware: ~30 Intel Xeon Dell R720xd nodes (~18 cores avg) + 3 AMD EPYC nodes;
# total 600 cores / 1200 threads / 2.5 TB RAM across ~33 nodes. Mostly Intel,
# so MKL is the BLAS of choice — compiler-agnostic and avoids the gnu7-only
# openblas modules on this cluster.
#
# Source from your interactive shell (`source projects/p2x2/ops/clusters/dirac.sh`)
# before building or submitting jobs.

cd ~/Code/macro_dr/macro_dr

module purge
module load gnu14            # gcc 14.3 — full C++20
module load mkl/2019.5.281   # BLAS + LAPACK via Intel MKL
module load gsl              # default gsl/2.7
module load ninja            # build-system used by CMake presets
module load cmake/3.30.8     # presets need ≥3.18

export PATH_MACRO=${HOME}/Code/macro_dr
export CLUSTER=dirac

# CMake's FindBLAS dispatch — sequential MKL (we use OMP at our level, not MKL's)
export BLA_VENDOR=Intel10_64lp_seq

# Scratch root for SLURM outputs (home quota is 160 GB; scratch is 4 TB BeeGFS,
# files >30 days auto-deleted). Submit lines should point WORKDIR here.
export SCRATCH_MACRO=/scratch/${USER}/macro_dr

if [ -z "${PARTITION}" ]; then
    export PARTITION=batch
fi
