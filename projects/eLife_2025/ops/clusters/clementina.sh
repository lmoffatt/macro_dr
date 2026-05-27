#!/bin/bash
# Cluster profile: Clementina XXI (ssh.clementinaxxi.org.ar), Servicio
# Meteorológico Nacional. OpenHPC + spack + Lmod, SLURM. Reachable only through
# the SMN VPN — connect first with:
#   sudo openfortivpn vpn.clementinaxxi.org.ar:10443 -u <user>
# then `ssh <user>@ssh.clementinaxxi.org.ar`.
#
# Hardware: 80 nodes cn[001-080], each 2x Intel Xeon CPU Max 9462 = 64 cores +
# 512 GB RAM + 64 GB HBM. cn[001-074] also carry 4x Intel GPU Max 1550 (Ponte
# Vecchio) — unused: macro_dr is CPU code, so we target the cpunode partition.
# The 512 GB nodes mean the figure_2 OOM that bound dirac (64-128 GB nodes) does
# not bind here, but keep the serial-combo strategy for core efficiency anyway.
#
# Source from your interactive shell before building or submitting:
#   source projects/eLife_2025/ops/clusters/clementina.sh

cd ~/Code/macro_dr/macro_dr

module purge
module load gcc/15.1.0       # newer than dirac's gnu14 — full C++23, known-good family
module load mkl/2025.3       # oneMKL: BLAS + LAPACK (sets MKLROOT)
module load gsl/2.8          # /data/shared/modules — compiler-agnostic C lib
module load cmake/3.24.2     # OpenHPC core module
# NOTE: no ninja module on this cluster → MACRODR_GENERATOR forces Unix Makefiles.

export PATH_MACRO=${HOME}/Code/macro_dr
export CLUSTER=clementina

# CMake FindBLAS: sequential MKL. We parallelize with OpenMP at the per-simulation
# level, so MKL must NOT spawn its own BLAS threads underneath (avoids nesting).
export BLA_VENDOR=Intel10_64lp_seq

# No ninja here → build_cluster.sh reads this and configures Unix Makefiles.
export MACRODR_GENERATOR="Unix Makefiles"

# Every SLURM job on Clementina must charge a project allocation. Project PCI 86
# → pci_86. build_cluster.sh and dispatch_figure_2.sh add --account=$ACCOUNT
# when this is set (dirac leaves it unset → no flag).
export ACCOUNT=pci_86

# Output goes under the shared IPAC dir on /data: /home has a 500 GB cap and the
# docs advise against running from it. Namespaced by user to avoid colliding
# with other pci_86 members. (No documented auto-purge, unlike dirac's scratch.)
export SCRATCH_MACRO=/data/contrib/pci_86/$(whoami)/macro_dr

# CPU partition (the GPU partition's Ponte Vecchio cards are unused). The build
# (~30 min) exceeds the 20-min `testing` cap, so build on cpunode; use
# PARTITION=testing only for <=20 min smoke runs of the binary.
if [ -z "${PARTITION}" ]; then
    export PARTITION=cpunode
fi
