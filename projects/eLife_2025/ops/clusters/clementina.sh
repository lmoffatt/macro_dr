#!/bin/bash
# Cluster profile: Clementina XXI (ssh.clementinaxxi.org.ar), Servicio
# Meteorológico Nacional. EasyBuild/Lmod software stages + SLURM. Reachable only
# through the SMN VPN — connect first with:
#   sudo openfortivpn vpn.clementinaxxi.org.ar:10443 -u <user>
# then `ssh <user>@ssh.clementinaxxi.org.ar`. Login node has no direct internet;
# git/pip go through the proxy 172.28.3.3:3128 (set http(s)_proxy in ~/.bashrc).
#
# Hardware: 80 nodes cn[001-080], each 2x Intel Xeon CPU Max 9462 = 64 cores +
# 512 GB RAM + 64 GB HBM. cn[001-074] also carry 4x Intel GPU Max 1550 (unused:
# macro_dr is CPU code → cpunode partition). The 512 GB nodes dissolve the
# figure_2 OOM that bound dirac, but keep the serial-combo strategy for core use.
#
# Toolchain note: the modern compiler lives behind the `stages/2026` EasyBuild
# stage and must be loaded first. MKL is chained to the Intel compiler module
# (one compiler at a time), so for the gcc15 build we use the stage's OpenBLAS
# (GCC-compatible, conflict-free) instead of MKL.
#
# Source from your interactive shell before building or submitting:
#   source projects/eLife_2025/ops/clusters/clementina.sh

cd ~/Code/macro_dr/macro_dr

module purge
module load stages/2026          # EasyBuild stack root — exposes the gcc15 toolchain
module load gcc/15.1.0           # newer than dirac's gnu14 — full C++23
module load openblas/0.3.29      # BLAS + LAPACK, GCC-built (MKL is gated behind icpx here)
module load cmake/3.31.11        # stage default; > dirac's 3.30.8
module load ninja/1.13.2         # present in the stage → keep Ninja generator (no Makefiles override)

export PATH_MACRO=${HOME}/Code/macro_dr
export CLUSTER=clementina

# GSL: the cluster's gsl/2.8 module is Intel-built (its modulefile loads
# intel/2024.0.0 and libgsl.so needs Intel runtime libs — __svml_*/libimf — so
# it won't link with gcc15), and there is no gcc15-native GSL in the stack. So
# we use a from-source GSL built with gcc15 into /data/contrib. Point the
# compiler, linker and CMake at it directly (no module). NOTE: GSL only backs an
# unused proportional-noise path in macro_dr, so it just needs to LINK; if the
# Bessel milestone ever calls GSL's special functions, rebuild it under
# -std=gnu17 (gcc15 defaults to C23, which miscompiles GSL 2.8 — values read 0).
export GSL_ROOT=/data/contrib/pci_86/opt/gsl-2.8
export GSL_ROOT_DIR="$GSL_ROOT"                            # CMake FindGSL
export CPATH="$GSL_ROOT/include:$CPATH"
export LIBRARY_PATH="$GSL_ROOT/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$GSL_ROOT/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$GSL_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH"

# CMake FindBLAS/FindLAPACK: OpenBLAS. Runtime keeps BLAS single-threaded
# (run_macroir.sh sets OPENBLAS_NUM_THREADS=1) because we parallelize at the
# per-simulation OpenMP level — a threaded BLAS underneath would oversubscribe.
export BLA_VENDOR=OpenBLAS

# Every SLURM job on Clementina must charge a project allocation. Project PCI 86
# → pci_86. build_cluster.sh and dispatch_figure_2.sh add --account=$ACCOUNT
# when set (dirac leaves it unset → no flag).
export ACCOUNT=pci_86

# Output goes under the shared IPAC dir on /data: /home has a 500 GB cap and the
# docs advise against running from it. Namespaced by user to avoid colliding
# with other pci_86 members.
export SCRATCH_MACRO=/data/contrib/pci_86/$(whoami)/macro_dr

# CPU partition (GPU partition's Ponte Vecchio cards unused). The build (~30 min)
# exceeds the 20-min `testing` cap, so build on cpunode; use PARTITION=testing
# only for <=20 min smoke runs of the binary.
if [ -z "${PARTITION}" ]; then
    export PARTITION=cpunode
fi
