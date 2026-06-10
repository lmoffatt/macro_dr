#!/bin/bash
# Dispatch figure_2 across Nchannels values — one SLURM job per Nch.
#
# Each job injects the figure's named parameters before the .macroir (the file
# behaves as a function of them). Injection order matters: axis_Nchanels must
# precede Num_ch, which references it.
#
# FILE CONTRACT — figure_2.macroir must NOT define these itself (else its own
# assignment clobbers the injection):
#   axis_Nchanels, Num_ch, n_simulations, filepath
# In particular, comment out the in-file `axis_Nchanels = axis(...)` line.
#
# Prereq: a built binary for the cluster (run once, when code changes):
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
#
# Usage (from the repo base):
#   projects/eLife_2025/ops/slurm/dispatch_figure_2.sh <cluster>     # e.g. dirac
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), CPUS, MEM,
# TIME, PARTITION, BIN. Example:
#   NCHS="10 100" N_SIMS="2048 1024" projects/.../dispatch_figure_2.sh dirac

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
SCRIPT="$HERE/../local/figure_2.macroir"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

# Load the cluster profile for CLUSTER/PARTITION/SCRATCH_MACRO and the repo cd.
# (module loads here are harmless — jobs get modules via run_macroir.sh on the node.)
PROFILE="$(readlink -f "$HERE/../clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[dispatch] no such cluster profile: $PROFILE" >&2; exit 1; }
# /etc/profile.d scripts (e.g. debuginfod.sh runs `cat /dev/null /etc/debuginfod/*.urls`
# on an empty glob → non-zero) and Lmod return non-zero; under `set -eo pipefail`
# that silently aborts before any job is submitted. Guard the env setup, restore after.
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$PROFILE"
set -e
# Passed to run_macroir.sh so it finds the profile from the sbatch spool copy
# (where its own $0 no longer points at the repo).
export MACRODR_PROFILE="$PROFILE"

# Default BIN to this cluster's latest build; override (export BIN=…) to pin one.
BIN="${BIN:-$(readlink -f "build/macrodr_cli-${CLUSTER}-current")}"
[ -x "$BIN" ] || {
    echo "[dispatch] binary not found: $BIN" >&2
    echo "           build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
    exit 1
}

# This run's grid. NCHS and N_SIMS are parallel arrays paired by index:
# job i runs NCHS[i] channels with N_SIMS[i] simulations.
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMS=(${N_SIMS:-1024 1024 1024 1024})
N_NOISE=(${N_NOISE:-0.1 0.1 0.1 0.1})
#N_ALGO=(${N_ALGO:-macro_IR macro_IRT })
N_ALGO=(${N_ALGO:-macro_IR macro_R macro_MR macro_NR macro_NMR})
# h_rel: relative step for the central-difference numerical Fisher. figure_2.macroir
# keeps h_rel_value commented (file contract) and expects it injected — same idiom as
# the local/debug dispatchers. Injected as a single-value axis (axis_h_fim) so the
# output matches them; pass H_RELS="1e-5 1e-6" to sweep h_rel inside each job.
H_RELS=(${H_RELS:-1e-5})


[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}

# Shared output dir on scratch; jobs write nch-distinct filenames into it.
WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/eLife_2025"
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"

# h_rel_value injection (figure_2.macroir references it but keeps it commented per
# the file contract). axis_h_fim must precede h_rel_value, which references it.
# Independent of the (algo, nch) loop, so build once.
join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }
axis_h_arg=$(printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= [%s])' "$(join_qcsv "${H_RELS[@]}")")
h_rel_arg=$( printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$(join_csv "${H_RELS[@]}")")

for j in "${!N_ALGO[@]}"; do

for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"
 
    algo="${N_ALGO[$j]}"
    
    # Map the algorithm label to its (recursive, averaging) settings. Mirrors the
    # commented-out indexed_bool_by/indexed_int_by rows in figure_2.macroir.
    case "$algo" in
        macro_NR)  recursive=false; averaging=0 ;  taylor=false ; micro=false ;;
        macro_R)   recursive=true;  averaging=0 ;taylor=false ; micro=false ;;
        macro_NMR) recursive=false; averaging=1 ;taylor=false ; micro=false ;;
        macro_MR)  recursive=true;  averaging=1 ;taylor=false ; micro=false ;;
        macro_IR)  recursive=true;  averaging=2 ;taylor=false ; micro=false ;;
        macro_IRT)  recursive=true;  averaging=2 ;taylor=true ; micro=false ;;
        micro_R)   recursive=true;  averaging=0 ;taylor=false ; micro=true ;;
        micro_MR)  recursive=true;  averaging=1 ;taylor=false ; micro=true ;;
        micro_IR)  recursive=true;  averaging=2 ;taylor=false ; micro=true ;;
        
        *) echo "[dispatch] unknown algorithm '$algo' (want macro_{NR,R,NMR,MR,IR, IRT})" >&2; exit 1 ;;
    esac

   case "$nnoise" in
        0.1)  vnoise=0.0001;;
        1)  vnoise=0.001;;
        10)  vnoise=0.01;;
        100)  vnoise=0.1;;
        *) echo "[dispatch] unknown noise level '$nnoise' (want 0.1, 1, 10, 100)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    # get_number(n=...) → size_t; a bare literal would be a double and
    # simulate's n_simulations expects unsigned long.
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_2_nch_%s_nsim_%s_%s_noise_%s"' "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    # algorithm — injected the same way as noise: the axis, plus the (recursive,
    # averaging) settings the label maps to. algorithm_axis must come first; the
    # other two reference it.
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    micro_arg=$( printf -- '--algo_micro_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$micro")

    # interval_in_tau — figure_2.macroir now has axis_interval + the experiment
    # step/sample sizes injection-ready (commented), so the dispatcher injects
    # them. Production = the full 7-interval grid (identical to the values that
    # used to be hardcoded in figure_2.macroir). axis_interval must precede the
    # exp_n_* (which reference it).
    axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= ["1","0.5","0.2","0.1","0.05","0.02","0.01"])')
    exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2,4,10,20,40,100,200])')
    exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')

    # MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop so the
    # per-simulation / bootstrap loops become the active OpenMP level. Without
    # it the combo loop stays parallel and memory ∝ concurrent combos → OOM
    # (cf. jobs 110042–110045, killed at ~129 GB). Load-bearing, not optional.
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-96G}" \
        --time="${TIME:-2-00:00:00}" \
        --job-name="f_${nch}c_${nsim}s" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",MACRODR_AXIS_SERIAL=1 \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$micro_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$axis_h_arg" "$h_rel_arg" \
        "$SCRIPT")

    echo "submitted fig2_nch_${nch}  job=${jobid}  n_sim=${nsim}  -> ${WORKDIR}/figures/data/figure_2_nch_${nch}_nsim_${nsim}_*"
done
done

echo
echo "watch: squeue -u \$(whoami)   logs: tail -f ${WORKDIR}/logs/slurm-<jobid>.out"
