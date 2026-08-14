#!/usr/bin/env bash
# =============================================================================
# Dispatch figure_3_mle_fisher_only / figure_3_mle_LSE_fisher_only — the analytic
# Gaussian Fisher against its finite-difference counterpart, over the WHOLE roster,
# in ONE lane at HEAD. Stages 2-5 only: no per-group cloud, no empirical capstone.
#
# WHY ONE DISPATCHER, TWO SCRIPTS. The two scripts differ only in the parameter block
# (six free parameters for macro/micro; unitary_current and Current_Noise FIXED for
# least squares, which are flat directions there). The DSL cannot inject a parameter's
# transformation, so the configuration has to live in the file, and the file is picked
# HERE from the algorithm's family. That closes the trap its predecessor documented in
# dispatch_figure_3_LSE_numfim.sh:144-147, where the case listed macro members the
# script could not legitimately run.
#
# WHY IT EXISTS. The Gaussian-vs-numerical comparison exists today only in 433ed13
# (macro NR/R/MR/IR + the superseded NMR), which the manuscript declares superseded and
# quotes nothing from; INR, VR, LSE and ILSE were never measured at all. See the header
# of ../local/figure_3_mle_fisher_only.macroir, including the two traps: the
# non-positivity of F_b at theta_sim is read off the SPECTRUM and never off
# Min_Eigenvalue, and per-replicate indefiniteness overstates it.
#
# WHAT MAKES THIS ONE CHEAP. Its predecessor ran Stages 1-6; Stage 1 (per-group MLE
# over the group_size axis 1/10/100) is what took days, and it feeds only Stage 6.
# Neither is needed here, so GROUP_SIZE / N_BOOT / MIN_GROUPS are NOT injected and must
# not be: the scripts do not define those names.
#
# SEED. SIM_SEED is injected. seed=0 in the DSL is the "draw from random_device"
# sentinel and the drawn value is never written anywhere, so a seed=0 ensemble cannot
# be regenerated — that is the hole 433ed13 has. One seed for the whole sweep also
# pairs the ensemble across algorithms within a cell, which is what makes the
# cross-member comparison a comparison.
#
# Injected names (the scripts must NOT define them): axis_Nchanels, Num_ch, axis_noise,
# current_noise, n_simulations, filepath, sim_seed, algorithm_axis,
# algo_recursive_approximation, algo_averaging_approximation, algo_taylor_approximation,
# algo_family_approximation, algo_variance_form, axis_interval, exp_n_step_/exp_n_samp_,
# axis_h_fim, h_rel_value, gn_max_iter.
#
# Tunables via env: NCHS and N_SIMS (parallel arrays, same length), N_NOISE, N_ALGO,
# H_RELS, SIM_SEED, GN_MAX_ITER, CPUS, MEM, TIME, PARTITION, BIN, RUN_DIR,
# DEPEND (job id to wait for — unset = no dependency).
#
# The defaults are the run this file was written for: the whole macro + least-squares
# roster, noise 0.1, four channel counts, ten thousand recordings.
#   projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
#
# Smoke first — one cheap job, so the parser talks instead of the cluster:
#   N_ALGO="macro_VR" NCHS="10" N_SIMS="8" N_NOISE="0.1" TIME=00:30:00 \
#     projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
#
# The micro arm is a separate sweep over the same lane (different channel counts):
#   N_ALGO="micro_IR micro_R" NCHS="5 10" N_SIMS="10000 10000" N_NOISE="0.1 0.1" \
#     projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
# =============================================================================
set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"          # .../eLife_2025/ops/slurm
WRAPPER="$HERE/run_macroir.sh"
SCRIPT_MACRO="$HERE/../local/figure_3_mle_fisher_only.macroir"
SCRIPT_LSE="$HERE/../local/figure_3_mle_LSE_fisher_only.macroir"
[ -f "$SCRIPT_MACRO" ] || { echo "[dispatch] missing script: $SCRIPT_MACRO" >&2; exit 1; }
[ -f "$SCRIPT_LSE" ]   || { echo "[dispatch] missing script: $SCRIPT_LSE" >&2; exit 1; }

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

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

# Per-commit output isolation. Ask the binary for its baked git hash (the same
# string it stamps as row 1 of every CSV) so two concurrently-running commit
# versions of macrodr_cli write to DISJOINT folders. RUN_DIR overrides the folder.
if ! commit="$("$BIN" --commit)"; then
    echo "[dispatch] could not query commit hash: '$BIN --commit' failed" >&2
    exit 1
fi
[ -n "$commit" ] || { echo "[dispatch] '$BIN --commit' returned empty" >&2; exit 1; }
run="${RUN_DIR:-$commit}"

# Optional SLURM job dependency: hold every dispatched job until another job
# finishes. Set DEPEND to either a bare job id (→ afterok) or a full SLURM
# dependency expression (passed verbatim).
DEP_SPEC=""
if [ -n "${DEPEND:-}" ]; then
    case "$DEPEND" in
        *[!0-9]*) DEP_SPEC="$DEPEND" ;;         # has a non-digit → full SLURM expr
        *)        DEP_SPEC="afterok:$DEPEND" ;; # bare job id → afterok
    esac
    echo "[dispatch] job dependency: --dependency=$DEP_SPEC"
fi

# This run's grid. NCHS and N_SIMS are parallel arrays paired by index.
NCHS=(${NCHS:-10 100 1000 10000})
N_SIMS=(${N_SIMS:-10000 10000 10000 10000})
N_NOISE=(${N_NOISE:-0.1 0.1 0.1 0.1})
N_ALGO=(${N_ALGO:-macro_NR macro_R macro_INR macro_MR macro_VR macro_IR nonlinearsqr nonlinearsqr_g})
# h_rel: relative step for the central-difference numerical Fisher (injected as a
# single-value axis_h_fim so the output matches figure_2). Sweep with H_RELS="1e-4 1e-5 1e-6"
# to show the comparison is not a step-size artifact; that multiplies Stage 3 only.
H_RELS=(${H_RELS:-1e-5})

# Explicit simulation seed — never 0 (see the header). Same value across the sweep, so
# every algorithm in a cell sees the same ensemble.
SIM_SEED="${SIM_SEED:-20260814}"
GN_MAX_ITER="${GN_MAX_ITER:-100}"      # GN iteration cap for the single joint fit

[ "${#NCHS[@]}" -eq "${#N_SIMS[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_SIMS (${#N_SIMS[@]} values) must be the same length" >&2
    exit 1
}
[ "${#NCHS[@]}" -eq "${#N_NOISE[@]}" ] || {
    echo "[dispatch] NCHS (${#NCHS[@]} values) and N_NOISE (${#N_NOISE[@]} values) must be the same length" >&2
    exit 1
}

# Shared output dir on scratch, under the per-commit (or RUN_DIR) folder.
case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/eLife_2025/$run" ;;
esac
mkdir -p "$WORKDIR/figures/data" "$WORKDIR/logs"
echo "[dispatch] commit=${commit}  run=${run}  WORKDIR=${WORKDIR}  (Fisher-only, stages 2-5)"
echo "[dispatch] seed=${SIM_SEED}  h_rel=[${H_RELS[*]}]  algos=[${N_ALGO[*]}]  cells=[${NCHS[*]}]  jobs=$(( ${#N_ALGO[@]} * ${#NCHS[@]} ))"

join_csv()  { local IFS=,; echo "$*"; }
join_qcsv() { local out=""; for v in "$@"; do [ -n "$out" ] && out+=","; out+="\"$v\""; done; echo "$out"; }

# Loop-invariant injections (built once). axis_h_fim must precede h_rel_value.
axis_h_arg=$(printf -- '--axis_h_fim = axis(name= "axis_h_fim", labels= [%s])' "$(join_qcsv "${H_RELS[@]}")")
h_rel_arg=$( printf -- '--h_rel_value = indexed_double_by(axis= axis_h_fim, values=[%s])' "$(join_csv "${H_RELS[@]}")")

# get_number(n=...) → size_t (a bare literal would be a double).
seed_arg=$(   printf -- '--sim_seed = get_number(n=%s)' "$SIM_SEED")
gnmaxit_arg=$(printf -- '--gn_max_iter = get_number(n=%s)' "$GN_MAX_ITER")

for j in "${!N_ALGO[@]}"; do
for i in "${!NCHS[@]}"; do
    nch="${NCHS[$i]}"
    nsim="${N_SIMS[$i]}"
    nnoise="${N_NOISE[$i]}"
    algo="${N_ALGO[$j]}"

    # Map the algorithm label to its (recursive, averaging, taylor, family,
    # variance_form) flags. family: 0=macro, 1=micro, 2=nonlinearsqr (LSE), and it
    # ALSO picks the script, because the two differ in which parameters are free.
    # variance_form: 0=total, 1=residual.
    #
    # macro_VR is macro_MR's mean and gain (averaging=1, start-conditioned, no
    # boundary cross-covariance) with macro_IR's residual interval variance. It is
    # the control that splits the MR->IR gap into a variance step and a gain step.
    case "$algo" in
        macro_NR)       recursive=false; averaging=0 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_R)        recursive=true;  averaging=0 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_INR)      recursive=false; averaging=1 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_MR)       recursive=true;  averaging=1 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_VR)       recursive=true;  averaging=1 ; taylor=false ; family=0 ; variance_form=1 ;;
        macro_IR)       recursive=true;  averaging=2 ; taylor=false ; family=0 ; variance_form=0 ;;
        macro_IRT)      recursive=true;  averaging=2 ; taylor=true  ; family=0 ; variance_form=0 ;;
        micro_R)        recursive=true;  averaging=0 ; taylor=false ; family=1 ; variance_form=0 ;;
        micro_MR)       recursive=true;  averaging=1 ; taylor=false ; family=1 ; variance_form=0 ;;
        micro_IR)       recursive=true;  averaging=2 ; taylor=false ; family=1 ; variance_form=0 ;;
        nonlinearsqr)   recursive=false; averaging=1 ; taylor=false ; family=2 ; variance_form=0 ;;
        nonlinearsqr_g) recursive=false; averaging=0 ; taylor=false ; family=2 ; variance_form=0 ;;
        *) echo "[dispatch] unknown algorithm '$algo' (want macro_{NR,R,INR,MR,VR,IR,IRT}, micro_{R,MR,IR}, nonlinearsqr, nonlinearsqr_g)" >&2; exit 1 ;;
    esac

    # The parameter configuration travels with the family: least squares fixes
    # unitary_current and Current_Noise, everything else fits all six.
    if [ "$family" -eq 2 ]; then SCRIPT="$SCRIPT_LSE"; else SCRIPT="$SCRIPT_MACRO"; fi

    case "$nnoise" in                    # label -> current_noise (vnoise = label / 1000)
        0.05) vnoise=0.00005;;
        0.1)  vnoise=0.0001;;
        0.2)  vnoise=0.0002;;
        0.3)  vnoise=0.0003;;
        0.5)  vnoise=0.0005;;
        1)    vnoise=0.001;;
        10)   vnoise=0.01;;
        100)  vnoise=0.1;;
        1000)  vnoise=1;;
        10000) vnoise=10;;
        100000) vnoise=100;;
        1000000) vnoise=1000;;
        *) echo "[dispatch] unknown noise level '$nnoise' (want 0.05, 0.1, 0.2, 0.3, 0.5, 1, 10, 100, 1000, 10000, 100000, 1000000; vnoise = label/1000)" >&2; exit 1 ;;
    esac

    # printf builds the injections so the DSL double-quotes need no shell escaping.
    axis_arg=$(printf -- '--axis_Nchanels = axis(name= "Num_ch", labels= ["%s"])' "$nch")
    num_arg=$( printf -- '--Num_ch = indexed_double_by(axis= axis_Nchanels, values=[%s])' "$nch")
    nsim_arg=$(printf -- '--n_simulations = get_number(n=%s)' "$nsim")
    fp_arg=$(  printf -- '--filepath = "figures/data/figure_3_fim_nch_%s_nsim_%s_%s_noise_%s"' "$nch" "$nsim" "$algo" "$nnoise")
    axis_noise_arg=$(printf -- '--axis_noise = axis(name= "noise_in_conductance_tau", labels= ["%s"])' "$nnoise")
    current_noise_arg=$(printf -- '--current_noise = indexed_double_by(axis= axis_noise, values=[%s])' "$vnoise")
    axis_algo_arg=$( printf -- '--algorithm_axis = axis(name= "algorithm", labels= ["%s"])' "$algo")
    recursive_arg=$( printf -- '--algo_recursive_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$recursive")
    averaging_arg=$( printf -- '--algo_averaging_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$averaging")
    taylor_arg=$( printf -- '--algo_taylor_approximation = indexed_bool_by(axis= algorithm_axis, values=[%s])' "$taylor")
    family_arg=$(printf -- '--algo_family_approximation = indexed_int_by(axis= algorithm_axis, values=[%s])' "$family")
    vform_arg=$( printf -- '--algo_variance_form = indexed_int_by(axis= algorithm_axis, values=[%s])' "$variance_form")

    # interval_in_tau grid (identical to figure_2). axis_interval precedes exp_n_*.
    # This is the axis the indefiniteness lives on: it switches on at the coarse end.
    axis_interval_arg=$(printf -- '--axis_interval = axis(name= "interval_in_tau", labels= ["1","0.5","0.2","0.1","0.05","0.02","0.01"])')
    exp_step_1_arg=$(printf -- '--exp_n_step_1 = indexed_size_by(axis= axis_interval, values=[2,4,10,20,40,100,200])')
    exp_samp_1_arg=$(printf -- '--exp_n_samp_1 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_2_arg=$(printf -- '--exp_n_step_2 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_2_arg=$(printf -- '--exp_n_samp_2 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')
    exp_step_3_arg=$(printf -- '--exp_n_step_3 = indexed_size_by(axis= axis_interval, values=[4,8,20,40,80,200,400])')
    exp_samp_3_arg=$(printf -- '--exp_n_samp_3 = indexed_size_by(axis= axis_interval, values=[500,250,100,50,25,10,5])')

    # MACRODR_AXIS_SERIAL=1 serializes the internal axis-combo loop so the
    # per-simulation / bootstrap loops become the active OpenMP level (load-bearing).
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        ${DEPEND:+--dependency="$DEP_SPEC"} \
        --cpus-per-task="${CPUS:-32}" \
        --mem="${MEM:-48G}" \
        --time="${TIME:-2-00:00:00}" \
        --job-name="f3fim_${algo}_${nch}c" \
        --output="$WORKDIR/logs/slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",MACRODR_AXIS_SERIAL=1 \
        "$WRAPPER" \
        "$axis_arg" "$num_arg" "$nsim_arg" "$fp_arg" "$seed_arg" \
        "$axis_noise_arg" "$current_noise_arg" \
        "$axis_algo_arg" "$recursive_arg" "$averaging_arg" "$taylor_arg" \
        "$family_arg" "$vform_arg" \
        "$axis_interval_arg" \
        "$exp_step_1_arg" "$exp_samp_1_arg" "$exp_step_2_arg" "$exp_samp_2_arg" \
        "$exp_step_3_arg" "$exp_samp_3_arg" \
        "$axis_h_arg" "$h_rel_arg" \
        "$gnmaxit_arg" \
        "$SCRIPT")

    echo "submitted f3fim_${algo}_nch_${nch}  job=${jobid}  n_sim=${nsim}  noise=${nnoise}  script=$(basename "$SCRIPT")  -> ${WORKDIR}/figures/data/figure_3_fim_nch_${nch}_nsim_${nsim}_${algo}_noise_${nnoise}_*"
done
done
