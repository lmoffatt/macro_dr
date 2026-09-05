#!/bin/bash
# SLURM payload for one Figure 1 job: the SAME three stages the local
# dispatcher runs, executed inside a single allocation:
#   1. figure_1_sim.macroir       simulate the replica recording
#   2. (stationary only)          nan-mask the equilibration segment (awk,
#                                 same line as the local dispatcher, so the
#                                 masked files are byte-comparable)
#   3. figure_1_evidence.macroir  paired evidence, both fits on the same data
#
# Everything arrives via environment (sbatch --export); no positional args.
# Required: CLUSTER, BIN, WORKDIR, MACRODR_PROFILE, SIM_SCRIPT, EVI_SCRIPT,
#   LABEL, PROT, TRUTH_MODEL, TRUTH_PAR, TEMPLATE, PRIOR_CCO, PRIOR_COC,
#   N1, N2, N3, NSAMP, AG2, AG3, SEED_SIM, SEED_CCO, SEED_COC,
#   SCOUTS, BETA_SIZE, MAX_ITER, ADAPT_EVERY, MAX_VALUES.
#
# Seeds are baked into LABEL by the dispatcher, so every output filename
# carries the seed and pairs with the local run of the same BASE_SEED.

#SBATCH --job-name=fig1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%j.out

set -eo pipefail

: "${CLUSTER:?}" ; : "${BIN:?}" ; : "${WORKDIR:?}" ; : "${MACRODR_PROFILE:?}"
: "${SIM_SCRIPT:?}" ; : "${EVI_SCRIPT:?}" ; : "${LABEL:?}" ; : "${PROT:?}"
: "${TRUTH_MODEL:?}" ; : "${TRUTH_PAR:?}" ; : "${TEMPLATE:?}"
: "${PRIOR_CCO:?}" ; : "${PRIOR_COC:?}"
: "${N1:?}" ; : "${N2:?}" ; : "${N3:?}" ; : "${NSAMP:?}" ; : "${AG2:?}" ; : "${AG3:?}"
: "${SEED_SIM:?}" ; : "${SEED_CCO:?}" ; : "${SEED_COC:?}"
: "${SCOUTS:?}" ; : "${BETA_SIZE:?}" ; : "${MAX_ITER:?}" ; : "${ADAPT_EVERY:?}"
: "${MAX_VALUES:?}"

# Same env dance as run_macroir.sh: drop inherited Lmod tracking, then source
# /etc/profile and the cluster profile with -e off (OpenHPC profile.d scripts
# and Lmod return non-zero on empty globs).
unset LOADEDMODULES _LMFILES_ 2>/dev/null || true
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$MACRODR_PROFILE"
set -e

# THREADS_PER_FIT overrides the allocation size: under run_figure_1_pack.sh
# several of these payloads share one whole-node allocation, and each must
# use its slice, not SLURM_CPUS_PER_TASK (= the whole node).
export OMP_NUM_THREADS="${THREADS_PER_FIT:-${SLURM_CPUS_PER_TASK:-16}}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

cd "$WORKDIR"
echo "[fig1-job] label=$LABEL prot=$PROT threads=$OMP_NUM_THREADS bin=$BIN"

# Drop the flood of benign range-canary warnings ([warn] check_g*_in_range,
# to_Probability*, and their torn fragments from concurrent stderr writes):
# measured 2026-09-05 at ~2300 lines/iteration under MCMC, 267 MB per
# 4-minute slice, i.e. tens of GB per job at campaign scale. Real errors do
# not carry these markers and pass through. The || true keeps an all-filtered
# (empty) stream from tripping pipefail; the binary's own exit status still
# propagates through the pipe.
filter_warns() {
    grep -v -e '\[warn\]' -e 'check_g' -e 'to_Probability' -e 'to_Covariance' \
            -e 'warn=' -e 'err=' -e 'excursion' -e 'P_ij' -e 'gsqr' || true
}

seg1_arg=$(printf -- '--exp_n_step_1 = get_number(n=%s)' "$N1")
seg2_arg=$(printf -- '--exp_n_step_2 = get_number(n=%s)' "$N2")
seg3_arg=$(printf -- '--exp_n_step_3 = get_number(n=%s)' "$N3")
nsamp_arg=$(printf -- '--exp_n_samp = get_number(n=%s)' "$NSAMP")
ag2_arg=$(printf -- '--agonist_2 = %s' "$AG2")
ag3_arg=$(printf -- '--agonist_3 = %s' "$AG3")

"$BIN" \
    "$seg1_arg" "$seg2_arg" "$seg3_arg" "$nsamp_arg" "$ag2_arg" "$ag3_arg" \
    "$(printf -- '--truth_model = "%s"' "$TRUTH_MODEL")" \
    "$(printf -- '--truth_par_file = "%s"' "$TRUTH_PAR")" \
    "$(printf -- '--observations_template = "%s"' "$TEMPLATE")" \
    "$(printf -- '--sim_prefix = "%s"' "$LABEL")" \
    "$(printf -- '--sim_seed = get_number(n=%s)' "$SEED_SIM")" \
    "$SIM_SCRIPT" 2>&1 | filter_warns

sim_csv="$(ls -t "${LABEL}"_*_simulation.csv 2>/dev/null | head -1)"
[ -n "$sim_csv" ] || { echo "[fig1-job] no simulation csv for $LABEL" >&2; exit 1; }
sim_csv="$(readlink -f "$sim_csv")"

if [ "$PROT" = "stationary" ]; then
    masked="$WORKDIR/data/${LABEL}_masked.csv"
    awk -F, -v a="$N1" -v b="$((N1 + N2))" \
        'NR==1{print; next} { if ($1+0>=a && $1+0<b) print $1",nan"; else print }' \
        "$sim_csv" > "$masked"
    obs="$masked"
else
    obs="$sim_csv"
fi

echo "[fig1-job] paired evidence on $(basename "$obs")"
"$BIN" \
    "$seg1_arg" "$seg2_arg" "$seg3_arg" "$nsamp_arg" "$ag2_arg" "$ag3_arg" \
    "$(printf -- '--observations_file = "%s"' "$obs")" \
    "$(printf -- '--prior_cco_file = "%s"' "$PRIOR_CCO")" \
    "$(printf -- '--prior_coc_file = "%s"' "$PRIOR_COC")" \
    "$(printf -- '--idname_cco = "%s_fit_CCO"' "$LABEL")" \
    "$(printf -- '--idname_coc = "%s_fit_COC"' "$LABEL")" \
    "$(printf -- '--num_scouts = get_number(n=%s)' "$SCOUTS")" \
    "$(printf -- '--beta_size = get_number(n=%s)' "$BETA_SIZE")" \
    "$(printf -- '--max_iter = get_number(n=%s)' "$MAX_ITER")" \
    "$(printf -- '--adapt_beta_every = get_number(n=%s)' "$ADAPT_EVERY")" \
    "$(printf -- '--max_values = get_number(n=%s)' "$MAX_VALUES")" \
    "$(printf -- '--seed_cco = get_number(n=%s)' "$SEED_CCO")" \
    "$(printf -- '--seed_coc = get_number(n=%s)' "$SEED_COC")" \
    "$EVI_SCRIPT" 2>&1 | filter_warns

echo "[fig1-job] done: $LABEL"
