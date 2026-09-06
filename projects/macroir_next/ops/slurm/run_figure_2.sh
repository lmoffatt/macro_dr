#!/bin/bash
# SLURM payload for one Figure 2 cell: sim (with acquisition filter) ->
# nan-mask (stationary) -> paired MEMBER evidence (canonical box member vs
# qdtf-Bessel member on the same recording). Mirrors run_figure_1.sh; the
# grid axes differ (truth in {box,bessel} via SIM_POLES, filter cutoff FC).
#
# Everything arrives via environment (sbatch --export or run_figure_2_pack).
# Required: CLUSTER, BIN, WORKDIR, MACRODR_PROFILE, SIM_SCRIPT, EVI_SCRIPT,
#   LABEL, PROT, SIM_POLES, FC, FIT_NPOLES, TRUTH_PAR, PRIOR, TEMPLATE,
#   N1, N2, N3, NSAMP, AG2, AG3, SEED_SIM, SEED_BOX, SEED_BES,
#   SCOUTS, BETA_SIZE, MAX_ITER, ADAPT_EVERY, MAX_VALUES.
# Optional: THREADS_PER_FIT (under the pack payload).

#SBATCH --job-name=fig2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%j.out

set -eo pipefail

: "${CLUSTER:?}" ; : "${BIN:?}" ; : "${WORKDIR:?}" ; : "${MACRODR_PROFILE:?}"
: "${SIM_SCRIPT:?}" ; : "${EVI_SCRIPT:?}" ; : "${LABEL:?}" ; : "${PROT:?}"
: "${SIM_POLES:?}" ; : "${FC:?}" ; : "${FIT_NPOLES:?}"
: "${TRUTH_PAR:?}" ; : "${PRIOR:?}" ; : "${TEMPLATE:?}"
: "${N1:?}" ; : "${N2:?}" ; : "${N3:?}" ; : "${NSAMP:?}" ; : "${AG2:?}" ; : "${AG3:?}"
: "${SEED_SIM:?}" ; : "${SEED_BOX:?}" ; : "${SEED_BES:?}"
: "${SCOUTS:?}" ; : "${BETA_SIZE:?}" ; : "${MAX_ITER:?}" ; : "${ADAPT_EVERY:?}"
: "${MAX_VALUES:?}"

unset LOADEDMODULES _LMFILES_ 2>/dev/null || true
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$MACRODR_PROFILE"
set -e

export OMP_NUM_THREADS="${THREADS_PER_FIT:-${SLURM_CPUS_PER_TASK:-16}}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

cd "$WORKDIR"
echo "[fig2-job] label=$LABEL prot=$PROT sim_poles=$SIM_POLES fc=$FC threads=$OMP_NUM_THREADS"

# Same benign-warning filter as run_figure_1.sh (2026-09-05).
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
    "$(printf -- '--truth_model = "scheme_CCO"')" \
    "$(printf -- '--truth_par_file = "%s"' "$TRUTH_PAR")" \
    "$(printf -- '--observations_template = "%s"' "$TEMPLATE")" \
    "$(printf -- '--sim_prefix = "%s"' "$LABEL")" \
    "$(printf -- '--sim_seed = get_number(n=%s)' "$SEED_SIM")" \
    "$(printf -- '--sim_filter_n_poles = get_number(n=%s)' "$SIM_POLES")" \
    "$(printf -- '--sim_filter_cutoff = %s' "$FC")" \
    "$SIM_SCRIPT" 2>&1 | filter_warns

sim_csv="$(ls -t "${LABEL}"_*_simulation.csv 2>/dev/null | head -1)"
[ -n "$sim_csv" ] || { echo "[fig2-job] no simulation csv for $LABEL" >&2; exit 1; }
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

echo "[fig2-job] paired member evidence on $(basename "$obs") (fit ${FIT_NPOLES}p @ ${FC} Hz)"
"$BIN" \
    "$seg1_arg" "$seg2_arg" "$seg3_arg" "$nsamp_arg" "$ag2_arg" "$ag3_arg" \
    "$(printf -- '--observations_file = "%s"' "$obs")" \
    "$(printf -- '--prior_file = "%s"' "$PRIOR")" \
    "$(printf -- '--fit_model = "scheme_CCO"')" \
    "$(printf -- '--idname_box = "%s_fit_box"' "$LABEL")" \
    "$(printf -- '--idname_bessel = "%s_fit_bessel"' "$LABEL")" \
    "$(printf -- '--num_scouts = get_number(n=%s)' "$SCOUTS")" \
    "$(printf -- '--beta_size = get_number(n=%s)' "$BETA_SIZE")" \
    "$(printf -- '--max_iter = get_number(n=%s)' "$MAX_ITER")" \
    "$(printf -- '--adapt_beta_every = get_number(n=%s)' "$ADAPT_EVERY")" \
    "$(printf -- '--max_values = get_number(n=%s)' "$MAX_VALUES")" \
    "$(printf -- '--seed_box = get_number(n=%s)' "$SEED_BOX")" \
    "$(printf -- '--seed_bessel = get_number(n=%s)' "$SEED_BES")" \
    "$(printf -- '--fit_filter_n_poles = get_number(n=%s)' "$FIT_NPOLES")" \
    "$(printf -- '--fit_filter_cutoff = %s' "$FC")" \
    "$EVI_SCRIPT" 2>&1 | filter_warns

echo "[fig2-job] done: $LABEL"
