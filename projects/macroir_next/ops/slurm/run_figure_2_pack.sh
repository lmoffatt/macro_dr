#!/bin/bash
# SLURM payload for one PACKED Figure 2 job: several (truth, fc, replica)
# cells run CONCURRENTLY inside a single whole-node allocation, each cell
# being the unchanged per-cell payload run_figure_2.sh. Same rationale and
# structure as run_figure_1_pack.sh (dirac caps a user at 4 RUNNING jobs).
#
# Required env: CLUSTER, BIN, WORKDIR, MACRODR_PROFILE, SIM_SCRIPT, EVI_SCRIPT,
#   RUN_ONE (path to run_figure_2.sh), MANIFEST (one cell per line:
#   LABEL PROT SIM_POLES FC TEMPLATE N1 N2 N3 AG2 AG3 SEED_SIM SEED_BOX
#   SEED_BES — no field may contain spaces), THREADS_PER_FIT, FIT_NPOLES,
#   TRUTH_PAR, PRIOR, NSAMP, SCOUTS, BETA_SIZE, MAX_ITER, ADAPT_EVERY,
#   MAX_VALUES.

#SBATCH --job-name=fig2pack
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%j.out

set -eo pipefail

: "${CLUSTER:?}" ; : "${BIN:?}" ; : "${WORKDIR:?}" ; : "${MACRODR_PROFILE:?}"
: "${SIM_SCRIPT:?}" ; : "${EVI_SCRIPT:?}" ; : "${RUN_ONE:?}" ; : "${MANIFEST:?}"
: "${THREADS_PER_FIT:?}" ; : "${FIT_NPOLES:?}" ; : "${TRUTH_PAR:?}" ; : "${PRIOR:?}"
: "${NSAMP:?}" ; : "${SCOUTS:?}" ; : "${BETA_SIZE:?}" ; : "${MAX_ITER:?}"
: "${ADAPT_EVERY:?}" ; : "${MAX_VALUES:?}"

[ -f "$MANIFEST" ] || { echo "[fig2-pack] manifest not found: $MANIFEST" >&2; exit 1; }
mkdir -p "$WORKDIR/logs"

echo "[fig2-pack] $(wc -l < "$MANIFEST") cells x ${THREADS_PER_FIT} threads on $(hostname)"

pids=()
labels=()
while read -r LABEL PROT SIM_POLES FC TEMPLATE N1 N2 N3 AG2 AG3 \
              SEED_SIM SEED_BOX SEED_BES; do
    [ -n "$LABEL" ] || continue
    env LABEL="$LABEL" PROT="$PROT" SIM_POLES="$SIM_POLES" FC="$FC" \
        TEMPLATE="$TEMPLATE" N1="$N1" N2="$N2" N3="$N3" AG2="$AG2" AG3="$AG3" \
        SEED_SIM="$SEED_SIM" SEED_BOX="$SEED_BOX" SEED_BES="$SEED_BES" \
        THREADS_PER_FIT="$THREADS_PER_FIT" \
        bash "$RUN_ONE" > "$WORKDIR/logs/${LABEL}_cell.out" 2>&1 &
    pids+=($!)
    labels+=("$LABEL")
done < "$MANIFEST"

fails=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        fails=$((fails + 1))
        echo "[fig2-pack] FAILED: ${labels[$i]} (see logs/${labels[$i]}_cell.out)" >&2
    fi
done

echo "[fig2-pack] done: $((${#pids[@]} - fails))/${#pids[@]} cells ok"
[ "$fails" -eq 0 ] || exit 1
