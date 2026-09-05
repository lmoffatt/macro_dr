#!/bin/bash
# SLURM payload for one PACKED Figure 1 job: several (truth, protocol,
# replica) pairs run CONCURRENTLY inside a single whole-node allocation.
#
# Why: dirac's QOS caps a user at 4 RUNNING jobs (MaxJobsPU=4, verified
# 2026-07-12), so one-pair-per-job wastes the quota: 4 jobs x 4 CPUs would
# keep 16 cores in flight on a 600-core cluster. Packing PACK pairs into a
# whole node turns the same 4 job slots into 4 x 32 = 128 cores, each pair at
# its measured-efficient thread count. Same answer eLife's figure_3 gave to
# this cap (few fat jobs, parallelism inside the allocation), realized here
# by running the UNCHANGED per-pair payload run_figure_1.sh in background,
# once per manifest line, and waiting.
#
# Required env: CLUSTER, BIN, WORKDIR, MACRODR_PROFILE, SIM_SCRIPT, EVI_SCRIPT,
#   RUN_ONE (path to run_figure_1.sh: $0 is spooled by sbatch, so the real
#   path must come from the dispatcher), MANIFEST (file, one pair per line:
#   LABEL PROT TRUTH_MODEL TRUTH_PAR TEMPLATE N1 N2 N3 AG2 AG3 SEED_SIM
#   SEED_CCO SEED_COC — no field may contain spaces), THREADS_PER_FIT,
#   PRIOR_CCO, PRIOR_COC, NSAMP, SCOUTS, BETA_SIZE, MAX_ITER, ADAPT_EVERY,
#   MAX_VALUES.
#
# Each pair writes its own log to $WORKDIR/logs/<label>_pair.out. If any pair
# fails the job exits non-zero (visible to afterok chains), but the finished
# pairs' outputs stay on disk.

#SBATCH --job-name=fig1pack
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%j.out

set -eo pipefail

: "${CLUSTER:?}" ; : "${BIN:?}" ; : "${WORKDIR:?}" ; : "${MACRODR_PROFILE:?}"
: "${SIM_SCRIPT:?}" ; : "${EVI_SCRIPT:?}" ; : "${RUN_ONE:?}" ; : "${MANIFEST:?}"
: "${THREADS_PER_FIT:?}" ; : "${PRIOR_CCO:?}" ; : "${PRIOR_COC:?}" ; : "${NSAMP:?}"
: "${SCOUTS:?}" ; : "${BETA_SIZE:?}" ; : "${MAX_ITER:?}" ; : "${ADAPT_EVERY:?}"
: "${MAX_VALUES:?}"

[ -f "$MANIFEST" ] || { echo "[fig1-pack] manifest not found: $MANIFEST" >&2; exit 1; }
mkdir -p "$WORKDIR/logs"

echo "[fig1-pack] $(wc -l < "$MANIFEST") pairs x ${THREADS_PER_FIT} threads on $(hostname)"

pids=()
labels=()
while read -r LABEL PROT TRUTH_MODEL TRUTH_PAR TEMPLATE N1 N2 N3 AG2 AG3 \
              SEED_SIM SEED_CCO SEED_COC; do
    [ -n "$LABEL" ] || continue
    env LABEL="$LABEL" PROT="$PROT" TRUTH_MODEL="$TRUTH_MODEL" \
        TRUTH_PAR="$TRUTH_PAR" TEMPLATE="$TEMPLATE" \
        N1="$N1" N2="$N2" N3="$N3" AG2="$AG2" AG3="$AG3" \
        SEED_SIM="$SEED_SIM" SEED_CCO="$SEED_CCO" SEED_COC="$SEED_COC" \
        THREADS_PER_FIT="$THREADS_PER_FIT" \
        bash "$RUN_ONE" > "$WORKDIR/logs/${LABEL}_pair.out" 2>&1 &
    pids+=($!)
    labels+=("$LABEL")
done < "$MANIFEST"

fails=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        fails=$((fails + 1))
        echo "[fig1-pack] FAILED: ${labels[$i]} (see logs/${labels[$i]}_pair.out)" >&2
    fi
done

echo "[fig1-pack] done: $((${#pids[@]} - fails))/${#pids[@]} pairs ok"
[ "$fails" -eq 0 ] || exit 1
