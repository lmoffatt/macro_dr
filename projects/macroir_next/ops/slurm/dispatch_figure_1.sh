#!/bin/bash
# SLURM dispatcher for Figure 1 on a cluster (e.g. dirac): same grid and same
# three-stage functionality as ops/local/dispatch_figure_1_local.sh, one SLURM
# job per (truth, protocol, replica). The payload ops/slurm/run_figure_1.sh
# runs sim -> nan-mask -> paired evidence inside the allocation.
#
# Seeds are deterministic from BASE_SEED and the job index, and the sim seed
# is baked into the label (hence into every output filename), IDENTICALLY to
# the local dispatcher: a local run and a dirac run with the same BASE_SEED
# pair up file by file. Byte-identical results additionally require the same
# thread count on both machines (one RNG stream per OMP thread: CPUS here must
# equal THREADS locally) and a numerically identical BLAS; otherwise expect
# float-level differences, not different draws.
#
# Cluster plumbing copied from projects/eLife_2025/ops/slurm/dispatch_figure_3.sh:
# profile sourcing, per-commit scratch isolation, DEPEND chaining, BIN pinning.
#
# Prereq: a built binary for the cluster:
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
# Usage (from the repo base):
#   projects/macroir_next/ops/slurm/dispatch_figure_1.sh <cluster>   # e.g. dirac
# Tunables via env: TRUTHS, PROTOCOLS, REPLICAS, NCH, TAU_MS, INTERVAL_IN_TAU,
# PRE_TAUS, PULSE_TAUS, POST_TAUS, GAP_TAUS, MEAS_TAUS, SCOUTS, BETA_SIZE,
# MAX_ITER, ADAPT_EVERY, MAX_VALUES, BASE_SEED, CPUS, MEM, TIME, PARTITION,
# ACCOUNT, BIN, DEPEND, RUN_DIR; and for the whole-node packing mode (the
# default, because dirac runs at most 4 jobs per user): PACK (pairs per node,
# default JOBS_PER_NODE from tuning.env else 8; PACK=1 = one pair per job),
# THREADS_PER_FIT (default CPUS_RECOMMENDED from tuning.env else 4),
# NODE_CPUS (default 32).

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../macroir_next/ops/slurm
PROJ="$(readlink -f "$HERE/../..")"                 # .../projects/macroir_next
PAYLOAD="$HERE/run_figure_1.sh"
PACK_PAYLOAD="$HERE/run_figure_1_pack.sh"
SIM_SCRIPT="$(readlink -f "$HERE/../local/figure_1_sim.macroir")"
EVI_SCRIPT="$(readlink -f "$HERE/../local/figure_1_evidence.macroir")"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

# Cluster profiles live with eLife_2025 (referenced, not duplicated).
PROFILE="$(readlink -f "$HERE/../../../eLife_2025/ops/clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[fig1] no such cluster profile: $PROFILE" >&2; exit 1; }
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$PROFILE"
set -e
export MACRODR_PROFILE="$PROFILE"

BIN="${BIN:-$(readlink -f "build/macrodr_cli-${CLUSTER}-current")}"
[ -x "$BIN" ] || {
    echo "[fig1] binary not found: $BIN" >&2
    echo "       build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
    exit 1
}

if ! commit="$("$BIN" --commit)"; then
    echo "[fig1] could not query commit hash: '$BIN --commit' failed" >&2
    exit 1
fi
[ -n "$commit" ] || { echo "[fig1] '$BIN --commit' returned empty" >&2; exit 1; }
run="${RUN_DIR:-$commit}"

DEP_SPEC=""
if [ -n "${DEPEND:-}" ]; then
    case "$DEPEND" in
        *[!0-9]*) DEP_SPEC="$DEPEND" ;;
        *)        DEP_SPEC="afterok:$DEPEND" ;;
    esac
    echo "[fig1] job dependency: --dependency=$DEP_SPEC"
fi

# ---- knobs (same defaults as the local dispatcher) --------------------------
TRUTHS=(${TRUTHS:-CCO COC})
PROTOCOLS=(${PROTOCOLS:-stationary episodic})
REPLICAS="${REPLICAS:-10}"

NCH="${NCH:-1000}"
TAU_MS="${TAU_MS:-16}"
INTERVAL_IN_TAU="${INTERVAL_IN_TAU:-0.1}"
PRE_TAUS="${PRE_TAUS:-1}"
PULSE_TAUS="${PULSE_TAUS:-5}"
POST_TAUS="${POST_TAUS:-5}"
GAP_TAUS="${GAP_TAUS:-20}"
MEAS_TAUS="${MEAS_TAUS:-10}"

SCOUTS="${SCOUTS:-32}"
# BETA_SIZE default is resolved AFTER the tuning hook below, so a discovered
# ladder size (BETA_EFF) can supply it; an explicit env value still wins.
BETA_SIZE="${BETA_SIZE:-}"
MAX_ITER="${MAX_ITER:-30000}"
ADAPT_EVERY="${ADAPT_EVERY:-256}"
MAX_VALUES="${MAX_VALUES:-128}"
BASE_SEED="${BASE_SEED:-910000}"

N_SAMP=$(python3 -c "print(max(1, round(50*$TAU_MS*$INTERVAL_IN_TAU)))")
IPT=$(python3 -c "print(max(1, round(1/$INTERVAL_IN_TAU)))")

# ---- scratch workdir, per commit --------------------------------------------
case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/macroir_next/$run" ;;
esac
mkdir -p "$WORKDIR/logs" "$WORKDIR/data"

# figure_0 tuning: if the probe/characterization ran for this commit, take its
# recommended CPUS as the default (explicit CPUS env still wins). See
# dispatch_figure_0.sh; run the probe before big campaigns.
TUNE_FILE="${TUNE_FILE:-$WORKDIR/figure_0/tuning.env}"
if [ ! -f "$TUNE_FILE" ] && [ -n "${DEPEND:-}" ]; then
    # sbatch freezes --cpus-per-task at submit time: chaining the campaign
    # behind a still-running figure_0 job means submitting WITHOUT its
    # measurement, and the dependency cannot fix that afterwards.
    echo "[fig1] WARNING: no tuning.env at $TUNE_FILE and DEPEND is set." >&2
    echo "[fig1]          These jobs will be submitted with CPUS=${CPUS:-32} (default)," >&2
    echo "[fig1]          ignoring whatever figure_0 measures. To use the measurement," >&2
    echo "[fig1]          wait for figure_0 to finish and dispatch WITHOUT DEPEND." >&2
fi
if [ -f "$TUNE_FILE" ]; then
    # shellcheck source=/dev/null
    source "$TUNE_FILE"
    CPUS="${CPUS:-${CPUS_RECOMMENDED:-32}}"
    # Start the campaign at the DISCOVERED ladder size: the probe measured
    # s/iter there, and starting where the ladder settles avoids the churn
    # (and the load mis-prediction) of growing from a small nominal count.
    BETA_SIZE="${BETA_SIZE:-${BETA_EFF:-16}}"
    echo "[fig1] tuning.env: CPUS=$CPUS beta_size=$BETA_SIZE jobs/node=${JOBS_PER_NODE:-?} s/iter=${SECONDS_PER_ITER:-?}"
    if [ -n "${SECONDS_PER_ITER:-}" ]; then
        est=$(python3 -c "print(round($SECONDS_PER_ITER*${MAX_ITER:-30000}/3600, 2))")
        est2=$(python3 -c "print(round(2*$SECONDS_PER_ITER*${MAX_ITER:-30000}/3600, 2))")
        echo "[fig1] predicted wall-clock: ~${est} h per fit, ~${est2} h per pair (= packed-job wall)"
    fi
fi
BETA_SIZE="${BETA_SIZE:-16}"

# ---- packing: several pairs per whole-node job -------------------------------
# dirac's QOS runs at most 4 jobs per user (MaxJobsPU=4, verified 2026-07-12),
# so one-pair-per-job caps the campaign at 4 x CPUS cores in flight. PACK > 1
# groups PACK pairs into ONE exclusive whole-node job (run_figure_1_pack.sh
# runs them concurrently, THREADS_PER_FIT threads each): 4 job slots then
# carry 4 whole nodes. Defaults come from figure_0's tuning (JOBS_PER_NODE,
# CPUS_RECOMMENDED); PACK=1 restores the one-pair-per-job submission.
PACK="${PACK:-${JOBS_PER_NODE:-8}}"
THREADS_PER_FIT="${THREADS_PER_FIT:-${CPUS_RECOMMENDED:-4}}"
NODE_CPUS="${NODE_CPUS:-32}"

echo "[fig1] commit=$commit  run=$run  WORKDIR=$WORKDIR"
echo "[fig1] tau=${TAU_MS}ms  interval=${INTERVAL_IN_TAU}tau  n_samp=$N_SAMP  intervals/tau=$IPT  NCH=$NCH"
if [ "$PACK" -gt 1 ]; then
    echo "[fig1] packing: $PACK pairs/job x $THREADS_PER_FIT threads/fit (whole node, $NODE_CPUS cpus)"
fi

# ---- N-channel variants of the model inputs (same generator as local) -------
python3 - "$PROJ" "$WORKDIR" "$NCH" << 'EOF'
import csv, math, sys, pathlib
proj, work, nch = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), float(sys.argv[3])
for name in ["scheme_CCO_par", "scheme_CCO_prior",
             "scheme_COC_twin_par", "scheme_COC_twin_prior"]:
    src = proj / "data" / "models" / f"{name}.csv"
    dst = work / "data" / f"{name}_N{int(nch)}.csv"
    rows = list(csv.reader(open(src)))
    for r in rows[1:]:
        if r[2] == "Num_ch_mean":
            r[4] = repr(nch)
            r[5] = repr(math.log10(nch))
    csv.writer(open(dst, "w", newline="")).writerows(rows)
EOF

PAR_CCO="$WORKDIR/data/scheme_CCO_par_N${NCH}.csv"
PRIOR_CCO="$WORKDIR/data/scheme_CCO_prior_N${NCH}.csv"
PAR_COC="$WORKDIR/data/scheme_COC_twin_par_N${NCH}.csv"
PRIOR_COC="$WORKDIR/data/scheme_COC_twin_prior_N${NCH}.csv"

# manifest plumbing for PACK > 1: cells accumulate into manifests/pack_<k>.txt
# and every PACK lines one whole-node job is submitted for that manifest.
mkdir -p "$WORKDIR/manifests"
pack_idx=0
pack_count=0
pack_manifest=""

submit_pack() {
    [ "$pack_count" -gt 0 ] || return 0
    pack_idx=$((pack_idx + 1))
    local jobid
    jobid=$(sbatch --parsable \
        --partition="${PARTITION:-batch}" \
        ${ACCOUNT:+--account="$ACCOUNT"} \
        ${DEPEND:+--dependency="$DEP_SPEC"} \
        --exclusive \
        --cpus-per-task="$NODE_CPUS" \
        --mem="${MEM:-0}" \
        --time="${TIME:-2-00:00:00}" \
        --job-name="f1pack_${pack_idx}" \
        --output="$WORKDIR/logs/pack_${pack_idx}_slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",SIM_SCRIPT="$SIM_SCRIPT",EVI_SCRIPT="$EVI_SCRIPT",RUN_ONE="$PAYLOAD",MANIFEST="$pack_manifest",THREADS_PER_FIT="$THREADS_PER_FIT",PRIOR_CCO="$PRIOR_CCO",PRIOR_COC="$PRIOR_COC",NSAMP="$N_SAMP",SCOUTS="$SCOUTS",BETA_SIZE="$BETA_SIZE",MAX_ITER="$MAX_ITER",ADAPT_EVERY="$ADAPT_EVERY",MAX_VALUES="$MAX_VALUES" \
        "$PACK_PAYLOAD")
    echo "[fig1] pack $pack_idx ($pack_count pairs) -> job $jobid"
    pack_count=0
    pack_manifest=""
}

job=0
for truth in "${TRUTHS[@]}"; do
for prot in "${PROTOCOLS[@]}"; do
for rep in $(seq 1 "$REPLICAS"); do
    job=$((job + 1))
    seed_sim=$((BASE_SEED + 10 * job + 1))
    seed_cco=$((BASE_SEED + 10 * job + 2))
    seed_coc=$((BASE_SEED + 10 * job + 3))
    label="fig1_${truth}_${prot}_rep${rep}_s${seed_sim}"

    case "$truth" in
        CCO) truth_model="scheme_CCO"; truth_par="$PAR_CCO" ;;
        COC) truth_model="scheme_COC"; truth_par="$PAR_COC" ;;
        *) echo "[fig1] unknown truth '$truth' (want CCO or COC)" >&2; exit 1 ;;
    esac

    n1=$((PRE_TAUS * IPT))
    case "$prot" in
        episodic)   n2=$((PULSE_TAUS * IPT)); n3=$((POST_TAUS * IPT)); ag2=10.0; ag3=0.0 ;;
        stationary) n2=$((GAP_TAUS  * IPT)); n3=$((MEAS_TAUS * IPT)); ag2=10.0; ag3=10.0 ;;
        *) echo "[fig1] unknown protocol '$prot'" >&2; exit 1 ;;
    esac
    total=$((n1 + n2 + n3))

    template="$WORKDIR/data/${label}_template.csv"
    python3 - "$template" "$total" << 'EOF'
import sys
path, n = sys.argv[1], int(sys.argv[2])
with open(path, "w") as f:
    f.write("i_step,patch_current\n")
    for i in range(n):
        f.write(f"{i},0\n")
EOF

    if [ "$PACK" -gt 1 ]; then
        if [ -z "$pack_manifest" ]; then
            pack_manifest="$WORKDIR/manifests/pack_$((pack_idx + 1)).txt"
            : > "$pack_manifest"
        fi
        echo "$label $prot $truth_model $truth_par $template $n1 $n2 $n3 $ag2 $ag3 $seed_sim $seed_cco $seed_coc" >> "$pack_manifest"
        pack_count=$((pack_count + 1))
        echo "[fig1] ($job) $label -> pack $((pack_idx + 1))"
        [ "$pack_count" -lt "$PACK" ] || submit_pack
    else
        jobid=$(sbatch --parsable \
            --partition="${PARTITION:-batch}" \
            ${ACCOUNT:+--account="$ACCOUNT"} \
            ${DEPEND:+--dependency="$DEP_SPEC"} \
            --cpus-per-task="${CPUS:-32}" \
            --mem="${MEM:-16G}" \
            --time="${TIME:-2-00:00:00}" \
            --job-name="f1_${truth}_${prot}_r${rep}" \
            --output="$WORKDIR/logs/${label}_slurm-%j.out" \
            --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",SIM_SCRIPT="$SIM_SCRIPT",EVI_SCRIPT="$EVI_SCRIPT",LABEL="$label",PROT="$prot",TRUTH_MODEL="$truth_model",TRUTH_PAR="$truth_par",TEMPLATE="$template",PRIOR_CCO="$PRIOR_CCO",PRIOR_COC="$PRIOR_COC",N1="$n1",N2="$n2",N3="$n3",NSAMP="$N_SAMP",AG2="$ag2",AG3="$ag3",SEED_SIM="$seed_sim",SEED_CCO="$seed_cco",SEED_COC="$seed_coc",SCOUTS="$SCOUTS",BETA_SIZE="$BETA_SIZE",MAX_ITER="$MAX_ITER",ADAPT_EVERY="$ADAPT_EVERY",MAX_VALUES="$MAX_VALUES" \
            "$PAYLOAD")
        echo "[fig1] ($job) $label -> job $jobid"
    fi
done
done
done
if [ "$PACK" -gt 1 ]; then
    submit_pack
    echo "[fig1] submitted $job pairs in $pack_idx whole-node jobs to $CLUSTER, outputs under $WORKDIR"
else
    echo "[fig1] submitted $job jobs to $CLUSTER, outputs under $WORKDIR"
fi
