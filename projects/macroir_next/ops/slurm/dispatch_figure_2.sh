#!/bin/bash
# SLURM dispatcher for Figure 2 on a cluster: the MEMBER-confusion matrix of
# the Bessel plan (M4). Grid {truth: box, bessel} x {FCS} x REPLICAS, one
# cell = sim (filtered truth) -> mask -> paired member evidence, packed into
# whole-node jobs like dispatch_figure_1.sh (dirac runs at most 4 jobs per
# user; PACK cells per exclusive node, THREADS_PER_FIT threads per fit).
#
# There is usually NO tuning.env for this campaign's commit (the figure_0
# probe measures the figure_1 workload, and the qdtf member costs more per
# iteration), so PACK/THREADS_PER_FIT default to figure_1's measured optimum
# (8 x 4 on dirac, 2026-09-05) and are meant to be overridden explicitly.
#
# NOTE on chaining behind a running figure_1 campaign: no DEPEND is needed to
# wait for resources — the QOS job cap serializes; submitted jobs start as
# figure_1 drains. If strict ordering is wanted use afterany (afterok would
# cancel the chain if any figure_1 pack reports a failed pair).
#
# Prereq: a built binary for the cluster AT THE COMMIT WITH THE QDTF WIRING:
#   projects/eLife_2025/ops/build_cluster.sh <cluster>
# Usage (from the repo base):
#   projects/macroir_next/ops/slurm/dispatch_figure_2.sh <cluster>
# Tunables via env: TRUTHS, FCS, REPLICAS, PROT, FIT_NPOLES, NCH, TAU_MS,
# INTERVAL_IN_TAU, PRE_TAUS, PULSE_TAUS, POST_TAUS, GAP_TAUS, MEAS_TAUS,
# SCOUTS, BETA_SIZE, MAX_ITER, ADAPT_EVERY, MAX_VALUES, BASE_SEED, PACK,
# THREADS_PER_FIT, NODE_CPUS, MEM, TIME, PARTITION, ACCOUNT, BIN, DEPEND,
# RUN_DIR.

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../macroir_next/ops/slurm
PROJ="$(readlink -f "$HERE/../..")"                 # .../projects/macroir_next
PAYLOAD="$HERE/run_figure_2.sh"
PACK_PAYLOAD="$HERE/run_figure_2_pack.sh"
SIM_SCRIPT="$(readlink -f "$HERE/../local/figure_2_sim.macroir")"
EVI_SCRIPT="$(readlink -f "$HERE/../local/figure_2_evidence.macroir")"

CLUSTER="${1:?Usage: $0 <cluster>   (e.g. dirac)}"

PROFILE="$(readlink -f "$HERE/../../../eLife_2025/ops/clusters/${CLUSTER}.sh")"
[ -f "$PROFILE" ] || { echo "[fig2] no such cluster profile: $PROFILE" >&2; exit 1; }
set +e
[ -f /etc/profile ] && source /etc/profile
# shellcheck source=/dev/null
source "$PROFILE"
set -e
export MACRODR_PROFILE="$PROFILE"

# Deferred-binary mode: with DEPEND set (e.g. afterok:<build job>) the
# dispatcher may run BEFORE the build finishes. The commit is already frozen
# by git at submission (build_cluster.sh tags the build with `git rev-parse
# --short HEAD`), so both the WORKDIR name and the future binary path are
# derivable without the binary. The jobs only touch BIN at RUN time, after
# the dependency released them. Requires a CLEAN tree (a dirty tree would
# stamp "<hash>-dirty" and the WORKDIR name would lie).
BIN_DEFAULT="build/macrodr_cli-${CLUSTER}-current"
if [ -n "${BIN:-}" ] || [ -x "$BIN_DEFAULT" ] || [ -z "${DEPEND:-}" ]; then
    BIN="${BIN:-$(readlink -f "$BIN_DEFAULT")}"
    [ -x "$BIN" ] || {
        echo "[fig2] binary not found: $BIN" >&2
        echo "       build first: projects/eLife_2025/ops/build_cluster.sh ${CLUSTER}" >&2
        echo "       or submit the build and chain: DEPEND=afterok:<build_jobid> $0 ${CLUSTER}" >&2
        exit 1
    }
    if ! commit="$("$BIN" --commit)"; then
        echo "[fig2] could not query commit hash: '$BIN --commit' failed" >&2
        exit 1
    fi
    [ -n "$commit" ] || { echo "[fig2] '$BIN --commit' returned empty" >&2; exit 1; }
else
    commit="$(git rev-parse --short HEAD)"
    if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
        echo "[fig2] ERROR: deferred-binary dispatch needs a CLEAN tree (the build" >&2
        echo "       would stamp ${commit}-dirty and the WORKDIR name would not match)." >&2
        exit 1
    fi
    BIN="$(pwd)/build/${CLUSTER}-${commit}/macrodr_cli"
    echo "[fig2] deferred binary: $BIN (does not exist yet; jobs run behind $DEPEND)"
    echo "[fig2] make sure the depended-on build job was submitted at THIS commit ($commit)."
fi
run="${RUN_DIR:-$commit}"

DEP_SPEC=""
if [ -n "${DEPEND:-}" ]; then
    case "$DEPEND" in
        *[!0-9]*) DEP_SPEC="$DEPEND" ;;
        *)        DEP_SPEC="afterany:$DEPEND" ;;
    esac
    echo "[fig2] job dependency: --dependency=$DEP_SPEC"
fi

# ---- knobs -------------------------------------------------------------------
TRUTHS=(${TRUTHS:-box bessel})
FCS=(${FCS:-600 200})
REPLICAS="${REPLICAS:-10}"
PROT="${PROT:-stationary}"
FIT_NPOLES="${FIT_NPOLES:-4}"

NCH="${NCH:-1000}"
TAU_MS="${TAU_MS:-16}"
INTERVAL_IN_TAU="${INTERVAL_IN_TAU:-0.1}"
PRE_TAUS="${PRE_TAUS:-1}"
PULSE_TAUS="${PULSE_TAUS:-5}"
POST_TAUS="${POST_TAUS:-5}"
GAP_TAUS="${GAP_TAUS:-20}"
MEAS_TAUS="${MEAS_TAUS:-10}"

SCOUTS="${SCOUTS:-32}"
BETA_SIZE="${BETA_SIZE:-16}"
MAX_ITER="${MAX_ITER:-30000}"
ADAPT_EVERY="${ADAPT_EVERY:-256}"
MAX_VALUES="${MAX_VALUES:-128}"
BASE_SEED="${BASE_SEED:-920000}"

PACK="${PACK:-8}"
THREADS_PER_FIT="${THREADS_PER_FIT:-4}"
NODE_CPUS="${NODE_CPUS:-32}"

N_SAMP=$(python3 -c "print(max(1, round(50*$TAU_MS*$INTERVAL_IN_TAU)))")
IPT=$(python3 -c "print(max(1, round(1/$INTERVAL_IN_TAU)))")

case "$run" in
    /*) WORKDIR="$run" ;;
    *)  WORKDIR="${SCRATCH_MACRO:-/scratch/$(whoami)/macro_dr}/macroir_next/$run" ;;
esac
mkdir -p "$WORKDIR/logs" "$WORKDIR/data" "$WORKDIR/manifests"

echo "[fig2] commit=$commit  run=$run  WORKDIR=$WORKDIR"
echo "[fig2] grid: truths=${TRUTHS[*]} fcs=${FCS[*]} replicas=$REPLICAS prot=$PROT"
if [ "$PACK" -gt 1 ]; then
    echo "[fig2] packing: $PACK cells/job x $THREADS_PER_FIT threads/fit (whole node, $NODE_CPUS cpus)"
fi

# ---- N-channel variants (scheme_CCO only: both arms fit the same model) -----
python3 - "$PROJ" "$WORKDIR" "$NCH" << 'EOF'
import csv, math, sys, pathlib
proj, work, nch = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), float(sys.argv[3])
for name in ["scheme_CCO_par", "scheme_CCO_prior"]:
    src = proj / "data" / "models" / f"{name}.csv"
    dst = work / "data" / f"{name}_N{int(nch)}.csv"
    rows = list(csv.reader(open(src)))
    for r in rows[1:]:
        if r[2] == "Num_ch_mean":
            r[4] = repr(nch)
            r[5] = repr(math.log10(nch))
    csv.writer(open(dst, "w", newline="")).writerows(rows)
EOF
TRUTH_PAR="$WORKDIR/data/scheme_CCO_par_N${NCH}.csv"
PRIOR="$WORKDIR/data/scheme_CCO_prior_N${NCH}.csv"

n1=$((PRE_TAUS * IPT))
case "$PROT" in
    episodic)   n2=$((PULSE_TAUS * IPT)); n3=$((POST_TAUS * IPT)); ag2=10.0; ag3=0.0 ;;
    stationary) n2=$((GAP_TAUS  * IPT)); n3=$((MEAS_TAUS * IPT)); ag2=10.0; ag3=10.0 ;;
    *) echo "[fig2] unknown protocol '$PROT'" >&2; exit 1 ;;
esac
total=$((n1 + n2 + n3))

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
        --job-name="f2pack_${pack_idx}" \
        --output="$WORKDIR/logs/f2pack_${pack_idx}_slurm-%j.out" \
        --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",SIM_SCRIPT="$SIM_SCRIPT",EVI_SCRIPT="$EVI_SCRIPT",RUN_ONE="$PAYLOAD",MANIFEST="$pack_manifest",THREADS_PER_FIT="$THREADS_PER_FIT",FIT_NPOLES="$FIT_NPOLES",TRUTH_PAR="$TRUTH_PAR",PRIOR="$PRIOR",NSAMP="$N_SAMP",SCOUTS="$SCOUTS",BETA_SIZE="$BETA_SIZE",MAX_ITER="$MAX_ITER",ADAPT_EVERY="$ADAPT_EVERY",MAX_VALUES="$MAX_VALUES" \
        "$PACK_PAYLOAD")
    echo "[fig2] pack $pack_idx ($pack_count cells) -> job $jobid"
    pack_count=0
    pack_manifest=""
}

job=0
for truth in "${TRUTHS[@]}"; do
for fc in "${FCS[@]}"; do
for rep in $(seq 1 "$REPLICAS"); do
    job=$((job + 1))
    seed_sim=$((BASE_SEED + 10 * job + 1))
    seed_box=$((BASE_SEED + 10 * job + 2))
    seed_bes=$((BASE_SEED + 10 * job + 3))
    label="fig2_${truth}_fc${fc}_rep${rep}_s${seed_sim}"

    case "$truth" in
        box)    sim_poles=0 ;;
        bessel) sim_poles="$FIT_NPOLES" ;;
        *) echo "[fig2] unknown truth '$truth' (want box or bessel)" >&2; exit 1 ;;
    esac

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
            pack_manifest="$WORKDIR/manifests/f2pack_$((pack_idx + 1)).txt"
            : > "$pack_manifest"
        fi
        echo "$label $PROT $sim_poles $fc $template $n1 $n2 $n3 $ag2 $ag3 $seed_sim $seed_box $seed_bes" >> "$pack_manifest"
        pack_count=$((pack_count + 1))
        echo "[fig2] ($job) $label -> pack $((pack_idx + 1))"
        [ "$pack_count" -lt "$PACK" ] || submit_pack
    else
        jobid=$(sbatch --parsable \
            --partition="${PARTITION:-batch}" \
            ${ACCOUNT:+--account="$ACCOUNT"} \
            ${DEPEND:+--dependency="$DEP_SPEC"} \
            --cpus-per-task="${CPUS:-32}" \
            --mem="${MEM:-16G}" \
            --time="${TIME:-2-00:00:00}" \
            --job-name="f2_${truth}_fc${fc}_r${rep}" \
            --output="$WORKDIR/logs/${label}_slurm-%j.out" \
            --export=ALL,CLUSTER="$CLUSTER",BIN="$BIN",WORKDIR="$WORKDIR",MACRODR_PROFILE="$PROFILE",SIM_SCRIPT="$SIM_SCRIPT",EVI_SCRIPT="$EVI_SCRIPT",LABEL="$label",PROT="$PROT",SIM_POLES="$sim_poles",FC="$fc",FIT_NPOLES="$FIT_NPOLES",TRUTH_PAR="$TRUTH_PAR",PRIOR="$PRIOR",TEMPLATE="$template",N1="$n1",N2="$n2",N3="$n3",NSAMP="$N_SAMP",AG2="$ag2",AG3="$ag3",SEED_SIM="$seed_sim",SEED_BOX="$seed_box",SEED_BES="$seed_bes",SCOUTS="$SCOUTS",BETA_SIZE="$BETA_SIZE",MAX_ITER="$MAX_ITER",ADAPT_EVERY="$ADAPT_EVERY",MAX_VALUES="$MAX_VALUES" \
            "$PAYLOAD")
        echo "[fig2] ($job) $label -> job $jobid"
    fi
done
done
done
if [ "$PACK" -gt 1 ]; then
    submit_pack
    echo "[fig2] submitted $job cells in $pack_idx whole-node jobs to $CLUSTER, outputs under $WORKDIR"
else
    echo "[fig2] submitted $job jobs to $CLUSTER, outputs under $WORKDIR"
fi
