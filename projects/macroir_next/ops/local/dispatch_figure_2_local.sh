#!/bin/bash
# Local dispatcher for Figure 2: the MEMBER-confusion matrix of the Bessel
# plan (M4). Grid: {truth: box, bessel} x {FCS: filter cutoffs} x REPLICAS,
# three stages per cell exactly like dispatch_figure_1_local.sh:
#   1. figure_2_sim.macroir       simulate the replica (sample_bessel: causal
#                                 filter riding the substep path; n_poles=0
#                                 for box truth, bit-identical to figure_1's)
#   2. (stationary) nan-mask      same awk as figure_1
#   3. figure_2_evidence.macroir  paired MEMBER fits on the same recording:
#                                 canonical box member vs qdtf-Bessel member,
#                                 same model, same prior, same data
#
# Cutoff choice: the box member's damage scales with fc*Delta. At the
# figure_1 protocol (tau=16 ms, interval 0.1 tau -> Delta=1.6 ms) the rig-like
# 10 kHz cutoff gives fc*Delta=16 (filter invisible); the defaults 600 and
# 200 Hz give fc*Delta = 1 and 0.32, where the filter bites. Rig-realistic
# cutoffs pair with native-rate sampling and faster kinetics — a later,
# bigger campaign; this lane's first job is to see the machinery WORK.
#
# Fit arm filter: always 4 poles at the cell's cutoff (rig metadata, never
# fitted); truth bessel simulates with the same (4, fc). Truth box ignores
# the cutoff (n_poles=0).
#
# Usage (from the repo base, after building gcc-release):
#   projects/macroir_next/ops/local/dispatch_figure_2_local.sh
# Smoke first: REPLICAS=1 (default) TRUTHS="bessel" FCS="600" MAX_ITER=200 ...
# Tunables via env: TRUTHS, FCS, REPLICAS, PROT, NCH, TAU_MS, INTERVAL_IN_TAU,
# PRE_TAUS, PULSE_TAUS, POST_TAUS, GAP_TAUS, MEAS_TAUS, SCOUTS, BETA_SIZE,
# MAX_ITER, ADAPT_EVERY, MAX_VALUES, BASE_SEED, THREADS, FIT_NPOLES, BIN,
# WORKDIR.

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../macroir_next/ops/local
BASE="$(readlink -f "$HERE/../../../..")"           # repo base
PROJ="$(readlink -f "$HERE/../..")"                 # .../projects/macroir_next
SIM_SCRIPT="$HERE/figure_2_sim.macroir"
EVI_SCRIPT="$HERE/figure_2_evidence.macroir"

BIN="${BIN:-$(readlink -f "$BASE/build/gcc-release/macrodr_cli")}"
[ -x "$BIN" ] || {
    echo "[fig2] binary not found: $BIN" >&2
    echo "       build first:  cmake --build --preset gcc-release" >&2
    exit 1
}
commit="$("$BIN" --commit)" || { echo "[fig2] '$BIN --commit' failed" >&2; exit 1; }
WORKDIR="${WORKDIR:-$PROJ/runs/$commit}"
mkdir -p "$WORKDIR/logs" "$WORKDIR/data"
WORKDIR="$(readlink -f "$WORKDIR")"

# ---- knobs -------------------------------------------------------------------
TRUTHS=(${TRUTHS:-box bessel})
FCS=(${FCS:-600 200})
REPLICAS="${REPLICAS:-1}"
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

N_SAMP=$(python3 -c "print(max(1, round(50*$TAU_MS*$INTERVAL_IN_TAU)))")
IPT=$(python3 -c "print(max(1, round(1/$INTERVAL_IN_TAU)))")

export OMP_NUM_THREADS="${THREADS:-$(nproc)}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${BLAS_THREADS:-1}"
export BLIS_NUM_THREADS="${BLAS_THREADS:-1}"

# ---- N-channel variants of the model inputs (same generator as figure_1) ----
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
PAR="$WORKDIR/data/scheme_CCO_par_N${NCH}.csv"
PRIOR="$WORKDIR/data/scheme_CCO_prior_N${NCH}.csv"

cd "$WORKDIR"

# Drop the flood of benign range-canary warnings; real errors pass through
# (same filter, same rationale as dispatch_figure_1_local.sh, 2026-09-05).
filter_warns() {
    grep -v -e '\[warn\]' -e 'check_g' -e 'to_Probability' -e 'to_Covariance' \
            -e 'warn=' -e 'err=' -e 'excursion' -e 'P_ij' -e 'gsqr' || true
}

echo "[fig2] bin=$BIN commit=$commit workdir=$WORKDIR threads=$OMP_NUM_THREADS"
echo "[fig2] grid: truths=${TRUTHS[*]} fcs=${FCS[*]} replicas=$REPLICAS prot=$PROT"

n1=$((PRE_TAUS * IPT))
case "$PROT" in
    episodic)   n2=$((PULSE_TAUS * IPT)); n3=$((POST_TAUS * IPT)); ag2=10.0; ag3=0.0 ;;
    stationary) n2=$((GAP_TAUS  * IPT)); n3=$((MEAS_TAUS * IPT)); ag2=10.0; ag3=10.0 ;;
    *) echo "[fig2] unknown protocol '$PROT'" >&2; exit 1 ;;
esac
total=$((n1 + n2 + n3))

job=0
for truth in "${TRUTHS[@]}"; do
for fc in "${FCS[@]}"; do
for rep in $(seq 1 "$REPLICAS"); do
    job=$((job + 1))
    seed_sim=$((BASE_SEED + 10 * job + 1))
    seed_box=$((BASE_SEED + 10 * job + 2))
    seed_bes=$((BASE_SEED + 10 * job + 3))
    label="fig2_${truth}_fc${fc}_rep${rep}_s${seed_sim}"
    log="$WORKDIR/logs/${label}.log"

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

    echo "[fig2] ($job) $label: sim (poles=$sim_poles fc=$fc)"
    "$BIN" \
        "$(printf -- '--exp_n_step_1 = get_number(n=%s)' "$n1")" \
        "$(printf -- '--exp_n_step_2 = get_number(n=%s)' "$n2")" \
        "$(printf -- '--exp_n_step_3 = get_number(n=%s)' "$n3")" \
        "$(printf -- '--exp_n_samp = get_number(n=%s)' "$N_SAMP")" \
        "$(printf -- '--agonist_2 = %s' "$ag2")" \
        "$(printf -- '--agonist_3 = %s' "$ag3")" \
        "$(printf -- '--truth_model = "scheme_CCO"')" \
        "$(printf -- '--truth_par_file = "%s"' "$PAR")" \
        "$(printf -- '--observations_template = "%s"' "$template")" \
        "$(printf -- '--sim_prefix = "%s"' "$label")" \
        "$(printf -- '--sim_seed = get_number(n=%s)' "$seed_sim")" \
        "$(printf -- '--sim_filter_n_poles = get_number(n=%s)' "$sim_poles")" \
        "$(printf -- '--sim_filter_cutoff = %s' "$fc")" \
        "$SIM_SCRIPT" 2>&1 | filter_warns | tee "$log"

    sim_csv="$(ls -t "${label}"_*_simulation.csv 2>/dev/null | head -1)"
    [ -n "$sim_csv" ] || { echo "[fig2] no simulation csv for $label" >&2; exit 1; }
    sim_csv="$(readlink -f "$sim_csv")"

    if [ "$PROT" = "stationary" ]; then
        masked="$WORKDIR/data/${label}_masked.csv"
        awk -F, -v a="$n1" -v b="$((n1 + n2))" \
            'NR==1{print; next} { if ($1+0>=a && $1+0<b) print $1",nan"; else print }' \
            "$sim_csv" > "$masked"
        obs="$masked"
    else
        obs="$sim_csv"
    fi

    echo "[fig2] ($job) $label: paired member evidence (fit filter ${FIT_NPOLES}p @ ${fc} Hz)"
    "$BIN" \
        "$(printf -- '--exp_n_step_1 = get_number(n=%s)' "$n1")" \
        "$(printf -- '--exp_n_step_2 = get_number(n=%s)' "$n2")" \
        "$(printf -- '--exp_n_step_3 = get_number(n=%s)' "$n3")" \
        "$(printf -- '--exp_n_samp = get_number(n=%s)' "$N_SAMP")" \
        "$(printf -- '--agonist_2 = %s' "$ag2")" \
        "$(printf -- '--agonist_3 = %s' "$ag3")" \
        "$(printf -- '--observations_file = "%s"' "$obs")" \
        "$(printf -- '--prior_file = "%s"' "$PRIOR")" \
        "$(printf -- '--fit_model = "scheme_CCO"')" \
        "$(printf -- '--idname_box = "%s_fit_box"' "$label")" \
        "$(printf -- '--idname_bessel = "%s_fit_bessel"' "$label")" \
        "$(printf -- '--num_scouts = get_number(n=%s)' "$SCOUTS")" \
        "$(printf -- '--beta_size = get_number(n=%s)' "$BETA_SIZE")" \
        "$(printf -- '--max_iter = get_number(n=%s)' "$MAX_ITER")" \
        "$(printf -- '--adapt_beta_every = get_number(n=%s)' "$ADAPT_EVERY")" \
        "$(printf -- '--max_values = get_number(n=%s)' "$MAX_VALUES")" \
        "$(printf -- '--seed_box = get_number(n=%s)' "$seed_box")" \
        "$(printf -- '--seed_bessel = get_number(n=%s)' "$seed_bes")" \
        "$(printf -- '--fit_filter_n_poles = get_number(n=%s)' "$FIT_NPOLES")" \
        "$(printf -- '--fit_filter_cutoff = %s' "$fc")" \
        "$EVI_SCRIPT" 2>&1 | filter_warns | tee -a "$log"

    echo "[fig2] ($job) $label done"
done
done
done
echo "[fig2] $job cells done, outputs under $WORKDIR"
