#!/bin/bash
# Local dispatcher for Figure 1: the evidence confusion matrix
#   {truth: CCO, COC-twin} x {protocol: episodic, stationary} x REPLICAS,
# each replica fitted with BOTH models (paired Delta logZ).
#
# Per job it runs THREE stages, sequentially:
#   1. figure_1_sim.macroir       simulate the replica recording
#   2. (stationary only)          nan-mask the equilibration segment in the csv
#      The nan samples are "no measurement": the likelihood propagates without
#      a Bayes update and contributes zero logL (qmodel isnan gating), so the
#      chain equilibrates INSIDE the scored record and the measured stretch
#      starts at true stationarity.
#   3. figure_1_evidence.macroir  paired evidence, both fits on the same data
#
# Same injection contract as projects/eLife_2025/ops/local/dispatch_figure_3_local.sh:
# names are passed as "--name = expression" argv entries before the script, and
# the .macroir files must not define them.
#
# Design defaults (override via env):
#   NCH=1000 channels; TAU_MS=100 (the recording's slow relaxation is ~6.2
#   nominal tau, 2026-09-05, so the segments are laid out on the slow scale);
#   INTERVAL_IN_TAU=0.1 (measurement step = 0.1 tau);
#   episodic  = PRE_TAUS 1 | PULSE_TAUS 5  (10 uM) | POST_TAUS 5 (0 uM)
#   stationary= PRE_TAUS 1 | GAP_TAUS  20  (10 uM, nan-masked) | MEAS_TAUS 10 (10 uM)
#   REPLICAS=10; MAX_ITER=30000 (~300 full-ladder score events at MAX_VALUES=128);
#   ADAPT_EVERY=256 and ADAPT_T0=3000 (phase-1 adaptation); BETA_SIZE=4, grown
#   by adjust_beta during phase 1; SCOUTS=32; and the drift-and-hold ladder
#   schedule PHASE1_END=5000, DRIFT=300, HOLD=3000, HOLD_BURNIN=100,
#   N_CYCLES=4, CYCLE_GAIN=0.3 (set_Ladder_schedule; theory in
#   legacy/parallel_tempering.h). Regression check of the schedule: N_CYCLES=0
#   PHASE1_END=$MAX_ITER HOLD_BURNIN=0 reproduces the pre-schedule run.
# Ladder criterion (injected into set_ThermoAlgorithm_dts / set_Ladder_schedule):
# EQUALIZER (deltaBeta_deltaL_vfm | Acceptance_vfm | Acceptance_fixed_vfm),
# DESIRED_ACC (target of the _fixed one), ADAPT_BETA_MIN (0 pins beta_min, 1
# adapts it). LADDER_COMBOS="eq/acc/min/tag ..." sweeps several criteria on
# the SAME recordings (seeds are per cell, not per submission; the tag goes
# into the label), e.g. the six-way sweep of 2026-09-07:
#   LADDER_COMBOS="deltaBeta_deltaL_vfm/0.25/0/mu-pin deltaBeta_deltaL_vfm/0.25/1/mu-free \
#                  Acceptance_vfm/0.25/0/acc-pin Acceptance_vfm/0.25/1/acc-free \
#                  Acceptance_fixed_vfm/0.234/0/fix-pin Acceptance_fixed_vfm/0.234/1/fix-free" \
#   TRUTHS=CCO PROTOCOLS=episodic REPLICAS=2 ...
#
# Prereq: a local build (the user compiles): build/gcc-release/macrodr_cli.
# Usage: projects/macroir_next/ops/local/dispatch_figure_1_local.sh
# Quick smoke: REPLICAS=1 TRUTHS="CCO" PROTOCOLS="stationary" MAX_ITER=200 ...

set -eo pipefail

HERE="$(dirname "$(readlink -f "$0")")"             # .../macroir_next/ops/local
BASE="$(readlink -f "$HERE/../../../..")"           # repo base (macro_dr)
PROJ="$(readlink -f "$HERE/../..")"                 # .../projects/macroir_next

SIM_SCRIPT="$(readlink -f "$HERE/figure_1_sim.macroir")"
EVI_SCRIPT="$(readlink -f "$HERE/figure_1_evidence.macroir")"

BIN="${BIN:-$(readlink -f "$BASE/build/gcc-release/macrodr_cli")}"
[ -x "$BIN" ] || {
    echo "[fig1] binary not found: $BIN" >&2
    echo "       build first:  cmake --build --preset gcc-release" >&2
    exit 1
}

# ---- knobs -----------------------------------------------------------------
TRUTHS=(${TRUTHS:-CCO COC})
PROTOCOLS=(${PROTOCOLS:-stationary episodic})
REPLICAS="${REPLICAS:-10}"

NCH="${NCH:-1000}"
TAU_MS="${TAU_MS:-100}"
INTERVAL_IN_TAU="${INTERVAL_IN_TAU:-0.1}"
PRE_TAUS="${PRE_TAUS:-1}"
PULSE_TAUS="${PULSE_TAUS:-5}"       # episodic 10 uM segment
POST_TAUS="${POST_TAUS:-5}"         # episodic 0 uM tail
GAP_TAUS="${GAP_TAUS:-20}"          # stationary unmeasured equilibration (10 uM)
MEAS_TAUS="${MEAS_TAUS:-10}"        # stationary measured stretch (10 uM)

SCOUTS="${SCOUTS:-32}"
# resolved after the tuning hook (BETA_EFF from figure_0 can supply it)
BETA_SIZE="${BETA_SIZE:-}"
MAX_ITER="${MAX_ITER:-30000}"
ADAPT_EVERY="${ADAPT_EVERY:-256}"
ADAPT_T0="${ADAPT_T0:-3000}"
PHASE1_END="${PHASE1_END:-5000}"
DRIFT="${DRIFT:-300}"
HOLD="${HOLD:-3000}"
HOLD_BURNIN="${HOLD_BURNIN:-100}"
N_CYCLES="${N_CYCLES:-4}"
CYCLE_GAIN="${CYCLE_GAIN:-0.3}"
EQUALIZER="${EQUALIZER:-deltaBeta_deltaL_vfm}"
DESIRED_ACC="${DESIRED_ACC:-0.25}"
ADAPT_BETA_MIN="${ADAPT_BETA_MIN:-0}"
LADDER_COMBOS=(${LADDER_COMBOS:-"${EQUALIZER}/${DESIRED_ACC}/${ADAPT_BETA_MIN}/"})
MAX_VALUES="${MAX_VALUES:-128}"
BASE_SEED="${BASE_SEED:-910000}"    # nonzero; 0 would mean random_device

# steps (measurement intervals) per tau, and sub-samples per step at 50 kHz:
# step duration = INTERVAL_IN_TAU * TAU_MS, n_samp = 50 samples/ms * step_ms.
N_SAMP=$(python3 -c "print(max(1, round(50*$TAU_MS*$INTERVAL_IN_TAU)))")
IPT=$(python3 -c "print(max(1, round(1/$INTERVAL_IN_TAU)))")   # intervals per tau

# ---- per-commit output isolation, under runs/ (hygiene rule) ---------------
if ! commit="$("$BIN" --commit)"; then
    echo "[fig1] could not query commit hash: '$BIN --commit' failed" >&2
    exit 1
fi
[ -n "$commit" ] || { echo "[fig1] '$BIN --commit' returned empty" >&2; exit 1; }
WORKDIR="${WORKDIR:-$PROJ/runs/$commit}"
mkdir -p "$WORKDIR/logs" "$WORKDIR/data"
WORKDIR="$(readlink -f "$WORKDIR")"

# figure_0 tuning: if the local probe ran for this commit, take its
# recommended thread count as the default (explicit THREADS env still wins).
TUNE_FILE="${TUNE_FILE:-$WORKDIR/figure_0/tuning.env}"
if [ -f "$TUNE_FILE" ]; then
    # shellcheck source=/dev/null
    source "$TUNE_FILE"
    THREADS="${THREADS:-${CPUS_RECOMMENDED:-$(nproc)}}"
    # the ladder starts SMALL and adjust_beta grows it in phase 1; BETA_EFF is ignored
    echo "[fig1] tuning.env: THREADS=$THREADS beta_size=${BETA_SIZE:-4} s/iter=${SECONDS_PER_ITER:-?}"
fi
BETA_SIZE="${BETA_SIZE:-4}"

export OMP_NUM_THREADS="${THREADS:-$(nproc)}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${BLAS_THREADS:-1}"
export BLIS_NUM_THREADS="${BLAS_THREADS:-1}"

# ---- N-channel variants of the model inputs --------------------------------
# The seeded par/prior csvs carry Num_ch_mean=5000; Figure 1 wants NCH.
# Generate variants once per NCH into WORKDIR/data (value + log10 column).
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

cd "$WORKDIR"

# Drop the flood of benign range-canary warnings ([warn] check_g*_in_range,
# to_Probability*, and their torn fragments from concurrent stderr writes):
# measured 2026-09-05 at ~2300 lines/iteration under MCMC, 267 MB per
# 4-minute slice — it bloats logs AND stalls the terminal through tee. Real
# errors do not carry these markers and pass through. The || true keeps an
# all-filtered (empty) stream from tripping pipefail.
filter_warns() {
    grep -v -e '\[warn\]' -e 'check_g' -e 'to_Probability' -e 'to_Covariance' \
            -e 'warn=' -e 'err=' -e 'excursion' -e 'P_ij' -e 'gsqr' || true
}

echo "[fig1] base=$BASE  bin=$BIN  commit=$commit"
echo "[fig1] cwd=$WORKDIR  threads=$OMP_NUM_THREADS"
echo "[fig1] tau=${TAU_MS}ms  interval=${INTERVAL_IN_TAU}tau  n_samp=$N_SAMP  intervals/tau=$IPT  NCH=$NCH"

job=0
cell=0
for truth in "${TRUTHS[@]}"; do
for prot in "${PROTOCOLS[@]}"; do
for rep in $(seq 1 "$REPLICAS"); do
    # Seeds are deterministic from BASE_SEED and the (truth, protocol,
    # replica) CELL index, and the sim seed goes INTO the label (hence into
    # every output filename), so a local run and a dirac run with the same
    # BASE_SEED pair up file by file for verification, and every ladder combo
    # below fits the SAME recording. Byte-identical results additionally
    # require the same THREADS on both machines (one RNG stream per OMP
    # thread) and a numerically identical BLAS; otherwise expect float-level
    # differences.
    cell=$((cell + 1))
    seed_sim=$((BASE_SEED + 10 * cell + 1))
    seed_cco=$((BASE_SEED + 10 * cell + 2))
    seed_coc=$((BASE_SEED + 10 * cell + 3))
for combo in "${LADDER_COMBOS[@]}"; do
    IFS='/' read -r eq acc bmin tag <<< "$combo"
    job=$((job + 1))
    label="fig1_${tag:+${tag}_}${truth}_${prot}_rep${rep}_s${seed_sim}"
    log="$WORKDIR/logs/${label}.log"

    case "$truth" in
        CCO) truth_model="scheme_CCO"; truth_par="$PAR_CCO" ;;
        COC) truth_model="scheme_COC"; truth_par="$PAR_COC" ;;
        *) echo "[fig1] unknown truth '$truth' (want CCO or COC)" >&2; exit 1 ;;
    esac

    # segments in measurement intervals; agonist_3 tells the protocols apart
    n1=$((PRE_TAUS * IPT))
    case "$prot" in
        episodic)   n2=$((PULSE_TAUS * IPT)); n3=$((POST_TAUS * IPT)); ag2=10.0; ag3=0.0 ;;
        stationary) n2=$((GAP_TAUS  * IPT)); n3=$((MEAS_TAUS * IPT)); ag2=10.0; ag3=10.0 ;;
        *) echo "[fig1] unknown protocol '$prot'" >&2; exit 1 ;;
    esac
    total=$((n1 + n2 + n3))

    # observations template: one row per measurement interval, zeros
    template="$WORKDIR/data/${label}_template.csv"
    python3 - "$template" "$total" << 'EOF'
import sys
path, n = sys.argv[1], int(sys.argv[2])
with open(path, "w") as f:
    f.write("i_step,patch_current\n")
    for i in range(n):
        f.write(f"{i},0\n")
EOF

    seg1_arg=$(printf -- '--exp_n_step_1 = get_number(n=%s)' "$n1")
    seg2_arg=$(printf -- '--exp_n_step_2 = get_number(n=%s)' "$n2")
    seg3_arg=$(printf -- '--exp_n_step_3 = get_number(n=%s)' "$n3")
    nsamp_arg=$(printf -- '--exp_n_samp = get_number(n=%s)' "$N_SAMP")
    ag2_arg=$(printf -- '--agonist_2 = %s' "$ag2")
    ag3_arg=$(printf -- '--agonist_3 = %s' "$ag3")

    echo "[fig1] ($job) $label: sim (n1=$n1 n2=$n2 n3=$n3 n_samp=$N_SAMP) -> $log"
    "$BIN" \
        "$seg1_arg" "$seg2_arg" "$seg3_arg" "$nsamp_arg" "$ag2_arg" "$ag3_arg" \
        "$(printf -- '--truth_model = "%s"' "$truth_model")" \
        "$(printf -- '--truth_par_file = "%s"' "$truth_par")" \
        "$(printf -- '--observations_template = "%s"' "$template")" \
        "$(printf -- '--sim_prefix = "%s"' "$label")" \
        "$(printf -- '--sim_seed = get_number(n=%s)' "$seed_sim")" \
        "$SIM_SCRIPT" 2>&1 | filter_warns | tee "$log"

    # the simulator names its output <prefix>_<model>_<time>_<seed>_simulation.csv
    sim_csv="$(ls -t "${label}"_*_simulation.csv 2>/dev/null | head -1)"
    [ -n "$sim_csv" ] || { echo "[fig1] no simulation csv for $label" >&2; exit 1; }
    sim_csv="$(readlink -f "$sim_csv")"

    if [ "$prot" = "stationary" ]; then
        # nan-mask the equilibration segment. awk on purpose (not python):
        # the dirac payload uses the SAME awk line, so the masked files are
        # byte-comparable between local and cluster runs.
        masked="$WORKDIR/data/${label}_masked.csv"
        awk -F, -v a="$n1" -v b="$((n1 + n2))" \
            'NR==1{print; next} { if ($1+0>=a && $1+0<b) print $1",nan"; else print }' \
            "$sim_csv" > "$masked"
        obs="$masked"
    else
        obs="$sim_csv"
    fi

    echo "[fig1] ($job) $label: paired evidence on $(basename "$obs")"
    "$BIN" \
        "$seg1_arg" "$seg2_arg" "$seg3_arg" "$nsamp_arg" "$ag2_arg" "$ag3_arg" \
        "$(printf -- '--observations_file = "%s"' "$obs")" \
        "$(printf -- '--prior_cco_file = "%s"' "$PRIOR_CCO")" \
        "$(printf -- '--prior_coc_file = "%s"' "$PRIOR_COC")" \
        "$(printf -- '--idname_cco = "%s_fit_CCO"' "$label")" \
        "$(printf -- '--idname_coc = "%s_fit_COC"' "$label")" \
        "$(printf -- '--num_scouts = get_number(n=%s)' "$SCOUTS")" \
        "$(printf -- '--beta_size = get_number(n=%s)' "$BETA_SIZE")" \
        "$(printf -- '--max_iter = get_number(n=%s)' "$MAX_ITER")" \
        "$(printf -- '--adapt_beta_every = get_number(n=%s)' "$ADAPT_EVERY")" \
        "$(printf -- '--adapt_beta_t0 = %s' "$ADAPT_T0")" \
        "$(printf -- '--phase1_end = get_number(n=%s)' "$PHASE1_END")" \
        "$(printf -- '--drift_iters = get_number(n=%s)' "$DRIFT")" \
        "$(printf -- '--hold_iters = get_number(n=%s)' "$HOLD")" \
        "$(printf -- '--hold_burnin = get_number(n=%s)' "$HOLD_BURNIN")" \
        "$(printf -- '--n_cycles = get_number(n=%s)' "$N_CYCLES")" \
        "$(printf -- '--cycle_gain = %s' "$CYCLE_GAIN")" \
        "$(printf -- '--adapt_beta_equalizer = "%s"' "$eq")" \
        "$(printf -- '--desired_acceptance = %s' "$acc")" \
        "$(printf -- '--adapt_beta_min = get_number(n=%s)' "$bmin")" \
        "$(printf -- '--max_values = get_number(n=%s)' "$MAX_VALUES")" \
        "$(printf -- '--seed_cco = get_number(n=%s)' "$seed_cco")" \
        "$(printf -- '--seed_coc = get_number(n=%s)' "$seed_coc")" \
        "$EVI_SCRIPT" 2>&1 | filter_warns | tee -a "$log"
done
done
done
done
echo "[fig1] done: $job jobs in $WORKDIR"
