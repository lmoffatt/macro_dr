#!/bin/bash
# figure_0 payload: parallelism characterization of the evidence machinery,
# runnable BOTH as a SLURM job (via dispatch_figure_0.sh) and directly on a
# local machine. Two modes:
#
#   MODE=full   the one-off characterization (the figure_0 baseline):
#               strong scaling over THREADS_LIST, then real packing tests
#               (P concurrent instances x TOTAL_CPUS/P threads), long slices.
#               Writes baseline.env next to tuning.env.
#   MODE=probe  the per-campaign preamble (~10 min): a few points around the
#               known optimum (from baseline.env if present), short slices.
#               Writes tuning.env and compares against baseline.env: a
#               deviation beyond TOL flags a changed node/build before the
#               campaign burns cluster hours.
#
# Method: every configuration runs the SAME reference job (the figure_1
# stationary CCO fit, adjust_beta=0) under `timeout SLICE`; iterations per
# second come from the iter/iter_time columns the binary already emits in
# <idname>_*__i_iter.csv (slope over the last 60% of iterations, so ladder
# warm-up and cache fill do not bias the rate). Killing at the timeout is
# harmless: only the csv is read.
#
# Timing-only caveat: RNG streams differ across thread counts, so numbers
# from different T are NOT comparable statistically, only in speed.
#
# Required env: BIN, WORKDIR, SIM_SCRIPT, EVI0_SCRIPT, PROJ.
# Optional: MODE (probe), TOTAL_CPUS (nproc), SLICE (180 full / 75 probe),
#   THREADS_LIST, PACK_LIST, NCH, TAU_MS, INTERVAL_IN_TAU, PRE_TAUS, GAP_TAUS,
#   MEAS_TAUS, SCOUTS, BETA_SIZE, ADAPT_EVERY, MAX_VALUES, SEED, TOL (0.15),
#   MACRODR_PROFILE + CLUSTER (sourced when set, i.e. under SLURM).

#SBATCH --job-name=fig0
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=1:30:00
#SBATCH --output=slurm-%j.out

set -eo pipefail

: "${BIN:?}" ; : "${WORKDIR:?}" ; : "${SIM_SCRIPT:?}" ; : "${EVI0_SCRIPT:?}" ; : "${PROJ:?}"

if [ -n "${MACRODR_PROFILE:-}" ]; then
    unset LOADEDMODULES _LMFILES_ 2>/dev/null || true
    set +e
    [ -f /etc/profile ] && source /etc/profile
    # shellcheck source=/dev/null
    source "$MACRODR_PROFILE"
    set -e
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

MODE="${MODE:-probe}"
TOTAL_CPUS="${TOTAL_CPUS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
TOL="${TOL:-0.15}"

NCH="${NCH:-1000}"
TAU_MS="${TAU_MS:-16}"
INTERVAL_IN_TAU="${INTERVAL_IN_TAU:-0.1}"
PRE_TAUS="${PRE_TAUS:-1}"
GAP_TAUS="${GAP_TAUS:-20}"
MEAS_TAUS="${MEAS_TAUS:-10}"
SCOUTS="${SCOUTS:-32}"
BETA_SIZE="${BETA_SIZE:-16}"
ADAPT_EVERY="${ADAPT_EVERY:-256}"
MAX_VALUES="${MAX_VALUES:-128}"
SEED="${SEED:-910001}"

N_SAMP=$(python3 -c "print(max(1, round(50*$TAU_MS*$INTERVAL_IN_TAU)))")
IPT=$(python3 -c "print(max(1, round(1/$INTERVAL_IN_TAU)))")
N1=$((PRE_TAUS * IPT)); N2=$((GAP_TAUS * IPT)); N3=$((MEAS_TAUS * IPT))

F0DIR="$WORKDIR/figure_0"
mkdir -p "$F0DIR"
cd "$F0DIR"

if [ "$MODE" = "full" ]; then
    SLICE="${SLICE:-180}"
    THREADS_LIST=(${THREADS_LIST:-1 2 4 8 16 32})
    PACK_LIST=(${PACK_LIST:-1 2 4 8})
else
    SLICE="${SLICE:-75}"
    # probe points around the known optimum (baseline.env), else around 8
    T_REF=8; P_REF=4
    if [ -f "$F0DIR/baseline.env" ]; then
        # shellcheck source=/dev/null
        source "$F0DIR/baseline.env"
        T_REF="${CPUS_RECOMMENDED:-8}"; P_REF="${JOBS_PER_NODE:-4}"
    fi
    t_lo=$((T_REF / 2)); [ "$t_lo" -lt 1 ] && t_lo=1
    t_hi=$((T_REF * 2)); [ "$t_hi" -gt "$TOTAL_CPUS" ] && t_hi=$TOTAL_CPUS
    THREADS_LIST=(${THREADS_LIST:-$t_lo $T_REF $t_hi})
    PACK_LIST=(${PACK_LIST:-$P_REF})
fi

# drop scaling points beyond the cores actually available (a local machine
# with 8 cores should not time T=16 or 32: oversubscription measures nothing)
_tl=()
for T in "${THREADS_LIST[@]}"; do
    [ "$T" -le "$TOTAL_CPUS" ] && _tl+=("$T")
done
THREADS_LIST=("${_tl[@]}")

echo "[fig0] mode=$MODE cpus=$TOTAL_CPUS slice=${SLICE}s bin=$BIN"
echo "[fig0] scaling T = ${THREADS_LIST[*]} ; packing P = ${PACK_LIST[*]}"

seg_args=()
seg_args+=("$(printf -- '--exp_n_step_1 = get_number(n=%s)' "$N1")")
seg_args+=("$(printf -- '--exp_n_step_2 = get_number(n=%s)' "$N2")")
seg_args+=("$(printf -- '--exp_n_step_3 = get_number(n=%s)' "$N3")")
seg_args+=("$(printf -- '--exp_n_samp = get_number(n=%s)' "$N_SAMP")")
seg_args+=("$(printf -- '--agonist_2 = 10.0')")
seg_args+=("$(printf -- '--agonist_3 = 10.0')")

# ---- reference recording: simulate once, reuse across every config ---------
OBS="$F0DIR/fig0_reference_masked.csv"
if [ ! -f "$OBS" ]; then
    # N-variant inputs (same generator as the figure_1 dispatchers)
    python3 - "$PROJ" "$F0DIR" "$NCH" << 'EOF'
import csv, math, sys, pathlib
proj, work, nch = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), float(sys.argv[3])
for name in ["scheme_CCO_par", "scheme_CCO_prior"]:
    src = proj / "data" / "models" / f"{name}.csv"
    dst = work / f"{name}_N{int(nch)}.csv"
    rows = list(csv.reader(open(src)))
    for r in rows[1:]:
        if r[2] == "Num_ch_mean":
            r[4] = repr(nch)
            r[5] = repr(math.log10(nch))
    csv.writer(open(dst, "w", newline="")).writerows(rows)
EOF
    total=$((N1 + N2 + N3))
    python3 - "$F0DIR/fig0_template.csv" "$total" << 'EOF'
import sys
path, n = sys.argv[1], int(sys.argv[2])
with open(path, "w") as f:
    f.write("i_step,patch_current\n")
    for i in range(n):
        f.write(f"{i},0\n")
EOF
    OMP_NUM_THREADS="$TOTAL_CPUS" "$BIN" \
        "${seg_args[@]}" \
        "$(printf -- '--truth_model = "scheme_CCO"')" \
        "$(printf -- '--truth_par_file = "%s"' "$F0DIR/scheme_CCO_par_N${NCH}.csv")" \
        "$(printf -- '--observations_template = "%s"' "$F0DIR/fig0_template.csv")" \
        "$(printf -- '--sim_prefix = "fig0_reference"')" \
        "$(printf -- '--sim_seed = get_number(n=%s)' "$SEED")" \
        "$SIM_SCRIPT" > "$F0DIR/sim.log" 2>&1
    sim_csv="$(ls -t fig0_reference_*_simulation.csv | head -1)"
    awk -F, -v a="$N1" -v b="$((N1 + N2))" \
        'NR==1{print; next} { if ($1+0>=a && $1+0<b) print $1",nan"; else print }' \
        "$sim_csv" > "$OBS"
fi
PRIOR="$F0DIR/scheme_CCO_prior_N${NCH}.csv"

# ---- role variants of the evidence script -----------------------------------
# adjust_beta is a bool in set_ThermoAlgorithm_dts and the DSL cannot pass a
# bool through an injected name (get_number returns size_t, an injected numeric
# literal compiles to double, and identifier arguments must match the parameter
# type exactly; only a literal IN argument position decodes to bool). So the
# canonical script carries the marker `adjust_beta_flag` and each role gets its
# own copy with the literal substituted.
EVI0_DISC="$F0DIR/figure_0_evidence_discover.macroir"
EVI0_TIME="$F0DIR/figure_0_evidence_timing.macroir"
sed 's/adjust_beta = adjust_beta_flag/adjust_beta = 1/' "$EVI0_SCRIPT" > "$EVI0_DISC"
sed 's/adjust_beta = adjust_beta_flag/adjust_beta = 0/' "$EVI0_SCRIPT" > "$EVI0_TIME"
if grep -v '^#' "$EVI0_DISC" | grep -q adjust_beta_flag; then
    echo "[fig0] ERROR: could not substitute the adjust_beta literal in $EVI0_SCRIPT" >&2
    echo "[fig0]        expected the marker line 'adjust_beta = adjust_beta_flag'" >&2
    exit 1
fi

# Drop the flood of benign range-canary warnings ([warn] check_g*_in_range,
# to_Probability*, and their torn fragments from concurrent stderr writes):
# measured 2026-09-05 at ~2300 lines/iteration under MCMC (267 MB in the
# 4-minute discovery slice alone). Real errors do not carry these markers
# and pass through. The || true keeps an all-filtered (empty) stream from
# tripping pipefail.
filter_warns() {
    grep -v -e '\[warn\]' -e 'check_g' -e 'to_Probability' -e 'to_Covariance' \
            -e 'warn=' -e 'err=' -e 'excursion' -e 'P_ij' -e 'gsqr' || true
}

# ---- one timed slice: run under timeout, read iterations/second ------------
# rate = slope of (iter, iter_time) over the last 60% of iterations.
read_rate() {  # $1 = idname; prints "<iters> <s_per_iter>" or "0 nan"
    local csv
    csv="$(ls -t "$1"_*__i_iter.csv 2>/dev/null | head -1)"
    [ -n "$csv" ] || { echo "0 nan"; return; }
    awk -F, 'NR>1 && !($1 in t) { t[$1]=$2; if ($1+0>m) m=$1+0 }
        END {
            lo = 0.4*m; i1 = -1; i2 = -1
            for (i in t) { iv = i+0
                if (iv >= lo) {
                    if (i1 < 0 || iv < i1) i1 = iv
                    if (iv > i2) i2 = iv } }
            if (i2 > i1) printf "%d %.6g\n", m, (t[i2]-t[i1])/(i2-i1)
            else printf "%d nan\n", m }' "$csv"
}

run_slice() {  # $1 = tag, $2 = threads, $3 = beta_size, $4 = adjust_flag, $5 = seconds
    local tag=$1 threads=$2 nbeta=$3 adj=$4 secs=$5
    local script="$EVI0_TIME"
    [ "$adj" = "1" ] && script="$EVI0_DISC"
    OMP_NUM_THREADS="$threads" timeout --signal=TERM "${secs}s" "$BIN" \
        "${seg_args[@]}" \
        "$(printf -- '--observations_file = "%s"' "$OBS")" \
        "$(printf -- '--prior_file = "%s"' "$PRIOR")" \
        "$(printf -- '--fit_model = "scheme_CCO"')" \
        "$(printf -- '--run_idname = "%s"' "$tag")" \
        "$(printf -- '--num_scouts = get_number(n=%s)' "$SCOUTS")" \
        "$(printf -- '--beta_size = get_number(n=%s)' "$nbeta")" \
        "$(printf -- '--max_iter = get_number(n=999999)')" \
        "$(printf -- '--adapt_beta_every = get_number(n=%s)' "$ADAPT_EVERY")" \
        "$(printf -- '--max_values = get_number(n=%s)' "$MAX_VALUES")" \
        "$(printf -- '--mcmc_seed = get_number(n=%s)' "$SEED")" \
        "$script" 2>&1 | filter_warns > "$F0DIR/${tag}.log" || true
}

read_nbeta() {  # $1 = idname; prints "<first_num_beta> <last_num_beta>"
    local csv
    csv="$(ls -t "$1"_*__i_iter.csv 2>/dev/null | head -1)"
    [ -n "$csv" ] || { echo "0 0"; return; }
    awk -F, 'NR==2{f=$4} NR>1{l=$4} END{print f+0, l+0}' "$csv"
}

table="$F0DIR/figure_0_table.txt"
: > "$table"
echo "config threads instances iters s_per_iter node_iters_per_s" | tee -a "$table"

# ---- stage 0: ladder discovery ----------------------------------------------
# The per-iteration load scales with n_beta x n_walkers and the parallel
# efficiency depends on that grain, so timing at a nominal beta_size measures
# the wrong machine. Run the CAMPAIGN configuration (adjust_beta = 1, starting
# at BETA_SIZE) until the timeout, read where the ladder settled, and time all
# the configurations below AT that count.
if [ -z "${DISC_SLICE:-}" ]; then
    if [ "$MODE" = "full" ]; then DISC_SLICE=240; else DISC_SLICE=120; fi
fi
run_slice "f0_discover" "$TOTAL_CPUS" "$BETA_SIZE" 1 "$DISC_SLICE"
# fail fast: no iter csv at all means the run died at startup (DSL/compile
# error, bad path...); every later slice would fail identically, so stop here
# with the actual error instead of cascading nan rates.
if ! ls f0_discover_*__i_iter.csv > /dev/null 2>&1; then
    echo "[fig0] ERROR: discovery produced no __i_iter.csv; the run died at startup." >&2
    echo "[fig0] last lines of $F0DIR/f0_discover.log:" >&2
    tail -n 30 "$F0DIR/f0_discover.log" >&2
    exit 1
fi
read -r nbeta_first nbeta_last <<< "$(read_nbeta f0_discover)"
BETA_EFF="$nbeta_last"
[ "$BETA_EFF" -ge 2 ] 2>/dev/null || BETA_EFF="$BETA_SIZE"
echo "[fig0] ladder discovery: started $BETA_SIZE, first seen $nbeta_first, settled $BETA_EFF"
if [ "$nbeta_first" != "$nbeta_last" ]; then
    echo "[fig0] NOTE: ladder still moving at the discovery timeout ($nbeta_first -> $nbeta_last);"
    echo "[fig0]       a longer DISC_SLICE would pin it better. Timing uses $BETA_EFF."
fi

declare -A RATE
# ---- strong scaling (at the settled ladder size) ----------------------------
for T in "${THREADS_LIST[@]}"; do
    tag="f0_scale_T${T}"
    run_slice "$tag" "$T" "$BETA_EFF" 0 "$SLICE"
    read -r iters spi <<< "$(read_rate "$tag")"
    nps=$(python3 -c "print('nan' if '$spi'=='nan' else round(1/$spi, 4))")
    RATE[$T]="$spi"
    echo "scale $T 1 $iters $spi $nps" | tee -a "$table"
done

# ---- packing: P concurrent instances x TOTAL_CPUS/P threads ----------------
best_nps="0"; best_cpus=""; best_jobs=""
for P in "${PACK_LIST[@]}"; do
    T=$((TOTAL_CPUS / P)); [ "$T" -lt 1 ] && continue
    pids=()
    for j in $(seq 1 "$P"); do
        run_slice "f0_pack_P${P}_j${j}" "$T" "$BETA_EFF" 0 "$SLICE" &
        pids+=($!)
    done
    wait "${pids[@]}"
    agg=0
    for j in $(seq 1 "$P"); do
        read -r iters spi <<< "$(read_rate "f0_pack_P${P}_j${j}")"
        agg=$(python3 -c "print($agg + (0 if '$spi'=='nan' else 1/$spi))")
    done
    agg=$(python3 -c "print(round($agg, 4))")
    echo "pack $T $P - - $agg" | tee -a "$table"
    if python3 -c "exit(0 if $agg > $best_nps else 1)"; then
        best_nps="$agg"; best_cpus="$T"; best_jobs="$P"
    fi
done

# fall back to the best scaling point if no packing row won
if [ -z "$best_cpus" ]; then
    for T in "${THREADS_LIST[@]}"; do
        spi="${RATE[$T]}"
        [ "$spi" = "nan" ] && continue
        nps=$(python3 -c "print(1/$spi)")
        if python3 -c "exit(0 if $nps > $best_nps else 1)"; then
            best_nps="$nps"; best_cpus="$T"; best_jobs=1
        fi
    done
fi
# if nothing produced a finite rate, stop with the first evidence log instead
# of interpolating empty variables into python (that was a SyntaxError crash)
if [ -z "$best_cpus" ] || [ -z "$best_jobs" ]; then
    echo "[fig0] ERROR: no configuration produced a finite rate." >&2
    first_log="$(ls -t "$F0DIR"/f0_*.log 2>/dev/null | tail -1)"
    if [ -n "$first_log" ]; then
        echo "[fig0] last lines of $first_log:" >&2
        tail -n 30 "$first_log" >&2
    fi
    exit 1
fi
spi_best=$(python3 -c "print(round($best_jobs/$best_nps, 6))")

# ---- full mode extra: work-vs-ladder linearity ------------------------------
# Two grain points at max threads confirm that s/iter scales ~linearly with
# n_beta, which is what lets a future campaign at another ladder size reuse
# this baseline by scaling instead of re-running full.
if [ "$MODE" = "full" ]; then
    for NB in $((BETA_EFF / 2)) $((BETA_EFF * 2)); do
        [ "$NB" -ge 2 ] || continue
        tag="f0_grain_NB${NB}"
        run_slice "$tag" "$TOTAL_CPUS" "$NB" 0 "$SLICE"
        read -r iters spi <<< "$(read_rate "$tag")"
        echo "grain $TOTAL_CPUS 1 $iters $spi (n_beta=$NB)" | tee -a "$table"
    done
fi

out="$F0DIR/tuning.env"
{
    echo "# figure_0 $MODE, $(date -Is), host $(hostname), commit $("$BIN" --commit)"
    echo "CPUS_RECOMMENDED=$best_cpus"
    echo "JOBS_PER_NODE=$best_jobs"
    echo "SECONDS_PER_ITER=$spi_best"
    echo "NODE_ITERS_PER_S=$best_nps"
    echo "BETA_EFF=$BETA_EFF"
} > "$out"
echo "[fig0] wrote $out: CPUS=$best_cpus JOBS_PER_NODE=$best_jobs node_iters/s=$best_nps"

if [ "$MODE" = "full" ]; then
    cp "$out" "$F0DIR/baseline.env"
    echo "[fig0] baseline recorded: $F0DIR/baseline.env"
elif [ -f "$F0DIR/baseline.env" ]; then
    base_nps="$(grep '^NODE_ITERS_PER_S=' "$F0DIR/baseline.env" | cut -d= -f2)"
    if [ -n "$base_nps" ] && python3 -c "exit(0 if abs($best_nps-$base_nps)/max($base_nps,1e-12) > $TOL else 1)"; then
        echo "[fig0] WARNING: node throughput $best_nps deviates >${TOL} from baseline $base_nps"
        echo "[fig0]          different node/build or a performance regression;"
        echo "[fig0]          consider rerunning MODE=full before the campaign."
    else
        echo "[fig0] within ${TOL} of baseline ($base_nps): good to dispatch."
    fi
fi
