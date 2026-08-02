#!/bin/bash
# What is queued / running / done, decoded back into the sweep coordinates
# (algorithm, N_ch, noise, n_sim) that the job name cannot carry.
#
# WHY THIS EXISTS. The dispatchers name a job f4<fam>_<Nch>_<noise> (e.g.
# f4G_1E4_1E5), which pins the family and the two decades but NOT the algorithm
# — so one grid cell run under three algorithms shows up as three identical
# names in squeue. The obvious fix, reading the submit arguments back from
# `scontrol show job`, does NOT work here: this Slurm (23.11) prints
#   Command=/…/run_macroir.sh
# with the argv dropped. The full submit line does survive in the accounting
# DB (sacct SubmitLine, Slurm >= 21.08), and the dispatcher's --filepath
# injection carries every coordinate in one token:
#   figures/data/figure_3_G_nch_<Nch>_nsim_<nsim>_<algo>_noise_<label>
# so that is what gets parsed back into columns below. group_size, the interval
# grid and WORKDIR come from the same line.
#
# NOISE is the dispatcher LABEL, not what the binary sees: current_noise =
# label/1000 (and the paper's dimensionless S = 0.1 * label). Printed as the
# label so it matches the job name and the output filenames.
#
# OUT counts the CSVs already sitting in that job's WORKDIR for that exact cell.
# 0 on a PENDING job = still to run; >0 on a PENDING job = a previous launch
# already produced it, so the queued one is a redo.
#
# Usage (from anywhere):
#   projects/eLife_2025/ops/slurm/queue_status.sh            # queue now (PD + R)
#   projects/eLife_2025/ops/slurm/queue_status.sh -a         # every state since -S
#   projects/eLife_2025/ops/slurm/queue_status.sh -a -S 2026-07-01
#   projects/eLife_2025/ops/slurm/queue_status.sh -n         # skip the file count (no scratch I/O)
#   projects/eLife_2025/ops/slurm/queue_status.sh -u otheruser

set -eo pipefail

U="$USER"
SINCE="$(date -d '60 days ago' +%F 2>/dev/null || echo 2000-01-01)"
STATES="PENDING,RUNNING"
COUNT_FILES=1

while [ $# -gt 0 ]; do
    case "$1" in
        -u) U="$2"; shift 2 ;;
        -S) SINCE="$2"; shift 2 ;;
        -a) STATES=""; shift ;;
        -n) COUNT_FILES=0; shift ;;
        -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "queue_status: unknown option '$1' (see -h)" >&2; exit 1 ;;
    esac
done

command -v sacct >/dev/null || { echo "queue_status: sacct not found — run this on the cluster" >&2; exit 1; }

shopt -s nullglob

# Reasons live only in squeue (sacct has no Reason column), so pull them first.
declare -A REASON=()
if command -v squeue >/dev/null; then
    while IFS='|' read -r jid rsn; do REASON[$jid]="$rsn"; done \
        < <(squeue -u "$U" -h -o "%i|%R" 2>/dev/null)
fi

hdr() { printf "%-8s %-3s %-11s %-10s %-13s %-7s %-9s %-6s %-4s %-9s %s\n" "$@"; }

# Launch-level knobs are constant within one dispatch; tallied and printed once
# at the bottom instead of repeated in every row.
declare -A KNOBS=()
declare -A PEND=()
n_rows=0

# The tallies below are built inside the row loop, so the loop cannot run in a
# pipe subshell (its variables would die with it) — rows land in a temp file and
# get sorted after.
TMP="$(mktemp "${TMPDIR:-/tmp}/queue_status.XXXXXX")"
trap 'rm -f "$TMP"' EXIT

{
while IFS='|' read -r jid jname state subtime elapsed submit; do
    # CANCELLED comes through as "CANCELLED by <uid>"; keep the first word so the
    # row keeps its column count (the sort keys are positional).
    state="${state%% *}"
    day="${subtime%%T*}"
    case "$state" in
        PENDING)        st=PD ;;
        RUNNING)        st=R  ;;
        COMPLETED)      st=CD ;;
        FAILED)         st=F  ;;
        TIMEOUT)        st=TO ;;
        CANCELLED*)     st=CA ;;
        OUT_OF_MEMORY)  st=OOM ;;
        NODE_FAIL)      st=NF ;;
        *)              st="${state:0:3}" ;;
    esac

    # --- the filepath injection: everything hangs off it -------------------
    algo="?"; nch="?"; noise="?"; nsim="?"; fp=""
    if [[ $submit == *figures/data/* ]]; then
        fp="figures/data/${submit#*figures/data/}"
        fp="${fp%%\"*}"                       # stop at the closing DSL quote
        b="${fp##*/}"                         # figure_3_G_nch_10_nsim_256_macro_IR_noise_100
        if [[ $b == *_nch_*_nsim_*_noise_* ]]; then
            r="${b#*_nch_}";    nch="${r%%_nsim_*}"
            r="${r#*_nsim_}";   nsim="${r%%_*}"
            r="${r#*_}";        algo="${r%%_noise_*}"
            noise="${r#*_noise_}"
        fi
    fi

    # --- launch knobs + workdir from the same line ------------------------
    gs="-"; iv="-"; wd=""
    if [[ $submit == *'"group_size", labels= ['* ]]; then
        gs="${submit#*'"group_size", labels= ['}"; gs="${gs%%]*}"; gs="${gs//\"/}"
    fi
    if [[ $submit == *'"interval_in_tau", labels= ['* ]]; then
        iv="${submit#*'"interval_in_tau", labels= ['}"; iv="${iv%%]*}"; iv="${iv//\"/}"
    fi
    if [[ $submit == *WORKDIR=* ]]; then
        wd="${submit#*WORKDIR=}"; wd="${wd%%,*}"; wd="${wd%% *}"
    fi
    run="${wd##*/}"; run="${run:--}"

    # --- how much output that exact cell already has ----------------------
    # Glob with a trailing "_" so noise 100 does not swallow noise 1000.
    n_out="-"
    if [ "$COUNT_FILES" = 1 ] && [ -n "$wd" ] && [ -n "$fp" ]; then
        outs=("$wd/${fp}"_*)
        n_out="${#outs[@]}"
    fi

    KNOBS["gs=[$gs]  intervals=[$iv]"]=$(( ${KNOBS["gs=[$gs]  intervals=[$iv]"]:-0} + 1 ))
    [ "$st" = PD ] && PEND[$algo]=$(( ${PEND[$algo]:-0} + 1 ))
    n_rows=$(( n_rows + 1 ))

    printf "%-8s %-3s %-11s %-10s %-13s %-7s %-9s %-6s %-4s %-9s %s\n" \
        "$jid" "$st" "$elapsed" "$day" "$algo" "$nch" "$noise" "$nsim" "$n_out" "$run" \
        "${REASON[$jid]:-$jname}"
done < <(sacct -u "$U" -X --noheader -P -S "$SINCE" ${STATES:+-s "$STATES"} \
             -o JobID,JobName,State,Submit,Elapsed,SubmitLine)
} > "$TMP"

hdr JOBID ST TIME SUBMIT ALGO NCH NOISE NSIM OUT RUN NODE/REASON
# -b is load-bearing: without it a field carries its leading padding, so a short
# TIME (more pad spaces) sorts ahead of a long one and the ALGO key is ignored.
# LC_ALL=C so "_" is not collated away and the order is the same everywhere.
LC_ALL=C sort -b -k5,5 -k6,6n -k7,7n "$TMP"

echo
echo "$n_rows jobs   (NOISE = dispatcher label; current_noise = label/1000.  OUT = CSVs already written for that cell)"
if [ ${#PEND[@]} -gt 0 ]; then
    echo "pendientes por algoritmo:"
    for a in "${!PEND[@]}"; do printf "   %-14s %s\n" "$a" "${PEND[$a]}"; done | sort
fi
if [ ${#KNOBS[@]} -gt 0 ]; then
    echo "knobs por tanda:"
    for k in "${!KNOBS[@]}"; do printf "   %-4s jobs   %s\n" "${KNOBS[$k]}" "$k"; done | sort -rn
fi
