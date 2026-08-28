#!/usr/bin/env bash
# Renders the Figure 4 set, skipping whatever is already up to date.
#
# There are two layers of not-recomputing and they do different jobs:
#
#   THIS SCRIPT decides whether to start R at all. A notebook is rebuilt when its PDF is missing or
#   older than its own .Rmd, than the two shared files, or than the newest battery file on disk.
#   When nothing has moved it starts no process and takes under a second.
#
#   THE CACHE in figure_4_common.R decides whether to parse the CSVs. Six notebooks read the same
#   gigabyte to draw a few thousand numbers; `cached()` keeps the small filtered products in
#   figures/.cache and reuses them until an input actually changes. Its fingerprint carries the
#   name, size and mtime of every data file AND of figure_4_common.R and figure_4_layout.R, so
#   editing the code that builds a product invalidates it too.
#
# So a data-only change re-renders everything but parses once; a code-only change in one notebook
# re-renders that one and parses nothing; and nothing changing does nothing.
#
#   ./render_figure_4_set.sh              rebuild what is stale
#   ./render_figure_4_set.sh -n           say what would be rebuilt, do nothing
#   ./render_figure_4_set.sh -f           rebuild everything
#   ./render_figure_4_set.sh figure_4     rebuild just this one, if stale (-f to force)
#
# One notebook per R process, deliberately: the roster is set through globals that persist for the
# session (FIG4_ALGOS and friends), so two notebooks in one session would have the second inherit
# the first one's roster.

set -u
cd "$(dirname "$0")"

DEPS=(figure_4_common.R figure_4_layout.R)

# notebook : the outputs it writes, comma separated
JOBS=(
  "figure_4:Figure_4.pdf"
  "figure_4_supplement_standard_error:Figure_4_supplement_standard_error.pdf"
  "figure_4_S1:Figure_4_supplement_sample_corr.pdf"
  # 2026-08-27: the recursive ladder moved from slot 3 to slot 2 with the manuscript renumber, so
  # the notebook that draws it is now figure_4_S2 and it writes Figure_4_S2.pdf. The output names
  # on the neighbouring lines are pre-2026-08-12 and were already stale before that swap.
  "figure_4_S2:Figure_4_S2.pdf"
  "figure_4_S3:Figure_4_S3.pdf"
  "figure_4_supplement_lag_kappa:Figure_4_supplement_lag_kappa.pdf"
  "figure_4_supplement_se_kappa:Figure_4_supplement_se_kappa.pdf"
  "figure_4_supplement_coverage:Figure_4_supplement_coverage.pdf"
  "figure_4_supplement_qq_grid:Figure_4_supplement_qq_grid.pdf,Figure_4_supplement_qq_pooled.pdf"
  "figure_4_supplement_other_parameters:Figure_4_supplement_other_parameters_bias.pdf,Figure_4_supplement_other_parameters_distortion.pdf"
)

DRY=0; FORCE=0; ONLY=""
for a in "$@"; do
  case "$a" in
    -n|--dry-run) DRY=1 ;;
    -f|--force)   FORCE=1 ;;
    -*)           echo "unknown flag $a" >&2; exit 2 ;;
    *)            ONLY="${a%.Rmd}" ;;
  esac
done

mtime() { [ -e "$1" ] && stat -c %Y "$1" || echo 0; }

# One sentinel for the data: the newest battery file anywhere on the search path. Cheaper than
# resolving each notebook's own cell list here, and conservative in the right direction — a new run
# for one member rebuilds the set, which is what you want anyway when a column fills in.
newest_data=$(find ../data -name '*_battery_*.csv' -printf '%T@\n' 2>/dev/null |
              sort -rn | head -1 | cut -d. -f1)
: "${newest_data:=0}"

dep_max=$newest_data
for d in "${DEPS[@]}"; do
  m=$(mtime "$d"); [ "$m" -gt "$dep_max" ] && dep_max=$m
done

todo=(); skipped=0
for job in "${JOBS[@]}"; do
  nb="${job%%:*}"; outs="${job#*:}"
  [ -n "$ONLY" ] && [ "$nb" != "$ONLY" ] && continue
  [ -f "$nb.Rmd" ] || { echo "  missing   $nb.Rmd"; continue; }

  need=$dep_max
  m=$(mtime "$nb.Rmd"); [ "$m" -gt "$need" ] && need=$m

  oldest=""
  IFS=, read -ra arr <<< "$outs"
  for o in "${arr[@]}"; do
    m=$(mtime "$o")
    [ -z "$oldest" ] && oldest=$m
    [ "$m" -lt "$oldest" ] && oldest=$m
  done

  if [ "$FORCE" -eq 1 ] || [ "$oldest" -eq 0 ] || [ "$oldest" -lt "$need" ]; then
    todo+=("$nb")
    if [ "$oldest" -eq 0 ]; then why="no output yet"; else why="inputs are newer"; fi
    echo "  BUILD     $nb   ($why)"
  else
    skipped=$((skipped + 1))
    echo "  up to date $nb"
  fi
done

[ "$DRY" -eq 1 ] && { echo; echo "dry run: ${#todo[@]} to build, $skipped up to date"; exit 0; }
[ "${#todo[@]}" -eq 0 ] && { echo; echo "nothing to do ($skipped up to date)"; exit 0; }

echo
fail=0
for nb in "${todo[@]}"; do
  echo "=== $nb"
  if ! Rscript -e "rmarkdown::render('$nb.Rmd', quiet = TRUE)"; then
    echo "!!! $nb FAILED" >&2; fail=$((fail + 1))
  fi
done

echo
echo "built ${#todo[@]}, failed $fail, skipped $skipped"
[ "$fail" -eq 0 ] || exit 1
