#!/usr/bin/env bash
# check.sh — the definition of done for the MacroIR paper (01_writing_plan.md, section 0).
#
# Derives the state of the manuscript from the repo. Nothing here is declared by hand;
# there is no ledger to keep in sync. Red is the normal state early on: read the FAILs
# as the task list.
#
#   ./check.sh            summary
#   ./check.sh -v         list every offending line
#
# Exit 0 only when the paper is submittable (all eight checks green).

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MS="$HERE/docs/manuscript-drafts"
TEX="$MS/elife_paper.tex"
BIB="$MS/biblio.bib"
SECTIONS="$MS/sections"
VERBOSE=0
[ "${1:-}" = "-v" ] && VERBOSE=1

PASS=0; FAIL=0; WARN=0
red()   { printf '\033[31mFAIL\033[0m  %s\n' "$1"; FAIL=$((FAIL+1)); }
green() { printf '\033[32mPASS\033[0m  %s\n' "$1"; PASS=$((PASS+1)); }
warn()  { printf '\033[33mWARN\033[0m  %s\n' "$1"; WARN=$((WARN+1)); }
detail(){ [ "$VERBOSE" = 1 ] && printf '        %s\n' "$1"; return 0; }

# Body = everything the author writes: the section files if T-1 has split them,
# otherwise the part of elife_paper.tex after \begin{document}. The preamble is exempt
# (\documentclass[9pt] is not a claim).
BODY_FILES=()
if [ -d "$SECTIONS" ] && compgen -G "$SECTIONS/*.tex" >/dev/null; then
  while IFS= read -r f; do BODY_FILES+=("$f"); done < <(find "$SECTIONS" -name '*.tex' | sort)
fi
BODY_TMP="$(mktemp)"; trap 'rm -f "$BODY_TMP"' EXIT
if [ ${#BODY_FILES[@]} -gt 0 ]; then
  for f in "${BODY_FILES[@]}"; do sed "s|^|${f}:|" "$f"; done > "$BODY_TMP"
else
  awk -v f="$TEX" 'p{print f":"FNR":"$0} /\\begin\{document\}/{p=1}' "$TEX" > "$BODY_TMP"
fi

echo "=== MacroIR paper — state of the manuscript ==="
echo

# --- 1. it compiles -----------------------------------------------------------------
if command -v latexmk >/dev/null 2>&1; then
  OUT="$(mktemp -d)"; trap 'rm -f "$BODY_TMP"; rm -rf "$OUT"' EXIT
  if (cd "$MS" && latexmk -pdf -outdir="$OUT" -interaction=nonstopmode -halt-on-error elife_paper.tex >/tmp/latexmk.log 2>&1); then
    green "1. compiles (latexmk -pdf)"
  else
    red "1. compiles (latexmk -pdf) — see /tmp/latexmk.log"
    detail "$(grep -m5 '^!' /tmp/latexmk.log || true)"
  fi
else
  warn "1. compiles — latexmk not installed"
fi

# --- 2. no holes --------------------------------------------------------------------
N_TODO=$(grep -c '\\todo{' "$BODY_TMP" || true)
if [ "$N_TODO" -eq 0 ]; then green "2. no \\todo{} remain"; else
  red "2. $N_TODO \\todo{} remain"
  while IFS= read -r l; do detail "$l"; done < <(grep -n '\\todo{' "$BODY_TMP" | head -20)
fi

# --- 3. every section written -------------------------------------------------------
# A section is unwritten if it holds only comments and blank lines (the current skeleton).
UNWRITTEN=0
while IFS= read -r sec; do
  name="${sec%%|*}"; n="${sec##*|}"
  if [ "$n" -eq 0 ]; then UNWRITTEN=$((UNWRITTEN+1)); detail "unwritten: $name"; fi
done < <(awk -F: '
  /\\(sub)?section\*?\{/ { if (cur!="") print cur"|"n; match($0,/\{[^}]*\}/); cur=substr($0,RSTART+1,RLENGTH-2); n=0; next }
  cur!="" { line=$0; sub(/^[^:]*:[0-9]*:/,"",line); gsub(/^[ \t]+/,"",line);
            if (line!="" && substr(line,1,1)!="%") n++ }
  END { if (cur!="") print cur"|"n }' "$BODY_TMP")
STALE=$(grep -cE '\b(TBD|XXX|FIXME)\b' "$BODY_TMP" || true)
if [ "$UNWRITTEN" -eq 0 ] && [ "$STALE" -eq 0 ]; then green "3. every section has prose"; else
  red "3. $UNWRITTEN section(s) hold only fill-hints; $STALE TBD/XXX/FIXME marker(s)"
fi

# --- 4. LINT-SRC: no number without a source ----------------------------------------
# A numeral in *prose* must carry "% src: <file>". Math is exempt (a formula is a derivation,
# not a measurement), as are cross-references, and section/figure numbering.
UNSOURCED="$(perl -ne '
  ($f,$rest) = /^([^:]*:\d*):(.*)$/ ? ($1,$2) : ($_,$_);
  $t = $rest;
  next if $t =~ /^\s*%/;                       # a whole-line comment is not prose
  $t =~ s/%.*$//;                              # strip trailing comment (incl. the % src: itself)
  $t =~ s/\$[^\$]*\$//g;                       # inline math
  $t =~ s/\\\[.*?\\\]//g;                      # display math
  $t =~ s/\\(ref|cite|label|includegraphics|input|eqref|autoref)\s*(\[[^\]]*\])?\{[^}]*\}//g;
  $t =~ s/\\(sub)?section\*?\{[^}]*\}//g;
  $t =~ s/(Figure|Fig\.?|Table|Section|Equation|Eq\.?)~?\s*\d+//gi;   # numbering, not claims
  next unless $t =~ /\d/;
  print "$f: $rest\n" unless $rest =~ /%\s*src:/;
' "$BODY_TMP")"
N_UNSRC=$(printf '%s' "$UNSOURCED" | grep -c . || true)
if [ "$N_UNSRC" -eq 0 ]; then green "4. LINT-SRC: every number carries % src:"; else
  red "4. LINT-SRC: $N_UNSRC line(s) carry a number with no % src:"
  while IFS= read -r l; do detail "$l"; done < <(printf '%s\n' "$UNSOURCED" | head -20)
fi

# --- 5. bibliography ----------------------------------------------------------------
N_CITE=$(grep -o '\\cite[tp]\?\s*{[^}]*}' "$BODY_TMP" | sed 's/.*{//;s/}//' | tr ',' '\n' | sed 's/ //g' | grep -c . || true)
if [ "$N_CITE" -eq 0 ]; then
  red "5. bibliography: zero \\cite in the manuscript (biblio.bib is still the P2X2 paper's references)"
else
  MISSING=0; UNCITED=0
  CITED="$(grep -o '\\cite[tp]\?\s*{[^}]*}' "$BODY_TMP" | sed 's/.*{//;s/}//' | tr ',' '\n' | sed 's/ //g' | sort -u | grep .)"
  KEYS="$(grep -o '^@[A-Za-z]*{[^,]*,' "$BIB" | sed 's/.*{//;s/,//' | sort -u)"
  while IFS= read -r k; do grep -qx "$k" <<<"$KEYS" || { MISSING=$((MISSING+1)); detail "cited, not in biblio.bib: $k"; }; done <<<"$CITED"
  while IFS= read -r k; do grep -qx "$k" <<<"$CITED" || UNCITED=$((UNCITED+1)); done <<<"$KEYS"
  if [ "$MISSING" -eq 0 ] && [ "$UNCITED" -eq 0 ]; then
    green "5. bibliography: $N_CITE citations, all keys resolve, no dead entries"
  else
    red "5. bibliography: $MISSING unresolved \\cite, $UNCITED uncited entry(ies) in biblio.bib"
  fi
fi

# --- 6. six figures, captioned ------------------------------------------------------
N_FIG=$(grep -c '\\includegraphics' "$BODY_TMP" || true)
N_CAP=$(grep -c '\\caption{' "$BODY_TMP" || true)
if [ "$N_FIG" -ge 6 ] && [ "$N_CAP" -ge "$N_FIG" ]; then
  green "6. figures: $N_FIG included, $N_CAP captioned"
else
  red "6. figures: $N_FIG of 6 included, $N_CAP captioned"
fi

# --- 7. front/back matter -----------------------------------------------------------
MISSING_MATTER=""
for s in "Data and code availability" "Author contributions" "Competing interests" "Funding" "Acknowledgements"; do
  n=$(awk -v s="$s" 'index($0,s){p=1;next} /\\section\*?\{/{p=0} p{l=$0; sub(/^[^:]*:[0-9]*:/,"",l); gsub(/^[ \t]+/,"",l); if(l!="" && substr(l,1,1)!="%") c++} END{print c+0}' "$BODY_TMP")
  [ "$n" -eq 0 ] && MISSING_MATTER="$MISSING_MATTER  $s"
done
if [ -z "$MISSING_MATTER" ]; then green "7. front/back matter complete"; else
  red "7. front/back matter empty:$MISSING_MATTER"
fi

# --- 8. word count ------------------------------------------------------------------
# Methods are excluded from eLife's count. The limit depends on D-1 (Research Article vs
# Tools & Resources); until D-1 is answered this reports, it does not judge.
if command -v detex >/dev/null 2>&1; then W=$(detex "$TEX" | wc -w); else W=$(sed 's/%.*//' "$BODY_TMP" | wc -w); fi
warn "8. word count ~$W (limit pending D-1: Research Article vs Tools & Resources)"

echo
printf 'pass %d   fail %d   warn %d\n' "$PASS" "$FAIL" "$WARN"
[ "$VERBOSE" = 0 ] && [ "$FAIL" -gt 0 ] && echo "(re-run with -v to list the offending lines)"
exit $(( FAIL > 0 ? 1 : 0 ))
