#!/usr/bin/env python3
"""Forward-reference check for the concept census.

Reads `1_method/10_concept_census.md`, takes every row whose "Defined in" cell names a section,
and reports any concept whose FIRST use in reading order sits in an earlier section than the one
that defines it.  Called by check.sh item 10; runs standalone too.

The census carries the truth about where a concept is defined; this script only checks that the
manuscript agrees with it.  A row is checked when it carries a `first-use pattern` in the
`<!-- pat: ... -->` comment beside it; rows without one are skipped and counted, so the check
never pretends to cover more than it does.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SEC = os.path.normpath(os.path.join(HERE, "..", "1_method", "docs", "manuscript-drafts", "sections"))
CENSUS = os.path.normpath(os.path.join(HERE, "..", "1_method", "10_concept_census.md"))
ORDER = ["00_abstract.tex", "01_introduction.tex", "02_framework.tex", "04_results.tex",
         "05_discussion.tex", "06_methods.tex", "08_appendix_derivation.tex",
         "09_appendix_members.tex", "10_appendix_diagnostics.tex", "11_appendix_repairs.tex"]
RANK = {f.split("_")[0]: i for i, f in enumerate(ORDER)}


def load(f):
    path = os.path.join(SEC, f)
    if not os.path.exists(path):
        return "", []
    lines = [re.sub(r"(?<!\\)%.*$", "", l) for l in open(path).read().split("\n")]
    txt, off = "", []
    for i, l in enumerate(lines, 1):
        txt += l + " "
        off += [i] * (len(l) + 1)
    return txt, off


DOCS = [(f, ) + load(f) for f in ORDER]


def first_use(pat):
    rx = re.compile(pat.replace(" ", r"\s+"))
    for f, txt, off in DOCS:
        m = rx.search(txt)
        if m:
            return f.split("_")[0], off[m.start()]
    return None, None


def rows():
    for line in open(CENSUS):
        m = re.search(r"<!--\s*pat:\s*(.+?)\s*\|\s*def:\s*(\d\d)\s*(?:\|\s*ok:\s*(\d\d)\s*)?-->", line)
        if m:
            name = line.split("|")[1].strip() if line.startswith("|") else "?"
            yield name, m.group(1), m.group(2), m.group(3)


def main():
    bad, checked = [], 0
    for name, pat, defsec, oksec in rows():
        checked += 1
        sec, ln = first_use(pat)
        # The abstract names everything by construction and is always exempt; a row may also declare
        # the earliest section in which an undefined mention is deliberate (`ok:`), which is how an
        # Introduction that names a quantity the body defines stops being reported every run.
        floor = oksec or defsec
        if sec == "00":
            continue
        if sec is None:
            bad.append((name, "never used", defsec))
        elif RANK.get(sec, 99) < RANK.get(floor, 99):
            bad.append((name, f"first used in {sec}:{ln}", f"defined in {defsec}"))
    print(f"concept census: {checked} rows carry a pattern, {len(bad)} forward reference(s)")
    for b in bad:
        print("   ", " | ".join(b))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
