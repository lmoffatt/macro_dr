#!/usr/bin/env python3
"""Measure ABSTRACT length of published eLife articles.

Same sample and same API as elife_length_survey.py, beside this file: the ids in
elife_length_survey_data.csv (413 VOR articles, Structural Biology and Molecular
Biophysics + Computational and Systems Biology, research-article and tools-resources).

Counted eLife's way: the abstract body only.  The structured-abstract section titles
("Background", "Methods", ...) are dropped, and eLife's separate one-sentence "impact
statement" is never part of the abstract and is not counted here.
"""
import csv, json, os, re, statistics, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
IDS = os.path.join(HERE, "elife_length_survey_data.csv")
CACHE = "/tmp/claude-1000/elife_abs/cache"
os.makedirs(CACHE, exist_ok=True)
UA = {"User-Agent": "abstract-survey/1.0 (research use)"}
TAG = re.compile(r"<[^>]+>")


def get(aid):
    path = os.path.join(CACHE, aid + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    for k in range(4):
        try:
            u = "https://api.elifesciences.org/articles/" + aid
            req = urllib.request.Request(u, headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                d = json.loads(r.read().decode())
            json.dump(d, open(path, "w"))
            return d
        except Exception:
            time.sleep(1.5 * (k + 1))
    return None


def text_of(blocks, out):
    """Walk eLife's content block tree, collecting paragraph text."""
    for b in blocks or []:
        if not isinstance(b, dict):
            continue
        t = b.get("type")
        if t == "paragraph":
            out.append(TAG.sub(" ", b.get("text", "")))
        elif t == "section":
            text_of(b.get("content"), out)          # title dropped on purpose
        elif t in ("list", "quote", "box"):
            text_of(b.get("content") or b.get("items"), out)
        elif t == "mathml":
            out.append("x")
    return out


def words(s):
    return len(s.split())


def main():
    ids, meta = [], {}
    for row in csv.DictReader(open(IDS)):
        ids.append(row["id"])
        meta[row["id"]] = row["type"]

    with ThreadPoolExecutor(max_workers=8) as ex:
        arts = list(ex.map(get, ids))

    rows, missing = [], 0
    for aid, d in zip(ids, arts):
        if not d:
            missing += 1
            continue
        abs_ = d.get("abstract")
        if not abs_:
            missing += 1
            continue
        n = words(" ".join(text_of(abs_.get("content"), [])))
        if n == 0:
            missing += 1
            continue
        rows.append((aid, meta[aid], d.get("published", "")[:10], n))

    out = os.path.join(HERE, "elife_abstract_survey_data.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "type", "published", "abstract_words"])
        w.writerows(rows)

    def report(label, sel):
        v = sorted(n for _, _, _, n in sel)
        if not v:
            return
        q = statistics.quantiles(v, n=100, method="inclusive")
        print(f"{label:28s} n={len(v):4d}  median={statistics.median(v):5.0f} "
              f" p25={q[24]:4.0f} p75={q[74]:4.0f} p90={q[89]:4.0f} "
              f" max={max(v):4d}")

    report("all", rows)
    report("research-article", [r for r in rows if r[1] == "research-article"])
    report("tools-resources", [r for r in rows if r[1] == "tools-resources"])
    v = sorted(n for _, _, _, n in rows)
    for thr in (150, 200, 220, 239, 250, 260, 300):
        over = sum(1 for n in v if n > thr)
        pct = 100.0 * sum(1 for n in v if n <= thr) / len(v)
        print(f"  over {thr:3d}: {over:4d} of {len(v)} ({100*over/len(v):4.1f}%)"
              f"   -> {thr} sits at percentile {pct:4.1f}")
    print("missing/unparsed:", missing, "| data:", out)


if __name__ == "__main__":
    main()
