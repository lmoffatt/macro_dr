#!/usr/bin/env python3
"""Measure main-text length of published eLife articles.

Main text is defined the way eLife's own advisory defines it: everything except
Materials and Methods, References, and figure legends.  Appendices are counted
separately, never inside the main text.
"""
import json, os, re, sys, time, urllib.request, collections

CACHE = "/tmp/claude-1000/elife_len/cache"
os.makedirs(CACHE, exist_ok=True)
UA = {"User-Agent": "length-survey/1.0 (research use)"}


def get(url, path=None, tries=4):
    if path and os.path.exists(path):
        return json.load(open(path))
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                d = json.loads(r.read().decode())
            if path:
                json.dump(d, open(path, "w"))
            return d
        except Exception as e:                      # transient API failures
            last = e
            time.sleep(1.5 * (k + 1))
    print("FAIL", url, last, file=sys.stderr)
    return None


def listing(subject, atype, want):
    """Most recent `want` VOR articles of one type in one subject."""
    out, page = [], 1
    while len(out) < want:
        u = (f"https://api.elifesciences.org/search?for=&subject[]={subject}"
             f"&type[]={atype}&per-page=100&page={page}&sort=date&order=desc")
        d = get(u)
        if not d or not d.get("items"):
            break
        for it in d["items"]:
            if it.get("status") == "vor" and it.get("type") == atype:
                out.append(it["id"])
        if len(d["items"]) < 100:
            break
        page += 1
    return out[:want]


TAG = re.compile(r"<[^>]+>")
CITE = re.compile(r"\s+")
METHODS = re.compile(
    r"^(materials?\s+and\s+methods?|methods?\s+and\s+materials?|methods?|"
    r"experimental\s+(procedures?|methods?|design)|star\s*.?\s*methods?|"
    r"materials?)\b", re.I)
# blocks whose text is a legend / caption / tabular content, excluded by the advisory
SKIP_BLOCK = {"figure", "image", "table", "video", "code", "mathml",
              "asset", "figure-supplement"}


def words(html):
    t = TAG.sub(" ", html or "")
    t = t.replace("&nbsp;", " ")
    t = CITE.sub(" ", t)
    return len([w for w in t.split() if any(c.isalnum() for c in w)])


def count(node, acc):
    """Walk, counting paragraph/list text only, never descending into assets."""
    if isinstance(node, list):
        for x in node:
            count(x, acc)
        return
    if not isinstance(node, dict):
        return
    t = node.get("type")
    if t in SKIP_BLOCK:
        acc["asset_blocks"] += 1
        if t in ("figure", "image"):
            acc["figures"] += 1
        if t == "table":
            acc["tables"] += 1
        return                                    # captions are excluded
    if t == "paragraph":
        acc["w"] += words(node.get("text", ""))
        return
    if t in ("list",):
        for item in node.get("items", []):
            if isinstance(item, str):
                acc["w"] += words(item)
            else:
                count(item, acc)
        return
    if t == "quote":
        count(node.get("text", []), acc)
        return
    if t in ("box", "section"):
        count(node.get("content", []), acc)
        return
    if "content" in node:
        count(node["content"], acc)


def measure(art):
    main = 0
    per_section = []
    methods_w = 0
    figs = tabs = 0
    for sec in art.get("body", []):
        title = (sec.get("title") or "").strip()
        acc = {"w": 0, "figures": 0, "tables": 0, "asset_blocks": 0}
        count(sec.get("content", []), acc)
        figs += acc["figures"]
        tabs += acc["tables"]
        if METHODS.match(title):
            methods_w += acc["w"]
        else:
            main += acc["w"]
            per_section.append((title, acc["w"]))
    app_w, app_n = 0, 0
    for ap in art.get("appendices", []) or []:
        acc = {"w": 0, "figures": 0, "tables": 0, "asset_blocks": 0}
        count(ap.get("content", []), acc)
        app_w += acc["w"]
        app_n += 1
    return dict(id=art["id"], type=art["type"], title=art["title"][:90],
                published=art.get("published", "")[:10],
                main=main, methods=methods_w, appendix_words=app_w,
                appendices=app_n, figures=figs, tables=tabs,
                sections=per_section,
                subjects=[s["id"] for s in art.get("subjects", [])])


if __name__ == "__main__":
    SUBJECTS = ["structural-biology-molecular-biophysics",
                "computational-systems-biology"]
    TYPES = ["research-article", "tools-resources"]
    WANT = {"research-article": 150, "tools-resources": 60}
    ids = collections.OrderedDict()
    for s in SUBJECTS:
        for t in TYPES:
            got = listing(s, t, WANT[t])
            print(f"{s:42s} {t:16s} {len(got)}", file=sys.stderr)
            for i in got:
                ids.setdefault(i, t)
    print(f"unique articles: {len(ids)}", file=sys.stderr)

    rows = []
    for n, (aid, t) in enumerate(ids.items(), 1):
        art = get(f"https://api.elifesciences.org/articles/{aid}",
                  os.path.join(CACHE, f"{aid}.json"))
        if not art or "body" not in art:
            continue
        rows.append(measure(art))
        if n % 25 == 0:
            print(f"  {n}/{len(ids)}", file=sys.stderr)
        time.sleep(0.06)
    json.dump(rows, open("/tmp/claude-1000/elife_len/rows.json", "w"))
    print(f"measured {len(rows)}", file=sys.stderr)
