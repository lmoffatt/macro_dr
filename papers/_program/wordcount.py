#!/usr/bin/env python3
"""Count a manuscript's main text the way eLife defines it.

Main text = every body section except Materials and Methods, with LaTeX comments,
figure and table captions, equations, and the back matter excluded.  Appendices are
counted separately and never inside the main text.

Usage:
    python3 wordcount.py                       # the default section set of paper 1
    python3 wordcount.py path/to/sections      # any directory of section .tex files
    python3 wordcount.py --by-subsection       # per-subsection breakdown

Targets it prints against come from elife_main_text_length_survey.md (n = 413).
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.normpath(os.path.join(
    HERE, "..", "1_method", "docs", "manuscript-drafts", "sections"))

# (file, counts toward main text?)
MAIN = ["01_introduction.tex", "02_theory_full.tex", "03_diagnostics.tex",
        "04_results.tex", "05_discussion.tex"]
EXCLUDED = ["00_abstract.tex", "06_methods.tex", "07_backmatter.tex",
            "08_appendix_derivation.tex", "09_appendix_members.tex"]

BENCH = [("median", 5209), ("p75", 6543), ("p90", 8017), ("p95", 9066),
         ("Munch e62714", 11515), ("longest of 413", 13300)]


def strip_comments(text):
    out = []
    for line in text.split("\n"):
        if line.lstrip().startswith("%"):
            continue
        # an unescaped % starts a comment
        out.append(re.sub(r"(?<!\\)%.*$", "", line))
    return "\n".join(out)


def brace_span(s, i):
    """index just past the group whose opening brace is at or after i"""
    while i < len(s) and s[i] not in "{[":
        i += 1
    if i >= len(s):
        return len(s), ""
    if s[i] == "[":                       # optional argument first
        d = 0
        while i < len(s):
            if s[i] == "[":
                d += 1
            elif s[i] == "]":
                d -= 1
                if d == 0:
                    i += 1
                    break
            i += 1
        while i < len(s) and s[i] != "{":
            i += 1
    d, j = 0, i
    while j < len(s):
        if s[j] == "{":
            d += 1
        elif s[j] == "}":
            d -= 1
            if d == 0:
                return j + 1, s[i + 1:j]
        j += 1
    return len(s), s[i + 1:]


def remove_macro_groups(text, macros):
    """delete \\macro[...]{...} and return (text, list of removed bodies)"""
    removed = []
    for mac in macros:
        while True:
            m = re.search(r"\\" + mac + r"\b", text)
            if not m:
                break
            end, body = brace_span(text, m.end())
            if mac == "figsupp":                       # \figsupp[alt]{caption}{art}
                end2, body2 = brace_span(text, end)
                body, end = body2 and body or body, end2
            removed.append(body)
            text = text[:m.start()] + " " + text[end:]
    return text, removed


ENVS_DROP = ["equation", "align", "equation*", "align*", "gather", "gather*",
             "tabular", "tabularx", "verbatim"]


def words(text):
    text = strip_comments(text)
    text, caps = remove_macro_groups(text, ["caption", "captionof", "figsupp", "label",
                                            "includegraphics", "src"])
    for env in ENVS_DROP:
        text = re.sub(r"\\begin\{" + re.escape(env) + r"\}.*?\\end\{"
                      + re.escape(env) + r"\}", " ", text, flags=re.S)
    text = re.sub(r"\\\[.*?\\\]", " ", text, flags=re.S)          # display math
    text = re.sub(r"\$[^$]*\$", " x ", text)                      # inline math -> 1 word
    text = re.sub(r"\\[a-zA-Z]+\*?", " ", text)                   # remaining macros
    text = re.sub(r"[{}~^_&\\]", " ", text)
    n = len([w for w in text.split() if any(c.isalnum() for c in w)])
    capn = sum(len([w for w in re.sub(r"\\[a-zA-Z]+\*?", " ", c).split()
                    if any(ch.isalnum() for ch in w)]) for c in caps)
    return n, capn


def subsections(path):
    text = strip_comments(open(path).read())
    parts, cur, buf = [], "(opening)", []
    for line in text.split("\n"):
        m = re.match(r"\s*\\(?:sub)?section\*?\{(.+?)\}", line)
        if m:
            parts.append((cur, "\n".join(buf)))
            cur, buf = m.group(1)[:58], []
        else:
            buf.append(line)
    parts.append((cur, "\n".join(buf)))
    return parts


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    d = args[0] if args else DEFAULT
    by_sub = "--by-subsection" in sys.argv
    total = capstot = 0
    print(f"{'file':28s} {'main':>7s} {'captions':>9s}")
    for f in MAIN:
        p = os.path.join(d, f)
        if not os.path.exists(p):
            continue
        n, c = words(open(p).read())
        total += n
        capstot += c
        print(f"{f:28s} {n:7d} {c:9d}")
        if by_sub:
            for title, body in subsections(p):
                sn, sc = words(body)
                if sn:
                    print(f"    {title[:50]:50s} {sn:6d}")
    print(f"{'-'*46}")
    print(f"{'MAIN TEXT (counted)':28s} {total:7d} {capstot:9d}")
    for f in EXCLUDED:
        p = os.path.join(d, f)
        if os.path.exists(p):
            n, c = words(open(p).read())
            print(f"{'  not counted: ' + f:28s} {n:7d} {c:9d}")
    print()
    for name, v in BENCH:
        print(f"   vs {name:16s} {v:6d}   {total / v:5.2f}x"
              + ("   OVER" if total > v else ""))


if __name__ == "__main__":
    main()
