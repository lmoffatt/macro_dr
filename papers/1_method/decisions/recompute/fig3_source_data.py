# Figure 3--source data 1 (2026-08-25, the digit-migration pass).
# The producer of every Figure 3 number is the caption-numbers chunk of
# projects/eLife_2025/figures/paper_both/figure_3.Rmd; its rendered output is figure_3.html.
# This script re-serializes that rendered table into a CSV, so the CSV is the producer's own
# output and no analysis logic is duplicated (drift-proof by construction). Re-knit the Rmd,
# re-run this, and the CSV follows.
import re, html, csv, os

SRC = "projects/eLife_2025/figures/paper_both/figure_3.html"
OUT_DIR = "projects/eLife_2025/figures/figure_3_source_data"
OUT = os.path.join(OUT_DIR, "figure_3_source_data_calibration.csv")

s = open(SRC, encoding="utf8", errors="replace").read()
pres = re.findall(r"<pre[^>]*>(.*?)</pre>", s, re.S)

rows = []
for p in pres:
    t = html.unescape(re.sub(r"<[^>]+>", "", p))
    if "D_bias" not in t or "F_accum" not in t:
        continue
    lines = [re.sub(r"^##\s*", "", l) for l in t.splitlines() if l.strip().startswith("##")]
    # The table prints in two blocks (wide data.frame wrap): first block carries
    # key..logL, second block logL_se F_fails G_fails. Parse the first block only;
    # logL_se comes from the second, matched by order. A <pre> holding the chunk's
    # R source also mentions the column names; only the OUTPUT pre has ## lines.
    head1 = next((i for i, l in enumerate(lines) if l.split()[:2] == ["key", "param"]), None)
    if head1 is None:
        continue
    block1, block2 = [], []
    for l in lines[head1 + 1:]:
        f = l.split()
        if not f:
            continue
        if f[0] == "logL_se":
            continue
        if len(f) >= 10 and not f[0] in ("TRUE", "FALSE"):
            block1.append(f)
        elif len(f) == 3:
            block2.append(f)
    for i, f in enumerate(block1):
        se = block2[i][0] if i < len(block2) else ""
        rows.append(dict(
            member=f[0], param=f[1], D_score_bias_frac=f[2], E_perinterval_frac=f[3],
            F_accumulated=f[4], F_lo=f[5], F_hi=f[6],
            G_score_acf_lag1=f[7], G_lo=f[8], G_hi=f[9], logL_mean=f[10], logL_se=se))
    break

assert rows, "caption-numbers table not found in figure_3.html"
os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("wrote", OUT, f"({len(rows)} rows)")
for r in rows[:3]:
    print(r)
