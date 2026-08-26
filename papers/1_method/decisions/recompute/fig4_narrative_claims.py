# fig4_narrative_claims.py -- measures the claims of the Figure 4 narrative, from the figure's own
# source data. Written 2026-08-26 to test, BEFORE rewriting the subsection, a reading of the plane
# that the current text does not make: that bias and distortion are organised by signal to noise
# and by the acquisition interval rather than member by member.
#
# Run from the repository root:  python3 papers/1_method/decisions/recompute/fig4_narrative_claims.py
#
# TRAPS OBSERVED
#   * m is the estimate; lo/hi are the interval. Dconf/Bconf are the CI EDGE (what the MAPS are
#     drawn from, figure_4_common.R:235-236) and must never be read as the estimate.
#   * THE TWO HALVES ARE ON DIFFERENT SCALES. In bias.csv m is the SIGNED bias in log10 units,
#     null at ZERO. In distortion.csv m is a variance ratio, null at ONE. Reading the bias as a
#     factor turns a 0.002 bias into "|log10| = 2.7" and every member looks catastrophic.
#   * the noise column is a LABEL: the dimensionless instrumental noise is S = 0.1 * label.
#   * distortion carries anchor (sim|pool) and comp (sample|corr|total); the body figure reads
#     anchor=pool, comp=total.
#   * nonlinearsqr_g is LSE (averaging=0), nonlinearsqr is ILSE (averaging=1).
import csv, math, statistics as st
from collections import defaultdict

SD = "projects/eLife_2025/figures/figure_4_source_data/"
NAME = {"macro_NR":"NR","macro_INR":"INR","macro_R":"R","macro_MR":"MR","macro_VR":"VR",
        "macro_IR":"IR","nonlinearsqr_g":"LSE","nonlinearsqr":"ILSE"}
ORDER = ["LSE","ILSE","NR","INR","R","MR","VR","IR"]
AVG   = {"INR","IR","ILSE"}
PARAMS = ["off","Num_ch_mean"]

def load(fn, null, **filt):
    """null = 0 for the bias (signed log10) and 1 for the distortion (a variance ratio)."""
    rows = []
    txt = open(SD+fn).read().split("\n",1)[1]
    for r in csv.DictReader(txt.splitlines()):
        if any(r[k]!=v for k,v in filt.items()): continue
        try: m = float(r["m"]); lo = float(r["lo"]); hi = float(r["hi"])
        except ValueError: continue
        if not math.isfinite(m): continue
        if null == 1 and m <= 0: continue          # greyed / non-positive cells
        # lab is the RAW noise column and S is the dimensionless noise, S = 0.1 * lab. Both are
        # carried because the roster split is on the label (IR alone was swept at 0.05, 0.2, 0.5)
        # while everything printed in the paper is in S.
        rows.append(dict(mem=NAME[r["algo"]], param=r["param"], m=m, lo=lo, hi=hi,
                         sig=(lo > null or hi < null), lab=float(r["noise"]),
                         S=0.1*float(r["noise"]), N=float(r["Num_ch"]),
                         d=float(r["interval_in_tau"])))
    return rows

B = load("figure_4_source_data_bias.csv", 0)
D = load("figure_4_source_data_distortion.csv", 1, anchor="pool", comp="total")
lg = abs        # the bias is already in log10 units
med = lambda v: st.median(v) if v else float("nan")

print("rows: bias %d, distortion %d (anchor=pool, comp=total)" % (len(B), len(D)))
print("params per member (bias):")
for m in ORDER:
    print("   %-5s %s" % (m, sorted({r['param'] for r in B if r['mem']==m})))

print("\n\n######## BIAS")
for p in PARAMS:
    print("\n--- |log10 bias| by member, parameter %s ---" % p)
    print("   %-5s %8s %8s %8s   %s" % ("mem","median","p90","max","cells with CI excluding 1"))
    for m in ORDER:
        rs = [r for r in B if r["mem"]==m and r["param"]==p]
        if not rs: continue
        v = sorted(lg(r["m"]) for r in rs)
        sig = [r for r in rs if r["sig"]]
        print("   %-5s %8.3f %8.3f %8.3f   %d/%d" %
              (m, med(v), v[int(.9*len(v))-1], v[-1], len(sig), len(rs)))

print("\n--- BIAS 1: where the three interval-averaged members exceed 0.05 in |log10| ---")
for m in ["INR","IR","ILSE"]:
    for p in PARAMS:
        rs = [r for r in B if r["mem"]==m and r["param"]==p]
        if not rs: continue
        bad = [r for r in rs if lg(r["m"])>0.05]
        if not bad: print("   %-5s %-13s none of %d" % (m,p,len(rs))); continue
        print("   %-5s %-13s %d of %d: N_ch %s | S %s | delta %s" %
              (m, p, len(bad), len(rs), sorted({int(r["N"]) for r in bad}),
               sorted({r["S"] for r in bad}), sorted({r["d"] for r in bad})))

print("\n--- BIAS 2 and 3: median |log10 bias| against the acquisition interval ---")
for p in PARAMS:
    print("   parameter %s" % p)
    for m in ORDER:
        rs = [r for r in B if r["mem"]==m and r["param"]==p]
        if not rs: continue
        out = ["d=%-5g %6.3f" % (d, med([lg(r["m"]) for r in rs if r["d"]==d]))
               for d in sorted({r["d"] for r in rs})]
        print("      %-5s %s %s" % (m, "AVG " if m in AVG else "    ", "  ".join(out)))

print("\n\n######## DISTORTION (factor m, 1 = calibrated)")
for p in PARAMS:
    print("\n--- distortion by member, parameter %s ---" % p)
    print("   %-5s %8s %8s %8s   %s" % ("mem","median","p90","max","outside [1/1.15, 1.15]"))
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p]
        if not rs: continue
        v = sorted(r["m"] for r in rs)
        out = [r for r in rs if r["m"]>1.15 or r["m"]<1/1.15]
        print("   %-5s %8.3f %8.3f %8.1f   %d/%d" %
              (m, med(v), v[int(.9*len(v))-1], v[-1], len(out), len(rs)))

print("\n--- DISTORTION 1: where IR leaves [1/1.15, 1.15] ---")
for p in PARAMS:
    bad = [r for r in D if r["mem"]=="IR" and r["param"]==p and (r["m"]>1.15 or r["m"]<1/1.15)]
    tot = len([r for r in D if r["mem"]=="IR" and r["param"]==p])
    byN = defaultdict(int); byS = defaultdict(int)
    for r in bad: byN[int(r["N"])] += 1; byS[r["S"]] += 1
    print("   %-13s %d/%d | by N_ch %s | by S %s" %
          (p, len(bad), tot, dict(sorted(byN.items())), {k: byS[k] for k in sorted(byS)}))

print("\n--- DISTORTION 2: distortion at the lowest SHARED noise (S=0.01), by channel count ---")
for p in PARAMS:
    print("   parameter %s" % p)
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p and abs(r["S"]-0.01)<1e-9]
        if not rs: continue
        out = ["N=%-6d %9.3f" % (int(N), med([r["m"] for r in rs if r["N"]==N]))
               for N in sorted({r["N"] for r in rs})]
        print("      %-5s %s" % (m, "  ".join(out)))

print("\n--- DISTORTION 3: median distortion against noise (median over delta and N_ch) ---")
for p in PARAMS:
    print("   parameter %s" % p)
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p]
        if not rs: continue
        out = ["%.3g:%.2f" % (S, med([r["m"] for r in rs if r["S"]==S]))
               for S in sorted({r["S"] for r in rs})]
        print("      %-5s %s" % (m, " ".join(out)))

print("\n\n######## SIGNAL TO NOISE: does the plane collapse on S/N_ch?")
print("Spread of log10(distortion) across channel counts at FIXED S/N_ch and fixed interval,")
print("against the spread at fixed S. Smaller first column means S/N_ch is the coordinate.")
for p in PARAMS:
    print("\n   parameter %s" % p)
    print("      %-5s %16s %14s" % ("mem","spread|S/N fixed","spread|S fixed"))
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p]
        if not rs: continue
        def spread(keyf):
            g = defaultdict(list)
            for r in rs: g[keyf(r)].append(math.log10(r["m"]))
            w = [max(v)-min(v) for v in g.values() if len(v)>1]
            return med(w)
        print("      %-5s %16.3f %14.3f" %
              (m, spread(lambda r: (round(math.log10(r["S"]/r["N"]),6), r["d"])),
                  spread(lambda r: (r["S"], r["d"]))))

print("\n--- the IR exception: at fixed S/N_ch and interval, does distortion fall as N_ch rises? ---")
for p in PARAMS:
    print("   parameter %s" % p)
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p]
        if not rs: continue
        pairs = defaultdict(dict)
        for r in rs: pairs[(round(math.log10(r["S"]/r["N"]),6), r["d"])][r["N"]] = r["m"]
        up = dn = 0
        for k,v in pairs.items():
            Ns = sorted(v)
            if len(Ns) < 2: continue
            # DISTANCE FROM CALIBRATION, not the raw factor: IR sits BELOW one on the channel
            # number, where a rising m is an improvement and the raw comparison inverts it.
            if abs(math.log10(v[Ns[-1]])) < abs(math.log10(v[Ns[0]])): dn += 1
            else: up += 1
        print("      %-5s improves with N_ch in %d of %d matched pairs" % (m, dn, dn+up))

print("\n--- the interval: median distortion by delta ---")
for p in PARAMS:
    print("   parameter %s" % p)
    for m in ORDER:
        rs = [r for r in D if r["mem"]==m and r["param"]==p]
        if not rs: continue
        out = ["d=%-5g %8.2f" % (d, med([r["m"] for r in rs if r["d"]==d]))
               for d in sorted({r["d"] for r in rs})]
        print("      %-5s %s %s" % (m, "AVG " if m in AVG else "    ", "  ".join(out)))




# ============================================================================================
# THE NUMBERS QUOTED IN THE BODY. Added 2026-08-26 when the subsection was rewritten; the block was
# rebuilt the same day, when a review of the printed text against this file found six defects. What
# changed and why:
#   * THE BASIS. The interval paragraph used to be computed on S/N_ch <= 1, a sub-grid the printed
#     text never declared, and on that sub-grid NR falls 43x and R 16x where the FULL grid gives 28x
#     and 8.5x. Everything below now runs on the full grid, and the paragraph was rewritten to quote
#     the bias as a FACTOR ON THE PARAMETER rather than as a ratio of |log10| values: NR going from
#     1.054 to 0.038 in log10 is a fall from a factor of 11.3 to 9 per cent, not "a factor of forty".
#   * THE ROSTER AND THE CELL SET. IR was swept at three noise labels nobody else was (0.05, 0.2,
#     0.5), so it carries 294 cells against everyone else's 210 and 17 of its 28 closing-rate
#     failures sit on ground no other member covers. Every head-to-head calibration count below is
#     therefore taken on the SHARED noise labels; the bias medians are quoted on the full grid
#     because they are identical on both (IR 1 per cent either way) and the others only have one.
#   * THE PARTIAL CORRECTIONS. "Neither partial correction improves on the member it corrects" was
#     false in the first moment: VR beats R on the channel number at all four channel counts and MR
#     at three. The true statement, emitted below, is that both shave the bias and pay for it in the
#     second moment as channels are added.
# ============================================================================================
SHARED_LABELS = {0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0}
BODY = ["LSE", "ILSE", "NR", "INR", "R", "IR"]
shared = lambda rs: [r for r in rs if r["lab"] in SHARED_LABELS]
pc = lambda l10: 100 * (10 ** l10 - 1)          # |log10 bias| -> per cent on the parameter

print("\n\n######## NUMBERS QUOTED IN THE SUBSECTION, paragraph by paragraph")

print("\n[verdict] median |bias| on the channel number, full grid then shared noise labels")
for m in ORDER:
    a = [lg(r["m"]) for r in B if r["mem"] == m and r["param"] == "Num_ch_mean"]
    s = [lg(r["m"]) for r in shared(B) if r["mem"] == m and r["param"] == "Num_ch_mean"]
    if a: print("   %-5s full %.4f (%2.0f%%) n=%-4d | shared %.4f (%2.0f%%) n=%d"
                % (m, med(a), pc(med(a)), len(a), med(s), pc(med(s)), len(s)))

print("\n[verdict] the honesty count, SHARED noise labels, both drawn parameters")
print("   the printed 'five per cent against R's twenty-nine and about seventy' is this column")
for m in BODY:
    rs = [r for r in shared(D) if r["mem"] == m and r["param"] in PARAMS]
    bad = [r["m"] for r in rs if r["m"] > 1.15 or r["m"] < 1 / 1.15]
    print("   %-5s outside the 15%% band %3d/%3d = %2.0f%%   median miss, as the factor the BAR is off by %.2f"
          % (m, len(bad), len(rs), 100 * len(bad) / len(rs), med(bad) ** 0.5))

print("\n[verdict] direction of the miss, shared labels, body roster")
rs = [r for r in shared(D) if r["mem"] in BODY and r["param"] in PARAMS]
o = sum(1 for r in rs if r["m"] > 1.15); u = sum(1 for r in rs if r["m"] < 1 / 1.15)
print("   too narrow %d, too wide %d, i.e. %.2f of the misses are too narrow" % (o, u, o / (o + u)))

print("\n[verdict] the two partial corrections against R, by channel count")
for p in PARAMS:
    print("   parameter %s" % p)
    for src, lab in ((B, "bias |log10|"), (D, "distortion  ")):
        for m in ("R", "MR", "VR", "IR"):
            rows = [r for r in src if r["mem"] == m and r["param"] == p]
            if not rows: continue
            print("      %-12s %-3s %s" % (lab, m, "  ".join(
                "N=%-6d %.4f" % (N, med([r["m"] if src is D else lg(r["m"])
                                         for r in rows if r["N"] == N]))
                for N in sorted({r["N"] for r in rows}))))

print("\n[channels] closing-rate distortion at S=0.01, ten channels to ten thousand")
for m in ORDER:
    rs = [r for r in D if r["mem"] == m and r["param"] == "off" and abs(r["S"] - 0.01) < 1e-9]
    if rs: print("   %-5s %s" % (m, "  ".join("%.3f" % med([r["m"] for r in rs if r["N"] == N])
                                              for N in sorted({r["N"] for r in rs}))))

print("\n[channels] where IR's error bar fails, and how much of that corner is IR-only ground")
for p in PARAMS:
    bad = [r for r in D if r["mem"] == "IR" and r["param"] == p and (r["m"] > 1.15 or r["m"] < 1/1.15)]
    byN = defaultdict(int)
    for r in bad: byN[int(r["N"])] += 1
    only = [r for r in bad if r["lab"] not in SHARED_LABELS]
    print("   %-13s %d fail | by N_ch %s | %d of them at the IR-only labels %s"
          % (p, len(bad), dict(sorted(byN.items())), len(only),
             sorted({r["lab"] for r in only})))

print("\n[noise] median closing-rate distortion against the noise level")
print("   the printed 'from 13 to 1 over five decades' and 'R starts at 1.5, inside after two'")
for m in ("LSE", "R", "IR"):
    rs = [r for r in D if r["mem"] == m and r["param"] == "off"]
    print("   %-5s %s" % (m, " ".join("%.3g:%.2f" % (S, med([r["m"] for r in rs if r["S"] == S]))
                                      for S in sorted({r["S"] for r in rs}))))

print("\n[noise] cells whose bias half-width exceeds 0.05 in log10, i.e. blind to a bias of an eighth")
g = defaultdict(lambda: [0, 0]); tot = bad = 0
for r in B:
    if r["mem"] not in BODY or r["param"] not in PARAMS: continue
    wide = (r["hi"] - r["lo"]) / 2 > 0.05
    tot += 1; bad += wide
    k = round(math.log10(r["S"] / r["N"])); g[k][0] += 1; g[k][1] += wide
print("   %d of %d body-roster cells (%.0f%%); 0.05 in log10 IS 12.2 per cent, hence 'an eighth'"
      % (bad, tot, 100 * bad / tot))
print("   by S/N_ch: " + " ".join("1e%d:%.0f%%" % (k, 100 * g[k][1] / g[k][0]) for k in sorted(g)))

print("\n[interval] median |bias| on the channel number, FULL grid, as a factor on the parameter")
for m in ("NR", "INR", "R", "IR"):
    rows = [r for r in B if r["mem"] == m and r["param"] == "Num_ch_mean"]
    per_d = {d: med([lg(r["m"]) for r in rows if r["d"] == d]) for d in sorted({r["d"] for r in rows})}
    print("   %-5s coarsest %.4f (x%.2f on the parameter) | finest %.4f (%.1f%%) | worst interval %.1f%%"
          % (m, per_d[1.0], 10 ** per_d[1.0], per_d[0.01], pc(per_d[0.01]), pc(max(per_d.values()))))

print("\n[interval] median |bias| on the closing rate by interval, the two least-squares arms")
for m in ("LSE", "ILSE"):
    rs = [r for r in B if r["mem"] == m and r["param"] == "off"]
    print("   %-5s %s" % (m, "  ".join("d=%-5g %.4f" % (d, med([lg(r["m"]) for r in rs if r["d"] == d]))
                                       for d in sorted({r["d"] for r in rs}))))
d1 = {m: med([lg(r["m"]) for r in B if r["mem"] == m and r["param"] == "off" and r["d"] == 1.0])
      for m in ("LSE", "ILSE")}
print("   ratio at the coarsest interval: %.1f  (in per cent on the parameter, %.2f%% against %.2f%%)"
      % (d1["LSE"] / d1["ILSE"], pc(d1["LSE"]), pc(d1["ILSE"])))

print("\n[qualifications] R at N_ch = 1e4, Delta = 1: the amplitude trade-off the body quotes")
for p in ("unitary_current", "Num_ch_mean"):
    v = [r["m"] for r in B if r["mem"] == "R" and r["param"] == p and r["N"] == 1e4 and r["d"] == 1.0]
    print("   %-16s median %+.3f in log10 over the noise sweep (n=%d)" % (p, med(v), len(v)))
print("   they cancel: the product i * N_ch moves by %+.3f in log10"
      % (med([r["m"] for r in B if r["mem"] == "R" and r["param"] == "unitary_current"
              and r["N"] == 1e4 and r["d"] == 1.0])
         + med([r["m"] for r in B if r["mem"] == "R" and r["param"] == "Num_ch_mean"
                and r["N"] == 1e4 and r["d"] == 1.0])))
print("\n[qualifications] WHAT IS ACTUALLY GREY IN THE BODY FIGURE. figure_4.Rmd renders through")
print("   figure_4_layout.R::blk(), which has NO unident tile layer and does not set grey_offscale;")
print("   every value is clamped to the ends of the scale. The kappa > 3e4 machinery lives in")
print("   figure_4_common.R::mapblock(), the superseded producer. So grey is (a) the grey85 panel")
print("   background below the shared noise floor, IR alone being swept at label 0.05, and (b) cells")
print("   with no finite value. On the two drawn parameters that is NR alone:")
for r in D:
    pass
raw = [r for r in csv.DictReader(open(SD + "figure_4_source_data_distortion.csv").read()
                                 .split("\n", 1)[1].splitlines())
       if r["anchor"] == "pool" and r["comp"] == "total" and r["param"] in PARAMS]
def nonpos(v):
    try: return float(v) <= 0
    except ValueError: return True
blank = [r for r in raw if nonpos(r.get("Dconf", ""))]
for r in blank:
    print("      %-14s %-13s N_ch=%-6s noise=%-8s delta=%-5s Dconf=%s"
          % (r["algo"], r["param"], r["Num_ch"], r["noise"], r["interval_in_tau"], r["Dconf"]))
