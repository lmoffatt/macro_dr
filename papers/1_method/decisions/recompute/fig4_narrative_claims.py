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
        rows.append(dict(mem=NAME[r["algo"]], param=r["param"], m=m, lo=lo, hi=hi,
                         sig=(lo > null or hi < null),
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
# THE NUMBERS QUOTED IN THE BODY. Added 2026-08-26 when the subsection was rewritten: everything
# the printed text asserts about Figure 4 is emitted here, in the order the paragraphs use it.
# ============================================================================================
print("\n\n######## NUMBERS QUOTED IN THE REWRITTEN SUBSECTION")

print("\n[verdict] median |bias| on the channel number, over the plane, log10 then per cent")
for m in ORDER:
    rs = [lg(r["m"]) for r in B if r["mem"]==m and r["param"]=="Num_ch_mean"]
    if rs: print("   %-5s %.3f  (%.0f%%)" % (m, med(rs), 100*(10**med(rs)-1)))

print("\n[verdict] direction of the distortion, both drawn parameters, all cells")
o = u = 0
for m in ORDER:
    c = [r["m"] for r in D if r["mem"]==m and r["param"] in PARAMS]
    o += sum(1 for x in c if x>1.15); u += sum(1 for x in c if x<1/1.15)
print("   too narrow %d, too wide %d, i.e. %.2f of the misses are too narrow" % (o,u,o/(o+u)))
print("   median miss among the failures, as the factor the error BAR is too narrow by:")
for m in ORDER:
    bad = [r["m"] for r in D if r["mem"]==m and r["param"] in PARAMS and (r["m"]>1.15 or r["m"]<1/1.15)]
    if bad: print("   %-5s variance %.2f -> bar %.2f" % (m, med(bad), med(bad)**0.5))

print("\n[channels] closing-rate distortion at S=0.01, ten channels to ten thousand")
for m in ORDER:
    rs = [r for r in D if r["mem"]==m and r["param"]=="off" and abs(r["S"]-0.01)<1e-9]
    if rs: print("   %-5s %s" % (m, "  ".join("%.3f"%med([r["m"] for r in rs if r["N"]==N])
                                              for N in sorted({r["N"] for r in rs}))))

print("\n[channels] where IR's error bar fails, by channel count and by noise")
for p in PARAMS:
    bad = [r for r in D if r["mem"]=="IR" and r["param"]==p and (r["m"]>1.15 or r["m"]<1/1.15)]
    byN = defaultdict(int); byS = defaultdict(int)
    for r in bad: byN[int(r["N"])] += 1; byS[r["S"]] += 1
    print("   %-13s %d cells | N_ch %s | S %s" % (p, len(bad), dict(sorted(byN.items())),
                                                  {k: byS[k] for k in sorted(byS)}))

print("\n[noise] median closing-rate distortion against the noise level")
for m in ("LSE","R","IR"):
    rs = [r for r in D if r["mem"]==m and r["param"]=="off"]
    print("   %-5s %s" % (m, " ".join("%.3g:%.2f" % (S, med([r["m"] for r in rs if r["S"]==S]))
                                      for S in sorted({r["S"] for r in rs}))))

print("\n[noise] cells whose bias half-width exceeds 0.05 in log10, i.e. blind to a 12% bias")
BODY = {"LSE","ILSE","NR","INR","R","IR"}
g = defaultdict(lambda: [0,0])
tot = bad = 0
for r in B:
    if r["mem"] not in BODY or r["param"] not in PARAMS: continue
    wide = (r["hi"]-r["lo"])/2 > 0.05
    tot += 1; bad += wide
    k = round(math.log10(r["S"]/r["N"]))
    g[k][0] += 1; g[k][1] += wide
print("   %d of %d body-roster cells (%.0f%%)" % (bad, tot, 100*bad/tot))
print("   by S/N_ch: " + " ".join("1e%d:%.0f%%" % (k, 100*g[k][1]/g[k][0]) for k in sorted(g)))

print("\n[interval] median |bias| on the channel number, coarsest against finest, S/N_ch <= 1")
for m in ORDER:
    for d in (1.0, 0.01):
        rs = [lg(r["m"]) for r in B if r["mem"]==m and r["param"]=="Num_ch_mean"
              and r["d"]==d and r["S"]/r["N"] <= 1]
        if rs: print("   %-5s d=%-5g %.3f  (%.1f%%)" % (m, d, med(rs), 100*(10**med(rs)-1)), end="")
    print()

print("\n[interval] median |bias| on the closing rate by interval, the two least-squares arms")
for m in ("LSE","ILSE"):
    rs = [r for r in B if r["mem"]==m and r["param"]=="off"]
    print("   %-5s %s" % (m, "  ".join("d=%-5g %.4f" % (d, med([lg(r["m"]) for r in rs if r["d"]==d]))
                                       for d in sorted({r["d"] for r in rs}))))
