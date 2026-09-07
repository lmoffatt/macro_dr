"""Rewrite the CCO truth + prior and regenerate its exact COC twin + prior.

Run from the repo base:
    python3 projects/macroir_next/data/models/make_cco_coc_twin.py

Truth: scheme_CCO with kon 6.73, koff 166, gating_off 45.3 and gating_on set
below. gating_on was 743 in the 2026-09 campaign (runs e743fd5), which puts
P_open at 0.826 at 10 uM and saturates the occupancy; 157 puts it at 0.50
(the design value of the eLife figures). The COC twin is exact at ONE agonist
concentration x, here the 10 uM of the episodic pulse and the stationary
protocol (theory/macroir/notes/cco_coc_discrimination_plan.md, section 1):

    k1 = kon*x, km1 = koff, k2 = gating_on, km2 = gating_off
    theta_1 < theta_2 = roots of s^2 - (k1 + km1 + k2) s + k1 k2 = 0
    a_1 = (theta_2 - k2)/(theta_2 - theta_1)
    on*x = theta_1, inactivating_off = theta_2, off = a_1 km2,
    inactivating_on = km2 - off

Priors: log10-normal, variance 2, centred on the truth (CCO) and on the twin
(COC), the shape used since 2026-09-01. The equilibrium equivalence (same
P_open, same closed-dwell mixture) is checked before anything is written.
"""
import csv, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
GATING_ON = 157.0
X_UM = 10.0
COMMON = [("unitary_current", 1.0), ("Current_Noise", 0.001),
          ("Current_Baseline", 1.0), ("Num_ch_mean", 5000.0)]

def twin(kon, koff, g_on, g_off, x):
    k1 = kon * x
    S, P = k1 + koff + g_on, k1 * g_on
    th1 = (S - math.sqrt(S * S - 4 * P)) / 2
    th2 = (S + math.sqrt(S * S - 4 * P)) / 2
    a1 = (th2 - g_on) / (th2 - th1)
    off = a1 * g_off
    return {"on": th1 / x, "off": off, "inactivating_on": g_off - off,
            "inactivating_off": th2}

def p_open_cco(kon, koff, g_on, g_off, x):
    r = kon * x / koff
    return r * g_on / g_off / (1 + r + r * g_on / g_off)

def p_open_coc(on, off, i_on, i_off, x):
    r = on * x / off
    return r / (1 + r + r * i_on / i_off)

def write(model, rates, with_variance):
    rows = [["model_name", "i_par", "parameter_name", "parameter_transformation",
             "parameter_value", "transformed_mean"] + (["transformed_variance"] if with_variance else [])]
    for i, (name, v) in enumerate(rates + COMMON):
        rows.append([model, i, name, "Log10", repr(float(v)), repr(math.log10(v))]
                    + ([2] if with_variance else []))
    suffix = "_prior.csv" if with_variance else "_par.csv"
    stem = model if model == "scheme_CCO" else "scheme_COC_twin"
    csv.writer(open(HERE / (stem + suffix), "w", newline="")).writerows(rows)

kon, koff, g_off = 6.73, 166.0, 45.3
tw = twin(kon, koff, GATING_ON, g_off, X_UM)
po_a = p_open_cco(kon, koff, GATING_ON, g_off, X_UM)
po_b = p_open_coc(tw["on"], tw["off"], tw["inactivating_on"], tw["inactivating_off"], X_UM)
assert min(tw.values()) > 0, tw
assert abs(po_a - po_b) < 1e-12, (po_a, po_b)
# closed-dwell mixture: both branches' exit rates and weights must coincide
kout_a, kout_b = g_off, tw["off"] + tw["inactivating_on"]
assert abs(kout_a - kout_b) < 1e-12
print(f"scheme_CCO gating_on={GATING_ON:g}: P_open({X_UM:g} uM) = {po_a:.4f}")
print("twin scheme_COC:", {k: round(v, 6) for k, v in tw.items()})
if "--check" in sys.argv:
    sys.exit(0)
write("scheme_CCO", [("kon", kon), ("koff", koff), ("gating_on", GATING_ON), ("gating_off", g_off)], False)
write("scheme_CCO", [("kon", kon), ("koff", koff), ("gating_on", GATING_ON), ("gating_off", g_off)], True)
coc = [("on", tw["on"]), ("off", tw["off"]), ("inactivating_on", tw["inactivating_on"]),
       ("inactivating_off", tw["inactivating_off"])]
write("scheme_COC", coc, False)
write("scheme_COC", coc, True)
print("written: scheme_CCO_{par,prior}.csv scheme_COC_twin_{par,prior}.csv")
