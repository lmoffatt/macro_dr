# Item (f) closure: Gaussian analytic Fisher vs finite-difference Fisher, full roster,
# from the fisher_only lane (a202e03, nsim 10000, noise 0.1, seed 20260814).
# Reads the FULL Likelihood_Gaussian_Fisher_Distortion matrix per interval
# (probit=mean, statistic=value, operation=probit) and eigendecomposes it here,
# NEVER the emitted spectrum/scalars: where F_b fails strict PSD the congruence
# helpers emit a ZERO matrix (jobs_to_be_run.md aviso 1, corrected 2026-08-15),
# detectable as an all-zero matrix -> counted as non-PSD, excluded from the median.
import csv, glob, math, collections
import numpy as np
FILES = sorted(glob.glob('projects/eLife_2025/figures/data/a202e03/figure_3_fim_nch_*_nsim_10000_*_noise_0.1_battery_sim.csv'))
MEMBERS = ['nonlinearsqr_g','nonlinearsqr','macro_NR','macro_INR','macro_R','macro_MR','macro_VR','macro_IR']
out = {}
for f in FILES:
    import re
    m = re.search(r'nch_(\d+)_nsim_10000_(.+)_noise_0.1_battery_sim', f)
    nch, algo = int(m.group(1)), m.group(2)
    if algo not in MEMBERS: continue
    mats = collections.defaultdict(dict)
    with open(f) as fh:
        fh.readline()
        for r in csv.DictReader(fh):
            if r['variable']!='Likelihood_Gaussian_Fisher_Distortion': continue
            if r['probit']!='mean' or r['statistic']!='value' or r['operation']!='probit': continue
            iv = float(r['interval_in_tau'])
            mats[iv][(int(r['value_row']), int(r['value_col']))] = float(r['value'])
    mags, bad = [], 0
    per_iv = {}
    for iv, d in sorted(mats.items()):
        n = max(k[0] for k in d)+1
        M = np.zeros((n,n))
        for (i,j),v in d.items(): M[i,j]=v
        if np.allclose(M,0): bad += 1; per_iv[iv]='nonPSD'; continue
        lam = np.linalg.eigvalsh((M+M.T)/2)
        if lam.min() <= 0: bad += 1; per_iv[iv]='nonpos'; continue
        mag = math.exp(np.mean(np.log(lam)))
        mags.append(mag); per_iv[iv]=round(mag,3)
    med = float(np.median(mags)) if mags else float('nan')
    out[(algo,nch)] = (med, min(mags) if mags else None, max(mags) if mags else None, bad, len(mats), n, per_iv)
print(f"{'member':16s}{'N_ch':>7s}  {'median':>7s} {'min':>7s} {'max':>7s}  nonPSD/ivs dim")
for algo in MEMBERS:
    for nch in (10,100,1000,10000):
        if (algo,nch) not in out: continue
        med,lo,hi,bad,niv,n,per = out[(algo,nch)]
        print(f'{algo:16s}{nch:7d}  {med:7.3f} {lo:7.3f} {hi:7.3f}  {bad:6d}/{niv:<3d} {n}x{n}')
print()
for algo in MEMBERS:
    for nch in (10,10000):
        if (algo,nch) in out: print(algo, nch, 'per-interval:', out[(algo,nch)][6])
