import csv, glob, os, re
DIRS=['/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data/1c2ae6f/','/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data/87889e6/']
PAR={0:'k_on',1:'k_off',2:'i',5:'N_ch'}
def rows(p):
    f=open(p); next(f)
    for r in csv.DictReader(f): yield r
    f.close()
res={}
files=[]
for _d in DIRS: files+=glob.glob(_d+'figure_3_G_nch_*_nsim_10000_macro_*_noise_*_battery_pool_G.csv')
for p in sorted(files):
    m=re.search(r'nch_(\d+)_nsim_10000_macro_([A-Za-z]+)_noise_([0-9.]+)_battery_pool_G',os.path.basename(p))
    if not m: continue
    nch,algo,noise=int(m.group(1)),m.group(2),m.group(3)
    for r in rows(p):
        if (r['variable']!='Likelihood_Gaussian_Information_Distortion' or r['probit']!='mean'
            or r['statistic']!='value' or r['calculus']!='primitive' or r['operation']!='probit'): continue
        if r['param_index']=='' or r['param_index']!=r['param_col']: continue
        i=int(r['param_index'])
        if i not in PAR: continue
        res.setdefault(algo,[]).append((float(r['value']),nch,noise,r['interval_in_tau'],PAR[i]))
for algo in ['IR','R','MR','NR','NMR']:
    v=res.get(algo)
    if not v: continue
    v.sort()
    lo,hi=v[0],v[-1]
    n_cells=len(set((x[1],x[2],x[3]) for x in v))
    within=sum(1 for x in v if 0.85<=x[0]<=1.15)
    print(f'{algo}: {len(v)} (param,cell) points over {n_cells} cells')
    print(f'   min D = {lo[0]:.3f}  at N_ch {lo[1]}, noise {lo[2]}, interval {lo[3]}, param {lo[4]}')
    print(f'   max D = {hi[0]:.3f}  at N_ch {hi[1]}, noise {hi[2]}, interval {hi[3]}, param {hi[4]}')
    print(f'   within 1 +- 0.15: {within}/{len(v)} = {100*within/len(v):.0f}%')
    print(f'   N_ch decades present: {sorted(set(x[1] for x in v))}')
    print(f'   noise levels present: {sorted(set(x[2] for x in v), key=float)}')
