import csv, collections
p="projects/eLife_2025/figures/figure_4_source_data/figure_4_source_data_distortion.csv"
rows=[]
with open(p) as f:
    f.readline()
    for r in csv.DictReader(f):
        if r['comp']!='total' or r['anchor']!='pool': continue
        try: r['m']=float(r['m']); r['Dconf']=float(r['Dconf'])
        except: continue
        rows.append(r)

def cell(r): return (r['noise'], r['Num_ch'], r['interval_in_tau'])
LO,HI=1/1.15,1.15
for key in ('m','Dconf'):
    print(f"\n=== counting on {key} ===")
    for algo in ('macro_IR','macro_R','macro_NR'):
        sub=[r for r in rows if r['algo']==algo]
        cells={cell(r) for r in sub}
        ok=sum(1 for r in sub if LO<=r[key]<=HI)
        print(f"{algo:11s} {ok:5d} of {len(sub):5d} = {100*ok/len(sub):5.1f}%   cells={len(cells)}")
# common sub-grid over the three
cs={}
for algo in ('macro_IR','macro_R','macro_NR'):
    cs[algo]={cell(r) for r in rows if r['algo']==algo}
common=cs['macro_IR']&cs['macro_R']&cs['macro_NR']
print(f"\ncommon cells over IR/R/NR: {len(common)}")
for key in ('m','Dconf'):
    print(f"--- on {key}, common sub-grid ---")
    for algo in ('macro_IR','macro_R','macro_NR'):
        sub=[r for r in rows if r['algo']==algo and cell(r) in common]
        ok=sum(1 for r in sub if LO<=r[key]<=HI)
        params=collections.Counter(r['param'] for r in sub)
        print(f"{algo:11s} {ok:5d} of {len(sub):5d} = {100*ok/len(sub):5.1f}%  params={dict(params)}")
