"""D-4 recompute on the FREEZE commit 1c2ae6f, Gaussian-Fisher anchor.
Mirrors D-4 section 7's recipe, with two deliberate changes:
  (a) commit 433ed13 -> 1c2ae6f (the freeze; 433ed13 is the numerical-Fisher demo)
  (b) variable Likelihood_Fisher_Covariance -> Gaussian_Fisher_Covariance,
      Likelihood_Information_Distortion -> Likelihood_Gaussian_Information_Distortion
      (the Gaussian Fisher is the program's declared anchor, _program/decisions.md sec.3)
"""
import csv, itertools, math, os, sys
import numpy as np

D = '/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data/1c2ae6f/'
PAR = {0: 'k_on', 1: 'k_off', 2: 'i', 5: 'N_ch'}     # the four the paper reports
PIDX = [0, 1, 2, 5]
GROUP = '10'
IV = '0.1'

def rows(path):
    f = open(path)
    next(f)                       # provenance line (git hash)
    for r in csv.DictReader(f):
        yield r
    f.close()

def emp_cov(algo, nch, noise):
    p = f'{D}figure_3_G_nch_{nch}_nsim_10000_macro_{algo}_noise_{noise}_mle_cloud_runs.csv'
    if not os.path.exists(p):
        return None
    runs = {}
    for r in rows(p):
        if (r['variable'] != 'Model_Parameters_Hat' or r['group_size'] != GROUP
                or r['interval_in_tau'] != IV):
            continue
        key = (r['simulation_index'], r['sample_index'], r['sub_index'], r['segment_index'])
        runs.setdefault(key, {})[int(r['param_index'])] = float(r['value'])
    M = np.array([[v[i] for i in PIDX] for v in runs.values() if all(i in v for i in PIDX)])
    return (np.cov(M, rowvar=False), M.shape[0]) if len(M) > 2 else None

def pool_matrix(algo, nch, noise, variable):
    p = f'{D}figure_3_G_nch_{nch}_nsim_10000_macro_{algo}_noise_{noise}_battery_pool_G.csv'
    if not os.path.exists(p):
        return None
    K = np.full((6, 6), np.nan)
    for r in rows(p):
        if (r['variable'] != variable or r['probit'] != 'mean' or r['statistic'] != 'value'
                or r['calculus'] != 'primitive' or r['operation'] != 'probit' or r['interval_in_tau'] != IV):
            continue
        if r['param_index'] == '' or r['param_col'] == '':
            continue
        K[int(r['param_index']), int(r['param_col'])] = float(r['value'])
    return K if not np.all(np.isnan(K)) else None

def diag_distortion(algo, nch, noise):
    p = f'{D}figure_3_G_nch_{nch}_nsim_10000_macro_{algo}_noise_{noise}_battery_pool_G.csv'
    if not os.path.exists(p):
        return {}
    out = {}
    for r in rows(p):
        if (r['variable'] != 'Likelihood_Gaussian_Information_Distortion'
                or r['probit'] != 'mean' or r['statistic'] != 'value'
                or r['interval_in_tau'] != IV or r['calculus'] != 'primitive' or r['operation'] != 'probit'):
            continue
        if r['param_index'] == '' or r['param_index'] != r['param_col']:
            continue
        i = int(r['param_index'])
        if i in PAR:
            out[PAR[i]] = float(r['value'])
    return out

def area(S):
    w = np.linalg.eigvalsh(S)
    return math.sqrt(abs(np.prod(w)))

print('=' * 78)
print('D-4 RECOMPUTE ON THE FREEZE COMMIT 1c2ae6f  (Gaussian-Fisher anchor)')
print('Headline cell: N_ch 100, noise label 0.1, interval 0.1 tau, group_size 10, nsim 10000')
print('=' * 78)

for algo in ['NR', 'R', 'MR', 'IR']:
    e = emp_cov(algo, 100, '0.1')
    F = pool_matrix(algo, 100, '0.1', 'Gaussian_Fisher_Covariance')
    C = pool_matrix(algo, 100, '0.1', 'Gaussian_Distortion_Corrected_Covariance')
    if e is None or F is None:
        print(f'\n{algo}: MISSING  (emp={e is not None}, fisher={F is not None})')
        continue
    E, n = e
    sub = np.ix_(PIDX, PIDX)
    Fm = F[sub] / 10.0                      # cov_scale = 1/group_size
    Cm = C[sub] / 10.0 if C is not None else None
    ratios, cratios = [], []
    for a, b in itertools.combinations(range(4), 2):
        ij = np.ix_([a, b], [a, b])
        ratios.append(area(E[ij]) / area(Fm[ij]))
        if Cm is not None:
            cratios.append(area(E[ij]) / area(Cm[ij]))
    g = math.exp(sum(math.log(x) for x in ratios) / len(ratios))
    line = (f'\n{algo:3s}  n_MLE={n:5d}   emp/Fisher area ratio: '
            f'geom-mean {g:.2f}  (per-pair {min(ratios):.2f}-{max(ratios):.2f})')
    if cratios:
        gc = math.exp(sum(math.log(x) for x in cratios) / len(cratios))
        line += f'\n     emp/corrected: geom-mean {gc:.2f} ({min(cratios):.2f}-{max(cratios):.2f})'
    print(line)
    print('     joint 4-param area ratio: %.2f' % (area(E) / area(Fm)))
    d = diag_distortion(algo, 100, '0.1')
    if d:
        print('     diagonal D: ' + '  '.join(f'{k} {v:.2f}' for k, v in d.items()))

print('\n' + '=' * 78)
print('N_ch TREND at noise 0.1, interval 0.1 tau  (diagonal D)')
print('=' * 78)
for algo in ['NR', 'R', 'MR', 'IR']:
    print(f'\n{algo}:')
    for nch in [10, 100, 1000, 10000]:
        d = diag_distortion(algo, nch, '0.1')
        if d:
            print(f'   N_ch {nch:6d}  ' + '  '.join(f'{k} {v:.2f}' for k, v in d.items()))
        else:
            print(f'   N_ch {nch:6d}  (no data)')
