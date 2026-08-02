"""The interval window x recursion factorial, on the CORRECTED INR run (1f7138b, 2026-07-31).

WHY IT EXISTS. decisions.md used to read the roster as a monotone ladder of cost, and its
INR entry was open: INR had been measured only as NMR, the build missing the N*ms
interval-variance term, where it was numerically indistinguishable from NR. With the term
restored the four members NR / INR / R / IR close a 2x2 (interval window x recursion) and
each margin turns out to move a DIFFERENT moment. This script measures both margins.

WHAT IT READS. battery_sim_G for the bias (anchored at theta_sim; at theta_pool the score
vanishes and the bias is zero by construction) and battery_pool_G for the three distortion
components, exactly the split figure_4_data.R uses.

THE LSE INDEX TRAP, and why this does not go through figures/figure_4_source_data. The macro
members carry the full parameter vector, 0=k_on 1=k_off 2=i 3=Current_Noise 4=Baseline
5=N_ch. The LSE run holds unitary_current and Current_Noise Fixed (ops/local/figure_4_LSE.macroir:54-55),
so its free vector RENUMBERS to 0=k_on 1=k_off 2=Baseline 3=N_ch. figure_4_data.R selects on
`param_index %in% PIDX` and labels by position in PIDX, so for the LSE column it prints
Baseline under the name unitary_current, N_ch under the name Current_Noise, and leaves the
N_ch panel empty. Known: figure_4_supplement_other_parameters.Rmd:32-44 excludes those rows,
and figures_build_plan.md section 4b records that the shared reader is still unfixed. Nothing
rendered is affected. Per-algorithm maps below, which is what lets this script report LSE's
real N_ch numbers at all.
"""
import csv, glob, os, re, statistics as st

ROOT = '/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data/'
# 1f7138b is the corrected INR run; the other three are the freeze.
DIRS = [ROOT + d + '/' for d in ('1f7138b', '1c2ae6f', '87889e6', '0ffbda7')]

MACRO_PAR = {0: 'k_on', 1: 'k_off', 2: 'i', 3: 'noise', 4: 'baseline', 5: 'N_ch'}
LSE_PAR = {0: 'k_on', 1: 'k_off', 2: 'baseline', 3: 'N_ch'}

BIAS = 'Probit_statistics_Gaussian_Distortion_Induced_Bias'
COMPS = {'total': 'Probit_statistics_Likelihood_Gaussian_Information_Distortion',
         'sample': 'Probit_statistics_Gaussian_Sample_Distortion',
         'corr': 'Probit_statistics_Likelihood_Correlation_Distortion'}

ORDER = ['nonlinearsqr', 'NR', 'INR', 'R', 'MR', 'VR', 'IR']
CELL = {'nonlinearsqr': 'off the lattice: fits the mean only',
        'NR': 'no window, no recursion', 'INR': 'window, no recursion',
        'R': 'no window, recursion', 'MR': 'one end, recursion',
        'VR': 'one end + residual variance, recursion', 'IR': 'both ends, recursion'}


def rows(p):
    with open(p) as f:
        next(f)                      # line 1 is the provenance stamp (the git hash)
        for r in csv.DictReader(f):
            yield r


def find(algo, nch, noise, kind):
    stem = ('figure_3_LSE_nch_%s_nsim_10000_nonlinearsqr_noise_%s_%s.csv'
            if algo == 'nonlinearsqr' else
            'figure_3_G_nch_%s_nsim_10000_macro_' + algo + '_noise_%s_%s.csv')
    for d in DIRS:
        hit = glob.glob(d + stem % (nch, noise, kind))
        if hit:
            return hit[0]
    return None


def collect(algo, nch, noise, par):
    """medians over the recording's intervals, for one (member, cell, parameter)"""
    pmap = LSE_PAR if algo == 'nonlinearsqr' else MACRO_PAR
    idx = {v: k for k, v in pmap.items()}.get(par)
    out = {'bias': [], 'total': [], 'sample': [], 'corr': []}
    if idx is None:
        return out
    fs = find(algo, nch, noise, 'battery_sim_G')
    if fs:
        for r in rows(fs):
            if (r['component_path'] == BIAS and r['param_index'] == str(idx)
                    and r['statistic'] == 'value' and r['probit'] == 'mean'):
                out['bias'].append(float(r['value']))
    fp = find(algo, nch, noise, 'battery_pool_G')
    if fp:
        for r in rows(fp):
            if (r['param_index'] == str(idx) and r['param_col'] == str(idx)
                    and r['statistic'] == 'value' and r['probit'] == 'mean'):
                for k, c in COMPS.items():
                    if r['component_path'] == c:
                        out[k].append(float(r['value']))
    return out


def med(v):
    return st.median(v) if v else float('nan')


if __name__ == '__main__':
    NOISE = '0.1'
    for par in ('k_off', 'N_ch'):
        print('\n=== %s, noise %s, median over the intervals ===' % (par, NOISE))
        print('%-13s %-38s %6s %8s %8s %8s %8s'
              % ('member', 'cell', 'N_ch', 'bias', 'total', 'sample', 'corr'))
        for algo in ORDER:
            for nch in ('10', '100', '1000', '10000'):
                c = collect(algo, nch, NOISE, par)
                if not c['total'] and not c['bias']:
                    continue
                print('%-13s %-38s %6s %8.3f %8.3f %8.3f %8.3f'
                      % (algo, CELL[algo], nch, med(c['bias']), med(c['total']),
                         med(c['sample']), med(c['corr'])))

    # The bias summary the roster claim rests on: pooled over EVERY cell and interval on disk,
    # because the amplitude bias is flat in N_ch and in noise and a single cell understates it.
    print('\n=== median |bias| over all cells and intervals ===')
    print('%-13s %12s %12s %12s' % ('member', 'k_off', 'i', 'N_ch'))
    for algo in ORDER:
        pmap = LSE_PAR if algo == 'nonlinearsqr' else MACRO_PAR
        line = []
        for par in ('k_off', 'i', 'N_ch'):
            v = []
            pat = ('figure_3_LSE_nch_*_nsim_10000_nonlinearsqr_noise_*_battery_sim_G.csv'
                   if algo == 'nonlinearsqr' else
                   'figure_3_G_nch_*_nsim_10000_macro_%s_noise_*_battery_sim_G.csv' % algo)
            idx = {b: a for a, b in pmap.items()}.get(par)
            if idx is not None:
                for d in DIRS:
                    for p in glob.glob(d + pat):
                        for r in rows(p):
                            if (r['component_path'] == BIAS and r['param_index'] == str(idx)
                                    and r['statistic'] == 'value' and r['probit'] == 'mean'):
                                v.append(abs(float(r['value'])))
            line.append('%12.3f' % med(v) if v else '%12s' % 'Fixed')
        print('%-13s %s' % (algo, ''.join(line)))
