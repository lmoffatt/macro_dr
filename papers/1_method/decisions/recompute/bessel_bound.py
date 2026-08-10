"""
How far the uniform acquisition window is from the real acquisition kernel, as a
function of one dimensionless group.

WHY THIS EXISTS. Theory, Methods and Discussion each concede that a real recording is
low-pass filtered (normally by a Bessel) before digitisation, while every likelihood in
the paper, and the simulator with them, uses the uniform average over the acquisition
window. All three concede it without bounding it, which leaves a referee free to read the
concession as "size unknown". This script bounds it. It computes no new physics: it is
the exact linear-filter algebra, evaluated numerically so the constants quoted in the
manuscript comments are measured rather than asserted.

WHAT IS COMPUTED. The real acquisition kernel is  boxcar(Delta) o Bessel(f_c).  The
likelihood implements the boxcar exactly and nothing else. For an input flat on the scale
of Delta (the instrumental noise, and also the gating fluctuation once its correlation
time is short against the window), the block average over one window has

    Var    = 2 S int_0^inf |H(f)|^2 sinc^2(pi f Delta) df
    Cov_k  = 2 S int_0^inf |H(f)|^2 sinc^2(pi f Delta) cos(2 pi k f Delta) df

with sinc(x) = sin(x)/x the boxcar transfer function and H the filter. Setting |H| = 1
gives Var = S/Delta and Cov_k = 0 exactly, which is the model's own premise, so the
departure of these two from (S/Delta, 0) IS the unmodelled part of the kernel.

RESULTS, at f_c*Delta >= 1 and to the precision printed (4-pole Bessel):

  * The whole departure is governed by f_c*Delta alone, and it is first order in its
    reciprocal:  variance deficit = 0.147/(f_c*Delta),  lag-1 correlation
    rho_1 = 0.075/(f_c*Delta).  Lag 2 and beyond are zero to four decimals, so the
    residue is a nearest-neighbour effect and nothing longer.
  * Var + 2*sum_k Cov_k = S/Delta to six decimals at every setting. The filter
    REDISTRIBUTES variance from the diagonal into the two neighbouring covariances and
    destroys none of it, which is forced by the filter having unit gain at DC.
  * A 4-pole Bessel, an 8-pole Bessel and the Gaussian of equal -3 dB frequency agree to
    the third decimal from f_c*Delta = 0.5 upward, so the number does not depend on the
    pole count and the Gaussian equivalent (sigma_f = sqrt(ln2)/(2 pi f_c) = 0.1325/f_c)
    can be quoted as the width.
  * The gating component is far better protected than the instrumental one. At the same
    f_c*Delta the deficit on a component of correlation time tau_c is second order in
    1/(f_c*Delta) rather than first, because a unit-DC-gain filter preserves the integral
    of the autocovariance and the block-average variance of a short-correlated component
    depends on that integral alone to leading order. At f_c*Delta = 5 the white term
    loses 2.94% and a gating term with tau_c = Delta loses 0.11%.
  * The one place the expansion does not apply is the interval holding a concentration or
    voltage jump: the filter's group delay plus its rise time displace and smear the
    transient by about 0.69/f_c for 4 poles, which is a fixed time and not a fraction of
    Delta.

Run: python3 bessel_bound.py     (~30 s; numpy + scipy only, no data files)
"""

import numpy as np
from scipy import signal

FC = 1.0  # everything is in units of the -3 dB frequency

# --------------------------------------------------------------------- filters
def bessel_ba(n, fc=FC):
    """n-pole analogue Bessel, normalised so that |H(fc)| = 1/sqrt(2)."""
    return signal.bessel(n, 2 * np.pi * fc, btype="low", analog=True, norm="mag")

def bessel_H2(n, fc=FC):
    b, a = bessel_ba(n, fc)
    def H2(f):
        _, h = signal.freqs(b, a, worN=2 * np.pi * np.atleast_1d(f))
        return np.abs(h) ** 2
    return H2

SIG_F = np.sqrt(np.log(2)) / (2 * np.pi * FC)   # Gaussian equivalent width, 0.1325/f_c

def gauss_H2(f):
    return np.exp(-4 * np.pi**2 * SIG_F**2 * np.atleast_1d(f) ** 2)

FILTERS = [("Bessel-4", bessel_H2(4)), ("Bessel-8", bessel_H2(8)), ("Gaussian", gauss_H2)]

# ------------------------------------------------------- block-average moments
def _simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-2:2].sum())

def moments(H2, D, kmax=8, S=1.0, fmax_fc=600.0, npts=2_000_001):
    """Var and Cov_1..Cov_kmax of consecutive block averages of filtered white noise."""
    f = np.linspace(0.0, fmax_fc * FC, npts)
    x = np.pi * f * D
    sc2 = np.ones_like(f)
    sc2[1:] = (np.sin(x[1:]) / x[1:]) ** 2
    w = H2(f) * sc2
    h = f[1] - f[0]
    var = 2 * S * _simpson(w, h)
    cov = np.array([2 * S * _simpson(w * np.cos(2 * np.pi * k * f * D), h)
                    for k in range(1, kmax + 1)])
    return var, cov

def block_var(psd, D, H2=None, fmax_fc=600.0, npts=2_000_001):
    """Block-average variance of a process of given PSD, with and without the filter."""
    f = np.linspace(0.0, fmax_fc * FC, npts)
    x = np.pi * f * D
    sc2 = np.ones_like(f)
    sc2[1:] = (np.sin(x[1:]) / x[1:]) ** 2
    y = psd(f) * sc2 * (H2(f) if H2 is not None else 1.0)
    return 2 * _simpson(y, f[1] - f[0])

# --------------------------------------------------------------------- report
def main():
    print("SANITY. |H(0)|^2 and |H(f_c)|^2 must read 1.000 and 0.500")
    for name, H2 in FILTERS:
        print(f"  {name:>9}: {float(H2(0.0)[0]):.6f}  {float(H2(FC)[0]):.6f}")
    print(f"  Gaussian equivalent width sigma_f*f_c = {SIG_F:.5f}"
          f"   (10-90 rise = {2*1.28155*SIG_F:.4f}/f_c)\n")

    print("1. THE DEPARTURE OF THE COMPOSITE KERNEL FROM THE UNIFORM WINDOW")
    print("   deficit = 1 - Var/(S/Delta); rho_k = Cov_k/Var; the last column is the")
    print("   conservation identity, which must read 1 exactly.")
    print(f"{'f_c*D':>6} {'filter':>9} | {'Var/(S/D)':>10} {'deficit':>9} {'rho_1':>9}"
          f" {'rho_2':>9} | {'Var+2sumCov':>12} | {'0.147/f_cD':>11} {'0.075/f_cD':>11}")
    print("-" * 112)
    store = {}
    for fcD in (0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        D = fcD / FC
        for name, H2 in FILTERS:
            var, cov = moments(H2, D)
            ideal = 1.0 / D
            deficit, rho1, rho2 = 1 - var / ideal, cov[0] / var, cov[1] / var
            print(f"{fcD:6.1f} {name:>9} | {var/ideal:10.4f} {deficit:9.4f} {rho1:9.4f}"
                  f" {rho2:9.4f} | {(var + 2*cov.sum())/ideal:12.6f} |"
                  f" {1.128*SIG_F/D:11.4f} {0.564*SIG_F/D:11.4f}")
            store[(fcD, name)] = (deficit, rho1, var, cov)
        print()

    print("2. WHAT IT COSTS A CONFIDENCE INTERVAL (4-pole Bessel).")
    print("   kappa is the effective-sample-size inflation the unmodelled correlation")
    print("   produces; the error bar scales as its square root.")
    print(f"{'f_c*D':>6} {'var deficit':>12} {'rho_1':>9} {'kappa':>9} {'CI error':>10}")
    for fcD in (0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        deficit, rho1, var, cov = store[(fcD, "Bessel-4")]
        kappa = 1 + 2 * cov.sum() / var
        print(f"{fcD:6.1f} {100*deficit:11.2f}% {100*rho1:8.2f}% {kappa:9.4f}"
              f" {100*(np.sqrt(max(kappa, 1e-12)) - 1):9.2f}%")
    print()

    print("3. THE GATING COMPONENT IS BETTER PROTECTED THAN THE INSTRUMENTAL ONE.")
    print("   Relative deficit in the block-average variance caused by the filter, for")
    print("   white noise and for an exponentially correlated component of time tau_c.")
    print(f"{'f_c*D':>6} {'white':>10} " + " ".join(f"{'tc/D=' + s:>10}"
                                                   for s in ("0.1", "0.3", "1", "3")))
    H2b = bessel_H2(4)
    for fcD in (1.0, 2.0, 5.0, 10.0):
        D = fcD / FC
        flat = lambda f: np.ones_like(f)
        row = [1 - block_var(flat, D, H2b) / block_var(flat, D)]
        for r in (0.1, 0.3, 1.0, 3.0):
            tc = r * D
            psd = lambda f, tc=tc: 2 * tc / (1 + (2 * np.pi * f * tc) ** 2)  # OU, R(0)=1
            row.append(1 - block_var(psd, D, H2b) / block_var(psd, D))
        print(f"{fcD:6.1f} " + " ".join(f"{100 * v:9.2f}%" for v in row))
    print()

    print("4. THE ONE PLACE THE EXPANSION DOES NOT APPLY: the interval holding a jump.")
    print("   Group delay and rise time are fixed times, not fractions of Delta.")
    print(f"{'poles':>6} {'group delay':>13} {'10-90 rise':>12} {'sum, in 1/f_c':>15}")
    for n in (2, 4, 6, 8):
        b, a = bessel_ba(n)
        w = np.linspace(1e-6, 2 * np.pi * 0.05 * FC, 400)
        _, h = signal.freqs(b, a, worN=w)
        gd = -np.gradient(np.unwrap(np.angle(h)), w)[0]
        t = np.linspace(0, 8 / FC, 400_001)
        _, y = signal.step((b, a), T=t)
        y = y / y[-1]
        rise = t[np.searchsorted(y, 0.90)] - t[np.searchsorted(y, 0.10)]
        print(f"{n:6d} {gd*FC:13.4f} {rise*FC:12.4f} {(gd + rise)*FC:15.4f}")


if __name__ == "__main__":
    main()
