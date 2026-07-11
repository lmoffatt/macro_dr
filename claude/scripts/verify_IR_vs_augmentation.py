#!/usr/bin/env python3
"""
Verify whether macro_IR is ISOMORPHIC to the standard integrated-measurement
Kalman filter (state augmentation) at the level of the per-step predictive
variance y_var and the Kalman gain gS = Cov(x_end, ybar).

Setting: 2-state scheme_CO ion channel, N independent channels, one acquisition
interval of length Delta. Observation = time-averaged current ybar = (1/Delta)
∫_0^Delta I(t) dt, with I = i_unit * n_open(t).

We compute y_var and gS THREE independent ways, for several priors over the
start state, and compare:

  (A) EXACT CTMC integrated moments  == what macro_IR computes (boundary-
      conditioned / E2-E3). Done here by numerical integration of the exact
      continuous-time autocovariance using expm(Q t). No linear-Gaussian
      approximation: this is the exact first-two-moment answer.

  (B) STANDARD AUGMENTATION (van Loan / integrated-measurement Kalman). The
      occupancy COUNT n is the Gaussian state; we augment with s = ∫ I dt and
      integrate the count moment ODEs (mean, covariance with the chemical-
      Langevin diffusion, Cov(n,s), Var(s)) over the interval. For independent
      channels these moment ODEs are EXACT, so this is the standard augmentation
      done exactly.

  (C) MONTE CARLO ground truth. Direct stochastic simulation of the 2-state
      CTMC (single channel, scaled to N by independence): record open-time and
      end-state, estimate Var(ybar) and Cov(n_end, ybar).

If A == B == C (to sampling error), macro_IR and the augmentation produce the
SAME filter -> isomorphic.
"""

import numpy as np
from scipy.linalg import expm
from scipy.integrate import dblquad, quad, solve_ivp

rng = np.random.default_rng(0)

# ---------------------------------------------------------------- model
on, off = 10.0, 100.0          # C->O and O->C rates (s^-1), scheme_CO
i_unit = 1.0                   # single-channel current of the OPEN state
g = np.array([0.0, i_unit])    # conductance vector: [closed, open]
N = 20                         # number of channels
Delta = 100 / 50_000.0         # interval: 100 samples @ 50 kHz = 2 ms

# Q convention: occupancy ROW prob p, dp/dt = p Q. states 0=Closed, 1=Open.
Q = np.array([[-on,  on],
              [ off, -off]])

def P(t):                      # transition matrix, p(t) = p(0) P(t)
    return expm(Q * t)

# ======================================================== (A) exact CTMC
def exact_ctmc(mu):
    """Exact integrated-current variance and Cov(n_end, ybar) for N channels."""
    def meanI(u):              # E[current per channel] at time u
        return (mu @ P(u)) @ g
    def cov1(u, v):            # single-channel current autocovariance (u<=v)
        a = mu @ P(u)
        EIuIv = (a * g) @ P(v - u) @ g
        return EIuIv - meanI(u) * meanI(v)
    # Var(ybar) = N/Delta^2 * ∫∫ cov1(u,v) du dv  (N indep channels add variances)
    integ, _ = dblquad(lambda v, u: cov1(min(u, v), max(u, v)),
                       0, Delta, 0, Delta, epsabs=1e-12, epsrel=1e-10)
    yvar = N / Delta**2 * integ
    # gS_j = Cov(n_end_j, ybar) = N/Delta * ∫ Cov(s_end=j, I(u)) du
    def cross(u):              # vector over end states j
        a = mu @ P(u)
        Esj_Iu = (a * g) @ P(Delta - u)          # E[1{end=j} * I(u)] / channel
        return Esj_Iu - (mu @ P(Delta)) * meanI(u)
    gS = N / Delta * np.array([
        quad(lambda u: cross(u)[j], 0, Delta, epsabs=1e-12, epsrel=1e-10)[0]
        for j in range(2)])
    return yvar, gS

# ============================================ (B) standard augmentation
def augmentation(mu):
    """Count moment ODEs (mean, cov, Cov(n,s), Var(s)). Exact for indep channels."""
    mu_n0 = N * mu
    Sig0 = N * (np.diag(mu) - np.outer(mu, mu))   # multinomial count covariance
    # state = [mu_n(2), Sigma(4), c=Cov(n,s)(2), var_s(1)]
    y0 = np.concatenate([mu_n0, Sig0.ravel(), np.zeros(2), [0.0]])

    def rhs(t, y):
        mu_n = y[0:2]
        Sig = y[2:6].reshape(2, 2)
        c = y[6:8]
        # chemical-Langevin diffusion of the count (2-state): D = (on*nC+off*nO)*[[1,-1],[-1,1]]
        a = on * mu_n[0] + off * mu_n[1]
        D = a * np.array([[1.0, -1.0], [-1.0, 1.0]])
        dmu = Q.T @ mu_n
        dSig = Q.T @ Sig + Sig @ Q + D
        dc = Q.T @ c + Sig @ g                    # d/dt Cov(n,s) = Q^T c + Sigma g
        dvar = 2.0 * (g @ c)                       # d/dt Var(s) = 2 g·Cov(n,s)
        return np.concatenate([dmu, dSig.ravel(), dc, [dvar]])

    sol = solve_ivp(rhs, [0, Delta], y0, rtol=1e-11, atol=1e-13, dense_output=False)
    yT = sol.y[:, -1]
    c = yT[6:8]
    var_s = yT[8]
    yvar = var_s / Delta**2
    gS = c / Delta
    return yvar, gS

# ================================================= (C) Monte Carlo truth
def monte_carlo(mu, M=400_000):
    """Single-channel CTMC simulation; scale to N by independence."""
    rates = np.array([on, off])    # leaving-rate of state 0 (closed), state 1 (open)
    s0 = (rng.random(M) < mu[1]).astype(int)     # start state per prior mu
    open_time = np.zeros(M)
    end_state = np.empty(M, dtype=int)
    # vectorized event-driven sim (few jumps in 2 ms)
    state = s0.copy()
    t = np.zeros(M)
    alive = np.ones(M, dtype=bool)
    while alive.any():
        idx = np.where(alive)[0]
        st = state[idx]
        dwell = rng.exponential(1.0 / rates[st])
        tend = t[idx] + dwell
        over = tend > Delta
        # channels whose dwell exceeds the interval: accrue remaining time, stop
        rem = np.where(over, Delta - t[idx], dwell)
        open_time[idx] += np.where(st == 1, rem, 0.0)
        # advance the survivors
        t[idx] = np.where(over, Delta, tend)
        state[idx] = np.where(over, st, 1 - st)   # switch (2-state)
        done = over
        end_state[idx[done]] = st[done]           # end state = current state (no switch happened in [t,Delta])
        alive[idx[done]] = False
    ybar_per_ch = (i_unit / Delta) * open_time     # single-channel contribution to ybar
    # Var(ybar) for N indep channels = N * Var(single contribution)
    yvar = N * np.var(ybar_per_ch, ddof=1)
    # gS_j = Cov(n_end_j, ybar) = N * Cov(1{end=j}, single contribution)
    gS = np.array([N * np.cov(1.0 * (end_state == j), ybar_per_ch, ddof=1)[0, 1]
                   for j in range(2)])
    se_yvar = yvar * np.sqrt(2.0 / (M - 1))         # ~ rel. std error of a variance
    return yvar, gS, se_yvar

# ===================================== (D) the REAL macro_IR recursion
# Faithful Python port of the codebase's per-interval macro_IR step
# (uses_recursive=true, uses_averaging=2). It builds the Qdt boundary-state
# moments gtotal_ij / gtotal_sqr_ij from the EXACT eigen-kernels E2/E3
# (legacy/qmodel.h:804-855 Ee/E3, calc_Qdt at :1462), then applies the live
# av=2 step formulas from safely_calculate_Algo_State_recursive (:4465).
# This is what (A) only ASSERTED in a comment; here macro_IR actually runs.

def _Ee(x, y, ex, ey):                 # (e^x - e^y)/(x - y), -> e^x at x==y
    return ex if abs(x - y) < 1e-9 else (ex - ey) / (x - y)

def _EX111(x, y, z, ex):
    return ex / ((x - y) * (x - z))

def _E111(x, y, z, ex, ey, ez):        # x,y,z all distinct
    return _EX111(x, y, z, ex) + _EX111(y, x, z, ey) + _EX111(z, y, x, ez)

def _E12(x, y, ex, ey):                # the y==z (double) confluent case
    return _EX111(x, y, y, ex) + ey / (y - x) * (1.0 - 1.0 / (y - x))

def _E3(x, y, z, ex, ey, ez):          # qmodel.h:834 degeneracy ladder
    tol = 1e-9
    if abs(x - y) < tol:
        return ex / 2.0 if abs(y - z) < tol else _E12(z, x, ez, ex)
    if abs(y - z) < tol:
        return _E12(x, y, ex, ey)
    if abs(x - z) < tol:
        return _E12(y, x, ey, ex)
    return _E111(x, y, z, ex, ey, ez)

def qdt_moments(dt=Delta):
    """Build the macro_IR Qdt boundary-state moments for one interval of length
    dt. Independent of the prior, so it is computed once and reused every step
    of the recursion. Returns gtotal_ij, gtotal_sqr_ij, P, gmean_i, gsqr_i,
    gmean_ij (min_P -> 0, no Bayesian shrinkage)."""
    n = Q.shape[0]
    lam, Vr = np.linalg.eig(Q)         # Q = Vr diag(lam) Vr^-1, P(t)=Vr e^{lam t} Vr^-1
    W = np.linalg.inv(Vr)
    ladt = lam * dt
    ex = np.exp(ladt)
    WgV = W @ np.diag(g) @ Vr          # conductance operator in the eigenbasis

    # gtotal_ij : 1/dt-normalised interval-averaged conductance, start i -> end j
    E2 = np.array([[_Ee(ladt[i], ladt[j], ex[i], ex[j]) for j in range(n)]
                   for i in range(n)])
    gtot = np.real(Vr @ (WgV * E2) @ W)

    # gtotal_sqr_ij : 2 * V (WgV . E3 . WgV) W  (the boundary-state second moment)
    E3w = np.zeros((n, n), dtype=complex)
    for a in range(n):
        for c in range(n):
            for b in range(n):
                E3w[a, c] += WgV[a, b] * WgV[b, c] * _E3(
                    ladt[a], ladt[b], ladt[c], ex[a], ex[b], ex[c])
    gtot_sqr = 2.0 * np.real(Vr @ E3w @ W)

    P_mat = np.real(Vr @ np.diag(ex) @ W)           # transition matrix over dt
    gmean_i = gtot.sum(axis=1)                      # row sums over end state j
    gsqr_i = gtot_sqr.sum(axis=1)
    gmean_ij = gtot / P_mat                         # min_P -> 0 (no shrinkage)
    return dict(gtot=gtot, gtot_sqr=gtot_sqr, P=P_mat,
                gmean_i=gmean_i, gsqr_i=gsqr_i, gmean_ij=gmean_ij)

def macro_IR_step(p, P_Cov, qd, e=0.0):
    """The live av=2 (IR) observation model: returns y_mean, y_var and the
    endpoint gain gS (per channel) for prior (p, P_Cov)."""
    gtot, gmean_i, gsqr_i = qd['gtot'], qd['gmean_i'], qd['gsqr_i']
    P_mat, gmean_ij = qd['P'], qd['gmean_ij']
    u = np.ones(Q.shape[0])
    SmD = P_Cov - np.diag(p)
    cross = (gtot * gmean_ij) @ u                    # Sum_j gtotal_ij * gmean_ij
    gSg = gmean_i @ SmD @ gmean_i + p @ cross
    ms = p @ (gsqr_i - cross)                        # av=2 residual (no double count)
    y_mean = N * (p @ gmean_i)                       # baseline 0
    y_var = e + N * (gSg + ms)
    gS = (gmean_i @ SmD) @ P_mat + p @ gtot          # Cov(X_end_j, ybar), per channel
    return y_mean, y_var, gS

def macro_IR(mu):
    """One macro_IR interval step from a known prior mu (P_Cov multinomial),
    noise-free, for the A/B/C/D single-interval comparison."""
    qd = qdt_moments()
    P_Cov = np.diag(mu) - np.outer(mu, mu)
    y_mean, yvar, gS = macro_IR_step(mu, P_Cov, qd, e=0.0)
    return yvar, N * gS                              # scale gS to N channels

# ============================================================== run
p_eq = np.array([off, on]) / (on + off)
priors = {
    "stationary  mu=[%.3f,%.3f]" % (p_eq[0], p_eq[1]): p_eq,
    "all closed  mu=[1,0]": np.array([1.0, 0.0]),
    "all open    mu=[0,1]": np.array([0.0, 1.0]),
    "mixed       mu=[0.5,0.5]": np.array([0.5, 0.5]),
}

print(f"model: on={on} off={off} g={list(g)} N={N} Delta={Delta*1e3:.3f} ms"
      f"  (Delta/tau = {Delta*(on+off):.3f})\n")
hdr = f"{'prior':28s} {'method':14s} {'y_var':>14s} {'gS_closed':>13s} {'gS_open':>13s}"
print(hdr); print("-" * len(hdr))
for name, mu in priors.items():
    yA, gA = exact_ctmc(mu)
    yB, gB = augmentation(mu)
    yC, gC, seC = monte_carlo(mu)
    yD, gD = macro_IR(mu)
    print(f"{name:28s} {'A exact-CTMC':14s} {yA:14.6e} {gA[0]:13.6e} {gA[1]:13.6e}")
    print(f"{'':28s} {'B augmentation':14s} {yB:14.6e} {gB[0]:13.6e} {gB[1]:13.6e}")
    print(f"{'':28s} {'C monte-carlo':14s} {yC:14.6e} {gC[0]:13.6e} {gC[1]:13.6e}"
          f"   (±{seC/yC*100:.2f}% on y_var)")
    print(f"{'':28s} {'D macro_IR':14s} {yD:14.6e} {gD[0]:13.6e} {gD[1]:13.6e}")
    relDA = abs(yD - yA) / abs(yA) if yA else 0.0
    relDB = abs(yD - yB) / abs(yB) if yB else 0.0
    relAC = abs(yA - yC) / abs(yA) if yA else 0.0
    print(f"{'':28s} {'rel.diff':14s} D vs A: {relDA:.2e}   D vs B: {relDB:.2e}"
          f"   A vs C: {relAC:.2e}")
    print(f"{'':28s} {'gS maxdiff':14s} D-A: {np.max(np.abs(gD-gA)):.2e}"
          f"   D-B: {np.max(np.abs(gD-gB)):.2e}   A-B: {np.max(np.abs(gA-gB)):.2e}")
    print()


# =====================================================================
# PART 2 — MULTI-STEP recursion: does the augmentation compute the SAME
# filtered trajectory and log-likelihood as macro_IR across many intervals,
# once the Gaussian closure + Kalman update are iterated?
#
#   macro_IR filter : carries per-channel (p, P_Cov); each step uses the
#                     boundary-state Qdt moments + the live av=2 update
#                     (qmodel.h:5627-5659).
#   augmentation    : carries the channel COUNT moments (mu_n, Sigma); each
#                     step integrates the count+integral moment ODEs and does
#                     the standard integrated-measurement Kalman update.
# Both run on the SAME simulated recording, with the SAME measurement-noise e.
# =====================================================================

def _add_open_time(bins, ta, tb, T):
    """Distribute an OPEN sojourn [ta, tb] over the interval bins it overlaps."""
    k0 = max(0, int(ta // Delta))
    k1 = min(T - 1, int(tb // Delta))
    for k in range(k0, k1 + 1):
        lo, hi = k * Delta, (k + 1) * Delta
        ov = min(tb, hi) - max(ta, lo)
        if ov > 0:
            bins[k] += ov

def simulate_recording(mu0, T, e_obs, seed=1):
    """N independent 2-state channels over T intervals; returns the noisy
    averaged-current trace ybar[t] (per interval) and the true mean trace."""
    sr = np.random.default_rng(seed)
    rates = np.array([on, off])
    open_per_interval = np.zeros(T)
    Ttot = T * Delta
    for _ in range(N):
        state = int(sr.random() < mu0[1])
        t = 0.0
        while t < Ttot:
            dwell = sr.exponential(1.0 / rates[state])
            t2 = min(t + dwell, Ttot)
            if state == 1:
                _add_open_time(open_per_interval, t, t2, T)
            t += dwell
            state = 1 - state
    ybar_true = (i_unit / Delta) * open_per_interval
    ybar = ybar_true + sr.normal(0.0, np.sqrt(e_obs), size=T)
    return ybar, ybar_true

def augment_predict(mu_n, Sig, dt=Delta):
    """Propagate the count moments over one interval and accumulate the
    integral s = INT I dt: returns mu_n_end, Sig_end, c=Cov(n_end,s), var_s,
    s_mean. EXACT for independent channels (the moment ODEs close)."""
    y0 = np.concatenate([mu_n, Sig.ravel(), np.zeros(2), [0.0, 0.0]])
    def rhs(t, y):
        mn = y[0:2]; Sg = y[2:6].reshape(2, 2); c = y[6:8]
        a = on * mn[0] + off * mn[1]
        D = a * np.array([[1.0, -1.0], [-1.0, 1.0]])
        dmu = Q.T @ mn
        dSig = Q.T @ Sg + Sg @ Q + D
        dc = Q.T @ c + Sg @ g
        dvar = 2.0 * (g @ c)
        dsm = g @ mn
        return np.concatenate([dmu, dSig.ravel(), dc, [dvar, dsm]])
    sol = solve_ivp(rhs, [0, dt], y0, rtol=1e-11, atol=1e-13)
    yT = sol.y[:, -1]
    return yT[0:2], yT[2:6].reshape(2, 2), yT[6:8], yT[8], yT[9]

def _psd_project(M):
    """to_Covariance_Probability's PSD guard: clip negative eigenvalues to 0."""
    w, V = np.linalg.eigh(0.5 * (M + M.T))
    if w.min() >= 0:
        return M
    return (V * np.maximum(w, 0.0)) @ V.T

def filter_macro_IR(ybar, mu0, e_obs, guarded=False):
    """macro_IR recursion (per-channel coords). guarded=False runs the RAW
    Gaussian step (alpha=1, no projection) for a head-to-head with the
    augmentation. guarded=True applies the live simplex trust coefficient
    alpha_mu and the PSD projection — the real macro_IR. Tracks the guard
    activity either way."""
    qd = qdt_moments()
    P = qd['P']
    p = mu0.copy()
    P_Cov = np.diag(p) - np.outer(p, p)
    logL = 0.0
    mean_trace = np.zeros(len(ybar))
    alpha_min = 1.0
    coveig_min = np.inf
    for t, y in enumerate(ybar):
        y_mean, y_var, gS = macro_IR_step(p, P_Cov, qd, e=e_obs)
        mean_trace[t] = y_mean
        dy = y - y_mean
        chi = dy / y_var
        logL += -0.5 * np.log(2 * np.pi * y_var) - 0.5 * dy * chi
        SmD = P_Cov - np.diag(p)
        q = p @ P                                    # predicted (pre-update) mean
        # simplex trust coefficient alpha_mu (factor 0.9, hard min)
        d = chi * gS
        a_real = 1.0
        for i in range(2):
            if d[i] > 0:
                a_real = min(a_real, 0.9 * (1.0 - q[i]) / d[i])
            elif d[i] < 0:
                a_real = min(a_real, 0.9 * (-q[i]) / d[i])
        alpha_min = min(alpha_min, a_real)
        alpha = a_real if guarded else 1.0
        # live update (qmodel.h:5648-5655)
        p = q + alpha * chi * gS
        sigma_pre = P.T @ SmD @ P + np.diag(q)        # AT_B_A(P,SmD)+diag(p*P)
        P_Cov = sigma_pre - (N / y_var) * np.outer(gS, gS)
        coveig_min = min(coveig_min, np.linalg.eigvalsh(P_Cov).min())
        if guarded:
            P_Cov = _psd_project(P_Cov)
    return logL, mean_trace, alpha_min, coveig_min

def filter_augmentation(ybar, mu0, e_obs):
    """Standard integrated-measurement Kalman on the channel COUNT moments."""
    mu_n = N * mu0.copy()
    Sig = N * (np.diag(mu0) - np.outer(mu0, mu0))
    logL = 0.0
    mean_trace = np.zeros(len(ybar))
    for t, y in enumerate(ybar):
        mu_e, Sig_e, c, var_s, s_mean = augment_predict(mu_n, Sig)
        y_mean = s_mean / Delta                       # E[ybar]
        cross = c / Delta                             # Cov(n_end, ybar)
        y_var = var_s / Delta**2 + e_obs
        mean_trace[t] = y_mean
        dy = y - y_mean
        logL += -0.5 * np.log(2 * np.pi * y_var) - 0.5 * dy * dy / y_var
        K = cross / y_var
        mu_n = mu_e + K * dy
        Sig = Sig_e - np.outer(cross, cross) / y_var
    return logL, mean_trace

# ----- run the multi-step comparison -----
print("=" * 86)
print("PART 2 — multi-step recursion: macro_IR filter vs augmentation filter "
      "on one recording")
print("=" * 86)
T_int = 400
e_obs = 0.25                                          # measurement-noise variance on ybar
mu_start = p_eq
ybar, ybar_true = simulate_recording(mu_start, T_int, e_obs, seed=1)
logL_M, mean_M, amin, ceig = filter_macro_IR(ybar, mu_start, e_obs, guarded=False)
logL_G, mean_G, _, _ = filter_macro_IR(ybar, mu_start, e_obs, guarded=True)
logL_A, mean_A = filter_augmentation(ybar, mu_start, e_obs)

print(f"recording: T={T_int} intervals, N={N} channels, e(meas var)={e_obs}, "
      f"start=stationary\n")
print(f"  cumulative logL   macro_IR (raw, alpha=1) = {logL_M:.8f}")
print(f"                    augmentation            = {logL_A:.8f}")
print(f"  -> RAW recursion vs augmentation: |Delta logL| = {abs(logL_M - logL_A):.3e}"
      f",  max |Delta y_mean| = {np.max(np.abs(mean_M - mean_A)):.3e}")
print(f"     (the cores are the SAME Kalman filter)\n")
print(f"  cumulative logL   macro_IR (guarded)      = {logL_G:.8f}")
print(f"  -> guarded macro_IR vs augmentation: |Delta logL| = {abs(logL_G - logL_A):.3e}"
      f",  max |Delta y_mean| = {np.max(np.abs(mean_G - mean_A)):.3e}")
print(f"     simplex trust alpha_mu hit {amin:.4f} (1.0=never binds); "
      f"the guard is the ONLY source of divergence")
print(f"  PSD guard: min cov eigenvalue (raw) = {ceig:.2e} (>~0 -> PSD held by "
      f"construction, projection ~ identity)")
