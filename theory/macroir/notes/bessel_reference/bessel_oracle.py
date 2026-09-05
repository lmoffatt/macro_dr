#!/usr/bin/env python3
# bessel_oracle.py — independent reference implementation of the MacroIR-Bessel
# augmented recursion, validated against exact Monte Carlo of the filtered
# process. Plays the role mr_vs_ir_boxes.py played for the MR/IR derivation:
# the formulas the C++ hot paths must reproduce, checked before any C++ exists.
#
# 2026-08-31. Plan: theory/macroir/notes/macroir_bessel_plan.md (M0/M2 math).
# Companion C++ kernels: legacy/acquisition_filter.h (phase-0 gates green).
#
# Design: everything is generic over a REAL linear observation system
# (A, b, read) driven by the ensemble current i(t) = Σ_c γ_{X_c(t)} + ξ(t):
#   * native read  (K+4): A = Bessel modal 2×2 blocks, read = filter output
#   * grouped read (K+5): same + one integrator row (the cumulator u,
#     legitimately reset at each window start), read = u/Δ
#   * box          (K+1): integrator alone on i(t), read = u/Δ — this IS the
#     current MacroIR member and is the regression anchor (gate O2)
# One code path covers the three; the collapse the plan demands is literal.
#
# Gates:
#   O1 sim-noise:   simulator's filtered-noise reads vs closed forms
#                   (same formulas the C++ gates validated vs quadrature).
#   O2 box anchor:  generic ensemble assembly vs the boundary-state (Σ^bnd)
#                   formulas of mr_vs_ir_from_macror.md, two windows deep.
#   O3 native MC:   window-1 predictive (ẑ, v) vs Monte Carlo.
#   O4 whiteness:   standardized residuals of the Bessel member ≈ white,
#                   while the box member on the same filtered data shows the
#                   positive lag-1 the plan predicts. Motivation numbers.
#   O5 grouped MC:  the K+5 read, window-1 predictive + whiteness.

import numpy as np

try:
    from scipy.linalg import expm as _expm
except Exception:  # scaling-and-squaring fallback, small matrices only

    def _expm(M):
        M = np.asarray(M, dtype=float)
        n = int(np.ceil(np.log2(max(1.0, np.linalg.norm(M, np.inf))))) + 1
        A = M / (2.0**n)
        E = np.eye(M.shape[0])
        P = np.eye(M.shape[0])
        for k in range(1, 19):
            P = P @ A / k
            E = E + P
        for _ in range(n):
            E = E @ E
        return E


rng = np.random.default_rng(20260831)

# ----------------------------------------------------------------------------
# Bessel filter (same construction as acquisition_filter.h)
# ----------------------------------------------------------------------------
def reverse_bessel_coeffs(n):
    prev, cur = [1.0], [1.0, 1.0]
    for k in range(2, n + 1):
        nxt = [0.0] * (k + 1)
        for i, c in enumerate(cur):
            nxt[i] += (2 * k - 1) * c
        for i, c in enumerate(prev):
            nxt[i + 2] += c
        prev, cur = cur, nxt
    return np.array(cur if n >= 1 else prev)


def make_bessel(order, fc_hz):
    a = reverse_bessel_coeffs(order)
    th0 = a[0]
    lo, hi = 0.0, 4.0 * order
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        m2 = abs(np.polyval(a[::-1], 1j * mid)) ** 2
        lo, hi = (mid, hi) if m2 < 2 * th0**2 else (lo, mid)
    w3 = 0.5 * (lo + hi)
    scale = 2 * np.pi * fc_hz / w3
    poles = np.roots(a[::-1]) * scale
    A = np.prod(-poles)
    residues = np.array(
        [A / np.prod([poles[k] - p for j, p in enumerate(poles) if j != k]) for k in range(len(poles))]
    )
    return poles, residues


def modal_real_system(poles, residues):
    """Real block-modal (u,v)=(2Re x, 2Im x) per conjugate pair; output = Σ u."""
    pairs = [(p, r) for p, r in zip(poles, residues) if p.imag > 0]
    m = 2 * len(pairs)
    A = np.zeros((m, m))
    b = np.zeros(m)
    c = np.zeros(m)
    for k, (p, r) in enumerate(pairs):
        s, w = -p.real, p.imag
        A[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = [[-s, -w], [w, -s]]
        b[2 * k] = 2 * r.real
        b[2 * k + 1] = 2 * r.imag
        c[2 * k] = 1.0
    return A, b, c


def add_cumulator(A, b, c):
    """Append du/dt = c·x (the K+5 row); reset per window is the caller's job."""
    m = A.shape[0]
    A2 = np.zeros((m + 1, m + 1))
    A2[:m, :m] = A
    A2[m, :m] = c
    b2 = np.concatenate([b, [0.0]])
    return A2, b2


# ----------------------------------------------------------------------------
# Exact propagators over a piece of length dt (Van Loan block tricks):
# Phi = e^{A dt};  psi = ∫ e^{As}b ds  (input loading, singular-A safe);
# Sig = S0 ∫ e^{As} b bᵀ e^{Aᵀs} ds  (noise Gramian).
# ----------------------------------------------------------------------------
def piece_ops(A, b, S0, dt, _cache={}):
    key = (A.tobytes(), b.tobytes(), S0, round(dt, 15))
    if key in _cache:
        return _cache[key]
    m = A.shape[0]
    M1 = np.zeros((m + 1, m + 1))
    M1[:m, :m] = A * dt
    M1[:m, m] = b * dt
    E1 = _expm(M1)
    Phi, psi = E1[:m, :m], E1[:m, m]
    M2 = np.zeros((2 * m, 2 * m))
    M2[:m, :m] = -A * dt
    M2[:m, m:] = np.outer(b, b) * dt
    M2[m:, m:] = A.T * dt
    E2 = _expm(M2)
    Sig = S0 * (E2[m:, m:].T @ E2[:m, m:])
    Sig = 0.5 * (Sig + Sig.T)
    w, V = np.linalg.eigh(Sig)
    L = V @ np.diag(np.sqrt(np.clip(w, 0.0, None)))
    _cache[key] = (Phi, psi, L)
    return Phi, psi, L


# ----------------------------------------------------------------------------
# Exact simulator: 2-state channels (C,O), γ = (0,1); obs system (A,b);
# native reads read x at window ends; grouped/box reset marked rows first.
# ----------------------------------------------------------------------------
def simulate(nsim, Nch, kco, koc, A, b, read, S0, dt_win, Twin, reset_rows=()):
    m = A.shape[0]
    p_open = kco / (kco + koc)
    Ttot = dt_win * Twin
    Z = np.empty((nsim, Twin))
    for s in range(nsim):
        # channel paths: jump times + state values per channel
        events = []  # (time, dG)
        G0 = 0.0
        for _ in range(Nch):
            st = 1 if rng.random() < p_open else 0
            G0 += st
            t = 0.0
            while True:
                rate = koc if st == 1 else kco
                t += rng.exponential(1.0 / rate)
                if t >= Ttot:
                    break
                events.append((t, -1.0 if st == 1 else 1.0))
                st = 1 - st
        events.sort()
        x = np.zeros(m)
        G = G0
        t0, ev = 0.0, 0
        for j in range(Twin):
            tend = (j + 1) * dt_win
            for r in reset_rows:
                x[r] = 0.0
            while ev < len(events) and events[ev][0] < tend:
                te, dG = events[ev]
                if te > t0:
                    Phi, psi, L = piece_ops(A, b, S0, te - t0)
                    x = Phi @ x + psi * G + L @ rng.standard_normal(m)
                    t0 = te
                G += dG
                ev += 1
            if tend > t0:
                Phi, psi, L = piece_ops(A, b, S0, tend - t0)
                x = Phi @ x + psi * G + L @ rng.standard_normal(m)
                t0 = tend
            Z[s, j] = read @ x
    return Z


# ----------------------------------------------------------------------------
# Window quantities by quadrature (the objects calc_Qdtf_eig will produce
# spectrally): P, pair-resolved firsts F1[i,j] (m-vectors), start-conditioned
# seconds S2[i] (m×m), old-state carry Phi, noise Gramian; γ general.
# ----------------------------------------------------------------------------
def window_quantities(Q, gamma, A, b, S0, dt, M=200):
    K = Q.shape[0]
    m = A.shape[0]
    h = dt / M
    s_nodes = np.linspace(0.0, dt, M + 1)
    Ps = np.array([_expm(Q * s) for s in s_nodes])          # P(s)
    Gs = np.array([_expm(A * (dt - s)) @ b for s in s_nodes])  # e^{A(Δ−s)}b
    wS = np.ones(M + 1)
    wS[1:-1:2], wS[2:-1:2] = 4.0, 2.0
    wS *= h / 3.0  # Simpson

    F1 = np.zeros((K, K, m))
    for i in range(K):
        for j in range(K):
            F1[i, j] = sum(
                wS[a] * sum(Ps[a, i, k] * gamma[k] * Ps[M - a, k, j] for k in range(K)) * Gs[a]
                for a in range(M + 1)
            )
    f1 = F1.sum(axis=1)  # start-conditioned firsts

    # seconds, trapezoid on the (s<t) triangle
    S2 = np.zeros((K, m, m))
    wT = np.ones(M + 1)
    wT[0] = wT[-1] = 0.5
    wT *= h
    for i in range(K):
        acc = np.zeros((m, m))
        for a in range(M + 1):
            row = np.zeros(m)
            pa = Ps[a, i, :]
            for bb in range(a, M + 1):
                lag = Ps[bb - a]
                val = sum(pa[k] * gamma[k] * sum(lag[k, l] * gamma[l] for l in range(K)) for k in range(K))
                acc += wT[a] * wT[bb] * val * (np.outer(Gs[a], Gs[bb]) + np.outer(Gs[bb], Gs[a])) * (
                    0.5 if bb == a else 1.0
                )
        S2[i] = acc
    Phi, _, _ = piece_ops(A, b, S0, dt)
    # noise Gramian (recompute; piece_ops caches chol only)
    m2 = np.zeros((2 * m, 2 * m))
    m2[:m, :m] = -A * dt
    m2[:m, m:] = np.outer(b, b) * dt
    m2[m:, m:] = A.T * dt
    E2 = _expm(m2)
    Sig = S0 * (E2[m:, m:].T @ E2[:m, m:])
    return dict(P=_expm(Q * dt), F1=F1, f1=f1, S2=S2, Phi=Phi, Sig=0.5 * (Sig + Sig.T))


# ----------------------------------------------------------------------------
# The augmented recursion (count scale). State: mu_N (K), S_NN (K×K),
# m_x (m), C_Nx (K×m), S_xx (m×m).
# ----------------------------------------------------------------------------
def predict(state, wq, reset_rows=()):
    mu, S, mx, C, X = state
    P, F1, f1, S2, Phi, Sig = wq["P"], wq["F1"], wq["f1"], wq["S2"], wq["Phi"], wq["Sig"]
    K, m = C.shape
    for r in reset_rows:  # legitimate recorder reset of the cumulator
        mx = mx.copy(); mx[r] = 0.0
        C = C.copy(); C[:, r] = 0.0
        X = X.copy(); X[r, :] = 0.0; X[:, r] = 0.0
    mu2 = P.T @ mu
    SmD = S - np.diag(mu)
    S2N = P.T @ SmD @ P + np.diag(mu2)                       # sigma_pre
    mx2 = Phi @ mx + np.einsum("i,im->m", mu, f1)
    # C_Nx: old-state carry + within-channel + between-channel (the gS shape)
    C2 = P.T @ C @ Phi.T
    for j in range(K):
        C2[j] += np.einsum("i,im->m", mu, F1[:, j, :]) - np.einsum("i,i,im->m", mu, P[:, j], f1)
        C2[j] += np.einsum("i,ik,km->m", P[:, j], S, f1)
    # S_xx
    M1 = np.einsum("ia,ib->ab", C, f1)                       # Cov(x0, Σw)
    Vw = np.einsum("i,iab->ab", mu, S2 - np.einsum("ia,ib->iab", f1, f1))
    Vw += np.einsum("ik,ia,kb->ab", S, f1, f1)
    X2 = Phi @ X @ Phi.T + Phi @ M1 + (Phi @ M1).T + Vw + Sig
    return (mu2, S2N, mx2, C2, X2)


def update(state, read, z, floor=0.0):
    mu, S, mx, C, X = state
    zhat = read @ mx
    v = read @ X @ read + floor
    kN = C @ read
    kx = X @ read
    d = z - zhat
    mu2 = mu + kN * d / v
    mx2 = mx + kx * d / v
    S2 = S - np.outer(kN, kN) / v
    C2 = C - np.outer(kN, kx) / v
    X2 = X - np.outer(kx, kx) / v
    ll = -0.5 * (np.log(2 * np.pi * v) + d * d / v)
    return (mu2, S2, mx2, C2, X2), zhat, v, ll


def run_member(Z, wq, read, init, reset_rows=()):
    nsim, Twin = Z.shape
    zh = np.empty_like(Z)
    vv = np.empty_like(Z)
    for s in range(nsim):
        st = init()
        for j in range(Twin):
            st = predict(st, wq, reset_rows)
            st, zh[s, j], vv[s, j], _ = update(st, read, Z[s, j])
    return zh, vv


# ----------------------------------------------------------------------------
# Boundary-state reference for the box anchor (mr_vs_ir_from_macror.md, count
# scale), pair-resolved Γ̄ and V̄ by quadrature.
# ----------------------------------------------------------------------------
def box_reference_blocks(wq, Q, gamma, S0, dt, mu, S, M=200):
    """Boundary-state (Σ^bnd) assembly of the three blocks for the box read,
    count scale, following mr_vs_ir_from_macror.md. Uses wq's F1 and the SAME
    quadrature grid/weights for the pair-resolved seconds, so any quadrature
    error cancels against the generic member and the gate tests pure algebra.

    Count-scale boundary covariance (verified against the multinomial law):
      Cov(B)_{(ij),(i'j')} = P_ij (S − diag μ)_{ii'} P_{i'j'}
                             + δ_{ii'} δ_{jj'} μ_i P_ij  (regenerated diagonal)
    """
    K = Q.shape[0]
    P, F1 = wq["P"], wq["F1"]
    h = dt / M
    s_nodes = np.linspace(0.0, dt, M + 1)
    Ps = np.array([_expm(Q * s) for s in s_nodes])
    wT = np.ones(M + 1)
    wT[0] = wT[-1] = 0.5
    wT *= h
    S2p = np.zeros((K, K))  # E[(∫γ)² · 1{end=j} | start=i], pair-resolved
    for i in range(K):
        for j in range(K):
            acc = 0.0
            for a in range(M + 1):
                for bb in range(a, M + 1):
                    val = 0.0
                    for k in range(K):
                        if gamma[k] == 0.0:
                            continue
                        for l in range(K):
                            if gamma[l] == 0.0:
                                continue
                            val += Ps[a, i, k] * gamma[k] * Ps[bb - a, k, l] * gamma[l] * Ps[M - bb, l, j]
                    acc += wT[a] * wT[bb] * 2.0 * val * (0.5 if bb == a else 1.0)
            S2p[i, j] = acc
    Gbar = np.where(P > 1e-300, F1[:, :, 0] / np.maximum(P, 1e-300) / dt, 0.0)
    Vbar = np.where(P > 1e-300, S2p / np.maximum(P, 1e-300) / dt**2 - Gbar**2, 0.0)
    mub = mu[:, None] * P
    SmD = S - np.diag(mu)
    yhat = np.sum(mub * Gbar)
    # Γ̄ᵀ Cov(B) Γ̄ by the contracted route (never forming K²×K²)
    var_state = np.einsum("ij,ik,kl,kl->", Gbar * P, SmD, P, Gbar) + np.sum(mub * Gbar**2)
    v = var_state + np.sum(mub * Vbar) + S0 / dt
    # cross block at the end frame: Cov(N_j(Δ), z)
    gbar0 = (P * Gbar).sum(axis=1)
    kappa = np.einsum("ij,ik,k->j", P, SmD, gbar0) + np.einsum("ij,ij,i->j", P, Gbar, mu)
    return yhat, v, kappa


# ----------------------------------------------------------------------------
# Closed forms for filtered-noise reads (ported from acquisition_filter.h)
# ----------------------------------------------------------------------------
def noise_read_var_cov(poles, residues, S0, dt, mmax=1):
    var = 0.0 + 0.0j
    covs = []
    for k, rk in enumerate(residues):
        for l, rl in enumerate(residues):
            p = poles[l]
            var += rk * rl / (-(poles[k] + p)) * (np.exp(p * dt) - 1 - p * dt) / p**2
    var = (2.0 * S0 / dt**2) * var.real
    for mm in range(1, mmax + 1):
        c = 0.0 + 0.0j
        for k, rk in enumerate(residues):
            for l, rl in enumerate(residues):
                p = poles[l]
                e1 = np.exp(p * dt) - 1
                c += rk * rl / (-(poles[k] + p)) * np.exp(p * dt * (mm - 1)) * e1 * e1 / p**2
        covs.append((S0 / dt**2) * c.real)
    return var, covs


# ============================================================================
def main():
    fc = 10e3
    poles, residues = make_bessel(4, fc)
    Af, bf, cf = modal_real_system(poles, residues)
    kco, koc = 100.0, 100.0
    Q = np.array([[-kco, kco], [koc, -koc]])
    gamma = np.array([0.0, 1.0])
    Nch = 10
    dt = 20e-6  # f_c·Δ = 0.2
    S0 = 2e-5   # S0/Δ = 1
    p_st = np.array([koc, kco]) / (kco + koc)

    def init_state(m):
        return lambda: (
            Nch * p_st.copy(),
            Nch * (np.diag(p_st) - np.outer(p_st, p_st)),
            np.zeros(m),
            np.zeros((2, m)),
            np.zeros((m, m)),
        )

    print("== O1: simulator noise vs closed forms (no channels, native read) ==")
    # Native read samples the filter output, so the references are the
    # stationary autocovariance R(0) and R(Δ) of the filtered noise (the
    # boxcar-read forms in noise_read_var_cov apply to the grouped read and
    # are validated against quadrature in test_acquisition_filter.cpp).
    Zn = simulate(3000, 0, kco, koc, Af, bf, cf, S0, dt, 60)
    Zn = Zn[:, 20:]  # discard filter transient from x(0)=0

    def R_noise(tau):
        s = 0.0 + 0.0j
        for k, rk in enumerate(residues):
            for l, rl in enumerate(residues):
                s += rk * rl * np.exp(poles[l] * tau) / (-(poles[k] + poles[l]))
        return (S0 * s).real

    var_th, cov1_th = R_noise(0.0), R_noise(dt)
    var_emp = Zn.var()
    cov1_emp = np.mean(Zn[:, 1:] * Zn[:, :-1]) - Zn.mean() ** 2
    n_eff = Zn.size
    print(f"  Var: emp {var_emp:.6g}  th {var_th:.6g}   rel err {abs(var_emp/var_th-1):.3%}")
    print(f"  Cov1: emp {cov1_emp:.6g}  th {cov1_th:.6g}  rel err {abs(cov1_emp/cov1_th-1):.3%}")
    ok1 = abs(var_emp / var_th - 1) < 5 / np.sqrt(n_eff) and abs(cov1_emp / cov1_th - 1) < 8 / np.sqrt(n_eff)
    print("  O1", "PASS" if ok1 else "FAIL")

    print("== O2: box anchor, generic assembly vs boundary-state formulas ==")
    A1 = np.zeros((1, 1)); b1 = np.ones(1); read1 = np.array([1.0 / dt])
    wq_box = window_quantities(Q, gamma, A1, b1, S0, dt)
    st = init_state(1)()
    ok2 = True
    for w in range(2):
        yr, vr, kr = box_reference_blocks(wq_box, Q, gamma, S0, dt, st[0], st[1])
        stp = predict(st, wq_box, reset_rows=(0,))
        _, zh, vv, _ = update(stp, read1, z=yr)  # z irrelevant: comparing prior blocks
        kN = stp[3] @ read1
        e1, e2, e3 = abs(zh - yr), abs(vv - vr), np.max(np.abs(kN - kr))
        print(f"  window {w}: |Δŷ|={e1:.2e}  |Δv|={e2:.2e}  |Δκ|={e3:.2e}")
        ok2 &= e1 < 1e-9 and e2 < 1e-9 and e3 < 1e-9
        st, _, _, _ = update(stp, read1, z=yr + 0.3 * np.sqrt(vv))
    print("  O2", "PASS" if ok2 else "FAIL")

    print("== O3/O4: native read (K+4), MC + whiteness ==")
    nsim, Twin = 3000, 100
    Z = simulate(nsim, Nch, kco, koc, Af, bf, cf, S0, dt, Twin)
    wq_f = window_quantities(Q, gamma, Af, bf, S0, dt)
    zh, vv = run_member(Z, wq_f, cf, init_state(Af.shape[0]))
    # window-1 predictive vs MC
    mu_emp, mu_th = Z[:, 0].mean(), zh[0, 0]
    va_emp, va_th = Z[:, 0].var(), vv[0, 0]
    se_mu = np.sqrt(va_emp / nsim)
    se_va = va_emp * np.sqrt(2.0 / (nsim - 1))
    print(f"  w1 mean: MC {mu_emp:.4f}  member {mu_th:.4f}  ({abs(mu_emp-mu_th)/se_mu:.1f} SE)")
    print(f"  w1 var:  MC {va_emp:.4f}  member {va_th:.4f}  ({abs(va_emp-va_th)/se_va:.1f} SE)")
    ok3 = abs(mu_emp - mu_th) < 4 * se_mu and abs(va_emp - va_th) < 4 * se_va
    print("  O3", "PASS" if ok3 else "FAIL")

    r = (Z - zh) / np.sqrt(vv)
    r = r[:, 10:]  # settle
    lag1 = np.mean(r[:, 1:] * r[:, :-1]) / np.mean(r * r)
    se_l = 1.0 / np.sqrt(r[:, 1:].size)
    # box member on the same filtered data
    zh_b, vv_b = run_member(Z, wq_box, read1, init_state(1), reset_rows=(0,))
    rb = (Z - zh_b) / np.sqrt(vv_b)
    rb = rb[:, 10:]
    lag1_b = np.mean(rb[:, 1:] * rb[:, :-1]) / np.mean(rb * rb)
    print(f"  residuals: mean {r.mean():+.4f}  var {r.var():.4f}  lag1 {lag1:+.4f} (SE {se_l:.4f})")
    print(f"  box member on filtered truth: var {rb.var():.4f}  lag1 {lag1_b:+.4f}")
    ok4 = abs(lag1) < 5 * se_l and abs(r.var() - 1.0) < 0.05 and lag1_b > 10 * se_l
    print("  O4", "PASS" if ok4 else "FAIL", "(box lag-1 is the motivation number)")

    print("== O5: grouped read (K+5, cumulator) MC ==")
    Ag, bg = add_cumulator(Af, bf, cf)
    mg = Ag.shape[0]
    readg = np.zeros(mg); readg[-1] = 1.0 / (5 * dt)
    Zg = simulate(1500, Nch, kco, koc, Ag, bg, readg, S0, 5 * dt, 30, reset_rows=(mg - 1,))
    wq_g = window_quantities(Q, gamma, Ag, bg, S0, 5 * dt)
    zh_g, vv_g = run_member(Zg, wq_g, readg, init_state(mg), reset_rows=(mg - 1,))
    mu_e, va_e = Zg[:, 0].mean(), Zg[:, 0].var()
    se_mu = np.sqrt(va_e / Zg.shape[0]); se_va = va_e * np.sqrt(2.0 / (Zg.shape[0] - 1))
    print(f"  w1 mean: MC {mu_e:.4f}  member {zh_g[0,0]:.4f}  ({abs(mu_e-zh_g[0,0])/se_mu:.1f} SE)")
    print(f"  w1 var:  MC {va_e:.4f}  member {vv_g[0,0]:.4f}  ({abs(va_e-vv_g[0,0])/se_va:.1f} SE)")
    rg = (Zg - zh_g) / np.sqrt(vv_g)
    rg = rg[:, 5:]
    lag1_g = np.mean(rg[:, 1:] * rg[:, :-1]) / np.mean(rg * rg)
    se_g = 1.0 / np.sqrt(rg[:, 1:].size)
    print(f"  residuals: var {rg.var():.4f}  lag1 {lag1_g:+.4f} (SE {se_g:.4f})")
    ok5 = (
        abs(mu_e - zh_g[0, 0]) < 4 * se_mu
        and abs(va_e - vv_g[0, 0]) < 4 * se_va
        and abs(lag1_g) < 5 * se_g
    )
    print("  O5", "PASS" if ok5 else "FAIL")

    print("ALL GATES:", "PASS" if all([ok1, ok2, ok3, ok4, ok5]) else "FAIL")


if __name__ == "__main__":
    main()
