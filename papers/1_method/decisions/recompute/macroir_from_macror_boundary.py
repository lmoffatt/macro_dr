# macroir_from_macror_boundary.py -- cited by 08_appendix_derivation.tex, which claims every printed
# MacroIR equation falls out of MacroR run on the boundary state and marginalized.
# Run from anywhere:  python3 papers/1_method/decisions/recompute/macroir_from_macror_boundary.py
# PROMOTED 2026-08-26 from tmp/, which is gitignored. Unchanged: it reads no files.
"""
Check Luciano's derivation of MacroIR:

  (1) take MacroR as given,
  (2) run it on the BOUNDARY state (the K^2 pairs (i0,it)) instead of the K states,
      with conductance Gamma_bar_{i0->it} and an extra "state-dependent noise"
      N_ch * sum_{i0,it} mu^bnd_{i0->it} * Vbar_{i0->it},
  (3) marginalise the result over the initial state i0.

Claim: every printed MacroIR equation of 02_theory.tex falls out, with no tilde
operator introduced anywhere.

Equations checked, by their manuscript labels:
  eq:mu-prop, eq:sig-prop, eq:ypred, eq:pred-var, eq:tilde, eq:vector-tilde,
  eq:mu-post, eq:sig-post, eq:mr-vr.

Also checks, by Monte Carlo, the boundary covariance the whole thing rests on.
"""

import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(20260806)

# ---------------------------------------------------------------- setup
K = 4
Nch = 137.0
Delta = 0.7
eps2 = 0.03          # instrumental noise over the window


def random_generator(K):
    Q = rng.uniform(0.05, 2.0, size=(K, K))
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def random_state(multinomial_start):
    """(mu0, Sigma0) per-channel: mu0 = E[N0]/Nch, Sigma0 = Cov(N0)/Nch."""
    mu0 = rng.dirichlet(np.ones(K))
    S = np.diag(mu0) - np.outer(mu0, mu0)          # multinomial, the exact case
    if not multinomial_start:
        # after one update the state is no longer multinomial: perturb, keeping
        # symmetry and the simplex constraint S @ 1 = 0
        A = rng.normal(size=(K, K)) * 0.02
        A = A + A.T
        A = A - np.outer(A.sum(1), np.ones(K)) / K
        A = A - np.outer(np.ones(K), A.sum(0)) / K
        S = S + A
    return mu0, S


Q = random_generator(K)
P = expm(Q * Delta)

# boundary-conditioned single-channel moments (arbitrary positive numbers: the
# identity below is algebraic and does not care how they were computed)
Gbar = rng.uniform(0.0, 1.0, size=(K, K))          # \overline{\Gamma}_{i0->it}
Vbar = rng.uniform(0.0, 0.3, size=(K, K))          # \overline{V}_{i0->it}

gamma_bar0 = (P * Gbar).sum(axis=1)                # eq:gammabar
G = Gbar * P                                       # eq:vector-tilde / eq:GH
H = Gbar**2 * P                                    # eq:tilde / eq:GH
sig2_gamma_bar0 = (P * Vbar).sum(axis=1)           # residual interval variance

idx = [(i0, it) for i0 in range(K) for it in range(K)]


def report(name, a, b):
    d = np.max(np.abs(np.asarray(a) - np.asarray(b)))
    scale = max(1.0, np.max(np.abs(np.asarray(b))))
    print(f"  {name:<58s} max|diff| = {d:.3e}   (rel {d/scale:.3e})")
    return d / scale


# ------------------------------------------------- MacroR, verbatim, on any space
def macror(mu, Sigma, gamma, extra_var, y_obs, Nch):
    """MacroR (Moffatt 2007) in the manuscript's per-channel normalisation:
    mu = E[N]/Nch, Sigma = Cov(N)/Nch.  extra_var is any observation noise that
    does not covary with the state."""
    y_pred = Nch * mu @ gamma
    var = extra_var + Nch * gamma @ Sigma @ gamma
    g = Sigma @ gamma                              # = Cov(N, y)/Nch
    delta = y_obs - y_pred
    mu_post = mu + g * delta / var
    Sigma_post = Sigma - Nch * np.outer(g, g) / var
    return y_pred, var, g, mu_post, Sigma_post


worst = 0.0
for multinomial_start in (True, False):
    tag = "multinomial start" if multinomial_start else "perturbed (post-update) start"
    print(f"\n=== {tag} ===")
    mu0, S0 = random_state(multinomial_start)
    SmD = S0 - np.diag(mu0)                        # the recurring Sigma0 - diag(mu0)

    # ---------------- step 2a: the boundary state's own (mu, Sigma), per channel
    mu_bnd = np.array([mu0[i0] * P[i0, it] for (i0, it) in idx])
    S_bnd = np.zeros((K * K, K * K))
    for a, (i0, it) in enumerate(idx):
        for b, (j0, jt) in enumerate(idx):
            S_bnd[a, b] = P[i0, it] * SmD[i0, j0] * P[j0, jt]
            if i0 == j0 and it == jt:
                S_bnd[a, b] += mu0[i0] * P[i0, it]

    # ---------------- step 2b: MacroR on that space, conductance = Gamma_bar,
    #                  plus the state-dependent noise as an extra variance
    gamma_bnd = Gbar.reshape(-1)
    state_dep_noise = Nch * (mu_bnd @ Vbar.reshape(-1))
    y_obs = 11.3

    y_pred_b, var_b, g_bnd, mu_bnd_post, S_bnd_post = macror(
        mu_bnd, S_bnd, gamma_bnd, eps2 + state_dep_noise, y_obs, Nch)

    # ---------------- step 3: marginalise over the initial state i0
    M = np.zeros((K, K * K))                       # marginaliser: sums over i0
    for a, (i0, it) in enumerate(idx):
        M[it, a] = 1.0
    mu_prop_m = M @ mu_bnd
    S_prop_m = M @ S_bnd @ M.T
    g_m = M @ g_bnd
    mu_post_m = M @ mu_bnd_post
    S_post_m = M @ S_bnd_post @ M.T

    # ---------------- the printed MacroIR equations of 02_theory.tex
    mu_prop = P.T @ mu0                                            # eq:mu-prop
    S_prop = P.T @ SmD @ P + np.diag(mu_prop)                      # eq:sig-prop
    y_pred = Nch * mu0 @ gamma_bar0                                # eq:ypred
    tilde = gamma_bar0 @ SmD @ gamma_bar0 + mu0 @ (H.sum(axis=1))  # eq:tilde
    var = eps2 + Nch * tilde + Nch * mu0 @ sig2_gamma_bar0         # eq:pred-var
    g = P.T @ SmD @ gamma_bar0 + G.T @ mu0                         # eq:vector-tilde
    delta = y_obs - y_pred
    mu_post = mu_prop + g * delta / var                            # eq:mu-post
    S_post = S_prop - Nch * np.outer(g, g) / var                   # eq:sig-post

    print(" boundary MacroR, marginalised over i0   vs   printed MacroIR equation")
    worst = max(worst, report("mu^prop            (eq:mu-prop)", mu_prop_m, mu_prop))
    worst = max(worst, report("Sigma^prop         (eq:sig-prop)", S_prop_m, S_prop))
    worst = max(worst, report("ybar^pred          (eq:ypred)", y_pred_b, y_pred))
    worst = max(worst, report("sigma^2_pred       (eq:pred-var + eq:tilde)", var_b, var))
    worst = max(worst, report("g                  (eq:vector-tilde)", g_m, g))
    worst = max(worst, report("mu^post            (eq:mu-post)", mu_post_m, mu_post))
    worst = max(worst, report("Sigma^post         (eq:sig-post)", S_post_m, S_post))

    # the bilinear tilde is just gamma^T Sigma gamma read on the boundary space
    worst = max(worst, report("tilde == gamma_bnd^T Sigma_bnd gamma_bnd",
                              gamma_bnd @ S_bnd @ gamma_bnd, tilde))

    # ---------------- the same recipe on the START state alone gives MR
    total_var_i0 = (P * (Vbar + (Gbar - gamma_bar0[:, None])**2)).sum(axis=1)
    y_pred_mr, var_mr, g_mr, _, _ = macror(
        mu0, S0, gamma_bar0, eps2 + Nch * mu0 @ total_var_i0, y_obs, Nch)
    worst = max(worst, report("MR and IR predict the same variance (eq:mr-vr)", var_mr, var_b))
    gap = np.max(np.abs(P.T @ g_mr - g))
    print(f"  {'MR gain vs IR gain (must NOT vanish)':<58s} max|diff| = {gap:.3e}"
          f"   [{'ok, they differ' if gap > 1e-6 else 'FAIL: identical'}]")

# ---------------------------------------------------------------- Monte Carlo
# the boundary covariance Sigma^bnd is the one input that is not a definition,
# so check it against simulated channels rather than trusting the algebra twice
print("\n=== Monte Carlo on the boundary covariance (multinomial start) ===")
mu0 = rng.dirichlet(np.ones(K))
SmD = np.diag(mu0) - np.outer(mu0, mu0) - np.diag(mu0)
Nch_i = 200
reps = 400000
start = rng.multinomial(Nch_i, mu0, size=reps)                     # N_0
bnd = np.zeros((reps, K, K))
for i0 in range(K):
    bnd[:, i0, :] = np.array([rng.multinomial(n, P[i0]) for n in start[:, i0]])
flat = bnd.reshape(reps, K * K)
emp = np.cov(flat.T) / Nch_i
theo = np.zeros((K * K, K * K))
for a, (i0, it) in enumerate(idx):
    for b, (j0, jt) in enumerate(idx):
        theo[a, b] = P[i0, it] * SmD[i0, j0] * P[j0, jt]
        if i0 == j0 and it == jt:
            theo[a, b] += mu0[i0] * P[i0, it]
err = np.max(np.abs(emp - theo))
C = Nch_i * theo                                                   # covariance of raw counts
se = np.sqrt((np.outer(np.diag(C), np.diag(C)) + C**2) / reps) / Nch_i
z = np.max(np.abs(emp - theo) / se)
print(f"  Sigma^bnd, empirical vs formula                    max|diff| = {err:.3e}")
print(f"  worst entry, in standard errors of the MC estimate            z = {z:.2f}"
      f"   [{'ok' if z < 5 else 'FAIL'}]")

print(f"\nworst relative discrepancy over all algebraic checks: {worst:.3e}")
