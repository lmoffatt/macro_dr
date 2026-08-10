"""MR contra IR: MacroR sobre el boundary state, y qué cambia si se marginaliza antes.

Acompaña a theory/macroir/notes/mr_vs_ir_from_macror.md. No lee datos: reconstruye
los objetos de la ventana con Van Loan (NO por la vía espectral, para no heredar
las ramas de coincidencia de E_2/E_3) y verifica el documento contra las fórmulas
que corren en legacy/qmodel.h:4557-4645.

Numpy + scipy, sin archivos, ~2 min (los Monte Carlo dominan).

El punto de partida es el BOUNDARY STATE: el indicador X_bnd sobre los K^2 pares
(i0, it) de estados en los dos extremos de la ventana. La corriente promediada al
intervalo no es lineal en la ocupación de ningún instante, y sí lo es en el par,
que es lo que permite correrle MacroR encima sin cambiar nada de MacroR.

Notación (la del manuscrito; el nombre del código entre corchetes):
    P            transición sobre la ventana                        [P]
    Gam_ij       media de la conductancia promediada al intervalo,
                 condicional al par                                 [gmean_ij]
    Vbar_ij      su varianza, condicional al par                    [gvar_ij]
    gbar0_i      Gam_ij promediado sobre it con pesos P             [gmean_i]
    mu0, S0      media y covarianza del prior al abrir la ventana   [P_mean, P_Cov]
    mu_bnd, S_bnd  media y covarianza de X_bnd sobre los K^2 pares  (no está en el código)
    D_ij         Gam_ij - gbar0_i                                   (no está en el código)

SECCIONES
  (0) El boundary state: mu_bnd y S_bnd contra Monte Carlo.
  (1) MacroR sobre el boundary state, literal, y su marginalización: reproduce
      exactamente lo que el código computa para IR (gSg, ms, gS en av=2).
  (2) MR: marginalizar it ANTES de condicionar. Reproduce av=1.
  (3) Las cuatro cajas y el bloque cruzado, escritos con D.
  (4) Var(corriente) contra Monte Carlo, y el tamaño de lo que MR reasigna.
"""
import numpy as np
from scipy.linalg import expm

K = 2
rng = np.random.default_rng(20260808)


# ----------------------------------------------------------------- la ventana
def window_objects(kon_A, koff, dt, g_unit=1.0):
    """P, Gam_ij, M2_ij, Vbar_ij por Van Loan."""
    Q = np.array([[-kon_A, kon_A], [koff, -koff]])
    G = np.diag([0.0, g_unit])
    P = expm(Q * dt)
    Z = np.zeros((K, K))
    E3 = expm(np.block([[Q, G, Z], [Z, Q, G], [Z, Z, Q]]) * dt)
    F1, F2 = E3[0:K, K:2 * K], E3[0:K, 2 * K:3 * K]
    Gam = F1 / (dt * P)
    M2 = 2 * F2 / (dt ** 2 * P)
    return P, Gam, M2, M2 - Gam ** 2


# ------------------------------------------------------ el estado generalizado
def boundary_state(P, mu0, S0):
    """mu_bnd (K^2) y S_bnd (K^2 x K^2) del indicador del par, un canal."""
    mu_bnd = (mu0[:, None] * P).reshape(-1)
    SmD = S0 - np.diag(mu0)                      # el prior sin su parte multinomial
    S_bnd = np.zeros((K * K, K * K))
    for i0 in range(K):
        for it in range(K):
            a = i0 * K + it
            for j0 in range(K):
                for jt in range(K):
                    b = j0 * K + jt
                    S_bnd[a, b] = P[i0, it] * SmD[i0, j0] * P[j0, jt]
    S_bnd += np.diag(mu_bnd)                      # el par que el canal efectivamente toma
    return mu_bnd, S_bnd


def macror_blocks(mu, S, cond, noise_per_state):
    """MacroR literal sobre CUALQUIER estado: los tres bloques, por canal.

    cond            conductancia por estado (lo que en MacroR es gamma)
    noise_per_state varianza de la corriente dado el estado (0 en MacroR original)
    """
    var_del_estado = cond @ S @ cond
    var_de_la_corriente = mu @ noise_per_state
    cov_estado_corriente = S @ cond
    media = mu @ cond
    return media, var_del_estado, var_de_la_corriente, cov_estado_corriente


# =============================================================== (0) y (1)
print("=" * 100)
print("(0)+(1)  MacroR sobre el boundary state; contra Monte Carlo y contra el código (av=2)")
print("=" * 100)
kon = koff = 100.
rates = np.array([kon, koff])
for dtn in [0.1, 1.0]:
    dt = dtn / koff
    mu0 = np.array([0.5, 0.5])
    S0 = np.diag(mu0) - np.outer(mu0, mu0)        # multinomial: simulable
    P, Gam, M2, Vbar = window_objects(kon, koff, dt)
    mu_bnd, S_bnd = boundary_state(P, mu0, S0)

    # --- Monte Carlo del par y de la corriente
    n = 300_000
    st0 = (rng.random(n) < 0.5).astype(int)
    st = st0.copy()
    t = np.zeros(n); acc = np.zeros(n); vivo = np.ones(n, bool)
    while vivo.any():
        idx = np.flatnonzero(vivo); s = st[idx]
        fin = np.minimum(t[idx] + rng.exponential(1.0 / rates[s]), dt)
        acc[idx] += (fin - t[idx]) * s
        t[idx] = fin
        salta = t[idx] < dt
        st[idx[salta]] = 1 - s[salta]
        vivo[idx[~salta]] = False
    A = acc / dt
    par = st0 * K + st
    onehot = np.zeros((n, K * K)); onehot[np.arange(n), par] = 1.0
    mc_mu = onehot.mean(0)
    mc_S = np.cov(onehot.T, ddof=1)
    mc_cov = np.array([np.cov(onehot[:, a], A, ddof=1)[0, 1] for a in range(K * K)])

    # --- MacroR sobre el boundary state, sin tocar MacroR
    Gv, Vv = Gam.reshape(-1), Vbar.reshape(-1)
    media, v_est, v_cor, cov_bnd = macror_blocks(mu_bnd, S_bnd, Gv, Vv)

    # --- marginalizar i0 (lo que la ventana siguiente no necesita)
    cov_marg = cov_bnd.reshape(K, K).sum(0)

    # --- lo que computa el código para av=2
    SmD = S0 - np.diag(mu0)
    gbar0 = (P * Gam).sum(1)
    gsqr_i = (P * M2).sum(1)
    Gm, H = Gam * P, (Gam ** 2) * P
    cod_gSg = gbar0 @ SmD @ gbar0 + mu0 @ H.sum(1)
    cod_ms = mu0 @ (gsqr_i - H.sum(1))
    cod_gS = P.T @ (SmD @ gbar0) + Gm.T @ mu0

    print(f"  dt*koff = {dtn}")
    print(f"    mu_bnd  contra MC : {np.abs(mu_bnd-mc_mu).max():.1e}"
          f"     S_bnd contra MC : {np.abs(S_bnd-mc_S).max():.1e}")
    print(f"    Cov(X_bnd, y) contra MC : {np.abs(cov_bnd-mc_cov).max():.1e}")
    print(f"    var del estado  : MacroR/bnd {v_est:.12f}  código {cod_gSg:.12f}"
          f"   dif {abs(v_est-cod_gSg):.1e}")
    print(f"    var de corriente: MacroR/bnd {v_cor:.12f}  código {cod_ms:.12f}"
          f"   dif {abs(v_cor-cod_ms):.1e}")
    print(f"    Cov marginalizada contra el código: {np.abs(cov_marg-cod_gS).max():.1e}")

# =============================================================== (2) y (3)
print()
print("=" * 100)
print("(2)+(3)  MR = marginalizar it ANTES de condicionar.  Las cuatro cajas con D.")
print("=" * 100)
CASOS = [(kon_, koff_, dtn, sc)
         for (kon_, koff_) in [(100., 100.), (300., 100.), (40., 250.), (1000., 10.)]
         for dtn in [0.01, 0.1, 1.0]
         for sc in [1.0, 0.31, 0.0]]
peor = {}
for kon_, koff_, dtn, sc in CASOS:
    dt = dtn / koff_
    p = kon_ / (kon_ + koff_)
    mu0 = np.array([1 - p, p])
    S0 = sc * (np.diag(mu0) - np.outer(mu0, mu0))
    P, Gam, M2, Vbar = window_objects(kon_, koff_, dt)
    mu_bnd, S_bnd = boundary_state(P, mu0, S0)
    Gv, Vv = Gam.reshape(-1), Vbar.reshape(-1)

    # IR: MacroR sobre el par, y marginalizar despues
    _, ir_est, ir_cor, ir_cov_bnd = macror_blocks(mu_bnd, S_bnd, Gv, Vv)
    ir_cov = ir_cov_bnd.reshape(K, K).sum(0)

    # MR: marginalizar it ANTES.  El estado vuelve a ser X_0 (mu0, S0); la
    # conductancia es gbar0 y el ruido por estado es la varianza TOTAL, que es
    # lo que la ley de la varianza total deja al promediar sobre it.
    gbar0 = (P * Gam).sum(1)
    var_tot = (P * (Vbar + (Gam - gbar0[:, None]) ** 2)).sum(1)
    _, mr_est, mr_cor, mr_cov0 = macror_blocks(mu0, S0, gbar0, var_tot)
    mr_cov = P.T @ mr_cov0                       # propagar al final de la ventana

    # la reescritura con D
    D = Gam - gbar0[:, None]
    prior = float(gbar0 @ S0 @ gbar0)
    ED2 = float((mu0[:, None] * P * D ** 2).sum())
    EVbar = float(mu0 @ (P * Vbar).sum(1))
    D_una_vez = (mu0[:, None] * P * D).sum(0)

    e = {
        "IR estado = prior + E[D2]": abs(ir_est - (prior + ED2)),
        "IR corriente = E[Vbar]": abs(ir_cor - EVbar),
        "MR estado = prior": abs(mr_est - prior),
        "MR corriente = E[Vbar] + E[D2]": abs(mr_cor - (EVbar + ED2)),
        "Var(corriente) total: MR = IR": abs((ir_est + ir_cor) - (mr_est + mr_cor)),
        "Cov_IR - Cov_MR = D contraido una vez": np.abs((ir_cov - mr_cov) - D_una_vez).max(),
        "ese vector suma cero sobre it": abs(D_una_vez.sum()),
        "D suma cero por FILA": np.abs((P * D).sum(1)).max(),
    }
    for k, v in e.items():
        peor[k] = max(peor.get(k, 0.0), float(v))
for k, v in peor.items():
    print(f"  {k:<42} error max = {v:.2e}   ({len(CASOS)} casos)")

# =============================================================== (4)
print()
print("=" * 100)
print("(4)  Var(corriente) por canal contra Monte Carlo, y el tamaño de lo reasignado")
print("=" * 100)
print(f"{'dt*koff':>9} {'MR':>14} {'IR':>14} {'MonteCarlo':>14} {'E[D2]/(prior+E[D2])':>22}")
for dtn in [0.01, 0.1, 0.5, 1.0, 2.0]:
    dt = dtn / 100.
    mu0 = np.array([0.5, 0.5])
    S0 = np.diag(mu0) - np.outer(mu0, mu0)
    P, Gam, M2, Vbar = window_objects(100., 100., dt)
    mu_bnd, S_bnd = boundary_state(P, mu0, S0)
    Gv, Vv = Gam.reshape(-1), Vbar.reshape(-1)
    _, ir_est, ir_cor, _ = macror_blocks(mu_bnd, S_bnd, Gv, Vv)
    gbar0 = (P * Gam).sum(1)
    var_tot = (P * (Vbar + (Gam - gbar0[:, None]) ** 2)).sum(1)
    _, mr_est, mr_cor, _ = macror_blocks(mu0, S0, gbar0, var_tot)
    D = Gam - gbar0[:, None]
    prior = float(gbar0 @ S0 @ gbar0)
    ED2 = float((mu0[:, None] * P * D ** 2).sum())
    n = 300_000
    st = (rng.random(n) < 0.5).astype(int)
    t = np.zeros(n); acc = np.zeros(n); vivo = np.ones(n, bool)
    while vivo.any():
        idx = np.flatnonzero(vivo); s = st[idx]
        fin = np.minimum(t[idx] + rng.exponential(1.0 / np.array([100., 100.])[s]), dt)
        acc[idx] += (fin - t[idx]) * s
        t[idx] = fin
        salta = t[idx] < dt
        st[idx[salta]] = 1 - s[salta]
        vivo[idx[~salta]] = False
    print(f"{dtn:>9} {mr_est+mr_cor:>14.9f} {ir_est+ir_cor:>14.9f} "
          f"{(acc/dt).var(ddof=1):>14.9f} {ED2/(prior+ED2):>22.4f}")
print("  con S0 = 0 la fracción es 1 en toda la columna: MR no puede corregir nada.")


# =============================================================== (5)
# EL TILDE: MacroIR como MacroR, literalmente.
#
# MacroR tiene DOS contracciones y nada mas:  gamma' Sigma gamma  y  Sigma gamma.
# El tilde marca "la misma contraccion, tomada sobre los K^2 pares". Si el tilde
# es de veras una sustitucion y no un operador nuevo, entonces las ecuaciones de
# MacroR valen VERBATIM y los tres miembros son tres lecturas del mismo par de
# contracciones:
#     R    sobre los K estados con gamma  (instantanea, marco de media ventana)
#     MR   sobre los K estados con gbar0  (media del intervalo dado el inicio)
#     IR   sobre los K^2 pares con Gam    (media del intervalo dado el par)
#
# (T1) las tres lecturas reproducen los tres miembros
# (T2) el EXCESO del tilde sobre la version sin tilde es D contraido tantas
#      veces como factores de conductancia tiene la contraccion:
#          escalar (dos gammas) -> E[D^2]
#          vector  (una gamma)  -> D contraido una vez, indice final abierto
# (T3) marginalizar i0 DESPUES de condicionar conmuta exactamente; es lo que
#      cierra la extension y devuelve un estado de K componentes a la ventana
#      siguiente.
print()
print("=" * 100)
print("(5)  El tilde: MacroIR = MacroR con las dos contracciones leidas sobre los pares")
print("=" * 100)

def marg_i0(K_):
    """matriz K x K^2 que suma el indice inicial (deja libre el final)."""
    M = np.zeros((K_, K_ * K_))
    for i0 in range(K_):
        for it in range(K_):
            M[it, i0 * K_ + it] = 1.0
    return M

M = marg_i0(K)
peor5 = {}
for kon_, koff_, dtn, sc in CASOS:
    dt = dtn / koff_
    p = kon_ / (kon_ + koff_)
    mu0 = np.array([1 - p, p])
    S0 = sc * (np.diag(mu0) - np.outer(mu0, mu0))
    P, Gam, M2, Vbar = window_objects(kon_, koff_, dt)
    mu_bnd, S_bnd = boundary_state(P, mu0, S0)
    Gv, Vv = Gam.reshape(-1), Vbar.reshape(-1)
    gbar0 = (P * Gam).sum(1)
    D = Gam - gbar0[:, None]

    # (T2) el exceso del tilde, en las dos contracciones
    esc_tilde = Gv @ S_bnd @ Gv                 # gamma' Sigma gamma  con tilde
    esc_plano = gbar0 @ S0 @ gbar0              # la misma, sin tilde
    vec_tilde = M @ (S_bnd @ Gv)                # Sigma gamma con tilde, marginalizado
    vec_plano = P.T @ (S0 @ gbar0)              # la misma, sin tilde, propagada
    ED2 = float((mu0[:, None] * P * D ** 2).sum())
    D1 = (mu0[:, None] * P * D).sum(0)

    # (T3) conmutacion: marginalizar i0 despues de actualizar
    eps2, N, delta = 1e-4, 100.0, 0.037
    s2 = eps2 + N * (esc_tilde + mu_bnd @ Vv)
    g_bnd = S_bnd @ Gv
    mu_post_bnd = mu_bnd + (g_bnd / s2) * delta
    S_post_bnd = S_bnd - (N / s2) * np.outer(g_bnd, g_bnd)
    mu_prop = P.T @ mu0
    S_prop = P.T @ (S0 - np.diag(mu0)) @ P + np.diag(mu_prop)
    g = M @ g_bnd
    mu_post = mu_prop + (g / s2) * delta
    S_post = S_prop - (N / s2) * np.outer(g, g)

    e5 = {
        "escalar: tilde - plano = E[D2]": abs((esc_tilde - esc_plano) - ED2),
        "vector : tilde - plano = D una vez": np.abs((vec_tilde - vec_plano) - D1).max(),
        "marginal de mu_bnd = mu_prop": np.abs(M @ mu_bnd - mu_prop).max(),
        "marginal de S_bnd  = S_prop": np.abs(M @ S_bnd @ M.T - S_prop).max(),
        "conmuta: marg(mu_post_bnd) = mu_post": np.abs(M @ mu_post_bnd - mu_post).max(),
        "conmuta: marg(S_post_bnd)  = S_post": np.abs(M @ S_post_bnd @ M.T - S_post).max(),
    }
    for k, v in e5.items():
        peor5[k] = max(peor5.get(k, 0.0), float(v))
for k, v in peor5.items():
    print(f"  {k:<40} error max = {v:.2e}   ({len(CASOS)} casos)")


# =============================================================== (6)
# R EN EL MARCO.  R no lleva termino de ruido por estado: la guardia del codigo
# es `variance::value && averaging::value>0` (qmodel.h:4586), asi que av=0 tiene
# solo las dos contracciones.  Lo que R SI tiene y las otras dos lecturas no es
# un ORIGEN DE TIEMPO propio: para av=0 el codigo relocaliza el prior a media
# ventana antes de nada (qmodel.h:4541-4555) y t_P pasa a ser P_half.
#
# La pregunta: eso es una lectura mas del marco (estado = ocupacion a media
# ventana) o es un caso aparte?
print()
print("=" * 100)
print("(6)  R en el marco: estado = ocupacion a media ventana")
print("=" * 100)
kon = koff = 100.
rates = np.array([kon, koff])
gamma = np.array([0.0, 1.0])
print(f"{'dt*koff':>8} {'Sig_mid vs MC':>14} {'contracc. vs código':>21} "
      f"{'media R/verdad':>15} {'var R/verdad':>13} {'cruz R/verdad':>14} {'cruz MR/verdad':>15}")
for dtn in [0.01, 0.1, 0.5, 1.0, 2.0]:
    dt = dtn / koff
    mu0 = np.array([0.5, 0.5])
    S0 = np.diag(mu0) - np.outer(mu0, mu0)
    P, Gam, M2, Vbar = window_objects(kon, koff, dt)
    Ph = expm(np.array([[-kon, kon], [koff, -koff]]) * dt / 2)

    # --- el estado de R: la ocupacion a media ventana, con SU gaussiana
    mu_mid = Ph.T @ mu0
    S_mid = Ph.T @ (S0 - np.diag(mu0)) @ Ph + np.diag(mu_mid)

    # --- MacroR sobre ese estado, ruido por estado CERO, y mapa al final por la
    #     mitad que falta
    media_R, var_est_R, var_cor_R, cov_R0 = macror_blocks(mu_mid, S_mid, gamma,
                                                          np.zeros(K))
    cov_R = Ph.T @ cov_R0

    # --- lo que computa el codigo en av=0 (qmodel.h:4578 y :4643, rama else)
    cod_gSg = gamma @ S_mid @ gamma
    cod_gS = Ph.T @ (S_mid @ gamma)
    cod_mean = mu_mid @ gamma

    # --- Monte Carlo: la ocupacion a media ventana es un estado de verdad?
    n = 300_000
    st = (rng.random(n) < 0.5).astype(int)
    t = np.zeros(n); vivo = np.ones(n, bool); mid = np.zeros(n, int)
    while vivo.any():
        idx = np.flatnonzero(vivo); s = st[idx]
        fin = np.minimum(t[idx] + rng.exponential(1.0 / rates[s]), dt / 2)
        t[idx] = fin
        salta = t[idx] < dt / 2
        st[idx[salta]] = 1 - s[salta]
        vivo[idx[~salta]] = False
    mid = st
    oh = np.zeros((n, K)); oh[np.arange(n), mid] = 1.0
    err_S = np.abs(np.cov(oh.T, ddof=1) - S_mid).max()

    # --- la verdad: los bloques exactos, que son los de IR
    mu_bnd, S_bnd = boundary_state(P, mu0, S0)
    Gv, Vv = Gam.reshape(-1), Vbar.reshape(-1)
    med_v, v_est_v, v_cor_v, cov_bnd = macror_blocks(mu_bnd, S_bnd, Gv, Vv)
    cov_v = cov_bnd.reshape(K, K).sum(0)
    gbar0 = (P * Gam).sum(1)
    var_tot = (P * (Vbar + (Gam - gbar0[:, None]) ** 2)).sum(1)
    _, _, _, mr0 = macror_blocks(mu0, S0, gbar0, var_tot)
    cov_MR = P.T @ mr0

    err_contr = max(abs(var_est_R - cod_gSg), np.abs(cov_R - cod_gS).max(),
                    abs(media_R - cod_mean))
    print(f"{dtn:>8} {err_S:>14.1e} {err_contr:>21.1e} "
          f"{media_R/med_v:>15.5f} {(var_est_R+var_cor_R)/(v_est_v+v_cor_v):>13.5f} "
          f"{cov_R[1]/cov_v[1]:>14.5f} {cov_MR[1]/cov_v[1]:>15.5f}")
print()
print("  media R/verdad : regla del punto medio contra la integral exacta")
print("  var R/verdad   : a R le falta TODA la varianza de intervalo (no tiene ms)")
print("  cruz R/verdad  : R queda cerca; MR se va por e^-(kon+koff)dt")
