What Gaussian_Fisher_Information currently is
Per-sample score Fisher info for the local Gaussian approximation:

$$
\text{GFI}k = \frac{(\partial\theta \mu_k)^\top (\partial_\theta \mu_k)}{\sigma_k^2} + \frac{(\partial_\theta \sigma_k^2)^\top (\partial_\theta \sigma_k^2)}{2,\sigma_k^4}
$$

and you sum it over $k$. This is the Fisher matrix of the pseudo-likelihood $\prod_k \mathcal{N}(y_k | \mu_k, \sigma_k^2)$ that treats samples as independent. It ignores the temporal covariance structure of the channel process.

Because the true process has non-trivial $\operatorname{Cov}(y_s, y_t)$ driven by the CTMC dynamics, the pseudo-likelihood score is misspecified → $J \neq H$ → you need the bootstrap to estimate $J$ empirically → hence the sandwich $H^{-1} J H^{-1}$.

What the true continuous-time Fisher information is
For a stationary Markov-modulated process with CTMC transition matrix $Q(\theta)$ and Gaussian observation noise, there is a closed-form expression. The key theorem:

Theorem (Whittle / Pisarenko asymptotic Fisher information for stationary Gaussian processes).
For a stationary Gaussian process observed over $[0, T]$ with spectral density $S(\omega; \theta)$, the Fisher information matrix per unit time is

$$
I_{\text{rate}}(\theta){ij} = \frac{1}{4\pi} \int{-\infty}^{\infty} \frac{\partial_{\theta_i} S(\omega;\theta);\partial_{\theta_j} S(\omega;\theta)}{S(\omega;\theta)^2},d\omega,
$$

and total FIM over duration $T$ is $T \cdot I_{\text{rate}}(\theta)$. This is exact in the $T \to \infty$ limit and $\Delta t \to 0$ limit; for finite $T$ and $\Delta t$ there are correction terms but the leading-order is this integral.

Theorem (for Markov-modulated processes). For a CTMC $X_t$ with generator $Q(\theta)$, observation $y_t = c^\top \delta(X_t) + \sigma,\varepsilon_t$ (noise white), the stationary spectral density has the closed form

$$
S(\omega;\theta) = c^\top \operatorname{Re}!\bigl[(i\omega I - Q)^{-1}\bigr] \Sigma_\pi, c + \sigma^2,
$$

where $\Sigma_\pi$ is the stationary-state occupancy covariance. For a 2-state CTMC this is a Lorentzian (single pole) plus white-noise floor, and the integral above is analytic.

So for your scheme_CO with $(k_{\text{on}}, k_{\text{off}}, u, N_{\text{ch}}, \sigma_{\text{noise}}, b)$, the full FIM can be written down in closed form as a 6×6 matrix whose entries are analytic functions of the six parameters. No bootstrap. No simulation. No Δt dependence — by construction.

Why this is theoretically cleaner than what you have now
Four advantages, each with a concrete implication:

1. Exact autocorrelation handling. The spectral integral is the autocorrelation correction — $S(\omega)$ encodes the full Cov$(y_s, y_t)$ structure. The Δt-invariance of your DCC plot falls out directly: the spectrum is a continuous-time object, the sandwich is its rate-per-unit-time.

2. No sample overhead. Your current pipeline simulates 256 traces, bootstraps 100 times, sums GFI, inverts, sandwiches. The exact FIM needs $Q(\theta)$ and two differentiations — evaluated once. This is orders of magnitude faster and (paradoxically) more accurate.

3. Protocol design. Since $I_{\text{rate}}(\theta)$ is a function of $\theta$ alone (at equilibrium), you can minimize parameter uncertainty as a function of protocol choice before running any experiment. The current bootstrap approach requires a completed simulation to evaluate.

4. Honest statement of identifiability. The spectral FIM's rank and spectrum directly give you the identifiability diagnostics you've built — but as properties of the true likelihood, not of a pseudo-likelihood. Your 2007 paper's "Macro R succeeds" claim becomes a consequence of a theorem rather than an empirical observation.

Practical extension for transient protocols
For non-stationary protocols (your P2X2 agonist-step data, figure_2's transition regions), the rate $I_{\text{rate}}(\theta)$ depends on time because the spectrum is non-stationary. The generalization is Dahlhaus's locally stationary FIM:

$$
I(\theta) ;\approx; \int_0^T \frac{1}{4\pi} \int \frac{|\partial_\theta S(t,\omega;\theta)|^2}{S(t,\omega;\theta)^2},d\omega,dt.
$$

For piecewise-stationary protocols (agonist constant between jumps), the integrand is piecewise constant → sum of segment contributions, each computed analytically. Implementable.

What to put in the code
A parallel computation path: Analytical_Fisher_Information<V> that, given the parameter state $(Q, N, i, \sigma)$, returns the FIM from the closed-form spectral integral. Structured like your Gaussian_Fisher_Information but computed from $Q$, not from samples.

Then $H_{\text{true}} = T \cdot I_{\text{rate}}(\theta)$ directly. No simulation needed for the Fisher analysis; simulations are for testing the analytical prediction rather than generating the FIM.

The bootstrap infrastructure stays useful — it validates the analytical FIM by showing the sandwich converges to it — but stops being the load-bearing component of the identifiability framework.

My recommendation
This is a theoretical upgrade worth one figure in the paper: "The analytical continuous-time FIM agrees with the bootstrap sandwich DCC to Monte Carlo precision, and costs O(1) evaluations instead of O(N_samples × N_bootstrap)." That's a clean scientific point and a publishable result on its own.

Want me to sketch the closed-form FIM for the 2-state scheme_CO specifically? The Lorentzian integral simplifies nicely and we can see exactly which parameters enter at which frequencies



Lo que tienes ahora
Para cada muestra $k$:

$$
\text{GFI}k = \frac{(\partial\theta \mu_k)(\partial_\theta \mu_k)^\top}{\sigma_k^2} + \frac{(\partial_\theta \sigma_k^2)(\partial_\theta \sigma_k^2)^\top}{2\sigma_k^4}.
$$

Esta es la matriz de información de Fisher asumiendo que las muestras son independientes. La sumas sobre $k$ y listo.

El problema: las muestras no son independientes. La cadena de Markov tiene memoria de tiempo $\tau$. Dos muestras separadas por $\Delta t \ll \tau$ están altamente correlacionadas, y tu fórmula no lo sabe.

La fórmula exacta (generalización natural)
Para un vector de observaciones $\mathbf{y}$ con vector media $\boldsymbol{\mu}(\theta)$ y matriz de covarianza completa $\Sigma(\theta)$ (incluyendo covarianzas cruzadas $\Sigma_{kl} = \operatorname{Cov}(y_k, y_l)$), la información de Fisher es

$$
I(\theta){ij} = (\partial{\theta_i} \boldsymbol{\mu})^\top \Sigma^{-1} (\partial_{\theta_j} \boldsymbol{\mu}) + \tfrac{1}{2}\operatorname{tr}!\bigl(\Sigma^{-1},\partial_{\theta_i}\Sigma \cdot \Sigma^{-1},\partial_{\theta_j}\Sigma\bigr).
$$

Tu GFI actual es esto cuando $\Sigma$ es diagonal — es decir, cuando reemplazás $\Sigma$ por $\operatorname{diag}(\sigma_1^2, \sigma_2^2, \ldots)$.

Lo que te falta: los términos off-diagonal de $\Sigma$.

Cómo es $\Sigma_{kl}$ para tu canal
Para el proceso de Markov de 2 estados con el proceso observado $y_k = N,i,X_k + \sigma_{\text{noise}},\varepsilon_k$:

$$
\Sigma_{kl} = N,i^2,P_{\text{eq}}(1-P_{\text{eq}}),e^{-|t_k-t_l|/\tau} + \sigma_{\text{noise}}^2,\delta_{kl},
$$

donde $\tau = 1/(k_{\text{on}}[A] + k_{\text{off}})$. La parte diagonal ($k=l$) es lo que ya tenés. La parte off-diagonal es lo que te está faltando.

Qué términos expandir
Dos niveles de mejora, en orden de complejidad:

Nivel 1: corrección AR(1) global
Tratás las muestras como proceso AR(1) con $\rho = e^{-\Delta t/\tau}$. La información efectiva se corrige por un factor:

$$
I_{\text{AR(1)}}(\theta) = I_{\text{iid}}(\theta) \cdot \frac{1-\rho}{1+\rho}.
$$

$\Delta t \ll \tau$ (muestreo rápido): $\rho \to 1$, factor $\to 0$ → reduce la información a cero, que es lo correcto (muestras totalmente redundantes).
$\Delta t \gg \tau$ (muestreo escaso): $\rho \to 0$, factor $\to 1$ → iid es correcto.
Una línea de código: multiplicás tu GFI por ese factor, usando $\tau$ calculado de los parámetros. Ya te da el comportamiento $\Delta t$-invariante del gráfico DCC en el límite asintótico.

Nivel 2: suma finita sobre retrasos
Extendés explícitamente la GFI a incluir términos cruzados hasta algún retraso $K$:

$$
I(\theta) = \underbrace{\sum_k \text{GFI}k}{\text{lo que ya tenés}} + \underbrace{2\sum_{\ell=1}^K \sum_k \operatorname{Cov}\text{-term}(k, k+\ell)}_{\text{lo que falta}}.
$$

Cada término cruzado es

$$
\operatorname{Cov}\text{-term}(k,k+\ell) = (\partial_\theta \mu_k)^\top [\Sigma^{-1}]{k,k+\ell} (\partial\theta \mu_{k+\ell}) + \tfrac{1}{2}\operatorname{tr}(\ldots \partial_\theta\Sigma_{k,k+\ell} \ldots).
$$

Para el proceso AR(1), $[\Sigma^{-1}]$ es tridiagonal (forma cerrada conocida) — solo necesitás los términos con $|k-l| \leq 1$, no hace falta ir hasta $K$ grande.

Nivel 3: matriz completa
Calculás $\Sigma(\theta)$ entera ($n \times n$), invertís, aplicás la fórmula de FIM gaussiana exacta. Caro si $n$ es grande, pero factible para $n$ moderados.

Recomendación concreta
Empezá con nivel 1: agregá el factor AR(1) a tu lambda de GFI:


auto lam_gfi = [](const auto& evo_i) {
    auto gfi_iid = sqr_X<true>(derivative(get<y_mean>(evo_i))()) 
                 * (1.0 / primitive(get<y_var>(evo_i))()) + ...;
    
    // corrección por autocorrelación AR(1):
    // rho = exp(-Δt/τ_kinetic), factor = (1-rho)/(1+rho)
    double tau_kin = estimate_tau_from_predictions(evo_i);
    double rho = std::exp(-dt / tau_kin);
    double ar1_factor = (1.0 - rho) / (1.0 + rho);
    
    return parameter_spd_payload(gfi_iid * ar1_factor, ...);
};
Esto te da el DCC $\Delta t$-invariante analíticamente, sin bootstrap, con una línea más. Si el resultado se parece al sandwich que calculás por bootstrap, tenés la corrección. Si no, pasás al nivel 2.

La clave conceptual: la covarianza temporal no está en los diagonales de $\Sigma$, está en los off-diagonals. Tu GFI actual los ignora. La mejora natural es incorporarlos, que se puede hacer progresivamente (AR(1) → retrasos finitos → matriz completa) según cuánta precisión querés.

¿Querés que te escriba el nivel 1 completo (cómo estimar $\tau$ efectivo del estado del filtro, cómo aplicarlo a cada GFI per-sample)? O preferís empezar con el nivel 2 que es más explícito?

Tenés razón, algo genérico. Pensemos estructuralmente.

La observación clave que hace esto genérico
La fórmula de Fisher correcta para un proceso de Markov oculto con observaciones Gaussianas, independientemente del esquema cinético, se factoriza según la descomposición de la verosimilitud en predicciones condicionales:

$$
\log L(\theta) ;=; \sum_k \log p(y_k \mid y_{1:k-1}, \theta).
$$

Esta es siempre exacta (teoría de la probabilidad condicional, no un supuesto). Cada $p(y_k \mid y_{1:k-1}, \theta)$ es la densidad predictiva del filtro: Gaussiana con media $\hat\mu_k$ y varianza $\hat\sigma_k^2$ condicionadas en la historia.

Si derivás Fisher de esta descomposición:

$$
I(\theta) ;=; \sum_k \left[\frac{(\partial_\theta \hat\mu_k)(\partial_\theta \hat\mu_k)^\top}{\hat\sigma_k^2} + \frac{(\partial_\theta \hat\sigma_k^2)(\partial_\theta \hat\sigma_k^2)^\top}{2\hat\sigma_k^4}\right].
$$

Es la misma forma de tu GFI — pero con $\hat\mu_k, \hat\sigma_k^2$ condicionales, no marginales.

No hay términos off-diagonal que expandir. La dependencia temporal ya está absorbida en los condicionales.

Qué diferencia hace
Algoritmo	$\mu_k, \sigma_k^2$ usados	¿GFI correcta?
Macro R / IR / IRT	Condicionales $\hat\mu_k, \hat\sigma_k^2$ del filtro recursivo	Sí, exacta — FIM real del proceso
Macro NR	Marginales $\mu_k^{\text{eq}}, \sigma_k^{2,\text{eq}}$ (ignoran historia)	No — es pseudo-FIM, necesita sandwich
Esto es genérico: aplica a cualquier $Q(\theta)$, cualquier número de estados, cualquier protocolo. No depende del esquema cinético concreto.

La teoría subyacente: descomposición innovations (Kolmogorov, Cramér). Para cualquier proceso Gaussiano estacionario o no, la verosimilitud exacta se factoriza en predicciones un-paso-adelante, y los "innovations" $y_k - \hat\mu_k$ son incorrelacionados entre sí bajo $\hat\Sigma_k = \hat\sigma_k^2$. Es decir, el filtro transforma un problema de $n$ muestras correlacionadas en $n$ muestras independientes con la varianza condicional correcta.

Lo que significa para tu pipeline
Ya tenés la mejora general si usás IR/MR/IRT en lugar de NR. No necesitás expandir más términos; necesitás que $\hat\mu_k, \hat\sigma_k^2$ en tu lambda de GFI sean los condicionales del filtro y no los marginales del modelo.

Revisar: en tu código, cuando corrés macro_IR, ¿qué está leyendo get<y_mean>(evo_i) y get<y_var>(evo_i) — las predicciones condicionales del filtro recursivo, o las medias/varianzas marginales del modelo?

Si son las condicionales del filtro, tu GFI con IR ya es la FIM exacta. La GFI con NR es la aproximación pseudo-likelihood, con todas las patologías conocidas.

Si son las marginales (el modelo evalúa $\mu, \sigma^2$ independientemente del estado del filtro), entonces falta la corrección — y la corrección genérica no es "expandir términos" sino usar las predicciones condicionales.

La propuesta estructural
La mejora de tu GFI no es algebraica (más términos) sino cambiar qué datos alimentás a la fórmula:

En vez de usar $\mu, \sigma^2$ del modelo marginal, usar $\hat\mu_k, \hat\sigma_k^2$ del filtro recursivo (Macro IR / Kalman bank para CTMC oculto).
Las derivadas $\partial_\theta \hat\mu_k, \partial_\theta \hat\sigma_k^2$ se propagan vía el filtro derivado (forward-mode autodiff a través del Kalman-like update).
La fórmula es la misma que ya tenés. No necesita cambios.
Lo que cambia es qué expresión de $\mu_k$ y $\sigma_k^2$ entra al lambda — las predictivas filtradas, no las marginales del modelo.

Verificación rápida que podés hacer
Comparar, para una celda bien condicionada, la GFI computada con IR contra el sandwich $H^{-1}{\text{NR}} J H^{-1}{\text{NR}}$ de bootstrap sobre NR. Deberían coincidir en el límite $T \to \infty$. Si coinciden, confirma que tu IR ya está computando la FIM correcta sin necesidad de bootstrap y sin expandir términos. La "mejora genérica" ya está dentro de tu código — solo no la estabas reconociendo como tal.

¿Podés verificar rápido en el código qué fuente tiene get<y_mean>(evo_i) para IR vs NR? Ese es el punto crítico.

podes repetir todo el texto de una forma que vea las formulas (mi cerebro no procesa latex)
Listo, mismas ideas sin LaTeX.

La observación clave
Para cualquier proceso de Markov oculto con observaciones gaussianas — cualquier esquema cinético, no solo 2-state — la log-verosimilitud exacta se escribe SIEMPRE como:


  log L(θ)  =  sum over k of  log p( y_k | y_1,...,y_{k-1}, θ )
Esto no es una aproximación — es la regla del producto de probabilidades condicionales. Cada término es la densidad predictiva un-paso-adelante del filtro: gaussiana con:


  media  μ̂_k  = E[y_k | y_1,...,y_{k-1}, θ]      ← predicción condicional
  varianza σ̂²_k = Var[y_k | y_1,...,y_{k-1}, θ]   ← varianza condicional
Punto clave: estos μ̂_k, σ̂²_k son diferentes de los μ, σ² marginales (las medias/varianzas de equilibrio del modelo). Los condicionales incorporan la historia y ya descontaron la correlación temporal con las muestras previas.

La fórmula de Fisher que sale de esa descomposición
Derivando la log-verosimilitud respecto a θ y cuadrando:


              ┌                                                          ┐
              │  (∂μ̂_k/∂θ)(∂μ̂_k/∂θ)ᵀ        (∂σ̂²_k/∂θ)(∂σ̂²_k/∂θ)ᵀ   │
  I(θ) = Σ_k  │  ─────────────────────  +  ───────────────────────────   │
              │         σ̂²_k                      2·(σ̂²_k)²              │
              └                                                          ┘
Es exactamente la forma de tu GFI, pero con μ̂_k, σ̂²_k condicionales en vez de marginales.

No hay términos off-diagonal que expandir. La dependencia temporal ya está absorbida en los condicionales. Lo que vos creías que había que "agregar" como corrección ya está incluido, si usás los condicionales correctos.

Qué diferencia hace algoritmo por algoritmo

Algoritmo     | μ_k, σ²_k usados                     | GFI correcta?
──────────────┼──────────────────────────────────────┼────────────────────
Macro R/IR/   | condicionales μ̂_k, σ̂²_k              | SÍ, exacta
IRT           | del filtro recursivo                 | (es la FIM real)
──────────────┼──────────────────────────────────────┼────────────────────
Macro NR      | marginales μ_k^eq, σ²_k^eq           | NO, es pseudo-FIM;
              | (ignoran historia)                   | necesita sandwich
Esto es genérico: aplica a cualquier Q(θ), cualquier número de estados, cualquier protocolo. No depende del esquema cinético concreto.

Teoría detrás
Se llama descomposición de innovations (Kolmogorov, Cramér). El filtro recursivo hace:


  innovation_k  =  y_k  −  μ̂_k
Teorema: las innovations {innovation_1, innovation_2, ...} son incorrelacionadas entre sí (bajo las varianzas condicionales σ̂²_k). Es decir:


  Cov(innovation_k, innovation_l)  =  σ̂²_k · δ_{k,l}
En otras palabras, el filtro transforma n muestras correlacionadas en n muestras "independientes" con las varianzas condicionales correctas. Por eso la FIM sale como suma simple sobre k — como si fueran iid — pero usando los valores condicionales.

Esto es un teorema general: vale para cualquier proceso Gaussiano (estacionario o no, cadenas de Markov, HMMs, difusiones). El filtro es el que convierte el problema correlacionado en independiente.

Qué significa esto para tu código
La mejora de tu GFI no es algebraica (agregar términos) sino qué alimentás a la fórmula:


  Actual (posiblemente NR):
    lam_gfi lee get<y_mean>(evo_i) y get<y_var>(evo_i)
    → ¿estos son marginales del modelo o condicionales del filtro?

  Correcto (IR/MR/IRT):
    get<y_mean>, get<y_var> deberían ser μ̂_k, σ̂²_k condicionales
    del filtro recursivo.
    Sus derivadas ∂μ̂_k/∂θ, ∂σ̂²_k/∂θ se propagan por autodiff
    a través del update recursivo (tipo Kalman).
La fórmula matemática que escribís es la misma; lo que cambia es cuál predicción (condicional vs marginal) va adentro.

Propuesta estructural
Verificar primero: cuando corrés Macro IR, ¿las y_mean y y_var que lee tu lambda son condicionales (output del filtro recursivo) o marginales (evaluación del modelo)?

Si son condicionales: tu GFI con IR ya es la FIM exacta. No falta nada. La "mejora" ya existe en tu código; solo no la estabas reconociendo como tal. El bootstrap sandwich con NR converge a esto en el límite T → ∞, pero IR lo da directamente sin bootstrap.

Si son marginales: la corrección genérica no es "expandir términos en una serie" — es usar las predicciones condicionales del filtro. Eso requiere cambiar qué datos alimentan el lambda.

Test de verificación que podés hacer
Comparar, para una celda bien condicionada (por ejemplo τ/Δt ≈ 1, sin degeneración):


  ¿ GFI(IR)  ≈  H⁻¹(NR) · J(NR) · H⁻¹(NR) ?
    |                 |
    sin bootstrap     con bootstrap (sandwich)
Si son iguales (en el límite T grande), confirma que tu IR ya computa la FIM correcta — la información exacta está en el filtro recursivo, no en un término nuevo a agregar. Y entonces entendés que NO hay que expandir más términos; la generalización está en el enfoque recursivo del filtro, no en la fórmula algebraica.

¿Podés revisar qué devuelven get<y_mean>(evo_i) y get<y_var>(evo_i) para el algoritmo IR vs NR? Ese es el punto crítico: si IR devuelve los condicionales del filtro (lo esperado), ya tenés la mejora. Si no, ahí está el bug/limitación a corregir.


Verdict: No bug — the ns-invariance is mathematically expected.
At fixed interval duration dt, varying ns (number of sub-samples per interval) does not change the per-interval likelihood, and this is correct, not a code error. The "perfect compensation mechanism" you sensed is real and falls straight out of the formula.

The one place ns enters: electronic noise
At legacy/qmodel.h:4432-4434 (non-recursive) and legacy/qmodel.h:4553-4555 (recursive):


auto e = get<Current_Noise>(m).value() * fs / get<number_of_samples>(t_Qdt).value()
       + get<Pink_Noise>(m).value()
       + get<Proportional_Noise>(m).value() * abs(y - r_y_mean());
Since dt = ns/fs (set in calc_Qdt* at e.g. legacy/qmodel.h:3371, legacy/qmodel.h:3391, legacy/qmodel.h:3416, legacy/qmodel.h:3477):

$$ e \;=\; \frac{S_e \, f_s}{n_s} + S_\mathrm{pink} + \alpha \, |y-\mu| \;=\; \frac{S_e}{\Delta t} + \ldots $$

At fixed dt, fs/ns = 1/dt — so e is independent of ns. This is just white-noise integration: averaging ns samples over a fixed window dt gives variance S_e·fs/ns = S_e/dt. Correct.

Channel noise has no ns-dependence at all
The channel contribution N·gSg (+ N·ms when variance::value) is built from gmean_i, gtotal_ij, gmean_ij, gvar_i at legacy/qmodel.h:4430 / legacy/qmodel.h:4539-4560. Those quantities come from Qdt, computed in calc_Qdt_eig legacy/qmodel.h:2511-2598. Reading that function carefully:

Everything is built from v_ladt = lambda * dt (qmodel.h:2519) and the integrals E2(λᵢdt,λⱼdt,…), E3(…) — dt only, no ns.
ns is only stored into the Qdt record at qmodel.h:2592 for later use by e. It does not enter any arithmetic producing gmean_i, gtotal_ij, gmean_ij, gsqr_i, gvar_i, etc.
The Taylor path is consistent: it sub-splits for numerical stability (qmodel.h:3063-3090), carries fractional sub_ns = ns·scale, but Qn_to_Qdt at qmodel.h:3310-3333 divides back out (r_gtotal_ij = PG_n/n, r_gtotal_sqr_ij = PGG_n · 2/n²), so the final Qdt quantities are again pure dt-functions.

The full likelihood formula (MacroR2, averaging=2, variance=V)

$$
\begin{aligned}
\mu_k  &= N \, \pi_k^{\top} \, \bar g(\Delta t) + b \\[2pt]
\sigma_k^2 &= \underbrace{\tfrac{S_e}{\Delta t} + S_{\text{pink}} + \alpha \, |\Delta y|}_{e(\Delta t)}
           \;+\; N \, g^{\top} \Sigma_k \, g
           \;+\; N \, \pi_k^{\top} \, \overline{\mathrm{var}}(g;\Delta t) \\[2pt]
\log L_k &= -\tfrac{1}{2} \log\!\bigl(2\pi \sigma_k^2\bigr) \;-\; \tfrac{1}{2} \, \frac{(y_k - \mu_k)^2}{\sigma_k^2}
\end{aligned}
$$

Every single term depends on dt only; ns cancels out. Therefore dlogL/dθ, d²logL/dθ², the Gaussian Fisher information per sample (src/core/likelihood.cpp:1499-1506)


gfi = sqr_X<true>(d y_mean/dθ) / y_var
    + sqr_X<true>(d y_var /dθ) / (2·y_var²);
and its sum H, its inverse Fisher_Covariance, and Information_Distortion_Matrix = W_H^{-1} · J · W_H^{-1} (src/core/likelihood.cpp:1614-1619) are by construction independent of ns at fixed dt. This is exactly the behaviour you observed.

Why this is the right physics (and the Poisson analogy)
MacroIR models the observation as the time-average of the continuous current over [0, dt]. That average is a single random variable whose first two moments are fully determined by dt and the model — sub-sampling the same continuous signal more finely within [0,dt] and then averaging gives you the same sufficient statistic. There is literally no new information in more sub-samples per interval, any more than 10-minute counts give extra info over the 1-hour count in your Poisson example.

The "compensation" you couldn't formulate is a tautology once you note that, with dt fixed:

white-noise contribution scales as 1/ns per sub-sample, but then you average ns of them → net variance ∝ 1/dt, no ns;
channel contribution never referenced ns to begin with.
What's not covered by this invariance (and where your MR vs IR comparison lives)
If you vary dt itself ("fracciones de tau"), Qdt changes and so does info-per-interval. For fixed total recording time T, number of intervals is T/dt, and the cumulative Fisher information scales roughly as T · I(dt)/dt. Whether that is ~constant across dt depends on how I(dt)/dt behaves for your model; for simple CO kinetics it's plausibly weakly dependent — that is also not a bug.
MacroR vs MacroIR differ only through averaging::value: MacroIR (=2) uses time-averaged gmean_i, gtotal_ij, gmean_ij from Qdt integrals; MacroR (=0 with instantaneous sample) uses g directly at the end of the step. So the MR→IR improvement shows up only when the interval is not short relative to τ — which matches what you're seeing in the cross-correlation plot.
Suggestions
The result is safe to put in the paper. You can state it as: "for models whose observation is the interval-averaged current, the per-interval Fisher information depends only on dt and the model, not on sub-sampling density, by sufficiency of the time average."
If you want me to write the short derivation block (~half a page, including the white-noise averaging step and the gtotal_ij-from-E2(λdt,λdt) line so any reviewer can verify), I can draft it into papers/macroir-elife-2025/docs/ — just say the word.
