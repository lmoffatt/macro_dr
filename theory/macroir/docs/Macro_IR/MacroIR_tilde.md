# MacroIR: Interval-Based Bayesian Updating Using Boundary-State Lifting

> **Revisada contra el código el 2026-08-09.** Lo que esta nota define y que no está escrito
> en ningún otro lado es el **operador**: el esquema lift–modulate–collapse de la §3 y la
> direccionalidad de la §3.4 (el tilde colapsa el índice **inicial** i₀; la contracción sobre
> i_t se obtiene transponiendo los pesos, sin operador hacia atrás). Eso es canónico y se usa
> tal cual.
>
> Las **ecuaciones** no lo eran. Tres estaban mal y quedaron corregidas contra
> `legacy/qmodel.h:4557-4645`, cada una marcada en su lugar: el signo de `diag(μ)` en la
> propagación de la covarianza (§3.2 y §4.1), la P que faltaba en el tilde vectorial (§3.3) y
> la P que sobraba en el escalar (§3.3). Dos más quedan **marcadas y sin tocar** en la §4.2,
> porque para reescribirlas hay que leer el bloque de actualización entero y no solo sus
> fórmulas. Verificación en `papers/1_method/decisions/recompute/mr_vs_ir_boxes.py`.
>
> Regla que gobierna esta revisión, y es la del repo: **el código productor decide, un
> documento de derivación no** (`02_theory.tex:756-764`). El paper existe justamente para
> mostrar que ese código reporta la incertidumbre que entrega, así que donde la nota y el
> código discrepan, la nota es la que está mal.
>
> Para la derivación primaria seguir usando `macroir_derivation.tex`; para qué separa a MR de
> IR, `../../notes/mr_vs_ir_from_macror.md`.

MacroIR performs Bayesian updating of a macroscopic state-occupancy distribution for an ensemble of Markov channels based on **interval-averaged** observations.  
The key difficulty is that an interval-averaged current depends on the entire hidden *trajectory* of each channel over \([0,t]\). MacroIR resolves this through:

\, - **Boundary states** \((i_0,i_t)\), encoding start–end microscopic states.
\, - A **unified lift–modulate–collapse transformation**, called the **tilde operator**, that embeds detailed boundary-state information into tractable macroscopic updates.

This document defines:

1. The interval problem and the boundary-state representation  
2. Predictive mean and variance  
3. The generalized tilde operator and its specializations  
4. The MacroIR interval update  
5. Conceptual summary


# 1. Interval Problem and Boundary-State Geometry

## 1.1 Macroscopic state variables

For an ensemble of \(N_{\text{ch}}\) independent channels with \(K\) microscopic states, let
\[
\boldsymbol{\mu}_0 \in \mathbb{R}^K,
\qquad
\boldsymbol{\Sigma}_0 \in \mathbb{R}^{K\times K}
\]
be the prior mean and covariance at time \(t=0\).

> **QUÉ OBJETO SON μ Y Σ, declarado acá porque no declararlo es la fuente de confusión más cara
> de toda esta familia (2026-08-09).** Con \(\mathbf{N}\) los conteos por estado sobre el
> ensamble, lo que el filtro arrastra es
> \[
> \boldsymbol{\mu} = \mathbb{E}[\mathbf{N}]/N_{\text{ch}},
> \qquad
> \boldsymbol{\Sigma} = \operatorname{Cov}(\mathbf{N})/N_{\text{ch}} .
> \]
>
> **Circulan dos descripciones y las dos son correctas, pero en momentos distintos.** Al
> inicializar, los canales son independientes y \(\Sigma\) coincide con la covarianza del
> indicador one-hot de **un** canal, \(\operatorname{diag}(\mu)-\mu\mu^{T}\), que es como la
> describe el código (`legacy/qmodel.h:867-874`, con su chequeo: estado determinístico \(\Rightarrow\)
> covarianza cero). **Después de la primera actualización dejan de coincidir**, porque condicionar
> sobre una observación compartida correlaciona los canales: desde ahí lo que se arrastra es
> \(\operatorname{Cov}(\mathbf{N})/N_{\text{ch}}\) y ningún canal tiene esa covarianza. La
> descripción de `macroir_derivation.tex:73-95` ("covarianzas totales divididas por
> \(N_{\text{ch}}\)") es la que sobrevive a la recursión; la del código describe la
> inicialización. Ninguna de las dos es el error: el error es no decir cuál momento describe cada
> una.
>
> **De ahí salen todos los \(N_{\text{ch}}\) de las ecuaciones**, y no de una convención de
> escala: lo observado es una cantidad poblacional y el estado se lleva por canal, así que el
> condicionamiento se hace sobre \(\operatorname{Cov}(\mathbf{N}) = N_{\text{ch}}\Sigma\) contra
> \(\operatorname{Var}(\bar y)\). El término de rango uno se lleva el conteo dos veces y volver a
> la escala por canal devuelve uno, dejando un factor; la corrección de la media se lo lleva una
> vez y lo devuelve, dejando ninguno.
>
> Verificado: la fórmula del downdate que corre coincide exactamente con
> \(\operatorname{Cov}(\mathbf{N}\mid \bar y)/N_{\text{ch}}\) y **no** con la covarianza posterior
> de un canal solo. Esto contesta el ítem que `02_theory.tex:911-915` dejó marcado como
> STILL OPEN.

Microscopic state evolution is Markovian with generator \(\mathbf{Q}\) and transition matrix:
\[
\mathbf{P}(t) = \exp(\mathbf{Q}t),
\qquad
P_{i_0\to i_t}(t) = [\mathbf{P}(t)]_{i_0,i_t}.
\]

## 1.2 Boundary states \((i_0,i_t)\)

A **boundary state** specifies a channel’s state at both ends of the interval \([0,t]\):
\[
(i_0,i_t) .
\]

For each boundary pair, precompute:
\[
\overline{\Gamma}_{i_0\to i_t}
\qquad
\text{(mean interval-averaged current)},
\]
\[
\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})
\qquad
\text{(interval-averaged current variance)}.
\]

These encode the microscopic model’s behavior over the interval.

To obtain a start-state–indexed conditional mean:
\[
(\overline{\gamma}_0)_{i_0} =
\sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}.
\]


# 2. Predictive Mean and Variance of the Interval-Averaged Observation

## 2.1 Predictive mean

The predicted macroscopic interval-averaged current is:
\[
\overline{y}^{\mathrm{pred}}_{0\to t}
\, =
N_{\text{ch}}\,
\boldsymbol{\mu}_0\cdot \overline{\boldsymbol{\gamma}}_0 .
\]

## 2.2 Predictive variance

Measurement noise:
\[
\epsilon^2_{0\to t} \qquad = \frac{\epsilon^2}{t} + \nu^2.
\]

Intrinsic channel variability consists of:

1. Uncertainty in the initial ensemble state,  
2. Variance of the interval current conditioned on boundary states.

Define:
\[
(\sigma^2_{\overline{\gamma}_0})_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,
\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t}).
\]

The full predictive variance is:
\[
\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}
\, =
\epsilon^2_{0\to t}
\, +
N_{\text{ch}}
\,\widetilde{\gamma^{T}\Sigma\gamma}
\, +
N_{\text{ch}}
\sum_{i_0}
(\boldsymbol{\mu}_0)_{i_0}\,
(\sigma^2_{\overline{\gamma}_0})_{i_0}.
\]

The central term \(\widetilde{\gamma^{T}\Sigma\gamma}\) is a **triple-tilde** scalar defined below.


# 3. The Tilde Operator: Generalized Lift–Modulate–Collapse

The **tilde operator** is a three-stage transformation applied to any state-space object \(X\):

1. **Lift**: reinterpret \(X\) in boundary-state coordinates.  
2. **Modulate**: apply interval-conditioned microscopic weights (e.g.  
   \(\overline{\Gamma}_{i_0\to i_t},\overline{\mathbf{V}},P(t)\)).  
3. **Collapse**: sum over the initial boundary index \(i_0\), yielding a macroscopic object.

Formally, a tilde always has the form:
\[
\widetilde{X}_{i_t}
\, =
\sum_{i_0}
\text{Lift}(X)_{i_0,i_t}
\cdot
\text{Modulate}(i_0,i_t).
\]

The operator depends on the interval:
\[
\widetilde{(\cdot)} \equiv \widetilde{(\cdot)}_{0\to t},
\]
though the interval is usually clear from context.

## 3.1 First-order tilde: mean and boundary expectations (special cases)

Mean propagation is a first-order tilde:
\[
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, =
\boldsymbol{\mu}_0\,\mathbf{P}(t)
\, =
\widetilde{\boldsymbol{\mu}_0}.
\]

Similarly, the start-state boundary-average
\[
(\overline{\gamma}_0)_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}
\]
is a *backward* first-order tilde (or forward tilde with time indices transposed).

These unify neatly under the lift–modulate–collapse schema.

## 3.2 Second-order tilde: covariance propagation (special case)

Covariance propagation follows exactly the same schema applied to a second-order object:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\, =
\mathbf{P}(t)^T\,
(\boldsymbol{\Sigma}_0 - \mathrm{diag}(\boldsymbol{\mu}_0))\,
\mathbf{P}(t)
\, +
\mathrm{diag}(\boldsymbol{\mu}^{\mathrm{prior}}(t)),
\]
which is a first-second–order tilde transformation. Los dos términos son las dos fuentes de
dispersión: el primero lleva la incertidumbre del prior a través de la ventana, una vez por
índice, y el segundo **regenera** la parte multinomial de un solo canal en el instante final,
que es la que el `- diag(μ₀)` de adentro sacó para no contarla dos veces.

> **CORREGIDO 2026-08-09: el segundo término iba con menos.** No es una convención distinta,
> es una errata, y aparecía dos veces (acá y en la §4.1, una sola errata copiada). Con el menos
> la matriz no es una covarianza: en el caso de dos estados simétrico da autovalores
> $[-1, -0.5]$ y traza $-1.5$ donde tiene que dar $0.5$. Con el más coincide exactamente con el
> margen sobre $i_0$ de $\Sigma^{\mathrm{bnd}}$ (1.5e-16) y con `legacy/qmodel.h:4655`
> (`sigma_pre = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P())`), que es la forma que corre.

## 3.3 True MacroIR tilde operators (second and third order)

For inference, we require tilde applied to second-order contractions involving the covariance.

### Double-tilde (vector)
For \(u\in\mathbb{R}^K\) with boundary weights \(\overline{\mathbf{W}}\):
\[
\widetilde{u^{T}\Sigma}
\, =
\overline{\mathbf{w}}_0^{T}\,
(\boldsymbol{\Sigma}_0 - \operatorname{diag}(\boldsymbol{\mu}_0))\,
\mathbf{P}(t)
\, +
(\overline{\mathbf{W}}\circ \mathbf{P}(t))^{T}\boldsymbol{\mu}_0,
\]
where
\[
(\overline{\mathbf{w}}_0)_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,\overline{W}_{i_0\to i_t}.
\]

> **CORREGIDO 2026-08-09: al primer término le faltaba la \(\mathbf{P}(t)\) del final.** Tiene
> que estar, y por una razón que es la del propio operador: el resultado se indexa por \(i_t\),
> así que el término que entra por el prior tiene que cruzar la ventana para llegar a ese
> índice. Sale de colapsar \(\Sigma^{\mathrm{bnd}}\overline{\mathbf{W}}\) sobre \(i_0\).
> Verificado contra `legacy/qmodel.h:4641`
> (`TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij()`) y contra Monte Carlo.
> **Cuidado al chequearlo**: con prior estacionario y multinomial las dos formas coinciden por
> accidente, porque \(\Sigma_0 - \operatorname{diag}(\mu_0) = -\mu_0\mu_0^{T}\) y
> \(\mathbf{P}^{T}\mu_0 = \mu_0\). Fuera de ese caso discrepan en los 12 de 12 casos probados,
> con cambio de signo incluido.

### Triple-tilde (scalar)
\[
\widetilde{u^{T}\Sigma w}
\, =
\overline{\mathbf{u}}_0^{T}
(\boldsymbol{\Sigma}_0 - \operatorname{diag}(\boldsymbol{\mu}_0))
\overline{\mathbf{w}}_0
+
\boldsymbol{\mu}_0^{T}
\left[
\overline{\mathbf{U}}
\circ
\overline{\mathbf{W}}
\circ
\mathbf{P}(t)
\right]\mathbf{1}.
\]

> **CORREGIDO 2026-08-09: la \(\mathbf{P}(t)\) aparecía dos veces.** El segundo término es
> \(\sum_{i_0,i_t}\mu^{\mathrm{bnd}}_{i_0 i_t}\,\overline{U}\,\overline{W}\) con
> \(\mu^{\mathrm{bnd}} = \mu_0 P\), o sea **una** P, no el producto de dos objetos que ya la
> llevan cada uno. Con la forma anterior el término da negativo, que para una contribución a una
> varianza es imposible: en dos estados a \(\Delta k_{\mathrm{off}} = 1\) daba \(-0.076\) contra
> \(+0.082\). Verificado contra `legacy/qmodel.h:4574`
> (`elemMult(t_gtotal_ij(), t_gmean_ij()) * u`, que es \(\overline{\Gamma}\circ\mathbf{P}\)
> multiplicado por \(\overline{\Gamma}\), una sola P).

For current-based likelihoods:
\, - Triple tilde: \(\widetilde{\gamma^{T}\Sigma\gamma}\).  
\, - Under TaylorIR: mixed forms \(\widetilde{\gamma^{T}\Sigma v}\).


## 3.4 Directionality

Tilde is defined as contraction over the **initial** boundary index \(i_0\).  
If contraction over \(i_t\) is needed, we apply tilde to:
\[
\overline{\mathbf{W}}^{T},
\qquad
\mathbf{P}(t)^{T}.
\]

No independent backward-tilde operator is necessary.


# 4. MacroIR Interval Update

## 4.1 Propagation step

Mean:
\[
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, =
\boldsymbol{\mu}_0 \mathbf{P}(t).
\]

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\ =
\mathbf{P}(t)^T
(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0)
\mathbf{P}(t)
\, +
\operatorname{diag}(\boldsymbol{\mu}^{\mathrm{prior}}(t)).
\]

> **CORREGIDO 2026-08-09**, misma errata que en la §3.2 y por las mismas razones. Era la misma
> fórmula copiada, así que el signo estaba mal en los dos lugares.

## 4.2 Measurement update

Let
\[
\widetilde{\gamma^{T}\Sigma} \in \mathbb{R}^K,
\qquad
\sigma^2 \qquad = \sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}.
\]

Define the innovation:
\[
\delta
\ =
\overline{y}^{\mathrm{obs}}_{0\to t}
\, -
\overline{y}^{\mathrm{pred}}_{0\to t}.
\]

Mean update:
\[
\boldsymbol{\mu}^{\mathrm{post}}(t)
\, =
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, +
\frac{1}{\sigma^2}
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\,
\widetilde{\gamma^{T}\Sigma}
\,\delta.
\]

Covariance update:
\[
\boldsymbol{\Sigma}^{\mathrm{post}}(t)
\, =
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\, -
\frac{1}{\sigma^2}
(\widetilde{\gamma^{T}\Sigma})^{T}
(\widetilde{\gamma^{T}\Sigma}).
\]

Esto es el condicionamiento gaussiano del conjunto (estado, corriente) sobre el valor grabado,
en el espacio macroscópico de K estados. No hace falta ninguna ecuación de medición del tipo
\(y = Hx + v\): los tres bloques del conjunto se calculan de la física del intervalo y Bayes
hace el resto.

> **DOS ECUACIONES DE ESTA SUBSECCIÓN QUEDAN MARCADAS Y SIN TOCAR, 2026-08-09.** Ninguna de las
> dos coincide con lo que corre, pero para reescribirlas hay que leer el bloque de actualización
> completo, con su coeficiente de confianza incluido, y no alcanza con comparar fórmulas. Las
> discrepancias, para que quien lo haga sepa qué busca:
>
> 1. **La media lleva un \(\Sigma^{\mathrm{prop}}\) de más.** Acá está escrito
>    \(\mu^{\mathrm{post}} = \mu^{\mathrm{prior}} + \sigma^{-2}\,\Sigma^{\mathrm{prop}}\,
>    \widetilde{\gamma^{T}\Sigma}\,\delta\), pero el tilde **ya contiene** la covarianza, así que
>    multiplicarla otra vez la cuenta dos veces. En Theory la corrección es
>    \(\mathbf{g}\,\delta/\sigma^2\) a secas (`02_theory.tex`, Eq. mu-post).
> 2. **Al downdate le falta el factor \(N_{\mathrm{ch}}\).** Theory lo lleva
>    (Eq. sig-post) y el código también (`XTX(gS)` pesado por `N / r_y_var()`). Sale de que el
>    estado se lleva en escala por canal mientras que lo observado es una cantidad poblacional,
>    o sea de que el condicionamiento se hace sobre \(\operatorname{Cov}(\mathbf{N})\) contra
>    \(\operatorname{Var}(\bar y)\); ver la declaración de convención de la §1.1, que es donde
>    esto se decide. **No decir que \(\mu\) y \(\Sigma\) son "de un canal"**: lo son al
>    inicializar y dejan de serlo en la primera actualización.
>
> Además el paso real está amortiguado por el coeficiente de confianza que mantiene a \(\mu\) en
> el símplex (`06_methods.tex`, los tres resguardos numéricos), que acá no figura.


# 5. Conceptual Summary

MacroIR reduces a trajectory-dependent microscopic observation to a macroscopic Bayesian update through the following structure:

1. **Boundary states** \((i_0,i_t)\) encode interval-conditioned behavior without explicit trajectory enumeration.  
2. **Boundary-conditioned matrices** \(\overline{\Gamma}\), \(\overline{V}\), and \(\mathbf{P}(t)\)  
   capture microscopic physics over \([0,t]\).  
3. The **tilde operator** is a unified lift–modulate–collapse transform:  
   - First-order tilde covers mean and boundary expectations,  
   - Second-order tilde covers covariance propagation,  
   - Second/third-order tilde contractions yield \(\widetilde{\gamma^{T}\Sigma}\) and  
     \(\widetilde{\gamma^{T}\Sigma\gamma}\) for inference.  
4. The update is the **Gaussian conditioning** of the pair (state, current) on the recorded
   value, in the \(K\)-dimensional macroscopic state space.  
5. This avoids constructing the full \(K^2\times K^2\) boundary covariance while preserving all necessary microscopic detail.

The tilde operator is the core mechanism that makes MacroIR efficient, expressive, and compatible with microscopic kinetics without explicit path integration.

Y el punto 5 no es una optimización: el operador **entra y sale** del espacio de K estados, con
el de \(K^2\) adentro suyo, así que el objeto grande no se construye por construcción y no por
astucia de implementación. Eso es lo que hace que el esquema de la §3 sea la definición y no una
descripción.

## 6. Lo que el orden de las tres etapas decide

Colapsar **después** de contraer sale gratis: es el marginal del posterior. Colapsar **antes**
cuesta, y cuesta exactamente lo que la corriente grabada habría dicho sobre el índice que se
colapsó. Los dos miembros de la familia que condicionan al intervalo son las dos maneras de
ordenar las mismas tres etapas:

    IR:   lift  →  modulate  →  collapse        colapsa i_0, al final
    MR:   lift  →  collapse  →  modulate        colapsa i_t, antes de contraer

MR ni siquiera necesita levantar la covarianza, porque levantarla y colapsarla sin contraer
devuelve \(\Sigma_0\). Las dos lecturas, su costo y su tamaño están en
`../../notes/mr_vs_ir_from_macror.md`.
