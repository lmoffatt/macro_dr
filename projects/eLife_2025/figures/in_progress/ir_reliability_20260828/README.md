# IR reliability: cotas de una variable, componentes, y los gusanos (2026-08-28)

Sesión de un día sobre "un número simple para decir cuándo IR es confiable". Todo sale de las
baterías `_G` ancladas en θ_sim (verdad) de `figures/data/{1c2ae6f,0ffbda7,87889e6,82b956f}`,
`figure_3_G_nch_*_nsim_10000_macro_IR_noise_*_battery_sim_G.csv`: 560 celdas IR,
N_ch ∈ {5..10⁴} × 12 ruidos × 7 intervalos. Cero corridas nuevas.

Convenciones: S̃ = 0.1·label (NOISE_AXIS_UNITS.md); Δ̂ = Δ·k_off = `interval_in_tau`;
matrices AGREGADAS (probit=mean, reensambladas 6×6, eigen) y no los escalares emitidos
(sesgo E[f]>f(E) cerca de I, ver memoria gaussian_fisher_distortion).

## Pipeline

Digests (leen las 560 baterías, ~1–2 min c/u, salida en `digest/`):

| script | produce | contenido |
|---|---|---|
| ir_collapse_probe.R | ir_collapse_cells.csv | GIDM total: ē, s, d_AI; zmax por celda (sim y pool) |
| ir_bias_cis.R | ir_bias_cis.csv | GDIB por parámetro con CI bootstrap + SE (diag GFC) |
| ir_bias_mahal.R | ir_bias_mahal.csv | √(b'Gb) con la GFC completa |
| ir_bias_whitened.R | ir_bias_whitened.csv | v = G^(1/2)·b, 6 componentes (raíz simétrica) |
| ir_components.R | ir_components_cells.csv | ē, s para total / sample / corr (**columna dAI = NA, bug tibble; usar 6(ē²+s²)**) |
| ir_r2_witness.R | ir_r2std_cells.csv | r̄²_std − 1, r̄_std (Sum_r2_std / n, n = 10/Δ̂) |
| ir_shape_witnesses.R | ir_score_JG_perparam.csv, ir_residual_shape.csv | J_pp, J_step_pp, G_pp por parámetro; τ_int, curtosis, memoria residual |
| ir_perdelta_all.R | ir_perdelta_all.csv | ajustes por intervalo (a, b, A, R²) para las 4 cantidades |
| ir_pair_boot.R | ir_pair_boot.csv | el par del paper (ē, s → m, a), la peor dirección max|log λ| y λmax/λmin de la matriz agregada, + cuantiles bootstrap EMITIDOS de ē, s², λmax, λmin (barras de error) |
| ir_dai_boot.R | ir_dai_boot.csv | cuantiles bootstrap emitidos del escalar afín (media sesgada arriba ×1.18 mediana; los anchos sirven) |
| ir_bias_Hb.R | ir_bias_Hb_se.csv | H·b por celda y SE delta de b'Hb desde los CI por parámetro |

Figuras (leen `digest/`, escriben `figures/`): `figure_ir_*.R`. Las que valen:

- `figure_ir_8panel_round.R` → Figure_IR_8panel_round: bias / total / sample / corr, media y desvío, cada una contra su variable redonda.
- `figure_ir_sumsq.R` → Figure_IR_sumsq: lo mismo con Σlog²λ y b'Gb (la elección final).
- `ir_perdelta_all.R` → Figure_IR_perdelta_all: amplitud y exponentes por intervalo, 4 cantidades.
- `ir_sample_perdelta.R` → Figure_IR_sample_perdelta: la componente de muestra, un ajuste por Δ̂.
- `figure_ir_ebar_anatomy.R` → Figure_IR_ebar_anatomy: por qué ē tiene dos ramas.
- `figure_ir_bias_whitened.R` → Figure_IR_bias_whitened: dirección del bias por parámetro.
- `figure_ir_sign_witness.R` (en ir_r2_witness.R) → Figure_IR_sign_witness: ē, Σv, r̄²−1.

El resto (`ir_collapse_fit/variant`, `ir_nls`, `ir_round_exponents`, `ir_sample_peak*`,
`ir_sample_floor`, `ir_sample_nls`, `figure_ir_reliability*`, `figure_ir_4panel/6panel/8panel*`)
son los pasos intermedios; se dejan porque documentan qué se probó y qué no funcionó.

## Hallazgos (estado al 2026-08-28)

VERIFICADOS celda por celda (envolventes con cero violaciones en 560 celdas):

1. **Cota de distorsión** (barras de error): s = sd log λ ≤ max(0.21·x^(−1/5), 0.045),
   x = N_ch·S̃·√Δ̂. Sin inflación neta: |ē| ≤ s en 560/560. En d_AI: ≤ 0.5·x^(−1/5).
   Factor peor-dirección en la barra ≤ exp(d_AI/2): x=1 → ×1.28, x=100 → ×1.10 (medidos: peor celda 1.49 en N=10, S̃=0.005, Δ̂=0.02; 1.16 con x≥1; 1.07 con x≥100).
   Exponentes libres (sobre piso): N −0.20, S̃ −0.24, Δ̂ −0.12; R² 0.82, dispersión ×1.14.
   EN EL PAR DEL PAPER (2026-08-28, es lo que va al suplemento; D_AI no se reporta): a = e^s ≤ exp(0.21·x^(−1/5)),
   |log m| ≤ log a en 560/560, y la PEOR DIRECCIÓN medida directamente: max|log λ_i| ≤ 0.4·x^(−1/5) (piso 0.069;
   0.38 ya da cero violaciones, 0.35 una). Factor en la información de la peor dirección: 1.5 a x=1 (1.22 en la
   barra), 1.17 a x=100, 1.78 en la peor celda (N=10, S̃=0.005, Δ̂=0.01; 1.33 en la barra). max|log λ|/s ≈ 1.74.
2. **Cota de bias**: b'Gb ≤ max(0.024·(N_ch·S̃^(3/4))^(−1), 0.003) por grabación
   (pendiente libre −0.94 ≈ −1 = TCL). Pooling: n* = 1/(b'Gb) ≈ 42·N_ch·S̃^(3/4) grabaciones.
   Per-parámetro: max_p |b_p|/SE_p ≤ max(0.8/√N_ch, piso). Ciego a Δ̂ (exponente ~0).
3. **Qué es N_ch·S̃**: = N_ch/v², v = i/σ_inst(τ) visibilidad de un canal. IR falla donde
   casi resolvés canales individuales y son pocos = frontera macro/micro. N·S̃·Δ̂ =
   saltos por muestra / visibilidad² (variable natural de la componente de correlación).
4. **Bias absoluto NO sirve**: G⁻¹s̄ explota donde i y N_ch pierden identificabilidad separada
   (errores relativos 100×–30000× a S̃ ≥ 10³, todos no significativos). Normalizar por
   Fisher (b'Gb) o por SE. |v| ≥ max_p|b_p|/SE_p siempre (Cauchy-Schwarz), mediana ×1.41.
5. **Dirección del bias rota con Δ̂**: Δ̂ ≤ 0.05 → baseline (+, 97%) y S (−); Δ̂ ≥ 0.5 →
   S (−, 97%), k_off (+, 90%), N_ch (+, 79%), i (−). k_on inmune. Norma casi independiente de Δ̂.
6. **ē tiene dos ramas con dos mecanismos** (testigos en residuos, misma batería):
   Δ̂ corto, ē>0: residuos leptocúrticos (κ +0.1..0.3) y autocorrelados (memoria +0.06..0.16),
   J/G > 0 en k_on, k_off, i, N_ch, mitad forma (J_step/G) mitad memoria (J/J_step).
   Δ̂ largo, ē<0: varianza sobre-predicha (r̄²−1 = −2..−7%), platicúrtico (κ −0.2..−0.4),
   cae en S y baseline, sin memoria. A x>100 los 6 parámetros calibrados a 0.3%.
   Matiz para Fig. 3: "score de IR blanco" vale fuera del rincón; adentro (N≤20, S̃≤0.02,
   Δ̂≤0.02) queda memoria hasta 0.2 en log-varianza en k_off.
7. **Componentes** (ē_total = ē_sample + ē_corr verificado a 1.7e-4):
   sample ∝ pico en Δ̂ = 0.1τ; corr ∝ N·S̃·Δ̂ (memoria muere linealmente con saltos por
   muestra). El √Δ̂ del total = media de Δ̂^0 y Δ̂^1. Exponentes redondos con R² −0.01:
   bias N·S̃^(3/4) pend −1/2; total N·S̃·√Δ̂ pend −1/5; sample N·S̃² pend −1/8 (CONTAMINADO,
   ver 9); corr N·S̃·Δ̂ pend −1/5.
8. **Por intervalo** (ir_perdelta_all): sólo el bias tiene exponentes constantes en Δ̂
   (a −0.9, b −0.6). Total constantes hasta 0.2τ (a −0.4, b −0.6), deriva a Δ̂ ≥ 0.5.
   Sample: exponentes ROTAN (a +0.15 → −0.85, b −1.32 → −0.24), amplitud pico 0.1τ.
   Corr: a −0.5 siempre, amplitud mínimo en 0.2τ (subida hacia τ sugerida, n=7–8).
   Ajustar por intervalo es la única forma de ver esto; R² 0.8–0.96 por intervalo.

RETRACTADOS en la misma sesión (no repetir):

- "Piso de bias 0.03–0.04 es real, techo de pooling ~350 grabaciones": NO. El piso es
  ruido de 10⁴ simulaciones (√(6/10⁴) = 0.024 en |v|, 0.015 en zmax; reparto uniforme
  entre 6 parámetros, signo aleatorio, plano en N), y las celdas "significativas" a N alto
  son el 25% = 1 − 0.95⁶ (prueba múltiple). Todos los pisos (s 0.045, ē 0.02, |v| 0.024)
  = resolución de n_sims = 10⁴.
- Ley global de la componente de muestra (N·S̃², −1/8, R² 0.52): contaminada por el piso,
  ver 9. Restar el piso tampoco (Luciano): ajustar por intervalo.

## Trampas

9. **El piso de la componente de muestra es ∝ Δ̂** (1.4e-3·Δ̂ en Σlog²λ, R² 0.994): J_step se
   estima con n_sims × n_pasos y n_pasos = T/Δ. Por eso en la cola los rojos quedan arriba de
   los azules. Cualquier estadístico construido con momentos POR PASO hereda un piso ∝ Δ̂;
   los construidos con totales por grabación (J, G, bias) tienen piso plano. Medir el piso
   por intervalo antes de leer colas.
10. `ir_components_cells.csv` columna `dAI` es NA (enmascaramiento de tibble: `tibble(s=s[2],
    dAI=s[3])` lee la columna recién creada). d_AI² = 6(ē²+s²) exacto.
11. Ancla: sim (θ_sim) para cotas de usuario; pool anula el bias (score≈0) y le borra la
    dependencia en ruido a la distorsión (R² 0.13 vs 0.61).
12. Envolventes con "cero violaciones sobre el piso": el piso hay que fijarlo como máximo de
    la región calibrada (x > 10⁴) por cantidad (y por intervalo para sample), y la constante
    ajustarla sobre celdas a > 1.5× piso, o las celdas rozando el piso a x grande la inflan.

## Qué iría al paper (orden de importancia, propuesto 2026-08-28)

1. Cotas 1 y 2 como figura de dos paneles (B y F de Figure_IR_6panel / o sumsq) + párrafo
   en Results; da número al [THRESHOLD-DJ] de la Discusión.
2. Una oración: N_ch·S̃ = N_ch/v², frontera macro/micro (hallazgo 3).
3. Methods, una oración: pisos = resolución de 10⁴ simulaciones, prueba múltiple.
4. Suplemento: anatomía de ē (6) + el matiz para la Fig. 3.
5. Suplemento: dirección del bias por intervalo (5); perfil de la componente de muestra (8).
No agregar: bias absoluto, techo de pooling, paneles de Σv (dependen de la base).

Pendiente antes de escribir 1: la frontera de región 0 de la Fig. 7 tiene pendiente −0.50 en
dlog(ruido)/dlog(N_ch) (texto vivo y caption 2026-08-26; el −0.56 que decía acá era el header
viejo de figure_7.Rmd; criterio: cruce por parámetro de 1.15 en θ_pool a Δ̂=0.1, N_ch 5–100);
la envolvente de s a Δ̂ fijo implica −0.8..−1 (curvas de nivel N·S̃ = const). Criterios y anclas
distintos, pero el paper no puede mostrar dos pendientes sin decir por qué. Además: el
[THRESHOLD-DJ] ya no existe en el tex vivo (resuelto 2026-08-04 declinando nombrar umbral; su
hueco es 05_discussion.tex:210–218); figure_6.Rmd:17 todavía apunta al marcador muerto.

## 2026-08-29: la envolvente de máximo se reemplazó por mediana + cuantil 95% (decisión de Luciano)

Por qué: la envolvente de cero violaciones es un máximo sobre 560 celdas ruidosas; la fijaba UNA celda
(N=20, S̃=1, Δ̂=1, x=20) más su ruido bootstrap, y la constante 0.21 que había puesto la inflaban celdas
apenas sobre el piso a x grande (con "sólo celdas > 1.5× piso" es 0.14; consciente del CI, 0.12). Además el
tope de la nube se aplana hacia el rincón (pendiente por década −0.15), así que una potencia tangente a x=20
quedaba ×1.7–2.6 arriba en x ≤ 0.01. Una forma saturante s ≤ 0.42(1+x/0.004)^(−1/5) pega, pero son dos
constantes y una meseta extrapolada.

Lo que va al paper (Figure 6—figure supplement 1, `paper_both/figure_6_S1.R`; en el par del paper, sin D_AI):
ley mediana con pendiente fija −1/5 (−1 para el bias) y cuantil 95% de v·x^(1/5) sobre celdas > 1.5× piso:

| | mediana | q95 | x=1 | x=100 | peor celda |
|---|---|---|---|---|---|
| a = e^s | exp(0.094·x^(−1/5)) | 0.118 | 1.10 (1.13) | 1.04 (1.05) | 1.38 |
| |ē| | 0.026 (n=32 sobre piso) | 0.037 | | | |
| peor dirección exp(max|log λ|) | exp(0.167·x^(−1/5)) | 0.212 | 1.18 (1.24); barra 1.09 (1.11) | 1.07 (1.09) | 1.78 (barra 1.33) |
| b'Hb | 0.0107/(N·S̃^(3/4)) | 0.0177 | n* = 93 (q95: 56) × N·S̃^(3/4) | | |

Barras: ē banda emitida; a y peor dirección semianchos emitidos de √Var(log λ) y log λmax centrados en el valor
de matriz (medias emitidas sesgadas ×1.18, anchos no); b'Hb ±1.96 SE delta (ir_bias_Hb.R). Digests nuevos:
ir_pair_boot.csv, ir_dai_boot.csv, ir_bias_Hb_se.csv. Tipografía: device `pdf` + Helvetica como las Figs 6 y 2 S1
(cairo mete Nimbus+DejaVu; el device png sustituye la cara: mirar SIEMPRE el PDF rasterizado).

## 2026-08-29: el piso binomial de Methods 405 NO está activo en el barrido (Luciano lo sospechó; verificado)

`is_Binomial_Approximation_valid(N, p, q, Np_min)` (legacy/qmodel.h:2693, Np_min = Binomial_magical_number =
5.0, models_simple.h:65) se consulta sólo en qmodel.h:7301 y :7363, que están en la rama `else` de
`if constexpr (!adaptive::value)` (bloque 6857–7232) de `log_Likelihood`. El barrido corre con
`adaptive_approximation = false` (ops/local/figure_3_mle_G.macroir:87), así que la prueba nunca se ejecuta y
la caída a la forma no recursiva nunca ocurre: las celdas N = 5 y N = 10 las puntuó IR recursivo. El punto
abierto "quién puntuó las celdas de 5 y 10 canales" se cierra: IR; la cima de la nube de Figure 6 S1 es
legítima. QUEDA PARA LUCIANO: Methods 405–407 ("a binomial-count floor of 5, below which ... the member falls
back to its non-recursive form") describe la constante como si actuara; la nota de verificación del
2026-08-13 verificó las oraciones de truncamiento (to_Probability en dc1d295), no esta compuerta. El paper no
menciona el flag adaptive en ningún lado vivo.
