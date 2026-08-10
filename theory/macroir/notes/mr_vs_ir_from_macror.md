---
date: 2026-08-08
status: derivación verificada; la reescritura de los sitios que dependen de ella está pendiente
scope: qué separa MacroMR de MacroIR, partiendo del boundary state, en la base de MacroR
verifica: papers/1_method/decisions/recompute/mr_vs_ir_boxes.py (36 casos a 1e-15, más Monte Carlo)
lee: legacy/qmodel.h:4557-4645 (camino vivo), :1743-1774 (los objetos de la ventana)
---

# MR contra IR, partiendo del boundary state

## Por qué existe este archivo

En el repo hay **una medición** canónica de la relación MR/IR
(`papers/1_method/figures_build_plan.md:195-245`, 2026-08-06) y **once paráfrasis** de ella.
No había derivación en ningún lado. La paráfrasis que se propagó, *"MR e IR asignan la misma
varianza y lo único que los separa es la actualización"*, es cierta de un total y falsa de todo
lo demás, y ya tuvo dos rondas de parches. Vuelve porque una medición no impide que una
paráfrasis derive. Una derivación sí.

**Vocabulario.** No se usa la palabra *gain* ni el marco de filtrado lineal. El objeto que esa
palabra nombra es una covarianza, y llamarlo ganancia lo hace sonar a parámetro del algoritmo
cuando es una propiedad del par de variables que se está condicionando. La traducción a ese
marco va al final y está sin escribir a propósito.

**El código no está mal.** Nada de acá reporta un error en el camino vivo. El código llega a
estos mismos números por la ruta contraída, que es la eficiente porque no arma nunca el objeto
de K²×K². Lo que no hace, ni tiene por qué, es exhibir el mecanismo. Eso es trabajo de esta nota.

## La base

MacroR lleva una gaussiana sobre las ocupaciones de un canal y la condiciona a lo grabado por
Bayes. Lo que condiciona es el conjunto (estado, corriente), y ese conjunto tiene **tres bloques**:

| bloque | qué es |
|---|---|
| Var(estado) | incertidumbre sobre en qué estado está el canal |
| Var(corriente) | varianza de lo que se graba: la que viene del estado, más la de la corriente dado el estado |
| Cov(estado, corriente) | cuánto informa lo grabado sobre el estado |

Eso es todo MacroR. Extenderlo al intervalo **no cambia ninguna de esas tres líneas**: cambia
sobre qué variable de estado se arman.

## El estado generalizado: por qué tiene que ser el par

La corriente promediada al intervalo no es lineal en la ocupación de ningún instante. No lo es
en el estado inicial X₀, no lo es en el final X_Δ. **Sí lo es en el par.**

> **X^bnd**: el indicador del par (i₀, i_t) de estados en los dos extremos de la ventana. Vale
> uno en el par que el canal recorrió y cero en los otros K²−1.

Sobre ese estado, un canal aporta en promedio Γ̄_{i₀→i_t}, así que la corriente vuelve a ser
lineal en la ocupación, y MacroR se le puede correr encima sin tocar nada de MacroR.

Su gaussiana sale del prior y del propagador:

- **μ^bnd_{(i₀,i_t)} = μ_{0,i₀} · P_{i₀i_t}**
- **Σ^bnd_{(i₀,i_t),(j₀,j_t)} = P_{i₀i_t} (Σ₀ − diag(μ₀))_{i₀j₀} P_{j₀j_t} + [diag(μ^bnd)]**

Dos fuentes de dispersión, y la resta no es un truco: el primer término lleva la incertidumbre
que ya tenía el prior a través de la ventana, una vez por índice; el segundo es la incertidumbre
que le queda a un canal que se sabe saliendo de i₀ sobre dónde aterriza, y es diagonal sobre los
pares porque el canal toma exactamente uno. Restar diag(μ₀) evita contar dos veces la parte
propia del prior, que el segundo término regenera a la resolución de los pares.

Verificado contra Monte Carlo del camino CTMC: μ^bnd y Σ^bnd a 7e-4 con n = 3·10⁵, que es error
de muestreo.

**Qué objeto son μ y Σ, y no decirlo genera la confusión más cara de esta familia.** Son
μ = E[N]/N_ch y Σ = Cov(N)/N_ch, con N los conteos sobre el ensamble. Coinciden con la media y la
covarianza del indicador de **un** canal al inicializar, y dejan de coincidir en la primera
actualización, porque condicionar sobre una observación compartida correlaciona los canales.
**No decir "por canal" a secas.** La declaración completa, con la verificación y con por qué
las dos descripciones que circulan son correctas en momentos distintos, está en
`../docs/Macro_IR/MacroIR_tilde.md` §1.1. Todo lo de abajo está en esa escala, y los factores
N_ch se omiten donde no cambian la comparación.

## MacroR sobre el boundary state, literal

Con Γ̄ como conductancia por estado y V̄ (la varianza de la corriente **dado el par**) cobrada
como ruido dependiente del estado:

| bloque | sobre X^bnd |
|---|---|
| media | N · μ^bndᵀ Γ̄ |
| Var(corriente) | ε² + N · **Γ̄ᵀ Σ^bnd Γ̄** + N · **μ^bndᵀ V̄** |
| Cov(estado, corriente) | N · **Σ^bnd Γ̄** |

Son las tres líneas de MacroR con Γ̄ donde iba γ y Σ^bnd donde iba Σ. Nada más.

**IR es eso, y después marginaliza i₀**, que es el índice que la ventana siguiente no necesita.

Verificado: correr MacroR así sobre los K² pares y marginalizar reproduce a 1e-16 lo que el
código computa para av=2, los tres bloques, incluido el cruzado.

## MR: la misma construcción, marginalizando antes

**MR marginaliza i_t primero y después corre MacroR.**

Marginalizar i_t de X^bnd devuelve exactamente X₀ con (μ₀, Σ₀): las sumas de μ^bnd y de Σ^bnd
sobre los índices finales colapsan a μ₀ y a Σ₀ − diag(μ₀) + diag(μ₀). Y por la ley de la
esperanza y de la varianza totales, sobre ese estado la conductancia es γ̄₀ y el ruido por estado
deja de ser V̄ y pasa a ser la varianza **total** por estado inicial:

> Var(corriente | i₀) = **E[V̄ | i₀] + Var_{i_t}[Γ̄ | i₀]**

O sea: al marginalizar, el spread del promedio entre estados finales se suma al ruido, que es
justamente por qué MR lleva la forma total y IR la residual. Después MacroR sobre X₀, y propagar
al final de la ventana, que es donde arranca la siguiente.

Verificado: reproduce av=1.

**Los dos marginalizan. La diferencia es cuál índice y cuándo.** IR marginaliza i₀ después de
condicionar; MR marginaliza i_t antes.

## El objeto que los separa

> **D_{i₀i_t} = Γ̄_{i₀→i_t} − γ̄₀_{i₀}**

el desvío de la conductancia condicional al par respecto de la condicional al inicio. No está en
el código; es la reescritura. Dos propiedades, y toda la asimetría sale de la primera:

1. **Suma cero por FILA**: Σ_{i_t} P_{i₀i_t} D_{i₀i_t} = 0 para cada i₀, por definición de γ̄₀.
2. **Perfil por COLUMNA no nulo**: a i_t fijo, sumando sobre i₀ con pesos μ₀P, no se anula. Ese
   perfil es la información sobre en qué estado terminó el canal.

Escribimos E[D²] = Σ μ_{0,i₀}P_{i₀i_t}D², y E[V̄] = Σ μ_{0,i₀}P_{i₀i_t}V̄_{i₀→i_t}.

## Bloque Var(corriente): las cuatro cajas

| | **la que viene del estado** | **la de la corriente dado el estado** |
|---|---|---|
| **IR** | γ̄₀ᵀΣ₀γ̄₀ **+ E[D²]** | E[V̄] |
| **MR** | γ̄₀ᵀΣ₀γ̄₀ | E[V̄] **+ E[D²]** |

El mismo número, E[D²], está en el bolsillo del estado en IR y en el de la corriente en MR.
La primera fila es Γ̄ᵀΣ^bndΓ̄ y μ^bndᵀV̄ desarmados; la segunda es lo que queda al marginalizar
i_t antes.

Los totales coinciden, y no por casualidad: los dos **son** la varianza verdadera del promedio
grabado dado (μ₀, Σ₀), y mover una varianza de un bolsillo al otro no cambia una suma.
Verificado contra Monte Carlo.

**Esto es lo que hay que decir, y no "tienen la misma varianza".** El segundo término es `gvar_ij`
en IR, y en MR es `gvar_ij` más lo que arrastró al marginalizar. Son objetos distintos. Que sumen
igual es consecuencia, no titular.

## Bloque Cov(estado, corriente)

> **Cov_IR(X_Δ, corriente) − Cov_MR(X_Δ, corriente) = ( Σ_{i₀} μ_{0,i₀} P_{i₀i_t} D_{i₀i_t} )_{i_t}**

o sea **D contraído una sola vez**, donde en la varianza aparecía contraído consigo mismo. No
depende de Σ₀: es el mismo vector con prior multinomial, encogido y con Σ₀ = 0.

Ese vector suma cero sobre i_t, que es por qué el bloque cruzado de IR pasa el canario
`to_Probability_displacement` (`qmodel.h:4650`): es un desplazamiento sobre el símplex.

## El mecanismo

Un solo objeto, D, dos contracciones, y se comportan distinto porque D tiene media cero por fila:

- **D²** se reasigna entero de bolsillo al marginalizar i_t, y el total de Var(corriente) no se
  entera. Por eso la media y la varianza predichas de MR e IR son idénticas desde el mismo estado.
- **D** desaparece, y no tiene dónde ir: un término de ruido no tiene primer momento. Pero su
  perfil por columna es lo que informa sobre X_Δ.

**MR marginaliza justo sobre el índice cuyo perfil llevaba la señal.** Conserva D², pierde D.
Como D es de media cero por fila, la pérdida no se ve en ninguno de los dos primeros momentos de
la predicción y aparece solo en el bloque cruzado, que es el único por el que el dato entra al
estado.

De ahí sale que MR predice el dato tan bien como IR y aprende de él peor. Y sale el caso extremo:
con Σ₀ = 0 el término de estado de MR es cero, su bloque cruzado es cero, y el intervalo se
desperdicia entero, mientras el de IR vale todo D. Es el ejemplo del canal cerrado con corriente
ya presente.

## Cuánto pesa

Fracción de la varianza-que-viene-del-estado que MR manda al ruido, E[D²]/(γ̄₀ᵀΣ₀γ̄₀ + E[D²]),
prior multinomial estacionario, k_on = k_off. Δ̃ = Δ·k_off.

| Δ̃ | 0.01 | 0.05 | 0.1 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|---|---|
| fracción | 0.010 | 0.048 | **0.091** | 0.316 | 0.432 | 0.491 |

Con Σ₀ = 0 la fracción es 1 en toda la fila. En la celda de referencia del manuscrito, Δ̃ = 0.1,
MR saca 9 puntos de su término de estado y los pone en el ruido.

**Predicción falsable, no medida todavía**: la distancia MR/IR debe crecer con Δ̃ y anularse
cuando Δ̃ → 0, porque E[D²] → 0. El plano de la Figura 4 tiene el eje para chequearlo.

## Qué se reescribe desde acá

La frase correcta, para los once sitios:

> MR le cobra al ruido de observación el spread del promedio de intervalo entre estados finales;
> IR lo lleva en el término que viene del estado. Las dos contabilidades suman la misma varianza
> del observable mientras comparten estado, y difieren en el término de intervalo, en la
> covarianza entre el estado y la corriente, y en todo lo que sigue a lo largo de una grabación.

Sitios que hoy dicen otra cosa, o la dicen sin el qualifier:

- `papers/1_method/figures_build_plan.md:195-245` (la fuente; dice "for any prior")
- `papers/_program/nomenclature.md:143-152`
- `papers/1_method/decisions.md:158`
- `papers/1_method/approach.md:165`
- `theory/macroir/notes/gvar_i_overcount_audit.md` (bloque de resolución)
- `theory/macroir/notes/README.md:17`
- `theory/macroir/notes/vr_variance_form_plan.md:20`
- manuscrito: `02_theory.tex:952-966`, el caption de la Tabla 1 en `:1009-1026`,
  `05_discussion.tex:114`, `08_appendix_members.tex:172-189`

**Aparte, y es la misma pasada**: la de-Kalman del 2026-08-06 sacó "Kalman correction" de la
actualización y dejó *the gain* en `02_theory.tex:670, 674, 807, 812, 959`, el caption de la
Tabla 1, `05_discussion.tex:114` y el Apéndice. El manuscrito ya define ese objeto correctamente
("the covariance between the state a channel is in and the current it carries", `:670`); falta
que lo llame así.

## Perímetro

**Verificado**: μ^bnd y Σ^bnd contra Monte Carlo; MacroR sobre el boundary state contra las
fórmulas del código para av=2, los tres bloques, a 1e-16; MR como marginalización previa contra
av=1; las ocho identidades en 36 casos (cuatro pares de tasas, tres ventanas, tres priors
incluido Σ₀ = 0) a 1e-15; y Var(corriente) contra Monte Carlo.
Script: `papers/1_method/decisions/recompute/mr_vs_ir_boxes.py`.

**Camino vivo del código**: `safely_calculate_Algo_State_recursive` (`legacy/qmodel.h:4519`),
llegado por `safely_calculate_Algo_State` (`:5964`). Los objetos de la ventana salen igual de
`calc_Qdt_eig` (`:1769-1774`) y de `calc_Qdtm_eig` (`:1666-1668`), así que el reparto no cambia
según qué flavor de Qdt le toque a cada miembro.

**No verificado, y es una pregunta abierta, no un hallazgo**: hay una familia gemela con otra
álgebra en `qmodel.h:3729-4078` (dentro de `#if 0`, rotulada dead code) y en
`micro_monoid.h:745-1088` (**sin fence**). Leí que en las dos, para av=2, el término de estado
lleva el término de frontera y el de corriente lee el **campo** `gvar_i`, que es la forma total
(`qmodel.h:1774`). Sobre el papel eso es el doble conteo que describe
`gvar_i_overcount_audit.md`. Pero: (a) no corrí nada de eso, (b) no sé si la copia de
`micro_monoid.h` es alcanzable (su helper `safely_calculate_Algo_Pmean_Pcov` no tiene definición
en el repo, o sea que no compilaría si se instanciara), y (c) `qmodel.h:1771-1773` dice
explícitamente que `micro_full` y `micro_monoid` **esperan** la forma total, así que puede ser la
correcta ahí por razones del camino micro que no estudié. Hay que mirarlo con el camino micro en
la mano antes de llamarlo bug.

**Tampoco verificado**: nada de esto corrió el binario. Es una reimplementación independiente en
Python, más Monte Carlo. Coincide con los números que ya estaban en los dumps de la Figura 1
(MR = IR = 1.048475625791748 mientras comparten prior) pero no es el binario. Y es K = 2; la
propiedad de suma cero por fila es general por definición de γ̄₀, así que el argumento lo es,
pero los números no.

## Traducción al marco de filtrado lineal

**NO ESCRIBIR TODAVÍA** (decisión de Luciano, 2026-08-08). Primero esta derivación queda firme en
la base de MacroR. La correspondencia con el filtrado lineal se agrega después, como marco
alternativo que explica lo mismo, y no como el marco en el que se explica.
