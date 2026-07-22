# Figura S1: score medio por paso s̄ₜ(t), 1000 recordings

Análisis de `projects/eLife_2025/figures/archive/supplement/figure_3_bias_pub.Rmd`, salida `Figure_S_score_mean.pdf`.

**Celda:** esquema C⇌O, N_ch = 100, Current_Noise = 1e-4 (la misma que Figura 2), interval_in_tau = 0.1, 1000 recordings independientes. Grilla de 100 pasos de tiempo, 40 de ellos bajo agonista (el "pulso").

**Cantidad graficada:** en cada paso t, el promedio sobre recordings del score por paso s_{t,k} = ∂logL_t/∂θ_k, para cuatro parámetros (k_on, k_off, i = corriente unitaria, N_ch). La banda es el intervalo de confianza (CI) bootstrap del 95 % sobre recordings; la línea nula está en 0. Layout: dos bloques con escala agrupada (NR/NMR | R/MR/IR), con y libre por parámetro, mismo tratamiento que Figura 3 B/C.

## Qué testea

Es el primer test de Bartlett resuelto en el tiempo. Para un likelihood bien especificado el score por paso es una diferencia de martingala de media cero, de modo que E[s_{t,k}] = 0 en cada t. El likelihood macro es una aproximación gaussiana por cierre de momentos de un proceso simulado de forma exacta, así que está mal especificado por construcción. Un score medio distinto de cero mide en qué paso y en qué parámetro la aproximación sesga el gradiente. Un punto cuenta como significativo cuando su CI no incluye al 0.

## Resultado principal

A 300 recordings la figura no discriminaba (CIs anchos). A 1000 los CIs se angostan alrededor de 1.8× y aparece la separación. Conteo de pasos significativos sobre los 40 del pulso:

| algoritmo | k_on | k_off | i | N_ch |
|---|---|---|---|---|
| NR  | 6  | 7  | 7  | 6  |
| NMR | 6  | 5  | 8  | 6  |
| R   | 24 | 16 | 38 | 24 |
| MR  | 26 | 16 | 38 | 27 |
| IR  | 2  | 5  | 2  | 4  |

IR se queda en 2 a 5 de 40, el nivel de falsos positivos esperado a 95 % (0.05 × 40 = 2). R y MR llegan a 38 de 40 en la corriente unitaria: sesgo sistemático en casi todo el pulso.

## Dos modos de fallo distintos

Las dos familias fallan de manera cualitativamente diferente, y eso se ve cruzando la tabla anterior con las magnitudes:

| max\|score medio\| | k_on | k_off | i | N_ch |
|---|---|---|---|---|
| NR  | 9.79 | 2.29 | 10.03 | 10.03 |
| NMR | 1.68 | 0.93 | 2.01  | 2.08  |
| R   | 0.63 | 0.41 | 1.14  | 0.79  |
| MR  | 0.48 | 0.51 | 1.38  | 0.91  |
| IR  | 0.35 | 0.30 | 0.71  | 0.66  |

**NR y NMR (no recursivos): fallo grande y localizado.** El max|score medio| llega a ~10 (k_on, i, N_ch en NR), pero solo 6 a 8 pasos del pulso resultan significativos. El sesgo se concentra en las transiciones (encendido y apagado del pulso), donde la predicción de media suave no puede seguir el cambio brusco. En la meseta del pulso el score vuelve cerca de 0.

**R y MR (recursivos sin medición integrada): fallo chico y persistente.** El max|score medio| es ~0.4 a 1.4, pero significativo a lo largo de casi todo el pulso (38/40 en i) y más allá. El filtro recursivo se sesga siempre que la dinámica dentro del intervalo no es trivial.

## Certificado de calibración de IR

Subir el poder de 300 a 1000 recordings no destapó sesgo en IR: durante el pulso queda en el nivel de azar (2 a 5 de 40), lo que confirma que está calibrado con poder de sobra. En magnitud IR es además el más chico (max ~0.3 a 0.7). En la corriente unitaria i, el parámetro más difícil, IR da 7 puntos significativos en total sobre 100, contra 70 en R y 78 en MR.

Nota: los primeros pasos pre-agonista dan score idénticamente 0 para todos los algoritmos (sin agonista, los canales cerrados no aportan gradiente). El transitorio visible en la figura es el del encendido, no la inicialización del filtro.

## MR peor que R

En conteo total sobre los 100 pasos, MR supera a R en sesgo de i: 78 contra 70 puntos significativos. La diferencia está en la cola posterior al pulso. MR arrastra el sesgo hacia la fase de decaimiento, R se recupera más rápido. Es la firma del condicionamiento al comienzo del intervalo: MR mistimea la actualización y el error se filtra al decaimiento. Coincide con el rol de MR como el arreglo simple que no alcanza, y con que IR (centrado en el intervalo vía medición integrada) es lo que restaura el gradiente insesgado.

## Enlace con las otras figuras

S1 es el compañero a nivel de inferencia de la Figura 3 panel B (media del residuo r̄(t), a nivel de dato). Junto con la autocorrelación del score (Figura 3 panel E) cierra los dos chequeos de Bartlett: media cero (S1) y blancura (ACF). Un score sesgado y al mismo tiempo correlacionado es la definición de F ≠ J, que es lo que alimenta la sobreconfianza de Fisher en la Figura 2 y la distorsión de la Figura 4.

## Estado y decisiones

- Construida y renderizada a 1000 recordings; escala agrupada en dos bloques, y libre por parámetro.
- Datos a Current_Noise = 1e-4, consistentes con Figura 2.
- Abierto: como ahora separa R/MR de IR de forma limpia, evaluar si sube de suplemento a figura principal.
- Caption pendiente: declarar n = 1000 recordings y que el sombreado gris es el pulso de agonista.
- Reproducibilidad: conteos calculados con `seed = 1`, B = 200 réplicas bootstrap. B no afecta la significancia (la fija n = 1000 recordings), solo suaviza la estimación del CI.
