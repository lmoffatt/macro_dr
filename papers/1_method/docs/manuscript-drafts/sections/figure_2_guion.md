# Guion de la Figura 2. Abierto 2026-08-27

Hermano de `figure_3_guion.md`, que existía y este no. Se abre ahora porque el pase de Figura 3 del
2026-08-27 dejó tres cosas que son de ESTA figura y vivían sólo en comentarios `%` o en el guion de
la otra. Mismo formato: lo verificado con su fuente, lo abierto marcado como abierto, y lo que se
sospecha dicho como sospecha.

**Quién la escribe:** `projects/eLife_2025/figures/paper_both/figure_2.Rmd`. Leyenda del manuscrito
en `04_results.tex` (bloque `fig:clouds`). La subsección es `04_results.tex:215-472`.

**Estado al 2026-08-27, después de la auditoría de Luciano (commit `728edf7`):** título nuevo, "The
window and the recursion repair different moments, and the correction is checked against the fits
themselves". **787 palabras en seis párrafos**: 155 / 262 / 84 / 60 / 167 / 41.

---

## Lo que sólo hace esta figura

Mide el RESULTADO en una celda, sobre estimaciones ajustadas. Figura 3 mide el MECANISMO en la misma
celda sin ajustar nada; Figura 4 mide el resultado en todo el plano. Y es la única que responde si
la corrección sándwich funciona, dividiendo la elipse empírica por la corregida. Eso no lo puede
hacer ninguna otra, porque hace falta la nube de ajustes.

---

## LAS DOS ANCLAS DENTRO DEL MISMO PANEL (verificado 2026-08-27)

Lo más importante que este guion tiene para decir hoy, y no está dicho en ninguna parte del
manuscrito.

`figure_2.Rmd`, cabecera líneas 24-25 y el cuerpo en 115-116 y 140-142:

- las **elipses** (Fisher y corregida) salen de `_battery_pool`, o sea **θ_pool**, con
  `cov_scale = 1/REP_GROUP`, `group_size = 10`, `noise = 0.1`;
- el **marcador de sesgo** (el círculo abierto, `truth + DIB`) sale de `_battery_sim`, o sea
  **θ_sim** (línea 65, `Distortion_Induced_Bias`; el vector se arma en 232).

Un panel, dos anclas. No es un error: cada momento se ancla donde corresponde (el sesgo se predice
en la verdad, la distorsión se reporta donde cae el estimador). Pero **no está dicho**, y el párrafo
2 dice hoy "every later figure reads them from the score at one anchor", que se lee como si ésta
usara una sola.

**Por qué importa ahora.** El pase de Figura 3 midió la misma partición en las dos anclas
(`figure_4_source_data_distortion.csv`, celda N_ch 100 / ruido 0,1 / intervalo 0,1, k_off,
`sample x corr = total`):

| miembro | @ theta_sim | @ theta_pool |
|---|---|---|
| R | 0,757 x 1,462 = 1,108 | 1,087 x 1,397 = 1,513 |
| MR | 0,572 x 2,007 = 1,149 | 1,002 x 1,827 = 1,816 |
| VR | 0,767 x 1,620 = 1,240 | 1,041 x 1,498 = 1,559 |
| NR | 1,445 x 19,64 = 14,61 | 1,019 x 19,47 = 17,32 |
| INR | 1,140 x 20,72 = 23,67 | 1,151 x 20,71 = 23,91 |
| IR | 1,115 x 1,011 = 1,130 | 1,113 x 1,010 = 1,127 |

**Los dos centrados leen igual en las dos anclas a la tercera cifra y todos los desplazados no.** O
sea que para IR e INR la mezcla de anclas de esta figura es inocua, y para R, MR, VR y NR los dos
marcadores del panel están en puntos distintos del espacio de parámetros. A DECIDIR: si la leyenda
lo dice en media línea o si se deja como está. Hoy la leyenda dice "Covariances are the
Gaussian-Fisher family at the pooled estimate" y del marcador no dice dónde se evalúa.

---

## ABIERTO, y es de esta figura: dos cómputos del mismo sesgo de primer orden

Traído textual de `figure_3_guion.md`, donde estaba archivado, porque el marcador es de Figura 2.

El "predicted (truth + bias)" es `truth + DIB`, con DIB =
`Probit_statistics_Likelihood_Distortion_Induced_Bias` emitido por el programa y leído de
`_battery_sim` (`figure_2.Rmd:232-259`). En log10, celda de la figura:

| | DIB, i | F⁻¹g de los digests, i | DIB, N_ch | F⁻¹g, N_ch | nube (medido), i |
|---|---|---|---|---|---|
| IR | +0,0014 | +0,0021 | −0,0038 | −0,0065 | sin resolver |
| NR | −0,0886 | −0,0683 | +0,1033 | +0,0820 | −0,103 |
| R | −0,1812 | −0,1307 | +0,2259 | +0,1461 | −0,156 |
| MR | −0,2865 | −0,1780 | +0,3241 | +0,1900 | −0,228 |

Mismos signos, misma estructura (todo en el par de amplitud, nada en las cinéticas, IR en cero),
magnitudes de los digests 25 a 40% más chicas. Dos cómputos del mismo objeto de primer orden sobre
baterías distintas (DIB: `battery_sim`, 10.000 sims; digests: 1.000 grabaciones, F reconstruida de
dm, dv e y_var). Candidatos SIN verificar: punto de evaluación, una F distinta (pooled contra
promedio de las por grabación), o un segundo término en la definición del DIB. **Nada del cuerpo
depende de cuál sea el correcto**: la prosa de Figura 2 cita el DIB, que es lo que el marcador
dibuja, y la de Figura 3 cita el suyo sin decir que sean el mismo número. Se cierra leyendo la
definición del DIB en el código.

Los tres pares que la prosa NO usa y no están en ningún otro lado (se retiraron al `%` el
2026-08-13): predicho contra medido en la corriente unitaria, −0,181 contra −0,156 para R y −0,089
contra −0,103 para NR; el marcador en el número de canales lee +0,324 MR, +0,226 R, +0,103 NR y
−0,004 IR.

---

## La oración del differenced check, que ahora tiene número

Párrafo 6, hoy: los dos brazos de mínimos cuadrados "are also the two arms whose analytic Fisher the
differenced check flags (Supplementary File~1), **at a size consistent with the residual measured
here**". El residuo medido acá es 1,6.

El lane `fisher_only` ya corrió y da el número (`fim_vs_gaussian_table.py` sobre
`figures/data/a202e03`, nsim 10.000, congruencia FIM numérico contra gaussiano analítico, media
geométrica de autovalores, >1 = el numérico excede al analítico). Medianas sobre los siete
intervalos a N_ch = 100, ruido 0,1:

| IR | INR | NR | VR | R | MR | LSE / ILSE |
|---|---|---|---|---|---|---|
| 1,000 | 0,999 | 0,991 | 0,935 | 0,874 | 0,872 | **0,784** |

Los dos brazos leen 0,784 en las cuatro décadas de canales, o sea que **el Fisher analítico excede
al numérico por 1/0,784 = 1,28**. Contra el 1,6 del cociente de elipses. Misma dirección, no el
mismo número. **La palabra "consistent" está cargando una comparación 1,28 contra 1,6 y no lo dice.**
A decidir: poner los dos números, o bajar "consistent" a "in the same direction". Lo que NO se puede
dejar es que el lector suponga que 1,6 está explicado.

Nota al pasar, del mismo lane: el analítico de NR está bien al 1% en esta celda (0,991) y se despega
recién en canales altos (1,162 a 1e3, 1,361 a 1e4) y con el signo contrario. O sea que ninguna de
las inversiones que reporta el cuerpo es artefacto del Fisher analítico.

---

## Trampas de dato

1. **`group_size`.** `_mle_cloud_runs.csv` mezcla `group_size` 10 (1.000 registros) y 100 (100). La
   figura del cuerpo filtra a 10 (`figure_2.Rmd:148`, `REP_GROUP`) y el suplemento 1 es la celda de
   `group_size` 100 a N_ch 10 y S̃ 0,005, la única nube del set con mil ajustes de ese tamaño. No
   promediar las dos.
2. **Dos archivos, dos anclas.** `_battery_pool` para las elipses, `_battery_sim` para el marcador.
   Ver arriba. `figure_2.Rmd` no usa `Dconf` en ningún lado (verificado, 0 ocurrencias), así que la
   trampa de "Dconf es el borde del IC, no el estimador" no la toca.
3. **La celda es una sola.** N_ch = 100, S̃ = 0,01, Δ̃ = 0,1. Todo lo que se afirme "en general" tiene
   que salir de Figura 4, no de acá.

---

## Presupuesto, si hace falta comprimir

787 palabras, y el desbalance es visible: el párrafo 3 tiene **262**, un tercio del bloque, y carga
cinco cosas (quién está centrado, la estructura de identificabilidad con su cita, el marcador
predicho, el orden en las cinéticas, y qué repara el promediado). Es el candidato natural si la
subsección tiene que achicarse; los párrafos 4, 5 y 7 ya están en el hueso (84, 60, 41).

---

## Lo que este guion NO tiene todavía

No hay una lectura panel por panel como la que `figure_3_guion.md` tiene fila por fila. Cuando se
haga, sale de la fuente y no del dibujo: el chunk de números de `figure_2.Rmd` y
`figure_2_source_data_*.csv`.
