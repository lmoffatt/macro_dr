# Jobs to be run (cluster)

Cola de corridas pendientes en cluster para las figuras del paper 1. Una entrada por job:
por qué es necesario, el comando exacto, qué produce, qué figura toca y qué editar después.
Cluster access es del autor. Precondición general: el binario de dirac debe hornear un commit
que tenga el algoritmo pedido (`projects/eLife_2025/ops/build_cluster.sh dirac`).

**Orden al 2026-08-14:** primero el Job 2 (barrido de roster, Fisher numérico contra Gaussiano),
que es el que cierra la deuda de `04_results.tex:35-62`. El Job 1 (micro) corre aparte y por la
misma lane. Los dos usan `dispatch_figure_3_fisher_only.sh`, que NO ajusta la nube por grupo: eso
era lo que hacía que estas corridas tardaran días y sólo alimentaba la capstone empírica, que ya
está medida en la lane gaussiana.

---

## Job 1 — micro numeric Fisher (sample distortion anclada en Fisher), N_ch 5 y 10

**Estado:** PENDIENTE (abierto 2026-07-22).

**Por qué es necesario.**
En `Figure_6_micro_macro_linear_sample.pdf` (Gaussian sample distortion), micro_IR se dispara en
la dirección de amplitud a pocos canales: en corriente (i) sube de 1.03 a ~2.6 entre Δ·k_off 0.01 y
0.2 a N_ch 5, y lo mismo en noise y baseline; macro_IR se queda en ~0.9-1.05. Los cinéticos
(k_on, k_off) de micro_IR están bien. Para saber si ese pico es no-Gaussianidad por-muestra **real**
del micro (el conteo discreto de canales abiertos hace fluctuar el score de amplitud más de lo que su
Fisher Gaussiano predice) o un **artefacto de ancla** (comparar J_s de micro contra el Fisher
Gaussiano G_b), hace falta la sample distortion anclada en el Fisher numérico
(`Probit_statistics_Likelihood_Sample_Distortion`).

Hoy no la tenemos para micro: micro se corrió sólo por el camino `_G` (ancla Gaussiana), donde toda
la familia de Fisher numérico queda vacía (sólo `bootstrap_count`, sin valores de matriz). El camino
numérico sí emite **ambas** familias en una sola corrida: `figure_3_mle.macroir` computa el Fisher
numérico y usa `likelihood_derivative_basic_diagnostics_paired`, cuyo preset `basic` emite en el mismo
Vector_Space `Likelihood_Sample_Distortion` (numérico) y `Gaussian_Sample_Distortion` (Gaussiano)
(verificado en `src/core/likelihood.cpp:3312-3313`). El dispatcher numérico ya mapea micro
(`dispatch_figure_3.sh:156-158`).

**Beneficio colateral.** micro_IR N_ch 10 nsim 10000 hoy tiene sólo el cloud (sin `battery_sim_G`),
por eso micro_IR se cae de la columna N_ch 10 en la figura. Esta corrida genera la batería completa
y lo recupera.

**Comando** (ACTUALIZADO 2026-08-14: pasa a la lane `fisher_only`, ver el Job 2).

```
N_ALGO="micro_IR micro_R" NCHS="5 10" N_SIMS="10000 10000" N_NOISE="0.1 0.1" H_RELS="1e-5" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
```

Parea NCHS/N_SIMS/N_NOISE por índice (dos celdas: N_ch 5 y 10, ruido 0.1, nsim 10000); N_ALGO es eje
aparte (micro_IR y micro_R). Total 4 jobs SLURM. `RUN_DIR=<folder>` si se quiere escribir en una
carpeta concreta (los archivos de esta lane llevan el stem `figure_3_fim_*`, así que no chocan ni con
el micro `_G` de `figures/data/87889e6` ni con los `figure_3_*` de 433ed13).

Ya no hay `GROUP_SIZE`: la sample distortion que este job busca sale de la batería apareada (Stage 5),
no de la nube, así que el eje de grupo se fue junto con el stage que costaba los días. Lo que se
pierde es la nube y la capstone empírica de micro, que este job nunca necesitó.

**Qué produce.** `_battery_sim` / `_battery_pool` (sin `_G`) con `Likelihood_Sample_Distortion` **y**
`Gaussian_Sample_Distortion` poblados, más `_mle_cloud` / `_pool` / `_empirical`.

**Figura y edición posterior.** `projects/eLife_2025/figures/in_progress/figure_6_micro_macro_linear.Rmd`.
Cuando esté el dato: agregar `Probit_statistics_Likelihood_Sample_Distortion` a `MATS` (o rehacer el
panel de sample para graficar Gaussian vs Fisher lado a lado), y apuntar la lectura de micro (y de la
celda macro N_ch 5) a la carpeta de esta corrida.

**Avisos.**
1. El Fisher por diferencias finitas es el potencialmente indefinido justo en la dirección de amplitud,
   que es donde micro_IR se dispara, así que `Likelihood_Sample_Distortion` ahí puede ser ruidoso.
   Igual, con las dos anclas se puede juzgar cuál es la razonable en vez de suponer.
2. Más pesada que la `_G` en un solo lugar: el stage de Fisher numérico son 2·n_params pasadas por
   registro a θ_sim y a θ_pool, o sea 24 con seis parámetros libres. A cambio, sacar la nube saca
   los miles de ajustes GN por celda, que era el costo real.

---

## Job 2 — Fisher numérico contra Gaussiano, roster completo, una sola lane

**Estado:** PENDIENTE (abierto 2026-08-14).

**Por qué es necesario.**
`04_results.tex:35-62`, ítem (f), lo declara deuda y aclara que es resultado, no limpieza. La
comparación entre el Fisher analítico Gaussiano y el de diferencias finitas existe hoy sólo en
`figures/data/433ed13` (macro NR/R/MR/IR más el NMR superseded), que el manuscrito declara superseded
y del que `04_results.tex:7` dice que no se cita nada; y para INR, VR, LSE e ILSE no existe en ningún
lado, porque el barrido `_G` es Gaussiano puro y nunca forma F_b (declara el componente
`Likelihood_Numerical_Fisher_Information` y emite cero filas de valor). Corriendo el roster entero en
un directorio de hash de HEAD la medición pasa a ser citable de una pieza, y de paso el NMR superseded
(el que no lleva el término N·ms) queda reemplazado por INR.

Segunda cosa que mide, y es una afirmación aparte: **el Fisher numérico en θ_sim puede no ser
positivo**. Eso ya está medido para macro (corrida de dirac del 2026-06-10, nsim 1024): el sesgo
genuino era sólo NR a muestreo grueso (λ_min ∝ −N_ch), R en interval_in_tau=1 y MR en bolsones, con
IR limpio en todas las celdas. Lo que falta son los miembros nuevos.

**Comando.** Los defaults del dispatcher SON esta corrida, así que:

```
projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
```

Roster `macro_NR macro_R macro_INR macro_MR macro_VR macro_IR nonlinearsqr nonlinearsqr_g`, N_ch
10/100/1000/10000, ruido 0.1, nsim 10000, h_rel 1e-5, semilla 20260814. Son 32 jobs SLURM. El eje de
intervalo (los siete valores de figure_2) se barre dentro de cada job, y es el eje donde la
indefinitud se prende, en el extremo grueso.

Antes, un smoke de un solo job para que hable el parser y no el cluster:

```
N_ALGO="macro_VR" NCHS="10" N_SIMS="8" N_NOISE="0.1" TIME=00:30:00 projects/eLife_2025/ops/slurm/dispatch_figure_3_fisher_only.sh dirac
```

VR es el que conviene para el smoke: es el único miembro que estrena la inyección de
`variance_form`, así que si algo del contrato de archivo está mal, falla ahí.

**Qué produce.** Por celda: `_pool`, `_battery_sim` y `_battery_pool` con stem `figure_3_fim_*`. Las
baterías son las apareadas, o sea que traen la familia anclada en F **y** la anclada en G en el mismo
archivo; ese apareo es la medición. No produce `_mle_cloud` ni `_empirical`.

**Figura y edición posterior.** `figures/archive/figure_4_gaussian_vs_numeric_fisher.Rmd`, que ya lee
`Likelihood_Gaussian_Fisher_Distortion`: agregar el hash nuevo a `FIG4_DATA_DIRS` (línea 103) y los
miembros nuevos a `FIG4_ALGOS` (línea 100). Después, la deuda de `04_results.tex:35-62` se cierra
eligiendo una de las tres formas que el archivo lista (oración con números, tabla suplementaria o
figura suplementaria).

**Avisos.**
1. **De dónde se lee la no positividad.** Del espectro
   (`Eigenvalue_Spectrum<Likelihood_Gaussian_Fisher_Distortion>`), que trae los autovalores negativos
   incluidos (`legacy/lapack_headers.h:3047-3057`, vía el fallback lenient de `compute_psd_decomp`).
   NUNCA del `Min_Eigenvalue<...>`: ese escalar sale de `compute_distortion_scalars`
   (`src/core/likelihood.cpp:591-595`), que se queda sólo con los autovalores por encima de la
   tolerancia de retención, así que nunca devuelve un negativo y se lee como un piso cerca de cero
   justo donde la matriz es indefinida. El comentario de `include/macrodr/cmd/likelihood.h:1012` llama
   a ese escalar "the key FIM_sim indefinite? readout" y no lo es. La inercia sí es legítima de leer
   sobre la distorsión: G_b^{-1/2}·F_b·G_b^{-1/2} es una congruencia por una simétrica definida
   positiva, así que por Sylvester tiene exactamente tantos negativos como F_b.
2. **Por réplica exagera.** λ_min es cóncava, así que E[λ_min(F_i)] ≤ λ_min(E[F_i]). El pipeline ya
   hace lo correcto (`compute_F_b` promedia sobre grabaciones antes de la congruencia), pero no hay
   que reportar nunca "el X% de las réplicas dio indefinida".
3. **Dos configuraciones de parámetros, a propósito.** Los miembros macro y micro corren con los seis
   parámetros libres; los de mínimos cuadrados con cuatro (`unitary_current` y `Current_Noise` fijos,
   direcciones planas ahí). Es la configuración en la que el paper cita cada uno, así que la pregunta
   por miembro es sana, pero una tabla cruzada de los escalares compara un espectro 6×6 con uno 4×4 y
   eso hay que decirlo en la caption.
4. **La semilla va explícita.** `SIM_SEED` se inyecta; `seed = 0` es el centinela que saca el valor de
   `random_device` (`legacy/mcmc.h:37-45`) y el valor resuelto no se escribe en ningún lado, que es el
   agujero de 433ed13. Una sola semilla para todo el barrido también aparea el ensemble entre
   algoritmos dentro de una celda, que es lo que hace que la comparación entre miembros sea comparación.
5. **El paso h.** Si el número de NR (la analítica sobreestimando 1.6) va a ir al texto, conviene una
   celda con `H_RELS="1e-4 1e-5 1e-6"` para mostrar que no es artefacto del paso. Multiplica sólo el
   stage 3.
