# Jobs to be run (cluster)

Cola de corridas pendientes en cluster para las figuras del paper 1. Una entrada por job:
por qué es necesario, el comando exacto, qué produce, qué figura toca y qué editar después.
Cluster access es del autor. Precondición general: el binario de dirac debe hornear un commit
que tenga el algoritmo pedido (`projects/eLife_2025/ops/build_cluster.sh dirac`).

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

**Comando.**

```
NCHS="5 10" N_SIMS="10000 10000" N_NOISE="0.1 0.1" GROUP_SIZE="10 100" \
N_ALGO="micro_IR micro_R" H_RELS="1e-5" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3.sh dirac
```

Parea NCHS/N_SIMS/N_NOISE por índice (dos celdas: N_ch 5 y 10, ruido 0.1, nsim 10000); N_ALGO es eje
aparte (micro_IR y micro_R); GROUP_SIZE se barre dentro de cada job. Total 4 jobs SLURM.
`H_RELS` por defecto ya es `1e-5`; se deja explícito. `RUN_DIR=<folder>` si se quiere escribir en una
carpeta concreta (los archivos numéricos NO llevan sufijo `_G`, así que no chocan con el micro `_G`
existente en `figures/data/87889e6`).

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
2. Más pesada que la `_G`: agrega el stage de Fisher numérico (2·n_params pasadas por registro a θ_sim
   y θ_pool). El lever de costo es GROUP_SIZE (grupo 10 = 1000 refits GN por celda; grupo 100 = 100).
   Si sólo interesa la sample distortion, `GROUP_SIZE="100"` alcanza y es mucho más barato.
