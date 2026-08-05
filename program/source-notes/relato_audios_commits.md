# MacroIR: audios y commits, semana a semana

Este documento pone las dos series una al lado de la otra: lo que Luciano dijo en los audios
y lo que efectivamente entró al repositorio. El relato armado solo con audios está en
[relato_semana_a_semana.md](relato_semana_a_semana.md) y sigue valiendo como lectura del
material hablado; acá el aporte es el pareo, y lo que se ve cuando las dos series no coinciden.

## Cómo se armó

- **Audios**: 288 archivos únicos (deduplicando `.mp3`/`.ogg`/`.opus` del mismo mensaje),
  repartidos en 100 días distintos, de los grupos "MacroIR" 1 a 13. Fecha tomada del nombre
  del archivo de WhatsApp.
- **Commits**: 295 desde el 20 de agosto de 2025 hasta el 2 de agosto de 2026, con fecha de autor.
- **Unidad**: semana ISO, de lunes a domingo. La ventana cubre 49 semanas.
- Tablas de trabajo: `tmp/semanas.md` (pareo completo), `tmp/audio_dates.txt`, `tmp/commits_window.txt`.

De las 49 semanas: 31 tienen audios y commits, 8 tienen solo audios, 5 tienen solo commits,
5 están vacías. Las semanas de solo audio son casi todas semanas de decisión; las de solo
commits son de ejecución o de escritura. La correlación numérica entre cantidad de audios y
cantidad de commits por semana es baja (r = 0.17), y no mejora corriendo los commits una
semana (r = 0.16). O sea que hablar mucho no predice commitear mucho: el vínculo es de
contenido, no de volumen.

## Una corrección al relato anterior

El relato de solo audios afirma que no existe un grupo "MacroIR 11" y que entre el 16 de enero
y el 22 de marzo de 2026 el paper estuvo parado. Las dos cosas hay que revisarlas, y la
evidencia es del propio repositorio.

La carpeta `Chat de WhatsApp con MacroIR 11/` existe, sin trackear en git, con 8 audios:
uno del 11 de febrero, dos del 13 de febrero y cinco del 6 de marzo (estos últimos sin
transcribir todavía). El del 13 de febrero empieza así:

> "Estuve varios días trabajando toda la semana con Macro IR, con el paper, trabajando
> fundamentalmente en la teoría de la deformación de la log-likelihood."

Y en ese mismo audio nombra la matriz de distorsión y la descompone en componente local por
sample y componente de correlación entre likelihoods sucesivas, que es exactamente la
descomposición que el relato anterior daba por nacida en marzo. En el segundo audio del 13 de
febrero ya está identificado el problema de las mediciones sin información cinética (antes de
que se abran los canales), la singularidad de la matriz en esas regiones, y la salida por el
test de Wald.

Del lado de los commits, el 10 de febrero entran dos cosas: `c29fd72` "boostrap apparently
working" y `1ee1a22` "Paper 2: add MacroIR eLife 2025 planning pack", diez documentos y 645
líneas con plan maestro, workboard, log de decisiones, storyboard de figuras y grilla de
experimentos. Y del 12 al 16 de marzo hay once commits más de bootstrap sin ningún audio
asociado.

La lectura correcta del hueco es otra: entre mediados de enero y fines de marzo hubo dos
semanas de trabajo real y bien documentado (10 al 13 de febrero, 12 al 16 de marzo) separadas
por silencios largos. No fueron dos meses parados, fueron dos meses discontinuos.

---

## Etapa 0. No hay paper, hay un programa (26 ago – 28 sep 2025)

**Semanas del 25/8 y del 1/9.** 27 audios, cero commits. Camina por la 9 de Julio pensando qué
publicar, aparece la idea de "una colaboración conmigo mismo", y queda fijado el test que va a
sostener el año: esperanza del score cero, covarianza del score igual a la FIM. El repositorio
no registra nada. Es la primera de las ocho semanas mudas en código, y ya marca el patrón.

**Semana del 8/9.** El 8 de septiembre escribe el primer esqueleto del paper, en cinco
secciones. El sábado 13 entra `0442a72`, "big big cleaning and directory re-organization":
1028 archivos tocados, 2.955.686 líneas borradas. La refundación del software empieza cinco
días después del primer esqueleto del paper, y no tiene nada que ver con él.

**Semanas del 15/9 y del 22/9.** 39 audios y 35 commits, y es el punto de máxima divergencia
entre los dos registros. Los audios hablan de DSL, de tests como tripleta (función,
postcondición, dominio), y el 26 de septiembre de Homotopy Type Theory. Los commits van todos
al mismo lado: CI en GitHub, medición de tiempos de compilación, limpieza de `grammar_typed.h`
sacando punteros crudos en tres fases, soporte JSON dentro del DSL. `include/macrodr` concentra
60 archivos tocados y `docs/perf` 25. Ninguna línea de ciencia. El 20 de septiembre menciona
al pasar que mandó la versión final del manuscrito anterior.

---

## Etapa 1. La decisión, y el código la sigue en tres días (29 sep – 31 oct 2025)

**Semana del 29/9.** El 29 de septiembre se queda sin tokens de Codex, se frena, y dice en voz
alta que trabaja mucho y no avanza. De ahí sale la decisión: los tests de FIM, sampling y
evidencia son el paper. El 2 de octubre define la secuencia concreta, correr muchas
simulaciones y calcular likelihood con score y FIM en cada una.

Los commits de esa misma semana son, en orden: `7a84394` "simulation works" (3/10),
`ebe866f` "likelihood command on simulation works" (4/10), `6d74b7f` "dlikelihood command
works" (5/10). Tres días, tres commits, y son literalmente los tres pasos que enumeró el 2 de
octubre. Es el pareo más limpio de todo el año.

**Semanas del 6/10 y del 13/10.** El 7 de octubre encuentra punteros colgando en el delta X.
El 9 entra `be42bc6` "it calculates derivative" junto con la remoción completa de la
normalización de Nelson del eigensystem. Después viene la pared: `84d1e8a` (15/10) "do not
compile but codex may screw up even further", `a06963a` (16/10) "compiles and the derivatives
are closer, but not completely right", `45507ed` (19/10) "changes in eigensystem interface,
derivative still not safe". Los mensajes de commit de esa semana son el diario de frustración
más explícito del repositorio.

El 17 de octubre llega el diagnóstico hablado: no se pueden derivar autovectores porque no son
una función sino un espacio, así que hay que derivar la matriz de probabilidad de transición y
las conductancias condicionales. El commit del 19 cambia la interfaz del eigensystem. Del
diagnóstico al test verde pasan siete días: `b091574` (22/10) "derivative test almost works,
bug on ATP step with nsamples==0 solved", `a18ffa7` (24/10) "derivative test successful",
`ca9177d` (25/10) "test on Qdtm derivative", `256cb15` (30/10) "test on derivatives passes".

**Semana del 27/10.** El 31 de octubre saca el resultado de la distancia KL entre prior y
posterior, y hace el diagnóstico del código como "un programa monstruoso". Los cinco commits
de ese día son de infraestructura de compilación en GitHub. El resultado teórico no toca el
repositorio.

---

## Etapa 2. El giro conceptual pasa fuera del repositorio (3 nov – 19 dic 2025)

**Semana del 3/11.** Seis audios, cero commits. Es la semana del segundo punto de inflexión, el
6 de noviembre: se da cuenta de que puede prescindir del MCMC y publicar el algoritmo apoyado
solo en la likelihood. Eso saca del paper la evidencia, el sampling y la confusion matrix. La
decisión que más contenido eliminó del paper no dejó una sola línea en git.

**Semana del 10/11.** El 10 arma la narrativa en cuatro puntos y enumera el roster completo de
algoritmos. El 12 confirma con Gemini que un filtro de Bessel después de un Kalman suma
estados, no los multiplica. Los commits llegan el 13: `0ffc5c4` crea el namespace `p2x2` y
`9068ca8` implementa los esquemas CO, CCO y COC. Vale la pena notar la fricción: los tres
esquemas extra son maquinaria de confusion matrix, la figura que había declarado caída siete
días antes. El código siguió caminando una semana hacia una figura ya abandonada.

**Semana del 17/11.** El 14 planifica los datos de la figura 1 y decide reemplazar ATP por
agonist. El 16 entran vectores y tuplas en el DSL (`08c0849`) y el experimento `eLife_2025`
queda definido (`6d2002e`). El 17 y 18, `6398363` y `65dd0b7` hacen el reemplazo ATP → agonist
y habilitan el guardado de substeps de la simulación, que es exactamente el dato de la figura 1.
Once commits en `projects/p2x2` y `include/macrodr`. Esta semana el pareo es uno a uno.

**Semana del 24/11.** Once audios, cero commits. Es la semana de la corrección de Taylor
encontrada con DeepSeek, del cálculo político sobre salami slicing y el techo de cristal en
eLife, de "los tres papers del buen humor", y del mecanismo de seguridad del simplex. Otra vez:
la semana de mayor densidad conceptual del bimestre no tiene código.

**Semana del 1/12.** El 2 de diciembre entra `a3e0a89`, "theoretical results". Es el commit más
engañoso del año. El mensaje habla de teoría y el grueso del diff son 41 documentos nuevos en
`docs/theoretical results` con las derivaciones de MacroIR, MacroTaylor y MacroTaylorIR. Pero
adentro también reescribe `legacy/qmodel.h` de punta a punta, 3056 líneas movidas, y en esa
reescritura el camino no recursivo pierde el término `N·ms` de la varianza de conductancia del
intervalo. El error se descubre y se corrige el 31 de julio de 2026 (`1f7138b`). Doscientos
cuarenta y un días.

**Semanas del 8/12 y del 15/12.** `write_csv` implementado y funcionando (13 y 14/12), bug de
la log-likelihood en diagnósticos resuelto (15/12), y `b5798a1` (18/12) "defines algorithm but
crashes". El 19 de diciembre graba el audio más honesto del período, "estoy totalmente
bloqueado, no sé bien qué poner", y se desatasca hablando: define qué mostrar en la figura 1 y
descubre que sin la corrección de varianza MacroIR sale cualquier cosa, y que para comparar
justo hay que centrar MacroR y MacroNR en el medio del intervalo. El commit del día anterior
dice que el programa define el algoritmo y crashea. Las dos series describen el mismo atasco
desde adentro y desde afuera.

---

## Etapa 3. La crisis y la opción nuclear (20 dic 2025 – 15 ene 2026)

**Semana del 22/12.** El 20 decide dejar micro IR afuera y encuentra el argumento pedagógico
del boundary state. El 23 diseña la figura 2 en el espacio (probabilidad al inicio,
probabilidad al final) y ese mismo día entra `c3493b9` "diagnostic of the meta-state". El 25,
en un audio de 16 minutos, tira abajo esa misma figura 2. El commit que la implementa tiene
dos días de vida útil.

**Semana del 29/12.** El 29 entra `5db2beb` "figure 2 bug in calc_P". El 30, en cuatro minutos,
se le cae el plan: no hay diferencia entre MacroMR y MacroIR en el test del gradiente, y
anuncia la opción nuclear, calcular la Fisher Information Matrix. El 2 de enero entra
`347be45` (manejo de `SymmetricMatrix` y `sqr_X`) y el 4 de enero `c2d1fa0` "Hessian and
covariance of Gradient implemented". Cinco días entre el anuncio de la opción nuclear y la
opción nuclear andando.

**Semanas del 5/1 y del 12/1.** Cinco audios, cero commits. El 11 de enero ordena las dos
estimaciones de la FIM y encuentra la distinción que motoriza todo el análisis posterior: la
covarianza de la suma de los scores contra la suma de las covarianzas individuales. El 14 y 15
nace el factor de expansión de la varianza. El objeto central del paper nace en una semana sin
código.

---

## Etapa 4. Los dos meses discontinuos (16 ene – 22 mar 2026)

**Semana del 9/2.** El 10 de febrero, `c29fd72` "boostrap apparently working" y `1ee1a22`, el
planning pack completo del paper para eLife. El 11 el audio discute cómo llamar al comando que
procesa el batch de simulaciones y qué mostrar en los paneles. El 13 nombra la matriz de
distorsión y la descompone en componente de sample y componente de correlación, y encuentra el
problema de las regiones sin información cinética. Es una semana de trabajo completa, con
teoría, código y planificación editorial, que el relato de solo audios no tenía.

**Semana del 2/3.** Cinco audios el 6 de marzo, sin transcribir en esta carpeta. Es el audio
donde dice "Macro IR me quedó totalmente fuera de mi mente por alguna razón".

**Semanas del 9/3 y del 16/3.** Once commits, cero audios: bootstrap implementado y compilando
(12/3), regularización de `idm_matrix` removida, "boost and analysis runs" y la implementación
de subspace (15/3). Es la imagen espejo de las semanas mudas: acá trabaja y no habla.

---

## Etapa 5. La maquinaria de validación (23 mar – 2 may 2026)

**Semana del 23/3.** Abre el grupo 12 con un roadmap seco, y el repositorio hace exactamente
eso: cinco commits seguidos llamados "project definitions and documents reorganization" 1 a 5,
más la limpieza del DSL (referencias, punteros crudos). 61 archivos en `program/source-notes`,
57 en `theory/macroir`, 42 en `papers/macroir-elife-2025`. Es la semana en que el repositorio
deja de ser un programa con documentos adentro y pasa a ser un paper con código adentro.
El 27, el audio de 23 minutos que empieza con "volví a casa y no pude hacer nada, me tiré en el
sillón".

**Semanas del 30/3 y del 6/4.** Los tipos indexados: `d9671aa` (30/3) implementado pero no
integrado, `ad685e9` (1/4) "half cooked refactoring that does not work", `567ff6b` el mismo día
"major refactoring now indexing is a reality", `9057b83` (8/4) "Indexed works on tuples, sets
and vectors, it took a lot". El 4 de abril habilita Taylor porque el gradiente le da bias
(`eae07dd`, "Macro_ITaylorR results can be saved").

**Semana del 6/4, cierre.** El sábado 11 de abril entra `1839232` "bug that produced Inf at low
Number channels solved", y ese mismo día graba: "ayer fue un día excepcional para MacroIR
porque finalmente pude analizar los datos". Los tres primeros hallazgos reales del paper
(el ruido instrumental no afecta los diagnósticos, el bias depende del número de canales, la
distorsión depende solo del largo del intervalo) llegan el día que se destraba el bug que
impedía correr con pocos canales.

**Semana del 13/4.** La semana más hablada del año: 23 audios, 16 de ellos el viernes 17. El 14
se angustia por la matriz de 600×600×1000 del bootstrap y decide bajarse, y al rato se da
cuenta de que el cálculo es barato. Ese mismo 14 entra `48ae58d` "cross correlation working".
El 16 llega el gráfico que buscaba, la cross-correlación del residuo estandarizado por
algoritmo y por lag, y el 18 y 19 entran el refactor del diagnóstico de derivadas
(`9b61728`) y los diagnósticos de no identificabilidad con guardas de NaN (`bdfd24e`).

**Semana del 20/4.** Un solo commit, `a6678b7` "working on figure 2 of paper". Los audios son
los de la paradoja confirmada (el error de los parámetros corregido por la covarianza del score
da constante) y la decisión de contarlo aunque no lo entienda. También la crisis de sentido del
22, "estaba pensando si MacroIR fue un fracaso".

**Semana del 27/4.** Doce audios, un commit: `30dd0e0` (29/4) "micro_R MR IR working", 35
archivos y 6277 líneas, con `tests/microir/test_micro_derivatives.cpp` y
`tests/math/test_micro_full.cpp` nuevos. El audio del mismo día es el cierre conceptual:
micro IR no tiene correlación temporal ni distorsión, así que lo que se ve en macro es producto
de la aproximación normal. El commit y el argumento entran juntos.

---

## Etapa 6. Cluster, escritura y el bug del simplex (4 may – 28 jun 2026)

**Semanas del 4/5 y del 11/5.** Acá pasa el episodio más importante del año del lado del
código, y es el que peor se ve desde los audios. El viernes 8 de mayo entra `b7ad53a`,
"Restore SymmetricMatrix full-storage", que cierra un bug abierto el 2 de enero: desde
`347be45`, `SymmetricMatrix::set` escribía un solo triángulo, y todo lo que no leía por
`operator()(i,j)` (entre otras cosas `Lapack_Full_Product`, que lee los dos) trataba media
matriz de covarianza como cero. El commit declara la ventana él mismo, "Bug window:
2026-01-02 → 2026-05-08", y también el alcance: "Pre-window manuscript figures clean;
in-window figure-2 derivative diagnostics suspect". El sesgo pegaba en los algoritmos
recursivos, o sea que durante cuatro meses la distorsión de MacroIR salió más alta que la
real, y quedaba tapado por la renormalización de `to_Probability`, así que nunca falló, solo
sesgó. El análisis completo está en [bugs_que_trabaron.md](bugs_que_trabaron.md).

Los audios no lo registran. El último antes del arreglo es del 6 de mayo y habla de bajarle
la complejidad cognitiva al loop de MacroR; el siguiente es del 18 y ya está en el abstract.
El arreglo cae en el medio de doce días de silencio. La única mención en todo el corpus es de
pasada el 25 de mayo, como ejemplo de tiempo perdido, dentro de una digresión sobre Luthier:
"uno pierde tanto tiempo con estas cosas, y date cuenta por ejemplo con eso, con el error con
las matrices simétricas".

El resto de las dos semanas: `a7d6ff9` monoide de micro_ir, `01483ad` (10/5) "CI should be
working again. MicroIR bug solved", y el 14 `a3241c0`, que arregla la fórmula de los momentos
condicionales de conductancia y jubila los `force_*_in_range` reemplazándolos por canarias.
Aparece por primera vez `docs/bibliography` con 25 archivos.

**Semana del 18/5.** Abre el grupo 13 el 19 "porque empiezo a escribir el manuscrito", y el
primer problema es el abstract, que describe MacroIR y no dice nada nuevo respecto de
Communications Biology. El 20 juega con la adjunción entre sampling y likelihood y él mismo la
descarta. Los commits del 21 son de cluster: slurm de eLife, `build_cluster.sh` acotando `-j`,
y el relajamiento de las cotas de la canaria del simplex a warn-only.

**Semana del 25/5.** Cincuenta y dos commits, cuatro audios. Es la semana más pesada del
repositorio en todo el año y casi toda es operación: dispatch de figura 2 por SLURM,
paralelización a nivel de simulación, Clementina lista con gcc15 y OpenBLAS, provenance por
sidecar `.binary`, y al final del domingo 31 el renombre conceptual de IDM a "Likelihood
Information Distortion" en teoría, docs y código, más el suplemento de la IDM posterior en tres
commits. 59 archivos en `projects/eLife_2025`. El audio del 31 cuenta la salida de la crisis
del Hessiano indefinido: medir en el óptimo y duplicar el análisis en versión likelihood y
versión posterior.

**Semana del 1/6.** El 2 y 3 la crisis vuelve por otro lado, dispersión grande de la matriz de
distorsión con 10.000 canales e intervalo 0.01, y sospecha de un bug de paralelismo. El 4
entran el Gauss-Newton genérico y el renombre a `Gaussian_Fisher_Information` (`a105078`), y
también `ad5d594` "audios, theory text etc". El 5, damping del GN y análisis MLE por grupo.

**Semana del 8/6.** El 9 encuentra el salto en la función trust coefficient y la reescribe como
un mínimo suave. El 10 cuenta el bug completo en doce minutos: el mismo alfa se usaba para
`P_min` y para `P_cov`, y en la región donde lo esperado coincidía con lo encontrado la derivada
caía en un régimen o en el otro. Ese mismo 10 de junio entran `361ccba`
"trust-coefficient FD-discontinuity fix + per-sample detailed Fisher diagnostic" (1614 líneas)
y `ac08e4d` con los scripts de caza del bug y el análisis en R. El bug del año se explica y se
arregla el mismo día.

**Semana del 15/6.** Dieciséis commits, casi todos de infraestructura de corrida: version stamp
por hash de git, dependencias entre jobs de SLURM, paralelización en dos niveles del MLE por
grupo, paralelización del bootstrap y de la simulación por recording, y `433ed13` con la
tolerancia adaptativa del Gauss-Newton. Los audios del 15 y el 18 cuentan el otro lado: con
pocos canales la covarianza empírica sale mucho más grande que la predicha, así que hay que
agrupar réplicas, y termina eligiendo grupos de 10, 100 y 1000.

**Semana del 22/6.** Ocho audios, cero commits. Es la semana del mapa de regímenes (multinomial
con pocos canales, poissoniano con intervalos cortos, gaussiano en el medio) y del audio de 32
minutos del 23 donde dicta el paper entero de punta a punta. La columna vertebral argumental
del paper se define, otra vez, en una semana sin un solo commit.

---

## Etapa 7. Tres papers durante tres días (29 jun – 2 ago 2026)

**Semana del 29/6.** El 2 y 3 de julio tiene cinco figuras y cinco suplementarias y el problema
de mostrar la exploración completa del espacio. El 4 y 5 entran `7766cd8` "paper figures 1-4
many 5 candidates", `b028c03` "gaussian fisher used for information and bias distortion",
`1c2ae6f` "gaussian also for corrected covariance" y el dispatch correspondiente. La decisión
hablada (definir la distorsión con la Fisher gaussiana en vez de la numérica, porque es más
estable) y su implementación caen en el mismo fin de semana. 185 archivos en
`projects/eLife_2025`.

**Semana del 6/7.** El 6 dispatch para ruido 0.05 a 0.5. El 11 (`019dbb6`) reorganización de
las figuras del paper, y ese mismo día el audio de 28 minutos que fija el alcance (dos estados,
no estacionario, sin micro IR, sin datos experimentales) y saca el hallazgo de que la Fisher
registra que uno deja de ganar información sobre el número de canales una vez que el número de
canales abiertos deja de subir. También ese día entran los audios del grupo 13 al repositorio
(`7239137`) y la bibliografía (`e9798cb`).

**Semana del 13/7.** Trece commits, cero audios. Es la semana de escritura pura: master list y
planning pack completo del paper, archivado de los drafts superseded dejando `elife_paper.tex`
como cabeza, mapa de prior art, argumento del score como martingala, plan de escritura, y al
final de la semana el plan de cuadrados mínimos no lineales (`67b9ce2`, `82b956f`). 81 archivos
en `papers/macroir-elife-2025`.

**Semana del 20/7.** El lunes 20 pasan dos cosas juntas. Por un lado LSE entra al programa:
`18ead3a` "nonlinearsqr runs", `a416eff` "figure 1 lse". Por otro, `dafbb5f` "decision: 3
papers" y `b6ff4e1` "restructure: three-paper program (method / map / micro)" crean
`papers/_program/` con `program.md`, `paper-2.md`, `paper-3.md` y el índice.

El jueves 23 entran `45102de` "paper both" y `2b3d972`, y el viernes 24 `02e068f` "paper both
working". La estructura de tres papers vivió tres días en el repositorio. El audio del 24 lo
cuenta desde adentro: estuvo a punto de dividir el paper en tres y lo volvió a juntar, ahora
con LSE incluido, porque "con least squares no tenés bias en los parámetros, pero el error que
te da está subestimado hasta cien veces".

**Semana del 27/7.** El 28 explica en el audio por qué se juntó todo: el primer paper era
invendible para eLife, porque MacroR no lo usa nadie y todo el mundo usa cuadrados mínimos. En
el mismo audio aparece el mapa de cinco regiones con las condiciones experimentales reales
dibujadas encima, y dicta la introducción completa.

Los commits del 29 son el paper escribiéndose: consolidación de los ocho planes de sección,
un brief por sección al lado de su `.tex`, el abstract reescrito en 219 palabras
(`5f1c921`), y `b861ad8`, "decisions: NMR is MacroINR, so the drop reason was false; reopen".
Ese commit reabre lo que el 31 se corrige en `1f7138b`, que restaura la varianza de conductancia
del intervalo en el camino no recursivo y renombra NMR a INR. Es el cierre del arco abierto el
2 de diciembre de 2025.

El audio del 1 de agosto narra esa corrección y saca de ella el resultado más limpio del paper:
la corrección de intervalo elimina el bias, la recursión elimina la inflación de la varianza, y
la recursión tiene que considerar los dos extremos. El 2 de agosto, `c675578` "LSE y ILSE" es el
último commit de la ventana, y los cinco audios de ese día dejan anotado el pendiente, la
anisotropía de la matriz de distorsión de cuadrados mínimos.

---

## Lo que se ve al parear, y no se ve en cada serie por separado

**Dónde vivió el trabajo.** El perfil de directorios tocados por semana traza tres fases
netas. De septiembre a enero manda `include/macrodr`, `src/cli`, `src/core` y `docs/perf`:
se está construyendo el programa. De marzo a junio manda `projects/eLife_2025`, hasta 59
archivos en una semana: se están produciendo los datos. De julio en adelante mandan
`papers/1_method`, `papers/_program` y `docs/bibliography`, hasta 297 archivos en la semana del
20 de julio: se está escribiendo. Las tres transiciones vienen precedidas por una semana de
audios sin commits.

**Las semanas mudas.** Las ocho semanas con audios y sin commits son las semanas del 25/8 y 1/9
(qué publicar), 3/11 (prescindir del MCMC), 24/11 (la corrección de Taylor y el cálculo del
salami), 5/1 y 12/1 (nace la matriz de distorsión), 2/3 (el audio del olvido) y 22/6 (el mapa de
regímenes y el dictado completo del paper). Cinco de los seis giros que hicieron el paper
ocurrieron en semanas sin un solo commit. La excepción es el del 28 de julio, que ocurre en una
semana de diez commits porque para entonces decidir y escribir eran la misma actividad.

**La velocidad de bajada.** Cuando la decisión hablada implica código, el código llega rápido y
la distancia es informativa. Tres días de la secuencia enumerada el 2 de octubre a los tres
commits que la implementan. Cinco días del anuncio de la opción nuclear al Hessiano andando.
Siete días del diagnóstico de los autovectores al test de derivadas verde. Cero días para el
bug del trust coefficient, para la cross-correlación y para micro IR. Lo que tarda meses no son
las implementaciones, son las decisiones de qué implementar.

**El código sigue caminando después de la decisión.** El 6 de noviembre declara caída la
confusion matrix y el 13 de noviembre entran tres esquemas cinéticos nuevos. El 23 de diciembre
entra el diagnóstico del meta-estado para la figura 2 y el 25 la figura 2 se cae. La inercia
entre decidir y dejar de hacer es de una semana o menos, pero existe y deja rastro.

**El arco largo.** El 2 de diciembre de 2025 un commit titulado "theoretical results" reescribe
`qmodel.h` entera y pierde un término de varianza en el camino no recursivo. El 31 de julio de
2026 otro commit lo restaura. Doscientos cuarenta y un días, y el error solo se encontró cuando
la escritura del paper obligó a explicar por qué un algoritmo del roster se comportaba como se
comportaba. Ninguna de las dos series lo muestra sola: el audio de agosto cuenta el hallazgo sin
fecha de origen, y el commit de diciembre no dice nada.

**El punto ciego del método.** Este documento se armó leyendo los audios y los títulos de los
commits, y esa combinación tiene un agujero exacto: el trabajo que no se habla y cuyo título de
commit es modesto desaparece de las dos series a la vez. El caso es el bug de las matrices
simétricas del 8 de mayo, que no está en ningún audio y cuyo título parece mantenimiento. Todo
lo que hacía falta estaba escrito en inglés claro en el cuerpo del commit, que en la primera
versión de este relato no se leyó.

El agujero es acotable. De los 295 commits de la ventana, 79 tienen cuerpo y 57 tienen más de
60 palabras, y todos son del 8 de mayo de 2026 en adelante: antes de esa fecha los mensajes son
de una línea ("bug solved?", "advances", "theoretical results"). O sea que hay dos regímenes
documentales y cada uno falla distinto. Hasta abril de 2026 el registro está en los audios y el
git no dice casi nada, con el caso extremo de `a3e0a89`, titulado "theoretical results",
escondiendo una reescritura completa del motor numérico. De mayo en adelante se invierte: los
commits explican y los audios se saltean semanas enteras. El primer commit con cuerpo largo de
todo el proyecto es, justamente, el que este relato había aplanado.

**Lo que el pareo no sostiene.** No hay relación de volumen entre las dos series. Las semanas de
más audios no son las de más commits, ni las anteriores a las de más commits. Quien busque en
esta historia un ritmo de "pensar una semana, ejecutar la siguiente" no lo va a encontrar en los
números; lo que hay es una alternancia irregular entre semanas de decisión sin código y semanas
de ejecución sin habla, con la mayoría de las semanas mezclando las dos cosas.
