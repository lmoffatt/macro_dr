# Los bugs que más trabaron el trabajo

Ordenados por cuánto costaron, no por orden cronológico. Cada entrada cruza el audio donde
Luciano cuenta el síntoma con el commit que lo cierra. Cuando las dos fuentes no coinciden,
está dicho.

Nota de vocabulario: FIM es la matriz de información de Fisher (Fisher Information Matrix);
IDM es la matriz de distorsión de la información, después renombrada Likelihood Information
Distortion; N_ch es el número de canales; Σ es la matriz de covarianza de las probabilidades
de estado, `P_Cov` en el código.

---

## 1. La matriz simétrica que guardaba medio triángulo (126 días, enero a mayo de 2026)

**Fechas**: se introduce el 2 de enero de 2026 (`347be45`, "better SymmetricMatrix handling and
sqr_X(SymmMatrix)") y se cierra el 8 de mayo de 2026 (`b7ad53a`, "Restore SymmetricMatrix
full-storage"). El propio commit del arreglo declara la ventana: "Bug window: 2026-01-02 →
2026-05-08".

**El bug**: el 2 de enero se adoptó una convención de almacenamiento canónico para matrices
simétricas, guardar y exponer solo el triángulo superior. `SymmetricMatrix::set` pasó de
escribir los dos triángulos a escribir uno solo, y `operator()(i,j)` pasó a plegar el índice
para que la lectura por par de índices siguiera funcionando. Internamente es consistente. El
problema es todo lo que no lee por `(i,j)`: el acceso crudo por `operator[]`, `zip`, `reduce`,
y sobre todo `Lapack_Full_Product`, que lee los dos triángulos del almacenamiento real. El
comentario que quedó en el código lo dice con todas las letras:

> "The base Matrix versions call Lapack_Full_Product, which reads BOTH triangles of its
> arguments. SymmetricMatrix only stores ONE triangle (set() writes one side; the other holds
> the construction value, i.e. zero). Without these overloads, TranspMult(M, Σ) returns a
> product where half of Σ is treated as zero, silently wrong, no compile error."

Es decir: durante cuatro meses, cualquier producto de la covarianza de estados que pasara por
esa ruta trató media matriz como cero.

**Qué impedía**: nada, y ahí está el daño. El programa corría, no crasheaba, no tiraba NaN, y
además el síntoma estaba tapado: el commit del arreglo dice que la desviación quedaba escondida
por la renormalización de `to_Probability`. Lo que salía era un sesgo, no una falla.

El sesgo golpeaba a los algoritmos recursivos, que son los que usan Σ en la ganancia. O sea que
pegaba justo donde el paper tenía su tesis: **la distorsión de MacroIR salía más alta que la
real**. Con media Σ leída como cero la ganancia queda mal, la actualización corrige de menos,
los residuos quedan correlacionados, y el desacuerdo entre el score y la Fisher (que es
exactamente lo que mide la matriz de distorsión) sale inflado. El algoritmo propio quedaba
peor de lo que es, medido con su propio instrumento.

Y todo eso corrió durante la etapa en que se estaba pensando el paper, no solo corriéndolo.
Adentro de la ventana caen los hallazgos de abril: que la distorsión no depende del número de
canales sino solo del largo del intervalo (11 de abril), la cross-correlación por algoritmo
(16 de abril), la paradoja de que la resolución de los parámetros no mejora con el sample rate
(21 de abril). El commit del arreglo es explícito sobre el alcance: "Pre-window manuscript
figures clean; in-window figure-2 derivative diagnostics suspect".

Conviene decir qué NO quedó contaminado. Todas las carpetas de datos que alimentan las figuras
del paper (`figures/data/`) llevan hash de junio y julio de 2026: `8fc274d` y `433ed13` del 16
y 17 de junio, `1c2ae6f` del 5 de julio, `82b956f`, `87889e6`, `0ffbda7`, `1f7138b` de julio.
Todas posteriores al arreglo. Lo que se perdió no son los números finales, son cuatro meses de
razonamiento sobre resultados que en parte eran artefacto.

**Qué lo destrabó**: volver a la conducta anterior a enero. `set` vuelve a espejar los dos
triángulos, con el comentario de que eso deja la matriz segura bajo cualquier ruta de acceso y
que las rutas conscientes de la simetría quedan como optimización de rendimiento, no como
requisito de corrección. Además se agregan sobrecargas de `TranspMult` y `MultTransp` que
enrutan por `Lapack_Sym_Product` (DSYMM), que solo mira el triángulo poblado; se alinea
`SymPosDefMatrix::operator-` con `operator+`, porque la iteración vieja sobre el cuadrado
completo escribía dos veces las off-diagonales una vez que el espejado entra en juego; y se
agrega la canaria `to_Probability_displacement`, que chequea que las filas sumen cero. Es
decir: se arregla la causa y se instala el detector que lo habría cazado en enero.

**Se llevó puesto un hallazgo del paper.** Esto se estableció el 2026-08-07, con los audios del
10 de mayo rescatados del grupo 12, dos días después del arreglo. El primero, de 33 segundos, es
una retractación:

> "Una de las paradojas que yo veía en los datos era que más información era menos. Y bueno, eso
> resulta que no es verdad. O sea, si yo tengo más resolución temporal, tengo igual o más
> resolución en el valor de los parámetros. Ese efecto paradójico era simplemente porque había un
> error de código, estaba mal de programa."

La paradoja es la del 21 de abril de 2026, medida en plena ventana del bug, sobre la que había
decidido escribir que no la entendía: "honestidad intelectual ante todo", "tenemos que hacer un
paper más o menos humano, mostrar que uno no entiende eso me parece que está bueno". Era un
artefacto. La versión correcta y más fina llega el 11 de julio, ya post-arreglo: aumentar la
resolución temporal no mejora la resolución de las constantes cinéticas, pero sí la del número
de canales y la conductancia.

Verificado el 2026-08-07: la afirmación paradójica **no aparece en ninguna sección del
manuscrito**. Sobrevivía solo en las notas de fuente, que quedaron corregidas.

**Y el diagnóstico lo cazó.** En el audio del mismo día a las 19:36, la observación que le da la
vuelta al asunto:

> "Una idea un tanto osada sería volver el error a una virtud, es decir, mostrar los datos, los
> análisis que yo había hecho con el algoritmo erróneo, mostrando que justamente no daban bien,
> en el sentido de que no se verificaba la matriz de distorsión de la información, daba un valor
> grande o cercano a dos. Entonces eso ponerlo como un indicador de que es una forma de encontrar
> errores."

O sea que la matriz de distorsión marcó el bug antes de que él supiera que había un bug, dando
alrededor de dos donde tenía que dar uno. Es el mejor argumento posible para la tesis del paper
(el diagnóstico es un test de validez) y viene de una experiencia propia y documentada. Que la
propuesta de contarlo esté sin usar es raro.

**El costo que no se ve en el código**: este es el bug que le dio origen a Luthier. En el audio
del 25 de mayo, hablando de reescribir MacroIR como un sistema que se testee a sí mismo:

> "Uno pierde tanto tiempo con estas cosas, y date cuenta por ejemplo con eso, con el error con
> las matrices simétricas. Yo creo que sí, que una conclusión es que tengo que lutierizarlo."

Es la lección de "el significado vive en los tipos" aprendida por el lado caro. El tipo
`SymmetricMatrix` cargaba una invariante (medio triángulo es el dato, el otro es basura) que no
todos sus consumidores respetaban, y el compilador no tenía cómo saberlo.

---

## 2. El coeficiente de confianza del simplex (dos episodios, abril y junio de 2026)

Un solo mecanismo con dos caras. Nació el 27 y 28 de noviembre de 2025, en un audio donde
Luciano lo describe como un factor que infla la varianza de medición cuando la probabilidad se
quiere salir de [0,1], "ad hoc pero con su elegancia dentro de lo ad hoc". Costó unas tres
semanas repartidas en dos episodios y terminó reescrito desde primeros principios.

### Episodio A: Inf con pocos canales

**Fecha**: se cierra el 11 de abril de 2026, commit `1839232`, "bug that produced Inf at low
Number channels solved".

**El bug**: el coeficiente (`alfa`) multiplicaba a la varianza dentro del logaritmo de la
log-likelihood, `log(2π·y_var·α)`, y lo mismo en la log-likelihood esperada y en la
normalización de ruido de Poisson. Con pocos canales el coeficiente entra en acción y se
achica, y el logaritmo se va a infinito.

**Qué impedía**: analizar el rincón de pocos canales, que es exactamente donde el paper tenía
algo que decir. Es el rincón multinomial del mapa de regímenes.

**Qué lo destrabó**: sacar el coeficiente del logaritmo. El diff es una remoción: `calculate_logL`
y `calculate_elogL` pierden el parámetro `alfa` en las cinco llamadas. Nunca tuvo que haber
estado ahí. Al día siguiente llega el primer hallazgo real del paper, que estaba esperando esto.

### Episodio B: el salto en la Fisher numérica

**Fechas**: síntoma el 31 de mayo, diagnóstico el 9 de junio, causa y arreglo el 10 de junio de
2026, commit `361ccba`. Diez días.

**El bug**: dos defectos acoplados. El mínimo estaba implementado con un `if` sobre el signo del
desplazamiento D, así que trataba distinto la dirección positiva de la negativa; en los puntos
donde lo esperado coincidía con lo medido, D vale cero y una perturbación mínima hacía que la
diferencia finita cayera de un lado en un régimen y del otro en el otro, produciendo un salto en
la derivada segunda. Y se usaba el mismo alfa para `P_mean` y para `P_cov`, por lo que el salto
aparecía en la covarianza y no en la media:

> "El salto se transportaba inmediatamente a la matriz de covarianza pero no a la matriz de
> P_min. [...] Si estabas en una región donde justo coincidía lo esperado con lo encontrado, el
> D va a ser cero, entonces una pequeña variabilidad te mueve para arriba o para abajo y
> calculás la derivada en un régimen o en el otro. [...] El error fue que se usaba el mismo
> alfa para P_min y para P_cov."

**Qué impedía**: confiar en cualquier número de la matriz de distorsión, con 10.000 canales e
intervalo 0.01, condiciones donde no debería pasar nada. El 2 de junio la hipótesis era mucho
peor que la realidad: "tengo miedo que sea algún error de paralelismo que escriban una thread
sobre la otra". Si hubiera sido eso, se caía toda la infraestructura de corrida.

**Qué lo destrabó**: dejar de mirar bootstraps y mirar réplicas individuales para cazar
outliers; construir comandos que detectaran dónde hay un salto en la derivada segunda, que es
lo que permitió localizar la función culpable; y después ir al pizarrón sin abrir ninguna IA y
reescribir la fórmula como un mínimo suave diferenciable, el log-sum-exp, "una hermosa
ecuación". El remate es que la corrección de `P_cov` no hacía falta: la sacó y el error
desapareció.

Notar que este episodio y el número 1 comparten víctima. Los dos corrompen Σ, los dos pegan en
los recursivos, y los dos son de la misma temporada. Entre enero y junio de 2026 la covarianza
de estados estuvo mal por dos razones independientes.

---

## 3. Las derivadas del eigensystem (octubre de 2025, unas tres semanas)

**Fechas**: aparece el 5 de octubre, se declara trabado el 17, se arregla entre el 19 y el 30
de octubre de 2025.

**El bug**: no es un error de programación sino un imposible matemático que el código intentaba
igual. Los autovectores de una matriz no son una función sino un espacio, así que su derivada
no está definida.

> "Los autovectores son indeterminados, no es una función autovectores, es un espacio, entonces
> sacar la derivada de eso es medio un quilombo, es medio como imposible."

**Qué impedía**: todo. En octubre de 2025 el paper entero era el test de que la esperanza del
score da cero y su covarianza da la FIM, y ese test necesita derivadas confiables. El 17 de
octubre dice que lleva "un par de semanas" trabado.

**Qué lo destrabó**: dejar de derivar lo que no se puede derivar, y derivar en cambio la matriz
de probabilidad de transición y las conductancias condicionales al estado inicial y final. Con
ChatGPT saca la fórmula de la exponencial de matriz por bloques con Padé. Descarta
backpropagation: "no es la idea de que sea más rápido, la idea es tener algo robusto y
confiable". El rastro en los commits es el diario de la pelea:

- `45507ed` (19/10) "changes in eigensystem interface, derivative still not safe"
- `b091574` (22/10) "derivative test almost works, bug on ATP step with nsamples==0 solved"
- `a18ffa7` (24/10) "derivative test successful"
- `256cb15` (30/10) "test on derivatives passes"

Hay un segundo destrabe, más chico y más astuto: saca de la comparación la derivada de la
conductancia condicional (que se divide por una P_ij chiquísima) y descarta los puntos donde la
diferencia finita cambia de signo entre la dirección positiva y la negativa. Es la misma
patología que en junio de 2026 le va a costar diez días, pero acá la esquiva en vez de
resolverla.

---

## 4. La fórmula de los momentos condicionales de conductancia (mayo de 2026)

**Fecha**: el cambio de código entra el 14 de mayo de 2026, commit `a3241c0`, seis días después
del arreglo de las matrices simétricas.

**El bug**: una ecuación mal. En el commit está dicho con precisión: "Full LTV
gvar_i = gsqr_i − gmean_i² at all sites (some paths were a partial sum before, missing the
cross-state mean-variance term)". Algunos caminos calculaban una suma parcial de la varianza de
conductancia y se comían el término de varianza entre estados finales. La carpeta de teoría que
acompaña el arreglo se llama, literalmente, `Gmean_ij_gvarij`.

En el raconto del 1 de agosto lo cuenta al revés, diciendo que usaba G_bar_i cuando tenía que
usar la suma en j de G_bar_ij. En ese mismo audio aclara que no se acuerda bien cómo fue, así
que para la dirección conviene el mensaje del commit. Lo que las dos fuentes coinciden en decir
es que había un error en esta familia de fórmulas y que era grave.

**Qué impedía**: nada visible. Contaminaba en silencio, como el número 1 y en la misma ventana.

> "Ahí vino toda una época donde, claro, andaba mal porque tenía mal las ecuaciones. Y yo no sé
> cómo sobrevivía eso."

**Qué lo destrabó**: implementar micro IR. El 29 de abril (`30dd0e0`) entra el algoritmo
microscópico exacto, que por construcción no tiene aproximación normal y por lo tanto tiene que
dar distorsión cero. Cuando micro IR dio limpio y macro no, el residuo dejó de poder atribuirse
a la aproximación. Micro IR entró al trabajo como argumento científico y terminó funcionando
como instrumento de debugging.

En el mismo commit se retiran los `force_gmean_in_range`, `force_gtotal_in_range` y
`force_gtotal_var_in_range`, reemplazados por regularización bayesiana con canarias.

---

## 4bis. El acumulador de la Fisher gaussiana era código muerto (marzo a junio de 2026)

**Fechas**: la Fisher gaussiana nace el 12 de marzo de 2026 (`6bfde5b`). El 4 de junio de 2026
(`a105078`) se descubre que no acumulaba nada. Casi tres meses.

**El bug**: en el camino de derivación automática, la guarda que decidía si acumular preguntaba
por un slot de Hessiano que no existía, así que la condición nunca se cumplía. El arreglo
cambia la guarda por una que pregunta por el componente correcto, y recién ahí empieza a sumar
`XXt(∂μ)/σ² + XXt(∂σ²)/(2σ⁴)`.

**Qué impedía**: nada visible, otra vez. La diferencia con los demás bugs silenciosos es que
acá el número no salía mal, no salía. Cualquier valor de Fisher gaussiana anterior al 4 de junio
de 2026 hay que tratarlo como no calculado, no como mal calculado.

**Qué lo destrabó**: el trabajo de armar el Gauss-Newton genérico, en el mismo commit. Al
necesitar la Fisher gaussiana como curvatura para optimizar, dejó de ser un diagnóstico que se
mira de reojo y pasó a ser una pieza de la que dependía la convergencia, y ahí se notó.

Vale la pena leerlo junto con el ancla de las figuras: la decisión de anclar la distorsión en la
Fisher gaussiana en vez de la numérica es del 4 y 5 de julio de 2026, exactamente un mes después
de que la gaussiana empezara a calcularse. No era una preferencia entre dos estimadores
disponibles, era una preferencia que recién se pudo evaluar cuando uno de los dos existió.

---

## 5. El término N·ms perdido en el camino no recursivo (241 días)

**Fechas**: se introduce el 2 de diciembre de 2025 (`a3e0a89`), se descubre y arregla el 31 de
julio de 2026 (`1f7138b`).

**El bug**: `a3e0a89`, un commit titulado "theoretical results" cuyo diff visible son 41
documentos de teoría, además reescribe `legacy/qmodel.h` entera (3056 líneas) y parte
`safely_calculate_y_mean_yvar_Pmean_PCov` en una copia recursiva y una no recursiva. La copia
no recursiva reimplementó la varianza sin el término `N·ms`, con ms = P_mean·gvar_i, la
varianza de la conductancia media del intervalo dada el estado inicial. La copia nueva aceptó e
ignoró el parámetro durante ocho meses.

**Qué impedía**: nada, y por eso duró tanto. El algoritmo corría y entraba en las figuras.
Simplemente no era el algoritmo que decía ser: todos los scripts pasaban
`variance_approximation = 1` creyendo que el término estaba.

**Qué lo destrabó**: escribir el paper. Para explicar por qué ese miembro del roster se
comportaba como se comportaba hubo que abrirlo. La verificación se hizo contra el fuente de la
sumisión a Communications Biology (`macro_dr_submission b4a0e28`), que sí lleva el término. El
commit del 29 de julio que reabre el tema lo dice sin vueltas: "decisions: NMR is MacroINR, so
the drop reason was false; reopen". El radio de daño es un solo algoritmo, y trae de arrastre un
renombre, MNR pasa a INR. Las figuras de INR quedan viejas hasta que se re-corra.

Vale como advertencia sobre commits grandes con mensaje inocente: un título de teoría
escondiendo una reescritura completa del motor numérico.

---

## 6. Los punteros colgantes del delta X (7 de octubre de 2025, nunca resuelto)

**El bug**: el sistema de derivadas guarda un puntero al delta X respecto del cual se deriva, y
la implementación a veces no lo inicializaba ni a null sino a cualquier cosa, con el resultado
de que aparecían matrices de una cantidad absurda de filas y columnas de la nada.

**Qué impedía**: la confianza en toda la maquinaria de derivación, justo cuando el paper
dependía de ella. En el mismo audio aparece un segundo síntoma que nunca cerró: una diferencia
entre la log-likelihood calculada directamente y la calculada por el algoritmo.

**Qué lo destrabó**: nada. Lo tapó.

> "Eso más o menos lo solucioné como para que no ocurra, pero no es una implementación muy
> robusta. [...] Tendría que tapar los agujeros que tiene, por lo menos para que no tenga bugs
> muy evidentes, pero no es una buena implementación."

---

## 7. La Fisher indefinida o singular (31 de mayo de 2026)

No es un bug, es una propiedad real, y va en esta lista porque trabó más que varios bugs y
porque cambió la metodología del paper.

**El bug que no es bug**: el Hessiano numérico medido fuera del máximo puede tener autovalores
negativos, cosa legítima para un Hessiano y no para una covarianza. Y las mediciones sin
información cinética (la corriente antes de que se abran los canales) no aportan curvatura, así
que la matriz queda singular y la distorsión indefinida. Esto último ya estaba identificado el
13 de febrero de 2026, tres meses y medio antes de que lo frenara.

**Qué impedía**: medir la matriz de distorsión de MacroR, que era el término de comparación.

**Qué lo destrabó**: dos decisiones caras. Primero, medir el Hessiano en el óptimo, lo que
obligó a agregar una etapa de optimización por máxima verosimilitud que no estaba en el plan:

> "Mi plan era, en principio, yo sampleo y mido. No optimizo. Y bueno, ahí tuve que optimizar."

Todo el aparato de Gauss-Newton, grupos de réplicas, anclas theta_pool y manejo de no
convergencia sale de acá. Segundo, el 4 y 5 de julio (`b028c03`, `1c2ae6f`), definir la
distorsión con la Fisher gaussiana en lugar de la numérica, que es semidefinida positiva por
construcción, al precio de volver a correr todo.

---

## 8. Las corridas del cluster (26 al 29 de mayo de 2026)

**Los bugs**: cuatro, todos de operación y todos caros en tiempo de reloj. `srun` fijaba el
trabajo OpenMP a un solo núcleo (`679cdcd`). El armado de módulos mataba el sbatch en dos
segundos (`673c44f`, resuelto guardando el setup con `set +e`). El TIME por defecto de dos días
excedía el límite de la partición y los trabajos eran rechazados (`55979a8`). Y 2000 muestras de
bootstrap sobre 16k samples daban un trabajo de 37 horas (`d28c0f3`, vuelto a 100).

**Qué impedían**: las corridas de producción, en el momento de máxima necesidad, con el
manuscrito ya empezado.

**Qué los destrabó**: puro oficio de cluster. Van en la lista porque en días perdidos pesan
tanto como un bug numérico, y porque no dejan huella en los audios: es trabajo invisible que
solo aparece en el log de git.

---

## 9. Encontrados y no arreglados

`Lapack_SymmPosDef_inv` devuelve basura fuera de la diagonal al invertir una matriz simétrica
definida positiva. El workaround es el patrón cholesky + inv(L) + XXT. La suite de tests de
álgebra lineal que lo cubriría está pendiente. Es de la misma familia que el número 1: una
rutina de LAPACK y una convención de triángulos que no se hablan.

En el volcado de derivadas de la likelihood, cada fila de evolución se escribe dos veces (una
copia de segmento en blanco y una de segmento cero). Eso produjo un "score al doble" que era
puro doble conteo, hasta que se dedujo por (sim, sample). El arreglo en C++ está diferido.

En el análisis en R, colapsar Dconf y Bconf al borde del intervalo de confianza
(`figure_4_common.R:200`) sirve para los mapas y genera máximos espurios en las líneas. Es un
ejemplo de una clase entera: bugs en la capa de lectura de los datos, donde no hay tests.

---

## Lo que se aprende del conjunto

**Los bugs caros no rompen, sesgan.** Los tres primeros de la lista por costo real (el triángulo
de la simétrica, la fórmula de gvar, el N·ms) no hicieron fallar nada. Corrieron, dieron
números, entraron en gráficos y en conclusiones. El costo no se mide en horas de debugging sino
en meses de razonamiento sobre resultados que eran en parte artefacto. El caso extremo es el
número 1: durante cuatro meses el instrumento de medición hacía ver al algoritmo propio peor de
lo que es.

**Lo que se rompe casi siempre es Σ.** La covarianza de las probabilidades de estado es el
objeto que concentra el daño: el triángulo faltante la corrompía en los productos, el
coeficiente de confianza le metía un salto, la fórmula de gvar le erraba al término entre
estados finales. No es casualidad. Σ es lo que distingue a los algoritmos recursivos de los que
no lo son, y por lo tanto es donde vive la contribución del paper, así que un error ahí sale
directamente como un resultado científico equivocado.

**La no diferenciabilidad es el otro hilo, y tiene tres episodios.** Los autovectores en octubre
de 2025; el 29 de diciembre de 2025, donde `5db2beb` (titulado "figure 2 bug in calc_P") elimina
la parte positiva suavizada de `to_Transition_Probability` y con eso reintroduce el quiebre en
el camino de derivadas; y el `if` sobre el signo de D en junio de 2026. El mismo parche,
descartar los puntos donde la diferencia finita cambia de signo, aparece como remedio en dos de
las tres ocasiones, con ocho meses de diferencia.

**Las invariantes que el tipo no puede defender son las que se violan.** `SymmetricMatrix`
significaba "medio triángulo es el dato", y bastó con que un consumidor legítimo (una rutina de
LAPACK de propósito general) no conociera esa convención para que media matriz pasara a valer
cero sin un solo warning de compilación. De ahí sale Luthier como proyecto, y de ahí sale la
regla de que el significado tiene que vivir en el tipo y no en un comentario.

**El algoritmo exacto es la mejor herramienta de debugging.** Micro IR entró como argumento
científico y destrabó el bug de ecuaciones, porque un algoritmo que tiene que dar identidad
exacta convierte cualquier desvío en evidencia. Patrón reusable: si hay una versión del cálculo
sin aproximación, aunque sea impracticablemente lenta, vale como oráculo.

**Los bugs silenciosos duran lo que tarda alguien en tener que explicarlos.** El término N·ms
sobrevivió 241 días porque nada fallaba. Lo mató la obligación de escribir en prosa por qué un
algoritmo se comportaba de cierta manera.

**Los detectores que se construyen para cazar un bug se quedan.** Las canarias de
desplazamiento y de probabilidad (mayo), los comandos para localizar saltos en la derivada
segunda (junio), el volcado por sample: todos nacieron de una cacería puntual y quedaron como
capa permanente. Es la parte del costo que se recupera.
