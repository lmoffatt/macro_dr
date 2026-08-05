# El paper contado desde sus figuras

Este relato se arma al revés que los otros dos. En vez de seguir el calendario, tira del hilo
desde lo que sobrevivió: cada figura publicada, hacia atrás hasta los datos que la alimentan,
la corrida que los produjo y los commits que hicieron ese número posible o correcto. Lo que no
llega a una figura no aparece, y eso es deliberado: es la única manera de separar el trabajo
que quedó del trabajo que se hizo.

Los otros dos documentos son [relato_audios_commits.md](relato_audios_commits.md), que es
cronológico y pareado, y [bugs_que_trabaron.md](bugs_que_trabaron.md), que ordena por costo.
Este es el que tiene el criterio de relevancia más duro.

Vocabulario: FIM es la matriz de información de Fisher; N_ch es el número de canales; SE es el
error estándar; θ_sim es el punto de simulación y θ_pool el máximo de verosimilitud conjunto;
LSE es cuadrados mínimos no lineales.

---

## 1. Lo que muestran las seis figuras

Puesto en fila, el argumento del paper es este.

**Figura 1** muestra un paso del filtro sobre una sola grabación de 20 canales, seis intervalos
de 2 ms a 50 kHz, en cuatro columnas: LSE, NR, R e IR. Cuatro filas: probabilidad de apertura
previa, corriente predicha con su dispersión y la innovación, posterior (en LSE y NR dice "no
update, open loop") y log-verosimilitud acumulada. Es la explicación del mecanismo.

**Figura 2** son nubes de máxima verosimilitud sobre 10.000 grabaciones, con tres elipses al
95%: la empírica, la de Fisher y la corregida tipo sándwich. Siete columnas, LSE, NR, INR, R,
MR, VR e IR. En la primera fila los factores de sobreconfianza son 15, 15, 13, 1,3, 2,0, 2,2 y
1,0. Ahí está la tesis en una línea: cuadrados mínimos y los no recursivos subestiman el error
entre trece y quince veces, IR da uno.

**Figura 3** pone a los siete miembros a puntuar el mismo ensemble de 1000 grabaciones
simuladas, sin ajustar nada, en siete filas: corriente predicha, residuo estandarizado al
cuadrado, información por intervalo, sesgo del score, varianza del score sobre la información,
lo mismo para las sumas, y la autocorrelación del score y del residuo. Es el mecanismo medido.

**Figura 4** es el plano de diseño, con el sesgo inducido por la distorsión y la distorsión de
información sobre una sola escala logarítmica de factor, en columnas anidadas por momento,
miembro y parámetro, y filas por número de canales. Es dónde falla cada uno. Tiene quince
suplementos.

**Figura 5** responde por cuál de los dos momentos descartados falla IR. El hallazgo es fino:
en la peor celda, con 10 canales y ruido 0,005, k_off sub-reporta el error estándar por un
factor 1,62 mientras N_ch lo sobre-reporta por 0,65, en la misma celda.

**Figura 6** es el mapa de recomendación, el plano de número de canales contra ruido partido en
cinco regiones por fronteras ajustadas como leyes de potencia, con las preparaciones
experimentales superpuestas. Es la figura que cierra el argumento editorial.

---

## 2. La cadena de cada figura

### Figura 1

Lee archivos sueltos, `figure_1_likelihood_diagnostic_{LSE,NR,R,IR}.csv` y
`figure_1_simulation.csv`, todos del 23 de julio de 2026 a las 14:38. No llevan fila de
procedencia, así que el hash del binario no se recupera del dato. Salen de
`ops/local/figure_1.macroir`, corrido a mano, entrada de ledger `run-20260723-143807`.

El commit que la hace posible en su forma actual es `a416eff` (20 de julio de 2026, "figure 1
lse"), que rutea el guard de `nonlinearsqr` en `calculate_mlikelihood_diagnostics`; el
comentario "ROUTED (was guarded)" sigue en `src/core/likelihood.cpp:913-917`. Sin ese commit no
hay columna LSE. Antes, el diagnóstico como tal nace en `1f917b9` y `b5798a1` (13 y 18 de
diciembre de 2025), que es cuando aparece `calc_likelihood_diagnostic`.

Verificación que vale la pena registrar: en los datos, la `y_var` de LSE es una constante única
(11,780) contra cinco valores distintos en NR y seis en IR. La banda plana que describe la
caption es real, no un artefacto de dibujo.

### Figura 2

Lee `_mle_cloud_runs.csv` para las nubes, `_battery_pool_G.csv` para las elipses ancladas en
θ_pool y `_battery_sim_G.csv` para el marcador de sesgo anclado en θ_sim. La resolución es por
camino de búsqueda con precedencia `1f7138b → 1c2ae6f → 0ffbda7 → 87889e6`, y el reparto
resultante es: INR sale de `1f7138b`, NR, R, MR e IR de `1c2ae6f`, y VR y LSE de `0ffbda7`. La
carpeta `87889e6` no se lee nunca, porque tiene `macro_NMR` y el roster pide `macro_INR`.

Dos corridas la alimentan: `ops/local/figure_3_mle_G.macroir` vía `dispatch_figure_3_G.sh` para
los seis macro, y `ops/local/figure_3_mle_LSE.macroir` vía `dispatch_figure_3_LSE.sh` para LSE.
Las dos celdas vacías de las filas 2 y 3 tienen explicación en el script: LSE fija la corriente
unitaria y el ruido, así que su nube solo trae cuatro parámetros.

El commit que la hace correcta es `1c2ae6f` (5 de julio de 2026), que completa la familia
gaussiana empezada el día anterior en `b028c03`. El directorio de datos se llama así porque el
binario estampa su propio hash: el ancla de la figura y el commit que la hizo válida son el
mismo objeto.

### Figura 3

Lee ocho digests `.rds` que `figure_3_digest.R` extrae por `awk` en streaming de los volcados
`figure_3_time_dlik_*.csv`, de alrededor de 1 GB cada uno. Los siete vivos son del 31 de julio
de 2026 entre las 17:41 y las 17:46, de una sola invocación de
`ops/local/figure_3_time.macroir`, ledger `run-20260731-174018`, semilla 20260722.

Un detalle que cambia lo que puede decir Methods: la información por intervalo no viene del
motor. El digest la reconstruye como `dm²/y_var + 0,5·dv²/y_var²`.

La duplicación por dos de cada fila de evolución, que en su momento produjo un score al doble,
está resuelta: el `awk` deduplica con `$4=="0"`, uno a uno y sin pérdida. El predecesor
archivado lo hacía promediando, que era más frágil. Y las derivadas nulas antes del agonista
están confirmadas en el dato: en las primeras muestras solo el ruido y la línea de base tienen
score no nulo, lo que sostiene la afirmación del segundo suplemento.

### Figura 4

Lee las mismas familias de CSV que la Figura 2, con el mismo mecanismo de camino de búsqueda,
declarando cuatro carpetas. INR viene en exclusiva de `1f7138b`, 25 celdas. La grilla no está
declarada en ningún lado: `autodetect_cells()` escanea el disco por expresión regular y acepta
la celda solo si existen las dos baterías.

La corrida es `ops/local/figure_4.macroir` más `figure_4_LSE.macroir`, despachadas por
`dispatch_figure_4.sh`, una tarea por algoritmo, número de canales y ruido. Las etapas son:
máxima verosimilitud conjunta por Gauss-Newton para obtener θ_pool, predicciones de la derivada
de la log-verosimilitud en θ_sim y en θ_pool, y batería gaussiana en ambos con bootstrap de 100
réplicas y lag máximo 10.

Hay una dependencia cruzada que no está declarada: la etapa de nube de máxima verosimilitud fue
eliminada de este pipeline el 23 de julio, así que los `_mle_cloud_runs` que leen cinco de sus
suplementos vienen del otro pipeline, el de `figure_3_mle_G.macroir`.

### Figura 5

Sale de `figure_3_mle_G.macroir` vía `dispatch_figure_3_G.sh` con 10.000 simulaciones, y sus 24
columnas dibujadas están todas en `1c2ae6f`, sin costura de procedencia. Los tres componentes
son la distorsión gaussiana total, la parte por muestra y la parte de correlación, todos sobre
la batería anclada en θ_pool. Dibuja la media bootstrap, no el borde del intervalo.

El `figure_5.Rmd` actual nació el 31 de julio de 2026, cuando el anterior pasó a llamarse
`figure_5_budget.Rmd`.

### Figura 6

IR sale de `figure_3_mle_G.macroir`, con 60 archivos en `1c2ae6f` y 20 en `0ffbda7` para los
ruidos altos. LSE sale de `figure_3_mle_LSE.macroir` vía `dispatch_figure_3_LSE.sh`, 32
archivos, todos en `0ffbda7`.

El linaje de LSE es corto y está bien documentado: `67b9ce2` (18 de julio) trae el plan, la
especificación y un oráculo en R; `82b956f` (19 de julio, "finishing nonlinearsqrfit") es la
implementación, con `family_type` de tres valores; `18ead3a` (20 de julio) agrega la Fisher de
LSE por intervalo y los scripts; `a416eff` (20 de julio) la ruta de diagnóstico. El 31 de julio
`1f7138b` toca la ruta LSE solo en comentarios, así que es numéricamente idéntica.

---

## 3. El linaje del motor

Siete piezas, y ninguna nació cuando uno esperaría.

**El motor de verosimilitud.** Las banderas `recursive` y `averaging` no entraron una por una:
nacieron juntas el 1 de septiembre de 2023, en `cafb6df`, como los seis puntos de una grilla de
dos por tres, con `averaging` ya como entero. La consolidación que importa es `7aac319` (7 de
junio de 2025), que pasa el valor del objeto al tipo y convierte el despacho en `if constexpr`.
Después se agregan ejes: Taylor en `ebe866f` (4 de octubre de 2025), la familia LSE en
`82b956f` (19 de julio de 2026) y la forma de varianza que da VR en `ef71be6` (21 de julio de
2026). O sea que de los ocho miembros del roster final, seis existían como coordenadas desde
2023 y dos llegaron en la última semana y media de producción.

**Las derivadas.** Nacen en septiembre de 2023. El test pasa en `a18ffa7` (24 de octubre de
2025), lo rompe `a3e0a89` en diciembre y lo reparan tres commits del 4 de diciembre. En
`01483ad` (10 de mayo de 2026) el test se mete adentro de las corridas, con un chequeo de
derivada de Clarke en el camino de corrección de varianza.

**La Fisher gaussiana y la numérica.** La gaussiana nace en `6bfde5b` (12 de marzo de 2026), la
numérica en `30dd0e0` (29 de abril de 2026). Y acá aparece un bug que no estaba en la lista:
`a105078` (4 de junio de 2026) descubre que **el acumulador gaussiano en el camino de derivación
automática era código muerto**, porque la guarda preguntaba por un slot que no existía. Recién
ahí empieza a acumular de verdad. Cualquier número de Fisher gaussiana anterior al 4 de junio de
2026 no es un número malo, es un número que no se calculó.

**La matriz de distorsión.** Nace completa el 12 de marzo de 2026 junto con el bootstrap y las
estadísticas de momentos. Se le quita la regularización el 15 de marzo, se le agrega la familia
de correlación el 14 de abril, y el 31 de mayo se renombra de IDM a Likelihood Information
Distortion en teoría, documentos y código a la vez, lo que cambia también los nombres de columna
de los CSV.

**La optimización.** `legacy/gauss_newton.h` nace el 4 de junio de 2026, o sea siete semanas
antes del final. Su historia es de ajustes forzados por los datos: la escalera de amortiguación
el 5 de junio, el retorno "stalled" el 14 de junio porque si no los grupos ruidosos no daban
máximo, y la tolerancia adaptativa en `433ed13` (17 de junio) después de medir que el 40% de los
grupos de 10 canales agotaba el tope de iteraciones y se comía el 86% del cómputo.

**El DSL de índices.** Nace el 30 de marzo de 2026 y se integra el 1 de abril. Es lo que permite
barrer la grilla completa dentro de un proceso, y sin eso el presupuesto de cómputo no cerraba.

**La procedencia de los CSV.** Es la pieza que explica una cosa que yo había leído mal. El sello
de commit entra primero como comentario dentro del CSV, rompe el test de conteo de líneas, se
muda a un sidecar `.binary` en `5913c56` (26 de mayo de 2026), y termina como fila 1 del propio
CSV en `5d4e828` (16 de junio de 2026), porque el sidecar no sobrevive a copiar o mover.

---

## 4. La ventana de 45 días, y por qué

Las nueve carpetas de datos que alimentan las figuras van del 16 de junio al 31 de julio de
2026. No hay una sola cifra publicada anterior. El ledger local de corridas tiene 751 entradas
con los picos en diciembre de 2025 (239), abril (179) y mayo (140), contra 37 en junio y 28 en
julio.

Es tentador leer eso como que el 91% del cómputo se tiró, y en parte es cierto, porque cada
arreglo grande invalidaba lo anterior. Pero hay una razón más aburrida y más importante: los
directorios nombrados por hash **no podían existir antes del 16 de junio**, que es cuando
`5d4e828` hizo al binario consciente de su propio commit. Antes de esa fecha la salida iba a
archivos sueltos sin sello. Así que la ventana de 45 días mide dos cosas a la vez, cuánto
sobrevivió y desde cuándo se puede saber qué sobrevivió.

Dentro de esa ventana, cada ancla es una decisión, no una corrida de rutina:

| Ancla | Fecha | Qué decidió |
|---|---|---|
| `8fc274d`, `dfa842d`, `5d9b43b` | 16-17 jun | corrección de la covarianza empírica y arreglos de dispatch |
| `433ed13` | 17 jun | batería con Fisher numérica, cinco miembros, tolerancia adaptativa |
| `1c2ae6f` | 5 jul | anclar en la Fisher gaussiana |
| `87889e6` | 11 jul | completar la batería gaussiana para los miembros que no eran IR |
| `82b956f` | 19 jul | entra LSE |
| `0ffbda7` | 21 jul | entra VR, y LSE pasa a dominar |
| `1f7138b` | 31 jul | INR corregido reemplaza a NMR |

---

## 5. Los commits cuyo título subestima el diff

Este es el resultado metodológico del ejercicio, y explica por qué un relato hecho con títulos
de commit falla. Son siete confirmados:

- `a3e0a89` (2 dic 2025) "theoretical results". Cuerpo vacío, 3300 líneas de código, `qmodel.h`
  reescrita entera, pierde el término `N·ms`.
- `5db2beb` (29 dic 2025) "figure 2 bug in calc_P". Además del arreglo, elimina la parte
  positiva suavizada en `to_Transition_Probability`, o sea reintroduce el quiebre en el camino de
  derivadas, y de paso mete el manejo de semillas (`SeedNumber`, `calc_seed`).
- `347be45` (2 ene 2026) "better SymmetricMatrix handling". Introduce el almacenamiento de medio
  triángulo. Ventana de daño hasta el 8 de mayo.
- `9b61728` (18 abr 2026). Crea `io/fingerprint.h`, un subsistema de identidad de datos que el
  asunto no menciona.
- `30dd0e0` (29 abr 2026) "micro_R MR IR working". Cuerpo vacío, 4576 líneas. Además de micro,
  crea la **Fisher numérica** y la distorsión gaussiana. Dos piezas macro bajo un asunto que dice
  micro.
- `a7d6ff9` (7 may 2026) "micro_ir monoid implementation". Cuerpo vacío, 8165 líneas. Extrae
  `qmodel_types.h` y reescribe `qmodel.h`: una reorganización del motor macro anunciada como una
  funcionalidad micro.
- `1c2ae6f` (5 jul 2026) "gaussian also for corrected covariance". Cuerpo vacío, y es la decisión
  de anclaje que le da nombre a la carpeta de datos de la mitad de las figuras.

El corte es el 8 de mayo de 2026. Antes, asunto de una línea y cuerpo vacío incluso para 8000
líneas. Después, cuerpos que declaran alcance, radio de daño y ventana temporal.

---

## 6. Lo que este ejercicio dejó para arreglar

Tres categorías, en orden de gravedad.

### Sustantivo

**La frontera de LSE en la Figura 6 está trazada con la Fisher gaussiana.** Los CSV de LSE
traen la columna de Fisher numérica vacía: el script omite esa etapa a propósito y el dispatcher
lo declara. No hay ninguna batería de LSE anclada en la Fisher numérica en disco. La variante
`figure_3_mle_LSE_numfisher.macroir` es del 1 de agosto de 2026, posterior al dibujo, y su
cabecera dice que sin ella no se pudo zanjar si el sándwich es correcto para cuadrados mínimos.
Esto contradice lo que el propio autor concluyó el 2 de agosto, que para LSE hay que usar la
numérica porque la gaussiana da mal.

### Texto que va publicado

Las captions de las Figuras 1 y 2 describen figuras distintas de las que existen. La de la 1
describe R, MR, VR e IR cuando el PDF tiene LSE, NR, R e IR. La de la 2 describe cuatro
algoritmos, dos filas y una parte B que ya se mudó al suplemento, cuando el PDF tiene siete
columnas y tres filas. Los números que cita sí coinciden.

`Figure_3_caption.md` dice que los volcados son del 22 de julio y que solo INR es del 31. Los
siete son del 31, de una sola corrida. La semilla sí es 20260722, así que el ensemble es el
mismo, pero la corrida no.

Las captions de la Figura 3 y `provenance.md` declaran motor `0ffbda7`; el sello real de los
volcados es `1f7138b-dirty`.

`decisions.md` define la Figura 5 como el presupuesto de información entre métodos. La figura
que hay es solo IR.

Si el texto dice que la Figura 6 lleva condiciones experimentales reales, hay que ajustarlo: son
rangos derivados de bibliografía, y el único dato experimental verdadero está apagado desde el
27 de julio por decisión del autor.

### Procedencia y reproducibilidad

`provenance.md` no tiene la corrida de la Figura 4, que es la figura central, ni la de LSE, ni
la re-corrida de INR. Su trampa 3 está invertida respecto del estado actual, y su afirmación de
que `1c2ae6f` alimenta solo notebooks exploratorios es falsa.

De tres carpetas-ancla (`0ffbda7`, `87889e6`, `1f7138b`) no hay entrada de ledger, así que no se
sabe qué máquina ni qué invocación las produjo.

`seed = 0` es el centinela de aleatorio y la semilla resuelta nunca se registra, así que ningún
ensemble se regenera bit a bit. El bootstrap sí es determinista.

Los CSV de la Figura 1 no llevan fila de procedencia.

`figure_5_budget.Rmd` rotula el label crudo de ruido como si fuera adimensional, un factor 10
fuera de la convención que las Figuras 5 y 6 sí respetan.

Hay una celda duplicada entre carpetas con valores distintos, `macro_IR` con 5 canales y label
0,1, que vale 1,6916 en `0ffbda7` y 1,7012 en `87889e6`. El chunk que dibuja el mapa la
deduplica; otro chunk no. Es el riesgo del camino de búsqueda materializándose.

---

## 7. Qué se ve desde acá que no se veía desde el calendario

**El paper se terminó de definir en las últimas tres semanas.** LSE entra el 19 de julio, VR el
21, INR se corrige el 31. Dos de los ocho miembros del roster y la corrección de un tercero
caben en doce días. El giro editorial de comparar contra lo que usa la gente no fue una idea que
se implementó después con calma: se implementó, se corrió y se dibujó en una semana y media.

**El commit que habilita la tesis es de una línea.** `a416eff`, "figure 1 lse", que rutea un
guard. Sin él no hay columna LSE en las Figuras 1 y 2, o sea que no hay comparación contra
cuadrados mínimos, o sea que no hay paper para eLife. En el relato cronológico es una línea más
de una semana ocupada.

**Las piezas de análisis son todas de los últimos cuatro meses.** La matriz de distorsión es del
12 de marzo de 2026, la Fisher numérica del 29 de abril, el Gauss-Newton del 4 de junio, el
anclaje gaussiano del 5 de julio. Los diez meses anteriores construyeron el motor de
verosimilitud y el lenguaje para manejarlo, que son condición necesaria y no aparecen en ninguna
figura.

**El instrumento de medida se terminó después que el objeto medido.** El algoritmo estaba desde
2023 como coordenadas de una grilla. Lo que faltaba, y lo que llevó el año, fue el aparato para
decidir si funciona.
