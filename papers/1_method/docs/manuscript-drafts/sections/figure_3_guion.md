# Guion de la Figura 3. Escrito 2026-08-12

Fuente de todos los números: el chunk `caption-numbers` de
`projects/eLife_2025/figures/paper_both/figure_3.Rmd`, leído en `figure_3.html` del 2026-08-05
(bloques 5, 9 y 11 del `<pre>`). Todo lo que sigue está impreso ahí y ninguno se lee del dibujo.

**Decisión del 2026-08-12: se queda en n = 1.000 grabaciones.** La opción de 10.000 quedó
postergada, no descartada; lo que hay que saber para reabrirla está al final.

---

# CUARTA PASADA, 2026-08-27. Auditoría de la narrativa y ocho reparaciones

No es una reescritura de eje: el eje de la tercera pasada (qué tipo de información da cada
variable) se queda. Lo que se reparó son defectos de contenido que las pasadas del 2026-08-25/26
dejaron al arreglar premisas sueltas sin releer el párrafo entero contra la tabla de la figura.

**HALLAZGO NUEVO, y cierra el ítem que este guion tenía abierto de otra forma: EL ANCLA.**
El 0,724 de la fila E (Ē al ancla `sim`) y el "1.09 to 0.998 for R" de la leyenda de
`fig:plane`--supp 1 (ancla `pool`) son el MISMO objeto en dos puntos de evaluación, y los dos
están bien. Medido de `figure_4_source_data_distortion.csv` (columna `m`, no `Dconf`), celda
Num_ch 100 / noise 0,1 / interval 0,1 / param_index 1, `sample x corr = total`:

| miembro | @ theta_sim | @ theta_pool |
|---|---|---|
| R | 0,757 x 1,462 = 1,108 | **1,087 x 1,397 = 1,513** |
| MR | 0,572 x 2,007 = 1,149 | 1,002 x 1,827 = 1,816 |
| VR | 0,767 x 1,620 = 1,240 | 1,041 x 1,498 = 1,559 |
| NR | 1,445 x 19,64 = 14,61 | 1,019 x 19,47 = 17,32 |
| INR | 1,140 x 20,72 = 23,67 | 1,151 x 20,71 = 23,91 |
| IR | 1,115 x 1,011 = 1,130 | 1,113 x 1,010 = 1,127 |
| LSE | 0,961 x 13,75 = 13,35 | 0,967 x 13,70 = 13,38 |

**El patrón ES el hallazgo:** los dos centrados leen igual en las dos anclas a la tercera cifra;
todos los desplazados no, y su factor por intervalo pasa de abajo de uno en la verdad a uno o más
en su propio óptimo. O sea: **la cancelación de R es una propiedad del punto de evaluación, no del
miembro.** En theta_pool, que es donde caen las estimaciones y donde lee la Figura 4B, R no cancela
nada: 1,51 con los dos factores arriba de uno. Refuerza el "sólo IR" en vez de debilitarlo, y es el
único lugar del paper donde la calibración aparente de R queda acotada. Escrito en el cuerpo y en
la leyenda de `fig:plane`--supp 1, con la tabla entera en un bloque `%` de `04_results.tex`.

**Las otras siete.**

1. **Universal falso, párrafo 4.** "per interval ... for every rung of the ladder, while accumulated
   it fails by an order of magnitude": el orden de magnitud es de LSE/ILSE, NR e INR y de nadie más
   (R 1,074/0,946, MR 1,188, VR 1,198/1,427, IR 1,080/0,971 en la source data 1). Es la misma falla
   que la auditoría del 2026-08-13, con el signo cambiado. El cuantificador ahora va con su roster.
2. **Dirección invertida, párrafo 5.** "a per-sample variance it under-reports by 28%": Ē = 0,724
   es CONSERVADOR por intervalo (barra 18% demasiado ancha), y como los dos factores se enunciaban
   como errores del mismo signo, la cancelación no se veía. Reescrito con los signos opuestos.
3. **La disociación central, que estaba medida y no dicha.** INR quinto en logL, 15 nats arriba de
   NR, y el peor de los ocho en J_T/F_T (20,84 contra 13,63). Es la respuesta a "¿para qué las filas
   B a G si ya tengo la logL?". Una oración al cierre del párrafo 1.
4. **La cláusula de LSE, dimensionada.** Declaraba el problema y no su tamaño: ahora dice que el
   perfil sólo puede quedar ARRIBA de la densidad, así que los 119 nats son cota inferior, y que el
   desvío medido es 2,2 nats (`fig3_tc.R`). La reserva pasó de reserva a argumento. La identidad KL
   queda acotada a los seis miembros de parámetros fijos.
5. **Tres escalas en una oración, párrafo 3.** "chance level" / "twice chance" / "four times it" /
   fracciones. Todas las fracciones ahora, con el 0,05 dicho una vez.
6. **`informative` se usaba y no se definía** en ningún lado del manuscrito (grep sobre las
   secciones y sobre Methods). Definido en el cuerpo al primer uso y en la leyenda, donde además
   "the information reaches zero" pasó a "reaches the numerical floor".
7. **`-0,004 ± 0,004` era el SE bootstrap** donde todo el resto del paper usa IC de 95%. Con ± SE el
   lector calcula [-0,0086, 0] y ve el cero excluido, que es lo contrario de lo que afirma la
   oración. Puesto el IC (-0,013 a +0,004).

**Además, sin agregar palabras:** el párrafo 2 recupera los NULOS de cada objeto y la línea de qué
filas están disponibles en un registro real (era el plan de la tercera pasada, la pasada del 08-25
lo había dejado en lista de filas); la ACF del residuo entra con sus números (0,85 / 0,51 / 0,12 /
-0,01), que vivían sólo en un comentario; y el párrafo 7 se reordenó para cerrar en el doble conteo
y no en un puntero a suplemento.

**Título, cuarta versión:** "The error bar fails for two reasons, one within a sample and one across
them, and they can cancel". El anterior decía "two independent reasons" y las dos fuentes no son
independientes: la figura las mide como los dos FACTORES de un producto y el hallazgo es que en R
tienen signos opuestos.

**PRESUPUESTO, y queda abierto.** 834 -> 1.041 palabras. El plan de la tercera pasada eran 892. Los
~200 son los tres resultados nuevos (la cota de LSE con su tamaño, la disociación logL/calibración,
el ancla) más los nulos del párrafo 2. El beat 7 ya no tiene la grasa que este guion le suponía (el
conteo de pasos informativos salió en una pasada anterior). Si hay que cerrar contra 892, los
candidatos son el ancla (~55, pero es el hallazgo), la cota de LSE (~30) o el mecanismo del washout
en el párrafo 7 (~33). Decisión de Luciano, no tomada.

---

## La escalera de logL se descompone en las otras dos filas (2026-08-12, sin escribir todavía)

Material nuevo, medido pero NO incorporado a la prosa. Script: `papers/1_method/decisions/recompute/fig3_kl.R`.

Tres identidades encadenadas. (1) E_p[log q] = −H(p) − KL(p‖q), y como el promedio es sobre
grabaciones del simulador exacto, las DIFERENCIAS de la escalera son diferencias de KL al proceso
verdadero, con −H(p) común a los ocho. **Corrige lo que se había dicho el mismo día, que la
likelihood sólo se podía ordenar y no testear: tiene criterio de verdad, salvo un corrimiento
compartido.** (2) Para una gaussiana con la varianza errada por un factor r, KL = ½(r − 1 − log r),
la pérdida de Stein, y ese r es la fila B. (3) Si el miembro asume intervalos independientes,
KL(p‖∏q_t) = Σ_t KL(p_t‖q_t) + [Σ_t H(p_t) − H(p)], y el segundo término es la correlación total de
la secuencia observada, propiedad del DATO y no del algoritmo.

Medido, costo marginal predicho desde la fila B sola y el resto del hueco contra IR:

| | logL | hueco | marginal (fila B) | dependencia |
|---|---|---|---|---|
| LSE / ILSE | −262,30 | 118,85 | 49,09 | **69,76** |
| NR | −230,15 | 86,70 | 14,81 | **71,89** |
| INR | −215,37 | 71,93 | 0,03 | **71,90** |
| R | −147,17 | 3,72 | 2,09 | 1,64 |
| MR | −153,87 | 10,43 | 4,51 | 5,92 |
| VR | −148,94 | 5,50 | 2,41 | 3,08 |
| IR | −143,45 | 0 | 0,05 | −0,05 |

Dos chequeos que no son ajuste: NR menos INR observado da 14,78 nats y la fila B sola lo predice en
14,78 a cuatro cifras (el término de dependencia se cancela entre dos miembros que ambos asumen
independencia, así que el residuo tiene que cargarlo entero), y los tres de lazo abierto dan 69,8,
71,9 y 71,9 de dependencia, tres algoritmos distintos al 3% sobre un número que la teoría dice que
es del proceso. Son **72 nats por grabación de cien muestras** de correlación temporal, que INR paga
enteros teniendo los marginales perfectos. Y el 81% de la penalidad marginal de NR (12,05 de 14,81)
sale de UN intervalo, el del escalón, que es el 28,45 de la fila B.

LO QUE NO DICE: la implicación es de una mano. Distorsión implica costo de verosimilitud; el
recíproco es falso, porque una media mal especificada cuesta KL y deja J = F intacto. La likelihood
ve los dos defectos y no los separa, que es para lo que sirve el resto de la figura. Y la moneda del
costo marginal es el cociente de varianzas del RESIDUO (fila B), no el J/F por parámetro (fila E):
R lee 0,782 en B y 0,724 en E, LSE lee 1,000 en B y 1,486 en E.

EL PISO, VERIFICADO CONTRA EL DATO y no ya inferido por resta (`papers/1_method/decisions/recompute/fig3_tc.R`). La forma testeable
no necesita H(p), que se cancela: **−logL_m = Σ_t H(p_t) + Σ_t KL(p_t‖q_t)** para todo miembro de
independencia. El primer sumando sale del dato crudo (cien marginales de la corriente observada,
mil muestras cada uno), el segundo de los momentos del residuo del miembro. Son caminos
independientes.

| | −logL | Σ H(p_t) | Σ KL(p_t‖q_t) | suma | dif |
|---|---|---|---|---|---|
| NR | 230,15 | 215,34 | 14,87 | 230,21 | −0,06 |
| INR | 215,37 | 215,34 | 0,09 | 215,42 | −0,05 |
| LSE | 262,30 | 215,34 | 49,17 | 264,51 | **−2,21** |

**La logL total de INR ES la suma de las entropías marginales verdaderas**, 215,37 contra 215,34
medido del dato: sus marginales son los verdaderos a 0,03 nats sobre cien intervalos y paga 72
igual. Controles: los marginales son gaussianos a 0,16 nats (Vasicek contra plug-in gaussiano), y
los ocho reportan exactamente −½log(2πv) − ½r² con su propia v, verificado a 1e-13. Y lo medido es
TC − KL(p‖q_IR), o sea 71,89 es COTA INFERIOR de la correlación total; lo que falta es la
divergencia propia de IR, que ninguna de estas cuentas ve.

**HALLAZGO, y es un problema del cuerpo del paper.** LSE falla el test por 2,21 nats porque **su
varianza predictiva varía entre grabaciones** (coeficiente de variación 0,478 en todos los
intervalos, contra 0,000 de NR e INR): estima una escala de ruido por grabación. Eso explica de una
el r̄² exactamente 1,0000, la pregunta que había quedado abierta sobre el ruido marginalizado, y
este desvío. La consecuencia: **la entrada de LSE en la escalera de costo no es el mismo tipo de
objeto que las otras siete**, que son densidades en parámetros fijos, mientras la suya es un perfil
sobre un parámetro ajustado de cada grabación. El orden no se cae (119 nats de hueco contra 2,2 de
discrepancia) pero la comparación tiene que declararse donde se anuncia la escalera. Para no
confundir: IR también tiene varianza variable entre grabaciones (cv 0,093) y ahí es condicionamiento
en el pasado, que es lo correcto; la de LSE es lazo abierto con un parámetro estimado.

CERRADO de los dos pendientes: el sesgo de media es despreciable (0,04 a 0,09 nats en todos los
miembros, separado en parte de varianza y parte de media), así que la fórmula de Stein aguanta. SIGUE
ABIERTO: buscar si la descomposición ya tiene nombre en la literatura, y explicar el residuo del
bloque recursivo (R 1,59, MR 5,87, VR 3,03), que no es media ni varianza del residuo; la sospecha es
Jensen, la dispersión del cociente de varianzas entre grabaciones dentro de un mismo intervalo, que
es justo lo que el condicionamiento en los dos extremos elimina.

---

## La figura fila por fila (agregado 2026-08-12, sobre la lectura de Luciano)

Una entrada por fila: qué afirma, el número que la carga, y la trampa. Todo verificado contra los
digests el 2026-08-12; lo que no está impreso hoy por el chunk `caption-numbers` está marcado.

**EL ORDEN, antes de las filas.** El orden de las verosimilitudes sigue al de las distorsiones con
DOS inversiones, no una. Distorsión de Figura 2 (razón de áreas del par cinético): IR 1,02 / R 1,32
/ MR 1,97 / VR 2,18 / INR 12,6 / LSE 14,7 / ILSE 14,8 / NR 14,8. logL: IR −143,45 / R −147,17 /
VR −148,94 / MR −153,87 / INR −215,37 / NR −230,15 / LSE e ILSE −262,30. Se invierten MR con VR y
NR con LSE (empate en distorsión, 32 nats en likelihood). **TRAMPA: depende de qué distorsión.** Con
el ratio acumulado de esta figura, INR es el PEOR de los ocho (20,84 contra 13,63 de NR) estando
quinto en likelihood y 15 nats por encima de NR. Esa inversión es la disociación central del paper,
así que el párrafo no puede abrir con "el orden se parece" sin decir dónde se rompe.

**A, la traza.** Afirma que la partición es 4+4 y no la escalera: LSE, ILSE, NR e INR predicen la
media determinista y R, MR, VR e IR siguen al registro. Es la lectura visual de los dos bloques de
escala que la figura ya usa. La anotación es la escalera de costo, span 118,85 nats.

**B, el residuo.** Afirma que sólo los dos miembros de intervalo tienen la varianza que dicen. Por
segmento (pulso = muestras 20 a 59, agonista leído de la columna y no supuesto):

| | 1-19 | s=20 | 21-30 | 31-59 | 60-100 | todo |
|---|---|---|---|---|---|---|
| LSE / ILSE | 0,010 | 0,27 | 1,65 | 1,91 | 0,69 | **1,0000** |
| NR | 0,997 | **28,45** | 1,36 | 0,96 | 0,64 | 1,153 |
| INR | 0,997 | 0,95 | 0,99 | 1,00 | 1,00 | 0,999 |
| R | 0,997 | 0,61 | 0,68 | 0,67 | 0,78 | 0,782 |
| MR | 0,997 | 0,95 | 0,63 | 0,58 | 0,63 | 0,690 |
| VR | 0,997 | **3,06** | 0,77 | 0,72 | 0,74 | 0,812 |
| IR | 0,997 | 0,95 | 1,01 | 1,00 | 1,00 | 0,999 |

DOS COSAS QUE NO ESTABAN DICHAS. El exceso de NR no está difuso al comienzo del pulso: está en UN
intervalo, el que recibe el escalón, y ahí vale 28,45 (VR falla en el mismo, 3,06). Es el fallo del
promediado de ventana en su forma pura, la ocupación cambia dentro de la ventana. Y LSE recorre un
rango de 190 veces, 0,010 antes del agonista contra 1,91 en el pulso, que promedia exactamente a
1,0000: el nivel de ruido marginalizado fija la media en uno y deja el curso temporal libre. Eso
cierra la pregunta que había quedado abierta sobre el 1,0000.

**C, la información.** Afirma dónde se mide cada parámetro y por cuánto tiempo. Al retirarse el
agonista la corriente ya no depende del número total de canales sino de los que quedan abiertos, que
un miembro que condicionó ya conoce, así que no aprende más ni sobre N_ch ni sobre **k_on**.
CUIDADO CON EL NOMBRE: la que sobrevive al lavado es k_off, la de cierre, justamente por la caída;
la que muere con el agonista es la de binding. Números: F(N_ch) de IR va de 0,36 a 2,4e−14 en nueve
intervalos mientras NR sólo baja de 400 a 91, y los pasos informativos son 80 para NR e INR contra
43 a 46 para los cuatro recursivos. F(k_off) SIGUE a la corriente sin ser proporcional a ella:
correlación log-log 0,96 con exponente 0,57.

**D, el sesgo.** Por intervalo, sólo IR está en cobertura nominal (0,0375 y 0,0652) e INR cerca
(0,10 y 0,11). Acumulado y con la matriz entera, INR e IR cubren el cero en los CUATRO parámetros
identificados y ningún otro lo hace. Trampa: "no distinto de cero" depende de las mil grabaciones.

**E, el ratio por intervalo.** Sólo los dos de intervalo están en uno (INR 1,096 y 1,006; IR 1,069 y
1,015); R, MR y VR reportan MENOS varianza que la que su Fisher indica (0,58 a 0,87). Faltan dos
columnas en la lectura: LSE está por ENCIMA (1,486) y NR es mixto (0,810 en k_off, 2,421 en N_ch).
El aumento durante el decaimiento es real y sistemático: IR va de 0,999 en el pulso a 1,187 y 1,310
en los últimos veinte intervalos, así que ni IR está plano en el tiempo.

**F, el acumulado.** Hasta 20,84. TRAMPA: con este estadístico solo, R también cubre el uno (1,074
[0,991, 1,17] y 0,946 [0,863, 1,03]), así que "sólo IR" no se sostiene sobre la fila F. Lo que
separa a IR es la factorización: IR tiene los dos factores en uno, R llega cancelando.

**G, la memoria.** Cero para IR, lenta para los de lazo abierto, rápida para los recursivos, y se
puede cuantificar: la ACF cae por debajo de 0,1 en el lag 9 a 12 (LSE, NR, INR), en el lag 2 (R, MR,
VR) y en el lag 1 (IR). Segunda serie, la del residuo, que es lo único de la figura computable en un
registro real: 0,846 / 0,510 / 0,634 / 0,118 / 0,217 / 0,168 / −0,010. TRAMPA GRANDE: la curva
dibujada no reconstruye la fila F. 1+2Σρ da 7,85 para NR donde la inflación medida es 16,83, porque
el score no es estacionario y el promedio de ACFs por grabación no pesa los intervalos donde la
correlación importa.

**Lo que ninguna fila dice por sí sola y hay que decir en prosa:** que LSE e ILSE son la misma
columna dos veces (la razón de que haya ocho), que E × G = F y que R llega a uno cancelando, que el
desplazamiento acumulado vive en el par de amplitud y no en las cinéticas, la diferencia entre
magnitud y significancia, que se dibujan dos parámetros y se miden cuatro, y que N_ch está excluido
de las dos columnas de mínimos cuadrados por degeneración y no por falta de datos.

---

## El trabajo, lo que sólo esta figura hace

Figura 2 mide el resultado en una celda, Figura 4 el resultado en todo el plano. Figura 3 es la
única que mide el MECANISMO: por qué hay sesgo y por qué hay distorsión, que son dos causas
distintas con dos reparaciones distintas. Lo que no sea eso pertenece a otra figura.

No le toca: volver a probar que MR y VR fallan (`04_results.tex:250-252`, con áreas de elipse), ni
argumentar memoria contra distorsión sobre el plano (Figura 5 y su suplemento 1, correlaciones de
rango +0,91 a +0,96), ni contar la escalera de costo como resultado (es contexto).

---

## La arquitectura: dos causas, dos filas, dos interruptores

| causa | fila | interruptor que la repara | las cuatro esquinas |
|---|---|---|---|
| mala estimación por muestra → **sesgo** | D | promediado de intervalo | NR 0,45 → INR 0,10; R 0,50 → IR 0,0375 |
| autocorrelación → **distorsión** | F | recursión | NR 13,6 → R 1,07; INR 20,8 → IR 1,08 |

Ese es el centro y hoy no está escrito así. El manuscrito cuenta el factorial con r̄²_std y J/F, o
sea con un testigo del dato y la cantidad correcta. Contarlo con D y F lo cuenta con las dos
cantidades que SON el sesgo y la distorsión, y las dos esquinas de cada fila salen de la misma
tabla impresa.

Los testigos quedan en su lugar: fila B es el testigo a nivel de dato (fracción de pasos cuyo IC
excluye el 1: INR 0,00 e IR 0,04 contra LSE/ILSE 0,98, R 0,78, MR 0,78, VR 0,79, NR 0,52) y fila G
es el testigo del mecanismo.

---

## La identidad que ordena el roster de filas

Var(Σ s_t)/Σ F_t = Ē × (1 + 2 Σ_{t<u} cov / Σ Var).

**La fila F es el producto de la fila E por un factor de inflación que la fila G mide.** No son
tres mediciones independientes, son dos factores y su producto, y la figura muestra los tres para
que el lector cierre la cuenta.

**MEDIDO 2026-08-12** (`papers/1_method/decisions/recompute/fig3_check.R`, reproduciendo el prep, la máscara y el bootstrap de la
notebook; el producto reproduce el acumulado a todos los dígitos, que es la identidad y sirve de
chequeo del código). k_off salvo donde dice N_ch:

| miembro | Ē (por intervalo, pesada) | inflación (correlación) | producto = J_T/F_T |
|---|---|---|---|
| LSE / ILSE | 1,486 | 10,07 | 14,96 |
| NR | 0,810 | 16,83 | 13,63 |
| NR, N_ch | 2,421 | 5,21 | 12,62 |
| INR | 1,096 | 19,02 | 20,84 |
| INR, N_ch | 1,006 | 11,76 | 11,83 |
| R | 0,724 | 1,483 | 1,074 |
| MR | 0,578 | 2,055 | 1,189 |
| VR | 0,722 | 1,659 | 1,198 |
| IR | 1,069 | 1,010 | 1,080 |

**El hallazgo no es la identidad, es lo que muestra: R llega a 1,07 CANCELANDO dos errores**, un
28% de sub-reporte por intervalo contra un 48% de inflación por correlación, y MR lo mismo con 0,58
contra 2,05. IR es el único con los dos factores en uno por separado. Eso explica de una lo que el
suplemento S1 reporta como rareza (la calibración por intervalo no es monótona en la escalera: INR
está mejor por intervalo que los tres recursivos que cuestan más), y cambia la lectura del
factorial: la recursión sola no "repara la distorsión", produce un sub-reporte por intervalo que
compensa su propia inflación.

**La versión AR(1) murió.** Predecir la inflación con sólo el ρ de lag 1 anda en los recursivos
(0,7% a 8%) y falla en los de lazo abierto (18% a 46%: NR N_ch predice 2,83 contra 5,21 medido),
porque su ACF no decae geométricamente. No se imprime. Lo que se imprime es la factorización, que
es exacta y más fuerte.

---

## El panel que falta, y la forma barata de ponerlo

La fila D es por intervalo. El sesgo del estimador lo determina E[Σ s_t]/√(Σ F_t), que no está en
ningún panel. Sin eso, "sólo IR tiene cero en todos lados" queda sin conectar con los marcadores
de Figura 2, y si los sesgos por intervalo se cancelan al sumar es una pregunta que el paper deja
abierta teniendo el dato en el disco.

**MEDIDO 2026-08-12, y CORREGIDO el mismo día.** Acumular el score está bien; dividir por √F_aa no.
El sesgo de primer orden es la ecuación VECTORIAL b = F⁻¹·E[Σ s_t] con F la matriz entera, y el
desplazamiento de un parámetro en unidades de SU error estándar marginal es b_a/√((F⁻¹)_aa).
Dividir el score acumulado de cada parámetro por la raíz de su propia diagonal trata a los otros
tres como conocidos, y acá no es una diferencia chica: en N_ch da −0,34 para R donde el marginal es
+0,82, signo incluido, porque N_ch y la corriente unitaria están fuertemente correlacionados.

z = b_a/√((F⁻¹)_aa), los cuatro parámetros identificados, IC bootstrap sobre grabaciones:

| miembro | k_on | k_off | i | N_ch |
|---|---|---|---|---|
| LSE | −0,005 | **+0,215** | +0,034 | (excluido) |
| ILSE | +0,016 | **+0,219** | +0,016 | (excluido) |
| NR | **−0,279** | **+0,538** | **−0,946** | **+0,919** |
| INR | +0,072 | +0,084 | +0,036 | −0,057 |
| R | **−0,167** | −0,055 | **−1,383** | **+0,816** |
| MR | **−0,176** | **−0,151** | **−1,956** | **+1,137** |
| VR | +0,050 | **−0,211** | **−1,091** | **+0,565** |
| IR | +0,035 | −0,011 | +0,022 | −0,036 |

En negrita, los que excluyen el cero. **Los dos miembros que promedian sobre la ventana (INR e IR)
cubren el cero en los cuatro y ningún otro lo hace**, aunque los intervalos de IR son tres a cinco
veces más angostos que los de INR (la anchura de INR es su propia distorsión: con la varianza de la
suma inflada 19 veces, su sesgo está mal determinado más que medido). Y el desplazamiento vive en
el par de amplitud, no en las cinéticas: R llega a −1,38 en la corriente unitaria y +0,82 en N_ch.

**Esto cierra el puente mecanismo→resultado.** En log10, R queda en +0,146 sobre N_ch e IR en
−0,0065, y la Figura 4 mide sobre estimaciones AJUSTADAS medianas de 0,099 y 0,110 para los dos que
no promedian contra 0,003 y 0,002 para sus contrapartes promediadas. El score predice, en orden, el
desplazamiento que se ve en las nubes. Y coincide con lo que la subsección de Figura 4 ya afirma:
que promediar es lo que saca el sesgo de amplitud.

**El panel sigue pendiente y ya no bloquea.** El número está en la prosa (beat 3) y el chunk que lo
imprime está en la notebook. Anotarlo en el panel D, como la fila A anota el logL, cuesta cero
espacio vertical y sigue siendo la manera de que se vea; queda como decisión de artefacto, no de
contenido.

---

## Los beats, con presupuesto

Hoy el bloque son **799 palabras de prosa** en siete párrafos, `04_results.tex:330-398`. Results
está en su piso (`08_length_plan.md`), así que esto es neutro: lo que entra sale de adentro.

**1. Encuadre. 85** (hoy 135). Qué mide y por qué una sola celda. La configuración de seis
parámetros de LSE se queda entera: es una trampa real y no está en ningún otro lado. La
justificación de mil contra diez mil pasa a cláusula.

**2. La escalera y el par que no se separa. 85** (hoy 151). Los 118,9 nats como contexto, una
oración. LSE e ILSE idénticos a tres cifras con la razón en media línea: con una sola varianza
constante el promediado no tiene sobre qué actuar. Es un resultado negativo y vale, pero no 151
palabras.

**3. EL CENTRO: el factorial en D y F. 190.** Las dos causas, las dos filas, las cuatro esquinas de
cada una, el interruptor que repara cada una. Cierra con la oración que hoy está dispersa: un
miembro que toma un eje y no el otro repara una falla y conserva la otra, así que la familia no se
ordena en una escala. Acá entra el sesgo acumulado si se anota el panel D.

**4. Significancia contra magnitud. 95.** Reemplaza a "Per interval every member reports its
information correctly", que la tabla de la propia figura contradice. Lo chico es la MAGNITUD
(mediana de log10(J_t/F_t) dentro de ±0,19, menos de un factor 1,6); en significancia la fila E
excluye el 1 en 95% de los pasos para LSE, 72,5% para R, 88,8% para MR, 100% para VR en N_ch. El
contraste queda 1,6 contra 20, más fuerte que "bien contra mal". Acá van las fracciones que hoy se
imprimen y no se reportan, **con n = 1.000 dicho al lado** (ver la nota sobre n al final) y con la
frase de que los pasos están correlacionados, así que son descriptivas y no un test.

**5. La factorización, y la cancelación de R. 120** (era 90; creció al medirse y se lo ganó). La
identidad, los dos factores de IR en uno por separado, y R llegando a 1,07 como producto de 0,72 y
1,48. Cierra con el caso puro: en los que no condicionan la falla es casi toda correlación, INR en
N_ch 1,01 contra 11,76. La predicción AR(1) no entra, murió al medirse.

**6. Milescu y la oferta al experimentalista. 90.** La ACF del score en lag 1 (IR −0,0043 ± 0,0043,
R 0,191, NR 0,864, LSE 0,856) es la correlación local que Milescu nombró, medida al nivel de la
inferencia y no del dato. Pegado y sin párrafo propio: la ACF del residuo es la misma memoria al
nivel del dato y la única cantidad de la figura computable sobre un registro real. Hoy está partido
entre este bloque y la leyenda.

**7. El washout. 120** (hoy 173). Con el calificador adelante y no atrás: es una propiedad de los
miembros que CONDICIONAN. N_ch y k_on al piso numérico en seis intervalos, k_off sobrevive por la
media y la corriente unitaria por la varianza; NR e INR retienen 15-16% de su información de pulso
sobre N_ch contra 0,5-0,8% de los recursivos. Y la evidencia barata que hoy no se usa: la máscara
de información deja el registro entero informativo para N_ch en NR e INR (denominador 80) y sólo el
pulso en los recursivos (43 en R, 46 en MR, 44 en VR, 46 en IR). Se lee de la tabla sin abrir la
figura.

**8. Dos cláusulas. 60.** El caveat de LSE en 35 (su Fisher gaussiano presume varianza constante y
el simulador la viola, así que el factor 15 mide esa ruptura y no una falla del mismo tipo que la
del miembro de lazo abierto), hoy 95. Y MR/VR en 25, confirmando en la moneda de esta figura lo que
Figura 2 ya dijo con elipses: acumulado en k_off, R 1,074, MR 1,188, VR 1,198, IR 1,080. No
interpolan, están del otro lado de R.

**Aritmética: 85+85+190+95+90+90+120+60 = 815** contra 799 de hoy. Los 16 salen del beat 2 si hay
que cerrar exacto.

---

## El título de la subsección

El actual, "Calibration holds interval by interval and fails only when the intervals accumulate",
afirma justo lo que el beat 4 corrige. Y el par interruptor-reparación ya es el título de la
subsección siguiente, la de Figura 4.

Recomendado: **"The per-interval error is small and systematic, and the accumulation multiplies it
twentyfold"**. Alternativa si se prefiere nombrar las dos causas: "Two failures with two causes,
the sample and the accumulation".

---

## La leyenda

`Figure_3_caption.md` se rehace, no se parchea: sigue en la versión de siete miembros sin ILSE y su
párrafo del factorial cita r̄²_std y J/F.

La leyenda del manuscrito (`04_results.tex:406`) está al día en roster pero es puramente
descriptiva. Con este guion absorbe, además de las instrucciones de lectura, las fracciones que
excluyen el nulo para B, D y E, que es la afirmación de significancia.

El párrafo de MR/VR de la leyenda vieja se acorta o se acota: hoy afirma que VR es peor que MR con
una historia mecanística, y la fila G da vuelta ese orden (MR 0,357 contra VR 0,256) aunque el
acumulado en N_ch (1,188 contra 1,427) y las elipses de Figura 2 lo sostengan.

---

## Lo que muere

1. La oración de que NR cambia el signo de la discrepancia y que un error por muestra no puede
   hacer eso, ~23 palabras. La subsume el beat 3 con los cuatro números.
2. Dos tercios del caveat de LSE, ~60.
3. La justificación de las mil grabaciones contra diez mil, ~35, a cláusula.
4. "The reading is interval by interval... two parameters to a page", ~30: es instrucción de
   lectura y pertenece a la leyenda del suplemento, que ya la tiene.

---

## Abierto, en orden de decisión

1. **¿Va la anotación del score acumulado en el panel D?** Es la única adición al artefacto y la
   que cierra el hueco. Barata.
2. **Rehacer el chequeo de la inflación** con Ē pesada por F y ACF completa, antes de que el beat 5
   exista.
3. **¿La información que NR e INR conservan después del lavado es la misma que ya tenían, contada
   de nuevo?** Si se sostiene, las filas C y F son un solo hecho y el beat 7 se pega al 3 en vez de
   vivir aparte.

---

## La nota sobre n, y cómo se reabre lo de 10.000

Postergado el 2026-08-12. Lo que se midió ese día, para no re-derivarlo:

- **Cómputo, no es problema.** La corrida de ocho miembros del 2026-08-05 escribió entre las
  23:13:00 y las 23:19:28. A 10× es del orden de una hora.
- **Disco, sí lo es.** Los dumps pesan 1,00 a 1,16 GB cada uno, o sea 10 a 11,6 GB por miembro a
  10.000 y ~90 GB los ocho, contra 36 GB libres en un disco al 98%.
- **Sale igual, miembro por miembro.** `seed = 20260722` está FIJA en
  `ops/local/figure_3_time.macroir:52`, así que corridas separadas siguen puntuando el ensemble
  idéntico, que es la premisa de la figura. Correr, digerir, borrar el CSV, siguiente: pico ~12 GB.
  `figure_3_digest.R` filtra con awk en streaming, así que la RAM no cambia. La mitad de cada dump
  es la copia duplicada que el digest tira con `segment_index == 0`: evitarla en el writer partiría
  los 90 GB al medio.
- **Un arreglo previo en la notebook.** `score_stats` hace `rowCumsum(Sb)` DENTRO del loop de 400
  bootstraps. Es innecesario: el cumsum es por fila, así que `rowCumsum(M[idx,])` es exactamente
  `rowCumsum(M)[idx,]`. Se calcula una vez afuera y adentro se indexa. Idéntico al bit, 400 veces
  menos trabajo, y a 10.000 grabaciones es la diferencia entre un knit y una tarde.
- **Lo que costaría, que no es cómputo.** Las bandas se angostan √10 y hay dos cantidades que
  dependen de eso. El ratio acumulado de IR es 1,080 con IC [0,991, 1,17], que cubre el 1 por poco;
  escalando la semi-amplitud por √10 queda cerca de [1,05, 1,11], que lo excluye. Igual R con 1,074
  y R en N_ch por abajo con 0,946. Y la fracción de pasos que excluye el nulo NO es invariante en
  n: mide detectabilidad, no tamaño, así que el 0,0375 y el 0,0652 de IR subirían, no porque IR
  empeore sino porque un sesgo de 0,02σ se vuelve resoluble. (El escalado por √10 es aritmética
  sobre el IC impreso, no una medición.)
- **Consecuencia editorial.** "Los dos recursivos están en uno" pasaría a "están en 1,07 y 1,08,
  significativamente conservadores en menos de un 10%, contra 13 a 21". El título de la subsección
  de Figura 2 y las frases que dicen que IR cubre el uno necesitarían el número al lado.

**Mientras se quede en 1.000:** las fracciones del beat 4 y de la leyenda se reportan con el n
dicho explícitamente, porque el estadístico depende de él, y la excepción de las mil grabaciones
sigue escrita en el encuadre (beat 1), en cláusula.

---

# APLICADO (2026-08-12). Entregado contra presupuesto

| beat | antes | presupuesto | ENTREGADO |
|---|---|---|---|
| 1 encuadre | 135 | 85 | 81 |
| 2 escalera + el par que no se separa | 151 | 85 | 119 |
| 3 el factorial en D y F, con el sesgo acumulado | — | 190 | 253 |
| 4 significancia contra magnitud | 136 | 95 | 112 |
| 5 la factorización y la cancelación de R | — | 120 | 126 |
| 6 Milescu y la oferta al experimentalista | 109 | 90 | 99 |
| 7 el washout | 173 | 120 | 182 |
| 8 dos cláusulas | 95 | 60 | 82 |
| | **799** | **845** | **1054** |

## ABIERTO: dos cómputos del mismo sesgo de primer orden que no coinciden (2026-08-12)

El marcador "predicted (truth + bias)" de Figura 2 es `truth + DIB`, donde DIB es el
`Probit_statistics_Likelihood_Distortion_Induced_Bias` que emite el programa, leído de `_battery_sim`
(`figure_2.Rmd:232-259`). En la celda de Figura 2 (N_ch = 100, interval 0,1, ruido 0,1 que es
S̃ = 0,01), en log10:

| | DIB, i | mío F⁻¹g, i | DIB, N_ch | mío, N_ch | nube (medido), i |
|---|---|---|---|---|---|
| IR | +0,0014 | +0,0021 | −0,0038 | −0,0065 | sin resolver |
| NR | −0,0886 | −0,0683 | +0,1033 | +0,0820 | −0,103 |
| R | −0,1812 | −0,1307 | +0,2259 | +0,1461 | −0,156 |
| MR | −0,2865 | −0,1780 | +0,3241 | +0,1900 | −0,228 |

**Mismos signos, misma estructura (todo en el par de amplitud, nada en las cinéticas, IR en cero),
magnitudes mías 25 a 40% más chicas.** Son dos cómputos del mismo objeto de primer orden sobre
baterías distintas: el DIB sale de `battery_sim` (10.000 simulaciones, otro commit de motor) y el mío
de los digests de Figura 3 (1.000 grabaciones), con F reconstruida de dm, dv e y_var. Candidatos a
explicar la brecha, ninguno verificado: el punto de evaluación (θ_sim contra θ_pool), una F distinta
(la del pooled contra el promedio de las por grabación), o un segundo término en la definición del
DIB. **Nada del cuerpo depende de cuál sea el correcto**: la prosa de Figura 2 cita el DIB, que es lo
que el marcador dibuja, y la de Figura 3 cita el suyo sin afirmar que sean el mismo número.
Corresponde cerrarlo leyendo la definición del DIB en el código.

Lo que sí se agregó a la subsección de Figura 2, que no decía nada del marcador: que es el
desplazamiento predicho por el score y la información sin ajustar nada, y que cae del mismo lado que
la nube y cerca (−0,287 contra −0,228 medido en MR, −0,181 contra −0,156 en R, −0,089 contra −0,103
en NR, +0,001 en IR donde la nube no resuelve desplazamiento). Es el único panel del paper donde el
mecanismo y el resultado están dibujados juntos.

---

## TERCERA REESCRITURA, 2026-08-12 (LA VIGENTE)

Luciano rechazó la segunda. **El eje no es la descomposición sino qué tipo de información da cada
variable**, y el bloque anticipa el análisis de Figura 4. Siete párrafos, 892 palabras.

| párrafo | palabras | qué dice |
|---|---|---|
| 1 encuadre | 82 | qué mide, y el plan: cercanía, desplazamiento, incertidumbre mal reportada |
| 2 la logL | 146 | mide CERCANÍA GLOBAL y nada más; sólo se conocen diferencias de distancia, no la distancia; por eso el resto puede testear donde la escalera sólo ordena; más la cláusula de LSE |
| 3 los tres objetos | 97 | residuo, gradiente y FIM: qué nulo tiene cada uno y cuál está disponible en un registro real |
| 4 el desplazamiento | 131 | media del score, por intervalo y acumulada con la matriz entera; predice el primer momento de fig:plane |
| 5 distorsión por muestra | 141 | r̄² por miembro, el intervalo del escalón (28,45), y LSE de 0,010 a 1,91 con media uno por construcción |
| 6 distorsión por correlación | 150 | INR con residuo calibrado y 20× de error; la factorización y la cancelación de R; Milescu; el residuo como oferta; y el puente explícito a fig:plane |
| 7 el FIM solo | 145 | dónde se mide cada parámetro; el lavado; el doble conteo |
| | **892** | |

**Fuera del cuerpo, por decisión de Luciano:** la descomposición de la logL de la familia no
recursiva. Es interesante y desvía del objetivo del paper. Quedó como nota invisible en
`04_results.tex`, en un bloque `%` inmediatamente después del párrafo de la logL, con la identidad,
los dos chequeos, el piso de 71,9 y los scripts. La versión larga sigue arriba en este guion.

Título, tercera versión: **"The uncertainty is misreported for two independent reasons, one within a
sample and one across them"**. Los dos anteriores quedan anotados en el `%` con la razón de su
retiro.

---

## SEGUNDA REESCRITURA, 2026-08-12 (superada el mismo día)

La de arriba quedó superada el mismo día. El bloque se rehízo sobre el eje nuevo, la descomposición
de la escalera, y con el orden de argumento y no de filas. **910 palabras en nueve párrafos, contra
1054 de la primera pasada y 799 del original.**

| párrafo | entregado | qué carga |
|---|---|---|
| 1 encuadre | 40 | qué mide y con qué |
| 2 por qué el logL es diagnóstico | 64 | E_p[log q] = −H(p) − KL, la escalera es una escalera de KL |
| 3 la descomposición | 136 | la identidad, los dos términos medidos, el piso de 71,9 que INR paga con marginales perfectos |
| 4 LSE | 98 | falla por 2,2; varianza por grabación (cv 0,478); su entrada no es el mismo objeto |
| 5 el intervalo del escalón | 70 | 12,1 de 14,8 nats en 1 de 100; r̄² 28,45 |
| 6 las dos fallas | 141 | b = F⁻¹E[Σs] y el ratio acumulado; el factorial; predice fig:plane |
| 7 E y la factorización | 111 | por intervalo casi bien; IR con los dos factores en uno, R cancelando |
| 8 la ACF y la oferta | 101 | Milescu al nivel de la inferencia; el residuo como lo computable |
| 9 el washout | 149 | resultado del experimento, no de los algoritmos; el doble conteo |
| | **910** | |

Lo que se cayó respecto de la primera pasada: el párrafo de la calificación del Fisher de mínimos
cuadrados (ya no hace falta porque el texto dejó de citar su ratio acumulado, y la reserva viaja en
una cláusula del párrafo 4), la cláusula de MR/VR (vive en la subsección de Figura 2 con sus áreas
de elipse), las fracciones que excluyen el nulo (viven en la leyenda) y el inventario de la
excepción de las mil grabaciones.

Título de la subsección, cambiado dos veces el mismo día y anotado en el `%`: quedó **"The
likelihood gap splits into what the marginals miss and what independence costs"**. Si alguna vez se
corta la descomposición, el título a restaurar es el anterior, que sigue escrito en el comentario.

**Siguen faltando 65 palabras contra el presupuesto de 845, y son 111 más que el original.** Es
mucho menos que las 255 de la primera pasada y ahora el bloque carga tres resultados que antes no
existían. Si hay que cerrarlas, el párrafo 9 es el que más da (149) y lo que sale sin perder nada es
el conteo de pasos informativos.

Dónde está el exceso y qué costaría sacarlo:

1. **Beat 7, +62.** Lleva tres cosas (el lavado, el mecanismo del condicionamiento, el conteo de
   pasos informativos) más dos citas de suplemento. El presupuesto de 120 estaba mal puesto. Lo que
   se puede sacar es el conteo de pasos informativos, ~45, que es evidencia nueva pero redundante
   con el 0,5-0,8% contra 15-16% que ya estaba.
2. **Beat 2, +34.** Ya perdió la oración del eje de ventana (queda como cláusula y anotada en el
   `%`). Lo que queda es la escalera y el par LSE/ILSE; bajar más es sacar la razón de por qué no se
   separan.
3. **Beat 4, +17 y beat 8, +22.** Prosa comprimible, ~30 entre los dos con esfuerzo.
4. **El candidato grande: beat 5 entero, 126.** Es el hallazgo más fuerte de la reescritura (R
   llega a uno cancelando dos errores) y sacarlo devuelve la sección a 875. No lo haría.

Mi orden si hay que cerrar: 1 primero (el conteo de pasos, que se puede mandar a la leyenda del
suplemento donde vive el dato), después 3. Con eso quedan ~925.

**Qué se tocó, además de la prosa:**

- `04_results.tex`: título de la subsección cambiado, con la nota de por qué; los ocho beats
  reescritos; la leyenda de `fig:time` absorbe las fracciones que excluyen el nulo y la
  factorización en la glosa de la fila F.
- `figures/paper_both/figure_3.Rmd`: bloque nuevo al final del chunk `caption-numbers` que calcula e
  imprime z_T, Ē e inflación con bootstrap y semilla propia (20260812). **Falta el re-knit**; el
  `.html` en disco es el del 2026-08-05. Los estimadores puntuales son deterministas, los bordes de
  los IC se van a mover en el último dígito.
- `Figure_3_caption.md`: reescrito entero sobre este guion, con la nota de qué cambió.
- Verificado de paso, contra los digests: r̄² por miembro (LSE e ILSE exactamente 1,0000, NR 1,1526,
  INR 0,9995, R 0,7816, MR 0,6904, VR 0,8118, IR 0,9989), la ACF del residuo en lag 1 (LSE e ILSE
  0,8457, NR 0,5100, R 0,1181, IR −0,0103, o sea el score lee 1,6 a 1,7 veces el residuo), y los
  logL con su error (LSE e ILSE −262,30 ± 0,73, IR −143,45 ± 0,24, span 118,85).
- Una inconsistencia resuelta al pasar: el área de elipse de VR figuraba como 2,19 en el caption
  viejo y 2,18 en `04_results.tex:252`. Se adoptó 2,18, el del manuscrito.

**Abierto, y anotado en el `%` del propio párrafo:** el r̄² de 1,0000 exacto de los dos brazos de
mínimos cuadrados sobre 80.000 residuos parece ser el nivel de ruido marginalizado fijando el
residuo estandarizado en uno por construcción. Está inferido del número, no leído de la definición
del miembro, así que no se afirma en el cuerpo.
