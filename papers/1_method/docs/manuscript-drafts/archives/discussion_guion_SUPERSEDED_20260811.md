# Guion de la Discussion, reconstruido desde el archivo y reescrito (2026-08-11)

**El trabajo de la Discussion:** hacer las cuatro cosas que ninguna otra sección hace. (1) Contestar
las tres objeciones predecibles del editor y del referee: en qué se diferencia de Münch, que eLife ya
publicó; si filtrar no es demasiado caro (Del Core); y por qué creerle a un estudio simulado de dos
estados. (2) Convertir lo medido en instrucciones que un experimentador pueda ejecutar. (3) Declarar
el perímetro, incluida la desviación del propio miembro calibrado. (4) Decir la apuesta: lo que una
aproximación erra sobre la información pasa al mecanismo que elige.

Lo que NO le toca: repetir Results, y volver a decir lo que Theory ya imprimió.

Presupuesto declarado en el encabezado del archivo: ~1500. Hoy: **3670**.

Los conteos de este documento son de prosa sola (comentarios `%` y `% src:` descontados), medidos el
2026-08-11 sobre `sections/05_discussion.tex`. Suman 3663 más 7 de los dos títulos.

---

## Como está hoy: veintitrés beats en tres bloques

Los bloques son del reordenamiento del 2026-08-11 (permutación pura, ninguna línea cambió).

### Bloque 1, lo que se midió y qué compra: 874 palabras

| beat | línea | palabras | qué dice |
|---|---|---|---|
| B1 | 49 | 30 | veredicto: calibrado sobre casi todo el plano, la desviación medida donde no |
| B2 | 55 | 39 | usalo en todas partes; el rincón de pocos canales; la reparación es el sandwich |
| B3 | 72 | 78 | qué lo compra: condicionar en los dos extremos; un extremo solo cuesta |
| B4 | 87-88 | 137 | la predicción VR registrada de antemano, más la calificación obligatoria |
| B5 | 90 | 62 | el diagnóstico de Milescu, ahora medido |
| B6 | 92 | 237 | el cociente gobierna todo; τ_int lee el régimen desde una grabación; el 2217 |
| B7 | 100, 108 | 150 | cómo gastar una grabación (washout), y la pregunta que le contesta a Del Core |
| B8 | 110 | 141 | el ordenamiento de los bordes del mapa; mapa de conceptos; concesión a LSE |

### Bloque 2, contra qué se compara: 1177 palabras

| beat | línea | palabras | qué dice |
|---|---|---|---|
| B9 | 125-127 | 112 | el dispositivo Kalman concedido; qué es específico de acá; el linaje |
| B10 | 138-139 | 347 | Münch: las tres separaciones, la ubicación, la generalización de la emisión |
| B11 | 183, 197 | 181 | Del Core y el costo, con la única cifra publicada; dónde no se justifica |
| B12 | 224 | 161 | el menú de reparaciones sin cambiar de verosimilitud |
| B13 | 226 | 87 | dónde la ruta clásica es fuerte, y hay que decirlo |
| B14 | 228-229 | 101 | la apuesta: el flip state tenía firma, el orden asimétrico no |
| B15 | 235, 245 | 188 | la matriz completa; volumen y tamaño muestral efectivo; qué no se dirime |

### Bloque 3, lo que este paper no dirime: 1612 palabras

| beat | línea | palabras | qué dice |
|---|---|---|---|
| B16 | 272-280 | 138 | preámbulo y los cuatro ítems que dependen del experimento del lector |
| B17 | 282-292 | 162 | el quinto: lo que el diagnóstico no puede ver, especificación compartida |
| B18 | 294-296 | 47 | qué entrega el paper |
| B19 | 298-300 | 198 | la desviación de dos lados del propio miembro, 0.645 a 1.701 |
| B20 | 322-323 | 185 | el modelo de ruido: gaussiano y blanco, no corrido |
| B21 | 343 | 370 | el kernel de Bessel, el bound, la forma operativa, el costo |
| B22 | 406-417 | 184 | el origen de la desviación no se conoce; cinco candidatos refutados; conjetura |
| B23 | 429-430 | 328 | el perímetro, la regla 4Np(1-p), los compañeros, la biblioteca |

### Los tres hechos que salen de la cuenta

1. **El 44% de la sección son limitaciones** (1612 de 3670), y otro 12% son concesiones dentro del
   bloque 2 (B11 "the cost is real", B13, la concesión de B9). El encabezado del archivo ya lo dice
   con otras palabras: "about 2270 words of caveat, limitation and concession before a single
   sentence of argument".
2. **La tesis del paper está en la posición 14**, a unas 2300 palabras del comienzo. En
   `abstract_guion.md` es el beat 1. Es la misma frase (`05_discussion.tex:228-229` es la forma
   comprimida de `01_introduction.tex:254-263`).
3. **El párrafo más grande de la sección no es Münch: es Bessel**, 370 contra 347, y está dentro del
   bloque de limitaciones.

### Dos hallazgos verificados hoy, ambos afectan la aritmética

**(a) B21 es la segunda vez que se cuenta el bound.** `02_theory_full.tex:94` ya imprime el grupo
adimensional, la razón de la ganancia unitaria en cero, el 3% a cinco veces la inversa del corte, el
primer orden para lo blanco contra segundo orden para lo que trae tiempo de correlación, el promedio
en bloques como remedio, el intervalo del salto como excepción y "left to later work". Lo único
propio de B21 son dos cosas: la cifra publicada del costo (más de diez veces para la comparación de
nueve esquemas) y el argumento de que la condición es más difícil justo donde el miembro calibrado
más conviene, y eso último ya está también en B23. Los tres sellos de `[BESSEL-BOUND] CLOSED` de este
archivo (líneas 324, 358, 462) están pegados con la puntuación rota, y las notas largas debajo de
ellos fueron escritas cuando la decisión estaba abierta.

**(b) Los punteros `% CUT` del archivo están corridos.** El pase de longitud del 2026-08-10 justificó
cada corte diciendo dónde sobrevive el contenido, y esos números ya no apuntan ahí. Verificado hoy,
viejo → real: 411-417 → **442** (el censo 93/70/31%); 193-198 → **218-224** (0.914 contra
0.947-0.949); 294-297 → **319-321** (la ACF del score); 307-311 → **330-331** (el washout);
639-643 → **758** (una a tres décadas); 405-407 → **424-425** (la generalidad más allá de los
canales); los pisos de alcanzabilidad → **764-766**. `abstract_guion.md` también trae el viejo
(433-435 por 442). Si el guion nuevo se ejecuta, los punteros hay que rehacerlos una sola vez y con
fecha, no repetirlos.

---

## Como debería ser: ocho beats

### 1. EL VERDICTO. 60 palabras
Fusión de B1 y B2. Calibrado sobre casi todo el plano medido, la desviación medida donde no, el
rincón de pocos canales nombrado, y la reparación es el sandwich que el diagnóstico ya entrega.
**Para el editor:** la sección abre con un resultado y no con un descargo.

### 2. QUÉ LO COMPRA. 90 palabras
B3 completo. Condicionar en los dos extremos; condicionar en el inicio solo no lo compra y en la
celda de referencia cuesta; toda la distancia está en el update. VR entra como una cláusula que
remite al suplemento, no como párrafo.
**Para el referee metodológico:** es el mecanismo, y es lo que hace que el resultado no sea una tabla.

### 3. LA APUESTA. 170 palabras
B14 (90) más lo defendible de B15 (80), subido de la posición 14 a la 3. El flip state dejaba una
firma que se señalaba en el registro; la activación asimétrica de las mismas grabaciones no deja
ninguna y descansa en razones de evidencia, así que lo que una aproximación erra sobre la información
pasa al mecanismo que elige. Cerca del máximo la misma matriz entra en la evidencia por el término de
volumen y el reescalado del tamaño muestral efectivo, derivados en otro lugar. Los tres niveles
separados, con las palabras de la corrección del 2026-08-10: sesgo de parámetros, calibración de la
información y validez de la evidencia son afirmaciones distintas, y la calibración es necesaria para
la tercera sin ser suficiente. Cierra con la frase que justifica el diseño simulado: ninguna cantidad
de datos grabados identifica el modelo generador.
**Para el editor:** es el beat que decide entre valuable e important, y hoy llega demasiado tarde para
decidir nada.

### 4. QUÉ PUEDE HACER UN EXPERIMENTADOR. 280 palabras
Comprime B6 + B7 + B8 (528) y sube la frase de la biblioteca de B23. Cuatro ítems, todos ya medidos:
- el cociente entre ruido instrumental y varianza de gating gobierna a toda la familia, así que la
  autocorrelación integrada del residuo estandarizado lee el régimen desde una sola grabación, sin
  ensemble y sin verdad conocida (90). **La advertencia va pegada y no puede faltar:** lee régimen y
  no certifica, `NR` erra hasta 2217× con residuo indistinguible de blanco (`04_results.tex:595`).
- cómo gastar una grabación: después de sacar el agonista no hay más información sobre N_ch ni k_on,
  y sí sobre k_off y la corriente unitaria (70).
- el ordenamiento de los bordes, con las dos cláusulas de alcance en una línea cada una: es un mapa
  de conceptos y no un diagrama de fases, y el brazo clásico es un ajuste de cuadrados mínimos y no
  análisis de fluctuaciones no estacionario (70).
- la biblioteca: corré el diagnóstico en tu propio esquema y en tu propio punto (50).

**Para el editor:** es la respuesta a "¿y qué hago con esto?", y es donde la biblioteca hace doble
trabajo, porque es también la respuesta a la objeción de dos estados.
**PENDIENTE, y cambia la frase:** hoy B23 dice que se pueden correr las dos identidades. Según
`abstract_guion.md` (2026-08-11, macroir `911828d`) el lazo ya existe, `mi.distortion(...)` en Python
y `macroir_distortion(...)` en R, y devuelve la matriz de distorsión, sus autovalores y el sesgo de
primer orden. La guarda tiene que viajar con la frase: la biblioteca implementa el miembro calibrado
y nada más, así que el lector audita SU punto con ESE miembro, y comparar peldaños necesita R y NR en
runtime.

### 5. CONTRA QUÉ. 310 palabras
B10 a 250, con las tres separaciones intactas (qué está dentro del veredicto, qué devuelve cada una,
qué cuesta cada una), la ubicación con el miembro recursivo instantáneo y la generalización de la
emisión que les fuerza su propia clausura. B9 dobla adentro en 60: el dispositivo está concedido, lo
específico es la realización en tiempo continuo y el mapa.
**Para el editor:** es el ítem más pedido de todas las lecturas externas, y el archivo tiene razón en
su propia nota: acortar una comparación es cómo una comparación se convierte en un desprecio. 250 es
el piso, no el objetivo.

### 6. EL COSTO Y LA RUTA CLÁSICA. 200 palabras
B11 (140) más B13 (60). El costo es real y ahora está cuantificado en unidades de incertidumbre mal
reportada; la única cifra es la publicada, doce a cincuenta y dos días por esquema, ofrecida como
escala y no como benchmark; dónde el ruido instrumental es grande los miembros baratos recuperan y
ahí el filtro no se justifica por calibración; y donde hay sweeps intercambiables la ruta clásica es
fuerte y hay que decirlo. El menú de reparaciones (B12) sale a Appendix 1.
**Para el referee:** las dos concesiones son las que compran credibilidad, y las dos son citables sin
que el paper pierda nada.

### 7. LO QUE ESTE PAPER NO DIRIME. 620 palabras, cuatro ítems
El preámbulo de los cuatro ítems del lector (B16, 138) se disuelve: tres de los cuatro ya se dicen
donde corresponde (costo y bootstrap en el beat 6, precisión en el 4), así que como enumeración es la
tercera vez.
- **el propio miembro calibrado** (200): la desviación de dos lados, 0.645 a 1.701, con el par de
  celdas; que no se vuelve segura por errar conservadoramente en promedio; dónde caen las dos
  excursiones y hasta dónde llega esa expectativa; que el origen no se conoce, que los cinco
  candidatos eran todos criterios sobre cumulantes por muestra y por eso ninguno lo localiza, el
  exponente −0.25 contra −0.5 y −1.0, y la conjetura en las palabras de Luciano. Fusión de B19 y B22,
  que hoy son 382 y dicen dos mitades de lo mismo.
- **el modelo de la medición** (180): una sola desviación vista dos veces, que es lo que la nota de la
  línea 327 del propio archivo ya dice. El ruido instrumental se supone gaussiano y blanco, las
  grabaciones se apartan de eso, no lo corrimos. El kernel uniforme es la primera aproximación al del
  filtro, y lo que se pierde está acotado y la condición se chequea de antemano con dos números que el
  experimentador eligió, remitiendo a Theory en vez de volver a derivarlo. Las dos cosas propias de
  B21 se quedan: la cifra del costo y que la condición es más difícil donde el miembro calibrado más
  conviene. Fusión de B20 y B21, que hoy son 555.
- **lo que el diagnóstico no puede ver** (120): dos implementaciones independientes, así que todo lo
  que está aguas abajo del punto donde se separan está testeado y la especificación que comparten
  aguas arriba no. Con el remedio general nombrado en una línea, no en cuatro.
- **el perímetro** (150): dos estados, p = 0.5, un solo salto, sólo verosimilitud; la regla
  4Np(1-p) con su residuo honesto (la asimetría binomial va como (1-2p) y nada acá la toca); los
  compañeros, con el kernel filtrado motivado en una cláusula en vez de igualado a los otros cuatro.

**Para el editor:** el título del bloque anuncia lo que trae, que es más visible que interleavear
caveats, y ninguno desaparece.

### 8. CIERRE. 40 palabras
La última oración de B18, que es la mejor que tiene la sección: el lector aporta su número de sweeps,
su rundown, su tamaño de modelo y su tolerancia y decide; lo que no puede medir por su cuenta, y es
lo que se entrega acá, es si el intervalo reportado es el intervalo que se entrega.

---

## Aritmética: lo planeado y lo entregado

**SUPERADO por el guion de 1500 al final de este documento (2026-08-11).** Esta tabla es el registro
de lo que se ejecutó y de por qué el presupuesto de 1830 estaba mal puesto; el plan vigente es el de
abajo.

Corregido el 2026-08-11 después de escribir la sección. La tabla anterior tenía dos defectos: B5
(Milescu, 62) no estaba asignado a ningún beat, y la frase de la biblioteca se contaba dos veces. El
beat 2 se queda con B5, que es el mecanismo y es lo que hace funcionar al beat 4.

| beat nuevo | de dónde | antes | presupuesto | ENTREGADO |
|---|---|---|---|---|
| 1 veredicto | B1+B2 | 69 | 60 | 69 |
| 2 qué lo compra, con Milescu | B3+B5+cláusula de B4 | 277 | 150 | 214 |
| 3 la apuesta | B14+B15 | 289 | 170 | 244 |
| 4 qué puede hacer | B6+B7+B8+biblioteca | 611 | 280 | 552 |
| 5 contra qué | B9+B10 | 459 | 310 | 447 |
| 6 costo y ruta clásica | B11+B13 (B12 al apéndice) | 429 | 200 | 367 |
| 7 no se dirime | B16+B17+B19..B23 | 1482 | 620 | 1174 |
| 8 cierre | B18 | 47 | 40 | 47 |
| | | **3670** | **1830** | **3114** |

El menú de reparaciones aterrizó: `sections/10_appendix_repairs.tex`, textual, con las tres cifras sin
verificar que estaban en su bloque de comentarios, cargado último en `elife_paper.tex` para que 08
siga siendo Appendix 1 y 09 Appendix 2 y no se renumere nada, y citado desde el beat 6 con una
cláusula de 31 palabras (que es la diferencia entre 3083 y 3114).

**El presupuesto no se cumplió y la razón no es de redacción: estaba mal puesto.** Cada beat entregado
está entre un 20% y un 90% por encima de lo presupuestado, y el patrón es sistemático, así que el
error está en la estimación y no en la ejecución. Dos casos lo muestran:

- **El beat 7 se presupuestó en 620 sumando cuatro ítems que suman 650**, y de esos cuatro tres son
  alcance declarado o límite medido donde cada cláusula dice una cosa distinta. El perímetro solo son
  245 palabras de las cuales no sobra ninguna: esquema, protocolo, datos simulados con su razón,
  inferencia sólo verosimilitud, el aislamiento de la clausura, la ambigüedad canales/aperturas, la
  regla 4Np(1-p) con su residuo, los cinco compañeros y el kernel filtrado motivado. Presupuestarlo en
  150 era decir que se puede declarar el mismo alcance en 150, y no se puede.
- **Münch quedó en 344 contra 250**, y eso fue deliberado: perdió once palabras de tejido conectivo y
  ninguna cláusula. La instrucción del archivo viejo es correcta y sigue en pie.

**Lo entregado: 3083, un 16% menos que el predecesor, con la estructura del guion completa y sin
tocar un solo caveat.** De donde salieron las 587: el menú al apéndice (−161), la fusión de Bessel con
el ruido y la remisión a Theory (−152), VR a una cláusula (−63), la enumeración de los cuatro ítems
(−83, de los cuales 68 volvieron al beat 6 donde se argumentan), la apuesta comprimida (−45), el
párrafo del cociente (−43) y ajustes menores.

**Para llegar a 1500 faltan 1583 palabras y no hay forma de conseguirlas sin decidir qué caveat
pierde.** Las cuatro opciones, con precio y con lo que cuesta cada una:

1. **Münch a 150** (−194). Precio: acortar una comparación es cómo una comparación se convierte en un
   desprecio, y es el ítem que toda lectura externa predijo que el editor iba a pedir.
2. **El bloque "no se dirime" a la mitad** (−590). Es la opción con más palabras y la única que las
   tiene. Precio: el perímetro, la desviación de dos lados del propio miembro, la conjetura o el
   modelo de la medición. Uno de esos cuatro tiene que salir entero, porque comprimirlos ya se hizo.
3. **El beat 4 a la mitad** (−276). Precio: es la única prosa accionable de la sección, y el washout
   se lleva la única cita de los suplementos 1 y 2 de la Figura 3 en todo el manuscrito, que habría
   que recolocar en Results.
4. **La apuesta y el puente a la evidencia a 120** (−124). Precio: es la tesis del paper y el beat 1
   del abstract. No.

Mi lectura: la sección no baja de ~2700 con todos los caveats puestos, y el 1500 era un presupuesto
heredado que nunca se confrontó con el contenido. La decisión que hay que tomar no es "cortar más",
es si el objetivo sigue siendo 1500.

## Las tres decisiones que no son de redacción

1. **El párrafo de VR** (B4, 137 palabras) sobre una predicción registrada de un miembro demotado al
   suplemento. Propuesta: una cláusula en el beat 2. Lo que se pierde es la demostración de que el
   criterio localiza un defecto derivado del álgebra antes de correr el miembro, que es un argumento
   sobre el diagnóstico y no sobre VR, y por eso puede vivir en una cláusula.
2. **El menú de reparaciones** (B12, 161). El encabezado del archivo ya lo nomina como el corte menos
   dañino, a Appendix 1, y observa que es un menú y no una concesión. Nada de eso se dice en otra
   parte del manuscrito, así que mover no es lo mismo que borrar.
3. **La enumeración de los cuatro ítems del lector** (B16, 138). Se disuelve en los beats 4 y 6.

## Deliberadamente afuera del guion nuevo, y de dónde salen si alguien los pide

- Los números del censo y del coverage, que ya se cortaron el 2026-08-10 y viven en
  `04_results.tex:442` y `:218-224`.
- La derivación del bound de Bessel: `02_theory_full.tex:94`, y acá sólo la condición y la remisión.
- La ACF del score por miembro: `04_results.tex:319-321`, y en `03_diagnostics.tex:5` la relación
  1.6-1.7× entre residuo y score.
- Los pisos de alcanzabilidad, 17 y 52 canales: `04_results.tex:764-766`.
- La generalidad más allá de los canales: `04_results.tex:424-425`. Si alguna vez se quiere como
  prosa, va en la Introduction.
- Los ratios de costo por evaluación (1.4× a 5×): no se re-derivaron de ningún artefacto primario de
  este repo, y el archivo tiene razón en dejarlos sin imprimir.

---

# EL GUION DE 1500 (2026-08-11). Qué cabeza rueda

## El criterio de la guillotina

Los dos pases anteriores fallaron por la misma razón: trataron "este caveat no puede salir del paper"
como si dijera "este caveat no puede salir de la Discussion". No es lo mismo, y ahí están las mil
seiscientas palabras. Tres reglas, en orden de aplicación:

1. **Segunda narración, muere.** Si el contenido está impreso en Theory, Results o Methods, en la
   Discussion queda el puntero y nada más. No es una concesión al lector, es que ya lo leyó.
2. **Caveat con casa propia, se muda.** Una limitación sobre la especificación del modelo pertenece a
   Methods, y una sobre el kernel de medición a Theory. Sale de la sección y **sigue en el paper**.
   Cada mudanza es un transplante que hay que ejecutar, y si no se ejecuta el paper pierde el
   contenido: están numeradas abajo.
3. **Lo que queda es lo que sólo la Discussion puede decir:** el veredicto, el mecanismo, la apuesta,
   la instrucción para el experimentador, la comparación con otros grupos, y el límite del propio
   miembro. Eso no se toca.

**Una constricción que ya no existe, verificada hoy.** El archivo viejo decía que el párrafo del
washout tenía la única cita de los suplementos 1 y 2 de la Figura 3 del manuscrito, así que no se
podía adelgazar sin recolocarla. Falso hoy: `04_results.tex:334` los cita, y `:570`, `:580`, `:593`
citan los de la Figura 5. Ningún párrafo de la Discussion es el único portador de una cita de
suplemento, así que todos se pueden cortar por contenido y no por bibliografía.

## La tabla, bloque por bloque

| bloque | entregado | nuevo | veredicto y qué se pierde |
|---|---|---|---|
| veredicto | 69 | 55 | ADELGAZA. Las dos oraciones en una |
| qué lo compra | 152 | 100 | ADELGAZA. VR de 74 a 25 palabras, con la calificación puesta |
| Milescu | 62 | 40 | ADELGAZA, y se pega al anterior: es el mismo mecanismo |
| la apuesta | 244 | 150 | ADELGAZA. Mueren las fórmulas de volumen y α\*, que están en la Introduction, y el recap de "cuál de los dos órdenes"; queda el puntero |
| el cociente | 194 | 100 | ADELGAZA. Muere el inventario de figuras; queda el cociente, τ_int y el 2217 |
| washout | 137 | 60 | ADELGAZA fuerte. Un resultado y una instrucción |
| el mapa | 137 | 75 | ADELGAZA. Muere la concesión del sesgo de LSE, textual en `04_results.tex:387-390` |
| la biblioteca | 84 | 40 | ADELGAZA |
| arte previo | 103 | 60 | ADELGAZA. Queda el dispositivo concedido, la realización y el linaje |
| **Münch** | 344 | 190 | **RUEDA LA MITAD.** Muere el recap de su método (cómo cuentan, la masa de 0 a 0.95, que lo propusieron como método general). Viven las TRES separaciones y la cláusula de ubicación |
| Del Core | 181 | 115 | ADELGAZA. Muere la oración de los grados de libertad (Results la tiene) |
| ruta clásica | 186 | 105 | ADELGAZA. Los tres ítems del lector a una cláusula; queda el puntero al apéndice |
| el propio miembro | 362 | 185 | ADELGAZA. Muere la explicación de por qué no se nombra un umbral en canales; queda la cláusula en el perímetro |
| **modelo de la medición** | 403 | 55 | **RUEDA LA CABEZA.** La mitad del ruido: Methods ya tiene la premisa (`06_methods.tex:41`, "The emission model is white instrumental noise only", con pink y proporcional en cero), así que acá queda sólo la consecuencia, que es lo único que Methods no dice. La mitad del kernel: está impresa en `02_theory_full.tex:94`, entera. Transplante T1 |
| **lo que el diagnóstico no puede ver** | 157 | 0 | **RUEDA ENTERA.** Se muda a Methods, al lado de la especificación de la que habla (`06_methods.tex:349` ya está en ese vecindario). En la Discussion queda una cláusula dentro del cierre. Transplante T2 |
| el perímetro | 245 | 115 | ADELGAZA a la mitad. La enumeración del alcance se muda a Methods, que ya tiene dos de los cuatro ítems (`:9` dos estados, `:64` un solo salto). Quedan la regla 4Np(1-p) con su residuo honesto y los compañeros con la razón del kernel filtrado. Transplante T3 |
| el cierre | 47 | 55 | CRECE 8. Absorbe la cláusula del diagnóstico |
| | **3114** | **1500** | |

## Los tres transplantes. Si no se ejecutan, el paper pierde contenido

**T1. La cifra del costo del kernel, a Theory.** `02_theory_full.tex:94` cierra con "Carrying a
general kernel through the same construction, with the coloured noise it brings, is left to later
work." Ahí va, en la misma oración: y se estimó que convolver el kernel del amplificador
explícitamente, a 10 kHz de Bessel con digitalización a 50 kHz, sube más de diez veces la demanda de
una comparación bayesiana de nueve esquemas (`moffatt2025bayesian`). **La guarda viaja con la cifra:**
lo estimado es la COMPARACIÓN DE NUEVE ESQUEMAS, no una evaluación de verosimilitud. Es la única cosa
del párrafo de Bessel que no está ya en Theory, junto con que el tamaño se calcula de antemano y no
depende del rig, que se queda en la Discussion como cláusula.

**T2. El límite del diagnóstico, a Methods.** Las dos implementaciones independientes, qué queda aguas
abajo del punto donde se separan, la especificación que comparten aguas arriba, y el remedio general
(verificar cada objeto contra los teoremas que debe satisfacer). ~120 palabras. Va donde Methods ya
dice que una especificación cuyos símbolos no se encuentran en el código no es una especificación.

**T3. La enumeración del alcance, a Methods.** Faltan dos de los cuatro ítems: que los datos son
simulados porque el diagnóstico necesita una verdad conocida, y que la inferencia es sólo
verosimilitud, sin prior ni cómputo de evidencia. ~35 palabras. Y **p = 0.5 tiene que quedar dicho en
Methods explícitamente**: hoy sale de los parámetros (k_on[A] = k_off = 100 s⁻¹ a 10 μM) y no está
escrito, y la Discussion lo necesita declarado para poder ofrecer la regla de transferencia.

## Las seis cabezas que ruedan del paper entero, no sólo de la sección

Estas no se mudan a ningún lado. Es lo que hay que firmar.

1. **El recap del método de Münch**, ~150 palabras. La comparación queda más seca: entra directo a las
   tres separaciones sin explicar antes cómo cuentan ellos. **Es la más discutible de las seis**, y la
   nota del archivo viejo ("acortar una comparación es cómo una comparación se convierte en un
   desprecio") apunta justamente acá. Contra-argumento: lo que se corta es el RESUMEN de su método, no
   ninguna de las separaciones, y un editor de eLife que publicó ese paper no necesita que le
   expliquen cómo cuentan.
2. **Las fórmulas del volumen y α\*** en la apuesta, ~45. La Introduction las dice en prosa ("a volume
   term and a rescaling of the effective sample size"). Se pierde la forma explícita.
3. **La concesión del sesgo de LSE** en el párrafo del mapa, ~30. Textual en `04_results.tex:387-390`,
   donde además está marcada como "worth making plainly".
4. **La oración de los grados de libertad** en Del Core, ~30. Results la tiene con número
   (`:552`, τ_int y las observaciones independientes efectivas).
5. **La explicación del umbral en canales**, ~60. Se pierde el razonamiento (que un umbral en canales
   estaría en la unidad equivocada); queda la ambigüedad declarada en el perímetro y los pisos de
   alcanzabilidad en Results.
6. **El recap de "cuál de los dos órdenes es el correcto"**, ~55 de sus 121. **NO puede ir a cero:**
   la nota del 2026-08-07 dice que es la oración que mantiene honesta la cita de P2X2, porque el
   abstract y la Introduction dicen que el orden se movió con la aproximación. Queda un puntero de
   doce palabras que dice que hace falta recalcular las evidencias y que acá no se hace.

## Qué NO se toca, y por qué

- **Las tres separaciones de Münch**, la cláusula de ubicación y la generalización de la emisión.
- **La apuesta** (flip state contra activación asimétrica) y la oración de que ninguna cantidad de
  datos grabados identifica el modelo generador, que es lo que hace del diseño simulado un requisito.
- **La desviación de dos lados del propio miembro** con sus dos números, el bound de no-normalidad,
  los cinco candidatos refutados y la conjetura en las palabras de Luciano.
- **El 2217**, que es la advertencia sobre la lectura barata.
- **La regla 4Np(1-p) con su residuo** (la asimetría binomial va como (1-2p)).
- **Las dos concesiones**: dónde los miembros baratos recuperan y dónde la ruta clásica es fuerte.

## Aritmética

55 + 100 + 40 + 150 + 100 + 60 + 75 + 40 + 60 + 190 + 115 + 105 + 185 + 55 + 0 + 115 + 55 = **1500**,
más los dos títulos. El margen es cero: cualquier cosa que se recupere sale de otra.
