# Guion de la Discussion. Plan vigente: 1500 palabras (2026-08-11)

**El guion viejo está archivado, no borrado:** `archives/discussion_guion_SUPERSEDED_20260811.md`,
393 líneas. Ahí viven, íntegros, (a) la reconstrucción de los veintitrés beats del predecesor de 3670
palabras con sus líneas verificadas, (b) el plan de ocho beats con presupuesto de 1830 y sus notas
"para el editor" beat por beat, y (c) la autopsia de por qué ese presupuesto estaba mal puesto, con la
tabla de planeado contra entregado. Nada de eso se repite acá. Lo que sigue vivo de ese documento y no
se puede perder de vista está apuntado en su lugar; el resto es historia y se lee allá.

La sección misma, en su versión de 3670 palabras, está en
`archives/05_discussion_SUPERSEDED_20260811.tex`.

---

## El trabajo de la Discussion

Hacer las cuatro cosas que ninguna otra sección hace.

1. Contestar las tres objeciones predecibles del editor y del referee: en qué se diferencia de Münch,
   que eLife ya publicó; si filtrar no es demasiado caro (Del Core); y por qué creerle a un estudio
   simulado de dos estados.
2. Convertir lo medido en instrucciones que un experimentador pueda ejecutar.
3. Declarar el perímetro, incluida la desviación del propio miembro calibrado.
4. Decir la apuesta: lo que una aproximación erra sobre la información pasa al mecanismo que elige.

Lo que NO le toca: repetir Results, y volver a decir lo que Theory ya imprimió.

---

## Estado entregado hoy: 3114 palabras, veintiún bloques

Prosa sola, comentarios `%` y `% src:` descontados, medido sobre `sections/05_discussion.tex` después
de la reescritura del 2026-08-11.

| línea | palabras | bloque |
|---|---|---|
| 73 | 69 | el veredicto |
| 75 | 152 | qué lo compra, con VR en cláusula |
| 76 | 62 | Milescu, ahora medido |
| 78 | 244 | la apuesta, y el puente a la evidencia |
| 80 | 194 | el cociente, τ_int desde una grabación, el 2217 |
| 82 | 52 | washout, el resultado |
| 83 | 85 | washout, la instrucción para el banco |
| 85 | 137 | el ordenamiento de los bordes del mapa |
| 86 | 84 | la biblioteca |
| 98 | 30 | el dispositivo Kalman concedido |
| 99 | 73 | qué es específico de acá, y el linaje |
| 101 | 344 | Münch |
| 131 | 181 | Del Core y el costo |
| 140 | 186 | la ruta clásica, con el puntero al apéndice |
| 144 | 79 | el propio miembro no está exento |
| 145 | 215 | dónde caen las excursiones; los cinco candidatos |
| 146 | 68 | la conjetura |
| 157 | 403 | el modelo de la medición, ruido y kernel fusionados |
| 184 | 157 | lo que el diagnóstico no puede ver |
| 186 | 245 | el perímetro |
| 197 | 47 | el cierre |
| | **3107** | más 7 de los dos títulos = **3114** |

---

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

---

## La tabla, bloque por bloque

| bloque | entregado | nuevo | veredicto y qué se pierde |
|---|---|---|---|
| veredicto | 69 | 55 | ADELGAZA. Las dos oraciones en una |
| qué lo compra | 152 | 100 | ADELGAZA. VR de 74 a 25 palabras, con la calificación puesta |
| Milescu | 62 | 40 | ADELGAZA, y se pega al anterior: es el mismo mecanismo |
| la apuesta | 244 | 150 | ADELGAZA. Mueren las fórmulas de volumen y α\*, que están en la Introduction, y el recap de "cuál de los dos órdenes"; queda el puntero |
| el cociente | 194 | 100 | ADELGAZA. Muere el inventario de figuras; queda el cociente, τ_int y el 2217 |
| washout (los dos) | 137 | 60 | ADELGAZA fuerte. Un resultado y una instrucción |
| el mapa | 137 | 75 | ADELGAZA. Muere la concesión del sesgo de LSE, textual en `04_results.tex:387-390` |
| la biblioteca | 84 | 40 | ADELGAZA |
| arte previo (los dos) | 103 | 60 | ADELGAZA. Queda el dispositivo concedido, la realización y el linaje |
| **Münch** | 344 | 190 | **RUEDA LA MITAD.** Muere el recap de su método (cómo cuentan, la masa de 0 a 0.95, que lo propusieron como método general). Viven las TRES separaciones y la cláusula de ubicación |
| Del Core | 181 | 115 | ADELGAZA. Muere la oración de los grados de libertad (Results la tiene) |
| ruta clásica | 186 | 105 | ADELGAZA. Los tres ítems del lector a una cláusula; queda el puntero al apéndice |
| el propio miembro (los tres) | 362 | 185 | ADELGAZA. Muere la explicación de por qué no se nombra un umbral en canales; queda la cláusula en el perímetro |
| **modelo de la medición** | 403 | 55 | **RUEDA LA CABEZA.** La mitad del ruido: Methods ya tiene la premisa (`06_methods.tex:41`, "The emission model is white instrumental noise only", con pink y proporcional en cero), así que acá queda sólo la consecuencia, que es lo único que Methods no dice. La mitad del kernel: está impresa en `02_theory_full.tex:94`, entera. Transplante T1 |
| **lo que el diagnóstico no puede ver** | 157 | 0 | **RUEDA ENTERA.** Se muda a Methods, al lado de la especificación de la que habla (`06_methods.tex:349` ya está en ese vecindario). En la Discussion queda una cláusula dentro del cierre. Transplante T2 |
| el perímetro | 245 | 115 | ADELGAZA a la mitad. La enumeración del alcance se muda a Methods, que ya tiene dos de los cuatro ítems (`:9` dos estados, `:64` un solo salto). Quedan la regla 4Np(1-p) con su residuo honesto y los compañeros con la razón del kernel filtrado. Transplante T3 |
| el cierre | 47 | 55 | CRECE 8. Absorbe la cláusula del diagnóstico |
| | **3114** | **1500** | |

---

## Los ocho beats, con el presupuesto nuevo

Las notas "para el editor" de cada beat, más largas, están en el guion archivado. Acá va lo que hay
que escribir y cuánto.

**1. EL VERDICTO. 55.** Calibrado sobre casi todo el plano medido, la desviación medida donde no, el
rincón de pocos canales nombrado, y la reparación es el sandwich que el diagnóstico ya entrega. La
sección abre con un resultado y no con un descargo.

**2. QUÉ LO COMPRA, CON EL MECANISMO. 140** (100 + 40). Condicionar en los dos extremos; un extremo
solo no lo compra y en la celda de referencia cuesta; toda la distancia está en el update. VR en una
cláusula de 25 palabras con su calificación (la varianza predictiva divide la ganancia, así que las
dos mitades no son una separación limpia). Y Milescu pegado: la correlación local que él nombró es lo
que aísla la descomposición, cada intervalo reporta bien y la falla está en cómo se acumulan.

**3. LA APUESTA. 150.** Sube a la tercera posición. El flip state dejaba una firma que se señalaba en
el registro; la activación asimétrica de las mismas grabaciones no deja ninguna y descansa en razones
de evidencia, así que lo que una aproximación erra sobre la información pasa al mecanismo que elige.
Una oración de que la misma matriz entra en la evidencia y el error pasa a toda comparación construida
sobre ella, sin las fórmulas. Un puntero de doce palabras a que hace falta recalcular las evidencias y
que acá no se hace. Cierra con la frase que justifica el diseño simulado: ninguna cantidad de datos
grabados identifica el modelo generador.

**4. QUÉ PUEDE HACER UN EXPERIMENTADOR. 275** (100 + 60 + 75 + 40). El cociente gobierna a toda la
familia y τ_int lo lee desde una sola grabación, sin ensemble y sin verdad conocida, con el 2217 pegado
como advertencia de que lee régimen y no certifica. Después del agonista no hay más información sobre
N_ch ni k_on, y sí sobre k_off y la corriente unitaria, con la instrucción para el banco en una
oración. El ordenamiento de los bordes con la cláusula de mapa de conceptos. Y la biblioteca: corré el
diagnóstico en tu esquema y en tu punto, con la guarda de que sólo está el miembro calibrado.

**5. CONTRA QUÉ. 250** (60 + 190). El dispositivo Kalman concedido con el 1e-8, la realización en
tiempo continuo y el rango de canales como lo específico, y el linaje. Münch entra directo a las tres
separaciones (qué está dentro del veredicto, qué devuelve cada una, qué cuesta cada una), la cláusula
de ubicación con el miembro recursivo instantáneo, y la generalización de la emisión que les fuerza su
propia clausura.

**6. EL COSTO Y LA RUTA CLÁSICA. 220** (115 + 105). El costo es real y ahora está cuantificado en
unidades de incertidumbre mal reportada; la única cifra es la publicada, ofrecida como escala y no
como benchmark; dónde el ruido instrumental es grande los baratos recuperan; el sesgo del primer
momento como cargo aparte. Y la ruta clásica es fuerte donde hay sweeps intercambiables, con los tres
ítems del lector en una cláusula y el puntero al Appendix de reparaciones.

**7. LO QUE ESTE PAPER NO DIRIME. 355** (185 + 55 + 115). El propio miembro no está exento, con los dos
números de la desviación de dos lados, el bound de no-normalidad, los cinco candidatos refutados y la
conjetura en las palabras de Luciano. El modelo de la medición en 55: la consecuencia del ruido no
gaussiano (que la distorsión pertenecería al modelo de ruido y el score no separa las dos), y una
cláusula de que el kernel uniforme se adopta a sabiendas y su residuo se calcula de antemano, con el
puntero a Theory. El perímetro con la regla 4Np(1-p), su residuo honesto y los compañeros.

**8. EL CIERRE. 55.** El lector aporta su número de sweeps, su rundown, su tamaño de modelo y su
tolerancia y decide; lo que no puede medir por su cuenta es si el intervalo reportado es el que se
entrega. Absorbe una cláusula sobre el límite del diagnóstico.

---

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

---

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

55 + 140 + 150 + 275 + 250 + 220 + 355 + 55 = **1500**, más los dos títulos. El margen es cero:
cualquier cosa que se recupere sale de otra.

## Deliberadamente afuera, y de dónde salen si alguien los pide

- Los números del censo y del coverage: `04_results.tex:442` y `:218-224`.
- La derivación del bound de Bessel: `02_theory_full.tex:94`.
- La ACF del score por miembro: `04_results.tex:319-321`, y la relación 1.6-1.7× entre residuo y
  score en `03_diagnostics.tex:5`.
- Los pisos de alcanzabilidad, 17 y 52 canales: `04_results.tex:764-766`.
- La generalidad más allá de los canales: `04_results.tex:424-425`. Si alguna vez se quiere como prosa,
  va en la Introduction.
- Los ratios de costo por evaluación (1.4× a 5×) y el piso de nsim: no se re-derivaron de ningún
  artefacto primario de este repo. Viven en los comentarios de `10_appendix_repairs.tex` y no se
  imprimen.

---

# APLICADO (2026-08-11). Entregado contra presupuesto

La sección está reescrita sobre este guion y los tres transplantes están ejecutados. **1852 palabras
contra 3670 de la mañana, un 50% menos**, y 352 por encima del 1500.

| bloque | antes | presupuesto | ENTREGADO |
|---|---|---|---|
| veredicto | 69 | 55 | 61 |
| qué lo compra | 152 | 100 | 99 |
| Milescu | 62 | 40 | 56 |
| la apuesta | 244 | 150 | 181 |
| el cociente | 194 | 100 | 102 |
| washout | 137 | 60 | 72 |
| el mapa + la biblioteca | 221 | 115 | 143 |
| arte previo | 103 | 60 | 90 |
| Münch | 344 | 190 | 228 |
| Del Core | 181 | 115 | 149 |
| ruta clásica | 186 | 105 | 122 |
| el propio miembro | 362 | 185 | 228 |
| modelo de la medición | 403 | 55 | 93 |
| el perímetro | 245 | 115 | 159 |
| el cierre | 47 | 55 | 62 |
| | **3670** | **1500** | **1852** |

**Los transplantes, ejecutados:**
- T1 en `02_theory_full.tex`, al cierre del párrafo del kernel: la estimación publicada de más de diez
  veces, con la guarda de que son nueve esquemas y no una evaluación.
- T2 en `06_methods.tex`, después de la oración que licencia al simulador como referencia: el párrafo
  entero de lo que la comparación no puede ver. En la Discussion queda la cláusula del cierre.
- T3 en `06_methods.tex`, en dos lugares: los dos ítems de alcance que faltaban, junto a la misma
  oración de licencia, y **p = 0.5 derivado de los parámetros** donde se dan los parámetros, con las
  dos propiedades sobre las que descansa la regla de transferencia (varianza máxima, asimetría nula).

**Cuentas del paper entero:** la Discussion baja 1818; Theory sube 36 y Methods 247, así que el
manuscrito baja 1535 netas y ninguna limitación salió del paper salvo las seis cabezas firmadas.

**Los 352 que faltan.** Seis bloques se los llevan casi enteros, y en cada uno lo que queda es
contenido y no prosa comprimible. Para cerrarlos hay que borrar algo de la lista "qué NO se toca", que
son cuatro candidatos y ninguno es indoloro:

1. **La conjetura del invariante**, 54. Es de Luciano y verbatim del audio; `approach.md` pide que la
   refutación de los cinco candidatos se vea en el paper. Borrarla deja el "no se conoce el origen"
   sin la única cosa positiva que se dice al respecto.
2. **La generalización de la emisión de Münch**, 45. Es la mitad generosa de la comparación (su
   generalización les fuerza una clausura de momentos, que es el mismo tipo de movimiento que las dos
   nuestras). Sin ella la comparación queda sólo en separaciones.
3. **La regla 4Np(1-p) con su residuo**, 47. Es la única respuesta accionable a la objeción de p = 0.5,
   y `00_abstract.tex` nota 1c la registra como parte de la respuesta a la objeción de dos estados.
4. **La ruta clásica entera**, 122. Es una concesión, y es la que compra credibilidad frente a un
   referee de electrofisiología clásica. Es la más grande y la más barata de justificar si se borra
   (Results no la contradice), y la más cara si el referee resulta ser ese.

Mi orden si hay que elegir: 4 primero (una concesión que el paper puede hacer en una cláusula dentro
de Del Core), después 2. Las dos primeras (la conjetura y la regla) no las tocaría: son lo único que
la sección dice y no está medido en Results, y las dos son respuestas a objeciones que ya aparecieron.
