### Orden de Archivos (Cronológico)

Basado en los metadatos de los nombres de archivo proporcionados, el orden lógico cronológico es el siguiente:

1. **PTT-20251117-WA0006.opus** (17 de noviembre)
2. **PTT-20251118-WA0006.opus** (18 de noviembre)
3. **PTT-20251119-WA0015.opus** (19 de noviembre)
4. **PTT-20251123-WA0008.opus** (23 de noviembre, secuencia 08)
5. **PTT-20251123-WA0010.opus** (23 de noviembre, secuencia 10)
6. **PTT-20251124-WA0000.opus** (24 de noviembre, secuencia 00)
7. **PTT-20251124-WA0001.opus** (24 de noviembre, secuencia 01)
8. **PTT-20251124-WA0002.opus** (24 de noviembre, secuencia 02)
9. **PTT-20251124-WA0003.opus** (24 de noviembre, secuencia 03)
10. **PTT-20251124-WA0004.opus** (24 de noviembre, secuencia 04)

---

### Transcripciones

#### 1. Archivo: PTT-20251117-WA0006.opus

Bueno acá abro un nuevo grupo Macro IR 9.
Sin saber bien qué voy a decir.
Porque bueno, ya más o menos está todo dicho en cuanto a los próximos pasos que tengo que hacer.
Entonces me estoy concentrando ahora en la figura uno.
Para la figura uno, eh lo único que era un poquito complejo era que yo quería mostrar
la evolución de este
de lo que sería como el número de canales activados
durante el tiempo.
Eso lo casi lo que lo podría
tomar de
de de
del programa como está ahora formulado que
que tiene la posibilidad de darme eso.
Podría tratar de hacer eso sí.
Y este...
Y bueno.
Y si hago eso entonces
yo podría este...
Sí, sin tener que
que eh...
correr el ruido dos veces, sí, podría.
Claro, porque el tema es el siguiente, yo quiero mostrar este cómo durante un un intervalo de tiempo
este los canales se abren y se cierran varias veces.
Entonces que que esos eventos se pierden y que se viola la condición
de markoviana no de de medición instantánea.
Bien, el tema es que para mostrar eso
este claro tendría que mostrar el ruido también, porque bueno las señales van a tener su ruido
y como el ruido aumenta con la frecuencia, entonces va a tener un ruido tan grande que no voy a poder mostrar lo que quiero mostrar
porque lo único va a haber es ruido.
Entonces por eso yo tengo que
hacer algo un tanto artificial si se quiere que es este mostrar
simplemente la
este el número de de
canales abiertos por unidad de tiempo que lo saco de la simulación.
No es cierto.
Bien.
Entonces este
bueno.
Está bien, entonces la idea es que
que nomás sin modificar demasiado el código
puedo puedo mostrar eso
y después lo que sí tengo que tener es
la este la nueva función que va a ser este promediar ¿no? promediar este
o resamplear, eso, resamplear.
Resamplear un una corriente este con un nuevo experimento ¿no? tengo el experimento original
mis observaciones, un nuevo experimento, entonces ahora resampleo.
Y el resampleo solamente es válido si el resampleo digamos este yo eh colapso mediciones.
Claro podría ser el resampleado en realidad este indicando las las mediciones que voy a
que voy a colapsar, eso está mejor porque
Pero o podría hacerlo este con un nuevo experimento y
y ese nuevo experimento lograr generar un
un algoritmo de de colapso de de mediciones
Podría ser también.
Eh...
los algoritmos estos de
que yo desarrollé para el para el fractioning.
Eh...
Bueno.
Tengo que resolver esa cosa, entonces si yo ya tengo eso, ya tengo
eh la figura uno en cuanto tengo la corriente, tengo la corriente promedio, y bueno después lo que tengo que tener son este
claro, las derivaciones para el likelihood que sería bueno el
la probabilidad eh digamos
sí.
Sí sí sí los p
eh a prior los p posterior.
El tema es si pongo el p prior y el posterior, ¿ves? ahí por ejemplo yo
claro tendría que generar un nuevo coso que me genere esos datos, el prior y el posterior.
Y no estaría mal tenerlos también ¿no?
Sí sí.
Es decir el prior y el posterior.
Eh sí, lamentablemente entonces tengo que
que generar un nuevo... es medio complicado pero bueno.
Sí me parece que es medio inevitable. Y sí para tener el prior y el posterior no
no hay mucha salvación.
Eh bien.
Sí prior y posterior y
sí.
La el la corriente esperada
claro la medida, la diferencia, todas esas cosas.
¿Eh sí?
Eh te voy a sacar todo eso
te voy a mostrar un poco cómo funciona.
Claro nuestro un algoritmo el algoritmo eh digamos Macro R o sea tengo que mostrar el Macro en R, el Macro R y el Macro IR el Macro
claro IR.
Y claro y después muestro en supplementary el Macro A
ATR.

#### 2. Archivo: PTT-20251118-WA0006.opus

Bueno hoy fue un día intenso con MacroIR.
Creo que ayer también.
Eh...
Bueno logré eh
generar los primeros resultados como para
el paper que es un resultado bastante pavo que es simplemente el crecimiento exponencial de la corriente
pero dentro de el eh de MacroIR ¿no? o sea generado desde comandos de MacroIR.
Yo puedo generar eh tanto la simulación este punto por punto de los
de los número de canales abiertos
y este y luego la la corriente resultante de ese intervalo.
Digamos eso como para mostrar
eh ilustrar cómo se este de dónde viene el problema ¿no? de que
digamos se abren y cierran canales durante el período de medición.
Bueno eh también digamos lo que hice es simplemente toda una estrategia de
eh template metaprogramming que es ponerle una cosa que la llamé simulation tag
en la cual puede derivar de distintos indicadores de cosas que tengo que poner o no en la
en el este ¿cómo se llama?
en la...
eh en el resultado final de la simulación ¿no? o sea si yo simulo y nomás tiro todo o si simulo y guardo toda la trayectoria de corriente o la trayectoria del número de
de canales en cada uno de los estados.
Entonces con eso bueno puedo decidir qué hacer y
y este y bueno y eventualmente este mostrar cosas más complejas este si quiero.
Bien.
Lo bueno es que esta misma estrategia la puedo usar para también este
guardarme, o sea guardar, o sea eh
no sería guardar sería sí este
retener
datos intermedios del algoritmo para este poder mostrarlo, por ejemplo el prior ¿no? el
porque yo estoy por ahora guardando el posterior pero estaría bueno también tener el prior y el posterior para mostrarlo en el gráfico.
Entonces eso también lo lo voy a implementar y con esas dos cosas este creo que voy a poder tener
eh a full la eh ¿cómo se dice?
el eh
eh la figura de las figuras del paper.
Eh o sea que eso estoy estoy bastante bien, creo que estoy encaminado.
En principio si no me distrajera demasiado DataBound tendría que poder solucionarlo esta semana por lo menos tener una figura ya armada
y este y otra figura la semana que viene y ya ya escribir definitivamente el paper.
Creo que creo que estamos estoy bastante cerca de eso.
Este sí.
Digamos como que que podría tenerlo para principios de diciembre este un manuscrito como para enviar a BioArchives y y luego ponerlo en en a consideración de eLife.
Eh sí lo que me va me va a matar es lo de lo de DataBound
que tengo que de alguna manera lograr porque yo estaba ya estaba convencido que
este Sorab se iba a la India y entonces me iba a dejar tranquilo durante dos meses, pero no, todavía está acá.
Este yo la verdad no entendí bien
cuánto tiempo se queda y cuánto se va
y bueno y eso un poco me me parte un poco al medio porque la verdad es que yo trabajo mejor cuando tengo una sola cosa en la cual preocuparme.
Si tengo varias un poco un poco me desarmo, no no funciono del todo bien.
Es como que
necesito mucha mucha maquinaria para funcionar.
Medio extraño eso que tengo, sí.
Pero bueno listo ya está este ahora eh
esto lo ahora tengo que pasar a DataBound a a ver qué qué hago para Sorab, a ver si lo puedo liquidar en una o dos horas y luego pasar de vuelta a MacroIR.
Voy a tener que hacer eso nomás.

#### 3. Archivo: PTT-20251119-WA0015.opus

Bueno, sigo acá con MacroIR.
Eh...
Bueno, hoy lo que hice fue eh
limpié un poco el código que hice ayer, un código de metaprogramación para poder este
agregar variables a las vector space este de acuerdo a lo que uno quiera. O sea,
las agrego con compile en compile time.
Y este, y entonces este,
digamos con un solo algoritmo
este, pu- genero todo tipo de combinaciones de
de qué quiero obtener, ¿no? En este caso particular era, bueno,
quería generar la, digamos la simulación de corrientes
y bueno y puedo guardar por ejemplo las corrientes este simuladas en cada etapa de
del sampleado eh de la variable multinormal de,
de multinominal,
que,
del digamos del
de la matriz de probabilidad P de transición
de estado, ¿no? O sea, de cada matriz de transición de estado  yo, yo sampleo
esa matriz y guardo el valor. Entonces con eso digamos este puedo guardar todo el recorrido y bueno, compararlo y ver,
qué sé yo, más que nada para debuggear sirve eso, pero bueno o para mostrar lo que pasa.
Este...
Bueno, y lo bueno es que claro, voy a usar la misma maquinaria para el Macro,
Macro R,
este para guardar el, justamente quería guardar el prior y el posterior, el prior además del posterior,
cosa de poder mostrar cómo funciona bien el
el algoritmo.
Este...
digamos es un overkill de alguna manera hacerlo así pero,
pero no sé me deja más tranquilo que queda una herramienta este flexible y poderosa.
Bueno, una vez que haga eso bueno ya ya ahora sí estoy en en la
el umbral de hacer la figura.
Estoy, me distraje bastante porque bueno salió el nuevo
el nuevo este versión de Gemini y decían que era muy bueno y
y sí, de hecho le ponés el el código de de mi repositorio de GitHub y me lo analiza todo, me lo lee todo.
Este es bastante bastante poderoso y bueno ya
ya lo usaré para,
bueno ahora compite directamente con ChatGPT a ver cuál de los dos es el
el mejor, ¿no?
Y no sabemos bueno qué sé yo, ya
ya eso es otra historia.
Y ahora lo que tengo que hacer, bueno tengo que ver, los, estoy medio cansado, el tema de DataBound
que tengo que generar las predicciones diarias que no sé por qué no
no no quiero predecirla, no, que no me animo a hacerlas. Bueno ahora las voy a hacer
este ahora voy a contar lo que voy a hacer en el próximo audio y listo.

#### 4. Archivo: PTT-20251123-WA0008.opus

Bueno, no sé hace cuánto que no grabo una puta idea nueva.
Estado bastante ocupado con MacroIR, o sea bastante compenetrado con
con escribir código.
Fui más o menos bien, qué sé yo, estoy bastante contento con eso.
Eh...
logré encontrar un error en en lo que era la
implementación de la corrección de Taylor.
Ahora con DeepSeek logré una
una eh
demostración de una fórmula que tiene mucho más sentido de la que tenía antes, que
que era incoherente porque
me daba que
que había un drift de la media sin
aún en ausencia de cambio en la
en el la diferencia entre la media esperada y la obteni-, la medida,
eso no podía ser.
Siempre me llamó la atención.
Bueno y eso no lo, eso no lo tiene, así que debo haber cometido algún error de de cuenta en el pasado.
Bueno todo esto en realidad es MacroIR.
Eh, no sé por qué lo digo acá en ideas.
Lo voy a pasar a MacroIR.

#### 5. Archivo: PTT-20251123-WA0010.opus

Macro IR nueve.
Estoy en Macro IR 9.
Ha sido una semana intensa para MacroIR.
Como decía,
eh todo mi desarrollo, no todo mi desarrollo porque la parte teórica del desarrollo de la Taylor expansion de el la la función de likelihood
estaba bien, pero me parece que en alguna cuenta le erré
y este ahora bueno este
con DeepSeek sacamos una una expresión que es bastante más compacta, digamos totalmente interpretable
de acá a la China, o sea se entiende bien.
Este y la verdad
eh que eso daría para un paper muy lindo este
eh para publicar en en eLife o algo así también.
Así que y este bueno y lo que me queda es extenderlo
para el caso de este intervalo.
Este y ya con eso, bueno no.
Lo extiendo para el caso de intervalos, agrego
el tema de la
de los filtros de Bessel
y este y después bueno implemento la back
backward gra- gradient
e implemento un un stochastic, no un stochastic descent no, sería algún
algo así, no sé, algún este MALA creo que es
alguno de esos este
algoritmos que no son el el este
Hamiltonian Monte Carlo porque Hamiltonian Monte Carlo este necesita de demasiados parámetros para tunear
y la idea es que no haya mucho que tunear
y que este el tuneo sea más bien por por la la escalera termodinámica, ¿no? esa es un poco la
la idea.
Este...
O sea que haga un un este un te...
eh un cla- no sería un... sí un parallel tempering
este pero no sé si con ensembles, podría ser con ensembles también la verdad.
Este...
Quizás con ensembles tenga buena un buen poder de diagnóstico
instantáneo.
Eh también cabe la posibilidad de
de plantear algún tipo de
se estaba pensando como un filtro de Kalman
eh como reemplazo a las este a al Metropolis Monte Carlo. Eso lo lo voy a preguntar ahora a ChatGPT.

#### 6. Archivo: PTT-20251124-WA0000.opus

Eh, una cosa que me impresionó mucho de el avance este de la expansión de Taylor
usando DeepSeek que encontré unas fórmulas que tienen sentido
es que
digamos el producto de la covarianza con el vector de
de conductancias diferenciales, llamémoslo conductancias diferenciales a la suma de las conductancias más
este peso por por la
eh la varianza de la conductancia.
Eh, ese vector
digamos si todo fuera idéntico, obviamente vos no
no tenés este cómo optimizar nada ¿no? O sea...
digamos eh
digamos el algoritmo debería tender a
a asignar a digamos aumentar la diferencial del grado de información de ese vector.
Bueno, eso por un lado. Y la otra cosa que estaba pensando es cómo usar eso para, o sea ver, digamos cómo ese vector cambia con el número de canales y con el
claro con el número de canales y con la eh el salto temporal.
Y con eso tratar de entender eh dos preguntas ¿no? o sea cuál es el salto temporal
y el número de canales ideal, porque vos podés por ejemplo,
o sea claramente tenés más información pero también hacés más cálculos
si vos subdividís, o sea si vos promediás corrientes o o tomás o no, ¿no?
Entonces, ¿cuánta información perdés promediando corrientes?
Eh y cómo, quizás
cómo promediar de una manera inteligente
en la cual vos digamos así puedas...
claro como separar la varianza de la media ¿no? o sea
la varianza de la media lo genial es que vos
con la media tenés la información total y con la varianza tenés la información que está dentro de la diferencia entre las
samples.
Entonces si yo lograra hacer algo así
con las mediciones, con la autocorrelación, podría acelerar mucho los
los cálculos.
Y lograr lo que siempre soñé que es
eh yo tengo una imagen como de una, como si fuera una fotografía ¿no? Que vos una fotografía
vos tenés un arreglo de
de
partículas de cloruro de plata, no sé qué son, listas para para ser este
estimuladas por la luz, las disponés y bueno cada
cada lugar es independiente del otro, ¿no? Solamente se estimula con la luz que llega
exactamente a un punto del espacio.
No tenés que deconvolucionar las cosas, esa es la brillantez de la fotografía.
Entonces, si lográramos un sistema
de que te lograra eh
deconvolucionar
digamos los este
eh...
claro tendrías un sistema lineal ¿no? Obviamente, porque lo lineal es lo único que
que permite separar las o las este contribuciones de cada punto independientemente de las otras.
Claro, si uno lograra hacer algo así entonces vos podrías
transformar digamos un registro en una una matriz de
de transiciones
y representarla.
Claro que sea en un en una sola en un solo paso.
¿Es eso posible?
En principio sabemos que no, que necesitás operaciones no lineales para eso.
Pero no sé, es interesante la pregunta. O sea, por ejemplo la existencia de las transformadas de Fourier
indica de que ese tipo de operaciones son posibles. Entonces la pregunta es cómo es que la transformada de Fourier
logra ser más que bueno tiene que ver con
con la
bueno eso con la separación de
senos y cosenos  pero además bueno con el la transformada rápida de Fourier
que cuya magia no no entiendo del todo.
Esa es mágico y habría que entenderlo bien.
Digamos cómo poder lograr tener esa magia en más lugares.
Eso es eso sería muy bueno, la única manera de avanzar en estas cosas.
Estos métodos intensivos son medio inútiles, no, están buenos para
para qué sé yo
publicar papers o divertirse pero
pero no logran avanzar nada el conocimiento bien, o tenés estas cosas que escalen como N o menos.
Bueno.
Encontrar esos métodos, desarrollarlos.
Esa es la única manera.
Sino todo lo demás es es medio al pedo me parece, pero bueno.
Nada, todo eso son reflexiones.
Pero sí, la la idea de dividir en
en este, o sea
o sea yo tengo qué sé yo tres registros y los sumo acá en mi mi coso, en mi bueno
dividirlos de vuelta y
y ver cómo combinar la suma de las de la división.
Quizás ahí hay haya algo, no sé, tengo que ver. Bien.
Sí, eso, hay algo ahí, obvio.
Hay cómo cómo hacerlo sistemático y cómo relacionar eso con la estructura de las
de de esos vectores, esas matrices, eso es es un desafío
interesante.

#### 7. Archivo: PTT-20251124-WA0001.opus

A lo que voy es que por ahí estudiando esas matrices se ve cuándo uno obtiene información separando y cuándo no y cuánta información, cuánto
cuánto más cómputo por cada extra bit que sacás, cada pregunta es cuánto cómputo sacás realmente.
Y cuánto, digamos, información nueva sacás por cada
mayor número de cómputos.
Y esa es la pregunta, como una especie de eficiencia informacional computacional, eso es lo que hay que entender.
O sea cuándo vale la pena separar y cuándo vale la pena promediar.
Claro, ¿cuándo promediás? Al promediar vos eliminás información, esa es la pregunta. Claro.
Ese es.
Ese es el tema.

#### 8. Archivo: PTT-20251124-WA0002.opus

Bueno, eso tiene una respuesta fácil, digamos en todo lo que son fluctuaciones aleatorias al promediar perdés la información, es obvio.
Las fluctuaciones fuera del equilibrio, no, no las este
no las perdés.

#### 9. Archivo: PTT-20251124-WA0003.opus

O sea las respuestas fuera del equilibrio las respuestas transitorias al al promediarlas las ves con mayor claridad.
Las respuestas las fluctuaciones
aleatorias, cuando estás en el equilibrio, ahí las perdés al promediar, eso es
creo que es bastante obvio.

#### 10. Archivo: PTT-20251124-WA0004.opus

Ahora y claramente en el equilibrio donde digamos la transformada de Fourier es este reina ¿no? Ahí sí usás transfourier y
y
podés este
sacás hasta el último detalle, o sea ahí sí
eh
sí, ahí ahí es donde
donde
ese método brilla.
Entonces sí, una especie de
de transfor- digamos de mmm
digamos sí eso que se llama transformada que se va moviendo, que no sé cómo se llama, eso tendría que estudiar un poco a ver si
cómo se puede entender.
Eh
ese sería el método.
Que tendríamos que tener de alguna manera.


**Orden de archivos:**

1. PTT-20251125-WA0004.opus

---

### Transcripción: PTT-20251125-WA0004.opus

Bueno acá vuelvo con macro IR.

Estuve como dos o tres días reescribiendo toda la teoría de macro IR y macro Taylor y macro Taylor IR.

Que son distintas aproximaciones a la función de likelihood de las corrientes macroscópicas.

Que tienen en cuenta, bueno por un lado, la variabilidad en la varianza entre estados, que esa es la de Taylor.

Y la integración basada en el intervalo de medición, que es la IR, y la combinación de ambas.

Bien.

Solo me quedaría incorporar el filtro de Bessel.

Y con eso ya tendríamos un cuasi gold standard por así llamarlo de macro R.

Digamos de ya, ya sería, ya me costaría imaginar algo mejor para agregar.

Bueno, lo que se podría agregar es una especie de meso R que es, digamos, un sistema que pueda hacer micro R para pocos canales y macro R para muchos, o sea un sistema adaptativo.

Esa sería la única cosa que podría mejorar o algún tipo de micro R estocástico, no sé.

Pero lo que es puramente macro R, digamos, creo que después de agregar el filtro de Bessel quedaría bien.

Y sería, este, sí.

Y solo me quedaría mejorar por el lado de los algoritmos de evidencia.

Que eso ya tuve una discusión con Gemini o Google o no me acuerdo cuál.

Que de cuáles eran de acuerdo al número que, bueno tendría que ir a o Hamilton Monte Carlo o alguno sin derivadas.

Este uno que, y quizás usando librerías de Python.

Eso me quedaría para después.

Entonces, con todas esas cosas más o menos me queda implícito un plan de publicación de papers.

Que tengo que ver un poco, o sea el tema es que todas estas cosas podrían ir en un solo paper grande.

Pero no tiene sentido hacerlo porque son demasiados avances conceptuales para un solo paper.

Y no me representan en mi currículum más que un paper.

En cambio, y no, digamos, y tengo un techo de cristal en cuanto a la categoría de revistas que las cuales puedo mandar.

Ese techo más o menos estaría en eLife.

Ya superar eLife, ponele ir a PNAS o a Nature Communications es imposible si no me alío con otra gente.

Entonces, lo mejor que puedo hacer es hacer un poco de salami slicing, sería no sé, algo más trozos de carne.

O sea de, no, digamos, de empaquetar en pequeños, en digamos, trozos digamos conceptualmente establecidos, no, bueno, macro IR, macro Taylor, macro Taylor IR.

Bessel.

Serían como cuatro papers y después bueno.

Evidencia usando algún método de mejor para grandes números.

Y este, y bueno y después ver.

Digamos lo que me quedaría ya a futuro es el tema de hacer, este.

Bueno, micro IR hacerlo.

Hacerlo, digamos, estocástico o algo así, cosa que que pueda ser computable.

Esa sería digamos otra, otra línea más de investigación posible.

Este, y con eso más o menos cubriríamos todos los huecos que hay en cuanto a el número de canales, este etcétera para para hacer la determinar la cinética.

Este, y bueno, no sé.

Probablemente, sí.

Eh.

Bueno por otro lado desde el punto de vista de la verificación de modelos, bueno, obviamente tenemos el método este de las del gradiente, la esperanza del gradiente sea cero con una varianza igual a la Fisher Information Matrix.

Que eso lo voy a meter ahora en este primer paper y, y digamos eso que sea como una especie.

Bueno, junto con el macro, macro IR, que son digamos las dos maneras de medir si la likelihood es correcta.

Y bueno después tenemos el tema este de múltiples canales.

Ah, y finalmente lo que me quedó en la cabeza es que este yo hago modelos conformacionales y bueno si yo acoplo digamos distintos modelos conformacionales, este el número de estados escala violentamente.

Y eso hace que que no funcione demasiado bien los algoritmos.

Y entonces ahí se plantea la la pregunta de si hay manera de de hacer una especie de de.

De, en lugar de establecer cada uno de los posibles estados en el canal.

La posibilidad de hacer distribuciones de estados, o sea con covarianza entre los estados, digamos de los distintos estadios conformacionales, cómo estarían este.

Digamos correlacionados con otros estadios.

Y de esa manera digamos modelar el acoplamiento alostérico, digamos de un punto de vista a lo estadístico.

Eso sería como un posible punto extra.

El tema es que bueno, es que en realidad vos no.

Eh.

Digamos macro R funciona porque tenés muchos canales y acá cada canal sería como.

Bueno en realidad sería una subunidad y tendría dos, tres subunidades nada más, o sea.

El número de subunidades es demasiado bajo.

Entonces estamos en una situación la cual va llegás a cien, doscientos estados fácilmente y.

Y ahí no sé.

Ahí digamos quizás lo mejor sea como una especie de de modelo donde vos solamente tenés registro de los estados más probables del canal y a la vez sea una especie de micro R.

Micro IR.

De los, digamos que.

Solamente tienen en cuenta.

Digamos una configuración posible de de los canales por vez.

Entonces sería como.

O sea yo lo que me imagino es, así como vos, digamos no es que vos explorás todas las combinaciones de parámetros posibles, sino que simplemente elegís una y punto.

Bueno, ahí lo mismo, o sea nosotros elegiríamos un estadio de digamos microscópico del canal y chau y bueno después veríamos otro, o sea.

Como que hacer un Metropolis sampling pero también dentro de los estadios del canal, entonces con eso la ventaja sería que bueno, no tenés que esparcir digamos que poner toda la matriz.

Este, del del canal, sino simplemente una matriz más reducida, no sé.

Esa es una idea digamos que que tengo que explorar.

Pero que sería un golazo para este poder trabajar en en estados más complejos todavía, o sea.

Donde a, ahí lo que importa es el número de.

De interacciones, que aumente digamos algo así como linealmente con el número de interacciones y no exponencialmente con el número de interacciones.

O sea, la idea ahí en el, en esto, sería romper el escalamiento lineal exponencial que tenemos con el número de interacciones a y llevarlo a una cosa lineal.

Que sería simplemente bueno, explorás algunas, no explorás todas.

Y bueno, de eso lo tomás como una sample.

Entonces, este nada, esa esa sería un poco una idea a futuro.

Digamos grande, como idea grande que no está desarrollada es, está esa.

Este, que sería bastante buena porque bueno, claro, es como una idea de que.

A ver yo me la imagino algo así como que vos.

Digamos hacés este, sería MCMC all the way down, ¿no?

O sea, que todo sea estocástico.

Y este y vos no, no probás todas las las combinaciones.

Sino que nada, este realmente fitteás hasta la última partícula de polvo.

Simulás eso, y bueno tomás este muestras así de todo eso y este.

Y bueno.

Y bueno, y que Dios te ayude.

Sería, pero bueno, es es un poco esa la idea.

Entonces.

Que es como un poco romper la unidad de canal, es como que todos, todos son al final partículas, ¿no?

O sea.

Todas son partículas este y bueno, una partícula forma un canal, muchas de esas forman digamos un conjunto, un ensamble de canales, etcétera.

Pero.

Pero digamos algo, algo así.

Y bueno y eventualmente claro lo ideal sería entrenarlo con técnicas de redes neuronales y que eso tenga una interpretación bayesiana.

Bueno ese es otro tema, ¿no?

O sea cómo hackear digamos el, los LL, no serían LLM, serían las eh, Deep DNN.

Que son las Deep Neural Networks.

Eh, desde el punto de vista este.

Eh, bayesiano, ¿no?

Si vos podés hacer una interpretación bayesiana de las DNN.

Eh, que claro ahora es eh.

Large Language Models, LLM, este es como que totalmente eclipsó a las DN, DNN.

Que eran lo que, lo que teníamos hasta hace dos años.

Eh, bien.

Sí, esa esa es una una pregunta importante, ¿no?

La digamos interpretación bayesiana de las DNN.

Eh, bueno, pero eso es otra historia.

Bueno, más o menos.

Esa sería digamos mi panorama con la cinética de canales, o sea.

Como vemos tengo tela para cortar un par de años en realidad.

Entonces, cuál es un poco mi mi mi futuro científico es bueno.

Es sacar estos papers por lo menos hasta.

Yo creo que digamos para poder tener una herramienta que sea útil este y que no me sienta incómodo y bueno sí, tengo que.

Bueno en realidad lo que tengo que hacer, que era lo que me había planteado antes, es estar un poco seguro de macro IR que funcione, ¿no?

O sea no hace falta necesariamente el macro IR con.

Macro TIR.

Con Bessel.

Pero bueno yo creo que eh, no es tanto laburo de código.

Hacerlo, así que yo creo que sí.

Que tengo que llegar a Bessel y tengo que llegar a un método eh.

Que no sea MC, porque MC digamos tiene una tasa de aceptación muy baja a a bajas, con muchas dimensiones a baja temperatura.

Entonces ahí claramente necesito algo que funcione mejor, algún método mejor.

Este, eso lo lo voy a explorar inmediatamente cuando termine con Bessel.

O sea que ese sería así como el, no sé, el quinto paper, no sé.

Este, la cagada es que eso recién ahí eh voy a poder este volver a tener datos experimentales, ¿no?

O sea ahí realmente yo recién ahí me sentiría cómodo haciendo este un análisis de datos experimentales, digamos ya con con digamos todas las herramientas bien.

Ah claro y ahí lo que me quedaría es este otra cosa para ver es el tema de los canales únicos, ¿no es cierto?

O sea.

Como, y bueno y también este de otras otras este tipo de de mediciones experimentales, qué sé yo.

Bueno difusión de moléculas únicas y motores moleculares, etcétera.

Todas esas cosas este cómo tratarlas.

Lo que pasa es que en realidad digamos la pregunta es.

Tiene sentido que yo me meta en todo eso.

Eh, y la respuesta es medio que no.

Porque eh yo creo que una vez que más o menos, o sea vamos a ser pragmáticos, o sea yo ahora eh.

Digamos ya tengo digamos la punta de cuál es mi mi.

Mi, mi futuro científico es bueno.

Es ser el Peter Turchin argentino, o sea.

Eh, es bien vendible, o sea es bien claro, o sea Peter Turchin es un biólogo que se dedicó a hacer modelos de la historia y bueno yo puedo hacer lo mismo, digamos.

A que que me que me la, me la chupen si no, o sea.

Lo puedo hacer.

Tranquilamente, sin drama.

Este.

Y bueno eh.

Y ese sería más o menos, yo creo que para para junio del año que viene yo en cuanto tenga digamos un poquito de respiro.

Yo ya tendría que turchinearme, ¿no?

O sea el proyecto de turchiniación mía este que por debajo sigue el proyecto de la HoTT, ¿no?

De de la Homotopy Type Theory pero eso como parte de la.

Del lenguaje.

Que indique eh.

Un lenguaje eh de programación que sea semánticamente este, que calcule semánticamente, etcétera, etcétera usando la Homotopy Type Theory como como una base.

Y bueno y la, el primer caso que yo quiero trabajar para eso es el caso de eh la las este.

¿Cómo se dice?

Base de datos de bases de datos, ¿no?

O sea tomar como un objeto una base de datos, es una fuente de información.

Y eso este.

Y claro y es la fuente de información de la fuente de información también, ¿no?

O sea es muy interesante eso.

Conceptualizar eso de una manera que sea rigurosa y este.

Y digamos rigurosa matemáticamente y que tenga digamos que tenga un mapeado semántico claro, eso para mí sería digamos una contribución este muy buena.

Eh, que digamos tiene que ver, que es en paralelo con digamos mi Peter Turchinización, ¿no?

O sea, yo haría eso este dentro de mi proyecto de bueno de modelar algo, o sea base de datos del Indec, lo que sea.

Este.

Y claro.

Claro estaba pensando.

Cómo es, por ejemplo base de datos del Indec vos lo que tenés que tener es algún tipo de, por ejemplo qué sé yo, anécdotas de cómo se tomaban los datos del Indec.

Y este y ver si eso es posible, no sé.

Eso es interesante.

Porque a mí lo que más me interesa es este claro que los datos sean falsos, ¿no?

O sea vos estás midiendo una cosa pero bueno qué qué tipo, o sea si te inventan puntos, etcétera, etcétera.

Yo no, digamos, o sea yo veo mucha eh motivación política debajo detrás de de datos este estadísticos.

Y bueno, digamos los datos son datos, o sea el tema es cómo vos los interpretás, lo podés interpretar como como justamente.

Eh, como un una forma de ver la intervención política sobre el Indec o o los deseos de cierta clase que maneja el Indec o lo que sea, digamos.

Este eso me parece que tiene eh tiene en sí mismo.

Un.

Una.

Digamos este, sería sería interesante elaborar teorías acerca de eso, ¿no?

Eh.

Bien.

Bien, eso es un punto para mí para mí súper claro que.

Que hay que manejar, entonces claro uno tendría que empezar con con cosas que este supuestamente son más creíbles o que entendés más cómo se generan, ¿no?

O sea quizás cuando sean menos creíbles entendés más cómo se generan porque eh puede ser la proyección, ¿no?

O sea.

Este qué sé yo, por ahí los datos por ejemplo estadísticos chinos que son los que ellos quieren hacer, bueno te está diciendo lo que ellos quieren hacer y eso es un dato más allá de ponerle si vos no creés en eso, ¿no?

Eh, bueno no sé.

En un punto digamos el valor de los datos siempre está dado por la la entropía y esas cosas que son son fáciles de medir y entonces ahí de alguna manera vos siempre algo tenés, o sea.

Las fuerzas de la naturaleza se se escapan siempre, ¿no?

No no hay manera de bloquearlas del todo.

Eh, eso es interesante.

O sea.

El tema es saber exactamente qué significa un dato, ¿no?

O sea esa es esa es la pregunta y y bueno.

Eh.

Y bueno y eso es lo que yo me me quiero empezar a meter eh ya a partir de no sé, abril del año que viene, una cosa así.

Cuando ya tenga el.

O sea ya tenga dos papers, yo tengo que tener dos papers este aprobados y un tercero enviado.

Ahí recién yo me puedo empezar a meter con con el tema de de modelos de de la Argentina.

Eso eso me queda súper claro que es así.

Digamos no no y eso más o menos supongo que será por abril.

Con un poco de suerte, quizás antes, pero bueno.

Y ahí no sé porque ahí bueno yo quisiera digamos terminar.

O sea la verdad como para sentirme, o sea vamos a pensar.

Si yo quiero abandonar canales.

Y que ya no no me moleste, lo que tengo que tener es una maquinaria aceitada que me permita obtener resultados de una manera este reliable y.

O sea tengo que hacer la máquina de hacer chorizos eh biológicos.

A ver cómo sería eso.

O sea la idea sería.

Digamos si es como en el proyecto de los siete locos que los siete locos iban a hacer la revolución con putas, bueno yo voy a hacer la revolución con eh con macro R.

Macro R me va a dar de comer porque voy a tener digamos una maquinaria que me permite fittear modelos cinéticos de lo que venga.

Este.

Y bueno ahí yo lo meto asiduamente, produzco, entrego y mi mi mente está en otro lado.

Entonces yo tengo que llegar a ese punto.

Y no estoy en ese punto porque por ahora los modelos son muy, tardan mucho, etcétera, etcétera.

O sea yo yo estoy este siendo demasiado exagerado en este en algunos puntos.

O sea yo lo que tengo que lograr optimizar con macro R.

Claro ahí está, esa es la cuestión.

Yo tengo que tener el tema de la.

La tasa de obtención de información estadística, ¿no?

O sea.

Eh, que eso pasa por bueno reducir el problema y este llegar digamos en un tiempo dado, o sea el al a a a obtener información este relevante en eh digamos dado bueno complejidad del modelo, complejidad de los datos y tiempo de cálculo, ¿no?

O sea o poder de cálculo.

Esos son los tres los tres puntos.

Entonces este.

Eh, claro.

Ahí ahí hay una cosa que.

Que eso eso entender digamos esa ese triángulo eh teóricamente eh sería eh sería digamos yo creo que ese es un punto muy importante que quizás sea el más el más importante que tenga que hacer con macro IR es ese.

Es como una especie de balanceo automático de la precisión, ¿no?

O sea.

Este.

Que claramente, o sea.

Eh, y bueno.

Digamos, o sea.

Y bueno.

Digamos si vos tenés muchísimos datos por ahí este claro cómo vos podés promediar sin este perder demasiada información, ¿no?

O sea, claro.

Eso es es el tema.

Claro, justamente yo lo que veo es cuando vos al promediar perdés información, ¿no?

Si vos tenés fluctuaciones aleatorias y las promediás, perdés algo, te queda todo determinista.

Entonces la pregunta es, vos podés este esas fluctuaciones aleatorias eh codificarlas con por ejemplo que yo algún modelo de digamos de unos cuatro o cinco canales abiertos y cerrados y eso después vos proyectás modelos más complejos en esos menos en esos eh sistemas este simplificados pero generalistas.

O sea.

Ahí hay un tema.

Ahí hay un tema muy interesante también.

Eh, que tengo que que es fundamental, ¿no?

Que es como una especie de eh digamos de cómo yo llevo de un este yo tengo un modelo eh alostérico complejo y todo y bueno yo eso lo reduzco a a un modelo más simple que sea trabajable, ¿no?

O sea.

Eh, cómo hacer eso.

Ese es un tema que que está bueno entender.

Eso por un lado.

Eh, que ese modelo más reducido sea el que digamos absorba toda la información de fluctuaciones este.

Y eh claro después bueno vos puedas irle agregando complejidad.

El tema es que claro obviamente si vos tenés más complejidad necesitás eh ver todos los datos.

O sea hay una parte en que en que necesitás todos los datos, pero por lo menos podés este explorar hasta llegar a ese punto donde vos estás entre dos o tres modelos y entonces ahí sí vos hacés un un testeo con una muestra más grande, algo así tendría que ser.

O sea como digamos no hacer, o sea la idea es esa, es como no hacer cálculos al pedo, ¿no?

O sea no digamos vos tenés que hacer los cálculos razonables, una especie de razonabilidad de los cálculos.

Y un balanceo automático de la razonabilidad de los cálculos, ese algoritmo digamos sería una especie de meta meta algoritmo.

Eh, creo que sería lo más lo más importante que podría hacer ahora.

Eh, lo más, sí.

Especie de meta meta algoritmo, o sea que vos elijas el algoritmo que más conviene en estas circunstancias, ¿no?

Que también es parecido a a la idea de las bases de datos, ¿no?

Bases de datos de bases de datos serían algoritmos de algoritmos, ¿no?

O sea yo ahora tengo que entrar en la meta en el meta mundo, ¿no?

O sea.

Eh, en un punto digamos, o sea una vez que termine de desarrollar estos métodos que, o sea de exponerlos, de mandarlos afuera y todo y tener digamos un.

Un cierto nombre o conocimiento en eso, lo que tengo que hacer es este agarrar los métodos que existen, las librerías que existen y hacer una especie de librería de librerías que este desde algún lenguaje de.

De computación vos puedas este generar programas que gene, que corran, ¿no?

O sea que podría ser en no sé, o en.

Habría que ver cuál es el que más conviene.

O sea obviamente ya teniendo un un scripting ability este vos tenés todo creo, ¿no?

O sea.

O sea en el fondo us, usando bash lo podrías hacer siempre.

Este.

Eh, sí.

O sea vos lo que tenés es una especie, o sea lo, los otros programas pasan a ser objetos de tu de tu programa, ¿no?

O sea.

Eso, esa idea eh, eso es lo que yo tendría que que hacer, ¿no es cierto?

Eh.

Sí.

Eh, y bueno, ahí habría que ver bien cómo se puede hacer eso.

Bueno se puede hacer digamos con un número más o menos reducido, pero no sé, habría que ver, o por ahí si tenés, si a vos parásitas ciertas herramientas eh.

Digamos todo es cuestión de tener un meta lenguaje eh.

Un lenguaje de de cómo instalar programas y y cómo eh correrlos y todo, entonces vos tenés este un programa de altísimo nivel que después elige el el programa que lo hace, ¿no?

Eh no está mal eso, me parece que eso es lo que tendría que que hacer, pero por lo menos digamos o sea empezar digamos o sea la idea es como hacer un lenguaje de lenguajes.

Eh pero que claro, que sea HoTT, ¿no?

O sea HoTT en el sentido de de de bueno que las funciones sean sean eh proposiciones, etcétera, etcétera.

O sea que la semántica esté bien bien eh puesta en ese lenguaje.

Y eh que la la ejecución digamos este se delegue en en otros en otros programas, ¿no?

Algo así.

Eh.

Eso es es es está muy bien.

Habría que desarrollarlo primero un poco desde el punto de vista teórico y después bueno implementarlo.

Pero.

Pero claro, la idea de eso es que entonces eso sí pueda elegir la herramienta que vos a, vos vas a usar para cada caso este y bueno, y no sé, con.

Algoritmo que sea.

Digamos.

Esa sería un poco el plan, ¿no?

Ese es un buen plan y ese plan, obviamente, ese programa de programas lo que usaría es para para hacerlo de Peter la la Peter Turchinización este de de mis este, sí, de eso.

Ser un proyecto de país.

O sea, yo lo que quiero hacer es eso, es que eso esté al servicio mío, bueno de la de la economía probablemente, lo que lo que yo decía es que digamos lo mejor es es predecir el el precio del dólar este como salida, ¿no?

Este.

Igual eso tengo que ver bien qué qué es lo que quiero hacer.

Pero sí.

Pero la idea es hacer como una especie de inteligencia este inteligencia humana, no sé.

Algo así.


**Orden de archivos:**

1. PTT-20251126-WA0007.opus

---

**Transcripción:**

Bueno, luego del audio más largo de mi historia, vamos a hablar de lo que se me ocurra aquí en la reserva mientras veo unas hermosas golondrinas. Que ya, ya hacen verano porque hay unas cuantas. Ahí se va.

Bueno. [silencio] ¿En qué punto estoy? Yo lo que estaba pensando es hacer con tres papers simultáneos de macro R, como si fueran Matrix 1, 2 y 3, y que serían, bueno, el de... la expansión por intervalo, la expansión de Taylor y la combinación de ambas.

La idea es, este... eh, bueno, mostrar los papers, o sea... yo creo que tengo que decirlo así abiertamente, que bueno, planteo los tres papers como uno a continuación del otro, eh... hago el "salami slicing" por una cuestión conceptual, porque son, este, conceptos complejos y que vale la pena explicarlos con... con, digamos, detenimiento cada uno.

Eh... el último, el tercero que tendría la combinación de los dos métodos, sería el que finalmente haga una... una comparación, este, clara de ellos. Este... los primeros, los otros dos serían para, este, mostrar, bueno, eh... que el problema que solucionan y co... cómo lo solucionan y qué es lo que... qué es lo que ganan.

Entonces, este, ¿en qué estado estoy como para hacer los... los tres papers del buen humor? Eh... tengo... eh... tengo, digamos, las fórmulas, las estuve trabajando de nuevo porque la parte del intervalo estaba bien, digamos, de hecho saqué un paper con eso, pero la parte de Taylor, la expansión de Taylor estaba mal. Estaba bien conceptualmente, pero no sé, había un error, algunos errores de cuenta medio extraños, me aparecían términos que no tenían mucho sentido.

Ahora tengo un, digamos, lo esperado, que es que si yo... no... si mi... eh... medición es igual a la esperada, no modifico mi prior. O sea, digamos, que sería lo elemental. Que en el caso de... de la... de la anterior, este, no... no era lo que ocurría, digamos, con... con la implementación que yo tenía, que se ve que algún error cometí en el momento.

O sea, gracias a DeepSeek hay una simplificación, que... un vector V corta que suma las contribuciones de la... media y de la... varianza de la conductancia. Este... y entonces con eso, este, digamos, la cosa... las cuentas se simplifican, también se entiende mejor, etcétera.

Entonces, en realidad lo que tendría que hacer es, este, bueno, presentar... este, los resultados de, o sea, de... del algoritmo, mostrar que en cierto... que... mejora... que funciona mejor que... que no hacerlo, digamos, ¿no? Que la... la aproximación esta es mejor que no hacerla. Y para eso, bueno, claro, para eso tengo que tener el micro R funcionando. Y la... varianza... la FIM... eh... la Fisher Information Matrix, etcétera.

Todo esto ya lo... lo sabía, no... no estoy aportando nada. Eh... nada nuevo. Eh... pero bueno. Tengo que llenar un poco el espacio el de... con algo, ¿no es cierto? Eh...

Entonces, ¿qué me queda? Ah, lo que me estaba quedando es, que es lo que anoche no podía dormir bien... eh... con la aproximación de Taylor, en realidad vos no... es una aproximación inexacta. Este, no es como la aproximación por intervalo que es... es exacta. Sino que esta es inexacta y lo que hacés es una aproximación de Newton, o sea... este... sucesiva y de... de maxima... maximum likelihood, ¿no? O sea, encontrás la... el punto de... de la... del gradiente que tiene... que es mínimo, que es cero. O sea, busco que el gradiente sea cero, eso lo declaro mi posterior... media... mean... y tomo el Hessiano en ese punto como el... este... el Hessiano.

En realidad, digamos, yo tendría dos Hessianos, ahí hay... hay una cuestión... que es que... podría tomar la media de los dos, no sé, eso habría que ver. Es una cuestión, una pregunta interesante. Porque si vos tenés... yo uso una... una aproximación, este, de Newton, o sea, yo no tengo que... tengo que el Hessiano no es constante, y el Hessiano no es constante porque depende, entre otras cosas, de... del delta, ¿no? De la diferencia entre lo observado y lo... este, y lo medido.

Eh... pero lo observado y lo medido en realidad es constante. La pregunta es... el Hessiano... es... este... claro, no, el Hessiano no... no cambia, creo. A ver, déjame pensar un poco. O sí cambia. Mmm. Esa es una muy buena pregunta. ¿Cambia el Hessiano o no cambia? Esa es la pregunta primera que tengo que hacer. Déjame ver. Mientras tanto, total... quemada.

A ver. El Hessiano es la... o sea... eh... vos tenés dos contribuciones... o sea... o sea, la pregunta es cómo entra P acá. O sea, que sería la... eh... claro, sí... eh... sí. Esa es la pregunta, cómo entra P. Que es la... Sí. Esa es la pregunta. Claro, y P... claro, entra... entra en Sigma... Vamos a... tengo que formularlo mejor. Tengo que pensar bien la... ahora la... sí... la expansión de Taylor. La tengo que pensar bien. Bueno, ese es el tema. La tengo que... que implementar bien. No lo puedo hacer acá en el aire.


**Orden de los archivos:**

He analizado los metadatos de los nombres de archivo para establecer el orden cronológico lógico:

1. **PTT-20251127-WA0000.opus** (27 de noviembre, secuencia 0000)
2. **PTT-20251127-WA0001.opus** (27 de noviembre, secuencia 0001)
3. **PTT-20251128-WA0004.opus** (28 de noviembre, secuencia 0004)

---

### Transcripción: PTT-20251127-WA0000.opus

Bueno, se me acaba de ocurrir un nuevo algoritmo del macro IR que sería macro IR doble recursivo, no sé cómo llamarlo, pero la idea es la siguiente. La idea tiene que ver con este... digamos, ¿cuál es el problema? El problema que tenemos es que si tenemos un impulso, en la respuesta a un impulso, que lo saca del equilibrio a un canal, entonces se empieza a relajar. Entonces, al relajarse, digamos, van decayendo las distintas constantes de tiempo, entonces vos necesitas tener mucha resolución al principio y luego se van equilibrando, ¿no? Primero las constantes de tiempo más rápidas se equilibran rápido, después las medianas en tiempo mediano, las largas en tiempo largo.

Entonces, lo ideal es tomar tiempos en una escala exponencial, es decir, que vos vas incrementando exponencialmente los tiempos. Primero un milisegundo, después dos, cuatro, ocho, así. Ahora bien, si vos tomas ocho milisegundos como el tercer punto, por ejemplo, o dieciséis, vos estás perdiendo la información que está entre el milisegundo cuatro y el milisegundo ocho. Entonces, ¿cómo vos podes usar esa información? Bueno, esa información ya no tiene tanto valor en cuanto a fluctuaciones determinísticas, sino como fluctuaciones estocásticas. Entonces vos podes ver, por ejemplo, en esos milisegundos eh cómo eh la digamos... cambia lo que sería la varianza en intervalos de un milisegundo, de dos milisegundos, de cuatro milisegundos. Eh digamos... porque fluctuaciones rápidas este van a aparecer más este... menos promediadas en intervalos de tiempo más cortos y en intervalos de tiempo más largos van a aparecer atenuadas y van a tener más preponderancia las fluctuaciones más lentas, ¿no? Entonces vos tendrías este lo que sería equivalente como a un análisis de frecuencia, este una especie de gráfico log de tiempo este de integración y vos vas a ver que digamos distintos componentes tienen distintos tiempos.

Entonces en principio debería aparecer la información este de esas fluctuaciones este explícita, ¿no es cierto? Ahora, ¿qué es lo que ocurre? Eh... uno podría predecir esos valores de varianza, porque ya lo hicimos, es la varianza de la media, pero claro, es la varianza de la media este digamos un poco más complicado, tomando este la varianza de la media... no sé, algo así... que sería eh vos este tu elemento es este intervalos de... vos tomas intervalos de un milisegundo, entonces tenes una media de esos intervalos de un milisegundo, y ahora eso vos este eh calculas la varianza de de esos intervalos de un milisegundo y este claro, y lo que tendrías que calcular es la varianza de esa varianza, ¿no? O sea cómo digamos cuáles son... porque después vos lo que queres hacer es comparar la varianza que vos medís contra la varianza esperada dado los parámetros del sistema.

Entonces todo eso digamos implica un cierto desarrollo este que yo lo intenté hacer en su momento y ya me daban este matrices de grado cuatro, o sea digamos tenía que integrar en grado cuatro, no sé qué. Pero yo creo que usando este método que consiste en eh... ay no me acuerdo cómo lo llamo yo, pero que vos vas llevando digamos la la media y la varianza por intervalos y después este vas este digamos eh recursivamente calculándola a medida que vas agregando intervalos, uno podría llegar a estimar esto.

Entonces, no sé, la idea sería poder digamos este tener para cada este para cada intervalo de medición, vos ya no tenes un número lineal de puntos, sino vas a tener un un logaritmo del número de puntos este de como como dato, ¿no? O sea, por ejemplo, si tenes treinta y dos puntos, por ejemplo, vas a tener este uno que es la media de esos treinta y dos puntos, otro que es este la varianza entre la primera y la segunda mitad, después la varianza entre el primer cuarto, el segundo cuarto, el tercer cuarto y el cuarto cuarto, después eh una varianza entre el primer octavo, segundo octavo, así, y hasta llegar digamos a los treinta y dos puntos. Es decir que vos este la información de N puntos la este codificas en logaritmo de N puntos, logaritmo base dos de N puntos. Eh y estás promediando digamos lo que sería la varianza. O sea que es una forma muy eficiente de extraer información.

Claro que vos perdes la información eh transitoria, ¿no es cierto? Entonces por eso digamos al principio vos este el primer punto digamos es como que lo lo analizas este distinto a los otros, ¿no? O sea... haber si vos tenes por ejemplo... si tenes dos puntos, bueno, elegís el primero y el segundo punto. Si tenes cuatro puntos, elegís el primero, el segundo y el cuarto punto, por ejemplo. Si tenes ocho puntos tenes el primero, el segundo, el cuarto y el octavo. Y además tenes este la varianza de eh claro, pero excluyendo el primero, ¿no? O sea tendrías a ver... a ver cómo sería... pensémoslo un poco.

Si vos tenes cuatro... tenes ocho puntos, ¿no es cierto? Tendrías el primero, el segundo, el cuarto... después tenes la varianza entre el segundo y el cuarto... y después tenes este la varianza entre el cuarto, el sexto eh y el octavo, y además la varianza entre el cuarto, quinto, cinco, sexto, sexto, séptimo, séptimo, octavo. Y así recursivamente este vos podrías tener digamos todos esos todas esas este variables que que fitear. Claro, y fitear... digamos vos esto eh... te quedaría que un problema de este claro... de N puntos te quedaría no sé un logaritmo de N más este dos veces, ponele, no sé, una cosa así.

Y eh... y casi que eso vos lo podrías fitear de una sin este sin necesidad de... no, a ver, tendrías... no, igual usarías un algoritmo recursivo, ¿no? Porque vos este... bueno podrías hacerlo, fitearlo de una, son aproximación, y después hacer el algoritmo recursivo en el cual digamos vos este tenes el posterior digamos de, bueno dado la primer medición, la segunda medición, la cuarta medición, etcétera, tenes ahí un un prior de eso y tendrías este que calcularías el posterior dado todas las otras observaciones. Eh... suena como un algoritmo súper súper eficiente eh que pierde un poco, pero no tanto de información comparado con el algoritmo que sería usar todos los puntos, que ahí tendrías todo. Está bien, y este bueno podría ser como una aproximación y después este... Claro, vos lo que digamos la idea es tener cajas, cajas chinas donde tenes aproximaciones más groseras afuera y y así a medida que vos vas este disminuyendo la temperatura, o sea en temperaturas altas tenes aproximaciones bien groseras y a medida que vas bajando la temperatura tenes aproximaciones más y más precisas, cosa de de no gastar cómputo al pedo. Esa es un poco la lógica que habría que implementar para tener digamos un algoritmo espantosamente eficiente.

---

### Transcripción: PTT-20251127-WA0001.opus

A ver si se me acaba de ocurrir una idea que soluciona un problema que hace mucho tiempo yo tenía en mente, ¿no? O sea, la pregunta es... eh si vos tenes una corriente en equilibrio, o sea un sistema marcoviano en equilibrio, o sea, si vos lo tenes fuera del equilibrio lo lo soltás a que vaya al equilibrio y ves un decaimiento exponencial. Ahora, si vos lo tenes en equilibrio no ves decaimiento porque estas cerca de la media, pero en realidad es mentira eso, porque nunca estas... estas cerca de la media, pero vos estas a una raíz cuadrada de N no sé qué de la media.

Entonces claro, si uno pudiera promediar todos los puntos que están más o menos lejos del equilibrio, o sea partiendo de ese punto ves toda la trayectoria que sigue, es como que vos estás promediando puntos que que se parecen y este como que está la probabilidad del estado de estar por abajo del equilibrio, por arriba, es es grande, o sea es este mayor que un medio, entonces dado eso, bueno, vos vas a tener un vas a tener un decaimiento exponencial o hacia un lado o hacia el otro, ¿no? hacia la media o o de arriba de la media hacia abajo.

Entonces, claro, uno podría directamente, digamos, agarrar los datos y haces una especie de convolución, pero una convolución solamente, o sea un promedio, ¿no? O sea, eh tomo digamos los datos, ponele de de un punto que está en el equilibrio, desde I uno hasta ponele I uno más N. Después veo I dos. Si I dos es mayor que la media, entonces lo pongo en el equipo mayores que la media, desde I dos hasta I dos más N. I tres, si es mayor que la media, bueno el equipo mayor de... si es que es menor que la media, entonces va al equipo menor que la media, lo sumo ahí desde I tres más N. Y así. Y así todos, digamos los voy sumando a uno o a otro, de acuerdo a que sean mayores o menores que la media.

Y entonces ahora voy a tener dos curvas que necesariamente empiezan de un valor que es mayor o menor que la media y va decayendo eh, o sea... eso suena bastante este, o sea es muy fácil de implementar y debería debería ocurrir que uno ve este claro, el tema es que eh eso sirve claro cuando vos tenes equilibrio, ¿no? Cuando vos vos este llegaste a un equilibrio, o sea tenes saturación de agonista y lo dejas ahí un tiempo. Mmm. Claro que es lo que no ocurre generalmente porque vos lo que tenes es una desensibilización, o sea, eso serviría, digamos en el caso de canales iónicos, solamente para casos donde estas cerca de la desensibilización.

Pero bueno, podes tener por arriba o por debajo de la media esperada también, ¿no? Eh en cuyo caso, bueno es un poco más complejo porque supuestamente no ha llegado al equilibrio, pero bueno podes también hacerlo así. Sí, no sé, es un digamos, claro, no no no tiene tanta aplicabilidad digamos esta esta herramienta de análisis, pero pero nada, era interesante. No, creo que digamos el lo anterior es más tiene una aplicabilidad más general que es más o menos, digamos tiene... no, no es parecido en nada, pero bueno. Pero la idea es este ir promediando tiempos cada vez más largos eh sí. Y ver este...

No, claro, pues eso lo que pensaba era la, claro la variabilidad lejos del eh claro si uno pudiera digamos eh simular un un cambio eh digamos como correr el algoritmo con alta resolución lejos de del inicio para ver si cambia algo del del sistema este al irse relajando, ¿no? O sea al estar lejos del del pulso inicial, si después el sistema es diferente, que debe serlo porque bueno entra entra en estados extraños, pero bueno eh nada eso. Pero sí, no, evidentemente lo tenes que digamos correr contra el el el estándar es que sean que incluya todos los puntos, ¿no? Eh hay que ver esto hasta qué punto uno puede hacer cosas este eh sucias y rápidas contra... Lo que pasa que lo lo que es este limpio y y no dan los tiempos para para correrlo, ¿no? Es es infinito. En fin eh bueno.

---

### Transcripción: PTT-20251128-WA0004.opus

Hace rato que no desgrabo. Vamos con macro IR. Bueno he avanzado muchísimo esta semana... eh digamos diseñé varios algoritmos nuevos, uno lo Git que consiste en hacer una transformación logit para trabajar con pocos canales. Y bueno, y después el el mixto que mezclaría un microscópico con un macroscópico, que ese bueno, no sé, está la idea no más, no no está implementado. Bueno en un audio anterior puse también otros más que que tenían que ver con intervalos este que se van incrementando y y podes ver... bueno es medio complicado porque no no me queda bien claro cómo funcionaría, ¿no? que es el de que vos tenes este simultáneamente eh varios eh intervalos de de tiempo y varianzas, etcétera. Eh y se lo que tendría que hacer es como ir corriendo sucesivamente este ir ajustando sucesivamente la los priors, ¿no? de las eh eso.

Bueno, en fin. Eh... o sea la idea sería bueno, vos partís de el prior y haces con con la media. Después haces el prior y haces con la varianza. Después si haces con, digamos el habiendo tenido el posterior de la media y varianza, haces el posterior de la eh varianza en intervalo dos. Y después el intervalo cuatro, etcétera. O sea como vas haciendo este posteriors de los posteriors, ¿no? O sea esa manera de de correrlo asesinaría los los ruidos, eso me parece que estaría bien. Habría que pensarlo un poco, pero en principio podría solucionar el problema y además simplificaría un poco todo.

Bueno, entonces, ¿en qué ando? Ahora bueno tengo que simplemente terminar mi mi paper de macro IR. Hace ya un par de semanas que estoy atascado porque atendiendo otras atendiendo otras cosas, bueno lo que hice es elaboré un documento donde muestro cómo eh cómo manejar que eso es un avance grande, es cómo manejar cuando el algoritmo se sale de de madres, es decir, cuando se sale del simplex de eh de los este del vector que define la eh la probabilidad de estar en cada estado, que es un vector de probabilidades, entonces las probabilidades no pueden ser menores que ceros o mayores que uno, y bueno si se sale de eso hay como un mecanismo que aumenta la la el ruido digamos de de medición para impedir que eso pase y bueno y eso se suma a la a al likelihood, ¿no? como como un ruido extra. Digamos me parece como una forma razonable, no no sé, es ad hoc digamos, pero bueno por lo menos como que se interpreta que que bueno si pasa eso tenes una medida que es medio chota, entonces bueno le crees menos al al ruido. Eh no deja de ser ad hoc digamos, pero bueno, no sé, tiene su elegancia dentro de lo ad hoc que es.

Este bueno ya con eso ya ahora me pondría a implementar macro IR. Eh yo creo que no voy a implementar, estoy digamos ahora la duda es por eso salí a pensar si iba a implementar macro Taylor IR o no. Este... la verdad es que no sé, eh podría serlo, o sea, porque ya lo tengo casi implementado el macro IR, podría hacerlo el macro Taylor, o sea la idea original era que el macro Taylor imitaba a Munch, por eso es que me metí en este quilombo. Este analizando las cosas no entiendo bien qué es lo que hace Munch, este la verdad tampoco es que es tan importante, entonces lo podría hacer simplemente es mostrar mi propio análisis que es hacer un uno de digamos el macro R, macro R que no sé cómo llamaría que sería no Taylor sino de digamos de media, macro MR, y después bueno él lo llamaba el macro de R que ahora lo llamo IR. Bueno el macro MR, que sería macro media eh recursive, que debería mostrar que con eso no alcanza, ¿no? y bueno necesitas el inter el IR. Después yo creo que necesitas el Taylor porque si no te quedan digamos este eh no, digamos el tema es que no tenes eh buena estimación del del número de canales y de la conductancia, eso te sale mal y eso me pareció en el paper. Así que yo creo que eso es un un efecto real. Este pero bueno quedaría para otro otro paper. En principio, pero bueno la pregunta no es eso sino si lo implemento ahora o no.

Eh... la verdad que no sé. Creo que no es tan difícil implementarlo eh pero es un poquito, pero no no sé si tanto. Eh yo lo que haría es primero haría el el Taylor no más, el Taylor, el el IR como está, y este y bueno y le quitaría la parte esta de eh que uso de digamos como que degrado el algoritmo de acuerdo a si algunas cosas digamos eso que safe eh el safety no sería cambiando el algoritmo sino aumentando este este este factor alfa, algo así de como inflación de varianza, eh digamos sería el mecanismo, o sea que mantendría siempre este doble intervalo recursivo, o sea no degrado el algoritmo. Este bueno lo mismo vale que en realidad para todos, ¿no? O sea eh todos los que serían eh macro R tendrían que tener ese mecanismo, o sea ese mecanismo es universal.

Este... y bueno y tendría digamos eh seis seis algoritmos que tratar que serían los macro R, macro MR, macro eh IR y los este eh NR, ¿no? los de cada uno, macro NR, macro MNR. Y creo que el eh eh no habría un un macro I I un macro INR porque sería el macro MNR, porque digamos siendo doble no no tiene un efecto, o sea sigue siendo el mismo algoritmo, o sea son cinco en realidad. Eh... hm. Sí, está bien porque no. Sí sería más o menos eh lo mismo. Eh... Sí, a los efectos de este de este gráfico no no el macro Taylor no aporta nada realmente, o sea el macro Taylor aportaría bueno después si yo muestro eso que la dirección este cambia con con el delta y la varianza cambia con el delta eh una sutileza un poco mayor. Eh... Sí, igual lo importante del macro Taylor es ese tema de estimar bien la la el número de canales y la conductancia, ¿no? O sea que lo estaría estimando mal en la conductancia.

Eh bien, eh... el macro Taylor no habría... bueno podría tener un eh... sí, o sea yo tengo de rango dos o de rango uno y tengo que sí implementar el rango uno por una cuestión eh si se quiere eh eh estética o pedagógica, pero sí porque sí, especialmente cuando cuando lo saqué para mostrar que bueno que eso es es necesario, es como un un partecita del paper, ¿viste? para ponerle un poquito de de juguito de de algo que decir. Sí. Entonces en principio iríamos con macro IR a full. Este no yo la verdad me parece que no voy a hacer nada de macro Taylor ahora, haría solamente macro IR y este y me concentraría en terminar esto y ver ver digamos y bueno si evidentemente hay problemas eh con con la varianza con la que yo creo que va a haber problemas porque ya los hubo eh eh hubo eh uh. Eh... entonces sí eh eh lo señalaría, ¿no? que el gradiente de la del número de canales y el gradiente de la de la conductancia va a dar medio mal. Eso es casi seguro.

Eh entonces y después bueno ahí lo lo mandaría con el macro Taylor porque eso es la solución. Mmm. Sí, no sé. Eh... Vamos a ver digamos un poco va va con las ganas, si tengo fuerzas para hacer macro Taylor yo lo haría digamos por una cuestión de que bueno lo tengo todo en la cabeza ahora y es más fácil así, si no después se me va a ir. Eso es un buen argumento para hacerlo, la verdad que sí yo creo que es un argumento bastante fuerte. Eh sí. Si es un argumento muy fuerte para tratar de terminar todo esto. No de, yo no haría el logit ese ya sería un poco too much y bueno menos el Bessel este igual lo tendría que hacer en algún momento pero bueno. No, yo creo que ahora eh bueno antes de fin de año yo tengo que mandar este este paper de macro IR eh solo y bueno y tener quizás el macro Taylor, después bueno macro el macro logit y el macro Bessel que estarían después, en realidad antes que el logit tendría que estar el macro Bessel, este y este igual en realidad ah, más importante que eso ahora me doy cuenta es el micro R, el micro R estoy en pelotas total este y ese es importante y era las razones por las cuales uso modelos tan simples. O sea creo que que más que pelotudez, sí o sea bueno igual tengo que hacer el el macro Taylor eh y después hacer el micro R y este y bueno sí.

Porque a ver la figura uno ¿qué es lo que tendría que tener? La figura uno tendría que tener este bueno es una figura para yo la había pensado eso ya ya creo lo había hablado en algún lado, serían varios algoritmos, este cómo funcionan, ¿no? entonces este y las distintas variables para dos dos mediciones, ¿no? dijimos y serían bueno bueno yo tengo la la medición de de I, o sea este tengo la la I predicha, la varianza de I, eh la la mu predicha eh la mu posterior, en el caso de del macro R común y para el otro tengo digamos este la mu cero y mu uno prior y posterior este que la mu uno posterior pasa a ser la mu cero prior. Eh claro que ahí yo lo llamo de otra manera, ¿no? que tiene esos esos nombrecitos que estaban muy bien. Este y claro eh la pregunta es si meto o sea lo que me quedaría acá nuevo serían bueno dos puntos no que serían el delta y la dirección.

Este... eh no estaría mal eh qué sé yo mostrarlo quizás eh... no sé, no estaría mal en un gráfico igual yo no sé eso si es pedagógico eh pero no sé si es para una revista digamos así buena así es ahí es un poco dudoso ¿no? digamos hacer un un diagrama de cómo proyecta en digamos vos lo que tendrías es digamos el el lo que hace macro R yo veo cuál es la la diferencia entre lo predicho y lo lo lo observado y yo digo es un número natural, un número real perdón, entre menos infinito y más infinito, y eso lo proyecto en un un eje de un vector de del número de estados eh que es el el vector este G que no es otra cosa que la conductancia proyectada sobre el la matriz de covarianzas que al proyectar la matriz la conductancia, la matriz de covarianza vos lo que te aseguras es que este entre otras cosas se mantiene eh el valor de P, ¿no? o sea el se mantiene este P eh o sea la que sería el vector de de probabilidades como probabilidades, ¿no? sale suma a uno y tenes suma a uno básicamente. Justamente y lo otro que es que para garantizar que esté entre cero y uno bueno eso lo hago por afuera digamos del algoritmo lo lo le pongo un un forzamiento ¿no? eh esa sería un poco la el asunto del forzamiento es importante porque justamente es lo que hace que uno pueda correr el algoritmo bien digamos sin preocuparse incluso con modelos que estén muy lejos de de la realidad, digamos igualmente van a te van a dar digamos un un resultado, digamos, ¿no? o sea no no se van a no te van a dar este NaN o infinito qué sé yo, sino te va a dar un valor. Este eso está bien. Digamos es una una buena cosa.

Este... Bien. Eh... Y bueno. Claro, el tema es hacer ma micro R, tengo que ver cómo está micro R que que eso funcione. Igual no es muy difícil micro R, no sé tengo que verlo un poco. Claro lo que es difícil es claro tengo que fabricar la la matriz Q de micro R, la Q I Q I matriz esa es medio quilombo. Eh sí eso es quilombo no sé. Lo hice o sea que debería funcionar en algún momento lo hice este claro un momento había hecho una que incluso construía lo que hacía falta o sea si vos partís de de un estadio donde están todos cerrados bueno un solo estado y bueno vas incorporando los estados a medida que que hacen falta, ¿no? o sea eso podría funcionar muy bien en el caso en que eh por ejemplo eh claro que vos tengas eh hagas un una especie de Montecarlo Markov chain de eso, cual vos no no estás este llevando la cuenta de todos los las posibilidades de toda la distribución sino simplemente seguís algunos este algunas partículas serían ¿no? o sea como un un este sí un seguimiento de partículas y no de de toda la distribución. Entonces este y haces estadística desde ahí ¿no?

Esa me parece que es una una buena una buena idea, bueno eso habría que es el tema ese para trabajar este el año que viene ¿no? de modelado digamos este así estocástico que no no planteas todos los metaestados sino digamos este haces una simulación. Claro es un poco la idea de si se quiere de dinámica de proteínas ¿no? de dinámica de proteínas vos no no contás considerando todos los ensambles de posibles distribuciones de las proteínas sino que simulas una una trayectoria individual ¿no? y rezas por que esa sea significativa. Eh digamos sería la misma la misma idea. Este... claro eso lo podes hacer muchas veces y y bueno y te da algo no sé. Sí. Eh... Bueno no sé eh... Voy a dar una vuelta más a ver si se me ocurre algo más y ya me vuelvo a eh a casa a terminar macro IR de una vez. A ver si puedo juntar un poco de fuerzas para ahora terminarlo. La verdad estoy con poca poca creatividad, bueno qué sé yo no hay mucho que hacer.

Eh a ver el algo state... el algo state eh qué tendría que agregar, bueno podría agregar tengo que agregar definitivamente el tema este del alfa, eso sin lugar a dudas tiene que estar incluido porque es un un detalle que justamente tiene que aparecer en el paper, digamos momento es en que eso o sea tiene que haber una figura de eso del del parámetro alfa. Eh con lo cual eso tiene que estar. Eh la dirección la incluyo en el algo state pero bueno no la incluyo no más para incluirla. Eh y el delta y el chi que sería o algo así también. El chi era importante por el tema ese de claro la distribución ¿no? de de los chis eh que es una cosa que lo lo mostré en en el paper de biophysical journal también así que eso podría andar también.

Eh pero bueno en principio eh claro, en realidad digamos a ver el chi lo puedo sacar pues tengo la sí la predicha la varianza y ahí lo puedo calcular digamos tengo que guardarlo ya puedo calcular después o sea es al pedo es al pedo guardarlo yo lo lo calculo digamos post post processing no no tiene sentido calcularlo antes salvo para verificar que que no haya metido la pata ponele pero no no es no parece ser necesario. Eh... No, sí, el único es ese ese ese famoso alfa sí. Eh con lo cual el macro algorithm dejaría de de existir, mirá vos fíjate ese detalle. Hm. Claro. Y bueno uno podría redefinir eh sí digamos que el macro R macro IR cuando no no tenes este digamos cuando no tenes información bueno es como que la likelihood es cero y no no te cambia nada es como que tenes varianza infinita ¿no? o sea no te cambia nada. Es decir que puede funcionar como un algoritmo normal.

Eh... no, sí quizás macro algoritmo debería volar y eso lo reemplaza este este parámetro alfa. Eh... Bueno, un tema que me queda es a ver si Q de Q de T eh claro ahí qué salvo de qué qué guardo de Q de T ¿no? Es un tema también. Eh... En principio claro ahí tenes el tema ¿no? si guardas cosas pesadas para guardar serían bueno eh la conductancia condicional entrada y salida y bueno la varianza y la la conductancia media digamos es se podría calcular pero bueno es medio engorroso al pedo guardar igual no hay más no son esas tres en principio. Está el G total pero yo no sé eso si hace falta eh o sea ahí me ahorra como una multiplicación pero no sé por ahí es más fácil eh sí si no andar arrastrando tanta tanta sí bueno no es tantos cálculos me parece. No, por ahí por ahí vuela voy a volar voy a quedarme con con lo mínimo. Eh esa es una buena pregunta.

He realizado la **Prueba de Señal** en los 6 archivos proporcionados. Todos contienen audio válido e inteligible.

A continuación, presento la **Lista de Orden Cronológico** basada en la fecha y secuencia numérica de los nombres de archivo, seguida de las **transcripciones verbatim**.

### Orden de Archivos

1. **PTT-20251130-WA0005.opus** (30 de Noviembre)
2. **PTT-20251203-WA0021.opus** (03 de Diciembre)
3. **PTT-20251219-WA0012.opus** (19 de Diciembre - Secuencia 1)
4. **PTT-20251219-WA0013.opus** (19 de Diciembre - Secuencia 2)
5. **PTT-20251219-WA0015.opus** (19 de Diciembre - Secuencia 3)
6. **PTT-20251219-WA0016.opus** (19 de Diciembre - Secuencia 4)

---

### Transcripciones

#### 1. Archivo: PTT-20251130-WA0005.opus

"Un cuento borgeano de la imposibilidad del... una refutación del pensamiento, el pensamiento no existe porque un pensamiento está constituido por la interacción de varias neuronas. ¿Unas cuántas? Y estas neuronas tienen que mandarse información unas a otras y eso lleva un tiempo. Entonces eh ocurre en en un espacio-tiempo finito, es decir no, la la información de una neurona a la otra tarda qué sé yo, ponele no sé cuántos milisegundos tardará, pero debe tardar del orden de 10 milisegundos, por lo menos. Y entonces este el pensamiento es, además que bueno sabemos que el pensamiento abarca varios segundos, pero pero no el pensamiento ocurre ocupa un lugar, un espacio-tiempo, un volumen de espacio-tiempo. Existimos en un volumen de espacio-tiempo. No en un instante. En un instante no existimos. Somos un volumen de espacio-tiempo, somos tiempo, claramente, sí. Eso es, es así, pero lo loco es que nos vamos desplazando en el tiempo, ¿no? Es un tiempo, claro el tiempo fluye, ¿no? Uno no puede estar en el mismo lugar, pero no puede estar en el mismo tiempo. Qué locura eso. Y es es algo que, claro es un proceso que fluye y fluye. Sin parar. No se detiene el tiempo. Qué loco. No podés congelar. Y ser es seguir siendo. No se es, se sigue siendo. No no se sigue siendo y y es un instante que desaparece en cada instante, esa idea, sí. Que uno es en un instante para dejar de serlo y seguir siendo en otro instante. Es medio es este vertiginoso cuando uno lo lo piensa. Igual ¿por qué es vertiginoso en realidad? Si es ineludible, digamos si somos eso, por qué digamos eso que somos, que es un volumen de espacio-tiempo, nos parece tan tan extraño, ¿no? Estamos acostumbrados a quizás por una cuestión de de leer libros, que uno lee un libro y una página y la página está ahí. Quizás, digamos los músicos tengan una noción más parecida porque la música existe en el espacio-tiempo. Existe en el tiempo más que nada. Más que en el espacio pero pero sí, no sé. Es interesante todo esto. Pero bueno, no sé. No hay no sé si hay mucho más acá tela para cortar. Tampoco sería que quería cortar una tela acá. ¿En qué en qué se gatillaría alguna idea de todas estas cosas, esta rumiación filosófica? ¿En qué se podría podría plasmar en algo que que nos cambie el espacio-tiempo? ¿Qué pensamiento nos podría llevar a a ser diferentes o o o no sé, o superar alguna limitación o o ver algo mejor o no sé? ¿Qué sería en en nuestra vida eso? ¿Qué pensamiento nos puede liberar? ¿Hacer entender algo? No sé. No sé. Creo que por acá. No hay mucho más que pensar."

#### 2. Archivo: PTT-20251203-WA0021.opus

"Bueno con Macro IR me fui a la mierda con un cambio de las clases de de cómo organizaba los eh las variables en Macro R y hice una clase de clases de... complicado un poco las cosas quizás, no sé. Estoy pensando en en cambiar patch state, que sea en vez de sea una clase que sea un este un sinónimo de un vector y después este pegar ese vector y y simplificar un poco la la forma de trabajar. Este eso lo lo voy a ejecutar ahora. Eh y bueno, eso no más."

#### 3. Archivo: PTT-20251219-WA0012.opus

"Bueno. Estoy en la fase final, va fase final... ya tengo los datos como para armar la figura uno del paper y estoy totalmente bloqueado, no sé bien qué poner. Eh... Y están tocando la puerta."

#### 4. Archivo: PTT-20251219-WA0013.opus

"Y estoy totalmente trabado y no sé qué carajo hacer y es como que no ni tengo ganas de hacer nada. Pero bueno tengo que hacerlo. ¿Qué es lo que tengo que hacer? ¿Qué es lo que tengo que hacer en este paper? Lo pensé hace tiempo y ya ni me acuerdo. Lo que tengo que mostrar son... digamos tengo que presentar el algoritmo Macro IR para lo cual tengo que mostrar el problema de de que digamos las mediciones eh eh individuales esconden el promedio temporal. Eh... después tengo que mostrar eh los elementos de Macro IR que son eh claro. El o sea Macro NR solamente no muestra digamos la el ¿cómo se llama? La probabilidad predicha y la covarianza se se calcula a partir de la de la probabilidad de estado. El vector probabilidad de estado, eso te determina una covarianza. En Macro R vos actualizás la covarianza y la y la probabilidad media en cada después de cada operación, de cada medición. Y finalmente... Bueno y después eh lo que podés decir es bueno, está bien pero vos estás suponiendo ahí, o sea vos ¿cómo podés tener en cuenta el hecho de que de que vos estás tirando una una probabilidad y una covarianza media, o sea que digamos que vos... que no se... que cómo se dice... que no está representada por por lo que ocurre durante el el intervalo, porque vos suponés que que es instantáneo? Entonces vas a tener eh digamos no estás teniendo en cuenta la el drift en la probabilidad durante el intervalo de medición. Entonces eh probablemente tengas estés subestimando el el ruido, la variabilidad. O sobreestimando también, puede ser, depende de de qué es lo que, digamos cómo lo promedies. Entonces eh digamos vos podés cambiar tu tu conductancia por la conductancia media dado que vos empezás la digamos que tu estado inicial es conocido, el estado inicial del del intervalo de medición vos vos suponés que sea uno y entonces este vos ajustás eh tu probabilidad de acuerdo a a la media a la corriente media esperada dado que empezás en este estado y después ajustás tu probabilidad a posterior y después lo lo actualizás. Y en el caso de la de de intervalo de intervalo recursivo vos este vos de alguna manera actualizás tu probabilidad a posteriori pero ya en el tiempo donde termina el intervalo. Entonces este estás aprovechando bien la información de de tu de tu intervalo de medición. Eh, bien. Entonces ¿cómo uno muestra eso? Bueno podés mostrar las probabilidades de los distintos estados, o sea te quedaría la figura la parte más de arriba te queda la corriente. Después te queda la... digamos lógicamente vos con la corriente vos este tenés la eh la media después, la covarianza, la I min, la I bar y ya está. Y ahí tenés todo. Claro vos podrías claro poner después la segunda covarianza. Mm, eso podría ser porque es como que te da una idea de de de cálculo ¿no? De la segunda covarianza, según dice... y después por ahí claro. Esa la siguiente medición. Puf puf puf puf puf puf. Claro. Sí está bien, eso me parece que estaría bien. Una idea de flujo de de de medición. Y puedo hacerlo así o puedo hacerlo con con la distribución de sí... ajá. Sí, está bien."

#### 5. Archivo: PTT-20251219-WA0015.opus

"Bueno empecé a hacer los gráficos, este tuve bastantes problemas porque claro me había un gráfico que me quedaba todo mal todo el tiempo que es el de Macro y R. Hasta que al final me di cuenta que si no le pongo el V corta, es decir la varianza, este la cosa anda muy mal. O sea que hay que incorporar la corrección de varianza sí o sí. Después el otro tema que tengo que ver es que los que digamos cuando el algoritmo es es cero, es decir average igual cero, tiene que que centrarlo en el medio de en el medio ¿cómo se llama? Del intervalo, sino es cualquier cosa. Eh, eso y creo que que con eso estaría más o menos y tendría eh sí. Eh... no sé si básicamente serían creo que serían tres algoritmos. Sería o o tres o cinco ¿no? A ver serían podría ser Macro NR, Macro R, Macro MR, Macro MNR y Macro IR. Y todos con la corrección de varianza, o sea sino no tiene sentido. Puedo hacerlo sin la corrección para mostrar que que el que la cosa anda muy mal, no sé, tengo que ver. Pero sí, podría ser. Este y bueno sí no no lo no ves la diferencia así en un solo en una sola este ¿cómo se dice? Gráfico. Necesito hacer un estudio más serio para mostrar la este que la diferencia de los algoritmos. Esto sería para presentarlos nada más."

#### 6. Archivo: PTT-20251219-WA0016.opus

"Bueno acá el tema entonces es si si re implemento Macro R y Macro NR de manera de que eh predigan el medio del intervalo. Yo creo que sí que tiene que ser eso porque lo otro no tiene mucho sentido. Entonces tengo que hacer esa implementación. Este y bueno. Y la pregunta es qué mierda mostrar. Entonces eh... eh... yo no sé, o sea tendría que mostrar quizás el... sí sería el el prior y posterior como distribución. Este podría poner el prior en términos de este sí de de probabilidad. Después el la la corriente esperada contra la corriente medida. La medida le puedo poner el la ¿cómo se dice? El error de medición. Y este y después con eso sacar la la posterior. Que la posterior este claro. Sí. Sí la pregunta es si pongo el el GS o no. O sea cuánto cuánto del algoritmo expongo en los gráficos ¿no? Eh eso es una buena pregunta. Eh no sé, no me queda muy claro. Eh yo creo que estaría bien exponerlo. Eh... o sea un poco la diferencia es cómo se calcula ese ¿no? El el GS. El gamma digamos eh gamma sigma ¿no? Que es la eh es la proyección de eh sobre el espacio de sobre el espacio de estados de de la dirección de digamos de de la diferencia entre lo medido y lo esperado, de cómo eso lo lo proyecto eso sobre el espacio de de estados ¿no? Eh... sí. Esa sería un poco la la variable. La pregunta es ¿qué es lo que quiero decir con eso? O sea yo lo que tengo que mostrar fundamentalmente es la necesidad de de este nuevo algoritmo. Después este mostrar que que el algoritmo funciona y lo importante es en qué rangos de parámetros funciona y eh y definir las condiciones donde eh los distintos eh algoritmos son este aceptables y condiciones de de aceptabilidad, ¿no es cierto? Había quedado en una cosa de ese estilo. Este yo eso lo tenía bastante hablado y manejado. No sé bien dónde está. Eh tendría que por ahí resumir todo pero medianamente, o sea la idea era esa ¿no? Era que tengo eh yo tenía definidas las figuras, tendría que buscar dónde está eso y centralizarlo todo. Entonces lo que voy a hacer ahora es eso, buscar toda la información que tengo de de cómo quiero armar las figuras y centralizarlo."

---

¿Deseas que extraiga algún punto específico de estas transcripciones?

