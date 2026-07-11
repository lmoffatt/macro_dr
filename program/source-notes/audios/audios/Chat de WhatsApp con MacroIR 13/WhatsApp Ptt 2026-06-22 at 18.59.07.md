**Archivo:** WhatsApp Ptt 2026-06-22 at 18.59.07.mp3
**Duración:** ~24:41
**Hablante:** Luciano

[00:00] Bueno, de vuelta a la 9 de julio ya es tarde, no sé.

[00:07] Vamos a jugar a Argentina, Messi metió dos goles, lleva cinco goles en el Mundial en dos partidos de Liria, con 39 años.

[00:17] En fin, estamos en momentos muy raros del universo.

[00:22] Y yo con este problema, que no termino de cerrarlo, pero creo que hoy hice un gran avance,

[00:31] porque yo me imagino así como un gráfico, como dos regiones,

[00:37] la región multinomial con pocos canales, la región paizoniana con intervalos muy cortos

[00:47] y después la zona gaussiana, que sería con muchos canales e intervalos no tan cortos

[00:54] y con más ruido instrumental.

[00:57] Yo creo que un poco esa es la definición del paper SS,

[01:06] que es esa región donde el ruido instrumental y los intervalos de medición

[01:14] son relativamente no tan cortos y bueno, o sea, como que lo podes considerar como un

[01:20] proceso gaussiano, o sea, como que el reloj es gaussiano, ¿no?

[01:24] En la determinación de la conductancia.

[01:28] Es el lugar donde este algoritmo es ideal, pero igualmente aún, digamos, en la región

[01:37] más poesioniana tampoco es que, digamos, la distorsión es muy grande, ¿no?

[01:43] no es una distorsión que es medible pero es pequeña

[01:50] no llega a ser factor 1 siquiera

[01:53] factor 1.2, 1.3

[01:59] pero los otros algoritmos no fallan garrafalmente

[02:08] Tengo que ver cómo se comportan los otros algoritmos con regímenes de mayor ruido instrumental.

[02:18] O sea, porque acá estoy usando un ruido instrumental que es 0,1 veces el kinetic rate.

[02:28] Estoy tomando, digamos, el error, el ruido instrumental es del orden de 0,1 veces la conductancia de un canal

[02:46] por el tiempo de vida de una transición.

[02:56] O sea que estoy considerando regímenes donde las señales del canal único son bien claras.

[03:05] Y eso puede no ser la situación.

[03:12] Vos podés tener cinéticas donde las aperturas y clausuras del canal sean más rápidas.

[03:19] Que no se pueda resolver fácilmente un canal abierto.

[03:25] Pero igual es una situación posible, el 0,1 es imposible.

[03:32] Con buenos sellos y buenas mediciones lo podrías obtener.

[03:39] Además que para ciertas constantes cinéticas están dentro de ese rango.

[03:47] Otras estarán más corridas hacia la derecha.

[03:53] que son del orden del ruido instrumental.

[04:04] Eso por un lado.

[04:08] Lo que quiero decir es que el método es bueno, no es perfecto, pero está muy confiable.

[04:22] confiable, mientras que los otros dejan bastante que desear, igual el macro R más o menos

[04:29] zafaría, digamos, es del orden de dos veces, que se yo, o sea, tiene un ruido bastante

[04:35] más grande, pero ya el otro, el no recursivo, no es, es muy malo, o sea, cientos de veces

[04:42] la inflación de la varianza es muy grande, la distorsión de la varianza es muy grande.

[04:52] O sea, no podés usarlo sin corrección, ¿no?

[04:57] Bien.

[04:59] Bueno, eso sería la parte esta de... ¿cómo se llama?

[05:10] de análisis

[05:13] entonces la pregunta es cuál es el mensaje del paper, el mensaje del paper es bueno yo

[05:20] presento el nuevo método de macro IR, qué problema resuelve

[05:27] digamos ilustro el problema, lo ilustro con estas mediciones, para

[05:38] el problema si tengo que mostrar cómo falla macro R y cómo fallaba macro NR y entonces

[05:47] después propongo una solución que es macro MR, que es la solución simple, que es considerar

[05:55] "open channel" o sea como que el canal tiene ruido y después consideras la siguiente

[06:12] solución que es el tal de macro IR. Claro ahora entiendo la diferencia entre macro IR

[06:19] macro MR es que macro R ya está centrado en el intervalo y macro MR está al comienzo

[06:26] el intervalo. Así que yo creo que un poco por eso es que anda tan mal también, ¿no?

[06:31] O sea, yo creo que en parte de eso es ser lo que lo jode. O sea, posiblemente, digamos,

[06:39] macro MR puesto centrado en el intervalo ande un poco mejor. Pero bueno, no tengo ganas

[06:45] hacer todas las cuentas

[06:48] con macro IR ya está una curiosidad seria

[06:55] entonces a ver cómo presento el paper esa es la pregunta

[06:59] la respuesta es bueno, tengo que preguntar cuál es el problema que se trata de resolver

[07:06] cuantificar el problema porque es una sospecha de un problema

[07:10] digamos, ya uno sabe que eso es así, pero bueno, vamos a cuantificar para ver si eso lo resuelvo.

[07:18] O sea, puedo tener una métrica de resolución del problema.

[07:22] Y bueno, entonces después lo resuelvo, o sea, propongo el nuevo método de macro IR.

[07:32] La RA ya estaba descrito ahí, pero bueno, lo presento de vuelta.

[07:39] Y lo que hago es muestro cómo el macro BIR soluciona esos problemas,

[07:46] pero no los soluciona del todo, es decir, que tiene un régimen donde anda bien

[07:50] y un régimen donde ya no anda tan bien.

[07:54] Y creo que ahí cierro el paper, o sea, no hay mucha más historia.

[08:06] La pregunta es, claro, yo, por ejemplo, que podría ser...

[08:11] O sea, a ver...

[08:13] O sea, yo tengo todos estos esquemas, ¿no?

[08:17] De, bueno, el experimento, digamos, simplificado y todo eso.

[08:23] Y que varío, digamos, la señal a ruido, varío el largo del intervalo, varío el número de canales.

[08:31] Son tres variables, o sea, es bastante extenso.

[08:36] con muchas simulaciones. La pregunta es, bueno, si lo muestro con todos los algoritmos o no.

[08:42] Yo creo que no está mal mostrar en todos los algoritmos,

[08:49] porque, bueno, el trabajo está hecho y, bueno, qué sé yo, o sea, ilustra un poco.

[08:56] Especialmente, bueno, yo ilustro con todo este desarrollo que hice.

[09:06] la distorsión de correlación, la distorsión geométrica o de sample,

[09:12] tengo que ver cuáles dos uso, y la distorsión de la información, que es la suma de las dos,

[09:16] o sea, la distorsión de la información es la que te la ves representada después en la

[09:24] en la correlación

[09:27] ¿cómo se llama?

[09:28] en la

[09:28] en la modificación

[09:34] de la matriz de covarianza

[09:37] la matriz de covarianza

[09:38] reducida e empírica

[09:40] esta la diferencia entre los dos

[09:43] está esta matriz de distorsión

[09:45] y bueno

[09:47] y eso te da una idea

[09:48] de qué tan buena es la aproximación

[09:50] de la Liklitz y cómo se puede

[09:53] mejorar

[09:54] y bueno

[10:00] entonces

[10:03] si volvemos al tema

[10:07] se plantea cuál es el problema

[10:09] el problema de sistemas

[10:12] estocásticos

[10:13] de sistemas estocásticos

[10:16] que modelan

[10:18] procesos biológicos

[10:20] como los canales

[10:23] pero pueden moderar muchos otros

[10:25] y como

[10:29] está el tema

[10:31] de que son idealmente

[10:33] instantáneos, pero la verdad

[10:35] es que siempre son integrados

[10:36] y que bueno, que no es lo mismo

[10:40] una cosa que la otra

[10:41] entonces

[10:41] uno puede decir, bueno, pero no importa

[10:45] pero no, bueno

[10:46] veamos las distribuciones

[10:49] o sea, hagamos unas pruebas

[10:51] muestro que dan mal

[10:53] digamos si vos no usás

[11:00] no recursivo

[11:02] o usás recursivo pero no

[11:04] no integrás en el tiempo

[11:07] te da mal

[11:08] entonces ahí se propone

[11:11] el nuevo modelo que supera eso

[11:13] y bueno

[11:16] y muestra que anda bien

[11:17] y cuáles son los límites

[11:19] y ahí después se dijo algo en la discusión

[11:23] acerca de cómo

[11:25] el paper

[11:27] digamos, cómo se podría

[11:29] corregir la evidencia

[11:31] y la likelihood de la evidencia

[11:33] y bueno, y listo

[11:35] y ya tengo

[11:37] digamos

[11:38] el macro IR

[11:41] puesto como el

[11:43] no es un gold standard, pero bueno

[11:45] sí, es un método

[11:47] de batalla

[11:48] práctico para procesos marcovianos integrados en el tiempo. Y que efectivamente cualquier

[12:01] proceso va a tener que ser usado por eso. O sea, no podés usar macro IR en ninguna circunstancia,

[12:08] no sirve, esa sirve, pero el macro IR es mejor.

[12:14] En cuanto, claro, ojo, es mejor en la estimación de la matriz de covarianzas

[12:23] y posiblemente en la comparación de modelos,

[12:25] pero puede ser que igual de otro más o menos sirva para determinar los parámetros,

[12:31] porque al final de cuentas este encuentro que es más o menos parecido,

[12:36] o la determinación de parámetros, no se nos parecía. Lo cual explicaría por qué Macro IR no ha sido desarrollado antes.

[12:45] Este problema no se lo atacó porque es un poco sutil.

[12:50] Ahí está, entonces hay un poco cierro que genera eso.

[12:55] Y bueno, eso sería todo el paper.

[13:02] bueno digo que ya fue usado para canalizión y pós etc

[13:10] está probado que sirve

[13:13] bueno ya está

[13:15] no hay mucho más nitida digamos

[13:19] entonces qué son los gráficos las figuras bueno tengo una figura que sería

[13:24] mostrar cómo funciona

[13:26] en realidad digamos yo empezaría ahora al revés

[13:29] saber

[13:32] yo estaría usando

[13:38] si son seis algoritmos los que estaría probando

[13:47] serían nr, nmr, r, mr, eir

[13:56] son cinco algoritmos

[14:06] bueno

[14:12] la pregunta es porque son todos esos gráficos son complejos son grandes pero bueno yo lo que

[14:25] mostraría el bias es uno y el otro es bueno todas las matrices de

[14:35] distorsión y la covarianza y después bueno podría mostrar también cómo se

[14:47] corrige la corrección de la covarianza

[14:51] eso está bien, es un linto gráfico y podría poner algún gráfico ilustrativo de los parámetros

[15:08] estimados con la covarianza estimada.

[15:14] creo que esos serían todos los gráficos

[15:19] sí

[15:28] entonces la pregunta es ¿hacen falta más gráficos? ¿hace falta por ejemplo poner

[15:32] yo la Lightly Hood, podrás poner la

[15:37] el RSTD, los residuos también los residuos yo creo que

[15:43] estaría bien porque eso también me da, eso es diagnóstico en cuanto a que fallan en otros,

[15:50] entonces también lo tendría que poner, o sea, los residuos sería, claro, sería el primero,

[15:57] el primer test, digamos, el segundo es el gradiente y el tercero es la matriz de covarianzas,

[16:07] Esa sería el primero es I, los residuos, los R cuadrados, o sea, los residuos, la normalidad

[16:18] de los residuos, y el primer test, el segundo test es la esperanza del gradiente, y después

[16:26] la otra es la covarianza del gradiente

[16:32] pero ya comparada con la derivación

[16:41] los tres test

[16:45] claro, vos decís, bueno, ¿cómo?

[16:52] yo comparo para todos los parámetros

[16:55] está bien eso comparado con los parámetros es un poco más pesado pero pero puedo hacerlo

[17:02] tengo la diferencia entre el ON y el OFF

[17:07] el ON tiene un doble, está expuesto al doble de RAID

[17:14] porque son la suma de los dos

[17:16] el otro es la mitad entonces eso es picaría

[17:24] que esté corrida la curva y después tengo los otros que son la visión especular

[17:31] el número de canales y la conductancia

[17:38] el producto de los dos tiene que dar la corriente

[17:46] entonces vos definiendo uno ya definís el otro

[17:51] ¿Y cómo definís el 1? Bueno, por el ruido, ¿no? O sea, claro, sí.

[18:00] No sé cómo lo defino, no sé, a ver.

[18:04] Pero bueno, el tema es que sí, está afectado por la...

[18:19] bueno por el hecho de que puedas tener

[18:24] digamos canales que quedan congelados y no se abren

[18:28] entonces eso lo ves como si hubiera

[18:34] menos canales creo, no sé, yo voy a pensar un poco eso a ver cómo

[18:39] cómo explico con los dedos digamos las direcciones de las curvas

[18:45] entonces la idea es como que

[18:49] que tengo digamos todos estos diagnósticos

[19:00] que son muy bien diagnósticos

[19:11] yo no mostraría más que eso, no tiene sentido mostrar otras cosas

[19:17] y bueno, con eso los caracterizo a los distintos algoritmos

[19:26] o sea, caracterización de los algoritmos

[19:29] de acuerdo a los distintos diagnósticos

[19:37] y bueno y eso sería todo el piper ya estoy exhausto con esto la verdad que sí ya estoy exhausto

[19:45] así que bueno preguntas que me quedan serían bueno

[19:52] le meto más canales por ejemplo meto canales de 1, 1, 2, 5, 20, 50 por ejemplo canales o no

[20:07] podría hacerlo para Macro IR, para nosotros me parece que no tiene sentido, seguro que no tiene sentido

[20:18] para ahí sí tendría

[20:24] para ver un poco bien en qué más definido es eso

[20:33] Pero igual, la verdad es que no sé mucho porque...

[20:38] O sea, quizás, donde es necesario entre 100 y 10 canales, sí.

[20:45] Quizás tenga que hacer.

[20:47] Entre 10 tengo que hacer 100, 50, 20.

[20:50] Si 50 y 20 sí me parece que son necesarios.

[20:54] Y además los compenso con otros.

[21:00] Así que eso yo creo que los voy a hacer.

[21:03] 50 y 20 son dos, es nada más

[21:06] no puedo hacer, no es caro

[21:09] si, eso lo voy a hacer

[21:13] por lo menos para hacer 0,1 de ruido

[21:18] podría hacer menos ruido también

[21:23] a ver si eso empeora

[21:33] podría ser, a ver qué pasa

[21:35] es mío, dedo por ahí pero no sé

[21:41] pero sí, podría ser, mucho menos ruido

[21:47] qué más podría hacer, bueno después tengo la parte de si quiero hacer

[21:59] digamos

[22:02] estacionario no veo aquí que sea súper necesario lo estacionario porque me parece que con esto ya

[22:11] ya digamos lo demuestro lo que es lo estacionario serviría para mostrar que

[22:18] bueno que macro y r y más resolución que macro r pero

[22:24] La verdad que, digamos, o sea, macro viernes es mejor, no sé si hace falta.

[22:32] Podría ser después si es necesario. Creo que no está bien.

[22:38] Estacionario no hace falta.

[22:41] Sí, quizás lo de resolver bien lo que sería así.

[22:48] O sea, tengo que tener resolución más fina hacia abajo y hacia arriba.

[22:54] y con eso quizás pueda definir bien las curvas estas

[23:00] y encontrar digamos una explicación que más o menos estaría

[23:08] el ruido poesoniano que aumentaría con el tiempo

[23:17] y el ruido gaussiano que aumenta con la frecuencia

[23:20] entonces bueno, uno compensa al otro y después queda al máximo y bueno, después el otro

[23:27] es la correlación por...

[23:31] porque bueno, el canal no se cierra, o sea que queda congelado, entonces ahí tenés la

[23:45] la distorsión por correlación

[23:47] practicada por eso

[23:50] la verdad es que no hace falta hacer más

[23:58] pero bueno, lo voy a hacer igual, total

[24:00] no pierdo nada

[24:01] diciéndolos, es tiempo que está ahí

[24:05] creo que está bien

[24:12] o sea, más o menos en lo que sería

[24:17] el mensaje del paper está escrito ya

[24:21] es bastante claro

[24:23] creo que ya lo tengo

[24:26] no sé si es un paper para Elias, la verdad

[24:35] quizás no lo sea

[24:36] pero bueno, es lo que es
