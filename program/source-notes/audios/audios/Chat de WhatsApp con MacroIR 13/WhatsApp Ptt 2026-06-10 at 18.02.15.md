**Archivo:** WhatsApp Ptt 2026-06-10 at 18.02.15.mp3
**Duración:** ~12:02
**Hablante:** Luciano

[00:00] Ok, 9 horas 9 de julio.

[00:02] Bueno, ayer fue un día excelente.

[00:06] Solucioné finalmente el bug este que aumentaba la variabilidad de macro IR en 10.000 canales,

[00:21] que era una cosa imposible, encontré la razón.

[00:24] El bug era bastante extraño, no sé si ya lo conté, lo que no.

[00:29] O sea, era bastante sutil, la verdad.

[00:33] O sea, porque lo que pasaba era que primero aparecía un salto en lo que llamaba Joltrass coefficient.

[00:46] Y ese salto se transportaba inmediatamente a la matriz de covarianza, pero no a la matriz de alpinín,

[00:56] es decir, a la matriz media de las prioridades de los estados,

[01:03] sino a la covarianza de esa matriz.

[01:05] Después sí aparecía en la PEMIN.

[01:09] Y era una cosa medio extrañísima.

[01:12] Y claro, agoté todas las posibilidades, o sea, estaba bien tras coeficiente,

[01:19] pero había una cosa que es que si vos tenías un IF,

[01:25] Es decir, trataba diferente los D, digamos, la dirección positiva que la negativa.

[01:35] Entonces, ¿qué significaba?

[01:37] Significaba que si vos estabas en una región donde justo coincidía lo esperado con lo encontrado,

[01:46] el D va a ser cero.

[01:48] Entonces una pequeña variabilidad, por ejemplo, en la corriente de baseline, te va a mover para arriba o para abajo.

[01:57] Entonces vos podés estar calculando en un régimen de positivo la derivada en una condición y en el de negativo en la otra, con lo cual ahí te aparece ese salto.

[02:12] y bueno, efectivamente coincidía

[02:14] justo las samples

[02:16] donde aparecía el salto

[02:19] en la derivada

[02:20] vos encontrabas que

[02:21] el valor

[02:23] de la Imin

[02:26] y la Ifit

[02:27] era el mismo, la diferencia era cero

[02:30] y claro

[02:32] y si eso ocurre, ¿qué es lo que pasaba?

[02:34] que cuando eso ocurre

[02:36] digamos

[02:37] te aparece la derivada

[02:40] del

[02:41] con llamada

[02:43] del

[02:44] transmin

[02:46] coefficient perturbada

[02:50] pero claro, no va a afectar

[02:52] al Pmin

[02:53] porque el Pmin no se corre nada porque la D es cero

[02:56] pero sí te afecta

[02:57] al PCOV porque estábamos

[03:00] usando el mismo alfa

[03:01] para las Pmin

[03:04] y las PCOV, que ese fue el error

[03:06] digamos, o sea, se usaba el mismo porque

[03:07] bueno, supuestamente indicaba no sé qué cosa de la variabilidad y qué sé yo, cosa que te mete la IA.

[03:15] Y bueno, le hice investigar a la IA, porque a mí parecía que era el pedo para macro IR,

[03:22] porque ya los cálculos eran correctos para la corrección de la matriz de covarianza,

[03:31] porque no dependía de la media ni nada

[03:34] y bueno efectivamente

[03:36] no hacía falta hacer la corrección de la PECOV

[03:39] entonces la volví la publicación de la PECOV

[03:41] y ya desapareció el error

[03:43] sí, hace falta una corrección de la PECOV

[03:46] para macro IRT

[03:47] pero macro IRT no está funcionando ahora

[03:50] así que bueno

[03:50] eso cuando lo haga funcionar

[03:52] así que bueno, solucioné el bug

[03:55] era bastante sutil

[03:56] y bueno, ahora ya estoy corriendo todo

[03:59] y viendo a ver lo que da

[04:01] y bueno, da cosas bastante

[04:03] al menos

[04:04] parece robustas

[04:06] bueno, una de las cosas

[04:09] que tengo que mirar con cuidado es

[04:11] claro, los

[04:13] autovalores negativos

[04:15] del

[04:16] gessiano numérico

[04:18] para todos

[04:21] los algoritmos y claro, si vos tenés un montón

[04:23] de gessianos

[04:25] con autovalores negativos

[04:27] aún para macro IR

[04:29] y hasta con 10.000 canales

[04:31] pero bueno

[04:33] si vos los promediás

[04:34] desaparecen casi todos

[04:37] salvo en los algoritmos

[04:39] que no funcionan

[04:40] el macro R para

[04:43] TAU1

[04:44] o el macro NR

[04:46] entonces bueno

[04:48] y no funciona para

[04:50] macro MR para 100 canales

[04:53] entonces

[04:56] bueno

[04:58] Bueno, con eso se cierra el tema de que puedo estimar efectivamente la matriz de distorsión

[05:06] de la información.

[05:07] Y ahora lo que me quedaría hacer, o sea, no voy a aplicar la parte posterior porque

[05:14] la verdad que no sé, no tiene mucho sentido, me parece.

[05:17] Ah, tendría que ver si son singulares, pero...

[05:22] Ah, eso es un directo.

[05:24] Sí, la matriz en Fishery Information, Rauciano, es singular.

[05:33] Si es singular, ahí sí tengo que ser seguro de usar los prior, que sea un postitio.

[05:41] Para el caso de Macro R.

[05:47] pero en principio había decidido no hacerlo, que iba a hacerlo todo "likely"

[05:53] y bueno, entonces yo lo que estaba evaluando era si me convenía o no

[06:01] hacer la matriz

[06:14] hacer encontrar los máximos locales y la verdad es que si todo me apunta que

[06:20] que me conviene porque entonces puedo hacer un test más y la verdad que se

[06:25] lo va a dar bastante robustez a todos es cierto

[06:29] queda un poco la pregunta de si apunto solo macro y r oa todos los algoritmos

[06:41] podría mostrar los casos en que los otros algoritmos fallan, por ejemplo,

[06:47] aunque el hierro falla en algunos puntos también.

[06:50] En realidad lo que habría que mostrar es los puntos donde los algoritmos no fallan,

[06:59] o sea, si los otros no fallan en algún punto, bueno, mostrarlo por ahí.

[07:02] Esa sería un poco la cosa.

[07:06] o ver en realidad, digamos, claro, quizás la regla de oro sería ver si, ahí está,

[07:15] si yo haciendo la corrección de la matriz de covarianza,

[07:20] estimo adecuadamente la covarianza de los parámetros estimados, ¿no?

[07:31] Esa quizás sería un poco la verdadera prueba, porque ya ni eso podés corregir,

[07:39] lo que vos podés corregir no queda mal, entonces bueno, ahí ya sería un algoritmo de poca utilidad,

[07:48] por lo menos para este caso de estimación.

[07:51] O sea, sí, estaríamos haciendo estimación clásica en este caso,

[07:57] no estaríamos haciendo lo bayesiano, es más simple, lo bayesiano lo puedes agregar después.

[08:09] Lo bayesiano lo necesitamos sí o sí para comparar modelos,

[08:13] pero no estaría ahora en la comparación de modelos, sino justamente en la estimación de la variabilidad de los parámetros,

[08:24] el error de los parámetros. Está bien, sí tendría que hacerlo así.

[08:30] Entonces ahí más o menos cerraría.

[08:33] Y bueno, ahí sí podría hacer un estudio de todos los algoritmos,

[08:39] centrado por ahí en eso, en los errores de los parámetros.

[08:52] que es en el fondo lo que podemos encontrar

[08:56] si el parámetro que obtenes es el verdadero, que no es, pues no es un bias

[09:01] y después si la varianza que vos estimas es la verdadera

[09:08] y ver hasta que punto

[09:11] yo creo que con eso estaríamos bien

[09:16] entonces tengo que hacer esa segunda etapa

[09:20] que es la de optimizar y bueno, si, voy a tener que hacer, de optimizar con distintos sizes,

[09:30] o no sé, bueno, al principio vamos a optimizar un solo parámetro de chau, o uno y, si,

[09:41] Sí, o puedo hacer de distintos niveles, sí.

[09:44] Y bueno...

[09:48] Y sí, tengo que poner como ciudadano primero a las Fisher,

[09:56] Fifth Information Matrix de distintos tipos,

[09:59] y no las inversas, porque eso...

[10:02] o las funciones debiladas, porque esa depende de si las funciones se comportan bien o no,

[10:07] y no siempre pasa.

[10:09] O sea que eso...

[10:11] Eso lo tengo que mostrar, o sea, una de las cosas que tengo que mostrar es que el porcentaje de las

[10:19] Digital Information Mentors son positivas definidas, o sea que no tienen muchos valores negativos

[10:29] ni cero, o sea que los resultados son mayores a cero.

[10:35] Sí, eso está bien.

[10:40] Y sí, claramente cuanto más...

[10:45] Porque si vos tenés pocos canales es como que lo que no es el modelo,

[10:51] lo que no es la parte de likelihood, sino que son las otras partes que tienen más información

[10:58] y esas son las que tienden a ser menos normales.

[11:03] Al final la mara de la Krihut es una destrucción normal.

[11:09] Es la idea.

[11:12] Lo que se aparta de la gaussiana son las otras cosas.

[11:16] A medida que tenés pocos canales, lo nuevo gaussiano pesa más.

[11:23] Tengo que ver ahora cómo hacer para que el taper no sea kilométrico.

[11:31] posiblemente pueda ser kilométrica

[11:33] en los suplementarios

[11:37] a todas las tablas

[11:38] y todo eso sí

[11:39] por eso no hay problema

[11:41] y bueno

[11:44] hay que concentrarme

[11:45] en que esta parte de acá

[11:48] cierre preguntas

[11:51] y respuestas

[11:53] o sea, tener preguntas y respuestas

[11:55] de todo, pero creo que tengo que hacer

[11:57] algunos cambios nomás

[11:58] al programa
