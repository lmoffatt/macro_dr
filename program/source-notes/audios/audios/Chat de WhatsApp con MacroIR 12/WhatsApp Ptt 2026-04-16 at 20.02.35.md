**Archivo:** WhatsApp Ptt 2026-04-16 at 20.02.35.mp3
**Duración:** ~23:05
**Hablante:** Luciano

[00:00] bueno bueno hoy fue un buen día logré finalmente el gráfico que estaba buscando

[00:11] y un poco bueno ya toda esta presión que tenía respecto de la cross correlation y todo eso se

[00:23] baja pues bueno ya sé lo que necesito. Bueno entonces en síntesis el gráfico

[00:30] consiste en la evolución de la cross correlation de distintas variables

[00:40] fundamentalmente sería del R, el residuo estandarizado y de la log likelihood

[00:49] para distintos lags, entonces para lag de 1, 2, etc. para distintos algoritmos.

[00:57] Entonces lo que se ve claramente es que el mejor algoritmo tiene, digamos, la correlación decae más rápido.

[01:07] O sea, el de Taylor es mejor que el IR y el IR es mejor que el MR y el NR es espantoso.

[01:17] Bueno, eso un poco es una de las figuras importantes del Pedro,

[01:24] porque muestra claramente, digamos, cómo, digamos, dónde está la eficacia de los algoritmos,

[01:33] que es en eso, en reducir la correlación temporal de este lado de la crejero.

[01:38] Bien, y eso después lo veo también en la Distortion, la Information Distortion Math,

[01:47] matrix y bueno y por otro lado está en paralelo el tema del bias que el bias tiene que ver

[01:56] con el número de canales y importa y aunque también importa el número de canales en esto

[02:01] de la distorsión porque digamos la distorsión ya sea más grande cuando vos tenés este

[02:06] número de canales más pequeño. Bien. Entonces en definitiva

[02:13] ya casi tendría el paper liquidado

[02:20] o sea, tengo que ahora reunir todas las figuras

[02:22] me quedaría por ver el

[02:27] bueno, tengo que, digamos ahora

[02:28] yo creo que hacer

[02:30] 100 y hasta por ahí quizás

[02:33] 1000, o sea, tener

[02:34] 1000 samples

[02:38] por tau, a ver que

[02:40] como da

[02:41] porque es importante ver

[02:44] las distintas

[02:45] cómo se distorsiona

[02:48] la

[02:50] la likelihood con

[02:52] distintos

[02:53] pasos de tiempo respecto

[02:56] del

[02:57] como se llama del tau

[03:00] y eventualmente estaría

[03:02] bueno meter dos taus

[03:04] uno largo y uno corto

[03:05] y ver eso como

[03:07] conviven los dos

[03:10] igual eso podría ser hasta

[03:16] te diría un segundo paper, no sé

[03:18] no sé

[03:23] me parece, o sea, ya como

[03:26] este paper, digamos

[03:27] hasta acá estamos muy bien

[03:29] y habría que ver si hacen falta

[03:32] dos estados o no

[03:33] podría dejarlo en un estado solo

[03:35] y después avanzar con

[03:37] dos estados para

[03:39] otro un paper subsiguiente

[03:43] o no no sé tengo que ver ver eso que hago exactamente por ahí por ahí no

[03:50] porque bueno digamos

[03:53] nada tengo que ver igual es fácil o sea ya lo tengo armado el otro modelo

[04:01] podría hacerlo pero

[04:05] No sé, me parece que este paper cerraría bastante bien con esta caracterización de la likelihood,

[04:13] o sea, de cómo se determina la validez de la likelihood, o sea, de una manera bastante, digamos, rigurosa.

[04:26] Sí, digamos, es rigurosa.

[04:33] me parece que está bien. Lo que no está del todo claro que sirva mucho es esto de la

[04:41] separación de la Distortion Matrix en la Sample Distortion Matrix y la Correlation Distortion Matrix

[04:52] porque claro, porque me queda que la Sample es menor, entonces eso un poco me confunde en

[05:02] en el sentido de, claro, que es como que yo estoy, no sé, me cuesta interpretarlo eso,

[05:10] o sea, me da 07, claro, y bueno, y eso un poco no lo entiendo del todo.

[05:18] O sea, eso es como que me da, a ver qué sería, no sé, lo que sería,

[05:27] la verdad que me cuesta un poco interpretarlo.

[05:32] Es como que ahí la...

[05:35] Claro.

[05:37] Sí, no sé.

[05:41] No sé bien cómo es.

[05:45] Así que la likelihood en realidad es más...

[05:53] O sea, yo la considero...

[05:59] más probable y es menos o algo así, no sé, como que la...

[06:07] algo así.

[06:08] Y eso no me queda bien claro por qué, bueno, sería porque justamente la...

[06:20] la...

[06:22] la...

[06:27] si, no se

[06:31] o sea, lo que pasa es que yo tengo dos likelihoods ahí y me confunde un poco

[06:38] vos tenes la likelihood, o sea tenes la probabilidad de los parámetros

[06:41] la probabilidad de medición

[06:44] y digamos que contribuye a que me

[06:47] me confunde un poco no tendría que pensarlo un poquito más eso

[06:54] si vosotros los parámetros y la prioridad de la medición no son dos grandes 22 espacios

[07:05] prioridades diferentes justamente que digamos a través de la prioridad de medición queremos

[07:09] centrar en tener la realidad de los parámetros. Ese mapa de una probabilidad a la otra es

[07:18] todo el modelado, es la idea esa. Cómo vos proyectas un espacio al otro. Entonces, claro,

[07:32] digamos esta distorsión

[07:34] te daría como que vos

[07:37] digamos que es lo que estaría pasando

[07:43] en

[07:43] en estos mapas

[07:47] o sea

[07:48] eso te tengo que pensar

[07:52] es como que una

[07:54] digamos o sea

[07:55] la función que yo obtengo

[07:58] comparada con la verdadera es diferente

[08:00] y es diferente de qué manera

[08:03] y en qué sentido

[08:06] eso lo tengo que pensar un poco

[08:11] yo creo que todo eso lo tengo que

[08:13] masticar un poco más

[08:15] eso está bien pensarlo así

[08:20] vos tenés dos probabilidades

[08:25] una que va

[08:25] del espacio

[08:27] o sea de un parámetro

[08:30] vos vas a

[08:34] o sea un punto en el espacio de parámetros

[08:37] vos te define

[08:40] ¿no es cierto?

[08:42] una función de probabilidad

[08:46] en el espacio de medición

[08:48] ¿no es cierto?

[08:49] ¿pero qué pasa?

[08:52] que vos en realidad

[08:54] tenés de movida

[08:57] un espacio

[08:58] del espacio de parámetros

[08:59] es un espacio de probabilidad y ese espacio de probabilidad digamos con la función del

[09:03] Alkyhut vos lo proyectas en un espacio en la probabilidad de mediciones y dado que vos

[09:15] tenés un punto de la medición vos ahí eso lo proyectas para atrás en el espacio de parámetros.

[09:27] Entonces esa es la idea de la... ¿cómo se llama?

[09:30] De esto, la deducción, no me acuerdo cómo se llama.

[09:36] Bueno, esto.

[09:38] Esto tiene un nombre que ya me olvidé.

[09:41] Bien.

[09:42] Entonces, yo lo que estoy viendo es...

[09:48] Yo lo que puedo simular con, digamos, cierta precisión es el proceso de...

[09:55] este dado un parámetro, generar digamos mediciones, o sea yo sampleo en el espacio de mediciones y

[10:08] entonces a partir del espacio de mediciones yo digamos tengo la probabilidad, la idea es que

[10:20] tengo dos funciones la función sampling y la función likelihood no es cierto entonces cuál es

[10:29] la idea de todo esto la idea es decir bueno la función de la iclijo de ésta representa a la

[10:37] función de sampling o no no es como puedo yo saber tal cosa bueno yo con la función de sampling y la

[10:50] función de likelihood lo que hago es yo tengo la likelihood y tengo la derivada de likelihood

[10:59] entonces tomo la propiedad de que claro que la esperanza de el gradiente o sea la experiencia

[11:16] la variante del gradiente es igual a la derivada segunda la esperanza de la

[11:22] vida segunda de la likelihood y este y eso de alguna manera te da

[11:30] eso te da digamos una identidad que si no está entonces eso significa que

[11:36] que digamos la la función del likelihood está digamos distorsionada

[11:43] ¿Distorsionada de qué manera? Bueno, probablemente que, por ejemplo, sea como equivalente a un grado de número de libertades mayor o menor, ¿no es cierto?

[12:02] En realidad la distorsión se da para distintos parámetros, los distintos parámetros tienen

[12:14] distintos grados de acoplamiento y bueno, claro, lo que son los parámetros que están

[12:22] relacionados con la cinética que son el número de canales, la conductancia y los parámetros

[12:31] cinéticos esos si están perturbados mientras que los parámetros son más clásicos

[12:38] que hice la corriente de base y el error de la medición y son no están correlacionados

[12:49] y se los estima bien sin problema

[12:57] Bien, el tema es, entonces después lo que puede pasar es, qué pasa si vos tenés dos parámetros

[13:04] con programas cinéticos de distinto, porque como esto está, digamos, como

[13:14] digamos, el error del parámetro depende de la relación entre el intervalo de medición y la

[13:25] constante cinética la constante de entonces 32 constantes cinéticas una va a estar más

[13:32] deformada que la otra y después bueno la deformación de los parámetros no cinéticos

[13:42] como número de canales y la conductancia va a ser la suma de la contribución de ambas

[13:47] parámetros cinéticos

[13:54] bien más o menos y las cosas intermedias y medio quilombo

[14:04] con lo cual cada vos lo que tiene que hacer básicamente es calcular la

[14:12] hacer la calculada historia matrix para calcular

[14:18] e

[14:24] el error de los parámetros corregidas con la distorsión corregida. La distorsión, claro,

[14:33] va a ser un factor aparentemente de 2 no más, o sea, no es un factor terrible, pero bueno,

[14:38] afecta, ¿no es cierto? Es como si tuviera la mitad de los datos.

[14:49] Pero bueno, no sé.

[14:50] Ya más o menos tengo ya una idea bastante definida de cómo es el algoritmo,

[15:03] o sea, en qué circunstancias funciona, en qué circunstancias funciona no tan bien,

[15:09] o sea, que tan bien y tan mal.

[15:15] Y bueno, claro, sí, la idea de que tendría que cada vez, como se dice, corregir esta matriz, es un temita, pero bueno, porque hay que simular todo eso, que se yo.

[15:38] Y bueno, va a haber que hacerlo.

[15:44] Ah, tengo que ver el tema este del número de...

[15:47] Sí, de samples dentro del intervalo de medición también.

[15:52] Es verdad eso, ¿eh?

[15:55] Claro, no, pero igual es...

[15:58] Yo eso creo que no afecta para nada.

[16:01] Voy a probar un par nomás, pero no.

[16:04] Eso no afecta a nada.

[16:10] y si bueno es porque la aproximación de la likelihood por

[16:17] digamos con una normal no es adecuada

[16:23] tendríamos que hacer digamos la distribución microscópica para tener

[16:30] una

[16:33] una deconvolución verdaderamente precisa.

[16:36] Pero bueno, tampoco ese es el objetivo.

[16:41] El objetivo es tener una medición que te permita discriminar cosas.

[16:46] O sea, acá lo jodido en realidad no es que el alcoholismo sea más o menos preciso,

[16:51] sino que no sabemos si dos cosas son diferentes o no.

[16:56] O sea, tenemos que tener una forma de saber, digamos, de estimar la probabilidad de los modelos.

[17:09] porque esto, lo que es malo de esto

[17:12] lo que un poco me pincha el globo

[17:14] es el tema de que

[17:18] la evidencia

[17:25] no es un valor absoluto

[17:30] está deformado

[17:33] está distorsionada la evidencia

[17:36] y va a estar distorsionada

[17:37] en distintos puntos de distinta manera entonces

[17:41] si se hace un poco difícil de validar las cosas entonces al final de cuentas lo lo único que cabe

[17:54] hacer es bueno es

[17:57] suponer bueno vos encontras varios modelos bueno

[18:07] los modelos alternativos que existían los datos y bueno vos obtenes los

[18:12] posterior de los parámetros de los modelos y bueno y ahora vos simulas a

[18:17] partir de esos posterior los resultados y y ves este

[18:23] hay que calcular realmente la

[18:31] O sea, corregís la evidencia para ambos.

[18:34] Lo que se podría hacer, eventualmente, es...

[18:39] digamos, tener un factor de corrección, poner que sea 2, qué sé yo, y chau.

[18:45] Y bueno.

[18:46] Pero, claro, el tema es que vos no sabés

[18:51] si distintos parámetros son deformados de distinta manera.

[19:01] Pero bueno, ahí no sé, ahí habría que ver.

[19:03] Yo creo que de algún modo me parece como que los modelos van a ser afectados por igual.

[19:12] O sea, que a lo sumo te cambia, digamos, te infla la métrica, pero la infla de una manera más o menos pareja.

[19:25] De todos modos, sí, no tenés seguridad, o sea, hay que hacer algún tipo de simulaciones para cada cosa,

[19:38] porque no es ciencia cierta esto, son aproximaciones y tenés que tener manera de cerrar el círculo,

[19:50] es decir, de generar los datos de la simulación, algo así.

[19:57] Así que bueno, eso es más o menos mi cosa acá, porque digamos, o sea,

[20:03] si vos tenés parámetros, los parámetros vos podés encontrarles el rango de...

[20:10] no el rango, digamos, lo que sería el error de cada parámetro, bueno,

[20:14] el error podrá ser mayor o menor, pero bueno, está por ahí.

[20:20] El problema es cuando vos comparás modelos alternativos, ahí como lo ves esa diferencia

[20:29] de error, vos podrías tener que un modelo sea más o menos favorecido por este tipo

[20:38] de errores de ensampleado y qué sé yo.

[20:42] Entonces, eso sería lo que habría que ver como para tener una seguridad más profunda de lo que uno mide.

[20:56] De todos modos, creo que ya está, el algoritmo estaría.

[21:03] y yo creo que la aplicación más polenta ahora sería, que ese sería el siguiente paper,

[21:11] es hacer modelos de pocos estados y tratar de fitear tramos de registros.

[21:22] O sea, hacer como una especie de microscopio de cinética.

[21:27] cinética una especie de así como uno hace la

[21:33] la transformada de Fourier en el tiempo algo así tener como una especie de

[21:39] cinética en el tiempo o sea como que vos tengas vos fites en dos tres cuatro

[21:44] estados con L en cada región y bueno tenés ahí los estados como estarían

[21:50] representados

[21:52] Yo creo que eso sería un buen avance porque sería como transformar un registro en algo

[22:03] un poco más interpretable para un electrofisiólogo.

[22:09] Eso creo que tendría como mucha aplicación porque si vos haces un modelo muy complejo

[22:15] con muchos estados, la gente no lo entiende.

[22:16] Pero uno es un modelo, dos o tres estados y sí, se entiende, digamos.

[22:20] vos tenés un cerrado, dos cerrados, tres cerrados, qué sé yo.

[22:22] No sé.

[22:24] Es algo que lo

[22:26] podés manejar, entonces

[22:29] yo creo que eso

[22:31] estaría bueno hacerlo, ¿ves?

[22:33] Habría que probarlo.

[22:35] Yo creo que eso es lo próximo

[22:39] que voy a hacer.

[22:40] Y bueno,

[22:45] y después, claro, está el otro

[22:47] paper también, que sería

[22:48] Bueno, voy a seguir trabajando en la cinética de P2X bien.

[22:55] Bueno, no sé, estoy un poco harto de todo, pero bueno, ya este Piper ya está prácticamente cocinado.
