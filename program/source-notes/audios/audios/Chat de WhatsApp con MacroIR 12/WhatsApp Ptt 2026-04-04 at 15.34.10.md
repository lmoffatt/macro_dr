**Archivo:** WhatsApp Ptt 2026-04-04 at 15.34.10.mp3
**Duración:** ~08:52
**Hablante:** Luciano

[00:00] Después de varios días de iteraciones, logré dos cosas, una me costó más que la otra,

[00:13] que es primero hacer como una especie, o sea, expandí mi lenguaje macro IR script

[00:20] De manera que puedo tener ejes, o sea, puedo correr mis diagnósticos de likelihood para varios modelos a la vez,

[00:34] varias y, por ejemplo, varias condiciones de lo que sea a la vez.

[00:39] bien, eso por un lado

[00:41] simplifica mucho

[00:44] el análisis

[00:46] porque bueno

[00:46] puedo probar cosas

[00:50] y verlas enseguida

[00:51] eso por un lado

[00:53] y por otro lado

[00:54] como me está dando

[00:57] algo medio mal

[00:59] en lo que es el gradiente

[01:01] me estaba dando un bias

[01:03] el gradiente no me da cero

[01:06] de algunos parámetros

[01:08] para los algoritmos estándares que estaba probando, habilité la opción nuclear que es el algoritmo Taylor,

[01:18] que hace la aproximación de Taylor de macro IR, que en realidad sería como el Extended Kalman Filter.

[01:29] Y bueno, veo que sí, mejora un poco eso, pero bueno, ahora las tengo que cuantificar bien.

[01:37] Y bueno, y vengo acá un poco porque, bueno, me...

[01:42] Creo que me duele un poco la panza.

[01:44] Tengo la panza medio hinchada.

[01:45] Tanto comer, no sé si ayer.

[01:48] Creo que comí mucho, tomé mucho alcohol y...

[01:52] Tengo miedo con la panza hinchada y...

[01:54] Y eso un poco no me deja pensar bien.

[01:57] Me molesta un poco.

[01:59] Pero bueno.

[02:02] La cuestión es que tengo que terminar este análisis.

[02:05] Mierda.

[02:07] y bueno, eso es muy simple

[02:10] no sé para qué vengo para acá

[02:12] para hacer algo que ya lo sé

[02:14] que es bueno, en lugar de analizarlo

[02:17] la evolución de las cosas

[02:19] la evolución

[02:21] está bien para ver cuando algo falla

[02:24] pero cuando falla

[02:26] sí medio caóticamente

[02:28] no sé si sirve mucho

[02:29] lo que tengo que ver más que nada

[02:32] es

[02:34] cómo

[02:35] el bias

[02:36] reacciona con

[02:40] disminuye con

[02:41] digamos, bueno, con el número de

[02:44] samples sería

[02:46] y también

[02:49] por otro lado con

[02:52] los

[02:54] parámetros que yo ya había dicho antes

[02:56] de lo que son tres dimensiones

[02:59] que eran

[02:59] las tengo que repetir

[03:02] pero bueno, las repito

[03:03] que son

[03:05] la time constant

[03:10] el

[03:12] el level noise

[03:15] y el número de canales

[03:17] y el número de

[03:20] estados

[03:21] tengo que analizar todas esas

[03:23] este

[03:25] variables

[03:27] para ver como

[03:29] como tengo el bias

[03:31] que bias tengo

[03:32] este

[03:34] a ver si

[03:36] en qué punto del bahías es máximo o mínimo lo que sea para cada uno de los algoritmos

[03:46] aparentemente ninguno de los algoritmos es bueno

[03:50] y bueno así que no sé pero con el todo probando ahora con el que sería el bahía

[03:59] es el otro

[04:02] el otro

[04:04] test que es el de

[04:06] el de

[04:07] la varianza

[04:09] el test de la varianza

[04:10] que ese

[04:13] es este

[04:14] no sé si es mejor o peor

[04:17] pero bueno, que se yo

[04:18] ese es un poco

[04:21] más este

[04:23] indica

[04:24] digamos, es un test

[04:27] de lo que sería

[04:28] la métrica de la Lightly Hood. Esto tiene que ver más con la posición, es decir, el

[04:38] bias y esto es con la métrica, o sea, si la métrica de Lightly Hood es confiable o

[04:43] no. Entonces tengo que ver eso. Después me quedaría el tema de la resolución, porque

[04:54] porque vos podes tener dos métodos que sean los dos con métrica y todo bien, pero bueno

[05:00] uno que tenga más resolución que la otra, o sea que uno que digamos asigne variabilidad

[05:07] al parámetro y no al ruido.

[05:11] Me voy a saber un poco la idea, y ver si eso pasa, o sea, bueno eso debería pasar con

[05:23] recursivo, el no recursivo debería ser bastante bueno en cuanto al bias, o sea no debería

[05:29] tener tanto bias, salvo por ahí por el tema de lo que no está pensado, que es para él

[05:36] se llama el ruido instrumental y el número de canales.

[05:58] digamos para eso no debería afectar tanto pero

[06:04] pero bueno yo creo que lo afecta por el tema de que eso está metido en cuanto al promediado

[06:20] de segmentos, quizás tendría que correr también algunas simulaciones sin simulaciones que

[06:31] simulen una precisión instantánea y ver ahí que el método sea perfecto. Eso debería

[06:48] hacer para poder descartar errores que yo no estoy considerando. Lo que vendría que

[06:58] ver ahora es a ver qué fuentes de errores puedo tener. Fundamentalmente una de las fuentes

[07:05] es la simulación. Porque yo estoy simulando al mismo tiempo, estoy testeando mutuamente

[07:14] simulación y la icliput y bueno puede ser que lo que esté fallando sea sea la simulación

[07:20] y ahora verificar un poco a ver si si eso podría hacerlo de una manera mejor quizás

[07:32] simulando lo con la fórmula que ya tengo no sea

[07:43] directamente no sabría que ver eso bien cómo lo podría hacer

[07:50] por ejemplo

[07:53] directamente

[07:56] tomarlo como si fuera un canal instantáneo

[08:02] no, eso no tiene sentido

[08:05] no, lamentablemente me parece que no

[08:10] si no tengo otro otro remedio que simularlo así o sea hacerlo con pasos muy pequeños

[08:17] ahora tengo que verificar un poco cómo es

[08:24] los pasos del algoritmo

[08:33] Porque yo me acuerdo que había...

[08:35] ¿Ves? Ah, sí.

[08:36] Podía haber errores si yo ponía un...

[08:39] Yo tengo un umbral.

[08:41] Ahí es donde jode.

[08:42] Yo tengo umbrales ahí la cosa...

[08:45] Se puede joder.

[08:46] Ahora me tengo que fijar ese tema.

[08:49] El tema de los umbrales.
