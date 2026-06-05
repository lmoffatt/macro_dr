**Archivo:** WhatsApp Ptt 2026-04-11 at 11.39.12.mp3
**Duración:** ~07:30
**Hablante:** Luciano

[00:01] Bueno, ayer fue un día excepcional para Macro IR porque finalmente pude analizar los datos, o sea, un avance real.

[00:09] Todo se lo debo a pensar un poco, porque yo había dispuesto unos gráficos de una manera, digamos, medio atontas y a locas,

[00:21] es decir, mostrando la relación de las variables relevantes,

[00:28] que son una variable que estima el error absoluto,

[00:35] que es el bias por la distorsión,

[00:40] después una medida que tiene que ver con también la distorsión del error,

[00:46] o sea de cómo el algoritmo me engaña en cuanto a su propia resolución.

[00:59] No, sería likelihood, o sea una estimación de la Information Distortion Matrix.

[01:06] Bueno, esas dos son.

[01:08] Ah, y después el error absoluto que se consigue.

[01:10] Son tres variables fundamentales, como la sensibilidad, el error y el hecho de que la corrección del algoritmo.

[01:22] Bueno, me puse a pensar que en realidad todo dependía del número de transiciones que ocurren en un intervalo.

[01:34] Entonces, ahí se unen los tres factores que estaba analizando, que es el número de canales, el largo del intervalo y el ruido, ¿no?

[01:45] Y el ruido, claro, no tendría rola y lo único que aportaría sería el largo del intervalo y el número de canales, ¿no?

[01:50] Porque, bueno, la cinética es constante.

[01:52] Y la cinética está, de alguna manera, representada por el ruido, ¿no?

[01:58] Porque, bueno, si vos tenés una cinética más rápida vas a tener más ruido que si tenés una cinética más lenta.

[02:04] porque el ruido es más o menos instrumental, es constante, no lo podes bajar mucho.

[02:08] Bueno, lo primero que encontré es que el ruido no tiene ningún efecto en ningún lado,

[02:16] y después el efecto es muy diferente en cuanto a la distorsión,

[02:23] es decir, que el gradiente no te dé cero, la distorsión de la media,

[02:30] que esa depende

[02:32] del número de canales

[02:34] o sea, depende

[02:36] digamos del número de transiciones, quiero decir

[02:38] y después

[02:39] el otro factor que es

[02:42] la distorsión

[02:45] de la likelihood

[02:46] increíblemente

[02:48] no depende del número de canales, solo depende

[02:50] del largo

[02:51] del

[02:53] intervalo de

[02:55] la integración

[02:58] sí, o sea

[03:03] da todo como

[03:04] todo y es independiente del número de canales y de ruido

[03:06] o sea, súper simple

[03:07] y se ve claramente que

[03:10] prácticamente

[03:12] la misma curva para

[03:14] los

[03:16] parámetros, digamos los parámetros

[03:18] que son cinéticos y los parámetros

[03:20] que son dimensionales

[03:22] o sea, el número de canales y

[03:24] la conductancia del canal

[03:26] esos cuatro siguen en el mismo patrón

[03:28] y después los que son los parámetros

[03:30] instrumentales que son

[03:32] el ruido instrumental

[03:33] y el bias instrumental

[03:37] es decir, la corriente de base

[03:38] esos no

[03:39] esos no tienen alteración

[03:44] no se comporta de otra manera

[03:47] entonces me queda

[03:50] un panorama muy claro

[03:52] donde aparentemente habría como

[03:57] como este efecto que yo había llamado

[04:00] la inflación en número de muestras

[04:06] o sea como que vos pensás que tenés más muestras

[04:09] de las que tenés, tendrías el doble

[04:12] de las muestras que vos creés que tenés

[04:17] respecto a las que tenés

[04:18] en el caso de la del muestreo con llamada a 100 y bueno y se ve un mínimo en la resolución

[04:30] también entonces ahora lo que me quedaría hacer es hacer la cross correlation temporal

[04:39] de una medición respecto de las mediciones sucesivas y el tema es que claro eso puede

[04:44] una matriz de 600 x 600 x 1000 es una cosa digamos más grande entonces debo ver cómo

[04:53] lo resuelvo eso de una manera más elegante. O sea, podría digamos hacer un vector de 600 nada más

[05:05] o tomar muestras aleatorias a distintos tiempos.

[05:11] Hay muchas palabras.

[05:12] No sé, tengo que ver un poco qué hago con eso.

[05:16] Eso es realmente mi problema ahora a nivel de código y computacional.

[05:25] O sea, cómo resuelvo cómo seguir.

[05:30] porque es como que si lo hago así, de una manera intuitiva pero sería una manera ingenua

[05:43] hacer la matriz de... lo que pasa es que no puedo hacer la cross correlation para todo

[05:52] lo podría hacer

[05:55] lo puedo hacer

[05:56] este

[05:58] pues claro, no es lo mismo al inicio

[06:03] que al final

[06:04] pero no sé

[06:07] tengo que ver bien cómo lo resuelvo

[06:09] es algo que un poco

[06:11] me tiene

[06:12] me tiene medio preocupado

[06:14] o sea, no

[06:16] lo que tengo que resolver es eso

[06:18] cómo hago esta co-crolación

[06:20] de la log likelihood durante un trace, un registro, o sea, las medidas sucesivas, cómo

[06:36] se correlacionan y cómo se correlacionan a los distintos tiempos.

[06:40] Tienes tres variables, tenes el tiempo que estás midiendo y la distancia de ese tiempo.

[06:51] Y bueno, entonces sí, lo pones así, lo podes graficar y te queda bien.

[07:04] Sí, lo que podría hacer es tomar intervalos semi-riolítnicos o algo así,

[07:13] como tomar la inicial contra los primeros 10, después tomar salteadito.

[07:21] Eso sería lo más razonable, yo creo, y tendría más o menos sentido.
