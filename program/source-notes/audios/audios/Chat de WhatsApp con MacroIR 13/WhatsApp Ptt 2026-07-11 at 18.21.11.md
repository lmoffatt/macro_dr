**Archivo:** WhatsApp Ptt 2026-07-11 at 18.21.11.mp3
**Duración:** ~27:54
**Hablante:** Luciano

[00:00] Bueno, sin un plan específico, salí a caminar un poco, estaba viendo las figuras, el tema es que, bueno, las figuras finales van a quedar con la Gaussian,

[00:18] como se llama la Gaussian cosa

[00:20] de

[00:21] el Gaussian Fishing Information Matrix

[00:25] en lugar del

[00:27] Numeric Information Matrix

[00:29] entonces

[00:30] pueden cambiar un poquito y como no

[00:32] todavía no corrí todas las

[00:35] las figuras

[00:39] de los otros

[00:41] algoritmos entonces no

[00:42] todavía no puedo terminarlas

[00:44] pero bueno

[00:45] con lo cual digamos que

[00:48] un poco congelado con finalizar las figuras. Entonces pensaba empezar en la escritura,

[00:55] pero claro, el tema es que con la escritura tengo que todavía ver el tema del micro R,

[01:03] micro IR, a ver si los incluyo o no. Entonces vamos a suponer que no los incluyo por un

[01:11] momento entonces creo que está bien con no incluirlos está bien el paper va a quedar bien

[01:19] y bueno entonces

[01:24] mi cerebro se va a derretir pero

[01:29] ya miles de veces he hecho esto y pensado uno y otra vez

[01:40] bueno cuál es el objetivo del paper cuál es el objetivo después de tanto tiempo pensando y pensando

[01:52] ya mi objetivo era

[01:57] curar o sea estar completamente seguro de cómo funcionaba macro y r y de hecho

[02:06] encontré algunas cosas que no estaban del todo bien, o sea que eso estuvo bien.

[02:11] Y bueno, la manera de estar completamente seguro es hacer un test, digamos, exigente,

[02:29] y un test que además muestre que macro IR funciona mal porque sabemos que tiene que funcionar mal

[02:36] entonces si encontramos la manera de que funcione mal en qué región entonces de alguna manera es

[02:43] como que eso te da cierta seguridad de algo o sea si vos un poco ves donde se aparta bueno por lo menos

[02:55] lo tenés agarrado por las astas, o sea, es una idea un tanto, no sé si es una idea muy científica,

[03:07] es una idea un poco emocional, pero bueno, ver un poco donde falla, la idea es que bueno,

[03:14] por lo menos estás detectando, porque tiene que fallar, o sea, es una aproximación, o sea que,

[03:20] Entonces, bueno, sabiendo que estás fallando, bueno, por lo menos, claro, tenés un poco de información acerca de cómo falla y la forma en qué falla tiene que ser coherente, ¿no?

[03:34] Algo así sería un poco la idea. Esto es medio complicado, pero sí, digamos, esa sería un poco la idea, ¿no?

[03:47] O sea, que encontrar que falla y que falla de una manera específica que tiene que ver con el tipo de aproximación que es.

[03:54] Ok.

[03:55] Entonces, ¿qué tipo de aproximación es?

[03:57] Bueno, se aproximan dos cosas.

[03:58] La primera cosa que se aproxima es que es macro, no micro.

[04:07] Entonces que aproximas una distribución multinomial, o sea una aproximación, digamos,

[04:18] sí, multinomial, no sé cómo llamarla a esta altura, microscópica, es una aproximación macroscópica.

[04:30] Bueno, necesito nombres para hacer una macroscópica.

[04:36] Y entonces, bueno, entre otras cosas vos lo que hacés es, en principio, porque la

[04:43] realidad microscópica no es que sea una multinomial, que sería fácil dentro de todo.

[04:49] No, es que es una cosa más complicada que una multinomial, porque es una multinomial

[04:55] con una covarianza arbitraria.

[04:58] Bueno, la media arbitraria puede tener la multinomial, pero la covarianza está un poco

[05:04] fija por la multinomial, entonces no es una multinomial.

[05:07] Y es más, todavía tiene una estructura de covarianza, yo creo que tiene más complejidad

[05:13] que un momento de segundo orden, tiene que tener un momento de todos los órdenes que

[05:19] se te canten probablemente.

[05:22] nada, vos tenés una distribución

[05:26] digamos, un tanto casi

[05:28] arbitraria por llamarlo

[05:29] lo ajustamos

[05:32] con una normal

[05:33] eso bueno, va por

[05:37] el teorema del central del límite

[05:39] y bueno, qué sé yo

[05:40] y eso

[05:41] tiene sentido

[05:45] para valores

[05:48] grandes

[05:49] entonces, al aumentar

[05:51] el número de canales es como que vos vas pasando de una una generalizada en multinomial ponele no

[05:59] sé cómo llamarla que eso es importante de definir cómo llamarla a una normal y por otro lado lo que

[06:06] tenemos es la la aproximación de la likelihood que también de nuevo la likelihood es una función tanto

[06:20] porque depende de cuántos, o sea, la likelihood del intervalo que empieza y termina, o sea,

[06:31] dado una distribución inicial de estados,

[06:44] nos da un estado inicial

[06:49] y un estado final

[06:53] cuál es la distribución esperada

[06:58] de corriente

[07:02] cuál es la corriente media esperada

[07:05] cuál es la distribución de esa corriente media

[07:08] y bueno, de vuelta

[07:10] yo lo aproximo con una distribución

[07:13] una gaussiana, una instrucción normal, pero en realidad es una instrucción compleja,

[07:19] arbitraria. Entonces, también ahí nosotros estamos viendo eso.

[07:25] Entonces, lo que indica eso es que de alguna manera al aumentar...

[07:32] Ah, y la otra cosa es que...

[07:36] que si yo aumento el número de canales, entonces yo voy a estar saltando para una medición

[07:55] de muchos estados, voy a estar haciendo la convolución de muchos estados, entonces la

[08:01] la composición de muchas cosas, generalmente te da una normal, es decir que a medida que

[08:11] yo tengo intervalos de medición más largos y mayor número de canales también se aproximan

[08:19] a la normal. Entonces es como que la idea es que esos dos bordes tienen que ir hacia

[08:29] normal. O sea que lo normal tendría que ser, digamos, a intervalos de medición largos

[08:36] y muchos canales. Eso debería ser lo ideal, ¿no? Y bueno, encuentro exactamente eso,

[08:45] salvo hay una cuestión que es un tanto, si se quiere, paradojal, que es que yo tengo

[08:57] más distorsión para intervalos de medición largos y pocos canales.

[09:04] Que para que si aumenta el número de intervalos, ahí mejora un poco la cosa.

[09:18] Que eso posiblemente sea porque estoy tomando más medición.

[09:25] o sea que al aumentar el número de mediciones también es como que se aproxima más a una normal

[09:34] probablemente sea eso el efecto

[09:38] pero bueno, más o menos lo que encontramos tiene sentido con lo que se espera

[09:49] No trato de hacer una relación cuantitativa porque se escapa un poco a mis capacidades

[10:02] neuronales actuales.

[10:03] Lo dejo como una pregunta abierta, o sea, exactamente cómo plantear una relación cuantitativa,

[10:14] de cuánto se aparta, porque además es como una especie de pozo sin fin,

[10:23] porque si yo encuentro una cosa predictiva de cuánto se aparta,

[10:28] entonces bueno, ahora lo mejoro macro IR y entonces pasaría a ser otro algoritmo.

[10:36] Y ese es un camino que no quiero ir, quiero frenarme a casa,

[10:42] y bueno está aquí entendí el resto quedará para otros investigadores que

[10:50] podría ser yo más adelante o vos, entonces eso está bien

[10:55] y bueno entonces la otra cosa es que bueno que macro IR y macro R coinciden

[11:02] cuando vos tenés 10 canales

[11:05] y cómo se dice

[11:12] intervalos de medición muy cortos lo cual es exactamente lo esperable porque si vos tenés

[11:18] intervalos de medición muy cortos entonces básicamente

[11:24] no hay mucha diferencia, vos no esperas que haya casi cambios en la conductancia durante una medición

[11:40] y qué más pero que justamente esa es la zona donde macro y r es más macho tonosa

[11:58] digamos después me queda ese tema no sea el ranking de los algoritmos bueno me queda que

[12:06] que hay un único ganador que es macro IR porque macro IR funciona donde macro IR funciona

[12:14] mal y funciona igual de mal que macro IR, bueno por definición, o sea que macro IR sería

[12:21] entre comillas un mal algoritmo comparado con macro IR y macro IR ni te cuento pero

[12:31] Bueno, macro NMR podría tener lugar en

[12:35] macro NMR

[12:39] en circunstancias donde la velocidad de acomplizar el modelo impida

[12:49] resolver digamos

[12:54] matices

[12:57] habría que ver un poco si la velocidad de macro porque macro en el real no ser

[13:09] recursivo digamos bueno podría ser altamente prelizable

[13:13] tendría suyo pero bueno

[13:20] Claro, la pregunta es, vos podrías usar macro...

[13:24] Claro, ahí está.

[13:27] Porque macro...

[13:29] O sea, una de las cosas que encuentro, a ver, es la siguiente.

[13:37] Que aumentar la resolución temporal no te otorga mayor resolución en las constantes cinéticas.

[13:50] Sí te la aporta en lo que es el número de canales y la conductancia.

[13:58] Para eso es indudable que tenés mayor resolución teniendo más mediciones.

[14:05] Eso no cabe duda, que sí.

[14:08] Pero no, en cuanto a las constantes cinéticas, no.

[14:13] básicamente más que dos veces la constante cinética no gana mucha resolución.

[14:22] Con lo cual, si uno quisiera tener un sistema rápido,

[14:28] quizás habría que ver de proponer macro MNR con grueso,

[14:37] porque ahí debería ser claro.

[14:40] en ese régimen funcione mejor que el macro MR. Tengo que ver bien ahí cómo es la diferencia

[14:50] entre esos dos, porque ese podría ser un algoritmo posible.

[14:58] Básicamente la conclusión es que hay que usar macro IR, otros algoritmos son inútiles.

[15:07] El Macro MR es un "Stromant" prácticamente, es una idea para mostrar que una idea simple no funciona.

[15:18] La verdad, no sé, podrían incluirlo.

[15:25] Quizás, no sé, la verdad que no lo sé.

[15:31] Podrían incluirlo.

[15:36] bien entonces

[15:39] qué me queda a ver

[15:43] o sea yo lo que dije es bueno primero eso que yo quería mostrar dónde falla y falla donde debe

[15:54] fallar y lo encontré que falla está bien en ese sentido igual la cantidad que falla es manejable

[16:02] para ciertos usos usar

[16:06] no sería invalidante

[16:08] después

[16:18] qué más

[16:36] la resolución de los digamos de los

[16:42] cosas también y bueno que podés ganar digamos si dividiendo o sea puede sacarlo un poco más el

[16:52] jugo teniendo midiendo varias corrientes de 10 canales de una sola de 100 con él pero tampoco es

[17:05] guay guay cuánto que sacas porque tenés un esfuerzo computacional importante

[17:10] que en principio digamos este

[17:15] medir todo junto digamos es más eficiente por menos computacionalmente

[17:24] igual esto es solo para condiciones

[17:31] digamos, no estacionarias, para condiciones estacionarias, ahí, digamos, claramente si vos promediás te quedás constante, o sea que no aplica esta idea.

[17:52] bueno está el tema este, una de las cosas que más me sorprendió es cómo

[17:59] como la

[18:02] si la matriz de Fisher

[18:08] registra el hecho de que

[18:13] vos no extraes más información acerca del número de canales por ejemplo

[18:19] luego de que cesa

[18:23] una vez que deja de subir el número de canales abiertos, cuando va bajando

[18:32] es como que no puedes ganar más información acerca del número de canales originales

[18:39] eso a mí un poco me huele a la cabeza

[18:43] un resultado muy fuerte, muy claro y habla un poco de la belleza o no sé cómo llamarlo,

[18:58] la agudeza de Néstor, la matriz de Fischer y los algoritmos recursivos, o sea que muestra

[19:08] que desaparece, se vuelve en cero

[19:14] la magnitud de la matriz de Fischer

[19:20] cuando resulta en el tiempo

[19:22] cuando llegas, pero vos no ganás más información

[19:27] acerca del número de canales originales

[19:31] bueno, sabés, va decayendo

[19:34] pero cuántos había originalmente

[19:37] dado que sabe que hay tantos ahora, ese es el tema, ese es el truco, dado que sabe que

[19:49] hablas, o sea no ganas más información. Eso me pareció como algo a mí me voló la

[20:03] la cabeza así que voy a escribirlo para el paper no creo que vaya el abstract porque

[20:09] es un poco demasiado

[20:11] y bueno

[20:15] y bueno entonces macro IR es un algoritmo que recomiendo para usar para corrientes macroscópicas

[20:26] o sea de todo tipo

[20:27] El tema obviamente es que yo lo probé con dos estados nada más y puedo con esto recomendarlo

[20:36] para todo tipo de estados y qué sé yo.

[20:38] Y la verdad que es un poco mucho.

[20:41] Pero bueno, por lo menos sé que ese, o sea, sé que los otros no funcionan, ¿no?

[20:51] O sea, este funciona con dos estados, otros no.

[20:54] si funciona con más estados

[20:58] en principio no lo puedo decir

[21:01] pero

[21:01] no sospecharía que sí

[21:04] la pregunta es

[21:06] ¿por qué no hago una prueba

[21:08] con una cosa compleja

[21:11] de muchos estados?

[21:12] la verdad

[21:15] es que quiero hacer yo un estudio

[21:17] un poquito más detallado

[21:19] de ver

[21:20] cómo unen

[21:22] digamos dos estados uno por ejemplo tengo

[21:27] cada un estado fijo y voy cambiando el otro

[21:33] ver eso como lo afecta yo quisiera hacer un estudio así

[21:41] digamos el principio digamos el hecho de que yo haya obtenido

[21:47] los resultados de comunicarios biologiques te hacen suponer que sí que

[21:52] podés usarlo tranquilamente pero bueno acá acá está un poco

[22:00] una zona de peligro en el sentido de bueno pueden pedir más

[22:07] más cosas o no yo la verdad prefiero hacer

[22:11] como se dice preciso con dos canales llegar digamos muy firme a esto y después

[22:18] es seguir con más que tirar una cosa así a la marchanta de muchos canales y decir, bueno,

[22:23] sí, funciona.

[22:24] Está bien.

[22:25] Sí, bueno.

[22:26] Lo que puedo decir es que, bueno, probablemente se descompongan, o sea, mi hipótesis es

[22:39] que, digamos, este, ¿cómo sería?

[22:46] el tema es

[22:48] si vos tenés

[22:49] constantes

[22:51] cinéticas de distinto orden

[22:54] conviviendo

[22:56] bueno

[22:58] ahí un poco

[22:59] el esquema este

[23:02] semilogarín

[23:05] no sé cómo llamarlo

[23:06] que vos haces

[23:08] o intervalos exponencialmente

[23:10] que se incrementen exponencialmente

[23:14] creo que es la salida lógica

[23:18] porque

[23:19] vos

[23:22] vas teniendo

[23:25] el mayor

[23:28] peso resolutivo

[23:30] a ver cómo sería

[23:35] vos le das

[23:38] digamos

[23:39] este rango dinámico

[23:40] a todos los

[23:43] con los más cortos resolves lo más rápido los más largos lentos y todos con la misma cantidad de

[23:51] puntos básicamente entonces es como un sistema un tanto democrático en el sentido de cuánto

[24:00] cuánto digamos este punto en la cina sacado ese tema lo que pasa es que yo lo podría simular

[24:08] tranquilamente y ver que es así y no lo estoy haciendo.

[24:13] Entonces, claro, yo tengo acá varias puntas para seguir.

[24:17] Entonces lo que voy a hacer es simplemente las voy a decir

[24:21] y nada, que harán para trabajar con los futuros.

[24:24] No puedo hacer todo en un solo trabajo.

[24:26] Es inhumano y está bien.

[24:34] Yo me quedo en esto tranquilo.

[24:37] bueno digo si es todo trabajo si un revivir me pide bueno lo hago pero en principio bueno

[24:46] entonces este eso sería eso sería todo o sea el mensaje más o menos es bastante simple

[24:54] que es que bueno macro y r es el único que sobrevivió a esto y sobrevive en esas

[25:03] circunstancias y en las otras más o menos no está tan mal, digamos, es una cosa razonable.

[25:10] Sí, que hice un sistema minimalista para estudiar las cosas con claridad, para entender bien,

[25:28] para separar lo que es estas variables, las tres o cuatro variables importantes

[25:35] que son el número de canales, el nivel de ruido y el largo intervalo de integración

[25:43] que el largo intervalo de integración es lo mismo que está expresado en términos de la constante sinélica.

[25:57] No trabajé con la prioridad de apertura, porque prefería trabajar con esas otras variables.

[26:13] Sí, podría haber trabajado, podría ser que una opening mucho más bajo, por ejemplo,

[26:22] y ver ahí la resolución del CAON sería más pobre, probablemente.

[26:30] O digamos, no sé, sí, más...

[26:35] Y eso se ve un poco porque, por ejemplo, el CAO...

[26:43] que el caón está un poco...

[26:48] a ver cómo es el tema...

[26:50] Claro, el caón que solamente está conviviendo con el caos,

[26:58] es decir, que ahí el sistema tiene una lambda, un autovalor que es el doble,

[27:03] porque suman los dos, claro, te queda mejor resuelto, ¿no?

[27:10] no sé, no es que mejor resuelto, si no, no tenés tanta distorsión, tenés menos distorsión,

[27:20] lo cual hace suponer que eso.

[27:23] Pero sí, no sé, sí, tendría que estudiar, qué sé yo, eso, tendría que estudiar miles de cosas,

[27:34] pero... no sé... es un poco mucho, o sea, ya está.

[27:43] Quise estudiar esto bien y dije que ya está.

[27:47] Bueno, ya está.
