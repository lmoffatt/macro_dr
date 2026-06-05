**Archivo:** WhatsApp Ptt 2026-04-17 at 14.55.48.mp3
**Duración:** ~14:48
**Hablante:** Luciano

[00:00] Pero bueno, supongamos que queremos analizar las contribuciones de cada medición a esto.

[00:10] Claro, bueno, después sí, tenés que ir a nivel de, sí, para cada medida, digamos, como es la correlación temporal.

[00:24] Ahí sí, no está bien.

[00:28] Claro, ahí me queda bien claro que bueno, no hay correlación cuando es puro ruido

[00:38] y después eso aumenta y baja.

[00:44] y bueno, seguramente

[00:46] se hace cero

[00:49] para

[00:49] para el, ¿cómo se llama?

[00:52] el, este

[00:56] sí, el off

[01:00] el on, por ejemplo, el on tendría que hacerse cero

[01:02] el off

[01:03] eso estaría bien

[01:10] claro, vos lo que tenés

[01:17] si con

[01:18] con la Fisher Information

[01:21] Metrics es cuánto

[01:23] contribuye, claro

[01:24] ahí sí, ahí sí es

[01:27] acá la Fisher

[01:29] justamente

[01:31] eso es una cosa que tiene

[01:32] pero no te cubre todo

[01:35] justamente

[01:39] Y ahí el quid, ¿no?

[01:43] Ah.

[01:45] ¿Por qué es así?

[01:56] Ahí hay algo raro.

[02:00] Bien.

[02:09] Bueno, entonces, nada.

[02:16] Sí, no sé, la Fisher Information Matrix, la verdad, que no tiene sentido guardarla

[02:29] por cada...

[02:30] No.

[02:31] Eso no aporta nada.

[02:37] Pues la...

[02:39] ¿Qué sería?

[02:41] Claro, sí lo que tengo es eso, ¿no?

[02:54] Y eso sí.

[02:55] Eso sería la...

[02:58] Eso es for...

[03:07] Por cada edición, digamos, cómo es la matriz de cross-correlación temporal, ¿no?

[03:25] Eso sí, eso lo puedo mostrar y sería, bueno, claro, digamos, los tiempos de OLAGS

[03:37] cross correlation, eso sí. Claro, bueno es el lag y también es la

[03:45] Kavlag, sería uno solo la matriz y el otro serían por lo menos 10 con L.

[03:53] Que creo que sí, que bueno para 100 se va a aguantar.

[04:04] veremos para más que pasa no para mí cuántas tienen que ser pero está bien eso es lo que

[04:14] tengo que hacer fundamentalmente es ahora hacer que él en la decir es si el momento

[04:27] tenga digamos un lag límite, vos le pones un lag máximo.

[04:35] Que también podría ser hasta programable.

[04:43] Vos lo calculas hasta que llegas, pero bueno a ver cómo es. Claro si lo vas calculando

[04:56] hasta que llegues a un punto, no sé.

[05:00] Luego lo que te queda es...

[05:06] que te queda una matriz y bueno...

[05:10] le metes el...

[05:14] el determinante de la matriz, no sé.

[05:18] Está bien, ahí lo hago con...

[05:22] hasta que un determinante de la matriz sea tanto.

[05:24] está bien creo que estaría bien es un poco peligroso mejor ponerle un número fijo y chao

[05:35] lo pongo por no por si dinámicamente no sea por comando totalmente lo voy a hacer así

[05:48] ese es el cambio, o sea que la Siri es Momentum

[05:55] te calculo los Momentum, las Cross Correlations

[06:03] con un lag máximo, que el lag máximo lo pones vos,

[06:11] lo manejas como un dato que se le entrega.

[06:17] Está bien.

[06:19] Ese es uno de los cambios, el otro cambio van a ser los determinantes.

[06:25] Y después el otro cambio es, obviamente, los K efectivos, ¿no?

[06:32] la multiplicación

[06:37] de los

[06:38] en efectivos

[06:40] el K efectivo

[06:41] el K efectivo sería como

[06:50] multiplicas K por el

[06:52] y ahí sí sería

[06:55] el K efectivo por el número de

[06:57] parámetros

[06:58] esa es la

[07:01] digamos si todo feria

[07:02] más o menos equivalente, esa es la, ¿cómo se llama?

[07:08] El factor que le va a la evidencia, ¿no?

[07:19] Porque creo que sí, que justamente en la evidencia vos le...

[07:24] Creo que hay una, creo que es la de...

[07:29] No me acuerdo como se llama las mediciones de los modelos, que se yo, que justamente

[07:38] le arrestan el número de parámetros de IO2, creo que es una cosa así.

[07:42] Entonces en este caso sería eso más o menos.

[07:46] Sí, está bien.

[07:51] Bien.

[07:53] Sí, con todo esto creo que lo podemos ir sacando y la idea sería, bueno, primero correr,

[08:03] no sé yo. Voy corriendo por trancos de, yo creo que lo mejor es eso, por trancos del

[08:12] El intervalo de medición.

[08:14] Eso es lo mejor. De ahí manejo los tiempos y todo.

[08:23] Cuanto tiempo se tarda cada cosa.

[08:26] Voy a hacer así.

[08:29] Y bueno.

[08:36] No sé, debe de ver cuántos niveles de ruido uso también.

[08:43] Ese es un tema.

[08:46] Porque estaba manejando creo que como cuatro niveles de ruido.

[08:50] Que no estaría mal, ¿no?

[08:53] Pero bueno.

[08:55] Claro, niveles de ruido sí.

[08:58] Ahí yo tendría que hacer unos pocos.

[09:02] No sé si hacer con dos.

[09:06] los niveles de ruido alcanza, puede ser, a ver si hay algún efecto, puede ser, dos o tres

[09:16] ponerme, y yo creo que con tres niveles está bien, que sería, bueno, digamos, bien por

[09:34] debajo y bien por arriba de...

[09:36] el ruido de gating o algo así.

[09:48] No, más que el ruido de gating sería la conductancia, ¿no?

[09:53] Sí.

[09:54] La conductancia...

[09:56] Sí, por eso, sí.

[09:58] Sí, está bien.

[10:02] Sigue limando 3 niveles de ruido, 5 o 7 niveles de canales y todos los niveles que pueda de

[10:13] intervalo de medición.

[10:16] Y ya está.

[10:18] Sí, va a ser un volumen de gatos más o menos grande, pero bueno, va a quedar algo muy bonito.

[10:32] Y creo que no hay otro.

[10:35] Ah, y después lo tengo que hacer todo eso mismo en el estacionario.

[10:43] Y ahí en el estacionario lo que puedo hacer es algo,

[10:48] que es ir diez veces tau, ¿no?

[10:55] O sea, es cuando vos estás muy por encima de tau, ¿qué pasa?

[11:01] Si eso se manifiesta, no, yo creo que algo se va a manifestar.

[11:06] Así que bueno.

[11:13] Claro que sería como los Mises de Benz, ¿no?

[11:15] O sea, los Mises de Benz son eventos que son mucho más rápidos que nuestra tasa de emisión.

[11:24] Bueno, por ahí los ves.

[11:26] Digamos, como una especie de fantasma estadístico.

[11:30] Sí, como pequeñas fluctuaciones, no sé.

[11:35] Eso es interesante verlo, a ver si se puede.

[11:43] Está bien, eso ya es el otro punto que sería el experimento figura 3 posiblemente, que

[11:53] sería el estacionario, este es con no estacionario, no me acuerdo cómo se llama, es pulsado,

[12:04] es un pulso, es como sube y baja, acá bueno, ese otro sería estacionario, y después el otro sería todo lo mismo, pero con dos estados, tres estados,

[12:20] o sea, no con dos estados sino con tres estados, a ver qué da, ¿qué sería la figura de cuatro quizás?

[12:34] que se hace estacionario o no estacionario.

[12:37] Y ahí con la FIRO 4K, vos lo que tendrías es si la...

[12:41] tenés otro dato más que es la distancia...

[12:44] bueno, ahí habría dos dimensiones más que ver,

[12:48] que serían la...

[12:51] digamos,

[12:59] el peso de cada estado,

[13:02] cada estado, cuánto contribuye cada uno y la diferencia constante de tiempo, mucho más

[13:20] rápido, etcétera.

[13:24] Creo que eso sería todo.

[13:30] Claro, ahí vos podrías ver también con un experimento en escalera, pero eso ya es

[13:44] Es otro tema, no sé.

[13:46] No sé, ahí digo, a ver.

[13:54] Bien.

[13:56] ¿Qué hago? Pues ahí nos vamos un poco a la mierda.

[14:02] Pero eso yo creo que sería para otro paper, que sería un poco analizar ese tipo de cosas.

[14:11] qué cosas se pueden diferenciar y qué cosas no etcétera

[14:16] y eso es importante

[14:24] eso es verdad que sí especialmente en contra las cosas que no se pueden diferenciar

[14:29] en fin

[14:41] El Rey de la vida.

[14:42] El Rey de la vida.

[14:44] El Rey de la vida.

[14:46] El Rey de la vida.
