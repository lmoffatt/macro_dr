**Archivo:** WhatsApp Ptt 2026-03-26 at 11.06.26.mp3
**Duración:** ~12:46
**Hablante:** Luciano

[00:00] Bueno, hace un par de días estoy trabajando full con Macro IR, con el nuevo Paper AnyLife.

[00:07] Estoy ahora caminando a las 9 de julio.

[00:10] ¿Por qué salí a caminar? Porque, bueno, escuchándome, digamos, se me prendió la lamparita.

[00:19] Y, digamos, se me ocurrió una forma de usar Macro IR que sea mucho más efectiva.

[00:28] y la misma consiste en fitear tramos cortos de los traces con un número muy reducido de estados,

[00:44] o sea, entre 2, 3 o hasta 4 estados, quizás con una topología hiperconectada,

[00:57] de manera que no haya ningún tipo de modelo que decidir,

[01:05] eso se me ocurre ahora, o no, o podés hacer con algunas rates que sean cero,

[01:14] o sea que sean muy grandes.

[01:20] Y bueno, nada, me parece como una forma muy eficiente de usar Macro IR,

[01:32] porque uno podría correrlo en distintas partes del trace y eso te daría una respuesta inmediata, básicamente, espero.

[01:44] Y te daría datos interesantes, como el número de canales, la probabilidad de la conductancia,

[01:53] y los R de opening and closing.

[02:00] Nada de eso.

[02:04] Bueno, después otra cosa que pensé es el tema de cómo analizo la eficiencia de macro IR

[02:13] comparado con macro MR, que es el que toma en un solo punto,

[02:25] actualiza no el boundary state sino el point state

[02:30] del canal y bueno y en principio las

[02:39] las variables o sea las meta variables no sé cómo llamarlo relevantes son bueno en principio

[02:51] una es quizás la más importante, es la tau, es decir, la constante de tiempo de la apertura,

[03:00] dividido el delta, que es el intervalo de medición, eso considerando una única tau.

[03:11] Después tenemos el otro importante, es el sigma dividido el épsilon,

[03:18] sigma sería el ruido

[03:20] el gating noise

[03:22] dividido el gaussian noise

[03:25] o sea el instrumental noise

[03:26] después

[03:28] esos son los dos

[03:30] digamos parámetros

[03:31] que claro que van de

[03:34] de cero a infinito

[03:36] y después tendríamos

[03:41] parámetros que son

[03:43] enteros

[03:45] que sería el número de canales

[03:46] y el número de estados

[03:49] o sea, bueno

[03:51] si irían de 1 a infinito

[03:53] digamos ambos

[03:54] y después

[03:56] cosas un poco más interesantes

[03:59] que es de donde me salió todo este tema

[04:00] es

[04:01] lo que sería

[04:03] yo pensé es la amplitud

[04:07] de

[04:08] digamos

[04:10] el rango

[04:12] cinético tomado como la tau

[04:14] máxima dividida a tau mínima

[04:17] o sea, o la diferencia

[04:19] de los logaritmos mejor

[04:20] del logaritmo de la tau máxima

[04:22] menos el logaritmo de la tau mínima

[04:25] o no sé, al revés, qué sé yo

[04:26] eso, y eventualmente uno podría

[04:29] pensar

[04:30] algún tipo de medición de entropía

[04:33] del espectro

[04:35] de

[04:36] espectro cinético

[04:39] que vos tengas

[04:40] para cada constante cinética

[04:43] vos vas a tener una proporción, digamos, un peso,

[04:49] y entonces vos ahí lo que tomás es básicamente una entropía,

[05:00] una función de entropía, pero lo interesante sería hacerla una densidad,

[05:08] tomar con la densidad de, ¿cómo se llama?

[05:12] de cinética

[05:18] con la cual, digamos

[05:20] vos de alguna manera

[05:23] incorporás

[05:25] en el tema

[05:29] en el esquema

[05:32] si las constates cinéticas están alejadas

[05:36] o tan lejos o cerca

[05:37] entonces para eso vos tendrías que tener una medición de la entropía que sea de alguna manera diferencial,

[05:51] o sea que esté sobre, no en el continuo, en el espectro continuo de las longitudes de onda,

[06:03] en este caso sería que serían las constantes cinéticas.

[06:07] Y después, bueno, tendrías las topologías también,

[06:10] sería que esos allá son cosas medio que no sé cómo se harían,

[06:16] que bueno, lo de topologías es lo que me di cuenta de que hacer esto,

[06:22] esto de poner, titear pocos estados está bueno,

[06:28] porque bueno, si yo tomo por ejemplo un estado de cuatro estados,

[06:31] dos cerrados y dos abiertos, lo que podría hacer es justamente ver las constantes cinéticas

[06:40] que estén cerca o lejos, etc. Por ejemplo, uno de tres vos podrías tener dos abiertos

[06:48] conectados con un cerrado y que uno de los abiertos sea un poco más rápido que el otro

[06:53] y la diferencia de rapidez, cómo lo podés detectar, etc.

[06:58] Eso es un poco la idea de ver qué es lo que uno puede ver con este sistema

[07:07] y qué es lo que no podés ver.

[07:09] Eso me parece que estaría bueno hacer ese tipo de experimentos numéricos

[07:17] porque son realmente, te dan una intuición de qué es lo que vos podés distinguir, ¿no?

[07:24] O sea, justamente, como para ver si uno le puede creer a estos experimentos alostéricos

[07:33] tan tan complejos, ¿no?

[07:34] Que es increíble que pueda dar, pero bueno, hay que tener un poco de base, ¿no?

[07:40] O sea, ahí está.

[07:42] El tema es que yo lo que quiero hacer es, con este tipo de análisis, es creíble la diferencia entre esos temas a los técnicos tan complejos y o no, o es un delirio.

[07:57] Entonces está bueno, eso está muy bueno.

[08:02] creo que eso le va a dar fuerza al paper porque es como que vos estás tratando de validar desde

[08:10] ese punto de vista también.

[08:15] Está bien.

[08:19] Bien.

[08:21] ¿Alguna cosa más?

[08:26] no bueno y la idea de que a una vida que era importante que me había quedado en la ducha

[08:33] es que vos no solamente fiteas digamos número de canal la conductancia las constantes cinéticas

[08:38] sino las prioridades iniciales las prioridades iniciales también son un dato que vos estás

[08:45] fiteando eso es muy importante también tenerlo en cuenta o sea ahí habría que justamente

[08:55] claro bueno bueno vos podés partir de un tal como hace eso claro unas

[09:02] probabilidades que sean hasta buena mira vos no se me había

[09:07] ocurrido podría aplicar siempre que vos tengas una especie de distribución

[09:13] totalmente

[09:16] agnóstica del estado del canal y bueno lo me dice una línea de base y ves si está

[09:24] cerrado, porque yo lo estoy tomando como que está cerrado de entrada, pero en realidad

[09:29] de ahí yo estoy metiendo por la ventana un dato que no tengo, mirá, y no es tan difícil

[09:36] plantear un esquema donde, claro, vos no sabes dónde está, ¿no?

[09:43] Claro, bueno, vos lo que tenés es que en realidad el canal está en un estado de equilibrio,

[09:50] pero claro, podría estar en un estado...

[09:53] Claro, no sé, ahí no sé cómo hacer para que cualquier estado del canal sea equiprobable, ¿no?

[10:03] Ahí justamente estamos medio jodidos porque la distribución normal no te deja hacer eso, porque es...

[10:11] Sí, porque no, no te deja.

[10:19] eso es un tema interesante, no sé, eso es un tema aparte, porque ahora yo lo que estoy

[10:25] haciendo es tomando como una distribución inicial hacia la Star, pero sí,

[10:33] bueno, en el caso de pocos canales, ¿cómo hacés para...? Justamente eso es un tema que no sé cómo resolverlo.

[10:45] Ah, porque si vos en realidad tenés pocos estados, podrías hacer un microscópico.

[10:50] Eso es verdad, bien, ojo.

[10:52] Ojo al piojo.

[10:54] No, pero microscópico no, porque son...

[10:59] No, no, estoy confundiendo, el número de canales sea reducido, no el número de estados.

[11:04] Bueno, sí, justamente está el tema de los microscópicos, yo pensé que lo podría tratar

[11:11] microscópicos

[11:12] al ser solamente dos estados

[11:15] y después me di cuenta

[11:17] que no, porque

[11:18] por una cosa que cada tanto vuelvo a descubrir

[11:22] vuelvo a olvidar

[11:23] de que no es reducible

[11:25] porque

[11:27] tenés

[11:28] una cantidad de estados posibles

[11:31] que es así medio

[11:32] fabulósica

[11:34] y entonces no sé por qué

[11:39] no podés, no me acuerdo bien, creo que era, no me acuerdo cómo es el tema, pero había

[11:54] algo que no podías hacer, no sé, eso lo tengo que escribir en algún lado, por qué

[11:59] no puedo porque el micro IR es un algoritmo que el micro IR es medio como que no es analítico

[12:16] el micro IR es analítico

[12:19] pero el micro IR

[12:20] no, tenés que usar

[12:22] MSM

[12:25] algo así

[12:27] una aproximación

[12:28] y no me acuerdo bien por qué mierda

[12:31] pero me acuerdo

[12:33] que sí, que era así

[12:34] yo digo chao

[12:36] puedo hacer micro IR

[12:38] pero no, no me acuerdo por qué era

[12:40] tengo que fijar bien
