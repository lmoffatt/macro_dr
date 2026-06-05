**Archivo:** WhatsApp Ptt 2026-03-27 at 09.53.55.mp3
**Duración:** ~17:04
**Hablante:** Luciano

[00:00] la 9 de julio

[00:05] bueno, ayer estuve trabajando en la figura

[00:11] de la, como llaman

[00:14] la optimización del prior, o sea la optimización, la actualización

[00:21] del prior, o sea la figura de las elipses de prior, posterior

[00:27] que representan de alguna manera el meta estado, el boundary state.

[00:36] Creo que es una figura importante, como estaba diciendo en un audio anterior,

[00:46] lo que demuestra esta figura comparándola en IR versus MR,

[00:53] es que en IR las elipses son más alargadas y en MR son más redondeadas,

[01:01] lo cual indica que es más balanceado el update entre el inicio y el final del intervalo.

[01:14] Pero de alguna manera la idea de la figura es dar una intuición

[01:21] de lo que significa este meta estado.

[01:25] ¿Qué es lo que tiene de malo la figura?

[01:30] Es que la figura justamente vos la podés hacer

[01:32] también con el argolito MR

[01:37] y en el MR no usás este meta estado

[01:40] sin embargo lo podés usar, lo cual puede confundir un poco.

[01:49] Y la verdad es que más allá de eso de que está alargado y redondeado,

[01:53] no hay otra cosa que muestre la diferencia entre un algoritmo y el otro.

[01:57] Entonces uno plantea la pregunta, ¿habrá alguna manera de demostrar el algoritmo mejor?

[02:12] Bueno, en la figura 1, la que muestro la evolución, la actualización por Markov y por Valles del estado del canal,

[02:33] Ahí sí se nota un poco más porque se ve que mientras la actualización por MR es básicamente paralela una línea a la otra, en la actualización por IR puede ser bastante sin paralelo.

[02:51] Es como que se ve ahí claro que la línea se mueve.

[03:00] Quizás, digamos, en el gráfico de las elipses, si yo uno los dos centros de las elipses,

[03:09] eso quizás se pueda ver un poco más, ahora que pienso.

[03:14] Probablemente se vea. Voy a hacer eso.

[03:17] Ese es un buen punto.

[03:19] Bueno, entonces, claro, con las elipses lo que podría hacer es

[03:27] no hacer solamente

[03:29] digamos

[03:31] presentar

[03:38] el estado

[03:38] el meta estado

[03:40] no solamente con el abierto

[03:42] abierto sino bueno hacer con

[03:44] cerrado cerrado, cerrado abierto

[03:46] abierto cerrado

[03:47] no sé

[03:50] cosa como que

[03:51] que dé una

[03:53] una idea un poco más

[03:56] general

[03:57] si eso no estaría mal

[04:04] eso podría funcionar

[04:10] como para mostrar el estado

[04:11] lo que pasa es que en realidad

[04:13] la pregunta es

[04:15] queremos que el lector pierda

[04:20] energía mental

[04:22] en algo que no conduce a nada

[04:24] porque

[04:28] a decir verdad

[04:32] estos

[04:32] digamos yo no

[04:35] las variables que uso

[04:38] para las elipses

[04:39] las uso para las elipses

[04:41] y nada más la correlación cruzada

[04:44] abierto cerrado no la uso

[04:46] luego, o sea me queda ahí nomás

[04:48] o sea para lo que es

[04:50] el algoritmo propiamente dicho

[04:52] esa elipse

[04:54] es un poco irrelevante

[04:56] entonces uno se pregunta

[04:59] ¿vale la pena ponerla?

[05:03] yo creo que la figura 1 sí

[05:07] es claro que sí

[05:08] está muy bien

[05:10] creo que por ahí sí es verdad

[05:16] la elipse podría servir como para mostrar definitivamente

[05:20] claro, lo que podría hacer

[05:24] es no poner las elipses

[05:25] para el caso de

[05:27] el MR, sino

[05:29] ahí poner realmente

[05:31] la actualización

[05:33] por

[05:34] por gaussianas, ¿no?

[05:37] que una gaussiana y la otra gaussiana

[05:39] o sea, la gaussiana en una sola dimensión

[05:41] también es una alternativa, o sea, eso tengo que

[05:47] considerar la alternativa de si hago elipse, elipse

[05:49] o hago Gaussian elipse

[05:51] y bueno, y después está el tema

[05:57] de cómo muestro todo

[05:59] el algoritmo, o sea el algoritmo

[06:01] ahora vamos a pensarlo, lo primero

[06:03] o sea si vamos a hacer miles de figuras

[06:05] la primera figura sería una figura

[06:07] mostrando cómo durante

[06:09] un intervalo de medición

[06:11] los canales se abren

[06:13] y cierran

[06:14] y que entonces vos estás tomando el premedio

[06:17] de un estado de

[06:19] que es el promedio de muchos estados

[06:22] o sea eso para plantear

[06:24] el problema me parece bien

[06:25] después entonces

[06:28] ahí planteo bueno

[06:30] cómo funciona el algoritmo que vos decís

[06:31] bueno

[06:32] claro

[06:35] tendría como primero

[06:46] el prior

[06:46] después del prior

[06:49] la predicted current

[06:51] comparada con la

[06:53] measured current

[06:54] después tengo

[06:56] el este

[06:59] la likelihood

[07:01] y el posterior

[07:03] y bueno

[07:04] y después tengo el next prior

[07:07] sería

[07:07] digamos eso sería uno

[07:16] y el otro sería que tengo la predicted current y tengo justamente al mismo tiempo que el

[07:28] posterior ya es el Markov Transition.

[07:46] Lo que podría hacer también es pulsar la elipse como para mostrar el posterior, o sea

[08:01] como el posterior yo formo la proyección de la elipse en el estado final del intervalo.

[08:15] tomar eso como el input y en el otro caso es tomo claro, podría ser alguna manera de

[08:31] de mostrar, bueno igual eso se muestra en ese gráfico mejor, en el gráfico de la figura

[08:37] 1. Bueno, claro y después que tengo, después lo que tengo que mostrar es el tema de, ahí

[08:54] tengo que remar

[08:57] una likelihood global de todo el corriente

[09:05] y claro digamos como yo sé

[09:09] qué

[09:11] que este método que está la equipo es precisa no sea cómo puedo comparar una

[09:22] la de click y con la otra, estoy simulando el mismo proceso, en un caso digamos, tomo

[09:30] la media inicial del intervalo, otra cosa tomo este boundary state.

[09:37] Bien, entonces lo que hago es, calculo lo que sería el gradiente, la esperanza del

[09:49] gradiente dado que hago muchas repeticiones muchas samples de la simulación para un conjunto

[10:01] de parámetros fijos y yo calculo el gradiente del score o sea el gradiente de la likelihood

[10:10] respecto de los parámetros del modelo y

[10:15] este

[10:18] y veo digamos si este la esperanza del score digamos

[10:26] está dentro de él haciendo un bootstrap hago hay un bootstrap y veo si la

[10:34] la esperanza del score es distinta de 0, está contenida dentro del 0

[10:41] bien

[10:43] y después el otro test que hago es veo si lo que sería

[10:54] las dos

[10:57] O sea, yo tendría dos aproximaciones a la covarianza de los parámetros que estarían dadas.

[11:06] Por un lado por el gesiano y por otro lado por el jacobiano.

[11:12] O sea, por el gesiano yo lo que uso es el gesiano de la distribución normal.

[11:24] Y tomo la esperanza de ese hegessiano y ese lo comparo con, digamos, ese sería el hegessiano teórico,

[11:33] el hegessiano, digamos, sampleado, que estaría dado por la covarianza del score.

[11:53] Entonces, ambos deberían ser iguales y si son diferentes, yo tengo ahí una especie de factor que te dice, bueno, cuánto estás, cuánto se aleja.

[12:18] O sea, a ver, vamos a ver dos cosas.

[12:23] Por lo tanto yo tengo, si el gradiente me da cero, en principio la cosa no sería sesgada,

[12:34] o sea, yo la esperanza de gradiente me da cero, digamos, vos no tendrías sesgo en la obtención de los parámetros.

[12:46] pero lo que podrías tener es un error

[12:50] en cuanto a la estimación de los grados de libertad

[12:54] de alguna manera, o sea, de lo que serían los estadísticos

[12:58] de la varianza de los, o sea, vos tendrías un error en la varianza de los parámetros

[13:02] con lo cual, digamos, todo lo que serían

[13:07] test estadísticos que te dicen que tan

[13:10] cerca estás o no de

[13:13] un valor, o sea, el error de los parámetros

[13:16] sería, digamos

[13:18] sería erróneo

[13:21] lo que sea

[13:22] entonces, claro

[13:25] vos podrías

[13:26] usar

[13:28] justamente

[13:31] la

[13:32] covarianza del gradiente

[13:37] para

[13:39] estimar

[13:40] el verdadero error de los parámetros

[13:43] eso en principio

[13:48] creo que sería válido

[13:51] aún cuando la likelihood

[13:53] sea errónea

[13:54] porque sale por una definición

[14:02] digamos, es la esperanza

[14:03] de los parámetros

[14:07] y entonces vos podrías estimar una especie de, digamos en realidad cuál es el problema

[14:25] de estimar la varianza del gradiente en el caso de una...

[14:42] O sea, a ver, yo tengo por ejemplo una, o sea, sampleo no ahora la, yo sampleo no la variable que estoy midiendo,

[15:07] sino que se amplía los parámetros con un Metropolis Monte Carlo.

[15:12] Ahí yo también puedo estimar la varianza de los parámetros.

[15:22] No, ahí puedo estimar la varianza del gradiente.

[15:31] este es el crediente y me va a dar cero porque justamente es lo que busca el Metrópolis

[15:44] Monte Carlo. Lo que pasa es que yo lo que no tengo ahí es la posibilidad de ver que

[15:56] que no sea sesgado el coso, claro, pero podría después hacer una simulación y ahí ver

[16:02] si el gradiente, digamos si hago muchas simulaciones, claro, yo tengo un punto por ejemplo, digo

[16:16] bueno, pum, este es mi dato, este es mi parámetro, entonces tomo ese parámetro con grader, hago

[16:24] lo simulo muchas veces y ahí veo cómo es la esperanza del gradiente, esperanza del score

[16:46] y ahí me tendría que dar cero, y si me da distinta de cero, bueno yo ahí tendría

[16:53] una medición de qué tan errado está el modelo.
